// Quick conformer search: alternates fresh distance-geometry embeddings with random torsion
// perturbations of the best structure found so far; every candidate is UFF-minimised and must keep
// all stereo specs.
import { Mol } from '../mol';
import { embedCoordinates } from './embed';
import { UFF, coordsOf, setCoords } from './uff';
import { lbfgs } from './minimize';
import { makeRng, Rng } from './rng';
import { stereoReport, branchAtoms, axisAngle, pos3 } from './stereo3d';

export interface ConformerOptions {
  /** Number of candidate conformers to generate (default 10, including the input geometry). */
  n?: number;
  seed?: number;
  /** Wall-clock budget (default 5000 ms); results are deterministic as long as the budget is not hit. */
  maxTimeMs?: number;
}

/** Rotatable bonds: acyclic single bonds between two non-terminal atoms, excluding amide C–N. */
function rotatableBonds(mol: Mol, ff: UFF): { bond: number; side: number[] }[] {
  const out: { bond: number; side: number[] }[] = [];
  const { bondOrder, rings } = ff.typing;
  mol.bonds.forEach((b, bi) => {
    if (bondOrder[bi] !== 1 || rings.bondInRing[bi]) return;
    if (mol.heavyDegree(b.a) < 2 || mol.heavyDegree(b.b) < 2) return;
    const sa = branchAtoms(mol, b.a, b.b);
    const sb = branchAtoms(mol, b.b, b.a);
    if (!sa || !sb) return;
    out.push({ bond: bi, side: sa.length <= sb.length ? sa : sb });
  });
  return out;
}

function twist(mol: Mol, bond: number, side: number[], angle: number): void {
  const b = mol.bonds[bond];
  // `side` was collected from the atom that is NOT the pivot; find which end it starts at
  const start = side[0];
  const pivot = start === b.b ? b.a : b.b;
  const p0 = pos3(mol, pivot), p1 = pos3(mol, start);
  const ax = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
  const l = Math.hypot(ax[0], ax[1], ax[2]) || 1;
  const R = axisAngle([ax[0] / l, ax[1] / l, ax[2] / l], angle);
  for (const a of side) {
    const at = mol.atoms[a];
    const v = [at.x - p1[0], at.y - p1[1], (at.z ?? 0) - p1[2]];
    at.x = p1[0] + R[0] * v[0] + R[1] * v[1] + R[2] * v[2];
    at.y = p1[1] + R[3] * v[0] + R[4] * v[1] + R[5] * v[2];
    at.z = p1[2] + R[6] * v[0] + R[7] * v[1] + R[8] * v[2];
  }
}

/** True if two atoms more than two bonds apart are closer than 1 Å (a move that would tangle the molecule). */
function hasClash(mol: Mol, nbrs: number[][]): boolean {
  const n = mol.atoms.length;
  const stamp = new Int32Array(n).fill(-1);
  for (let i = 0; i < n; i++) {
    stamp[i] = i;
    for (const j of nbrs[i]) {
      stamp[j] = i;
      for (const k of nbrs[j]) stamp[k] = i;
    }
    const a = mol.atoms[i];
    for (let j = i + 1; j < n; j++) {
      if (stamp[j] === i) continue;
      const b = mol.atoms[j];
      const dx = a.x - b.x, dy = a.y - b.y, dz = (a.z ?? 0) - (b.z ?? 0);
      if (dx * dx + dy * dy + dz * dz < 1) return true;
    }
  }
  return false;
}

function minimise(mol: Mol, maxIter = 2000, gradTol = 1e-3): number {
  const ff = new UFF(mol);
  const x = coordsOf(mol);
  const r = lbfgs((xx, g) => ff.energy(xx, g), x, { maxIter, gradTol });
  setCoords(mol, x);
  return r.value;
}

/**
 * Searches for low-energy conformers of an H-explicit 3D molecule (e.g. from embed3D).
 * The input is not modified. `energies` lists the UFF energies (kcal/mol) of all accepted candidates,
 * sorted ascending; `best` is a copy of the molecule in the lowest-energy conformation.
 */
export function conformerSearch(mol3d: Mol, opts: ConformerOptions = {}): { best: Mol; energies: number[] } {
  const total = Math.max(1, opts.n ?? 10);
  const seed = opts.seed ?? 1;
  const deadline = Date.now() + (opts.maxTimeMs ?? 5000);
  const rng: Rng = makeRng(seed * 31 + 5);

  let best = mol3d.clone();
  let bestE = minimise(best);
  const energies = [bestE];
  const ff0 = new UFF(best);
  const rot = rotatableBonds(best, ff0);
  const nbrs = ff0.typing.nbrs;

  for (let k = 1; k < total && Date.now() < deadline; k++) {
    let cand = best.clone();
    if (k % 2 === 1 || rot.length === 0) {
      // fresh embedding (samples ring conformations and global shape)
      embedCoordinates(cand, { seed: seed + 101 * k, maxAttempts: 3 });
    } else {
      // perturb 1–3 rotatable bonds of the current best; reject moves that pass atoms through each other
      let ok = false;
      for (let tries = 0; tries < 10 && !ok; tries++) {
        cand = best.clone();
        const m = 1 + Math.floor(rng() * Math.min(3, rot.length));
        for (let t = 0; t < m; t++) {
          const r = rot[Math.floor(rng() * rot.length)];
          twist(cand, r.bond, r.side, (Math.PI / 3) * (1 + Math.floor(rng() * 5)));
        }
        ok = !hasClash(cand, nbrs);
      }
      if (!ok) continue;
    }
    const e = minimise(cand, 1500, 5e-3); // screening-level minimisation
    if (!Number.isFinite(e) || !stereoReport(cand).allOK) continue;
    energies.push(e);
    if (e < bestE - 1e-6) {
      bestE = e;
      best = cand;
    }
  }
  // polish the winner to the default convergence threshold
  bestE = Math.min(bestE, minimise(best));
  energies.sort((a, b) => a - b);
  energies[0] = Math.min(energies[0], bestE);
  return { best, energies };
}
