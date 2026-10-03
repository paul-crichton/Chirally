// Stereochemistry of 3D coordinates: verification against TetraSpec / DbSpec and local repairs.
import { Mol, TetraSpec, DbSpec } from '../mol';
import { ccwFromPositions } from '../stereo2d';

export type V3 = [number, number, number];

const sub = (a: V3, b: V3): V3 => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const add = (a: V3, b: V3): V3 => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const scale = (a: V3, s: number): V3 => [a[0] * s, a[1] * s, a[2] * s];
const dot = (a: V3, b: V3) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a: V3, b: V3): V3 => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const norm = (a: V3) => Math.hypot(a[0], a[1], a[2]);
const unit = (a: V3): V3 => {
  const l = norm(a) || 1;
  return [a[0] / l, a[1] / l, a[2] / l];
};

export const pos3 = (mol: Mol, i: number): V3 => [mol.atoms[i].x, mol.atoms[i].y, mol.atoms[i].z ?? 0];

/** Dihedral angle a–b–c–d in degrees (−180…180). */
export function dihedral(p0: V3, p1: V3, p2: V3, p3: V3): number {
  const b0 = sub(p0, p1), b1 = sub(p2, p1), b2 = sub(p3, p2);
  const b1u = unit(b1);
  const v = sub(b0, scale(b1u, dot(b0, b1u)));
  const w = sub(b2, scale(b1u, dot(b2, b1u)));
  const xx = dot(v, w);
  const yy = dot(cross(b1u, v), w);
  return (Math.atan2(yy, xx) * 180) / Math.PI;
}

/** Bond angle a–b–c in degrees. */
export function angle(a: V3, b: V3, c: V3): number {
  const u = sub(a, b), v = sub(c, b);
  const cs = dot(u, v) / ((norm(u) || 1) * (norm(v) || 1));
  return (Math.acos(Math.max(-1, Math.min(1, cs))) * 180) / Math.PI;
}

/**
 * Positions of a spec's four neighbours. A -1 entry (lone pair / implicit H) becomes a virtual point
 * opposite the other neighbours.
 */
function specPositions(mol: Mol, spec: TetraSpec): V3[] | null {
  const c = pos3(mol, spec.center);
  const real = spec.nbrs.filter((v) => v >= 0).map((v) => pos3(mol, v));
  if (real.length < 3) return null;
  let virt: V3 | null = null;
  if (real.length === 3) {
    const m = scale(add(add(real[0], real[1]), real[2]), 1 / 3);
    virt = add(c, unit(sub(c, m)));
  }
  return spec.nbrs.map((v) => (v >= 0 ? pos3(mol, v) : virt!));
}

/** Signed chirality volume of a spec in the current 3D coordinates (< 0 ⇔ counter-clockwise ⇔ '@'). */
export function tetraVolume(mol: Mol, spec: TetraSpec): number {
  const p = specPositions(mol, spec);
  if (!p) return 0;
  const a = sub(p[1], p[0]), b = sub(p[2], p[0]), c = sub(p[3], p[0]);
  return dot(a, cross(b, c));
}

/** True if the 3D coordinates realise the spec (coordinates must be right-handed, y up). */
export function tetraSatisfied(mol: Mol, spec: TetraSpec): boolean {
  const p = specPositions(mol, spec);
  if (!p) return false;
  if (Math.abs(tetraVolume(mol, spec)) < 1e-3) return false;
  return ccwFromPositions(p[0], p[1], p[2], p[3]) === spec.ccw;
}

/**
 * True if a four-coordinate centre lies inside the tetrahedron of its neighbours (a geometry that a force
 * field relaxes without inverting). Centres with a lone pair / fewer neighbours always pass.
 */
export function centreInside(mol: Mol, spec: TetraSpec): boolean {
  if (spec.nbrs.some((v) => v < 0)) return true;
  const p = spec.nbrs.map((v) => pos3(mol, v));
  const c = pos3(mol, spec.center);
  const vol = (q: V3[]) => dot(sub(q[1], q[0]), cross(sub(q[2], q[0]), sub(q[3], q[0])));
  const V = vol(p);
  for (let k = 0; k < 4; k++) {
    const q = p.slice();
    q[k] = c;
    if (vol(q) * Math.sign(V) <= 0) return false;
  }
  return true;
}

/** True if the double-bond spec is realised (|dihedral a–A=B–b| < 90° ⇔ cis). */
export function dbSatisfied(mol: Mol, spec: DbSpec): boolean {
  const b = mol.bonds[spec.bond];
  // spec.a is attached to b.a, spec.b to b.b (tolerate swapped attachment)
  let ea = b.a, eb = b.b;
  if (mol.bondBetween(spec.a, ea) < 0) [ea, eb] = [eb, ea];
  const d = dihedral(pos3(mol, spec.a), pos3(mol, ea), pos3(mol, eb), pos3(mol, spec.b));
  return Math.abs(d) < 90 === spec.cis;
}

export interface StereoReport {
  tetra: { spec: TetraSpec; ok: boolean }[];
  db: { spec: DbSpec; ok: boolean }[];
  allOK: boolean;
}

/** Checks every tetrahedral and double-bond spec of an H-explicit 3D molecule. */
export function stereoReport(mol: Mol): StereoReport {
  const tetra = mol.tetra.map((spec) => ({ spec, ok: tetraSatisfied(mol, spec) }));
  const db = mol.dbStereo.map((spec) => ({ spec, ok: dbSatisfied(mol, spec) }));
  return { tetra, db, allOK: tetra.every((t) => t.ok) && db.every((d) => d.ok) };
}

/** Atoms reachable from `start` without passing through `block` (null if this returns to `block`, i.e. a ring). */
export function branchAtoms(mol: Mol, block: number, start: number, limit = Infinity): number[] | null {
  const seen = new Set<number>([block, start]);
  const out = [start];
  for (let q = 0; q < out.length; q++) {
    for (const w of mol.neighbors(out[q])) {
      if (w === block) {
        if (out[q] !== start) return null;
        continue;
      }
      if (!seen.has(w)) {
        seen.add(w);
        out.push(w);
        if (out.length > limit) return null;
      }
    }
  }
  return out;
}

/** Rotation matrix (row-major 3×3) turning unit vector u onto unit vector v. */
function rotationBetween(u: V3, v: V3): number[] {
  const c = dot(u, v);
  const ax = cross(u, v);
  const s = norm(ax);
  if (s < 1e-9) {
    if (c > 0) return [1, 0, 0, 0, 1, 0, 0, 0, 1];
    // 180°: rotate about any axis perpendicular to u
    const p = unit(Math.abs(u[0]) < 0.9 ? cross(u, [1, 0, 0]) : cross(u, [0, 1, 0]));
    return axisAngle(p, Math.PI);
  }
  return axisAngle(scale(ax, 1 / s), Math.atan2(s, c));
}

/** Rodrigues rotation matrix about a unit axis. */
export function axisAngle(k: V3, th: number): number[] {
  const c = Math.cos(th), s = Math.sin(th), t = 1 - c;
  const [x, y, z] = k;
  return [t * x * x + c, t * x * y - s * z, t * x * z + s * y, t * x * y + s * z, t * y * y + c, t * y * z - s * x, t * x * z - s * y, t * y * z + s * x, t * z * z + c];
}

const mulMV = (m: number[], v: V3): V3 => [m[0] * v[0] + m[1] * v[1] + m[2] * v[2], m[3] * v[0] + m[4] * v[1] + m[5] * v[2], m[6] * v[0] + m[7] * v[1] + m[8] * v[2]];

const setPos = (mol: Mol, i: number, p: V3) => {
  mol.atoms[i].x = p[0];
  mol.atoms[i].y = p[1];
  mol.atoms[i].z = p[2];
};

/**
 * Inverts a wrong stereocentre by an "umbrella flip": the centre is reflected through the plane of three
 * neighbours and the fourth neighbour's (acyclic) branch is moved rigidly to the opposite side.
 * Returns false if no suitable branch exists (all four bonds in rings).
 */
export function invertCenter(mol: Mol, spec: TetraSpec): boolean {
  const ci = spec.center;
  const c = pos3(mol, ci);
  const nb = mol.neighbors(ci);
  // Candidate movers: lone pair (no atom), else the smallest acyclic branch (H preferred).
  let mover = -1;
  let branch: number[] | null = null;
  if (nb.length === 3) mover = -2; // lone pair
  else {
    let best = Infinity;
    for (const j of nb) {
      const br = branchAtoms(mol, ci, j, 60);
      if (br && br.length < best) {
        best = br.length;
        mover = j;
        branch = br;
      }
    }
  }
  if (mover === -1) return false;
  const others = nb.filter((j) => j !== mover);
  if (others.length !== 3) return false;
  const [pa, pb, pc] = others.map((j) => pos3(mol, j));
  const nrm = unit(cross(sub(pb, pa), sub(pc, pa)));
  const c2 = sub(c, scale(nrm, 2 * dot(sub(c, pa), nrm)));
  const centroid = scale(add(add(pa, pb), pc), 1 / 3);
  setPos(mol, ci, c2);
  if (branch) {
    const u = unit(sub(pos3(mol, mover), c));
    let v = sub(c2, centroid);
    if (norm(v) < 1e-6) v = scale(nrm, dot(sub(c2, c), nrm) >= 0 ? 1 : -1);
    const R = rotationBetween(u, unit(v));
    for (const a of branch) setPos(mol, a, add(c2, mulMV(R, sub(pos3(mol, a), c))));
  }
  return true;
}

/**
 * Fixes a wrong double-bond configuration by rotating the substituents of one end by 180° about the bond.
 * Returns false when both ends are part of rings.
 */
export function flipDoubleBond(mol: Mol, spec: DbSpec): boolean {
  const b = mol.bonds[spec.bond];
  for (const [end, other] of [[b.b, b.a], [b.a, b.b]]) {
    const moving: number[] = [];
    let ok = true;
    for (const s of mol.neighbors(end)) {
      if (s === other) continue;
      const br = branchAtoms(mol, end, s);
      if (!br) {
        ok = false;
        break;
      }
      moving.push(...br);
    }
    if (!ok) continue;
    const pe = pos3(mol, end), po = pos3(mol, other);
    const R = axisAngle(unit(sub(pe, po)), Math.PI);
    for (const a of moving) setPos(mol, a, add(pe, mulMV(R, sub(pos3(mol, a), pe))));
    return true;
  }
  return false;
}
