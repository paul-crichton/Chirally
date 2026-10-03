import { describe, it, expect } from 'vitest';
import { Mol, TetraSpec, DbSpec } from '../src/chem/mol';
import { parseSmiles } from '../src/chem/smiles';
import { perceiveStereo2D, ccwFromPositions } from '../src/chem/stereo2d';
import {
  embed3D, optimizeGeometry, conformerSearch, generate3D, toXYZ, UFF, coordsOf, stereoReport, addExplicitHydrogens,
  typeUFF, uffEnergy,
} from '../src/chem/3d/index';

type V3 = [number, number, number];

const SMILES: Record<string, string> = {
  ethane: 'CC',
  ethene: 'C=C',
  ethyne: 'C#C',
  benzene: 'c1ccccc1',
  cyclohexane: 'C1CCCCC1',
  water: 'O',
  ammonia: 'N',
  methanol: 'CO',
  aceticAcid: 'CC(=O)O',
  lAlanine: 'N[C@@H](C)C(=O)O',
  rButanol: 'CC[C@@H](C)O',
  sButanol: 'CC[C@H](C)O',
  eButene: 'C/C=C/C',
  zButene: 'C/C=C\\C',
  caffeine: 'Cn1cnc2c1c(=O)n(C)c(=O)n2C',
  ibuprofen: 'CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O',
  cholesterol: 'C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C',
  sulfanilamide: 'Nc1ccc(cc1)S(N)(=O)=O',
  trimethylPhosphate: 'COP(=O)(OC)OC',
};
/** Erythromycin A: 51 heavy atoms, 118 atoms with H, 18 stereocentres, 14-membered macrolactone. */
const ERYTHROMYCIN =
  'CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O';

// ── geometry helpers ──────────────────────────────────────────────────────────
const P = (m: Mol, i: number): V3 => [m.atoms[i].x, m.atoms[i].y, m.atoms[i].z ?? 0];
const sub = (a: V3, b: V3): V3 => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dot = (a: V3, b: V3) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a: V3, b: V3): V3 => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const len = (a: V3) => Math.hypot(a[0], a[1], a[2]);
const dist = (m: Mol, i: number, j: number) => len(sub(P(m, i), P(m, j)));
const angle = (m: Mol, i: number, j: number, k: number) => {
  const u = sub(P(m, i), P(m, j)), v = sub(P(m, k), P(m, j));
  return (Math.acos(Math.max(-1, Math.min(1, dot(u, v) / (len(u) * len(v))))) * 180) / Math.PI;
};
const dihedral = (m: Mol, i: number, j: number, k: number, l: number) => {
  const b0 = sub(P(m, j), P(m, i)), b1 = sub(P(m, k), P(m, j)), b2 = sub(P(m, l), P(m, k));
  const n1 = cross(b0, b1), n2 = cross(b1, b2);
  const y = dot(cross(n1, n2), b1) / len(b1);
  return (Math.atan2(y, dot(n1, n2)) * 180) / Math.PI;
};
/** Largest distance of the given atoms from their least-squares plane (via the smallest-variance normal). */
function planarity(m: Mol, atoms: number[]): number {
  const c = atoms.reduce<V3>((s, i) => [s[0] + P(m, i)[0] / atoms.length, s[1] + P(m, i)[1] / atoms.length, s[2] + P(m, i)[2] / atoms.length], [0, 0, 0]);
  // try normals from all atom pairs' cross products and keep the best (robust for small rings)
  let best = Infinity;
  for (let a = 0; a < atoms.length; a++) {
    for (let b = a + 1; b < atoms.length; b++) {
      const nrm = cross(sub(P(m, atoms[a]), c), sub(P(m, atoms[b]), c));
      const l = len(nrm);
      if (l < 1e-3) continue;
      const dev = Math.max(...atoms.map((i) => Math.abs(dot(sub(P(m, i), c), nrm)) / l));
      best = Math.min(best, dev);
    }
  }
  return best;
}
/** All bond lengths between elements el1/el2 with the given Mol bond order. */
const bondLengths = (m: Mol, el1: string, el2: string, order?: number) =>
  m.bonds
    .filter((b) => (order === undefined || b.order === order) && [m.atoms[b.a].el, m.atoms[b.b].el].sort().join() === [el1, el2].sort().join())
    .map((b) => dist(m, b.a, b.b));
const heavy = (m: Mol) => m.atoms.map((_, i) => i).filter((i) => m.atoms[i].el !== 'H');

/**
 * Independent stereo check of an INPUT spec (heavy-atom indices, -1 = implicit H / lone pair) against the
 * 3D output: the implicit H is the centre's neighbour not listed in the spec.
 */
function inputTetraOK(out: Mol, spec: TetraSpec): boolean {
  const listed = spec.nbrs.filter((v) => v >= 0);
  const extra = out.neighbors(spec.center).filter((v) => !listed.includes(v));
  const c = P(out, spec.center);
  const pos = spec.nbrs.map((v): V3 => {
    if (v >= 0) return P(out, v);
    if (extra.length === 1) return P(out, extra[0]);
    const mean = listed.reduce<V3>((s, i) => [s[0] + P(out, i)[0] / 3, s[1] + P(out, i)[1] / 3, s[2] + P(out, i)[2] / 3], [0, 0, 0]);
    const d = sub(c, mean);
    return [c[0] + d[0], c[1] + d[1], c[2] + d[2]];
  });
  return ccwFromPositions(pos[0], pos[1], pos[2], pos[3]) === spec.ccw;
}
function inputDbOK(out: Mol, spec: DbSpec): boolean {
  const b = out.bonds[spec.bond];
  const [ea, eb] = out.bondBetween(spec.a, b.a) >= 0 ? [b.a, b.b] : [b.b, b.a];
  return Math.abs(dihedral(out, spec.a, ea, eb, spec.b)) < 90 === spec.cis;
}

function build(name: string, seed = 1): { input: Mol; out: Mol; result: ReturnType<typeof optimizeGeometry> } {
  const input = parseSmiles(SMILES[name]);
  const out = embed3D(input, { seed });
  const result = optimizeGeometry(out);
  return { input, out, result };
}

// ── tests ─────────────────────────────────────────────────────────────────────
describe('explicit hydrogens', () => {
  it('appends H atoms and re-points implicit-H stereo entries', () => {
    const input = parseSmiles(SMILES.lAlanine);
    const m = addExplicitHydrogens(input);
    expect(m.atoms.length).toBe(13);
    expect(m.atoms.slice(0, input.atoms.length).map((a) => a.el)).toEqual(input.atoms.map((a) => a.el));
    const t = m.tetra[0];
    expect(t.nbrs).not.toContain(-1);
    const hIdx = t.nbrs[input.tetra[0].nbrs.indexOf(-1)];
    expect(m.atoms[hIdx].el).toBe('H');
    expect(m.bondBetween(hIdx, t.center)).toBeGreaterThanOrEqual(0);
  });
  it('keeps a lone-pair entry (-1) on centres without hydrogen', () => {
    const m = addExplicitHydrogens(parseSmiles('C[S@](=O)c1ccccc1'));
    expect(m.tetra[0].nbrs).toContain(-1);
  });
});

describe('UFF typing', () => {
  it('assigns hybridisation-specific types', () => {
    const lab = (smi: string) => typeUFF(addExplicitHydrogens(parseSmiles(smi))).labels;
    expect(lab('CC=CC#N').slice(0, 5)).toEqual(['C_3', 'C_2', 'C_2', 'C_1', 'N_1']);
    expect(lab('c1ccncc1').slice(0, 6)).toEqual(['C_R', 'C_R', 'C_R', 'N_R', 'C_R', 'C_R']);
    expect(lab('CC(=O)NC').slice(0, 5)).toEqual(['C_3', 'C_2', 'O_2', 'N_R', 'C_3']);
    expect(lab('CS(=O)(=O)N')[1]).toBe('S_3+6');
    expect(lab('COP(=O)(OC)OC')[2]).toBe('P_3+5');
    expect(lab('CSC')[1]).toBe('S_3+2');
    expect(lab('O')[0]).toBe('O_3');
    expect(lab('ClCBr').filter((l) => l !== 'H_')).toEqual(['Cl', 'C_3', 'Br']);
  });
});

describe('UFF energy and gradient', () => {
  it('analytic gradient matches finite differences', () => {
    for (const smi of [SMILES.caffeine, SMILES.sulfanilamide, SMILES.trimethylPhosphate, 'CC(=O)Nc1ccc(O)cc1', 'C#CC=CC1CC1', 'CP(C)C', 'OO']) {
      const m = embed3D(parseSmiles(smi), { seed: 3 });
      for (const electrostatics of [false, true]) {
        const ff = new UFF(m, { electrostatics });
        const x = coordsOf(m);
        for (let i = 0; i < x.length; i++) x[i] += 0.05 * Math.sin(1.7 * i); // off-equilibrium
        const g = new Float64Array(x.length);
        ff.energy(x, g);
        const h = 1e-5;
        for (let i = 0; i < x.length; i++) {
          const xp = x.slice(), xm = x.slice();
          xp[i] += h;
          xm[i] -= h;
          const num = (ff.energy(xp) - ff.energy(xm)) / (2 * h);
          expect(Math.abs(num - g[i])).toBeLessThan(1e-4 * Math.max(1, Math.abs(num)));
        }
      }
    }
  });

  it('optimisation decreases the energy monotonically and converges', () => {
    const m = embed3D(parseSmiles(SMILES.ibuprofen), { seed: 2 });
    const e0 = uffEnergy(m).energy;
    const trace: number[] = [];
    const r = optimizeGeometry(m, { onProgress: (_, e) => trace.push(e) });
    expect(r.converged).toBe(true);
    expect(r.rmsGrad).toBeLessThan(1e-3);
    expect(r.energy).toBeLessThan(e0);
    for (let k = 1; k < trace.length; k++) expect(trace[k]).toBeLessThanOrEqual(trace[k - 1] + 1e-9);
    const sum = Object.values(r.terms).reduce((s, v) => s + v, 0);
    expect(sum).toBeCloseTo(r.energy, 8);
    expect(r.terms.elec).toBe(0);
  });

  it('optional electrostatics contributes a Coulomb term', () => {
    const m = embed3D(parseSmiles(SMILES.aceticAcid));
    const r = optimizeGeometry(m, { electrostatics: true });
    expect(r.terms.elec).not.toBe(0);
    expect(r.converged).toBe(true);
  });
});

describe('optimised geometries', () => {
  it('ethane: C–C, C–H, tetrahedral angles, staggered', () => {
    const { out, result } = build('ethane');
    expect(result.converged).toBe(true);
    expect(bondLengths(out, 'C', 'C')[0]).toBeCloseTo(1.53, 1);
    for (const d of bondLengths(out, 'C', 'H')) expect(Math.abs(d - 1.09)).toBeLessThan(0.03);
    const angs: number[] = [];
    for (const c of [0, 1]) {
      const nb = out.neighbors(c);
      for (let a = 0; a < nb.length; a++) for (let b = a + 1; b < nb.length; b++) angs.push(angle(out, nb[a], c, nb[b]));
    }
    for (const a of angs) expect(Math.abs(a - 109.5)).toBeLessThan(3);
    expect(Math.abs(angs.reduce((s, a) => s + a, 0) / angs.length - 109.5)).toBeLessThan(1);
    const h0 = out.neighbors(0).find((v) => out.atoms[v].el === 'H')!;
    const h1 = out.neighbors(1).find((v) => out.atoms[v].el === 'H')!;
    const phi = Math.abs(dihedral(out, h0, 0, 1, h1)) % 120;
    expect(Math.abs(phi - 60)).toBeLessThan(5);
  });

  it('ethene: C=C and planarity', () => {
    const { out } = build('ethene');
    expect(Math.abs(bondLengths(out, 'C', 'C', 2)[0] - 1.33)).toBeLessThan(0.03);
    expect(planarity(out, out.atoms.map((_, i) => i))).toBeLessThan(0.01);
  });

  it('ethyne: C≡C and linearity', () => {
    const { out } = build('ethyne');
    expect(Math.abs(bondLengths(out, 'C', 'C', 3)[0] - 1.2)).toBeLessThan(0.03);
    for (const c of [0, 1]) {
      const [a, b] = out.neighbors(c);
      expect(angle(out, a, c, b)).toBeGreaterThan(178);
    }
  });

  it('benzene: aromatic C–C ≈ 1.39 Å, planar, 120°', () => {
    const { out } = build('benzene');
    for (const d of bondLengths(out, 'C', 'C')) expect(Math.abs(d - 1.39)).toBeLessThan(0.03);
    expect(planarity(out, out.atoms.map((_, i) => i))).toBeLessThan(0.01);
    for (let i = 0; i < 6; i++) expect(Math.abs(angle(out, (i + 5) % 6, i, (i + 1) % 6) - 120)).toBeLessThan(1);
  });

  it('cyclohexane: chair after conformer search', () => {
    const m = embed3D(parseSmiles(SMILES.cyclohexane), { seed: 2 });
    const { best, energies } = conformerSearch(m, { n: 10, seed: 1 });
    expect(energies.length).toBeGreaterThan(1);
    expect(energies[0]).toBeLessThanOrEqual(energies[energies.length - 1]);
    const tors = [0, 1, 2, 3, 4, 5].map((k) => dihedral(best, k, (k + 1) % 6, (k + 2) % 6, (k + 3) % 6));
    for (let k = 0; k < 6; k++) {
      expect(Math.abs(tors[k])).toBeGreaterThan(45);
      expect(Math.abs(tors[k])).toBeLessThan(65);
      expect(Math.sign(tors[k])).toBe(-Math.sign(tors[(k + 1) % 6])); // alternating = chair
    }
    for (const d of bondLengths(best, 'C', 'C')) expect(Math.abs(d - 1.53)).toBeLessThan(0.03);
  });

  it('water and ammonia angles', () => {
    const w = build('water').out;
    expect(Math.abs(angle(w, 1, 0, 2) - 104.5)).toBeLessThan(2);
    for (const d of bondLengths(w, 'O', 'H')) expect(Math.abs(d - 0.97)).toBeLessThan(0.04);
    const a = build('ammonia').out;
    for (const [i, j] of [[1, 2], [1, 3], [2, 3]]) expect(Math.abs(angle(a, i, 0, j) - 106.7)).toBeLessThan(3);
    expect(planarity(a, [1, 2, 3, 0])).toBeGreaterThan(0.2); // pyramidal
  });

  it('methanol and acetic acid', () => {
    const me = build('methanol').out;
    expect(Math.abs(bondLengths(me, 'C', 'O')[0] - 1.42)).toBeLessThan(0.04);
    const ac = build('aceticAcid').out;
    expect(Math.abs(bondLengths(ac, 'C', 'O', 2)[0] - 1.21)).toBeLessThan(0.03);
    expect(Math.abs(bondLengths(ac, 'C', 'O', 1)[0] - 1.36)).toBeLessThan(0.05);
    expect(Math.abs(bondLengths(ac, 'C', 'C')[0] - 1.5)).toBeLessThan(0.04);
    expect(Math.abs(angle(ac, 2, 1, 3) - 122)).toBeLessThan(5);
    expect(planarity(ac, [0, 1, 2, 3])).toBeLessThan(0.02); // carboxyl group planar
  });

  it('caffeine: planar ring system and carbonyls', () => {
    const { out, result } = build('caffeine');
    expect(result.converged).toBe(true);
    const ringAtoms = heavy(out).filter((i) => out.neighbors(i).filter((j) => out.atoms[j].el !== 'H').length >= 2 && out.atoms[i].el !== 'O');
    const ring = ringAtoms.filter((i) => !(out.atoms[i].el === 'C' && out.neighbors(i).filter((j) => out.atoms[j].el === 'H').length === 3));
    expect(planarity(out, ring)).toBeLessThan(0.05);
    for (const d of bondLengths(out, 'C', 'O', 2)) expect(Math.abs(d - 1.22)).toBeLessThan(0.03);
  });

  it('sulfonamide and phosphate ester: tetrahedral S and P', () => {
    const s = build('sulfanilamide').out;
    const S = s.atoms.findIndex((a) => a.el === 'S');
    const sn = s.neighbors(S);
    for (let a = 0; a < 4; a++) for (let b = a + 1; b < 4; b++) expect(Math.abs(angle(s, sn[a], S, sn[b]) - 109.5)).toBeLessThan(10);
    // UFF's S=O rest length (≈1.50 Å) is known to be longer than experiment (1.43 Å)
    for (const d of bondLengths(s, 'S', 'O', 2)) expect(Math.abs(d - 1.48)).toBeLessThan(0.06);
    const ringC = heavy(s).filter((i) => s.atoms[i].el === 'C');
    expect(planarity(s, ringC)).toBeLessThan(0.02);
    const p = build('trimethylPhosphate').out;
    const Pi = p.atoms.findIndex((a) => a.el === 'P');
    const pn = p.neighbors(Pi);
    for (let a = 0; a < 4; a++) for (let b = a + 1; b < 4; b++) expect(Math.abs(angle(p, pn[a], Pi, pn[b]) - 109.5)).toBeLessThan(10);
    expect(Math.max(...bondLengths(p, 'P', 'O', 2))).toBeLessThan(Math.min(...bondLengths(p, 'P', 'O', 1)));
  });

  it('(E)- and (Z)-2-butene keep their configuration', () => {
    const e = build('eButene').out;
    const z = build('zButene').out;
    expect(Math.abs(dihedral(e, 0, 1, 2, 3))).toBeGreaterThan(170);
    expect(Math.abs(dihedral(z, 0, 1, 2, 3))).toBeLessThan(10);
    expect(Math.abs(bondLengths(e, 'C', 'C', 2)[0] - 1.33)).toBeLessThan(0.03);
  });
});

describe('stereo preservation', () => {
  const stereoMols = ['lAlanine', 'rButanol', 'sButanol', 'eButene', 'zButene', 'ibuprofen', 'cholesterol'];
  it('every tetrahedral and double-bond spec survives embedding + optimisation', () => {
    for (const name of stereoMols) {
      for (const seed of [1, 2, 3]) {
        const { input, out } = build(name, seed);
        expect(input.tetra.length + input.dbStereo.length, name).toBeGreaterThan(0);
        for (const t of input.tetra) expect(inputTetraOK(out, t), `${name} seed ${seed} centre ${t.center}`).toBe(true);
        for (const d of input.dbStereo) expect(inputDbOK(out, d), `${name} seed ${seed} bond ${d.bond}`).toBe(true);
        expect(stereoReport(out).allOK).toBe(true);
      }
    }
  });

  it('cholesterol keeps all 8 stereocentres', () => {
    const { input, out } = build('cholesterol', 5);
    expect(input.tetra.length).toBe(8);
    expect(input.tetra.every((t) => inputTetraOK(out, t))).toBe(true);
  });

  it('(R)- and (S)-2-butanol come out as mirror images', () => {
    const r = build('rButanol').out, s = build('sButanol').out;
    const vol = (m: Mol) => {
      const c = 2, nb = [1, 3, 4]; // CH2, CH3, O
      const v = nb.map((i) => sub(P(m, i), P(m, c))) as [V3, V3, V3];
      return dot(v[0], cross(v[1], v[2]));
    };
    expect(Math.sign(vol(r))).toBe(-Math.sign(vol(s)));
  });

  it('honours specs derived from a 2D wedge drawing (perceiveStereo2D)', () => {
    // L-alanine: centre (0,0); N lower-left, CH3 up (wedge), COOH lower-right; y is DOWN.
    const m = new Mol();
    const c = m.addAtom({ el: 'C', x: 0, y: 0 });
    const n = m.addAtom({ el: 'N', x: -0.87, y: 0.5 });
    const me = m.addAtom({ el: 'C', x: 0, y: -1 });
    const cc = m.addAtom({ el: 'C', x: 0.87, y: 0.5 });
    m.addAtom({ el: 'O', x: 1.74, y: 0 });
    m.addAtom({ el: 'O', x: 0.87, y: 1.5 });
    m.addBond(c, n);
    m.addBond(c, me, 1, 'wedge');
    m.addBond(c, cc);
    m.addBond(cc, 4, 2);
    m.addBond(cc, 5);
    perceiveStereo2D(m);
    expect(m.tetra.length).toBe(1);
    const out = embed3D(m);
    optimizeGeometry(out);
    expect(inputTetraOK(out, m.tetra[0])).toBe(true);
    // same absolute configuration as the SMILES route
    const ref = build('lAlanine').out; // atoms: N0 C1 C2(Me) C3(OOH)
    const vol = (mm: Mol, centre: number, a: number, b: number, d: number) => {
      const v = [a, b, d].map((i) => sub(P(mm, i), P(mm, centre))) as [V3, V3, V3];
      return dot(v[0], cross(v[1], v[2]));
    };
    expect(Math.sign(vol(out, c, n, me, cc))).toBe(Math.sign(vol(ref, 1, 0, 2, 3)));
  });

  it('keeps sulfoxide (lone-pair) stereo', () => {
    const input = parseSmiles('C[S@](=O)c1ccccc1');
    const out = embed3D(input, { seed: 4 });
    optimizeGeometry(out);
    expect(inputTetraOK(out, input.tetra[0])).toBe(true);
  });
});

describe('determinism, components, export and performance', () => {
  it('is deterministic for a given seed', () => {
    const a = embed3D(parseSmiles(SMILES.ibuprofen), { seed: 42 });
    const b = embed3D(parseSmiles(SMILES.ibuprofen), { seed: 42 });
    const c = embed3D(parseSmiles(SMILES.ibuprofen), { seed: 43 });
    expect(coordsOf(a)).toEqual(coordsOf(b));
    expect(coordsOf(a)).not.toEqual(coordsOf(c));
    const ra = optimizeGeometry(a), rb = optimizeGeometry(b);
    expect(ra.energy).toBe(rb.energy);
    const g1 = generate3D(parseSmiles(SMILES.cholesterol), { seed: 7 });
    const g2 = generate3D(parseSmiles(SMILES.cholesterol), { seed: 7 });
    expect(coordsOf(g1.mol)).toEqual(coordsOf(g2.mol));
  });

  it('places disconnected components apart', () => {
    const out = embed3D(parseSmiles('C[N+](C)(C)C.[Cl-]'));
    optimizeGeometry(out);
    const cl = out.atoms.findIndex((a) => a.el === 'Cl');
    for (let i = 0; i < out.atoms.length; i++) if (i !== cl) expect(dist(out, i, cl)).toBeGreaterThan(2);
  });

  it('generate3D returns an optimised, stereo-correct structure', () => {
    const { mol, result } = generate3D(parseSmiles(SMILES.cyclohexane), { conformers: 6 });
    expect(result.converged).toBe(true);
    const tors = [0, 1, 2, 3, 4, 5].map((k) => dihedral(mol, k, (k + 1) % 6, (k + 2) % 6, (k + 3) % 6));
    expect(tors.every((t, k) => Math.sign(t) === -Math.sign(tors[(k + 1) % 6]))).toBe(true);
  });

  it('writes XYZ', () => {
    const { out } = build('water');
    const xyz = toXYZ(out, 'water\nmolecule');
    const lines = xyz.trim().split('\n');
    expect(lines[0]).toBe('3');
    expect(lines[1]).toBe('water molecule');
    expect(lines.length).toBe(5);
    expect(lines[2]).toMatch(/^O\s+-?\d+\.\d{6}\s+-?\d+\.\d{6}\s+-?\d+\.\d{6}$/);
  });

  it('optimises a 50-heavy-atom molecule (≈120 atoms) in under a second, keeping 18 stereocentres', () => {
    const input = parseSmiles(ERYTHROMYCIN);
    expect(input.atoms.length).toBe(51);
    expect(input.tetra.length).toBe(18);
    const m = embed3D(input, { seed: 1 });
    const t0 = performance.now();
    const r = optimizeGeometry(m);
    const dt = performance.now() - t0;
    expect(dt).toBeLessThan(1000);
    expect(r.converged).toBe(true);
    for (const t of input.tetra) expect(inputTetraOK(m, t)).toBe(true);
  });
});
