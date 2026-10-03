// Gasteiger–Marsili (PEOE) partial charges, used by the optional electrostatic term.
import { Mol } from '../mol';

// Orbital electronegativity parameters χ = a + b·q + c·q²  (Gasteiger & Marsili, Tetrahedron 1980)
const PARAMS: Record<string, [number, number, number]> = {
  H: [7.17, 6.24, -0.56],
  C3: [7.98, 9.18, 1.88], C2: [8.79, 9.32, 1.51], C1: [10.39, 9.45, 0.73],
  N3: [11.54, 10.82, 1.36], N2: [12.87, 11.15, 0.85], N1: [15.68, 11.7, -0.27],
  O3: [14.18, 12.92, 1.39], O2: [17.07, 13.79, 0.47],
  S3: [10.14, 9.13, 1.38], S2: [10.88, 9.49, 1.33],
  F: [14.66, 13.85, 2.31], Cl: [11.0, 9.69, 1.35], Br: [10.08, 8.47, 1.16], I: [9.9, 7.96, 0.96],
  P: [8.9, 8.24, 0.96], B: [5.98, 6.82, 1.605], Si: [7.3, 6.567, 0.657],
};

function paramsFor(el: string, hyb: number): [number, number, number] {
  if (PARAMS[el]) return PARAMS[el];
  const h = hyb === 1 ? 1 : hyb === 2 ? 2 : 3;
  return PARAMS[el + h] ?? PARAMS[el + '3'] ?? PARAMS.C3;
}

/**
 * Computes Gasteiger charges for an H-explicit molecule.
 * `hyb[i]` (1 = sp, 2 = sp2, 3 = sp3) selects the parameter set; `nbrs` lists covalent neighbours.
 */
export function gasteigerCharges(mol: Mol, hyb: ArrayLike<number>, nbrs: number[][], iterations = 6): number[] {
  const n = mol.atoms.length;
  const prm = mol.atoms.map((a, i) => paramsFor(a.el, hyb[i]));
  const chiPlus = prm.map((p, i) => (mol.atoms[i].el === 'H' ? 20.02 : p[0] + p[1] + p[2]));
  const q = mol.atoms.map((a) => a.charge);
  const chi = new Float64Array(n);
  for (let it = 1; it <= iterations; it++) {
    for (let i = 0; i < n; i++) chi[i] = prm[i][0] + prm[i][1] * q[i] + prm[i][2] * q[i] * q[i];
    const damp = Math.pow(0.5, it);
    const dq = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      for (const j of nbrs[i]) {
        if (j < i) continue;
        const diff = chi[j] - chi[i];
        const delta = (diff / (diff > 0 ? chiPlus[i] : chiPlus[j])) * damp;
        dq[i] += delta;
        dq[j] -= delta;
      }
    }
    for (let i = 0; i < n; i++) q[i] += dq[i];
  }
  return q;
}
