// UFF atom typing, hybridisation and bond orders for an H-explicit molecule.
import { Mol } from '../mol';
import { atomicNumber } from '../elements';
import { perceiveRings, RingInfo } from '../rings';
import { perceiveAromaticity } from '../aromaticity';
import { UFF_PARAMS, UFFAtomParams, defaultLabelFor, genericParams } from './uffParams';

export interface UFFTyping {
  /** UFF label per atom ('C_3', 'N_R', …; generic fallbacks end in '_g'). */
  labels: string[];
  params: UFFAtomParams[];
  /** 1 = sp (linear), 2 = sp2 (trigonal), 3 = sp3, 0 = no torsional treatment (metals, noble gases, …). */
  hyb: Int8Array;
  /** UFF bond order per Mol bond: 1, 1.41 (amide), 1.5 (aromatic), 2, 3; 0 = not covalent (ignored). */
  bondOrder: Float64Array;
  aromaticAtom: boolean[];
  rings: RingInfo;
  /** Covalent neighbours (bonds with bondOrder > 0). */
  nbrs: number[][];
}

const GROUP6 = new Set([8, 16, 34, 52, 84]);
const SP3_MAIN_GROUP = new Set(['F', 'Cl', 'Br', 'I', 'Si', 'Ge', 'Sn', 'Pb', 'As', 'Se', 'Te', 'Sb', 'Bi', 'Al', 'Ga', 'In']);
export const isGroup6 = (el: string): boolean => GROUP6.has(atomicNumber(el));

/** True if a bond is covalent for force-field purposes (excludes H-bonds and zero-order bonds). */
export function isCovalent(order: number, style: string): boolean {
  return order > 0 && style !== 'hbond';
}

/** Assigns UFF atom types. The molecule should carry explicit hydrogens (see addExplicitHydrogens). */
export function typeUFF(mol: Mol): UFFTyping {
  const n = mol.atoms.length;
  const rings = perceiveRings(mol);
  const aro = perceiveAromaticity(mol, rings);
  const nbrs: number[][] = mol.atoms.map(() => []);
  const isArBond = mol.bonds.map((b, bi) => aro.bonds[bi] || b.order === 1.5);
  mol.bonds.forEach((b) => {
    if (!isCovalent(b.order, b.style)) return;
    nbrs[b.a].push(b.b);
    nbrs[b.b].push(b.a);
  });

  // Bond-order statistics per atom
  const nDouble = new Int8Array(n), nTriple = new Int8Array(n), nArom = new Int8Array(n);
  const vsum = new Float64Array(n);
  mol.bonds.forEach((b, bi) => {
    if (!isCovalent(b.order, b.style)) return;
    for (const a of [b.a, b.b]) {
      if (isArBond[bi]) { nArom[a]++; vsum[a] += 1.5; }
      else if (b.order >= 3) { nTriple[a]++; vsum[a] += 3; }
      else if (b.order === 2) { nDouble[a]++; vsum[a] += 2; }
      else vsum[a] += b.order;
    }
  });
  const hasPi = (j: number) => nDouble[j] + nTriple[j] + nArom[j] > 0;
  // Conjugation for lone-pair donors: bonded to a C/N/O/B atom carrying a multiple or aromatic bond.
  const conjugated = (i: number) =>
    nbrs[i].some((j) => ['C', 'N', 'O', 'B'].includes(mol.atoms[j].el) && hasPi(j));

  const labels: string[] = new Array(n);
  const hyb = new Int8Array(n);
  for (let i = 0; i < n; i++) {
    const a = mol.atoms[i];
    const deg = nbrs[i].length;
    const arom = aro.atoms[i] || nArom[i] >= 2;
    let label = '';
    let h = 3;
    switch (a.el) {
      case 'H':
        label = deg >= 2 ? 'H_b' : 'H_';
        h = 0;
        break;
      case 'C':
        if (nTriple[i] || nDouble[i] >= 2) { label = 'C_1'; h = 1; }
        else if (arom) { label = 'C_R'; h = 2; }
        else if (nDouble[i] || nArom[i]) { label = 'C_2'; h = 2; }
        else if (deg === 3 && (a.charge > 0 || (a.radical ?? 0) > 0)) { label = 'C_2'; h = 2; } // planar cation/radical
        else label = 'C_3';
        break;
      case 'N':
        if (nTriple[i] || nDouble[i] >= 2) { label = 'N_1'; h = 1; }
        else if (arom) { label = 'N_R'; h = 2; }
        else if (nDouble[i] || nArom[i]) { label = deg >= 3 ? 'N_R' : 'N_2'; h = 2; }
        else if (deg === 3 && a.charge === 0 && conjugated(i)) { label = 'N_R'; h = 2; } // amide, aniline, enamine
        else label = 'N_3';
        break;
      case 'O':
        if (nTriple[i]) { label = 'O_1'; h = 1; }
        else if (nDouble[i]) { label = 'O_2'; h = 2; }
        else if (arom) { label = 'O_R'; h = 2; }
        else if (deg === 1 && a.charge < 0 && conjugated(i)) { label = 'O_2'; h = 2; } // carboxylate, nitro O-
        else if (deg === 2 && conjugated(i)) { label = 'O_R'; h = 2; } // ester, acid, phenol, enol ether
        else label = 'O_3';
        break;
      case 'S':
        if (arom) { label = 'S_R'; h = 2; }
        else if (deg === 1 && nDouble[i]) { label = 'S_2'; h = 2; }
        else if (deg >= 4 || vsum[i] >= 6) label = 'S_3+6';
        else if (deg === 3 || vsum[i] >= 4) label = 'S_3+4';
        else label = 'S_3+2';
        break;
      case 'P':
        if (deg >= 4 || vsum[i] >= 5) label = a.charge > 0 && !nDouble[i] ? 'P_3+q' : 'P_3+5';
        else label = 'P_3+3';
        break;
      case 'B':
        if (deg >= 4) label = 'B_3';
        else { label = 'B_2'; h = 2; }
        break;
      default:
        label = defaultLabelFor(a.el) ?? '';
        // Main-group elements are treated as sp3; metals and noble gases get no torsions.
        h = SP3_MAIN_GROUP.has(a.el) ? 3 : 0;
    }
    labels[i] = label;
    hyb[i] = h;
  }

  const params = labels.map((l, i) => {
    const p = l ? UFF_PARAMS.get(l) : undefined;
    if (p) return p;
    const g = genericParams(mol.atoms[i].el);
    labels[i] = g.label;
    return g;
  });

  // UFF bond orders: aromatic 1.5, amide C–N 1.41
  const bondOrder = new Float64Array(mol.bonds.length);
  const isCarbonylC = (c: number) =>
    mol.atoms[c].el === 'C' && hyb[c] === 2 &&
    nbrs[c].some((o) => mol.atoms[o].el === 'O' && labels[o] === 'O_2' && mol.bonds[mol.bondBetween(c, o)].order === 2);
  mol.bonds.forEach((b, bi) => {
    if (!isCovalent(b.order, b.style)) { bondOrder[bi] = 0; return; }
    if (isArBond[bi]) { bondOrder[bi] = 1.5; return; }
    if (b.order === 1) {
      const [c, nn] = mol.atoms[b.a].el === 'C' ? [b.a, b.b] : [b.b, b.a];
      if (mol.atoms[nn].el === 'N' && (labels[nn] === 'N_R' || labels[nn] === 'N_2') && isCarbonylC(c)) {
        bondOrder[bi] = 1.41;
        return;
      }
    }
    bondOrder[bi] = Math.max(1, Math.min(4, b.order));
  });

  return { labels, params, hyb, bondOrder, aromaticAtom: aro.atoms, rings, nbrs };
}

/** UFF natural bond length (Å) with bond-order (Pauling) and electronegativity (O'Keefe–Breese) corrections. */
export function uffBondLength(order: number, p1: UFFAtomParams, p2: UFFAtomParams): number {
  const ri = p1.r1, rj = p2.r1;
  const rBO = -0.1332 * (ri + rj) * Math.log(order);
  const sx = Math.sqrt(p1.Xi) - Math.sqrt(p2.Xi);
  const rEN = (ri * rj * sx * sx) / (p1.Xi * ri + p2.Xi * rj);
  return ri + rj + rBO - rEN;
}
