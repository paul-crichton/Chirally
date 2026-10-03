// Normalised heavy-atom graph used by the IUPAC namer.
//
// Built from an (abbreviation-expanded) Mol: ordinary hydrogens are folded into per-atom
// counts, delocalised (order 1.5) bonds are kekulized, and ring information is perceived.
// The namer never looks at coordinates.
import { Mol } from '../mol';
import { implicitH } from '../valence';
import { kekulize } from '../aromaticity';
import { perceiveRings, RingInfo } from '../rings';

export const HALOGENS = new Set(['F', 'Cl', 'Br', 'I']);

export class NGraph {
  /** Heavy-atom molecule (Kekulé bonds, explicit hCount on every atom). */
  mol: Mol;
  n: number;
  el: string[] = [];
  charge: number[] = [];
  h: number[] = [];
  iso: (number | undefined)[] = [];
  /** Index of the atom in the source Mol. */
  src: number[] = [];
  /** Source bond index of each bond (-1 if none). */
  srcBond: number[] = [];
  nb: number[][] = [];
  rings: RingInfo;
  inRing: boolean[];
  /** Ring-system index of each atom (-1 if acyclic). */
  sysOf: number[];
  /** Atoms of each ring system. */
  sysAtoms: number[][];
  private bondMap = new Map<number, number>();

  constructor(mol: Mol, src: number[], srcBond: number[]) {
    this.mol = mol;
    this.n = mol.atoms.length;
    this.src = src;
    this.srcBond = srcBond;
    for (let i = 0; i < this.n; i++) {
      const a = mol.atoms[i];
      this.el.push(a.el);
      this.charge.push(a.charge);
      this.h.push(a.hCount ?? implicitH(mol, i));
      this.iso.push(a.isotope);
      this.nb.push(mol.neighbors(i));
    }
    mol.bonds.forEach((b, bi) => {
      this.bondMap.set(b.a * 100003 + b.b, bi);
      this.bondMap.set(b.b * 100003 + b.a, bi);
    });
    this.rings = perceiveRings(mol);
    this.inRing = this.rings.inRing;
    this.sysOf = new Array(this.n).fill(-1);
    this.sysAtoms = [];
    // Ring systems: group SSSR systems (shared atoms) – rings.ts already provides them.
    this.rings.systems.forEach((sys, k) => {
      const set = new Set<number>();
      for (const r of sys) for (const a of this.rings.rings[r]) set.add(a);
      const atoms = [...set].sort((p, q) => p - q);
      for (const a of atoms) this.sysOf[a] = k;
      this.sysAtoms.push(atoms);
    });
  }

  bond(i: number, j: number): number {
    return this.bondMap.get(i * 100003 + j) ?? -1;
  }

  /** Bond order between i and j (0 if not bonded). */
  order(i: number, j: number): number {
    const b = this.bond(i, j);
    return b < 0 ? 0 : this.mol.bonds[b].order;
  }

  degree(i: number): number {
    return this.nb[i].length;
  }

  /** Neighbours of i other than those in `exclude`. */
  others(i: number, ...exclude: number[]): number[] {
    return this.nb[i].filter((j) => !exclude.includes(j));
  }

  /** Neighbour of i reached by a double bond to an atom of element `el` that has no other neighbours. */
  terminalDouble(i: number, el: string): number {
    for (const j of this.nb[i]) if (this.el[j] === el && this.order(i, j) === 2 && this.degree(j) === 1) return j;
    return -1;
  }

  countTerminalDouble(i: number, el: string): number {
    let c = 0;
    for (const j of this.nb[i]) if (this.el[j] === el && this.order(i, j) === 2 && this.degree(j) === 1 && this.charge[j] === 0) c++;
    return c;
  }
}

export interface BuildGraphResult {
  graph: NGraph | null;
  error?: string;
  warnings: string[];
}

/**
 * Builds the naming graph from a Mol whose abbreviations have already been expanded.
 * Hydrogen atoms bonded to one heavy atom are folded into H counts.
 */
export function buildGraph(em: Mol): BuildGraphResult {
  const warnings: string[] = [];
  const n = em.atoms.length;
  const keep: number[] = [];
  const fold = new Map<number, number>(); // heavy atom → folded explicit H
  const isFolded = new Array(n).fill(false);
  for (let i = 0; i < n; i++) {
    const a = em.atoms[i];
    if (a.el === 'H' && !a.charge && em.degree(i) === 1) {
      const j = em.neighbors(i)[0];
      const b = em.bonds[em.adj[i][0]];
      if (em.atoms[j].el !== 'H' && b.order === 1) {
        if (a.isotope) warnings.push('isotopic hydrogen treated as ordinary hydrogen');
        isFolded[i] = true;
        fold.set(j, (fold.get(j) ?? 0) + 1);
        continue;
      }
    }
    keep.push(i);
  }
  const map = new Array(n).fill(-1);
  const m = new Mol();
  for (const i of keep) {
    const a = em.atoms[i];
    if (a.el === 'R' || a.el === '*' || a.alias) {
      return { graph: null, error: 'structure contains pseudo atoms (R/*)', warnings };
    }
    if (a.radical) return { graph: null, error: 'radicals are not supported', warnings };
    map[i] = m.atoms.length;
    // total H: implicit (computed on the source mol) + folded explicit H atoms
    let h = implicitH(em, i);
    h += fold.get(i) ?? 0;
    m.atoms.push({ id: m.atoms.length + 1, el: a.el, x: 0, y: 0, charge: a.charge, isotope: a.isotope, hCount: h });
  }
  const srcBond: number[] = [];
  const aroAtoms = new Array(m.atoms.length).fill(false);
  const aroBonds: boolean[] = [];
  em.bonds.forEach((b, bi) => {
    const p = map[b.a], q = map[b.b];
    if (p < 0 || q < 0) return;
    if (b.style === 'hbond' || b.order === 0) return;
    const idx = m.bonds.length;
    m.bonds.push({ id: idx + 1, a: p, b: q, order: b.order, style: 'plain' });
    srcBond.push(bi);
    if (b.order === 1.5) {
      aroAtoms[p] = aroAtoms[q] = true;
      aroBonds[idx] = true;
    } else aroBonds[idx] = false;
  });
  m.invalidate();
  if (aroBonds.some((x) => x)) {
    const ok = kekulize(m, aroAtoms, aroBonds);
    if (!ok || m.bonds.some((b) => b.order === 1.5)) {
      return { graph: null, error: 'could not assign a Kekulé structure to the aromatic bonds', warnings };
    }
  }
  if (m.bonds.some((b) => b.order !== 1 && b.order !== 2 && b.order !== 3)) {
    return { graph: null, error: 'unsupported bond order', warnings };
  }
  if (m.atoms.some((a) => a.isotope)) warnings.push('isotope labels are not expressed in the name');
  const src = keep.slice();
  return { graph: new NGraph(m, src, srcBond), warnings };
}
