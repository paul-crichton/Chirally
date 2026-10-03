// Core molecular graph used by every chemistry algorithm.
//
// Coordinate convention: x to the right, y DOWN (screen convention), z toward the viewer.
// One model unit = one standard bond length (BOND_LENGTH = 1).
// File formats with y-up conventions (MOL, CDXML…) flip y on import/export.

export const BOND_LENGTH = 1;

export type BondStyle =
  | 'plain'
  | 'wedge'   // solid wedge, narrow end at bond.a (stereocentre), wide end at bond.b
  | 'hash'    // hashed wedge, narrow end at bond.a
  | 'wavy'    // unknown stereo
  | 'bold'    // thick bond (perspective)
  | 'dashed'  // dashed line (partial / perspective bond)
  | 'hollow'  // hollow wedge
  | 'dative'  // coordinate bond drawn as an arrow from a to b
  | 'hbond'   // dotted hydrogen bond (order 0)
  | 'crossed'; // crossed double bond (unknown E/Z)

export interface Atom {
  id: number;
  /** Element symbol ('C', 'N', …). Generic/pseudo atoms use 'R' (with alias) or '*'. */
  el: string;
  x: number;
  y: number;
  z?: number;
  charge: number;
  /** Mass number; undefined = natural abundance. */
  isotope?: number;
  /** Explicit hydrogen count override; undefined = computed from valence. */
  hCount?: number;
  /** Number of unpaired electrons (1 = radical, 2 = carbene/diradical). */
  radical?: number;
  /** Contracted abbreviation label (e.g. 'OMe', 'Ph', 'CO2H'); see abbreviations.ts. */
  abbrev?: string;
  /** Display text for pseudo atoms (R, R1, Ar, X, …). */
  alias?: string;
  /** Draw this atom's lone pairs. */
  lonePairs?: boolean;
  /** Reaction atom-atom mapping number. */
  map?: number;
  color?: string;
}

export interface Bond {
  id: number;
  a: number; // atom index (in Mol) / atom id (in Document)
  b: number;
  /** 1, 2, 3; 1.5 = aromatic/delocalised; 0 = zero-order (H-bond, ionic). */
  order: number;
  style: BondStyle;
  /** Placement of the second line of a double bond. */
  dbPos?: 'auto' | 'left' | 'right' | 'center';
  color?: string;
}

/**
 * Tetrahedral stereo specification in SMILES semantics.
 * Looking from nbrs[0] toward the centre, nbrs[1..3] are arranged counter-clockwise when ccw = true
 * (SMILES '@'), clockwise when false ('@@'). -1 denotes an implicit hydrogen or lone pair.
 */
export interface TetraSpec {
  center: number;
  nbrs: [number, number, number, number];
  ccw: boolean;
}

/** Double-bond stereo: atom `a` (neighbour of bond's first atom) and atom `b` (neighbour of bond's second atom) are cis or trans. */
export interface DbSpec {
  bond: number;
  a: number;
  b: number;
  cis: boolean;
}

export class Mol {
  atoms: Atom[] = [];
  bonds: Bond[] = [];
  tetra: TetraSpec[] = [];
  dbStereo: DbSpec[] = [];
  /** Free-form properties (e.g. SD file data fields, title). */
  props: Record<string, string> = {};
  name = '';
  private _adj: number[][] | null = null;

  get adj(): number[][] {
    if (!this._adj) {
      const adj: number[][] = this.atoms.map(() => []);
      this.bonds.forEach((b, i) => {
        adj[b.a].push(i);
        adj[b.b].push(i);
      });
      this._adj = adj;
    }
    return this._adj;
  }

  /** Must be called after structural edits made directly on atoms/bonds arrays. */
  invalidate(): void {
    this._adj = null;
  }

  addAtom(a: Partial<Atom> & { el: string }): number {
    const idx = this.atoms.length;
    this.atoms.push({ id: idx + 1, x: 0, y: 0, charge: 0, ...a } as Atom);
    this._adj = null;
    return idx;
  }

  addBond(a: number, b: number, order = 1, style: BondStyle = 'plain'): number {
    const idx = this.bonds.length;
    this.bonds.push({ id: idx + 1, a, b, order, style });
    this._adj = null;
    return idx;
  }

  neighbors(i: number): number[] {
    return this.adj[i].map((bi) => this.other(bi, i));
  }

  other(bi: number, ai: number): number {
    const b = this.bonds[bi];
    return b.a === ai ? b.b : b.a;
  }

  bondBetween(i: number, j: number): number {
    for (const bi of this.adj[i]) {
      if (this.other(bi, i) === j) return bi;
    }
    return -1;
  }

  degree(i: number): number {
    return this.adj[i].length;
  }

  heavyDegree(i: number): number {
    let n = 0;
    for (const j of this.neighbors(i)) if (this.atoms[j].el !== 'H') n++;
    return n;
  }

  clone(): Mol {
    const m = new Mol();
    m.atoms = this.atoms.map((a) => ({ ...a }));
    m.bonds = this.bonds.map((b) => ({ ...b }));
    m.tetra = this.tetra.map((t) => ({ ...t, nbrs: [...t.nbrs] as TetraSpec['nbrs'] }));
    m.dbStereo = this.dbStereo.map((d) => ({ ...d }));
    m.props = { ...this.props };
    m.name = this.name;
    return m;
  }

  /** Connected components as lists of atom indices. */
  components(): number[][] {
    const seen = new Array(this.atoms.length).fill(false);
    const out: number[][] = [];
    for (let s = 0; s < this.atoms.length; s++) {
      if (seen[s]) continue;
      const comp: number[] = [];
      const stack = [s];
      seen[s] = true;
      while (stack.length) {
        const v = stack.pop()!;
        comp.push(v);
        for (const w of this.neighbors(v)) {
          if (!seen[w]) {
            seen[w] = true;
            stack.push(w);
          }
        }
      }
      out.push(comp.sort((a, b) => a - b));
    }
    return out;
  }

  /**
   * Returns a new Mol containing only the given atoms (and bonds among them).
   * `map[oldIndex] = newIndex` (or -1). Stereo specs are carried over when fully contained.
   */
  subset(atomIdx: Iterable<number>): { mol: Mol; map: number[] } {
    const keep = new Set(atomIdx);
    const map = new Array(this.atoms.length).fill(-1);
    const m = new Mol();
    for (let i = 0; i < this.atoms.length; i++) {
      if (keep.has(i)) {
        map[i] = m.atoms.length;
        m.atoms.push({ ...this.atoms[i] });
      }
    }
    const bmap = new Array(this.bonds.length).fill(-1);
    this.bonds.forEach((b, bi) => {
      if (map[b.a] >= 0 && map[b.b] >= 0) {
        bmap[bi] = m.bonds.length;
        m.bonds.push({ ...b, a: map[b.a], b: map[b.b] });
      }
    });
    for (const t of this.tetra) {
      if (map[t.center] < 0) continue;
      if (t.nbrs.some((n) => n >= 0 && map[n] < 0)) continue;
      m.tetra.push({ center: map[t.center], nbrs: t.nbrs.map((n) => (n < 0 ? -1 : map[n])) as TetraSpec['nbrs'], ccw: t.ccw });
    }
    for (const d of this.dbStereo) {
      if (bmap[d.bond] < 0 || map[d.a] < 0 || map[d.b] < 0) continue;
      m.dbStereo.push({ bond: bmap[d.bond], a: map[d.a], b: map[d.b], cis: d.cis });
    }
    m.props = { ...this.props };
    m.name = this.name;
    return { mol: m, map };
  }

  /** Removes atoms (and incident bonds) in place; returns old→new index map. */
  removeAtoms(remove: Iterable<number>): number[] {
    const rm = new Set(remove);
    const keep: number[] = [];
    for (let i = 0; i < this.atoms.length; i++) if (!rm.has(i)) keep.push(i);
    const { mol, map } = this.subset(keep);
    this.atoms = mol.atoms;
    this.bonds = mol.bonds;
    this.tetra = mol.tetra;
    this.dbStereo = mol.dbStereo;
    this._adj = null;
    return map;
  }

  /** Appends another molecule; returns index offset of the appended atoms. */
  append(other: Mol): number {
    const off = this.atoms.length;
    const boff = this.bonds.length;
    let maxId = this.atoms.reduce((m, a) => Math.max(m, a.id), 0);
    for (const a of other.atoms) this.atoms.push({ ...a, id: ++maxId });
    let maxBid = this.bonds.reduce((m, b) => Math.max(m, b.id), 0);
    for (const b of other.bonds) this.bonds.push({ ...b, a: b.a + off, b: b.b + off, id: ++maxBid });
    for (const t of other.tetra)
      this.tetra.push({ center: t.center + off, nbrs: t.nbrs.map((n) => (n < 0 ? -1 : n + off)) as TetraSpec['nbrs'], ccw: t.ccw });
    for (const d of other.dbStereo) this.dbStereo.push({ bond: d.bond + boff, a: d.a + off, b: d.b + off, cis: d.cis });
    this._adj = null;
    return off;
  }

  bbox(): { minX: number; minY: number; maxX: number; maxY: number } {
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    for (const a of this.atoms) {
      minX = Math.min(minX, a.x); maxX = Math.max(maxX, a.x);
      minY = Math.min(minY, a.y); maxY = Math.max(maxY, a.y);
    }
    if (!this.atoms.length) return { minX: 0, minY: 0, maxX: 0, maxY: 0 };
    return { minX, minY, maxX, maxY };
  }

  translate(dx: number, dy: number): void {
    for (const a of this.atoms) {
      a.x += dx;
      a.y += dy;
    }
  }
}

/** Parity of the permutation that maps array `from` onto array `to` (same elements). true = even. */
export function permutationParityEven(from: readonly number[], to: readonly number[]): boolean {
  const perm = to.map((v) => from.indexOf(v));
  if (perm.some((p) => p < 0)) throw new Error('permutationParity: arrays differ');
  let even = true;
  const seen = new Array(perm.length).fill(false);
  for (let i = 0; i < perm.length; i++) {
    if (seen[i]) continue;
    let len = 0;
    let j = i;
    while (!seen[j]) {
      seen[j] = true;
      j = perm[j];
      len++;
    }
    if (len % 2 === 0) even = !even;
  }
  return even;
}

/**
 * Re-expresses a tetrahedral spec with a different neighbour ordering.
 * Returns the `ccw` flag valid for `order` (which must be a permutation of spec.nbrs).
 */
export function tetraCcwForOrder(spec: TetraSpec, order: readonly number[]): boolean {
  return permutationParityEven(spec.nbrs, order) ? spec.ccw : !spec.ccw;
}
