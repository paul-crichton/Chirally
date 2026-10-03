// CIP stereodescriptors (R/S, r/s, E/Z) from mol.tetra / mol.dbStereo.
//
// Hierarchical-digraph implementation of the IUPAC 2013 rules (P-92) in the machine-oriented form
// of Hanson, Musacchio, Mayfield, Vainio, Yerin, Redkin, J. Chem. Inf. Model. 2018, 58, 1755
// (the algorithm behind RDKit's new CIP labeller):
//  • The digraph is rooted at the stereocentre (or a double-bond atom) and expanded lazily. Multiple
//    bonds give duplicate atoms on both ends — except for the root's own multiple bonds (S=O, P=O),
//    which are treated as single; ring closures give a duplicate of the revisited atom. Duplicates have
//    only phantom substituents (atomic number 0), as does a lone pair on an S/P/N centre. Implicit
//    hydrogens are ordinary atoms (Z = 1).
//  • Mancude rings are Kekulé-independent: aromatic bonds count as single and every atom needing a
//    π partner gets one duplicate whose atomic number is the mean over its possible double-bond
//    partners (P-92.1.4.4), so pyridine C2 gets 6.5 whichever Kekulé structure was drawn.
//  • Branches are compared sphere by sphere, visiting branches in the order of their own rank, and
//    each rule is exhausted over the whole digraph before the next: 1a atomic number; 1b ring
//    duplicates first, nearer duplicated atoms first; 2 mass (only if an isotope is specified);
//    3 Z > E; 4a chiral > pseudoasymmetric (or E/Z) > none; 4b like > unlike; 4c r > s; 5 R > S.
//  • Stereo units inside branches carry auxiliary descriptors computed in the same digraph, re-rooted
//    at that unit (deepest first: a unit only sees descriptors of farther nodes). This yields e.g.
//    (1r,4r) for trans-tranexamic acid and the pseudoasymmetric C3 of 2,3,4-trihydroxyglutaric acid.
//  • Rules 4b/5 at a root compare descriptor sequences of the branches relative to a reference
//    (R-like/S-like); a decision by rule 5 whose R- and S-referenced comparisons disagree marks the
//    centre pseudoasymmetric → lower-case r/s.
//  • Centres with two ligands that stay tied are not stereocentres and receive no label.
import { Mol, TetraSpec, DbSpec, tetraCcwForOrder } from './mol';
import { element } from './elements';
import { implicitH } from './valence';
import { perceiveRings } from './rings';
import { perceiveAromaticity } from './aromaticity';

export type CenterLabel = 'R' | 'S' | 'r' | 's';
export type BondLabel = 'E' | 'Z';

export interface CIPResult {
  /** atom index → descriptor */
  centers: Map<number, CenterLabel>;
  /** bond index → descriptor */
  bonds: Map<number, BondLabel>;
}

// Descriptor codes used for auxiliary descriptors
const NONE = 0, D_R = 1, D_S = 2, D_r = 3, D_s = 4, D_E = 5, D_Z = 6, D_SEQCIS = 7, D_SEQTRANS = 8;
type Desc = number;
/** R-like / S-like reference of a descriptor (PairList semantics): 1 = R, 2 = S, 0 = none. */
const refOf = (d: Desc) => (d === D_R || d === D_SEQCIS ? 1 : d === D_S || d === D_SEQTRANS ? 2 : 0);

/** Digraph node. Traversals use (node, from) handles so that the same digraph can be re-rooted. */
interface Node {
  /** Atom index (real atom, or the atom a duplicate stands for); -1 for implicit H / phantom / mancude duplicate. */
  atom: number;
  /** TetraSpec ligand slot: atom index, -1 (implicit H / lone pair), -2 (multiple-bond or mancude duplicate). */
  slot: number;
  z: number;
  mass: number;
  /** An isotope mass number was specified (rule 2 only applies then). */
  labelled: boolean;
  dup: boolean;
  ringDup: boolean;
  /** Depth (distance from the original root); for ring duplicates the depth of the duplicated atom. */
  dist: number;
  depth: number;
  parent: Node | null;
  /** Bond to the parent (-1 for none). */
  bond: number;
  /** Placeholder for the other atom of a double bond above an E/Z root (never traversed). */
  virtual?: boolean;
  kids?: Node[];
  sortCache?: Map<Node | null, Map<number, Node[]>>;
  tetraAux?: Desc;
  dbAux?: Desc;
}

// Node-level rules (BFS comparison); rules 4b and 5 only act at a root (see rootCompare).
const R1A = 0, R1B = 1, R2 = 2, R3 = 3, R4A = 4, R4C = 5;
const NODE_RULES = [R1A, R1B, R2, R3, R4A, R4C];
const MAX_LEVEL = NODE_RULES.length - 1;

/** Digraph size limits: per stereo unit, and per assignCIP call (bounds pathological cages). */
const NODE_BUDGET = 20000;
const TOTAL_BUDGET = 250000;
class BudgetExceeded extends Error {}

interface AdjEntry {
  nbr: number;
  bond: number;
  /** Bond order used for duplication: aromatic (mancude) bonds count as 1. */
  order: number;
}

class CIPContext {
  readonly n: number;
  readonly z: number[];
  readonly weight: number[];
  readonly isoMass: number[];
  readonly hcount: number[];
  readonly adj: AdjEntry[][];
  readonly mancudeZ: (number | undefined)[];
  readonly tetra = new Map<number, TetraSpec>();
  readonly db = new Map<number, DbSpec>();
  /** Atom → stereo double bonds it belongs to. */
  readonly dbOfAtom = new Map<number, DbSpec>();
  private nodes = 0;
  private totalNodes = 0;
  private origin: Node | null = null;
  /** Auxiliary descriptors are visible only on nodes deeper than this. */
  private horizon = 0;

  constructor(readonly mol: Mol) {
    const n = (this.n = mol.atoms.length);
    this.z = mol.atoms.map((a) => (a.abbrev ? 0 : element(a.el)?.z ?? 0));
    this.weight = mol.atoms.map((a) => (a.abbrev ? 0 : element(a.el)?.mass ?? 0));
    this.isoMass = mol.atoms.map((a) => a.isotope ?? 0);
    this.hcount = mol.atoms.map((a, i) => (a.abbrev || !element(a.el) ? 0 : implicitH(mol, i)));
    const rings = perceiveRings(mol);
    const aro = perceiveAromaticity(mol, rings);
    const mancudeBond = mol.bonds.map((b, i) => aro.bonds[i] || b.order === 1.5);
    this.adj = mol.atoms.map(() => []);
    mol.bonds.forEach((b, bi) => {
      if (b.order === 0 || b.style === 'hbond') return;
      const order = mancudeBond[bi] ? 1 : Math.max(1, Math.round(b.order));
      this.adj[b.a].push({ nbr: b.b, bond: bi, order });
      this.adj[b.b].push({ nbr: b.a, bond: bi, order });
    });
    this.mancudeZ = new Array(n).fill(undefined);
    this.perceiveMancude(mancudeBond);
    for (const t of mol.tetra) this.tetra.set(t.center, t);
    for (const d of mol.dbStereo) {
      const b = mol.bonds[d.bond];
      if (!b || b.order !== 2) continue;
      this.db.set(d.bond, d);
      this.dbOfAtom.set(b.a, d);
      this.dbOfAtom.set(b.b, d);
    }
  }

  /**
   * Atoms carrying a π bond inside an aromatic system get one duplicate whose atomic number is the
   * mean over the partners they are double-bonded to in at least one Kekulé structure.
   */
  private perceiveMancude(mancudeBond: boolean[]): void {
    const mol = this.mol;
    const need = new Array(this.n).fill(false);
    for (let i = 0; i < this.n; i++) {
      let dbl = 0, delocal = 0;
      for (const bi of mol.adj[i]) {
        if (!mancudeBond[bi]) continue;
        if (mol.bonds[bi].order === 2) dbl++;
        if (mol.bonds[bi].order === 1.5) delocal++;
      }
      need[i] = dbl === 1 || delocal >= 2;
    }
    const partners = (i: number) =>
      mol.adj[i].filter((bi) => mancudeBond[bi] && need[mol.other(bi, i)]).map((bi) => mol.other(bi, i));
    const seen = new Array(this.n).fill(false);
    for (let s = 0; s < this.n; s++) {
      if (!need[s] || seen[s]) continue;
      const comp: number[] = [];
      const st = [s];
      seen[s] = true;
      while (st.length) {
        const v = st.pop()!;
        comp.push(v);
        for (const w of partners(v)) if (!seen[w]) { seen[w] = true; st.push(w); }
      }
      comp.sort((a, b) => a - b);
      // enumerate Kekulé structures (perfect matchings) and record which partners occur
      const occurs = new Map<number, Set<number>>();
      const match = new Map<number, number>();
      let found = 0;
      let steps = 0;
      const rec = (): void => {
        if (found >= 2048 || ++steps > 100000) return;
        const v = comp.find((a) => !match.has(a));
        if (v === undefined) {
          found++;
          for (const [a, b] of match) {
            if (!occurs.has(a)) occurs.set(a, new Set());
            occurs.get(a)!.add(b);
          }
          return;
        }
        for (const w of partners(v)) {
          if (match.has(w)) continue;
          match.set(v, w);
          match.set(w, v);
          rec();
          match.delete(v);
          match.delete(w);
        }
      };
      rec();
      for (const a of comp) {
        const p = found ? [...(occurs.get(a) ?? [])] : partners(a);
        if (p.length) this.mancudeZ[a] = p.reduce((t, b) => t + this.z[b], 0) / p.length;
      }
    }
  }

  // ─────────────────────────── digraph construction ───────────────────────────

  private make(init: Omit<Node, 'kids' | 'sortCache' | 'tetraAux' | 'dbAux'>): Node {
    if (++this.nodes > NODE_BUDGET || ++this.totalNodes > TOTAL_BUDGET) throw new BudgetExceeded();
    return init;
  }

  private atomNode(atom: number, parent: Node | null, bond: number, depth: number): Node {
    return this.make({
      atom, slot: atom, z: this.z[atom], mass: this.isoMass[atom] || this.weight[atom], labelled: this.isoMass[atom] > 0,
      dup: false, ringDup: false, dist: depth, depth, parent, bond,
    });
  }

  private dupNode(atom: number, parent: Node, ring: boolean, dist: number, slot: number): Node {
    return this.make({
      atom, slot, z: this.z[atom], mass: this.weight[atom], labelled: false,
      dup: true, ringDup: ring, dist, depth: parent.depth + 1, parent, bond: -1,
    });
  }

  private hNode(parent: Node): Node {
    const d = parent.depth + 1;
    return this.make({ atom: -1, slot: -1, z: 1, mass: element('H')!.mass, labelled: false, dup: false, ringDup: false, dist: d, depth: d, parent, bond: -1 });
  }

  private lonePair(parent: Node): Node {
    return this.make({ atom: -1, slot: -1, z: 0, mass: 0, labelled: false, dup: true, ringDup: false, dist: parent.depth + 1, depth: parent.depth + 1, parent, bond: -1 });
  }

  /** Starts a new digraph rooted at `atom`; `via` (other atom of a double bond) counts as visited. */
  private newRoot(atom: number, via = -1): Node {
    this.nodes = 0;
    let parent: Node | null = null;
    if (via >= 0) {
      parent = this.atomNode(via, null, -1, -1);
      parent.virtual = true;
    }
    const root = this.atomNode(atom, parent, -1, 0);
    this.origin = root;
    this.horizon = 0;
    return root;
  }

  /** Depth of the (non-duplicate) ancestor carrying `atom`, or null if it is not on the path. */
  private ancestorDepth(node: Node, atom: number): number | null {
    for (let p: Node | null = node; p; p = p.parent) if (p.atom === atom && !p.dup) return p.depth;
    return null;
  }

  private kidsOf(node: Node): Node[] {
    if (node.kids) return node.kids;
    const kids: Node[] = [];
    node.kids = kids;
    if (node.dup || node.atom < 0) return kids;
    const x = node.atom;
    for (const e of this.adj[x]) {
      if (e.bond === node.bond) {
        // bond back to the parent: only its multiple-bond duplicates (never towards the root)
        if (node.parent !== this.origin) for (let k = 1; k < e.order; k++) kids.push(this.dupNode(e.nbr, node, false, node.depth - 1, -2));
        continue;
      }
      const anc = this.ancestorDepth(node, e.nbr);
      if (anc !== null) {
        kids.push(this.dupNode(e.nbr, node, true, anc, e.nbr));
        for (let k = 1; k < e.order; k++) kids.push(this.dupNode(e.nbr, node, false, anc, -2));
      } else {
        kids.push(this.atomNode(e.nbr, node, e.bond, node.depth + 1));
        for (let k = 1; k < e.order; k++) kids.push(this.dupNode(e.nbr, node, false, node.depth + 1, -2));
      }
    }
    const mz = this.mancudeZ[x];
    if (mz !== undefined) {
      kids.push(this.make({ atom: -1, slot: -2, z: mz, mass: 0, labelled: false, dup: true, ringDup: false, dist: node.depth + 1, depth: node.depth + 1, parent: node, bond: -1 }));
    }
    for (let k = 0; k < this.hcount[x]; k++) kids.push(this.hNode(node));
    return kids;
  }

  /** Neighbours of `n` in the digraph seen from `from` (kids plus parent, minus where we came from). */
  private children(n: Node, from: Node | null): Node[] {
    const out = this.kidsOf(n).filter((k) => k !== from);
    const p = n.parent;
    if (p && p !== from && !p.virtual) out.push(p);
    return out;
  }

  private sortedChildren(n: Node, from: Node | null, level: number): Node[] {
    if (!n.sortCache) n.sortCache = new Map();
    let byKey = n.sortCache.get(from);
    if (!byKey) {
      byKey = new Map();
      n.sortCache.set(from, byKey);
    }
    // descriptor-based rules depend on which auxiliary descriptors are visible
    const key = level * 100000 + (level >= R3 ? this.horizon + 1 : 0);
    let s = byKey.get(key);
    if (!s) {
      s = this.children(n, from).sort((a, b) => this.cmpBranch(b, a, n, level));
      byKey.set(key, s);
    }
    return s;
  }

  // ─────────────────────────── rules ───────────────────────────

  private visible(n: Node): boolean {
    return n.depth > this.horizon && !n.dup && n.atom >= 0;
  }

  private prop(n: Node, rule: number): number {
    switch (rule) {
      case R1A: return n.z;
      case R1B: return n.ringDup ? 1e6 - n.dist : 0;
      case R3: {
        if (!this.visible(n)) return 0;
        const d = this.dbAux(n);
        return d === D_Z ? 2 : d === D_E ? 1 : 0;
      }
      case R4A: {
        if (!this.visible(n)) return 0;
        const t = this.tetraAux(n);
        if (t === D_R || t === D_S) return 2;
        if (t === D_r || t === D_s) return 1;
        const d = this.dbAux(n);
        if (d === D_SEQCIS || d === D_SEQTRANS) return 2;
        if (d === D_E || d === D_Z) return 1;
        return 0;
      }
      case R4C: {
        if (!this.visible(n)) return 0;
        const t = this.tetraAux(n);
        return t === D_r ? 2 : t === D_s ? 1 : 0;
      }
    }
    return 0;
  }

  private cmpProp(a: Node | undefined, b: Node | undefined, rule: number): number {
    if (rule === R2) {
      // Rule 2 only distinguishes when an isotope is specified on either side
      if (!(a?.labelled || b?.labelled)) return 0;
      const ma = a ? a.mass : 0, mb = b ? b.mass : 0;
      return ma === mb ? 0 : ma > mb ? 1 : -1;
    }
    const pa = a ? this.prop(a, rule) : 0;
    const pb = b ? this.prop(b, rule) : 0;
    return pa === pb ? 0 : pa > pb ? 1 : -1;
  }

  /** Breadth-first comparison of branches a (entered from fa) and b (from fb) under one rule. */
  private cmpRule(a: Node, fa: Node | null, b: Node, fb: Node | null, rule: number): number {
    let c = this.cmpProp(a, b, rule);
    if (c) return c;
    let qa: [Node, Node | null][] = [[a, fa]];
    let qb: [Node, Node | null][] = [[b, fb]];
    while (qa.length) {
      const na: [Node, Node | null][] = [];
      const nb: [Node, Node | null][] = [];
      for (let i = 0; i < qa.length; i++) {
        const [x, fx] = qa[i];
        const [y, fy] = qb[i];
        const kx = this.sortedChildren(x, fx, rule);
        const ky = this.sortedChildren(y, fy, rule);
        const len = Math.max(kx.length, ky.length);
        for (let j = 0; j < len; j++) {
          c = this.cmpProp(kx[j], ky[j], rule);
          if (c) return c;
        }
        if (kx.length !== ky.length) return kx.length > ky.length ? 1 : -1;
        for (let j = 0; j < kx.length; j++) {
          if (kx[j].atom < 0 && ky[j].atom < 0) continue; // leaves (H, phantoms, mancude duplicates)
          na.push([kx[j], x]);
          nb.push([ky[j], y]);
        }
      }
      qa = na;
      qb = nb;
    }
    return 0;
  }

  /** Hierarchical comparison of two children of `from`, applying node rules 0..level in turn. */
  private cmpBranch(a: Node, b: Node, from: Node | null, level: number): number {
    for (let r = 0; r <= level; r++) {
      const c = this.cmpRule(a, from, b, from, NODE_RULES[r]);
      if (c) return c;
    }
    return 0;
  }

  // ─────────────────────────── rules 4b and 5 (root only) ───────────────────────────

  /** Visible auxiliary descriptors of a branch in hierarchical order (levels of tied groups). */
  private levels(start: Node, from: Node | null): Node[][] {
    const out: Node[][] = [];
    let q: [Node, Node | null][] = [[start, from]];
    while (q.length) {
      out.push(q.map((e) => e[0]));
      const next: [Node, Node | null][] = [];
      for (const [n, f] of q) for (const k of this.sortedChildren(n, f, MAX_LEVEL)) if (!k.dup && k.atom >= 0) next.push([k, n]);
      q = next;
    }
    return out;
  }

  private descOf(n: Node): Desc {
    if (!this.visible(n)) return NONE;
    const t = this.tetraAux(n);
    return t !== NONE ? t : this.dbAux(n);
  }

  /** Reference descriptor(s) of a branch: majority R-like/S-like of the first descriptor-bearing level. */
  private references(lv: Node[][]): number[] {
    for (const level of lv) {
      let r = 0, s = 0;
      for (const n of level) {
        const ref = refOf(this.descOf(n));
        if (ref === 1) r++;
        else if (ref === 2) s++;
      }
      if (r + s === 0) continue;
      return r > s ? [1] : s > r ? [2] : [1, 2];
    }
    return [];
  }

  /** like(1)/unlike(0) sequence of a branch's descriptors relative to `ref` (hierarchical order). */
  private likeSeq(lv: Node[][], ref: number): number[] {
    const out: number[] = [];
    for (const level of lv)
      for (const n of level) {
        const r = refOf(this.descOf(n));
        if (r) out.push(r === ref ? 1 : 0);
      }
    return out;
  }

  private static cmpSeq(a: number[], b: number[]): number {
    const len = Math.min(a.length, b.length);
    for (let k = 0; k < len; k++) if (a[k] !== b[k]) return a[k] > b[k] ? 1 : -1;
    return 0;
  }

  private rule4b(a: Node, b: Node, root: Node): number {
    const la = this.levels(a, root), lb = this.levels(b, root);
    const refA = this.references(la), refB = this.references(lb);
    if (!refA.length || !refB.length) return 0;
    if (refA.length === 1 && refB.length === 1) {
      // positional like/unlike comparison of the two branches against their own reference
      const fa = la.flat(), fb = lb.flat();
      const len = Math.min(fa.length, fb.length);
      for (let k = 0; k < len; k++) {
        const likeA = refOf(this.descOf(fa[k])) === refA[0];
        const likeB = refOf(this.descOf(fb[k])) === refB[0];
        if (likeA !== likeB) return likeA ? 1 : -1;
      }
      return 0;
    }
    // several equally possible references: compare the best like/unlike lists
    const best = (lv: Node[][], refs: number[]) =>
      refs.map((r) => this.likeSeq(lv, r)).sort((p, q) => -CIPContext.cmpSeq(p, q))[0] ?? [];
    return CIPContext.cmpSeq(best(la, refA), best(lb, refB));
  }

  /** Rule 5: ±1 decided, ±2 decided with enantiomorphic ligands (pseudoasymmetry). */
  private rule5(a: Node, b: Node, root: Node): number {
    const la = this.levels(a, root), lb = this.levels(b, root);
    const cmpR = CIPContext.cmpSeq(this.likeSeq(la, 1), this.likeSeq(lb, 1));
    const cmpS = CIPContext.cmpSeq(this.likeSeq(la, 2), this.likeSeq(lb, 2));
    if (cmpR < 0) return cmpS < 0 ? -1 : -2;
    if (cmpR > 0) return cmpS > 0 ? 1 : 2;
    return 0;
  }

  /** Full comparison of two ligands of `root`; |result| = 2 flags a rule-5 pseudoasymmetric decision. */
  private rootCompare(a: Node, b: Node, root: Node): number {
    for (const r of [R1A, R1B, R2, R3, R4A]) {
      const c = this.cmpRule(a, root, b, root, r);
      if (c) return c;
    }
    let c = this.rule4b(a, b, root);
    if (c) return c;
    c = this.cmpRule(a, root, b, root, R4C);
    if (c) return c;
    return this.rule5(a, b, root);
  }

  /** Sorts ligands of `root` (highest first); null when two of them are equivalent. */
  private rank(ligs: Node[], root: Node): { order: Node[]; pseudo: boolean } | null {
    const order = [...ligs];
    let pseudo = false;
    for (let i = 1; i < order.length; i++) {
      for (let j = i; j > 0; j--) {
        const c = this.rootCompare(order[j - 1], order[j], root);
        if (c === 0) return null;
        if (Math.abs(c) === 2) pseudo = true;
        if (c > 0) break;
        [order[j - 1], order[j]] = [order[j], order[j - 1]];
      }
    }
    return { order, pseudo };
  }

  // ─────────────────────────── auxiliary descriptors ───────────────────────────

  /** Runs f with the digraph re-rooted at `n` (only deeper descriptors visible). */
  private withHorizon<T>(n: Node, f: () => T): T {
    const saved = this.horizon;
    this.horizon = n.depth;
    try {
      return f();
    } finally {
      this.horizon = saved;
    }
  }

  private tetraAux(n: Node): Desc {
    if (n.tetraAux !== undefined) return n.tetraAux;
    n.tetraAux = NONE;
    const spec = this.tetra.get(n.atom);
    if (!spec || n.dup || n.atom < 0 || !n.parent || n.parent.virtual) return NONE;
    const ligs = this.children(n, null).filter((k) => k.slot !== -2);
    if (ligs.length === 3 && spec.nbrs.includes(-1)) ligs.push(this.lonePair(n));
    if (ligs.length !== 4) return NONE;
    const res = this.withHorizon(n, () => this.rank(ligs, n));
    if (!res) return NONE;
    const slots = res.order.map((k) => (k === n.parent ? k.atom : k.slot));
    if (!sameSet(slots, spec.nbrs)) return NONE;
    const isR = tetraCcwForOrder(spec, [slots[3], slots[0], slots[1], slots[2]]);
    n.tetraAux = res.pseudo ? (isR ? D_r : D_s) : isR ? D_R : D_S;
    return n.tetraAux;
  }

  /** In-digraph E/Z of the stereo double bond whose nearer atom is `n`. */
  private dbAux(n: Node): Desc {
    if (n.dbAux !== undefined) return n.dbAux;
    n.dbAux = NONE;
    const spec = this.dbOfAtom.get(n.atom);
    if (!spec || n.dup || n.atom < 0 || !n.parent || n.parent.virtual) return NONE;
    const b = this.mol.bonds[spec.bond];
    const partner = b.a === n.atom ? b.b : b.a;
    const far = this.kidsOf(n).find((k) => k.atom === partner && !k.dup);
    if (!far) return NONE; // n is the far atom (or the bond closes a ring): descriptor lives elsewhere
    const res = this.withHorizon(n, () => this.ezFromRoots(n, far, spec));
    n.dbAux = res;
    return res;
  }

  /** E/Z from two digraph roots (r1 carries atom a1 of the bond, r2 the other atom). */
  private ezFromRoots(r1: Node, r2: Node, spec: DbSpec): Desc {
    const b = this.mol.bonds[spec.bond];
    const ligs1 = this.children(r1, null).filter((k) => k.atom !== r2.atom && k.slot !== -2);
    const ligs2 = this.children(r2, null).filter((k) => k.atom !== r1.atom && k.slot !== -2);
    const best = (ligs: Node[], root: Node): { top: Node; pseudo: boolean } | null => {
      if (ligs.length === 1) return { top: ligs[0], pseudo: false };
      if (ligs.length !== 2) return null;
      const r = this.rank(ligs, root);
      return r ? { top: r.order[0], pseudo: r.pseudo } : null;
    };
    const h1 = best(ligs1, r1);
    const h2 = best(ligs2, r2);
    if (!h1 || !h2) return NONE;
    const atomOf = (k: Node, root: Node) => (k === root.parent ? k.atom : k.slot);
    const carrier1 = r1.atom === b.a ? spec.a : spec.b;
    const carrier2 = r2.atom === b.a ? spec.a : spec.b;
    let cis = spec.cis;
    if (atomOf(h1.top, r1) !== carrier1) cis = !cis;
    if (atomOf(h2.top, r2) !== carrier2) cis = !cis;
    if (h1.pseudo !== h2.pseudo) return cis ? D_SEQCIS : D_SEQTRANS;
    return cis ? D_Z : D_E;
  }

  // ─────────────────────────── public entry points ───────────────────────────

  labelCenter(spec: TetraSpec): CenterLabel | null {
    const c = spec.center;
    if (c < 0 || c >= this.n) return null;
    const real = spec.nbrs.filter((x) => x >= 0);
    if (new Set(real).size !== real.length) return null;
    const root = this.newRoot(c);
    const ligs: Node[] = [];
    let hLeft = this.hcount[c];
    for (const nb of spec.nbrs) {
      if (nb >= 0) {
        const e = this.adj[c].find((x) => x.nbr === nb);
        if (!e) return null;
        ligs.push(this.atomNode(nb, root, e.bond, 1));
      } else if (hLeft > 0) {
        hLeft--;
        ligs.push(this.hNode(root));
      } else {
        ligs.push(this.lonePair(root));
      }
    }
    root.kids = ligs;
    const r = this.rank(ligs, root);
    if (!r) return null;
    const slots = r.order.map((k) => k.slot);
    const isR = tetraCcwForOrder(spec, [slots[3], slots[0], slots[1], slots[2]]);
    if (r.pseudo) return isR ? 'r' : 's';
    return isR ? 'R' : 'S';
  }

  labelBond(spec: DbSpec): BondLabel | null {
    const b = this.mol.bonds[spec.bond];
    if (!b || b.order !== 2) return null;
    if (!this.adj[b.a].some((e) => e.nbr === spec.a) || !this.adj[b.b].some((e) => e.nbr === spec.b)) return null;
    // each end is ranked in its own digraph; the other double-bond atom counts as visited
    const end = (x: number, y: number): Node | null => {
      const root = this.newRoot(x, y);
      const ligs: Node[] = [];
      for (const e of this.adj[x]) if (e.nbr !== y) ligs.push(this.atomNode(e.nbr, root, e.bond, 1));
      for (let k = 0; k < this.hcount[x]; k++) ligs.push(this.hNode(root));
      root.kids = ligs;
      if (ligs.length === 1) return ligs[0];
      if (ligs.length !== 2) return null;
      const r = this.rank(ligs, root);
      return r ? r.order[0] : null;
    };
    const h1 = end(b.a, b.b);
    const h2 = end(b.b, b.a);
    if (!h1 || !h2) return null;
    let cis = spec.cis;
    if (h1.slot !== spec.a) cis = !cis;
    if (h2.slot !== spec.b) cis = !cis;
    return cis ? 'Z' : 'E';
  }
}

function sameSet(a: readonly number[], b: readonly number[]): boolean {
  if (a.length !== b.length) return false;
  const s = [...a].sort((p, q) => p - q);
  const t = [...b].sort((p, q) => p - q);
  return s.every((v, i) => v === t[i]);
}

/**
 * Assigns CIP descriptors from mol.tetra / mol.dbStereo (populate them first via parseSmiles
 * or perceiveStereo2D). Implicit hydrogens are handled; abbreviations must be expanded by the caller.
 * Never throws: units that cannot be analysed (or exceed the digraph size budget) stay unlabelled.
 */
export function assignCIP(mol: Mol): CIPResult {
  const result: CIPResult = { centers: new Map(), bonds: new Map() };
  if (!mol.tetra.length && !mol.dbStereo.length) return result;
  let ctx: CIPContext;
  try {
    ctx = new CIPContext(mol);
  } catch {
    return result;
  }
  for (const t of ctx.tetra.values()) {
    try {
      const l = ctx.labelCenter(t);
      if (l) result.centers.set(t.center, l);
    } catch {
      // budget exceeded or malformed spec: leave unlabelled
    }
  }
  for (const d of ctx.db.values()) {
    try {
      const l = ctx.labelBond(d);
      if (l) result.bonds.set(d.bond, l);
    } catch {
      // leave unlabelled
    }
  }
  return result;
}

/** @internal Averaged duplicate atomic numbers of mancude atoms (exposed for tests). */
export function _mancudeDuplicateZ(mol: Mol): (number | undefined)[] {
  return new CIPContext(mol).mancudeZ;
}
