// Minimal SMARTS engine used for descriptor atom typing (Wildman–Crippen, RDKit HBD/HBA,
// rotatable bonds) and functional-group perception.
//
// Supported: bracket atoms with element symbols (aliphatic upper-case / aromatic lower-case),
// #n, *, A, a, H<n> (total H), D<n> (explicit degree), X<n> (total connections), v<n> (valence),
// R / R<n> (SSSR ring count), r / r<n> (SSSR ring size), x<n> (ring bonds), charges (+, -, +2, --, +0),
// isotopes, recursive $(…), logical ! & , ; ; bonds - = # : ~ @ with the same logic operators;
// branches and ring closures. Chirality, atom maps and hybridisation (^) are parsed and ignored.
//
// Matching runs on a MatchGraph: a flat, aromaticity-perceived view of a Mol in which hydrogens
// are either folded into H counts ('none', like an RDKit molecule without explicit H) or all
// present as atoms ('all', like RDKit after AddHs).
import { element } from '../elements';
import { Mol } from '../mol';
import { implicitH, bondValence } from '../valence';
import { perceiveRings, RingInfo } from '../rings';
import { perceiveAromaticity } from '../aromaticity';

// ───────────────────────────── Match graph ─────────────────────────────

export interface MatchGraph {
  n: number;
  /** Atomic number (0 for pseudo atoms). */
  z: number[];
  arom: boolean[];
  charge: number[];
  /** Mass number, 0 = unspecified. */
  isotope: number[];
  /** Total attached hydrogens (implicit + explicit H neighbours). */
  hTotal: number[];
  /** Hydrogens not represented as graph atoms. */
  hImplicit: number[];
  /** Number of graph neighbours (explicit degree). */
  degree: number[];
  /** Kekulé bond-order sum + hydrogens (SMARTS v). */
  valence: number[];
  /** Number of SSSR rings containing the atom. */
  ringCount: number[];
  /** Sizes of the SSSR rings containing the atom. */
  ringSizes: number[][];
  nbr: number[][];
  /** Bond index (into the bond arrays below) parallel to nbr. */
  nbrBond: number[][];
  bondA: number[];
  bondB: number[];
  bondOrder: number[];
  bondArom: boolean[];
  bondRing: boolean[];
  /** Index of the originating Mol atom; for added hydrogens the heavy atom they hang on. */
  src: number[];
  /** True for hydrogens created from implicit H counts. */
  addedH: boolean[];
  /** Graph index of each Mol atom (-1 when the atom was folded into an H count). */
  fromMol: number[];
  rings: RingInfo;
  /** Memo for recursive SMARTS: pattern id → per-atom result (0 unknown, 1 true, 2 false). */
  recCache: Map<number, Uint8Array>;
}

/** True for an "ordinary" explicit hydrogen that RDKit would remove (neutral, no isotope, one single bond). */
function foldableH(mol: Mol, i: number): boolean {
  const a = mol.atoms[i];
  if (a.el !== 'H' || a.isotope || a.charge || mol.degree(i) !== 1) return false;
  const j = mol.neighbors(i)[0];
  return mol.atoms[j].el !== 'H' && mol.bonds[mol.adj[i][0]].order === 1;
}

/**
 * Builds a MatchGraph. Aromaticity is perceived with perceiveAromaticity (the Mol should be in a
 * Kekulé form or use 1.5 bond orders). Bonds of order 0 / H-bonds are ignored.
 */
export function buildMatchGraph(mol: Mol, hMode: 'none' | 'all', rings?: RingInfo): MatchGraph {
  const ri = rings ?? perceiveRings(mol);
  const aro = perceiveAromaticity(mol, ri);
  const nA = mol.atoms.length;
  const fromMol = new Array(nA).fill(-1);
  const g: MatchGraph = {
    n: 0, z: [], arom: [], charge: [], isotope: [], hTotal: [], hImplicit: [], degree: [], valence: [],
    ringCount: [], ringSizes: [], nbr: [], nbrBond: [], bondA: [], bondB: [], bondOrder: [], bondArom: [],
    bondRing: [], src: [], addedH: [], fromMol, rings: ri, recCache: new Map(),
  };
  const folded = new Array(nA).fill(false);
  if (hMode === 'none') for (let i = 0; i < nA; i++) folded[i] = foldableH(mol, i);

  const addAtom = (src: number, added: boolean): number => {
    const k = g.n++;
    const a = mol.atoms[src];
    const z = added ? 1 : element(a.el)?.z ?? 0;
    g.z.push(z);
    g.arom.push(!added && aro.atoms[src]);
    g.charge.push(added ? 0 : a.charge);
    g.isotope.push(added ? 0 : a.isotope ?? 0);
    g.hTotal.push(0);
    g.hImplicit.push(0);
    g.degree.push(0);
    g.valence.push(0);
    g.ringCount.push(added ? 0 : ri.atomRings[src].length);
    g.ringSizes.push(added ? [] : ri.atomRings[src].map((r) => ri.rings[r].length));
    g.nbr.push([]);
    g.nbrBond.push([]);
    g.src.push(src);
    g.addedH.push(added);
    return k;
  };
  const addBond = (a: number, b: number, order: number, arom: boolean, ring: boolean) => {
    const k = g.bondA.length;
    g.bondA.push(a);
    g.bondB.push(b);
    g.bondOrder.push(order);
    g.bondArom.push(arom);
    g.bondRing.push(ring);
    g.nbr[a].push(b);
    g.nbrBond[a].push(k);
    g.nbr[b].push(a);
    g.nbrBond[b].push(k);
    g.degree[a]++;
    g.degree[b]++;
  };

  for (let i = 0; i < nA; i++) if (!folded[i]) fromMol[i] = addAtom(i, false);
  mol.bonds.forEach((b, bi) => {
    if (b.order === 0 || b.style === 'hbond') return;
    const ga = fromMol[b.a], gb = fromMol[b.b];
    if (ga < 0 || gb < 0) return;
    const arom = aro.bonds[bi] || b.order === 1.5;
    addBond(ga, gb, b.order === 1.5 ? 1 : Math.round(b.order), arom, ri.bondInRing[bi]);
  });
  // Hydrogen bookkeeping and valence
  for (let i = 0; i < nA; i++) {
    const k = fromMol[i];
    if (k < 0) continue;
    const a = mol.atoms[i];
    const pseudo = !element(a.el) || !!a.abbrev;
    const ih = pseudo ? 0 : implicitH(mol, i);
    let explicitH = 0;
    let foldedH = 0;
    let bsum = 0;
    for (const bi of mol.adj[i]) {
      const b = mol.bonds[bi];
      bsum += bondValence(b.order, b.style);
      const j = mol.other(bi, i);
      if (mol.atoms[j].el === 'H') {
        explicitH++;
        if (folded[j]) foldedH++;
      }
    }
    g.hTotal[k] = ih + explicitH;
    g.valence[k] = Math.round(bsum + ih);
    if (hMode === 'all') {
      for (let h = 0; h < ih; h++) {
        const hk = addAtom(i, true);
        addBond(k, hk, 1, false, false);
        g.hTotal[hk] = 0;
        g.valence[hk] = 1;
      }
      g.hImplicit[k] = 0;
    } else {
      g.hImplicit[k] = ih + foldedH;
    }
  }
  return g;
}

// ───────────────────────────── SMARTS compilation ─────────────────────────────

type AtomTest = (g: MatchGraph, i: number) => boolean;
type BondTest = (g: MatchGraph, b: number) => boolean;

export interface SmartsPattern {
  source: string;
  atoms: AtomTest[];
  /** Tree parent of each pattern atom (-1 for the first atom). */
  parent: number[];
  parentBond: BondTest[];
  /** Ring-closure constraints checked once both ends are mapped (stored on the later atom). */
  closures: { other: number; test: BondTest }[][];
}

let recId = 0;
const cache = new Map<string, SmartsPattern>();

const anyBondDefault: BondTest = (g, b) => g.bondArom[b] || g.bondOrder[b] === 1;

class Parser {
  pos = 0;
  constructor(public s: string) {}
  peek(o = 0): string {
    return this.s[this.pos + o] ?? '';
  }
  fail(msg: string): never {
    throw new Error(`SMARTS ${msg} at ${this.pos} in "${this.s}"`);
  }
  num(): number | null {
    const m = /^\d+/.exec(this.s.slice(this.pos));
    if (!m) return null;
    this.pos += m[0].length;
    return +m[0];
  }
}

/** Generic precedence parser shared by atom and bond expressions: ';' < ',' < '&'/implicit < '!'. */
function parseLogic(p: Parser, end: (c: string) => boolean, prim: () => (g: MatchGraph, i: number) => boolean): (g: MatchGraph, i: number) => boolean {
  const unary = (): ((g: MatchGraph, i: number) => boolean) => {
    if (p.peek() === '!') {
      p.pos++;
      const f = unary();
      return (g, i) => !f(g, i);
    }
    return prim();
  };
  const high = () => {
    const parts = [unary()];
    for (;;) {
      const c = p.peek();
      if (c === '&') {
        p.pos++;
        parts.push(unary());
      } else if (c && c !== ',' && c !== ';' && !end(c)) parts.push(unary());
      else break;
    }
    return parts.length === 1 ? parts[0] : (g: MatchGraph, i: number) => parts.every((f) => f(g, i));
  };
  const or = () => {
    const parts = [high()];
    while (p.peek() === ',') {
      p.pos++;
      parts.push(high());
    }
    return parts.length === 1 ? parts[0] : (g: MatchGraph, i: number) => parts.some((f) => f(g, i));
  };
  const parts = [or()];
  while (p.peek() === ';') {
    p.pos++;
    parts.push(or());
  }
  return parts.length === 1 ? parts[0] : (g, i) => parts.every((f) => f(g, i));
}

const AROMATIC_SYMBOLS = ['se', 'as', 'te', 'c', 'n', 'o', 's', 'p', 'b'];

function elementTest(z: number, aromatic: boolean | null): AtomTest {
  if (aromatic === null) return (g, i) => g.z[i] === z;
  return aromatic ? (g, i) => g.z[i] === z && g.arom[i] : (g, i) => g.z[i] === z && !g.arom[i];
}

/** Parses one atom primitive inside brackets. */
function atomPrimitive(p: Parser): AtomTest {
  const c = p.peek();
  if (c === '$' && p.peek(1) === '(') {
    // recursive SMARTS: find matching parenthesis
    let depth = 0;
    const start = p.pos + 2;
    let k = p.pos + 1;
    for (; k < p.s.length; k++) {
      if (p.s[k] === '(') depth++;
      else if (p.s[k] === ')') {
        depth--;
        if (depth === 0) break;
      }
    }
    if (depth !== 0) p.fail('unbalanced $(');
    const inner = compileSmarts(p.s.slice(start, k));
    p.pos = k + 1;
    const id = ++recId;
    return (g, i) => {
      let memo = g.recCache.get(id);
      if (!memo) {
        memo = new Uint8Array(g.n);
        g.recCache.set(id, memo);
      }
      if (memo[i] === 0) memo[i] = matchAt(inner, g, i) ? 1 : 2;
      return memo[i] === 1;
    };
  }
  if (c === '#') {
    p.pos++;
    const z = p.num();
    if (z === null) p.fail('expected atomic number');
    return elementTest(z, null);
  }
  if (c >= '0' && c <= '9') {
    const iso = p.num()!;
    return (g, i) => g.isotope[i] === iso;
  }
  if (c === '*') {
    p.pos++;
    return () => true;
  }
  if (c === '+' || c === '-') {
    const sign = c === '+' ? 1 : -1;
    p.pos++;
    let n = p.num();
    if (n === null) {
      n = 1;
      while (p.peek() === c) {
        n++;
        p.pos++;
      }
    }
    const want = sign * n;
    return (g, i) => g.charge[i] === want;
  }
  if (c === '@') {
    // chirality is not used for matching
    p.pos++;
    if (p.peek() === '@') p.pos++;
    return () => true;
  }
  if (c === '^') {
    p.pos++;
    p.num();
    return () => true;
  }
  if (c === ':') {
    p.pos++;
    p.num();
    return () => true;
  }
  if (c >= 'A' && c <= 'Z') {
    const two = c + p.peek(1);
    if (p.peek(1) >= 'a' && p.peek(1) <= 'z' && element(two)) {
      p.pos += 2;
      return elementTest(element(two)!.z, false);
    }
    if ('HDXRA'.includes(c)) {
      p.pos++;
      const n = p.num();
      switch (c) {
        case 'H': { const v = n ?? 1; return (g, i) => g.hTotal[i] === v; }
        case 'D': { const v = n ?? 1; return (g, i) => g.degree[i] === v; }
        case 'X': { const v = n ?? 1; return (g, i) => g.degree[i] + g.hImplicit[i] === v; }
        case 'R': return n === null ? (g, i) => g.ringCount[i] > 0 : (g, i) => g.ringCount[i] === n;
        case 'A': return (g, i) => !g.arom[i] && g.z[i] > 0;
      }
    }
    const e = element(c);
    if (!e) p.fail(`unknown element ${c}`);
    p.pos++;
    return elementTest(e.z, false);
  }
  if (c >= 'a' && c <= 'z') {
    for (const sym of AROMATIC_SYMBOLS) {
      if (p.s.startsWith(sym, p.pos)) {
        p.pos += sym.length;
        const z = element(sym[0].toUpperCase() + sym.slice(1))!.z;
        return elementTest(z, true);
      }
    }
    p.pos++;
    if (c === 'a') return (g, i) => g.arom[i];
    const n = p.num();
    switch (c) {
      case 'v': { const v = n ?? 1; return (g, i) => g.valence[i] === v; }
      case 'r': return n === null ? (g, i) => g.ringCount[i] > 0 : (g, i) => g.ringSizes[i].includes(n);
      case 'x': {
        return (g, i) => {
          let k = 0;
          for (const b of g.nbrBond[i]) if (g.bondRing[b]) k++;
          return n === null ? k > 0 : k === n;
        };
      }
      case 'h': { const v = n ?? 1; return (g, i) => g.hImplicit[i] === v; }
    }
  }
  return p.fail(`unexpected '${c}'`);
}

function bondPrimitive(p: Parser): BondTest {
  const c = p.peek();
  p.pos++;
  switch (c) {
    case '-': case '/': case '\\': return (g, b) => !g.bondArom[b] && g.bondOrder[b] === 1;
    case '=': return (g, b) => !g.bondArom[b] && g.bondOrder[b] === 2;
    case '#': return (g, b) => !g.bondArom[b] && g.bondOrder[b] === 3;
    case ':': return (g, b) => g.bondArom[b];
    case '~': return () => true;
    case '@': return (g, b) => g.bondRing[b];
  }
  return p.fail(`unexpected bond '${c}'`);
}

const BOND_CHARS = '-=#:~@!&,;/\\';

/** Compiles (and caches) a SMARTS string. Throws on syntax errors. */
export function compileSmarts(s: string): SmartsPattern {
  const hit = cache.get(s);
  if (hit) return hit;
  const p = new Parser(s);
  const pat: SmartsPattern = { source: s, atoms: [], parent: [], parentBond: [], closures: [] };
  const stack: number[] = [];
  let prev = -1;
  let pendingBond: BondTest | null = null;
  const ringOpen = new Map<number, { atom: number; bond: BondTest | null }>();

  const addAtom = (test: AtomTest) => {
    const k = pat.atoms.length;
    pat.atoms.push(test);
    pat.parent.push(prev);
    pat.parentBond.push(pendingBond ?? anyBondDefault);
    pat.closures.push([]);
    if (prev < 0 && k > 0) p.fail('disconnected patterns are not supported');
    pendingBond = null;
    prev = k;
  };

  while (p.pos < s.length) {
    const c = p.peek();
    if (c === '(') {
      stack.push(prev);
      p.pos++;
    } else if (c === ')') {
      prev = stack.pop() ?? p.fail('unmatched )');
      p.pos++;
    } else if (BOND_CHARS.includes(c)) {
      pendingBond = parseLogic(p, (ch) => !BOND_CHARS.includes(ch), () => bondPrimitive(p));
    } else if ((c >= '0' && c <= '9') || c === '%') {
      let num: number;
      if (c === '%') {
        p.pos++;
        num = p.num() ?? p.fail('bad ring number');
      } else {
        num = +c;
        p.pos++;
      }
      const open = ringOpen.get(num);
      if (open) {
        pat.closures[prev].push({ other: open.atom, test: pendingBond ?? open.bond ?? anyBondDefault });
        ringOpen.delete(num);
      } else {
        ringOpen.set(num, { atom: prev, bond: pendingBond });
      }
      pendingBond = null;
    } else if (c === '[') {
      p.pos++;
      const test = parseLogic(p, (ch) => ch === ']', () => atomPrimitive(p));
      if (p.peek() !== ']') p.fail('expected ]');
      p.pos++;
      addAtom(test);
    } else {
      // unbracketed atom (organic subset, aromatic symbols, *, A, a)
      const two = s.slice(p.pos, p.pos + 2);
      if (two === 'Cl' || two === 'Br') {
        p.pos += 2;
        addAtom(elementTest(element(two)!.z, false));
      } else if ('BCNOPSFI'.includes(c)) {
        p.pos++;
        addAtom(elementTest(element(c)!.z, false));
      } else if ('bcnops'.includes(c)) {
        p.pos++;
        addAtom(elementTest(element(c.toUpperCase())!.z, true));
      } else if (c === '*') {
        p.pos++;
        addAtom(() => true);
      } else if (c === 'A') {
        p.pos++;
        addAtom((g, i) => !g.arom[i] && g.z[i] > 0);
      } else if (c === 'a') {
        p.pos++;
        addAtom((g, i) => g.arom[i]);
      } else p.fail(`unexpected '${c}'`);
    }
  }
  if (ringOpen.size) p.fail('unclosed ring');
  cache.set(s, pat);
  return pat;
}

// ───────────────────────────── Matching ─────────────────────────────

/**
 * Backtracking match with pattern atom 0 mapped to graph atom `root`.
 * `onMatch` receives each complete mapping; returning true stops the search.
 */
function search(pat: SmartsPattern, g: MatchGraph, root: number, onMatch: (map: number[]) => boolean): boolean {
  if (!pat.atoms[0](g, root)) return false;
  const n = pat.atoms.length;
  const map = new Array(n).fill(-1);
  const used = new Uint8Array(g.n);
  map[0] = root;
  used[root] = 1;
  const rec = (k: number): boolean => {
    if (k === n) return onMatch(map);
    const par = map[pat.parent[k]];
    const nb = g.nbr[par];
    const nbB = g.nbrBond[par];
    for (let t = 0; t < nb.length; t++) {
      const cand = nb[t];
      if (used[cand]) continue;
      if (!pat.parentBond[k](g, nbB[t])) continue;
      if (!pat.atoms[k](g, cand)) continue;
      let ok = true;
      for (const cl of pat.closures[k]) {
        const o = map[cl.other];
        const bi = bondIndex(g, cand, o);
        if (bi < 0 || !cl.test(g, bi)) {
          ok = false;
          break;
        }
      }
      if (!ok) continue;
      map[k] = cand;
      used[cand] = 1;
      if (rec(k + 1)) return true;
      used[cand] = 0;
      map[k] = -1;
    }
    return false;
  };
  // ring closures on atom 0 are impossible (it is always the first atom), so start at k = 1
  return rec(1);
}

function bondIndex(g: MatchGraph, a: number, b: number): number {
  const nb = g.nbr[a];
  for (let t = 0; t < nb.length; t++) if (nb[t] === b) return g.nbrBond[a][t];
  return -1;
}

/** True if the pattern matches with its first atom on graph atom i. */
export function matchAt(pat: SmartsPattern, g: MatchGraph, i: number): boolean {
  return search(pat, g, i, () => true);
}

/**
 * All matches as arrays of graph atom indices (pattern order). With `unique` (default true),
 * matches covering the same atom set are reported once (RDKit uniquify semantics).
 */
export function findMatches(pat: SmartsPattern, g: MatchGraph, unique = true, limit = 10000): number[][] {
  const out: number[][] = [];
  const seen = new Set<string>();
  for (let i = 0; i < g.n && out.length < limit; i++) {
    search(pat, g, i, (map) => {
      if (unique) {
        const key = [...map].sort((a, b) => a - b).join(',');
        if (seen.has(key)) return false;
        seen.add(key);
      }
      out.push([...map]);
      return out.length >= limit;
    });
  }
  return out;
}
