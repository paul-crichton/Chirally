// Ring-system parents: fused templates, Hantzsch–Widman monocycles, carbocycles,
// von Baeyer bicycles and monospiro systems. Each parent comes with the list of admissible
// numberings; the namer picks the best one by the lowest-locant rules.
import { NGraph } from './graph';
import { compiledTemplates, locantValue } from './templates';
import { alkaneStem, multiplier } from './numerals';

export interface Loc {
  label: string;
  value: number;
}
export type Numbering = Map<number, Loc>;

export interface RingParent {
  kind: 'ring';
  sys: number;
  atoms: number[];
  atomSet: Set<number>;
  /**
   * 'mancude': fixed name; hydrogenation expressed by hydro prefixes / indicated hydrogen.
   * 'ene': stem + ane/ene/yne endings (cycloalkanes, von Baeyer, spiro, replacement names).
   * 'benzene': fully aromatic benzene ring.
   */
  style: 'mancude' | 'ene' | 'benzene';
  /** Name of the mancude parent (style 'mancude'/'benzene'). */
  name: string;
  /** Name used when the ring has no endocyclic double bond (HW saturated names); '' if none. */
  satName: string;
  /** Stem for style 'ene' (e.g. "cyclohex", "7-oxabicyclo[2.2.1]hept"). */
  stem: string;
  numberings: Numbering[];
  /** Homogeneous carbocyclic monocycle (benzene / cycloalkane): locant omission rules apply. */
  carbocycleMono: boolean;
  nRings: number;
  /** sp2-capable atoms of the mancude parent. */
  A: Set<number>;
  needsIH: boolean;
  validIH: Set<number>;
  /** Seniority descriptors (P-44.2). */
  hetCount: Map<string, number>;
}

const HETERO_ORDER = ['O', 'S', 'Se', 'Te', 'N', 'P', 'As', 'Sb', 'Bi', 'Si', 'Ge', 'Sn', 'Pb', 'B', 'Hg'];
const REPL_PREFIX: Record<string, string> = {
  O: 'oxa', S: 'thia', Se: 'selena', Te: 'tellura', N: 'aza', P: 'phospha', As: 'arsa', Sb: 'stiba',
  Bi: 'bisma', Si: 'sila', Ge: 'germa', Sn: 'stanna', Pb: 'plumba', B: 'bora', Hg: 'mercura',
};

export function heteroRank(el: string): number {
  const i = HETERO_ORDER.indexOf(el);
  return i < 0 ? 99 : i;
}

function mkLoc(label: string): Loc {
  return { label, value: locantValue(label) };
}

/** Can the atom take part in a ring double bond of the mancude parent? */
export function sp2Capable(g: NGraph, i: number): boolean {
  const el = g.el[i];
  if (el === 'C' || el === 'N' || el === 'P' || el === 'As' || el === 'Sb' || el === 'B' || el === 'Si' || el === 'Ge') return true;
  if ((el === 'O' || el === 'S' || el === 'Se' || el === 'Te') && g.charge[i] === 1) return true;
  return false;
}

/** Backtracking perfect-matching test on a small graph. */
export function hasPerfectMatching(verts: number[], adj: (v: number) => number[]): boolean {
  if (verts.length % 2) return false;
  const set = new Set(verts);
  const matched = new Set<number>();
  let steps = 0;
  const rec = (): boolean => {
    if (++steps > 100000) return false;
    let v = -1;
    for (const x of verts) if (!matched.has(x)) { v = x; break; }
    if (v < 0) return true;
    matched.add(v);
    for (const w of adj(v)) {
      if (!set.has(w) || matched.has(w)) continue;
      matched.add(w);
      if (rec()) return true;
      matched.delete(w);
    }
    matched.delete(v);
    return false;
  };
  return rec();
}

function sysAdjacency(g: NGraph, atoms: number[]): Map<number, number[]> {
  const set = new Set(atoms);
  const m = new Map<number, number[]>();
  for (const a of atoms) m.set(a, g.nb[a].filter((b) => set.has(b)));
  return m;
}

/** Ordered cycle of a monocyclic ring system. */
function cycleOrder(adj: Map<number, number[]>, atoms: number[]): number[] {
  const out = [atoms[0]];
  let prev = -1;
  let cur = atoms[0];
  while (out.length < atoms.length) {
    const nx = adj.get(cur)!.find((x) => x !== prev && (out.length < 2 || x !== out[out.length - 2]) && !out.includes(x));
    if (nx === undefined) break;
    prev = cur;
    cur = nx;
    out.push(cur);
  }
  return out;
}

/** All 2n numberings of a ring given as an ordered cycle. */
function cycleNumberings(cyc: number[]): Numbering[] {
  const n = cyc.length;
  const out: Numbering[] = [];
  for (let s = 0; s < n; s++) {
    for (const dir of [1, -1]) {
      const m: Numbering = new Map();
      for (let k = 0; k < n; k++) m.set(cyc[(s + dir * k + n * 2) % n], mkLoc(String(k + 1)));
      out.push(m);
    }
  }
  return out;
}

function cmpArr(a: number[], b: number[]): number {
  for (let i = 0; i < Math.min(a.length, b.length); i++) if (a[i] !== b[i]) return a[i] - b[i];
  return a.length - b.length;
}

/** Keeps the numberings that give heteroatoms the lowest locants (as a set, then by seniority O > S > … > N …). */
function filterHetero(g: NGraph, atoms: number[], nums: Numbering[]): Numbering[] {
  const het = atoms.filter((a) => g.el[a] !== 'C');
  if (!het.length) return nums;
  const key = (m: Numbering): number[] => {
    const all = het.map((a) => m.get(a)!.value).sort((p, q) => p - q);
    const bySen: number[] = [];
    for (const el of HETERO_ORDER) {
      const l = het.filter((a) => g.el[a] === el).map((a) => m.get(a)!.value).sort((p, q) => p - q);
      bySen.push(...l);
    }
    return [...all, ...bySen];
  };
  let best: number[] | null = null;
  let out: Numbering[] = [];
  for (const m of nums) {
    const k = key(m);
    const c = best ? cmpArr(k, best) : -1;
    if (c < 0) {
      best = k;
      out = [m];
    } else if (c === 0) out.push(m);
  }
  return out;
}

/** "1,4-dioxa", "7-oxa", "3-oxa-8-aza": replacement prefixes for a fixed numbering. */
function replacementPrefixes(g: NGraph, atoms: number[], m: Numbering, omitSingleLocant: boolean): string {
  const het = atoms.filter((a) => g.el[a] !== 'C');
  if (!het.length) return '';
  const parts: string[] = [];
  for (const el of HETERO_ORDER) {
    const locs = het.filter((a) => g.el[a] === el).map((a) => m.get(a)!).sort((p, q) => p.value - q.value);
    if (!locs.length) continue;
    const pre = multiplier(locs.length) + REPL_PREFIX[el];
    parts.push((omitSingleLocant && het.length === 1 ? '' : locs.map((l) => l.label).join(',') + '-') + pre);
  }
  return parts.join('-');
}

// ───────────── Hantzsch–Widman names ─────────────

const HW_RETAINED: Record<string, string> = {
  azole: 'pyrrole', '1,3-diazole': 'imidazole', '1,2-diazole': 'pyrazole', oxole: 'furan', thiole: 'thiophene',
  selenole: 'selenophene', azine: 'pyridine', '1,2-diazine': 'pyridazine', '1,3-diazine': 'pyrimidine',
  '1,4-diazine': 'pyrazine', oxine: 'pyran', thiine: 'thiopyran',
  azolidine: 'pyrrolidine', '1,3-diazolidine': 'imidazolidine', '1,2-diazolidine': 'pyrazolidine',
  azinane: 'piperidine', '1,4-diazinane': 'piperazine', '1,4-oxazinane': 'morpholine', '1,4-thiazinane': 'thiomorpholine',
};

const GROUP_6A = new Set(['O', 'S', 'Se', 'Te', 'Bi', 'Hg']);
const GROUP_6B = new Set(['N', 'Si', 'Ge', 'Sn', 'Pb']);

function hwStem(size: number, saturated: boolean, hasN: boolean, last: string): string | null {
  switch (size) {
    case 3: return saturated ? (hasN ? 'iridine' : 'irane') : hasN ? 'irine' : 'irene';
    case 4: return saturated ? (hasN ? 'etidine' : 'etane') : 'ete';
    case 5: return saturated ? (hasN ? 'olidine' : 'olane') : 'ole';
    case 6:
      if (GROUP_6A.has(last)) return saturated ? 'ane' : 'ine';
      if (GROUP_6B.has(last)) return saturated ? 'inane' : 'ine';
      return saturated ? 'inane' : 'inine';
    case 7: return saturated ? 'epane' : 'epine';
    case 8: return saturated ? 'ocane' : 'ocine';
    case 9: return saturated ? 'onane' : 'onine';
    case 10: return saturated ? 'ecane' : 'ecine';
  }
  return null;
}

function hwName(g: NGraph, atoms: number[], m: Numbering, saturated: boolean): string | null {
  const het = atoms.filter((a) => g.el[a] !== 'C');
  const size = atoms.length;
  const els = HETERO_ORDER.filter((el) => het.some((a) => g.el[a] === el));
  if (els.length !== new Set(het.map((a) => g.el[a])).size) return null; // unknown heteroatom
  let pre = '';
  const locs: string[] = [];
  // joins name parts eliding a final "a" before a vowel ("oxa"+"aza" → "oxaza", "tetra"+"aza" → "tetraza")
  const join = (a: string, b: string) => (a.endsWith('a') && /^[aeiou]/.test(b) ? a.slice(0, -1) + b : a + b);
  for (const el of els) {
    const l = het.filter((a) => g.el[a] === el).map((a) => m.get(a)!).sort((p, q) => p.value - q.value);
    locs.push(...l.map((x) => x.label));
    pre = join(pre, join(multiplier(l.length), REPL_PREFIX[el]));
  }
  const last = els[els.length - 1];
  const stem = hwStem(size, saturated, het.some((a) => g.el[a] === 'N'), last);
  if (!stem) return null;
  // elide the final "a" of the prefix block before a vowel
  const name = join(pre, stem);
  const sameEl = els.length === 1;
  const omitLoc = het.length === 1 || (sameEl && het.length === size - 1);
  const full = omitLoc ? name : locs.join(',') + '-' + name;
  return HW_RETAINED[full] ?? full;
}

// ───────────── template isomorphism ─────────────

function templateIsomorphisms(g: NGraph, atoms: number[], adj: Map<number, number[]>, tpl: ReturnType<typeof compiledTemplates>[number]): Map<number, number>[] {
  const nT = tpl.labels.length;
  if (nT !== atoms.length) return [];
  // template BFS order
  const order: number[] = [0];
  const seen = new Array(nT).fill(false);
  seen[0] = true;
  for (let k = 0; k < order.length; k++) for (const w of tpl.adj[order[k]]) if (!seen[w]) { seen[w] = true; order.push(w); }
  if (order.length !== nT) return [];
  const tMap = new Array(nT).fill(-1); // template → molecule
  const used = new Set<number>();
  const results: Map<number, number>[] = [];
  const elOk = (a: number, t: number) => g.el[a] === tpl.els[t] && adj.get(a)!.length === tpl.adj[t].length;
  const rec = (k: number) => {
    if (results.length > 500) return;
    if (k === nT) {
      const m = new Map<number, number>();
      for (let t = 0; t < nT; t++) m.set(tMap[t], t);
      results.push(m);
      return;
    }
    const t = order[k];
    let cands: number[];
    if (k === 0) cands = atoms;
    else {
      const tp = tpl.adj[t].find((x) => tMap[x] >= 0)!;
      cands = adj.get(tMap[tp])!;
    }
    for (const a of cands) {
      if (used.has(a) || !elOk(a, t)) continue;
      let ok = true;
      for (const tn of tpl.adj[t]) {
        if (tMap[tn] >= 0 && !adj.get(a)!.includes(tMap[tn])) { ok = false; break; }
      }
      if (!ok) continue;
      tMap[t] = a;
      used.add(a);
      rec(k + 1);
      used.delete(a);
      tMap[t] = -1;
    }
  };
  rec(0);
  return results;
}

// ───────────── von Baeyer / spiro ─────────────

function vonBaeyer(g: NGraph, atoms: number[], adj: Map<number, number[]>): { stem: string; numberings: Numbering[]; fused: boolean } | null {
  const bh = atoms.filter((a) => adj.get(a)!.length === 3);
  if (bh.length !== 2 || atoms.some((a) => adj.get(a)!.length > 3)) return null;
  const [b1, b2] = bh;
  // bridges from b1 to b2 (atom lists, excluding bridgeheads, ordered from b1)
  const bridges: number[][] = [];
  for (const start of adj.get(b1)!) {
    const path: number[] = [];
    let prev = b1, cur = start;
    while (cur !== b2) {
      path.push(cur);
      const nx = adj.get(cur)!.find((x) => x !== prev);
      if (nx === undefined) return null;
      prev = cur;
      cur = nx;
      if (path.length > atoms.length) return null;
    }
    bridges.push(path);
  }
  if (bridges.length !== 3) return null;
  const sorted = [...bridges].sort((p, q) => q.length - p.length);
  const [a, b, c] = sorted.map((x) => x.length);
  const total = a + b + c + 2;
  const nums: Numbering[] = [];
  const perms = (arr: number[][]): number[][][] => {
    if (arr.length <= 1) return [arr];
    const out: number[][][] = [];
    arr.forEach((x, i) => {
      for (const rest of perms(arr.filter((_, j) => j !== i))) out.push([x, ...rest]);
    });
    return out;
  };
  for (const p of perms(bridges)) {
    if (!(p[0].length >= p[1].length && p[1].length >= p[2].length)) continue;
    for (const first of [b1, b2]) {
      const second = first === b1 ? b2 : b1;
      const fromFirst = (br: number[]) => (first === b1 ? br : [...br].reverse());
      const m: Numbering = new Map();
      let k = 1;
      m.set(first, mkLoc(String(k++)));
      for (const x of fromFirst(p[0])) m.set(x, mkLoc(String(k++)));
      m.set(second, mkLoc(String(k++)));
      for (const x of [...fromFirst(p[1])].reverse()) m.set(x, mkLoc(String(k++)));
      for (const x of fromFirst(p[2])) m.set(x, mkLoc(String(k++)));
      nums.push(m);
    }
  }
  const filtered = filterHetero(g, atoms, nums);
  const pre = replacementPrefixes(g, atoms, filtered[0], false);
  return { stem: pre + `bicyclo[${a}.${b}.${c}]` + alkaneStem(total), numberings: filtered, fused: c === 0 };
}

function spiro(g: NGraph, atoms: number[], adj: Map<number, number[]>): { stem: string; numberings: Numbering[] } | null {
  const sp = atoms.filter((a) => adj.get(a)!.length === 4);
  if (sp.length !== 1 || atoms.some((a) => adj.get(a)!.length !== 2 && a !== sp[0])) return null;
  const s = sp[0];
  // walk the two rings
  const rings: number[][] = [];
  const nbs = [...adj.get(s)!];
  const used = new Set<number>();
  for (const st of nbs) {
    if (used.has(st)) continue;
    const ring: number[] = [];
    let prev = s, cur = st;
    while (cur !== s) {
      ring.push(cur);
      used.add(cur);
      const nx = adj.get(cur)!.find((x) => x !== prev);
      if (nx === undefined) return null;
      prev = cur;
      cur = nx;
    }
    rings.push(ring);
  }
  if (rings.length !== 2) return null;
  rings.sort((p, q) => p.length - q.length);
  const [x, y] = [rings[0].length, rings[1].length];
  const nums: Numbering[] = [];
  const orders = x === y ? [[0, 1], [1, 0]] : [[0, 1]];
  for (const [i1, i2] of orders) {
    for (const d1 of [false, true]) {
      for (const d2 of [false, true]) {
        const r1 = d1 ? [...rings[i1]].reverse() : rings[i1];
        const r2 = d2 ? [...rings[i2]].reverse() : rings[i2];
        const m: Numbering = new Map();
        let k = 1;
        for (const a of r1) m.set(a, mkLoc(String(k++)));
        m.set(s, mkLoc(String(k++)));
        for (const a of r2) m.set(a, mkLoc(String(k++)));
        nums.push(m);
      }
    }
  }
  const filtered = filterHetero(g, atoms, nums);
  const pre = replacementPrefixes(g, atoms, filtered[0], false);
  return { stem: pre + `spiro[${x}.${y}]` + alkaneStem(x + y + 1), numberings: filtered };
}

// ───────────── main entry ─────────────

export interface RingParentResult {
  parent: RingParent | null;
  error?: string;
}

export function buildRingParent(g: NGraph, sys: number): RingParentResult {
  const atoms = g.sysAtoms[sys];
  const atomSet = new Set(atoms);
  const adj = sysAdjacency(g, atoms);
  let nEdges = 0;
  for (const a of atoms) nEdges += adj.get(a)!.length;
  nEdges /= 2;
  const nRings = nEdges - atoms.length + 1;
  const hetCount = new Map<string, number>();
  for (const a of atoms) if (g.el[a] !== 'C') hetCount.set(g.el[a], (hetCount.get(g.el[a]) ?? 0) + 1);
  const D = atoms.filter((a) => g.nb[a].some((b) => atomSet.has(b) && g.order(a, b) === 2));

  const base = (): RingParent => ({
    kind: 'ring', sys, atoms, atomSet, style: 'mancude', name: '', satName: '', stem: '', numberings: [],
    carbocycleMono: false, nRings, A: new Set(), needsIH: false, validIH: new Set(), hetCount,
  });
  const fillMancude = (rp: RingParent) => {
    const A = atoms.filter((a) => sp2Capable(g, a));
    const Aset = new Set(A);
    const nonIso = A.filter((a) => adj.get(a)!.some((b) => Aset.has(b)));
    rp.A = Aset;
    const f = (v: number) => adj.get(v)!;
    if (!hasPerfectMatching(nonIso, f)) {
      rp.needsIH = true;
      for (const p of nonIso) if (hasPerfectMatching(nonIso.filter((x) => x !== p), f)) rp.validIH.add(p);
    }
  };

  if (nRings === 1) {
    const cyc = cycleOrder(adj, atoms);
    if (cyc.length !== atoms.length) return { parent: null, error: 'ring perception failed' };
    const nums = cycleNumberings(cyc);
    const rp = base();
    if (!hetCount.size) {
      rp.carbocycleMono = true;
      rp.numberings = nums;
      if (atoms.length === 6 && D.length === 6) {
        rp.style = 'benzene';
        rp.name = 'benzene';
      } else {
        rp.style = 'ene';
        rp.stem = 'cyclo' + alkaneStem(atoms.length);
      }
      return { parent: rp };
    }
    const filtered = filterHetero(g, atoms, nums);
    if (atoms.length <= 10) {
      const mName = hwName(g, atoms, filtered[0], false);
      const sName = hwName(g, atoms, filtered[0], true);
      if (!mName || !sName) return { parent: null, error: 'unsupported heteroatom in ring' };
      rp.style = 'mancude';
      rp.name = mName;
      rp.satName = sName;
      rp.numberings = filtered;
      fillMancude(rp);
      return { parent: rp };
    }
    // large heteromonocycles: replacement nomenclature ("1,4,7,10-tetraoxacyclododecane")
    if ([...hetCount.keys()].some((el) => !REPL_PREFIX[el])) return { parent: null, error: 'unsupported heteroatom in ring' };
    rp.style = 'ene';
    rp.numberings = filtered;
    rp.stem = replacementPrefixes(g, atoms, filtered[0], true) + 'cyclo' + alkaneStem(atoms.length);
    return { parent: rp };
  }

  // fused / bridged templates
  for (const tpl of compiledTemplates()) {
    if (tpl.labels.length !== atoms.length || tpl.nBonds !== nEdges) continue;
    const isos = templateIsomorphisms(g, atoms, adj, tpl);
    if (!isos.length) continue;
    const rp = base();
    rp.numberings = isos.map((iso) => {
      const m: Numbering = new Map();
      for (const [a, t] of iso) m.set(a, mkLoc(tpl.labels[t]));
      return m;
    });
    if (tpl.def.saturated) {
      if (D.length) return { parent: null, error: 'unsaturated ' + tpl.def.name + ' not supported' };
      rp.style = 'ene';
      rp.stem = tpl.def.name.replace(/ane$/, '');
      return { parent: rp };
    }
    rp.style = 'mancude';
    rp.name = tpl.def.name;
    fillMancude(rp);
    return { parent: rp };
  }

  if (nRings === 2) {
    const sp = spiro(g, atoms, adj);
    if (sp) {
      if ([...hetCount.keys()].some((el) => !REPL_PREFIX[el])) return { parent: null, error: 'unsupported heteroatom in ring' };
      const rp = base();
      rp.style = 'ene';
      rp.stem = sp.stem;
      rp.numberings = sp.numberings;
      return { parent: rp };
    }
    const vb = vonBaeyer(g, atoms, adj);
    if (vb) {
      // ortho-fused bicycles with two rings of ≥5 members require fusion nomenclature (P-52.2.4.3)
      const sizes = vb.stem.match(/\[(\d+)\.(\d+)\.(\d+)\]/);
      if (vb.fused && sizes && +sizes[2] + 2 >= 5) return { parent: null, error: 'fused ring system not supported (no template)' };
      if ([...hetCount.keys()].some((el) => !REPL_PREFIX[el])) return { parent: null, error: 'unsupported heteroatom in ring' };
      const rp = base();
      rp.style = 'ene';
      rp.stem = vb.stem;
      rp.numberings = vb.numberings;
      return { parent: rp };
    }
  }
  return { parent: null, error: 'ring system not supported' };
}
