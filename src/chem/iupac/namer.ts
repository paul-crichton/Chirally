// Core substitutive naming engine (IUPAC 2013 recommendations, PubChem/LexiChem-compatible style).
//
// Overview
//  1. Atoms are classified: skeletal atoms (ring atoms and acyclic carbons) can belong to a parent
//     hydride; "group carbons" (–COOH, –CHO, –CN, –CONR2, –COOR, –COX) may only terminate chains
//     or be expressed as -carboxylic acid etc.; acyclic heteroatoms always belong to groups.
//  2. The principal characteristic group class is chosen by seniority (P-41).
//  3. Parent candidates (ring systems and acyclic chains) are enumerated and the senior one is
//     selected (P-44, 2013 rules: principal groups first, then rings > chains, then length …).
//  4. The parent is numbered by the lowest-locant rules (P-31.1.4) and everything attached to it is
//     named recursively as substituent prefixes.
import { NGraph, HALOGENS } from './graph';
import { buildRingParent, RingParent, Numbering, Loc } from './ringnames';
import { alkaneStem, multiplier, multiplierComplex } from './numerals';
import { SubName, PrefixItem, formatPrefixes, elide, enclose, alphaKey, startsAmbiguous } from './assemble';

export class NamingError extends Error {}

/** Principal characteristic group classes – lower value = more senior (P-41). */
export const CLS = {
  ANION_COO: 10,
  ANION_SO3: 11,
  ANION_INORG: 12,
  ANION_O: 13,
  CATION: 20,
  ACID: 30,
  CARBAMIC: 31,
  SULFONIC: 32,
  SULFINIC: 33,
  ESTER: 40,
  SULFONATE_ESTER: 41,
  CARBAMATE: 42,
  ESTER_INORG: 43,
  ACYL_HALIDE: 50,
  SULFONYL_HALIDE: 51,
  AMIDE: 60,
  UREA: 61,
  SULFONAMIDE: 62,
  NITRILE: 70,
  ALDEHYDE: 80,
  KETONE: 90,
  THIONE: 91,
  ALCOHOL: 100,
  THIOL: 101,
  AMINE: 110,
  IMINE: 120,
};

type GCType = 'acid' | 'carboxylate' | 'ester' | 'amide' | 'acylhalide' | 'aldehyde' | 'nitrile' | 'carbonic' | 'bad';

type SuffixKind =
  | 'acid' | 'carboxylate' | 'ester' | 'acylhalide' | 'amide' | 'nitrile' | 'aldehyde'
  | 'one' | 'thione' | 'ol' | 'olate' | 'thiol' | 'amine' | 'aminium' | 'imine'
  | 'sulfonic' | 'sulfonate' | 'sulfonamide' | 'sulfonylhalide' | 'sulfinic' | 'ium'
  | 'yl' | 'ylidene' | 'ylidyne' | 'oyl' | 'carbonyl';

const GC_KINDS = new Set<SuffixKind>(['acid', 'carboxylate', 'ester', 'acylhalide', 'amide', 'nitrile', 'aldehyde']);

interface GroupInfo {
  cls: number;
  kind: SuffixKind;
  gc: number;
  /** Nitrogen whose other substituents become N-prefixes (-1 if none). */
  n: number;
  consumed: number[];
  /** Ester: oxygen and alkyl attachment atom. */
  esterO: number;
  esterR: number;
  halide: string;
}

const HALO_PREFIX: Record<string, string> = { F: 'fluoro', Cl: 'chloro', Br: 'bromo', I: 'iodo' };
const HALIDE: Record<string, string> = { F: 'fluoride', Cl: 'chloride', Br: 'bromide', I: 'iodide' };

interface Parent {
  kind: 'chain' | 'ring';
  atoms: number[];
  ring?: RingParent;
}

interface FV {
  atom: number;
  from: number;
  order: number;
  kind: 'yl' | 'oyl' | 'carbonyl';
  /** For 'carbonyl': the C=O carbon between ring atom and the attachment point. */
  carbonylC?: number;
}

interface BuildOpts {
  cls: number | null;
  /** Chain parents: principal group carbons are included at chain ends ("-oic acid"). */
  inclMode?: boolean;
  fv?: FV;
  /** Substituent naming context: 'oxy' disables benzyl/tert-butyl (LexiChem style). */
  ctx?: string;
}

export interface BuildResult {
  text: string;
  /** Parent + suffix part only (for functional class names). */
  suffixLocs: number[];
  unsatLocs: number[];
  nPrefixes: number;
  prefixLocs: number[];
  alpha: string[];
  locants: Map<number, string>;
  /** Ester alkyl groups with the locant of their acid group. */
  esterAlkyls: { sub: SubName; loc: string }[];
  nSuffix: number;
  hasPrefixes: boolean;
  /** Chain, free valence at C1, order 1 (oxy contraction). */
  oxyContract: boolean;
  phenyl: boolean;
  phenylLike: boolean;
}

interface Item {
  atom: number;
  sub: SubName;
}

interface SuffixItem {
  atom: number;
  kind: SuffixKind;
  halide?: string;
}

export interface CIPInfo {
  centers: Map<number, string>;
  /** key "i,j" (i<j) → E/Z */
  bonds: Map<string, string>;
}

function cmpNum(a: number[], b: number[]): number {
  for (let i = 0; i < Math.min(a.length, b.length); i++) if (a[i] !== b[i]) return a[i] - b[i];
  return a.length - b.length;
}

function cmpStr(a: string[], b: string[]): number {
  for (let i = 0; i < Math.min(a.length, b.length); i++) if (a[i] !== b[i]) return a[i] < b[i] ? -1 : 1;
  return a.length - b.length;
}

export class Namer {
  g: NGraph;
  warnings: string[] = [];
  cip: CIPInfo | null;
  stereoUsed = new Set<string>();
  gcType: (GCType | null)[];
  skel: boolean[];
  private ringCache = new Map<number, RingParent>();
  private subCache = new Map<string, SubName>();
  private deadline: number;

  constructor(g: NGraph, cip: CIPInfo | null = null) {
    this.g = g;
    this.cip = cip;
    this.deadline = Date.now() + 2000;
    this.gcType = new Array(g.n).fill(null);
    this.skel = new Array(g.n).fill(false);
    for (let i = 0; i < g.n; i++) {
      if (g.inRing[i]) this.skel[i] = true;
      else if (g.el[i] === 'C') {
        const t = this.classifyCarbon(i);
        this.gcType[i] = t;
        if (!t) this.skel[i] = true;
      }
    }
  }

  private check(): void {
    if (Date.now() > this.deadline) throw new NamingError('naming took too long');
  }

  // ───────────────────────── atom classification ─────────────────────────

  /** Classifies an acyclic carbon: group carbon type or null for skeletal carbons. */
  private classifyCarbon(c: number): GCType | null {
    const g = this.g;
    let dO = -1, dS = -1, dN = -1, tN = -1;
    for (const j of g.nb[c]) {
      const o = g.order(c, j);
      if (o === 2 && g.el[j] === 'O' && g.degree(j) === 1) dO = j;
      else if (o === 2 && g.el[j] === 'S' && g.degree(j) === 1) dS = j;
      else if (o === 2 && g.el[j] === 'N') dN = j;
      else if (o === 3 && g.el[j] === 'N' && g.degree(j) === 1) tN = j;
    }
    const het = g.nb[c].filter((j) => g.el[j] !== 'C' && !g.inRing[j] && g.order(c, j) === 1);
    if (tN >= 0) {
      if (het.length) return 'bad'; // cyanate / thiocyanate / cyanamide carbon (named from the heteroatom)
      return 'nitrile';
    }
    if (dO >= 0) {
      if (het.length === 0) return g.h[c] >= 1 ? 'aldehyde' : null;
      if (het.length === 2) return 'carbonic';
      const x = het[0];
      if (g.el[x] === 'O') {
        if (g.charge[x] === -1) return 'carboxylate';
        const o = g.others(x, c);
        if (o.length === 0) return 'acid';
        if (o.length === 1 && (g.el[o[0]] === 'C' || g.inRing[o[0]]) && g.charge[x] === 0) {
          // anhydrides (O–C(=O) on both sides) are not supported
          if (!g.inRing[o[0]] && this.isCarbonylC(o[0]) && g.el[o[0]] === 'C') return 'bad';
          return 'ester';
        }
        return 'bad';
      }
      if (g.el[x] === 'N') return g.charge[x] === 0 ? 'amide' : 'bad';
      if (HALOGENS.has(g.el[x])) return 'acylhalide';
      return 'bad';
    }
    if ((dS >= 0 || dN >= 0) && het.length) return 'bad';
    return null;
  }

  private isCarbonylC(c: number): boolean {
    return this.g.terminalDouble(c, 'O') >= 0;
  }

  /** Suffix information for a group carbon. */
  gcInfo(gc: number): GroupInfo | null {
    const g = this.g;
    const t = this.gcType[gc];
    if (!t || t === 'bad' || t === 'carbonic') return null;
    const base = (cls: number, kind: SuffixKind, consumed: number[]): GroupInfo => ({ cls, kind, gc, n: -1, consumed, esterO: -1, esterR: -1, halide: '' });
    if (t === 'nitrile') {
      const n = g.nb[gc].find((j) => g.el[j] === 'N' && g.order(gc, j) === 3)!;
      return base(CLS.NITRILE, 'nitrile', [n]);
    }
    const o = g.terminalDouble(gc, 'O');
    if (t === 'aldehyde') return base(CLS.ALDEHYDE, 'aldehyde', [o]);
    const x = g.nb[gc].find((j) => g.el[j] !== 'C' && !g.inRing[j] && g.order(gc, j) === 1)!;
    switch (t) {
      case 'acid': return base(CLS.ACID, 'acid', [o, x]);
      case 'carboxylate': return base(CLS.ANION_COO, 'carboxylate', [o, x]);
      case 'ester': {
        const gi = base(CLS.ESTER, 'ester', [o, x]);
        gi.esterO = x;
        gi.esterR = g.others(x, gc)[0];
        return gi;
      }
      case 'amide': {
        const gi = base(CLS.AMIDE, 'amide', [o, x]);
        gi.n = x;
        return gi;
      }
      case 'acylhalide': {
        const gi = base(CLS.ACYL_HALIDE, 'acylhalide', [o, x]);
        gi.halide = g.el[x];
        return gi;
      }
    }
    return null;
  }

  /** Is N an amine-type nitrogen (only C/H neighbours, no acyl/sulfonyl/hetero attached)? */
  private isAmineN(x: number): boolean {
    const g = this.g;
    if (g.el[x] !== 'N' || g.inRing[x]) return false;
    for (const j of g.nb[x]) {
      if (g.order(x, j) !== 1) return false;
      if (g.el[j] !== 'C' && !g.inRing[j]) return false;
      if (this.gcType[j]) return false;
      if (g.el[j] === 'C' && !g.inRing[j] && (g.terminalDouble(j, 'S') >= 0 || g.nb[j].some((k) => g.el[k] === 'N' && g.order(j, k) === 2))) return false;
    }
    return true;
  }

  /** Characteristic group formed by neighbour x of parent atom s (suffix candidate), or null. */
  attachedGroup(s: number, x: number): GroupInfo | null {
    const g = this.g;
    if (this.gcType[x]) {
      const gi = this.gcInfo(x);
      return gi;
    }
    if (g.inRing[x] || this.skel[x]) return null;
    const ord = g.order(s, x);
    const base = (cls: number, kind: SuffixKind, consumed: number[], n = -1): GroupInfo => ({ cls, kind, gc: -1, n, consumed, esterO: -1, esterR: -1, halide: '' });
    const sC = g.el[s] === 'C';
    if (g.el[x] === 'O') {
      if (!sC) return null;
      if (ord === 2 && g.degree(x) === 1) return base(CLS.KETONE, 'one', [x]);
      if (ord === 1 && g.degree(x) === 1) {
        if (g.charge[x] === -1) return base(CLS.ANION_O, 'olate', [x]);
        if (g.charge[x] === 0 && g.h[x] === 1) return base(CLS.ALCOHOL, 'ol', [x]);
      }
      return null;
    }
    if (g.el[x] === 'S') {
      if (!sC) return null;
      if (ord === 2 && g.degree(x) === 1) return base(CLS.THIONE, 'thione', [x]);
      if (ord !== 1 || g.charge[x] !== 0) return null;
      if (g.degree(x) === 1 && g.h[x] === 1) return base(CLS.THIOL, 'thiol', [x]);
      const oxo = g.nb[x].filter((j) => g.el[j] === 'O' && g.order(x, j) === 2 && g.degree(j) === 1);
      const rest = g.others(x, s, ...oxo);
      if (oxo.length === 2 && rest.length === 1) {
        const y = rest[0];
        if (g.el[y] === 'O' && g.degree(y) === 1) {
          if (g.charge[y] === -1) return base(CLS.ANION_SO3, 'sulfonate', [x, ...oxo, y]);
          if (g.h[y] === 1) return base(CLS.SULFONIC, 'sulfonic', [x, ...oxo, y]);
        }
        if (g.el[y] === 'O' && g.degree(y) === 2 && g.charge[y] === 0) {
          const r = g.others(y, x)[0];
          if (this.skel[r] || g.inRing[r]) {
            const gi = base(CLS.SULFONATE_ESTER, 'sulfonate', [x, ...oxo, y]);
            gi.esterO = y;
            gi.esterR = r;
            return gi;
          }
        }
        if (g.el[y] === 'N' && !g.inRing[y] && g.charge[y] === 0) return base(CLS.SULFONAMIDE, 'sulfonamide', [x, ...oxo, y], y);
        if (HALOGENS.has(g.el[y])) {
          const gi = base(CLS.SULFONYL_HALIDE, 'sulfonylhalide', [x, ...oxo, y]);
          gi.halide = g.el[y];
          return gi;
        }
      }
      if (oxo.length === 1 && rest.length === 1) {
        const y = rest[0];
        if (g.el[y] === 'O' && g.degree(y) === 1 && g.h[y] === 1) return base(CLS.SULFINIC, 'sulfinic', [x, ...oxo, y]);
      }
      return null;
    }
    if (g.el[x] === 'N') {
      if (!sC) {
        // N-amino on a ring nitrogen (e.g. piperidin-1-amine)
        if (!(g.inRing[s] && g.el[s] === 'N' && ord === 1 && g.degree(x) === 1 && g.charge[x] === 0)) return null;
        return base(CLS.AMINE, 'amine', [x], x);
      }
      if (ord === 2) {
        if (g.charge[x] !== 0) return null;
        const o = g.others(x, s);
        if (o.length > 1) return null;
        return base(CLS.IMINE, 'imine', [x], x);
      }
      if (ord !== 1) return null;
      if (g.charge[x] === 1 && g.nb[x].every((j) => g.order(x, j) === 1 && (g.el[j] === 'C' || g.inRing[j]) && !this.gcType[j])) {
        return base(CLS.CATION, 'aminium', [x], x);
      }
      if (g.charge[x] === 0 && this.isAmineN(x)) return base(CLS.AMINE, 'amine', [x], x);
      return null;
    }
    return null;
  }

  // ───────────────────────── carbonic acid derivatives / inorganic esters ─────────────────────────

  /** Class of a carbonic-acid derivative carbon (urea, carbamate, carbonate, …). */
  carbonicClass(c: number): number {
    const g = this.g;
    const het = g.nb[c].filter((j) => g.el[j] !== 'C' && !g.inRing[j] && g.order(c, j) === 1);
    const els = het.map((j) => g.el[j]).sort().join(',');
    const oKinds = het.filter((j) => g.el[j] === 'O').map((j) => (g.charge[j] === -1 ? 'm' : g.degree(j) === 1 ? 'h' : 'r'));
    if (els === 'N,N') return CLS.UREA;
    if (els === 'N,O') return oKinds[0] === 'm' ? CLS.ANION_COO : oKinds[0] === 'h' ? CLS.CARBAMIC : CLS.CARBAMATE;
    if (els === 'O,O') {
      if (oKinds.includes('m')) return CLS.ANION_COO;
      if (oKinds.includes('r')) return CLS.ESTER_INORG;
      return CLS.ACID;
    }
    if (/^(Br|Cl|F|I),N$|^N,(Br|Cl|F|I)$/.test(els)) return CLS.ACYL_HALIDE;
    if (/^(Br|Cl|F|I),O$|^O,(Br|Cl|F|I)$/.test(els)) return oKinds[0] === 'r' ? CLS.ESTER_INORG : 999;
    return 999;
  }

  /** Inorganic oxoacid ester centre (sulfate, phosphate, nitrate) – returns class or 999. */
  inorganicCentre(z: number): { cls: number; kind: string; ors: number[] } | null {
    const g = this.g;
    if (g.inRing[z]) return null;
    const el = g.el[z];
    const oxo = g.nb[z].filter((j) => g.el[j] === 'O' && g.degree(j) === 1 && (g.order(z, j) === 2 || (el === 'N' && g.charge[j] === -1)));
    const sing = g.nb[z].filter((j) => !oxo.includes(j));
    if (sing.some((j) => g.el[j] !== 'O' || g.order(z, j) !== 1)) return null;
    let kind = '';
    if (el === 'S' && g.charge[z] === 0 && oxo.length === 2 && sing.length === 2) kind = 'sulfate';
    else if (el === 'P' && g.charge[z] === 0 && oxo.length === 1 && sing.length === 3) kind = 'phosphate';
    else if (el === 'N' && ((g.charge[z] === 0 && oxo.length === 2) || (g.charge[z] === 1 && oxo.length === 2)) && sing.length === 1) kind = 'nitrate';
    if (!kind) return null;
    const esterified = sing.filter((o) => g.degree(o) === 2);
    if (!esterified.length) return null;
    const ok = esterified.every((o) => {
      const r = g.others(o, z)[0];
      return this.skel[r];
    });
    if (!ok) return null;
    const anion = sing.some((o) => g.charge[o] === -1 && g.degree(o) === 1);
    return { cls: anion ? CLS.ANION_INORG : CLS.ESTER_INORG, kind, ors: esterified };
  }

  // ───────────────────────── ring parents ─────────────────────────

  ringParent(sys: number): RingParent {
    let rp = this.ringCache.get(sys);
    if (!rp) {
      const r = buildRingParent(this.g, sys);
      if (!r.parent) throw new NamingError(r.error ?? 'unsupported ring system');
      rp = r.parent;
      this.ringCache.set(sys, rp);
    }
    return rp;
  }

  // ───────────────────────── top level ─────────────────────────

  /** Names the whole (single-component) graph. */
  nameAll(): { name: string; locants: Map<number, string> } {
    const g = this.g;
    if (this.gcType.some((t) => t === 'bad')) throw new NamingError('functional group not supported (thio/imido acids, anhydrides, cyanates)');
    for (let i = 0; i < g.n; i++) {
      if (this.skel[i] && g.charge[i] !== 0 && !(g.el[i] === 'N' && g.charge[i] === 1 && g.inRing[i])) {
        throw new NamingError('charged skeletal atoms are not supported');
      }
    }
    // principal class
    const classes = this.principalClasses();
    for (const cls of classes) {
      const r = this.nameWithClass(cls);
      if (r) return r;
    }
    const r = this.nameWithClass(null);
    if (!r) throw new NamingError('no parent hydride found');
    return r;
  }

  /** Candidate principal classes present in the molecule, most senior first. */
  principalClasses(): number[] {
    const g = this.g;
    const set = new Set<number>();
    for (let s = 0; s < g.n; s++) {
      if (this.gcType[s] && this.gcType[s] !== 'carbonic') {
        const gi = this.gcInfo(s);
        if (gi) set.add(gi.cls);
      }
      if (this.gcType[s] === 'carbonic') {
        const c = this.carbonicClass(s);
        if (c < 999) set.add(c);
      }
      if (!this.skel[s]) {
        const ic = this.inorganicCentre(s);
        if (ic) set.add(ic.cls);
        continue;
      }
      if (g.inRing[s] && g.el[s] === 'N' && g.charge[s] === 1) set.add(CLS.CATION);
      for (const x of g.nb[s]) {
        if (this.skel[x] || this.gcType[x]) continue;
        const gi = this.attachedGroup(s, x);
        if (gi) set.add(gi.cls);
      }
    }
    return [...set].sort((a, b) => a - b);
  }

  private nameWithClass(cls: number | null): { name: string; locants: Map<number, string> } | null {
    if (cls === CLS.UREA || cls === CLS.CARBAMIC || cls === CLS.CARBAMATE) return this.nameCarbonic(cls);
    if (cls === CLS.ESTER_INORG || cls === CLS.ANION_INORG) {
      const r = this.nameInorganicEster(cls);
      if (r) return r;
      return this.nameCarbonic(cls);
    }
    if (cls === CLS.ACYL_HALIDE) {
      const r = this.nameCarbonic(cls);
      if (r) return r;
    }
    if (cls === CLS.ANION_COO || cls === CLS.ACID) {
      const r = this.nameCarbonic(cls);
      if (r) return r;
    }
    const cand = this.selectParent(cls, null);
    if (!cand) return null;
    const res = cand.res;
    if (cls === CLS.ESTER || cls === CLS.SULFONATE_ESTER) {
      return { name: this.esterName(res), locants: res.locants };
    }
    return { name: res.text, locants: res.locants };
  }

  private esterName(res: BuildResult): string {
    const alk = res.esterAlkyls;
    const texts = new Set(alk.map((a) => a.sub.text));
    let pre: string;
    if (texts.size === 1) {
      const s = alk[0].sub;
      const n = alk.length;
      if (n === 1) pre = s.text;
      else if (s.compound) pre = multiplierComplex(n) + enclose(s.text);
      else pre = multiplier(n) + (startsAmbiguous(s.text) ? enclose(s.text) : s.text);
    } else {
      const sorted = [...alk].sort((a, b) => (alphaKey(a.sub.text) < alphaKey(b.sub.text) ? -1 : 1));
      pre = sorted.map((a) => (a.loc && alk.length > 1 ? a.loc + '-' : '') + a.sub.text).join(' ');
    }
    return pre + ' ' + res.text;
  }

  // ───────────────────────── parent selection ─────────────────────────

  /** Atoms reachable from `start` without crossing the bond start–from (whole component if from < 0). */
  region(start: number, from: number): Set<number> {
    const g = this.g;
    const seen = new Set<number>([start]);
    const st = [start];
    while (st.length) {
      const v = st.pop()!;
      for (const w of g.nb[v]) {
        if (v === start && w === from) continue;
        if (!seen.has(w)) {
          seen.add(w);
          st.push(w);
        }
      }
    }
    return seen;
  }

  /** Can group carbon gc terminate a chain for the given principal class / context? */
  private chainExtendable(gc: number, cls: number | null, main: boolean): boolean {
    const t = this.gcType[gc];
    if (!t || t === 'bad' || t === 'carbonic') return false;
    const gi = this.gcInfo(gc)!;
    if (cls !== null && gi.cls === cls) return true;
    if (t === 'nitrile') return false;
    if (t === 'acid' || t === 'carboxylate') return main && cls !== null && cls < gi.cls;
    return true; // aldehyde, ester, amide, acyl halide → oxo + prefix
  }

  /**
   * Enumerates chain candidates within `reg`. If `through` ≥ 0 the chain must contain it
   * (substituent free valence); if `startAt` ≥ 0 the chain must start at it (acyl groups).
   */
  private chainCandidates(reg: Set<number>, cls: number | null, main: boolean, through: number, startAt: number): number[][] {
    const g = this.g;
    const core = [...reg].filter((a) => this.skel[a] && !g.inRing[a] && g.el[a] === 'C');
    const coreSet = new Set(core);
    const cadj = new Map<number, number[]>();
    for (const a of core) cadj.set(a, g.nb[a].filter((b) => coreSet.has(b)));
    // terminal group carbons
    const ext = new Map<number, number[]>();
    for (const a of core) {
      const l = g.nb[a].filter((b) => reg.has(b) && this.gcType[b] && this.chainExtendable(b, cls, main));
      if (l.length) ext.set(a, l);
    }
    const chains: number[][] = [];
    const pathBetween = (u: number, v: number): number[] => {
      if (u === v) return [u];
      const prev = new Map<number, number>([[u, -1]]);
      const q = [u];
      for (let k = 0; k < q.length; k++) {
        const x = q[k];
        if (x === v) break;
        for (const y of cadj.get(x)!) if (!prev.has(y)) { prev.set(y, x); q.push(y); }
      }
      if (!prev.has(v)) return [];
      const p: number[] = [];
      for (let x = v; x !== -1; x = prev.get(x)!) p.push(x);
      return p.reverse();
    };
    type End = { u: number; g: number };
    const ends: End[] = [];
    for (const a of core) {
      if (cadj.get(a)!.length <= 1 || a === through || a === startAt) ends.push({ u: a, g: -1 });
      for (const gc of ext.get(a) ?? []) ends.push({ u: a, g: gc });
    }
    // standalone group carbons (formic acid, formamide, …)
    for (const a of reg) {
      if (this.gcType[a] && this.chainExtendable(a, cls, main) && !g.nb[a].some((b) => coreSet.has(b) || g.inRing[b])) {
        if ((through < 0 || through === a) && (startAt < 0 || startAt === a)) chains.push([a]);
        // two directly bonded group carbons (oxalic acid, glyoxal, oxamide …)
        for (const b of g.nb[a]) {
          if (!reg.has(b) || !this.gcType[b] || !this.chainExtendable(b, cls, main)) continue;
          const ch = [a, b];
          if (through >= 0 && !ch.includes(through)) continue;
          if (startAt >= 0 && a !== startAt && b !== startAt) continue;
          chains.push(ch);
        }
      }
    }
    if (startAt >= 0 && this.gcType[startAt]) {
      // acyl chain starting at a group carbon is handled by the caller
    }
    for (let i = 0; i < ends.length; i++) {
      for (let j = i; j < ends.length; j++) {
        const e1 = ends[i], e2 = ends[j];
        if (i === j && e1.g >= 0) continue;
        const p = pathBetween(e1.u, e2.u);
        if (!p.length) continue;
        if (i === j && p.length !== 1) continue;
        if (e1.u === e2.u && e1.g >= 0 && e2.g >= 0 && e1.g === e2.g) continue;
        const chain = [...(e1.g >= 0 ? [e1.g] : []), ...p, ...(e2.g >= 0 ? [e2.g] : [])];
        if (i === j && e1.g < 0 && cadj.get(e1.u)!.length > 0 && e1.u !== through && e1.u !== startAt) continue;
        if (through >= 0 && !chain.includes(through)) continue;
        if (startAt >= 0 && chain[0] !== startAt && chain[chain.length - 1] !== startAt) continue;
        chains.push(chain);
      }
    }
    // dedupe (reverse duplicates)
    const seen = new Set<string>();
    const out: number[][] = [];
    for (const c of chains) {
      const k1 = c.join(','), k2 = [...c].reverse().join(',');
      if (seen.has(k1) || seen.has(k2)) continue;
      seen.add(k1);
      out.push(startAt >= 0 && c[0] !== startAt ? [...c].reverse() : c);
    }
    return out;
  }

  /** Number of principal groups a chain can express as suffixes, and whether group carbons are included. */
  private chainCount(chain: number[], cls: number | null): { n: number; incl: boolean } {
    if (cls === null) return { n: 0, incl: false };
    const g = this.g;
    const inC = new Set(chain);
    let incl = 0, att = 0, other = 0;
    for (const s of chain) {
      if (this.gcType[s]) {
        const gi = this.gcInfo(s);
        if (gi && gi.cls === cls) incl++;
        continue;
      }
      for (const x of g.nb[s]) {
        if (inC.has(x)) continue;
        const gi = this.attachedGroup(s, x);
        if (!gi || gi.cls !== cls) continue;
        if (gi.gc >= 0) att++;
        else other++;
      }
    }
    if (incl > 0) return { n: incl + other, incl: true };
    return { n: att + other, incl: false };
  }

  private ringCount(rp: RingParent, cls: number | null): number {
    if (cls === null) return 0;
    const g = this.g;
    let n = 0;
    for (const s of rp.atoms) {
      if (cls === CLS.CATION && g.el[s] === 'N' && g.charge[s] === 1) n++;
      for (const x of g.nb[s]) {
        if (rp.atomSet.has(x)) continue;
        const gi = this.attachedGroup(s, x);
        if (gi && gi.cls === cls) n++;
      }
    }
    return n;
  }

  private multipleBonds(chain: number[]): { mb: number; db: number } {
    let mb = 0, db = 0;
    for (let i = 0; i + 1 < chain.length; i++) {
      const o = this.g.order(chain[i], chain[i + 1]);
      if (o >= 2) mb++;
      if (o === 2) db++;
    }
    return { mb, db };
  }

  /** Ring seniority key (P-44.2): larger is more senior. */
  private ringKey(rp: RingParent): number[] {
    const g = this.g;
    const het = rp.hetCount;
    const hetTotal = [...het.values()].reduce((a, b) => a + b, 0);
    const order = ['N', 'O', 'S', 'Se', 'Te', 'P', 'As', 'Si', 'B'];
    let D = 0;
    for (const a of rp.atoms) for (const b of g.nb[a]) if (b > a && rp.atomSet.has(b) && g.order(a, b) === 2) D++;
    return [het.has('N') ? 1 : 0, hetTotal ? 1 : 0, rp.nRings, rp.atoms.length, hetTotal, het.size, ...order.map((e) => het.get(e) ?? 0), D];
  }

  /**
   * Selects the senior parent for the principal class `cls`. With `fv`, the parent must contain the
   * free-valence atom (substituent naming).
   */
  selectParent(cls: number | null, fv: FV | null, ctx?: string): { parent: Parent; res: BuildResult } | null {
    this.check();
    const g = this.g;
    const reg = fv ? this.region(fv.kind === 'carbonyl' ? fv.carbonylC! : fv.atom, fv.from) : this.region(0, -1);
    type Cand = { parent: Parent; n: number; key: number[]; incl: boolean };
    const cands: Cand[] = [];
    // ring candidates
    if (!fv || fv.kind !== 'oyl') {
      const systems = new Set<number>();
      for (const a of reg) if (g.sysOf[a] >= 0) systems.add(g.sysOf[a]);
      for (const sys of systems) {
        if (fv && g.sysOf[fv.atom] !== sys) continue;
        const rp = this.ringParent(sys);
        cands.push({ parent: { kind: 'ring', atoms: rp.atoms, ring: rp }, n: this.ringCount(rp, cls), key: [1, ...this.ringKey(rp)], incl: false });
      }
    }
    // chain candidates
    if (!fv || (!g.inRing[fv.atom] && fv.kind !== 'carbonyl')) {
      const chains = this.chainCandidates(reg, cls, !fv, fv && fv.kind === 'yl' ? fv.atom : -1, fv && fv.kind === 'oyl' ? fv.atom : -1);
      for (const ch of chains) {
        const { n, incl } = this.chainCount(ch, cls);
        const { mb, db } = this.multipleBonds(ch);
        cands.push({ parent: { kind: 'chain', atoms: ch }, n, key: [0, ch.length, mb, db], incl });
      }
    }
    if (!cands.length) return null;
    let maxN = Math.max(...cands.map((c) => c.n));
    if (cls !== null && !fv && maxN === 0) return null;
    let best = cands.filter((c) => c.n === maxN);
    // rings before chains, then ring seniority / chain length and unsaturation
    best.sort((a, b) => -cmpNum(a.key, b.key));
    best = best.filter((c) => cmpNum(c.key, best[0].key) === 0);
    let top: { parent: Parent; res: BuildResult } | null = null;
    for (const c of best) {
      const res = this.buildName(c.parent, { cls, inclMode: c.incl, fv: fv ?? undefined, ctx });
      if (!top || this.compareResults(res, top.res) < 0) top = { parent: c.parent, res };
    }
    void maxN;
    maxN = 0;
    return top;
  }

  private compareResults(a: BuildResult, b: BuildResult): number {
    return (
      cmpNum(a.suffixLocs, b.suffixLocs) ||
      cmpNum(a.unsatLocs, b.unsatLocs) ||
      b.nPrefixes - a.nPrefixes ||
      cmpNum(a.prefixLocs, b.prefixLocs) ||
      cmpStr(a.alpha, b.alpha) ||
      (a.text < b.text ? -1 : a.text > b.text ? 1 : 0)
    );
  }

  // ───────────────────────── building a name on a parent ─────────────────────────

  buildName(P: Parent, opts: BuildOpts): BuildResult {
    this.check();
    const g = this.g;
    const inP = new Set(P.atoms);
    const consumed = new Set<number>();
    const suffixes: SuffixItem[] = [];
    const prefixes: Item[] = [];
    const nSubs: { nAtom: number; host: number; sub: SubName }[] = [];
    const esterAlk: { host: number; sub: SubName }[] = [];
    const fv = opts.fv;
    const cls = opts.cls;

    const addNSubs = (nAtom: number, host: number, exclude: number[]) => {
      for (const y of g.nb[nAtom]) {
        if (exclude.includes(y) || consumed.has(y)) continue;
        nSubs.push({ nAtom, host, sub: this.nameSubstituent(y, nAtom, g.order(nAtom, y), '') });
      }
    };
    const takeGroup = (s: number, gi: GroupInfo, inChain: boolean) => {
      suffixes.push({ atom: s, kind: inChain || !GC_KINDS.has(gi.kind) ? gi.kind : gi.kind, halide: gi.halide });
      (suffixes[suffixes.length - 1] as SuffixItem & { attached?: boolean });
      for (const c of gi.consumed) consumed.add(c);
      if (gi.gc >= 0 && !inChain) consumed.add(gi.gc);
      if (gi.n >= 0) addNSubs(gi.n, s, [s, gi.gc, ...gi.consumed]);
      if (gi.esterR >= 0) {
        const host = s;
        esterAlk.push({ host, sub: this.nameSubstituent(gi.esterR, gi.esterO, 1, 'ester') });
        consumed.add(gi.esterR);
      }
    };

    // free valence
    if (fv) {
      if (fv.kind === 'yl') suffixes.push({ atom: fv.atom, kind: fv.order === 2 ? 'ylidene' : fv.order === 3 ? 'ylidyne' : 'yl' });
      else if (fv.kind === 'oyl') {
        suffixes.push({ atom: fv.atom, kind: 'oyl' });
        consumed.add(g.terminalDouble(fv.atom, 'O'));
      } else {
        suffixes.push({ atom: fv.atom, kind: 'carbonyl' });
        consumed.add(fv.carbonylC!);
        consumed.add(g.terminalDouble(fv.carbonylC!, 'O'));
      }
    }
    const attachedSuffix = new Set<number>(); // hosts of attached group-carbon suffixes (for -carboxylic acid)
    for (const s of P.atoms) {
      // group carbon inside a chain
      if (this.gcType[s] && P.kind === 'chain') {
        const gi = this.gcInfo(s);
        if (gi && cls !== null && gi.cls === cls && opts.inclMode) takeGroup(s, gi, true);
      }
      if (P.kind === 'ring' && cls === CLS.CATION && g.el[s] === 'N' && g.charge[s] === 1) suffixes.push({ atom: s, kind: 'ium' });
      for (const x of g.nb[s]) {
        if (inP.has(x) || consumed.has(x)) continue;
        if (fv && s === fv.atom && x === fv.from) continue;
        if (fv && fv.kind === 'carbonyl' && x === fv.carbonylC) continue;
        if (cls !== null && !(this.gcType[s] && P.kind === 'chain')) {
          const gi = this.attachedGroup(s, x);
          if (gi && gi.cls === cls && !(P.kind === 'chain' && opts.inclMode && gi.gc >= 0)) {
            takeGroup(s, gi, false);
            if (gi.gc >= 0) attachedSuffix.add(s);
            continue;
          }
        }
        if (consumed.has(x)) continue;
        prefixes.push({ atom: s, sub: this.nameSubstituent(x, s, g.order(s, x), '') });
      }
    }
    // mark attached-group-carbon suffix kinds
    const isAttachedGC = (it: SuffixItem) => GC_KINDS.has(it.kind) && !(P.kind === 'chain' && this.gcType[it.atom]);

    // ── numbering ──
    const numberings: Numbering[] = P.kind === 'chain' ? this.chainNumberings(P.atoms) : P.ring!.numberings;
    const rp = P.ring;
    // ring hydro data (style mancude)
    let S: number[] = [], K = new Set<number>(), Dset = new Set<number>();
    let useSat = false;
    if (rp && rp.style === 'mancude') {
      for (const a of rp.atoms) if (g.nb[a].some((b) => rp.atomSet.has(b) && g.order(a, b) === 2)) Dset.add(a);
      S = rp.atoms.filter((a) => rp.A.has(a) && !Dset.has(a));
      for (const a of S) if (g.nb[a].some((b) => !rp.atomSet.has(b) && g.order(a, b) === 2)) K.add(a);
      useSat = !!rp.satName && Dset.size === 0;
    }
    const unsatBonds: { a: number; b: number; order: number }[] = [];
    if (P.kind === 'chain') {
      for (let i = 0; i + 1 < P.atoms.length; i++) {
        const o = g.order(P.atoms[i], P.atoms[i + 1]);
        if (o >= 2) unsatBonds.push({ a: P.atoms[i], b: P.atoms[i + 1], order: o });
      }
    } else if (rp!.style === 'ene') {
      for (const a of rp!.atoms) for (const b of g.nb[a]) if (b > a && rp!.atomSet.has(b) && g.order(a, b) >= 2) unsatBonds.push({ a, b, order: g.order(a, b) });
    }
    const fvAtom = fv ? fv.atom : -1;
    // locant of a multiple bond: lower locant; ring closure 1–n of a monocycle → n; otherwise "a(b)"
    const nRingAtoms = rp ? rp.atoms.length : 0;
    const bondLoc = (N: Numbering, a: number, b: number): Loc => {
      const la = N.get(a)!, lb = N.get(b)!;
      const [lo, hi] = la.value <= lb.value ? [la, lb] : [lb, la];
      if (!rp || /^\d+$/.test(lo.label) && /^\d+$/.test(hi.label) && hi.value - lo.value === 1) return lo;
      if (rp.nRings === 1 && lo.value === 1 && hi.value === nRingAtoms) return hi;
      if (/^\d+$/.test(lo.label) && /^\d+$/.test(hi.label) && hi.value - lo.value === 1) return lo;
      return { label: lo.label + '(' + hi.label + ')', value: lo.value + 0.001 };
    };
    // forced (implied) saturated positions for the hydro display
    const hEff = (a: number) => g.h[a];
    const forced = new Set<number>();
    if (S.length) {
      for (const a of S) if (K.has(a) || (hEff(a) === 0)) forced.add(a);
      for (const a of S) {
        if (forced.has(a)) continue;
        const partner = g.nb[a].some((q) => rp!.atomSet.has(q) && rp!.A.has(q) && !K.has(q) && (Dset.has(q) || hEff(q) > 0));
        if (!partner) forced.add(a);
      }
    }
    const E = S.filter((a) => !forced.has(a));

    let bestN: Numbering | null = null;
    let bestKey: number[][] | null = null;
    const prefixAll = [...prefixes.map((p) => ({ atom: p.atom, sub: p.sub }))];
    const alphaOrder = [...new Set(prefixAll.map((p) => p.sub.text))].sort((a, b) => (alphaKey(a) < alphaKey(b) ? -1 : alphaKey(a) > alphaKey(b) ? 1 : a < b ? -1 : 1));
    for (const N of numberings) {
      const L = (a: number) => N.get(a)!.value;
      const key: number[][] = [];
      // (a) indicated hydrogen of the parent
      let ihAtom = -1;
      if (rp && rp.style === 'mancude' && rp.needsIH && !useSat) {
        let best = Infinity;
        for (const a of S) if (rp.validIH.has(a) && L(a) < best) { best = L(a); ihAtom = a; }
        key.push([best]);
      }
      // (b) principal groups and free valence
      key.push(suffixes.map((s) => L(s.atom)).sort((p, q) => p - q));
      // (c) hydro prefixes / ene-yne endings
      if (rp && rp.style === 'mancude' && !useSat) key.push(S.filter((a) => a !== ihAtom && !K.has(a)).map(L).sort((p, q) => p - q));
      else key.push(unsatBonds.map((b) => bondLoc(N, b.a, b.b).value).sort((p, q) => p - q));
      key.push(unsatBonds.filter((b) => b.order === 2).map((b) => bondLoc(N, b.a, b.b).value).sort((p, q) => p - q));
      // (d) detachable prefixes together
      key.push(prefixAll.map((p) => L(p.atom)).sort((p, q) => p - q));
      // (e) in order of citation
      for (const t of alphaOrder) key.push(prefixAll.filter((p) => p.sub.text === t).map((p) => L(p.atom)).sort((p, q) => p - q));
      // (f) stereo: Z before E, R before S (CIP-based lowest locants) – approximated by locants of R/Z
      if (this.cip) {
        const rz: number[] = [];
        for (const a of P.atoms) {
          const d = this.cip.centers.get(a);
          if (d === 'R' || d === 'r') rz.push(L(a));
        }
        key.push(rz.sort((p, q) => p - q));
      }
      if (!bestKey || this.cmpKeys(key, bestKey) < 0) {
        bestKey = key;
        bestN = N;
      }
    }
    const N = bestN!;
    const L = (a: number) => N.get(a)!;
    const locants = new Map<number, string>();
    for (const a of P.atoms) locants.set(a, L(a).label);

    // ── assemble ──
    const nSuffix = suffixes.length;
    const suffixKind = suffixes.length ? suffixes[0].kind : null;
    const isChain = P.kind === 'chain';
    const chainLen = isChain ? P.atoms.length : 0;
    let parentCore = '';
    let hydroPart = '';
    // unsaturation text for chains / ene rings
    const unsatText = (stem: string, isRing: boolean): string => {
      const enesL = unsatBonds.filter((b) => b.order === 2).map((b) => bondLoc(N, b.a, b.b)).sort((p, q) => p.value - q.value);
      const ynesL = unsatBonds.filter((b) => b.order === 3).map((b) => bondLoc(N, b.a, b.b)).sort((p, q) => p.value - q.value);
      const enes = enesL.map((l) => l.value);
      const ynes = ynesL.map((l) => l.value);
      const labOf = new Map<number, string>();
      for (const l of [...enesL, ...ynesL]) labOf.set(l.value, l.label);
      const lab = (v: number) => labOf.get(v) ?? String(v);
      if (!enes.length && !ynes.length) return stem + 'ane';
      const omitLoc = (isChain && chainLen <= 2) || (isRing && rp?.carbocycleMono && enes.length === 1 && !ynes.length && enes[0] === 1);
      let s = stem;
      const first = enes.length ? enes.length : ynes.length;
      if (first > 1) s += 'a';
      if (enes.length) s += (omitLoc ? '' : '-' + enes.map(lab).join(',') + '-') + multiplier(enes.length) + 'ene';
      if (ynes.length) {
        if (enes.length) s = s.slice(0, -1); // "en" before "yne"
        s += (omitLoc ? '' : '-' + ynes.map(lab).join(',') + '-') + multiplier(ynes.length) + 'yne';
      }
      return s;
    };
    if (isChain) parentCore = unsatText(alkaneStem(chainLen), false);
    else if (rp!.style === 'ene') parentCore = unsatText(rp!.stem, true);
    else if (rp!.style === 'benzene') parentCore = 'benzene';
    else if (useSat) parentCore = rp!.satName;
    else {
      // hydro prefixes and indicated hydrogen (PubChem-style display)
      const eSorted = [...E].sort((a, b) => L(a).value - L(b).value);
      let ind = -1;
      let hydro = eSorted;
      if (eSorted.length % 2 === 1) {
        ind = eSorted[0];
        hydro = eSorted.slice(1);
      }
      parentCore = (ind >= 0 ? L(ind).label + 'H-' : '') + rp!.name;
      if (hydro.length) hydroPart = hydro.map((a) => L(a).label).join(',') + '-' + multiplier(hydro.length) + 'hydro' + (/^\d/.test(parentCore) ? '-' : '');
    }

    // suffix text
    let suffixText = '';
    let retainedDone = false;
    const sufLocs = suffixes.map((s) => L(s.atom)).sort((p, q) => p.value - q.value);
    const omitAllLocants = isChain && chainLen === 1;
    let omitSuffixLoc = omitAllLocants;
    const totalItems = suffixes.length + prefixes.length + nSubs.length;
    if (isChain && chainLen === 2) {
      if (totalItems === 1 && !nSubs.length) omitSuffixLoc = true;
      else if (suffixes.length === 1) omitSuffixLoc = true;
    }
    if (rp && rp.carbocycleMono && rp.style !== 'ene' && suffixes.length === 1 && (suffixKind === 'yl' || rp.style === 'benzene')) omitSuffixLoc = true;
    if (rp && rp.carbocycleMono && suffixes.length === 1 && prefixes.length === 0 && unsatBonds.length === 0) omitSuffixLoc = true;
    if (rp && rp.carbocycleMono && rp.style === 'ene' && suffixKind === 'yl' && unsatBonds.length === 0) omitSuffixLoc = true;
    const locStr = (locs: Loc[]) => locs.map((l) => l.label).join(',');
    const halideName = (h?: string) => HALIDE[h ?? 'Cl'];
    if (suffixKind) {
      const n = suffixes.length;
      const mult = multiplier(n);
      const inChainGC = isChain && GC_KINDS.has(suffixKind) && !isAttachedGC(suffixes[0]);
      if (inChainGC || suffixKind === 'oyl') {
        const t: Record<string, string> = {
          acid: 'oic acid', carboxylate: 'oate', ester: 'oate', amide: 'amide', nitrile: 'nitrile', aldehyde: 'al', oyl: 'oyl',
          acylhalide: 'oyl ' + halideName(suffixes[0].halide),
        };
        let body = t[suffixKind];
        if (n > 1) {
          body = suffixKind === 'acylhalide' ? mult + 'oyl ' + mult + halideName(suffixes[0].halide) : mult + body;
        }
        suffixText = body;
      } else {
        const t: Record<string, string> = {
          acid: 'carboxylic acid', carboxylate: 'carboxylate', ester: 'carboxylate', amide: 'carboxamide', nitrile: 'carbonitrile',
          aldehyde: 'carbaldehyde', one: 'one', thione: 'thione', ol: 'ol', olate: 'olate', thiol: 'thiol', amine: 'amine',
          aminium: 'aminium', imine: 'imine', sulfonic: 'sulfonic acid', sulfonate: 'sulfonate', sulfonamide: 'sulfonamide',
          sulfonylhalide: 'sulfonyl ' + halideName(suffixes[0].halide), sulfinic: 'sulfinic acid', ium: 'ium', yl: 'yl', ylidene: 'ylidene',
          ylidyne: 'ylidyne', carbonyl: 'carbonyl',
        };
        let body = t[suffixKind];
        if (suffixKind === 'acylhalide') body = 'carbonyl ' + halideName(suffixes[0].halide);
        if (n > 1) {
          if (suffixKind === 'acylhalide') body = mult + 'carbonyl ' + mult + halideName(suffixes[0].halide);
          else if (suffixKind === 'sulfonylhalide') body = mult + 'sulfonyl ' + mult + halideName(suffixes[0].halide);
          else body = (mult.endsWith('a') && /^o/.test(body) ? mult.slice(0, -1) : mult) + body; // "tetrol", "pentol"
        }
        suffixText = (omitSuffixLoc ? '' : '-' + locStr(sufLocs) + '-') + body;
      }
    }

    // retained names
    let fullParent = '';
    const benz = rp && rp.style === 'benzene';
    if (benz && suffixes.length === 1 && !fv) {
      const map: Record<string, string> = {
        acid: 'benzoic acid', carboxylate: 'benzoate', ester: 'benzoate', amide: 'benzamide', nitrile: 'benzonitrile',
        aldehyde: 'benzaldehyde', ol: 'phenol', olate: 'phenolate', amine: 'aniline', aminium: 'anilinium',
        acylhalide: 'benzoyl ' + halideName(suffixes[0].halide),
      };
      if (map[suffixKind!]) {
        fullParent = map[suffixKind!];
        retainedDone = true;
      }
    }
    if (benz && fv && fv.kind === 'carbonyl') {
      fullParent = 'benzoyl';
      retainedDone = true;
    }
    if (benz && suffixes.length === 2 && suffixKind === 'acid' && !fv && prefixes.length + nSubs.length >= 0) {
      const v = sufLocs.map((l) => l.value);
      const nm = v[1] === 2 ? 'phthalic acid' : v[1] === 3 ? 'isophthalic acid' : v[1] === 4 ? 'terephthalic acid' : '';
      if (nm && v[0] === 1) {
        fullParent = nm;
        retainedDone = true;
      }
    }
    if (isChain && chainLen <= 2 && suffixKind && (GC_KINDS.has(suffixKind) || suffixKind === 'oyl') && !isAttachedGC(suffixes[0])) {
      const c1: Record<string, string> = {
        acid: 'formic acid', carboxylate: 'formate', ester: 'formate', amide: 'formamide', nitrile: 'formonitrile', aldehyde: 'formaldehyde', oyl: 'formyl',
      };
      const c2: Record<string, string> = {
        acid: 'acetic acid', carboxylate: 'acetate', ester: 'acetate', amide: 'acetamide', nitrile: 'acetonitrile', aldehyde: 'acetaldehyde', oyl: 'acetyl',
        acylhalide: 'acetyl ' + halideName(suffixes[0].halide),
      };
      if (chainLen === 1 && suffixes.length === 1 && c1[suffixKind]) { fullParent = c1[suffixKind]; retainedDone = true; }
      if (chainLen === 2 && suffixes.length === 1 && c2[suffixKind] && unsatBonds.length === 0) { fullParent = c2[suffixKind]; retainedDone = true; }
      if (chainLen === 2 && suffixes.length === 2 && suffixKind === 'acid') { fullParent = 'oxalic acid'; retainedDone = true; }
      if (chainLen === 2 && suffixes.length === 2 && (suffixKind === 'ester' || suffixKind === 'carboxylate')) { fullParent = 'oxalate'; retainedDone = true; }
      if (chainLen === 2 && suffixes.length === 2 && suffixKind === 'amide') { fullParent = 'oxamide'; retainedDone = true; }
      if (chainLen === 2 && suffixes.length === 2 && suffixKind === 'aldehyde') { fullParent = 'oxaldehyde'; retainedDone = true; }
      if (chainLen === 2 && suffixes.length === 2 && suffixKind === 'nitrile') { fullParent = 'oxalonitrile'; retainedDone = true; }
      if (chainLen === 2 && suffixes.length === 2 && suffixKind === 'acylhalide') { fullParent = 'oxalyl di' + halideName(suffixes[0].halide); retainedDone = true; }
    }
    if (!retainedDone) {
      if (fv && fv.kind === 'yl' && benz && fv.order === 1) fullParent = 'phenyl';
      else if (fv && fv.kind === 'yl' && isChain && L(fv.atom).value === 1) {
        // free valence at C1: "propyl", "prop-2-enyl", "methylidene"
        const sfx = suffixText.replace(/^-1-/, '');
        fullParent = unsatBonds.length ? elide(parentCore, sfx) : alkaneStem(chainLen) + sfx;
      } else if (fv && fv.kind === 'yl' && rp?.carbocycleMono && rp.style === 'ene' && unsatBonds.length === 0) {
        fullParent = rp.stem + suffixText.replace(/^-1-/, '');
      } else fullParent = elide(parentCore, suffixText);
    }

    // prefixes
    const locantless = omitAllLocants;
    let omitPrefixLoc = omitAllLocants;
    if (isChain && chainLen === 2 && totalItems === 1) omitPrefixLoc = true;
    if (rp && rp.carbocycleMono && totalItems === 1 && unsatBonds.length === 0 && !fv) omitPrefixLoc = true;
    const items: PrefixItem[] = [];
    for (const p of prefixes) items.push({ sub: p.sub, locant: omitPrefixLoc ? '' : L(p.atom).label, value: L(p.atom).value });
    // N-substituents: N, N', N'' by order of the host locant
    const nAtoms = [...new Set(nSubs.map((x) => x.nAtom))];
    const nHostLoc = (nA: number) => L(nSubs.find((x) => x.nAtom === nA)!.host).value;
    const allN = suffixes.length; // all suffix nitrogens count for primes
    void allN;
    const suffixNs: number[] = [];
    for (const s of suffixes) {
      // nitrogen of each suffix group, ordered by locant
      void s;
    }
    nAtoms.sort((p, q) => nHostLoc(p) - nHostLoc(q));
    const nGroupCount = suffixes.filter((s) => ['amine', 'amide', 'aminium', 'imine', 'sulfonamide'].includes(s.kind)).length;
    for (const x of nSubs) {
      let label = 'N';
      if (nGroupCount > 1) {
        // primes follow the order of the suffix locants
        const hostLocs = suffixes.filter((s) => ['amine', 'amide', 'aminium', 'imine', 'sulfonamide'].includes(s.kind)).map((s) => L(s.atom).value).sort((p, q) => p - q);
        const idx = hostLocs.indexOf(L(x.host).value);
        label = 'N' + "'".repeat(Math.max(0, idx));
      }
      items.push({ sub: x.sub, locant: label, value: -1 });
      suffixNs.push(x.nAtom);
    }
    let prefixText = formatPrefixes(items, locantless);

    // special retained substituent names
    let text = '';
    const onlyPrefix = prefixes.length === 1 && !nSubs.length ? prefixes[0].sub : null;
    if (!fv && benz && suffixes.length === 0 && prefixes.length === 1 && onlyPrefix) {
      if (onlyPrefix.text === 'methyl') { text = 'toluene'; prefixText = ''; fullParent = ''; }
      else if (onlyPrefix.text === 'ethenyl') { text = 'styrene'; prefixText = ''; fullParent = ''; }
      else if (onlyPrefix.text === 'methoxy') { text = 'anisole'; prefixText = ''; fullParent = ''; }
      else if (onlyPrefix.phenyl) { text = '1,1\'-biphenyl'; prefixText = ''; fullParent = ''; }
    }
    if (fv && fv.kind === 'yl' && opts.ctx !== 'oxy' && isChain) {
      // benzyl / benzylidene / tert-butyl (PubChem style; not used inside "…oxy")
      if (chainLen === 1 && prefixes.length === 1 && prefixes[0].sub.phenyl && !nSubs.length) {
        if (fv.order === 1) { text = 'benzyl'; prefixText = ''; fullParent = ''; }
        else if (fv.order === 2) { text = 'benzylidene'; prefixText = ''; fullParent = ''; }
      }
      if (chainLen === 3 && fv.order === 1 && L(fv.atom).value === 2 && prefixes.length === 1 && prefixes[0].sub.text === 'methyl' && prefixes[0].atom === fv.atom && unsatBonds.length === 0) {
        text = 'tert-butyl'; prefixText = ''; fullParent = '';
      }
    }
    if (!text) {
      const rest = hydroPart + fullParent;
      text = prefixText + (prefixText && /^\d/.test(rest) && /[a-z)\]}]$/.test(prefixText) ? '-' : '') + rest;
    }

    // stereodescriptors for this parent
    const stereo = this.stereoPrefix(P, L, unsatBonds.length > 0 || true);
    if (stereo) text = stereo + text;

    // comparison data
    const prefixLocs = prefixes.map((p) => L(p.atom).value).sort((p, q) => p - q);
    const alpha = items.map((i) => alphaKey(i.sub.text)).sort();
    const unsatLocs = rp && rp.style === 'mancude' && !useSat ? S.filter((a) => !K.has(a)).map((a) => L(a).value).sort((p, q) => p - q) : unsatBonds.map((b) => bondLoc(N, b.a, b.b).value).sort((p, q) => p - q);
    const esterAlkyls = esterAlk.map((e) => ({ sub: e.sub, loc: L(e.host).label }));
    const hasPrefixes = items.length > 0 || hydroPart !== '';
    return {
      text,
      suffixLocs: sufLocs.map((l) => l.value),
      unsatLocs,
      nPrefixes: items.length,
      prefixLocs,
      alpha,
      locants,
      esterAlkyls,
      nSuffix,
      hasPrefixes,
      oxyContract: !!fv && fv.kind === 'yl' && fv.order === 1 && ((isChain && L(fv.atom).value === 1) || !!benz) && text.endsWith('yl') && text !== 'benzyl',
      phenyl: !!benz && !!fv && fv.kind === 'yl' && fv.order === 1 && items.length === 0,
      phenylLike: !!benz && !!fv && fv.kind === 'yl' && fv.order === 1 && !stereo,
    };
  }

  private cmpKeys(a: number[][], b: number[][]): number {
    for (let i = 0; i < Math.min(a.length, b.length); i++) {
      const c = cmpNum(a[i], b[i]);
      if (c) return c;
    }
    return 0;
  }

  private chainNumberings(atoms: number[]): Numbering[] {
    const mk = (arr: number[]) => {
      const m: Numbering = new Map();
      arr.forEach((a, i) => m.set(a, { label: String(i + 1), value: i + 1 }));
      return m;
    };
    if (atoms.length === 1) return [mk(atoms)];
    return [mk(atoms), mk([...atoms].reverse())];
  }

  // ───────────────────────── stereodescriptors ─────────────────────────

  private stereoPrefix(P: Parent, L: (a: number) => Loc, _x: boolean): string {
    if (!this.cip) return '';
    const g = this.g;
    const inP = new Set(P.atoms);
    const items: { v: number; text: string; db: boolean }[] = [];
    for (const a of P.atoms) {
      const d = this.cip.centers.get(a);
      if (d && !this.stereoUsed.has('a' + a)) {
        this.stereoUsed.add('a' + a);
        items.push({ v: L(a).value, text: L(a).label + d, db: false });
      }
    }
    for (const [k, d] of this.cip.bonds) {
      const [i, j] = k.split(',').map(Number);
      if (this.stereoUsed.has('b' + k)) continue;
      const ii = inP.has(i), jj = inP.has(j);
      if (!ii && !jj) continue;
      // double bond between a parent atom and a substituent atom is cited with the parent locant
      if (!(ii && jj)) {
        const other = ii ? j : i;
        if (g.el[other] !== 'C' && g.el[other] !== 'N') continue;
      }
      this.stereoUsed.add('b' + k);
      const v = Math.min(ii ? L(i).value : Infinity, jj ? L(j).value : Infinity);
      const lab = ii && (!jj || L(i).value <= L(j).value) ? L(i).label : L(j).label;
      items.push({ v, text: lab + d, db: true });
    }
    if (!items.length) return '';
    const dbs = items.filter((x) => x.db);
    let parts: string[];
    if (dbs.length === 1) {
      // PubChem style: a single E/Z descriptor is cited without locant, first
      const rest = items.filter((x) => !x.db).sort((a, b) => a.v - b.v);
      parts = [dbs[0].text.replace(/^[\d]+[a-z]?'*/, ''), ...rest.map((x) => x.text)];
    } else parts = items.sort((a, b) => a.v - b.v).map((x) => x.text);
    return '(' + parts.join(',') + ')-';
  }

  // ───────────────────────── substituent prefixes ─────────────────────────

  /** Name of the substituent rooted at x, attached to parent atom `from` by a bond of order `order`. */
  nameSubstituent(x: number, from: number, order: number, ctx: string): SubName {
    const key = x + '|' + from + '|' + order + '|' + ctx;
    let s = this.subCache.get(key);
    if (!s) {
      this.check();
      s = this.nameSubstituentRaw(x, from, order, ctx);
      this.subCache.set(key, s);
    }
    return s;
  }

  private nameSubstituentRaw(x: number, from: number, order: number, ctx: string): SubName {
    const g = this.g;
    if (this.gcType[x]) return this.nameGCSub(x, from);
    if (this.skel[x]) return this.nameSkeletalSub(x, from, order, ctx);
    return this.nameHetSub(x, from, order);
  }

  /** Skeletal substituent (alkyl, aryl, acyl, …). */
  private nameSkeletalSub(x: number, from: number, order: number, ctx: string): SubName {
    const g = this.g;
    if (g.charge[x] !== 0) throw new NamingError('charged substituent atoms are not supported');
    // acyl: acyclic C(=O) attached by a single bond
    if (!g.inRing[x] && g.el[x] === 'C' && order === 1 && g.terminalDouble(x, 'O') >= 0) return this.acylName(x, from);
    if (!g.inRing[x] && g.el[x] === 'C' && order === 1 && (g.terminalDouble(x, 'S') >= 0)) throw new NamingError('thioacyl groups not supported');
    if (!g.inRing[x] && g.el[x] !== 'C') throw new NamingError('unsupported acyclic skeletal atom');
    const sel = this.selectParent(null, { atom: x, from, order, kind: 'yl' }, ctx);
    if (!sel) throw new NamingError('cannot name substituent');
    const r = sel.res;
    const substituted = r.hasPrefixes && r.text !== 'tert-butyl' && r.text !== 'benzyl' && r.text !== 'benzylidene';
    return {
      text: r.text,
      compound: substituted,
      enclose: substituted,
      oxyContract: r.oxyContract,
      phenyl: r.phenyl,
      phenylLike: r.phenylLike,
    };
  }

  /** Acyl group C(=O)–R attached through the carbonyl carbon c. */
  private acylName(c: number, from: number): SubName {
    const g = this.g;
    const r = g.nb[c].find((j) => j !== from && (this.skel[j] || g.inRing[j]) && g.order(c, j) === 1);
    if (r === undefined) {
      // no carbon substituent: formyl (H), or heteroatom → handled by group-carbon naming
      return { text: 'formyl', compound: false, enclose: false };
    }
    let res: BuildResult;
    if (g.inRing[r]) {
      const sel = this.selectParent(null, { atom: r, from: c, order: 1, kind: 'carbonyl', carbonylC: c });
      if (!sel) throw new NamingError('cannot name acyl group');
      res = sel.res;
    } else {
      const sel = this.selectParent(null, { atom: c, from, order: 1, kind: 'oyl' });
      if (!sel) throw new NamingError('cannot name acyl group');
      res = sel.res;
    }
    return { text: res.text, compound: res.hasPrefixes, enclose: false };
  }

  /** Substituent rooted at a group carbon. */
  private nameGCSub(gc: number, from: number): SubName {
    const g = this.g;
    const t = this.gcType[gc]!;
    const fromIsGroupHet = !this.skel[from] && !g.inRing[from] && g.el[from] !== 'C';
    if (t === 'carbonic') {
      // from one heteroatom: name by the other one
      const other = g.nb[gc].find((j) => j !== from && g.el[j] !== 'C' && !g.inRing[j] && g.order(gc, j) === 1)!;
      return this.carbonylWith(other, gc);
    }
    if (fromIsGroupHet) {
      // gc reached from its own heteroatom (ester O / amide N): acyl group
      return this.acylName(gc, from);
    }
    switch (t) {
      case 'acid': return { text: 'carboxy', compound: false, enclose: false };
      case 'carboxylate': return { text: 'carboxylato', compound: false, enclose: false };
      case 'aldehyde': return { text: 'formyl', compound: false, enclose: false };
      case 'nitrile': return { text: 'cyano', compound: false, enclose: false };
    }
    const x = g.nb[gc].find((j) => j !== from && g.el[j] !== 'C' && !g.inRing[j] && g.order(gc, j) === 1)!;
    return this.carbonylWith(x, gc);
  }

  /** "-C(=O)-X" named from the carbon side: carboxy, alkoxycarbonyl, carbamoyl, carbonochloridoyl. */
  private carbonylWith(x: number, gc: number): SubName {
    const g = this.g;
    if (g.el[x] === 'O') {
      if (g.charge[x] === -1) return { text: 'carboxylato', compound: false, enclose: false };
      const o = g.others(x, gc);
      if (!o.length) return { text: 'carboxy', compound: false, enclose: false };
      const r = this.nameSubstituent(o[0], x, 1, 'oxy');
      return { text: this.oxyName(r) + 'carbonyl', compound: true, enclose: false };
    }
    if (g.el[x] === 'N') {
      const subs = g.others(x, gc).map((y) => this.nameSubstituent(y, x, g.order(x, y), ''));
      if (!subs.length) return { text: 'carbamoyl', compound: false, enclose: false };
      return { text: this.combineSubs(subs) + 'carbamoyl', compound: true, enclose: true };
    }
    if (HALOGENS.has(g.el[x])) return { text: 'carbono' + HALIDE[g.el[x]].replace(/ide$/, 'id') + 'oyl', compound: false, enclose: false };
    throw new NamingError('unsupported carbonyl substituent');
  }

  /** "methoxy", "phenoxy", "propan-2-yloxy", "(2-methylpropan-2-yl)oxy". */
  oxyName(r: SubName): string {
    if (r.oxyContract) return r.text.slice(0, -2) + 'oxy';
    return (r.compound || startsAmbiguous(r.text) ? enclose(r.text) : r.text) + 'oxy';
  }

  /** Joins substituents of an amino/carbamoyl/silyl centre: "dimethyl", "ethyl(methyl)", "bis(2-hydroxyethyl)". */
  combineSubs(subs: SubName[]): string {
    const items: PrefixItem[] = subs.map((s) => ({ sub: s, locant: '', value: 0 }));
    return formatPrefixes(items, true);
  }

  private withTail(r: SubName, tail: string): string {
    return (r.compound || startsAmbiguous(r.text) ? enclose(r.text) : r.text) + tail;
  }

  /** Heteroatom-rooted substituents. */
  private nameHetSub(x: number, from: number, order: number): SubName {
    const g = this.g;
    const el = g.el[x];
    const simple = (text: string): SubName => ({ text, compound: false, enclose: false });
    const others = g.others(x, from);
    if (HALOGENS.has(el)) {
      if (others.length || g.charge[x]) throw new NamingError('hypervalent halogen not supported');
      return simple(HALO_PREFIX[el]);
    }
    if (el === 'O') {
      if (order === 2) return simple('oxo');
      if (order !== 1) throw new NamingError('unsupported oxygen substituent');
      if (!others.length) return simple(g.charge[x] === -1 ? 'oxido' : 'hydroxy');
      if (others.length !== 1 || g.charge[x] !== 0) throw new NamingError('unsupported oxygen substituent');
      const y = others[0];
      if (this.gcType[y]) {
        const a = this.nameGCSub(y, x);
        return { text: this.withTail(a, 'oxy'), compound: a.compound, enclose: false };
      }
      if (this.skel[y]) {
        const r = this.nameSubstituent(y, x, 1, 'oxy');
        return { text: this.oxyName(r), compound: r.compound, enclose: false };
      }
      if (g.el[y] === 'O') {
        const o2 = g.others(y, x);
        if (!o2.length && g.h[y] === 1) return simple('hydroperoxy');
        if (o2.length === 1 && this.skel[o2[0]]) {
          const r = this.nameSubstituent(o2[0], y, 1, '');
          return { text: this.withTail(r, 'peroxy'), compound: true, enclose: false };
        }
        throw new NamingError('unsupported peroxide');
      }
      const h = this.nameHetSub(y, x, g.order(x, y));
      return { text: this.withTail(h, 'oxy'), compound: h.compound, enclose: false };
    }
    if (el === 'N') return this.nameNSub(x, from, order, others);
    if (el === 'S') {
      if (order === 2) {
        if (others.length) throw new NamingError('unsupported sulfur substituent');
        return simple('sulfanylidene');
      }
      if (g.charge[x] !== 0) throw new NamingError('charged sulfur not supported');
      const oxo = others.filter((j) => g.el[j] === 'O' && g.order(x, j) === 2 && g.degree(j) === 1);
      const rest = others.filter((j) => !oxo.includes(j));
      if (oxo.length === 0) {
        if (!rest.length) return simple('sulfanyl');
        if (rest.length === 1) {
          const y = rest[0];
          if (g.el[y] === 'C' && g.nb[y].some((k) => g.el[k] === 'N' && g.order(y, k) === 3)) return simple('thiocyanato');
          if (this.skel[y]) {
            const r = this.nameSubstituent(y, x, 1, '');
            return { text: this.withTail(r, 'sulfanyl'), compound: true, enclose: false };
          }
          if (this.gcType[y]) {
            const a = this.nameGCSub(y, x);
            return { text: this.withTail(a, 'sulfanyl'), compound: true, enclose: false };
          }
          if (g.el[y] === 'S') {
            const r2 = g.others(y, x);
            if (r2.length === 1 && this.skel[r2[0]]) {
              const r = this.nameSubstituent(r2[0], y, 1, '');
              return { text: this.withTail(r, 'disulfanyl'), compound: true, enclose: false };
            }
          }
        }
        throw new NamingError('unsupported sulfur substituent');
      }
      if (rest.length !== 1) throw new NamingError('unsupported sulfur substituent');
      const y = rest[0];
      const tail = oxo.length === 1 ? 'sulfinyl' : oxo.length === 2 ? 'sulfonyl' : '';
      if (!tail) throw new NamingError('unsupported sulfur substituent');
      if (g.el[y] === 'O' && g.degree(y) === 1) {
        if (oxo.length === 2) return simple(g.charge[y] === -1 ? 'sulfonato' : 'sulfo');
        return simple(g.charge[y] === -1 ? 'sulfinato' : 'sulfino');
      }
      if (g.el[y] === 'N' && !g.inRing[y] && oxo.length === 2) {
        const subs = g.others(y, x).map((z) => this.nameSubstituent(z, y, g.order(y, z), ''));
        if (!subs.length) return simple('sulfamoyl');
        return { text: this.combineSubs(subs) + 'sulfamoyl', compound: true, enclose: true };
      }
      if (HALOGENS.has(g.el[y]) && oxo.length === 2) return simple(HALO_PREFIX[g.el[y]] + 'sulfonyl');
      if (g.el[y] === 'O' && g.degree(y) === 2) {
        const r = this.nameSubstituent(g.others(y, x)[0], y, 1, 'oxy');
        return { text: this.oxyName(r) + tail, compound: true, enclose: false };
      }
      if (this.skel[y] || g.inRing[y]) {
        const r = this.nameSubstituent(y, x, 1, '');
        return { text: this.withTail(r, tail), compound: true, enclose: false };
      }
      throw new NamingError('unsupported sulfur substituent');
    }
    if (el === 'P') {
      const oxo = others.filter((j) => g.el[j] === 'O' && g.order(x, j) === 2 && g.degree(j) === 1);
      const rest = others.filter((j) => !oxo.includes(j));
      if (oxo.length === 1 && rest.length === 2 && rest.every((j) => g.el[j] === 'O')) {
        if (rest.every((j) => g.degree(j) === 1 && g.h[j] === 1)) return simple('phosphono');
        if (rest.every((j) => g.degree(j) === 2)) {
          const subs = rest.map((j) => this.nameSubstituent(g.others(j, x)[0], j, 1, 'oxy'));
          const names = subs.map((s) => ({ ...s, text: this.oxyName(s), compound: s.compound }));
          return { text: this.combineSubs(names) + 'phosphoryl', compound: true, enclose: false };
        }
      }
      throw new NamingError('unsupported phosphorus substituent');
    }
    if (el === 'Si') {
      if (order !== 1 || g.charge[x] !== 0) throw new NamingError('unsupported silicon substituent');
      const subs = others.map((y) => this.nameSubstituent(y, x, g.order(x, y), ''));
      if (!subs.length) return simple('silyl');
      return { text: this.combineSubs(subs) + 'silyl', compound: true, enclose: subs.length > 1 && new Set(subs.map((s) => s.text)).size > 1 };
    }
    throw new NamingError('substituent element ' + el + ' not supported');
  }

  private nameNSub(x: number, from: number, order: number, others: number[]): SubName {
    const g = this.g;
    const simple = (text: string): SubName => ({ text, compound: false, enclose: false });
    const oDbl = others.filter((j) => g.el[j] === 'O' && g.degree(j) === 1 && (g.order(x, j) === 2 || g.charge[j] === -1));
    // nitro / nitroso
    if (order === 1 && oDbl.length === 2 && others.length === 2) return simple('nitro');
    if (order === 1 && oDbl.length === 1 && others.length === 1 && g.order(x, oDbl[0]) === 2) return simple('nitroso');
    // azido
    if (order === 1 && others.length === 1 && g.el[others[0]] === 'N' && g.order(x, others[0]) === 2) {
      const n2 = others[0];
      const n3 = g.others(n2, x);
      if (n3.length === 1 && g.el[n3[0]] === 'N' && g.order(n2, n3[0]) === 2 && g.degree(n3[0]) === 1) return simple('azido');
      if (n3.length === 0 || (n3.length === 1)) {
        // diazenyl (azo) – not supported
      }
    }
    // isocyanato / isothiocyanato / isocyano
    if (order === 1 && others.length === 1 && g.el[others[0]] === 'C' && g.order(x, others[0]) === 2) {
      const c = others[0];
      if (g.terminalDouble(c, 'O') >= 0 && g.degree(c) === 2) return simple('isocyanato');
      if (g.terminalDouble(c, 'S') >= 0 && g.degree(c) === 2) return simple('isothiocyanato');
    }
    if (order === 1 && others.length === 1 && g.el[others[0]] === 'C' && g.order(x, others[0]) === 3 && g.degree(others[0]) === 1) return simple('isocyano');
    if (order === 2) {
      if (g.charge[x] !== 0) throw new NamingError('unsupported imino substituent');
      if (!others.length) return simple('imino');
      if (others.length === 1) {
        const y = others[0];
        if (g.el[y] === 'O' && g.degree(y) === 1 && g.h[y] === 1) return { text: 'hydroxyimino', compound: true, enclose: false };
        const r = this.nameSubstituent(y, x, g.order(x, y), '');
        return { text: this.withTail(r, 'imino'), compound: true, enclose: false };
      }
      throw new NamingError('unsupported imino substituent');
    }
    if (order !== 1) throw new NamingError('unsupported nitrogen substituent');
    if (g.charge[x] === 1) {
      const subs = others.map((y) => this.nameSubstituent(y, x, g.order(x, y), ''));
      if (!subs.length) return simple('azaniumyl');
      return { text: this.combineSubs(subs) + 'azaniumyl', compound: true, enclose: true };
    }
    if (g.charge[x] !== 0) throw new NamingError('unsupported nitrogen substituent');
    if (!others.length) return simple('amino');
    if (others.some((y) => g.order(x, y) !== 1)) throw new NamingError('unsupported nitrogen substituent');
    // hydrazinyl
    if (others.some((y) => g.el[y] === 'N')) {
      if (others.length === 1 && g.others(others[0], x).length === 0 && g.order(x, others[0]) === 1) return simple('hydrazinyl');
      throw new NamingError('hydrazines not supported');
    }
    const subs = others.map((y) => ({ y, s: this.nameSubstituent(y, x, 1, '') }));
    // acylamino with retained contracted forms
    if (subs.length === 1) {
      const t = subs[0].s.text;
      if (t === 'acetyl') return simple('acetamido');
      if (t === 'formyl') return simple('formamido');
      if (t === 'benzoyl') return simple('benzamido');
    }
    // anilino
    const ph = subs.find((p) => p.s.phenylLike);
    if (ph) {
      const rest = subs.filter((p) => p !== ph).map((p) => p.s);
      const an = ph.s.text.replace(/phenyl$/, 'anilino');
      if (!rest.length) return { text: an, compound: ph.s.compound, enclose: false };
      const items: PrefixItem[] = rest.map((s) => ({ sub: s, locant: 'N', value: 0 }));
      return { text: formatPrefixes(items, false) + (startsAmbiguous(an) ? '-' + enclose(an) : '-' + an).replace(/^-/, ''), compound: true, enclose: true };
    }
    return { text: this.combineSubs(subs.map((p) => p.s)) + 'amino', compound: true, enclose: true };
  }

  // ───────────────────────── functional parents ─────────────────────────

  /** Urea, carbamic acid, carbamates, carbonates, carbonic acid, carbamoyl halides, chloroformates. */
  private nameCarbonic(cls: number): { name: string; locants: Map<number, string> } | null {
    const g = this.g;
    const cs: number[] = [];
    for (let i = 0; i < g.n; i++) if (this.gcType[i] === 'carbonic' && this.carbonicClass(i) === cls) cs.push(i);
    if (!cs.length) return null;
    const c = cs[0];
    const het = g.nb[c].filter((j) => g.el[j] !== 'C' && !g.inRing[j] && g.order(c, j) === 1);
    const locants = new Map<number, string>();
    const nSubsOf = (n: number) => g.others(n, c).map((y) => this.nameSubstituent(y, n, g.order(n, y), ''));
    if (cls === CLS.UREA) {
      const [n1, n2] = het;
      const s1 = nSubsOf(n1), s2 = nSubsOf(n2);
      // lowest locants: the nitrogen with more substituents (then alphabetically first) gets 1
      let a = s1, b = s2, na = n1, nb = n2;
      const firstKey = (l: SubName[]) => l.map((s) => alphaKey(s.text)).sort()[0] ?? '~';
      if (s2.length > s1.length || (s2.length === s1.length && firstKey(s2) < firstKey(s1))) {
        a = s2; b = s1; na = n2; nb = n1;
      }
      locants.set(na, '1');
      locants.set(c, '2');
      locants.set(nb, '3');
      const total = a.length + b.length;
      const items: PrefixItem[] = [
        ...a.map((s) => ({ sub: s, locant: total === 1 ? '' : '1', value: 1 })),
        ...b.map((s) => ({ sub: s, locant: total === 1 ? '' : '3', value: 3 })),
      ];
      return { name: formatPrefixes(items, false) + 'urea', locants };
    }
    const nAt = het.find((j) => g.el[j] === 'N');
    const oAts = het.filter((j) => g.el[j] === 'O');
    const xAt = het.find((j) => HALOGENS.has(g.el[j]));
    const nPrefix = (n: number) => {
      const subs = nSubsOf(n);
      return formatPrefixes(subs.map((s) => ({ sub: s, locant: 'N', value: 0 })), false);
    };
    const alkyl = (o: number) => this.nameSubstituent(g.others(o, c)[0], o, 1, 'ester');
    const alkText = (list: SubName[]) => {
      const texts = new Set(list.map((s) => s.text));
      if (texts.size === 1 && list.length > 1) return (list[0].compound ? multiplierComplex(list.length) + enclose(list[0].text) : multiplier(list.length) + list[0].text);
      return list.map((s) => s.text).sort((p, q) => (alphaKey(p) < alphaKey(q) ? -1 : 1)).join(' ');
    };
    if (cls === CLS.CARBAMIC && nAt !== undefined) return { name: nPrefix(nAt) + 'carbamic acid', locants };
    if (cls === CLS.CARBAMATE && nAt !== undefined) return { name: alkyl(oAts[0]).text + ' ' + nPrefix(nAt) + 'carbamate', locants };
    if (cls === CLS.ANION_COO && nAt !== undefined) return { name: nPrefix(nAt) + 'carbamate', locants };
    if (cls === CLS.ACYL_HALIDE && nAt !== undefined && xAt !== undefined) return { name: nPrefix(nAt) + 'carbamoyl ' + HALIDE[g.el[xAt]], locants };
    if (cls === CLS.ESTER_INORG && xAt !== undefined && oAts.length === 1) return { name: alkyl(oAts[0]).text + ' carbono' + HALIDE[g.el[xAt]].replace(/ide$/, 'idate'), locants };
    if (oAts.length === 2) {
      const rs = oAts.filter((o) => g.degree(o) === 2).map(alkyl);
      const hs = oAts.filter((o) => g.degree(o) === 1 && g.charge[o] === 0).length;
      if (cls === CLS.ACID && !rs.length) return { name: 'carbonic acid', locants };
      const parts: string[] = [];
      if (rs.length) parts.push(alkText(rs));
      if (hs) parts.push(multiplier(hs) + 'hydrogen');
      return { name: parts.join(' ') + ' carbonate', locants };
    }
    return null;
  }

  /** Sulfate / phosphate / nitrate esters: "dodecyl hydrogen sulfate", "trimethyl phosphate". */
  private nameInorganicEster(cls: number): { name: string; locants: Map<number, string> } | null {
    const g = this.g;
    let best: { z: number; info: NonNullable<ReturnType<Namer['inorganicCentre']>> } | null = null;
    for (let z = 0; z < g.n; z++) {
      const info = this.inorganicCentre(z);
      if (!info || info.cls !== cls) continue;
      if (!best) best = { z, info };
      else {
        // prefer the centre whose alkyl attaches to the least substituted carbon
        const deg = (b: typeof best) => Math.min(...b!.info.ors.map((o) => g.degree(g.others(o, b!.z)[0])));
        if (deg({ z, info }) < deg(best)) best = { z, info };
      }
    }
    if (!best) return null;
    const { z, info } = best;
    const sing = g.nb[z].filter((j) => g.el[j] === 'O' && g.order(z, j) === 1 && !(g.el[z] === 'N' && g.charge[j] === -1 && g.degree(j) === 1 && g.nb[z].filter((k) => g.order(z, k) === 2).length < 2));
    const esterO = info.ors;
    const subs = esterO.map((o) => this.nameSubstituent(g.others(o, z)[0], o, 1, 'ester'));
    let hCount = 0;
    if (info.kind !== 'nitrate') hCount = sing.filter((o) => g.degree(o) === 1 && g.charge[o] === 0).length;
    const texts = new Set(subs.map((s) => s.text));
    let alk: string;
    if (texts.size === 1 && subs.length > 1) alk = subs[0].compound ? multiplierComplex(subs.length) + enclose(subs[0].text) : multiplier(subs.length) + subs[0].text;
    else alk = subs.map((s) => s.text).sort((p, q) => (alphaKey(p) < alphaKey(q) ? -1 : 1)).join(' ');
    const parts = [alk];
    if (hCount) parts.push(multiplier(hCount) + 'hydrogen');
    parts.push(info.kind);
    return { name: parts.join(' '), locants: new Map() };
  }
}
