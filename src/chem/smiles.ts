// SMILES reader and (canonical) writer.
import { BY_SYMBOL, element } from './elements';
import { Mol, TetraSpec, DbSpec, tetraCcwForOrder } from './mol';
import { kekulize, perceiveAromaticity } from './aromaticity';
import { implicitH, bondOrderSum, chargedValences } from './valence';
import { canonicalRanks } from './canon';
import { perceiveRings } from './rings';

export class SmilesError extends Error {
  constructor(msg: string, public pos: number) {
    super(`${msg} (at position ${pos + 1})`);
  }
}

const ORGANIC = new Set(['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', '*']);
const AROMATIC_ORGANIC = new Set(['b', 'c', 'n', 'o', 'p', 's']);
const AROMATIC_BRACKET = new Set(['b', 'c', 'n', 'o', 'p', 's', 'se', 'as', 'te']);

interface PendingBond {
  order: number;
  dir?: '/' | '\\';
  aromatic?: boolean;
  explicit: boolean;
}

/**
 * Parses a SMILES string (OpenSMILES subset incl. stereo, isotopes, charges, ring closures, '.').
 * Atoms receive no coordinates (x = y = 0); call a layout routine afterwards.
 * A reaction SMILES (with '>') is rejected — use parseReactionSmiles.
 */
export function parseSmiles(input: string): Mol {
  const s = input.trim().split(/\s+/)[0];
  const mol = new Mol();
  const aromaticAtom: boolean[] = [];
  const aromaticBond: boolean[] = [];
  const organicAtom: boolean[] = [];
  const nbrOrder: (number | { ring: number })[][] = [];
  const chiral: (null | '@' | '@@')[] = [];
  const bracketH: boolean[] = [];
  const hadPrev: boolean[] = [];
  const bondDirs: { bond: number; from: number; to: number; dir: '/' | '\\' }[] = [];
  const ringOpen = new Map<number, { atom: number; bond: PendingBond | null; slot: { ring: number } }>();
  const stack: number[] = [];
  let prev = -1;
  let pending: PendingBond | null = null;
  let i = 0;

  const addBondTo = (a: number, b: number, pb: PendingBond | null) => {
    let order = pb?.order ?? 1;
    let aro = false;
    if (!pb || !pb.explicit) {
      if (aromaticAtom[a] && aromaticAtom[b]) aro = true;
    } else if (pb.aromatic) aro = true;
    if (aro) order = 1;
    const bi = mol.addBond(a, b, order);
    aromaticBond[bi] = aro;
    if (pb?.dir) bondDirs.push({ bond: bi, from: a, to: b, dir: pb.dir });
    return bi;
  };

  const newAtom = (el: string, aro: boolean, organic: boolean) => {
    const idx = mol.addAtom({ el });
    aromaticAtom[idx] = aro;
    organicAtom[idx] = organic;
    nbrOrder[idx] = [];
    chiral[idx] = null;
    bracketH[idx] = false;
    hadPrev[idx] = prev >= 0;
    if (prev >= 0) {
      addBondTo(prev, idx, pending);
      nbrOrder[prev].push(idx);
      nbrOrder[idx].push(prev);
    }
    pending = null;
    prev = idx;
    return idx;
  };

  while (i < s.length) {
    const c = s[i];
    if (c === '(') {
      if (prev < 0) throw new SmilesError('Branch without preceding atom', i);
      stack.push(prev);
      i++;
    } else if (c === ')') {
      if (!stack.length) throw new SmilesError('Unmatched )', i);
      prev = stack.pop()!;
      pending = null;
      i++;
    } else if (c === '.') {
      prev = -1;
      pending = null;
      i++;
    } else if ('-=#$:/\\'.includes(c)) {
      const map: Record<string, PendingBond> = {
        '-': { order: 1, explicit: true },
        '=': { order: 2, explicit: true },
        '#': { order: 3, explicit: true },
        '$': { order: 4, explicit: true },
        ':': { order: 1, aromatic: true, explicit: true },
        '/': { order: 1, dir: '/', explicit: true },
        '\\': { order: 1, dir: '\\', explicit: true },
      };
      pending = map[c];
      i++;
    } else if (c === '%' || (c >= '0' && c <= '9')) {
      let num: number;
      if (c === '%') {
        if (s[i + 1] === '(') {
          const end = s.indexOf(')', i);
          num = parseInt(s.slice(i + 2, end), 10);
          i = end + 1;
        } else {
          num = parseInt(s.slice(i + 1, i + 3), 10);
          i += 3;
        }
      } else {
        num = +c;
        i++;
      }
      if (prev < 0) throw new SmilesError('Ring closure without atom', i - 1);
      const open = ringOpen.get(num);
      if (open) {
        const pb = pending ?? open.bond;
        let bi: number;
        if (pending && pending.dir) {
          bi = addBondTo(prev, open.atom, pending);
        } else if (open.bond && open.bond.dir) {
          bi = addBondTo(open.atom, prev, open.bond);
        } else {
          bi = addBondTo(open.atom, prev, pb);
        }
        void bi;
        const slotIdx = nbrOrder[open.atom].indexOf(open.slot);
        nbrOrder[open.atom][slotIdx] = prev;
        nbrOrder[prev].push(open.atom);
        ringOpen.delete(num);
      } else {
        const slot = { ring: num };
        nbrOrder[prev].push(slot);
        ringOpen.set(num, { atom: prev, bond: pending, slot });
      }
      pending = null;
    } else if (c === '[') {
      const end = s.indexOf(']', i);
      if (end < 0) throw new SmilesError('Unclosed [', i);
      const body = s.slice(i + 1, end);
      const m = /^(\d+)?([A-Z][a-z]?|[a-z][a-z]?|\*)(@@|@(?:TH[12]|AL[12]|SP[123]|TB\d+|OH\d+)?)?(H\d*)?([+-]+\d*|[+-]\d+)?(?::(\d+))?$/.exec(body);
      if (!m) throw new SmilesError(`Cannot parse bracket atom [${body}]`, i);
      let sym = m[2];
      let aro = false;
      if (sym !== '*' && sym[0] === sym[0].toLowerCase()) {
        if (!AROMATIC_BRACKET.has(sym)) throw new SmilesError(`Unknown aromatic atom ${sym}`, i);
        aro = true;
        sym = sym[0].toUpperCase() + sym.slice(1);
      }
      if (sym !== '*' && !BY_SYMBOL.has(sym)) {
        // e.g. [Cl] is fine; [Xx] not
        throw new SmilesError(`Unknown element ${sym}`, i);
      }
      const idx = newAtom(sym === '*' ? '*' : sym, aro, false);
      const atom = mol.atoms[idx];
      if (m[1]) atom.isotope = +m[1];
      if (m[3]) chiral[idx] = m[3].startsWith('@@') ? '@@' : m[3] === '@' || m[3] === '@TH1' ? '@' : m[3] === '@TH2' ? '@@' : null;
      const hTok = m[4];
      atom.hCount = hTok ? (hTok.length > 1 ? +hTok.slice(1) : 1) : 0;
      bracketH[idx] = !!hTok && atom.hCount > 0;
      if (m[5]) {
        const ch = m[5];
        const sign = ch[0] === '+' ? 1 : -1;
        if (/^[+-]\d+$/.test(ch)) atom.charge = sign * +ch.slice(1);
        else atom.charge = sign * ch.length;
      }
      if (m[6]) atom.map = +m[6];
      i = end + 1;
    } else {
      // organic subset
      let sym = '';
      if (c === 'C' && s[i + 1] === 'l') sym = 'Cl';
      else if (c === 'B' && s[i + 1] === 'r') sym = 'Br';
      else sym = c;
      if (ORGANIC.has(sym)) {
        newAtom(sym, false, true);
        i += sym.length;
      } else if (AROMATIC_ORGANIC.has(c)) {
        newAtom(c.toUpperCase(), true, true);
        i++;
      } else {
        throw new SmilesError(`Unexpected character '${c}'`, i);
      }
    }
  }
  if (stack.length) throw new SmilesError('Unclosed branch', s.length - 1);
  if (ringOpen.size) throw new SmilesError(`Unclosed ring ${[...ringOpen.keys()][0]}`, s.length - 1);

  // Implicit hydrogens of organic-subset atoms
  for (let a = 0; a < mol.atoms.length; a++) {
    if (!organicAtom[a] || mol.atoms[a].el === '*') {
      if (mol.atoms[a].el === '*') mol.atoms[a].hCount = 0;
      continue;
    }
    let sum = 0;
    for (const bi of mol.adj[a]) sum += aromaticBond[bi] ? 1 : mol.bonds[bi].order;
    if (aromaticAtom[a]) sum += 1;
    const vals = chargedValences(mol.atoms[a].el, 0);
    const v = vals.find((x) => x >= sum);
    mol.atoms[a].hCount = v === undefined ? 0 : v - sum;
  }

  if (aromaticBond.some((x) => x)) {
    kekulize(mol, aromaticAtom, aromaticBond);
  }

  // Tetrahedral stereo
  for (let a = 0; a < mol.atoms.length; a++) {
    const ch = chiral[a];
    if (!ch) continue;
    const order = nbrOrder[a].map((x) => (typeof x === 'number' ? x : -1));
    const h = mol.atoms[a].hCount ?? 0;
    if (h === 1 || (order.length === 3 && h === 0)) {
      // implicit H (or lone pair) goes right after the preceding atom, or first if there is none
      const pos = hadPrev[a] ? 1 : 0;
      order.splice(pos, 0, -1);
    }
    if (order.length !== 4) continue;
    mol.tetra.push({ center: a, nbrs: order as TetraSpec['nbrs'], ccw: ch === '@' });
  }

  // Double-bond stereo from directional bonds
  if (bondDirs.length) {
    const dirOf = (center: number, nbr: number): 'up' | 'down' | null => {
      for (const d of bondDirs) {
        if (d.from === center && d.to === nbr) return d.dir === '/' ? 'up' : 'down';
        if (d.to === center && d.from === nbr) return d.dir === '/' ? 'down' : 'up';
      }
      return null;
    };
    mol.bonds.forEach((b, bi) => {
      if (b.order !== 2) return;
      let ra = -1, da: 'up' | 'down' | null = null;
      for (const n of mol.neighbors(b.a)) {
        if (n === b.b) continue;
        const d = dirOf(b.a, n);
        if (d) { ra = n; da = d; break; }
      }
      let rb = -1, db: 'up' | 'down' | null = null;
      for (const n of mol.neighbors(b.b)) {
        if (n === b.a) continue;
        const d = dirOf(b.b, n);
        if (d) { rb = n; db = d; break; }
      }
      if (ra >= 0 && rb >= 0) mol.dbStereo.push({ bond: bi, a: ra, b: rb, cis: da === db });
    });
  }

  // Drop explicit H counts that equal the computed implicit value (keeps the editor flexible)
  for (let a = 0; a < mol.atoms.length; a++) {
    const at = mol.atoms[a];
    const hc = at.hCount;
    if (hc === undefined) continue;
    delete at.hCount;
    if (implicitH(mol, a) !== hc) at.hCount = hc;
  }
  return mol;
}

/** Parses "reactants>agents>products"; each part may contain several '.'-separated molecules. */
export function parseReactionSmiles(s: string): { reactants: Mol[]; agents: Mol[]; products: Mol[] } {
  const parts = s.trim().split('>');
  if (parts.length !== 3) throw new SmilesError('Reaction SMILES needs exactly two ">"', 0);
  const split = (p: string) => (p ? p.split('.').filter(Boolean).map((x) => parseSmiles(x)) : []);
  return { reactants: split(parts[0]), agents: split(parts[1]), products: split(parts[2]) };
}

// ───────────────────────────── Writer ─────────────────────────────

export interface SmilesWriteOptions {
  /** Write aromatic atoms in lowercase (default true). */
  aromatic?: boolean;
  /** Include stereo descriptors (default true). Requires mol.tetra / mol.dbStereo to be populated. */
  stereo?: boolean;
  /** Canonical atom ordering (default true). */
  canonical?: boolean;
  /** Include atom-map numbers. */
  maps?: boolean;
}

/**
 * Writes SMILES. Explicit hydrogen atoms are written as bracket atoms only when they carry
 * isotope/charge or bridge two atoms; otherwise callers should suppress them first (see suppressHydrogens).
 */
export function writeSmiles(input: Mol, opts: SmilesWriteOptions = {}): string {
  const mol = input;
  const useArom = opts.aromatic !== false;
  const useStereo = opts.stereo !== false;
  const n = mol.atoms.length;
  if (!n) return '';
  const rings = perceiveRings(mol);
  const aro = useArom ? perceiveAromaticity(mol, rings) : { atoms: new Array(n).fill(false), bonds: new Array(mol.bonds.length).fill(false) };
  const ranks = opts.canonical === false ? mol.atoms.map((_, i) => i) : canonicalRanks(mol, { rings, aromatic: aro.atoms });

  const tetraByCenter = new Map<number, TetraSpec>();
  if (useStereo) for (const t of mol.tetra) tetraByCenter.set(t.center, t);
  const dbByBond = new Map<number, DbSpec>();
  if (useStereo) for (const d of mol.dbStereo) dbByBond.set(d.bond, d);

  // DFS spanning forest
  const visited = new Array(n).fill(false);
  const parent = new Array(n).fill(-1);
  const children: number[][] = mol.atoms.map(() => []);
  const ringClosures: { a: number; b: number; bond: number }[] = [];
  const closureSeen = new Set<number>();
  const roots: number[] = [];
  const sortedAtoms = mol.atoms.map((_, i) => i).sort((p, q) => ranks[p] - ranks[q]);
  const dfsOrder: number[] = [];

  const dfs = (root: number) => {
    visited[root] = true;
    const rec = (v: number) => {
      dfsOrder.push(v);
      const nb = mol.adj[v]
        .map((bi) => ({ bi, w: mol.other(bi, v) }))
        .sort((p, q) => ranks[p.w] - ranks[q.w]);
      for (const { bi, w } of nb) {
        if (w === parent[v]) continue;
        if (visited[w]) {
          if (!closureSeen.has(bi)) {
            closureSeen.add(bi);
            ringClosures.push({ a: w, b: v, bond: bi });
          }
          continue;
        }
        visited[w] = true;
        parent[w] = v;
        children[v].push(w);
        rec(w);
      }
    };
    rec(root);
  };

  for (const a of sortedAtoms) {
    if (visited[a]) continue;
    // prefer starting from a terminal atom of the lowest rank in this component
    roots.push(a);
    dfs(a);
  }

  // Ring-closure bookkeeping: for each atom, list of closures it opens/closes in output order
  const closuresAt: Map<number, { partner: number; bond: number }[]> = new Map();
  // Opening atom = the one visited first in DFS order
  const pos = new Array(n).fill(0);
  dfsOrder.forEach((v, k) => (pos[v] = k));
  for (const rc of ringClosures) {
    const [first, second] = pos[rc.a] < pos[rc.b] ? [rc.a, rc.b] : [rc.b, rc.a];
    if (!closuresAt.has(first)) closuresAt.set(first, []);
    if (!closuresAt.has(second)) closuresAt.set(second, []);
    closuresAt.get(first)!.push({ partner: second, bond: rc.bond });
    closuresAt.get(second)!.push({ partner: first, bond: rc.bond });
  }
  // order closures at each atom: closing ones (partner earlier) first by partner position, then openings by partner rank
  for (const [v, list] of closuresAt) {
    list.sort((p, q) => {
      const pc = pos[p.partner] < pos[v] ? 0 : 1;
      const qc = pos[q.partner] < pos[v] ? 0 : 1;
      if (pc !== qc) return pc - qc;
      return pos[p.partner] - pos[q.partner];
    });
  }

  // Determine bond directions for stereo double bonds
  const bondDir = new Map<number, { from: number; dir: '/' | '\\' }>();
  if (useStereo && dbByBond.size) {
    assignBondDirections(mol, dbByBond, bondDir, pos);
  }

  const ringNum = new Map<number, number>(); // bond -> ring digit
  const used = new Set<number>();
  const nextDigit = () => {
    let d = 1;
    while (used.has(d)) d++;
    used.add(d);
    return d;
  };

  const bondSymbol = (bi: number, from: number, to: number): string => {
    const b = mol.bonds[bi];
    const dir = bondDir.get(bi);
    if (dir) {
      const d = dir.from === from ? dir.dir : dir.dir === '/' ? '\\' : '/';
      return d;
    }
    if (aro.bonds[bi]) return '';
    if (b.order === 2) return '=';
    if (b.order === 3) return '#';
    if (b.order === 4) return '$';
    if (b.order === 1.5) return ':';
    if (useArom && aro.atoms[from] && aro.atoms[to]) return '-';
    return '';
  };

  const atomSymbol = (v: number, outOrder: number[]): string => {
    const a = mol.atoms[v];
    const isAro = aro.atoms[v];
    let sym = a.el === 'R' ? '*' : a.el;
    const hTotal = a.el === '*' || a.el === 'R' ? 0 : implicitH(mol, v);
    let chir = '';
    const t = tetraByCenter.get(v);
    if (t) {
      // output neighbour order with -1 for implicit H at its position
      const order = [...outOrder];
      const hasH = t.nbrs.includes(-1);
      if (hasH) {
        const p = parent[v] >= 0 ? 1 : 0;
        order.splice(p, 0, -1);
      }
      if (order.length === 4 && order.every((x) => t.nbrs.includes(x))) {
        chir = tetraCcwForOrder(t, order) ? '@' : '@@';
      }
    }
    const organicOK =
      ORGANIC.has(sym) && !a.isotope && a.charge === 0 && !chir && !a.map && !(a.radical ?? 0) &&
      hTotal === defaultSmilesH(mol, v, isAro, aro.bonds);
    if (isAro) sym = sym.toLowerCase();
    if (organicOK && (!isAro || AROMATIC_ORGANIC.has(sym))) return sym;
    let out = '[';
    if (a.isotope) out += a.isotope;
    out += sym;
    out += chir;
    if (hTotal > 0) out += 'H' + (hTotal > 1 ? hTotal : '');
    if (a.charge) out += (a.charge > 0 ? '+' : '-') + (Math.abs(a.charge) > 1 ? Math.abs(a.charge) : '');
    if (a.map && opts.maps !== false) out += ':' + a.map;
    out += ']';
    return out;
  };

  const write = (v: number): string => {
    const closures = closuresAt.get(v) ?? [];
    const outOrder: number[] = [];
    if (parent[v] >= 0) outOrder.push(parent[v]);
    let ringStr = '';
    for (const c of closures) {
      outOrder.push(c.partner);
      let d = ringNum.get(c.bond);
      let bs = '';
      if (d === undefined) {
        d = nextDigit();
        ringNum.set(c.bond, d);
        bs = bondSymbol(c.bond, v, c.partner);
      } else {
        used.delete(d);
        ringNum.delete(c.bond);
        bs = bondSymbol(c.bond, v, c.partner);
        // avoid writing the bond symbol twice for non-directional bonds
        if (!bondDir.has(c.bond)) bs = '';
      }
      ringStr += bs + (d > 9 ? '%' + d : String(d));
    }
    for (const w of children[v]) outOrder.push(w);
    let out = atomSymbol(v, outOrder) + ringStr;
    const ch = children[v];
    for (let k = 0; k < ch.length; k++) {
      const w = ch[k];
      const bi = mol.bondBetween(v, w);
      const sub = bondSymbol(bi, v, w) + write(w);
      out += k < ch.length - 1 ? '(' + sub + ')' : sub;
    }
    return out;
  };

  return roots.map((r) => write(r)).join('.');
}

/** Implicit H count SMILES would assume for an organic-subset atom written without brackets. */
function defaultSmilesH(mol: Mol, v: number, aromatic: boolean, aroBonds: boolean[]): number {
  const a = mol.atoms[v];
  if (!ORGANIC.has(a.el) || a.el === '*') return -1;
  let sum = 0;
  for (const bi of mol.adj[v]) {
    const b = mol.bonds[bi];
    if (b.style === 'hbond' || b.style === 'dative') continue;
    // aromatic bonds are implicit and count 1 each; the π bond adds one more below
    sum += aromatic && (aroBonds[bi] || b.order === 1.5) ? 1 : b.order;
  }
  if (aromatic) sum += 1;
  const vals = chargedValences(a.el, 0);
  const val = vals.find((x) => x >= sum - 1e-6);
  return val === undefined ? 0 : Math.round(val - sum);
}

function assignBondDirections(
  mol: Mol,
  dbByBond: Map<number, DbSpec>,
  bondDir: Map<number, { from: number; dir: '/' | '\\' }>,
  pos: number[],
): void {
  // For each stereo double bond pick one single bond at each end and assign '/' or '\'
  // relative to the double-bond atom ("up"/"down"). Conflicts with already-assigned bonds are resolved
  // by deriving the needed direction from them.
  const sideOf = new Map<string, 'up' | 'down'>(); // key "center,nbr"
  const setSide = (center: number, nbr: number, side: 'up' | 'down') => {
    const bi = mol.bondBetween(center, nbr);
    // Express as direction written from the atom that comes first in output order
    const [from, to] = pos[center] < pos[nbr] ? [center, nbr] : [nbr, center];
    // if from === center: '/' means nbr is up. If from === nbr: '/' means nbr is down (center up from nbr).
    let dir: '/' | '\\';
    if (from === center) dir = side === 'up' ? '/' : '\\';
    else dir = side === 'up' ? '\\' : '/';
    bondDir.set(bi, { from, dir });
    void to;
    sideOf.set(center + ',' + nbr, side);
  };
  const getSide = (center: number, nbr: number): 'up' | 'down' | null => {
    const bi = mol.bondBetween(center, nbr);
    const d = bondDir.get(bi);
    if (!d) return null;
    if (d.from === center) return d.dir === '/' ? 'up' : 'down';
    return d.dir === '/' ? 'down' : 'up';
  };
  const singleNbrs = (center: number, exclude: number) =>
    mol.neighbors(center).filter((w) => w !== exclude && mol.bonds[mol.bondBetween(center, w)].order === 1);

  for (const [bi, d] of dbByBond) {
    const b = mol.bonds[bi];
    const nbA = singleNbrs(b.a, b.b);
    const nbB = singleNbrs(b.b, b.a);
    if (!nbA.length || !nbB.length) continue;
    // reference neighbours
    let refA = d.a, refB = d.b;
    let cis = d.cis;
    if (!nbA.includes(refA)) { refA = nbA[0]; cis = !cis; }
    if (!nbB.includes(refB)) { refB = nbB[0]; cis = !cis; }
    // existing assignment at A side?
    let sideA = getSide(b.a, refA);
    if (!sideA) {
      // try deriving from other neighbour of A
      const otherA = nbA.find((w) => w !== refA);
      const so = otherA !== undefined ? getSide(b.a, otherA) : null;
      sideA = so ? (so === 'up' ? 'down' : 'up') : 'up';
      setSide(b.a, refA, sideA);
    }
    const wantB: 'up' | 'down' = cis ? sideA : sideA === 'up' ? 'down' : 'up';
    const existingB = getSide(b.b, refB);
    if (!existingB) setSide(b.b, refB, wantB);
    // also mark the other neighbours for consistency
    for (const w of nbA) if (w !== refA && !getSide(b.a, w)) setSide(b.a, w, sideA === 'up' ? 'down' : 'up');
    for (const w of nbB) if (w !== refB && !getSide(b.b, w)) setSide(b.b, w, wantB === 'up' ? 'down' : 'up');
  }
}

/**
 * Removes "ordinary" explicit hydrogen atoms (neutral, natural isotope, single bond to one heavy atom)
 * and folds them into the heavy atom's hCount. Stereo specs that referenced them get -1.
 */
export function suppressHydrogens(mol: Mol): Mol {
  const m = mol.clone();
  const remove: number[] = [];
  for (let i = 0; i < m.atoms.length; i++) {
    const a = m.atoms[i];
    if (a.el !== 'H' || a.isotope || a.charge || m.degree(i) !== 1) continue;
    const j = m.neighbors(i)[0];
    if (m.atoms[j].el === 'H') continue;
    const b = m.bonds[m.adj[i][0]];
    if (b.order !== 1) continue;
    remove.push(i);
  }
  if (!remove.length) return m;
  const rm = new Set(remove);
  // fix heavy-atom H counts
  const addH = new Map<number, number>();
  for (const h of remove) {
    const j = m.neighbors(h)[0];
    addH.set(j, (addH.get(j) ?? 0) + 1);
  }
  const before = new Map<number, number>();
  for (const [j] of addH) before.set(j, implicitH(m, j));
  for (const t of m.tetra) t.nbrs = t.nbrs.map((x) => (rm.has(x) ? -1 : x)) as TetraSpec['nbrs'];
  for (const d of m.dbStereo) {
    if (rm.has(d.a) || rm.has(d.b)) {
      const b = m.bonds[d.bond];
      const swap = (end: number, ref: number) => {
        const alt = m.neighbors(end).find((w) => w !== ref && w !== (end === b.a ? b.b : b.a) && !rm.has(w));
        return alt;
      };
      if (rm.has(d.a)) {
        const alt = swap(b.a, d.a);
        if (alt === undefined) { d.bond = -1; continue; }
        d.a = alt; d.cis = !d.cis;
      }
      if (rm.has(d.b)) {
        const alt = swap(b.b, d.b);
        if (alt === undefined) { d.bond = -1; continue; }
        d.b = alt; d.cis = !d.cis;
      }
    }
  }
  m.dbStereo = m.dbStereo.filter((d) => d.bond >= 0);
  const map = m.removeAtoms(remove);
  for (const [j, k] of addH) {
    const nj = map[j];
    const a = m.atoms[nj];
    const want = (a.hCount ?? before.get(j)!) + k;
    delete a.hCount;
    if (implicitH(m, nj) !== want) a.hCount = want;
  }
  return m;
}

/** Converts an element symbol from SMILES casing. Exposed for tests. */
export function isKnownElement(sym: string): boolean {
  return !!element(sym);
}
