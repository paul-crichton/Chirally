// Arrow-pushing engine: applies curved (electron-flow) arrows to produce the next intermediate,
// with full electron bookkeeping (lone pairs, bonds, formal charges, radicals) and octet checks.
import { ChemDoc, CurvedArrowObj, Anchor } from './types';
import { docToMol, adjacency, fragmentOf } from './document';
import { Mol } from '../chem/mol';
import { implicitH, nonBondingElectrons, bondOrderSum, octetLimit } from '../chem/valence';
import { valenceElectrons, element } from '../chem/elements';

export interface MechanismWarning {
  message: string;
  atomIds?: number[];
}

export interface MechanismResult {
  /** Product structure; atoms keep the document ids of the atoms they came from. */
  product: Mol;
  /** Document atom ids of the reactant fragments used. */
  reactantAtomIds: number[];
  warnings: MechanismWarning[];
  /** true when σ-connectivity is unchanged (resonance structures → ↔ arrow). */
  resonance: boolean;
  /** Human-readable summary of what happened, e.g. "formed C–O, broke C–Br". */
  summary: string[];
}

function atomName(mol: Mol, i: number): string {
  const a = mol.atoms[i];
  return a.abbrev ?? a.alias ?? a.el;
}

/** Groups curved arrows into independent "steps": arrows sharing fragments belong together. */
export function arrowGroups(doc: ChemDoc): { arrows: number[]; atoms: number[]; maxX: number }[] {
  const adj = adjacency(doc);
  const fragOf = new Map<number, number>();
  const frags: number[][] = [];
  const getFrag = (atomId: number) => {
    let f = fragOf.get(atomId);
    if (f === undefined) {
      const comp = fragmentOf(doc, atomId, adj);
      f = frags.length;
      frags.push(comp);
      for (const a of comp) fragOf.set(a, f);
    }
    return f;
  };
  const anchorAtoms = (an: Anchor): number[] => {
    if (an.type === 'atom') return [an.id];
    if (an.type === 'bond') {
      const b = doc.bonds.get(an.id);
      return b ? [b.a, b.b] : [];
    }
    if (an.type === 'between') return [an.a, an.b];
    return [];
  };
  // union-find over fragments
  const parent: number[] = [];
  const find = (x: number): number => (parent[x] === x ? x : (parent[x] = find(parent[x])));
  const arrowFrags = new Map<number, number[]>();
  for (const c of doc.curved.values()) {
    const fs = [...anchorAtoms(c.from), ...anchorAtoms(c.to)].map(getFrag);
    for (const f of fs) if (parent[f] === undefined) parent[f] = f;
    for (let k = 1; k < fs.length; k++) {
      const a = find(fs[0]), b = find(fs[k]);
      if (a !== b) parent[a] = b;
    }
    arrowFrags.set(c.id, fs);
  }
  const groups = new Map<number, { arrows: number[]; atoms: Set<number> }>();
  for (const [cid, fs] of arrowFrags) {
    if (!fs.length) continue;
    const root = find(fs[0]);
    if (!groups.has(root)) groups.set(root, { arrows: [], atoms: new Set() });
    const g = groups.get(root)!;
    g.arrows.push(cid);
    for (const f of fs) for (const a of frags[f]) g.atoms.add(a);
  }
  return [...groups.values()].map((g) => {
    let maxX = -Infinity;
    for (const a of g.atoms) maxX = Math.max(maxX, doc.atoms.get(a)!.x);
    return { arrows: g.arrows, atoms: [...g.atoms], maxX };
  });
}

/**
 * Applies the given curved arrows. Each arrow moves 2 electrons (or 1 for fishhooks) from its source
 * (atom lone pair or bond) to its target (atom, bond, or the space between two atoms).
 */
export function applyArrows(doc: ChemDoc, arrowIds: number[]): MechanismResult {
  const warnings: MechanismWarning[] = [];
  const summary: string[] = [];
  const arrows = arrowIds.map((id) => doc.curved.get(id)).filter(Boolean) as CurvedArrowObj[];
  // involved fragments
  const adj = adjacency(doc);
  const involved = new Set<number>();
  const addFrag = (id: number) => {
    if (involved.has(id) || !doc.atoms.has(id)) return;
    for (const a of fragmentOf(doc, id, adj)) involved.add(a);
  };
  for (const c of arrows) {
    for (const an of [c.from, c.to]) {
      if (an.type === 'atom') addFrag(an.id);
      else if (an.type === 'bond') {
        const b = doc.bonds.get(an.id);
        if (b) addFrag(b.a);
      } else if (an.type === 'between') {
        addFrag(an.a);
        addFrag(an.b);
      }
    }
  }
  const { mol, index } = docToMol(doc, involved);
  const n = mol.atoms.length;
  const H = mol.atoms.map((_, i) => implicitH(mol, i));
  const N = mol.atoms.map((_, i) => nonBondingElectrons(mol, i));
  const N0 = [...N];
  // bond electrons keyed by "i,j" (i<j)
  const key = (i: number, j: number) => (i < j ? `${i},${j}` : `${j},${i}`);
  const be = new Map<string, number>();
  const bondStyle = new Map<string, { style: string; id: number }>();
  for (const b of mol.bonds) {
    const k = key(b.a, b.b);
    be.set(k, Math.round(b.order * 2));
    bondStyle.set(k, { style: b.style, id: b.id });
  }
  const before = new Map(be);
  const idx = (id: number) => index.get(id);
  const nearer = (pair: [number, number], target: number): number => {
    const t = mol.atoms[target];
    const [a, b] = pair;
    const da = Math.hypot(mol.atoms[a].x - t.x, mol.atoms[a].y - t.y);
    const db = Math.hypot(mol.atoms[b].x - t.x, mol.atoms[b].y - t.y);
    return da <= db ? a : b;
  };

  for (const c of arrows) {
    const k = c.electrons;
    // ── source ──
    let srcAtom = -1;
    let srcPair: [number, number] | null = null;
    if (c.from.type === 'atom') {
      const i = idx(c.from.id);
      if (i === undefined) continue;
      srcAtom = i;
      if (N[i] < k) warnings.push({ message: `${atomName(mol, i)} has no ${k === 2 ? 'lone pair' : 'non-bonding electron'} to donate`, atomIds: [mol.atoms[i].id] });
      N[i] -= k;
    } else if (c.from.type === 'bond') {
      const b = doc.bonds.get(c.from.id);
      const i = b ? idx(b.a) : undefined, j = b ? idx(b.b) : undefined;
      if (i === undefined || j === undefined) continue;
      srcPair = [i, j];
      const kk = key(i, j);
      be.set(kk, (be.get(kk) ?? 0) - k);
      if ((be.get(kk) ?? 0) < 0) warnings.push({ message: 'An arrow starts from a bond with no electrons left', atomIds: [mol.atoms[i].id, mol.atoms[j].id] });
    } else {
      warnings.push({ message: 'An electron-pushing arrow must start at an atom (lone pair) or a bond' });
      continue;
    }
    // ── target ──
    const addBondE = (i: number, j: number) => {
      const kk = key(i, j);
      be.set(kk, (be.get(kk) ?? 0) + k);
    };
    if (c.to.type === 'atom') {
      const t = idx(c.to.id);
      if (t === undefined) continue;
      if (srcAtom >= 0) {
        if (t === srcAtom) {
          N[t] += k;
          continue;
        }
        addBondE(srcAtom, t);
      } else if (srcPair) {
        if (t === srcPair[0] || t === srcPair[1]) N[t] += k;
        else addBondE(nearer(srcPair, t), t);
      }
    } else if (c.to.type === 'bond') {
      const b = doc.bonds.get(c.to.id);
      const i = b ? idx(b.a) : undefined, j = b ? idx(b.b) : undefined;
      if (i === undefined || j === undefined) continue;
      if (srcAtom >= 0 && srcAtom !== i && srcAtom !== j) {
        // lone pair pushed toward a remote bond: interpret as forming a bond to the nearer atom
        addBondE(nearer([i, j], srcAtom), srcAtom);
        warnings.push({ message: 'Arrow from a lone pair to a remote bond interpreted as bond formation to the nearer atom' });
      } else addBondE(i, j);
    } else if (c.to.type === 'between') {
      const i = idx(c.to.a), j = idx(c.to.b);
      if (i === undefined || j === undefined) continue;
      addBondE(i, j);
    } else {
      warnings.push({ message: 'An arrow ends in empty space — point it at an atom, a bond or between two atoms' });
      // electrons are lost → give them back to the source to keep counts consistent
      if (srcAtom >= 0) N[srcAtom] += k;
      else if (srcPair) addBondE(srcPair[0], srcPair[1]);
    }
  }

  // ── build product ──
  const product = new Mol();
  for (let i = 0; i < n; i++) product.atoms.push({ ...mol.atoms[i] });
  const pairs: [number, number, number][] = [];
  for (const [kk, e] of be) {
    const [i, j] = kk.split(',').map(Number);
    let ee = e;
    if (ee % 2 !== 0) {
      warnings.push({ message: `Odd number of electrons left between ${atomName(mol, i)} and ${atomName(mol, j)} — check the fishhook arrows`, atomIds: [mol.atoms[i].id, mol.atoms[j].id] });
      ee -= 1;
      N[i] += 1; // keep the electron on one atom (radical)
    }
    if (ee <= 0) continue;
    let order = ee / 2;
    if (order > 3) {
      warnings.push({ message: `Bond order ${order} between ${atomName(mol, i)} and ${atomName(mol, j)}`, atomIds: [mol.atoms[i].id, mol.atoms[j].id] });
      order = 3;
    }
    pairs.push([i, j, order]);
  }
  for (const [i, j, order] of pairs) {
    const kk = key(i, j);
    const old = bondStyle.get(kk);
    const bi = product.addBond(i, j, order, 'plain');
    if (old) {
      product.bonds[bi].id = old.id;
      const st = old.style;
      // keep stereo/visual styles on unchanged single bonds
      if (order === 1 && (before.get(kk) ?? 0) === 2) product.bonds[bi].style = st as any;
    } else product.bonds[bi].id = -1;
    // summary
    const prev = (before.get(kk) ?? 0) / 2;
    if (prev === 0) summary.push(`formed ${atomName(mol, i)}–${atomName(mol, j)} bond`);
    else if (order !== prev) summary.push(`${atomName(mol, i)}–${atomName(mol, j)}: bond order ${prev} → ${order}`);
  }
  for (const [kk, e] of before) {
    if (e > 0 && !pairs.some(([i, j]) => key(i, j) === kk)) {
      const [i, j] = kk.split(',').map(Number);
      summary.push(`broke ${atomName(mol, i)}–${atomName(mol, j)} bond`);
    }
  }
  product.invalidate();

  // formal charges, radicals and fixed hydrogen counts
  for (let i = 0; i < n; i++) {
    const a = product.atoms[i];
    if (a.abbrev || !element(a.el)) continue;
    const V = valenceElectrons(a.el);
    const B = bondOrderSum(product, i) + H[i];
    const Ni = Math.max(0, N[i]);
    if (N[i] < 0) warnings.push({ message: `${atomName(mol, i)} gave away more electrons than it had`, atomIds: [a.id] });
    a.charge = Math.round(V - Ni - B);
    const rad = Ni % 2;
    // keep an existing diradical/carbene designation when the electron count is unchanged
    if (rad) a.radical = 1;
    else if (N0[i] === N[i] && mol.atoms[i].radical === 2) a.radical = 2;
    else delete a.radical;
    // fix H counts so they don't drift
    delete a.hCount;
    if (implicitH(product, i) !== H[i]) a.hCount = H[i];
    // octet check
    const count = 2 * B + Ni;
    const lim = octetLimit(a.el);
    if (count > lim) warnings.push({ message: `${atomName(mol, i)} would have ${count} valence electrons (limit ${lim})`, atomIds: [a.id] });
  }

  // resonance: σ framework unchanged?
  const sig = (m: Map<string, number>) => [...m.entries()].filter(([, e]) => e > 0).map(([kk]) => kk).sort().join('|');
  const resonance = sig(before) === sig(new Map(pairs.map(([i, j, o]) => [key(i, j), o * 2])));
  if (!arrows.length) warnings.push({ message: 'No electron-pushing arrows to apply' });
  return { product, reactantAtomIds: [...involved], warnings, resonance, summary };
}

/** Atoms whose electron count violates the octet/duet rule in the current drawing. */
export function octetViolations(doc: ChemDoc): { atomId: number; count: number; limit: number }[] {
  const { mol } = docToMol(doc);
  const out: { atomId: number; count: number; limit: number }[] = [];
  for (let i = 0; i < mol.atoms.length; i++) {
    const a = mol.atoms[i];
    if (a.abbrev || !element(a.el)) continue;
    const e = element(a.el)!;
    if (e.group < 13 && a.el !== 'H') continue;
    const B = bondOrderSum(mol, i) + implicitH(mol, i);
    const count = 2 * B + nonBondingElectrons(mol, i);
    const lim = octetLimit(a.el);
    if (count > lim) out.push({ atomId: a.id, count, limit: lim });
  }
  return out;
}
