// Arrow-pushing engine: applies curved (electron-flow) arrows to produce the next intermediate,
// with full electron bookkeeping (lone pairs, bonds, formal charges, radicals) and octet checks.
import { ChemDoc, CurvedArrowObj, Anchor, ArrowObj } from './types';
import { docToMol, adjacency, fragmentOf, docBounds } from './document';
import { Mol } from '../chem/mol';
import { clean2D } from '../chem/clean2d';
import { implicitH, nonBondingElectrons, bondOrderSum, octetLimit } from '../chem/valence';
import { valenceElectrons, element } from '../chem/elements';

export interface MechanismWarning {
  message: string;
  atomIds?: number[];
  /** 'error': the product is not a valid structure (it is previewed, not inserted); 'warning': inserted, but check it. */
  level: 'error' | 'warning';
}

export interface MechanismResult {
  /** Product structure; atoms keep the document ids of the atoms they came from. */
  product: Mol;
  /** Document atom ids of the reactant fragments used. */
  reactantAtomIds: number[];
  warnings: MechanismWarning[];
  /** true for a valid step that only moves π / lone-pair electrons (σ-connectivity unchanged → ↔ arrow). */
  resonance: boolean;
  /** No error-level warnings: the product is a valid next intermediate. */
  ok: boolean;
  /** The arrows changed at least one bond or electron count. */
  changed: boolean;
  /** Human-readable summary of what happened, e.g. "formed C–O, broke C–Br". */
  summary: string[];
}

function atomName(mol: Mol, i: number): string {
  const a = mol.atoms[i];
  return a.abbrev ?? a.alias ?? a.el;
}

/** Groups curved arrows into independent "steps": arrows sharing fragments belong together. */
export interface ArrowGroup {
  arrows: number[];
  atoms: number[];
  maxX: number;
  /** Every arrow in the group already produced a step whose reaction/resonance arrow is still in the drawing. */
  applied: boolean;
}

/** True when the curved arrow was applied and the arrow drawn for that step still exists. */
export function isApplied(doc: ChemDoc, c: CurvedArrowObj): boolean {
  return c.step !== undefined && doc.arrows.has(c.step);
}

export function arrowGroups(doc: ChemDoc): ArrowGroup[] {
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
    for (const a of g.atoms) maxX = Math.max(maxX, doc.atoms.get(a)?.x ?? -Infinity);
    return { arrows: g.arrows, atoms: [...g.atoms], maxX, applied: g.arrows.every((id) => isApplied(doc, doc.curved.get(id)!)) };
  });
}

/**
 * Applies the given curved arrows. Each arrow moves 2 electrons (or 1 for fishhooks) from its source
 * (atom lone pair or bond) to its target (atom, bond, or the space between two atoms).
 */
export function applyArrows(doc: ChemDoc, arrowIds: number[]): MechanismResult {
  const warnings: MechanismWarning[] = [];
  const error = (message: string, atomIds?: number[]) => warnings.push({ message, atomIds, level: 'error' });
  const warn = (message: string, atomIds?: number[]) => warnings.push({ message, atomIds, level: 'warning' });
  const summary: string[] = [];
  const arrows = arrowIds.map((id) => doc.curved.get(id)).filter(Boolean) as CurvedArrowObj[];
  const missing = () => error('An arrow is attached to an atom or bond that no longer exists');
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
    const addBondE = (i: number, j: number) => {
      const kk = key(i, j);
      be.set(kk, (be.get(kk) ?? 0) + k);
    };
    // ── source ──
    let srcAtom = -1;
    let srcPair: [number, number] | null = null;
    if (c.from.type === 'atom') {
      const i = idx(c.from.id);
      if (i === undefined) {
        missing();
        continue;
      }
      srcAtom = i;
      if (N[i] < k) error(`${atomName(mol, i)} has no ${k === 2 ? 'lone pair' : 'non-bonding electron'} to donate`, [mol.atoms[i].id]);
      N[i] -= k;
    } else if (c.from.type === 'bond') {
      const b = doc.bonds.get(c.from.id);
      const i = b ? idx(b.a) : undefined, j = b ? idx(b.b) : undefined;
      if (i === undefined || j === undefined) {
        missing();
        continue;
      }
      srcPair = [i, j];
      const kk = key(i, j);
      be.set(kk, (be.get(kk) ?? 0) - k);
      if ((be.get(kk) ?? 0) < 0) error('An arrow starts from a bond with no electrons left', [mol.atoms[i].id, mol.atoms[j].id]);
    } else {
      error('An electron-pushing arrow must start at an atom (lone pair) or a bond');
      continue;
    }
    // ── target ──
    // electrons go back to where they came from when the target is unusable
    const giveBack = () => {
      if (srcAtom >= 0) N[srcAtom] += k;
      else if (srcPair) addBondE(srcPair[0], srcPair[1]);
    };
    if (c.to.type === 'atom') {
      const t = idx(c.to.id);
      if (t === undefined) {
        missing();
        giveBack();
        continue;
      }
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
      if (i === undefined || j === undefined) {
        missing();
        giveBack();
        continue;
      }
      if (srcAtom >= 0 && srcAtom !== i && srcAtom !== j) {
        // lone pair pushed toward a remote bond: interpret as forming a bond to the nearer atom
        addBondE(nearer([i, j], srcAtom), srcAtom);
        warn('Arrow from a lone pair to a remote bond interpreted as bond formation to the nearer atom');
      } else addBondE(i, j);
    } else if (c.to.type === 'between') {
      const i = idx(c.to.a), j = idx(c.to.b);
      if (i === undefined || j === undefined) {
        missing();
        giveBack();
        continue;
      }
      addBondE(i, j);
    } else {
      error('An arrow ends in empty space — point it at an atom, a bond or between two atoms');
      // electrons are lost → give them back to the source to keep counts consistent
      giveBack();
    }
  }

  const changed =
    N.some((v, i) => v !== N0[i]) || [...new Set([...before.keys(), ...be.keys()])].some((kk) => Math.max(0, be.get(kk) ?? 0) !== (before.get(kk) ?? 0));

  // ── build product ──
  const product = new Mol();
  for (let i = 0; i < n; i++) product.atoms.push({ ...mol.atoms[i] });
  const pairs: [number, number, number][] = [];
  for (const [kk, e] of be) {
    const [i, j] = kk.split(',').map(Number);
    let ee = e;
    if (ee % 2 !== 0) {
      error(`Odd number of electrons left between ${atomName(mol, i)} and ${atomName(mol, j)} — check the fishhook arrows`, [mol.atoms[i].id, mol.atoms[j].id]);
      ee -= 1;
      N[i] += 1; // keep the electron on one atom (radical)
    }
    if (ee <= 0) continue;
    let order = ee / 2;
    if (order > 3) {
      error(`Bond order ${order} between ${atomName(mol, i)} and ${atomName(mol, j)}`, [mol.atoms[i].id, mol.atoms[j].id]);
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
    if (N[i] < 0) error(`${atomName(mol, i)} gave away more electrons than it had`, [a.id]);
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
    if (count > lim) error(`${atomName(mol, i)} would have ${count} valence electrons (limit ${lim})`, [a.id]);
  }

  // resonance: σ framework unchanged?
  const sig = (m: Map<string, number>) => [...m.entries()].filter(([, e]) => e > 0).map(([kk]) => kk).sort().join('|');
  if (!arrows.length) error('No electron-pushing arrows to apply');
  const ok = !warnings.some((w) => w.level === 'error');
  // a failed or empty step is never a resonance structure
  const resonance = ok && changed && sig(before) === sig(new Map(pairs.map(([i, j, o]) => [key(i, j), o * 2])));
  return { product, reactantAtomIds: [...involved], warnings, resonance, ok, changed, summary };
}

export interface StepPlacement {
  /** The product, moved to where it goes in the drawing (atoms keep the ids of the atoms they came from). */
  mol: Mol;
  /** The reaction (→) or resonance (↔) arrow drawn between reactants and product. */
  arrow: Omit<ArrowObj, 'id' | 'type'>;
}

/** Where a step's product and arrow go: to the right of the reactants and of anything else drawn in that row. */
export function placeStep(doc: ChemDoc, r: MechanismResult): StepPlacement | null {
  const ids = new Set(r.reactantAtomIds);
  const b = docBounds(doc, ids);
  if (!b || !r.product.atoms.length) return null;
  const mol = r.product.clone();
  const inBand = (y1: number, y2: number) => y2 >= b.minY - 1.5 && y1 <= b.maxY + 1.5;
  let bandMax = b.maxX;
  for (const a of doc.atoms.values()) if (!ids.has(a.id) && a.x > b.minX && inBand(a.y, a.y)) bandMax = Math.max(bandMax, a.x);
  for (const o of doc.arrows.values()) if (Math.max(o.x1, o.x2) > b.minX && inBand(Math.min(o.y1, o.y2), Math.max(o.y1, o.y2))) bandMax = Math.max(bandMax, o.x1, o.x2);
  for (const t of doc.texts.values()) if (t.x > b.maxX && inBand(t.y, t.y)) bandMax = Math.max(bandMax, t.x + 1);
  const startX = bandMax + 0.8;
  const arrowGap = 3.4;
  // tidy bonds that were created between separate fragments
  const formedLong = mol.bonds.some((bd) => {
    const A = mol.atoms[bd.a], B = mol.atoms[bd.b];
    return Math.hypot(A.x - B.x, A.y - B.y) > 1.6;
  });
  if (formedLong) {
    try {
      clean2D(mol);
    } catch {
      /* keep raw geometry */
    }
  }
  const pb = mol.bbox();
  const midY = (b.minY + b.maxY) / 2;
  mol.translate(startX + arrowGap - pb.minX + 0.2, midY - (pb.minY + pb.maxY) / 2);
  return { mol, arrow: { kind: r.resonance ? 'resonance' : 'reaction', x1: startX, y1: midY, x2: startX + arrowGap - 0.8, y2: midY } };
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
