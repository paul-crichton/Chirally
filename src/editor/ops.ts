// Pure document operations used by the editor tools (no UI state).
import { ChemDoc, DocAtom, DocBond } from '../doc/types';
import {
  addAtom, addBond, adjacency, bondBetween, neighborsOf, removeAtom, removeBond, insertMol, pruneCurved, docToMol,
} from '../doc/document';
import { implicitH } from '../chem/valence';
import { freeAngles } from '../render/scene';
import { BondStyle, Mol } from '../chem/mol';
import { Pt, add, sub, mul, norm, len, dist, fromAngle, angleOf, perp, rotate, lerp } from '../render/geom';
import { element, ISOTOPE_ALIASES } from '../chem/elements';
import { ABBREVIATIONS, findAbbreviation, parseCondensedLabel, reverseLabel, labelOffsetForNet, labelDisplayCharge, chargeSuffix } from '../chem/abbreviations';

export const MERGE_TOL = 0.25;

export interface Selection {
  atoms: Set<number>;
  bonds: Set<number>;
  objects: Set<number>;
}

export const emptySelection = (): Selection => ({ atoms: new Set(), bonds: new Set(), objects: new Set() });

export function selectionSize(s: Selection): number {
  return s.atoms.size + s.bonds.size + s.objects.size;
}

/** Directions (unit vectors) of the bonds around an atom. */
export function bondDirs(doc: ChemDoc, atomId: number, adj = adjacency(doc)): Pt[] {
  const a = doc.atoms.get(atomId)!;
  return neighborsOf(doc, atomId, adj).map((n) => norm(sub(doc.atoms.get(n)!, a)));
}

/**
 * Ideal direction for a new bond from `atomId` (ChemDraw-like): zig-zag continuation for chains,
 * bisector of the largest gap otherwise, linear for triple bonds/allenes.
 */
export function idealBondDirection(doc: ChemDoc, atomId: number, newOrder = 1, adj = adjacency(doc)): Pt {
  const a = doc.atoms.get(atomId)!;
  const nbrs = neighborsOf(doc, atomId, adj);
  if (nbrs.length === 0) return fromAngle(-Math.PI / 6);
  if (nbrs.length === 1) {
    const n = doc.atoms.get(nbrs[0])!;
    const d = norm(sub(a, n)); // continuing direction
    const bond = bondBetween(doc, atomId, nbrs[0])!;
    // linear for triple bonds or cumulated double bonds
    if (bond.order === 3 || newOrder === 3 || (bond.order === 2 && newOrder === 2)) return d;
    // zig-zag: rotate ±60° from the continuing direction, choose the side that is trans to n's other neighbour
    const cands = [rotateVec(d, Math.PI / 3), rotateVec(d, -Math.PI / 3)];
    const nn = neighborsOf(doc, nbrs[0], adj).filter((x) => x !== atomId);
    if (nn.length) {
      const ref = doc.atoms.get(nn[0])!;
      // choose candidate farthest from ref (trans zigzag)
      const p0 = add(a, cands[0]);
      const p1 = add(a, cands[1]);
      return dist(p0, ref) >= dist(p1, ref) ? cands[0] : cands[1];
    }
    // prefer the upward/right-going candidate for a pleasing default
    const score = (v: Pt) => -v.y * 0.6 + v.x * 0.4;
    return score(cands[0]) >= score(cands[1]) ? cands[0] : cands[1];
  }
  // largest angular gap
  const angs = nbrs.map((nid) => angleOf(sub(doc.atoms.get(nid)!, a))).sort((p, q) => p - q);
  let best = 0, bestMid = 0;
  for (let i = 0; i < angs.length; i++) {
    const a0 = angs[i];
    const a1 = i + 1 < angs.length ? angs[i + 1] : angs[0] + 2 * Math.PI;
    if (a1 - a0 > best) {
      best = a1 - a0;
      bestMid = (a0 + a1) / 2;
    }
  }
  // with 3 neighbours in a "T" use the gap bisector; with two neighbours at 120°, gap bisector gives 120° again
  return fromAngle(bestMid);
}

function rotateVec(v: Pt, ang: number): Pt {
  const c = Math.cos(ang), s = Math.sin(ang);
  return { x: v.x * c - v.y * s, y: v.x * s + v.y * c };
}

/** Finds an atom within `tol` of p (excluding ids in `exclude`). */
export function atomNear(doc: ChemDoc, p: Pt, tol = MERGE_TOL, exclude?: Set<number>): DocAtom | null {
  let best: DocAtom | null = null;
  let bd = tol;
  for (const a of doc.atoms.values()) {
    if (exclude?.has(a.id)) continue;
    const d = Math.hypot(a.x - p.x, a.y - p.y);
    if (d < bd) {
      bd = d;
      best = a;
    }
  }
  return best;
}

/** Adds a bond from atom `fromId` to point p (merging with an existing atom there). */
export function bondToPoint(doc: ChemDoc, fromId: number, p: Pt, order: number, style: BondStyle, el = 'C'): { atom: DocAtom; bond: DocBond } {
  const existing = atomNear(doc, p, MERGE_TOL, new Set([fromId]));
  const target = existing ?? addAtom(doc, { el, x: p.x, y: p.y });
  let bond = bondBetween(doc, fromId, target.id);
  if (bond) {
    bond.order = order;
    bond.style = style;
    if (style === 'wedge' || style === 'hash' || style === 'hollow' || style === 'dative') {
      bond.a = fromId;
      bond.b = target.id;
    }
  } else bond = addBond(doc, fromId, target.id, order, style);
  return { atom: target, bond };
}

/** ChemDraw-style click on a bond with the plain bond tool: 1 → 2 → 3 → 1. */
export function cycleBondOrder(b: DocBond): void {
  if (b.style !== 'plain' && b.style !== 'bold' && b.style !== 'dashed') {
    b.style = 'plain';
    b.order = 1;
    return;
  }
  b.order = b.order === 1 || b.order === 1.5 ? 2 : b.order === 2 ? 3 : 1;
  delete b.dbPos;
}

/** Applies a bond type to an existing bond (ChemDraw semantics: the single-bond tool cycles 1 → 2 → 3). */
export function applyBondType(b: DocBond, order: number, style: BondStyle): void {
  if (order === 1 && style === 'plain') {
    if (b.style === 'plain' && b.order !== 0) cycleBondOrder(b);
    else {
      b.style = 'plain';
      b.order = 1;
    }
    return;
  }
  if (b.order === order && b.style === style) {
    if (style === 'wedge' || style === 'hash' || style === 'hollow' || style === 'dative') {
      const t = b.a;
      b.a = b.b;
      b.b = t;
      return;
    }
    if (order === 2 && style === 'plain') {
      // cycle double-bond placement: auto → right → left → center → auto
      const seq: DocBond['dbPos'][] = ['auto', 'right', 'left', 'center'];
      const k = seq.indexOf(b.dbPos ?? 'auto');
      b.dbPos = seq[(k + 1) % seq.length];
      return;
    }
    return;
  }
  b.order = order;
  b.style = style;
  if (order !== 2) delete b.dbPos;
}

/** Parses typed atom label text into atom properties. Returns null if not understood as chemistry. */
const PSEUDO_LABEL = /^(R\d*'*|R[a-z]|Ar|X\d*|Y|Z|Nu|Nuc|Nu:|B:|Base|E|LG|PG|Pg|M|L|Hal|A|Q|G\d*)$/;

/** "+", "2-", "−", "3+" → signed charge. */
function chargeOf(s: string): number {
  const sgn = /[\-−]/.test(s) ? -1 : 1;
  const num = s.replace(/[+\-−]/, '');
  return sgn * (num ? +num : 1);
}

/** The abbreviation a label names: exact spellings first, also reversed (MeS is SMe, not mesityl), then any case. */
function abbreviationKey(t: string): string | null {
  if (ABBREVIATIONS[t]) return t;
  const rev = reverseLabel(t);
  if (rev !== t && ABBREVIATIONS[rev]) return rev;
  return findAbbreviation(t) ?? (rev !== t ? findAbbreviation(rev) : null);
}

/**
 * A group label with a typed net charge (`net`), or with no sign (`net` undefined: the group as defined, so
 * "COO" is carboxylate). Abbreviation atoms store their charge as an offset from the group's own charge.
 */
function groupLabel(t: string, net: number | undefined, allowGeneric: boolean): Partial<DocAtom> | null {
  const key = abbreviationKey(t === 'NC' && (net ?? 0) < 0 ? 'CN' : t); // N≡C⁻ written either way round is cyanide
  if (key) return { abbrev: key, el: 'C', charge: net === undefined ? 0 : labelOffsetForNet(key, net), alias: undefined, hCount: undefined, isotope: undefined };
  if (parseCondensedLabel(t)) return { abbrev: t, el: 'C', charge: net ?? 0, alias: undefined, hCount: undefined, isotope: undefined };
  if (allowGeneric && (PSEUDO_LABEL.test(t) || t.length <= 6)) return { el: 'R', alias: t, abbrev: undefined, charge: net ?? 0, hCount: undefined, isotope: undefined };
  return null;
}

export function parseAtomLabel(text: string): Partial<DocAtom> | null {
  const t = text.trim().replace(/−/g, '-');
  if (!t) return null;
  // ions typed with their sign come first: "N3-" is azide, not nitride
  if (ABBREVIATIONS[t]?.category === 'ion') return { abbrev: t, el: 'C', charge: 0, alias: undefined, hCount: undefined, isotope: undefined };
  // element with optional H count and charge: "N", "NH2", "OH", "N+", "O-", "NH3+", "13C", "Fe2+", "D"
  const m = /^(\d+)?([A-Z][a-z]?)(?:H(\d*))?(\d*[+\-]|[+\-]\d*)?$/.exec(t);
  if (m) {
    let el = m[2];
    let isotope = m[1] ? +m[1] : undefined;
    if (ISOTOPE_ALIASES[el] && !m[1]) {
      isotope = ISOTOPE_ALIASES[el].isotope;
      el = ISOTOPE_ALIASES[el].el;
    }
    if (element(el)) {
      const out: Partial<DocAtom> = { el, isotope, abbrev: undefined, alias: undefined };
      if (m[3] !== undefined) out.hCount = m[3] === '' ? 1 : +m[3];
      out.charge = m[4] ? chargeOf(m[4]) : 0;
      // a halogen cation typed without H is the bare cation (Br⁺), not H₂Br⁺
      if (out.charge > 0 && m[3] === undefined && element(el)!.group === 17) out.hCount = 0;
      return out;
    }
  }
  // formulas written hydrogens first: H2O, H3O+, HO-, H2N-, HCl
  const f = /^H(\d*)([A-Z][a-z]?)(\d*[+\-]|[+\-]\d*)?$/.exec(t);
  if (f && f[2] !== 'H' && element(f[2])) {
    return { el: f[2], hCount: f[1] ? +f[1] : 1, charge: f[3] ? chargeOf(f[3]) : 0, abbrev: undefined, alias: undefined, isotope: undefined };
  }
  // a label spelled exactly (ions such as NO2+ or BF4- include their sign), or a condensed formula
  const exact = groupLabel(t, undefined, false);
  if (exact) return exact;
  // a trailing sign is the group's net charge: CN-, MeO-, COO-, Ph3P+, SMe2+, RS-, E+, R1-
  const s1 = /^(.+)([+\-])$/.exec(t);
  if (s1) {
    const g = groupLabel(s1[1], s1[2] === '+' ? 1 : -1, true);
    if (g) return g;
  }
  // a charge size written after the formula ("2-") only when nothing else fits
  const s2 = /^(.+?)(\d+)([+\-])$/.exec(t);
  if (s2) {
    const g = groupLabel(s2[1], chargeOf(s2[2] + s2[3]), true);
    if (g) return g;
  }
  // generic pseudo atoms: R, R1, R', Ar, X, Y, Z, Nu, E, LG, Pg…
  if (PSEUDO_LABEL.test(t) || t.length <= 6) {
    return { el: 'R', alias: t, abbrev: undefined, charge: 0, hCount: undefined };
  }
  return null;
}

/** The text that re-creates an atom's label when typed (used to start editing it): "OMe-", "E+", "COO-", "NH3+". */
export function labelText(a: DocAtom, bonded: boolean): string {
  const text = a.abbrev ?? a.alias;
  if (!text) return '';
  return text + chargeSuffix(a.abbrev ? labelDisplayCharge(a.abbrev, a.charge, bonded) : a.charge || 0);
}

export function setAtomLabel(doc: ChemDoc, atomId: number, text: string): boolean {
  const a = doc.atoms.get(atomId);
  if (!a) return false;
  const p = parseAtomLabel(text);
  if (!p) return false;
  for (const [k, v] of Object.entries(p)) {
    if (v === undefined) delete (a as any)[k];
    else (a as any)[k] = v;
  }
  return true;
}

/**
 * Draws the implicit hydrogens of the given atoms as H atoms (so curved arrows can start or end on a C–H
 * bond or its H). Returns the number of H atoms added.
 */
export function drawHydrogens(doc: ChemDoc, atomIds: Iterable<number>): number {
  let added = 0;
  const adj = adjacency(doc);
  for (const id of atomIds) {
    const a = doc.atoms.get(id);
    if (!a || a.abbrev || a.alias || !element(a.el) || a.el === 'H') continue;
    const nbrs = neighborsOf(doc, id, adj);
    const { mol, index } = docToMol(doc, [id, ...nbrs]);
    const k = implicitH(mol, index.get(id)!);
    if (k < 1) continue;
    const occupied = nbrs.map((n) => angleOf(sub(doc.atoms.get(n)!, a)));
    const hs: { atom: number; bond: number }[] = [];
    for (const ang of freeAngles(occupied, k, occupied.length ? occupied[0] + Math.PI : -Math.PI / 2)) {
      const h = addAtom(doc, { el: 'H', x: a.x + Math.cos(ang) * 0.8, y: a.y + Math.sin(ang) * 0.8 });
      hs.push({ atom: h.id, bond: addBond(doc, id, h.id).id });
      added++;
    }
    // every hydrogen is drawn now: an explicit count goes to 0 (it would otherwise add them again)
    if (a.hCount !== undefined) a.hCount = 0;
    // arrows that started or ended on the label's H now use the first drawn H (its bond as a source)
    for (const c of doc.curved.values()) {
      if (c.from.type === 'atom' && c.from.id === id && c.from.h) c.from = { type: 'bond', id: hs[0].bond };
      if (c.to.type === 'atom' && c.to.id === id && c.to.h) c.to = { type: 'atom', id: hs[0].atom };
    }
  }
  return added;
}

export function setElement(doc: ChemDoc, atomId: number, el: string): void {
  const a = doc.atoms.get(atomId);
  if (!a) return;
  a.el = el;
  delete a.abbrev;
  delete a.alias;
  delete a.hCount;
  delete a.isotope;
}

/** Deletes selected atoms/bonds/objects. Carbon atoms left without any bond by bond deletion are removed. */
export function deleteSelection(doc: ChemDoc, sel: Selection): void {
  const touched = new Set<number>();
  for (const bid of sel.bonds) {
    const b = doc.bonds.get(bid);
    if (!b) continue;
    touched.add(b.a);
    touched.add(b.b);
    removeBond(doc, bid);
  }
  for (const aid of sel.atoms) {
    for (const b of doc.bonds.values()) if (b.a === aid || b.b === aid) { touched.add(b.a); touched.add(b.b); }
    removeAtom(doc, aid);
  }
  for (const oid of sel.objects) {
    doc.arrows.delete(oid);
    doc.curved.delete(oid);
    doc.texts.delete(oid);
    doc.shapes.delete(oid);
  }
  const adj = adjacency(doc);
  for (const id of touched) {
    const a = doc.atoms.get(id);
    if (!a) continue;
    if ((adj.get(id)?.length ?? 0) === 0 && a.el === 'C' && !a.abbrev && !a.alias && !a.charge && !a.isotope) removeAtom(doc, id);
  }
  pruneCurved(doc);
}

/** Merges atom `dropId` into `keepId` (bonds are re-attached; duplicate bonds removed). */
export function mergeAtoms(doc: ChemDoc, keepId: number, dropId: number): void {
  if (keepId === dropId) return;
  for (const b of [...doc.bonds.values()]) {
    if (b.a !== dropId && b.b !== dropId) continue;
    const other = b.a === dropId ? b.b : b.a;
    if (other === keepId) {
      doc.bonds.delete(b.id);
      continue;
    }
    const dup = bondBetween(doc, keepId, other);
    if (dup) {
      if (b.order > dup.order) dup.order = b.order;
      doc.bonds.delete(b.id);
      continue;
    }
    if (b.a === dropId) b.a = keepId;
    else b.b = keepId;
  }
  // keep the more informative atom properties
  const keep = doc.atoms.get(keepId)!, drop = doc.atoms.get(dropId)!;
  if (keep.el === 'C' && !keep.abbrev && !keep.alias && (drop.el !== 'C' || drop.abbrev || drop.alias)) {
    keep.el = drop.el;
    if (drop.abbrev) keep.abbrev = drop.abbrev;
    if (drop.alias) keep.alias = drop.alias;
    keep.charge = drop.charge;
  }
  // re-anchor curved arrows
  for (const c of doc.curved.values()) {
    for (const an of [c.from, c.to]) {
      if (an.type === 'atom' && an.id === dropId) an.id = keepId;
      if (an.type === 'between') {
        if (an.a === dropId) an.a = keepId;
        if (an.b === dropId) an.b = keepId;
      }
    }
  }
  doc.atoms.delete(dropId);
  pruneCurved(doc);
}

/** After moving `moved` atoms, merge any that landed on a non-moved atom. Returns number of merges. */
export function fuseOverlaps(doc: ChemDoc, moved: Set<number>, tol = MERGE_TOL): number {
  let n = 0;
  for (const id of [...moved]) {
    const a = doc.atoms.get(id);
    if (!a) continue;
    const t = atomNear(doc, a, tol, moved);
    if (t) {
      mergeAtoms(doc, t.id, id);
      n++;
    }
  }
  return n;
}

/** Finds a pair (moved atom, target atom) that would merge if dropped now — used for drop highlights. */
export function pendingMerges(doc: ChemDoc, moved: Set<number>, tol = MERGE_TOL): [number, number][] {
  const out: [number, number][] = [];
  for (const id of moved) {
    const a = doc.atoms.get(id);
    if (!a) continue;
    const t = atomNear(doc, a, tol, moved);
    if (t) out.push([id, t.id]);
  }
  return out;
}

// ───────────── rings & templates ─────────────

export function ringRadius(n: number, L = 1): number {
  return L / (2 * Math.sin(Math.PI / n));
}

/** Places a new ring centred at c (first vertex at the top). */
export function placeRing(doc: ChemDoc, c: Pt, n: number, aromatic: boolean, startAngle = -Math.PI / 2): number[] {
  const R = ringRadius(n);
  const ids: number[] = [];
  for (let k = 0; k < n; k++) {
    const p = add(c, fromAngle(startAngle + (2 * Math.PI * k) / n, R));
    const ex = atomNear(doc, p, MERGE_TOL);
    ids.push(ex ? ex.id : addAtom(doc, { el: 'C', x: p.x, y: p.y }).id);
  }
  closeRingBonds(doc, ids, aromatic);
  return ids;
}

/** Creates ring bonds among consecutive ids, choosing alternating double bonds if aromatic. */
function closeRingBonds(doc: ChemDoc, ids: number[], aromatic: boolean, fixedFirst = false): void {
  const n = ids.length;
  const edges: [number, number][] = [];
  for (let k = 0; k < n; k++) edges.push([ids[k], ids[(k + 1) % n]]);
  for (const [a, b] of edges) if (a !== b && !bondBetween(doc, a, b)) addBond(doc, a, b, 1);
  if (!aromatic) return;
  const hasDouble = (id: number, chosen: Set<number>) => {
    if (chosen.has(id)) return true;
    for (const b of doc.bonds.values()) if ((b.a === id || b.b === id) && b.order >= 2) return true;
    return false;
  };
  let bestSet: [number, number][] = [];
  for (let off = 0; off < 2; off++) {
    const chosen = new Set<number>();
    const set: [number, number][] = [];
    for (let k = 0; k < n; k++) {
      const e = edges[(k + off) % n];
      if (fixedFirst && (k + off) % n === 0) continue;
      const bd = bondBetween(doc, e[0], e[1]);
      if (!bd || bd.order !== 1) continue;
      if (hasDouble(e[0], chosen) || hasDouble(e[1], chosen)) continue;
      chosen.add(e[0]);
      chosen.add(e[1]);
      set.push(e);
    }
    if (set.length > bestSet.length) bestSet = set;
  }
  for (const [a, b] of bestSet) bondBetween(doc, a, b)!.order = 2;
}

/** Fuses a new n-membered ring onto bond `bondId`, on the less crowded side. */
export function fuseRingOnBond(doc: ChemDoc, bondId: number, n: number, aromatic: boolean, sideHint?: Pt): number[] {
  const b = doc.bonds.get(bondId)!;
  const A = doc.atoms.get(b.a)!, B = doc.atoms.get(b.b)!;
  const L = dist(A, B) || 1;
  const mid = lerp(A, B, 0.5);
  const u = norm(sub(B, A));
  const pv = perp(u);
  let side: number;
  if (sideHint) {
    side = (sideHint.x - mid.x) * pv.x + (sideHint.y - mid.y) * pv.y >= 0 ? 1 : -1;
  } else {
    // count neighbours of A and B on each side
    let s = 0;
    const adj = adjacency(doc);
    for (const id of [...neighborsOf(doc, A.id, adj), ...neighborsOf(doc, B.id, adj)]) {
      if (id === A.id || id === B.id) continue;
      const p = doc.atoms.get(id)!;
      s += Math.sign((p.x - mid.x) * pv.x + (p.y - mid.y) * pv.y);
    }
    side = s > 0 ? -1 : 1;
  }
  const apo = L / (2 * Math.tan(Math.PI / n));
  const C = add(mid, mul(pv, apo * side));
  const R = ringRadius(n, L);
  const thA = angleOf(sub(A, C));
  const thB = angleOf(sub(B, C));
  let step = (2 * Math.PI) / n;
  // direction such that A → B is one step
  const diff = Math.atan2(Math.sin(thB - thA), Math.cos(thB - thA));
  if (diff < 0) step = -step;
  const ids = [A.id, B.id];
  for (let k = 2; k < n; k++) {
    const p = add(C, fromAngle(thA + step * k, R));
    const ex = atomNear(doc, p, MERGE_TOL * 1.2);
    ids.push(ex ? ex.id : addAtom(doc, { el: 'C', x: p.x, y: p.y }).id);
  }
  closeRingBonds(doc, ids, aromatic, false);
  return ids;
}

/** Attaches a ring to `atomId` through a new single bond (substituent), or spiro-fuses when `spiro`. */
export function attachRingToAtom(doc: ChemDoc, atomId: number, n: number, aromatic: boolean, spiro = false, dirOverride?: Pt): number[] {
  const a = doc.atoms.get(atomId)!;
  const adj = adjacency(doc);
  const deg = adj.get(atomId)?.length ?? 0;
  const d = dirOverride ?? (deg === 0 ? { x: 0, y: -1 } : idealBondDirection(doc, atomId, 1, adj));
  const R = ringRadius(n);
  if (deg === 0 || spiro) {
    // the atom itself becomes a ring vertex
    const C = add(a, mul(d, R));
    const th0 = angleOf(sub(a, C));
    const ids = [atomId];
    for (let k = 1; k < n; k++) {
      const p = add(C, fromAngle(th0 + (2 * Math.PI * k) / n, R));
      const ex = atomNear(doc, p, MERGE_TOL);
      ids.push(ex ? ex.id : addAtom(doc, { el: 'C', x: p.x, y: p.y }).id);
    }
    closeRingBonds(doc, ids, aromatic);
    return ids;
  }
  const first = add(a, d);
  const C = add(first, mul(d, R));
  const th0 = angleOf(sub(first, C));
  const ids: number[] = [];
  for (let k = 0; k < n; k++) {
    const p = add(C, fromAngle(th0 + (2 * Math.PI * k) / n, R));
    const ex = atomNear(doc, p, MERGE_TOL);
    ids.push(ex ? ex.id : addAtom(doc, { el: 'C', x: p.x, y: p.y }).id);
  }
  addBond(doc, atomId, ids[0], 1);
  closeRingBonds(doc, ids, aromatic);
  return ids;
}

/** Cyclohexane chair (bond length ≈ 1, opposite bonds parallel), centred at the origin. */
export const CHAIR: Pt[] = [
  { x: -1.2, y: 0.4 }, { x: -0.75, y: -0.45 }, { x: 0.2, y: -0.2 }, { x: 1.2, y: -0.4 }, { x: 0.75, y: 0.45 }, { x: -0.2, y: 0.2 },
];

export function placeChair(doc: ChemDoc, c: Pt): number[] {
  const ids = CHAIR.map((p) => {
    const q = { x: c.x + p.x, y: c.y + p.y };
    const ex = atomNear(doc, q, MERGE_TOL);
    return ex ? ex.id : addAtom(doc, { el: 'C', x: q.x, y: q.y }).id;
  });
  for (let k = 0; k < 6; k++) if (!bondBetween(doc, ids[k], ids[(k + 1) % 6])) addBond(doc, ids[k], ids[(k + 1) % 6], 1);
  return ids;
}

/** Inserts a template molecule centred at p. Returns new atom ids. */
export function insertTemplateAt(doc: ChemDoc, mol: Mol, p: Pt): number[] {
  const bb = mol.bbox();
  const cx = (bb.minX + bb.maxX) / 2, cy = (bb.minY + bb.maxY) / 2;
  const { atomIds } = insertMol(doc, mol, p.x - cx, p.y - cy);
  fuseOverlaps(doc, new Set(atomIds));
  return atomIds.filter((id) => doc.atoms.has(id));
}

/**
 * Attaches template atom `attach` (default 0) to `atomId` via a new bond, rotating the template so it
 * points along the ideal direction.
 */
export function attachTemplate(doc: ChemDoc, mol: Mol, atomId: number, attach = 0): number[] {
  const a = doc.atoms.get(atomId)!;
  const d = idealBondDirection(doc, atomId);
  const m = mol.clone();
  // template orientation: vector from attach atom to template centroid
  const cen = m.atoms.reduce((s, at) => add(s, at), { x: 0, y: 0 });
  const cc = mul(cen, 1 / m.atoms.length);
  const at = m.atoms[attach];
  const v = m.atoms.length > 1 ? norm(sub(cc, at)) : d;
  const ang = angleOf(d) - angleOf(v);
  const target = add(a, d);
  for (const x of m.atoms) {
    const r = rotate(x, at, ang);
    x.x = r.x - at.x + target.x;
    x.y = r.y - at.y + target.y;
  }
  const { atomIds } = insertMol(doc, m);
  addBond(doc, atomId, atomIds[attach], 1);
  fuseOverlaps(doc, new Set(atomIds));
  return atomIds.filter((id) => doc.atoms.has(id));
}

// ───────────── transforms ─────────────

export function selectionCenter(doc: ChemDoc, sel: Selection): Pt | null {
  const pts: Pt[] = [];
  for (const id of selectedAtomIds(doc, sel)) pts.push(doc.atoms.get(id)!);
  for (const id of sel.objects) {
    const o = doc.arrows.get(id) ?? doc.shapes.get(id);
    if (o) {
      pts.push({ x: o.x1, y: o.y1 }, { x: o.x2, y: o.y2 });
      continue;
    }
    const t = doc.texts.get(id);
    if (t) pts.push(t);
  }
  if (!pts.length) return null;
  let x1 = Infinity, y1 = Infinity, x2 = -Infinity, y2 = -Infinity;
  for (const p of pts) {
    x1 = Math.min(x1, p.x); x2 = Math.max(x2, p.x);
    y1 = Math.min(y1, p.y); y2 = Math.max(y2, p.y);
  }
  return { x: (x1 + x2) / 2, y: (y1 + y2) / 2 };
}

/** Atoms affected by a selection: selected atoms plus endpoints of selected bonds. */
export function selectedAtomIds(doc: ChemDoc, sel: Selection): Set<number> {
  const s = new Set<number>(sel.atoms);
  for (const bid of sel.bonds) {
    const b = doc.bonds.get(bid);
    if (b) {
      s.add(b.a);
      s.add(b.b);
    }
  }
  return s;
}

export function transformSelection(doc: ChemDoc, sel: Selection, fn: (p: Pt) => Pt): void {
  for (const id of selectedAtomIds(doc, sel)) {
    const a = doc.atoms.get(id)!;
    const p = fn(a);
    a.x = p.x;
    a.y = p.y;
  }
  for (const id of sel.objects) {
    const ar = doc.arrows.get(id) ?? doc.shapes.get(id);
    if (ar) {
      const p1 = fn({ x: ar.x1, y: ar.y1 }), p2 = fn({ x: ar.x2, y: ar.y2 });
      ar.x1 = p1.x; ar.y1 = p1.y; ar.x2 = p2.x; ar.y2 = p2.y;
      continue;
    }
    const t = doc.texts.get(id);
    if (t) {
      const p = fn(t);
      t.x = p.x;
      t.y = p.y;
      continue;
    }
    const c = doc.curved.get(id);
    if (c) {
      for (const an of [c.from, c.to]) {
        if (an.type === 'point') {
          const p = fn(an);
          an.x = p.x;
          an.y = p.y;
        }
      }
    }
  }
  // free-point curved arrows attached to moved atoms are moved implicitly through anchors
}

export function moveSelection(doc: ChemDoc, sel: Selection, dx: number, dy: number): void {
  transformSelection(doc, sel, (p) => ({ x: p.x + dx, y: p.y + dy }));
}

export function rotateSelection(doc: ChemDoc, sel: Selection, center: Pt, ang: number): void {
  transformSelection(doc, sel, (p) => rotate(p, center, ang));
}

/**
 * Flips the selection about a vertical ('h') or horizontal ('v') axis through its centre.
 * With keepStereo (default) wedges/hashes are swapped so the molecule is unchanged (a 3D rotation);
 * without, the result is the mirror image (enantiomer).
 */
export function flipSelection(doc: ChemDoc, sel: Selection, axis: 'h' | 'v', keepStereo = true): void {
  const c = selectionCenter(doc, sel);
  if (!c) return;
  transformSelection(doc, sel, (p) => (axis === 'h' ? { x: 2 * c.x - p.x, y: p.y } : { x: p.x, y: 2 * c.y - p.y }));
  const atoms = selectedAtomIds(doc, sel);
  for (const b of doc.bonds.values()) {
    if (!atoms.has(b.a) || !atoms.has(b.b)) continue;
    if (keepStereo) {
      if (b.style === 'wedge') b.style = 'hash';
      else if (b.style === 'hash') b.style = 'wedge';
    }
    if (b.dbPos === 'left') b.dbPos = 'right';
    else if (b.dbPos === 'right') b.dbPos = 'left';
  }
  // curved arrows inside the selection: mirror their bulge
  for (const cv of doc.curved.values()) {
    const inSel = (an: typeof cv.from) =>
      (an.type === 'atom' && atoms.has(an.id)) ||
      (an.type === 'bond' && sel.bonds.has(an.id)) ||
      (an.type === 'between' && atoms.has(an.a) && atoms.has(an.b)) ||
      (an.type === 'point' && sel.objects.has(cv.id));
    if (inSel(cv.from) && inSel(cv.to)) {
      cv.c1.h = -cv.c1.h;
      cv.c2.h = -cv.c2.h;
    }
  }
}

/** Scales the selection about its centre. */
export function scaleSelection(doc: ChemDoc, sel: Selection, s: number): void {
  const c = selectionCenter(doc, sel);
  if (!c) return;
  transformSelection(doc, sel, (p) => ({ x: c.x + (p.x - c.x) * s, y: c.y + (p.y - c.y) * s }));
}

/** Zig-zag chain points from p0 in direction `ang`, `n` bonds. `firstUp` decides the initial zig. */
export function chainPoints(p0: Pt, ang: number, n: number, firstUp: boolean): Pt[] {
  const pts: Pt[] = [];
  let p = p0;
  for (let k = 0; k < n; k++) {
    const s = (k % 2 === 0) === firstUp ? -1 : 1;
    p = add(p, fromAngle(ang + (s * Math.PI) / 6, 1));
    pts.push(p);
  }
  return pts;
}

/** Normalises bond lengths of a pasted/imported structure: returns median bond length. */
export function medianBondLength(doc: ChemDoc, atomIds?: Set<number>): number {
  const ls: number[] = [];
  for (const b of doc.bonds.values()) {
    if (atomIds && !atomIds.has(b.a)) continue;
    const A = doc.atoms.get(b.a)!, B = doc.atoms.get(b.b)!;
    ls.push(dist(A, B));
  }
  if (!ls.length) return 1;
  ls.sort((p, q) => p - q);
  return ls[Math.floor(ls.length / 2)];
}

export { len };
