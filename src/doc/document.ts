import { Mol, TetraSpec } from '../chem/mol';
import { ChemDoc, DocAtom, DocBond, STYLE_PRESETS, DocStyle, ArrowObj, CurvedArrowObj, TextObj, ShapeObj } from './types';

export const FILE_FORMAT = 'chirally';
export const FILE_VERSION = 1;

export function createDoc(style: DocStyle = STYLE_PRESETS.Chirally): ChemDoc {
  return {
    version: 1,
    atoms: new Map(),
    bonds: new Map(),
    arrows: new Map(),
    curved: new Map(),
    texts: new Map(),
    shapes: new Map(),
    nextId: 1,
    style: { ...style },
    meta: { title: 'Untitled', created: new Date().toISOString() },
  };
}

export function cloneDoc(d: ChemDoc): ChemDoc {
  const cloneMap = <T extends object>(m: Map<number, T>): Map<number, T> => {
    const out = new Map<number, T>();
    for (const [k, v] of m) out.set(k, structuredCloneLite(v));
    return out;
  };
  return {
    version: 1,
    atoms: cloneMap(d.atoms),
    bonds: cloneMap(d.bonds),
    arrows: cloneMap(d.arrows),
    curved: cloneMap(d.curved),
    texts: cloneMap(d.texts),
    shapes: cloneMap(d.shapes),
    nextId: d.nextId,
    style: { ...d.style },
    meta: { ...d.meta },
  };
}

function structuredCloneLite<T>(v: T): T {
  // objects in the document are shallow except curved-arrow anchors / control points
  const o: any = { ...(v as any) };
  for (const k of Object.keys(o)) {
    const x = o[k];
    if (x && typeof x === 'object') o[k] = Array.isArray(x) ? [...x] : { ...x };
  }
  return o;
}

export function newId(doc: ChemDoc): number {
  return doc.nextId++;
}

export function addAtom(doc: ChemDoc, a: Partial<DocAtom> & { el: string; x: number; y: number }): DocAtom {
  const atom: DocAtom = { charge: 0, ...a, id: newId(doc) } as DocAtom;
  doc.atoms.set(atom.id, atom);
  return atom;
}

export function addBond(doc: ChemDoc, a: number, b: number, order = 1, style: DocBond['style'] = 'plain'): DocBond {
  const bond: DocBond = { id: newId(doc), a, b, order, style };
  doc.bonds.set(bond.id, bond);
  return bond;
}

export function bondBetween(doc: ChemDoc, a: number, b: number): DocBond | undefined {
  for (const bd of doc.bonds.values()) {
    if ((bd.a === a && bd.b === b) || (bd.a === b && bd.b === a)) return bd;
  }
  return undefined;
}

/** Adjacency map atomId → bond ids. Build once per operation for speed. */
export function adjacency(doc: ChemDoc): Map<number, number[]> {
  const adj = new Map<number, number[]>();
  for (const id of doc.atoms.keys()) adj.set(id, []);
  for (const b of doc.bonds.values()) {
    adj.get(b.a)?.push(b.id);
    adj.get(b.b)?.push(b.id);
  }
  return adj;
}

export function neighborsOf(doc: ChemDoc, atomId: number, adj = adjacency(doc)): number[] {
  return (adj.get(atomId) ?? []).map((bid) => {
    const b = doc.bonds.get(bid)!;
    return b.a === atomId ? b.b : b.a;
  });
}

export function removeAtom(doc: ChemDoc, id: number): void {
  doc.atoms.delete(id);
  for (const [bid, b] of [...doc.bonds]) if (b.a === id || b.b === id) doc.bonds.delete(bid);
  pruneCurved(doc);
}

export function removeBond(doc: ChemDoc, id: number): void {
  doc.bonds.delete(id);
  pruneCurved(doc);
}

/** Removes curved arrows whose anchors no longer exist. */
export function pruneCurved(doc: ChemDoc): void {
  const ok = (an: CurvedArrowObj['from']) => {
    if (an.type === 'atom') return doc.atoms.has(an.id);
    if (an.type === 'bond') return doc.bonds.has(an.id);
    if (an.type === 'between') return doc.atoms.has(an.a) && doc.atoms.has(an.b);
    return true;
  };
  for (const [id, c] of [...doc.curved]) if (!ok(c.from) || !ok(c.to)) doc.curved.delete(id);
}

/** Connected fragments as lists of atom ids. */
export function fragments(doc: ChemDoc, adj = adjacency(doc)): number[][] {
  const seen = new Set<number>();
  const out: number[][] = [];
  for (const id of doc.atoms.keys()) {
    if (seen.has(id)) continue;
    const comp: number[] = [];
    const st = [id];
    seen.add(id);
    while (st.length) {
      const v = st.pop()!;
      comp.push(v);
      for (const w of neighborsOf(doc, v, adj)) {
        if (!seen.has(w)) {
          seen.add(w);
          st.push(w);
        }
      }
    }
    out.push(comp);
  }
  return out;
}

export function fragmentOf(doc: ChemDoc, atomId: number, adj = adjacency(doc)): number[] {
  const seen = new Set<number>([atomId]);
  const st = [atomId];
  while (st.length) {
    const v = st.pop()!;
    for (const w of neighborsOf(doc, v, adj)) {
      if (!seen.has(w)) {
        seen.add(w);
        st.push(w);
      }
    }
  }
  return [...seen];
}

/**
 * Builds a Mol from document atoms (all, or the given ids). Mol atoms keep their document `id`;
 * Mol bond `id`s are document bond ids. `index` maps atom id → Mol index.
 */
export function docToMol(doc: ChemDoc, atomIds?: Iterable<number>): { mol: Mol; index: Map<number, number> } {
  const mol = new Mol();
  const index = new Map<number, number>();
  const ids = atomIds ? [...atomIds] : [...doc.atoms.keys()];
  for (const id of ids) {
    const a = doc.atoms.get(id);
    if (!a) continue;
    index.set(id, mol.atoms.length);
    mol.atoms.push({ ...a });
  }
  for (const b of doc.bonds.values()) {
    const ia = index.get(b.a);
    const ib = index.get(b.b);
    if (ia === undefined || ib === undefined) continue;
    mol.bonds.push({ ...b, a: ia, b: ib });
  }
  mol.invalidate();
  return { mol, index };
}

/**
 * Inserts a Mol into the document (coordinates taken as-is). Returns the new atom ids in Mol order.
 * Stereo specs are not stored in documents (stereo lives in wedges/geometry).
 */
export function insertMol(doc: ChemDoc, mol: Mol, dx = 0, dy = 0): { atomIds: number[]; bondIds: number[] } {
  const atomIds: number[] = [];
  const bondIds: number[] = [];
  for (const a of mol.atoms) {
    const { id: _id, ...rest } = a;
    void _id;
    const na = addAtom(doc, { ...rest, x: a.x + dx, y: a.y + dy });
    atomIds.push(na.id);
  }
  for (const b of mol.bonds) {
    const nb = addBond(doc, atomIds[b.a], atomIds[b.b], b.order, b.style);
    if (b.dbPos) nb.dbPos = b.dbPos;
    if (b.color) nb.color = b.color;
    bondIds.push(nb.id);
  }
  return { atomIds, bondIds };
}

export function isDocEmpty(doc: ChemDoc): boolean {
  return !doc.atoms.size && !doc.arrows.size && !doc.texts.size && !doc.shapes.size && !doc.curved.size;
}

// ───────────── serialization (native .chirally JSON) ─────────────

export interface SerializedDoc {
  format: typeof FILE_FORMAT;
  version: number;
  meta: ChemDoc['meta'];
  style: DocStyle;
  atoms: DocAtom[];
  bonds: DocBond[];
  arrows: ArrowObj[];
  curved: CurvedArrowObj[];
  texts: TextObj[];
  shapes: ShapeObj[];
}

const round = (v: number) => Math.round(v * 10000) / 10000;

export function serializeDoc(doc: ChemDoc): SerializedDoc {
  return {
    format: FILE_FORMAT,
    version: FILE_VERSION,
    meta: { ...doc.meta, modified: new Date().toISOString() },
    style: { ...doc.style },
    atoms: [...doc.atoms.values()].map((a) => ({ ...a, x: round(a.x), y: round(a.y) })),
    bonds: [...doc.bonds.values()].map((b) => ({ ...b })),
    arrows: [...doc.arrows.values()],
    curved: [...doc.curved.values()],
    texts: [...doc.texts.values()],
    shapes: [...doc.shapes.values()],
  };
}

export function deserializeDoc(data: unknown): ChemDoc {
  const d = data as Partial<SerializedDoc>;
  if (!d || d.format !== FILE_FORMAT) throw new Error('Not a Chirally document');
  const doc = createDoc({ ...STYLE_PRESETS.Chirally, ...(d.style ?? {}) });
  doc.meta = { title: 'Untitled', ...(d.meta ?? {}) };
  let maxId = 0;
  const put = <T extends { id: number }>(m: Map<number, T>, list?: T[]) => {
    for (const o of list ?? []) {
      m.set(o.id, { ...o });
      maxId = Math.max(maxId, o.id);
    }
  };
  put(doc.atoms, d.atoms);
  put(doc.bonds, d.bonds);
  put(doc.arrows, d.arrows);
  put(doc.curved, d.curved);
  put(doc.texts, d.texts);
  put(doc.shapes, d.shapes);
  // drop dangling bonds
  for (const [id, b] of [...doc.bonds]) if (!doc.atoms.has(b.a) || !doc.atoms.has(b.b)) doc.bonds.delete(id);
  pruneCurved(doc);
  doc.nextId = maxId + 1;
  return doc;
}

export function docToJSON(doc: ChemDoc, pretty = false): string {
  return JSON.stringify(serializeDoc(doc), null, pretty ? 1 : 0);
}

export function docFromJSON(text: string): ChemDoc {
  return deserializeDoc(JSON.parse(text));
}

/** Bounding box of all (or selected) chemistry + objects. */
export function docBounds(doc: ChemDoc, ids?: Set<number>): { minX: number; minY: number; maxX: number; maxY: number } | null {
  let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
  const inc = (x: number, y: number) => {
    minX = Math.min(minX, x); maxX = Math.max(maxX, x);
    minY = Math.min(minY, y); maxY = Math.max(maxY, y);
  };
  const want = (id: number) => !ids || ids.has(id);
  for (const a of doc.atoms.values()) if (want(a.id)) inc(a.x, a.y);
  for (const o of doc.arrows.values()) if (want(o.id)) { inc(o.x1, o.y1); inc(o.x2, o.y2); }
  for (const o of doc.texts.values()) if (want(o.id)) inc(o.x, o.y);
  for (const o of doc.shapes.values()) if (want(o.id)) { inc(o.x1, o.y1); inc(o.x2, o.y2); }
  if (minX === Infinity) return null;
  return { minX, minY, maxX, maxY };
}

/** Stereo specs are not persisted; helper for code that wants to keep TetraSpec typing handy. */
export type { TetraSpec };
