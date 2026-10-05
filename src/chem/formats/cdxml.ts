// ChemDraw XML (CDXML) reader/writer for whole canvas documents.
//
// CDXML coordinates are points with y DOWN (same orientation as Chirally); one model unit is the
// document's BondLength (default 14.4 pt). Reading turns fragments into atoms/bonds, <t> into text
// objects, <arrow>/<graphic>/<curve> into arrows, curved arrows and shapes. Unknown elements are
// skipped. Writing produces a CDXML file ChemDraw can open (nodes, bonds with wedge Display
// attributes, nickname nodes with expanded inner fragments, text, arrows, curved arrows, shapes).
import { Atom, BondStyle, Mol } from '../mol';
import { atomicNumber } from '../elements';
import { implicitH } from '../valence';
import { abbreviationMol, expandAbbreviations, chargeSuffix, labelDisplayCharge } from '../abbreviations';
import { ChemDoc, DocBond, ArrowObj, ArrowKind, CurvedArrowObj, TextObj, ShapeObj, ShapeKind, DocStyle, STYLE_PRESETS, Anchor } from '../../doc/types';
import { createDoc, addAtom, addBond, newId, docToMol, fragments, adjacency, neighborsOf, docBounds, implicitHydrogenPoint } from '../../doc/document';
import { parseXml, XmlElement, childElements, firstChild, localName, textContent, attrs, escapeText } from './xml';
import { FormatError, guard, elementByNumber, setAbbrev, dropDefaultHCounts, splitLabelCharge } from './common';

const DEFAULT_BOND_LENGTH_PT = 14.4;

type Pt = { x: number; y: number };

const num = (s: string | undefined): number => (s === undefined || s.trim() === '' ? NaN : Number(s));
const int = (s: string | undefined): number => {
  const v = num(s);
  return Number.isFinite(v) ? Math.round(v) : NaN;
};
const nums = (s: string | undefined): number[] =>
  s === undefined ? [] : s.trim().split(/[\s,]+/).filter(Boolean).map(Number);

function parsePoint(s: string | undefined): Pt | null {
  const v = nums(s);
  return v.length >= 2 && Number.isFinite(v[0]) && Number.isFinite(v[1]) ? { x: v[0], y: v[1] } : null;
}

function parseBox(s: string | undefined): [number, number, number, number] | null {
  const v = nums(s);
  return v.length >= 4 && v.slice(0, 4).every(Number.isFinite) ? [v[0], v[1], v[2], v[3]] : null;
}

const normBox = (b: [number, number, number, number]): [number, number, number, number] => [
  Math.min(b[0], b[2]), Math.min(b[1], b[3]), Math.max(b[0], b[2]), Math.max(b[1], b[3]),
];

function median(v: number[]): number {
  const s = [...v].sort((a, b) => a - b);
  const m = s.length >> 1;
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
}

// ───────────────────────────── reading ─────────────────────────────

interface Collected {
  fragments: XmlElement[];
  texts: XmlElement[];
  arrows: XmlElement[];
  graphics: XmlElement[];
  curves: XmlElement[];
  steps: XmlElement[];
}

const LEAF_ELEMENTS = new Set(['colortable', 'fonttable', 'n', 'b', 's', 'objecttag', 'annotation', 'chemicalproperty', 'represent']);

/** Gathers drawable objects from pages/groups/schemes (nested fragments inside nodes excluded). */
function collect(el: XmlElement, out: Collected): void {
  for (const c of childElements(el)) {
    const name = localName(c.name);
    if (name === 'fragment') out.fragments.push(c);
    else if (name === 't') out.texts.push(c);
    else if (name === 'arrow') out.arrows.push(c);
    else if (name === 'graphic') out.graphics.push(c);
    else if (name === 'curve') out.curves.push(c);
    else if (name === 'step') out.steps.push(c);
    else if (!LEAF_ELEMENTS.has(name)) collect(c, out);
  }
}

interface ReadCtx {
  doc: ChemDoc;
  unit: number;
  labelSize: number;
  colors: string[];
  nodeIds: Map<string, number>;
  bondKeys: Set<string>;
  textIds: Map<string, TextObj>;
  arrowIds: Map<string, ArrowObj>;
}

/** Reads a CDXML document into a Chirally document. */
export function readCDXML(xml: string): ChemDoc {
  return guard('CDXML', () => {
    const root = parseXml(xml, 'CDXML');
    if (localName(root.name) !== 'CDXML') throw new FormatError('CDXML', `root element is <${root.name}>, expected <CDXML>`);
    const items: Collected = { fragments: [], texts: [], arrows: [], graphics: [], curves: [], steps: [] };
    collect(root, items);
    const unit = chooseUnit(root, items.fragments);
    const doc = createDoc(styleFromRoot(root, unit));
    doc.meta.title = root.attrs.Name?.replace(/\.cdxml$/i, '') || 'Untitled';
    const ctx: ReadCtx = {
      doc, unit,
      labelSize: num(root.attrs.LabelSize) > 0 ? num(root.attrs.LabelSize) : 10,
      colors: parseColorTable(root),
      nodeIds: new Map(), bondKeys: new Set(), textIds: new Map(), arrowIds: new Map(),
    };
    for (const f of items.fragments) readFragment(ctx, f, items);
    fixHydrogens(doc);
    for (const t of items.texts) readText(ctx, t);
    for (const a of items.arrows) readArrow(ctx, a);
    for (const g of items.graphics) readGraphic(ctx, g);
    for (const c of items.curves) readCurve(ctx, c);
    for (const s of items.steps) applyStep(ctx, s);
    return doc;
  });
}

/** BondLength attribute, unless the drawn bonds clearly use a different length. */
function chooseUnit(root: XmlElement, frags: XmlElement[]): number {
  const lengths: number[] = [];
  const visit = (f: XmlElement) => {
    const pos = new Map<string, Pt>();
    for (const n of childElements(f, 'n')) {
      const p = parsePoint(n.attrs.p);
      if (p && n.attrs.id) pos.set(n.attrs.id, p);
    }
    for (const b of childElements(f, 'b')) {
      const p = pos.get(b.attrs.B), q = pos.get(b.attrs.E);
      if (p && q) {
        const d = Math.hypot(p.x - q.x, p.y - q.y);
        if (d > 1e-3) lengths.push(d);
      }
    }
    for (const sub of childElements(f, 'fragment')) visit(sub);
  };
  frags.forEach(visit);
  const bl = num(root.attrs.BondLength);
  const med = lengths.length ? median(lengths) : NaN;
  if (bl > 0) return !lengths.length || Math.abs(med - bl) / bl < 0.25 ? bl : med;
  return lengths.length ? med : DEFAULT_BOND_LENGTH_PT;
}

function styleFromRoot(root: XmlElement, unit: number): DocStyle {
  const s: DocStyle = { ...STYLE_PRESETS['ACS 1996'], name: 'CDXML', bondLengthPt: unit };
  const a = root.attrs;
  const rel = (v: string | undefined) => (num(v) > 0 ? num(v) / unit : undefined);
  s.lineWidth = rel(a.LineWidth) ?? s.lineWidth;
  s.boldWidth = rel(a.BoldWidth) ?? s.boldWidth;
  s.hashSpacing = rel(a.HashSpacing) ?? s.hashSpacing;
  s.labelMargin = rel(a.MarginWidth) ?? s.labelMargin;
  s.fontSize = rel(a.LabelSize) ?? s.fontSize;
  if (num(a.BondSpacing) > 0) s.bondSpacing = num(a.BondSpacing) / 100;
  if (a.ShowNonTerminalCarbonLabels === 'yes') s.showCarbons = 'all';
  else if (a.ShowTerminalCarbonLabels === 'yes') s.showCarbons = 'terminal';
  if (a.HideImplicitHydrogens === 'yes') s.showImplicitH = false;
  return s;
}

/** Colour table: indices 0/1 are black/white, table entries start at index 2. */
function parseColorTable(root: XmlElement): string[] {
  const ct = firstChild(root, 'colortable');
  const hex = (v: string | undefined) => Math.max(0, Math.min(255, Math.round((num(v) || 0) * 255))).toString(16).padStart(2, '0');
  const list = ['#000000', '#ffffff'];
  if (ct) for (const c of childElements(ct, 'color')) list.push('#' + hex(c.attrs.r) + hex(c.attrs.g) + hex(c.attrs.b));
  return list;
}

function colorOf(ctx: ReadCtx, attr: string | undefined): string | undefined {
  const i = int(attr);
  if (!Number.isFinite(i)) return undefined;
  const c = ctx.colors[i];
  return c && c !== '#000000' ? c : undefined;
}

interface Run {
  text: string;
  face: number;
  size: number;
  color?: string;
}

function readRuns(ctx: ReadCtx | null, t: XmlElement): Run[] {
  const ss = childElements(t, 's');
  if (!ss.length) return [{ text: textContent(t), face: 0, size: NaN }];
  return ss.map((s) => ({ text: textContent(s), face: int(s.attrs.face) || 0, size: num(s.attrs.size), color: ctx ? colorOf(ctx, s.attrs.color) : undefined }));
}

const PSEUDO_NODE_TYPES = new Set(['GenericNickname', 'Unspecified', 'Unknown', 'AnonymousAlternativeGroup', 'NamedAlternativeGroup', 'LinkNode']);

function readFragment(ctx: ReadCtx, frag: XmlElement, items: Collected): void {
  const { doc, unit } = ctx;
  const pendingNet = new Map<number, number>();
  for (const n of childElements(frag, 'n')) {
    const type = n.attrs.NodeType ?? 'Element';
    if (type === 'ExternalConnectionPoint' || type === 'MultiAttachment' || type === 'VariableAttachment') continue;
    const bb = parseBox(n.attrs.BoundingBox);
    const p = parsePoint(n.attrs.p) ?? (bb ? { x: (bb[0] + bb[2]) / 2, y: (bb[1] + bb[3]) / 2 } : { x: 0, y: 0 });
    const tEl = firstChild(n, 't');
    const label = tEl ? readRuns(null, tEl).map((r) => r.text).join('').trim() : '';
    const z = int(n.attrs.Element);
    const atom: Partial<Atom> & { el: string; x: number; y: number } = {
      el: Number.isFinite(z) ? elementByNumber(z) ?? 'C' : 'C',
      x: p.x / unit,
      y: p.y / unit,
      charge: int(n.attrs.Charge) || 0,
    };
    const iso = int(n.attrs.Isotope);
    if (iso > 0) atom.isotope = iso;
    const rad = n.attrs.Radical;
    if (rad === 'Doublet') atom.radical = 1;
    else if (rad === 'Singlet' || rad === 'Triplet') atom.radical = 2;
    const nh = int(n.attrs.NumHydrogens);
    if (nh >= 0) atom.hCount = nh;
    const color = colorOf(ctx, n.attrs.color);
    if (color) atom.color = color;
    // a label written with its charge ("OMe-", "PPh3+", "E+"): the sign is the group's net charge
    const signed = splitLabelCharge(label);
    let net: number | undefined;
    if (type === 'Nickname' || type === 'Fragment') {
      // right-justified labels are drawn reversed ("MeO" for OMe attached on its right)
      const reversed = tEl?.attrs.LabelJustification === 'Right' || tEl?.attrs.LabelAlignment === 'Right';
      if (label && setAbbrev(atom as Atom, label, reversed)) {
        // the label as written is a known group or ion
      } else if (signed && setAbbrev(atom as Atom, signed[0], reversed)) net = signed[1];
      else makePseudo(atom, label || '?');
    } else if (PSEUDO_NODE_TYPES.has(type)) {
      const name = n.attrs.GenericNickname || (signed ? signed[0] : label) || 'R';
      makePseudo(atom, name);
      if (signed && label !== name) atom.charge = signed[1];
    }
    const a = addAtom(doc, atom);
    if (n.attrs.id) ctx.nodeIds.set(n.attrs.id, a.id);
    if (net !== undefined) pendingNet.set(a.id, net);
  }
  for (const b of childElements(frag, 'b')) {
    let a = ctx.nodeIds.get(b.attrs.B), e = ctx.nodeIds.get(b.attrs.E);
    if (a === undefined || e === undefined || a === e) continue;
    const key = a < e ? `${a},${e}` : `${e},${a}`;
    if (ctx.bondKeys.has(key)) continue;
    ctx.bondKeys.add(key);
    let order = 1;
    let style: BondStyle = 'plain';
    const ord = (b.attrs.Order ?? '1').trim().split(/\s+/)[0];
    if (ord === 'dative') style = 'dative';
    else if (ord === 'hydrogen') { order = 0; style = 'hbond'; }
    else if (ord === 'ionic') order = 0;
    else {
      const v = Number(ord);
      if (v === 1.5) order = 1.5;
      else if (v === 0.5) { order = 1; style = 'dashed'; }
      else if (Number.isFinite(v)) order = Math.max(0, Math.min(4, Math.floor(v)));
    }
    if (style === 'plain') {
      let swap = false;
      switch (b.attrs.Display) {
        case 'WedgeBegin': style = 'wedge'; break;
        case 'WedgeEnd': style = 'wedge'; swap = true; break;
        case 'WedgedHashBegin': style = 'hash'; break;
        case 'WedgedHashEnd': style = 'hash'; swap = true; break;
        case 'HollowWedgeBegin': style = 'hollow'; break;
        case 'HollowWedgeEnd': style = 'hollow'; swap = true; break;
        case 'Wavy': case 'WavyWedgeBegin': case 'WavyWedgeEnd': style = 'wavy'; break;
        case 'Bold': style = 'bold'; break;
        case 'Dash': case 'Hash': case 'Dot': case 'DashDot': style = 'dashed'; break;
      }
      if (order !== 1 && (style === 'wedge' || style === 'hash' || style === 'hollow' || style === 'wavy')) style = 'plain';
      if (swap) [a, e] = [e, a];
    }
    const bond = addBond(doc, a, e, order, style);
    const dp = b.attrs.DoublePosition;
    if (dp === 'Left' || dp === 'Right' || dp === 'Center') bond.dbPos = dp.toLowerCase() as 'left' | 'right' | 'center';
    const color = colorOf(ctx, b.attrs.color);
    if (color) bond.color = color;
  }
  // the offset behind a signed label depends on whether it is bonded (a free PPh3 is the neutral molecule)
  if (pendingNet.size) {
    const adj = adjacency(doc);
    for (const [aid, q] of pendingNet) {
      const at = doc.atoms.get(aid)!;
      at.charge = q - labelDisplayCharge(at.abbrev!, 0, neighborsOf(doc, aid, adj).length > 0);
    }
  }
  for (const sub of childElements(frag, 'fragment')) readFragment(ctx, sub, items);
  for (const t of childElements(frag, 't')) items.texts.push(t);
}

function makePseudo(atom: Partial<Atom>, label: string): void {
  atom.el = 'R';
  atom.alias = label;
  atom.charge = 0;
  delete atom.abbrev;
  delete atom.hCount;
  delete atom.isotope;
}

/** Keeps NumHydrogens only where it differs from the computed implicit-H count. */
function fixHydrogens(doc: ChemDoc): void {
  if (![...doc.atoms.values()].some((a) => a.hCount !== undefined)) return;
  const { mol } = docToMol(doc);
  dropDefaultHCounts(mol);
  for (const a of mol.atoms) if (a.hCount === undefined) delete doc.atoms.get(a.id)!.hCount;
}

const wrapMarkup = (t: string, c: string) => (t.length === 1 ? c + t : `${c}{${t}}`);

/** Formula-style text → markup: digits after element symbols subscript, trailing signs superscript. */
function formulaToMarkup(t: string): string {
  let charge = '';
  const m = /^(.*[A-Za-z0-9)\]])([+\-−]+)$/.exec(t);
  if (m) {
    t = m[1];
    charge = m[2];
  }
  t = t.replace(/([A-Za-z)\]])(\d+)/g, (_, a: string, d: string) => a + wrapMarkup(d, '_'));
  return t + (charge ? wrapMarkup(charge, '^') : '');
}

/** CDXML style runs → Chirally text markup (+ whole-text flags). */
function runsToMarkup(runs: Run[]): { text: string; formula?: boolean; bold?: boolean; italic?: boolean } {
  const visible = runs.filter((r) => r.text.trim());
  const join = runs.map((r) => r.text).join('').replace(/\r\n?/g, '\n');
  if (visible.length && visible.every((r) => (r.face & 96) === 96)) return { text: join, formula: true };
  const allBold = visible.length > 0 && visible.every((r) => r.face & 1);
  const allItalic = visible.length > 0 && visible.every((r) => r.face & 2);
  let s = '';
  for (const r of runs) {
    let t = r.text.replace(/\r\n?/g, '\n');
    if (t.trim()) {
      if ((r.face & 96) === 96) t = formulaToMarkup(t);
      else if (r.face & 32) t = wrapMarkup(t, '_');
      else if (r.face & 64) t = wrapMarkup(t, '^');
      if (r.face & 2 && !allItalic) t = `*${t}*`;
      if (r.face & 1 && !allBold) t = `**${t}**`;
    }
    s += t;
  }
  const out: { text: string; bold?: boolean; italic?: boolean } = { text: s };
  if (allBold) out.bold = true;
  if (allItalic) out.italic = true;
  return out;
}

function readText(ctx: ReadCtx, t: XmlElement): void {
  const runs = readRuns(ctx, t);
  const m = runsToMarkup(runs);
  if (!m.text.trim()) return;
  const bb = parseBox(t.attrs.BoundingBox);
  const p = parsePoint(t.attrs.p) ?? (bb ? { x: Math.min(bb[0], bb[2]), y: Math.max(bb[1], bb[3]) } : null);
  if (!p) return;
  const obj: TextObj = { id: newId(ctx.doc), type: 'text', x: p.x / ctx.unit, y: p.y / ctx.unit, text: m.text };
  if (m.formula) obj.formula = true;
  if (m.bold) obj.bold = true;
  if (m.italic) obj.italic = true;
  const just = t.attrs.Justification ?? t.attrs.CaptionJustification;
  if (just === 'Center') obj.align = 'center';
  else if (just === 'Right') obj.align = 'right';
  const size = runs.find((r) => Number.isFinite(r.size))?.size;
  if (size && Math.abs(size / ctx.labelSize - 1) > 0.02) obj.size = Math.round((size / ctx.labelSize) * 1000) / 1000;
  const color = runs.find((r) => r.color)?.color ?? colorOf(ctx, t.attrs.color);
  if (color) obj.color = color;
  ctx.doc.texts.set(obj.id, obj);
  if (t.attrs.id) ctx.textIds.set(t.attrs.id, obj);
}

const hasHead = (v: string | undefined) => !!v && v !== 'None';

function straightKind(a: Record<string, string>): { kind: ArrowKind; reverse: boolean } {
  const h = hasHead(a.ArrowheadHead), t = hasHead(a.ArrowheadTail);
  const spacing = num(a.ArrowShaftSpacing) > 0;
  const reverse = !h && t;
  if (hasHead(a.NoGo)) return { kind: 'noGo', reverse };
  if (a.ArrowheadType === 'Hollow' || (spacing && h !== t)) return { kind: 'retro', reverse };
  if (spacing || (a.ArrowheadHead?.startsWith('Half') && a.ArrowheadTail?.startsWith('Half'))) {
    const ratio = num(a.ArrowEquilibriumRatio);
    return { kind: Number.isFinite(ratio) && ratio !== 1 && ratio !== 100 ? 'unbalancedEq' : 'equilibrium', reverse: false };
  }
  if (h && t) return { kind: 'resonance', reverse: false };
  if (!h && !t) return { kind: 'line', reverse: false };
  return { kind: /dash/i.test(a.LineType ?? '') ? 'dashed' : 'reaction', reverse };
}

function addArrow(ctx: ReadCtx, tail: Pt, head: Pt, kind: ArrowKind, el: XmlElement): ArrowObj {
  const u = ctx.unit;
  const obj: ArrowObj = { id: newId(ctx.doc), type: 'arrow', kind, x1: tail.x / u, y1: tail.y / u, x2: head.x / u, y2: head.y / u };
  const color = colorOf(ctx, el.attrs.color);
  if (color) obj.color = color;
  ctx.doc.arrows.set(obj.id, obj);
  if (el.attrs.id) ctx.arrowIds.set(el.attrs.id, obj);
  return obj;
}

function readArrow(ctx: ReadCtx, a: XmlElement): void {
  let head = parsePoint(a.attrs.Head3D);
  let tail = parsePoint(a.attrs.Tail3D);
  const center = parsePoint(a.attrs.Center3D);
  const ang = num(a.attrs.AngularSize);
  if (!head || !tail) {
    const bb = parseBox(a.attrs.BoundingBox);
    if (!bb) return;
    const [x1, y1, x2, y2] = normBox(bb);
    if (x2 - x1 >= y2 - y1) {
      tail = { x: x1, y: (y1 + y2) / 2 };
      head = { x: x2, y: (y1 + y2) / 2 };
    } else {
      tail = { x: (x1 + x2) / 2, y: y1 };
      head = { x: (x1 + x2) / 2, y: y2 };
    }
  }
  if (center && Number.isFinite(ang) && Math.abs(ang) > 0.5) {
    readArcArrow(ctx, a, tail, head, center, Math.abs(ang));
    return;
  }
  const { kind, reverse } = straightKind(a.attrs);
  if (reverse) addArrow(ctx, head, tail, kind, a);
  else addArrow(ctx, tail, head, kind, a);
}

/** Converts a cubic Bézier (absolute points) to a curved-arrow object with point anchors. */
function addCurved(ctx: ReadCtx, S: Pt, c1: Pt, c2: Pt, E: Pt, electrons: 1 | 2, el: XmlElement): void {
  const u = ctx.unit;
  const sx = S.x / u, sy = S.y / u, ex = E.x / u, ey = E.y / u;
  const len = Math.hypot(ex - sx, ey - sy);
  if (len < 1e-6) return;
  const ux = (ex - sx) / len, uy = (ey - sy) / len;
  const vx = -uy, vy = ux; // u rotated +90°
  const rel = (c: Pt) => {
    const dx = c.x / u - sx, dy = c.y / u - sy;
    return { t: (dx * ux + dy * uy) / len, h: (dx * vx + dy * vy) / len };
  };
  const obj: CurvedArrowObj = {
    id: newId(ctx.doc), type: 'curved', electrons,
    from: { type: 'point', x: sx, y: sy }, to: { type: 'point', x: ex, y: ey },
    c1: rel(c1), c2: rel(c2),
  };
  const color = colorOf(ctx, el.attrs.color);
  if (color) obj.color = color;
  ctx.doc.curved.set(obj.id, obj);
}

const electronsOf = (head: string | undefined): 1 | 2 => (head && head.startsWith('Half') ? 1 : 2);

/**
 * Arc arrow (ChemDraw >= 13 curved arrow): tail -> head along the ellipse given by Center3D,
 * MajorAxisEnd3D and MinorAxisEnd3D, sweeping |AngularSize| degrees. The sweep direction is taken
 * from the geometry (which way round reaches the head); for ambiguous half-circles the sign of
 * AngularSize in the major->minor axis frame decides.
 */
function readArcArrow(ctx: ReadCtx, a: XmlElement, tail: Pt, head: Pt, C: Pt, angleDeg: number): void {
  let T = tail, H = head;
  let headType = a.attrs.ArrowheadHead;
  let sign = num(a.attrs.AngularSize) < 0 ? -1 : 1;
  if (!hasHead(headType) && hasHead(a.attrs.ArrowheadTail)) {
    [T, H] = [H, T];
    headType = a.attrs.ArrowheadTail;
    sign = -sign;
  }
  const theta = (Math.min(angleDeg, 359) * Math.PI) / 180;
  const maj = parsePoint(a.attrs.MajorAxisEnd3D), min = parsePoint(a.attrs.MinorAxisEnd3D);
  // ellipse frame: E(phi) = C + cos(phi) * ax + sin(phi) * bx
  let ax = { x: T.x - C.x, y: T.y - C.y };
  let bx = { x: -ax.y, y: ax.x };
  if (maj && min) {
    const a2 = { x: maj.x - C.x, y: maj.y - C.y }, b2 = { x: min.x - C.x, y: min.y - C.y };
    const det = a2.x * b2.y - a2.y * b2.x;
    if (Math.abs(det) > 1e-6 * Math.hypot(a2.x, a2.y) * Math.hypot(b2.x, b2.y)) {
      ax = a2;
      bx = b2;
    }
  }
  const det = ax.x * bx.y - ax.y * bx.x;
  if (Math.abs(det) < 1e-9) return;
  const param = (p: Pt) => {
    const vx = p.x - C.x, vy = p.y - C.y;
    return Math.atan2((ax.x * vy - ax.y * vx) / det, (vx * bx.y - vy * bx.x) / det);
  };
  const TWO_PI = 2 * Math.PI;
  const norm = (v: number) => ((v % TWO_PI) + TWO_PI) % TWO_PI;
  const pT = param(T), pH = param(H);
  const errPos = Math.abs(norm(pH - pT) - theta), errNeg = Math.abs(norm(pT - pH) - theta);
  const dir = Math.abs(errPos - errNeg) < 0.05 ? sign : errPos < errNeg ? 1 : -1;
  const delta = dir * theta;
  const dE = (phi: number) => ({ x: -Math.sin(phi) * ax.x + Math.cos(phi) * bx.x, y: -Math.sin(phi) * ax.y + Math.cos(phi) * bx.y });
  const k = (4 / 3) * Math.tan(delta / 4);
  const d1 = dE(pT), d2 = dE(pT + delta);
  const c1 = { x: T.x + k * d1.x, y: T.y + k * d1.y };
  const c2 = { x: H.x - k * d2.x, y: H.y - k * d2.y };
  addCurved(ctx, T, c1, c2, H, electronsOf(headType), a);
}

function readCurve(ctx: ReadCtx, c: XmlElement): void {
  const v = nums(c.attrs.CurvePoints);
  const P: Pt[] = [];
  for (let i = 0; i + 1 < v.length; i += 2) if (Number.isFinite(v[i]) && Number.isFinite(v[i + 1])) P.push({ x: v[i], y: v[i + 1] });
  const n = P.length;
  if (n < 4) return;
  const h = c.attrs.ArrowheadHead, t = c.attrs.ArrowheadTail;
  if (!hasHead(h) && !hasHead(t)) return; // decorative curve
  let S: Pt, c1: Pt, c2: Pt, E: Pt;
  if (n % 3 === 0) [S, c1, c2, E] = [P[1], P[2], P[n - 3], P[n - 2]]; // in/anchor/out triples
  else [S, c1, c2, E] = [P[0], P[1], P[n - 2], P[n - 1]];
  if (!hasHead(h)) addCurved(ctx, E, c2, c1, S, electronsOf(t), c);
  else addCurved(ctx, S, c1, c2, E, electronsOf(h), c);
}

const SYMBOL_TEXT: Record<string, string> = {
  Plus: '+', CirclePlus: '+', Minus: '−', CircleMinus: '−', Radical: '•', RadicalCation: '•+', RadicalAnion: '•−',
  Dagger: '†', DoubleDagger: '‡', Electron: '•',
};

function readGraphic(ctx: ReadCtx, g: XmlElement): void {
  if (g.attrs.SupersededBy) return; // an <arrow> element carries this object
  const u = ctx.unit;
  const type = g.attrs.GraphicType;
  const bb = parseBox(g.attrs.BoundingBox);
  const color = colorOf(ctx, g.attrs.color);
  const shape = (kind: ShapeKind, box: [number, number, number, number], extra: Partial<ShapeObj> = {}) => {
    const obj: ShapeObj = { id: newId(ctx.doc), type: 'shape', kind, x1: box[0] / u, y1: box[1] / u, x2: box[2] / u, y2: box[3] / u, ...extra };
    if (color) obj.color = color;
    ctx.doc.shapes.set(obj.id, obj);
  };
  const ellipseBox = (): { box: [number, number, number, number]; angle?: number } | null => {
    const c = parsePoint(g.attrs.Center3D), maj = parsePoint(g.attrs.MajorAxisEnd3D), min = parsePoint(g.attrs.MinorAxisEnd3D);
    if (c && maj) {
      const rx = Math.hypot(maj.x - c.x, maj.y - c.y);
      const ry = min ? Math.hypot(min.x - c.x, min.y - c.y) : rx;
      const ang = Math.atan2(maj.y - c.y, maj.x - c.x);
      const box: [number, number, number, number] = [c.x - rx, c.y - ry, c.x + rx, c.y + ry];
      return Math.abs(Math.sin(ang)) > 1e-3 ? { box, angle: ang } : { box };
    }
    return bb ? { box: normBox(bb) } : null;
  };
  switch (type) {
    case 'Line': {
      if (!bb) return;
      const head = { x: bb[0], y: bb[1] }, tail = { x: bb[2], y: bb[3] };
      const at = g.attrs.ArrowType;
      if (!at || at === 'NoHead') {
        shape('line', bb, /dash/i.test(g.attrs.LineType ?? '') ? { dashed: true } : {});
        return;
      }
      const kinds: Record<string, ArrowKind> = { Resonance: 'resonance', Equilibrium: 'equilibrium', Hollow: 'retro', RetroSynthetic: 'retro' };
      const kind = kinds[at] ?? (/dash/i.test(g.attrs.LineType ?? '') ? 'dashed' : 'reaction');
      addArrow(ctx, tail, head, kind, g);
      return;
    }
    case 'Rectangle': {
      if (!bb) return;
      const rt = g.attrs.RectangleType ?? '';
      shape(/RoundEdge/.test(rt) ? 'roundRect' : 'rect', normBox(bb), /Dash/.test(rt) ? { dashed: true } : {});
      return;
    }
    case 'Oval': {
      const e = ellipseBox();
      if (e) shape('ellipse', e.box, { ...(e.angle !== undefined ? { angle: e.angle } : {}), ...(/Dash/.test(g.attrs.OvalType ?? '') ? { dashed: true } : {}) });
      return;
    }
    case 'Orbital': {
      const e = ellipseBox();
      if (e) shape(/^s/i.test(g.attrs.OrbitalType ?? '') ? 'orbitalS' : 'orbitalP', e.box, e.angle !== undefined ? { angle: e.angle } : {});
      return;
    }
    case 'Bracket': {
      if (!bb) return;
      const bt = g.attrs.BracketType ?? 'SquarePair';
      shape(/Round/.test(bt) ? 'paren' : /Curly/.test(bt) ? 'brace' : 'bracket', normBox(bb));
      return;
    }
    case 'Symbol': {
      const text = SYMBOL_TEXT[g.attrs.SymbolType ?? ''];
      if (!text || !bb) return;
      const nb = normBox(bb);
      const x = (nb[0] + nb[2]) / 2 / u, y = (nb[1] + nb[3]) / 2 / u;
      // charge/radical symbols attached to atoms duplicate the atom's own properties
      for (const a of ctx.doc.atoms.values()) {
        if (Math.hypot(a.x - x, a.y - y) > 0.8) continue;
        if ((text === '+' && a.charge > 0) || (text === '−' && a.charge < 0) || (text === '•' && a.radical)) return;
      }
      const obj: TextObj = { id: newId(ctx.doc), type: 'text', x, y, text, align: 'center' };
      if (color) obj.color = color;
      ctx.doc.texts.set(obj.id, obj);
      return;
    }
  }
}

/** <step>: texts listed above/below a reaction arrow become the arrow's conditions. */
function applyStep(ctx: ReadCtx, s: XmlElement): void {
  const ids = (k: string) => (s.attrs[k] ?? '').trim().split(/\s+/).filter(Boolean);
  const arrow = ids('ReactionStepArrows').map((id) => ctx.arrowIds.get(id)).find(Boolean);
  if (!arrow) return;
  const take = (list: string[]) => {
    const parts: string[] = [];
    for (const id of list) {
      const t = ctx.textIds.get(id);
      if (!t || !ctx.doc.texts.has(t.id)) continue;
      parts.push(t.text);
      ctx.doc.texts.delete(t.id);
    }
    return parts.join('\n');
  };
  const above = take(ids('ReactionStepObjectsAboveArrow'));
  const below = take(ids('ReactionStepObjectsBelowArrow'));
  if (above) arrow.above = above;
  if (below) arrow.below = below;
}

// ───────────────────────────── writing ─────────────────────────────

const f2 = (v: number) => {
  const s = (Math.round(v * 100) / 100).toFixed(2).replace(/\.?0+$/, '');
  return s === '-0' ? '0' : s;
};

/** Parses Chirally text markup into CDXML style runs (face bits: 1 bold, 2 italic, 32 sub, 64 sup). */
function markupToRuns(text: string, baseFace: number): { text: string; face: number }[] {
  const runs: { text: string; face: number }[] = [];
  let bold = false, italic = false;
  const push = (t: string, face: number) => {
    if (!t) return;
    const last = runs[runs.length - 1];
    if (last && last.face === face) last.text += t;
    else runs.push({ text: t, face });
  };
  const cur = () => baseFace | (bold ? 1 : 0) | (italic ? 2 : 0);
  let i = 0;
  while (i < text.length) {
    const c = text[i];
    if (c === '*' && text[i + 1] === '*') {
      bold = !bold;
      i += 2;
    } else if (c === '*') {
      italic = !italic;
      i++;
    } else if ((c === '_' || c === '^') && i + 1 < text.length) {
      let body: string;
      if (text[i + 1] === '{') {
        const end = text.indexOf('}', i + 2);
        body = end < 0 ? text.slice(i + 2) : text.slice(i + 2, end);
        i = end < 0 ? text.length : end + 1;
      } else {
        body = text[i + 1];
        i += 2;
      }
      push(body, cur() | (c === '_' ? 32 : 64));
    } else {
      push(c, cur());
      i++;
    }
  }
  return runs;
}

interface WriteCtx {
  BL: number;
  ox: number;
  oy: number;
  labelSize: number;
  nextId: number;
  z: number;
  colors: string[];
}

/** Writes a Chirally document as CDXML. */
export function writeCDXML(doc: ChemDoc): string {
  const BL = doc.style?.bondLengthPt > 0 ? doc.style.bondLengthPt : DEFAULT_BOND_LENGTH_PT;
  const bounds = docBounds(doc) ?? { minX: 0, minY: 0, maxX: 0, maxY: 0 };
  const margin = 54;
  const ctx: WriteCtx = {
    BL,
    ox: margin - bounds.minX * BL,
    oy: margin - bounds.minY * BL,
    labelSize: Math.round((doc.style?.fontSize > 0 ? doc.style.fontSize * BL : 10) * 100) / 100,
    nextId: 1,
    z: 1,
    colors: [],
  };
  const X = (x: number) => f2(x * BL + ctx.ox);
  const Y = (y: number) => f2(y * BL + ctx.oy);
  const P = (x: number, y: number) => `${X(x)} ${Y(y)}`;
  const id = () => ctx.nextId++;
  const Z = () => ctx.z++;
  const colorIdx = (c: string | undefined): number | undefined => {
    if (!c) return undefined;
    const hex = normalizeColor(c);
    if (!hex || hex === '#000000') return undefined;
    let k = ctx.colors.indexOf(hex);
    if (k < 0) k = ctx.colors.push(hex) - 1;
    return k + 4; // 0 black, 1 white (reserved); table: 2 white, 3 black, then custom colours
  };
  const body: string[] = [];
  const steps: string[] = [];

  // ── chemistry ──
  const { mol, index } = docToMol(doc);
  const adj = adjacency(doc);
  const nodeId = new Map<number, number>();
  const fs = ctx.labelSize;
  const label = (x: number, y: number, runs: { text: string; face: number }[], color?: number) => {
    const w = runs.reduce((s, r) => s + r.text.length, 0) * fs * 0.6;
    const lx = x * BL + ctx.ox - fs * 0.33, ly = y * BL + ctx.oy + fs * 0.35;
    return `<t${attrs({ p: `${f2(lx)} ${f2(ly)}`, BoundingBox: `${f2(lx)} ${f2(ly - fs * 0.75)} ${f2(lx + w)} ${f2(ly + fs * 0.1)}`, LabelJustification: 'Left', LabelAlignment: 'Left' })}>` +
      runs.map((r) => `<s${attrs({ font: 3, size: fs, color: color ?? 0, face: r.face })}>${escapeText(r.text)}</s>`).join('') + '</t>';
  };
  const bondsByFrag = new Map<number, DocBond[]>();
  // fragments and their atoms in document order (keeps atom order stable across a round trip)
  const order = new Map<number, number>();
  [...doc.atoms.keys()].forEach((aid, k) => order.set(aid, k));
  const fragList = fragments(doc, adj).map((f) => f.sort((p, q) => order.get(p)! - order.get(q)!));
  fragList.sort((p, q) => order.get(p[0])! - order.get(q[0])!);
  const fragOf = new Map<number, number>();
  fragList.forEach((f, k) => f.forEach((aid) => fragOf.set(aid, k)));
  for (const b of doc.bonds.values()) {
    const k = fragOf.get(b.a);
    if (k === undefined || !doc.atoms.has(b.b)) continue;
    if (!bondsByFrag.has(k)) bondsByFrag.set(k, []);
    bondsByFrag.get(k)!.push(b);
  }
  fragList.forEach((frag, fk) => {
    const fragId = id();
    const xs = frag.map((a) => doc.atoms.get(a)!.x), ys = frag.map((a) => doc.atoms.get(a)!.y);
    const fbb = `${X(Math.min(...xs))} ${Y(Math.min(...ys))} ${X(Math.max(...xs))} ${Y(Math.max(...ys))}`;
    const parts: string[] = [];
    for (const aid of frag) nodeId.set(aid, id());
    for (const aid of frag) {
      const a = doc.atoms.get(aid)!;
      const mi = index.get(aid)!;
      const col = colorIdx(a.color);
      const common = { id: nodeId.get(aid), p: P(a.x, a.y), Z: Z() };
      if (a.abbrev) {
        const nb = neighborsOf(doc, aid, adj).map((n) => doc.atoms.get(n)!)[0];
        const inner = nicknameFragment(a, nb ? { x: nb.x, y: nb.y } : { x: a.x - 1, y: a.y }, P, id, Z, label);
        if (inner) {
          const text = a.abbrev + chargeSuffix(labelDisplayCharge(a.abbrev, a.charge, !!nb));
          parts.push(`<n${attrs({ ...common, NodeType: 'Nickname', NeedsClean: 'yes', color: col })}>${inner}${label(a.x, a.y, [{ text, face: 96 }], col)}</n>`);
          continue;
        }
      }
      if (a.abbrev || a.el === 'R' || a.el === '*' || !atomicNumber(a.el)) {
        const name = a.abbrev ?? a.alias ?? (a.el === '*' ? 'A' : 'R');
        const generic = /^(R\d*'*|X|Ar|Ak|Hal|M|Q|A|G\d*)$/.test(name);
        const net = a.abbrev ? labelDisplayCharge(a.abbrev, a.charge, neighborsOf(doc, aid, adj).length > 0) : a.charge || 0;
        const text = name + chargeSuffix(net); // E+, Nu-, OMe- keep their sign
        parts.push(`<n${attrs({ ...common, NodeType: generic ? 'GenericNickname' : 'Unspecified', GenericNickname: generic ? name : undefined, color: col })}>${label(a.x, a.y, [{ text, face: 0 }], col)}</n>`);
        continue;
      }
      const h = implicitH(mol, mi);
      const showLabel = a.el !== 'C' || a.charge !== 0 || a.isotope !== undefined || (a.radical ?? 0) > 0 || doc.style?.showCarbons === 'all';
      const na = {
        ...common,
        Element: a.el !== 'C' ? atomicNumber(a.el) : undefined,
        NumHydrogens: showLabel ? h : undefined,
        Charge: a.charge || undefined,
        Isotope: a.isotope,
        Radical: a.radical === 1 ? 'Doublet' : a.radical === 2 ? 'Triplet' : undefined,
        color: col,
      };
      if (!showLabel) {
        parts.push(`<n${attrs(na)}/>`);
        continue;
      }
      const runs: { text: string; face: number }[] = [];
      if (a.isotope) runs.push({ text: String(a.isotope), face: 64 });
      const ch = a.charge ? (Math.abs(a.charge) > 1 ? Math.abs(a.charge) : '') + (a.charge > 0 ? '+' : '-') : '';
      runs.push({ text: a.el + (h > 0 ? 'H' + (h > 1 ? h : '') : '') + ch, face: 96 });
      parts.push(`<n${attrs(na)}>${label(a.x, a.y, runs, col)}</n>`);
    }
    for (const b of bondsByFrag.get(fk) ?? []) {
      parts.push(`<b${attrs({ id: id(), Z: Z(), B: nodeId.get(b.a), E: nodeId.get(b.b), ...bondAttrs(b.order, b.style, b.dbPos), color: colorIdx(b.color) })}/>`);
    }
    body.push(`<fragment${attrs({ id: fragId, BoundingBox: fbb, Z: Z() })}>${parts.join('')}</fragment>`);
  });

  // ── text ──
  const textXml = (x: number, y: number, runs: { text: string; face: number }[], size: number, align: string | undefined, color: string | undefined, tid = id()) => {
    const lines = runs.map((r) => r.text).join('').split('\n');
    const w = Math.max(...lines.map((l) => l.length)) * size * 0.55;
    const left = align === 'center' ? x - w / 2 : align === 'right' ? x - w : x;
    const bbox = `${f2(left)} ${f2(y - size * 0.8)} ${f2(left + w)} ${f2(y + size * 0.25 + (lines.length - 1) * size * 1.2)}`;
    const col = colorIdx(color) ?? 0;
    return {
      id: tid,
      xml: `<t${attrs({ id: tid, p: `${f2(x)} ${f2(y)}`, BoundingBox: bbox, Z: Z(), Justification: align === 'center' ? 'Center' : align === 'right' ? 'Right' : undefined, InterpretChemically: 'no' })}>` +
        runs.map((r) => `<s${attrs({ font: 3, size: f2(size), color: col, face: r.face })}>${escapeText(r.text)}</s>`).join('') + '</t>',
    };
  };
  for (const t of doc.texts.values()) {
    const base = (t.bold ? 1 : 0) | (t.italic ? 2 : 0);
    const runs = t.formula ? [{ text: t.text, face: 96 | base }] : markupToRuns(t.text, base);
    if (!runs.length) continue;
    body.push(textXml(t.x * BL + ctx.ox, t.y * BL + ctx.oy, runs, ctx.labelSize * (t.size ?? 1), t.align, t.color).xml);
  }

  // ── arrows ──
  const arrowCommon = { FillType: 'None', ArrowheadType: 'Solid', HeadSize: 1000, ArrowheadCenterSize: 875, ArrowheadWidth: 250 };
  for (const a of doc.arrows.values()) {
    const aid = id();
    const k = a.kind;
    const extra: Record<string, string | number | undefined> = { ArrowheadHead: 'Full' };
    if (k === 'equilibrium' || k === 'unbalancedEq') Object.assign(extra, { ArrowheadHead: 'HalfLeft', ArrowheadTail: 'HalfLeft', ArrowShaftSpacing: 400, ArrowEquilibriumRatio: k === 'unbalancedEq' ? 200 : undefined });
    else if (k === 'retro') Object.assign(extra, { ArrowheadType: 'Hollow', ArrowShaftSpacing: 800 });
    else if (k === 'resonance') extra.ArrowheadTail = 'Full';
    else if (k === 'dashed') extra.LineType = 'Dashed';
    else if (k === 'noGo') extra.NoGo = 'Cross';
    else if (k === 'line') extra.ArrowheadHead = undefined;
    const bb = `${X(Math.min(a.x1, a.x2))} ${Y(Math.min(a.y1, a.y2) - 0.1)} ${X(Math.max(a.x1, a.x2))} ${Y(Math.max(a.y1, a.y2) + 0.1)}`;
    body.push(`<arrow${attrs({ id: aid, BoundingBox: bb, Z: Z(), ...arrowCommon, ...extra, Head3D: `${P(a.x2, a.y2)} 0`, Tail3D: `${P(a.x1, a.y1)} 0`, color: colorIdx(a.color) })}/>`);
    // reaction conditions as captions + a <step> linking them to the arrow
    const mx = ((a.x1 + a.x2) / 2) * BL + ctx.ox;
    const my = ((a.y1 + a.y2) / 2) * BL + ctx.oy;
    const above: number[] = [], below: number[] = [];
    if (a.above) {
      const nLines = a.above.split('\n').length;
      const t = textXml(mx, my - fs * 0.6 - (nLines - 1) * fs * 1.2, markupToRuns(a.above, 0), fs, 'center', a.color);
      body.push(t.xml);
      above.push(t.id);
    }
    if (a.below) {
      const t = textXml(mx, my + fs * 1.3, markupToRuns(a.below, 0), fs, 'center', a.color);
      body.push(t.xml);
      below.push(t.id);
    }
    if (above.length || below.length) {
      steps.push(`<step${attrs({ id: id(), ReactionStepArrows: aid, ReactionStepObjectsAboveArrow: above.join(' ') || undefined, ReactionStepObjectsBelowArrow: below.join(' ') || undefined })}/>`);
    }
  }

  // ── curved arrows (written as ChemDraw arc arrows) ──
  for (const c of doc.curved.values()) {
    const S = anchorPoint(doc, c.from), E = anchorPoint(doc, c.to);
    if (!S || !E) continue;
    const len = Math.hypot(E.x - S.x, E.y - S.y);
    if (len < 1e-6) continue;
    const ux = (E.x - S.x) / len, uy = (E.y - S.y) / len;
    const cp = (q: { t: number; h: number }) => ({ x: S.x + ux * q.t * len - uy * q.h * len, y: S.y + uy * q.t * len + ux * q.h * len });
    const c1 = cp(c.c1), c2 = cp(c.c2);
    const M = { x: (S.x + 3 * c1.x + 3 * c2.x + E.x) / 8, y: (S.y + 3 * c1.y + 3 * c2.y + E.y) / 8 };
    const arc = circleThrough(S, M, E);
    const head = c.electrons === 1 ? 'HalfLeft' : 'Full';
    const pts = [S, c1, c2, E];
    const bb = `${X(Math.min(...pts.map((p) => p.x)))} ${Y(Math.min(...pts.map((p) => p.y)))} ${X(Math.max(...pts.map((p) => p.x)))} ${Y(Math.max(...pts.map((p) => p.y)))}`;
    const common = { id: id(), BoundingBox: bb, Z: Z(), ...arrowCommon, ArrowheadHead: head, Head3D: `${P(E.x, E.y)} 0`, Tail3D: `${P(S.x, S.y)} 0`, color: colorIdx(c.color) };
    if (!arc) {
      body.push(`<arrow${attrs(common)}/>`);
      continue;
    }
    // Major axis points at the tail; minor axis is the direction of travel (sweep tail → head)
    const vT = { x: S.x - arc.c.x, y: S.y - arc.c.y };
    const perp = { x: -vT.y * arc.dir, y: vT.x * arc.dir };
    body.push(`<arrow${attrs({
      ...common,
      AngularSize: f2(arc.sweep * 180 / Math.PI),
      Center3D: `${P(arc.c.x, arc.c.y)} 0`,
      MajorAxisEnd3D: `${P(S.x, S.y)} 0`,
      MinorAxisEnd3D: `${P(arc.c.x + perp.x, arc.c.y + perp.y)} 0`,
    })}/>`);
  }

  // ── shapes ──
  for (const s of doc.shapes.values()) {
    const x1 = Math.min(s.x1, s.x2), y1 = Math.min(s.y1, s.y2), x2 = Math.max(s.x1, s.x2), y2 = Math.max(s.y1, s.y2);
    const box = `${X(x1)} ${Y(y1)} ${X(x2)} ${Y(y2)}`;
    const col = colorIdx(s.color);
    const g = (a: Record<string, string | number | undefined>) => body.push(`<graphic${attrs({ id: id(), BoundingBox: box, Z: Z(), ...a, color: col })}/>`);
    const cx = (x1 + x2) / 2, cy = (y1 + y2) / 2, rx = (x2 - x1) / 2, ry = (y2 - y1) / 2;
    const ang = s.angle ?? 0;
    const axes = {
      Center3D: `${P(cx, cy)} 0`,
      MajorAxisEnd3D: `${P(cx + rx * Math.cos(ang), cy + rx * Math.sin(ang))} 0`,
      MinorAxisEnd3D: `${P(cx - ry * Math.sin(ang), cy + ry * Math.cos(ang))} 0`,
    };
    switch (s.kind) {
      case 'rect': g({ GraphicType: 'Rectangle', RectangleType: s.dashed ? 'Dashed' : undefined }); break;
      case 'roundRect': g({ GraphicType: 'Rectangle', RectangleType: s.dashed ? 'RoundEdge Dashed' : 'RoundEdge' }); break;
      case 'ellipse': g({ GraphicType: 'Oval', OvalType: s.dashed ? 'Dashed' : undefined, ...axes }); break;
      case 'line': body.push(`<graphic${attrs({ id: id(), BoundingBox: `${P(s.x1, s.y1)} ${P(s.x2, s.y2)}`, Z: Z(), GraphicType: 'Line', LineType: s.dashed ? 'Dashed' : undefined, color: col })}/>`); break;
      case 'bracket': case 'tsBracket': g({ GraphicType: 'Bracket', BracketType: 'SquarePair' }); break;
      case 'paren': g({ GraphicType: 'Bracket', BracketType: 'RoundPair' }); break;
      case 'brace': g({ GraphicType: 'Bracket', BracketType: 'CurlyPair' }); break;
      case 'orbitalS': g({ GraphicType: 'Orbital', OrbitalType: 's', ...axes }); break;
      case 'orbitalP': g({ GraphicType: 'Orbital', OrbitalType: 'p', ...axes }); break;
    }
    const lab = s.label ?? (s.kind === 'tsBracket' ? '‡' : undefined);
    if (lab) body.push(textXml(x2 * BL + ctx.ox + 2, y1 * BL + ctx.oy + fs * 0.7, markupToRuns(lab, 0), fs, undefined, s.color).xml);
  }
  if (steps.length) body.push(`<scheme${attrs({ id: id() })}>${steps.join('')}</scheme>`);

  // ── document ──
  const W = Math.max(540, (bounds.maxX - bounds.minX) * BL + 2 * margin);
  const H = Math.max(720, (bounds.maxY - bounds.minY) * BL + 2 * margin);
  const st = doc.style;
  const rootAttrs = attrs({
    CreationProgram: 'Chirally',
    Name: doc.meta?.title || undefined,
    BoundingBox: `0 0 ${f2(W)} ${f2(H)}`,
    WindowPosition: '0 0',
    WindowSize: '0 0',
    FractionalWidths: 'yes',
    InterpretChemically: 'yes',
    ShowAtomQuery: 'yes',
    ShowAtomStereo: 'no',
    ShowAtomEnhancedStereo: 'yes',
    ShowAtomNumber: 'no',
    ShowBondQuery: 'yes',
    ShowBondRxn: 'yes',
    ShowBondStereo: 'no',
    ShowTerminalCarbonLabels: st?.showCarbons === 'terminal' || st?.showCarbons === 'all' ? 'yes' : 'no',
    ShowNonTerminalCarbonLabels: st?.showCarbons === 'all' ? 'yes' : 'no',
    HideImplicitHydrogens: st && st.showImplicitH === false ? 'yes' : 'no',
    LabelFont: 3,
    LabelSize: f2(ctx.labelSize),
    LabelFace: 96,
    CaptionFont: 3,
    CaptionSize: f2(ctx.labelSize),
    HashSpacing: f2((st?.hashSpacing ?? 2.5 / 14.4) * BL),
    MarginWidth: f2((st?.labelMargin ?? 1.6 / 14.4) * BL),
    LineWidth: f2((st?.lineWidth ?? 0.6 / 14.4) * BL),
    BoldWidth: f2((st?.boldWidth ?? 2 / 14.4) * BL),
    BondLength: f2(BL),
    BondSpacing: f2((st?.bondSpacing ?? 0.18) * 100),
    ChainAngle: 120,
    LabelJustification: 'Auto',
    CaptionJustification: 'Left',
    color: 0,
    bgcolor: 1,
  });
  const rgb = (hex: string) => {
    const v = parseInt(hex.slice(1), 16);
    return attrs({ r: f2(((v >> 16) & 255) / 255), g: f2(((v >> 8) & 255) / 255), b: f2((v & 255) / 255) });
  };
  const colortable = ['<color r="1" g="1" b="1"/>', '<color r="0" g="0" b="0"/>', ...ctx.colors.map((c) => `<color${rgb(c)}/>`)].join('\n');
  return [
    '<?xml version="1.0" encoding="UTF-8" ?>',
    '<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >',
    `<CDXML${rootAttrs}>`,
    `<colortable>\n${colortable}\n</colortable>`,
    '<fonttable>\n<font id="3" charset="iso-8859-1" name="Arial"/>\n</fonttable>',
    `<page${attrs({ id: id(), BoundingBox: `0 0 ${f2(W)} ${f2(H)}`, HeaderPosition: 36, FooterPosition: 36, PrintTrimMarks: 'yes', HeightPages: Math.ceil(H / 720), WidthPages: Math.ceil(W / 540) })}>`,
    ...body,
    '</page>',
    '</CDXML>',
    '',
  ].join('\n');
}

function bondAttrs(order: number, style: BondStyle, dbPos: string | undefined): Record<string, string | undefined> {
  const out: Record<string, string | undefined> = {};
  if (style === 'dative') out.Order = 'dative';
  else if (style === 'hbond' || order === 0) {
    out.Order = 'hydrogen';
    out.Display = 'Dash';
  } else if (order !== 1) out.Order = order === 1.5 ? '1.5' : String(order);
  const disp: Partial<Record<BondStyle, string>> = {
    wedge: 'WedgeBegin', hash: 'WedgedHashBegin', hollow: 'HollowWedgeBegin', wavy: 'Wavy', bold: 'Bold', dashed: 'Dash',
  };
  if (disp[style]) out.Display = disp[style];
  if (dbPos === 'left' || dbPos === 'right' || dbPos === 'center') out.DoublePosition = dbPos[0].toUpperCase() + dbPos.slice(1);
  return out;
}

/** Inner fragment of a nickname node: the expanded group plus an ExternalConnectionPoint. */
function nicknameFragment(
  a: Atom, nb: Pt,
  P: (x: number, y: number) => string, id: () => number, Z: () => number,
  label: (x: number, y: number, runs: { text: string; face: number }[]) => string,
): string | null {
  if (!a.abbrev || !abbreviationMol(a.abbrev)) return null;
  const tmp = new Mol();
  const n0 = tmp.addAtom({ el: 'C', x: nb.x, y: nb.y });
  const n1 = tmp.addAtom({ el: a.el, abbrev: a.abbrev, charge: a.charge, x: a.x, y: a.y });
  tmp.addBond(n0, n1);
  const ex = expandAbbreviations(tmp);
  const ids = ex.atoms.map(() => id());
  const parts: string[] = [];
  ex.atoms.forEach((g, i) => {
    if (i === 0) {
      parts.push(`<n${attrs({ id: ids[0], p: P(g.x, g.y), Z: Z(), NodeType: 'ExternalConnectionPoint' })}/>`);
      return;
    }
    const h = implicitH(ex, i);
    const lab = g.el !== 'C' || g.charge !== 0;
    const na = { id: ids[i], p: P(g.x, g.y), Z: Z(), Element: g.el !== 'C' ? atomicNumber(g.el) : undefined, NumHydrogens: lab ? h : undefined, Charge: g.charge || undefined };
    if (!lab) parts.push(`<n${attrs(na)}/>`);
    else {
      const ch = g.charge ? (Math.abs(g.charge) > 1 ? Math.abs(g.charge) : '') + (g.charge > 0 ? '+' : '-') : '';
      parts.push(`<n${attrs(na)}>${label(g.x, g.y, [{ text: g.el + (h > 0 ? 'H' + (h > 1 ? h : '') : '') + ch, face: 96 }])}</n>`);
    }
  });
  for (const b of ex.bonds) parts.push(`<b${attrs({ id: id(), Z: Z(), B: ids[b.a], E: ids[b.b], ...bondAttrs(b.order, b.style, undefined) })}/>`);
  return `<fragment${attrs({ id: id() })}>${parts.join('')}</fragment>`;
}

/** Model-space position of a curved-arrow anchor. */
function anchorPoint(doc: ChemDoc, an: Anchor): Pt | null {
  switch (an.type) {
    case 'point': return { x: an.x, y: an.y };
    case 'atom': {
      if (an.h) return implicitHydrogenPoint(doc, an.id); // on the atom's implicit H
      const a = doc.atoms.get(an.id);
      return a ? { x: a.x, y: a.y } : null;
    }
    case 'bond': {
      const b = doc.bonds.get(an.id);
      const p = b && doc.atoms.get(b.a), q = b && doc.atoms.get(b.b);
      return p && q ? { x: (p.x + q.x) / 2, y: (p.y + q.y) / 2 } : null;
    }
    case 'between': {
      const p = doc.atoms.get(an.a), q = doc.atoms.get(an.b);
      return p && q ? { x: (p.x + q.x) / 2, y: (p.y + q.y) / 2 } : null;
    }
  }
  return null;
}

/**
 * Circle through start S, mid-curve M and end E: centre, sweep angle (radians) from S to E via M,
 * and rotation direction (+1 = from x toward y axis). Null when the points are collinear.
 */
function circleThrough(S: Pt, M: Pt, E: Pt): { c: Pt; sweep: number; dir: 1 | -1 } | null {
  const d = 2 * (S.x * (M.y - E.y) + M.x * (E.y - S.y) + E.x * (S.y - M.y));
  const span = Math.hypot(E.x - S.x, E.y - S.y);
  if (Math.abs(d) < 1e-9 * Math.max(1, span * span) || Math.abs(d) / (span * span) < 1e-3) return null;
  const s2 = S.x * S.x + S.y * S.y, m2 = M.x * M.x + M.y * M.y, e2 = E.x * E.x + E.y * E.y;
  const c = {
    x: (s2 * (M.y - E.y) + m2 * (E.y - S.y) + e2 * (S.y - M.y)) / d,
    y: (s2 * (E.x - M.x) + m2 * (S.x - E.x) + e2 * (M.x - S.x)) / d,
  };
  const ang = (p: Pt) => Math.atan2(p.y - c.y, p.x - c.x);
  const norm = (a: number) => ((a % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
  const aS = ang(S);
  const toM = norm(ang(M) - aS), toE = norm(ang(E) - aS);
  // positive direction reaches M before E → sweep toE; otherwise go the other way round
  if (toM <= toE) return { c, sweep: toE, dir: 1 };
  return { c, sweep: 2 * Math.PI - toE, dir: -1 };
}

function normalizeColor(c: string): string | null {
  const s = c.trim().toLowerCase();
  if (/^#[0-9a-f]{6}$/.test(s)) return s;
  if (/^#[0-9a-f]{3}$/.test(s)) return '#' + s[1] + s[1] + s[2] + s[2] + s[3] + s[3];
  const m = /^rgba?\((\d+)[,\s]+(\d+)[,\s]+(\d+)/.exec(s);
  if (m) return '#' + [m[1], m[2], m[3]].map((v) => Math.min(255, +v).toString(16).padStart(2, '0')).join('');
  const named: Record<string, string> = { black: '#000000', white: '#ffffff', red: '#ff0000', blue: '#0000ff', green: '#008000' };
  return named[s] ?? null;
}
