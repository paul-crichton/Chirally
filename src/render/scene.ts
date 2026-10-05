// Converts a ChemDoc into a list of drawing primitives (model units). Used by the canvas view,
// SVG/PNG export and thumbnails, so every visual rule lives here.
import { ChemDoc, DocStyle, ArrowObj, CurvedArrowObj, TextObj, ShapeObj } from '../doc/types';
import { docToMol } from '../doc/document';
import { Mol } from '../chem/mol';
import { implicitH, hasValenceError, lonePairCount } from '../chem/valence';
import { perceiveRings, RingInfo } from '../chem/rings';
import { perceiveAromaticity } from '../chem/aromaticity';
import { reverseLabel } from '../chem/abbreviations';
import { labelColor, element } from '../chem/elements';
import { Box, Pt, clipSegment, norm, perp, sub, add, mul, len, dot, lerp, angleOf, fromAngle, cross, bezierTangent } from './geom';
import { measureText, richLines, TextRun, SCRIPT_SCALE, SUB_SHIFT, SUP_SHIFT, formulaRuns } from './text';

export type PathSeg =
  | ['M', number, number]
  | ['L', number, number]
  | ['Q', number, number, number, number]
  | ['C', number, number, number, number, number, number]
  | ['Z'];

export type Prim =
  | { k: 'line'; x1: number; y1: number; x2: number; y2: number; w: number; color: string; dash?: number[]; cap?: 'butt' | 'round' | 'square' }
  | { k: 'poly'; pts: number[]; fill?: string; stroke?: string; w?: number; closed?: boolean; dash?: number[] }
  | { k: 'path'; d: PathSeg[]; fill?: string; stroke?: string; w?: number; dash?: number[]; cap?: 'butt' | 'round' }
  | { k: 'circle'; x: number; y: number; r: number; fill?: string; stroke?: string; w?: number; dash?: number[] }
  | { k: 'text'; x: number; y: number; text: string; size: number; family: string; color: string; bold?: boolean; italic?: boolean };

export interface SceneOptions {
  /** Default drawing colour (bonds, carbon, text). */
  ink: string;
  /** Dark background → brighter heteroatom colours. */
  dark?: boolean;
  /** atom id → CIP label shown when style.showStereoLabels */
  stereoLabels?: Map<number, string>;
  /** bond id → E/Z label */
  bondLabels?: Map<number, string>;
  /** atom id → small annotation (atom numbers / IUPAC locants) */
  atomNotes?: Map<number, string>;
  /** atoms whose label is hidden (being edited) */
  hiddenLabels?: Set<number>;
  /** Objects not to draw (being edited in place) */
  hiddenObjects?: Set<number>;
  /** atom id → highlight colour (e.g. electron-count warnings) */
  atomHalo?: Map<number, string>;
  /** Draw valence error boxes */
  showErrors?: boolean;
}

export interface Scene {
  prims: Prim[];
  /** atom id → label boxes (for hit testing / clipping); absent when the atom has no label */
  labelBoxes: Map<number, Box[]>;
  /** object id → bounding box (texts, arrow captions) */
  objectBoxes: Map<number, Box>;
  mol: Mol;
  index: Map<number, number>;
  rings: RingInfo;
  /** atom id → total implicit/explicit H count */
  hCount: Map<number, number>;
  bounds: Box | null;
}

const CAP = 0.72; // cap height as fraction of font size

const DARK_COLORS: Record<string, string> = {
  N: '#7b93ff', O: '#ff5a5a', S: '#e8c840', P: '#ff9a3c', F: '#7bdc7b', Cl: '#5fd35f', Br: '#d9786e',
  I: '#c77dff', B: '#ffb0b0', Si: '#e0c090', Na: '#c79bff', K: '#b88bff', H: '#d0d0d0',
};

export function atomInkColor(el: string, ink: string, style: DocStyle, dark?: boolean): string {
  if (!style.colorAtoms || el === 'C' || !element(el)) return ink;
  if (dark) return DARK_COLORS[el] ?? (el === 'C' ? ink : labelColor(el));
  const c = labelColor(el);
  return c === '#000000' ? ink : c;
}

interface PosRun {
  x: number;
  y: number;
  text: string;
  size: number;
  bold?: boolean;
  italic?: boolean;
}

export function buildScene(doc: ChemDoc, opts: SceneOptions): Scene {
  const st = doc.style;
  const prims: Prim[] = [];
  const { mol, index } = docToMol(doc);
  const rings = perceiveRings(mol);
  const n = mol.atoms.length;
  const fam = st.fontFamily;
  const fs = st.fontSize;
  const hc: number[] = new Array(n);
  const hCount = new Map<number, number>();
  for (let i = 0; i < n; i++) {
    hc[i] = implicitH(mol, i);
    hCount.set(mol.atoms[i].id, hc[i]);
  }
  const pos: Pt[] = mol.atoms.map((a) => ({ x: a.x, y: a.y }));
  const nbrDirs: Pt[][] = mol.atoms.map((_, i) => mol.neighbors(i).map((j) => norm(sub(pos[j], pos[i]))));

  // ───── atom labels ─────
  const labelBoxes = new Map<number, Box[]>();
  const boxesByIdx: (Box[] | null)[] = new Array(n).fill(null);
  const labelRuns: { runs: PosRun[]; color: string; i: number }[] = [];
  const hSideOf: (string | null)[] = new Array(n).fill(null);

  for (let i = 0; i < n; i++) {
    const a = mol.atoms[i];
    if (opts.hiddenLabels?.has(a.id)) continue;
    const deg = mol.degree(i);
    if (!labelVisible(mol, i, deg, st, hc[i])) continue;
    const color = a.color ?? (a.abbrev || a.alias ? opts.ink : atomInkColor(a.el, opts.ink, st, opts.dark));
    const lab = layoutLabel(mol, i, nbrDirs[i], hc[i], st);
    hSideOf[i] = lab.hSide;
    boxesByIdx[i] = lab.boxes;
    labelBoxes.set(a.id, lab.boxes);
    labelRuns.push({ runs: lab.runs, color, i });
  }

  // ───── bonds ─────
  let aro: boolean[] | null = null;
  let aroBonds: boolean[] | null = null;
  if (st.aromaticCircles) {
    const ar = perceiveAromaticity(mol, rings);
    aro = ar.atoms;
    aroBonds = ar.bonds;
  }
  const w = st.lineWidth;
  const sp = st.bondSpacing;
  mol.bonds.forEach((b, bi) => {
    const color = b.color ?? opts.ink;
    const A = pos[b.a], B = pos[b.b];
    const bxA = boxesByIdx[b.a], bxB = boxesByIdx[b.b];
    const order = aroBonds && aroBonds[bi] && (b.order === 1 || b.order === 2) ? 1 : b.order;
    const style = b.style;
    const u = norm(sub(B, A));
    const L = len(sub(B, A));
    if (L < 1e-6) return;
    const pv = perp(u);
    const line = (p: Pt, q: Pt, extra: Partial<Extract<Prim, { k: 'line' }>> = {}) =>
      prims.push({ k: 'line', x1: p.x, y1: p.y, x2: q.x, y2: q.y, w, color, cap: 'round', ...extra });
    const clipped = clipSegment(A, B, bxA, bxB);

    if (style === 'wedge' || style === 'hollow') {
      if (!clipped) return;
      const [p, q] = clipped;
      const ww = st.wedgeWidth / 2;
      const tq = len(sub(q, A)) / L; // width proportional to distance from narrow end
      const wq = Math.max(w / 2, ww * tq);
      const wp = w / 2;
      const pts = [p.x + pv.x * wp, p.y + pv.y * wp, q.x + pv.x * wq, q.y + pv.y * wq, q.x - pv.x * wq, q.y - pv.y * wq, p.x - pv.x * wp, p.y - pv.y * wp];
      if (style === 'wedge') prims.push({ k: 'poly', pts, fill: color, stroke: color, w: w * 0.4, closed: true });
      else prims.push({ k: 'poly', pts, stroke: color, w, closed: true });
      return;
    }
    if (style === 'hash') {
      if (!clipped) return;
      const [p, q] = clipped;
      const ww = st.wedgeWidth / 2;
      const segLen = len(sub(q, p));
      const count = Math.max(3, Math.round(segLen / st.hashSpacing));
      for (let k = 1; k <= count; k++) {
        const t = k / count;
        const c = lerp(p, q, t);
        const dA = len(sub(c, A)) / L;
        const hw = Math.max(w * 0.6, ww * dA);
        line(add(c, mul(pv, hw)), sub(c, mul(pv, hw)), { w: w * 0.9, cap: 'butt' });
      }
      return;
    }
    if (style === 'wavy') {
      if (!clipped) return;
      const [p, q] = clipped;
      const segLen = len(sub(q, p));
      const waves = Math.max(2, Math.round(segLen / 0.16));
      const amp = 0.07;
      const d: PathSeg[] = [['M', p.x, p.y]];
      for (let k = 0; k < waves; k++) {
        const t0 = k / waves, t1 = (k + 1) / waves;
        const mid = lerp(p, q, (t0 + t1) / 2);
        const side = k % 2 === 0 ? 1 : -1;
        const c = add(mid, mul(pv, amp * 2 * side));
        const e = lerp(p, q, t1);
        d.push(['Q', c.x, c.y, e.x, e.y]);
      }
      prims.push({ k: 'path', d, stroke: color, w, cap: 'round' });
      return;
    }
    if (style === 'hbond' || order === 0) {
      if (!clipped) return;
      line(clipped[0], clipped[1], { dash: [w * 0.6, w * 2.6] });
      return;
    }
    if (style === 'dative') {
      if (!clipped) return;
      const [p, q] = clipped;
      const headLen = 0.22;
      const tip = q;
      const base = sub(tip, mul(u, headLen));
      line(p, base);
      prims.push(arrowHeadPoly(tip, u, headLen, 0.09, color, 'full'));
      return;
    }
    if (order === 1 && style !== 'crossed') {
      if (!clipped) return;
      if (style === 'bold') line(clipped[0], clipped[1], { w: st.boldWidth, cap: 'butt' });
      else if (style === 'dashed') line(clipped[0], clipped[1], { dash: [0.09, 0.07], cap: 'butt' });
      else line(clipped[0], clipped[1]);
      return;
    }

    // multiple bonds
    if (order === 3 || order === 4) {
      if (clipped) line(clipped[0], clipped[1]);
      const offs = order === 3 ? [sp, -sp] : [sp, -sp, 2 * sp];
      for (const o of offs) {
        const pa = add(A, mul(pv, o)), pb = add(B, mul(pv, o));
        const c = clipSegment(pa, pb, bxA, bxB);
        if (!c) continue;
        // shorten slightly at unlabeled ends
        const s0 = bxA ? c[0] : add(c[0], mul(u, 0.04));
        const s1 = bxB ? c[1] : sub(c[1], mul(u, 0.04));
        line(s0, s1);
      }
      return;
    }

    if (style === 'crossed') {
      const h = sp / 2;
      const c1 = clipSegment(add(A, mul(pv, h)), sub(B, mul(pv, h)), bxA, bxB);
      const c2 = clipSegment(sub(A, mul(pv, h)), add(B, mul(pv, h)), bxA, bxB);
      if (c1) line(c1[0], c1[1]);
      if (c2) line(c2[0], c2[1]);
      return;
    }

    // double (2) or delocalised (1.5)
    const placement = doubleBondPlacement(mol, bi, rings, pos, b.dbPos ?? 'auto', boxesByIdx);
    if (placement.mode === 'center' && order === 2) {
      const h = sp / 2;
      const c1 = clipSegment(add(A, mul(pv, h)), add(B, mul(pv, h)), bxA, bxB);
      const c2 = clipSegment(sub(A, mul(pv, h)), sub(B, mul(pv, h)), bxA, bxB);
      if (c1) line(c1[0], c1[1], style === 'bold' ? { w: st.boldWidth } : {});
      if (c2) line(c2[0], c2[1]);
      return;
    }
    const side = placement.side;
    if (clipped) {
      if (style === 'bold') line(clipped[0], clipped[1], { w: st.boldWidth, cap: 'butt' });
      else if (style === 'dashed') line(clipped[0], clipped[1], { dash: [0.09, 0.07], cap: 'butt' });
      else line(clipped[0], clipped[1]);
    }
    const off = mul(pv, sp * side);
    let pa = add(A, off), pb = add(B, off);
    const trimA = bxA ? 0 : innerTrim(nbrDirs[b.a], u, side, sp, mol.atoms[b.a].id, pv);
    const trimB = bxB ? 0 : innerTrim(nbrDirs[b.b], mul(u, -1), -side, sp, mol.atoms[b.b].id, mul(pv, -1));
    pa = add(pa, mul(u, trimA));
    pb = sub(pb, mul(u, trimB));
    const c = clipSegment(pa, pb, bxA, bxB);
    if (c) line(c[0], c[1], order === 1.5 ? { dash: [0.08, 0.06], cap: 'butt' } : {});
  });

  // aromatic circles
  if (aro && st.aromaticCircles) {
    const ar = perceiveAromaticity(mol, rings);
    rings.rings.forEach((r, k) => {
      if (!r.every((i) => ar.atoms[i]) || !rings.ringBonds[k].every((bi) => ar.bonds[bi])) return;
      const c = r.reduce((s, i) => add(s, pos[i]), { x: 0, y: 0 });
      const cc = mul(c, 1 / r.length);
      const apo = r.reduce((s, i) => s + len(sub(pos[i], cc)), 0) / r.length * Math.cos(Math.PI / r.length);
      prims.push({ k: 'circle', x: cc.x, y: cc.y, r: apo * 0.62, stroke: opts.ink, w });
    });
  }

  // halos (behind labels)
  if (opts.atomHalo) {
    for (const [id, col] of opts.atomHalo) {
      const i = index.get(id);
      if (i === undefined) continue;
      prims.unshift({ k: 'circle', x: pos[i].x, y: pos[i].y, r: 0.32, fill: col });
    }
  }

  // ───── labels ─────
  for (const L of labelRuns) {
    for (const r of L.runs) prims.push({ k: 'text', x: r.x, y: r.y, text: r.text, size: r.size, family: fam, color: L.color, bold: r.bold, italic: r.italic });
  }

  // ───── atom decorations: charges / radicals on unlabeled atoms, lone pairs, stereo labels, notes, errors ─────
  for (let i = 0; i < n; i++) {
    const a = mol.atoms[i];
    const labeled = !!boxesByIdx[i];
    const occupied = nbrDirs[i].map(angleOf);
    const hs = hSideOf[i];
    if (hs === 'right') occupied.push(0);
    if (hs === 'left') occupied.push(Math.PI);
    if (hs === 'down') occupied.push(Math.PI / 2);
    if (hs === 'up') occupied.push(-Math.PI / 2);
    const color = a.color ?? atomInkColor(a.el, opts.ink, st, opts.dark);
    const rad = labeled ? labelRadius(boxesByIdx[i]!, pos[i]) : 0;
    // superscript charge / isotope occupy the upper corners of a label
    if (labeled && a.charge) occupied.push(-Math.PI / 4);
    if (labeled && a.isotope) occupied.push((-3 * Math.PI) / 4);
    if (!labeled && a.charge) {
      const ang = freeAngles(occupied, 1, -Math.PI / 4)[0];
      occupied.push(ang);
      const p = add(pos[i], fromAngle(ang, 0.36));
      const txt = chargeText(a.charge);
      const size = fs * 0.8;
      const tw = measureText(txt, size, fam);
      prims.push({ k: 'text', x: p.x - tw / 2, y: p.y + size * CAP / 2, text: txt, size, family: fam, color: opts.ink });
    }
    const radicals = a.radical ?? 0;
    if (radicals) {
      const ang = labeled ? -Math.PI / 2 : freeAngles(occupied, 1, -Math.PI / 2)[0];
      if (!labeled) occupied.push(ang);
      const r0 = labeled ? rad + 0.1 : 0.26;
      const c = add(pos[i], fromAngle(ang, r0));
      const t = perp(fromAngle(ang));
      if (radicals === 1) prims.push({ k: 'circle', x: c.x, y: c.y, r: 0.045, fill: opts.ink });
      else {
        prims.push({ k: 'circle', x: c.x + t.x * 0.06, y: c.y + t.y * 0.06, r: 0.045, fill: opts.ink });
        prims.push({ k: 'circle', x: c.x - t.x * 0.06, y: c.y - t.y * 0.06, r: 0.045, fill: opts.ink });
      }
    }
    const showLP = a.lonePairs || (st.showLonePairs && a.el !== 'C' && a.el !== 'H');
    if (showLP && !a.abbrev) {
      const nlp = lonePairCount(mol, i);
      if (nlp > 0) {
        const angs = freeAngles(occupied, nlp, -Math.PI / 2);
        for (const ang of angs) {
          occupied.push(ang);
          const r0 = labeled ? rad + 0.09 : 0.24;
          const c = add(pos[i], fromAngle(ang, r0));
          const t = perp(fromAngle(ang));
          for (const s of [1, -1]) prims.push({ k: 'circle', x: c.x + t.x * 0.065 * s, y: c.y + t.y * 0.065 * s, r: 0.038, fill: color });
        }
      }
    }
    const stereo = st.showStereoLabels ? opts.stereoLabels?.get(a.id) : undefined;
    if (stereo) {
      const ang = freeAngles(occupied, 1, Math.PI / 2)[0];
      occupied.push(ang);
      const p = add(pos[i], fromAngle(ang, (labeled ? rad : 0.1) + 0.32));
      const txt = `(${stereo})`;
      const size = fs * 0.62;
      const tw = measureText(txt, size, fam, false, true);
      prims.push({ k: 'text', x: p.x - tw / 2, y: p.y + size * CAP / 2, text: txt, size, family: fam, color: '#9b3fb5', italic: true });
    }
    const note = opts.atomNotes?.get(a.id);
    if (note) {
      const ang = freeAngles(occupied, 1, (3 * Math.PI) / 4)[0];
      const p = add(pos[i], fromAngle(ang, (labeled ? rad : 0.05) + 0.27));
      const size = fs * 0.55;
      const tw = measureText(note, size, fam);
      prims.push({ k: 'text', x: p.x - tw / 2, y: p.y + size * CAP / 2, text: note, size, family: fam, color: '#1c7ed6' });
    }
    if (opts.showErrors !== false && st.showValenceErrors && hasValenceError(mol, i)) {
      if (labeled) {
        const bx = unionBox(boxesByIdx[i]!);
        prims.push({ k: 'poly', pts: [bx.x1, bx.y1, bx.x2, bx.y1, bx.x2, bx.y2, bx.x1, bx.y2], stroke: '#e03131', w: w * 0.8, closed: true });
      } else prims.push({ k: 'circle', x: pos[i].x, y: pos[i].y, r: 0.2, stroke: '#e03131', w: w * 0.8 });
    }
  }

  // E/Z labels on bonds
  if (st.showStereoLabels && opts.bondLabels) {
    for (const [bid, lab] of opts.bondLabels) {
      const b = doc.bonds.get(bid);
      if (!b) continue;
      const A = doc.atoms.get(b.a)!, B = doc.atoms.get(b.b)!;
      const mid = lerp(A, B, 0.5);
      const pv = perp(norm(sub(B, A)));
      const p = add(mid, mul(pv, -0.35));
      const txt = `(${lab})`;
      const size = fs * 0.62;
      const tw = measureText(txt, size, fam, false, true);
      prims.push({ k: 'text', x: p.x - tw / 2, y: p.y + size * CAP / 2, text: txt, size, family: fam, color: '#9b3fb5', italic: true });
    }
  }

  // ───── objects ─────
  const objectBoxes = new Map<number, Box>();
  for (const s of doc.shapes.values()) {
    if (opts.hiddenObjects?.has(s.id)) continue;
    shapePrims(s, st, opts.ink, prims);
    objectBoxes.set(s.id, { x1: Math.min(s.x1, s.x2), y1: Math.min(s.y1, s.y2), x2: Math.max(s.x1, s.x2), y2: Math.max(s.y1, s.y2) });
  }
  for (const ar of doc.arrows.values()) {
    if (opts.hiddenObjects?.has(ar.id)) continue;
    const bx = arrowPrims(ar, st, ar.color ?? opts.ink, prims);
    objectBoxes.set(ar.id, bx);
  }
  for (const c of doc.curved.values()) {
    if (opts.hiddenObjects?.has(c.id)) continue;
    const g = curvedGeometry(doc, c, labelBoxes);
    if (!g) continue;
    curvedPrims(g, c, st, c.color ?? '#d6336c', prims);
  }
  for (const t of doc.texts.values()) {
    if (opts.hiddenObjects?.has(t.id)) continue;
    const bx = textPrims(t, st, t.color ?? opts.ink, prims);
    objectBoxes.set(t.id, bx);
  }

  return { prims, labelBoxes, objectBoxes, mol, index, rings, hCount, bounds: primBounds(prims) };
}

// ───────────────────────── labels ─────────────────────────

function labelVisible(mol: Mol, i: number, deg: number, st: DocStyle, h: number): boolean {
  const a = mol.atoms[i];
  if (a.abbrev || a.alias) return true;
  if (a.el !== 'C') return true;
  if (deg === 0) return true;
  if (a.isotope) return true;
  if (a.map) return false;
  if (st.showCarbons === 'all') return true;
  if (st.showCarbons === 'terminal' && deg === 1) return true;
  if (a.hCount !== undefined) {
    // explicit H count that differs from what the skeleton implies → show
    const save = a.hCount;
    delete a.hCount;
    const def = implicitH(mol, i);
    a.hCount = save;
    if (def !== h) return true;
  }
  return false;
}

export function chargeText(c: number): string {
  if (!c) return '';
  return (Math.abs(c) > 1 ? String(Math.abs(c)) : '') + (c > 0 ? '+' : '−');
}

/** Choose where implicit hydrogens go relative to the symbol. */
function chooseHSide(dirs: Pt[], el: string): 'left' | 'right' | 'up' | 'down' {
  if (!dirs.length) return ['O', 'S', 'Se', 'Te', 'F', 'Cl', 'Br', 'I'].includes(el) ? 'left' : 'right';
  let right = false, left = false, up = false, down = false;
  for (const d of dirs) {
    if (d.x > 0.35) right = true;
    if (d.x < -0.35) left = true;
    if (d.y < -0.6) up = true;
    if (d.y > 0.6) down = true;
  }
  if (dirs.length === 1) return dirs[0].x > 0.1 ? 'left' : 'right';
  if (!right) return 'right';
  if (!left) return 'left';
  if (!down) return 'down';
  if (!up) return 'up';
  return 'right';
}

interface LabelLayout {
  runs: PosRun[];
  boxes: Box[];
  hSide: 'left' | 'right' | 'up' | 'down' | null;
}

function tokenRuns(text: string): TextRun[] {
  // digits following letters/brackets become subscripts; trailing +/- superscript
  return formulaRuns(text);
}

function layoutLabel(mol: Mol, i: number, dirs: Pt[], h: number, st: DocStyle): LabelLayout {
  const a = mol.atoms[i];
  const fs = st.fontSize;
  const fam = st.fontFamily;
  const m = st.labelMargin;
  const capH = fs * CAP;
  const base = a.y + capH / 2;
  const runs: PosRun[] = [];
  const boxes: Box[] = [];

  if (a.abbrev || a.alias) {
    let text = a.abbrev ?? a.alias!;
    const fromRight = dirs.length === 1 && dirs[0].x > 0.3;
    if (a.abbrev && fromRight) text = reverseLabel(text);
    const tr = tokenRuns(text);
    // measure runs; anchor = first (or last when reversed) capital-letter token
    const widths = tr.map((r) => measureText(r.text, r.sub || r.sup ? fs * SCRIPT_SCALE : fs, fam));
    const total = widths.reduce((s, x) => s + x, 0);
    // anchor width: first element-like token
    const firstTok = /^[A-Z][a-z]?/.exec(tr[0]?.text ?? '')?.[0] ?? tr[0]?.text ?? '';
    const anchorW = measureText(firstTok, fs, fam);
    let x0: number;
    if (a.abbrev && fromRight) {
      const lastRun = tr[tr.length - 1];
      const lastTok = /[A-Z][a-z]?$/.exec(lastRun.text)?.[0] ?? lastRun.text;
      const lw = measureText(lastTok, fs, fam);
      x0 = a.x - (total - lw / 2);
      if (lastRun.sub || lastRun.sup) x0 = a.x - total + widths[widths.length - 1] / 2;
    } else if (dirs.length === 0 || a.alias) {
      x0 = a.x - (dirs.length === 0 ? total / 2 : anchorW / 2);
    } else {
      x0 = a.x - anchorW / 2;
    }
    let x = x0;
    tr.forEach((r, k) => {
      const size = r.sub || r.sup ? fs * SCRIPT_SCALE : fs;
      const y = r.sub ? base + fs * SUB_SHIFT : r.sup ? base - fs * SUP_SHIFT : base;
      runs.push({ x, y, text: r.text, size });
      x += widths[k];
    });
    if (a.charge) {
      // charge after the label: OMe⁻, E⁺
      const txt = chargeText(a.charge);
      const size = fs * SCRIPT_SCALE;
      runs.push({ x, y: base - fs * SUP_SHIFT, text: txt, size });
      x += measureText(txt, size, fam);
    }
    boxes.push({ x1: x0 - m, y1: a.y - capH / 2 - m, x2: x + m, y2: a.y + capH / 2 + m * 1.6 });
    return { runs, boxes, hSide: null };
  }

  const sym = a.el === '*' ? '*' : a.el;
  const wMain = measureText(sym, fs, fam);
  const mainL = a.x - wMain / 2;
  const mainR = a.x + wMain / 2;
  runs.push({ x: mainL, y: base, text: sym, size: fs });
  boxes.push({
    x1: mainL - m, y1: a.y - capH / 2 - m, x2: mainR + m, y2: a.y + capH / 2 + m,
    shape: sym.length === 1 ? 'ellipse' : 'rect',
  });
  // enlarge ellipse a bit so diagonal bonds keep the same gap as horizontal ones
  if (sym.length === 1) {
    const bx = boxes[0];
    const cx = a.x, cy = a.y;
    const rx = (bx.x2 - bx.x1) / 2 * 1.18, ry = (bx.y2 - bx.y1) / 2 * 1.18;
    boxes[0] = { x1: cx - rx, y1: cy - ry, x2: cx + rx, y2: cy + ry, shape: 'ellipse' };
  }

  let rightEdge = mainR;
  let leftEdge = mainL;
  let hSide: LabelLayout['hSide'] = null;
  if (h > 0 && st.showImplicitH) {
    hSide = chooseHSide(dirs, a.el);
    const wH = measureText('H', fs, fam);
    const cnt = h > 1 ? String(h) : '';
    const wC = cnt ? measureText(cnt, fs * SCRIPT_SCALE, fam) : 0;
    if (hSide === 'right') {
      runs.push({ x: mainR, y: base, text: 'H', size: fs });
      if (cnt) runs.push({ x: mainR + wH, y: base + fs * SUB_SHIFT, text: cnt, size: fs * SCRIPT_SCALE });
      boxes.push({ x1: mainR - m * 0.3, y1: a.y - capH / 2 - m, x2: mainR + wH + wC + m, y2: a.y + capH / 2 + m * (cnt ? 2 : 1) });
      rightEdge = mainR + wH + wC;
    } else if (hSide === 'left') {
      const x0 = mainL - wH - wC;
      runs.push({ x: x0, y: base, text: 'H', size: fs });
      if (cnt) runs.push({ x: x0 + wH, y: base + fs * SUB_SHIFT, text: cnt, size: fs * SCRIPT_SCALE });
      boxes.push({ x1: x0 - m, y1: a.y - capH / 2 - m, x2: mainL + m * 0.3, y2: a.y + capH / 2 + m * (cnt ? 2 : 1) });
      leftEdge = x0;
    } else {
      const dy = (hSide === 'down' ? 1 : -1) * fs * 1.05;
      const x0 = a.x - wH / 2;
      const b2 = base + dy;
      runs.push({ x: x0, y: b2, text: 'H', size: fs });
      if (cnt) runs.push({ x: x0 + wH, y: b2 + fs * SUB_SHIFT, text: cnt, size: fs * SCRIPT_SCALE });
      boxes.push({ x1: x0 - m, y1: a.y + dy - capH / 2 - m, x2: x0 + wH + wC + m, y2: a.y + dy + capH / 2 + m });
    }
  }
  if (a.charge) {
    const txt = chargeText(a.charge);
    const size = fs * SCRIPT_SCALE;
    const x = hSide === 'left' || hSide === 'up' || hSide === 'down' || !hSide ? mainR : rightEdge;
    runs.push({ x, y: base - fs * SUP_SHIFT, text: txt, size });
    const tw = measureText(txt, size, fam);
    boxes.push({ x1: x - m * 0.3, y1: a.y - capH / 2 - fs * 0.45 - m, x2: x + tw + m * 0.5, y2: a.y });
  }
  if (a.isotope) {
    const txt = String(a.isotope);
    const size = fs * SCRIPT_SCALE;
    const tw = measureText(txt, size, fam);
    const x = leftEdge - tw;
    runs.push({ x, y: base - fs * SUP_SHIFT, text: txt, size });
    boxes.push({ x1: x - m * 0.5, y1: a.y - capH / 2 - fs * 0.45 - m, x2: leftEdge, y2: a.y });
  }
  return { runs, boxes, hSide };
}

function unionBox(bs: Box[]): Box {
  return bs.reduce(
    (u, b) => ({ x1: Math.min(u.x1, b.x1), y1: Math.min(u.y1, b.y1), x2: Math.max(u.x2, b.x2), y2: Math.max(u.y2, b.y2) }),
    { x1: Infinity, y1: Infinity, x2: -Infinity, y2: -Infinity },
  );
}

function labelRadius(bs: Box[], c: Pt): number {
  const b = bs[0];
  return Math.max(c.x - b.x1, b.x2 - c.x, c.y - b.y1, b.y2 - c.y) * 0.85;
}

/** Spreads `count` directions into the largest free angular gaps around an atom. */
export function freeAngles(occupied: number[], count: number, preferred: number): number[] {
  const out: number[] = [];
  const occ = occupied.map((a) => ((a % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI));
  for (let k = 0; k < count; k++) {
    if (!occ.length) {
      const a = preferred + (k * 2 * Math.PI) / Math.max(count, 1);
      out.push(a);
      occ.push(((a % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI));
      continue;
    }
    const s = [...occ].sort((p, q) => p - q);
    let bestGap = -1, bestMid = preferred;
    for (let j = 0; j < s.length; j++) {
      const a0 = s[j];
      const a1 = j + 1 < s.length ? s[j + 1] : s[0] + 2 * Math.PI;
      const gap = a1 - a0;
      const mid = a0 + gap / 2;
      // tie-break toward the preferred direction
      const pref = Math.cos(mid - preferred) * 0.01;
      if (gap + pref > bestGap) {
        bestGap = gap + pref;
        bestMid = mid;
      }
    }
    out.push(bestMid);
    occ.push(((bestMid % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI));
  }
  return out;
}

// ───────────────────────── bonds ─────────────────────────

function doubleBondPlacement(
  mol: Mol,
  bi: number,
  rings: RingInfo,
  pos: Pt[],
  dbPos: 'auto' | 'left' | 'right' | 'center',
  boxes: (Box[] | null)[],
): { mode: 'center' | 'side'; side: number } {
  if (dbPos === 'center') return { mode: 'center', side: 1 };
  if (dbPos === 'left') return { mode: 'side', side: -1 };
  if (dbPos === 'right') return { mode: 'side', side: 1 };
  const b = mol.bonds[bi];
  const A = pos[b.a], B = pos[b.b];
  const u = norm(sub(B, A));
  const pv = perp(u);
  const mid = lerp(A, B, 0.5);
  const rs = rings.bondRings[bi];
  if (rs.length) {
    // prefer the smallest ring with the most double bonds
    let best = rs[0];
    let bestScore = -Infinity;
    for (const r of rs) {
      const dbl = rings.ringBonds[r].filter((x) => mol.bonds[x].order === 2).length;
      const score = dbl * 10 - rings.rings[r].length;
      if (score > bestScore) {
        bestScore = score;
        best = r;
      }
    }
    const ring = rings.rings[best];
    const c = mul(ring.reduce((s, i) => add(s, pos[i]), { x: 0, y: 0 }), 1 / ring.length);
    return { mode: 'side', side: dot(pv, sub(c, mid)) >= 0 ? 1 : -1 };
  }
  const nA = mol.neighbors(b.a).filter((j) => j !== b.b);
  const nB = mol.neighbors(b.b).filter((j) => j !== b.a);
  if (nA.length === 0 && nB.length === 0) return { mode: 'center', side: 1 };
  if ((nA.length === 0 && nB.length === 2) || (nB.length === 0 && nA.length === 2)) return { mode: 'center', side: 1 };
  // terminal labeled atom (C=O) with a single other substituent on the far end → still offset
  if ((nA.length === 0 && boxes[b.a]) || (nB.length === 0 && boxes[b.b])) {
    if (nA.length >= 2 || nB.length >= 2) return { mode: 'center', side: 1 };
  }
  let score = 0;
  for (const j of nA) score += Math.sign(dot(pv, sub(pos[j], A)));
  for (const j of nB) score += Math.sign(dot(pv, sub(pos[j], B)));
  if (score === 0) {
    const ref = nA[0] ?? nB[0];
    const base = nA.length ? A : B;
    score = Math.sign(dot(pv, sub(pos[ref], base))) || 1;
  }
  return { mode: 'side', side: score > 0 ? 1 : -1 };
}

/** Distance to shorten the inner line of an offset double bond at the end whose bond direction is `u`. */
function innerTrim(dirs: Pt[], u: Pt, side: number, sp: number, _id: number, pv: Pt): number {
  // neighbour bonds on the offset side
  let best = Infinity;
  for (const d of dirs) {
    if (dot(d, u) > 0.999) continue; // the bond itself
    const s = dot(d, mul(pv, side));
    if (s <= 0.05) continue;
    const ang = Math.acos(Math.max(-1, Math.min(1, dot(d, u))));
    best = Math.min(best, ang);
  }
  if (best === Infinity) return sp * 0.5;
  const t = sp / Math.tan(best / 2);
  return Math.max(0.04, Math.min(0.4, t));
}

// ───────────────────────── arrows ─────────────────────────

export function arrowHeadPoly(tip: Pt, u: Pt, length: number, halfW: number, color: string, kind: 'full' | 'halfLeft' | 'halfRight'): Prim {
  const pv = perp(u);
  const back = sub(tip, mul(u, length));
  const notch = sub(tip, mul(u, length * 0.75));
  const l = add(back, mul(pv, -halfW));
  const r = add(back, mul(pv, halfW));
  let pts: number[];
  if (kind === 'full') pts = [tip.x, tip.y, l.x, l.y, notch.x, notch.y, r.x, r.y];
  else if (kind === 'halfLeft') pts = [tip.x, tip.y, l.x, l.y, notch.x, notch.y];
  else pts = [tip.x, tip.y, notch.x, notch.y, r.x, r.y];
  return { k: 'poly', pts, fill: color, stroke: color, w: 0.01, closed: true };
}

const HEAD_LEN = 0.3;
const HEAD_HW = 0.11;

export function arrowPrims(ar: ArrowObj, st: DocStyle, color: string, prims: Prim[]): Box {
  const P = { x: ar.x1, y: ar.y1 }, Q = { x: ar.x2, y: ar.y2 };
  const L = len(sub(Q, P));
  const u = L > 1e-6 ? norm(sub(Q, P)) : { x: 1, y: 0 };
  const pv = perp(u);
  const w = st.lineWidth;
  const ln = (a: Pt, b: Pt, dash?: number[]) => prims.push({ k: 'line', x1: a.x, y1: a.y, x2: b.x, y2: b.y, w, color, dash, cap: 'butt' });
  switch (ar.kind) {
    case 'reaction':
    case 'dashed':
    case 'noGo': {
      ln(P, sub(Q, mul(u, HEAD_LEN * 0.7)), ar.kind === 'dashed' ? [0.12, 0.08] : undefined);
      prims.push(arrowHeadPoly(Q, u, HEAD_LEN, HEAD_HW, color, 'full'));
      if (ar.kind === 'noGo') {
        const m = lerp(P, Q, 0.5);
        const s = 0.18;
        const d1 = { x: (u.x + pv.x) * s, y: (u.y + pv.y) * s };
        const d2 = { x: (u.x - pv.x) * s, y: (u.y - pv.y) * s };
        ln(sub(m, d1), add(m, d1));
        ln(sub(m, d2), add(m, d2));
      }
      break;
    }
    case 'resonance': {
      ln(add(P, mul(u, HEAD_LEN * 0.7)), sub(Q, mul(u, HEAD_LEN * 0.7)));
      prims.push(arrowHeadPoly(Q, u, HEAD_LEN, HEAD_HW, color, 'full'));
      prims.push(arrowHeadPoly(P, mul(u, -1), HEAD_LEN, HEAD_HW, color, 'full'));
      break;
    }
    case 'equilibrium':
    case 'unbalancedEq': {
      const o = 0.07;
      const top = mul(pv, -o), bot = mul(pv, o);
      const shrink = ar.kind === 'unbalancedEq' ? L * 0.2 : 0;
      ln(add(P, top), sub(add(Q, top), mul(u, HEAD_LEN * 0.6)));
      prims.push(arrowHeadPoly(add(Q, top), u, HEAD_LEN, HEAD_HW * 1.1, color, 'halfLeft'));
      const p2 = add(add(P, mul(u, shrink)), bot), q2 = add(sub(Q, mul(u, shrink)), bot);
      ln(add(p2, mul(u, HEAD_LEN * 0.6)), q2);
      prims.push(arrowHeadPoly(p2, mul(u, -1), HEAD_LEN, HEAD_HW * 1.1, color, 'halfLeft'));
      break;
    }
    case 'retro': {
      const o = 0.065;
      const headBack = HEAD_LEN * 0.9;
      ln(add(P, mul(pv, -o)), add(sub(Q, mul(u, headBack)), mul(pv, -o)));
      ln(add(P, mul(pv, o)), add(sub(Q, mul(u, headBack)), mul(pv, o)));
      const b = sub(Q, mul(u, HEAD_LEN * 1.1));
      const l = add(b, mul(pv, -HEAD_HW * 1.6)), r = add(b, mul(pv, HEAD_HW * 1.6));
      prims.push({ k: 'poly', pts: [l.x, l.y, Q.x, Q.y, r.x, r.y], stroke: color, w, closed: false });
      break;
    }
    case 'line':
    default:
      ln(P, Q);
  }
  let box: Box = { x1: Math.min(P.x, Q.x) - 0.15, y1: Math.min(P.y, Q.y) - 0.15, x2: Math.max(P.x, Q.x) + 0.15, y2: Math.max(P.y, Q.y) + 0.15 };
  // captions
  const mid = lerp(P, Q, 0.5);
  const size = st.fontSize * 0.85;
  // "above" = the side facing screen-up for left-to-right arrows
  let upSide = mul(pv, -1);
  if (upSide.y > 0.01 || (Math.abs(upSide.y) < 0.01 && upSide.x > 0)) upSide = mul(upSide, -1);
  const vertical = Math.abs(u.x) < 0.5;
  const place = (text: string, sideSign: number) => {
    const lines = text.split('\n');
    const lh = size * 1.2;
    const blockH = lines.length * lh;
    if (vertical) {
      // text to the right (above) or left (below) of vertical arrows
      const sx = sideSign > 0 ? 1 : -1;
      const x = mid.x + sx * 0.2;
      let y = mid.y - blockH / 2 + size * CAP;
      for (const l of lines) {
        const bx = textLine(l, x, y, size, st.fontFamily, color, sx > 0 ? 'left' : 'right', prims, true);
        box = unionBox([box, bx]);
        y += lh;
      }
      return;
    }
    const gap = 0.16;
    const start = add(mid, mul(upSide, sideSign * gap));
    let y = sideSign > 0 ? start.y - blockH + size * CAP + (lh - size) * 0 : start.y + size * CAP;
    for (const l of lines) {
      const bx = textLine(l, start.x, y, size, st.fontFamily, color, 'center', prims, true);
      box = unionBox([box, bx]);
      y += lh;
    }
  };
  if (ar.above) place(ar.above, 1);
  if (ar.below) place(ar.below, -1);
  return box;
}

/** Draws one line of auto-formatted text; returns its box. */
function textLine(text: string, x: number, baseline: number, size: number, fam: string, color: string, align: 'left' | 'center' | 'right', prims: Prim[], formula: boolean): Box {
  const runs = richLines(text, false)[0];
  // inside reaction captions, auto-subscript formula-like words (H2O, CH2Cl2) but keep plain words
  const expanded: TextRun[] = [];
  for (const r of runs) {
    if (formula && !r.sub && !r.sup) {
      for (const part of r.text.split(/(\s+|,|;)/)) {
        if (/^[A-Z][A-Za-z0-9()]*\d[A-Za-z0-9()]*[+\-]?$/.test(part)) expanded.push(...formulaRuns(part));
        else if (part) expanded.push({ ...r, text: part });
      }
    } else expanded.push(r);
  }
  let total = 0;
  const ws = expanded.map((r) => {
    const s = r.sub || r.sup ? size * SCRIPT_SCALE : size;
    const wv = measureText(r.text, s, fam, r.bold, r.italic);
    total += wv;
    return wv;
  });
  let cx = align === 'left' ? x : align === 'right' ? x - total : x - total / 2;
  const x0 = cx;
  expanded.forEach((r, k) => {
    const s = r.sub || r.sup ? size * SCRIPT_SCALE : size;
    const y = r.sub ? baseline + size * SUB_SHIFT : r.sup ? baseline - size * SUP_SHIFT : baseline;
    if (r.text.trim()) prims.push({ k: 'text', x: cx, y, text: r.text, size: s, family: fam, color, bold: r.bold, italic: r.italic });
    cx += ws[k];
  });
  return { x1: x0 - 0.05, y1: baseline - size * 0.9, x2: cx + 0.05, y2: baseline + size * 0.3 };
}

export function textPrims(t: TextObj, st: DocStyle, color: string, prims: Prim[]): Box {
  const size = st.fontSize * (t.size ?? 1);
  const lines = richLines(t.text || ' ', !!t.formula);
  const lh = size * 1.25;
  let y = t.y;
  let box: Box | null = null;
  for (const runs of lines) {
    let total = 0;
    const ws = runs.map((r) => {
      const s = r.sub || r.sup ? size * SCRIPT_SCALE : size;
      const wv = measureText(r.text, s, st.fontFamily, r.bold || t.bold, r.italic || t.italic);
      total += wv;
      return wv;
    });
    const align = t.align ?? 'left';
    let cx = align === 'left' ? t.x : align === 'right' ? t.x - total : t.x - total / 2;
    const x0 = cx;
    runs.forEach((r, k) => {
      const s = r.sub || r.sup ? size * SCRIPT_SCALE : size;
      const yy = r.sub ? y + size * SUB_SHIFT : r.sup ? y - size * SUP_SHIFT : y;
      prims.push({ k: 'text', x: cx, y: yy, text: r.text, size: s, family: st.fontFamily, color, bold: r.bold || t.bold, italic: r.italic || t.italic });
      cx += ws[k];
    });
    const lb: Box = { x1: x0 - 0.06, y1: y - size * 0.9, x2: Math.max(cx, x0 + size * 0.3) + 0.06, y2: y + size * 0.3 };
    box = box ? unionBox([box, lb]) : lb;
    y += lh;
  }
  return box!;
}

// ───────────────────────── curved arrows ─────────────────────────

export interface CurvedGeom {
  p0: Pt;
  c1: Pt;
  c2: Pt;
  p3: Pt;
  rawStart: Pt;
  rawEnd: Pt;
}

export function anchorPoint(doc: ChemDoc, an: CurvedArrowObj['from']): Pt | null {
  switch (an.type) {
    case 'atom': {
      const a = doc.atoms.get(an.id);
      return a ? { x: a.x, y: a.y } : null;
    }
    case 'bond': {
      const b = doc.bonds.get(an.id);
      if (!b) return null;
      const A = doc.atoms.get(b.a), B = doc.atoms.get(b.b);
      if (!A || !B) return null;
      return lerp(A, B, 0.5);
    }
    case 'between': {
      const A = doc.atoms.get(an.a), B = doc.atoms.get(an.b);
      if (!A || !B) return null;
      return lerp(A, B, 0.5);
    }
    case 'point':
      return { x: an.x, y: an.y };
  }
}

export function curvedGeometry(doc: ChemDoc, c: CurvedArrowObj, labelBoxes?: Map<number, Box[]>): CurvedGeom | null {
  const S = anchorPoint(doc, c.from);
  const E = anchorPoint(doc, c.to);
  if (!S || !E) return null;
  let d = sub(E, S);
  let L = len(d);
  if (L < 1e-6) {
    d = { x: 1, y: 0 };
    L = 0.6;
  }
  const u = mul(d, 1 / len(d));
  const v = perp(u);
  // short arrows (bond → own atom) still need a visible arc
  const Lh = Math.max(L, 0.9);
  const cp1 = add(S, add(mul(u, c.c1.t * L), mul(v, c.c1.h * Lh)));
  const cp2 = add(S, add(mul(u, c.c2.t * L), mul(v, c.c2.h * Lh)));
  const trim = (an: CurvedArrowObj['from'], p: Pt, toward: Pt, isEnd: boolean): Pt => {
    const dir = norm(sub(toward, p));
    let r = 0;
    if (an.type === 'atom') {
      const bx = labelBoxes?.get(an.id);
      r = bx ? Math.max(0.3, (bx[0].x2 - bx[0].x1) / 2 + 0.08) : isEnd ? 0.24 : 0.2;
    } else if (an.type === 'bond') r = 0.08;
    else if (an.type === 'between') r = 0.05;
    return add(p, mul(dir, r));
  };
  const p0 = trim(c.from, S, cp1, false);
  const p3 = trim(c.to, E, cp2, true);
  return { p0, c1: cp1, c2: cp2, p3, rawStart: S, rawEnd: E };
}

function curvedPrims(g: CurvedGeom, c: CurvedArrowObj, st: DocStyle, color: string, prims: Prim[]): void {
  const w = st.lineWidth * 1.05;
  const tan = norm(bezierTangent(g.p0, g.c1, g.c2, g.p3, 1));
  const headLen = 0.24;
  // stop the shaft a little before the tip so the head looks crisp
  const shaftEnd = sub(g.p3, mul(tan, headLen * 0.6));
  const c2 = add(g.c2, sub(shaftEnd, g.p3));
  prims.push({ k: 'path', d: [['M', g.p0.x, g.p0.y], ['C', g.c1.x, g.c1.y, c2.x, c2.y, shaftEnd.x, shaftEnd.y]], stroke: color, w, cap: 'round' });
  if (c.electrons === 2) prims.push(arrowHeadPoly(g.p3, tan, headLen, 0.1, color, 'full'));
  else {
    // fishhook: barb on the outside of the curve
    const bend = cross(sub(g.c2, g.p0), sub(g.p3, g.p0));
    prims.push(arrowHeadPoly(g.p3, tan, headLen, 0.11, color, bend > 0 ? 'halfRight' : 'halfLeft'));
  }
}

// ───────────────────────── shapes ─────────────────────────

function shapePrims(s: ShapeObj, st: DocStyle, ink: string, prims: Prim[]): void {
  const color = s.color ?? ink;
  const w = st.lineWidth;
  const x1 = Math.min(s.x1, s.x2), x2 = Math.max(s.x1, s.x2);
  const y1 = Math.min(s.y1, s.y2), y2 = Math.max(s.y1, s.y2);
  const dash = s.dashed ? [0.1, 0.07] : undefined;
  const tick = Math.min(0.2, (x2 - x1) * 0.15);
  switch (s.kind) {
    case 'rect':
      prims.push({ k: 'poly', pts: [x1, y1, x2, y1, x2, y2, x1, y2], stroke: color, fill: s.fill, w, closed: true, dash });
      break;
    case 'roundRect': {
      const r = Math.min(0.25, (x2 - x1) / 4, (y2 - y1) / 4);
      prims.push({
        k: 'path',
        d: [
          ['M', x1 + r, y1], ['L', x2 - r, y1], ['Q', x2, y1, x2, y1 + r], ['L', x2, y2 - r], ['Q', x2, y2, x2 - r, y2],
          ['L', x1 + r, y2], ['Q', x1, y2, x1, y2 - r], ['L', x1, y1 + r], ['Q', x1, y1, x1 + r, y1], ['Z'],
        ],
        stroke: color, fill: s.fill, w, dash,
      });
      break;
    }
    case 'ellipse': {
      const cx = (x1 + x2) / 2, cy = (y1 + y2) / 2, rx = (x2 - x1) / 2, ry = (y2 - y1) / 2;
      const k = 0.5523;
      prims.push({
        k: 'path',
        d: [
          ['M', cx + rx, cy],
          ['C', cx + rx, cy + ry * k, cx + rx * k, cy + ry, cx, cy + ry],
          ['C', cx - rx * k, cy + ry, cx - rx, cy + ry * k, cx - rx, cy],
          ['C', cx - rx, cy - ry * k, cx - rx * k, cy - ry, cx, cy - ry],
          ['C', cx + rx * k, cy - ry, cx + rx, cy - ry * k, cx + rx, cy],
          ['Z'],
        ],
        stroke: color, fill: s.fill, w, dash,
      });
      break;
    }
    case 'line':
      prims.push({ k: 'line', x1: s.x1, y1: s.y1, x2: s.x2, y2: s.y2, w, color, dash });
      break;
    case 'bracket':
    case 'tsBracket':
      prims.push({ k: 'poly', pts: [x1 + tick, y1, x1, y1, x1, y2, x1 + tick, y2], stroke: color, w, closed: false });
      prims.push({ k: 'poly', pts: [x2 - tick, y1, x2, y1, x2, y2, x2 - tick, y2], stroke: color, w, closed: false });
      break;
    case 'paren': {
      const bow = Math.min(0.25, (x2 - x1) * 0.12);
      prims.push({ k: 'path', d: [['M', x1 + bow, y1], ['Q', x1 - bow, (y1 + y2) / 2, x1 + bow, y2]], stroke: color, w });
      prims.push({ k: 'path', d: [['M', x2 - bow, y1], ['Q', x2 + bow, (y1 + y2) / 2, x2 - bow, y2]], stroke: color, w });
      break;
    }
    case 'brace': {
      const bw = Math.min(0.18, (x2 - x1) * 0.1);
      const my = (y1 + y2) / 2;
      prims.push({ k: 'path', d: [['M', x1 + bw, y1], ['Q', x1, y1, x1, y1 + bw], ['L', x1, my - bw], ['Q', x1, my, x1 - bw, my], ['Q', x1, my, x1, my + bw], ['L', x1, y2 - bw], ['Q', x1, y2, x1 + bw, y2]], stroke: color, w });
      prims.push({ k: 'path', d: [['M', x2 - bw, y1], ['Q', x2, y1, x2, y1 + bw], ['L', x2, my - bw], ['Q', x2, my, x2 + bw, my], ['Q', x2, my, x2, my + bw], ['L', x2, y2 - bw], ['Q', x2, y2, x2 - bw, y2]], stroke: color, w });
      break;
    }
    case 'orbitalP': {
      const cx = (x1 + x2) / 2, cy = (y1 + y2) / 2;
      const hw = (x2 - x1) / 2, hh = (y2 - y1) / 2;
      const lobe = (dir: number, fill?: string) => {
        prims.push({
          k: 'path',
          d: [['M', cx, cy], ['C', cx - hw * 1.1, cy + dir * hh * 0.35, cx - hw, cy + dir * hh, cx, cy + dir * hh], ['C', cx + hw, cy + dir * hh, cx + hw * 1.1, cy + dir * hh * 0.35, cx, cy], ['Z']],
          stroke: color, fill, w,
        });
      };
      lobe(-1, s.fill ?? 'rgba(120,120,120,0.45)');
      lobe(1);
      break;
    }
    case 'orbitalS': {
      const cx = (x1 + x2) / 2, cy = (y1 + y2) / 2;
      prims.push({ k: 'circle', x: cx, y: cy, r: Math.min(x2 - x1, y2 - y1) / 2, stroke: color, fill: s.fill, w });
      break;
    }
  }
  const label = s.kind === 'tsBracket' ? s.label ?? '‡' : s.label;
  if (label) {
    const size = st.fontSize * 0.9;
    prims.push({ k: 'text', x: x2 + 0.06, y: y1 + size * 0.55, text: label, size, family: st.fontFamily, color });
  }
}

export function primBounds(prims: Prim[]): Box | null {
  let x1 = Infinity, y1 = Infinity, x2 = -Infinity, y2 = -Infinity;
  const inc = (x: number, y: number, pad = 0) => {
    x1 = Math.min(x1, x - pad); y1 = Math.min(y1, y - pad);
    x2 = Math.max(x2, x + pad); y2 = Math.max(y2, y + pad);
  };
  for (const p of prims) {
    switch (p.k) {
      case 'line': inc(p.x1, p.y1, p.w); inc(p.x2, p.y2, p.w); break;
      case 'poly': for (let i = 0; i < p.pts.length; i += 2) inc(p.pts[i], p.pts[i + 1]); break;
      case 'path': for (const s of p.d) if (s[0] !== 'Z') { const n = s.length; inc(s[n - 2] as number, s[n - 1] as number); if (s[0] === 'C' || s[0] === 'Q') inc(s[1] as number, s[2] as number); } break;
      case 'circle': inc(p.x, p.y, p.r); break;
      case 'text': inc(p.x, p.y - p.size * 0.8); inc(p.x + measureText(p.text, p.size, p.family, p.bold, p.italic), p.y + p.size * 0.25); break;
    }
  }
  if (x1 === Infinity) return null;
  return { x1, y1, x2, y2 };
}
