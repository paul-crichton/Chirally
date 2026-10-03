// Interactive tools (ChemDraw-style). Each tool translates pointer gestures into document operations.
import type { Editor } from './editor';
import { Hit, PEvent, Tool, ToolId } from './types';
import {
  addAtom, addBond, adjacency, bondBetween, cloneDoc, neighborsOf,
} from '../doc/document';
import {
  idealBondDirection, bondToPoint, applyBondType, atomNear, deleteSelection, emptySelection, fuseOverlaps, pendingMerges,
  placeRing, fuseRingOnBond, attachRingToAtom, placeChair, insertTemplateAt, attachTemplate, moveSelection,
  rotateSelection, selectedAtomIds, chainPoints, setElement, Selection, MERGE_TOL, selectionCenter,
} from './ops';
import { Pt, add, sub, mul, norm, dist, fromAngle, angleOf, snapAngle, pointInPolygon, perp, lerp, len } from '../render/geom';
import { ChemDoc, CurvedArrowObj, Anchor } from '../doc/types';
import { curvedGeometry } from '../render/scene';

const DEG15 = Math.PI / 12;

// ───────────────────────── helpers ─────────────────────────

function ghostColor(ed: Editor): string {
  return ed.theme.dark ? 'rgba(120,190,255,0.75)' : 'rgba(25,113,194,0.6)';
}

function drawGhost(ed: Editor, ctx: CanvasRenderingContext2D, segs: [Pt, Pt][], dots: Pt[] = [], dashed = false): void {
  ctx.save();
  ctx.strokeStyle = ghostColor(ed);
  ctx.fillStyle = ghostColor(ed);
  ctx.lineWidth = Math.max(1.5, ed.doc.style.lineWidth * ed.view.scale);
  ctx.lineCap = 'round';
  if (dashed) ctx.setLineDash([5, 4]);
  for (const [a, b] of segs) {
    const p = ed.toScreen(a), q = ed.toScreen(b);
    ctx.beginPath();
    ctx.moveTo(p.x, p.y);
    ctx.lineTo(q.x, q.y);
    ctx.stroke();
  }
  for (const d of dots) {
    const p = ed.toScreen(d);
    ctx.beginPath();
    ctx.arc(p.x, p.y, 3, 0, Math.PI * 2);
    ctx.fill();
  }
  ctx.restore();
}

/** Draws atoms/bonds that would be added by an operation (diff between docs). */
function drawDocDiff(ed: Editor, ctx: CanvasRenderingContext2D, before: ChemDoc, after: ChemDoc): void {
  const segs: [Pt, Pt][] = [];
  for (const b of after.bonds.values()) {
    const old = before.bonds.get(b.id);
    if (old && old.order === b.order) continue;
    const A = after.atoms.get(b.a)!, B = after.atoms.get(b.b)!;
    segs.push([A, B]);
    if (b.order >= 2) {
      // draw an inner line for double bonds in the ghost
      const u = norm(sub(B, A));
      const pv = perp(u);
      // offset toward the centroid of new atoms
      const off = mul(pv, 0.18);
      segs.push([add(lerp(A, B, 0.15), off), add(lerp(A, B, 0.85), off)]);
    }
  }
  drawGhost(ed, ctx, segs);
}

function snapToGrid(ed: Editor, p: Pt): Pt {
  if (!ed.grid.snap) return p;
  const g = ed.grid.size;
  return { x: Math.round(p.x / g) * g, y: Math.round(p.y / g) * g };
}

/** End point for a bond dragged from `from` towards the cursor: snaps to atoms, else 15° angle + unit length. */
function snapBondEnd(ed: Editor, from: Pt, cursor: Pt, fromAtom: number | undefined, alt: boolean): { p: Pt; atom?: number } {
  const ex = atomNear(ed.doc, cursor, 0.38, fromAtom !== undefined ? new Set([fromAtom]) : undefined);
  if (ex) return { p: { x: ex.x, y: ex.y }, atom: ex.id };
  const d = sub(cursor, from);
  if (len(d) < 1e-6) return { p: add(from, fromAngle(-Math.PI / 6)) };
  let ang = angleOf(d);
  if (!alt) ang = snapAngle(ang, DEG15);
  const p = add(from, fromAngle(ang, 1));
  const ex2 = atomNear(ed.doc, p, MERGE_TOL, fromAtom !== undefined ? new Set([fromAtom]) : undefined);
  if (ex2) return { p: { x: ex2.x, y: ex2.y }, atom: ex2.id };
  return { p };
}

function anchorOfHit(hit: Hit): Anchor | null {
  if (!hit) return null;
  if (hit.kind === 'atom') return { type: 'atom', id: hit.id };
  if (hit.kind === 'bond') return { type: 'bond', id: hit.id };
  return null;
}

// ───────────────────────── Select / Lasso ─────────────────────────

type SelMode = 'none' | 'move' | 'rotate' | 'marquee' | 'lasso' | 'arrowEnd' | 'curvedCp' | 'shapeCorner';

class SelectTool implements Tool {
  id: ToolId;
  cursor = 'default';
  private mode: SelMode = 'none';
  private last: Pt = { x: 0, y: 0 };
  private start: Pt = { x: 0, y: 0 };
  private rect: { x1: number; y1: number; x2: number; y2: number } | null = null;
  private lasso: Pt[] = [];
  private rotCenter: Pt = { x: 0, y: 0 };
  private rotStart = 0;
  private rotApplied = 0;
  private target: Hit = null;
  private guides: { x?: number; y?: number } = {};
  private merges: [number, number][] = [];
  private moved = new Set<number>();
  private additive = false;
  private baseSel: Selection = emptySelection();
  private moveOrigin = new Map<number, Pt>();
  private totalDelta: Pt = { x: 0, y: 0 };

  constructor(private ed: Editor, lasso: boolean) {
    this.id = lasso ? 'lasso' : 'select';
  }

  hint(): string {
    return 'Click to select · drag to move · Shift-click to add · double-click an atom to edit its label · drag empty space to marquee-select';
  }

  down(e: PEvent): void {
    const ed = this.ed;
    this.start = { x: e.x, y: e.y };
    this.last = { x: e.x, y: e.y };
    this.guides = {};
    this.merges = [];
    this.totalDelta = { x: 0, y: 0 };
    const h = e.hit;
    this.target = h;
    if (h?.kind === 'rotate') {
      this.mode = 'rotate';
      this.rotCenter = selectionCenter(ed.doc, ed.sel) ?? { x: e.x, y: e.y };
      this.rotStart = angleOf(sub({ x: e.x, y: e.y }, this.rotCenter));
      this.rotApplied = 0;
      return;
    }
    if (h && h.kind === 'arrow' && h.part !== 'body' && ed.sel.objects.has(h.id)) {
      this.mode = 'arrowEnd';
      return;
    }
    if (h && h.kind === 'curved' && h.part !== 'body' && ed.sel.objects.has(h.id)) {
      this.mode = 'curvedCp';
      return;
    }
    if (h && h.kind === 'shape' && h.part !== 'body' && ed.sel.objects.has(h.id)) {
      this.mode = 'shapeCorner';
      return;
    }
    if (h) {
      const sel = ed.sel;
      const isSel = this.isSelected(h);
      if (e.shift) {
        const s = cloneSel(sel);
        this.toggle(s, h);
        ed.setSelection(s);
        this.mode = 'none';
        return;
      }
      if (!isSel) {
        const s = emptySelection();
        this.toggle(s, h);
        ed.setSelection(s);
      }
      this.mode = 'move';
      this.moved = selectedAtomIds(ed.doc, ed.sel);
      this.moveOrigin.clear();
      for (const id of this.moved) {
        const a = ed.doc.atoms.get(id)!;
        this.moveOrigin.set(id, { x: a.x, y: a.y });
      }
      return;
    }
    this.additive = e.shift;
    this.baseSel = e.shift ? cloneSel(ed.sel) : emptySelection();
    if (!e.shift) ed.clearSelection();
    this.mode = this.id === 'lasso' || e.alt ? 'lasso' : 'marquee';
    this.rect = { x1: e.x, y1: e.y, x2: e.x, y2: e.y };
    this.lasso = [{ x: e.x, y: e.y }];
  }

  private isSelected(h: Hit): boolean {
    const s = this.ed.sel;
    if (!h) return false;
    if (h.kind === 'atom') return s.atoms.has(h.id);
    if (h.kind === 'bond') return s.bonds.has(h.id);
    if (h.kind === 'rotate') return false;
    return s.objects.has(h.id);
  }

  private toggle(s: Selection, h: Hit): void {
    if (!h || h.kind === 'rotate') return;
    const flip = (set: Set<number>, id: number) => (set.has(id) ? set.delete(id) : set.add(id));
    if (h.kind === 'atom') flip(s.atoms, h.id);
    else if (h.kind === 'bond') {
      const b = this.ed.doc.bonds.get(h.id)!;
      if (s.bonds.has(h.id)) {
        s.bonds.delete(h.id);
      } else {
        s.bonds.add(h.id);
        s.atoms.add(b.a);
        s.atoms.add(b.b);
      }
    } else flip(s.objects, h.id);
  }

  drag(e: PEvent): void {
    const ed = this.ed;
    const p = { x: e.x, y: e.y };
    switch (this.mode) {
      case 'move': {
        ed.begin();
        // total displacement from the drag start, then snapping
        let dx = p.x - this.start.x, dy = p.y - this.start.y;
        const snapped = this.snapMove(dx, dy, e.alt);
        dx = snapped.dx;
        dy = snapped.dy;
        const ddx = dx - this.totalDelta.x, ddy = dy - this.totalDelta.y;
        moveSelection(ed.doc, ed.sel, ddx, ddy);
        this.totalDelta = { x: dx, y: dy };
        this.merges = pendingMerges(ed.doc, this.moved);
        ed.touch();
        break;
      }
      case 'rotate': {
        ed.begin();
        let ang = angleOf(sub(p, this.rotCenter)) - this.rotStart;
        if (e.shift) ang = snapAngle(ang, DEG15);
        rotateSelection(ed.doc, ed.sel, this.rotCenter, ang - this.rotApplied);
        this.rotApplied = ang;
        ed.setStatus(`Rotate ${Math.round((ang * 180) / Math.PI)}° (hold Shift for 15° steps)`);
        ed.touch();
        break;
      }
      case 'arrowEnd': {
        const h = this.target as Extract<Hit, { kind: 'arrow' }>;
        const a = ed.doc.arrows.get(h.id);
        if (!a) return;
        ed.begin();
        const other = h.part === 'end' ? { x: a.x1, y: a.y1 } : { x: a.x2, y: a.y2 };
        let q = p;
        if (!e.alt) {
          const d = sub(p, other);
          q = add(other, fromAngle(snapAngle(angleOf(d), DEG15), len(d)));
        }
        if (h.part === 'end') { a.x2 = q.x; a.y2 = q.y; } else { a.x1 = q.x; a.y1 = q.y; }
        ed.touch();
        break;
      }
      case 'curvedCp': {
        const h = this.target as Extract<Hit, { kind: 'curved' }>;
        const c = ed.doc.curved.get(h.id);
        if (!c) return;
        ed.begin();
        const g = curvedGeometry(ed.doc, c);
        if (!g) return;
        const S = g.rawStart, E = g.rawEnd;
        let d = sub(E, S);
        let L = len(d);
        if (L < 1e-6) { d = { x: 1, y: 0 }; L = 0.6; }
        const u = mul(d, 1 / len(d));
        const v = perp(u);
        const rel = sub(p, S);
        const t = (rel.x * u.x + rel.y * u.y) / L, hh = (rel.x * v.x + rel.y * v.y) / Math.max(L, 0.9);
        if (h.part === 'c1') c.c1 = { t, h: hh };
        else c.c2 = { t, h: hh };
        ed.touch();
        break;
      }
      case 'shapeCorner': {
        const h = this.target as Extract<Hit, { kind: 'shape' }>;
        const s = ed.doc.shapes.get(h.id);
        if (!s) return;
        ed.begin();
        if (h.part === 'p2') { s.x2 = p.x; s.y2 = p.y; } else { s.x1 = p.x; s.y1 = p.y; }
        ed.touch();
        break;
      }
      case 'marquee': {
        this.rect!.x2 = p.x;
        this.rect!.y2 = p.y;
        this.liveSelect();
        break;
      }
      case 'lasso': {
        this.lasso.push(p);
        this.liveSelect();
        break;
      }
    }
  }

  /** Snapping while moving: merge targets, terminal bond angle/length, alignment guides, grid. */
  private snapMove(dx: number, dy: number, alt: boolean): { dx: number; dy: number } {
    const ed = this.ed;
    this.guides = {};
    if (alt) return { dx, dy };
    const moved = this.moved;
    // single atom: terminal atoms snap to unit length & 15° around their neighbour
    if (moved.size === 1 && this.ed.sel.objects.size === 0) {
      const id = [...moved][0];
      const o = this.moveOrigin.get(id)!;
      const target = { x: o.x + dx, y: o.y + dy };
      const near = atomNear(ed.doc, target, MERGE_TOL * 1.4, moved);
      if (near) return { dx: near.x - o.x, dy: near.y - o.y };
      const adj = adjacency(ed.doc);
      const nb = neighborsOf(ed.doc, id, adj);
      if (nb.length === 1) {
        const n = ed.doc.atoms.get(nb[0])!;
        const d = sub(target, n);
        const L = len(d);
        if (Math.abs(L - 1) < 0.35) {
          const q = add(n, fromAngle(snapAngle(angleOf(d), DEG15), 1));
          return { dx: q.x - o.x, dy: q.y - o.y };
        }
      }
    }
    // alignment guides against static atoms
    const tol = 6 / ed.view.scale;
    let bestX: { d: number; corr: number; x: number } | null = null;
    let bestY: { d: number; corr: number; y: number } | null = null;
    const statics: Pt[] = [];
    for (const a of ed.doc.atoms.values()) if (!moved.has(a.id)) statics.push(a);
    if (statics.length < 2000 && moved.size <= 60) {
      for (const id of moved) {
        const o = this.moveOrigin.get(id)!;
        const px = o.x + dx, py = o.y + dy;
        for (const s of statics) {
          const ddx = s.x - px, ddy = s.y - py;
          if (Math.abs(ddx) < tol && (!bestX || Math.abs(ddx) < bestX.d)) bestX = { d: Math.abs(ddx), corr: ddx, x: s.x };
          if (Math.abs(ddy) < tol && (!bestY || Math.abs(ddy) < bestY.d)) bestY = { d: Math.abs(ddy), corr: ddy, y: s.y };
        }
      }
    }
    if (bestX) { dx += bestX.corr; this.guides.x = bestX.x; }
    if (bestY) { dy += bestY.corr; this.guides.y = bestY.y; }
    if (ed.grid.snap && !bestX && !bestY) {
      const g = ed.grid.size;
      dx = Math.round(dx / g) * g;
      dy = Math.round(dy / g) * g;
    }
    return { dx, dy };
  }

  private liveSelect(): void {
    const ed = this.ed;
    const s = cloneSel(this.baseSel);
    const inside = (p: Pt) => {
      if (this.mode === 'lasso') return this.lasso.length > 2 && pointInPolygon(p, this.lasso);
      const r = this.rect!;
      return p.x >= Math.min(r.x1, r.x2) && p.x <= Math.max(r.x1, r.x2) && p.y >= Math.min(r.y1, r.y2) && p.y <= Math.max(r.y1, r.y2);
    };
    for (const a of ed.doc.atoms.values()) if (inside(a)) s.atoms.add(a.id);
    for (const b of ed.doc.bonds.values()) if (s.atoms.has(b.a) && s.atoms.has(b.b)) s.bonds.add(b.id);
    const scene = ed.getScene();
    for (const m of [ed.doc.arrows, ed.doc.texts, ed.doc.shapes, ed.doc.curved] as Map<number, unknown>[]) {
      for (const id of m.keys()) {
        const bx = ed.objectBox(id, scene);
        if (bx && inside({ x: bx.x1, y: bx.y1 }) && inside({ x: bx.x2, y: bx.y2 })) s.objects.add(id);
      }
    }
    ed.setSelection(s);
  }

  up(e: PEvent & { moved?: boolean }): void {
    const ed = this.ed;
    if (this.mode === 'move') {
      if (e.moved && ed.history.inProgress) {
        const n = fuseOverlaps(ed.doc, this.moved);
        ed.commit(n ? 'Move & join' : 'Move');
        if (n) ed.setStatus(`Joined ${n} atom${n > 1 ? 's' : ''}`);
      }
    } else if (this.mode === 'rotate' || this.mode === 'arrowEnd' || this.mode === 'curvedCp' || this.mode === 'shapeCorner') {
      ed.commit(this.mode === 'rotate' ? 'Rotate' : 'Edit');
    } else if ((this.mode === 'marquee' || this.mode === 'lasso') && !e.moved) {
      if (!this.additive) ed.clearSelection();
    }
    this.mode = 'none';
    this.rect = null;
    this.lasso = [];
    this.guides = {};
    this.merges = [];
    ed.requestRender();
  }

  cancel(): void {
    this.mode = 'none';
    this.rect = null;
    this.lasso = [];
    this.guides = {};
    this.merges = [];
  }

  dblclick(e: PEvent): void {
    const ed = this.ed;
    const h = e.hit;
    if (!h) return;
    if (h.kind === 'atom') ed.editAtomLabel(h.id);
    else if (h.kind === 'bond') {
      const b = ed.doc.bonds.get(h.id)!;
      ed.selectFragments([b.a]);
    } else if (h.kind === 'text') editTextObject(ed, h.id);
    else if (h.kind === 'arrow') editArrowCaption(ed, h.id, { x: e.x, y: e.y });
  }

  overlay(ctx: CanvasRenderingContext2D): void {
    const ed = this.ed;
    ctx.save();
    if (this.mode === 'marquee' && this.rect) {
      const p1 = ed.toScreen({ x: this.rect.x1, y: this.rect.y1 }), p2 = ed.toScreen({ x: this.rect.x2, y: this.rect.y2 });
      ctx.fillStyle = ed.theme.dark ? 'rgba(77,171,247,0.12)' : 'rgba(25,113,194,0.08)';
      ctx.strokeStyle = ed.theme.accent;
      ctx.setLineDash([4, 3]);
      ctx.fillRect(p1.x, p1.y, p2.x - p1.x, p2.y - p1.y);
      ctx.strokeRect(p1.x, p1.y, p2.x - p1.x, p2.y - p1.y);
    }
    if (this.mode === 'lasso' && this.lasso.length > 1) {
      ctx.strokeStyle = ed.theme.accent;
      ctx.fillStyle = ed.theme.dark ? 'rgba(77,171,247,0.12)' : 'rgba(25,113,194,0.08)';
      ctx.setLineDash([4, 3]);
      ctx.beginPath();
      this.lasso.forEach((p, i) => {
        const q = ed.toScreen(p);
        if (i === 0) ctx.moveTo(q.x, q.y);
        else ctx.lineTo(q.x, q.y);
      });
      ctx.closePath();
      ctx.fill();
      ctx.stroke();
    }
    ctx.setLineDash([3, 3]);
    ctx.strokeStyle = '#e8590c';
    ctx.lineWidth = 1;
    if (this.guides.x !== undefined) {
      const x = ed.toScreen({ x: this.guides.x, y: 0 }).x;
      ctx.beginPath();
      ctx.moveTo(x, 0);
      ctx.lineTo(x, ed.canvas.clientHeight);
      ctx.stroke();
    }
    if (this.guides.y !== undefined) {
      const y = ed.toScreen({ x: 0, y: this.guides.y }).y;
      ctx.beginPath();
      ctx.moveTo(0, y);
      ctx.lineTo(ed.canvas.clientWidth, y);
      ctx.stroke();
    }
    ctx.setLineDash([]);
    for (const [, t] of this.merges) {
      const a = ed.doc.atoms.get(t);
      if (!a) continue;
      const p = ed.toScreen(a);
      ctx.strokeStyle = '#2f9e44';
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.arc(p.x, p.y, Math.max(8, 0.3 * ed.view.scale), 0, Math.PI * 2);
      ctx.stroke();
    }
    ctx.restore();
  }
}

function cloneSel(s: Selection): Selection {
  return { atoms: new Set(s.atoms), bonds: new Set(s.bonds), objects: new Set(s.objects) };
}

export function editTextObject(ed: Editor, id: number): void {
  const t = ed.doc.texts.get(id);
  if (!t) return;
  ed.hiddenObjects.add(id);
  ed.rev++;
  ed.requestRender();
  const size = ed.doc.style.fontSize * (t.size ?? 1);
  ed.openInlineEditor(
    { x: t.x, y: t.y - size * 0.85 },
    t.text,
    (text) => {
      ed.hiddenObjects.delete(id);
      ed.begin();
      if (!text.trim()) ed.doc.texts.delete(id);
      else ed.doc.texts.get(id)!.text = text;
      ed.commit('Edit text');
    },
    {
      multiline: true,
      width: 220,
      placeholder: 'Text (H_2O, ^+ superscript, **bold**) — Ctrl+Enter to finish',
      onCancel: () => {
        ed.hiddenObjects.delete(id);
        if (!ed.doc.texts.get(id)?.text) {
          ed.doc.texts.delete(id);
        }
        ed.touch();
      },
    },
  );
}

/** Inline editing of reagents (above) or conditions (below) of a reaction arrow. */
export function editArrowCaption(ed: Editor, id: number, at: Pt): void {
  const a = ed.doc.arrows.get(id);
  if (!a) return;
  const mid = { x: (a.x1 + a.x2) / 2, y: (a.y1 + a.y2) / 2 };
  const u = norm({ x: a.x2 - a.x1, y: a.y2 - a.y1 });
  let up = perp(u);
  if (up.y > 0) up = mul(up, -1);
  const above = (at.x - mid.x) * up.x + (at.y - mid.y) * up.y >= 0;
  const key = above ? 'above' : 'below';
  const pos = add(mid, mul(up, above ? 0.75 : -0.15));
  ed.openInlineEditor(
    { x: pos.x - 1.2, y: pos.y - 0.2 },
    (a[key] ?? '').replace(/\n/g, ' / '),
    (text) => {
      ed.begin();
      const v = text.split(' / ').join('\n').trim();
      if (v) a[key] = v;
      else delete a[key];
      ed.commit(above ? 'Edit reagents' : 'Edit conditions');
    },
    { width: 180, placeholder: above ? 'Reagents (above)' : 'Conditions (below)' },
  );
}

// ───────────────────────── Eraser ─────────────────────────

class EraseTool implements Tool {
  id: ToolId = 'erase';
  cursor = 'cell';
  private active = false;
  constructor(private ed: Editor) {}
  hint(): string {
    return 'Click or drag over atoms, bonds and objects to erase them';
  }
  private erase(h: Hit): void {
    if (!h || h.kind === 'rotate') return;
    const ed = this.ed;
    ed.begin();
    const s = emptySelection();
    if (h.kind === 'atom') s.atoms.add(h.id);
    else if (h.kind === 'bond') s.bonds.add(h.id);
    else s.objects.add(h.id);
    deleteSelection(ed.doc, s);
    ed.touch();
  }
  down(e: PEvent): void {
    this.active = true;
    this.erase(e.hit);
  }
  drag(e: PEvent): void {
    if (this.active) this.erase(e.hit);
  }
  up(): void {
    this.active = false;
    this.ed.commit('Erase');
  }
  cancel(): void {
    this.active = false;
  }
}

// ───────────────────────── Bond ─────────────────────────

class BondTool implements Tool {
  id: ToolId = 'bond';
  cursor = 'crosshair';
  private startAtom: number | undefined;
  private startPoint: Pt | null = null;
  private startBond: number | undefined;
  private preview: { from: Pt; to: Pt } | null = null;
  private hoverGhost: [Pt, Pt][] = [];
  constructor(private ed: Editor) {}
  hint(): string {
    return 'Click an atom to sprout a bond · drag to draw at 15° steps (Alt = free) · click a bond to change it · hover + 1/2/3, W, H, O, N…';
  }
  down(e: PEvent): void {
    this.startAtom = e.hit?.kind === 'atom' ? e.hit.id : undefined;
    this.startBond = e.hit?.kind === 'bond' ? e.hit.id : undefined;
    this.startPoint = this.startAtom !== undefined ? { ...this.ed.doc.atoms.get(this.startAtom)! } : snapToGrid(this.ed, { x: e.x, y: e.y });
    this.preview = null;
    this.hoverGhost = [];
  }
  drag(e: PEvent): void {
    if (!this.startPoint) return;
    if (this.startBond !== undefined && this.startAtom === undefined) {
      // dragging from a bond's middle: start a new bond from its nearer atom
      const b = this.ed.doc.bonds.get(this.startBond)!;
      const A = this.ed.doc.atoms.get(b.a)!, B = this.ed.doc.atoms.get(b.b)!;
      this.startAtom = dist(A, this.startPoint) < dist(B, this.startPoint) ? A.id : B.id;
      this.startPoint = { ...this.ed.doc.atoms.get(this.startAtom)! };
      this.startBond = undefined;
    }
    const end = snapBondEnd(this.ed, this.startPoint, { x: e.x, y: e.y }, this.startAtom, e.alt);
    this.preview = { from: this.startPoint, to: end.p };
  }
  up(e: PEvent & { moved?: boolean }): void {
    const ed = this.ed;
    const { order, style } = ed.settings.bond;
    if (!e.moved || !this.preview) {
      if (this.startBond !== undefined) {
        ed.mutate('Change bond', (d) => applyBondType(d.bonds.get(this.startBond!)!, order, style));
      } else if (this.startAtom !== undefined) {
        ed.mutate('Add bond', (d) => {
          const a = d.atoms.get(this.startAtom!)!;
          const dir = idealBondDirection(d, a.id, order);
          bondToPoint(d, a.id, add(a, dir), order, style);
        });
      } else if (this.startPoint) {
        ed.mutate('Add bond', (d) => {
          const a = addAtom(d, { el: 'C', x: this.startPoint!.x, y: this.startPoint!.y });
          bondToPoint(d, a.id, add(a, fromAngle(-Math.PI / 6)), order, style);
        });
      }
    } else {
      const { from, to } = this.preview;
      ed.mutate('Add bond', (d) => {
        let sid = this.startAtom;
        if (sid === undefined) {
          const ex = atomNear(d, from, MERGE_TOL);
          sid = ex ? ex.id : addAtom(d, { el: 'C', x: from.x, y: from.y }).id;
        }
        bondToPoint(d, sid, to, order, style);
      });
    }
    this.reset();
  }
  private reset(): void {
    this.startAtom = undefined;
    this.startBond = undefined;
    this.startPoint = null;
    this.preview = null;
  }
  cancel(): void {
    this.reset();
  }
  hover(e: PEvent): void {
    const ed = this.ed;
    this.hoverGhost = [];
    if (e.touch) return;
    if (e.hit?.kind === 'atom') {
      const a = ed.doc.atoms.get(e.hit.id)!;
      const dir = idealBondDirection(ed.doc, a.id, ed.settings.bond.order);
      this.hoverGhost = [[a, add(a, dir)]];
    } else if (!e.hit) {
      const p = snapToGrid(ed, { x: e.x, y: e.y });
      this.hoverGhost = [[p, add(p, fromAngle(-Math.PI / 6))]];
    }
  }
  overlay(ctx: CanvasRenderingContext2D): void {
    if (this.preview) drawGhost(this.ed, ctx, [[this.preview.from, this.preview.to]], [this.preview.to]);
    else if (this.hoverGhost.length && this.ed.mouseInside) drawGhost(this.ed, ctx, this.hoverGhost, [], true);
  }
}

// ───────────────────────── Chain ─────────────────────────

class ChainTool implements Tool {
  id: ToolId = 'chain';
  cursor = 'crosshair';
  private startAtom: number | undefined;
  private start: Pt | null = null;
  private pts: Pt[] = [];
  private firstUp = true;
  constructor(private ed: Editor) {}
  hint(): string {
    return 'Drag to draw a zig-zag carbon chain (the bond count is shown) · Alt = free angle';
  }
  down(e: PEvent): void {
    this.startAtom = e.hit?.kind === 'atom' ? e.hit.id : undefined;
    this.start = this.startAtom !== undefined ? { ...this.ed.doc.atoms.get(this.startAtom)! } : snapToGrid(this.ed, { x: e.x, y: e.y });
    this.pts = [];
  }
  drag(e: PEvent): void {
    if (!this.start) return;
    const d = sub({ x: e.x, y: e.y }, this.start);
    let ang = angleOf(d);
    if (!e.alt) ang = snapAngle(ang, Math.PI / 6);
    const proj = len(d);
    const n = Math.max(1, Math.round(proj / 0.866));
    // choose the initial zig so the chain continues an existing bond in zig-zag fashion
    let firstUp = true;
    if (this.startAtom !== undefined) {
      const ideal = idealBondDirection(this.ed.doc, this.startAtom);
      const up = add(this.start, fromAngle(ang - Math.PI / 6));
      const dn = add(this.start, fromAngle(ang + Math.PI / 6));
      const target = add(this.start, ideal);
      firstUp = dist(up, target) <= dist(dn, target);
    } else firstUp = Math.cos(ang) >= 0;
    this.firstUp = firstUp;
    this.pts = chainPoints(this.start, ang, n, firstUp);
    this.ed.setStatus(`Chain: ${n} bond${n > 1 ? 's' : ''}`);
  }
  up(e: PEvent & { moved?: boolean }): void {
    const ed = this.ed;
    if (!this.start) return;
    const start = this.start;
    const pts = e.moved ? this.pts : [];
    ed.mutate('Add chain', (d) => {
      let prev = this.startAtom ?? (atomNear(d, start, MERGE_TOL)?.id ?? addAtom(d, { el: 'C', x: start.x, y: start.y }).id);
      if (!pts.length) {
        const a = d.atoms.get(prev)!;
        bondToPoint(d, prev, add(a, idealBondDirection(d, prev)), 1, 'plain');
        return;
      }
      for (const p of pts) prev = bondToPoint(d, prev, p, 1, 'plain').atom.id;
    });
    this.start = null;
    this.pts = [];
  }
  cancel(): void {
    this.start = null;
    this.pts = [];
  }
  overlay(ctx: CanvasRenderingContext2D): void {
    if (!this.start || !this.pts.length) return;
    const segs: [Pt, Pt][] = [];
    let prev = this.start;
    for (const p of this.pts) {
      segs.push([prev, p]);
      prev = p;
    }
    drawGhost(this.ed, ctx, segs, this.pts);
    const last = this.ed.toScreen(this.pts[this.pts.length - 1]);
    ctx.save();
    ctx.font = 'bold 12px system-ui, sans-serif';
    ctx.fillStyle = this.ed.theme.accent;
    ctx.fillText(String(this.pts.length), last.x + 8, last.y - 8);
    ctx.restore();
  }
}

// ───────────────────────── Ring / Template ─────────────────────────

class RingTool implements Tool {
  id: ToolId;
  cursor = 'crosshair';
  private ghost: { before: ChemDoc; after: ChemDoc } | null = null;
  private dragCenter: Pt | null = null;
  private dragAngle = -Math.PI / 2;
  private downHit: Hit = null;
  private downPoint: Pt = { x: 0, y: 0 };
  constructor(private ed: Editor, private template: boolean) {
    this.id = template ? 'template' : 'ring';
  }
  hint(): string {
    return this.template
      ? 'Click to place the template · click an atom to attach it · drag to rotate'
      : 'Click empty space to place a ring · click a bond to fuse · click an atom to attach (Shift = spiro) · drag to rotate';
  }
  /** Applies the ring/template action to `d`. */
  private apply(d: ChemDoc, hit: Hit, p: Pt, shift: boolean, angle?: number): void {
    const ed = this.ed;
    if (this.template) {
      const t = ed.settings.template;
      if (!t) return;
      let mol = t.mol;
      if (angle !== undefined) {
        mol = mol.clone();
        const bb = mol.bbox();
        const c = { x: (bb.minX + bb.maxX) / 2, y: (bb.minY + bb.maxY) / 2 };
        const rot = angle + Math.PI / 2;
        for (const a of mol.atoms) {
          const dx = a.x - c.x, dy = a.y - c.y;
          a.x = c.x + dx * Math.cos(rot) - dy * Math.sin(rot);
          a.y = c.y + dx * Math.sin(rot) + dy * Math.cos(rot);
        }
      }
      if (hit?.kind === 'atom') attachTemplate(d, mol, hit.id);
      else insertTemplateAt(d, mol, p);
      return;
    }
    const r = ed.settings.ring;
    if ('chair' in r) {
      placeChair(d, hit?.kind === 'atom' ? add(d.atoms.get(hit.id)!, { x: 1.2, y: -0.4 }) : p);
      return;
    }
    if (hit?.kind === 'bond') fuseRingOnBond(d, hit.id, r.size, r.aromatic, p);
    else if (hit?.kind === 'atom') {
      const a = d.atoms.get(hit.id)!;
      const dir = angle !== undefined ? fromAngle(angle) : undefined;
      void a;
      attachRingToAtom(d, hit.id, r.size, r.aromatic, shift, dir);
    } else placeRing(d, p, r.size, r.aromatic, angle ?? (r.size % 2 === 0 ? -Math.PI / 2 : -Math.PI / 2));
  }
  private preview(hit: Hit, p: Pt, shift: boolean, angle?: number): void {
    const before = this.ed.doc;
    const after = cloneDoc(before);
    try {
      this.apply(after, hit, p, shift, angle);
      this.ghost = { before, after };
    } catch {
      this.ghost = null;
    }
  }
  hover(e: PEvent): void {
    if (e.touch) {
      this.ghost = null;
      return;
    }
    this.preview(e.hit, snapToGrid(this.ed, { x: e.x, y: e.y }), e.shift);
  }
  down(e: PEvent): void {
    this.downHit = e.hit;
    this.downPoint = snapToGrid(this.ed, { x: e.x, y: e.y });
    this.dragCenter = null;
  }
  drag(e: PEvent): void {
    const p = { x: e.x, y: e.y };
    if (this.downHit?.kind === 'bond') {
      this.preview(this.downHit, p, e.shift);
      return;
    }
    const origin = this.downHit?.kind === 'atom' ? this.ed.doc.atoms.get(this.downHit.id)! : this.downPoint;
    this.dragCenter = origin;
    let ang = angleOf(sub(p, origin));
    if (!e.alt) ang = snapAngle(ang, DEG15);
    this.dragAngle = ang;
    this.preview(this.downHit, this.downPoint, e.shift, ang);
  }
  up(e: PEvent & { moved?: boolean }): void {
    const ed = this.ed;
    const hit = this.downHit;
    const p = e.moved && hit?.kind === 'bond' ? { x: e.x, y: e.y } : this.downPoint;
    const angle = e.moved && hit?.kind !== 'bond' ? this.dragAngle : undefined;
    ed.mutate(this.template ? 'Insert template' : 'Add ring', (d) => this.apply(d, hit, p, e.shift, angle));
    this.ghost = null;
    this.downHit = null;
  }
  cancel(): void {
    this.ghost = null;
    this.downHit = null;
  }
  overlay(ctx: CanvasRenderingContext2D): void {
    if (this.ghost && (this.ed.mouseInside || this.ed.isDragging)) drawDocDiff(this.ed, ctx, this.ghost.before, this.ghost.after);
  }
}

// ───────────────────────── Atom ─────────────────────────

class AtomTool implements Tool {
  id: ToolId = 'atom';
  cursor = 'crosshair';
  private startAtom: number | undefined;
  private preview: { from: Pt; to: Pt } | null = null;
  constructor(private ed: Editor) {}
  hint(): string {
    return `Click to place or change atoms to ${this.ed.settings.atom.label ?? this.ed.settings.atom.el} · drag from an atom to add a bonded atom`;
  }
  private applyTo(d: ChemDoc, id: number): void {
    const s = this.ed.settings.atom;
    if (s.label) this.ed.applyAtomLabel(id, s.label);
    else setElement(d, id, s.el);
  }
  down(e: PEvent): void {
    this.startAtom = e.hit?.kind === 'atom' ? e.hit.id : undefined;
    this.preview = null;
  }
  drag(e: PEvent): void {
    if (this.startAtom === undefined) return;
    const a = this.ed.doc.atoms.get(this.startAtom)!;
    const end = snapBondEnd(this.ed, a, { x: e.x, y: e.y }, this.startAtom, e.alt);
    this.preview = { from: a, to: end.p };
  }
  up(e: PEvent & { moved?: boolean }): void {
    const ed = this.ed;
    if (e.moved && this.preview && this.startAtom !== undefined) {
      const to = this.preview.to;
      ed.mutate('Add atom', (d) => {
        const r = bondToPoint(d, this.startAtom!, to, 1, 'plain');
        this.applyTo(d, r.atom.id);
      });
    } else if (e.hit?.kind === 'atom') {
      ed.mutate('Change atom', (d) => this.applyTo(d, (e.hit as { id: number }).id));
    } else if (e.hit?.kind === 'bond') {
      // apply to the nearer end of the bond
      const b = ed.doc.bonds.get(e.hit.id)!;
      const A = ed.doc.atoms.get(b.a)!, B = ed.doc.atoms.get(b.b)!;
      const id = dist(A, e) < dist(B, e) ? A.id : B.id;
      ed.mutate('Change atom', (d) => this.applyTo(d, id));
    } else {
      const p = snapToGrid(ed, { x: e.x, y: e.y });
      ed.mutate('Add atom', (d) => {
        const a = addAtom(d, { el: 'C', x: p.x, y: p.y });
        this.applyTo(d, a.id);
      });
    }
    this.startAtom = undefined;
    this.preview = null;
  }
  cancel(): void {
    this.startAtom = undefined;
    this.preview = null;
  }
  overlay(ctx: CanvasRenderingContext2D): void {
    if (this.preview) drawGhost(this.ed, ctx, [[this.preview.from, this.preview.to]], [this.preview.to]);
  }
}

// ───────────────────────── simple atom-click tools ─────────────────────────

class AtomClickTool implements Tool {
  cursor = 'pointer';
  constructor(private ed: Editor, public id: ToolId, private label: string, private fn: (ed: Editor, atomId: number) => void, private hintText: string) {}
  hint(): string {
    return this.hintText;
  }
  up(e: PEvent & { moved?: boolean }): void {
    if (e.moved) return;
    let id: number | undefined;
    if (e.hit?.kind === 'atom') id = e.hit.id;
    if (id === undefined) return;
    this.ed.mutate(this.label, () => this.fn(this.ed, id!));
  }
}

// ───────────────────────── Text / Plus ─────────────────────────

class TextTool implements Tool {
  id: ToolId = 'text';
  cursor = 'text';
  constructor(private ed: Editor) {}
  hint(): string {
    return 'Click to add text (H_2O subscripts, ^+ superscripts, **bold**, *italic*) · click an atom to type a label (OMe, CO2H, NH2…)';
  }
  up(e: PEvent & { moved?: boolean }): void {
    if (e.moved) return;
    const ed = this.ed;
    const h = e.hit;
    if (h?.kind === 'atom') return ed.editAtomLabel(h.id);
    if (h?.kind === 'text') return editTextObject(ed, h.id);
    if (h?.kind === 'arrow') return editArrowCaption(ed, h.id, { x: e.x, y: e.y });
    const p = snapToGrid(ed, { x: e.x, y: e.y });
    ed.begin();
    const id = ed.doc.nextId++;
    ed.doc.texts.set(id, { id, type: 'text', x: p.x, y: p.y, text: '' });
    ed.commit('Add text');
    editTextObject(ed, id);
  }
}

class PlusTool implements Tool {
  id: ToolId = 'plus';
  cursor = 'copy';
  constructor(private ed: Editor) {}
  hint(): string {
    return 'Click to place a “+” between reactants or products';
  }
  up(e: PEvent & { moved?: boolean }): void {
    if (e.moved) return;
    const ed = this.ed;
    const size = 1.4;
    ed.mutate('Add plus', (d) => {
      const id = d.nextId++;
      d.texts.set(id, { id, type: 'text', x: e.x, y: e.y + d.style.fontSize * size * 0.36, text: '+', size, align: 'center' });
    });
  }
}

// ───────────────────────── Arrows ─────────────────────────

class ArrowTool implements Tool {
  id: ToolId = 'arrow';
  cursor = 'crosshair';
  private from: Pt | null = null;
  private to: Pt | null = null;
  constructor(private ed: Editor) {}
  hint(): string {
    return 'Drag to draw a reaction arrow (15° steps, Alt = free) · click to add a default arrow · double-click above/below an arrow for reagents/conditions';
  }
  down(e: PEvent): void {
    if (e.hit?.kind === 'arrow') {
      this.from = null;
      return;
    }
    this.from = snapToGrid(this.ed, { x: e.x, y: e.y });
    this.to = null;
  }
  drag(e: PEvent): void {
    if (!this.from) return;
    const d = sub({ x: e.x, y: e.y }, this.from);
    let ang = angleOf(d);
    if (!e.alt) ang = snapAngle(ang, DEG15);
    this.to = add(this.from, fromAngle(ang, Math.max(0.5, len(d))));
  }
  up(e: PEvent & { moved?: boolean }): void {
    const ed = this.ed;
    if (!this.from) {
      if (!e.moved && e.hit?.kind === 'arrow') {
        // click on an existing arrow cycles its kind to the current one
        const id = e.hit.id;
        ed.mutate('Change arrow', (d) => (d.arrows.get(id)!.kind = ed.settings.arrow));
      }
      return;
    }
    const from = this.from;
    const to = e.moved && this.to ? this.to : add(from, { x: 2.6, y: 0 });
    ed.mutate('Add arrow', (d) => {
      const id = d.nextId++;
      d.arrows.set(id, { id, type: 'arrow', kind: ed.settings.arrow, x1: from.x, y1: from.y, x2: to.x, y2: to.y });
    });
    this.from = null;
    this.to = null;
  }
  cancel(): void {
    this.from = null;
    this.to = null;
  }
  overlay(ctx: CanvasRenderingContext2D): void {
    if (this.from && this.to) drawGhost(this.ed, ctx, [[this.from, this.to]], [this.to]);
  }
}

class CurvedTool implements Tool {
  id: ToolId = 'curved';
  cursor = 'crosshair';
  private from: Anchor | null = null;
  private fromPt: Pt | null = null;
  private to: Anchor | null = null;
  private toPt: Pt | null = null;
  private flip = false;
  constructor(private ed: Editor) {}
  hint(): string {
    return `Drag from a lone pair (atom) or bond to the electron destination (atom, bond, or between two atoms) · Shift flips the curve · ${this.ed.settings.curved === 2 ? 'electron pair' : 'single electron (fishhook)'}`;
  }
  down(e: PEvent): void {
    const an = anchorOfHit(e.hit);
    if (e.hit?.kind === 'curved') {
      // click on an existing curved arrow → handled in up (flip)
      this.from = null;
      return;
    }
    this.from = an ?? { type: 'point', x: e.x, y: e.y };
    this.fromPt = { x: e.x, y: e.y };
    this.to = null;
    this.toPt = null;
  }
  drag(e: PEvent): void {
    if (!this.from) return;
    this.flip = e.shift;
    const p = { x: e.x, y: e.y };
    this.toPt = p;
    this.to = this.targetAnchor(e);
  }
  private targetAnchor(e: PEvent): Anchor {
    const ed = this.ed;
    const h = ed.hitTest({ x: e.x, y: e.y }, e.touch);
    if (h?.kind === 'atom' && !(this.from?.type === 'atom' && this.from.id === h.id)) return { type: 'atom', id: h.id };
    if (h?.kind === 'bond' && !(this.from?.type === 'bond' && this.from.id === h.id)) return { type: 'bond', id: h.id };
    // midpoint between two non-bonded atoms → "between" (new bond target)
    let best: { a: number; b: number; d: number } | null = null;
    const atoms = [...ed.doc.atoms.values()];
    if (atoms.length < 400) {
      for (let i = 0; i < atoms.length; i++) {
        for (let j = i + 1; j < atoms.length; j++) {
          const A = atoms[i], B = atoms[j];
          const L = dist(A, B);
          if (L > 3.2 || L < 0.6) continue;
          if (bondBetween(ed.doc, A.id, B.id)) continue;
          const m = lerp(A, B, 0.5);
          const dd = dist(m, e);
          if (dd < 0.28 && (!best || dd < best.d)) best = { a: A.id, b: B.id, d: dd };
        }
      }
    }
    if (best) return { type: 'between', a: best.a, b: best.b };
    return { type: 'point', x: e.x, y: e.y };
  }
  /** Bulge to the side away from the molecule's bulk for a clean look. */
  private defaultBend(from: Anchor, to: Anchor): number {
    const ed = this.ed;
    const fake: CurvedArrowObj = { id: -1, type: 'curved', electrons: 2, from, to, c1: { t: 0.3, h: -0.5 }, c2: { t: 0.7, h: -0.5 } };
    const g = curvedGeometry(ed.doc, fake);
    if (!g) return -0.5;
    const mid = lerp(g.rawStart, g.rawEnd, 0.5);
    const u = norm(sub(g.rawEnd, g.rawStart));
    const v = perp(u);
    let s = 0;
    for (const a of ed.doc.atoms.values()) {
      const d = dist(a, mid);
      if (d > 2.5) continue;
      s += Math.sign((a.x - mid.x) * v.x + (a.y - mid.y) * v.y) / (0.5 + d);
    }
    return s > 0 ? -0.5 : 0.5;
  }
  up(e: PEvent & { moved?: boolean }): void {
    const ed = this.ed;
    if (!this.from) {
      if (!e.moved && e.hit?.kind === 'curved') {
        const id = e.hit.id;
        ed.mutate('Flip curved arrow', (d) => {
          const c = d.curved.get(id)!;
          c.c1.h = -c.c1.h;
          c.c2.h = -c.c2.h;
        });
      }
      return;
    }
    if (!e.moved || !this.to) {
      this.from = null;
      return;
    }
    const from = this.from, to = this.to;
    let h = this.defaultBend(from, to);
    if (this.flip) h = -h;
    ed.mutate('Add electron-pushing arrow', (d) => {
      const id = d.nextId++;
      d.curved.set(id, { id, type: 'curved', electrons: ed.settings.curved, from, to, c1: { t: 0.25, h }, c2: { t: 0.75, h } });
    });
    this.from = null;
    this.to = null;
  }
  cancel(): void {
    this.from = null;
    this.to = null;
  }
  overlay(ctx: CanvasRenderingContext2D): void {
    if (!this.from || !this.to) return;
    const ed = this.ed;
    let h = this.defaultBend(this.from, this.to);
    if (this.flip) h = -h;
    const fake: CurvedArrowObj = { id: -1, type: 'curved', electrons: ed.settings.curved, from: this.from, to: this.to, c1: { t: 0.25, h }, c2: { t: 0.75, h } };
    const g = curvedGeometry(ed.doc, fake, ed.getScene().labelBoxes);
    if (!g) return;
    ctx.save();
    ctx.strokeStyle = ghostColor(ed);
    ctx.lineWidth = 2;
    const p0 = ed.toScreen(g.p0), c1 = ed.toScreen(g.c1), c2 = ed.toScreen(g.c2), p3 = ed.toScreen(g.p3);
    ctx.beginPath();
    ctx.moveTo(p0.x, p0.y);
    ctx.bezierCurveTo(c1.x, c1.y, c2.x, c2.y, p3.x, p3.y);
    ctx.stroke();
    // target highlight
    if (this.to.type !== 'point') {
      ctx.fillStyle = 'rgba(214,51,108,0.25)';
      ctx.beginPath();
      ctx.arc(p3.x, p3.y, 9, 0, Math.PI * 2);
      ctx.fill();
    }
    ctx.restore();
  }
}

// ───────────────────────── Shapes ─────────────────────────

class ShapeTool implements Tool {
  id: ToolId = 'shape';
  cursor = 'crosshair';
  private from: Pt | null = null;
  private to: Pt | null = null;
  constructor(private ed: Editor) {}
  hint(): string {
    return 'Drag to draw brackets, transition-state brackets, boxes, ellipses or orbitals';
  }
  down(e: PEvent): void {
    this.from = { x: e.x, y: e.y };
    this.to = null;
  }
  drag(e: PEvent): void {
    if (!this.from) return;
    let to = { x: e.x, y: e.y };
    if (e.shift) {
      const s = Math.max(Math.abs(to.x - this.from.x), Math.abs(to.y - this.from.y));
      to = { x: this.from.x + Math.sign(to.x - this.from.x) * s, y: this.from.y + Math.sign(to.y - this.from.y) * s };
    }
    this.to = to;
  }
  up(e: PEvent & { moved?: boolean }): void {
    const ed = this.ed;
    if (!this.from) return;
    const kind = ed.settings.shape;
    let p1 = this.from, p2 = this.to ?? this.from;
    if (!e.moved || !this.to) {
      const hw = kind === 'orbitalP' ? 0.35 : kind === 'orbitalS' ? 0.4 : 1.2;
      const hh = kind === 'orbitalP' ? 0.8 : kind === 'orbitalS' ? 0.4 : 1;
      p1 = { x: this.from.x - hw, y: this.from.y - hh };
      p2 = { x: this.from.x + hw, y: this.from.y + hh };
    }
    ed.mutate('Add shape', (d) => {
      const id = d.nextId++;
      d.shapes.set(id, { id, type: 'shape', kind, x1: p1.x, y1: p1.y, x2: p2.x, y2: p2.y });
    });
    this.from = null;
    this.to = null;
  }
  cancel(): void {
    this.from = null;
    this.to = null;
  }
  overlay(ctx: CanvasRenderingContext2D): void {
    if (!this.from || !this.to) return;
    const a = this.from, b = this.to;
    drawGhost(this.ed, ctx, [[a, { x: b.x, y: a.y }], [{ x: b.x, y: a.y }, b], [b, { x: a.x, y: b.y }], [{ x: a.x, y: b.y }, a]], [], true);
  }
}

class PanTool implements Tool {
  id: ToolId = 'pan';
  cursor = 'grab';
  hint(): string {
    return 'Drag to pan · pinch or Ctrl+wheel to zoom';
  }
}

// ───────────────────────── factory ─────────────────────────

export function createTools(ed: Editor): Record<ToolId, Tool> {
  return {
    select: new SelectTool(ed, false),
    lasso: new SelectTool(ed, true),
    erase: new EraseTool(ed),
    bond: new BondTool(ed),
    chain: new ChainTool(ed),
    ring: new RingTool(ed, false),
    template: new RingTool(ed, true),
    atom: new AtomTool(ed),
    charge: new AtomClickTool(ed, 'charge', 'Change charge', (e, id) => {
      const a = e.doc.atoms.get(id)!;
      a.charge = (a.charge ?? 0) + e.settings.charge;
    }, 'Click an atom to add a positive or negative formal charge'),
    radical: new AtomClickTool(ed, 'radical', 'Radical', (e, id) => {
      const a = e.doc.atoms.get(id)!;
      a.radical = ((a.radical ?? 0) + 1) % 3;
      if (!a.radical) delete a.radical;
    }, 'Click an atom to cycle: radical (•) → diradical/carbene (••) → none'),
    lonepair: new AtomClickTool(ed, 'lonepair', 'Lone pairs', (e, id) => {
      const a = e.doc.atoms.get(id)!;
      a.lonePairs = !a.lonePairs;
      if (!a.lonePairs) delete a.lonePairs;
    }, 'Click an atom to show/hide its lone pairs (computed from the electron count)'),
    text: new TextTool(ed),
    arrow: new ArrowTool(ed),
    curved: new CurvedTool(ed),
    shape: new ShapeTool(ed),
    plus: new PlusTool(ed),
    pan: new PanTool(),
  };
}

export { addBond };
