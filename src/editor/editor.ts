// The interactive canvas editor: view transform, input dispatch, selection, history, rendering.
import { ChemDoc, CurvedArrowObj } from '../doc/types';
import { createDoc, cloneDoc, adjacency, docToMol, fragmentOf, insertMol, serializeDoc, deserializeDoc, docBounds } from '../doc/document';
import { buildScene, Scene, curvedGeometry } from '../render/scene';
import { drawPrims, ViewTransform } from '../render/draw';
import { Pt, distToSegment, distToBezier, dist, Box } from '../render/geom';
import { History } from './history';
import { Selection, emptySelection, selectedAtomIds, selectionSize, deleteSelection, setAtomLabel, moveSelection, selectionCenter } from './ops';
import { Hit, PEvent, Tool, ToolId, ToolSettings } from './types';
import { createTools } from './tools';
import { implicitH } from '../chem/valence';
import { Mol } from '../chem/mol';

export interface Theme {
  ink: string;
  bg: string;
  dark: boolean;
  accent: string;
  hover: string;
  selection: string;
  grid: string;
}

export const LIGHT_THEME: Theme = {
  ink: '#111111', bg: '#ffffff', dark: false, accent: '#1971c2', hover: 'rgba(25,113,194,0.22)', selection: 'rgba(25,113,194,0.28)', grid: '#e9ecef',
};
export const DARK_THEME: Theme = {
  ink: '#e6e6e6', bg: '#1b1d21', dark: true, accent: '#4dabf7', hover: 'rgba(77,171,247,0.25)', selection: 'rgba(77,171,247,0.3)', grid: '#2a2d33',
};

export interface Annotations {
  stereo?: Map<number, string>;
  bondLabels?: Map<number, string>;
  notes?: Map<number, string>;
  halo?: Map<number, string>;
}

type Listener = () => void;

export class Editor {
  readonly container: HTMLElement;
  readonly canvas: HTMLCanvasElement;
  readonly ctx: CanvasRenderingContext2D;
  doc: ChemDoc;
  rev = 0;
  history = new History();
  view: ViewTransform = { scale: 42, ox: 0, oy: 0 };
  sel: Selection = emptySelection();
  hover: Hit = null;
  toolId: ToolId = 'bond';
  tools: Record<ToolId, Tool>;
  settings: ToolSettings = {
    bond: { order: 1, style: 'plain' },
    ring: { size: 6, aromatic: true },
    template: null,
    atom: { el: 'N' },
    charge: 1,
    arrow: 'reaction',
    curved: 2,
    shape: 'bracket',
  };
  theme: Theme = LIGHT_THEME;
  annotations: Annotations = {};
  grid = { show: false, snap: false, size: 0.5 };
  wheelZooms = true;
  /** Last pointer position in model coordinates. */
  mouse: Pt = { x: 0, y: 0 };
  mouseInside = false;
  /** Hidden labels/objects while editing inline. */
  hiddenLabels = new Set<number>();
  hiddenObjects = new Set<number>();
  readonly isTouchDevice: boolean;

  private scene: Scene | null = null;
  private sceneKey = '';
  private renderQueued = false;
  private listeners: Record<string, Listener[]> = {};
  private pointers = new Map<number, { sx: number; sy: number }>();
  private pinch: { d: number; cx: number; cy: number; scale: number; ox: number; oy: number } | null = null;
  private down: { e: PEvent; moved: boolean; pointerId: number; pan: boolean; time: number } | null = null;
  private spaceHeld = false;
  private longPressTimer: number | null = null;
  private lastTap: { t: number; sx: number; sy: number } | null = null;
  private inlineEditor: HTMLElement | null = null;
  private dpr = 1;
  status = '';

  constructor(container: HTMLElement) {
    this.container = container;
    this.canvas = document.createElement('canvas');
    this.canvas.className = 'cw-canvas';
    this.canvas.tabIndex = 0;
    this.canvas.setAttribute('aria-label', 'Chemical structure drawing canvas');
    container.appendChild(this.canvas);
    this.ctx = this.canvas.getContext('2d')!;
    this.doc = createDoc();
    this.isTouchDevice = matchMedia('(pointer: coarse)').matches;
    this.tools = createTools(this);
    this.bindEvents();
    const ro = new ResizeObserver(() => this.resize());
    ro.observe(container);
    this.resize();
    this.view.ox = container.clientWidth / 2;
    this.view.oy = container.clientHeight / 2;
  }

  // ───────────── events ─────────────

  on(evt: 'change' | 'commit' | 'selection' | 'tool' | 'view' | 'status' | 'hover', fn: Listener): void {
    (this.listeners[evt] ??= []).push(fn);
  }
  contextMenuHandler: ((sx: number, sy: number, hit: Hit) => void) | null = null;
  emit(evt: string): void {
    for (const fn of this.listeners[evt] ?? []) fn();
  }

  setStatus(msg: string): void {
    this.status = msg;
    this.emit('status');
  }

  // ───────────── document & history ─────────────

  setDoc(doc: ChemDoc, resetHistory = true): void {
    this.doc = doc;
    this.sel = emptySelection();
    if (resetHistory) this.history.clear();
    this.touch();
    this.emit('selection');
  }

  /** Marks the document as modified (re-render + listeners). */
  touch(): void {
    this.rev++;
    this.requestRender();
    this.emit('change');
  }

  begin(): void {
    this.history.begin(this.doc);
  }

  commit(label: string): void {
    if (!this.history.inProgress) return;
    this.history.commit(label);
    this.doc.meta.modified = new Date().toISOString();
    this.cleanSelection();
    this.touch();
    this.emit('commit');
  }

  /** Abort an in-progress change and restore the document. */
  cancelChange(): void {
    const prev = this.history.cancel();
    if (prev) {
      this.doc = prev;
      this.cleanSelection();
      this.touch();
    }
  }

  /** Convenience: run a mutation as one undoable step. */
  mutate(label: string, fn: (doc: ChemDoc) => void): void {
    this.begin();
    fn(this.doc);
    this.commit(label);
  }

  undo(): void {
    if (this.history.inProgress) this.cancelChange();
    const d = this.history.undo(this.doc);
    if (d) {
      this.doc = d;
      this.cleanSelection();
      this.touch();
      this.emit('commit');
      this.setStatus('Undo');
    }
  }

  redo(): void {
    const d = this.history.redo(this.doc);
    if (d) {
      this.doc = d;
      this.cleanSelection();
      this.touch();
      this.emit('commit');
      this.setStatus('Redo');
    }
  }

  private cleanSelection(): void {
    const s = this.sel;
    for (const id of [...s.atoms]) if (!this.doc.atoms.has(id)) s.atoms.delete(id);
    for (const id of [...s.bonds]) if (!this.doc.bonds.has(id)) s.bonds.delete(id);
    for (const id of [...s.objects])
      if (!this.doc.arrows.has(id) && !this.doc.texts.has(id) && !this.doc.shapes.has(id) && !this.doc.curved.has(id)) s.objects.delete(id);
    this.emit('selection');
  }

  // ───────────── selection ─────────────

  setSelection(sel: Selection): void {
    this.sel = sel;
    this.requestRender();
    this.emit('selection');
  }

  clearSelection(): void {
    if (!selectionSize(this.sel)) return;
    this.setSelection(emptySelection());
  }

  selectAll(): void {
    const s = emptySelection();
    for (const id of this.doc.atoms.keys()) s.atoms.add(id);
    for (const id of this.doc.bonds.keys()) s.bonds.add(id);
    for (const m of [this.doc.arrows, this.doc.texts, this.doc.shapes, this.doc.curved]) for (const id of m.keys()) s.objects.add(id);
    this.setSelection(s);
  }

  /** Selects whole fragments containing any selected atom. */
  selectFragments(atomIds: number[]): void {
    const s = emptySelection();
    const adj = adjacency(this.doc);
    for (const id of atomIds) for (const a of fragmentOf(this.doc, id, adj)) s.atoms.add(a);
    for (const b of this.doc.bonds.values()) if (s.atoms.has(b.a) && s.atoms.has(b.b)) s.bonds.add(b.id);
    this.setSelection(s);
  }

  hasSelection(): boolean {
    return selectionSize(this.sel) > 0;
  }

  /** Atoms to analyse: selected atoms (whole fragments) or everything. */
  analysisAtomIds(): number[] {
    const ids = selectedAtomIds(this.doc, this.sel);
    if (!ids.size) return [...this.doc.atoms.keys()];
    return [...ids];
  }

  deleteSelected(): void {
    if (!this.hasSelection()) return;
    this.mutate('Delete', (d) => deleteSelection(d, this.sel));
    this.clearSelection();
  }

  // ───────────── tools ─────────────

  setTool(id: ToolId): void {
    this.tool.cancel?.();
    this.toolId = id;
    this.canvas.style.cursor = this.tool.cursor ?? 'crosshair';
    this.requestRender();
    this.emit('tool');
    this.setStatus(this.tool.hint?.() ?? '');
  }

  get tool(): Tool {
    return this.tools[this.toolId];
  }

  // ───────────── view ─────────────

  toModel(sx: number, sy: number): Pt {
    return { x: (sx - this.view.ox) / this.view.scale, y: (sy - this.view.oy) / this.view.scale };
  }

  toScreen(p: Pt): Pt {
    return { x: p.x * this.view.scale + this.view.ox, y: p.y * this.view.scale + this.view.oy };
  }

  zoomAt(factor: number, sx: number, sy: number): void {
    const ns = Math.max(6, Math.min(400, this.view.scale * factor));
    const f = ns / this.view.scale;
    this.view.ox = sx - (sx - this.view.ox) * f;
    this.view.oy = sy - (sy - this.view.oy) * f;
    this.view.scale = ns;
    this.requestRender();
    this.emit('view');
  }

  zoomBy(factor: number): void {
    this.zoomAt(factor, this.canvas.clientWidth / 2, this.canvas.clientHeight / 2);
  }

  setZoomPercent(p: number): void {
    this.zoomBy((42 * p) / 100 / this.view.scale);
  }

  get zoomPercent(): number {
    return Math.round((this.view.scale / 42) * 100);
  }

  /** Fit the drawing (or a box) into view. */
  fitToContent(box?: Box | null, maxScale = 60): void {
    const b = box ?? this.getScene().bounds;
    const W = this.canvas.clientWidth, H = this.canvas.clientHeight;
    if (!b || !W || !H) {
      this.view = { scale: 42, ox: W / 2, oy: H / 2 };
      this.requestRender();
      this.emit('view');
      return;
    }
    const pad = 1.2;
    const bw = b.x2 - b.x1 + 2 * pad, bh = b.y2 - b.y1 + 2 * pad;
    const s = Math.min(maxScale, W / bw, H / bh);
    this.view.scale = Math.max(6, s);
    this.view.ox = W / 2 - ((b.x1 + b.x2) / 2) * this.view.scale;
    this.view.oy = H / 2 - ((b.y1 + b.y2) / 2) * this.view.scale;
    this.requestRender();
    this.emit('view');
  }

  /** Model point at the visible centre of the canvas. */
  viewCenter(): Pt {
    return this.toModel(this.canvas.clientWidth / 2, this.canvas.clientHeight / 2);
  }

  /** A free spot to the right of existing content for inserting new structures. */
  freeSpot(width = 3): Pt {
    const b = docBounds(this.doc);
    if (!b) return this.viewCenter();
    const vc = this.viewCenter();
    // if the view centre is empty enough use it
    let crowded = false;
    for (const a of this.doc.atoms.values()) if (Math.abs(a.x - vc.x) < width && Math.abs(a.y - vc.y) < 2) { crowded = true; break; }
    if (!crowded) return vc;
    return { x: b.maxX + width / 2 + 1.5, y: (b.minY + b.maxY) / 2 };
  }

  private resize(): void {
    const W = this.container.clientWidth, H = this.container.clientHeight;
    this.dpr = window.devicePixelRatio || 1;
    this.canvas.width = Math.max(1, Math.round(W * this.dpr));
    this.canvas.height = Math.max(1, Math.round(H * this.dpr));
    this.canvas.style.width = W + 'px';
    this.canvas.style.height = H + 'px';
    this.requestRender();
  }

  // ───────────── scene / rendering ─────────────

  getScene(): Scene {
    const key = `${this.rev}|${this.theme.ink}|${this.hiddenLabels.size}|${this.hiddenObjects.size}|${this.annotationsRev}`;
    if (!this.scene || key !== this.sceneKey) {
      this.scene = buildScene(this.doc, {
        ink: this.theme.ink,
        dark: this.theme.dark,
        stereoLabels: this.annotations.stereo,
        bondLabels: this.annotations.bondLabels,
        atomNotes: this.annotations.notes,
        atomHalo: this.annotations.halo,
        hiddenLabels: this.hiddenLabels,
        hiddenObjects: this.hiddenObjects,
      });
      this.sceneKey = key;
    }
    return this.scene;
  }

  private annotationsRev = 0;
  setAnnotations(a: Annotations): void {
    this.annotations = a;
    this.annotationsRev++;
    this.requestRender();
  }

  requestRender(): void {
    if (this.renderQueued) return;
    this.renderQueued = true;
    requestAnimationFrame(() => {
      this.renderQueued = false;
      this.render();
    });
  }

  render(): void {
    const ctx = this.ctx;
    const W = this.canvas.clientWidth, H = this.canvas.clientHeight;
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    ctx.fillStyle = this.theme.bg;
    ctx.fillRect(0, 0, W, H);
    if (this.grid.show) this.drawGrid(ctx, W, H);
    const scene = this.getScene();
    this.drawSelectionUnder(ctx, scene);
    this.drawHover(ctx, scene);
    drawPrims(ctx, scene.prims, this.view);
    this.tool.overlay?.(ctx);
    this.drawSelectionBox(ctx);
    if (this.emptyHint && !this.doc.atoms.size && !this.doc.arrows.size && !this.doc.texts.size && !this.doc.shapes.size) this.drawEmptyHint(ctx, W, H);
  }

  /** Lines shown in the middle of an empty canvas (set by the app). */
  emptyHint: string[] | null = null;

  private drawEmptyHint(ctx: CanvasRenderingContext2D, W: number, H: number): void {
    const lines = this.emptyHint!;
    ctx.save();
    ctx.textAlign = 'center';
    ctx.fillStyle = this.theme.dark ? 'rgba(230,230,230,0.38)' : 'rgba(30,35,40,0.36)';
    const size = W < 500 ? 13 : 15;
    lines.forEach((l, i) => {
      ctx.font = `${i === 0 ? '600 ' : ''}${i === 0 ? size + 3 : size}px system-ui, -apple-system, sans-serif`;
      ctx.fillText(l, W / 2, H / 2 - ((lines.length - 1) * 26) / 2 + i * 26);
    });
    ctx.restore();
  }

  private drawGrid(ctx: CanvasRenderingContext2D, W: number, H: number): void {
    const step = this.grid.size * this.view.scale;
    if (step < 6) return;
    ctx.fillStyle = this.theme.grid;
    const x0 = ((this.view.ox % step) + step) % step;
    const y0 = ((this.view.oy % step) + step) % step;
    for (let x = x0; x < W; x += step) for (let y = y0; y < H; y += step) ctx.fillRect(x - 1, y - 1, 2, 2);
  }

  private drawSelectionUnder(ctx: CanvasRenderingContext2D, scene: Scene): void {
    if (!selectionSize(this.sel)) return;
    const s = this.view.scale;
    ctx.save();
    ctx.strokeStyle = this.theme.selection;
    ctx.fillStyle = this.theme.selection;
    ctx.lineCap = 'round';
    ctx.lineWidth = Math.max(6, 0.24 * s);
    for (const bid of this.sel.bonds) {
      const b = this.doc.bonds.get(bid);
      if (!b) continue;
      const A = this.toScreen(this.doc.atoms.get(b.a)!), B = this.toScreen(this.doc.atoms.get(b.b)!);
      ctx.beginPath();
      ctx.moveTo(A.x, A.y);
      ctx.lineTo(B.x, B.y);
      ctx.stroke();
    }
    for (const aid of this.sel.atoms) {
      const a = this.doc.atoms.get(aid);
      if (!a) continue;
      const p = this.toScreen(a);
      ctx.beginPath();
      ctx.arc(p.x, p.y, Math.max(5, 0.2 * s), 0, Math.PI * 2);
      ctx.fill();
    }
    for (const oid of this.sel.objects) {
      const bx = this.objectBox(oid, scene);
      if (!bx) continue;
      const p1 = this.toScreen({ x: bx.x1, y: bx.y1 }), p2 = this.toScreen({ x: bx.x2, y: bx.y2 });
      ctx.fillRect(p1.x, p1.y, p2.x - p1.x, p2.y - p1.y);
    }
    ctx.restore();
  }

  objectBox(id: number, scene = this.getScene()): Box | null {
    const c = this.doc.curved.get(id);
    if (c) {
      const g = curvedGeometry(this.doc, c, scene.labelBoxes);
      if (!g) return null;
      const xs = [g.p0.x, g.c1.x, g.c2.x, g.p3.x], ys = [g.p0.y, g.c1.y, g.c2.y, g.p3.y];
      return { x1: Math.min(...xs), y1: Math.min(...ys), x2: Math.max(...xs), y2: Math.max(...ys) };
    }
    return scene.objectBoxes.get(id) ?? null;
  }

  private drawHover(ctx: CanvasRenderingContext2D, scene: Scene): void {
    const h = this.hover;
    if (!h || this.down?.moved) return;
    const s = this.view.scale;
    ctx.save();
    ctx.fillStyle = this.theme.hover;
    ctx.strokeStyle = this.theme.hover;
    if (h.kind === 'atom') {
      const a = this.doc.atoms.get(h.id);
      if (a) {
        const p = this.toScreen(a);
        const boxes = scene.labelBoxes.get(h.id);
        let r = Math.max(7, 0.26 * s);
        if (boxes) r = Math.max(r, ((boxes[0].x2 - boxes[0].x1) / 2) * s * 0.95);
        ctx.beginPath();
        ctx.arc(p.x, p.y, r, 0, Math.PI * 2);
        ctx.fill();
      }
    } else if (h.kind === 'bond') {
      const b = this.doc.bonds.get(h.id);
      if (b) {
        const A = this.toScreen(this.doc.atoms.get(b.a)!), B = this.toScreen(this.doc.atoms.get(b.b)!);
        ctx.lineWidth = Math.max(8, 0.28 * s);
        ctx.lineCap = 'round';
        ctx.beginPath();
        ctx.moveTo(A.x, A.y);
        ctx.lineTo(B.x, B.y);
        ctx.stroke();
      }
    } else if (h.kind === 'curved' || h.kind === 'arrow' || h.kind === 'text' || h.kind === 'shape') {
      const bx = this.objectBox(h.id, scene);
      if (bx) {
        const p1 = this.toScreen({ x: bx.x1, y: bx.y1 }), p2 = this.toScreen({ x: bx.x2, y: bx.y2 });
        ctx.fillRect(p1.x - 2, p1.y - 2, p2.x - p1.x + 4, p2.y - p1.y + 4);
      }
    }
    ctx.restore();
  }

  /** Screen-space bounding box of the selection (for the floating toolbar and rotate handle). */
  selectionScreenBox(): { x1: number; y1: number; x2: number; y2: number } | null {
    if (!selectionSize(this.sel)) return null;
    let x1 = Infinity, y1 = Infinity, x2 = -Infinity, y2 = -Infinity;
    const inc = (p: Pt) => {
      const q = this.toScreen(p);
      x1 = Math.min(x1, q.x); x2 = Math.max(x2, q.x);
      y1 = Math.min(y1, q.y); y2 = Math.max(y2, q.y);
    };
    for (const id of selectedAtomIds(this.doc, this.sel)) inc(this.doc.atoms.get(id)!);
    const scene = this.getScene();
    for (const id of this.sel.objects) {
      const b = this.objectBox(id, scene);
      if (b) {
        inc({ x: b.x1, y: b.y1 });
        inc({ x: b.x2, y: b.y2 });
      }
    }
    if (x1 === Infinity) return null;
    const pad = 10;
    return { x1: x1 - pad, y1: y1 - pad, x2: x2 + pad, y2: y2 + pad };
  }

  rotateHandle(): Pt | null {
    if (this.toolId !== 'select' && this.toolId !== 'lasso') return null;
    const atoms = selectedAtomIds(this.doc, this.sel);
    if (atoms.size + this.sel.objects.size < 2) return null;
    const b = this.selectionScreenBox();
    if (!b) return null;
    return { x: (b.x1 + b.x2) / 2, y: b.y1 - 22 };
  }

  private drawSelectionBox(ctx: CanvasRenderingContext2D): void {
    const b = this.selectionScreenBox();
    if (!b) return;
    ctx.save();
    ctx.strokeStyle = this.theme.accent;
    ctx.lineWidth = 1;
    ctx.setLineDash([4, 3]);
    ctx.strokeRect(b.x1, b.y1, b.x2 - b.x1, b.y2 - b.y1);
    ctx.setLineDash([]);
    const h = this.rotateHandle();
    if (h) {
      ctx.beginPath();
      ctx.moveTo((b.x1 + b.x2) / 2, b.y1);
      ctx.lineTo(h.x, h.y + 6);
      ctx.stroke();
      ctx.beginPath();
      ctx.arc(h.x, h.y, 6, 0, Math.PI * 2);
      ctx.fillStyle = this.theme.bg;
      ctx.fill();
      ctx.stroke();
      // curved arrow glyph
      ctx.beginPath();
      ctx.arc(h.x, h.y, 3, -Math.PI * 0.9, Math.PI * 0.5);
      ctx.stroke();
    }
    ctx.restore();
  }

  // ───────────── hit testing ─────────────

  hitTest(p: Pt, touch = false, opts: { atoms?: boolean; bonds?: boolean; objects?: boolean } = {}): Hit {
    const scene = this.getScene();
    const s = this.view.scale;
    const tol = (touch ? 18 : 8) / s;
    // rotate handle
    const rh = this.rotateHandle();
    if (rh) {
      const sp = this.toScreen(p);
      if (Math.hypot(sp.x - rh.x, sp.y - rh.y) < (touch ? 16 : 9)) return { kind: 'rotate' };
    }
    // control points of a selected curved arrow
    if (opts.objects !== false) {
      for (const id of this.sel.objects) {
        const c = this.doc.curved.get(id);
        if (!c) continue;
        const g = curvedGeometry(this.doc, c, scene.labelBoxes);
        if (!g) continue;
        if (dist(p, g.c1) < tol * 1.3) return { kind: 'curved', id, part: 'c1' };
        if (dist(p, g.c2) < tol * 1.3) return { kind: 'curved', id, part: 'c2' };
      }
      for (const id of this.sel.objects) {
        const a = this.doc.arrows.get(id);
        if (a) {
          if (dist(p, { x: a.x2, y: a.y2 }) < tol * 1.3) return { kind: 'arrow', id, part: 'end' };
          if (dist(p, { x: a.x1, y: a.y1 }) < tol * 1.3) return { kind: 'arrow', id, part: 'start' };
        }
        const sh = this.doc.shapes.get(id);
        if (sh) {
          if (dist(p, { x: sh.x2, y: sh.y2 }) < tol * 1.3) return { kind: 'shape', id, part: 'p2' };
          if (dist(p, { x: sh.x1, y: sh.y1 }) < tol * 1.3) return { kind: 'shape', id, part: 'p1' };
        }
      }
    }
    // atoms
    if (opts.atoms !== false) {
      let best: number | null = null;
      let bd = Math.max(0.24, tol * 1.2);
      for (const a of this.doc.atoms.values()) {
        const d = Math.hypot(a.x - p.x, a.y - p.y);
        let dd = d;
        const boxes = scene.labelBoxes.get(a.id);
        if (boxes) {
          for (const bx of boxes) if (p.x >= bx.x1 && p.x <= bx.x2 && p.y >= bx.y1 && p.y <= bx.y2) dd = Math.min(dd, 0.05);
        }
        if (dd < bd) {
          bd = dd;
          best = a.id;
        }
      }
      if (best !== null) return { kind: 'atom', id: best };
    }
    if (opts.bonds !== false) {
      let best: number | null = null;
      let bd = Math.max(0.14, tol);
      for (const b of this.doc.bonds.values()) {
        const A = this.doc.atoms.get(b.a), B = this.doc.atoms.get(b.b);
        if (!A || !B) continue;
        const d = distToSegment(p, A, B);
        if (d < bd) {
          bd = d;
          best = b.id;
        }
      }
      if (best !== null) return { kind: 'bond', id: best };
    }
    if (opts.objects !== false) {
      for (const c of this.doc.curved.values()) {
        const g = curvedGeometry(this.doc, c, scene.labelBoxes);
        if (g && distToBezier(p, g.p0, g.c1, g.c2, g.p3) < Math.max(0.12, tol)) return { kind: 'curved', id: c.id, part: 'body' };
      }
      for (const a of this.doc.arrows.values()) {
        if (distToSegment(p, { x: a.x1, y: a.y1 }, { x: a.x2, y: a.y2 }) < Math.max(0.15, tol)) return { kind: 'arrow', id: a.id, part: 'body' };
        const bx = scene.objectBoxes.get(a.id);
        if (bx && (a.above || a.below) && p.x >= bx.x1 && p.x <= bx.x2 && p.y >= bx.y1 && p.y <= bx.y2) return { kind: 'arrow', id: a.id, part: 'body' };
      }
      for (const t of this.doc.texts.values()) {
        const bx = scene.objectBoxes.get(t.id);
        if (bx && p.x >= bx.x1 - tol && p.x <= bx.x2 + tol && p.y >= bx.y1 - tol && p.y <= bx.y2 + tol) return { kind: 'text', id: t.id };
      }
      for (const sh of this.doc.shapes.values()) {
        const x1 = Math.min(sh.x1, sh.x2), x2 = Math.max(sh.x1, sh.x2), y1 = Math.min(sh.y1, sh.y2), y2 = Math.max(sh.y1, sh.y2);
        const inside = p.x >= x1 - tol && p.x <= x2 + tol && p.y >= y1 - tol && p.y <= y2 + tol;
        if (!inside) continue;
        const nearEdge = Math.min(Math.abs(p.x - x1), Math.abs(p.x - x2), Math.abs(p.y - y1), Math.abs(p.y - y2)) < Math.max(0.15, tol);
        if (nearEdge || sh.fill || sh.kind === 'orbitalP' || sh.kind === 'orbitalS') return { kind: 'shape', id: sh.id, part: 'body' };
      }
    }
    return null;
  }

  // ───────────── input ─────────────

  private makeEvent(ev: PointerEvent | MouseEvent): PEvent {
    const r = this.canvas.getBoundingClientRect();
    const sx = ev.clientX - r.left, sy = ev.clientY - r.top;
    const m = this.toModel(sx, sy);
    const touch = (ev as PointerEvent).pointerType === 'touch' || (ev as PointerEvent).pointerType === 'pen';
    return {
      x: m.x, y: m.y, sx, sy,
      shift: ev.shiftKey, alt: ev.altKey, mod: ev.ctrlKey || ev.metaKey,
      button: ev.button, touch,
      hit: this.hitTest(m, touch),
      start: this.down?.e ? { x: this.down.e.x, y: this.down.e.y } : undefined,
    };
  }

  private bindEvents(): void {
    const c = this.canvas;
    c.addEventListener('pointerdown', (ev) => this.onPointerDown(ev));
    c.addEventListener('pointermove', (ev) => this.onPointerMove(ev));
    c.addEventListener('pointerup', (ev) => this.onPointerUp(ev));
    c.addEventListener('pointercancel', (ev) => this.onPointerCancel(ev));
    c.addEventListener('pointerleave', () => {
      this.mouseInside = false;
      if (!this.down) {
        this.hover = null;
        this.requestRender();
      }
    });
    c.addEventListener('pointerenter', () => (this.mouseInside = true));
    c.addEventListener('dblclick', (ev) => {
      const e = this.makeEvent(ev);
      this.tool.dblclick?.(e);
    });
    c.addEventListener('contextmenu', (ev) => {
      ev.preventDefault();
      const e = this.makeEvent(ev);
      this.contextMenuHandler?.(e.sx, e.sy, e.hit);
    });
    c.addEventListener('wheel', (ev) => this.onWheel(ev), { passive: false });
    window.addEventListener('keydown', (ev) => {
      if (ev.key === ' ' && this.isCanvasContext(ev)) {
        if (!this.spaceHeld) {
          this.spaceHeld = true;
          this.canvas.style.cursor = 'grab';
        }
        ev.preventDefault();
      }
    });
    window.addEventListener('keyup', (ev) => {
      if (ev.key === ' ') {
        this.spaceHeld = false;
        this.canvas.style.cursor = this.tool.cursor ?? 'crosshair';
      }
    });
  }

  /** True when keyboard events should go to the canvas (not to an input field). */
  isCanvasContext(ev: KeyboardEvent): boolean {
    const t = ev.target as HTMLElement | null;
    if (!t) return true;
    if (t.isContentEditable) return false;
    const tag = t.tagName;
    if (tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT') return false;
    return !this.inlineEditor;
  }

  private onPointerDown(ev: PointerEvent): void {
    this.canvas.focus({ preventScroll: true });
    this.commitInlineEditor();
    this.pointers.set(ev.pointerId, { sx: ev.clientX, sy: ev.clientY });
    try {
      this.canvas.setPointerCapture(ev.pointerId);
    } catch {
      /* ignore */
    }
    if (this.pointers.size === 2) {
      // start pinch: abort any tool action
      this.clearLongPress();
      if (this.down && this.down.moved) {
        this.tool.cancel?.();
        this.cancelChange();
      }
      this.down = null;
      const [p1, p2] = [...this.pointers.values()];
      const r = this.canvas.getBoundingClientRect();
      this.pinch = {
        d: Math.hypot(p1.sx - p2.sx, p1.sy - p2.sy),
        cx: (p1.sx + p2.sx) / 2 - r.left,
        cy: (p1.sy + p2.sy) / 2 - r.top,
        scale: this.view.scale, ox: this.view.ox, oy: this.view.oy,
      };
      return;
    }
    if (this.pointers.size > 2) return;
    const e = this.makeEvent(ev);
    this.mouse = { x: e.x, y: e.y };
    if (ev.button === 2) return;
    const pan = ev.button === 1 || this.spaceHeld || this.toolId === 'pan';
    this.down = { e, moved: false, pointerId: ev.pointerId, pan, time: performance.now() };
    if (pan) {
      this.canvas.style.cursor = 'grabbing';
      return;
    }
    if (e.touch) {
      this.clearLongPress();
      this.longPressTimer = window.setTimeout(() => {
        if (this.down && !this.down.moved) {
          const d = this.down;
          this.down = null;
          this.tool.cancel?.();
          this.cancelChange();
          this.contextMenuHandler?.(d.e.sx, d.e.sy, d.e.hit);
        }
      }, 550);
    }
    this.tool.down?.(e);
    this.requestRender();
  }

  private clearLongPress(): void {
    if (this.longPressTimer !== null) {
      clearTimeout(this.longPressTimer);
      this.longPressTimer = null;
    }
  }

  private onPointerMove(ev: PointerEvent): void {
    if (this.pointers.has(ev.pointerId)) this.pointers.set(ev.pointerId, { sx: ev.clientX, sy: ev.clientY });
    if (this.pinch && this.pointers.size >= 2) {
      const [p1, p2] = [...this.pointers.values()];
      const r = this.canvas.getBoundingClientRect();
      const d = Math.hypot(p1.sx - p2.sx, p1.sy - p2.sy);
      const cx = (p1.sx + p2.sx) / 2 - r.left, cy = (p1.sy + p2.sy) / 2 - r.top;
      const f = Math.max(6, Math.min(400, (this.pinch.scale * d) / Math.max(1, this.pinch.d))) / this.pinch.scale;
      this.view.scale = this.pinch.scale * f;
      this.view.ox = cx - (this.pinch.cx - this.pinch.ox) * f;
      this.view.oy = cy - (this.pinch.cy - this.pinch.oy) * f;
      this.requestRender();
      this.emit('view');
      return;
    }
    const e = this.makeEvent(ev);
    this.mouse = { x: e.x, y: e.y };
    if (this.down && this.down.pointerId === ev.pointerId) {
      const d0 = this.down.e;
      const movedPx = Math.hypot(e.sx - d0.sx, e.sy - d0.sy);
      if (!this.down.moved && movedPx > (e.touch ? 8 : 3)) {
        this.down.moved = true;
        this.clearLongPress();
      }
      if (this.down.pan) {
        this.view.ox += ev.movementX || 0;
        this.view.oy += ev.movementY || 0;
        if (!ev.movementX && !ev.movementY) {
          // movementX unsupported (some touch browsers): use absolute delta
          this.view.ox = this.view.ox + (e.sx - (this.lastPan?.sx ?? e.sx));
          this.view.oy = this.view.oy + (e.sy - (this.lastPan?.sy ?? e.sy));
        }
        this.lastPan = { sx: e.sx, sy: e.sy };
        this.requestRender();
        this.emit('view');
        return;
      }
      if (this.down.moved) {
        e.start = { x: d0.x, y: d0.y };
        this.tool.drag?.(e);
        this.requestRender();
      }
      return;
    }
    // hover
    const prev = this.hover;
    this.hover = e.hit;
    this.tool.hover?.(e);
    if (JSON.stringify(prev) !== JSON.stringify(this.hover)) this.emit('hover');
    this.requestRender();
  }
  private lastPan: { sx: number; sy: number } | null = null;

  private onPointerUp(ev: PointerEvent): void {
    this.pointers.delete(ev.pointerId);
    this.clearLongPress();
    if (this.pinch) {
      if (this.pointers.size < 2) this.pinch = null;
      return;
    }
    const d = this.down;
    if (!d || d.pointerId !== ev.pointerId) return;
    this.down = null;
    this.lastPan = null;
    if (d.pan) {
      this.canvas.style.cursor = this.spaceHeld ? 'grab' : this.tool.cursor ?? 'crosshair';
      return;
    }
    const e = this.makeEvent(ev);
    e.start = { x: d.e.x, y: d.e.y };
    if (!d.moved) {
      // treat as click at the down position (avoids jitter)
      e.x = d.e.x;
      e.y = d.e.y;
      e.hit = d.e.hit;
    }
    (e as PEvent & { moved?: boolean }).moved = d.moved;
    this.tool.up?.(e);
    // double tap detection for touch
    if (e.touch && !d.moved) {
      const now = performance.now();
      if (this.lastTap && now - this.lastTap.t < 320 && Math.hypot(e.sx - this.lastTap.sx, e.sy - this.lastTap.sy) < 14) {
        this.tool.dblclick?.(e);
        this.lastTap = null;
      } else this.lastTap = { t: now, sx: e.sx, sy: e.sy };
    }
    this.hover = this.hitTest({ x: e.x, y: e.y }, e.touch);
    this.requestRender();
  }

  private onPointerCancel(ev: PointerEvent): void {
    this.pointers.delete(ev.pointerId);
    this.clearLongPress();
    if (this.pointers.size < 2) this.pinch = null;
    if (this.down && this.down.pointerId === ev.pointerId) {
      this.down = null;
      this.tool.cancel?.();
      this.cancelChange();
    }
  }

  private onWheel(ev: WheelEvent): void {
    ev.preventDefault();
    const r = this.canvas.getBoundingClientRect();
    const sx = ev.clientX - r.left, sy = ev.clientY - r.top;
    const zoom = ev.ctrlKey || ev.metaKey || (this.wheelZooms && !ev.shiftKey);
    if (zoom) {
      const k = ev.deltaMode === 1 ? 0.05 : 0.0022;
      this.zoomAt(Math.exp(-ev.deltaY * k), sx, sy);
    } else {
      this.view.ox -= ev.shiftKey && !ev.deltaX ? ev.deltaY : ev.deltaX;
      this.view.oy -= ev.shiftKey && !ev.deltaX ? 0 : ev.deltaY;
      this.requestRender();
      this.emit('view');
    }
  }

  get isDragging(): boolean {
    return !!this.down?.moved;
  }

  // ───────────── inline editing ─────────────

  /** Shows an inline text input over the canvas. `onDone(text)` is called on commit (not on cancel). */
  openInlineEditor(at: Pt, initial: string, onDone: (text: string) => void, opts: { multiline?: boolean; width?: number; placeholder?: string; onCancel?: () => void } = {}): void {
    this.commitInlineEditor();
    const sp = this.toScreen(at);
    const el = document.createElement(opts.multiline ? 'textarea' : 'input') as HTMLInputElement | HTMLTextAreaElement;
    el.className = 'cw-inline-editor';
    el.value = initial;
    el.placeholder = opts.placeholder ?? '';
    el.spellcheck = false;
    el.setAttribute('autocapitalize', 'off');
    el.setAttribute('autocomplete', 'off');
    const fontPx = Math.max(13, this.doc.style.fontSize * this.view.scale);
    el.style.fontSize = fontPx + 'px';
    el.style.left = sp.x + 'px';
    el.style.top = sp.y + 'px';
    if (opts.width) el.style.width = opts.width + 'px';
    let done = false;
    const finish = (commit: boolean) => {
      if (done) return;
      done = true;
      const v = el.value;
      el.remove();
      this.inlineEditor = null;
      this.canvas.focus({ preventScroll: true });
      if (commit) onDone(v);
      else opts.onCancel?.();
    };
    (el as any)._finish = finish;
    el.addEventListener('keydown', (ev) => {
      const ke = ev as KeyboardEvent;
      ke.stopPropagation();
      if (ke.key === 'Escape') finish(false);
      else if (ke.key === 'Enter' && (!opts.multiline || ke.ctrlKey || ke.metaKey || !ke.shiftKey && el.tagName === 'INPUT')) {
        ke.preventDefault();
        finish(true);
      }
    });
    el.addEventListener('blur', () => setTimeout(() => finish(true), 0));
    this.container.appendChild(el);
    this.inlineEditor = el;
    el.focus();
    el.select();
  }

  commitInlineEditor(): void {
    const el = this.inlineEditor as any;
    if (el && el._finish) el._finish(true);
  }

  get isEditingInline(): boolean {
    return !!this.inlineEditor;
  }

  /** Inline label editing for an atom (ChemDraw: double-click / Enter on hovered atom). */
  editAtomLabel(atomId: number, initial?: string): void {
    const a = this.doc.atoms.get(atomId);
    if (!a) return;
    const cur = initial ?? (a.abbrev ?? a.alias ?? (a.el === 'C' && !a.charge ? '' : labelTextOf(this, atomId)));
    this.hiddenLabels.add(atomId);
    this.rev++;
    this.requestRender();
    this.openInlineEditor(
      a,
      cur,
      (text) => {
        this.hiddenLabels.delete(atomId);
        if (!text.trim()) {
          this.touch();
          return;
        }
        this.begin();
        if (!this.applyAtomLabel(atomId, text)) {
          this.cancelChange();
          this.setStatus(`Could not interpret label “${text}”`);
          this.touch();
          return;
        }
        this.commit('Edit label');
      },
      {
        placeholder: 'N, OH, CO2H, OMe, Ph…',
        onCancel: () => {
          this.hiddenLabels.delete(atomId);
          this.touch();
        },
      },
    );
  }

  /** Applies typed label text; normalises H counts equal to the default. */
  applyAtomLabel(atomId: number, text: string): boolean {
    if (!setAtomLabel(this.doc, atomId, text)) return false;
    const a = this.doc.atoms.get(atomId)!;
    if (a.hCount !== undefined) {
      const save = a.hCount;
      delete a.hCount;
      const { mol, index } = docToMol(this.doc, fragmentOf(this.doc, atomId));
      if (implicitH(mol, index.get(atomId)!) !== save) a.hCount = save;
    }
    return true;
  }

  // ───────────── keyboard ─────────────

  /** Handles a keydown; returns true if consumed. Hover hotkeys follow ChemDraw/Ketcher conventions. */
  handleKey(ev: KeyboardEvent): boolean {
    return handleEditorKey(this, ev);
  }

  // ───────────── clipboard ─────────────

  private internalClipboard: { json: string; text: string } | null = null;

  /** Serialises the selection (or everything) to a partial document JSON string. */
  selectionToJSON(): string {
    const ids = selectedAtomIds(this.doc, this.sel);
    const d = createDoc(this.doc.style);
    const all = !this.hasSelection();
    for (const a of this.doc.atoms.values()) if (all || ids.has(a.id)) d.atoms.set(a.id, { ...a });
    for (const b of this.doc.bonds.values()) if (d.atoms.has(b.a) && d.atoms.has(b.b)) d.bonds.set(b.id, { ...b });
    for (const m of ['arrows', 'texts', 'shapes', 'curved'] as const) {
      for (const [id, o] of this.doc[m] as Map<number, any>) if (all || this.sel.objects.has(id)) (d[m] as Map<number, any>).set(id, JSON.parse(JSON.stringify(o)));
    }
    // curved arrows whose anchors are inside the selection
    for (const [id, c] of this.doc.curved) {
      const ok = (an: CurvedArrowObj['from']) => (an.type === 'atom' ? d.atoms.has(an.id) : an.type === 'bond' ? d.bonds.has(an.id) : an.type === 'between' ? d.atoms.has(an.a) && d.atoms.has(an.b) : true);
      if (ok(c.from) && ok(c.to) && (c.from.type !== 'point' || c.to.type !== 'point')) d.curved.set(id, JSON.parse(JSON.stringify(c)));
    }
    d.nextId = this.doc.nextId;
    return JSON.stringify(serializeDoc(d));
  }

  setInternalClipboard(json: string, text: string): void {
    this.internalClipboard = { json, text };
  }

  getInternalClipboard(): { json: string; text: string } | null {
    return this.internalClipboard;
  }

  /** Pastes a serialized ChemWrite fragment, offset to `at` (or slightly shifted). Returns new selection. */
  pasteDocJSON(json: string, at?: Pt): void {
    const frag = deserializeDoc(JSON.parse(json));
    const b = docBounds(frag);
    if (!b) return;
    const cx = (b.minX + b.maxX) / 2, cy = (b.minY + b.maxY) / 2;
    const target = at ?? { x: cx + 0.6, y: cy + 0.6 };
    const dx = target.x - cx, dy = target.y - cy;
    this.begin();
    const idMap = new Map<number, number>();
    const nid = () => this.doc.nextId++;
    const sel = emptySelection();
    for (const a of frag.atoms.values()) {
      const id = nid();
      idMap.set(a.id, id);
      this.doc.atoms.set(id, { ...a, id, x: a.x + dx, y: a.y + dy });
      sel.atoms.add(id);
    }
    for (const bd of frag.bonds.values()) {
      const id = nid();
      idMap.set(bd.id, id);
      this.doc.bonds.set(id, { ...bd, id, a: idMap.get(bd.a)!, b: idMap.get(bd.b)! });
      sel.bonds.add(id);
    }
    for (const m of ['arrows', 'texts', 'shapes'] as const) {
      for (const o of (frag[m] as Map<number, any>).values()) {
        const id = nid();
        const c = { ...o, id };
        if ('x1' in c) {
          c.x1 += dx; c.x2 += dx; c.y1 += dy; c.y2 += dy;
        } else {
          c.x += dx; c.y += dy;
        }
        (this.doc[m] as Map<number, any>).set(id, c);
        sel.objects.add(id);
      }
    }
    for (const c of frag.curved.values()) {
      const id = nid();
      const remap = (an: CurvedArrowObj['from']): CurvedArrowObj['from'] => {
        if (an.type === 'atom') return { type: 'atom', id: idMap.get(an.id)! };
        if (an.type === 'bond') return { type: 'bond', id: idMap.get(an.id)! };
        if (an.type === 'between') return { type: 'between', a: idMap.get(an.a)!, b: idMap.get(an.b)! };
        return { type: 'point', x: an.x + dx, y: an.y + dy };
      };
      this.doc.curved.set(id, { ...c, id, from: remap(c.from), to: remap(c.to), c1: { ...c.c1 }, c2: { ...c.c2 } });
      sel.objects.add(id);
    }
    this.commit('Paste');
    this.setSelection(sel);
  }

  /** Inserts a laid-out Mol at a point (centred) as one undo step and selects it. */
  insertMolecule(mol: Mol, at?: Pt, label = 'Insert structure'): number[] {
    if (!mol.atoms.length) return [];
    const bb = mol.bbox();
    const w = bb.maxX - bb.minX;
    const p = at ?? this.freeSpot(w);
    this.begin();
    const { atomIds, bondIds } = insertMol(this.doc, mol, p.x - (bb.minX + bb.maxX) / 2, p.y - (bb.minY + bb.maxY) / 2);
    this.commit(label);
    const sel = emptySelection();
    atomIds.forEach((id) => sel.atoms.add(id));
    bondIds.forEach((id) => sel.bonds.add(id));
    this.setSelection(sel);
    return atomIds;
  }

  nudgeSelection(dx: number, dy: number): void {
    if (!this.hasSelection()) return;
    this.mutate('Nudge', (d) => moveSelection(d, this.sel, dx, dy));
  }

  selectionCenter(): Pt | null {
    return selectionCenter(this.doc, this.sel);
  }

  snapshotDoc(): ChemDoc {
    return cloneDoc(this.doc);
  }
}

/** Text shown in the inline editor for an atom label. */
export function labelTextOf(ed: Editor, atomId: number): string {
  const a = ed.doc.atoms.get(atomId)!;
  if (a.abbrev) return a.abbrev;
  if (a.alias) return a.alias;
  let s = (a.isotope ? a.isotope : '') + a.el;
  const scene = ed.getScene();
  const h = scene.hCount.get(atomId) ?? 0;
  if (a.el !== 'C' && h) s += 'H' + (h > 1 ? h : '');
  if (a.charge) s += (Math.abs(a.charge) > 1 ? Math.abs(a.charge) : '') + (a.charge > 0 ? '+' : '-');
  return s;
}

// keyboard handling lives in its own module to keep this file manageable
import { handleEditorKey } from './keys';
