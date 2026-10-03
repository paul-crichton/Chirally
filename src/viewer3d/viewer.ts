// Interactive 3D molecule viewer drawn with the Canvas 2D API (no WebGL).
// Perspective projection, depth-sorted shaded spheres and cylinders (painter's algorithm), depth cueing,
// trackball rotation with inertia, zoom, pan, touch gestures, keyboard control, measurement mode and
// PNG export. Colours follow the CPK scheme of elements.ts; the background and foreground follow the
// CSS variables --viewer-bg / --viewer-fg set on the container (with light/dark fallbacks).
import { Mol } from '../chem/mol';
import { element } from '../chem/elements';

export type Style3D = 'ballstick' | 'spacefill' | 'sticks' | 'wireframe';

type RGB = [number, number, number];

/** van der Waals radii (Å) for space-filling display (Bondi / Alvarez); default 2.0. */
const VDW: Record<string, number> = {
  H: 1.1, He: 1.4, Li: 1.82, B: 1.92, C: 1.7, N: 1.55, O: 1.52, F: 1.47, Ne: 1.54, Na: 2.27, Mg: 1.73, Al: 1.84,
  Si: 2.1, P: 1.8, S: 1.8, Cl: 1.75, Ar: 1.88, K: 2.75, Ca: 2.31, Fe: 2.0, Cu: 1.4, Zn: 1.39, Ga: 1.87, Ge: 2.11,
  As: 1.85, Se: 1.9, Br: 1.85, Kr: 2.02, Pd: 1.63, Ag: 1.72, Sn: 2.17, Sb: 2.06, Te: 2.06, I: 1.98, Xe: 2.16,
  Pt: 1.75, Au: 1.66, Hg: 1.55, Pb: 2.02,
};
const vdw = (el: string) => VDW[el] ?? 2.0;

const ACCENT: RGB = [255, 152, 0];
const FOG_MAX = 0.55; // fraction of background colour mixed into the farthest atoms
const SPRITE = 128; // sphere sprite resolution (px)

const clamp = (v: number, lo: number, hi: number) => Math.max(lo, Math.min(hi, v));
const mix = (a: RGB, b: RGB, t: number): RGB => [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t];
const css = (c: RGB, alpha = 1) =>
  alpha >= 1 ? `rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0})` : `rgba(${c[0] | 0},${c[1] | 0},${c[2] | 0},${alpha})`;
const WHITE: RGB = [255, 255, 255];
const BLACK: RGB = [0, 0, 0];
const luminance = (c: RGB) => (0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2]) / 255;

function hexToRgb(hex: string): RGB {
  const m = /^#?([0-9a-f]{6})$/i.exec(hex.trim());
  if (!m) return [255, 20, 147];
  const v = parseInt(m[1], 16);
  return [(v >> 16) & 255, (v >> 8) & 255, v & 255];
}

interface VAtom {
  x: number; y: number; z: number; // Å, centred on the molecule
  el: string;
  color: RGB;
  hex: string;
  isH: boolean;
  label: string;
}

interface VBond {
  a: number;
  b: number;
  order: number;
  /** Unit vector perpendicular to the bond, in the plane of the substituents (for multiple bonds). */
  perp: [number, number, number];
  dashed: boolean;
}

/** Draw-list entry for the painter's algorithm. */
interface Item {
  depth: number;
  kind: 0 | 1; // 0 = sphere, 1 = half bond
  i: number; // atom index (sphere) or bond index (half bond)
  from: number; // half bond: atom the half starts at
}

type Mat3 = Float64Array;
const identity = (): Mat3 => Float64Array.of(1, 0, 0, 0, 1, 0, 0, 0, 1);
function axisAngle(x: number, y: number, z: number, th: number): Mat3 {
  const c = Math.cos(th), s = Math.sin(th), t = 1 - c;
  return Float64Array.of(t * x * x + c, t * x * y - s * z, t * x * z + s * y, t * x * y + s * z, t * y * y + c, t * y * z - s * x, t * x * z - s * y, t * y * z + s * x, t * z * z + c);
}
function mul(a: Mat3, b: Mat3): Mat3 {
  const r = new Float64Array(9);
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) r[i * 3 + j] = a[i * 3] * b[j] + a[i * 3 + 1] * b[3 + j] + a[i * 3 + 2] * b[6 + j];
  return r;
}
/**
 * Proper rotation (det = +1) whose rows are the principal axes of a point cloud, largest spread first:
 * the molecule is shown face-on with its long axis horizontal. Never a reflection, so chirality is kept.
 */
function principalFrame(pts: { x: number; y: number; z: number }[]): Mat3 {
  const a = new Float64Array(9);
  for (const p of pts) {
    const v = [p.x, p.y, p.z];
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) a[i * 3 + j] += v[i] * v[j];
  }
  const V = identity();
  for (let sweep = 0; sweep < 30; sweep++) {
    const off = a[1] * a[1] + a[2] * a[2] + a[5] * a[5];
    if (off < 1e-14) break;
    for (const [p, q] of [[0, 1], [0, 2], [1, 2]]) {
      const apq = a[p * 3 + q];
      if (Math.abs(apq) < 1e-14) continue;
      const th = (a[q * 3 + q] - a[p * 3 + p]) / (2 * apq);
      const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1));
      const c = 1 / Math.sqrt(t * t + 1), sn = t * c;
      for (let k = 0; k < 3; k++) {
        const akp = a[k * 3 + p], akq = a[k * 3 + q];
        a[k * 3 + p] = c * akp - sn * akq;
        a[k * 3 + q] = sn * akp + c * akq;
      }
      for (let k = 0; k < 3; k++) {
        const apk = a[p * 3 + k], aqk = a[q * 3 + k];
        a[p * 3 + k] = c * apk - sn * aqk;
        a[q * 3 + k] = sn * apk + c * aqk;
      }
      for (let k = 0; k < 3; k++) {
        const vkp = V[k * 3 + p], vkq = V[k * 3 + q];
        V[k * 3 + p] = c * vkp - sn * vkq;
        V[k * 3 + q] = sn * vkp + c * vkq;
      }
    }
  }
  const order = [0, 1, 2].sort((i, j) => a[j * 3 + j] - a[i * 3 + i]);
  const e1 = [V[order[0]], V[3 + order[0]], V[6 + order[0]]];
  const e2 = [V[order[1]], V[3 + order[1]], V[6 + order[1]]];
  const e3 = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]];
  const m = Float64Array.of(e1[0], e1[1], e1[2], e2[0], e2[1], e2[2], e3[0], e3[1], e3[2]);
  return m.every(Number.isFinite) ? orthonormalise(m) : identity();
}

/** Re-orthonormalise a rotation matrix (prevents drift after many incremental rotations). */
function orthonormalise(m: Mat3): Mat3 {
  const r0 = [m[0], m[1], m[2]], r1 = [m[3], m[4], m[5]];
  const n0 = Math.hypot(r0[0], r0[1], r0[2]);
  for (let k = 0; k < 3; k++) r0[k] /= n0;
  const d = r0[0] * r1[0] + r0[1] * r1[1] + r0[2] * r1[2];
  for (let k = 0; k < 3; k++) r1[k] -= d * r0[k];
  const n1 = Math.hypot(r1[0], r1[1], r1[2]);
  for (let k = 0; k < 3; k++) r1[k] /= n1;
  const r2 = [r0[1] * r1[2] - r0[2] * r1[1], r0[2] * r1[0] - r0[0] * r1[2], r0[0] * r1[1] - r0[1] * r1[0]];
  return Float64Array.of(r0[0], r0[1], r0[2], r1[0], r1[1], r1[2], r2[0], r2[1], r2[2]);
}

export class Viewer3D {
  /** Called with a human-readable result whenever a measurement (distance/angle/dihedral) completes. */
  onMeasure?: (text: string) => void;

  private container: HTMLElement;
  private canvas: HTMLCanvasElement;
  private ctx: CanvasRenderingContext2D;
  private atoms: VAtom[] = [];
  private bonds: VBond[] = [];
  private atomBonds: number[][] = [];
  private radius = 5;
  private signature = '';

  // view state
  private rot: Mat3 = identity();
  /** Home orientation (principal axes of the current molecule). */
  private home: Mat3 = identity();
  private rotCount = 0;
  private zoom = 1;
  private panX = 0;
  private panY = 0;
  private style: Style3D = 'ballstick';
  private showLabels = false;
  private showH = true;
  private autoRotate = false;
  private measureMode = false;
  private selection: number[] = [];
  private measureText = '';

  // projection of the last frame (CSS px), used for picking
  private sx = new Float64Array(0);
  private sy = new Float64Array(0);
  private sz = new Float64Array(0);
  private sr = new Float64Array(0);
  private vis = new Uint8Array(0);

  // rendering / animation
  private cssW = 300;
  private cssH = 150;
  private dpr = 1;
  private raf = 0;
  private lastT = 0;
  private velX = 0; // inertia (px per ms of equivalent drag)
  private velY = 0;
  private sprites = new Map<string, HTMLCanvasElement>();
  private bg: RGB = WHITE;
  private fg: RGB = [31, 35, 40];
  private bgCss = '#ffffff';
  private themeKey = '';
  private themeCheckedAt = -1e9;
  private parseCtx: CanvasRenderingContext2D | null = null;

  // input
  private pointers = new Map<number, { x: number; y: number }>();
  private mode: 'rotate' | 'pan' | 'pinch' | null = null;
  private down = { x: 0, y: 0, t: 0, moved: false, button: 0, type: 'mouse' };
  private samples: { t: number; dx: number; dy: number }[] = [];
  private pinch = { dist: 0, mx: 0, my: 0 };
  private lastTap = 0;
  private ro: ResizeObserver | null = null;
  private cleanup: (() => void)[] = [];
  private disposed = false;

  constructor(container: HTMLElement) {
    this.container = container;
    const canvas = document.createElement('canvas');
    canvas.className = 'viewer3d-canvas';
    Object.assign(canvas.style, {
      display: 'block', width: '100%', height: '100%', touchAction: 'none', outline: 'none',
      userSelect: 'none', webkitUserSelect: 'none', cursor: 'grab',
    } as Partial<CSSStyleDeclaration>);
    canvas.setAttribute('role', 'img');
    canvas.setAttribute('aria-label', '3D molecular model');
    // keyboard focus: use the container if it is focusable, else the canvas itself
    if (!container.hasAttribute('tabindex')) canvas.tabIndex = 0;
    container.appendChild(canvas);
    this.canvas = canvas;
    const ctx = canvas.getContext('2d');
    if (!ctx) throw new Error('Canvas 2D context unavailable');
    this.ctx = ctx;

    this.listen(canvas, 'pointerdown', (e) => this.onPointerDown(e as PointerEvent));
    this.listen(canvas, 'pointermove', (e) => this.onPointerMove(e as PointerEvent));
    this.listen(canvas, 'pointerup', (e) => this.onPointerUp(e as PointerEvent));
    this.listen(canvas, 'pointercancel', (e) => this.onPointerUp(e as PointerEvent, true));
    this.listen(canvas, 'wheel', (e) => this.onWheel(e as WheelEvent), { passive: false });
    this.listen(canvas, 'dblclick', (e) => {
      e.preventDefault();
      this.resetView();
    });
    this.listen(canvas, 'contextmenu', (e) => e.preventDefault());
    // keep the page from scrolling / zooming while a finger is on the model (older iOS ignores touch-action)
    this.listen(canvas, 'touchmove', (e) => e.preventDefault(), { passive: false });
    this.listen(container, 'keydown', (e) => this.onKey(e as KeyboardEvent));

    if (typeof ResizeObserver !== 'undefined') {
      this.ro = new ResizeObserver(() => this.updateSize());
      this.ro.observe(container);
    } else {
      this.listen(window, 'resize', () => this.updateSize());
    }
    if (typeof matchMedia === 'function') {
      const mq = matchMedia('(prefers-color-scheme: dark)');
      const onTheme = () => {
        this.themeCheckedAt = -1e9;
        this.requestRender();
      };
      mq.addEventListener?.('change', onTheme);
      this.cleanup.push(() => mq.removeEventListener?.('change', onTheme));
    }
    this.updateSize();
  }

  // ───────────────────────────── public API ─────────────────────────────

  /**
   * Displays a 3D molecule (Å, right-handed, y up). Passing the same molecule again with new coordinates
   * (e.g. during an optimisation) keeps the current view; a different molecule resets it.
   */
  setMolecule(mol3d: Mol | null): void {
    const sig = mol3d ? mol3d.atoms.map((a) => a.el).join('') + '|' + mol3d.bonds.length : '';
    const keepView = sig !== '' && sig === this.signature;
    this.signature = sig;
    this.atoms = [];
    this.bonds = [];
    this.atomBonds = [];
    if (mol3d && mol3d.atoms.length) {
      const n = mol3d.atoms.length;
      let cx = 0, cy = 0, cz = 0;
      for (const a of mol3d.atoms) {
        cx += a.x;
        cy += a.y;
        cz += a.z ?? 0;
      }
      cx /= n;
      cy /= n;
      cz /= n;
      let r = 0;
      this.atoms = mol3d.atoms.map((a, i) => {
        const hex = a.el === 'H' ? '#FFFFFF' : element(a.el)?.color ?? '#FF1493'; // CPK (pseudo atoms: pink)
        const v: VAtom = {
          x: a.x - cx, y: a.y - cy, z: (a.z ?? 0) - cz, el: a.el, hex, color: hexToRgb(hex), isH: a.el === 'H',
          label: `${a.el === 'R' || a.el === '*' ? a.alias ?? a.el : a.el}${i + 1}`,
        };
        r = Math.max(r, Math.hypot(v.x, v.y, v.z));
        return v;
      });
      this.radius = Math.max(r, 1.5);
      this.atomBonds = this.atoms.map(() => []);
      mol3d.bonds.forEach((b) => {
        if (b.a === b.b || b.a < 0 || b.b < 0 || b.a >= n || b.b >= n) return;
        const bi = this.bonds.length;
        this.bonds.push({ a: b.a, b: b.b, order: b.order, perp: [0, 0, 0], dashed: b.order === 0 || b.style === 'hbond' });
        this.atomBonds[b.a].push(bi);
        this.atomBonds[b.b].push(bi);
      });
      for (const b of this.bonds) b.perp = this.bondPerp(b);
    }
    if (!keepView) this.home = this.atoms.length > 2 ? principalFrame(this.atoms) : identity();
    this.sx = new Float64Array(this.atoms.length);
    this.sy = new Float64Array(this.atoms.length);
    this.sz = new Float64Array(this.atoms.length);
    this.sr = new Float64Array(this.atoms.length);
    this.vis = new Uint8Array(this.atoms.length);
    if (!keepView) {
      this.selection = [];
      this.measureText = '';
      this.rot = this.home.slice();
      this.zoom = 1;
      this.panX = this.panY = 0;
      this.velX = this.velY = 0;
    } else {
      this.selection = this.selection.filter((i) => i < this.atoms.length);
      this.updateMeasurement(false);
    }
    this.requestRender();
  }

  setStyle(s: Style3D): void {
    this.style = s;
    this.requestRender();
  }

  /** Shows element symbols with atom numbers (1-based, matching the molecule's atom order). */
  setShowLabels(on: boolean): void {
    this.showLabels = on;
    this.requestRender();
  }

  setShowHydrogens(on: boolean): void {
    this.showH = on;
    if (!on) {
      const before = this.selection.length;
      this.selection = this.selection.filter((i) => !this.atoms[i].isH);
      if (this.selection.length !== before) this.updateMeasurement(false);
    }
    this.requestRender();
  }

  setAutoRotate(on: boolean): void {
    this.autoRotate = on;
    this.requestRender();
  }

  resetView(): void {
    this.rot = this.home.slice();
    this.zoom = 1;
    this.panX = this.panY = 0;
    this.velX = this.velY = 0;
    this.requestRender();
  }

  /** Measurement mode: click 2 atoms → distance, 3 → angle, 4 → dihedral; the result is shown in-canvas. */
  setMeasureMode(on: boolean): void {
    this.measureMode = on;
    this.selection = [];
    this.measureText = '';
    this.canvas.style.cursor = on ? 'crosshair' : 'grab';
    this.requestRender();
  }

  /** Renders the current view at `scale` × the on-screen CSS size and returns a PNG data URL. */
  toPNG(scale = 2): string {
    const c = document.createElement('canvas');
    c.width = Math.max(1, Math.round(this.cssW * scale));
    c.height = Math.max(1, Math.round(this.cssH * scale));
    const ctx = c.getContext('2d');
    if (!ctx) return '';
    this.readTheme(true);
    ctx.setTransform(scale, 0, 0, scale, 0, 0);
    this.draw(ctx, this.cssW, this.cssH, false);
    return c.toDataURL('image/png');
  }

  dispose(): void {
    if (this.disposed) return;
    this.disposed = true;
    if (this.raf) cancelAnimationFrame(this.raf);
    this.raf = 0;
    this.ro?.disconnect();
    for (const f of this.cleanup) f();
    this.cleanup = [];
    this.sprites.clear();
    this.canvas.remove();
  }

  // ───────────────────────────── setup helpers ─────────────────────────────

  private listen(target: EventTarget, type: string, fn: (e: Event) => void, opts?: AddEventListenerOptions): void {
    target.addEventListener(type, fn, opts);
    this.cleanup.push(() => target.removeEventListener(type, fn, opts));
  }

  /** Direction in which the extra lines of a multiple bond are offset: in the plane of a neighbour. */
  private bondPerp(b: VBond): [number, number, number] {
    const A = this.atoms[b.a], B = this.atoms[b.b];
    let ux = B.x - A.x, uy = B.y - A.y, uz = B.z - A.z;
    const l = Math.hypot(ux, uy, uz) || 1;
    ux /= l; uy /= l; uz /= l;
    const tryRef = (center: number, other: number): [number, number, number] | null => {
      for (const bj of this.atomBonds[center]) {
        const nb = this.bonds[bj];
        const o = nb.a === center ? nb.b : nb.a;
        if (o === other) continue;
        const C = this.atoms[o], P = this.atoms[center];
        let vx = C.x - P.x, vy = C.y - P.y, vz = C.z - P.z;
        const d = vx * ux + vy * uy + vz * uz;
        vx -= d * ux; vy -= d * uy; vz -= d * uz;
        const lv = Math.hypot(vx, vy, vz);
        if (lv > 1e-3) return [vx / lv, vy / lv, vz / lv];
      }
      return null;
    };
    const p = tryRef(b.a, b.b) ?? tryRef(b.b, b.a);
    if (p) return p;
    // linear environment: any perpendicular
    const ref = Math.abs(ux) < 0.9 ? [1, 0, 0] : [0, 1, 0];
    let px = uy * ref[2] - uz * ref[1], py = uz * ref[0] - ux * ref[2], pz = ux * ref[1] - uy * ref[0];
    const lp = Math.hypot(px, py, pz) || 1;
    px /= lp; py /= lp; pz /= lp;
    return [px, py, pz];
  }

  private updateSize(): void {
    if (this.disposed) return;
    const rect = this.canvas.getBoundingClientRect();
    const w = Math.max(1, Math.round(rect.width || this.container.clientWidth || 300));
    const h = Math.max(1, Math.round(rect.height || this.container.clientHeight || 150));
    const dpr = typeof window !== 'undefined' && window.devicePixelRatio ? window.devicePixelRatio : 1;
    if (w === this.cssW && h === this.cssH && dpr === this.dpr && this.canvas.width) {
      this.requestRender();
      return;
    }
    this.cssW = w;
    this.cssH = h;
    this.dpr = dpr;
    this.canvas.width = Math.round(w * dpr);
    this.canvas.height = Math.round(h * dpr);
    this.requestRender();
  }

  /** Reads --viewer-bg / --viewer-fg (throttled); clears cached sprites when the background changes. */
  private readTheme(force = false): void {
    const now = typeof performance !== 'undefined' ? performance.now() : Date.now();
    if (!force && now - this.themeCheckedAt < 500) return;
    this.themeCheckedAt = now;
    const cs = getComputedStyle(this.container);
    const dark = typeof matchMedia === 'function' && matchMedia('(prefers-color-scheme: dark)').matches;
    const bg = cs.getPropertyValue('--viewer-bg').trim() || (dark ? '#1b1d21' : '#ffffff');
    const fg = cs.getPropertyValue('--viewer-fg').trim() || (dark ? '#e6e6e6' : '#1f2328');
    const key = bg + '|' + fg;
    if (key === this.themeKey) return;
    this.themeKey = key;
    this.bg = this.parseColor(bg, WHITE);
    this.fg = this.parseColor(fg, BLACK);
    this.bgCss = css(this.bg);
    this.sprites.clear();
  }

  /** Converts any CSS colour to RGB using the canvas colour parser. */
  private parseColor(c: string, fallback: RGB): RGB {
    if (!this.parseCtx) {
      const pc = document.createElement('canvas');
      pc.width = pc.height = 1;
      this.parseCtx = pc.getContext('2d', { willReadFrequently: true } as CanvasRenderingContext2DSettings);
    }
    const ctx = this.parseCtx;
    if (!ctx) return fallback;
    ctx.clearRect(0, 0, 1, 1);
    ctx.fillStyle = '#000';
    ctx.fillStyle = c;
    ctx.fillRect(0, 0, 1, 1);
    const d = ctx.getImageData(0, 0, 1, 1).data;
    if (d[3] === 0) return fallback;
    return [d[0], d[1], d[2]];
  }

  // ───────────────────────────── animation ─────────────────────────────

  private requestRender(): void {
    if (this.disposed || this.raf) return;
    this.raf = requestAnimationFrame((t) => this.tick(t));
  }

  private tick(t: number): void {
    this.raf = 0;
    const dt = this.lastT ? clamp(t - this.lastT, 1, 50) : 16;
    this.lastT = t;
    let animating = false;
    if (this.autoRotate && this.pointers.size === 0) {
      this.rotateAbout(0, 1, 0, 0.0006 * dt);
      animating = true;
    }
    if (this.velX || this.velY) {
      if (this.pointers.size === 0) {
        this.rotateByDrag(this.velX * dt, this.velY * dt);
        const decay = Math.pow(0.94, dt / 16.7);
        this.velX *= decay;
        this.velY *= decay;
        if (Math.hypot(this.velX, this.velY) < 0.004) this.velX = this.velY = 0;
        else animating = true;
      }
    }
    if (this.dpr !== (window.devicePixelRatio || 1)) this.updateSize();
    this.render();
    if (animating) this.requestRender();
    else this.lastT = 0;
  }

  private render(): void {
    this.readTheme();
    const ctx = this.ctx;
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    this.draw(ctx, this.cssW, this.cssH, true);
  }

  private rotateAbout(x: number, y: number, z: number, angle: number): void {
    this.rot = mul(axisAngle(x, y, z, angle), this.rot);
  }

  /** Trackball: a drag of (dx, dy) CSS px rotates about the in-screen axis perpendicular to it. */
  private rotateByDrag(dx: number, dy: number): void {
    const l = Math.hypot(dx, dy);
    if (l < 1e-9) return;
    // screen y points down, view y points up: drag right → rotate about +y, drag down → about +x
    this.rotateAbout(dy / l, dx / l, 0, l * 0.01);
    if (++this.rotCount % 64 === 0) this.rot = orthonormalise(this.rot);
  }

  // ───────────────────────────── drawing ─────────────────────────────

  private atomRadius(a: VAtom): number {
    switch (this.style) {
      case 'spacefill': return vdw(a.el);
      case 'ballstick': return Math.max(0.2, 0.22 * vdw(a.el));
      case 'sticks': return 0.16;
      default: return 0.08;
    }
  }

  /** Projects all atoms; returns px-per-Å at the centre plane. */
  private project(w: number, h: number): number {
    const R = this.rot;
    const D = Math.max(this.radius, 3) * 4; // camera distance (Å)
    const pad = this.style === 'spacefill' ? 1.8 : 0.8;
    const s0 = (0.5 * Math.min(w, h) * this.zoom) / (this.radius + pad);
    const cx = w / 2 + this.panX, cy = h / 2 + this.panY;
    for (let i = 0; i < this.atoms.length; i++) {
      const a = this.atoms[i];
      const x = R[0] * a.x + R[1] * a.y + R[2] * a.z;
      const y = R[3] * a.x + R[4] * a.y + R[5] * a.z;
      const z = R[6] * a.x + R[7] * a.y + R[8] * a.z;
      const p = D / Math.max(D - z, 0.1 * D);
      this.sx[i] = cx + x * s0 * p;
      this.sy[i] = cy - y * s0 * p;
      this.sz[i] = z;
      this.sr[i] = this.atomRadius(a) * s0 * p;
      this.vis[i] = this.showH || !a.isH ? 1 : 0;
    }
    return s0;
  }

  /** Screen position of a point given in molecule coordinates (Å, centred). */
  private projectPoint(x0: number, y0: number, z0: number, w: number, h: number, s0: number): [number, number, number, number] {
    const R = this.rot;
    const D = Math.max(this.radius, 3) * 4;
    const x = R[0] * x0 + R[1] * y0 + R[2] * z0;
    const y = R[3] * x0 + R[4] * y0 + R[5] * z0;
    const z = R[6] * x0 + R[7] * y0 + R[8] * z0;
    const p = D / Math.max(D - z, 0.1 * D);
    return [w / 2 + this.panX + x * s0 * p, h / 2 + this.panY - y * s0 * p, z, p];
  }

  private draw(ctx: CanvasRenderingContext2D, w: number, h: number, interactive: boolean): void {
    ctx.save();
    ctx.fillStyle = this.bgCss;
    ctx.fillRect(0, 0, w, h);
    if (!this.atoms.length) {
      ctx.fillStyle = css(this.fg, 0.55);
      ctx.font = '13px system-ui, sans-serif';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillText('No 3D structure', w / 2, h / 2);
      ctx.restore();
      return;
    }
    const s0 = this.project(w, h);
    const n = this.atoms.length;
    let zmin = Infinity, zmax = -Infinity;
    for (let i = 0; i < n; i++) {
      if (!this.vis[i]) continue;
      zmin = Math.min(zmin, this.sz[i]);
      zmax = Math.max(zmax, this.sz[i]);
    }
    const zr = Math.max(zmax - zmin, 1e-6);
    const fogOf = (z: number) => (FOG_MAX * (zmax - z)) / zr;

    // build the painter's list
    const items: Item[] = [];
    const wire = this.style === 'wireframe';
    for (let i = 0; i < n; i++) {
      if (!this.vis[i]) continue;
      // wireframe draws no spheres, except for atoms without any visible bond (ions, water O in hidden-H mode…)
      if (wire && this.atomBonds[i].some((bi) => this.vis[this.bonds[bi].a] && this.vis[this.bonds[bi].b])) continue;
      items.push({ depth: this.sz[i], kind: 0, i, from: i });
    }
    if (this.style !== 'spacefill') {
      this.bonds.forEach((b, bi) => {
        if (!this.vis[b.a] || !this.vis[b.b]) return;
        for (const from of [b.a, b.b]) {
          const to = from === b.a ? b.b : b.a;
          // depth of the half's midpoint (from the sphere surface to the bond midpoint)
          items.push({ depth: (3 * this.sz[from] + this.sz[to]) / 4 + 1e-6, kind: 1, i: bi, from });
        }
      });
    }
    items.sort((p, q) => p.depth - q.depth);

    const labelSize = this.style === 'wireframe' || this.style === 'sticks';
    for (const it of items) {
      if (it.kind === 0) {
        const i = it.i;
        const r = wire ? Math.max(this.sr[i], 3) : this.sr[i];
        if (r < 0.3) continue;
        const f = fogOf(this.sz[i]);
        ctx.drawImage(this.sprite(this.atoms[i].color, Math.round(f * 20)), this.sx[i] - r, this.sy[i] - r, 2 * r, 2 * r);
        if (this.showLabels && !labelSize) this.drawLabel(ctx, i, f, false);
      } else {
        this.drawHalfBond(ctx, it.i, it.from, s0, w, h, fogOf);
      }
    }
    if (this.showLabels && labelSize) {
      // thin styles: labels on top, nearest last
      const order = [...Array(n).keys()].filter((i) => this.vis[i]).sort((p, q) => this.sz[p] - this.sz[q]);
      for (const i of order) this.drawLabel(ctx, i, fogOf(this.sz[i]), true);
    }
    this.drawOverlay(ctx, w, interactive);
    ctx.restore();
  }

  private drawHalfBond(ctx: CanvasRenderingContext2D, bi: number, from: number, s0: number, w: number, h: number, fogOf: (z: number) => number): void {
    const b = this.bonds[bi];
    const to = from === b.a ? b.b : b.a;
    const A = this.atoms[from], B = this.atoms[to];
    const dx = B.x - A.x, dy = B.y - A.y, dz = B.z - A.z;
    const L = Math.hypot(dx, dy, dz) || 1;
    const ux = dx / L, uy = dy / L, uz = dz / L;
    const style = this.style;
    // start at the sphere surface (slightly inside) so bonds emerge from the right place in perspective
    const r0 = style === 'ballstick' ? 0.85 * this.atomRadius(A) : 0;
    const mx = (A.x + B.x) / 2, my = (A.y + B.y) / 2, mz = (A.z + B.z) / 2;
    const multi = style === 'ballstick' || style === 'wireframe';
    const order = b.dashed ? 1 : b.order >= 3 ? 3 : b.order >= 2 ? 2 : b.order === 1.5 ? 2 : 1;
    const lines = multi ? order : 1;
    const bondR = style === 'sticks' ? 0.16 : style === 'wireframe' ? 0 : lines === 1 ? 0.11 : 0.075;
    const sep = style === 'wireframe' ? 0.16 : 0.19;
    const color = A.color;
    const f = fogOf((this.sz[from] + this.sz[to]) / 2);
    const base = mix(color, this.bg, f);
    for (let k = 0; k < lines; k++) {
      const off = (k - (lines - 1) / 2) * sep;
      const ox = b.perp[0] * off, oy = b.perp[1] * off, oz = b.perp[2] * off;
      const p1 = this.projectPoint(A.x + ux * r0 + ox, A.y + uy * r0 + oy, A.z + uz * r0 + oz, w, h, s0);
      const p2 = this.projectPoint(mx + ox, my + oy, mz + oz, w, h, s0);
      const aromaticSecond = b.order === 1.5 && k === lines - 1 && lines > 1;
      if (style === 'wireframe') {
        // thin lines need contrast: tint colours close to the background (white H on white) towards the foreground
        const c = Math.abs(luminance(color) - luminance(this.bg)) < 0.3 ? mix(color, this.fg, 0.55) : color;
        ctx.strokeStyle = css(mix(c, this.bg, f));
        ctx.lineWidth = Math.max(1, 1.6 * ((p1[3] + p2[3]) / 2));
        ctx.lineCap = 'round';
        ctx.setLineDash(b.dashed || aromaticSecond ? [3, 3] : []);
        ctx.beginPath();
        ctx.moveTo(p1[0], p1[1]);
        ctx.lineTo(p2[0], p2[1]);
        ctx.stroke();
        ctx.setLineDash([]);
        continue;
      }
      const width = Math.max(1, 2 * bondR * s0 * ((p1[3] + p2[3]) / 2) * (aromaticSecond ? 0.7 : 1));
      if (b.dashed) {
        ctx.strokeStyle = css(base, 0.8);
        ctx.lineWidth = Math.max(1, width * 0.35);
        ctx.setLineDash([4, 4]);
        ctx.beginPath();
        ctx.moveTo(p1[0], p1[1]);
        ctx.lineTo(p2[0], p2[1]);
        ctx.stroke();
        ctx.setLineDash([]);
        continue;
      }
      this.cylinder(ctx, p1[0], p1[1], p2[0], p2[1], width, base, style === 'sticks');
    }
  }

  /** Shaded cylinder: a thick line stroked with a gradient across its width (light from the upper left). */
  private cylinder(ctx: CanvasRenderingContext2D, x1: number, y1: number, x2: number, y2: number, width: number, c: RGB, round: boolean): void {
    const dx = x2 - x1, dy = y2 - y1;
    const L = Math.hypot(dx, dy);
    if (L < 0.25 && !round) return;
    let nx = L > 1e-6 ? -dy / L : 0, ny = L > 1e-6 ? dx / L : 1;
    if (nx + ny > 0) {
      nx = -nx;
      ny = -ny;
    }
    const hw = width / 2;
    const g = ctx.createLinearGradient(x1 + nx * hw, y1 + ny * hw, x1 - nx * hw, y1 - ny * hw);
    g.addColorStop(0, css(mix(c, BLACK, 0.35)));
    g.addColorStop(0.28, css(mix(c, WHITE, 0.45)));
    g.addColorStop(0.55, css(c));
    g.addColorStop(1, css(mix(c, BLACK, 0.5)));
    ctx.strokeStyle = g;
    ctx.lineWidth = width;
    ctx.lineCap = round ? 'round' : 'butt';
    ctx.beginPath();
    ctx.moveTo(x1, y1);
    ctx.lineTo(x2, y2);
    ctx.stroke();
  }

  /** Cached radial-gradient sphere image for a colour and fog level (0…20 → 0…1 of FOG_MAX). */
  private sprite(color: RGB, fogLevel: number): HTMLCanvasElement {
    const key = `${color.join(',')}|${fogLevel}`;
    let c = this.sprites.get(key);
    if (c) return c;
    c = document.createElement('canvas');
    c.width = c.height = SPRITE;
    const g = c.getContext('2d')!;
    const f = fogLevel / 20;
    const fogged = (x: RGB) => css(mix(x, this.bg, f));
    const R = SPRITE / 2;
    const grad = g.createRadialGradient(R * 0.7, R * 0.62, R * 0.06, R, R, R);
    const light = luminance(color) > 0.85;
    grad.addColorStop(0, fogged(mix(color, WHITE, 0.75)));
    grad.addColorStop(0.3, fogged(mix(color, WHITE, 0.2)));
    grad.addColorStop(0.75, fogged(light ? mix(color, BLACK, 0.12) : color));
    grad.addColorStop(1, fogged(mix(color, BLACK, light ? 0.38 : 0.5)));
    g.fillStyle = grad;
    g.beginPath();
    g.arc(R, R, R - 1, 0, Math.PI * 2);
    g.fill();
    if (light) {
      // subtle outline so white hydrogens stay visible on light backgrounds
      g.strokeStyle = css(mix([110, 110, 110], this.bg, f), 0.55);
      g.lineWidth = SPRITE * 0.025;
      g.beginPath();
      g.arc(R, R, R - 1 - g.lineWidth / 2, 0, Math.PI * 2);
      g.stroke();
    }
    if (this.sprites.size > 800) this.sprites.clear();
    this.sprites.set(key, c);
    return c;
  }

  private drawLabel(ctx: CanvasRenderingContext2D, i: number, fog: number, halo: boolean): void {
    const a = this.atoms[i];
    const r = this.sr[i];
    const size = halo ? 11 : clamp(r * 0.75, 8, 16);
    ctx.font = `600 ${size}px system-ui, sans-serif`;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    const x = this.sx[i], y = this.sy[i];
    if (halo) {
      ctx.lineWidth = 3;
      ctx.strokeStyle = css(this.bg, 0.85);
      ctx.strokeText(a.label, x, y);
      ctx.fillStyle = css(mix(this.fg, this.bg, fog));
    } else {
      ctx.fillStyle = luminance(a.color) > 0.55 ? css(mix([20, 20, 20], this.bg, fog)) : css(mix(WHITE, this.bg, fog));
    }
    ctx.fillText(a.label, x, y);
  }

  private drawOverlay(ctx: CanvasRenderingContext2D, w: number, interactive: boolean): void {
    if (!this.measureMode && !this.selection.length) return;
    const acc = css(ACCENT);
    // dashed lines between selected atoms
    if (this.selection.length > 1) {
      ctx.strokeStyle = acc;
      ctx.lineWidth = 1.5;
      ctx.setLineDash([5, 4]);
      ctx.beginPath();
      this.selection.forEach((i, k) => (k ? ctx.lineTo(this.sx[i], this.sy[i]) : ctx.moveTo(this.sx[i], this.sy[i])));
      ctx.stroke();
      ctx.setLineDash([]);
    }
    this.selection.forEach((i, k) => {
      const r = Math.max(this.sr[i], 5) + 3;
      ctx.strokeStyle = acc;
      ctx.lineWidth = 2.5;
      ctx.beginPath();
      ctx.arc(this.sx[i], this.sy[i], r, 0, Math.PI * 2);
      ctx.stroke();
      ctx.fillStyle = acc;
      ctx.font = '600 11px system-ui, sans-serif';
      ctx.textAlign = 'left';
      ctx.textBaseline = 'bottom';
      ctx.fillText(String(k + 1), this.sx[i] + r * 0.75, this.sy[i] - r * 0.75);
    });
    const text = this.measureText || (this.measureMode && interactive ? 'Measure: click 2 atoms (distance), 3 (angle) or 4 (dihedral)' : '');
    if (!text) return;
    ctx.font = '600 12px system-ui, sans-serif';
    const tw = Math.min(ctx.measureText(text).width, w - 24);
    const bx = 8, by = 8, bh = 24;
    ctx.fillStyle = css(this.bg, 0.85);
    ctx.strokeStyle = css(this.measureText ? ACCENT : this.fg, this.measureText ? 0.9 : 0.25);
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.roundRect?.(bx, by, tw + 16, bh, 6);
    if (!ctx.roundRect) ctx.rect(bx, by, tw + 16, bh);
    ctx.fill();
    ctx.stroke();
    ctx.fillStyle = css(this.fg, this.measureText ? 1 : 0.7);
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';
    ctx.fillText(text, bx + 8, by + bh / 2, w - 24);
  }

  // ───────────────────────────── measurement ─────────────────────────────

  private pick(x: number, y: number): number {
    let best = -1, bestZ = -Infinity;
    for (let i = 0; i < this.atoms.length; i++) {
      if (!this.vis[i]) continue;
      const r = Math.max(this.sr[i], 7);
      const dx = x - this.sx[i], dy = y - this.sy[i];
      if (dx * dx + dy * dy <= r * r && this.sz[i] > bestZ) {
        bestZ = this.sz[i];
        best = i;
      }
    }
    return best;
  }

  private handleClick(x: number, y: number): void {
    if (!this.measureMode) return;
    const i = this.pick(x, y);
    if (i < 0) {
      this.selection = [];
      this.measureText = '';
      this.requestRender();
      return;
    }
    const at = this.selection.indexOf(i);
    if (at >= 0) this.selection.splice(at, 1);
    else {
      if (this.selection.length >= 4) this.selection = [];
      this.selection.push(i);
    }
    this.updateMeasurement(true);
    this.requestRender();
  }

  private updateMeasurement(notify: boolean): void {
    const s = this.selection.map((i) => this.atoms[i]);
    const names = this.selection.map((i) => this.atoms[i].label).join('–');
    const sub = (a: VAtom, b: VAtom) => [a.x - b.x, a.y - b.y, a.z - b.z];
    const dot = (u: number[], v: number[]) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    const cross = (u: number[], v: number[]) => [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    const deg = (r: number) => (r * 180) / Math.PI;
    let text = '';
    if (s.length === 2) {
      text = `Distance ${names}: ${Math.hypot(...sub(s[0], s[1])).toFixed(3)} Å`;
    } else if (s.length === 3) {
      const u = sub(s[0], s[1]), v = sub(s[2], s[1]);
      const c = dot(u, v) / (Math.hypot(...u) * Math.hypot(...v) || 1);
      text = `Angle ${names}: ${deg(Math.acos(clamp(c, -1, 1))).toFixed(1)}°`;
    } else if (s.length === 4) {
      const b0 = sub(s[1], s[0]), b1 = sub(s[2], s[1]), b2 = sub(s[3], s[2]);
      const n1 = cross(b0, b1), n2 = cross(b1, b2);
      const yv = dot(cross(n1, n2), b1) / (Math.hypot(...b1) || 1);
      text = `Dihedral ${names}: ${deg(Math.atan2(yv, dot(n1, n2))).toFixed(1)}°`;
    }
    const changed = text !== this.measureText;
    this.measureText = text;
    if (notify && text && changed) this.onMeasure?.(text);
  }

  // ───────────────────────────── input ─────────────────────────────

  private localXY(e: { clientX: number; clientY: number }): [number, number] {
    const r = this.canvas.getBoundingClientRect();
    return [e.clientX - r.left, e.clientY - r.top];
  }

  private onPointerDown(e: PointerEvent): void {
    const focusTarget = this.container.hasAttribute('tabindex') ? this.container : this.canvas;
    focusTarget.focus?.({ preventScroll: true });
    try {
      this.canvas.setPointerCapture(e.pointerId);
    } catch {
      /* ignore (synthetic events) */
    }
    const [x, y] = this.localXY(e);
    this.pointers.set(e.pointerId, { x, y });
    this.velX = this.velY = 0;
    if (this.pointers.size === 1) {
      this.mode = e.button === 2 || e.button === 1 || e.shiftKey ? 'pan' : 'rotate';
      this.down = { x, y, t: performance.now(), moved: false, button: e.button, type: e.pointerType };
      this.samples = [];
      if (this.mode === 'rotate' && !this.measureMode) this.canvas.style.cursor = 'grabbing';
    } else if (this.pointers.size === 2) {
      this.mode = 'pinch';
      this.down.moved = true;
      this.startPinch();
    }
    e.preventDefault();
  }

  private startPinch(): void {
    const [p, q] = [...this.pointers.values()];
    this.pinch = { dist: Math.hypot(p.x - q.x, p.y - q.y) || 1, mx: (p.x + q.x) / 2, my: (p.y + q.y) / 2 };
  }

  private onPointerMove(e: PointerEvent): void {
    const p = this.pointers.get(e.pointerId);
    if (!p) return;
    const [x, y] = this.localXY(e);
    const dx = x - p.x, dy = y - p.y;
    p.x = x;
    p.y = y;
    if (Math.hypot(x - this.down.x, y - this.down.y) > 4) this.down.moved = true;
    if (this.mode === 'rotate' && this.pointers.size === 1) {
      this.rotateByDrag(dx, dy);
      const t = performance.now();
      this.samples.push({ t, dx, dy });
      while (this.samples.length && t - this.samples[0].t > 100) this.samples.shift();
    } else if (this.mode === 'pan') {
      this.panX += dx;
      this.panY += dy;
    } else if (this.mode === 'pinch' && this.pointers.size >= 2) {
      const [a, b] = [...this.pointers.values()];
      const dist = Math.hypot(a.x - b.x, a.y - b.y) || 1;
      const mx = (a.x + b.x) / 2, my = (a.y + b.y) / 2;
      this.zoomAround(dist / this.pinch.dist, mx, my);
      this.panX += mx - this.pinch.mx;
      this.panY += my - this.pinch.my;
      this.pinch = { dist, mx, my };
    } else return;
    this.requestRender();
  }

  private onPointerUp(e: PointerEvent, cancelled = false): void {
    if (!this.pointers.has(e.pointerId)) return;
    this.pointers.delete(e.pointerId);
    try {
      this.canvas.releasePointerCapture(e.pointerId);
    } catch {
      /* ignore */
    }
    if (this.pointers.size === 1 && this.mode === 'pinch') {
      // continue rotating with the remaining finger
      const [rest] = [...this.pointers.values()];
      this.mode = 'rotate';
      this.down = { ...this.down, x: rest.x, y: rest.y, moved: true };
      this.samples = [];
      return;
    }
    if (this.pointers.size > 0) return;
    const [x, y] = this.localXY(e);
    const now = performance.now();
    if (!cancelled && this.mode === 'rotate' && !this.down.moved && now - this.down.t < 500 && this.down.button === 0) {
      // a tap / click
      if (this.down.type !== 'mouse' && now - this.lastTap < 320 && !this.measureMode) {
        this.resetView();
        this.lastTap = 0;
      } else {
        this.lastTap = now;
        this.handleClick(x, y);
      }
    } else if (!cancelled && this.mode === 'rotate' && this.down.moved) {
      // inertia from the last ~80 ms of movement
      const recent = this.samples.filter((s) => now - s.t < 80);
      if (recent.length >= 2) {
        const span = Math.max(16, now - recent[0].t);
        const vx = recent.reduce((s, r) => s + r.dx, 0) / span;
        const vy = recent.reduce((s, r) => s + r.dy, 0) / span;
        if (Math.hypot(vx, vy) > 0.05) {
          this.velX = clamp(vx, -3, 3);
          this.velY = clamp(vy, -3, 3);
        }
      }
    }
    this.mode = null;
    if (!this.measureMode) this.canvas.style.cursor = 'grab';
    this.requestRender();
  }

  private zoomAround(factor: number, x: number, y: number): void {
    const z = clamp(this.zoom * factor, 0.1, 40);
    const f = z / this.zoom;
    // keep the point under the cursor fixed
    const cx = this.cssW / 2 + this.panX, cy = this.cssH / 2 + this.panY;
    this.panX += (x - cx) * (1 - f);
    this.panY += (y - cy) * (1 - f);
    this.zoom = z;
  }

  private onWheel(e: WheelEvent): void {
    e.preventDefault();
    const unit = e.deltaMode === 1 ? 33 : e.deltaMode === 2 ? 400 : 1;
    const [x, y] = this.localXY(e);
    this.zoomAround(Math.exp(-e.deltaY * unit * 0.0015), x, y);
    this.requestRender();
  }

  private onKey(e: KeyboardEvent): void {
    if (e.ctrlKey || e.metaKey || e.altKey) return;
    const step = (e.shiftKey ? 15 : 5) * (Math.PI / 180);
    switch (e.key) {
      case 'ArrowLeft': this.rotateAbout(0, 1, 0, -step); break;
      case 'ArrowRight': this.rotateAbout(0, 1, 0, step); break;
      case 'ArrowUp': this.rotateAbout(1, 0, 0, -step); break;
      case 'ArrowDown': this.rotateAbout(1, 0, 0, step); break;
      case '+': case '=': this.zoomAround(1.15, this.cssW / 2 + this.panX, this.cssH / 2 + this.panY); break;
      case '-': case '_': this.zoomAround(1 / 1.15, this.cssW / 2 + this.panX, this.cssH / 2 + this.panY); break;
      case 'Home': case '0': this.resetView(); break;
      case 'Escape':
        if (!this.selection.length) return;
        this.selection = [];
        this.measureText = '';
        break;
      default: return;
    }
    e.preventDefault();
    this.velX = this.velY = 0;
    this.requestRender();
  }
}
