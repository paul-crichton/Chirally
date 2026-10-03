// Small 2D geometry toolkit shared by the layout / clean-up code.
// Coordinates follow the Mol convention (x right, y down); angles are plain atan2 angles in that frame.

export interface Pt {
  x: number;
  y: number;
}

export const TAU = Math.PI * 2;
export const DEG = Math.PI / 180;

/** Normalises an angle to [0, 2π). */
export function normAngle(a: number): number {
  a %= TAU;
  return a < 0 ? a + TAU : a;
}

/** Smallest absolute difference between two angles (0..π). */
export function angleDiff(a: number, b: number): number {
  const d = Math.abs(normAngle(a) - normAngle(b));
  return d > Math.PI ? TAU - d : d;
}

export function dist(a: Pt, b: Pt): number {
  return Math.hypot(a.x - b.x, a.y - b.y);
}

export function dist2(a: Pt, b: Pt): number {
  const dx = a.x - b.x, dy = a.y - b.y;
  return dx * dx + dy * dy;
}

/** z-component of (b − a) × (c − a); its sign tells on which side of line a→b point c lies. */
export function cross3(a: Pt, b: Pt, c: Pt): number {
  return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

export function polar(origin: Pt, angle: number, r = 1): Pt {
  return { x: origin.x + r * Math.cos(angle), y: origin.y + r * Math.sin(angle) };
}

export function angleTo(from: Pt, to: Pt): number {
  return Math.atan2(to.y - from.y, to.x - from.x);
}

/** Circumradius of a regular polygon with n unit edges. */
export function polygonRadius(n: number): number {
  return 1 / (2 * Math.sin(Math.PI / n));
}

/** True if open segments p1p2 and p3p4 properly intersect (shared endpoints do not count). */
export function segmentsCross(p1: Pt, p2: Pt, p3: Pt, p4: Pt, eps = 1e-6): boolean {
  const d1 = cross3(p3, p4, p1);
  const d2 = cross3(p3, p4, p2);
  const d3 = cross3(p1, p2, p3);
  const d4 = cross3(p1, p2, p4);
  return ((d1 > eps && d2 < -eps) || (d1 < -eps && d2 > eps)) && ((d3 > eps && d4 < -eps) || (d3 < -eps && d4 > eps));
}

/** Distance from point p to segment ab. */
export function pointSegDist(p: Pt, a: Pt, b: Pt): number {
  const dx = b.x - a.x, dy = b.y - a.y;
  const l2 = dx * dx + dy * dy;
  let t = l2 > 0 ? ((p.x - a.x) * dx + (p.y - a.y) * dy) / l2 : 0;
  t = Math.max(0, Math.min(1, t));
  return Math.hypot(p.x - (a.x + t * dx), p.y - (a.y + t * dy));
}

/** Rigid 2D transform: p' = R·(reflect ? mirror(p) : p) + t, mirror = (x, −y). */
export interface Xform {
  cos: number;
  sin: number;
  tx: number;
  ty: number;
  reflect: boolean;
}

export function applyXform(t: Xform, p: Pt): Pt {
  const y = t.reflect ? -p.y : p.y;
  return { x: t.cos * p.x - t.sin * y + t.tx, y: t.sin * p.x + t.cos * y + t.ty };
}

/**
 * Least-squares rigid alignment (Kabsch in 2D) mapping `src` onto `dst`.
 * With allowReflection the mirrored fit is also tried and kept when strictly better.
 * Optional per-point weights.
 */
export function fitXform(src: Pt[], dst: Pt[], allowReflection = false, w?: number[]): { xf: Xform; rms: number } {
  const best = fitOne(src, dst, false, w);
  if (!allowReflection) return best;
  const alt = fitOne(src, dst, true, w);
  return alt.rms < best.rms - 1e-9 ? alt : best;
}

function fitOne(src: Pt[], dst: Pt[], reflect: boolean, w?: number[]): { xf: Xform; rms: number } {
  const n = src.length;
  let sw = 0, sx = 0, sy = 0, dx = 0, dy = 0;
  for (let i = 0; i < n; i++) {
    const wi = w ? w[i] : 1;
    const py = reflect ? -src[i].y : src[i].y;
    sw += wi;
    sx += wi * src[i].x; sy += wi * py;
    dx += wi * dst[i].x; dy += wi * dst[i].y;
  }
  if (sw <= 0) return { xf: { cos: 1, sin: 0, tx: 0, ty: 0, reflect }, rms: 0 };
  sx /= sw; sy /= sw; dx /= sw; dy /= sw;
  let a = 0, b = 0;
  for (let i = 0; i < n; i++) {
    const wi = w ? w[i] : 1;
    const px = src[i].x - sx, py = (reflect ? -src[i].y : src[i].y) - sy;
    const qx = dst[i].x - dx, qy = dst[i].y - dy;
    a += wi * (px * qx + py * qy);
    b += wi * (px * qy - py * qx);
  }
  const ang = n > 1 && (a !== 0 || b !== 0) ? Math.atan2(b, a) : 0;
  const c = Math.cos(ang), s = Math.sin(ang);
  const xf: Xform = { cos: c, sin: s, tx: dx - (c * sx - s * sy), ty: dy - (s * sx + c * sy), reflect };
  let err = 0;
  for (let i = 0; i < n; i++) {
    const p = applyXform(xf, src[i]);
    err += (w ? w[i] : 1) * dist2(p, dst[i]);
  }
  return { xf, rms: Math.sqrt(err / sw) };
}

/** Deterministic pseudo-random generator (mulberry32). */
export function rng(seed: number): () => number {
  let s = seed >>> 0;
  return () => {
    s = (s + 0x6d2b79f5) >>> 0;
    let t = s;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** Uniform spatial hash for neighbour queries on a flat [x0, y0, x1, y1, …] array. */
export class Grid {
  private cells = new Map<number, number[]>();
  constructor(private xy: Float64Array, private cell: number, items: Iterable<number>) {
    for (const i of items) this.insert(i);
  }
  private key(cx: number, cy: number): number {
    return (cx + 32768) * 65536 + (cy + 32768);
  }
  insert(i: number): void {
    const cx = Math.floor(this.xy[2 * i] / this.cell), cy = Math.floor(this.xy[2 * i + 1] / this.cell);
    const k = this.key(cx, cy);
    let c = this.cells.get(k);
    if (!c) this.cells.set(k, (c = []));
    c.push(i);
  }
  /** Calls fn for every item within `r` (≤ cell size) of (x, y). */
  near(x: number, y: number, r: number, fn: (j: number, d2: number) => void): void {
    const cx = Math.floor(x / this.cell), cy = Math.floor(y / this.cell);
    const r2 = r * r;
    const span = Math.ceil(r / this.cell);
    for (let gx = cx - span; gx <= cx + span; gx++) {
      for (let gy = cy - span; gy <= cy + span; gy++) {
        const c = this.cells.get(this.key(gx, gy));
        if (!c) continue;
        for (const j of c) {
          const dx = this.xy[2 * j] - x, dy = this.xy[2 * j + 1] - y;
          const d2 = dx * dx + dy * dy;
          if (d2 <= r2) fn(j, d2);
        }
      }
    }
  }
}
