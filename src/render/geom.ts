// Small 2D geometry helpers (model units, y down).

export interface Pt {
  x: number;
  y: number;
}

export const add = (a: Pt, b: Pt): Pt => ({ x: a.x + b.x, y: a.y + b.y });
export const sub = (a: Pt, b: Pt): Pt => ({ x: a.x - b.x, y: a.y - b.y });
export const mul = (a: Pt, s: number): Pt => ({ x: a.x * s, y: a.y * s });
export const len = (a: Pt): number => Math.hypot(a.x, a.y);
export const dist = (a: Pt, b: Pt): number => Math.hypot(a.x - b.x, a.y - b.y);
export const norm = (a: Pt): Pt => {
  const l = Math.hypot(a.x, a.y) || 1;
  return { x: a.x / l, y: a.y / l };
};
/** Perpendicular (rotated +90° in screen coordinates: (1,0) → (0,1)). */
export const perp = (a: Pt): Pt => ({ x: -a.y, y: a.x });
export const dot = (a: Pt, b: Pt): number => a.x * b.x + a.y * b.y;
export const cross = (a: Pt, b: Pt): number => a.x * b.y - a.y * b.x;
export const lerp = (a: Pt, b: Pt, t: number): Pt => ({ x: a.x + (b.x - a.x) * t, y: a.y + (b.y - a.y) * t });
export const angleOf = (a: Pt): number => Math.atan2(a.y, a.x);
export const fromAngle = (t: number, r = 1): Pt => ({ x: Math.cos(t) * r, y: Math.sin(t) * r });
export const rotate = (p: Pt, c: Pt, ang: number): Pt => {
  const s = Math.sin(ang), co = Math.cos(ang);
  const dx = p.x - c.x, dy = p.y - c.y;
  return { x: c.x + dx * co - dy * s, y: c.y + dx * s + dy * co };
};

export function distToSegment(p: Pt, a: Pt, b: Pt): number {
  const dx = b.x - a.x, dy = b.y - a.y;
  const l2 = dx * dx + dy * dy;
  let t = l2 ? ((p.x - a.x) * dx + (p.y - a.y) * dy) / l2 : 0;
  t = Math.max(0, Math.min(1, t));
  return Math.hypot(p.x - (a.x + t * dx), p.y - (a.y + t * dy));
}

export interface Box {
  x1: number;
  y1: number;
  x2: number;
  y2: number;
  /** Clip shape: 'rect' (default) or 'ellipse' inscribed in the box. */
  shape?: 'rect' | 'ellipse';
  /** 'h': the box around an atom label's implicit hydrogens. */
  role?: 'h';
}

/** Parameter t ∈ [0,1] at which the segment p→q leaves the box (0 if p is outside). */
export function exitParam(p: Pt, q: Pt, b: Box): number {
  const dx = q.x - p.x, dy = q.y - p.y;
  if (b.shape === 'ellipse') {
    const cx = (b.x1 + b.x2) / 2, cy = (b.y1 + b.y2) / 2;
    const rx = (b.x2 - b.x1) / 2, ry = (b.y2 - b.y1) / 2;
    if (rx <= 0 || ry <= 0) return 0;
    // solve ((px + t dx - cx)/rx)^2 + ((py + t dy - cy)/ry)^2 = 1
    const ax = (p.x - cx) / rx, ay = (p.y - cy) / ry;
    const bx = dx / rx, by = dy / ry;
    const A = bx * bx + by * by;
    const B = 2 * (ax * bx + ay * by);
    const C = ax * ax + ay * ay - 1;
    if (C > 0) return 0; // p outside
    const disc = B * B - 4 * A * C;
    if (disc < 0 || A === 0) return 0;
    const t = (-B + Math.sqrt(disc)) / (2 * A);
    return Math.max(0, Math.min(1, t));
  }
  // Liang–Barsky
  let t0 = 0, t1 = 1;
  const clip = (pp: number, qq: number) => {
    if (pp === 0) return qq >= 0;
    const r = qq / pp;
    if (pp < 0) {
      if (r > t1) return false;
      if (r > t0) t0 = r;
    } else {
      if (r < t0) return false;
      if (r < t1) t1 = r;
    }
    return true;
  };
  if (
    clip(-dx, p.x - b.x1) &&
    clip(dx, b.x2 - p.x) &&
    clip(-dy, p.y - b.y1) &&
    clip(dy, b.y2 - p.y)
  ) {
    if (t0 <= 1e-9) return Math.min(1, t1);
  }
  return 0;
}

/** Shortens segment a→b so it starts outside boxesA and ends outside boxesB. Returns null if nothing is left. */
export function clipSegment(a: Pt, b: Pt, boxesA: Box[] | null, boxesB: Box[] | null): [Pt, Pt] | null {
  let t0 = 0;
  let t1 = 1;
  if (boxesA) for (const bx of boxesA) t0 = Math.max(t0, exitParam(a, b, bx));
  if (boxesB) for (const bx of boxesB) t1 = Math.min(t1, 1 - exitParam(b, a, bx));
  if (t1 - t0 < 0.02) return null;
  return [lerp(a, b, t0), lerp(a, b, t1)];
}

export function bezierPoint(p0: Pt, p1: Pt, p2: Pt, p3: Pt, t: number): Pt {
  const u = 1 - t;
  const a = u * u * u, b = 3 * u * u * t, c = 3 * u * t * t, d = t * t * t;
  return { x: a * p0.x + b * p1.x + c * p2.x + d * p3.x, y: a * p0.y + b * p1.y + c * p2.y + d * p3.y };
}

export function bezierTangent(p0: Pt, p1: Pt, p2: Pt, p3: Pt, t: number): Pt {
  const u = 1 - t;
  return {
    x: 3 * u * u * (p1.x - p0.x) + 6 * u * t * (p2.x - p1.x) + 3 * t * t * (p3.x - p2.x),
    y: 3 * u * u * (p1.y - p0.y) + 6 * u * t * (p2.y - p1.y) + 3 * t * t * (p3.y - p2.y),
  };
}

export function distToBezier(p: Pt, p0: Pt, p1: Pt, p2: Pt, p3: Pt, steps = 24): number {
  let best = Infinity;
  let prev = p0;
  for (let i = 1; i <= steps; i++) {
    const q = bezierPoint(p0, p1, p2, p3, i / steps);
    best = Math.min(best, distToSegment(p, prev, q));
    prev = q;
  }
  return best;
}

export function pointInPolygon(p: Pt, poly: Pt[]): boolean {
  let inside = false;
  for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
    const a = poly[i], b = poly[j];
    if (a.y > p.y !== b.y > p.y && p.x < ((b.x - a.x) * (p.y - a.y)) / (b.y - a.y) + a.x) inside = !inside;
  }
  return inside;
}

/** Snap an angle (radians) to the nearest multiple of `step`. */
export function snapAngle(a: number, step: number): number {
  return Math.round(a / step) * step;
}
