// Iterative distance-constraint relaxation ("position based dynamics" style).
// Each constraint pulls/pushes its two atoms toward a target distance; immovable atoms act as anchors.
// Optional soft repulsion keeps non-bonded atoms apart. Deterministic (fixed iteration order).
import { Grid } from './geom';

export interface DistCon {
  i: number;
  j: number;
  /** Target distance. */
  d: number;
  /** Stiffness 0..1 (fraction of the error corrected per pass). */
  w: number;
  /** Inequality: only act when the atoms are closer than d. */
  min?: boolean;
}

export interface RefineOptions {
  iterations?: number;
  /** movable[i] = false keeps atom i fixed. Default: all movable. */
  movable?: ArrayLike<boolean | number>;
  /** Atoms taking part in repulsion (default: all atoms that appear in constraints). */
  atoms?: number[];
  /** Minimum distance for non-excluded pairs; 0 disables repulsion. */
  repelDist?: number;
  repelWeight?: number;
  /** Pair keys (min*n+max) exempt from repulsion – typically 1-2 and 1-3 pairs. */
  exclude?: Set<number>;
}

export const pairKey = (i: number, j: number, n: number): number => (i < j ? i * n + j : j * n + i);

/** Relaxes coordinates in xy (flat [x0,y0,x1,y1,…], n atoms) toward the constraints. */
export function relax(xy: Float64Array, n: number, cons: DistCon[], opt: RefineOptions = {}): void {
  const iters = opt.iterations ?? 100;
  const mov = opt.movable;
  const canMove = (i: number) => (mov ? !!mov[i] : true);
  const repel = opt.repelDist ?? 0;
  const repW = opt.repelWeight ?? 0.5;
  let atoms = opt.atoms;
  if (!atoms) {
    const s = new Set<number>();
    for (const c of cons) { s.add(c.i); s.add(c.j); }
    atoms = [...s].sort((a, b) => a - b);
  }
  let pairs: number[] = [];
  const rebuild = () => {
    pairs = [];
    if (repel <= 0) return;
    const reach = repel * 1.25;
    const g = new Grid(xy, reach, atoms!);
    for (const i of atoms!) {
      g.near(xy[2 * i], xy[2 * i + 1], reach, (j) => {
        if (j <= i) return;
        if (!canMove(i) && !canMove(j)) return;
        if (opt.exclude && opt.exclude.has(pairKey(i, j, n))) return;
        pairs.push(i, j);
      });
    }
  };
  const project = (i: number, j: number, d: number, w: number, minOnly: boolean) => {
    const mi = canMove(i) ? 1 : 0, mj = canMove(j) ? 1 : 0;
    if (!mi && !mj) return;
    let dx = xy[2 * j] - xy[2 * i];
    let dy = xy[2 * j + 1] - xy[2 * i + 1];
    let r = Math.hypot(dx, dy);
    if (minOnly && r >= d) return;
    if (r < 1e-9) {
      // coincident atoms: separate along a deterministic pseudo-random direction
      const a = ((i * 7919 + j * 104729) % 360) * (Math.PI / 180);
      dx = Math.cos(a) * 1e-3; dy = Math.sin(a) * 1e-3; r = 1e-3;
    }
    const k = ((r - d) / r) * w / (mi + mj);
    if (mi) { xy[2 * i] += dx * k; xy[2 * i + 1] += dy * k; }
    if (mj) { xy[2 * j] -= dx * k; xy[2 * j + 1] -= dy * k; }
  };
  for (let it = 0; it < iters; it++) {
    if (repel > 0 && it % 8 === 0) rebuild();
    for (const c of cons) project(c.i, c.j, c.d, c.w, !!c.min);
    for (let p = 0; p < pairs.length; p += 2) project(pairs[p], pairs[p + 1], repel, repW, true);
  }
}
