// Limited-memory BFGS minimiser with a backtracking (Armijo + interpolation) line search.

/** Objective: returns f(x) and writes ∇f(x) into g (g must be fully overwritten). */
export type Objective = (x: Float64Array, g: Float64Array) => number;

export interface MinimizeOptions {
  maxIter?: number;
  /** Convergence threshold on the RMS gradient component. */
  gradTol?: number;
  /** Maximum displacement of any single point per step (same units as x). */
  maxStep?: number;
  /** Coordinates per point (3 for xyz, 4 for 4D embedding) — used for step capping. */
  dim?: number;
  /** Number of correction pairs kept. */
  memory?: number;
  onProgress?: (iteration: number, value: number) => void;
}

export interface MinimizeResult {
  value: number;
  iterations: number;
  converged: boolean;
  rmsGrad: number;
}

const dot = (a: Float64Array, b: Float64Array): number => {
  let s = 0;
  for (let i = 0; i < a.length; i++) s += a[i] * b[i];
  return s;
};

/** Minimises `f` starting from `x` (updated in place). */
export function lbfgs(f: Objective, x: Float64Array, opts: MinimizeOptions = {}): MinimizeResult {
  const n = x.length;
  const maxIter = opts.maxIter ?? 1000;
  const gradTol = opts.gradTol ?? 1e-3;
  const maxStep = opts.maxStep ?? 0.3;
  const dim = opts.dim ?? 3;
  const m = opts.memory ?? 8;
  if (n === 0) return { value: f(x, new Float64Array(0)), iterations: 0, converged: true, rmsGrad: 0 };

  const S: Float64Array[] = [], Y: Float64Array[] = [];
  const rho: number[] = [];
  const alpha = new Float64Array(m);
  let g = new Float64Array(n);
  let gn = new Float64Array(n);
  const xn = new Float64Array(n);
  const d = new Float64Array(n);
  let fx = f(x, g);
  let rms = Math.sqrt(dot(g, g) / n);
  let it = 0;
  let converged = false;
  let stall = 0;

  const resetMemory = () => {
    S.length = 0;
    Y.length = 0;
    rho.length = 0;
  };

  while (it < maxIter) {
    if (!(rms > gradTol)) {
      converged = Number.isFinite(rms);
      break;
    }
    // Two-loop recursion: d = -H·g
    d.set(g);
    const k = S.length;
    for (let i = k - 1; i >= 0; i--) {
      alpha[i] = rho[i] * dot(S[i], d);
      const y = Y[i];
      for (let j = 0; j < n; j++) d[j] -= alpha[i] * y[j];
    }
    const gamma = k ? dot(S[k - 1], Y[k - 1]) / dot(Y[k - 1], Y[k - 1]) : 1;
    for (let j = 0; j < n; j++) d[j] *= gamma;
    for (let i = 0; i < k; i++) {
      const beta = rho[i] * dot(Y[i], d);
      const s = S[i];
      for (let j = 0; j < n; j++) d[j] += (alpha[i] - beta) * s[j];
    }
    for (let j = 0; j < n; j++) d[j] = -d[j];
    let gd = dot(g, d);
    if (!(gd < 0)) {
      // Not a descent direction: fall back to steepest descent.
      resetMemory();
      for (let j = 0; j < n; j++) d[j] = -g[j];
      gd = -dot(g, g);
    }
    // Cap the largest per-point displacement.
    let maxD = 0;
    for (let p = 0; p < n; p += dim) {
      let s2 = 0;
      for (let c = 0; c < dim && p + c < n; c++) s2 += d[p + c] * d[p + c];
      if (s2 > maxD) maxD = s2;
    }
    maxD = Math.sqrt(maxD);
    if (maxD > maxStep) {
      const sc = maxStep / maxD;
      for (let j = 0; j < n; j++) d[j] *= sc;
      gd *= sc;
    }

    // Backtracking line search with quadratic interpolation.
    let step = 1;
    let fn = Infinity;
    let ok = false;
    for (let ls = 0; ls < 30; ls++) {
      for (let j = 0; j < n; j++) xn[j] = x[j] + step * d[j];
      fn = f(xn, gn);
      if (Number.isFinite(fn) && fn <= fx + 1e-4 * step * gd) {
        ok = true;
        break;
      }
      if (!Number.isFinite(fn)) {
        step *= 0.1;
        continue;
      }
      const denom = 2 * (fn - fx - gd * step);
      let ns = denom > 0 ? (-gd * step * step) / denom : 0.5 * step;
      ns = Math.min(Math.max(ns, 0.1 * step), 0.5 * step);
      step = ns;
    }
    it++;
    if (!ok) {
      if (S.length) {
        resetMemory();
        continue;
      }
      break; // steepest descent cannot make progress: we are at (numerical) convergence
    }
    // L-BFGS update
    const s = new Float64Array(n), y = new Float64Array(n);
    for (let j = 0; j < n; j++) {
      s[j] = xn[j] - x[j];
      y[j] = gn[j] - g[j];
    }
    const sy = dot(s, y);
    if (sy > 1e-12) {
      if (S.length === m) {
        S.shift();
        Y.shift();
        rho.shift();
      }
      S.push(s);
      Y.push(y);
      rho.push(1 / sy);
    }
    x.set(xn);
    const tmp = g;
    g = gn;
    gn = tmp;
    const df = fx - fn;
    fx = fn;
    rms = Math.sqrt(dot(g, g) / n);
    opts.onProgress?.(it, fx);
    // Stop when the function no longer changes (flat valley at machine precision).
    if (df <= 1e-13 * Math.max(1, Math.abs(fx))) {
      if (++stall >= 10) break;
    } else stall = 0;
  }
  if (!converged) converged = rms <= gradTol;
  return { value: fx, iterations: it, converged, rmsGrad: rms };
}
