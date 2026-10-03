// Distance-geometry embedding ("DG-lite"):
//   1. bounds matrix from topology and ideal UFF geometry (1-2, 1-3, 1-4 incl. E/Z and ring planarity,
//      van der Waals lower bounds for everything else), triangle smoothing;
//   2. random distance matrix within the bounds → metric-matrix (MDS) embedding in 4D;
//   3. refinement of a distance-violation error function with chiral-volume and sp2-planarity terms,
//      first in 4D (lets atoms pass through each other), then squeezing out the 4th dimension, then 3D;
//   4. verification of every tetrahedral / double-bond spec, local repair (umbrella flip, 180° twist)
//      or a retry with a new random distance matrix.
import { Mol } from '../mol';
import { addExplicitHydrogens } from './hydrogens';
import { typeUFF, UFFTyping, uffBondLength } from './uffTyping';
import { makeRng, Rng, gaussian } from './rng';
import { lbfgs } from './minimize';
import { tetraSatisfied, dbSatisfied, invertCenter, flipDoubleBond, centreInside } from './stereo3d';

const BIG = 1000;
/** Approximate van der Waals radii (Å, Bondi) used for non-bonded lower bounds. */
const VDW: Record<string, number> = {
  H: 1.1, C: 1.7, N: 1.55, O: 1.52, F: 1.47, P: 1.8, S: 1.8, Cl: 1.75, Br: 1.85, I: 1.98, B: 1.92, Si: 2.1, Se: 1.9,
};
const vdwRadius = (el: string) => VDW[el] ?? 1.9;
const RAD = Math.PI / 180;

export interface EmbedOptions {
  /** Random seed (default 1): the same seed always gives the same coordinates. */
  seed?: number;
  /** Maximum number of random restarts per connected component (default 10). */
  maxAttempts?: number;
}

/** Returns a NEW Mol with explicit hydrogens and 3D coordinates (Å, right-handed, y up). */
export function embed3D(mol: Mol, opts: EmbedOptions = {}): Mol {
  const m = addExplicitHydrogens(mol);
  embedCoordinates(m, opts);
  return m;
}

/**
 * Computes 3D coordinates for an H-explicit molecule in place (atoms[i].x/y/z).
 * Returns the number of stereo specs that could not be satisfied (0 on success).
 */
export function embedCoordinates(m: Mol, opts: EmbedOptions = {}): number {
  const seed = opts.seed ?? 1;
  const attempts = Math.max(1, opts.maxAttempts ?? 10);
  let failures = 0;
  let xCursor = 0;
  const comps = m.components().sort((a, b) => b.length - a.length);
  for (const comp of comps) {
    let coords: number[][];
    if (comp.length === 1) coords = [[0, 0, 0]];
    else {
      const { mol: sub } = m.subset(comp);
      failures += embedComponent(sub, seed, attempts);
      coords = sub.atoms.map((a) => [a.x, a.y, a.z ?? 0]);
    }
    // centre the component, then place it to the right of the previous one
    const c = [0, 1, 2].map((d) => coords.reduce((s, p) => s + p[d], 0) / coords.length);
    let minX = Infinity, maxX = -Infinity;
    for (const p of coords) {
      for (let d = 0; d < 3; d++) p[d] -= c[d];
      minX = Math.min(minX, p[0]);
      maxX = Math.max(maxX, p[0]);
    }
    const shift = xCursor === 0 && comp === comps[0] ? 0 : xCursor - minX + 3.5;
    comp.forEach((ai, k) => {
      const at = m.atoms[ai];
      at.x = coords[k][0] + shift;
      at.y = coords[k][1];
      at.z = coords[k][2];
    });
    xCursor = maxX + shift;
  }
  // centre everything at the origin
  if (m.atoms.length) {
    const c = [0, 0, 0];
    for (const a of m.atoms) { c[0] += a.x; c[1] += a.y; c[2] += a.z ?? 0; }
    for (const a of m.atoms) { a.x -= c[0] / m.atoms.length; a.y -= c[1] / m.atoms.length; a.z = (a.z ?? 0) - c[2] / m.atoms.length; }
  }
  return failures;
}

// ───────────────────────────── geometry targets ─────────────────────────────

function ringAngle(size: number, planar: boolean, theta0: number): number {
  if (size === 3) return 60;
  if (size === 4) return planar ? 90 : 88.5;
  if (size === 5) return planar ? 108 : 104.5;
  if (planar) return size === 6 ? 120 : Math.min(128, (180 * (size - 2)) / size);
  if (theta0 < 100) return theta0 + 6; // S, Se, P … in saturated rings
  return size === 6 ? 111 : 114;
}

/** Ideal valence angles (degrees) per centre and neighbour pair, corrected for ring membership. */
function idealAngles(mol: Mol, ty: UFFTyping): { get: (j: number, i: number, k: number) => number; strained: Uint8Array; loose: Uint8Array } {
  const n = mol.atoms.length;
  const ringPairs: Map<number, number>[] = mol.atoms.map(() => new Map());
  const key = (a: number, b: number) => (a < b ? a * n + b : b * n + a);
  for (const r of ty.rings.rings) {
    for (let t = 0; t < r.length; t++) {
      const j = r[t], a = r[(t + r.length - 1) % r.length], b = r[(t + 1) % r.length];
      const k = key(a, b);
      const prev = ringPairs[j].get(k);
      if (prev === undefined || r.length < prev) ringPairs[j].set(k, r.length);
    }
  }
  const tables: Map<number, number>[] = [];
  const strained = new Uint8Array(n);
  const loose = new Uint8Array(n);
  for (let j = 0; j < n; j++) {
    const nb = ty.nbrs[j];
    const deg = nb.length;
    const th0 = ty.params[j].theta0;
    const planar = ty.hyb[j] === 2 && deg <= 3;
    const linear = ty.hyb[j] === 1 && deg <= 2;
    const m = new Map<number, number>();
    tables.push(m);
    if (deg < 2) continue;
    const pairs: [number, number][] = [];
    for (let p = 0; p < deg; p++) for (let q = p + 1; q < deg; q++) pairs.push([nb[p], nb[q]]);
    if (linear) {
      for (const [a, b] of pairs) m.set(key(a, b), 180);
      continue;
    }
    if (deg > 4) {
      loose[j] = 1;
      for (const [a, b] of pairs) m.set(key(a, b), 90);
      continue;
    }
    let sumFixed = 0, sumDeficit = 0;
    const base = deg === 4 ? 109.47 : th0;
    for (const [a, b] of pairs) {
      const size = ringPairs[j].get(key(a, b));
      if (size === undefined) continue;
      const ang = ringAngle(size, planar, th0);
      if (size <= 5) strained[j] = 1;
      m.set(key(a, b), ang);
      sumFixed += ang;
      sumDeficit += base - ang;
    }
    const rest = pairs.length - m.size;
    if (rest > 0) {
      let ang: number;
      if (planar && deg === 3) ang = (360 - sumFixed) / rest;
      else if (m.size > 0) ang = base + (sumDeficit / rest) * (deg === 4 ? 1 : 0.5);
      else ang = base;
      ang = Math.max(60, Math.min(180, ang));
      for (const [a, b] of pairs) if (!m.has(key(a, b))) m.set(key(a, b), ang);
    }
  }
  return {
    get: (j, i, k) => tables[j].get(key(i, k)) ?? ty.params[j].theta0,
    strained,
    loose,
  };
}

/** Distance between the ends of i–j–k–l at dihedral φ (radians). */
function dist14(r1: number, r2: number, r3: number, th1: number, th2: number, phi: number): number {
  const ix = r1 * Math.cos(th1), iy = r1 * Math.sin(th1);
  const lx = r2 - r3 * Math.cos(th2), lr = r3 * Math.sin(th2);
  return Math.hypot(lx - ix, lr * Math.cos(phi) - iy, lr * Math.sin(phi));
}

export interface Bounds {
  n: number;
  lo: Float64Array;
  up: Float64Array;
  /** Topological (bond-count) distances, 255 = unreachable. */
  topo: Uint8Array;
}

/** Builds and triangle-smooths the distance bounds matrix of a connected, H-explicit molecule. */
export function buildBounds(mol: Mol, ty: UFFTyping): Bounds {
  const n = mol.atoms.length;
  const nb = ty.nbrs;
  const lo = new Float64Array(n * n);
  const up = new Float64Array(n * n).fill(BIG);
  const set = new Uint8Array(n * n);
  const topo = new Uint8Array(n * n).fill(255);
  // BFS topological distances
  const queue = new Int32Array(n);
  for (let s = 0; s < n; s++) {
    topo[s * n + s] = 0;
    let head = 0, tail = 0;
    queue[tail++] = s;
    while (head < tail) {
      const v = queue[head++];
      const dv = topo[s * n + v];
      if (dv >= 254) continue;
      for (const w of nb[v]) {
        if (topo[s * n + w] === 255) {
          topo[s * n + w] = dv + 1;
          queue[tail++] = w;
        }
      }
    }
  }
  const setRange = (i: number, j: number, l: number, u: number) => {
    const ij = i * n + j, ji = j * n + i;
    if (set[ij]) {
      const nl = Math.max(lo[ij], l), nu = Math.min(up[ij], u);
      if (nl <= nu) { l = nl; u = nu; }
      else { l = Math.min(lo[ij], l); u = Math.max(up[ij], u); }
    }
    lo[ij] = lo[ji] = l;
    up[ij] = up[ji] = u;
    set[ij] = set[ji] = 1;
  };

  // 1-2: UFF rest lengths
  const r0 = new Map<number, number>();
  const rk = (i: number, j: number) => (i < j ? i * n + j : j * n + i);
  mol.bonds.forEach((b, bi) => {
    const bo = ty.bondOrder[bi];
    if (bo <= 0 || b.a === b.b) return;
    const r = uffBondLength(bo, ty.params[b.a], ty.params[b.b]);
    r0.set(rk(b.a, b.b), r);
    setRange(b.a, b.b, r - 0.01, r + 0.01);
  });
  const rest = (i: number, j: number) => r0.get(rk(i, j)) ?? ty.params[i].r1 + ty.params[j].r1;

  // 1-3: law of cosines with ideal angles
  const ang = idealAngles(mol, ty);
  for (let j = 0; j < n; j++) {
    const list = nb[j];
    for (let p = 0; p < list.length; p++) {
      for (let q = p + 1; q < list.length; q++) {
        const i = list[p], k = list[q];
        if (topo[i * n + k] !== 2) continue;
        const r1 = rest(i, j), r2 = rest(j, k);
        if (ang.loose[j]) {
          setRange(i, k, Math.sqrt(r1 * r1 + r2 * r2 - 2 * r1 * r2 * Math.cos(80 * RAD)), r1 + r2);
          continue;
        }
        const th = ang.get(j, i, k) * RAD;
        const d = Math.sqrt(Math.max(0.01, r1 * r1 + r2 * r2 - 2 * r1 * r2 * Math.cos(th)));
        const tol = ang.strained[j] ? 0.08 : 0.04;
        setRange(i, k, d - tol, d + tol);
      }
    }
  }

  // 1-4: torsional ranges; fixed for stereo double bonds and planar ring bonds
  const dbSpecByBond = new Map(mol.dbStereo.map((d) => [d.bond, d]));
  const { rings, bondRings } = ty.rings;
  mol.bonds.forEach((b, bi) => {
    if (ty.bondOrder[bi] <= 0) return;
    const j = b.a, k = b.b;
    const planarBond = ty.hyb[j] === 2 && ty.hyb[k] === 2;
    const multiple = ty.bondOrder[bi] >= 1.5;
    let ring: number[] | null = null;
    for (const r of bondRings[bi]) if (!ring || rings[r].length < ring.length) ring = rings[r];
    const spec = dbSpecByBond.get(bi);
    for (const i of nb[j]) {
      if (i === k) continue;
      for (const l of nb[k]) {
        if (l === j || l === i || topo[i * n + l] !== 3) continue;
        const r1 = rest(i, j), r2 = rest(j, k), r3 = rest(k, l);
        const th1 = ang.get(j, i, k) * RAD, th2 = ang.get(k, j, l) * RAD;
        const at = (phiDeg: number) => dist14(r1, r2, r3, th1, th2, phiDeg * RAD);
        const dc = at(0), dt = at(180);
        if (spec && ty.bondOrder[bi] >= 2) {
          // stereo double bond: substituents on the same side as the spec'd pair are cis/trans accordingly
          const sideI = i === spec.a || i === spec.b ? 0 : 1;
          const sideL = l === spec.a || l === spec.b ? 0 : 1;
          const cis = sideI === sideL ? spec.cis : !spec.cis;
          const d = cis ? dc : dt;
          setRange(i, l, d - 0.04, d + 0.04);
        } else if (planarBond && ring && ring.length <= 8) {
          const cis = ring.includes(i) === ring.includes(l);
          if (multiple) {
            const d = cis ? dc : dt;
            setRange(i, l, d - 0.05, d + 0.05);
          } else if (cis) setRange(i, l, dc - 0.05, at(35) + 0.05);
          else setRange(i, l, at(145) - 0.05, dt + 0.05);
        } else if (ring && ring.includes(i) && ring.includes(l) && ring.length <= 7) {
          setRange(i, l, dc - 0.05, at(ring.length === 6 ? 65 : 90) + 0.05);
        } else {
          setRange(i, l, dc - 0.05, dt + 0.05);
        }
      }
    }
  });

  // everything else: van der Waals lower bounds
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      if (set[i * n + j]) continue;
      const t = topo[i * n + j];
      const s = t === 4 ? 0.7 : 0.8;
      const l = s * (vdwRadius(mol.atoms[i].el) + vdwRadius(mol.atoms[j].el));
      lo[i * n + j] = lo[j * n + i] = l;
    }
  }
  triangleSmooth(lo, up, n);
  return { n, lo, up, topo };
}

/** Floyd–Warshall triangle-inequality smoothing of upper and lower bounds (in place). */
export function triangleSmooth(lo: Float64Array, up: Float64Array, n: number): number {
  let inconsistent = 0;
  for (let k = 0; k < n; k++) {
    const kn = k * n;
    for (let i = 0; i < n; i++) {
      if (i === k) continue;
      const in_ = i * n;
      const uik = up[in_ + k], lik = lo[in_ + k];
      for (let j = i + 1; j < n; j++) {
        if (j === k) continue;
        const ukj = up[kn + j], lkj = lo[kn + j];
        const ij = in_ + j;
        let u = up[ij], l = lo[ij];
        if (u > uik + ukj) u = uik + ukj;
        if (l < lik - ukj) l = lik - ukj;
        else if (l < lkj - uik) l = lkj - uik;
        if (l > u) {
          inconsistent++;
          l = u;
        }
        up[ij] = up[j * n + i] = u;
        lo[ij] = lo[j * n + i] = l;
      }
    }
  }
  return inconsistent;
}

// ───────────────────────────── metric embedding ─────────────────────────────

/** Jacobi eigen-decomposition of a small symmetric matrix (row-major p×p). */
function jacobi(a: Float64Array, p: number): { vals: number[]; vecs: Float64Array } {
  const v = new Float64Array(p * p);
  for (let i = 0; i < p; i++) v[i * p + i] = 1;
  for (let sweep = 0; sweep < 60; sweep++) {
    let off = 0;
    for (let i = 0; i < p; i++) for (let j = i + 1; j < p; j++) off += a[i * p + j] * a[i * p + j];
    if (off < 1e-18) break;
    for (let r = 0; r < p; r++) {
      for (let s = r + 1; s < p; s++) {
        const ars = a[r * p + s];
        if (Math.abs(ars) < 1e-14) continue;
        const theta = (a[s * p + s] - a[r * p + r]) / (2 * ars);
        const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        const c = 1 / Math.sqrt(t * t + 1), sn = t * c;
        for (let k = 0; k < p; k++) {
          const akr = a[k * p + r], aks = a[k * p + s];
          a[k * p + r] = c * akr - sn * aks;
          a[k * p + s] = sn * akr + c * aks;
        }
        for (let k = 0; k < p; k++) {
          const ark = a[r * p + k], ask = a[s * p + k];
          a[r * p + k] = c * ark - sn * ask;
          a[s * p + k] = sn * ark + c * ask;
        }
        for (let k = 0; k < p; k++) {
          const vkr = v[k * p + r], vks = v[k * p + s];
          v[k * p + r] = c * vkr - sn * vks;
          v[k * p + s] = sn * vkr + c * vks;
        }
      }
    }
  }
  return { vals: Array.from({ length: p }, (_, i) => a[i * p + i]), vecs: v };
}

/** Leading eigenpairs of a symmetric n×n matrix by subspace iteration + Rayleigh–Ritz. */
function topEigen(B: Float64Array, n: number, want: number, rng: Rng): { vals: number[]; vecs: Float64Array[] } {
  const p = Math.min(n, want + 2);
  let V: Float64Array[] = Array.from({ length: p }, () => Float64Array.from({ length: n }, () => rng() - 0.5));
  const orthonormalise = (M: Float64Array[]) => {
    for (let c = 0; c < M.length; c++) {
      for (let d = 0; d < c; d++) {
        let s = 0;
        for (let i = 0; i < n; i++) s += M[c][i] * M[d][i];
        for (let i = 0; i < n; i++) M[c][i] -= s * M[d][i];
      }
      let l = 0;
      for (let i = 0; i < n; i++) l += M[c][i] * M[c][i];
      l = Math.sqrt(l) || 1;
      for (let i = 0; i < n; i++) M[c][i] /= l;
    }
  };
  const mul = (v: Float64Array) => {
    const out = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      let s = 0;
      const row = i * n;
      for (let j = 0; j < n; j++) s += B[row + j] * v[j];
      out[i] = s;
    }
    return out;
  };
  orthonormalise(V);
  for (let it = 0; it < 40; it++) {
    V = V.map(mul);
    orthonormalise(V);
  }
  const BV = V.map(mul);
  const H = new Float64Array(p * p);
  for (let a = 0; a < p; a++) for (let b = 0; b < p; b++) {
    let s = 0;
    for (let i = 0; i < n; i++) s += V[a][i] * BV[b][i];
    H[a * p + b] = s;
  }
  const { vals, vecs } = jacobi(H, p);
  const order = vals.map((_, i) => i).sort((a, b) => vals[b] - vals[a]);
  return {
    vals: order.map((i) => vals[i]),
    vecs: order.map((c) => {
      const out = new Float64Array(n);
      for (let a = 0; a < p; a++) {
        const w = vecs[a * p + c];
        for (let i = 0; i < n; i++) out[i] += w * V[a][i];
      }
      return out;
    }),
  };
}

/** Random distance matrix within the bounds, embedded in 4D through the metric matrix. */
function initialCoords(bd: Bounds, rng: Rng): Float64Array {
  const { n, lo, up } = bd;
  const D2 = new Float64Array(n * n);
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const l = lo[i * n + j], u = Math.min(up[i * n + j], l + 8);
      const d = l + rng() * (u - l);
      D2[i * n + j] = D2[j * n + i] = d * d;
    }
  }
  const row = new Float64Array(n);
  let tot = 0;
  for (let i = 0; i < n; i++) {
    let s = 0;
    for (let j = 0; j < n; j++) s += D2[i * n + j];
    row[i] = s / n;
    tot += row[i];
  }
  tot /= n;
  const B = D2; // reuse storage: B = −½ J D² J
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) B[i * n + j] = -0.5 * (D2[i * n + j] - row[i] - row[j] + tot);
  const { vals, vecs } = topEigen(B, n, 4, rng);
  const x = new Float64Array(n * 4);
  for (let d = 0; d < 4; d++) {
    const s = d < vals.length && vals[d] > 1e-6 ? Math.sqrt(vals[d]) : 0;
    for (let i = 0; i < n; i++) x[i * 4 + d] = (s ? s * vecs[d][i] : 0) + 0.3 * gaussian(rng);
  }
  return x;
}

// ───────────────────────────── error function ─────────────────────────────

interface DGTerms {
  n: number;
  pi: Int32Array; pj: Int32Array; lo2: Float64Array; up2: Float64Array;
  /** chiral constraints: 4 atom indices each, desired sign (−1 = '@' / ccw) and minimal |volume| */
  ch: Int32Array; chSign: Float64Array; chMin: Float64Array;
  /** planarity: centre + 3 neighbours */
  pl: Int32Array;
}

function dgTerms(mol: Mol, ty: UFFTyping, bd: Bounds): DGTerms {
  const n = bd.n;
  const pi: number[] = [], pj: number[] = [], lo2: number[] = [], up2: number[] = [];
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const l = bd.lo[i * n + j], u = bd.up[i * n + j];
      if (l <= 0 && u >= BIG) continue;
      pi.push(i); pj.push(j); lo2.push(l * l); up2.push(Math.min(u, BIG) ** 2);
    }
  }
  const ch: number[] = [], chSign: number[] = [], chMin: number[] = [];
  const centres = new Set<number>();
  for (const t of mol.tetra) {
    const lp = t.nbrs.filter((v) => v < 0).length;
    if (lp > 1) continue;
    centres.add(t.center);
    const sign = t.ccw ? -1 : 1;
    if (lp) {
      // a lone pair is represented by the centre itself (same side of the plane of the other three)
      ch.push(...t.nbrs.map((v) => (v < 0 ? t.center : v)));
      chSign.push(sign);
      chMin.push(0.4);
      continue;
    }
    // neighbour volume: at least half of the ideal tetrahedral value
    const vIdeal = idealTetraVolume(t.nbrs.map((v) => restLen(bd, t.center, v)));
    ch.push(...t.nbrs);
    chSign.push(sign);
    chMin.push(0.5 * vIdeal);
    // the centre must lie inside the neighbour tetrahedron: replacing any vertex by the centre keeps the sign
    for (let k = 0; k < 4; k++) {
      ch.push(...t.nbrs.map((v, q) => (q === k ? t.center : v)));
      chSign.push(sign);
      chMin.push(0.1 * vIdeal);
    }
  }
  const pl: number[] = [];
  for (let j = 0; j < n; j++) {
    if (ty.hyb[j] === 2 && ty.nbrs[j].length === 3 && !centres.has(j)) pl.push(j, ...ty.nbrs[j]);
  }
  return {
    n,
    pi: Int32Array.from(pi), pj: Int32Array.from(pj), lo2: Float64Array.from(lo2), up2: Float64Array.from(up2),
    ch: Int32Array.from(ch), chSign: Float64Array.from(chSign), chMin: Float64Array.from(chMin),
    pl: Int32Array.from(pl),
  };
}

/** Rest length of bond i–j taken from the (1-2) bounds. */
function restLen(bd: Bounds, i: number, j: number): number {
  const k = i * bd.n + j;
  return 0.5 * (bd.lo[k] + bd.up[k]);
}

/** |chiral volume| of a regular tetrahedron of neighbours at the given distances from the centre. */
function idealTetraVolume(r: number[]): number {
  const s = 1 / Math.sqrt(3);
  const u = [[s, s, s], [s, -s, -s], [-s, s, -s], [-s, -s, s]];
  const p = u.map((v, k) => v.map((c) => c * r[k]));
  const a = p[1].map((c, d) => c - p[0][d]), b = p[2].map((c, d) => c - p[0][d]), c = p[3].map((c2, d) => c2 - p[0][d]);
  return Math.abs(a[0] * (b[1] * c[2] - b[2] * c[1]) + a[1] * (b[2] * c[0] - b[0] * c[2]) + a[2] * (b[0] * c[1] - b[1] * c[0]));
}

/** DG error function in `dim` (3 or 4) dimensions; `w4` penalises the 4th coordinate. */
function dgObjective(T: DGTerms, dim: number, w4: number, wChiral = 1, wPlanar = 0.2) {
  return (x: Float64Array, g: Float64Array): number => {
    g.fill(0);
    let e = 0;
    const { pi, pj, lo2, up2 } = T;
    for (let p = 0; p < pi.length; p++) {
      const i = pi[p] * dim, j = pj[p] * dim;
      let d2 = 0;
      for (let c = 0; c < dim; c++) {
        const t = x[i + c] - x[j + c];
        d2 += t * t;
      }
      let dE: number;
      if (d2 > up2[p]) {
        const v = d2 / up2[p] - 1;
        e += v * v;
        dE = (2 * v) / up2[p];
      } else if (d2 < lo2[p]) {
        const s = lo2[p] + d2;
        const v = (2 * lo2[p]) / s - 1;
        e += v * v;
        dE = (2 * v * -2 * lo2[p]) / (s * s);
      } else continue;
      for (let c = 0; c < dim; c++) {
        const t = 2 * dE * (x[i + c] - x[j + c]);
        g[i + c] += t;
        g[j + c] -= t;
      }
    }
    // chiral volumes V = (p1 − p0)·((p2 − p0) × (p3 − p0)) on xyz
    const { ch, chSign, chMin } = T;
    for (let t = 0; t < chSign.length; t++) {
      const i0 = ch[4 * t] * dim, i1 = ch[4 * t + 1] * dim, i2 = ch[4 * t + 2] * dim, i3 = ch[4 * t + 3] * dim;
      const ax = x[i1] - x[i0], ay = x[i1 + 1] - x[i0 + 1], az = x[i1 + 2] - x[i0 + 2];
      const bx = x[i2] - x[i0], by = x[i2 + 1] - x[i0 + 1], bz = x[i2 + 2] - x[i0 + 2];
      const cx = x[i3] - x[i0], cy = x[i3 + 1] - x[i0 + 1], cz = x[i3 + 2] - x[i0 + 2];
      const bcx = by * cz - bz * cy, bcy = bz * cx - bx * cz, bcz = bx * cy - by * cx;
      const V = ax * bcx + ay * bcy + az * bcz;
      const target = chSign[t] * chMin[t];
      const viol = chSign[t] < 0 ? V - target : target - V; // > 0 when violated
      if (viol <= 0) continue;
      e += wChiral * viol * viol;
      const dV = 2 * wChiral * viol * (chSign[t] < 0 ? 1 : -1);
      // ∂V/∂p1 = b × c ; ∂V/∂p2 = c × a ; ∂V/∂p3 = a × b ; ∂V/∂p0 = −(sum)
      const cax = cy * az - cz * ay, cay = cz * ax - cx * az, caz = cx * ay - cy * ax;
      const abx = ay * bz - az * by, aby = az * bx - ax * bz, abz = ax * by - ay * bx;
      g[i1] += dV * bcx; g[i1 + 1] += dV * bcy; g[i1 + 2] += dV * bcz;
      g[i2] += dV * cax; g[i2 + 1] += dV * cay; g[i2 + 2] += dV * caz;
      g[i3] += dV * abx; g[i3 + 1] += dV * aby; g[i3 + 2] += dV * abz;
      g[i0] -= dV * (bcx + cax + abx); g[i0 + 1] -= dV * (bcy + cay + aby); g[i0 + 2] -= dV * (bcz + caz + abz);
    }
    // sp2 planarity: volume spanned by the three neighbour vectors → 0
    const pl = T.pl;
    for (let t = 0; t < pl.length; t += 4) {
      const i0 = pl[t] * dim, i1 = pl[t + 1] * dim, i2 = pl[t + 2] * dim, i3 = pl[t + 3] * dim;
      const ax = x[i1] - x[i0], ay = x[i1 + 1] - x[i0 + 1], az = x[i1 + 2] - x[i0 + 2];
      const bx = x[i2] - x[i0], by = x[i2 + 1] - x[i0 + 1], bz = x[i2 + 2] - x[i0 + 2];
      const cx = x[i3] - x[i0], cy = x[i3 + 1] - x[i0 + 1], cz = x[i3 + 2] - x[i0 + 2];
      const bcx = by * cz - bz * cy, bcy = bz * cx - bx * cz, bcz = bx * cy - by * cx;
      const V = ax * bcx + ay * bcy + az * bcz;
      e += wPlanar * V * V;
      const dV = 2 * wPlanar * V;
      const cax = cy * az - cz * ay, cay = cz * ax - cx * az, caz = cx * ay - cy * ax;
      const abx = ay * bz - az * by, aby = az * bx - ax * bz, abz = ax * by - ay * bx;
      g[i1] += dV * bcx; g[i1 + 1] += dV * bcy; g[i1 + 2] += dV * bcz;
      g[i2] += dV * cax; g[i2 + 1] += dV * cay; g[i2 + 2] += dV * caz;
      g[i3] += dV * abx; g[i3 + 1] += dV * aby; g[i3 + 2] += dV * abz;
      g[i0] -= dV * (bcx + cax + abx); g[i0 + 1] -= dV * (bcy + cay + aby); g[i0 + 2] -= dV * (bcz + caz + abz);
    }
    if (dim === 4 && w4 > 0) {
      for (let i = 3; i < x.length; i += 4) {
        e += w4 * x[i] * x[i];
        g[i] += 2 * w4 * x[i];
      }
    }
    return e;
  };
}

/** Embeds one connected, H-explicit molecule in place; returns the number of unsatisfied stereo specs. */
function embedComponent(mol: Mol, seed: number, attempts: number): number {
  const n = mol.atoms.length;
  const ty = typeUFF(mol);
  const bd = buildBounds(mol, ty);
  const T = dgTerms(mol, ty, bd);
  const f3 = dgObjective(T, 3, 0);
  let best: { x: Float64Array; bad: number; err: number } | null = null;

  for (let attempt = 0; attempt < attempts; attempt++) {
    const rng = makeRng(seed * 1009 + attempt * 7919 + 17);
    const x4 = initialCoords(bd, rng);
    lbfgs(dgObjective(T, 4, 0), x4, { dim: 4, maxIter: 1000, gradTol: 1e-4, maxStep: 0.5 });
    lbfgs(dgObjective(T, 4, 0.2), x4, { dim: 4, maxIter: 400, gradTol: 1e-4, maxStep: 0.5 });
    lbfgs(dgObjective(T, 4, 5), x4, { dim: 4, maxIter: 400, gradTol: 1e-4, maxStep: 0.5 });
    const x3 = new Float64Array(n * 3);
    for (let i = 0; i < n; i++) for (let c = 0; c < 3; c++) x3[3 * i + c] = x4[4 * i + c];
    lbfgs(f3, x3, { dim: 3, maxIter: 600, gradTol: 1e-4, maxStep: 0.5 });
    write(mol, x3);

    // verify stereo; repair locally if needed and re-refine
    let bad = countBad(mol);
    if (bad > 0) {
      for (const t of mol.tetra) if (!tetraSatisfied(mol, t)) invertCenter(mol, t);
      for (const d of mol.dbStereo) if (!dbSatisfied(mol, d)) flipDoubleBond(mol, d);
      read(mol, x3);
      lbfgs(f3, x3, { dim: 3, maxIter: 600, gradTol: 1e-4, maxStep: 0.5 });
      write(mol, x3);
      bad = countBad(mol);
    }
    const err = maxLocalViolation(bd, x3);
    if (!best || bad < best.bad || (bad === best.bad && err < best.err)) best = { x: x3.slice(), bad, err };
    if (bad === 0 && err < 0.25) break;
  }
  write(mol, best!.x);
  return best!.bad;
}

/** Largest violation (Å) of the 1-2, 1-3 and 1-4 distance bounds — detects tangled embeddings. */
function maxLocalViolation(bd: Bounds, x: Float64Array): number {
  const { n, lo, up, topo } = bd;
  let worst = 0;
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const t = topo[i * n + j];
      if (t > 3) continue;
      const d = Math.hypot(x[3 * i] - x[3 * j], x[3 * i + 1] - x[3 * j + 1], x[3 * i + 2] - x[3 * j + 2]);
      worst = Math.max(worst, d - up[i * n + j], lo[i * n + j] - d);
    }
  }
  return worst;
}

function countBad(mol: Mol): number {
  let bad = 0;
  for (const t of mol.tetra) if (!tetraSatisfied(mol, t) || !centreInside(mol, t)) bad++;
  for (const d of mol.dbStereo) if (!dbSatisfied(mol, d)) bad++;
  return bad;
}

function write(mol: Mol, x: Float64Array): void {
  mol.atoms.forEach((a, i) => {
    a.x = x[3 * i];
    a.y = x[3 * i + 1];
    a.z = x[3 * i + 2];
  });
}

function read(mol: Mol, x: Float64Array): void {
  mol.atoms.forEach((a, i) => {
    x[3 * i] = a.x;
    x[3 * i + 1] = a.y;
    x[3 * i + 2] = a.z ?? 0;
  });
}
