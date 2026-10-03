// Perspective-style depiction for cages: embed the ring system in 3D by stress majorization on
// topological distances, then pick the viewing direction whose projection is easiest to read
// (no coincident atoms, no atoms sitting on bonds, uniform bond lengths, few crossings).
import { Mol } from '../mol';
import { Pt, dist, segmentsCross, pointSegDist } from './geom';

type V3 = [number, number, number];

export function embed3DProjected(mol: Mol, atoms: number[], bonds: number[]): Map<number, Pt> | null {
  const n = atoms.length;
  if (n < 4 || n > 120) return null;
  const local = new Map<number, number>();
  atoms.forEach((a, i) => local.set(a, i));
  const adj: number[][] = atoms.map(() => []);
  const edges: [number, number][] = [];
  for (const bi of bonds) {
    const b = mol.bonds[bi];
    const i = local.get(b.a), j = local.get(b.b);
    if (i === undefined || j === undefined) continue;
    adj[i].push(j); adj[j].push(i);
    edges.push([i, j]);
  }
  // topological distances
  const sp: number[][] = [];
  for (let s = 0; s < n; s++) {
    const d = new Array(n).fill(Infinity);
    d[s] = 0;
    const q = [s];
    for (let k = 0; k < q.length; k++) for (const w of adj[q[k]]) if (d[w] === Infinity) { d[w] = d[q[k]] + 1; q.push(w); }
    sp.push(d);
  }
  const target = (k: number) => (k === 1 ? 1 : k === 2 ? 1.63 : 2.3 + 0.75 * (k - 3));
  const D: number[][] = sp.map((row) => row.map((k) => (k === 0 ? 0 : target(Math.min(k, 30)))));

  // classical MDS initialisation (power iteration with deflation)
  const B: number[][] = [];
  const rowMean = D.map((r) => r.reduce((s, v) => s + v * v, 0) / n);
  const allMean = rowMean.reduce((s, v) => s + v, 0) / n;
  for (let i = 0; i < n; i++) {
    B.push([]);
    for (let j = 0; j < n; j++) B[i].push(-0.5 * (D[i][j] * D[i][j] - rowMean[i] - rowMean[j] + allMean));
  }
  const X: V3[] = atoms.map(() => [0, 0, 0]);
  for (let c = 0; c < 3; c++) {
    let v = atoms.map((_, i) => Math.sin(i * 1.7 + c * 2.3) + 0.01 * (i % 7));
    let lambda = 0;
    for (let it = 0; it < 60; it++) {
      const w = new Array(n).fill(0);
      for (let i = 0; i < n; i++) { let s = 0; for (let j = 0; j < n; j++) s += B[i][j] * v[j]; w[i] = s; }
      const norm = Math.hypot(...w) || 1;
      lambda = norm;
      v = w.map((x) => x / norm);
    }
    // Rayleigh quotient for sign of eigenvalue
    let num = 0;
    for (let i = 0; i < n; i++) { let s = 0; for (let j = 0; j < n; j++) s += B[i][j] * v[j]; num += v[i] * s; }
    lambda = num;
    const sc = Math.sqrt(Math.max(lambda, 1e-6));
    for (let i = 0; i < n; i++) X[i][c] = v[i] * sc;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) B[i][j] -= lambda * v[i] * v[j];
  }
  // stress majorization (localized updates), weights 1/d²
  for (let it = 0; it < 150; it++) {
    for (let i = 0; i < n; i++) {
      let sx = 0, sy = 0, sz = 0, sw = 0;
      for (let j = 0; j < n; j++) {
        if (i === j) continue;
        const dij = D[i][j];
        const w = 1 / (dij * dij);
        const dx = X[i][0] - X[j][0], dy = X[i][1] - X[j][1], dz = X[i][2] - X[j][2];
        const r = Math.hypot(dx, dy, dz) || 1e-6;
        sx += w * (X[j][0] + (dij * dx) / r);
        sy += w * (X[j][1] + (dij * dy) / r);
        sz += w * (X[j][2] + (dij * dz) / r);
        sw += w;
      }
      X[i] = [sx / sw, sy / sw, sz / sw];
    }
  }

  // choose the viewing direction
  const bonded = new Set(edges.map(([i, j]) => (i < j ? i * n + j : j * n + i)));
  let best: Pt[] | null = null, bestS = Infinity;
  const NDIR = n > 40 ? 120 : 240;
  for (let k = 0; k < NDIR; k++) {
    // Fibonacci hemisphere
    const z = 1 - (k + 0.5) / NDIR;
    const r = Math.sqrt(1 - z * z);
    const phi = k * 2.399963229728653;
    const v: V3 = [r * Math.cos(phi), r * Math.sin(phi), z];
    const e1 = normalize(Math.abs(v[0]) < 0.9 ? cross(v, [1, 0, 0]) : cross(v, [0, 1, 0]));
    const e2 = cross(v, e1);
    const P: Pt[] = X.map((p) => ({ x: dot(p, e1), y: dot(p, e2) }));
    const s = projectionScore(P, edges, bonded, n);
    if (s < bestS) { bestS = s; best = P; }
  }
  if (!best) return null;
  const lens = edges.map(([i, j]) => dist(best![i], best![j])).sort((a, b) => a - b);
  const med = lens[Math.floor(lens.length / 2)] || 1;
  const out = new Map<number, Pt>();
  atoms.forEach((a, i) => out.set(a, { x: best![i].x / med, y: best![i].y / med }));
  return out;
}

function projectionScore(P: Pt[], edges: [number, number][], bonded: Set<number>, n: number): number {
  const lens = edges.map(([i, j]) => dist(P[i], P[j]));
  const sorted = [...lens].sort((a, b) => a - b);
  const med = sorted[Math.floor(sorted.length / 2)] || 1;
  const minBond = sorted[0] / med;
  let varSum = 0;
  for (const l of lens) varSum += (l / med - 1) ** 2;
  let minNB = Infinity;
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      if (bonded.has(i * n + j)) continue;
      minNB = Math.min(minNB, dist(P[i], P[j]) / med);
    }
  }
  let minAB = Infinity, crossings = 0;
  for (let e = 0; e < edges.length; e++) {
    const [a, b] = edges[e];
    for (let i = 0; i < n; i++) {
      if (i === a || i === b) continue;
      minAB = Math.min(minAB, pointSegDist(P[i], P[a], P[b]) / med);
    }
    for (let f = e + 1; f < edges.length; f++) {
      const [c, d] = edges[f];
      if (a === c || a === d || b === c || b === d) continue;
      if (segmentsCross(P[a], P[b], P[c], P[d])) crossings++;
    }
  }
  return (
    6 * Math.max(0, 0.7 - minBond) +
    8 * Math.max(0, 0.6 - minNB) +
    6 * Math.max(0, 0.3 - minAB) +
    0.4 * crossings +
    (0.5 * varSum) / lens.length
  );
}

const dot = (a: V3, b: V3) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a: V3, b: V3): V3 => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const normalize = (a: V3): V3 => {
  const l = Math.hypot(a[0], a[1], a[2]) || 1;
  return [a[0] / l, a[1] / l, a[2] / l];
};
