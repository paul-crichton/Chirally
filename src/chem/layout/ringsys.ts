// Layout of one ring system (fused, spiro, bridged rings, macrocycles, cages) in local coordinates.
//
// Strategy (in order):
//  1. single ring → regular polygon (macrocycles with trans ring double bonds get a turn sequence);
//  2. cage/bridged skeleton templates matched by graph isomorphism (norbornane, adamantane, cubane…);
//  3. ring-by-ring polygon fusion (ortho-, peri-fused and spiro systems are exact);
//  4. for bridged systems: arcs for the bridges + constraint relaxation, compared with a
//     perspective projection of a 3D stress embedding; the best scoring candidate wins.
import { Mol, DbSpec } from '../mol';
import { RingInfo } from '../rings';
import { Pt, TAU, dist, polygonRadius, segmentsCross, angleTo, normAngle, cross3 } from './geom';
import { relax, DistCon, pairKey } from './refine';
import { matchCageTemplate } from './cages';
import { embed3DProjected } from './embed3d';

export interface SysLayout {
  pos: Map<number, Pt>;
  /** Penalty: 0 = perfect regular depiction. */
  score: number;
  method: 'polygon' | 'macrocycle' | 'fused' | 'template' | 'bridged' | 'projected';
}

interface SysCtx {
  mol: Mol;
  info: RingInfo;
  rings: number[];
  atoms: number[];
  atomSet: Set<number>;
  bonds: number[];
}

/** Lays out ring system `sysIdx` of `info`. */
export function layoutRingSystem(mol: Mol, info: RingInfo, sysIdx: number): SysLayout {
  const rings = info.systems[sysIdx];
  const atomSet = new Set<number>();
  for (const r of rings) for (const a of info.rings[r]) atomSet.add(a);
  const atoms = [...atomSet].sort((a, b) => a - b);
  const bonds: number[] = [];
  mol.bonds.forEach((b, bi) => {
    if (atomSet.has(b.a) && atomSet.has(b.b)) bonds.push(bi);
  });
  const ctx: SysCtx = { mol, info, rings, atoms, atomSet, bonds };

  if (rings.length === 1) return singleRing(ctx);

  const tmpl = matchCageTemplate(mol, atoms, bonds);
  if (tmpl) return { pos: tmpl, score: 0, method: 'template' };

  // Ring-by-ring fusion, trying a few start rings (best heuristic first).
  const starts = [...rings].sort((p, q) => startRingKey(ctx, q) - startRingKey(ctx, p) || p - q);
  const cands: { pos: Map<number, Pt>; bridged: boolean; score: number }[] = [];
  for (const s of starts.slice(0, 6)) {
    const f = fusedPlacement(ctx, s);
    const score = scoreLayout(ctx, f.pos);
    if (score < 0.05) return { pos: f.pos, score, method: 'fused' };
    cands.push({ ...f, score });
  }
  cands.sort((p, q) => p.score - q.score);
  const best = cands[0];
  if (!best.bridged && best.score < 0.5) return { pos: best.pos, score: best.score, method: 'fused' };

  // Bridged / strained: relax the best fusion candidates; a 3D projection is tried when no
  // planar candidate is acceptable.
  const out: SysLayout[] = [{ pos: best.pos, score: best.score, method: best.bridged ? 'bridged' : 'fused' }];
  for (const c of cands.slice(0, 2)) {
    const p = relaxSystem(ctx, c.pos);
    out.push({ pos: p, score: scoreLayout(ctx, p), method: 'bridged' });
  }
  for (const env of envelopeCandidates(ctx).slice(0, 2)) {
    const p = relaxSystem(ctx, env);
    out.push({ pos: p, score: scoreLayout(ctx, p), method: 'bridged' });
  }
  out.sort((p, q) => p.score - q.score);
  if (out[0].score > 4) {
    const proj = embed3DProjected(mol, atoms, bonds);
    if (proj) {
      const p = relaxSystem(ctx, proj, 0.3);
      // perspective drawings tolerate crossings; mildly favour planar solutions
      out.push({ pos: p, score: scoreLayout(ctx, p) + 0.5, method: 'projected' });
    }
  }
  out.sort((p, q) => p.score - q.score);
  if (debugHook) debugHook(out.map((o) => `${o.method}:${o.score.toFixed(2)}`).join(' ') + ` | fused best ${best.score.toFixed(2)} bridged=${best.bridged}`);
  return out[0];
}

/** Diagnostics hook (tests/tools only). */
export let debugHook: ((msg: string) => void) | null = null;
export function setLayoutDebugHook(f: ((msg: string) => void) | null): void {
  debugHook = f;
}

/** Preference for the first ring: many fused neighbours, then larger rings. */
function startRingKey(ctx: SysCtx, r: number): number {
  const ring = ctx.info.rings[r];
  let nb = 0;
  for (const q of ctx.rings) {
    if (q === r) continue;
    let shared = 0;
    for (const a of ctx.info.rings[q]) if (ring.includes(a)) shared++;
    if (shared >= 2) nb++;
  }
  return nb * 100 + Math.min(ring.length, 12);
}

// ───────────────────────────── single rings ─────────────────────────────

function regularPolygon(ring: number[], center: Pt, startAngle: number, dir: 1 | -1, edge = 1): Map<number, Pt> {
  const n = ring.length;
  const R = polygonRadius(n) * edge;
  const m = new Map<number, Pt>();
  ring.forEach((a, k) => {
    const t = startAngle + (dir * k * TAU) / n;
    m.set(a, { x: center.x + R * Math.cos(t), y: center.y + R * Math.sin(t) });
  });
  return m;
}

function singleRing(ctx: SysCtx): SysLayout {
  const ring = ctx.info.rings[ctx.rings[0]];
  const n = ring.length;
  if (n >= 8) {
    const turns = macrocycleTurns(ctx.mol, ring);
    if (turns) {
      const pos = turnPolygon(ctx, ring, turns);
      return { pos, score: scoreLayout(ctx, pos), method: 'macrocycle' };
    }
  }
  // flat bottom edge for odd rings, vertical side edges for even rings (rotated later anyway)
  const start = Math.PI / 2 + Math.PI / n;
  return { pos: regularPolygon(ring, { x: 0, y: 0 }, start, 1), score: 0, method: 'polygon' };
}

/**
 * For rings ≥ 8 with stereo double bonds that must be drawn trans, returns a turn sequence
 * (+1 convex, −1 reflex) per ring atom satisfying every cis/trans spec; null if a plain polygon works.
 */
function macrocycleTurns(mol: Mol, ring: number[]): number[] | null {
  const n = ring.length;
  // relation between consecutive turns: 0 free, 1 equal (cis), -1 opposite (trans)
  const rel = new Array(n).fill(0);
  let anyTrans = false;
  for (let k = 0; k < n; k++) {
    const a = ring[k], b = ring[(k + 1) % n];
    const bi = mol.bondBetween(a, b);
    const spec = mol.dbStereo.find((d) => d.bond === bi);
    if (!spec) continue;
    const cis = ringCis(mol, spec, ring[(k - 1 + n) % n], a, b, ring[(k + 2) % n]);
    if (cis === null) continue;
    rel[k] = cis ? 1 : -1;
    if (!cis) anyTrans = true;
  }
  if (!anyTrans) return null;
  // Segments joined by constrained bonds are solved together; free bonds decouple them.
  const turns = new Array(n).fill(0);
  let startK = rel.findIndex((r) => r === 0);
  if (startK < 0) startK = 0; // fully constrained ring: parity decides
  for (let s = 0; s < n; s++) {
    const k0 = (startK + 1 + s) % n; // atom index after a free bond
    if (turns[k0] !== 0) continue;
    // walk the segment
    const seg: number[] = [k0];
    const val: number[] = [1];
    let k = k0;
    while (rel[k] !== 0 && seg.length < n) {
      const nk = (k + 1) % n;
      if (nk === k0) break;
      val.push(val[val.length - 1] * rel[k]);
      seg.push(nk);
      k = nk;
    }
    const plus = val.filter((v) => v > 0).length;
    const sign = plus * 2 >= val.length ? 1 : -1;
    seg.forEach((idx, j) => (turns[idx] = val[j] * sign));
  }
  return turns.map((t) => (t === 0 ? 1 : t));
}

/** Whether the ring neighbours p (of a) and q (of b) are cis across double bond a=b, from spec d. */
function ringCis(mol: Mol, d: DbSpec, p: number, a: number, b: number, q: number): boolean | null {
  const bond = mol.bonds[d.bond];
  // orient spec refs: ra is the ref on a's side
  let ra = d.a, rb = d.b;
  if (bond.a !== a) { ra = d.b; rb = d.a; }
  if (!mol.neighbors(a).includes(ra) || !mol.neighbors(b).includes(rb)) {
    [ra, rb] = [rb, ra];
    if (!mol.neighbors(a).includes(ra) || !mol.neighbors(b).includes(rb)) return null;
  }
  let cis = d.cis;
  if (ra !== p) cis = !cis;
  if (rb !== q) cis = !cis;
  return cis;
}

/** Closed polygon following the turn sequence (equal edges), relaxed to close exactly. */
function turnPolygon(ctx: SysCtx, ring: number[], turns: number[]): Map<number, Pt> {
  const n = ring.length;
  const P = turns.filter((t) => t > 0).length;
  const M = n - P;
  const reflex = Math.PI / 3;
  let convex = (TAU + M * reflex) / P;
  if (convex > (2 * Math.PI) / 3) convex = (2 * Math.PI) / 3;
  const pts: Pt[] = [];
  let x = 0, y = 0, h = 0;
  for (let k = 0; k < n; k++) {
    pts.push({ x, y });
    const t = turns[(k + 1) % n] > 0 ? convex : -reflex;
    x += Math.cos(h); y += Math.sin(h);
    h += t;
  }
  // distribute closure error linearly
  const ex = x, ey = y;
  pts.forEach((p, k) => { p.x -= (ex * k) / n; p.y -= (ey * k) / n; });
  const pos = new Map<number, Pt>();
  ring.forEach((a, k) => pos.set(a, pts[k]));
  // relax with exact 1-2 and turn-dependent 1-3 targets
  const idx = ring;
  const xy = new Float64Array(ctx.mol.atoms.length * 2);
  for (const a of idx) { const p = pos.get(a)!; xy[2 * a] = p.x; xy[2 * a + 1] = p.y; }
  const cons: DistCon[] = [];
  for (let k = 0; k < n; k++) {
    const a = ring[k], b = ring[(k + 1) % n], c = ring[(k + 2) % n];
    cons.push({ i: a, j: b, d: 1, w: 1 });
    const turn = turns[(k + 1) % n] > 0 ? convex : reflex;
    const interior = Math.PI - turn;
    cons.push({ i: a, j: c, d: 2 * Math.sin(interior / 2), w: 0.5 });
  }
  const nAll = ctx.mol.atoms.length;
  const excl = new Set<number>();
  for (const c of cons) excl.add(pairKey(c.i, c.j, nAll));
  relax(xy, nAll, cons, { iterations: 300, atoms: idx, repelDist: 1.2, repelWeight: 0.3, exclude: excl });
  for (const a of idx) pos.set(a, { x: xy[2 * a], y: xy[2 * a + 1] });
  return pos;
}

// ───────────────────────────── fused systems ─────────────────────────────

interface FusedResult {
  pos: Map<number, Pt>;
  bridged: boolean;
}

function fusedPlacement(ctx: SysCtx, start: number): FusedResult {
  const { info } = ctx;
  const pos = new Map<number, Pt>();
  const startRing = info.rings[start];
  const n0 = startRing.length;
  regularPolygon(startRing, { x: 0, y: 0 }, Math.PI / 2 + Math.PI / n0, 1).forEach((p, a) => pos.set(a, p));
  const done = new Set<number>([start]);
  let bridged = false;
  while (done.size < ctx.rings.length) {
    // next ring: most placed atoms, preferring one contiguous shared run
    let bestR = -1, bestKey = -Infinity;
    for (const r of ctx.rings) {
      if (done.has(r)) continue;
      const ring = info.rings[r];
      const cnt = ring.filter((a) => pos.has(a)).length;
      if (!cnt) continue;
      const runs = placedRuns(ring, pos);
      const key = cnt * 10 + (runs.length === 1 ? 5 : 0) + (cnt === 2 ? 2 : 0) - ring.length * 0.01;
      if (key > bestKey) { bestKey = key; bestR = r; }
    }
    if (bestR < 0) break; // disconnected (should not happen within a system)
    done.add(bestR);
    const ring = info.rings[bestR];
    if (ring.every((a) => pos.has(a))) continue;
    if (placeRing(ctx, ring, pos)) bridged = true;
  }
  return { pos, bridged };
}

/** Maximal cyclic runs of placed atoms, as index lists into `ring`. */
function placedRuns(ring: number[], pos: Map<number, Pt>): number[][] {
  const n = ring.length;
  const placed = ring.map((a) => pos.has(a));
  if (placed.every(Boolean)) return [ring.map((_, i) => i)];
  // start just after an unplaced atom
  let s = placed.findIndex((p) => !p);
  const runs: number[][] = [];
  let cur: number[] = [];
  for (let k = 1; k <= n; k++) {
    const i = (s + k) % n;
    if (placed[i]) cur.push(i);
    else if (cur.length) { runs.push(cur); cur = []; }
  }
  if (cur.length) runs.push(cur);
  return runs;
}

/** Places the unplaced atoms of `ring`; returns true if a non-regular (bridged) placement was needed. */
function placeRing(ctx: SysCtx, ring: number[], pos: Map<number, Pt>): boolean {
  const n = ring.length;
  const runs = placedRuns(ring, pos);
  const placedCount = runs.reduce((s, r) => s + r.length, 0);

  if (placedCount === 1) {
    // spiro: polygon in the free exterior direction of the shared atom
    const s = ring[runs[0][0]];
    const ps = pos.get(s)!;
    const nb = ctx.mol.neighbors(s).filter((j) => pos.has(j));
    const dirAng = freeDirection(ps, nb.map((j) => pos.get(j)!));
    const R = polygonRadius(n);
    const c = { x: ps.x + R * Math.cos(dirAng), y: ps.y + R * Math.sin(dirAng) };
    const k0 = runs[0][0];
    const order = ring.map((_, k) => ring[(k0 + k) % n]);
    const cands = [1, -1].map((d) => regularPolygon(order, c, dirAng + Math.PI, d as 1 | -1));
    const best = pickCandidate(ctx, cands, pos, order.slice(1));
    best.forEach((p, a) => { if (!pos.has(a)) pos.set(a, p); });
    return false;
  }

  if (runs.length === 1) {
    const run = runs[0].map((i) => ring[i]); // e2 … e1 in ring order
    const free: number[] = [];
    for (let k = 1; k <= n - run.length; k++) free.push(ring[(runs[0][runs[0].length - 1] + k) % n]);
    // regular polygon candidates consistent with the run
    const regs = regularFromRun(ring, runs[0], pos);
    if (regs.length) {
      const best = pickCandidate(ctx, regs, pos, free);
      if (candidateClash(ctx, best, pos, free) < 1) {
        best.forEach((p, a) => { if (!pos.has(a)) pos.set(a, p); });
        return false;
      }
    }
    placeArc(ctx, run[run.length - 1], run[0], free, pos);
    return true;
  }

  // several separate placed segments: bridge each gap with an arc
  for (let r = 0; r < runs.length; r++) {
    const cur = runs[r], next = runs[(r + 1) % runs.length];
    const e1 = ring[cur[cur.length - 1]], e2 = ring[next[0]];
    const free: number[] = [];
    for (let k = cur[cur.length - 1] + 1; ; k++) {
      const i = k % n;
      if (i === next[0]) break;
      free.push(ring[i]);
    }
    if (free.length) placeArc(ctx, e1, e2, free, pos);
  }
  return true;
}

/** Direction (angle) of the widest gap among neighbour positions around p. */
export function freeDirection(p: Pt, nbrs: Pt[]): number {
  if (!nbrs.length) return 0;
  if (nbrs.length === 1) return angleTo(nbrs[0], p);
  const angs = nbrs.map((q) => normAngle(angleTo(p, q))).sort((a, b) => a - b);
  let best = 0, bestGap = -1;
  for (let k = 0; k < angs.length; k++) {
    const a = angs[k], b = k + 1 < angs.length ? angs[k + 1] : angs[0] + TAU;
    if (b - a > bestGap) { bestGap = b - a; best = a + (b - a) / 2; }
  }
  return best;
}

/** Regular polygons through the placed run (both orientations), keeping only consistent ones. */
function regularFromRun(ring: number[], runIdx: number[], pos: Map<number, Pt>): Map<number, Pt>[] {
  const n = ring.length;
  const a0 = ring[runIdx[0]], a1 = ring[runIdx[1]];
  const p0 = pos.get(a0)!, p1 = pos.get(a1)!;
  const e = dist(p0, p1);
  if (e < 1e-6) return [];
  const out: Map<number, Pt>[] = [];
  const order = ring.map((_, k) => ring[(runIdx[0] + k) % n]); // run first, in ring order
  const R = polygonRadius(n) * e;
  const apo = (e / 2) / Math.tan(Math.PI / n);
  const mx = (p0.x + p1.x) / 2, my = (p0.y + p1.y) / 2;
  const nx = -(p1.y - p0.y) / e, ny = (p1.x - p0.x) / e;
  for (const side of [1, -1]) {
    const c = { x: mx + side * nx * apo, y: my + side * ny * apo };
    const t0 = angleTo(c, p0);
    let step = normAngle(angleTo(c, p1) - t0);
    if (step > Math.PI) step -= TAU;
    const m = new Map<number, Pt>();
    order.forEach((a, k) => m.set(a, { x: c.x + R * Math.cos(t0 + k * step), y: c.y + R * Math.sin(t0 + k * step) }));
    let ok = true;
    for (let k = 2; k < runIdx.length; k++) {
      if (dist(m.get(order[k])!, pos.get(order[k])!) > 0.2 * e) { ok = false; break; }
    }
    if (ok) out.push(m);
  }
  return out;
}

/** Places `free` atoms on a circular arc of unit chords from e1 to e2 (the side with fewer clashes). */
function placeArc(ctx: SysCtx, e1: number, e2: number, free: number[], pos: Map<number, Pt>): void {
  const P = pos.get(e1)!, Q = pos.get(e2)!;
  const N = free.length + 1;
  const L = dist(P, Q);
  const cands: Map<number, Pt>[] = [];
  if (L >= N - 1e-6 || L < 1e-6) {
    const m = new Map<number, Pt>();
    const dx = L < 1e-6 ? 1 : (Q.x - P.x), dy = L < 1e-6 ? 0 : (Q.y - P.y);
    free.forEach((a, k) => m.set(a, { x: P.x + (dx * (k + 1)) / N, y: P.y + (dy * (k + 1)) / N }));
    cands.push(m);
  } else {
    // solve sin(Nα/2)/sin(α/2) = L for α in (0, 2π/N)
    let lo = 1e-9, hi = TAU / N - 1e-9;
    for (let it = 0; it < 60; it++) {
      const mid = (lo + hi) / 2;
      const v = Math.sin((N * mid) / 2) / Math.sin(mid / 2);
      if (v > L) lo = mid; else hi = mid;
    }
    const alpha = (lo + hi) / 2;
    const rho = 1 / (2 * Math.sin(alpha / 2));
    const mx = (P.x + Q.x) / 2, my = (P.y + Q.y) / 2;
    const ux = (Q.x - P.x) / L, uy = (Q.y - P.y) / L;
    for (const side of [1, -1]) {
      const nx = -uy * side, ny = ux * side; // bulge direction
      const h = rho * Math.cos((N * alpha) / 2);
      const c = { x: mx - nx * h, y: my - ny * h };
      const t0 = angleTo(c, P);
      // choose the rotation sense whose arc midpoint lies on the bulge side
      const tm1 = t0 + (N * alpha) / 2, tm2 = t0 - (N * alpha) / 2;
      const s1 = (c.x + rho * Math.cos(tm1) - mx) * nx + (c.y + rho * Math.sin(tm1) - my) * ny;
      const s2 = (c.x + rho * Math.cos(tm2) - mx) * nx + (c.y + rho * Math.sin(tm2) - my) * ny;
      const dir = s1 >= s2 ? 1 : -1;
      const m = new Map<number, Pt>();
      free.forEach((a, k) => {
        const t = t0 + dir * (k + 1) * alpha;
        m.set(a, { x: c.x + rho * Math.cos(t), y: c.y + rho * Math.sin(t) });
      });
      cands.push(m);
    }
  }
  const best = pickCandidate(ctx, cands, pos, free, [e1, e2]);
  best.forEach((p, a) => { if (!pos.has(a)) pos.set(a, p); });
}

/** Clash/crossing penalty of candidate positions (for `fresh` atoms) against already placed atoms. */
function candidateClash(ctx: SysCtx, cand: Map<number, Pt>, pos: Map<number, Pt>, fresh: number[], ends: number[] = []): number {
  const { mol } = ctx;
  let s = 0;
  for (const a of fresh) {
    const p = cand.get(a)!;
    for (const [b, q] of pos) {
      if (b === a) continue;
      const d = dist(p, q);
      if (d < 1.0 && mol.bondBetween(a, b) < 0) s += (1.0 - d) * (1.0 - d) * (d < 0.5 ? 20 : 4);
    }
  }
  // crossings of new bonds with placed bonds
  const get = (a: number) => cand.get(a) ?? pos.get(a);
  const freshSet = new Set(fresh);
  const newBonds: [number, number][] = [];
  const chain = ends.length === 2 ? [ends[0], ...fresh, ends[1]] : fresh;
  for (const bi of ctx.bonds) {
    const b = mol.bonds[bi];
    if ((freshSet.has(b.a) || freshSet.has(b.b)) && get(b.a) && get(b.b)) newBonds.push([b.a, b.b]);
  }
  void chain;
  for (const bi of ctx.bonds) {
    const b = mol.bonds[bi];
    if (freshSet.has(b.a) || freshSet.has(b.b)) continue;
    const pa = pos.get(b.a), pb = pos.get(b.b);
    if (!pa || !pb) continue;
    for (const [x, y] of newBonds) {
      if (x === b.a || x === b.b || y === b.a || y === b.b) continue;
      if (segmentsCross(get(x)!, get(y)!, pa, pb)) s += 3;
    }
  }
  return s;
}

function pickCandidate(ctx: SysCtx, cands: Map<number, Pt>[], pos: Map<number, Pt>, fresh: number[], ends: number[] = []): Map<number, Pt> {
  // tie-break: farther from the centroid of placed atoms
  let cx = 0, cy = 0;
  for (const p of pos.values()) { cx += p.x; cy += p.y; }
  cx /= pos.size || 1; cy /= pos.size || 1;
  let best = cands[0], bestS = Infinity;
  for (const c of cands) {
    let fx = 0, fy = 0;
    for (const a of fresh) { const p = c.get(a)!; fx += p.x; fy += p.y; }
    fx /= fresh.length || 1; fy /= fresh.length || 1;
    const s = candidateClash(ctx, c, pos, fresh, ends) - 0.01 * Math.hypot(fx - cx, fy - cy);
    if (s < bestS) { bestS = s; best = c; }
  }
  return best;
}

// ───────────────────────────── envelope + Tutte ─────────────────────────────

/**
 * Bridged-system candidates: a long peripheral cycle (found among sums of SSSR rings) drawn as a
 * convex polygon, the remaining atoms placed at the barycentre of their neighbours (Tutte).
 * For planar skeletons whose cycle bounds the outer face this yields crossing-free drawings.
 */
function envelopeCandidates(ctx: SysCtx): Map<number, Pt>[] {
  const { mol, info } = ctx;
  const R = ctx.rings.length;
  if (R < 2 || R > 12) return [];
  const nb = ctx.bonds.length;
  const edgeIdx = new Map<number, number>();
  ctx.bonds.forEach((bi, k) => edgeIdx.set(bi, k));
  const words = Math.ceil(nb / 32);
  const ringVec = ctx.rings.map((r) => {
    const v = new Uint32Array(words);
    for (const bi of info.ringBonds[r]) { const k = edgeIdx.get(bi)!; v[k >>> 5] ^= 1 << (k & 31); }
    return v;
  });
  const cycles = new Map<string, number[]>();
  const v = new Uint32Array(words);
  for (let mask = 1; mask < 1 << R; mask++) {
    if ((mask & (mask - 1)) === 0) continue; // single rings are not envelopes
    v.fill(0);
    for (let r = 0; r < R; r++) if (mask & (1 << r)) for (let w = 0; w < words; w++) v[w] ^= ringVec[r][w];
    const nbr = new Map<number, number[]>();
    const es: number[] = [];
    let ok = true;
    for (let k = 0; k < nb && ok; k++) {
      if (!(v[k >>> 5] & (1 << (k & 31)))) continue;
      es.push(k);
      const b = mol.bonds[ctx.bonds[k]];
      for (const [x, y] of [[b.a, b.b], [b.b, b.a]]) {
        let l = nbr.get(x);
        if (!l) nbr.set(x, (l = []));
        l.push(y);
        if (l.length > 2) ok = false;
      }
    }
    if (!ok || es.length < 6) continue;
    const key = es.join(',');
    if (cycles.has(key)) continue;
    // walk the cycle
    const startA = nbr.keys().next().value as number;
    const order = [startA];
    let prev = -1, cur = startA;
    for (;;) {
      const l = nbr.get(cur)!;
      if (l.length !== 2) { ok = false; break; }
      const nx = l[0] === prev ? l[1] : l[0];
      if (nx === startA) break;
      order.push(nx);
      prev = cur; cur = nx;
      if (order.length > es.length) { ok = false; break; }
    }
    if (!ok || order.length !== es.length) continue;
    cycles.set(key, order);
  }
  const sorted = [...cycles.values()].sort((a, b) => b.length - a.length || a[0] - b[0]);
  const out: Map<number, Pt>[] = [];
  for (const cyc of sorted.slice(0, 3)) {
    const m = cyc.length;
    const xy = new Map<number, Pt>();
    regularPolygon(cyc, { x: 0, y: 0 }, Math.PI / 2, 1).forEach((p, a) => xy.set(a, p));
    const inner = ctx.atoms.filter((a) => !xy.has(a));
    inner.forEach((a, k) => xy.set(a, { x: 0.01 * Math.cos(k * 2.4), y: 0.01 * Math.sin(k * 2.4) }));
    const nbrs = new Map<number, number[]>();
    for (const a of inner) nbrs.set(a, mol.neighbors(a).filter((j) => ctx.atomSet.has(j)));
    for (let it = 0; it < 300; it++) {
      for (const a of inner) {
        const l = nbrs.get(a)!;
        let x = 0, y = 0;
        for (const j of l) { const p = xy.get(j)!; x += p.x; y += p.y; }
        xy.set(a, { x: x / l.length, y: y / l.length });
      }
    }
    // shrink the polygon toward a bond length compatible with the inner bridges
    void m;
    out.push(xy);
  }
  return out;
}

// ───────────────────────────── relaxation & scoring ─────────────────────────────

/**
 * Distance constraints describing an ideal ring system: unit bonds, regular small rings (≤ 8,
 * all pairwise distances), 120° angles at macrocycle atoms that belong to no small ring, and
 * a minimum separation for every pair of neighbours of an atom (keeps angles from collapsing).
 */
export function ringSystemConstraints(mol: Mol, info: RingInfo, rings: number[], bonds: number[], wScale = 1): DistCon[] {
  const cons: DistCon[] = [];
  const n = mol.atoms.length;
  const have = new Set<number>();
  const add = (c: DistCon) => {
    const k = pairKey(c.i, c.j, n);
    if (have.has(k)) return;
    have.add(k);
    cons.push(c);
  };
  for (const bi of bonds) add({ i: mol.bonds[bi].a, j: mol.bonds[bi].b, d: 1, w: 0.9 * wScale });
  const inSmall = new Set<number>();
  for (const r of rings) if (info.rings[r].length <= 8) for (const a of info.rings[r]) inSmall.add(a);
  const sysAtoms = new Set<number>();
  for (const r of rings) for (const a of info.rings[r]) sysAtoms.add(a);
  for (const r of rings) {
    const ring = info.rings[r];
    const m = ring.length;
    if (m <= 8) {
      const R = polygonRadius(m);
      for (let k = 2; k <= Math.floor(m / 2); k++) {
        const d = 2 * R * Math.sin((k * Math.PI) / m);
        for (let i = 0; i < m; i++) {
          if (k * 2 === m && i >= m / 2) break;
          add({ i: ring[i], j: ring[(i + k) % m], d, w: (k === 2 ? 0.5 : 0.25) * wScale });
        }
      }
    } else {
      for (let i = 0; i < m; i++) {
        const c = ring[(i + 1) % m];
        if (inSmall.has(c)) continue;
        add({ i: ring[i], j: ring[(i + 2) % m], d: Math.sqrt(3), w: 0.3 * wScale });
      }
    }
  }
  for (const a of sysAtoms) {
    const nb = mol.neighbors(a).filter((j) => sysAtoms.has(j));
    for (let i = 0; i < nb.length; i++) for (let j = i + 1; j < nb.length; j++) add({ i: nb[i], j: nb[j], d: 1.25, w: 0.5 * wScale, min: true });
  }
  return cons;
}

/** Reflects ring atoms that fold their (small) ring inside-out across the line of their ring neighbours. */
function unfoldRings(ctx: SysCtx, xy: Float64Array): void {
  const { mol, info } = ctx;
  const P = (i: number): Pt => ({ x: xy[2 * i], y: xy[2 * i + 1] });
  for (let pass = 0; pass < 3; pass++) {
    let changed = false;
    for (const r of ctx.rings) {
      const ring = info.rings[r];
      const m = ring.length;
      if (m > 8) continue;
      let area = 0;
      for (let k = 0; k < m; k++) { const p = P(ring[k]), q = P(ring[(k + 1) % m]); area += p.x * q.y - q.x * p.y; }
      for (let k = 0; k < m; k++) {
        const a = ring[(k - 1 + m) % m], c = ring[k], b = ring[(k + 1) % m];
        const t = cross3(P(a), P(c), P(b));
        if (t * area >= 0) continue;
        const sysDeg = mol.neighbors(c).filter((j) => ctx.atomSet.has(j)).length;
        if (sysDeg !== 2) continue;
        // reflect c across line a-b
        const pa = P(a), pb = P(b), pc = P(c);
        const dx = pb.x - pa.x, dy = pb.y - pa.y, l2 = dx * dx + dy * dy || 1;
        const u = ((pc.x - pa.x) * dx + (pc.y - pa.y) * dy) / l2;
        const fx = pa.x + u * dx, fy = pa.y + u * dy;
        xy[2 * c] = 2 * fx - pc.x; xy[2 * c + 1] = 2 * fy - pc.y;
        changed = true;
      }
    }
    if (!changed) break;
  }
}

function relaxSystem(ctx: SysCtx, start: Map<number, Pt>, wScale = 1): Map<number, Pt> {
  const { mol, info } = ctx;
  const nAll = mol.atoms.length;
  const xy = new Float64Array(nAll * 2);
  for (const [a, p] of start) { xy[2 * a] = p.x; xy[2 * a + 1] = p.y; }
  unfoldRings(ctx, xy);
  const cons = ringSystemConstraints(mol, info, ctx.rings, ctx.bonds, wScale);
  const excl = new Set<number>();
  for (const c of cons) excl.add(pairKey(c.i, c.j, nAll));
  relax(xy, nAll, cons, { iterations: 160, atoms: ctx.atoms, repelDist: 0.95, repelWeight: 0.5, exclude: excl });
  const out = new Map<number, Pt>();
  for (const a of ctx.atoms) out.set(a, { x: xy[2 * a], y: xy[2 * a + 1] });
  return out;
}

/** Layout penalty: bond lengths, small-ring regularity/folding, non-bonded clashes, crossings. */
function scoreLayout(ctx: SysCtx, pos: Map<number, Pt>): number {
  const { mol, info } = ctx;
  let s = 0;
  for (const bi of ctx.bonds) {
    const b = mol.bonds[bi];
    const d = dist(pos.get(b.a)!, pos.get(b.b)!);
    s += 2 * (d - 1) * (d - 1);
    if (d < 0.75 || d > 1.3) s += 1;
  }
  for (const r of ctx.rings) {
    const ring = info.rings[r];
    const m = ring.length;
    if (m > 8) continue;
    const ideal = 2 * polygonRadius(m) * Math.sin((2 * Math.PI) / m);
    let area = 0;
    for (let k = 0; k < m; k++) { const p = pos.get(ring[k])!, q = pos.get(ring[(k + 1) % m])!; area += p.x * q.y - q.x * p.y; }
    for (let k = 0; k < m; k++) {
      const a = pos.get(ring[k])!, c = pos.get(ring[(k + 1) % m])!, b = pos.get(ring[(k + 2) % m])!;
      const d = dist(a, b);
      s += 0.5 * (d - ideal) * (d - ideal);
      if (cross3(a, c, b) * area < 0) s += 2; // folded vertex
    }
  }
  const at = ctx.atoms;
  for (let i = 0; i < at.length; i++) {
    for (let j = i + 1; j < at.length; j++) {
      if (mol.bondBetween(at[i], at[j]) >= 0) continue;
      const d = dist(pos.get(at[i])!, pos.get(at[j])!);
      if (d < 0.8) s += (0.8 - d) * 10 + (d < 0.45 ? 20 : 0);
    }
  }
  const bl = ctx.bonds.map((bi) => mol.bonds[bi]);
  for (let i = 0; i < bl.length; i++) {
    for (let j = i + 1; j < bl.length; j++) {
      const p = bl[i], q = bl[j];
      if (p.a === q.a || p.a === q.b || p.b === q.a || p.b === q.b) continue;
      if (segmentsCross(pos.get(p.a)!, pos.get(p.b)!, pos.get(q.a)!, pos.get(q.b)!)) s += 3;
    }
  }
  return s;
}

