// Clash detection and resolution for 2D depictions.
// Clashes are non-bonded atom pairs that are too close and bonds that cross. They are removed by
// mirroring or rotating the smaller side of acyclic single bonds on the path between the clashing
// atoms (greedy, best improvement first), followed by a gentle relaxation if anything severe remains.
import { Mol } from '../mol';
import { RingInfo } from '../rings';
import { Grid, segmentsCross } from './geom';
import { relax, DistCon, pairKey } from './refine';

export interface OverlapOptions {
  /** Atoms of the component(s) being resolved. */
  atoms: number[];
  /** movable[i] false → atom i never moves. */
  movable: ArrayLike<boolean | number>;
  /** Returns false if the new coordinates of `moved` atoms violate double-bond stereo. */
  stereoOK: (moved: number[]) => boolean;
  /** Ring systems (atom lists) kept rigid during the final relaxation. */
  ringSystems?: number[][];
}

const CLASH = 0.9; // non-bonded pairs closer than this are penalised
const TRIGGER = 0.75; // …and closer than this trigger a resolution attempt
const SEVERE = 0.55;

function pen(d: number): number {
  if (d >= CLASH) return 0;
  const e = CLASH - d;
  return e * e * (d < SEVERE ? 12 : 3);
}

export function resolveOverlaps(mol: Mol, xy: Float64Array, info: RingInfo, opt: OverlapOptions): void {
  const n = mol.atoms.length;
  const atoms = opt.atoms;
  if (atoms.length < 4) return;
  const inSet = new Uint8Array(n);
  for (const a of atoms) inSet[a] = 1;
  const compBonds: number[] = [];
  mol.bonds.forEach((b, bi) => { if (inSet[b.a] && inSet[b.b]) compBonds.push(bi); });
  const stereoBonds = new Set(mol.dbStereo.map((d) => d.bond));
  const rotatable = (bi: number) => {
    const b = mol.bonds[bi];
    return !info.bondInRing[bi] && b.order !== 3 && !stereoBonds.has(bi) && b.order !== 0;
  };
  // side cache: key = bond*2 + (0: side of b.a, 1: side of b.b)
  const sideCache = new Map<number, number[]>();
  const side = (bi: number, end: number): number[] => {
    const b = mol.bonds[bi];
    const key = bi * 2 + (end === b.a ? 0 : 1);
    let s = sideCache.get(key);
    if (s) return s;
    const other = end === b.a ? b.b : b.a;
    const seen = new Set<number>([end, other]);
    const st = [end];
    s = [end];
    while (st.length) {
      const v = st.pop()!;
      for (const w of mol.neighbors(v)) {
        if (seen.has(w)) continue;
        seen.add(w); s.push(w); st.push(w);
      }
    }
    sideCache.set(key, s);
    return s;
  };

  const bondedOr13 = (i: number, j: number) => mol.bondBetween(i, j) >= 0;

  const P = (i: number) => ({ x: xy[2 * i], y: xy[2 * i + 1] });

  // Spatial indexes of the current (static) state, rebuilt once per round.
  let atomGrid: Grid;
  let midXY = new Float64Array(compBonds.length * 2);
  let bondGrid: Grid;
  const buildGrids = () => {
    atomGrid = new Grid(xy, 1.0, atoms);
    midXY = new Float64Array(compBonds.length * 2);
    compBonds.forEach((bi, k) => {
      const b = mol.bonds[bi];
      midXY[2 * k] = (xy[2 * b.a] + xy[2 * b.b]) / 2;
      midXY[2 * k + 1] = (xy[2 * b.a + 1] + xy[2 * b.b + 1]) / 2;
    });
    bondGrid = new Grid(midXY, 1.6, compBonds.map((_, k) => k));
  };
  const bondsOf = new Map<number, number[]>(); // atom -> indices into compBonds
  compBonds.forEach((bi, k) => {
    const b = mol.bonds[bi];
    for (const a of [b.a, b.b]) {
      let l = bondsOf.get(a);
      if (!l) bondsOf.set(a, (l = []));
      l.push(k);
    }
  });

  /** Penalty of moved atoms (set M, at their current coordinates) against the static rest. */
  const localScore = (moved: number[], mset: Set<number>): number => {
    let s = 0;
    for (const i of moved) {
      atomGrid.near(xy[2 * i], xy[2 * i + 1], CLASH, (j, d2) => {
        if (mset.has(j) || j === i || bondedOr13(i, j)) return;
        s += pen(Math.sqrt(d2));
      });
    }
    // crossings between moved bonds and static bonds (static midpoints from the grid)
    const seen = new Set<number>();
    for (const i of moved) {
      for (const k of bondsOf.get(i) ?? []) {
        if (seen.has(k)) continue;
        seen.add(k);
        const b = mol.bonds[compBonds[k]];
        const pa = P(b.a), pb = P(b.b);
        bondGrid.near((pa.x + pb.x) / 2, (pa.y + pb.y) / 2, 1.6, (k2) => {
          const c = mol.bonds[compBonds[k2]];
          if (mset.has(c.a) || mset.has(c.b)) return;
          if (c.a === b.a || c.a === b.b || c.b === b.a || c.b === b.b) return;
          if (segmentsCross(pa, pb, P(c.a), P(c.b))) s += 4;
        });
      }
    }
    return s;
  };

  const findClashes = (): { i: number; j: number; sev: number }[] => {
    const out: { i: number; j: number; sev: number }[] = [];
    buildGrids();
    for (const i of atoms) {
      atomGrid.near(xy[2 * i], xy[2 * i + 1], CLASH, (j, d2) => {
        if (j <= i || bondedOr13(i, j)) return;
        out.push({ i, j, sev: pen(Math.sqrt(d2)) });
      });
    }
    // crossings (bond midpoints within ~1.6 of each other)
    compBonds.forEach((bi, k) => {
      const b = mol.bonds[bi];
      bondGrid.near(midXY[2 * k], midXY[2 * k + 1], 1.6, (k2) => {
        if (k2 <= k) return;
        const c = mol.bonds[compBonds[k2]];
        if (c.a === b.a || c.a === b.b || c.b === b.a || c.b === b.b) return;
        if (segmentsCross(P(b.a), P(b.b), P(c.a), P(c.b))) out.push({ i: b.a, j: c.a, sev: 4 });
      });
    });
    return out.sort((p, q) => q.sev - p.sev);
  };

  const pathBonds = (s: number, t: number): number[] => {
    const prev = new Map<number, number>();
    prev.set(s, -1);
    const q = [s];
    for (let k = 0; k < q.length && !prev.has(t); k++) {
      const v = q[k];
      for (const bi of mol.adj[v]) {
        const w = mol.other(bi, v);
        if (prev.has(w) || !inSet[w]) continue;
        prev.set(w, bi);
        q.push(w);
      }
    }
    if (!prev.has(t)) return [];
    const out: number[] = [];
    let c = t;
    while (c !== s) {
      const bi = prev.get(c)!;
      out.push(bi);
      c = mol.other(bi, c);
    }
    return out;
  };

  const apply = (moved: number[], f: (x: number, y: number) => [number, number]) => {
    for (const i of moved) {
      const [x, y] = f(xy[2 * i], xy[2 * i + 1]);
      xy[2 * i] = x; xy[2 * i + 1] = y;
    }
  };

  const tried = new Set<string>();
  // deterministic work budget (atoms moved over all trial moves) keeps huge, crowded molecules fast
  let budget = Math.max(80000, 120 * atoms.length);
  for (let round = 0; round < 40 && budget > 0; round++) {
    const clashes = findClashes().filter((c) => c.sev >= pen(TRIGGER));
    if (!clashes.length) break;
    let bestGain = 0.01, bestMove: { moved: number[]; to: Float64Array } | null = null;
    const cands = new Set<number>();
    for (const cl of clashes.slice(0, 8)) for (const bi of pathBonds(cl.i, cl.j)) if (rotatable(bi)) cands.add(bi);
    {
      for (const bi of cands) {
        if (budget <= 0) break;
        const b = mol.bonds[bi];
        // move the smaller side (the other one only if the smaller contains fixed atoms)
        const sA = side(bi, b.a), sB = side(bi, b.b);
        const order: [number, number][] = sA.length <= sB.length ? [[b.b, b.a], [b.a, b.b]] : [[b.a, b.b], [b.b, b.a]];
        let usedSide = false;
        for (const [pivot, end] of order) {
          if (usedSide) break;
          const moved = side(bi, end);
          if (moved.length > atoms.length - 1) continue;
          if (moved.some((i) => !opt.movable[i])) continue;
          // the larger side only moves when the smaller one is anchored, and never when huge
          if (moved.length * 2 > atoms.length + 2 && order[0][1] !== end && moved.length > 12) continue;
          usedSide = true;
          const mset = new Set(moved);
          const before = localScore(moved, mset);
          if (before <= 0) continue;
          const save = new Float64Array(moved.length * 2);
          moved.forEach((i, k) => { save[2 * k] = xy[2 * i]; save[2 * k + 1] = xy[2 * i + 1]; });
          const px = xy[2 * pivot], py = xy[2 * pivot + 1];
          const ex = xy[2 * end], ey = xy[2 * end + 1];
          const moves: { f: (x: number, y: number) => [number, number]; cost: number }[] = [];
          // mirror across the bond axis (pointless for a single terminal atom)
          if (moved.length > 1) {
            const dx = ex - px, dy = ey - py;
            const l2 = dx * dx + dy * dy || 1;
            moves.push({
              f: (x, y) => {
                const t = ((x - px) * dx + (y - py) * dy) / l2;
                const fx = px + t * dx, fy = py + t * dy;
                return [2 * fx - x, 2 * fy - y];
              },
              cost: 0,
            });
          }
          for (const deg of [30, -30, 60, -60, 90, -90]) {
            const a = (deg * Math.PI) / 180, c = Math.cos(a), s = Math.sin(a);
            moves.push({ f: (x, y) => [px + c * (x - px) - s * (y - py), py + s * (x - px) + c * (y - py)], cost: 0.15 * Math.abs(deg) / 30 });
          }
          budget -= moved.length * (moves.length + 1);
          for (const mv of moves) {
            apply(moved, mv.f);
            const ok = opt.stereoOK(moved);
            const after = ok ? localScore(moved, mset) + mv.cost : Infinity;
            const gain = before - after - moved.length * 1e-4;
            if (gain > bestGain) {
              const to = new Float64Array(moved.length * 2);
              moved.forEach((i, k) => { to[2 * k] = xy[2 * i]; to[2 * k + 1] = xy[2 * i + 1]; });
              const sig = moved[0] + ':' + to.slice(0, 4).join(',');
              if (!tried.has(sig)) { bestGain = gain; bestMove = { moved, to }; }
            }
            moved.forEach((i, k) => { xy[2 * i] = save[2 * k]; xy[2 * i + 1] = save[2 * k + 1]; });
          }
        }
      }
    }
    if (!bestMove) break;
    bestMove.moved.forEach((i, k) => { xy[2 * i] = bestMove!.to[2 * k]; xy[2 * i + 1] = bestMove!.to[2 * k + 1]; });
    tried.add(bestMove.moved[0] + ':' + bestMove.to.slice(0, 4).join(','));
  }

  // remaining severe contacts: relax angles a little, keeping ring systems rigid
  const severe = findClashes().filter((c) => c.sev >= pen(SEVERE));
  if (severe.length) relaxComponent(mol, xy, atoms, opt, 60);
}

/** Gentle relaxation: unit bonds, current angles and rigid ring systems, plus repulsion. */
export function relaxComponent(mol: Mol, xy: Float64Array, atoms: number[], opt: OverlapOptions, iterations: number): void {
  const n = mol.atoms.length;
  const inSet = new Uint8Array(n);
  for (const a of atoms) inSet[a] = 1;
  const cons: DistCon[] = [];
  const excl = new Set<number>();
  const d = (i: number, j: number) => Math.hypot(xy[2 * i] - xy[2 * j], xy[2 * i + 1] - xy[2 * j + 1]);
  for (const b of mol.bonds) {
    if (!inSet[b.a] || !inSet[b.b]) continue;
    cons.push({ i: b.a, j: b.b, d: 1, w: 0.8 });
    excl.add(pairKey(b.a, b.b, n));
  }
  for (const a of atoms) {
    const nb = mol.neighbors(a);
    for (let i = 0; i < nb.length; i++) {
      for (let j = i + 1; j < nb.length; j++) {
        // keep current angles, but never below ~85° (collapsed angles are what we are fixing)
        cons.push({ i: nb[i], j: nb[j], d: Math.max(d(nb[i], nb[j]), 1.35), w: 0.25 });
        excl.add(pairKey(nb[i], nb[j], n));
      }
    }
  }
  for (const sys of opt.ringSystems ?? []) {
    if (sys.length > 40) continue;
    for (let i = 0; i < sys.length; i++) {
      for (let j = i + 1; j < sys.length; j++) {
        cons.push({ i: sys[i], j: sys[j], d: d(sys[i], sys[j]), w: 0.3 });
        excl.add(pairKey(sys[i], sys[j], n));
      }
    }
  }
  relax(xy, n, cons, { iterations, movable: opt.movable, atoms, repelDist: 1.0, repelWeight: 0.5, exclude: excl });
}
