// Hand-tuned depictions of common bridged/cage ring skeletons (ChemDraw-like perspective drawings).
// A ring system whose skeleton graph (elements and bond orders ignored) is isomorphic to a template
// takes the template coordinates. Coordinates: y down, bond length ≈ 1.
import { Mol } from '../mol';
import { Pt } from './geom';

interface CageTemplate {
  name: string;
  /** Skeleton as an edge list over atoms 0..n-1. */
  edges: [number, number][];
  xy: [number, number][];
}

const TEMPLATES: CageTemplate[] = [];

/** Registers a template from edges + coordinates (used by the table below). */
function T(name: string, edges: [number, number][], xy: [number, number][]): void {
  TEMPLATES.push({ name, edges, xy });
}

// ── template table ──────────────────────────────────────────────────────────
// Edge lists follow the atom numbering of the xy arrays. Designed from projections of 3D
// conformers, then idealised by hand (bond lengths 0.65–1.2, no atom closer than 0.5 to another).

// bicyclo[2.2.1]heptane (norbornane, camphor, borneol…): 0,1 bridgeheads, 2 = C7, 3-4 back, 5-6 front bridge
T('norbornane',
  [[0, 2], [0, 4], [0, 6], [1, 2], [1, 3], [1, 5], [3, 4], [5, 6]],
  [[-1.0, -0.1], [0.75, -0.1], [-0.15, -0.85], [1.15, 0.6], [0.1, 0.3], [0.3, 1.05], [-0.7, 0.8]]);

// bicyclo[2.2.2]octane (quinuclidine, DABCO, 1,8-cineole): 0,1 bridgeheads; bridges 2-6 (top), 5-7 (middle), 3-4 (bottom)
T('bicyclo[2.2.2]octane',
  [[0, 2], [0, 3], [0, 5], [1, 4], [1, 6], [1, 7], [2, 6], [3, 4], [5, 7]],
  [[-1.2, 0], [1.2, 0], [-0.5, -0.85], [-0.5, 0.8], [0.5, 0.8], [-0.45, 0.2], [0.5, -0.85], [0.45, 0.2]]);

// bicyclo[3.2.1]octane (tropane, cocaine): 1,2 bridgeheads; 0 one-atom bridge; 3-4 two-atom; 5-7-6 three-atom
T('bicyclo[3.2.1]octane',
  [[0, 1], [0, 2], [1, 3], [1, 5], [2, 4], [2, 6], [3, 4], [5, 7], [6, 7]],
  [[0, -0.75], [-0.95, -0.1], [0.95, -0.1], [-0.45, 0.35], [0.45, 0.35], [-0.95, 0.85], [0.95, 0.85], [0, 1.35]]);

// adamantane: 0-3 bridgeheads, 4-9 methylenes
T('adamantane',
  [[0, 4], [0, 5], [0, 6], [1, 4], [1, 7], [1, 8], [2, 5], [2, 7], [2, 9], [3, 6], [3, 8], [3, 9]],
  [[-0.77, 0.84], [0.37, -0.77], [0.99, 0.57], [-0.58, -0.64], [-0.41, 0.07], [0.21, 1.41], [-1.35, 0.2], [1.35, -0.2], [-0.21, -1.41], [0.41, -0.07]]);

// cubane: vertex i = bits (x, y, depth); back face offset up-right (classic perspective cube)
T('cubane',
  [[0, 1], [0, 2], [0, 4], [1, 3], [1, 5], [2, 3], [2, 6], [3, 7], [4, 5], [4, 6], [5, 7], [6, 7]],
  [[-0.78, -0.28], [0.22, -0.28], [-0.78, 0.72], [0.22, 0.72], [-0.23, -0.73], [0.77, -0.73], [-0.23, 0.27], [0.77, 0.27]]);

// bicyclo[3.1.1]heptane (pinane, α/β-pinene): 0,1 bridgeheads; 2 outer and 3 inner one-atom bridges; 4-5-6
T('bicyclo[3.1.1]heptane',
  [[0, 2], [1, 2], [0, 3], [1, 3], [0, 4], [4, 5], [5, 6], [6, 1]],
  [[0.85, -0.6], [-0.85, -0.6], [0, -1.5], [0, -0.95], [0.9, 0.45], [0, 1.0], [-0.9, 0.45]]);

// bicyclo[1.1.1]pentane: 0,1 bridgeheads; 2,3,4 bridges
T('bicyclo[1.1.1]pentane',
  [[0, 2], [0, 3], [0, 4], [1, 2], [1, 3], [1, 4]],
  [[0, -0.75], [0, 0.75], [-0.85, 0], [0.85, 0], [0.3, 0.05]]);

/**
 * Ring-system atoms → template coordinates, or null when no template matches.
 * Among the template's symmetry-equivalent mappings, the one that puts substituted atoms on the
 * most exposed template positions wins.
 */
export function matchCageTemplate(mol: Mol, atoms: number[], bonds: number[]): Map<number, Pt> | null {
  const n = atoms.length;
  const cand = TEMPLATES.filter((t) => t.xy.length === n && t.edges.length === bonds.length);
  if (!cand.length) return null;
  const local = new Map<number, number>();
  atoms.forEach((a, i) => local.set(a, i));
  const gAdj: Set<number>[] = atoms.map(() => new Set());
  for (const bi of bonds) {
    const b = mol.bonds[bi];
    const i = local.get(b.a)!, j = local.get(b.b)!;
    gAdj[i].add(j); gAdj[j].add(i);
  }
  // substituent load of each system atom
  const load = atoms.map((a) => {
    let w = 0;
    for (const j of mol.neighbors(a)) if (!local.has(j)) w += mol.atoms[j].el === 'H' ? 0.2 : 1;
    return w;
  });
  for (const t of cand) {
    const tAdj: Set<number>[] = t.xy.map(() => new Set());
    for (const [i, j] of t.edges) { tAdj[i].add(j); tAdj[j].add(i); }
    const expo = exposure(t, tAdj);
    let best: number[] | null = null, bestS = -Infinity;
    isomorphisms(tAdj, gAdj, (map) => {
      let s = 0;
      map.forEach((g, ti) => (s += load[g] * expo[ti]));
      if (s > bestS + 1e-9) { bestS = s; best = map.slice(); }
    });
    if (!best) continue;
    const out = new Map<number, Pt>();
    (best as number[]).forEach((g, ti) => out.set(atoms[g], { x: t.xy[ti][0], y: t.xy[ti][1] }));
    return out;
  }
  return null;
}

/**
 * How much room each template atom has for substituents: clearance of a probe point one bond
 * length out along the bisector of the atom's widest free gap (distance to other atoms/bonds).
 */
function exposure(t: CageTemplate, adj: Set<number>[]): number[] {
  const n = t.xy.length;
  const segDist = (px: number, py: number, a: number[], b: number[]) => {
    const dx = b[0] - a[0], dy = b[1] - a[1];
    const l2 = dx * dx + dy * dy || 1;
    const u = Math.max(0, Math.min(1, ((px - a[0]) * dx + (py - a[1]) * dy) / l2));
    return Math.hypot(px - a[0] - u * dx, py - a[1] - u * dy);
  };
  return t.xy.map(([x, y], i) => {
    const angs = [...adj[i]].map((j) => Math.atan2(t.xy[j][1] - y, t.xy[j][0] - x)).sort((a, b) => a - b);
    let gap = 0, mid = 0;
    for (let k = 0; k < angs.length; k++) {
      const a = angs[k], b = k + 1 < angs.length ? angs[k + 1] : angs[0] + 2 * Math.PI;
      if (b - a > gap) { gap = b - a; mid = (a + b) / 2; }
    }
    const px = x + Math.cos(mid), py = y + Math.sin(mid);
    let clear = 1.2;
    for (let j = 0; j < n; j++) if (j !== i) clear = Math.min(clear, Math.hypot(px - t.xy[j][0], py - t.xy[j][1]));
    for (const [a, b] of t.edges) if (a !== i && b !== i) clear = Math.min(clear, segDist(px, py, t.xy[a], t.xy[b]));
    return clear + (0.2 * gap) / Math.PI;
  });
}

/** Enumerates graph isomorphisms (template→graph index maps) by backtracking; capped. */
function isomorphisms(tAdj: Set<number>[], gAdj: Set<number>[], visit: (map: number[]) => void, cap = 500): void {
  const n = tAdj.length;
  const tDeg = tAdj.map((s) => s.size).sort().join(',');
  const gDeg = gAdj.map((s) => s.size).sort().join(',');
  if (tDeg !== gDeg) return;
  // BFS order over template atoms starting from max degree
  const order: number[] = [];
  const seen = new Array(n).fill(false);
  let start = 0;
  for (let i = 1; i < n; i++) if (tAdj[i].size > tAdj[start].size) start = i;
  const q = [start];
  seen[start] = true;
  while (q.length) {
    const v = q.shift()!;
    order.push(v);
    for (const w of tAdj[v]) if (!seen[w]) { seen[w] = true; q.push(w); }
  }
  if (order.length !== n) return;
  const map = new Array(n).fill(-1);
  const used = new Array(n).fill(false);
  let found = 0;
  const rec = (k: number): void => {
    if (found >= cap) return;
    if (k === n) { found++; visit(map); return; }
    const t = order[k];
    for (let g = 0; g < n; g++) {
      if (used[g] || gAdj[g].size !== tAdj[t].size) continue;
      let ok = true;
      for (let j = 0; j < k && ok; j++) {
        const tj = order[j];
        if (tAdj[t].has(tj) !== gAdj[g].has(map[tj])) ok = false;
      }
      if (!ok) continue;
      map[t] = g; used[g] = true;
      rec(k + 1);
      used[g] = false;
    }
    map[t] = -1;
  };
  rec(0);
}
