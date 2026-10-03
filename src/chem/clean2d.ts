// "Clean Structure": regularises an existing (possibly hand-drawn) depiction while keeping it
// recognisably the same drawing.
//
//  1. stereo as currently drawn is perceived first (E/Z from coordinates, tetrahedral from wedges);
//  2. the drawing is rescaled to unit bond length; every ring system is snapped to its ideal
//     layout (best-fit rotation/translation onto the current ring atoms);
//  3. a distance-constraint relaxation (unit bonds, ideal angles in the current cyclic neighbour
//     order, rigid ring systems, cis/trans 1-4 distances, repulsion) starts from those coordinates;
//  4. clashes are resolved, the result is aligned back onto the original (rotation + translation,
//     never a reflection), E/Z is repaired if needed and wedges are re-derived when the drawn
//     configuration changed. Very messy input falls back to a fresh layout aligned to the original.
import { Mol, DbSpec, TetraSpec, tetraCcwForOrder } from './mol';
import { perceiveRings, RingInfo, smallestRingSizeOfBond } from './rings';
import { perceiveStereo2D, assignWedgesFromSpecs } from './stereo2d';
import { layoutRingSystem } from './layout/ringsys';
import { layoutMol } from './layout2d';
import { Pt, fitXform, applyXform, normAngle, cross3, segmentsCross, TAU } from './layout/geom';
import { relax, DistCon, pairKey } from './layout/refine';
import { resolveOverlaps } from './layout/overlap';

export interface CleanOptions {
  /** Only these atoms move; all others act as anchors. */
  atoms?: Set<number>;
  /** Relaxation passes (default 200). */
  iterations?: number;
}

export function clean2D(mol: Mol, opts: CleanOptions = {}): void {
  const n = mol.atoms.length;
  if (n < 2 || !mol.bonds.length) return;
  const movable = new Uint8Array(n);
  for (let i = 0; i < n; i++) movable[i] = !opts.atoms || opts.atoms.has(i) ? 1 : 0;
  const allMove = !opts.atoms || mol.atoms.every((_, i) => movable[i]);
  if (!movable.some((v) => v)) return;

  const info = perceiveRings(mol);
  const orig = mol.atoms.map((a) => ({ x: a.x, y: a.y }));

  // 1. stereo as drawn
  const drawn = mol.clone();
  const { tetra: drawnTetra, db: drawnDb } = perceiveStereo2D(drawn, { rings: info });
  const savedTetra = mol.tetra, savedDb = mol.dbStereo;

  const xy = new Float64Array(2 * n);
  mol.atoms.forEach((a, i) => { xy[2 * i] = a.x; xy[2 * i + 1] = a.y; });

  // degenerate input (no usable coordinates): plain layout
  const lens = mol.bonds.map((b) => Math.hypot(xy[2 * b.a] - xy[2 * b.b], xy[2 * b.a + 1] - xy[2 * b.b + 1])).sort((p, q) => p - q);
  const med = lens[Math.floor(lens.length / 2)];
  if (!(med > 1e-3)) {
    if (allMove) {
      mol.dbStereo = drawnDb.length ? drawnDb : savedDb;
      try { layoutMol(mol); } finally { mol.dbStereo = savedDb; }
    } else {
      const fixed = new Set<number>();
      for (let i = 0; i < n; i++) if (!movable[i]) fixed.add(i);
      layoutMol(mol, { fixed });
    }
    return;
  }

  // 2. scale to unit bonds (only when everything moves)
  if (allMove && Math.abs(med - 1) > 0.02) {
    let cx = 0, cy = 0;
    for (let i = 0; i < n; i++) { cx += xy[2 * i]; cy += xy[2 * i + 1]; }
    cx /= n; cy /= n;
    for (let i = 0; i < n; i++) { xy[2 * i] = cx + (xy[2 * i] - cx) / med; xy[2 * i + 1] = cy + (xy[2 * i + 1] - cy) / med; }
  }

  // ideal ring-system geometry, fitted onto the current drawing (macrocycles honour drawn E/Z)
  mol.dbStereo = drawnDb;
  const ringIdeal = new Map<number, Pt>(); // atom -> snapped ideal position
  const sysOf = new Int32Array(n).fill(-1);
  const sysAtomsList: number[][] = [];
  try {
    info.systems.forEach((_, s) => {
      const L = layoutRingSystem(mol, info, s);
      const atoms = [...L.pos.keys()];
      sysAtomsList.push(atoms);
      for (const a of atoms) sysOf[a] = s;
      const src = atoms.map((a) => L.pos.get(a)!);
      const dst = atoms.map((a) => ({ x: xy[2 * a], y: xy[2 * a + 1] }));
      // anchors (immovable atoms) dominate the fit
      const w = atoms.map((a) => (movable[a] ? 1 : 25));
      const { xf } = fitXform(src, dst, true, w);
      atoms.forEach((a, k) => ringIdeal.set(a, applyXform(xf, src[k])));
    });
  } finally {
    mol.dbStereo = savedDb;
  }
  for (const [a, p] of ringIdeal) if (movable[a]) { xy[2 * a] = p.x; xy[2 * a + 1] = p.y; }

  // 3. constraints
  const cons = buildConstraints(mol, info, xy, ringIdeal, sysOf, sysAtomsList);
  const excl = new Set<number>();
  for (const c of cons) excl.add(pairKey(c.i, c.j, n));
  const all = mol.atoms.map((_, i) => i);
  relax(xy, n, cons, { iterations: opts.iterations ?? 200, movable, atoms: all, repelDist: 1.0, repelWeight: 0.4, exclude: excl });

  // 4. clashes, then a short polish
  const dbActive = drawnDb.filter((d) => isActiveDb(info, d));
  const stereoOK = (moved: number[]) => {
    const ms = new Set(moved);
    for (const d of dbActive) {
      const b = mol.bonds[d.bond];
      if (![b.a, b.b, d.a, d.b].some((a) => ms.has(a))) continue;
      if (!dbOK(mol, xy, d)) return false;
    }
    return true;
  };
  for (const comp of mol.components()) {
    if (comp.length < 4) continue;
    resolveOverlaps(mol, xy, info, { atoms: comp, movable: movable, stereoOK, ringSystems: sysAtomsList });
  }
  relax(xy, n, cons, { iterations: 40, movable, atoms: all, repelDist: 1.0, repelWeight: 0.4, exclude: excl });

  // E/Z repair (mirror the smaller side of an acyclic double bond)
  for (const d of dbActive) {
    if (dbOK(mol, xy, d)) continue;
    if (info.bondInRing[d.bond]) continue;
    mirrorSide(mol, xy, d.bond, movable);
  }

  // messy input: compare with a fresh layout fitted onto the original drawing
  let result = xy;
  if (allMove) {
    alignOnto(xy, orig, false);
    const q = quality(mol, xy);
    if (q > 1.5) {
      const alt = mol.clone();
      alt.dbStereo = drawnDb;
      alt.tetra = drawnTetra;
      layoutMol(alt);
      const axy = new Float64Array(2 * n);
      alt.atoms.forEach((a, i) => { axy[2 * i] = a.x; axy[2 * i + 1] = a.y; });
      alignOnto(axy, orig, false);
      if (quality(mol, axy) + 0.5 < q) result = axy;
    }
  }
  for (let i = 0; i < n; i++) {
    if (!movable[i]) continue;
    mol.atoms[i].x = result[2 * i];
    mol.atoms[i].y = result[2 * i + 1];
  }

  // tetrahedral stereo: re-derive wedges for centres whose drawn configuration changed
  if (drawnTetra.length) {
    const now = mol.clone();
    const { tetra: nowTetra } = perceiveStereo2D(now, { rings: info });
    const changed = drawnTetra.filter((t) => !sameConfig(t, nowTetra.find((u) => u.center === t.center)));
    if (changed.length) {
      mol.tetra = changed;
      try { assignWedgesFromSpecs(mol, info); } finally { mol.tetra = savedTetra; }
    }
  }
  mol.tetra = savedTetra;
  mol.dbStereo = savedDb;
}

// ───────────────────────────── constraints ─────────────────────────────

function buildConstraints(mol: Mol, info: RingInfo, xy: Float64Array, ringIdeal: Map<number, Pt>, sysOf: Int32Array, systems: number[][]): DistCon[] {
  const n = mol.atoms.length;
  const cons: DistCon[] = [];
  const have = new Set<number>();
  const add = (i: number, j: number, d: number, w: number) => {
    const k = pairKey(i, j, n);
    if (have.has(k)) return;
    have.add(k);
    cons.push({ i, j, d, w });
  };
  for (const b of mol.bonds) add(b.a, b.b, 1, 0.9);
  // rigid ring systems (all pairs for normal sizes, near pairs for very large systems)
  for (const atoms of systems) {
    for (let p = 0; p < atoms.length; p++) {
      for (let q = p + 1; q < atoms.length; q++) {
        const a = ringIdeal.get(atoms[p])!, b = ringIdeal.get(atoms[q])!;
        const d = Math.hypot(a.x - b.x, a.y - b.y);
        if (atoms.length > 40 && d > 3.2) continue;
        add(atoms[p], atoms[q], d, d < 1.8 ? 0.6 : 0.35);
      }
    }
  }
  const ang = (c: number, j: number) => Math.atan2(xy[2 * j + 1] - xy[2 * c + 1], xy[2 * j] - xy[2 * c]);
  // ideal angles around every atom, keeping the current cyclic order of neighbours
  for (let c = 0; c < n; c++) {
    const nb = mol.neighbors(c);
    if (nb.length < 2) continue;
    const target = new Map<number, number>(); // neighbour -> target angle
    const s = sysOf[c];
    const ringNb = s >= 0 ? nb.filter((j) => sysOf[j] === s) : [];
    const exo = nb.filter((j) => !ringNb.includes(j));
    if (ringNb.length >= 2 && exo.length) {
      // substituents of ring atoms: evenly inside the widest exterior gap of the ideal ring geometry
      const cp = ringIdeal.get(c)!;
      const ra = ringNb.map((j) => ({ j, a: normAngle(Math.atan2(ringIdeal.get(j)!.y - cp.y, ringIdeal.get(j)!.x - cp.x)) })).sort((p, q) => p.a - q.a);
      let best = { start: 0, size: 0 }, bestS = -Infinity;
      for (let k = 0; k < ra.length; k++) {
        const p = ra[k], q = ra[(k + 1) % ra.length];
        let size = q.a - p.a;
        if (size <= 0) size += TAU;
        const interior = ringsShare(info, c, p.j, q.j);
        const sc = size - (interior ? Math.PI : 0);
        if (sc > bestS) { bestS = sc; best = { start: p.a, size }; }
      }
      // order substituents by their current angle measured from the gap start
      const cur = exo.map((j) => ({ j, a: normAngle(ang(c, j) - best.start) })).sort((p, q) => p.a - q.a);
      const m = cur.length;
      const step = Math.min(best.size / (m + 1), (2 * Math.PI) / 3);
      const mid = best.start + best.size / 2;
      cur.forEach((e, k) => target.set(e.j, mid + (k - (m - 1) / 2) * step));
      for (const r of ra) target.set(r.j, r.a);
    } else if (ringNb.length < 2) {
      // chain atom: 120° (trigonal / zig-zag), 180° (sp), 90° (cross) or even spacing
      const sorted = nb.map((j) => ({ j, a: ang(c, j) })).sort((p, q) => p.a - q.a);
      const d = nb.length;
      const linear = d === 2 && isLinear(mol, c);
      const steps = d === 2 ? (linear ? Math.PI : (2 * Math.PI) / 3) : d === 3 ? (2 * Math.PI) / 3 : d === 4 ? Math.PI / 2 : TAU / d;
      sorted.forEach((e, k) => target.set(e.j, sorted[0].a + k * steps));
    } else continue;
    const list = [...target.keys()];
    for (let p = 0; p < list.length; p++) {
      for (let q = p + 1; q < list.length; q++) {
        const i = list[p], j = list[q];
        if (sysOf[i] >= 0 && sysOf[i] === sysOf[j] && sysOf[c] === sysOf[i]) continue; // ring-internal
        let da = Math.abs(normAngle(target.get(i)! - target.get(j)!));
        if (da > Math.PI) da = TAU - da;
        add(i, j, 2 * Math.sin(da / 2), 0.5);
      }
    }
  }
  // cis/trans 1-4 distances across double bonds, as currently drawn
  mol.bonds.forEach((b, bi) => {
    if (b.order !== 2) return;
    const rs = smallestRingSizeOfBond(info, bi);
    if (rs > 0 && rs < 8) return;
    const A = { x: xy[2 * b.a], y: xy[2 * b.a + 1] }, B = { x: xy[2 * b.b], y: xy[2 * b.b + 1] };
    for (const x of mol.neighbors(b.a)) {
      if (x === b.b) continue;
      for (const y of mol.neighbors(b.b)) {
        if (y === b.a) continue;
        const sx = cross3(A, B, { x: xy[2 * x], y: xy[2 * x + 1] });
        const sy = cross3(A, B, { x: xy[2 * y], y: xy[2 * y + 1] });
        if (Math.abs(sx) < 1e-3 || Math.abs(sy) < 1e-3) continue;
        add(x, y, sx * sy > 0 ? 2.0 : Math.sqrt(7), 0.3);
      }
    }
  });
  return cons;
}

function ringsShare(info: RingInfo, c: number, p: number, q: number): boolean {
  for (const r of info.atomRings[c]) {
    const ring = info.rings[r];
    const k = ring.indexOf(c);
    const a = ring[(k + 1) % ring.length], b = ring[(k - 1 + ring.length) % ring.length];
    if ((a === p && b === q) || (a === q && b === p)) return true;
  }
  return false;
}

function isLinear(mol: Mol, u: number): boolean {
  const adj = mol.adj[u];
  if (adj.length !== 2) return false;
  const o1 = mol.bonds[adj[0]].order, o2 = mol.bonds[adj[1]].order;
  return o1 === 3 || o2 === 3 || (o1 === 2 && o2 === 2);
}

// ───────────────────────────── stereo helpers ─────────────────────────────

function isActiveDb(info: RingInfo, d: DbSpec): boolean {
  const rs = smallestRingSizeOfBond(info, d.bond);
  return rs === 0 || rs >= 8;
}

function dbOK(mol: Mol, xy: Float64Array, d: DbSpec): boolean {
  const b = mol.bonds[d.bond];
  let ra = d.a, rb = d.b;
  if (mol.bondBetween(b.a, ra) < 0) { ra = d.b; rb = d.a; }
  const P = (i: number): Pt => ({ x: xy[2 * i], y: xy[2 * i + 1] });
  const sa = cross3(P(b.a), P(b.b), P(ra));
  const sb = cross3(P(b.a), P(b.b), P(rb));
  if (Math.abs(sa) < 1e-4 || Math.abs(sb) < 1e-4) return false;
  return sa * sb > 0 === d.cis;
}

/** Mirrors the smaller side of bond bi across the bond axis (if that side is movable). */
function mirrorSide(mol: Mol, xy: Float64Array, bi: number, movable: Uint8Array): void {
  const b = mol.bonds[bi];
  const side = (start: number, other: number) => {
    const seen = new Set([start, other]);
    const st = [start];
    const out: number[] = [];
    while (st.length) {
      const v = st.pop()!;
      for (const w of mol.neighbors(v)) if (!seen.has(w)) { seen.add(w); out.push(w); st.push(w); }
    }
    return out;
  };
  const sa = side(b.a, b.b), sb = side(b.b, b.a);
  const cands = [sa, sb].filter((s) => s.every((i) => movable[i]) && !s.includes(b.a) && !s.includes(b.b)).sort((p, q) => p.length - q.length);
  if (!cands.length) return;
  const px = xy[2 * b.a], py = xy[2 * b.a + 1];
  const dx = xy[2 * b.b] - px, dy = xy[2 * b.b + 1] - py;
  const l2 = dx * dx + dy * dy || 1;
  for (const i of cands[0]) {
    const t = ((xy[2 * i] - px) * dx + (xy[2 * i + 1] - py) * dy) / l2;
    const fx = px + t * dx, fy = py + t * dy;
    xy[2 * i] = 2 * fx - xy[2 * i];
    xy[2 * i + 1] = 2 * fy - xy[2 * i + 1];
  }
}

function sameConfig(a: TetraSpec, b: TetraSpec | undefined): boolean {
  if (!b) return false;
  const sa = [...a.nbrs].sort((p, q) => p - q).join(','), sb = [...b.nbrs].sort((p, q) => p - q).join(',');
  if (sa !== sb) return false;
  return tetraCcwForOrder(a, b.nbrs) === b.ccw;
}

// ───────────────────────────── alignment & quality ─────────────────────────────

/** Rigidly moves xy onto the reference points (rotation + translation, optional reflection). */
function alignOnto(xy: Float64Array, ref: Pt[], allowReflection: boolean): void {
  const n = ref.length;
  const src: Pt[] = [], dst: Pt[] = [];
  for (let i = 0; i < n; i++) { src.push({ x: xy[2 * i], y: xy[2 * i + 1] }); dst.push(ref[i]); }
  // normalise reference scale for the fit only
  const { xf } = fitXform(src, dst, allowReflection);
  for (let i = 0; i < n; i++) {
    const p = applyXform(xf, src[i]);
    xy[2 * i] = p.x; xy[2 * i + 1] = p.y;
  }
}

/** Depiction problems: bond length deviations, close contacts, crossings. */
function quality(mol: Mol, xy: Float64Array): number {
  let q = 0;
  const P = (i: number): Pt => ({ x: xy[2 * i], y: xy[2 * i + 1] });
  for (const b of mol.bonds) {
    const d = Math.hypot(xy[2 * b.a] - xy[2 * b.b], xy[2 * b.a + 1] - xy[2 * b.b + 1]);
    if (d < 0.8 || d > 1.25) q += 0.5;
  }
  const n = mol.atoms.length;
  if (n <= 400) {
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      if (mol.bondBetween(i, j) >= 0) continue;
      if (Math.hypot(xy[2 * i] - xy[2 * j], xy[2 * i + 1] - xy[2 * j + 1]) < 0.5) q += 1;
    }
    const bl = mol.bonds;
    for (let i = 0; i < bl.length; i++) for (let j = i + 1; j < bl.length; j++) {
      const p = bl[i], r = bl[j];
      if (p.a === r.a || p.a === r.b || p.b === r.a || p.b === r.b) continue;
      if (segmentsCross(P(p.a), P(p.b), P(r.a), P(r.b))) q += 1;
    }
  }
  return q;
}
