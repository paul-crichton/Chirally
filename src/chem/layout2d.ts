// 2D coordinate generation ("structure diagram generation") for depiction.
//
// Pipeline per connected component:
//   1. perceive rings; lay out every ring system as a rigid block (see layout/ringsys.ts);
//   2. start from the largest ring system (or one end of the longest chain) and grow outwards
//      breadth-first: substituent directions are spread evenly around each atom (ring substituents
//      radially, chains as trans zig-zags, sp centres linear, 4-coordinate chain atoms as a cross),
//      ring blocks are attached through their exterior bisector; E/Z specs choose sides/flips;
//   3. resolve clashes by flipping/rotating around acyclic single bonds (layout/overlap.ts);
//   4. orient (long axis horizontal, bonds snapped to the 30°/90°/150° lattice);
// then components are arranged left-to-right (counter-ions next to their partner) and wedges are
// drawn for mol.tetra. Bond length is 1, y points down. Deterministic.
import { Mol, DbSpec } from './mol';
import { perceiveRings, RingInfo, smallestRingSizeOfBond } from './rings';
import { assignWedgesFromSpecs } from './stereo2d';
import { layoutRingSystem, SysLayout } from './layout/ringsys';
import { Pt, TAU, normAngle, angleTo, cross3, fitXform, applyXform } from './layout/geom';
import { resolveOverlaps, relaxComponent } from './layout/overlap';

export interface LayoutOptions {
  /** Atoms that keep their current coordinates; everything else is placed around them. */
  fixed?: Set<number>;
}

const GAP = 1.5; // spacing between components

/** Assigns 2D depiction coordinates to all atoms (or only those not in opts.fixed). */
export function layoutMol(mol: Mol, opts: LayoutOptions = {}): void {
  const n = mol.atoms.length;
  if (!n) return;
  const fixed = opts.fixed && opts.fixed.size ? new Set([...opts.fixed].filter((i) => i >= 0 && i < n)) : null;
  if (fixed && fixed.size === n) return;
  const info = perceiveRings(mol);
  const L = new Layout(mol, info, fixed);
  L.run();
  for (let i = 0; i < n; i++) {
    if (fixed && fixed.has(i)) continue;
    mol.atoms[i].x = L.xy[2 * i];
    mol.atoms[i].y = L.xy[2 * i + 1];
  }
  if (mol.tetra.length) drawWedges(mol, info, fixed);
}

/** Draws wedges for mol.tetra, leaving stereocentres among fixed atoms untouched. */
function drawWedges(mol: Mol, info: RingInfo, fixed: Set<number> | null): void {
  if (!fixed) {
    assignWedgesFromSpecs(mol, info);
    return;
  }
  const all = mol.tetra;
  const keep = new Map<number, { a: number; b: number; style: string }>();
  mol.bonds.forEach((b, bi) => {
    if ((b.style === 'wedge' || b.style === 'hash') && fixed.has(b.a)) keep.set(bi, { a: b.a, b: b.b, style: b.style });
  });
  mol.tetra = all.filter((t) => !fixed.has(t.center));
  try {
    if (mol.tetra.length) assignWedgesFromSpecs(mol, info);
  } finally {
    mol.tetra = all;
  }
  // restore wedges of fixed centres that a new centre may have taken over
  for (const [bi, k] of keep) {
    const b = mol.bonds[bi];
    b.a = k.a; b.b = k.b; b.style = k.style as typeof b.style;
  }
}

/** Stereo double-bond specs that the depiction must honour (acyclic or in rings ≥ 8). */
function activeDbSpecs(mol: Mol, info: RingInfo): DbSpec[] {
  return mol.dbStereo.filter((d) => {
    const b = mol.bonds[d.bond];
    if (!b || b.order !== 2) return false;
    const rs = smallestRingSizeOfBond(info, d.bond);
    return rs === 0 || rs >= 8;
  });
}

/** True when the drawn geometry agrees with the spec (false also for degenerate/linear geometry). */
function dbSatisfied(mol: Mol, xy: Float64Array, d: DbSpec): boolean {
  const b = mol.bonds[d.bond];
  const P = (i: number): Pt => ({ x: xy[2 * i], y: xy[2 * i + 1] });
  // orient refs so that ra hangs on b.a
  let ra = d.a, rb = d.b;
  if (mol.bondBetween(b.a, ra) < 0) { ra = d.b; rb = d.a; }
  const sa = cross3(P(b.a), P(b.b), P(ra));
  const sb = cross3(P(b.a), P(b.b), P(rb));
  if (Math.abs(sa) < 1e-4 || Math.abs(sb) < 1e-4) return false;
  return sa * sb > 0 === d.cis;
}

class Layout {
  readonly n: number;
  readonly xy: Float64Array;
  private placed: Uint8Array;
  private parent: Int32Array;
  private atomSys: Int32Array;
  private sysLayouts = new Map<number, SysLayout>();
  private sysPlaced: boolean[];
  private specs: DbSpec[];
  private specsByAtom = new Map<number, DbSpec[]>();
  private grid = new Map<number, number[]>();
  private branchMemo = new Map<number, number>();

  constructor(private mol: Mol, private info: RingInfo, private fixed: Set<number> | null) {
    const n = (this.n = mol.atoms.length);
    this.xy = new Float64Array(2 * n);
    this.placed = new Uint8Array(n);
    this.parent = new Int32Array(n).fill(-1);
    this.atomSys = new Int32Array(n).fill(-1);
    info.systems.forEach((sys, s) => {
      for (const r of sys) for (const a of info.rings[r]) this.atomSys[a] = s;
    });
    this.sysPlaced = info.systems.map(() => false);
    this.specs = activeDbSpecs(mol, info);
    for (const d of this.specs) {
      const b = mol.bonds[d.bond];
      for (const a of new Set([b.a, b.b, d.a, d.b])) {
        let l = this.specsByAtom.get(a);
        if (!l) this.specsByAtom.set(a, (l = []));
        l.push(d);
      }
    }
    if (fixed) {
      for (const i of fixed) {
        this.xy[2 * i] = mol.atoms[i].x;
        this.xy[2 * i + 1] = mol.atoms[i].y;
      }
    }
  }

  run(): void {
    const comps = this.mol.components();
    const free: number[][] = [];
    const anchored: number[][] = [];
    for (const c of comps) {
      if (this.fixed && c.some((a) => this.fixed!.has(a))) anchored.push(c);
      else free.push(c);
    }
    for (const c of anchored) this.layoutAnchored(c);
    for (const c of free) {
      this.layoutComponent(c);
      this.orient(c);
    }
    this.arrange(free, anchored);
  }

  // ───────────────────────── placement helpers ─────────────────────────

  private P(i: number): Pt {
    return { x: this.xy[2 * i], y: this.xy[2 * i + 1] };
  }

  private cellKey(x: number, y: number): number {
    return (Math.floor(x / 1.2) + 32768) * 65536 + (Math.floor(y / 1.2) + 32768);
  }

  private setPos(i: number, p: Pt): void {
    this.xy[2 * i] = p.x;
    this.xy[2 * i + 1] = p.y;
    if (!this.placed[i]) {
      this.placed[i] = 1;
      const k = this.cellKey(p.x, p.y);
      let c = this.grid.get(k);
      if (!c) this.grid.set(k, (c = []));
      c.push(i);
    }
  }

  private rebuildGrid(): void {
    this.grid.clear();
    for (let i = 0; i < this.n; i++) {
      if (!this.placed[i]) continue;
      const k = this.cellKey(this.xy[2 * i], this.xy[2 * i + 1]);
      let c = this.grid.get(k);
      if (!c) this.grid.set(k, (c = []));
      c.push(i);
    }
  }

  /** Crowding penalty of putting atom `self` at p (ignores `self`, its bonded neighbours and `skip`). */
  private crowd(p: Pt, self: number, skip?: Set<number>): number {
    let s = 0;
    const cx = Math.floor(p.x / 1.2), cy = Math.floor(p.y / 1.2);
    for (let gx = cx - 1; gx <= cx + 1; gx++) {
      for (let gy = cy - 1; gy <= cy + 1; gy++) {
        const c = this.grid.get((gx + 32768) * 65536 + (gy + 32768));
        if (!c) continue;
        for (const j of c) {
          if (j === self || (skip && skip.has(j))) continue;
          const d = Math.hypot(this.xy[2 * j] - p.x, this.xy[2 * j + 1] - p.y);
          if (d < 1.05 && this.mol.bondBetween(self, j) < 0) s += (1.05 - d) * (1.05 - d) * (d < 0.5 ? 10 : 1);
        }
      }
    }
    return s;
  }

  /** Number of violated E/Z specs among those touching `atoms` whose atoms are all placed. */
  private stereoViolations(atoms: number[]): number {
    if (!this.specs.length) return 0;
    let v = 0;
    const seen = new Set<DbSpec>();
    for (const a of atoms) {
      for (const d of this.specsByAtom.get(a) ?? []) {
        if (seen.has(d)) continue;
        seen.add(d);
        const b = this.mol.bonds[d.bond];
        if (!this.placed[b.a] || !this.placed[b.b] || !this.placed[d.a] || !this.placed[d.b]) continue;
        if (!dbSatisfied(this.mol, this.xy, d)) v++;
      }
    }
    return v;
  }

  /** Heavy-atom size of the branch hanging from `from` through `to`. */
  private branchWeight(from: number, to: number): number {
    const key = from * this.n + to;
    const m = this.branchMemo.get(key);
    if (m !== undefined) return m;
    const seen = new Set<number>([from, to]);
    const st = [to];
    let w = 0;
    while (st.length) {
      const v = st.pop()!;
      w += this.mol.atoms[v].el === 'H' ? 0.01 : 1;
      for (const x of this.mol.neighbors(v)) if (!seen.has(x)) { seen.add(x); st.push(x); }
      if (w > 5000) break;
    }
    this.branchMemo.set(key, w);
    return w;
  }

  private isLinear(u: number): boolean {
    const adj = this.mol.adj[u];
    if (adj.length !== 2) return false;
    const o1 = this.mol.bonds[adj[0]].order, o2 = this.mol.bonds[adj[1]].order;
    return o1 === 3 || o2 === 3 || (o1 === 2 && o2 === 2);
  }

  private sysLayout(s: number): SysLayout {
    let L = this.sysLayouts.get(s);
    if (!L) {
      L = layoutRingSystem(this.mol, this.info, s);
      this.sysLayouts.set(s, L);
    }
    return L;
  }

  /** True if p–c–q are consecutive atoms of one SSSR ring (gap p..q at c is a ring interior). */
  private interior(c: number, p: number, q: number): boolean {
    for (const r of this.info.atomRings[c]) {
      const ring = this.info.rings[r];
      const k = ring.indexOf(c);
      const a = ring[(k + 1) % ring.length], b = ring[(k - 1 + ring.length) % ring.length];
      if ((a === p && b === q) || (a === q && b === p)) return true;
    }
    return false;
  }

  /**
   * Best free angular gap at atom c given neighbour positions; prefers ring exteriors,
   * then wide gaps, then gaps facing away from `away` (e.g. the ring-system centroid).
   */
  private bestGap(c: number, cp: Pt, nbrs: number[], pos: (i: number) => Pt, away?: Pt): { start: number; size: number } {
    // ring interiors are only meaningful for planar (polygon/fused) ring-system layouts
    const sys = this.atomSys[c];
    const planar = sys < 0 || ['polygon', 'macrocycle', 'fused'].includes(this.sysLayout(sys).method);
    const items = nbrs.map((j) => ({ j, a: normAngle(angleTo(cp, pos(j))) })).sort((p, q) => p.a - q.a);
    let best = { start: 0, size: TAU }, bestS = -Infinity;
    const awayAng = away && Math.hypot(cp.x - away.x, cp.y - away.y) > 1e-3 ? angleTo(away, cp) : null;
    for (let k = 0; k < items.length; k++) {
      const p = items[k], q = items[(k + 1) % items.length];
      let size = q.a - p.a;
      if (size <= 0) size += TAU;
      if (items.length === 1) size = TAU;
      let s = size;
      if (planar && items.length > 1 && this.interior(c, p.j, q.j)) s -= Math.PI;
      if (awayAng !== null) s += 0.35 * Math.cos(p.a + size / 2 - awayAng);
      if (s > bestS + 1e-9) { bestS = s; best = { start: p.a, size }; }
    }
    return best;
  }

  // ───────────────────────── component layout ─────────────────────────

  private layoutComponent(comp: number[]): void {
    const { mol, info } = this;
    if (comp.length === 1) {
      this.setPos(comp[0], { x: 0, y: 0 });
      return;
    }
    // root: the largest ring system, else one end of the longest chain
    const systems = new Set<number>();
    for (const a of comp) if (this.atomSys[a] >= 0) systems.add(this.atomSys[a]);
    const queue: number[] = [];
    if (systems.size) {
      let best = -1, bestSize = -1;
      for (const s of [...systems].sort((a, b) => a - b)) {
        const atomsN = new Set(info.systems[s].flatMap((r) => info.rings[r])).size;
        if (atomsN > bestSize) { bestSize = atomsN; best = s; }
      }
      const L = this.sysLayout(best);
      // centre the block at the origin
      let cx = 0, cy = 0;
      for (const p of L.pos.values()) { cx += p.x; cy += p.y; }
      cx /= L.pos.size; cy /= L.pos.size;
      const blockAtoms = [...L.pos.keys()].sort((a, b) => a - b);
      for (const a of blockAtoms) this.setPos(a, { x: L.pos.get(a)!.x - cx, y: L.pos.get(a)!.y - cy });
      this.sysPlaced[best] = true;
      queue.push(...blockAtoms);
    } else {
      const root = longestPathEnd(mol, comp);
      this.setPos(root, { x: 0, y: 0 });
      queue.push(root);
    }
    this.grow(queue);
    this.resolve(comp);
  }

  /** Layout of a component containing fixed atoms. */
  private layoutAnchored(comp: number[]): void {
    const fixed = this.fixed!;
    for (const a of comp) if (fixed.has(a)) this.setPos(a, { x: this.mol.atoms[a].x, y: this.mol.atoms[a].y });
    // ring systems with fixed atoms are fitted onto them
    const systems = new Set<number>();
    for (const a of comp) if (this.atomSys[a] >= 0) systems.add(this.atomSys[a]);
    for (const s of [...systems].sort((a, b) => a - b)) {
      const L = this.sysLayout(s);
      const atoms = [...L.pos.keys()];
      const fx = atoms.filter((a) => fixed.has(a));
      if (!fx.length) continue;
      this.sysPlaced[s] = true;
      if (fx.length === atoms.length) continue;
      if (fx.length >= 2) {
        const src = fx.map((a) => L.pos.get(a)!);
        const dst = fx.map((a) => this.P(a));
        const { xf } = fitXform(src, dst, true);
        for (const a of atoms) if (!fixed.has(a)) this.setPos(a, applyXform(xf, L.pos.get(a)!));
      } else {
        const f = fx[0];
        const outside = this.mol.neighbors(f).filter((j) => this.placed[j] && !L.pos.has(j));
        const fp = this.P(f);
        const gap = outside.length ? this.bestGap(f, fp, outside, (i) => this.P(i)) : { start: 0, size: TAU };
        const into = gap.start + gap.size / 2; // the ring goes into the widest free gap
        this.placeBlockAt(s, f, fp, into + Math.PI, new Set([f]));
      }
    }
    const queue = comp.filter((a) => this.placed[a]).sort((a, b) => a - b);
    this.grow(queue);
    this.resolve(comp);
    // bonds between fixed and re-placed parts may be stretched: polish the free atoms only
    let bad = false;
    for (const b of this.mol.bonds) {
      if (!comp.includes(b.a)) continue;
      const d = Math.hypot(this.xy[2 * b.a] - this.xy[2 * b.b], this.xy[2 * b.a + 1] - this.xy[2 * b.b + 1]);
      if ((d < 0.8 || d > 1.25) && !(fixed.has(b.a) && fixed.has(b.b))) bad = true;
    }
    if (bad) {
      const movable = new Uint8Array(this.n);
      for (const a of comp) movable[a] = fixed.has(a) ? 0 : 1;
      relaxComponent(this.mol, this.xy, comp, { atoms: comp, movable: movable, stereoOK: () => true, ringSystems: this.ringSystemsOf(comp) }, 80);
    }
  }

  private ringSystemsOf(comp: number[]): number[][] {
    const set = new Set(comp);
    return this.info.systems
      .map((sys) => [...new Set(sys.flatMap((r) => this.info.rings[r]))])
      .filter((atoms) => atoms.length && set.has(atoms[0]));
  }

  /** Breadth-first growth from already placed atoms. */
  private grow(queue: number[]): void {
    for (let qi = 0; qi < queue.length; qi++) {
      const u = queue[qi];
      const placedNow = this.placeChildren(u);
      for (const a of placedNow) queue.push(a);
    }
  }

  /** Places all unplaced neighbours of u; returns newly placed atoms. */
  private placeChildren(u: number): number[] {
    const { mol } = this;
    const nb = mol.neighbors(u);
    const placedNb = nb.filter((j) => this.placed[j]);
    let kids = nb.filter((j) => !this.placed[j]);
    if (!kids.length) return [];
    // heaviest branch first
    kids = kids
      .map((c) => ({ c, w: this.branchWeight(u, c) }))
      .sort((p, q) => q.w - p.w || p.c - q.c)
      .map((x) => x.c);
    const up = this.P(u);
    const assignments = this.candidateAngles(u, up, placedNb, kids);

    // evaluate assignments: stereo first, then crowding, then preference order
    let best = assignments[0], bestS = Infinity;
    const skip = new Set<number>([u]);
    assignments.forEach((angs, rank) => {
      const pts = angs.map((a) => ({ x: up.x + Math.cos(a), y: up.y + Math.sin(a) }));
      // tentative placement (atom positions only; blocks are flipped later)
      const tmp: number[] = [];
      kids.forEach((c, k) => {
        if (!this.placed[c]) { this.xy[2 * c] = pts[k].x; this.xy[2 * c + 1] = pts[k].y; this.placed[c] = 1; tmp.push(c); }
      });
      const viol = this.stereoViolations(kids);
      for (const c of tmp) this.placed[c] = 0;
      let cr = 0;
      kids.forEach((c, k) => (cr += this.crowd(pts[k], c, skip)));
      const s = viol * 1000 + cr * 10 + rank * 0.01;
      if (s < bestS) { bestS = s; best = angs; }
    });

    const out: number[] = [];
    // simple atoms first so that block flips can see them
    const blockKids: { c: number; ang: number }[] = [];
    kids.forEach((c, k) => {
      const s = this.atomSys[c];
      if (s >= 0 && !this.sysPlaced[s]) blockKids.push({ c, ang: best[k] });
      else {
        this.setPos(c, { x: up.x + Math.cos(best[k]), y: up.y + Math.sin(best[k]) });
        this.parent[c] = u;
        out.push(c);
      }
    });
    for (const { c, ang } of blockKids) {
      const s = this.atomSys[c];
      if (this.sysPlaced[s]) continue;
      const cp = { x: up.x + Math.cos(ang), y: up.y + Math.sin(ang) };
      const placedAtoms = this.placeBlockAt(s, c, cp, ang + Math.PI, new Set([u]));
      this.parent[c] = u;
      out.push(...placedAtoms);
    }
    return out;
  }

  /**
   * Places ring system s so that atom c sits at cp with its exterior direction pointing along
   * `extAngle` (towards the attachment). Chooses the mirror image by stereo and crowding.
   */
  private placeBlockAt(s: number, c: number, cp: Pt, extAngle: number, skip: Set<number>): number[] {
    const L = this.sysLayout(s);
    const local = L.pos;
    const pc = local.get(c)!;
    const ringNb = this.mol.neighbors(c).filter((j) => local.has(j));
    let cx = 0, cy = 0;
    for (const p of local.values()) { cx += p.x; cy += p.y; }
    cx /= local.size; cy /= local.size;
    const g = this.bestGap(c, pc, ringNb, (i) => local.get(i)!, { x: cx, y: cy });
    const ext = g.start + g.size / 2;
    const atoms = [...local.keys()].filter((a) => !this.placed[a] || a === c).sort((a, b) => a - b);
    let bestPts: Map<number, Pt> | null = null, bestS = Infinity;
    for (const flip of [false, true]) {
      const rot = extAngle - (flip ? -ext : ext);
      const cr = Math.cos(rot), sr = Math.sin(rot);
      const pts = new Map<number, Pt>();
      for (const a of atoms) {
        const p = local.get(a)!;
        const dx = p.x - pc.x;
        const dy = flip ? -(p.y - pc.y) : p.y - pc.y;
        pts.set(a, { x: cp.x + cr * dx - sr * dy, y: cp.y + sr * dx + cr * dy });
      }
      // tentative
      const tmp: number[] = [];
      for (const [a, p] of pts) {
        if (this.placed[a]) continue;
        this.xy[2 * a] = p.x; this.xy[2 * a + 1] = p.y; this.placed[a] = 1; tmp.push(a);
      }
      const viol = this.stereoViolations(atoms);
      for (const a of tmp) this.placed[a] = 0;
      let crw = 0;
      const sk = new Set([...skip, ...atoms]);
      for (const [a, p] of pts) crw += this.crowd(p, a, sk);
      const sc = viol * 1000 + crw * 10 + (flip ? 0.01 : 0);
      if (sc < bestS) { bestS = sc; bestPts = pts; }
    }
    for (const [a, p] of bestPts!) if (!this.placed[a]) this.setPos(a, p);
    this.sysPlaced[s] = true;
    return atoms;
  }

  /** Candidate direction assignments (angles aligned with `kids`), best preference first. */
  private candidateAngles(u: number, up: Pt, placedNb: number[], kids: number[]): number[][] {
    const m = kids.length;
    const out: number[][] = [];
    if (placedNb.length === 0) {
      if (m === 1) return [[-Math.PI / 6]];
      if (m === 2 && this.isLinear(u)) return [[0, Math.PI]];
      if (m === 2) return [[-Math.PI / 6, (7 * Math.PI) / 6], [(7 * Math.PI) / 6, -Math.PI / 6]];
      const start = m === 4 ? 0 : -Math.PI / 2;
      const base = kids.map((_, k) => start + (k * TAU) / m);
      return [base];
    }
    if (placedNb.length === 1) {
      const a = placedNb[0];
      const th = angleTo(up, this.P(a));
      if (this.isLinear(u) && m === 1) return [[th + Math.PI]];
      const transSign = this.transSign(u, a, up);
      if (m === 1) {
        const opts = [th + (2 * Math.PI) / 3, th - (2 * Math.PI) / 3];
        // trans option first
        if (transSign !== 0 && this.sideOf(a, u, opts[1]) === transSign) opts.reverse();
        return opts.map((x) => [x]);
      }
      if (m === 2) {
        let A = th + (2 * Math.PI) / 3, B = th - (2 * Math.PI) / 3;
        if (transSign !== 0 && this.sideOf(a, u, B) === transSign) [A, B] = [B, A];
        return [[A, B], [B, A]];
      }
      if (m === 3) {
        let p90 = th + Math.PI / 2, m90 = th - Math.PI / 2;
        if (transSign !== 0 && this.sideOf(a, u, m90) === transSign) [p90, m90] = [m90, p90];
        const dirs = [p90, th + Math.PI, m90];
        for (const perm of [[0, 1, 2], [1, 0, 2], [0, 2, 1], [2, 1, 0], [1, 2, 0], [2, 0, 1]]) out.push(perm.map((k) => dirs[k]));
        return out;
      }
      const base = kids.map((_, k) => th + ((k + 1) * TAU) / (m + 1));
      return [base, [...base].reverse()];
    }
    // two or more placed neighbours: spread inside the best free gap
    const s = this.atomSys[u];
    let away: Pt | undefined;
    if (s >= 0) {
      let cx = 0, cy = 0, k = 0;
      for (const r of this.info.systems[s]) for (const a of this.info.rings[r]) if (this.placed[a]) { cx += this.xy[2 * a]; cy += this.xy[2 * a + 1]; k++; }
      if (k) away = { x: cx / k, y: cy / k };
    }
    const g = this.bestGap(u, up, placedNb, (i) => this.P(i), away);
    let base: number[];
    if (m === 1) base = [g.start + g.size / 2];
    else {
      // keep the substituents ≤ 120° apart from each other when the gap is wide
      const step = Math.min(g.size / (m + 1), (2 * Math.PI) / 3);
      const mid = g.start + g.size / 2;
      base = kids.map((_, k) => mid + (k - (m - 1) / 2) * step);
    }
    return [base, [...base].reverse()];
  }

  /** Side (−1/0/+1) of line a→u on which a "trans" continuation from u must lie. */
  private transSign(u: number, a: number, up: Pt): number {
    // reference: the atom before a (walking back over collinear sp atoms)
    let ref = this.parent[a];
    let guard = 0;
    while (ref >= 0 && Math.abs(cross3(this.P(a), up, this.P(ref))) < 1e-6 && guard++ < 50) ref = this.parent[ref];
    if (ref < 0 || ref === u) return 0;
    const s = cross3(this.P(a), up, this.P(ref));
    return s > 0 ? -1 : s < 0 ? 1 : 0;
  }

  private sideOf(a: number, u: number, ang: number): number {
    const up = this.P(u);
    const s = cross3(this.P(a), up, { x: up.x + Math.cos(ang), y: up.y + Math.sin(ang) });
    return s > 0 ? 1 : s < 0 ? -1 : 0;
  }

  // ───────────────────────── clean-up, orientation, arrangement ─────────────────────────

  private resolve(comp: number[]): void {
    const movable = new Uint8Array(this.n);
    for (const a of comp) movable[a] = this.fixed && this.fixed.has(a) ? 0 : 1;
    const stereoOK = (moved: number[]) => {
      const seen = new Set<DbSpec>();
      for (const a of moved) {
        for (const d of this.specsByAtom.get(a) ?? []) {
          if (seen.has(d)) continue;
          seen.add(d);
          if (!dbSatisfied(this.mol, this.xy, d)) return false;
        }
      }
      return true;
    };
    resolveOverlaps(this.mol, this.xy, this.info, {
      atoms: comp,
      movable: movable,
      stereoOK,
      ringSystems: this.ringSystemsOf(comp),
    });
    this.rebuildGrid();
  }

  /**
   * Unsubstituted single ring: heteroatoms at the bottom, otherwise odd rings with a vertex up,
   * six-membered rings with vertical sides and other even rings with a flat top. Returns the
   * rotation to apply, or null if the component is not a lone ring.
   */
  private loneRingRotation(comp: number[], cx: number, cy: number): number | null {
    if (this.info.rings.length === 0 || comp.some((a) => !this.info.inRing[a])) return null;
    const ringsHere = this.info.rings.filter((r) => comp.includes(r[0]));
    if (ringsHere.length !== 1) return null;
    const ring = ringsHere[0];
    const n = ring.length;
    const xy = this.xy;
    const hetero = ring.filter((a) => this.mol.atoms[a].el !== 'C');
    const ref = hetero.length ? hetero : [ring[0]];
    let hx = 0, hy = 0;
    for (const a of ref) { hx += xy[2 * a] - cx; hy += xy[2 * a + 1] - cy; }
    let rot: number;
    if (hetero.length && Math.hypot(hx, hy) > 1e-6) {
      // heteroatom centroid straight down, then snap so that a vertex/edge is symmetric
      rot = Math.PI / 2 - Math.atan2(hy, hx);
      const step = TAU / n;
      const a0 = Math.atan2(xy[2 * hetero[0] + 1] - cy, xy[2 * hetero[0]] - cx) + rot;
      const k = Math.round((a0 - Math.PI / 2) / (step / 2));
      rot -= a0 - Math.PI / 2 - (k * step) / 2;
    } else {
      const a0 = Math.atan2(xy[2 * ring[0] + 1] - cy, xy[2 * ring[0]] - cx);
      // vertex directions wanted: top vertex (odd, 6), or flat top edge (4, 8, …)
      const want = n % 2 === 1 || n === 6 ? -Math.PI / 2 : -Math.PI / 2 + Math.PI / n;
      rot = want - a0;
    }
    return rot;
  }

  /** Index of the largest ring system in the component (−1 if acyclic). */
  private coreSystem(comp: number[]): number {
    let best = -1, bestN = 0;
    const seen = new Set<number>();
    for (const a of comp) {
      const s = this.atomSys[a];
      if (s < 0 || seen.has(s)) continue;
      seen.add(s);
      const nAt = new Set(this.info.systems[s].flatMap((r) => this.info.rings[r])).size;
      if (nAt > bestN) { bestN = nAt; best = s; }
    }
    return best;
  }

  /** For a steroid-like (6-6-6-5 fused, 17 atoms) ring system, the atoms of the five-membered ring. */
  private steroidFiveRing(s: number): number[] | null {
    const rings = this.info.systems[s].map((r) => this.info.rings[r]);
    if (rings.length !== 4) return null;
    const sizes = rings.map((r) => r.length).sort();
    if (sizes.join() !== '5,6,6,6') return null;
    if (new Set(rings.flat()).size !== 17) return null;
    return rings.find((r) => r.length === 5)!;
  }

  /** Rotates a free component: long axis horizontal, bonds snapped to the hexagonal lattice. */
  private orient(comp: number[]): void {
    if (comp.length < 2) return;
    const xy = this.xy;
    // the core (largest) ring system decides a few special cases
    const core = this.coreSystem(comp);
    const coreMethod = core >= 0 ? this.sysLayout(core).method : null;
    const keepRotation = coreMethod === 'template' || coreMethod === 'projected';
    let cx = 0, cy = 0;
    for (const a of comp) { cx += xy[2 * a]; cy += xy[2 * a + 1]; }
    cx /= comp.length; cy /= comp.length;
    let sxx = 0, syy = 0, sxy = 0;
    for (const a of comp) {
      const dx = xy[2 * a] - cx, dy = xy[2 * a + 1] - cy;
      sxx += dx * dx; syy += dy * dy; sxy += dx * dy;
    }
    let phi = 0.5 * Math.atan2(2 * sxy, sxx - syy);
    const tr = (sxx + syy) / 2, det = Math.sqrt(((sxx - syy) / 2) ** 2 + sxy * sxy);
    const l1 = tr + det, l2 = tr - det;
    let ratio = l2 > 1e-9 ? l1 / l2 : 1e9;
    if (ratio < 1.3) {
      // nearly isotropic (small or round molecules): align the longest chain instead
      const [pa, pb] = longestPath(this.mol, comp);
      if (pa !== pb) {
        phi = Math.atan2(xy[2 * pb + 1] - xy[2 * pa + 1], xy[2 * pb] - xy[2 * pa]);
        ratio = 1.6;
      }
    }
    const bondAngles: number[] = [];
    const inComp = new Set(comp);
    for (const b of this.mol.bonds) {
      if (!inComp.has(b.a)) continue;
      bondAngles.push(Math.atan2(xy[2 * b.b + 1] - xy[2 * b.a + 1], xy[2 * b.b] - xy[2 * b.a]));
    }
    // search a rotation near −phi maximising alignment with 30° + k·60° directions
    // (the lattice repeats every 60°, so ±30° covers every alignment; elongated molecules
    // keep a mild preference for their principal axis)
    let bestRot = -phi, bestS = -Infinity;
    for (let d = -30; d <= 30; d += 0.5) {
      const rot = -phi + (d * Math.PI) / 180;
      let s = 0;
      for (const t of bondAngles) s -= Math.cos(6 * (t + rot));
      s /= bondAngles.length || 1;
      s -= ratio > 1.15 ? 0.12 * Math.min(1, (ratio - 1.15) * 2) * (d / 30) ** 2 : 0;
      if (s > bestS + 1e-9) { bestS = s; bestRot = rot; }
    }
    if (keepRotation) bestRot = 0; // perspective cage drawings keep their designed orientation
    const loneRing = this.loneRingRotation(comp, cx, cy);
    if (loneRing !== null) bestRot = loneRing;
    const c = Math.cos(bestRot), s = Math.sin(bestRot);
    for (const a of comp) {
      const dx = xy[2 * a] - cx, dy = xy[2 * a + 1] - cy;
      xy[2 * a] = c * dx - s * dy;
      xy[2 * a + 1] = s * dx + c * dy;
    }
    // re-centre on the bounding box
    let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    for (const a of comp) {
      minX = Math.min(minX, xy[2 * a]); maxX = Math.max(maxX, xy[2 * a]);
      minY = Math.min(minY, xy[2 * a + 1]); maxY = Math.max(maxY, xy[2 * a + 1]);
    }
    const mx = (minX + maxX) / 2, my = (minY + maxY) / 2;
    for (const a of comp) { xy[2 * a] -= mx; xy[2 * a + 1] -= my; }
    // mirror choices (E/Z is mirror invariant; wedges are assigned afterwards):
    // ring systems on the left (chains to the right), otherwise the first atom on the left;
    // heteroatoms preferably pointing up.
    const first = comp[0];
    const w = Math.max(maxX - minX, 1e-6), h = Math.max(maxY - minY, 1e-6);
    let flipX = false, flipY = false;
    let rx = 0, rn = 0;
    for (const a of comp) if (this.info.inRing[a]) { rx += xy[2 * a]; rn++; }
    if (rn && rn < comp.length && Math.abs(rx / rn) > 0.08 * w) flipX = rx / rn > 0;
    else if (Math.abs(xy[2 * first]) > 0.05 * w) flipX = xy[2 * first] > 0;
    // exocyclic heteroatoms (C=O, OH, NH2) up, ring heteroatoms down (indole NH, pyridine N…)
    let hy = 0, hn = 0;
    for (const a of comp) {
      const el = this.mol.atoms[a].el;
      if (el === 'C' || el === 'H') continue;
      const wgt = this.info.inRing[a] ? -1.5 : 1;
      hy += wgt * xy[2 * a + 1];
      hn += Math.abs(wgt);
    }
    if (hn && hy / hn > 0.02 * h) flipY = true;
    else if (!hn && xy[2 * first + 1] < -0.02 * h) flipY = true;
    if (keepRotation) flipY = false;
    const ringD = core >= 0 ? this.steroidFiveRing(core) : null;
    if (ringD) {
      // steroids: conventional orientation with ring D at the upper right
      const sysAtoms = [...this.sysLayout(core).pos.keys()];
      const cen = (l: number[]) => {
        let x = 0, y = 0;
        for (const a of l) { x += xy[2 * a]; y += xy[2 * a + 1]; }
        return { x: x / l.length, y: y / l.length };
      };
      const cs = cen(sysAtoms), cd = cen(ringD);
      flipX = cd.x < cs.x;
      flipY = cd.y > cs.y;
    }
    if (flipX || flipY) {
      for (const a of comp) {
        if (flipX) xy[2 * a] = -xy[2 * a];
        if (flipY) xy[2 * a + 1] = -xy[2 * a + 1];
      }
    }
  }

  /** Arranges free components left to right next to the anchored ones; counter-ions beside partners. */
  private arrange(free: number[][], anchored: number[][]): void {
    const xy = this.xy;
    const bbox = (comp: number[]) => {
      let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
      for (const a of comp) {
        minX = Math.min(minX, xy[2 * a]); maxX = Math.max(maxX, xy[2 * a]);
        minY = Math.min(minY, xy[2 * a + 1]); maxY = Math.max(maxY, xy[2 * a + 1]);
      }
      return { minX, maxX, minY, maxY };
    };
    const move = (comp: number[], dx: number, dy: number) => {
      for (const a of comp) { xy[2 * a] += dx; xy[2 * a + 1] += dy; }
    };
    // counter-ions: single (heavy) atom ions with a partner of opposite charge
    const isIon = (comp: number[]) => {
      const heavy = comp.filter((a) => this.mol.atoms[a].el !== 'H');
      return heavy.length === 1 && this.mol.atoms[heavy[0]].charge !== 0;
    };
    const ions = free.filter((c) => isIon(c) && free.length + anchored.length > 1);
    let rowComps = free.filter((c) => !ions.includes(c));
    if (!rowComps.length && ions.length) {
      rowComps = [ions[0]];
      ions.shift();
    }
    // row layout
    let x = 0, cy = 0;
    if (anchored.length) {
      const all = anchored.flat();
      const bb = bbox(all);
      x = bb.maxX + GAP;
      cy = (bb.minY + bb.maxY) / 2;
    }
    const rowStart = x;
    for (const c of rowComps) {
      const bb = bbox(c);
      move(c, x - bb.minX, cy - (bb.minY + bb.maxY) / 2);
      x += bb.maxX - bb.minX + GAP;
    }
    if (!anchored.length && rowComps.length > 1) {
      // centre the row on the origin
      const width = x - GAP - rowStart;
      for (const c of rowComps) move(c, -width / 2, 0);
    } else if (!anchored.length && rowComps.length === 1) {
      const bb = bbox(rowComps[0]);
      move(rowComps[0], -(bb.minX + bb.maxX) / 2, -(bb.minY + bb.maxY) / 2);
    }
    // counter-ions
    const placedAtoms = [...rowComps.flat(), ...anchored.flat()];
    const used = new Set<number>();
    for (const ion of ions) {
      const heavy = ion.find((a) => this.mol.atoms[a].el !== 'H')!;
      const q = this.mol.atoms[heavy].charge;
      const partners = placedAtoms.filter((a) => this.mol.atoms[a].charge * q < 0);
      partners.sort((a, b) => (used.has(a) ? 1 : 0) - (used.has(b) ? 1 : 0) || this.mol.degree(a) - this.mol.degree(b) || a - b);
      let done = false;
      for (const p of partners) {
        const pp = this.P(p);
        const nb = this.mol.neighbors(p);
        const g = nb.length ? this.bestGap(p, pp, nb, (i) => this.P(i)) : { start: -Math.PI / 2, size: TAU };
        const base = nb.length ? g.start + g.size / 2 : 0;
        // prefer "O⁻ Na⁺" side by side (ion to the right), then the free direction of the partner
        const angs = [0, ...[0, 30, -30, 60, -60, 90, -90, 180].map((o) => base + (o * Math.PI) / 180)];
        for (const ang of angs) {
          const target = { x: pp.x + 1.5 * Math.cos(ang), y: pp.y + 1.5 * Math.sin(ang) };
          let ok = true;
          for (const a of placedAtoms) {
            if (a === p) continue;
            if (Math.hypot(xy[2 * a] - target.x, xy[2 * a + 1] - target.y) < 1.3) { ok = false; break; }
          }
          if (!ok) continue;
          const hp = this.P(heavy);
          move(ion, target.x - hp.x, target.y - hp.y);
          placedAtoms.push(...ion);
          used.add(p);
          done = true;
          break;
        }
        if (done) break;
      }
      if (!done) {
        const all = bbox(placedAtoms.length ? placedAtoms : ion);
        const bb = bbox(ion);
        move(ion, all.maxX + GAP - bb.minX, (all.minY + all.maxY) / 2 - (bb.minY + bb.maxY) / 2);
        placedAtoms.push(...ion);
      }
    }
  }
}

/** One end of the longest (heavy-atom) path of an acyclic component (double BFS). */
function longestPathEnd(mol: Mol, comp: number[]): number {
  const [a, b] = longestPath(mol, comp);
  return Math.min(a, b);
}

/** Both ends of a longest heavy-atom path (double BFS; exact for trees, a heuristic otherwise). */
function longestPath(mol: Mol, comp: number[]): [number, number] {
  const heavy = comp.filter((a) => mol.atoms[a].el !== 'H');
  const start = heavy.length ? heavy[0] : comp[0];
  const far = (s: number): number => {
    const d = new Map<number, number>([[s, 0]]);
    const q = [s];
    let last = s;
    for (let k = 0; k < q.length; k++) {
      const v = q[k];
      for (const w of mol.neighbors(v)) {
        if (d.has(w) || (mol.atoms[w].el === 'H' && heavy.length)) continue;
        d.set(w, d.get(v)! + 1);
        q.push(w);
        if (d.get(w)! > d.get(last)! || (d.get(w)! === d.get(last)! && w < last)) last = w;
      }
    }
    return last;
  };
  const a = far(start);
  const b = far(a);
  return [a, b];
}
