// Stereochemistry from 2D drawings (wedges + coordinates) and the reverse (wedges from specs).
import { Mol, TetraSpec, DbSpec } from './mol';
import { implicitH } from './valence';
import { perceiveRings, RingInfo, smallestRingSizeOfBond } from './rings';
import { symmetryClasses } from './canon';

type V3 = [number, number, number];

const sub = (a: V3, b: V3): V3 => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const cross = (a: V3, b: V3): V3 => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const dot = (a: V3, b: V3) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

/**
 * Signed handedness of four neighbour positions (right-handed, y up).
 * Returns true if, viewed from p0 toward the centre, p1→p2→p3 run counter-clockwise.
 */
export function ccwFromPositions(p0: V3, p1: V3, p2: V3, p3: V3): boolean {
  const s = dot(sub(p1, p0), cross(sub(p2, p0), sub(p3, p0)));
  return s < 0;
}

const STEREO_ELEMENTS = new Set(['C', 'Si', 'Ge', 'Sn', 'N', 'P', 'As', 'S', 'Se', 'B']);

/** Atoms that can be tetrahedral stereocentres (ignores symmetry). */
export function isPotentialStereocenter(mol: Mol, i: number): boolean {
  const a = mol.atoms[i];
  if (!STEREO_ELEMENTS.has(a.el) || a.abbrev) return false;
  const h = implicitH(mol, i);
  const deg = mol.degree(i);
  const total = deg + h;
  for (const bi of mol.adj[i]) {
    const o = mol.bonds[bi].order;
    if (o !== 1) {
      // sulfoxides / phosphine oxides with S=O / P=O are stereogenic
      if (!((a.el === 'S' || a.el === 'P' || a.el === 'Se') && o === 2)) return false;
    }
  }
  if (h > 1) return false;
  if (total === 4) return true;
  if (total === 3 && (a.el === 'S' || a.el === 'Se' || a.el === 'P' || (a.el === 'N' && a.charge === 0 && false))) return true;
  return false;
}

/**
 * Derives tetrahedral and double-bond stereo specs from wedge/hash bonds and 2D coordinates.
 * Wedges are interpreted with their narrow end (bond.a) at the stereocentre.
 * Results are written to mol.tetra / mol.dbStereo (replacing previous content) and returned.
 */
export function perceiveStereo2D(mol: Mol, opts: { rings?: RingInfo; requireAsymmetry?: boolean } = {}): { tetra: TetraSpec[]; db: DbSpec[] } {
  const rings = opts.rings ?? perceiveRings(mol);
  const tetra: TetraSpec[] = [];
  const db: DbSpec[] = [];
  const sym = opts.requireAsymmetry === false ? null : symmetryClasses(mol);

  for (let i = 0; i < mol.atoms.length; i++) {
    let hasStereoBond = false;
    for (const bi of mol.adj[i]) {
      const b = mol.bonds[bi];
      if ((b.style === 'wedge' || b.style === 'hash') && b.a === i) hasStereoBond = true;
    }
    if (!hasStereoBond) continue;
    if (!isPotentialStereocenter(mol, i)) continue;
    const spec = tetraFromWedges(mol, i);
    if (spec) tetra.push(spec);
  }

  // Double bonds
  mol.bonds.forEach((b, bi) => {
    if (b.order !== 2 || b.style === 'crossed') return;
    const ringSize = smallestRingSizeOfBond(rings, bi);
    if (ringSize > 0 && ringSize < 8) return;
    const ends = [b.a, b.b];
    const refs: number[] = [];
    for (const e of ends) {
      const other = e === b.a ? b.b : b.a;
      const nb = mol.neighbors(e).filter((w) => w !== other);
      const h = implicitH(mol, e);
      if (nb.length === 0 || nb.length + h > 2) return;
      // wavy bond attached → unspecified
      for (const bj of mol.adj[e]) if (mol.bonds[bj].style === 'wavy') return;
      if (nb.length === 2 && sym && sym[nb[0]] === sym[nb[1]]) return;
      if (nb.length === 1 && h === 1 && mol.atoms[nb[0]].el === 'H') return;
      if (nb.length === 2 && h === 0) {
        // two substituents: if one is explicit H and other not, fine
      }
      refs.push(nb[0]);
    }
    if (refs.length !== 2) return;
    const A = mol.atoms[b.a], B = mol.atoms[b.b];
    const ra = mol.atoms[refs[0]], rb = mol.atoms[refs[1]];
    const dx = B.x - A.x, dy = B.y - A.y;
    const sa = dx * (ra.y - A.y) - dy * (ra.x - A.x);
    const sb = dx * (rb.y - A.y) - dy * (rb.x - A.x);
    if (Math.abs(sa) < 1e-3 || Math.abs(sb) < 1e-3) return; // linear – undefined
    db.push({ bond: bi, a: refs[0], b: refs[1], cis: sa * sb > 0 });
  });

  mol.tetra = tetra;
  mol.dbStereo = db;
  return { tetra, db };
}

/** Builds a TetraSpec for centre i from wedge/hash bonds and coordinates. */
export function tetraFromWedges(mol: Mol, i: number): TetraSpec | null {
  const c = mol.atoms[i];
  const nbrs = mol.neighbors(i);
  const pos: V3[] = [];
  for (const j of nbrs) {
    const bi = mol.bondBetween(i, j);
    const b = mol.bonds[bi];
    let z = 0;
    if (b.a === i) {
      if (b.style === 'wedge') z = 1;
      else if (b.style === 'hash') z = -1;
    }
    const a = mol.atoms[j];
    let vx = a.x - c.x;
    let vy = -(a.y - c.y); // flip to y-up
    const len = Math.hypot(vx, vy) || 1;
    vx /= len;
    vy /= len;
    pos.push([vx, vy, z * 0.8]);
  }
  const order = [...nbrs];
  if (nbrs.length === 3) {
    // implicit H / lone pair: opposite to the sum of the others, z opposite of mean
    const sx = pos[0][0] + pos[1][0] + pos[2][0];
    const sy = pos[0][1] + pos[1][1] + pos[2][1];
    const sz = pos[0][2] + pos[1][2] + pos[2][2];
    let hx = -sx, hy = -sy;
    const l = Math.hypot(hx, hy);
    if (l < 1e-3) {
      hx = 0; hy = 0;
    } else {
      hx /= l; hy /= l;
    }
    // if all three are in plane and one wedge sets z, H goes opposite
    const hz = sz === 0 ? 0 : -Math.sign(sz) * 0.8;
    let hp: V3 = [hx * 0.5, hy * 0.5, hz];
    if (l < 1e-3 && hz === 0) return null;
    if (l < 1e-3) hp = [0, 0, hz];
    pos.push(hp);
    order.push(-1);
  } else if (nbrs.length !== 4) {
    return null;
  } else {
    // Four neighbours: if no z info differentiates, undefined
    if (pos.every((p) => p[2] === 0)) return null;
    // Perspective fix: when a single wedge is drawn and the other three are in-plane,
    // keep as is – the signed volume handles it.
  }
  const ccw = ccwFromPositions(pos[0], pos[1], pos[2], pos[3]);
  const vol = Math.abs(dot(sub(pos[1], pos[0]), cross(sub(pos[2], pos[0]), sub(pos[3], pos[0]))));
  if (vol < 1e-4) return null;
  return { center: i, nbrs: order as TetraSpec['nbrs'], ccw };
}

/**
 * Sets wedge/hash bonds so that the 2D drawing reproduces mol.tetra.
 * Existing stereo bond styles at those centres are cleared first.
 */
export function assignWedgesFromSpecs(mol: Mol, rings?: RingInfo): void {
  const info = rings ?? perceiveRings(mol);
  const centers = new Set(mol.tetra.map((t) => t.center));
  // clear existing wedges touching stereocentres
  for (const b of mol.bonds) {
    if ((b.style === 'wedge' || b.style === 'hash') && (centers.has(b.a) || centers.has(b.b))) b.style = 'plain';
  }
  const usedBond = new Set<number>();
  for (const t of mol.tetra) {
    const i = t.center;
    // candidate bonds ranked by preference
    const cands = mol.adj[i]
      .filter((bi) => mol.bonds[bi].order === 1 && !usedBond.has(bi))
      .map((bi) => {
        const j = mol.other(bi, i);
        let score = 0;
        if (!info.bondInRing[bi]) score += 4;
        if (mol.atoms[j].el === 'H') score += 3;
        if (mol.degree(j) === 1) score += 2;
        if (centers.has(j)) score -= 5; // avoid wedges shared by two centres
        if (mol.atoms[j].el !== 'C') score += 0.5;
        return { bi, score };
      })
      .sort((p, q) => q.score - p.score);
    if (!cands.length) continue;
    const bi = cands[0].bi;
    const b = mol.bonds[bi];
    if (b.a !== i) {
      const tmp = b.a;
      b.a = b.b;
      b.b = tmp;
    }
    b.style = 'wedge';
    const spec = tetraFromWedges(mol, i);
    if (!spec) {
      b.style = 'plain';
      continue;
    }
    const order = spec.nbrs;
    // compare chirality: map t to spec ordering
    const want = sameOrderCcw(t, order);
    if (want === null) {
      b.style = 'plain';
      continue;
    }
    if (want !== spec.ccw) b.style = 'hash';
    usedBond.add(bi);
  }
}

function sameOrderCcw(t: TetraSpec, order: readonly number[]): boolean | null {
  // both contain the same set of neighbours (with -1 for implicit H)
  const a = [...t.nbrs].sort((p, q) => p - q).join(',');
  const b = [...order].sort((p, q) => p - q).join(',');
  if (a !== b) return null;
  const perm = order.map((v) => t.nbrs.indexOf(v));
  let even = true;
  const seen = new Array(4).fill(false);
  for (let k = 0; k < 4; k++) {
    if (seen[k]) continue;
    let len = 0;
    let j = k;
    while (!seen[j]) {
      seen[j] = true;
      j = perm[j];
      len++;
    }
    if (len % 2 === 0) even = !even;
  }
  return even ? t.ccw : !t.ccw;
}

/**
 * Places double-bond substituents so the drawing matches mol.dbStereo (used by layout).
 * Returns the specs that are currently NOT satisfied by the coordinates.
 */
export function unsatisfiedDbStereo(mol: Mol): DbSpec[] {
  const bad: DbSpec[] = [];
  for (const d of mol.dbStereo) {
    const b = mol.bonds[d.bond];
    const A = mol.atoms[b.a], B = mol.atoms[b.b];
    const ra = mol.atoms[d.a], rb = mol.atoms[d.b];
    const dx = B.x - A.x, dy = B.y - A.y;
    const sa = dx * (ra.y - A.y) - dy * (ra.x - A.x);
    const sb = dx * (rb.y - A.y) - dy * (rb.x - A.x);
    const cis = sa * sb > 0;
    if (cis !== d.cis) bad.push(d);
  }
  return bad;
}
