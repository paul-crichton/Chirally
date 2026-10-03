// UFF force field: energy and analytic gradient (kcal/mol, Å).
// Functional forms follow Rappé et al. (JACS 1992) as implemented in common open-source UFF codes:
// harmonic bond stretch, cosine-Fourier angle bend, cosine torsions, Wilson out-of-plane inversion,
// Lennard-Jones 12-6 van der Waals for pairs separated by more than two bonds, optional Coulomb.
import { Mol } from '../mol';
import { typeUFF, UFFTyping, uffBondLength, isGroup6 } from './uffTyping';
import { gasteigerCharges } from './charges';

export type TermName = 'stretch' | 'bend' | 'torsion' | 'oop' | 'vdw' | 'elec';
export type TermEnergies = Record<TermName, number>;

export interface UFFOptions {
  /** Include Coulomb electrostatics with Gasteiger charges (default false, as in standard UFF). */
  electrostatics?: boolean;
  /** Relative dielectric constant (default 1, or 4 with distanceDependent). */
  dielectric?: number;
  /** Distance-dependent dielectric ε = dielectric·r (default false). */
  distanceDependent?: boolean;
  /** Partial charges overriding the Gasteiger charges. */
  charges?: number[];
  /** Non-bonded pairs farther apart than this (Å) when the force field is set up are ignored (default 12). */
  nonBondedCutoff?: number;
}

const G = 332.06; // UFF's energy-conversion constant (kcal·Å/mol/e²)
const COULOMB = 332.0637;

export const emptyTerms = (): TermEnergies => ({ stretch: 0, bend: 0, torsion: 0, oop: 0, vdw: 0, elec: 0 });

/** Flat [x0, y0, z0, x1, …] coordinate array from a molecule. */
export function coordsOf(mol: Mol): Float64Array {
  const x = new Float64Array(mol.atoms.length * 3);
  mol.atoms.forEach((a, i) => {
    x[3 * i] = a.x;
    x[3 * i + 1] = a.y;
    x[3 * i + 2] = a.z ?? 0;
  });
  return x;
}

export function setCoords(mol: Mol, x: Float64Array): void {
  mol.atoms.forEach((a, i) => {
    a.x = x[3 * i];
    a.y = x[3 * i + 1];
    a.z = x[3 * i + 2];
  });
}

/** UFF angle force constant (kcal/mol/rad²). */
function angleForceConstant(theta0: number, r12: number, r23: number, z1: number, z3: number): number {
  const c0 = Math.cos(theta0);
  const r13sq = r12 * r12 + r23 * r23 - 2 * r12 * r23 * c0;
  const r13 = Math.sqrt(r13sq);
  const beta = (2 * G) / (r12 * r23);
  const pre = (beta * z1 * z3) / Math.pow(r13, 5);
  const rTerm = r12 * r23;
  return pre * rTerm * (3 * rTerm * (1 - c0 * c0) - r13sq * c0);
}

/** Coordination-specific angle expansion: 0 = general Fourier, 1 = linear, 3 = trigonal, 4 = square/octahedral. */
function angleOrder(theta0deg: number, degree: number): number {
  if (degree > 4) return 4;
  if (Math.abs(theta0deg - 180) < 1) return 1;
  if (Math.abs(theta0deg - 120) < 1) return 3;
  if (Math.abs(theta0deg - 90) < 1) return 4;
  return 0;
}

export class UFF {
  readonly n: number;
  readonly typing: UFFTyping;
  // bond stretch
  private bI: Int32Array; private bJ: Int32Array; private bK: Float64Array; private bR: Float64Array;
  // angle bend
  private aI: Int32Array; private aJ: Int32Array; private aK: Int32Array;
  private aF: Float64Array; private aC0: Float64Array; private aC1: Float64Array; private aC2: Float64Array; private aO: Int8Array;
  // torsion
  private tI: Int32Array; private tJ: Int32Array; private tK: Int32Array; private tL: Int32Array;
  private tV: Float64Array; private tC: Float64Array; private tN: Int8Array;
  // inversion (oI, oK in plane with centre oJ; oL out of plane)
  private oI: Int32Array; private oJ: Int32Array; private oK: Int32Array; private oL: Int32Array;
  private oF: Float64Array; private oC0: Float64Array; private oC1: Float64Array; private oC2: Float64Array;
  // non-bonded
  private vI: Int32Array; private vJ: Int32Array; private vX2: Float64Array; private vD: Float64Array;
  private vQ: Float64Array | null; private distDep: boolean;
  /** Natural bond lengths per Mol bond (0 for non-covalent bonds). */
  readonly restLengths: Float64Array;

  constructor(mol: Mol, opts: UFFOptions = {}, x0?: Float64Array) {
    const n = (this.n = mol.atoms.length);
    const ty = (this.typing = typeUFF(mol));
    const { params, hyb, nbrs, bondOrder } = ty;
    const x = x0 ?? coordsOf(mol);

    // ── bonds
    const bonds: number[][] = [];
    const rest = (this.restLengths = new Float64Array(mol.bonds.length));
    const r0Between = new Map<number, number>();
    const key = (i: number, j: number) => (i < j ? i * n + j : j * n + i);
    mol.bonds.forEach((b, bi) => {
      const bo = bondOrder[bi];
      if (bo <= 0 || b.a === b.b) return;
      const r0 = uffBondLength(bo, params[b.a], params[b.b]);
      const k = (2 * G * params[b.a].Z1 * params[b.b].Z1) / (r0 * r0 * r0);
      rest[bi] = r0;
      r0Between.set(key(b.a, b.b), r0);
      bonds.push([b.a, b.b, k, r0]);
    });
    this.bI = Int32Array.from(bonds, (t) => t[0]);
    this.bJ = Int32Array.from(bonds, (t) => t[1]);
    this.bK = Float64Array.from(bonds, (t) => t[2]);
    this.bR = Float64Array.from(bonds, (t) => t[3]);
    const r0Of = (i: number, j: number) => r0Between.get(key(i, j)) ?? params[i].r1 + params[j].r1;

    // ── angles
    const angles: number[][] = [];
    for (let j = 0; j < n; j++) {
      const nb = nbrs[j];
      const deg = nb.length;
      if (deg < 2) continue;
      let theta0deg = params[j].theta0;
      if (deg > 4) theta0deg = 90;
      const order = angleOrder(theta0deg, deg);
      const th0 = (theta0deg * Math.PI) / 180;
      for (let p = 0; p < deg; p++) {
        for (let q = p + 1; q < deg; q++) {
          const i = nb[p], k = nb[q];
          const kf = angleForceConstant(th0, r0Of(i, j), r0Of(j, k), params[i].Z1, params[k].Z1);
          const c = Math.cos(th0);
          const s2 = Math.max(1 - c * c, 1e-8);
          const C2 = 1 / (4 * s2);
          const C1 = -4 * C2 * c;
          const C0 = C2 * (2 * c * c + 1);
          angles.push([i, j, k, kf, C0, C1, C2, order]);
        }
      }
    }
    this.aI = Int32Array.from(angles, (t) => t[0]);
    this.aJ = Int32Array.from(angles, (t) => t[1]);
    this.aK = Int32Array.from(angles, (t) => t[2]);
    this.aF = Float64Array.from(angles, (t) => t[3]);
    this.aC0 = Float64Array.from(angles, (t) => t[4]);
    this.aC1 = Float64Array.from(angles, (t) => t[5]);
    this.aC2 = Float64Array.from(angles, (t) => t[6]);
    this.aO = Int8Array.from(angles, (t) => t[7]);

    // ── torsions
    const tors: number[][] = [];
    const eq17 = (bo: number, j: number, k: number) => 5 * Math.sqrt(params[j].Uj * params[k].Uj) * (1 + 4.18 * Math.log(bo));
    mol.bonds.forEach((b, bi) => {
      const bo = bondOrder[bi];
      if (bo <= 0) return;
      const j = b.a, k = b.b;
      const hj = hyb[j], hk = hyb[k];
      if ((hj !== 2 && hj !== 3) || (hk !== 2 && hk !== 3)) return;
      const nj = nbrs[j].length - 1, nk = nbrs[k].length - 1;
      if (nj < 1 || nk < 1) return;
      const elJ = mol.atoms[j].el, elK = mol.atoms[k].el;
      for (const i of nbrs[j]) {
        if (i === k) continue;
        for (const l of nbrs[k]) {
          if (l === j || l === i) continue;
          let V: number, order: number, cosTerm: number;
          if (hj === 3 && hk === 3) {
            V = Math.sqrt(params[j].Vi * params[k].Vi);
            order = 3;
            cosTerm = -1; // φ0 = 60°
            if (bo === 1 && isGroup6(elJ) && isGroup6(elK)) {
              const vj = elJ === 'O' ? 2 : 6.8, vk = elK === 'O' ? 2 : 6.8;
              V = Math.sqrt(vj * vk);
              order = 2;
              cosTerm = -1; // φ0 = 90°
            }
          } else if (hj === 2 && hk === 2) {
            V = eq17(bo, j, k);
            order = 2;
            cosTerm = 1; // φ0 = 180°
          } else {
            // sp2–sp3
            V = 1;
            order = 6;
            cosTerm = 1; // φ0 = 0°
            if (bo === 1) {
              const sp3 = hj === 3 ? j : k, sp2 = hj === 3 ? k : j;
              const endOnSp2 = sp2 === j ? i : l;
              if (isGroup6(mol.atoms[sp3].el) && !isGroup6(mol.atoms[sp2].el)) {
                V = eq17(bo, j, k);
                order = 2;
                cosTerm = -1; // φ0 = 90°
              } else if (hyb[endOnSp2] === 2) {
                V = 2;
                order = 3;
                cosTerm = -1; // φ0 = 180° (propene-like)
              }
            }
          }
          V /= nj * nk;
          if (V > 1e-6) tors.push([i, j, k, l, V, cosTerm, order]);
        }
      }
    });
    this.tI = Int32Array.from(tors, (t) => t[0]);
    this.tJ = Int32Array.from(tors, (t) => t[1]);
    this.tK = Int32Array.from(tors, (t) => t[2]);
    this.tL = Int32Array.from(tors, (t) => t[3]);
    this.tV = Float64Array.from(tors, (t) => t[4]);
    this.tC = Float64Array.from(tors, (t) => t[5]);
    this.tN = Int8Array.from(tors, (t) => t[6]);

    // ── inversions (sp2 C/N/O planarity; pyramidal group-15 centres)
    const invs: number[][] = [];
    const W0: Record<string, number> = { P: 84.4339, As: 86.9735, Sb: 87.7047, Bi: 90 };
    for (let j = 0; j < n; j++) {
      if (nbrs[j].length !== 3) continue;
      const el = mol.atoms[j].el;
      const lab = ty.labels[j];
      let K: number, C0: number, C1: number, C2: number;
      if (lab === 'C_2' || lab === 'C_R' || lab === 'N_2' || lab === 'N_R' || lab === 'O_2' || lab === 'O_R') {
        C0 = 1; C1 = -1; C2 = 0;
        const boundToO2 = el === 'C' && nbrs[j].some((o) => ty.labels[o] === 'O_2');
        K = boundToO2 ? 50 : 6;
      } else if (W0[el] !== undefined && hyb[j] === 3) {
        const w0 = (W0[el] * Math.PI) / 180;
        C2 = 1;
        C1 = -4 * Math.cos(w0);
        C0 = -(C1 * Math.cos(w0) + C2 * Math.cos(2 * w0));
        K = 22 / (C0 + C1 + C2);
      } else continue;
      K /= 3;
      const [a, b, c] = nbrs[j];
      invs.push([a, j, b, c, K, C0, C1, C2], [a, j, c, b, K, C0, C1, C2], [b, j, c, a, K, C0, C1, C2]);
    }
    this.oI = Int32Array.from(invs, (t) => t[0]);
    this.oJ = Int32Array.from(invs, (t) => t[1]);
    this.oK = Int32Array.from(invs, (t) => t[2]);
    this.oL = Int32Array.from(invs, (t) => t[3]);
    this.oF = Float64Array.from(invs, (t) => t[4]);
    this.oC0 = Float64Array.from(invs, (t) => t[5]);
    this.oC1 = Float64Array.from(invs, (t) => t[6]);
    this.oC2 = Float64Array.from(invs, (t) => t[7]);

    // ── non-bonded pairs: separated by ≥ 3 bonds and within the cutoff
    const cutoff = opts.nonBondedCutoff ?? 12;
    const cut2 = cutoff * cutoff;
    const charges = opts.electrostatics ? opts.charges ?? gasteigerCharges(mol, hyb, nbrs) : null;
    const eps = opts.dielectric ?? (opts.distanceDependent ? 4 : 1);
    this.distDep = !!opts.distanceDependent;
    const pi: number[] = [], pj: number[] = [], px2: number[] = [], pd: number[] = [], pq: number[] = [];
    const stamp = new Int32Array(n).fill(-1);
    for (let i = 0; i < n; i++) {
      stamp[i] = i;
      for (const j of nbrs[i]) {
        stamp[j] = i;
        for (const k of nbrs[j]) stamp[k] = i;
      }
      for (let j = i + 1; j < n; j++) {
        if (stamp[j] === i) continue;
        const dx = x[3 * i] - x[3 * j], dy = x[3 * i + 1] - x[3 * j + 1], dz = x[3 * i + 2] - x[3 * j + 2];
        if (dx * dx + dy * dy + dz * dz > cut2) continue;
        pi.push(i);
        pj.push(j);
        px2.push(params[i].x1 * params[j].x1); // x_ij² with x_ij = √(x_i x_j)
        pd.push(Math.sqrt(params[i].D1 * params[j].D1));
        if (charges) pq.push((COULOMB * charges[i] * charges[j]) / eps);
      }
    }
    this.vI = Int32Array.from(pi);
    this.vJ = Int32Array.from(pj);
    this.vX2 = Float64Array.from(px2);
    this.vD = Float64Array.from(pd);
    this.vQ = charges ? Float64Array.from(pq) : null;
  }

  /** Numbers of terms of each kind (for diagnostics). */
  get counts(): Record<string, number> {
    return { bonds: this.bI.length, angles: this.aI.length, torsions: this.tI.length, inversions: this.oI.length, nonBonded: this.vI.length };
  }

  /**
   * Total energy (kcal/mol) at coordinates x (length 3n). If `grad` is given it is overwritten with ∂E/∂x;
   * if `terms` is given it receives the per-term breakdown.
   */
  energy(x: Float64Array, grad?: Float64Array | null, terms?: TermEnergies): number {
    const g = grad ?? null;
    if (g) g.fill(0);
    let eS = 0, eB = 0, eT = 0, eO = 0, eV = 0, eE = 0;

    // ── bond stretch: E = ½ k (r − r0)²
    for (let t = 0; t < this.bI.length; t++) {
      const i = 3 * this.bI[t], j = 3 * this.bJ[t];
      const dx = x[i] - x[j], dy = x[i + 1] - x[j + 1], dz = x[i + 2] - x[j + 2];
      const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
      const dr = r - this.bR[t];
      eS += 0.5 * this.bK[t] * dr * dr;
      if (g) {
        const f = (this.bK[t] * dr) / Math.max(r, 1e-10);
        g[i] += f * dx; g[i + 1] += f * dy; g[i + 2] += f * dz;
        g[j] -= f * dx; g[j + 1] -= f * dy; g[j + 2] -= f * dz;
      }
    }

    // ── angle bend: general E = K (C0 + C1 cosθ + C2 cos2θ); special E = K/n² (1 − cos nθ)
    for (let t = 0; t < this.aI.length; t++) {
      const i = 3 * this.aI[t], j = 3 * this.aJ[t], k = 3 * this.aK[t];
      const ax = x[i] - x[j], ay = x[i + 1] - x[j + 1], az = x[i + 2] - x[j + 2];
      const bx = x[k] - x[j], by = x[k + 1] - x[j + 1], bz = x[k + 2] - x[j + 2];
      const ra2 = ax * ax + ay * ay + az * az, rb2 = bx * bx + by * by + bz * bz;
      if (ra2 < 1e-16 || rb2 < 1e-16) continue;
      const ra = Math.sqrt(ra2), rb = Math.sqrt(rb2);
      let c = (ax * bx + ay * by + az * bz) / (ra * rb);
      if (c > 1) c = 1; else if (c < -1) c = -1;
      const K = this.aF[t];
      let e: number, dEdc: number;
      switch (this.aO[t]) {
        case 1: e = K * (1 + c); dEdc = K; break;
        case 2: e = (K / 4) * (2 - 2 * c * c); dEdc = -K * c; break;
        case 3: e = (K / 9) * (1 - (4 * c * c * c - 3 * c)); dEdc = (-K / 9) * (12 * c * c - 3); break;
        case 4: e = (K / 16) * (1 - (8 * c * c * c * c - 8 * c * c + 1)); dEdc = (-K / 16) * (32 * c * c * c - 16 * c); break;
        default:
          e = K * (this.aC0[t] + this.aC1[t] * c + this.aC2[t] * (2 * c * c - 1));
          dEdc = K * (this.aC1[t] + 4 * this.aC2[t] * c);
      }
      eB += e;
      if (g) {
        const inv = 1 / (ra * rb);
        const ca = c / ra2, cb = c / rb2;
        const gix = dEdc * (bx * inv - ca * ax), giy = dEdc * (by * inv - ca * ay), giz = dEdc * (bz * inv - ca * az);
        const gkx = dEdc * (ax * inv - cb * bx), gky = dEdc * (ay * inv - cb * by), gkz = dEdc * (az * inv - cb * bz);
        g[i] += gix; g[i + 1] += giy; g[i + 2] += giz;
        g[k] += gkx; g[k + 1] += gky; g[k + 2] += gkz;
        g[j] -= gix + gkx; g[j + 1] -= giy + gky; g[j + 2] -= giz + gkz;
      }
    }

    // ── torsion: E = ½ V (1 − cos(nφ0) cos(nφ)), cos(nφ) via Chebyshev polynomials of cosφ
    for (let t = 0; t < this.tI.length; t++) {
      const i = 3 * this.tI[t], j = 3 * this.tJ[t], k = 3 * this.tK[t], l = 3 * this.tL[t];
      const r1x = x[j] - x[i], r1y = x[j + 1] - x[i + 1], r1z = x[j + 2] - x[i + 2];
      const r2x = x[k] - x[j], r2y = x[k + 1] - x[j + 1], r2z = x[k + 2] - x[j + 2];
      const r3x = x[l] - x[k], r3y = x[l + 1] - x[k + 1], r3z = x[l + 2] - x[k + 2];
      // t = r1 × r2, u = r2 × r3
      const tx = r1y * r2z - r1z * r2y, ty = r1z * r2x - r1x * r2z, tz = r1x * r2y - r1y * r2x;
      const ux = r2y * r3z - r2z * r3y, uy = r2z * r3x - r2x * r3z, uz = r2x * r3y - r2y * r3x;
      const tl2 = tx * tx + ty * ty + tz * tz, ul2 = ux * ux + uy * uy + uz * uz;
      if (tl2 < 1e-12 || ul2 < 1e-12) continue;
      const tl = Math.sqrt(tl2), ul = Math.sqrt(ul2);
      let c = (tx * ux + ty * uy + tz * uz) / (tl * ul);
      if (c > 1) c = 1; else if (c < -1) c = -1;
      let T: number, dT: number;
      switch (this.tN[t]) {
        case 1: T = c; dT = 1; break;
        case 2: T = 2 * c * c - 1; dT = 4 * c; break;
        case 3: T = 4 * c * c * c - 3 * c; dT = 12 * c * c - 3; break;
        default: {
          // n = 6
          const c2 = c * c;
          T = ((32 * c2 - 48) * c2 + 18) * c2 - 1;
          dT = ((192 * c2 - 192) * c2 + 36) * c;
        }
      }
      const V = this.tV[t], ct = this.tC[t];
      eT += 0.5 * V * (1 - ct * T);
      if (g) {
        const dEdc = -0.5 * V * ct * dT;
        // ∂c/∂t and ∂c/∂u
        const gtx = (ux / ul - (c * tx) / tl) / tl, gty = (uy / ul - (c * ty) / tl) / tl, gtz = (uz / ul - (c * tz) / tl) / tl;
        const gux = (tx / tl - (c * ux) / ul) / ul, guy = (ty / tl - (c * uy) / ul) / ul, guz = (tz / tl - (c * uz) / ul) / ul;
        // ∂c/∂r1 = r2 × gt ; ∂c/∂r2 = gt × r1 + r3 × gu ; ∂c/∂r3 = gu × r2
        const d1x = r2y * gtz - r2z * gty, d1y = r2z * gtx - r2x * gtz, d1z = r2x * gty - r2y * gtx;
        const d2x = gty * r1z - gtz * r1y + (r3y * guz - r3z * guy);
        const d2y = gtz * r1x - gtx * r1z + (r3z * gux - r3x * guz);
        const d2z = gtx * r1y - gty * r1x + (r3x * guy - r3y * gux);
        const d3x = guy * r2z - guz * r2y, d3y = guz * r2x - gux * r2z, d3z = gux * r2y - guy * r2x;
        g[i] -= dEdc * d1x; g[i + 1] -= dEdc * d1y; g[i + 2] -= dEdc * d1z;
        g[j] += dEdc * (d1x - d2x); g[j + 1] += dEdc * (d1y - d2y); g[j + 2] += dEdc * (d1z - d2z);
        g[k] += dEdc * (d2x - d3x); g[k + 1] += dEdc * (d2y - d3y); g[k + 2] += dEdc * (d2z - d3z);
        g[l] += dEdc * d3x; g[l + 1] += dEdc * d3y; g[l + 2] += dEdc * d3z;
      }
    }

    // ── inversion: E = K (C0 + C1 cosω + C2 cos2ω), ω = angle between J→L and plane IJK
    for (let t = 0; t < this.oI.length; t++) {
      const i = 3 * this.oI[t], j = 3 * this.oJ[t], k = 3 * this.oK[t], l = 3 * this.oL[t];
      const ax = x[i] - x[j], ay = x[i + 1] - x[j + 1], az = x[i + 2] - x[j + 2];
      const bx = x[k] - x[j], by = x[k + 1] - x[j + 1], bz = x[k + 2] - x[j + 2];
      const cx = x[l] - x[j], cy = x[l + 1] - x[j + 1], cz = x[l + 2] - x[j + 2];
      const mx = ay * bz - az * by, my = az * bx - ax * bz, mz = ax * by - ay * bx; // m = a × b
      const ml2 = mx * mx + my * my + mz * mz, cl2 = cx * cx + cy * cy + cz * cz;
      if (ml2 < 1e-16 || cl2 < 1e-16) continue;
      const ml = Math.sqrt(ml2), cl = Math.sqrt(cl2);
      let f = (mx * cx + my * cy + mz * cz) / (ml * cl); // sin ω
      if (f > 1) f = 1; else if (f < -1) f = -1;
      const s = Math.max(Math.sqrt(1 - f * f), 1e-8); // cos ω
      const K = this.oF[t], C1 = this.oC1[t], C2 = this.oC2[t];
      eO += K * (this.oC0[t] + C1 * s + C2 * (2 * s * s - 1));
      if (g) {
        const dEdf = (K * (C1 + 4 * C2 * s) * -f) / s;
        const im = 1 / (ml * cl);
        const gmx = cx * im - (f * mx) / ml2, gmy = cy * im - (f * my) / ml2, gmz = cz * im - (f * mz) / ml2;
        const gcx = mx * im - (f * cx) / cl2, gcy = my * im - (f * cy) / cl2, gcz = mz * im - (f * cz) / cl2;
        // ∂f/∂a = b × gm ; ∂f/∂b = gm × a
        const dax = by * gmz - bz * gmy, day = bz * gmx - bx * gmz, daz = bx * gmy - by * gmx;
        const dbx = gmy * az - gmz * ay, dby = gmz * ax - gmx * az, dbz = gmx * ay - gmy * ax;
        g[i] += dEdf * dax; g[i + 1] += dEdf * day; g[i + 2] += dEdf * daz;
        g[k] += dEdf * dbx; g[k + 1] += dEdf * dby; g[k + 2] += dEdf * dbz;
        g[l] += dEdf * gcx; g[l + 1] += dEdf * gcy; g[l + 2] += dEdf * gcz;
        g[j] -= dEdf * (dax + dbx + gcx); g[j + 1] -= dEdf * (day + dby + gcy); g[j + 2] -= dEdf * (daz + dbz + gcz);
      }
    }

    // ── van der Waals (LJ 12-6) and optional Coulomb
    const vQ = this.vQ;
    for (let t = 0; t < this.vI.length; t++) {
      const i = 3 * this.vI[t], j = 3 * this.vJ[t];
      const dx = x[i] - x[j], dy = x[i + 1] - x[j + 1], dz = x[i + 2] - x[j + 2];
      const r2 = Math.max(dx * dx + dy * dy + dz * dz, 1e-6);
      const q2 = this.vX2[t] / r2;
      const q6 = q2 * q2 * q2;
      const q12 = q6 * q6;
      const D = this.vD[t];
      eV += D * (q12 - 2 * q6);
      let dEdrOverR = (12 * D * (q6 - q12)) / r2;
      if (vQ) {
        const r = Math.sqrt(r2);
        if (this.distDep) {
          eE += vQ[t] / r2;
          dEdrOverR += (-2 * vQ[t]) / (r2 * r2);
        } else {
          eE += vQ[t] / r;
          dEdrOverR += -vQ[t] / (r2 * r);
        }
      }
      if (g) {
        g[i] += dEdrOverR * dx; g[i + 1] += dEdrOverR * dy; g[i + 2] += dEdrOverR * dz;
        g[j] -= dEdrOverR * dx; g[j + 1] -= dEdrOverR * dy; g[j + 2] -= dEdrOverR * dz;
      }
    }

    if (terms) {
      terms.stretch = eS; terms.bend = eB; terms.torsion = eT; terms.oop = eO; terms.vdw = eV; terms.elec = eE;
    }
    return eS + eB + eT + eO + eV + eE;
  }
}
