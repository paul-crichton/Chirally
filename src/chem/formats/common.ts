// Helpers shared by the file-format readers/writers: errors, coordinate normalisation,
// aromatic-bond kekulisation, alias/superatom handling, stereo perception.
import { Atom, Mol, TetraSpec, DbSpec } from '../mol';
import { BY_SYMBOL, ELEMENTS } from '../elements';
import { implicitH } from '../valence';
import { kekulize } from '../aromaticity';
import { findAbbreviation, parseCondensedLabel, abbreviationMol, reverseLabel } from '../abbreviations';
import { perceiveStereo2D, isPotentialStereocenter, ccwFromPositions } from '../stereo2d';
import { perceiveRings, smallestRingSizeOfBond } from '../rings';
import { symmetryClasses } from '../canon';

/** Error thrown by all readers for malformed/unsupported input. */
export class FormatError extends Error {
  constructor(public format: string, message: string) {
    super(`${format}: ${message}`);
    this.name = 'FormatError';
  }
}

/** Runs a reader, converting unexpected internal errors into a descriptive FormatError. */
export function guard<T>(format: string, fn: () => T): T {
  try {
    return fn();
  } catch (e) {
    if (e instanceof FormatError) throw e;
    throw new FormatError(format, `unreadable input (${e instanceof Error ? e.message : String(e)})`);
  }
}

/** Splits text into lines, normalising CRLF/CR and stripping a UTF-8 BOM. */
export function splitLines(text: string): string[] {
  if (typeof text !== 'string') throw new FormatError('Input', 'expected text');
  const t = text.charCodeAt(0) === 0xfeff ? text.slice(1) : text;
  return t.replace(/\r\n?/g, '\n').split('\n');
}

/** Standard bond length (Å) assumed when writing 2D coordinates. */
export const DEFAULT_BOND_LENGTH_ANGSTROM = 1.5;

/**
 * True if the molecule carries usable coordinates, i.e. it is not a multi-atom structure whose atoms
 * all sit at the same point (files written without a layout). Callers should lay out molecules for
 * which this returns false.
 */
export function hasCoordinates(mol: Mol): boolean {
  const n = mol.atoms.length;
  if (n === 0) return false;
  if (n === 1) return true;
  const a0 = mol.atoms[0];
  for (const a of mol.atoms) {
    if (Math.abs(a.x - a0.x) > 1e-6 || Math.abs(a.y - a0.y) > 1e-6 || Math.abs((a.z ?? 0) - (a0.z ?? 0)) > 1e-6) return true;
  }
  return false;
}

export function is3D(mol: Mol): boolean {
  return mol.props.dim === '3D';
}

function median(values: number[]): number {
  const v = [...values].sort((a, b) => a - b);
  const m = v.length >> 1;
  return v.length % 2 ? v[m] : (v[m - 1] + v[m]) / 2;
}

/**
 * Converts raw file coordinates (Å, y up) of one or more molecules that share a coordinate frame
 * (e.g. the components of a reaction) to Chirally model coordinates:
 *  - 3D (z has spread): y flipped, no scaling, props.dim = '3D';
 *  - 2D: z dropped, y flipped, scaled so the median bond length becomes 1;
 *  - no coordinates (all atoms coincident): left untouched (callers lay them out).
 */
export function normalizeImportedCoords(mols: Mol[], fallbackBondLength = DEFAULT_BOND_LENGTH_ANGSTROM): void {
  let zmin = Infinity, zmax = -Infinity;
  for (const m of mols) for (const a of m.atoms) {
    const z = a.z ?? 0;
    if (z < zmin) zmin = z;
    if (z > zmax) zmax = z;
  }
  if (zmax - zmin > 1e-4) {
    for (const m of mols) {
      for (const a of m.atoms) {
        a.y = -a.y + 0; // + 0 avoids -0
        a.z = a.z ?? 0;
      }
      m.props.dim = '3D';
    }
    return;
  }
  for (const m of mols) for (const a of m.atoms) delete a.z;
  const all = new Mol();
  all.atoms = mols.flatMap((m) => m.atoms);
  if (all.atoms.length < 2 || !hasCoordinates(all)) return;
  const lengths: number[] = [];
  for (const m of mols) {
    for (const b of m.bonds) {
      const p = m.atoms[b.a], q = m.atoms[b.b];
      const d = Math.hypot(p.x - q.x, p.y - q.y);
      if (d > 1e-6) lengths.push(d);
    }
  }
  const scale = 1 / (lengths.length ? median(lengths) : fallbackBondLength);
  for (const m of mols) for (const a of m.atoms) {
    a.x = a.x * scale + 0;
    a.y = -a.y * scale + 0;
  }
}

/** Implicit H count the valence model would assign if the atom had no explicit hCount. */
export function defaultImplicitH(mol: Mol, i: number): number {
  const a = mol.atoms[i];
  const saved = a.hCount;
  delete a.hCount;
  const h = implicitH(mol, i);
  if (saved !== undefined) a.hCount = saved;
  return h;
}

/** Drops explicit hCount values that equal the computed default (keeps documents tidy). */
export function dropDefaultHCounts(mol: Mol): void {
  for (let i = 0; i < mol.atoms.length; i++) {
    const a = mol.atoms[i];
    if (a.hCount === undefined) continue;
    if (a.abbrev || a.el === 'R' || a.el === '*') {
      delete a.hCount;
      continue;
    }
    if (defaultImplicitH(mol, i) === a.hCount) delete a.hCount;
  }
}

/**
 * Turns bonds flagged aromatic (order 1.5 in the file) into alternating single/double bonds.
 * Leaves them at 1.5 if no Kekulé structure exists.
 */
export function kekulizeFlagged(mol: Mol, aromaticBonds: boolean[]): void {
  if (!aromaticBonds.some(Boolean)) return;
  const atoms = new Array(mol.atoms.length).fill(false);
  aromaticBonds.forEach((f, bi) => {
    if (f) {
      atoms[mol.bonds[bi].a] = true;
      atoms[mol.bonds[bi].b] = true;
    }
  });
  mol.invalidate();
  try {
    kekulize(mol, atoms, aromaticBonds);
  } catch {
    aromaticBonds.forEach((f, bi) => f && (mol.bonds[bi].order = 1.5));
  }
}

/** Normalises an element symbol written in any case ('CL' → 'Cl'); returns null if unknown. */
export function normalizeElement(sym: string): string | null {
  if (BY_SYMBOL.has(sym)) return sym;
  if (!sym) return null;
  const s = sym[0].toUpperCase() + sym.slice(1).toLowerCase();
  return BY_SYMBOL.has(s) ? s : null;
}

export function elementByNumber(z: number): string | null {
  return ELEMENTS[z]?.symbol ?? null;
}

/**
 * Applies a display alias (molfile 'A' record, CDXML label…) to an atom:
 *  - a label that just spells the atom's own element with hydrogens ("OH" on O) is ignored;
 *  - a known abbreviation / condensed label ("OMe", "CH2CH2OH") makes it an abbreviation atom;
 *  - anything else makes it a pseudo atom (el 'R') displaying the alias.
 */
export function applyAlias(a: Atom, alias: string): void {
  const label = alias.trim();
  if (!label) return;
  // "OH", "NH2" (or reversed "HO", "H2N") on the matching element are plain element labels
  for (const l of [label, reverseLabel(label)]) {
    const own = /^([A-Z][a-z]?)(?:H\d*)?$/.exec(l);
    if (own && own[1] === a.el && BY_SYMBOL.has(a.el)) return;
  }
  if (/^R\d*$|^R'+$/.test(label) || !setAbbrev(a, label)) {
    a.el = 'R';
    a.alias = label;
    a.charge = 0;
    delete a.isotope;
    delete a.hCount;
    delete a.abbrev;
  }
}

/**
 * Abbreviation key for a label: a known abbreviation (case-insensitive) or a condensed formula
 * such as "CH2OH". Labels drawn right-to-left ("MeO", "HO2C") are also recognised; with
 * `preferReversed` (label known to be right-justified) the reversed reading wins.
 */
export function resolveAbbrev(label: string, preferReversed = false): string | null {
  const l = label.trim();
  if (!l) return null;
  const rev = reverseLabel(l);
  const order = preferReversed ? [rev, l] : [l, rev];
  for (const c of order) {
    const k = findAbbreviation(c);
    if (k) return k;
  }
  for (const c of order) if (parseCondensedLabel(c)) return c;
  return null;
}

/** Makes `a` an abbreviation atom if `label` is a known/condensed abbreviation. */
export function setAbbrev(a: Atom, label: string, preferReversed = false): boolean {
  const key = resolveAbbrev(label, preferReversed);
  if (!key) return false;
  const g = abbreviationMol(key);
  a.abbrev = key;
  if (g && g.atoms.length) a.el = g.atoms[0].el;
  a.charge = 0;
  delete a.isotope;
  delete a.hCount;
  delete a.radical;
  delete a.alias;
  return true;
}

export interface SupGroup {
  /** Atom indices (0-based) belonging to the superatom. */
  atoms: number[];
  label: string;
}

/**
 * Contracts superatom S-groups (molfile SUP) whose label is a known abbreviation into a single
 * abbreviation atom, provided the group has exactly one attachment atom and its heavy-atom
 * composition and net charge match the abbreviation's structure. Other groups stay expanded.
 */
export function contractSuperatoms(mol: Mol, groups: SupGroup[]): void {
  const remove = new Set<number>();
  for (const g of groups) {
    const set = new Set(g.atoms.filter((i) => i >= 0 && i < mol.atoms.length));
    if (!set.size || [...set].some((i) => remove.has(i))) continue;
    const key = resolveAbbrev(g.label);
    if (!key) continue;
    const ref = abbreviationMol(key);
    if (!ref || !ref.atoms.length) continue;
    // attachment atoms = group atoms bonded to the outside
    const attach = new Set<number>();
    for (const b of mol.bonds) {
      const ia = set.has(b.a), ib = set.has(b.b);
      if (ia !== ib) attach.add(ia ? b.a : b.b);
    }
    if (attach.size !== 1) continue;
    const att = [...attach][0];
    if (mol.atoms[att].el !== ref.atoms[0].el) continue;
    const heavy = (els: string[]) => els.filter((e) => e !== 'H').sort().join(',');
    const groupEls = [...set].map((i) => mol.atoms[i].el);
    if (heavy(groupEls) !== heavy(ref.atoms.map((a) => a.el))) continue;
    const charge = [...set].reduce((s, i) => s + mol.atoms[i].charge, 0);
    if (charge !== ref.atoms.reduce((s, a) => s + a.charge, 0)) continue;
    setAbbrev(mol.atoms[att], key);
    for (const i of set) if (i !== att) remove.add(i);
  }
  if (remove.size) mol.removeAtoms(remove);
}

/**
 * Size limits for automatic stereo perception on import. Ring perception and symmetry classes
 * grow steeply with size (≈0.3 s at these limits), so very large or ring-rich structures are
 * imported without stereo specs (wedges/coordinates are kept, so callers can perceive later).
 */
export const STEREO_MAX_ATOMS = 1000;
export const STEREO_MAX_RINGS = 60;

/** Fills mol.tetra / mol.dbStereo from wedges (2D) or geometry (3D). Best effort; never throws. */
export function perceiveStereo(mol: Mol): void {
  if (mol.atoms.length < 2 || !hasCoordinates(mol)) return;
  if (mol.atoms.length > STEREO_MAX_ATOMS) return;
  const threeD = is3D(mol);
  // nothing that could carry 2D stereo: skip the (costly) analysis
  if (!threeD && !mol.bonds.some((b) => b.style === 'wedge' || b.style === 'hash' || (b.order === 2 && b.style !== 'crossed'))) return;
  if (mol.bonds.length - mol.atoms.length + mol.components().length > STEREO_MAX_RINGS) return;
  try {
    if (threeD) perceiveStereo3D(mol);
    else perceiveStereo2D(mol);
  } catch {
    mol.tetra = [];
    mol.dbStereo = [];
  }
}

type V3 = [number, number, number];

/** Stereo specs from 3D coordinates (model frame y down → converted to right-handed y up). */
export function perceiveStereo3D(mol: Mol): void {
  const pos = (i: number): V3 => {
    const a = mol.atoms[i];
    return [a.x, -a.y, a.z ?? 0];
  };
  const sym = symmetryClasses(mol);
  const rings = perceiveRings(mol);
  const tetra: TetraSpec[] = [];
  for (let i = 0; i < mol.atoms.length; i++) {
    if (!isPotentialStereocenter(mol, i)) continue;
    const nbrs = mol.neighbors(i);
    if (nbrs.length < 3 || nbrs.length > 4) continue;
    const h = implicitH(mol, i);
    const classes = nbrs.map((j) => sym[j]);
    if (new Set(classes).size !== classes.length) continue;
    if (h > 0 && nbrs.length === 4) continue;
    const c = pos(i);
    const vec: V3[] = nbrs.map((j) => {
      const p = pos(j);
      const v: V3 = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
      const l = Math.hypot(v[0], v[1], v[2]) || 1;
      return [v[0] / l, v[1] / l, v[2] / l];
    });
    const order = [...nbrs];
    if (vec.length === 3) {
      vec.push([-(vec[0][0] + vec[1][0] + vec[2][0]), -(vec[0][1] + vec[1][1] + vec[2][1]), -(vec[0][2] + vec[1][2] + vec[2][2])]);
      order.push(-1);
    }
    // reject (near-)planar centres
    const [p0, p1, p2, p3] = vec;
    const d1: V3 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    const d2: V3 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
    const d3: V3 = [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]];
    const vol = d1[0] * (d2[1] * d3[2] - d2[2] * d3[1]) - d1[1] * (d2[0] * d3[2] - d2[2] * d3[0]) + d1[2] * (d2[0] * d3[1] - d2[1] * d3[0]);
    if (Math.abs(vol) < 0.05) continue;
    tetra.push({ center: i, nbrs: order as TetraSpec['nbrs'], ccw: ccwFromPositions(p0, p1, p2, p3) });
  }
  const db: DbSpec[] = [];
  mol.bonds.forEach((b, bi) => {
    if (b.order !== 2 || b.style === 'crossed') return;
    const rs = smallestRingSizeOfBond(rings, bi);
    if (rs > 0 && rs < 8) return;
    const refs: number[] = [];
    for (const e of [b.a, b.b]) {
      const other = e === b.a ? b.b : b.a;
      const nb = mol.neighbors(e).filter((w) => w !== other);
      if (nb.length === 0 || nb.length + implicitH(mol, e) > 2) return;
      if (nb.length === 2 && sym[nb[0]] === sym[nb[1]]) return;
      refs.push(nb[0]);
    }
    const A = pos(b.a), B = pos(b.b), ra = pos(refs[0]), rb = pos(refs[1]);
    const n: V3 = [B[0] - A[0], B[1] - A[1], B[2] - A[2]];
    const nl = Math.hypot(n[0], n[1], n[2]) || 1;
    const u: V3 = [n[0] / nl, n[1] / nl, n[2] / nl];
    const perp = (v: V3): V3 => {
      const d = v[0] * u[0] + v[1] * u[1] + v[2] * u[2];
      return [v[0] - d * u[0], v[1] - d * u[1], v[2] - d * u[2]];
    };
    const va = perp([ra[0] - A[0], ra[1] - A[1], ra[2] - A[2]]);
    const vb = perp([rb[0] - B[0], rb[1] - B[1], rb[2] - B[2]]);
    const dot = va[0] * vb[0] + va[1] * vb[1] + va[2] * vb[2];
    const norm = Math.hypot(...va) * Math.hypot(...vb);
    if (norm < 1e-6 || Math.abs(dot) / norm < 0.2) return; // near-perpendicular: undefined
    db.push({ bond: bi, a: refs[0], b: refs[1], cis: dot > 0 });
  });
  mol.tetra = tetra;
  mol.dbStereo = db;
}

/** Number formatting helpers for fixed-width writers. */
export const pad = (v: string | number, w: number) => String(v).padStart(w);
export function fixed(v: number, digits: number): string {
  const s = (Number.isFinite(v) ? v : 0).toFixed(digits);
  return /^-0\.?0*$/.test(s) ? s.slice(1) : s;
}
