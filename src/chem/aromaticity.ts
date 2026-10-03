import { Mol } from './mol';
import { perceiveRings, RingInfo } from './rings';
import { implicitH, chargedValences } from './valence';

export interface AromaticityResult {
  atoms: boolean[];
  bonds: boolean[];
  rings: RingInfo;
}

const LONE_PAIR_DONORS = new Set(['N', 'O', 'S', 'Se', 'P', 'Te', 'As']);
const ELECTRONEGATIVE = new Set(['O', 'N', 'S', 'Se']);

/**
 * π-electron contribution of atom i to a ring belonging to ring system `sysAtoms`.
 * Returns -1 if the atom cannot be part of an aromatic ring.
 */
function piElectrons(mol: Mol, i: number, sysAtoms: Set<number>): number {
  const a = mol.atoms[i];
  if (a.abbrev) return -1;
  let dblInSys = 0;
  let dblExoEN = 0;
  let dblExoOther = 0;
  let triple = 0;
  let aromaticBonds = 0;
  for (const bi of mol.adj[i]) {
    const b = mol.bonds[bi];
    const j = mol.other(bi, i);
    if (b.order === 2) {
      if (sysAtoms.has(j)) dblInSys++;
      else if (ELECTRONEGATIVE.has(mol.atoms[j].el)) dblExoEN++;
      else dblExoOther++;
    } else if (b.order === 3) triple++;
    else if (b.order === 1.5) aromaticBonds++;
  }
  if (triple) return -1;
  if (dblExoOther) return -1;
  if (dblInSys > 1) return -1; // cumulated
  if (aromaticBonds >= 2 && dblInSys === 0) {
    // already flagged aromatic (drawn with delocalised bonds) – treat like a double-bonded atom
    if (a.el === 'C' || a.el === 'N' || a.el === 'B') return 1;
    return 2;
  }
  if (dblInSys === 1) return 1;
  if (dblExoEN === 1) return 0; // e.g. C=O in pyridone/quinone
  // no double bond: lone pair donor, carbanion, carbocation, boron
  const deg = mol.degree(i) + implicitH(mol, i);
  if (a.el === 'C') {
    if (a.charge === -1) return 2;
    if (a.charge === 1) return 0;
    return -1;
  }
  if (a.el === 'B' && a.charge === 0) return 0;
  if (LONE_PAIR_DONORS.has(a.el)) {
    if (a.charge > 0) return -1;
    if (a.el === 'N' || a.el === 'P' || a.el === 'As') return deg <= 3 ? 2 : -1;
    return deg <= 2 || a.charge < 0 ? 2 : -1;
  }
  return -1;
}

function isHuckel(n: number): boolean {
  return n >= 2 && (n - 2) % 4 === 0;
}

/** Perceives aromatic atoms and bonds from a Kekulé (or partially delocalised) structure. */
export function perceiveAromaticity(mol: Mol, rings?: RingInfo): AromaticityResult {
  const info = rings ?? perceiveRings(mol);
  const atoms = new Array(mol.atoms.length).fill(false);
  const bonds = new Array(mol.bonds.length).fill(false);

  for (const system of info.systems) {
    const sysAtoms = new Set<number>();
    for (const r of system) for (const a of info.rings[r]) sysAtoms.add(a);
    const pi = new Map<number, number>();
    for (const a of sysAtoms) pi.set(a, piElectrons(mol, a, sysAtoms));
    const ringOK = (atomList: number[]) => atomList.every((a) => pi.get(a)! >= 0);
    const ringSum = (atomList: number[]) => atomList.reduce((s, a) => s + pi.get(a)!, 0);
    const aromaticRings = new Set<number>();
    for (const r of system) {
      const ring = info.rings[r];
      if (ring.length > 24) continue;
      if (ringOK(ring) && isHuckel(ringSum(ring))) aromaticRings.add(r);
    }
    // Fused pairs (envelopes) for systems such as azulene
    if (system.length > 1) {
      for (let p = 0; p < system.length; p++) {
        for (let q = p + 1; q < system.length; q++) {
          const r1 = system[p], r2 = system[q];
          if (aromaticRings.has(r1) && aromaticRings.has(r2)) continue;
          const shared = info.ringBonds[r1].filter((b) => info.ringBonds[r2].includes(b));
          if (shared.length !== 1) continue;
          const env = new Set([...info.rings[r1], ...info.rings[r2]]);
          const list = [...env];
          if (ringOK(list) && isHuckel(ringSum(list))) {
            aromaticRings.add(r1);
            aromaticRings.add(r2);
          }
        }
      }
    }
    for (const r of aromaticRings) {
      for (const a of info.rings[r]) atoms[a] = true;
      for (const b of info.ringBonds[r]) bonds[b] = true;
    }
  }
  return { atoms, bonds, rings: info };
}

/**
 * Assigns alternating single/double bonds to bonds flagged aromatic (e.g. after SMILES parsing).
 * `aromaticAtoms[i]` marks atoms written in lowercase; their implicit-H counts must already be fixed
 * (hCount) or derivable. Returns false if no valid Kekulé structure exists.
 */
export function kekulize(mol: Mol, aromaticAtoms: boolean[], aromaticBonds: boolean[]): boolean {
  const n = mol.atoms.length;
  // Determine which aromatic atoms need a double bond
  const need = new Array(n).fill(false);
  for (let i = 0; i < n; i++) {
    if (!aromaticAtoms[i]) continue;
    const a = mol.atoms[i];
    let sigma = 0;
    let exoDouble = 0;
    let aroCount = 0;
    for (const bi of mol.adj[i]) {
      const b = mol.bonds[bi];
      if (aromaticBonds[bi]) {
        aroCount++;
        sigma += 1;
      } else {
        sigma += b.order;
        if (b.order === 2) exoDouble++;
      }
    }
    const h = a.hCount ?? 0;
    const vals = chargedValences(a.el, a.charge);
    const used = sigma + h + (a.radical ?? 0);
    let v = vals.find((x) => x >= used) ?? vals[vals.length - 1] ?? 0;
    // A pi bond is required if one more bond fits the smallest adequate valence
    // (pyridine N: 3 - 2 = 1 → needs; pyrrole [nH]: 3 - 3 = 0 → no)
    if (aroCount >= 2 && v - used >= 1) {
      // For atoms like S that allow valence 4, prefer the lowest valence (no double bond) when used == 2
      if ((a.el === 'S' || a.el === 'Se' || a.el === 'O' || a.el === 'Te') && used >= 2) need[i] = false;
      else need[i] = true;
    }
    if (exoDouble) need[i] = false;
  }

  const candBonds: number[] = [];
  aromaticBonds.forEach((f, bi) => {
    if (f) {
      mol.bonds[bi].order = 1;
      const b = mol.bonds[bi];
      if (need[b.a] && need[b.b]) candBonds.push(bi);
    }
  });

  const matched = new Array(n).fill(-1);
  const order: number[] = [];
  for (let i = 0; i < n; i++) if (need[i]) order.push(i);
  if (order.length % 2 === 1) {
    // try to drop an atom that can live without a double bond (e.g. pyrrole-type N written as 'n')
    return fallbackKekulize(mol, aromaticAtoms, aromaticBonds, need, candBonds);
  }
  const opts = (i: number) => mol.adj[i].filter((bi) => aromaticBonds[bi] && need[mol.other(bi, i)]);
  const ok = matchAll(mol, order, matched, opts);
  if (!ok) return fallbackKekulize(mol, aromaticAtoms, aromaticBonds, need, candBonds);
  applyMatching(mol, matched);
  return true;
}

function applyMatching(mol: Mol, matched: number[]): void {
  for (let i = 0; i < matched.length; i++) {
    const j = matched[i];
    if (j > i) {
      const bi = mol.bondBetween(i, j);
      if (bi >= 0) mol.bonds[bi].order = 2;
    }
  }
}

function matchAll(mol: Mol, atoms: number[], matched: number[], opts: (i: number) => number[]): boolean {
  // Order atoms by number of options (most constrained first) and backtrack.
  let steps = 0;
  const LIMIT = 200000;
  const rec = (): boolean => {
    if (++steps > LIMIT) return false;
    let best = -1;
    let bestOpts: number[] = [];
    for (const i of atoms) {
      if (matched[i] >= 0) continue;
      const o = opts(i).filter((bi) => matched[mol.other(bi, i)] < 0);
      if (o.length === 0) return false;
      if (best < 0 || o.length < bestOpts.length) {
        best = i;
        bestOpts = o;
        if (o.length === 1) break;
      }
    }
    if (best < 0) return true;
    for (const bi of bestOpts) {
      const j = mol.other(bi, best);
      matched[best] = j;
      matched[j] = best;
      if (rec()) return true;
      matched[best] = -1;
      matched[j] = -1;
    }
    return false;
  };
  return rec();
}

function fallbackKekulize(
  mol: Mol,
  aromaticAtoms: boolean[],
  aromaticBonds: boolean[],
  need: boolean[],
  _cand: number[],
): boolean {
  // Candidate atoms that may become pyrrole-like ([nH], [pH]) if a perfect matching fails.
  const flexible = [];
  for (let i = 0; i < mol.atoms.length; i++) {
    const a = mol.atoms[i];
    if (need[i] && (a.el === 'N' || a.el === 'P') && a.charge === 0 && mol.degree(i) + (a.hCount ?? 0) === 2) flexible.push(i);
  }
  const base = need.slice();
  const tryWith = (drop: number[]): boolean => {
    const nd = base.slice();
    for (const d of drop) nd[d] = false;
    const atoms = [];
    for (let i = 0; i < nd.length; i++) if (nd[i]) atoms.push(i);
    if (atoms.length % 2) return false;
    const matched = new Array(mol.atoms.length).fill(-1);
    const opts = (i: number) => mol.adj[i].filter((bi) => aromaticBonds[bi] && nd[mol.other(bi, i)]);
    if (!matchAll(mol, atoms, matched, opts)) return false;
    applyMatching(mol, matched);
    for (const d of drop) mol.atoms[d].hCount = (mol.atoms[d].hCount ?? 0) + 1;
    return true;
  };
  for (const f of flexible) if (tryWith([f])) return true;
  for (let p = 0; p < flexible.length; p++)
    for (let q = p + 1; q < flexible.length; q++) if (tryWith([flexible[p], flexible[q]])) return true;
  // give up: leave as single bonds but mark them delocalised
  aromaticBonds.forEach((f, bi) => {
    if (f) mol.bonds[bi].order = 1.5;
  });
  void aromaticAtoms;
  return false;
}
