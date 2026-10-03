import { allowedValences, element, IMPLICIT_H_ELEMENTS, valenceElectrons } from './elements';
import { Mol } from './mol';

/** Bond-order contribution of a bond to valence (hydrogen/dative/zero-order bonds count 0). */
export function bondValence(order: number, style: string): number {
  if (style === 'hbond' || style === 'dative') return 0;
  return order;
}

/** Sum of explicit bond orders around atom i. */
export function bondOrderSum(mol: Mol, i: number): number {
  let s = 0;
  for (const bi of mol.adj[i]) {
    const b = mol.bonds[bi];
    s += bondValence(b.order, b.style);
  }
  return s;
}

/**
 * Valences allowed for an element carrying the given formal charge, using the isoelectronic rule
 * for main-group elements (N+ behaves like C, O- like F, C- like N, B- like C, …).
 */
export function chargedValences(el: string, charge: number): number[] {
  const e = element(el);
  if (!e) return [];
  if (el === 'H') return charge === 0 ? [1] : [0];
  const g = e.group;
  if (g >= 13 && g <= 18) {
    // Isoelectronic rule: C+ (eff. group 13) → 3, N+ (14) → 4, O- (17) → 1, B- (14) → 4 …
    const eff = g - charge;
    if (eff < 13) return [Math.max(0, eff - 10)];
    if (eff > 18) return [0];
    // Period-3+ elements keep their expanded valences even when charged.
    return allowedValences(eff, e.period);
  }
  if (g === 1 || g === 2) return [Math.max(0, g - Math.abs(charge))];
  return [];
}

/** Number of implicit hydrogens on atom i (0 when an explicit hCount is set — use totalH for the total). */
export function implicitH(mol: Mol, i: number): number {
  const a = mol.atoms[i];
  if (a.hCount !== undefined && a.hCount !== null) return a.hCount;
  if (a.abbrev) return 0;
  if (!IMPLICIT_H_ELEMENTS.has(a.el)) return 0;
  const used = bondOrderSum(mol, i) + (a.radical ?? 0);
  const vals = chargedValences(a.el, a.charge);
  if (!vals.length) return 0;
  // Aromatic bond orders (1.5) can leave fractional sums: round to nearest half-bond downward.
  const u = Math.round(used * 2) / 2;
  for (const v of vals) {
    if (v >= u - 1e-6) return Math.max(0, Math.floor(v - u + 1e-6));
  }
  return 0;
}

/** Number of explicit hydrogen atoms bonded to atom i. */
export function explicitHNeighbors(mol: Mol, i: number): number {
  let n = 0;
  for (const j of mol.neighbors(i)) if (mol.atoms[j].el === 'H') n++;
  return n;
}

/** Total hydrogens attached to atom i (implicit/explicit count + H atoms in the graph). */
export function totalH(mol: Mol, i: number): number {
  return implicitH(mol, i) + explicitHNeighbors(mol, i);
}

/** Total valence currently used (bonds + implicit H + radical electrons). */
export function usedValence(mol: Mol, i: number): number {
  return bondOrderSum(mol, i) + implicitH(mol, i);
}

/** True if the atom exceeds every allowed valence for its element/charge (drawn as a red warning). */
export function hasValenceError(mol: Mol, i: number): boolean {
  const a = mol.atoms[i];
  if (a.abbrev || a.el === 'R' || a.el === '*' || a.alias) return false;
  const e = element(a.el);
  if (!e) return false;
  const vals = chargedValences(a.el, a.charge);
  if (!vals.length) return false; // metals: anything goes
  const used = bondOrderSum(mol, i) + (a.hCount ?? 0) + (a.radical ?? 0);
  const max = Math.max(...vals);
  return used > max + 1e-6;
}

/**
 * Electron bookkeeping for one atom: non-bonding electrons N = V − FC − B,
 * where V = valence electrons, FC = formal charge, B = number of bonds (incl. H).
 */
export function nonBondingElectrons(mol: Mol, i: number): number {
  const a = mol.atoms[i];
  if (a.abbrev || !element(a.el)) return 0;
  const V = valenceElectrons(a.el);
  const B = bondOrderSum(mol, i) + implicitH(mol, i);
  return Math.max(0, Math.round(V - a.charge - B));
}

/** Lone pairs on atom i (non-bonding electrons not used for radicals, halved). */
export function lonePairCount(mol: Mol, i: number): number {
  const n = nonBondingElectrons(mol, i) - (mol.atoms[i].radical ?? 0);
  return Math.max(0, Math.floor(n / 2));
}

/** Electrons around the atom counted for the octet rule (2 per bond + non-bonding). */
export function octetCount(mol: Mol, i: number): number {
  const B = bondOrderSum(mol, i) + implicitH(mol, i);
  return 2 * B + nonBondingElectrons(mol, i);
}

/** Maximum electron count allowed by the octet/duet rule (period 2 strict, period ≥3 may expand). */
export function octetLimit(el: string): number {
  const e = element(el);
  if (!e) return 99;
  if (e.period === 1) return 2;
  if (e.period === 2) return 8;
  return 12;
}
