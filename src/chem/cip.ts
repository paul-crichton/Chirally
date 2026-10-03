// CIP stereodescriptors (R/S, E/Z). STUB — implemented by the properties/CIP module owner.
import { Mol } from './mol';

export type CenterLabel = 'R' | 'S' | 'r' | 's';
export type BondLabel = 'E' | 'Z';

export interface CIPResult {
  /** atom index → descriptor */
  centers: Map<number, CenterLabel>;
  /** bond index → descriptor */
  bonds: Map<number, BondLabel>;
}

/**
 * Assigns CIP descriptors from mol.tetra / mol.dbStereo (populate them first via parseSmiles
 * or perceiveStereo2D). Implicit hydrogens are handled; abbreviations must be expanded by the caller.
 */
export function assignCIP(mol: Mol): CIPResult {
  void mol;
  return { centers: new Map(), bonds: new Map() };
}
