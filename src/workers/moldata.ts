// Plain-data transport of Mol objects across worker boundaries.
import { Mol } from '../chem/mol';

export interface MolData {
  atoms: Mol['atoms'];
  bonds: Mol['bonds'];
  tetra: Mol['tetra'];
  dbStereo: Mol['dbStereo'];
}

export function toData(m: Mol): MolData {
  return { atoms: m.atoms, bonds: m.bonds, tetra: m.tetra, dbStereo: m.dbStereo };
}

export function fromData(d: MolData, dim3 = true): Mol {
  const m = new Mol();
  m.atoms = d.atoms.map((a) => ({ ...a }));
  m.bonds = d.bonds.map((b) => ({ ...b }));
  m.tetra = d.tetra.map((t) => ({ ...t, nbrs: [...t.nbrs] as Mol['tetra'][0]['nbrs'] }));
  m.dbStereo = d.dbStereo.map((x) => ({ ...x }));
  if (dim3) m.props = { dim: '3D' };
  return m;
}
