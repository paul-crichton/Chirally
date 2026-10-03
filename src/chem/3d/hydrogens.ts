// Implicit → explicit hydrogens, keeping stereo specs consistent.
import { Mol, TetraSpec } from '../mol';
import { implicitH } from '../valence';

/**
 * Returns a copy of `mol` in which every implicit hydrogen is an explicit H atom.
 * New H atoms are appended after the existing atoms, so original atom (and bond) indices are unchanged.
 * Tetrahedral specs whose neighbour list contains -1 (implicit H) are re-pointed to the new H atom;
 * a -1 on a centre that has no hydrogen (i.e. a lone pair, e.g. sulfoxides) is kept as -1.
 */
export function addExplicitHydrogens(mol: Mol): Mol {
  const out = mol.clone();
  const n = mol.atoms.length;
  const counts = mol.atoms.map((_, i) => implicitH(mol, i));
  const firstH = new Array<number>(n).fill(-1);
  let maxId = out.atoms.reduce((m, a) => Math.max(m, a.id), 0);
  let maxBid = out.bonds.reduce((m, b) => Math.max(m, b.id), 0);
  for (let i = 0; i < n; i++) {
    const a = out.atoms[i];
    const h = counts[i];
    for (let k = 0; k < h; k++) {
      const hi = out.atoms.length;
      out.atoms.push({ id: ++maxId, el: 'H', x: a.x, y: a.y, z: a.z ?? 0, charge: 0 });
      out.bonds.push({ id: ++maxBid, a: i, b: hi, order: 1, style: 'plain' });
      if (k === 0) firstH[i] = hi;
    }
    // All hydrogens are now explicit: pin the count so valence code does not add more.
    if (h > 0 || a.hCount !== undefined) a.hCount = 0;
  }
  out.invalidate();
  out.tetra = out.tetra.map((t) => ({
    center: t.center,
    nbrs: t.nbrs.map((v) => (v < 0 && firstH[t.center] >= 0 ? firstH[t.center] : v)) as TetraSpec['nbrs'],
    ccw: t.ccw,
  }));
  return out;
}
