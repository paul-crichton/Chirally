// XYZ coordinate files: atom count, comment line, then "symbol x y z" (Å) per atom.
// Bonds are perceived from covalent radii; bond orders for C/N/O/S are assigned heuristically.
// Like molfiles, y is flipped into Chirally's y-down frame (z kept); props.dim = '3D'.
import { Mol } from '../mol';
import { element } from '../elements';
import { FormatError, guard, splitLines, normalizeElement, elementByNumber, perceiveStereo, fixed } from './common';

/** Distance tolerance (Å) added to the sum of covalent radii when perceiving bonds. */
export const BOND_TOLERANCE = 0.45;

/** Reads the first frame of an XYZ file. */
export function readXYZ(text: string): Mol {
  return guard('XYZ', () => {
    const lines = splitLines(text);
    let k = 0;
    while (k < lines.length && !lines[k].trim()) k++;
    const n = /^\s*(\d+)\s*$/.exec(lines[k] ?? '');
    if (!n) throw new FormatError('XYZ', 'first line must contain the number of atoms');
    const count = parseInt(n[1], 10);
    const avail = Math.max(0, lines.length - (k + 2));
    if (count > avail) throw new FormatError('XYZ', `file announces ${count} atoms but has only ${avail} coordinate lines`);
    const mol = new Mol();
    mol.name = (lines[k + 1] ?? '').trim();
    for (let i = 0; i < count; i++) {
      const line = lines[k + 2 + i] ?? '';
      const t = line.trim().split(/\s+/);
      if (t.length < 4) throw new FormatError('XYZ', `cannot parse atom line ${k + 3 + i}: "${line}"`);
      const x = Number(t[1]), y = Number(t[2]), z = Number(t[3]);
      if (![x, y, z].every(Number.isFinite)) throw new FormatError('XYZ', `invalid coordinates on line ${k + 3 + i}`);
      let el: string | null;
      if (/^\d+$/.test(t[0])) el = elementByNumber(parseInt(t[0], 10));
      else el = normalizeElement(t[0].replace(/[\d_.:-].*$/, '')) ?? normalizeElement(t[0].slice(0, 1));
      if (t[0] === 'D' || t[0] === 'T') mol.addAtom({ el: 'H', isotope: t[0] === 'D' ? 2 : 3, x, y: -y + 0, z });
      else mol.addAtom({ el: el ?? '*', alias: el ? undefined : t[0], x, y: -y + 0, z });
    }
    for (const a of mol.atoms) if (a.alias === undefined) delete a.alias;
    mol.props.dim = '3D';
    perceiveBonds(mol);
    assignBondOrders(mol);
    perceiveStereo(mol);
    return mol;
  });
}

/** Writes a molecule as XYZ (Å; 2D molecules are written with 1.5 Å bonds and z = 0). */
export function writeXYZ(mol: Mol, comment = mol.name): string {
  const threeD = mol.props.dim === '3D' || mol.atoms.some((a) => a.z);
  const s = threeD ? 1 : 1.5;
  const out = [String(mol.atoms.length), comment.replace(/[\r\n]+/g, ' ')];
  for (const a of mol.atoms) {
    const sym = element(a.el) ? a.el : 'X';
    out.push(`${sym.padEnd(2)} ${fixed(a.x * s, 5).padStart(11)} ${fixed(-a.y * s, 5).padStart(11)} ${fixed((a.z ?? 0) * s, 5).padStart(11)}`);
  }
  return out.join('\n') + '\n';
}

/**
 * Adds single bonds between atoms closer than the sum of their covalent radii + tolerance.
 * Uses a uniform grid so large structures stay fast. Hydrogens keep only their nearest partner.
 */
export function perceiveBonds(mol: Mol, tolerance = BOND_TOLERANCE): void {
  const n = mol.atoms.length;
  const rad = mol.atoms.map((a) => element(a.el)?.covRadius || 0.75);
  const maxR = Math.max(0.5, ...rad);
  const cell = 2 * maxR + tolerance;
  const grid = new Map<string, number[]>();
  const key = (i: number, j: number, k: number) => `${i},${j},${k}`;
  const cellOf = (i: number) => {
    const a = mol.atoms[i];
    return [Math.floor(a.x / cell), Math.floor(a.y / cell), Math.floor((a.z ?? 0) / cell)];
  };
  for (let i = 0; i < n; i++) {
    const c = cellOf(i);
    const kk = key(c[0], c[1], c[2]);
    let list = grid.get(kk);
    if (!list) grid.set(kk, (list = []));
    list.push(i);
  }
  const cand: { i: number; j: number; d: number }[] = [];
  for (let i = 0; i < n; i++) {
    const [cx, cy, cz] = cellOf(i);
    const a = mol.atoms[i];
    for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
      const list = grid.get(key(cx + dx, cy + dy, cz + dz));
      if (!list) continue;
      for (const j of list) {
        if (j <= i) continue;
        const b = mol.atoms[j];
        const d = Math.hypot(a.x - b.x, a.y - b.y, (a.z ?? 0) - (b.z ?? 0));
        if (d > 0.4 && d < rad[i] + rad[j] + tolerance) cand.push({ i, j, d });
      }
    }
  }
  // shortest contacts first; hydrogens bond once
  cand.sort((p, q) => p.d - q.d);
  const hBonded = new Array(n).fill(false);
  for (const { i, j } of cand) {
    const hi = mol.atoms[i].el === 'H', hj = mol.atoms[j].el === 'H';
    if ((hi && hBonded[i]) || (hj && hBonded[j])) continue;
    mol.addBond(i, j, 1);
    if (hi) hBonded[i] = true;
    if (hj) hBonded[j] = true;
  }
}

const TARGET_VALENCE: Record<string, number> = { C: 4, N: 3, O: 2, S: 2, P: 3, B: 3, Si: 4 };

/**
 * Heuristic bond orders for a structure with explicit hydrogens: atoms whose neighbour count is
 * below their usual valence receive multiple bonds. Forced assignments (atoms with a single
 * unsaturated neighbour) are made first; remaining conjugated systems are resolved by matching.
 */
export function assignBondOrders(mol: Mol): void {
  const n = mol.atoms.length;
  const free = mol.atoms.map((a, i) => {
    const v = TARGET_VALENCE[a.el];
    if (v === undefined || a.charge) return 0;
    if (a.el === 'S' && mol.degree(i) > 2) return 0;
    return Math.max(0, v - mol.degree(i));
  });
  const unsatNbrs = (i: number) => mol.adj[i].filter((bi) => free[mol.other(bi, i)] > 0);
  const raise = (bi: number, by: number) => {
    const b = mol.bonds[bi];
    b.order += by;
    free[b.a] -= by;
    free[b.b] -= by;
  };
  const propagate = () => {
    let changed = true;
    while (changed) {
      changed = false;
      for (let i = 0; i < n; i++) {
        if (free[i] <= 0) continue;
        const nb = unsatNbrs(i);
        if (nb.length !== 1) continue;
        const bi = nb[0];
        const j = mol.other(bi, i);
        const by = Math.min(free[i], free[j], 3 - mol.bonds[bi].order);
        if (by > 0) {
          raise(bi, by);
          changed = true;
        }
      }
    }
  };
  propagate();
  // remaining unsaturated atoms (conjugated rings/chains): perfect matching by bounded
  // backtracking, falling back to a greedy most-constrained-first matching
  const atoms: number[] = [];
  for (let i = 0; i < n; i++) if (free[i] > 0 && unsatNbrs(i).length) atoms.push(i);
  if (atoms.length) {
    const partners = (i: number, matched: number[]) =>
      unsatNbrs(i).map((bi) => mol.other(bi, i)).filter((j) => matched[j] === -1);
    const perfect = (): number[] | null => {
      const matched = new Array(n).fill(-1);
      let steps = 0;
      const rec = (): boolean => {
        if (++steps > 50000) return false;
        let best = -1;
        let opts: number[] = [];
        for (const i of atoms) {
          if (matched[i] !== -1) continue;
          const o = partners(i, matched);
          if (!o.length) return false;
          if (best < 0 || o.length < opts.length) {
            best = i;
            opts = o;
          }
        }
        if (best < 0) return true;
        for (const j of opts) {
          matched[best] = j;
          matched[j] = best;
          if (rec()) return true;
          matched[best] = -1;
          matched[j] = -1;
        }
        return false;
      };
      return atoms.length % 2 === 0 && rec() ? matched : null;
    };
    const greedy = (): number[] => {
      const matched = new Array(n).fill(-1);
      for (;;) {
        let best = -1;
        let opts: number[] = [];
        for (const i of atoms) {
          if (matched[i] !== -1) continue;
          const o = partners(i, matched);
          if (!o.length) continue;
          if (best < 0 || o.length < opts.length) {
            best = i;
            opts = o;
          }
        }
        if (best < 0) return matched;
        matched[best] = opts[0];
        matched[opts[0]] = best;
      }
    };
    const matched = perfect() ?? greedy();
    for (let i = 0; i < n; i++) {
      const j = matched[i];
      if (j > i) raise(mol.bondBetween(i, j), 1);
    }
    propagate();
  }
  mol.invalidate();
}
