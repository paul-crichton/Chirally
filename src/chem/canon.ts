import { atomicNumber } from './elements';
import { Mol } from './mol';
import { totalH } from './valence';
import { RingInfo, perceiveRings } from './rings';

/**
 * Canonical atom ranks (0..n-1, all distinct) using iterative refinement of atom invariants
 * (Weininger/Morgan style) with deterministic tie-breaking.
 * If `breakTies` is false the symmetry classes are returned instead (equal ranks = equivalent atoms).
 */
export function canonicalRanks(mol: Mol, opts: { breakTies?: boolean; rings?: RingInfo; aromatic?: boolean[] } = {}): number[] {
  const n = mol.atoms.length;
  if (n === 0) return [];
  const rings = opts.rings ?? perceiveRings(mol);
  const inv: string[] = mol.atoms.map((a, i) => {
    const z = a.abbrev ? 200 : a.el === 'R' || a.el === '*' ? 0 : atomicNumber(a.el);
    return [
      String(mol.degree(i)).padStart(2, '0'),
      String(z).padStart(3, '0'),
      String(a.isotope ?? 0).padStart(3, '0'),
      String(a.charge + 50).padStart(3, '0'),
      String(totalH(mol, i)).padStart(2, '0'),
      rings.inRing[i] ? '1' : '0',
      opts.aromatic?.[i] ? '1' : '0',
      a.abbrev ?? a.alias ?? '',
      String(a.radical ?? 0),
    ].join('|');
  });
  let ranks = rankStrings(inv);
  ranks = refine(mol, ranks);
  if (opts.breakTies === false) return ranks;
  // tie breaking
  for (;;) {
    const counts = new Map<number, number>();
    for (const r of ranks) counts.set(r, (counts.get(r) ?? 0) + 1);
    let tied = -1;
    for (const [r, c] of [...counts.entries()].sort((p, q) => p[0] - q[0])) {
      if (c > 1) {
        tied = r;
        break;
      }
    }
    if (tied < 0) break;
    const idx = ranks.indexOf(tied);
    // Double all ranks and lower one member of the tied class by one
    ranks = ranks.map((r) => r * 2 + 1);
    ranks[idx] -= 1;
    ranks = refine(mol, rankNumbers(ranks));
  }
  return ranks;
}

function bondCode(order: number): number {
  return order === 1.5 ? 4 : Math.round(order);
}

function refine(mol: Mol, ranks: number[]): number[] {
  const n = ranks.length;
  // neighbour lists (atom index, bond code) computed once
  const nbr: number[][] = new Array(n);
  const code: number[][] = new Array(n);
  for (let i = 0; i < n; i++) {
    const adj = mol.adj[i];
    nbr[i] = adj.map((bi) => mol.other(bi, i));
    code[i] = adj.map((bi) => bondCode(mol.bonds[bi].order));
  }
  const keys: Float64Array[] = new Array(n);
  const idx = new Array(n);
  let classes = new Set(ranks).size;
  const cmp = (p: number, q: number): number => {
    if (ranks[p] !== ranks[q]) return ranks[p] - ranks[q];
    const a = keys[p], b = keys[q];
    if (a.length !== b.length) return a.length - b.length;
    for (let k = 0; k < a.length; k++) if (a[k] !== b[k]) return a[k] - b[k];
    return 0;
  };
  for (let iter = 0; iter < n + 2; iter++) {
    for (let i = 0; i < n; i++) {
      const ni = nbr[i], ci = code[i];
      const kv = new Float64Array(ni.length);
      for (let k = 0; k < ni.length; k++) kv[k] = ranks[ni[k]] * 8 + ci[k];
      kv.sort();
      keys[i] = kv;
      idx[i] = i;
    }
    idx.sort(cmp);
    const next = new Array(n).fill(0);
    let r = 0;
    for (let k = 0; k < n; k++) {
      if (k > 0 && cmp(idx[k - 1], idx[k]) !== 0) r = k;
      next[idx[k]] = r;
    }
    let nc = 0;
    for (let k = 0; k < n; k++) if (k === 0 || next[idx[k]] !== next[idx[k - 1]]) nc++;
    ranks = next;
    if (nc === classes) break;
    classes = nc;
  }
  return ranks;
}

function rankStrings(v: string[]): number[] {
  const sorted = [...new Set(v)].sort();
  const m = new Map(sorted.map((s, i) => [s, i]));
  // dense → convert to "competition" ranks so ties stay tied and ranks fit 0..n-1
  return rankNumbers(v.map((s) => m.get(s)!));
}

function rankNumbers(v: number[]): number[] {
  const idx = v.map((_, i) => i).sort((p, q) => v[p] - v[q]);
  const out = new Array(v.length).fill(0);
  let r = 0;
  for (let k = 0; k < idx.length; k++) {
    if (k > 0 && v[idx[k]] !== v[idx[k - 1]]) r = k;
    out[idx[k]] = r;
  }
  return out;
}

/** Symmetry classes: atoms with equal values are topologically equivalent (graph automorphism approximation). */
export function symmetryClasses(mol: Mol): number[] {
  return canonicalRanks(mol, { breakTies: false });
}
