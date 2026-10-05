// Layout of a mechanism step's product. Molecules that the step joins are brought together without being
// distorted (the smaller one is moved rigidly next to its new partner, then only the joint is tidied), rings
// the step closes are cleaned, and the separate species of the product are set out in a row.
import { Mol } from '../chem/mol';
import { clean2D } from '../chem/clean2d';
import { implicitH } from '../chem/valence';

const TAU = 2 * Math.PI;
const wrap = (a: number) => ((a % TAU) + TAU) % TAU;
const angDist = (a: number, b: number) => Math.min(wrap(a - b), wrap(b - a));

/** Direction for a new bond at an atom whose other bonds point along `occupied`, as close to `preferred` as fits. */
function freeDirection(occupied: number[], preferred: number): number {
  if (!occupied.length) return preferred;
  if (occupied.length === 1) {
    // zig-zag: 120° from the existing bond, on the side the new partner was drawn
    const a = occupied[0] + TAU / 3, b = occupied[0] - TAU / 3;
    return angDist(a, preferred) <= angDist(b, preferred) ? a : b;
  }
  const s = occupied.map(wrap).sort((p, q) => p - q);
  let best = preferred, bestGap = -1;
  s.forEach((a0, k) => {
    const a1 = k + 1 < s.length ? s[k + 1] : s[0] + TAU;
    const gap = a1 - a0, mid = a0 + gap / 2;
    const score = gap - 0.01 * angDist(mid, preferred);
    if (score > bestGap) (bestGap = score), (best = mid);
  });
  return best;
}

/** Moves piece `q` (which contains v) rigidly so that v sits one bond length from u, its bulk pointing away from u. */
function attach(mol: Mol, q: number[], u: number, v: number, isPlaced: (atom: number) => boolean): void {
  const U = mol.atoms[u], V = mol.atoms[v];
  const occupied = mol
    .neighbors(u)
    .filter((n) => n !== v && isPlaced(n))
    .map((n) => Math.atan2(mol.atoms[n].y - U.y, mol.atoms[n].x - U.x));
  const dir = freeDirection(occupied, Math.atan2(V.y - U.y, V.x - U.x));
  if (q.length > 1) {
    let cx = 0, cy = 0, k = 0;
    for (const x of q) if (x !== v) (cx += mol.atoms[x].x), (cy += mol.atoms[x].y), k++;
    const rot = dir - Math.atan2(cy / k - V.y, cx / k - V.x);
    const c = Math.cos(rot), s = Math.sin(rot);
    const vx = V.x, vy = V.y;
    for (const x of q) {
      const a = mol.atoms[x];
      const dx = a.x - vx, dy = a.y - vy;
      a.x = vx + dx * c - dy * s;
      a.y = vy + dx * s + dy * c;
    }
  }
  const tx = U.x + Math.cos(dir) - V.x, ty = U.y + Math.sin(dir) - V.y;
  for (const x of q) (mol.atoms[x].x += tx), (mol.atoms[x].y += ty);
}

/** Runs clean2D on one fragment (only `movable` atoms move when given) and copies the result back. */
function cleanFragment(mol: Mol, frag: number[], movable: Set<number> | null): void {
  const sub = new Mol();
  const at = new Map<number, number>();
  for (const i of frag) {
    at.set(i, sub.atoms.length);
    sub.atoms.push({ ...mol.atoms[i] });
  }
  const bondOf: number[] = [];
  mol.bonds.forEach((b, bi) => {
    if (!at.has(b.a) || !at.has(b.b)) return;
    sub.bonds.push({ ...b, a: at.get(b.a)!, b: at.get(b.b)! });
    bondOf.push(bi);
  });
  sub.invalidate();
  try {
    clean2D(sub, movable ? { atoms: new Set([...movable].filter((i) => at.has(i)).map((i) => at.get(i)!)) } : {});
  } catch {
    return; // keep the geometry as it is
  }
  for (const i of frag) {
    const a = sub.atoms[at.get(i)!];
    mol.atoms[i].x = a.x;
    mol.atoms[i].y = a.y;
  }
  // clean2D may re-derive wedges to keep the configuration
  sub.bonds.forEach((b, k) => {
    const o = mol.bonds[bondOf[k]];
    o.style = b.style;
    o.a = frag[b.a];
    o.b = frag[b.b];
  });
  mol.invalidate();
}

/**
 * Tidies a product in place. Bonds the step formed are recognised by id -1. Within each fragment, the pieces
 * those bonds joined are brought together around the largest one; a fragment where the step closed a ring is
 * cleaned as a whole, otherwise only the moved pieces and their joints are.
 */
export function tidyProduct(mol: Mol): void {
  for (const frag of mol.components()) {
    const inFrag = new Set(frag);
    const isNew = (bi: number) => mol.bonds[bi].id === -1;
    const formed = mol.bonds.map((_, bi) => bi).filter((bi) => isNew(bi) && inFrag.has(mol.bonds[bi].a));
    if (!formed.length) continue; // only bonds broken or shifted: keep the drawing
    // the pieces the fragment falls into without the new bonds: the molecules that came together
    const piece = new Map<number, number>();
    const pieces: number[][] = [];
    for (const s of frag) {
      if (piece.has(s)) continue;
      const p: number[] = [];
      const stack = [s];
      piece.set(s, pieces.length);
      while (stack.length) {
        const x = stack.pop()!;
        p.push(x);
        for (const bi of mol.adj[x]) {
          if (isNew(bi)) continue;
          const y = mol.other(bi, x);
          if (!piece.has(y)) (piece.set(y, pieces.length), stack.push(y));
        }
      }
      pieces.push(p);
    }
    let ring = false;
    const joins = new Map<string, number>();
    for (const bi of formed) {
      const pa = piece.get(mol.bonds[bi].a)!, pb = piece.get(mol.bonds[bi].b)!;
      if (pa === pb) ring = true;
      else {
        const k = pa < pb ? `${pa},${pb}` : `${pb},${pa}`;
        joins.set(k, (joins.get(k) ?? 0) + 1);
      }
    }
    if ([...joins.values()].some((c) => c > 1)) ring = true; // two new bonds between the same pieces close a ring
    // bring the pieces together along the new bonds, starting from the largest
    const largest = pieces.reduce((best, p, i) => (p.length > pieces[best].length ? i : best), 0);
    const placed = new Set([largest]);
    const moved = new Set<number>();
    for (let progress = true; progress; ) {
      progress = false;
      for (const bi of formed) {
        const b = mol.bonds[bi];
        const pa = piece.get(b.a)!, pb = piece.get(b.b)!;
        if (placed.has(pa) === placed.has(pb)) continue;
        const [u, v] = placed.has(pa) ? [b.a, b.b] : [b.b, b.a];
        const q = piece.get(v)!;
        attach(mol, pieces[q], u, v, (x) => placed.has(piece.get(x)!));
        placed.add(q);
        for (const x of pieces[q]) moved.add(x);
        moved.add(u);
        progress = true;
      }
    }
    cleanFragment(mol, frag, ring ? null : moved);
  }
}

/**
 * Rough horizontal extent of an atom's label (left, right) in bond lengths: nothing for a plain carbon in a
 * chain; symbol, hydrogens (which may sit on either side) and charge otherwise; the text of a label atom.
 */
function labelExtent(mol: Mol, i: number): [number, number] {
  const a = mol.atoms[i];
  const text = a.abbrev ?? a.alias;
  if (text) {
    const w = 0.34 * text.length + (a.charge ? 0.25 : 0);
    return [w, w];
  }
  const bare = a.el === 'C' && !a.charge && !a.radical && !a.isotope && mol.neighbors(i).length > 0;
  if (bare) return [0.12, 0.12];
  const h = implicitH(mol, i);
  const sym = 0.2 * a.el.length + 0.1;
  const hw = h ? 0.36 + (h > 1 ? 0.2 : 0) : 0;
  return [sym + hw, sym + hw + (a.charge ? 0.28 : 0)];
}

export interface ProductRow {
  /** Centres of the "+" signs between the species. */
  plus: { x: number; y: number }[];
  /** Right edge of the row. */
  maxX: number;
}

/**
 * Sets the fragments of `mol` out left to right from x0, centred on y: the largest species first, the others
 * (leaving groups, by-products) after it in the order they were drawn, with room for a "+" between them.
 */
export function arrangeRow(mol: Mol, x0: number, y: number): ProductRow {
  const frags = mol.components();
  const heavy = (f: number[]) => f.filter((i) => mol.atoms[i].el !== 'H').length;
  const cx = (f: number[]) => f.reduce((s, i) => s + mol.atoms[i].x, 0) / f.length;
  const main = frags.reduce((best, f) => (heavy(f) > heavy(best) || (heavy(f) === heavy(best) && cx(f) < cx(best)) ? f : best), frags[0]);
  const ordered = [main, ...frags.filter((f) => f !== main).sort((a, b) => cx(a) - cx(b))];
  const plus: { x: number; y: number }[] = [];
  let x = x0;
  ordered.forEach((f, k) => {
    // extent including atom labels, so a "+" never lands on an "H₂O" or "Br⁻"
    const ext = f.map((i) => labelExtent(mol, i));
    const minX = Math.min(...f.map((i, q) => mol.atoms[i].x - ext[q][0]));
    const maxX = Math.max(...f.map((i, q) => mol.atoms[i].x + ext[q][1]));
    const ys = f.map((i) => mol.atoms[i].y);
    const midY = (Math.min(...ys) + Math.max(...ys)) / 2;
    if (k > 0) {
      plus.push({ x: x + 0.6, y });
      x += 1.2;
    }
    const dx = x - minX, dy = y - midY;
    for (const i of f) (mol.atoms[i].x += dx), (mol.atoms[i].y += dy);
    x = maxX + dx;
  });
  return { plus, maxX: x };
}
