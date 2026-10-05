// Arrow-pushing engine: applies curved (electron-flow) arrows to produce the next intermediate,
// with full electron bookkeeping (lone pairs, bonds, formal charges, radicals) and octet checks.
import { ChemDoc, CurvedArrowObj, Anchor, ArrowObj } from './types';
import { docToMol, adjacency, fragmentOf, docBounds } from './document';
import { Mol } from '../chem/mol';
import { clean2D } from '../chem/clean2d';
import { implicitH, nonBondingElectrons, bondOrderSum, octetLimit, bondValence } from '../chem/valence';
import { kekulize } from '../chem/aromaticity';
import { abbreviationMol, expandAbbreviations } from '../chem/abbreviations';
import { valenceElectrons, element } from '../chem/elements';
import { writeSmiles, suppressHydrogens } from '../chem/smiles';

export interface MechanismWarning {
  message: string;
  atomIds?: number[];
  /** 'error': the product is not a valid structure (it is previewed, not inserted); 'warning': inserted, but check it. */
  level: 'error' | 'warning';
}

export interface MechanismResult {
  /** Product structure; atoms keep the document ids of the atoms they came from. */
  product: Mol;
  /** Document atom ids of the reactant fragments used. */
  reactantAtomIds: number[];
  warnings: MechanismWarning[];
  /** true for a valid step that only moves π / lone-pair electrons (σ-connectivity unchanged → ↔ arrow). */
  resonance: boolean;
  /** No error-level warnings: the product is a valid next intermediate. */
  ok: boolean;
  /** The arrows changed at least one bond or electron count. */
  changed: boolean;
  /** Human-readable summary of what happened, e.g. "formed C–O, broke C–Br". */
  summary: string[];
}

function atomName(mol: Mol, i: number): string {
  const a = mol.atoms[i];
  return a.abbrev ?? a.alias ?? a.el;
}

/** Groups curved arrows into independent "steps": arrows sharing fragments belong together. */
export interface ArrowGroup {
  arrows: number[];
  atoms: number[];
  maxX: number;
  /** Every arrow in the group already produced a step whose reaction/resonance arrow is still in the drawing. */
  applied: boolean;
}

/** True when the curved arrow was applied and the arrow drawn for that step still exists. */
export function isApplied(doc: ChemDoc, c: CurvedArrowObj): boolean {
  return c.step !== undefined && doc.arrows.has(c.step);
}

export function arrowGroups(doc: ChemDoc): ArrowGroup[] {
  const adj = adjacency(doc);
  const fragOf = new Map<number, number>();
  const frags: number[][] = [];
  const getFrag = (atomId: number) => {
    let f = fragOf.get(atomId);
    if (f === undefined) {
      const comp = fragmentOf(doc, atomId, adj);
      f = frags.length;
      frags.push(comp);
      for (const a of comp) fragOf.set(a, f);
    }
    return f;
  };
  const anchorAtoms = (an: Anchor): number[] => {
    if (an.type === 'atom') return [an.id];
    if (an.type === 'bond') {
      const b = doc.bonds.get(an.id);
      return b ? [b.a, b.b] : [];
    }
    if (an.type === 'between') return [an.a, an.b];
    return [];
  };
  // union-find over fragments
  const parent: number[] = [];
  const find = (x: number): number => (parent[x] === x ? x : (parent[x] = find(parent[x])));
  const arrowFrags = new Map<number, number[]>();
  for (const c of doc.curved.values()) {
    const fs = [...anchorAtoms(c.from), ...anchorAtoms(c.to)].map(getFrag);
    for (const f of fs) if (parent[f] === undefined) parent[f] = f;
    for (let k = 1; k < fs.length; k++) {
      const a = find(fs[0]), b = find(fs[k]);
      if (a !== b) parent[a] = b;
    }
    arrowFrags.set(c.id, fs);
  }
  const groups = new Map<number, { arrows: number[]; atoms: Set<number> }>();
  for (const [cid, fs] of arrowFrags) {
    if (!fs.length) continue;
    const root = find(fs[0]);
    if (!groups.has(root)) groups.set(root, { arrows: [], atoms: new Set() });
    const g = groups.get(root)!;
    g.arrows.push(cid);
    for (const f of fs) for (const a of frags[f]) g.atoms.add(a);
  }
  return [...groups.values()].map((g) => {
    let maxX = -Infinity;
    for (const a of g.atoms) maxX = Math.max(maxX, doc.atoms.get(a)?.x ?? -Infinity);
    return { arrows: g.arrows, atoms: [...g.atoms], maxX, applied: g.arrows.every((id) => isApplied(doc, doc.curved.get(id)!)) };
  });
}

/** Where an arrow takes its electrons from. */
type Source = { atom: number } | { bond: [number, number] };

/** Where an arrow puts its electrons. */
type Destination =
  | { t: 'lp'; atom: number } // non-bonding electrons on an atom
  | { t: 'bond'; a: number; b: number } // into the a–b bond (formed if absent)
  | { t: 'pick'; base: number; ends: [number, number]; prefer: number; weight: number; nearest: number } // a bond from `base` to one of `ends`
  | { t: 'back' }; // unusable target: the electrons stay where they were

interface ArrowPlan {
  k: number;
  src: Source;
  dest: Destination;
}

interface Outcome {
  product: Mol;
  warnings: MechanismWarning[];
  summary: string[];
  changed: boolean;
  resonance: boolean;
  errors: number;
}

/**
 * Applies the given curved arrows. Each arrow moves 2 electrons (or 1 for fishhooks) from its source
 * (atom lone pair or bond) to its target (atom, bond, or the space between two atoms).
 *
 * Some drawings do not say which atom gets a new bond: an arrow from a bond to an atom outside it, or from a
 * lone pair to a bond the atom is not part of. Every reading is tried; readings that break the octet rule
 * lose, then the chemically expected one wins (1,2-shifts, hydrogen transfer, the polarity of the bond,
 * Markovnikov/Michael selectivity), and only then the atom drawn nearer.
 */
export function applyArrows(doc: ChemDoc, arrowIds: number[]): MechanismResult {
  const fixed: MechanismWarning[] = [];
  const error = (message: string, atomIds?: number[]) => fixed.push({ message, atomIds, level: 'error' });
  const arrows = arrowIds.map((id) => doc.curved.get(id)).filter(Boolean) as CurvedArrowObj[];
  const missing = () => error('An arrow is attached to an atom or bond that no longer exists');
  // involved fragments
  const adj = adjacency(doc);
  const involved = new Set<number>();
  const addFrag = (id: number) => {
    if (involved.has(id) || !doc.atoms.has(id)) return;
    for (const a of fragmentOf(doc, id, adj)) involved.add(a);
  };
  for (const c of arrows) {
    for (const an of [c.from, c.to]) {
      if (an.type === 'atom') addFrag(an.id);
      else if (an.type === 'bond') {
        const b = doc.bonds.get(an.id);
        if (b) addFrag(b.a);
      } else if (an.type === 'between') {
        addFrag(an.a);
        addFrag(an.b);
      }
    }
  }
  const { mol: drawn, index } = docToMol(doc, involved);
  // shorthand labels (OMe, CO2H, Ph, …) take part as the atoms they stand for; the attachment atom keeps the
  // label's index (so arrows still find it) and the label groups the arrows leave alone are labels again afterwards
  const mol = expandAbbreviations(drawn);
  const labelGroups: { at: number; label: string; members: number[]; attachCharge: number }[] = [];
  drawn.atoms.forEach((a, i) => {
    const g = a.abbrev ? abbreviationMol(a.abbrev) : null;
    if (!g) return;
    const members = [i];
    for (let k = drawn.atoms.length; k < mol.atoms.length; k++) if (mol.atoms[k].id === a.id) members.push(k);
    labelGroups.push({ at: i, label: a.abbrev!, members, attachCharge: g.atoms[0].charge });
  });
  // an arrow at the H of an atom's label acts on one of its implicit hydrogens: that H becomes a real atom for the
  // step (-1: the atom has none) and is folded back into the implicit H count of whichever atom holds it afterwards
  const hydrogenOf = new Map<number, number>();
  for (const c of arrows)
    for (const an of [c.from, c.to]) {
      const i = an.type === 'atom' && an.h ? index.get(an.id) : undefined;
      if (i === undefined || hydrogenOf.has(i)) continue;
      const a = mol.atoms[i];
      const hc = implicitH(mol, i);
      if (hc < 1) {
        error(`${atomName(mol, i)} has no hydrogen to move`, [a.id]);
        hydrogenOf.set(i, -1);
        continue;
      }
      let dx = 0, dy = 0;
      for (const nb of mol.neighbors(i)) (dx += a.x - mol.atoms[nb].x), (dy += a.y - mol.atoms[nb].y);
      const l = Math.hypot(dx, dy);
      a.hCount = hc - 1;
      const hi = mol.atoms.length;
      mol.atoms.push({ id: a.id, el: 'H', charge: 0, x: a.x + (l > 1e-6 ? dx / l : 0) * 0.9, y: a.y + (l > 1e-6 ? dy / l : -1) * 0.9 });
      mol.bonds.push({ id: -1, a: i, b: hi, order: 1, style: 'plain' });
      mol.invalidate();
      hydrogenOf.set(i, hi);
    }
  const n = mol.atoms.length;
  // rings drawn with delocalised bonds: the bookkeeping needs a Kekulé structure, with the bonds that arrows
  // start from as double bonds where possible; rings no arrow touches are drawn delocalised again afterwards
  const aromatic = mol.bonds.map((b) => b.order === 1.5);
  const aromaticSystems: { keys: string[]; atoms: number[] }[] = [];
  if (aromatic.some(Boolean)) {
    const fromBonds = arrows.flatMap((c) => (c.from.type === 'bond' ? [c.from.id] : [])).map((id) => mol.bonds.findIndex((b) => b.id === id));
    if (!kekuliseForArrows(mol, aromatic, fromBonds))
      error('This ring could not be given alternating double bonds — draw it in Kekulé form', mol.atoms.filter((_, i) => mol.adj[i].some((bi) => aromatic[bi])).map((a) => a.id));
    // connected sets of delocalised bonds
    const seen = new Set<number>();
    mol.bonds.forEach((_, start) => {
      if (!aromatic[start] || seen.has(start)) return;
      const stack = [start], keys: string[] = [], atoms = new Set<number>();
      seen.add(start);
      while (stack.length) {
        const bi = stack.pop()!;
        const b = mol.bonds[bi];
        keys.push(b.a < b.b ? `${b.a},${b.b}` : `${b.b},${b.a}`);
        for (const x of [b.a, b.b]) {
          atoms.add(x);
          for (const nb of mol.adj[x]) if (aromatic[nb] && !seen.has(nb)) (seen.add(nb), stack.push(nb));
        }
      }
      aromaticSystems.push({ keys, atoms: [...atoms] });
    });
  }
  const H = mol.atoms.map((_, i) => implicitH(mol, i));
  const N0 = mol.atoms.map((_, i) => nonBondingElectrons(mol, i));
  // generic atoms (E⁺, Nu⁻, R) that arrows touch count as one-bond atoms: E⁺ has no electrons to give, Nu⁻ a lone pair
  const pseudoTouched = new Set<number>();
  for (const c of arrows)
    for (const an of [c.from, c.to]) {
      const b = an.type === 'bond' ? doc.bonds.get(an.id) : undefined;
      const ids = an.type === 'atom' ? [an.id] : b ? [b.a, b.b] : an.type === 'between' ? [an.a, an.b] : [];
      for (const id of ids) {
        const i = index.get(id);
        if (i !== undefined && !mol.atoms[i].abbrev && !element(mol.atoms[i].el)) pseudoTouched.add(i);
      }
    }
  for (const i of pseudoTouched) N0[i] = Math.max(0, 1 - mol.atoms[i].charge - bondOrderSum(mol, i));
  const charge0 = mol.atoms.reduce((s, a) => s + (a.charge || 0), 0);
  // bond electrons keyed by "i,j" (i<j)
  const key = (i: number, j: number) => (i < j ? `${i},${j}` : `${j},${i}`);
  const before = new Map<string, number>();
  const bondStyle = new Map<string, { style: string; id: number }>();
  // dative and hydrogen bonds hold no electrons of their own (the donor keeps its lone pair, as in valence
  // counting): they are copied to the product unless an arrow breaks them
  const zeroBonds = new Map<string, Mol['bonds'][number]>();
  for (const b of mol.bonds) {
    const k = key(b.a, b.b);
    if (bondValence(b.order, b.style) === 0) {
      zeroBonds.set(k, b);
      continue;
    }
    before.set(k, Math.round(b.order * 2));
    bondStyle.set(k, { style: b.style, id: b.id });
  }
  const brokenZero = new Set<string>();
  const idx = (id: number) => index.get(id);
  const nearer = (pair: [number, number], target: number): number => {
    const t = mol.atoms[target];
    const [a, b] = pair;
    const da = Math.hypot(mol.atoms[a].x - t.x, mol.atoms[a].y - t.y);
    const db = Math.hypot(mol.atoms[b].x - t.x, mol.atoms[b].y - t.y);
    return da <= db ? a : b;
  };

  // ── chemistry used to resolve ambiguous arrows ──
  const order0 = (i: number, j: number) => (before.get(key(i, j)) ?? 0) / 2;
  const en = (i: number) => element(mol.atoms[i].el)?.en || 2.5;
  const isH = (i: number) => mol.atoms[i].el === 'H';
  const others = (i: number, except: number) => mol.neighbors(i).filter((m) => m !== except);
  /** How well atom p carries a positive charge or an unpaired electron: substitution, adjacent lone pairs and π bonds. */
  const cationStability = (p: number, except: number) => {
    let s = 0;
    for (const m of others(p, except)) {
      if (!isH(m)) s += 1;
      if (['N', 'O', 'S'].includes(mol.atoms[m].el) && N0[m] >= 2) s += 2;
      if (others(m, p).some((q) => order0(m, q) >= 2)) s += 1.5;
    }
    return s;
  };
  /** Atom p is next to an acceptor π bond (C=O, C≡N, N=O, …), so it can take a negative charge. */
  const besideAcceptor = (p: number, except: number) =>
    others(p, except).some((m) => others(m, p).some((q) => order0(m, q) >= 2 && en(q) > en(m) + 0.3));

  /** Bond i–j pushes its electrons toward atom t (not in the bond): which end bonds to t? → [preferred end, weight]. */
  const bondToAtomEnd = (i: number, j: number, t: number): [number, number] => {
    const bi = order0(i, t) > 0, bj = order0(j, t) > 0;
    const pi = order0(i, j) >= 2;
    if (bi !== bj) {
      const bonded = bi ? i : j, other = bi ? j : i;
      // π electrons shift toward a neighbour (conjugation); a σ bond next to t migrates its far end (1,2-shift)
      return [pi ? bonded : other, 10];
    }
    if (isH(i) !== isH(j)) return [isH(i) ? i : j, 10]; // hydrogen is what transfers
    const d = en(i) - en(j);
    if (Math.abs(d) >= 0.3) return [d > 0 ? i : j, 5]; // the electrons sit on the more electronegative end
    if (pi) {
      // Markovnikov: bond at the end that leaves the more stable cation (or radical) behind
      const si = cationStability(j, i), sj = cationStability(i, j);
      if (si !== sj) return [si > sj ? i : j, 3];
    }
    return [nearer([i, j], t), 0];
  };
  /** A lone pair on s attacks bond i–j (s not in it): which end does s bond to? → [preferred end, weight]. */
  const atomToBondEnd = (s: number, i: number, j: number): [number, number] => {
    const bi = order0(s, i) > 0, bj = order0(s, j) > 0;
    if (bi !== bj) return [bi ? i : j, 10]; // a π bond to the neighbour
    const d = en(i) - en(j);
    if (Math.abs(d) >= 0.3) return [d < 0 ? i : j, 5]; // attack the electrophilic (less electronegative) end
    const ai = besideAcceptor(j, i), aj = besideAcceptor(i, j);
    if (ai !== aj) return [ai ? i : j, 3]; // Michael: attack the end away from the acceptor
    return [nearer([i, j], s), 0];
  };

  // ── read the arrows ──
  const plans: ArrowPlan[] = [];
  for (const c of arrows) {
    const k = c.electrons;
    let src: Source;
    if (c.from.type === 'atom') {
      const i = idx(c.from.id);
      if (i === undefined) {
        missing();
        continue;
      }
      if (c.from.h) {
        // from the H of a label: the X–H bond
        const hi = hydrogenOf.get(i)!;
        if (hi < 0) continue;
        src = { bond: [i, hi] };
      } else src = { atom: i };
    } else if (c.from.type === 'bond') {
      const b = doc.bonds.get(c.from.id);
      const i = b ? idx(b.a) : undefined, j = b ? idx(b.b) : undefined;
      if (i === undefined || j === undefined) {
        missing();
        continue;
      }
      src = { bond: [i, j] };
      if (zeroBonds.has(key(i, j))) {
        // breaking a dative or hydrogen bond moves no electrons
        const t = c.to.type === 'atom' ? idx(c.to.id) : undefined;
        if (t === i || t === j) brokenZero.add(key(i, j));
        else error('An arrow from a dative or hydrogen bond can only break it — point it at one of its two atoms', [mol.atoms[i].id, mol.atoms[j].id]);
        continue;
      }
    } else {
      error('An electron-pushing arrow must start at an atom (lone pair) or a bond');
      continue;
    }
    const s = 'atom' in src ? src.atom : -1;
    const pair = 'bond' in src ? src.bond : null;
    let dest: Destination = { t: 'back' };
    if (c.to.type === 'atom') {
      let t = idx(c.to.id);
      if (t !== undefined && c.to.h) t = (hydrogenOf.get(t) ?? -1) >= 0 ? hydrogenOf.get(t) : undefined; // the H itself
      if (t === undefined) {
        if (!c.to.h) missing();
      } else if (s >= 0) {
        if (t === s) fixed.push({ message: 'An arrow starts and ends on the same atom — it moves no electrons', atomIds: [mol.atoms[s].id], level: 'warning' });
        dest = t === s ? { t: 'lp', atom: t } : { t: 'bond', a: s, b: t };
      } else if (pair) {
        if (t === pair[0] || t === pair[1]) dest = { t: 'lp', atom: t };
        else {
          const [prefer, weight] = bondToAtomEnd(pair[0], pair[1], t);
          dest = { t: 'pick', base: t, ends: pair, prefer, weight, nearest: nearer(pair, t) };
        }
      }
    } else if (c.to.type === 'bond') {
      const b = doc.bonds.get(c.to.id);
      const i = b ? idx(b.a) : undefined, j = b ? idx(b.b) : undefined;
      if (i === undefined || j === undefined) missing();
      else if (s >= 0 && s !== i && s !== j) {
        const [prefer, weight] = atomToBondEnd(s, i, j);
        dest = { t: 'pick', base: s, ends: [i, j], prefer, weight, nearest: nearer([i, j], s) };
      } else dest = { t: 'bond', a: i, b: j };
    } else if (c.to.type === 'between') {
      const i = idx(c.to.a), j = idx(c.to.b);
      if (i === undefined || j === undefined) missing();
      else if (s >= 0 && s !== i && s !== j)
        error(`An arrow from ${atomName(mol, s)} points between two other atoms — point it at the atom that ${atomName(mol, s)} bonds to`, [mol.atoms[s].id, mol.atoms[i].id, mol.atoms[j].id]);
      else if (pair && ![i, j].some((x) => x === pair[0] || x === pair[1]))
        error('An arrow from a bond points between two atoms that are not part of that bond', [mol.atoms[i].id, mol.atoms[j].id]);
      else dest = { t: 'bond', a: i, b: j };
    } else error('An arrow ends in empty space — point it at an atom, a bond or between two atoms');
    plans.push({ k, src, dest });
  }
  if (!arrows.length) error('No electron-pushing arrows to apply');

  const fishhooks = plans.some((p) => p.k === 1);

  // ── electron bookkeeping for one reading of the arrows ──
  const sig = (m: Map<string, number>) => [...m.entries()].filter(([, e]) => e > 0).map(([kk]) => kk).sort().join('|');
  const simulate = (picks: number[]): Outcome => {
    const warnings: MechanismWarning[] = [];
    const err = (message: string, atomIds?: number[]) => warnings.push({ message, atomIds, level: 'error' });
    const summary: string[] = [];
    const N = [...N0];
    const be = new Map(before);
    const addBondE = (i: number, j: number, k: number) => be.set(key(i, j), (be.get(key(i, j)) ?? 0) + k);
    let p = 0;
    for (const { k, src, dest } of plans) {
      if ('atom' in src) {
        const i = src.atom;
        if (N[i] < k) err(`${atomName(mol, i)} has no ${k === 2 ? 'lone pair' : 'non-bonding electron'} to donate`, [mol.atoms[i].id]);
        N[i] -= k;
      } else {
        const [i, j] = src.bond;
        addBondE(i, j, -k);
        if ((be.get(key(i, j)) ?? 0) < 0) err('An arrow starts from a bond with no electrons left', [mol.atoms[i].id, mol.atoms[j].id]);
      }
      if (dest.t === 'lp') N[dest.atom] += k;
      else if (dest.t === 'bond') addBondE(dest.a, dest.b, k);
      else if (dest.t === 'pick') addBondE(dest.base, picks[p++], k);
      else if ('atom' in src) N[src.atom] += k;
      else addBondE(src.bond[0], src.bond[1], k);
    }
    const changed =
      brokenZero.size > 0 ||
      N.some((v, i) => v !== N0[i]) ||
      [...new Set([...before.keys(), ...be.keys()])].some((kk) => Math.max(0, be.get(kk) ?? 0) !== (before.get(kk) ?? 0));

    // build the product
    const product = new Mol();
    for (let i = 0; i < n; i++) product.atoms.push({ ...mol.atoms[i] });
    const pairs: [number, number, number][] = [];
    for (const [kk, e] of be) {
      const [i, j] = kk.split(',').map(Number);
      let ee = e;
      if (ee % 2 !== 0) {
        err(`Odd number of electrons left between ${atomName(mol, i)} and ${atomName(mol, j)} — check the ${fishhooks ? 'fishhook arrows' : 'arrows'}`, [mol.atoms[i].id, mol.atoms[j].id]);
        ee -= 1;
        // the stray electron stays on the atom that gained fewer (after a lone fishhook both ends are radicals)
        N[N[i] - N0[i] <= N[j] - N0[j] ? i : j] += 1;
      }
      if (ee <= 0) continue;
      let order = ee / 2;
      if (order > 3) {
        err(`Bond order ${order} between ${atomName(mol, i)} and ${atomName(mol, j)}`, [mol.atoms[i].id, mol.atoms[j].id]);
        order = 3;
      }
      pairs.push([i, j, order]);
    }
    for (const [i, j, order] of pairs) {
      const kk = key(i, j);
      const old = bondStyle.get(kk);
      const bi = product.addBond(i, j, order, 'plain');
      if (old) {
        product.bonds[bi].id = old.id;
        // keep stereo/visual styles on unchanged single bonds
        if (order === 1 && (before.get(kk) ?? 0) === 2) product.bonds[bi].style = old.style as Mol['bonds'][number]['style'];
      } else product.bonds[bi].id = -1;
      const prev = (before.get(kk) ?? 0) / 2;
      if (prev === 0) summary.push(`formed ${atomName(mol, i)}–${atomName(mol, j)} bond`);
      else if (order !== prev) summary.push(`${atomName(mol, i)}–${atomName(mol, j)}: bond order ${prev} → ${order}`);
    }
    for (const [kk, e] of before) {
      if (e > 0 && !pairs.some(([i, j]) => key(i, j) === kk)) {
        const [i, j] = kk.split(',').map(Number);
        summary.push(`broke ${atomName(mol, i)}–${atomName(mol, j)} bond`);
      }
    }
    for (const [kk, b] of zeroBonds) {
      if (brokenZero.has(kk)) {
        summary.push(`broke ${atomName(mol, b.a)}–${atomName(mol, b.b)} ${b.style === 'dative' ? 'dative' : 'hydrogen'} bond`);
        continue;
      }
      const bi = product.addBond(b.a, b.b, b.order, b.style);
      product.bonds[bi].id = b.id;
    }
    product.invalidate();

    // formal charges, radicals and fixed hydrogen counts
    for (let i = 0; i < n; i++) {
      const a = product.atoms[i];
      const generic = pseudoTouched.has(i);
      if (a.abbrev || (!element(a.el) && !generic)) continue;
      const V = generic ? 1 : valenceElectrons(a.el);
      const B = bondOrderSum(product, i) + H[i];
      const Ni = Math.max(0, N[i]);
      if (N[i] < 0) err(`${atomName(mol, i)} gave away more electrons than it had`, [a.id]);
      a.charge = Math.round(V - Ni - B);
      const rad = Ni % 2;
      // keep an existing diradical/carbene designation when the electron count is unchanged
      if (rad) a.radical = 1;
      else if (N0[i] === N[i] && mol.atoms[i].radical === 2) a.radical = 2;
      else delete a.radical;
      const count = 2 * B + Ni;
      const lim = generic ? 2 : octetLimit(a.el);
      if (count > lim) err(`${atomName(mol, i)} would have ${count} valence electrons (limit ${lim})`, [a.id]);
    }
    // electrons are conserved, so the total charge must be too
    const charge1 = product.atoms.reduce((s, a) => s + (a.charge || 0), 0);
    if (charge1 !== charge0) err(`The total charge would change from ${signed(charge0)} to ${signed(charge1)} — electrons went missing or appeared`);
    // rings drawn delocalised that the arrows left alone are drawn delocalised again
    for (const sys of aromaticSystems) {
      if (sys.atoms.some((i) => N[i] !== N0[i]) || sys.keys.some((kk) => (be.get(kk) ?? 0) !== (before.get(kk) ?? 0))) continue;
      for (const kk of sys.keys) {
        const [i, j] = kk.split(',').map(Number);
        const bi = product.bondBetween(i, j);
        if (bi >= 0) product.bonds[bi].order = 1.5;
      }
    }
    // fix H counts so they don't drift
    for (let i = 0; i < n; i++) {
      const a = product.atoms[i];
      if (a.abbrev || !element(a.el)) continue;
      delete a.hCount;
      if (implicitH(product, i) !== H[i]) a.hCount = H[i];
    }
    // hydrogens taken from labels go back to being implicit H of the atom that holds them now
    const drop = new Set<number>();
    for (const hi of hydrogenOf.values()) {
      if (hi < 0) continue;
      const nb = product.neighbors(hi);
      const h = product.atoms[hi];
      if (nb.length !== 1 || h.charge || h.radical || product.bonds[product.bondBetween(hi, nb[0])].order !== 1) continue;
      const holder = product.atoms[nb[0]];
      if (holder.abbrev || !element(holder.el)) continue;
      holder.hCount = (holder.hCount ?? implicitH(product, nb[0])) + 1;
      drop.add(hi);
    }
    // label groups the arrows left unchanged become labels again (with the charge their attachment atom gained)
    const expanded: string[] = [];
    for (const g of labelGroups) {
      const inGroup = new Set(g.members);
      const inner = g.members.slice(1);
      const pairKeys = [...new Set([...before.keys(), ...be.keys()])].filter((kk) => kk.split(',').every((x) => inGroup.has(+x)));
      const intact =
        inner.every((k) => N[k] === N0[k]) &&
        pairKeys.every((kk) => Math.max(0, be.get(kk) ?? 0) === (before.get(kk) ?? 0)) &&
        product.bonds.every((b) => !(inner.includes(b.a) && !inGroup.has(b.b)) && !(inner.includes(b.b) && !inGroup.has(b.a)));
      if (!intact) {
        expanded.push(g.label);
        continue;
      }
      const a = product.atoms[g.at];
      a.charge -= g.attachCharge;
      a.el = drawn.atoms[g.at].el;
      a.abbrev = g.label;
      delete a.hCount;
      inner.forEach((k) => drop.add(k));
    }
    for (const label of expanded) warnings.push({ message: `The ${label} label is drawn in full because the arrows change atoms inside it`, level: 'warning' });
    const errors = warnings.filter((w) => w.level === 'error').length + fixed.filter((w) => w.level === 'error').length;
    // a failed or empty step is never a resonance structure
    const resonance = errors === 0 && changed && sig(before) === sig(new Map(pairs.map(([i, j, o]) => [key(i, j), o * 2])));
    const out = drop.size ? withoutAtoms(product, drop) : product;
    // keep H counts implicit wherever the valence rules give the same number
    out.atoms.forEach((a, i) => {
      if (a.hCount === undefined || a.abbrev || !element(a.el)) return;
      const hc = a.hCount;
      delete a.hCount;
      if (implicitH(out, i) !== hc) a.hCount = hc;
    });
    return { product: out, warnings, summary, changed, resonance, errors };
  };

  // ── choose the reading of ambiguous arrows ──
  const pickPlans = plans.map((p) => p.dest).filter((d): d is Extract<Destination, { t: 'pick' }> => d.t === 'pick');
  const preferred = pickPlans.map((d) => d.prefer);
  const alt = (d: Extract<Destination, { t: 'pick' }>, end: number) => (end === d.ends[0] ? d.ends[1] : d.ends[0]);
  const cost = (picks: number[]) => picks.reduce((s, e, q) => s + (e !== pickPlans[q].prefer ? pickPlans[q].weight : 0) + (e !== pickPlans[q].nearest ? 0.5 : 0), 0);
  let bestPicks = preferred;
  let best = simulate(preferred);
  if (pickPlans.length && pickPlans.length <= 8) {
    let bestScore = best.errors * 1000 + cost(preferred);
    for (let mask = 1; mask < 1 << pickPlans.length; mask++) {
      const picks = preferred.map((e, q) => (mask & (1 << q) ? alt(pickPlans[q], e) : e));
      const o = simulate(picks);
      const score = o.errors * 1000 + cost(picks);
      if (score < bestScore) {
        bestScore = score;
        best = o;
        bestPicks = picks;
      }
    }
  }
  const warnings = [...fixed, ...best.warnings];
  // an arrow that only the drawing's geometry decided, and whose other reading is just as valid
  pickPlans.forEach((d, q) => {
    if (d.weight > 0 || best.errors) return;
    const picks = [...bestPicks];
    picks[q] = alt(d, picks[q]);
    const other = simulate(picks);
    const a = smilesOf(other.product);
    if (other.errors || (a !== null && a === smilesOf(best.product))) return;
    const [x, y] = [bestPicks[q], picks[q]];
    warnings.push({
      message: `An arrow could make a ${atomName(mol, d.base)}–${atomName(mol, x)} or a ${atomName(mol, d.base)}–${atomName(mol, y)} bond; the nearer atom was used. Point the arrow at an atom (or between two atoms) to choose.`,
      atomIds: [mol.atoms[d.base].id, mol.atoms[x].id, mol.atoms[y].id],
      level: 'warning',
    });
  });
  const labelIds = new Map(labelGroups.map((g) => [drawn.atoms[g.at].id, g.label]));
  const inLabel = [...new Set(warnings.filter((w) => w.level === 'error').flatMap((w) => w.atomIds ?? []).filter((id) => labelIds.has(id)))];
  for (const id of inLabel)
    warnings.push({ message: `Arrows can only reach the ${labelIds.get(id)} label as a whole — expand it to draw arrows to the atoms inside`, atomIds: [id], level: 'warning' });
  const ok = !warnings.some((w) => w.level === 'error');
  return { product: best.product, reactantAtomIds: [...involved], warnings, resonance: ok && best.resonance, ok, changed: best.changed, summary: best.summary };
}

/**
 * Gives delocalised (order 1.5) bonds a Kekulé structure in place, trying first to make the bonds in `prefer`
 * (bonds that arrows start from) double. Implicit hydrogen counts are fixed first so they do not change.
 */
function kekuliseForArrows(mol: Mol, aromatic: boolean[], prefer: number[]): boolean {
  const atoms = mol.atoms.map(() => false);
  mol.bonds.forEach((b, bi) => {
    if (aromatic[bi]) atoms[b.a] = atoms[b.b] = true;
  });
  atoms.forEach((f, i) => {
    if (f && mol.atoms[i].hCount === undefined) mol.atoms[i].hCount = implicitH(mol, i);
  });
  const attempt = (doubles: number[]): boolean => {
    for (const bi of doubles) mol.bonds[bi].order = 2;
    if (kekulize(mol, atoms, aromatic.map((f, bi) => f && !doubles.includes(bi)))) return true;
    aromatic.forEach((f, bi) => f && (mol.bonds[bi].order = 1.5));
    return false;
  };
  const wanted = [...new Set(prefer.filter((bi) => bi >= 0 && aromatic[bi]))];
  return (wanted.length > 0 && attempt(wanted)) || attempt([]);
}

const signed = (c: number) => (c > 0 ? `+${c}` : String(c));

/** A copy of `m` without the given atoms (and their bonds). */
function withoutAtoms(m: Mol, drop: Set<number>): Mol {
  const out = new Mol();
  const map = new Map<number, number>();
  m.atoms.forEach((a, i) => {
    if (drop.has(i)) return;
    map.set(i, out.atoms.length);
    out.atoms.push(a);
  });
  for (const b of m.bonds) if (map.has(b.a) && map.has(b.b)) out.bonds.push({ ...b, a: map.get(b.a)!, b: map.get(b.b)! });
  out.invalidate();
  return out;
}

/** Canonical SMILES without explicit hydrogens (so drawn and implicit H compare equal), or null if it cannot be written. */
function smilesOf(m: Mol): string | null {
  try {
    return writeSmiles(suppressHydrogens(m));
  } catch {
    return null;
  }
}

export interface StepPlacement {
  /** The product, moved to where it goes in the drawing (atoms keep the ids of the atoms they came from). */
  mol: Mol;
  /** The reaction (→) or resonance (↔) arrow drawn between reactants and product. */
  arrow: Omit<ArrowObj, 'id' | 'type'>;
}

/** Where a step's product and arrow go: to the right of the reactants and of anything else drawn in that row. */
export function placeStep(doc: ChemDoc, r: MechanismResult): StepPlacement | null {
  const ids = new Set(r.reactantAtomIds);
  const b = docBounds(doc, ids);
  if (!b || !r.product.atoms.length) return null;
  const mol = r.product.clone();
  const inBand = (y1: number, y2: number) => y2 >= b.minY - 1.5 && y1 <= b.maxY + 1.5;
  let bandMax = b.maxX;
  for (const a of doc.atoms.values()) if (!ids.has(a.id) && a.x > b.minX && inBand(a.y, a.y)) bandMax = Math.max(bandMax, a.x);
  for (const o of doc.arrows.values()) if (Math.max(o.x1, o.x2) > b.minX && inBand(Math.min(o.y1, o.y2), Math.max(o.y1, o.y2))) bandMax = Math.max(bandMax, o.x1, o.x2);
  for (const t of doc.texts.values()) if (t.x > b.maxX && inBand(t.y, t.y)) bandMax = Math.max(bandMax, t.x + 1);
  const startX = bandMax + 0.8;
  const arrowGap = 3.4;
  // tidy bonds that were created between separate fragments
  const formedLong = mol.bonds.some((bd) => {
    const A = mol.atoms[bd.a], B = mol.atoms[bd.b];
    return Math.hypot(A.x - B.x, A.y - B.y) > 1.6;
  });
  if (formedLong) {
    try {
      clean2D(mol);
    } catch {
      /* keep raw geometry */
    }
  }
  const pb = mol.bbox();
  const midY = (b.minY + b.maxY) / 2;
  mol.translate(startX + arrowGap - pb.minX + 0.2, midY - (pb.minY + pb.maxY) / 2);
  return { mol, arrow: { kind: r.resonance ? 'resonance' : 'reaction', x1: startX, y1: midY, x2: startX + arrowGap - 0.8, y2: midY } };
}

/** Atoms whose electron count violates the octet/duet rule in the current drawing. */
export function octetViolations(doc: ChemDoc): { atomId: number; count: number; limit: number }[] {
  const { mol } = docToMol(doc);
  const out: { atomId: number; count: number; limit: number }[] = [];
  for (let i = 0; i < mol.atoms.length; i++) {
    const a = mol.atoms[i];
    if (a.abbrev || !element(a.el)) continue;
    const e = element(a.el)!;
    if (e.group < 13 && a.el !== 'H') continue;
    const B = bondOrderSum(mol, i) + implicitH(mol, i);
    const count = 2 * B + nonBondingElectrons(mol, i);
    const lim = octetLimit(a.el);
    if (count > lim) out.push({ atomId: a.id, count, limit: lim });
  }
  return out;
}
