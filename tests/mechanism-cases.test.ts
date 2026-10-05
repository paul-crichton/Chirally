// Regression suite for the arrow-pushing engine (src/doc/mechanism.ts), built from a probing exercise
// (polar, aromatic/radical/pericyclic and edge-case probes) whose findings were re-run by an
// independent verifier.
//
//   textbook mechanisms   cases the probes found correct. They must keep passing.
//   the blocks after it   cases the verifier confirmed as wrong in the first engine. Each asserts the
//                         chemically correct outcome; all of them pass since the engine was fixed. Mark a
//                         new known bug it.fails(...) with the correct expectation until it is fixed.
//
// Geometry matters. When an arrow goes from a bond to an atom outside that bond, or from a lone pair to a
// remote bond, the engine resolves the end by chemistry and only then by distance, so several cases check
// that the drawing's geometry does not decide. Coordinates are therefore copied from the probes verbatim
// (1 = one bond length, y points down). Never "tidy" them.
//
// Chained steps are applied the way the app does it (MechanismPanel.applyStep): the product is laid out,
// inserted to the right in the same document, and the next arrows are drawn on that copy.
import { describe, it, expect } from 'vitest';
import { createDoc, addAtom, addBond, insertMol, docBounds, docToMol, bondBetween } from '../src/doc/document';
import { applyArrows, arrowGroups, placeStep, MechanismResult } from '../src/doc/mechanism';
import { parseSmiles, writeSmiles, suppressHydrogens } from '../src/chem/smiles';
import { layoutMol } from '../src/chem/layout2d';
import { perceiveStereo2D } from '../src/chem/stereo2d';
import { assignCIP } from '../src/chem/cip';
import { expandAbbreviations } from '../src/chem/abbreviations';
import { setAtomLabel } from '../src/editor/ops';
import type { ChemDoc, Anchor, DocAtom } from '../src/doc/types';
import type { Mol, BondStyle } from '../src/chem/mol';

// ───────────────────────────── helpers ─────────────────────────────

type AtomSpec = [el: string, x: number, y: number, charge?: number, extra?: Partial<DocAtom>];
type BondSpec = [a: string, b: string, order?: number, style?: BondStyle];
type XY = readonly [number, number];

/** A drawing: a document plus human-readable atom labels. */
class Scene {
  readonly doc: ChemDoc = createDoc();
  /** label → document atom id. After step() the labels point at the inserted product copy. */
  readonly L: Record<string, number> = {};
  /** Offset of the latest product copy, so atom() can place atoms in the product's own frame. */
  private ox = 0;
  private oy = 0;

  atom(label: string, el: string, x: number, y: number, extra: Partial<DocAtom> = {}): number {
    const a = addAtom(this.doc, { charge: 0, ...extra, el, x: x + this.ox, y: y + this.oy });
    this.L[label] = a.id;
    return a.id;
  }

  /** Adds atoms (in key order) and then bonds, like the probes' build(). */
  build(atoms: Record<string, AtomSpec>, bonds: BondSpec[] = []): this {
    for (const [label, [el, x, y, charge, extra]] of Object.entries(atoms)) this.atom(label, el, x, y, { ...extra, charge: charge ?? 0 });
    for (const [a, b, order, style] of bonds) this.bond(a, b, order, style);
    return this;
  }

  /** Adds an atom at an offset from an existing one (optionally bonded to `bondTo`, as bond bondTo→new). */
  near(label: string, el: string, ref: string, dx: number, dy: number, opts: { bondTo?: string; charge?: number; radical?: number } = {}): number {
    const r = this.at(ref);
    const a = addAtom(this.doc, { el, x: r.x + dx, y: r.y + dy, charge: opts.charge ?? 0, ...(opts.radical ? { radical: opts.radical } : {}) });
    this.L[label] = a.id;
    if (opts.bondTo) addBond(this.doc, this.id(opts.bondTo), a.id, 1);
    return a.id;
  }

  /** An atom whose label was typed in the editor (same path as the app: setAtomLabel). */
  typed(label: string, text: string, x: number, y: number): number {
    const id = this.atom(label, 'C', x, y);
    if (!setAtomLabel(this.doc, id, text)) throw new Error(`label not understood: ${text}`);
    return id;
  }

  bond(a: string, b: string, order = 1, style: BondStyle = 'plain'): number {
    return addBond(this.doc, this.id(a), this.id(b), order, style).id;
  }

  /**
   * Inserts a SMILES; atom i is labelled `${prefix}${i}`. Coordinates are either given explicitly (one
   * [x, y] per atom, SMILES order) or laid out with layoutMol and translated exactly like the edge probes'
   * put(smiles, x0, y0). Use explicit coordinates whenever the outcome depends on geometry.
   */
  smiles(prefix: string, smiles: string, at: readonly XY[] | { x0: number; y0?: number }): string[] {
    const m = parseSmiles(smiles);
    if ('x0' in at) {
      layoutMol(m);
      const bb = m.bbox();
      m.translate(at.x0 - bb.minX, (at.y0 ?? 0) - (bb.minY + bb.maxY) / 2);
    } else {
      if (at.length !== m.atoms.length) throw new Error(`${smiles}: ${m.atoms.length} atoms, ${at.length} coordinates`);
      m.atoms.forEach((a, i) => { a.x = at[i][0]; a.y = at[i][1]; });
    }
    const { atomIds } = insertMol(this.doc, m, this.ox, this.oy);
    return atomIds.map((id, i) => {
      this.L[prefix + i] = id;
      return prefix + i;
    });
  }

  id(label: string): number {
    const id = this.L[label];
    if (id === undefined) throw new Error(`no atom labelled ${label}`);
    return id;
  }
  at(label: string): DocAtom {
    return this.doc.atoms.get(this.id(label))!;
  }
  hasBond(a: string, b: string): boolean {
    return !!bondBetween(this.doc, this.id(a), this.id(b));
  }

  // anchors
  A(l: string): Anchor {
    return { type: 'atom', id: this.id(l) };
  }
  B(l1: string, l2: string): Anchor {
    const b = bondBetween(this.doc, this.id(l1), this.id(l2));
    if (!b) throw new Error(`no bond ${l1}-${l2}`);
    return { type: 'bond', id: b.id };
  }
  BT(l1: string, l2: string): Anchor {
    return { type: 'between', a: this.id(l1), b: this.id(l2) };
  }

  /** Curved arrow: 2 electrons (full head) or 1 (fishhook). */
  arrow(from: Anchor, to: Anchor, electrons: 1 | 2 = 2): number {
    const id = this.doc.nextId++;
    this.doc.curved.set(id, { id, type: 'curved', electrons, from, to, c1: { t: 0.3, h: -0.4 }, c2: { t: 0.7, h: -0.4 } });
    return id;
  }
  fish(from: Anchor, to: Anchor): number {
    return this.arrow(from, to, 1);
  }

  /**
   * Applies the arrows like MechanismPanel.applyStep (src/app/panels/mechanism.ts): run the engine, place the
   * product with placeStep (the same DOM-free layout the panel uses), insert it into the same document, and
   * move every label onto the copy so the next step's arrows are drawn on it. r.product becomes the placed
   * product, so layout assertions measure what the app draws.
   */
  step(arrowIds: number[]): MechanismResult {
    const doc = this.doc;
    const r = applyArrows(doc, arrowIds);
    const p = placeStep(doc, r);
    if (!p) return r;
    // offset of the main species, so atom() keeps placing atoms in the product's own frame
    const main = p.mol.components().reduce((a, b) => (b.length > a.length ? b : a));
    const dx = main.reduce((t, i) => t + p.mol.atoms[i].x - r.product.atoms[i].x, 0) / main.length;
    const dy = main.reduce((t, i) => t + p.mol.atoms[i].y - r.product.atoms[i].y, 0) / main.length;
    const arrowId = doc.nextId++;
    doc.arrows.set(arrowId, { id: arrowId, type: 'arrow', ...p.arrow });
    const { atomIds } = insertMol(doc, p.mol);
    const copyOf = new Map(p.mol.atoms.map((a, i) => [a.id, atomIds[i]] as const));
    for (const [label, id] of Object.entries(this.L)) {
      const c = copyOf.get(id);
      if (c !== undefined) this.L[label] = c;
    }
    this.ox = dx;
    this.oy = dy;
    r.product = p.mol;
    return r;
  }
}

const KEKULE = [2, 1, 2, 1, 2, 1]; // C1=C2, C3=C4, C5=C6 – what the ring tool draws
const AROMATIC = [1.5, 1.5, 1.5, 1.5, 1.5, 1.5]; // the aromatic/delocalised bond tool (keys 4 / a)

/** Regular hexagon C1..C6 starting at the top, clockwise on screen (as in both probes). */
function hexagon(sc: Scene, prefix: string, cx: number, cy: number, orders: number[]): void {
  for (let k = 0; k < 6; k++) {
    const t = -Math.PI / 2 + (Math.PI / 3) * k;
    sc.atom(`${prefix}${k + 1}`, 'C', cx + Math.cos(t), cy + Math.sin(t));
  }
  for (let k = 0; k < 6; k++) sc.bond(`${prefix}${k + 1}`, `${prefix}${((k + 1) % 6) + 1}`, orders[k]);
}

/** Canonical form of a dot-separated SMILES: every fragment re-parsed and re-written, fragments sorted. */
function canon(smiles: string): string {
  return smiles.split('.').filter(Boolean).map((f) => writeSmiles(parseSmiles(f))).sort().join('.');
}

/**
 * The chemistry a structure stands for: abbreviations expanded and ordinary explicit H atoms folded into
 * their heavy atom (suppressHydrogens), as the probes did. expandAbbreviations keeps the formal charge drawn
 * on an abbreviation atom (e.g. "OMe" with −1) on its attachment atom, so a charge an abbreviation gains or
 * loses shows up in the SMILES.
 */
function chemMol(m: Mol): Mol {
  return suppressHydrogens(expandAbbreviations(m));
}

function smilesOf(m: Mol): string {
  return canon(writeSmiles(chemMol(m)));
}

function expectProduct(r: MechanismResult, smiles: string): void {
  expect(smilesOf(r.product)).toBe(canon(smiles));
}

const netCharge = (m: Mol) => m.atoms.reduce((s, a) => s + (a.charge || 0), 0);

// All assertions about warnings go through these helpers, so they can be adapted in one place.
// Both errors and advisory warnings count: a clean textbook step produces neither.
const errorsOf = (r: MechanismResult) => r.warnings.map((w) => w.message);
function expectClean(r: MechanismResult) {
  expect(errorsOf(r)).toEqual([]);
}
/** The step is reported as a problem rather than silently offered as a valid intermediate. */
function expectFlagged(r: MechanismResult) {
  expect(errorsOf(r).length).toBeGreaterThan(0);
}

/** CIP labels from wedges + 2D coordinates: centre atom id → R/S, plus E/Z labels of stereo double bonds. */
function cip(m: Mol): { centers: Map<number, string>; bonds: string[] } {
  const ex = expandAbbreviations(m.clone());
  perceiveStereo2D(ex);
  const res = assignCIP(ex);
  const centers = new Map<number, string>();
  for (const [c, l] of res.centers) centers.set(ex.atoms[c].id, l);
  return { centers, bonds: [...res.bonds.values()] };
}

/** Closest pair of non-bonded atoms, either within one fragment or between different fragments. */
function closestContact(m: Mol, within: 'same fragment' | 'different fragments'): number {
  const comp = new Map<number, number>();
  m.components().forEach((c, k) => c.forEach((i) => comp.set(i, k)));
  let best = Infinity;
  for (let i = 0; i < m.atoms.length; i++) {
    for (let j = i + 1; j < m.atoms.length; j++) {
      if (m.bondBetween(i, j) >= 0) continue;
      if ((comp.get(i) === comp.get(j)) !== (within === 'same fragment')) continue;
      best = Math.min(best, Math.hypot(m.atoms[i].x - m.atoms[j].x, m.atoms[i].y - m.atoms[j].y));
    }
  }
  return best;
}

const bondLengths = (m: Mol) => m.bonds.map((b) => Math.hypot(m.atoms[b.a].x - m.atoms[b.b].x, m.atoms[b.a].y - m.atoms[b.b].y));

// ───────────────────────────── shared drawings ─────────────────────────────

/** Acetone C=O pointing up: Cc (0,0), Oc (0,−1); `me1Right` puts Me1 at lower right (polar probe variants). */
function acetone(sc: Scene, me1Right = false): void {
  sc.build(
    { Cc: ['C', 0, 0], Oc: ['O', 0, -1], Me1: ['C', me1Right ? 0.866 : -0.866, 0.5], Me2: ['C', me1Right ? -0.866 : 0.866, 0.5] },
    [['Cc', 'Oc', 2], ['Cc', 'Me1'], ['Cc', 'Me2']],
  );
}

/** Methyl acetate with hydroxide drawn up-left (polar probe). */
function methylAcetateAndHydroxide(sc: Scene): void {
  sc.build(
    { MeAc: ['C', 0, 0], Cc: ['C', 0.866, -0.5], Oc: ['O', 0.866, -1.5], Oe: ['O', 1.732, 0], MeO: ['C', 2.598, -0.5] },
    [['MeAc', 'Cc'], ['Cc', 'Oc', 2], ['Cc', 'Oe'], ['Oe', 'MeO']],
  );
  sc.build({ Ohyd: ['O', -0.8, -2.3, -1] });
}

function tertButylCation(sc: Scene): void {
  sc.build({ Ct: ['C', 0, 0, 1], Me1: ['C', -1, 0], Me2: ['C', 0, -1], Me3: ['C', 0, 1] }, [['Ct', 'Me1'], ['Ct', 'Me2'], ['Ct', 'Me3']]);
}

/** Kekulé (or 1.5-bond) benzene with an explicit H on C1 and Br2 drawn upper-left (nearer C1) or right (nearer C2). */
function benzeneBr2(sc: Scene, orders: number[], side: 'left' | 'right'): void {
  hexagon(sc, 'C', 0, 0, orders);
  sc.atom('H1', 'H', 0, -2);
  sc.bond('C1', 'H1');
  if (side === 'left') {
    sc.atom('Bra', 'Br', -1.2, -2.2);
    sc.atom('Brb', 'Br', -2.2, -2.8);
  } else {
    sc.atom('Bra', 'Br', 2.0, -1.0);
    sc.atom('Brb', 'Br', 3.0, -1.0);
  }
  sc.bond('Bra', 'Brb');
}

/** Step 1 of electrophilic bromination: the Wheland (arenium) ion with Br on C1, C2 cationic. */
function arenium(): { sc: Scene; r: MechanismResult } {
  const sc = new Scene();
  benzeneBr2(sc, KEKULE, 'left');
  const r = sc.step([sc.arrow(sc.B('C1', 'C2'), sc.A('Bra')), sc.arrow(sc.B('Bra', 'Brb'), sc.A('Brb'))]);
  return { sc, r };
}

/** 3-methylbutan-2-yl cation, C2+ and an explicit H on C3 (aromatic-radical-pericyclic probe). */
function methylbutylCation(sc: Scene): void {
  sc.build(
    { C1: ['C', 0, 0], C2: ['C', 0.866, -0.5, 1], C3: ['C', 1.732, 0], C4: ['C', 2.598, -0.5], C5: ['C', 2.232, 0.866], H3: ['H', 1.232, 0.866] },
    [['C1', 'C2'], ['C2', 'C3'], ['C3', 'C4'], ['C3', 'C5'], ['C3', 'H3']],
  );
}

/** Neopentyl cation Cp+ with three methyls on Cq. */
function neopentylCation(sc: Scene): void {
  sc.build(
    { Cp: ['C', 0, 0, 1], Cq: ['C', 0.866, -0.5], Me1: ['C', 1.732, 0], Me2: ['C', 0.866, -1.5], Me3: ['C', 1.4, 0.35] },
    [['Cp', 'Cq'], ['Cq', 'Me1'], ['Cq', 'Me2'], ['Cq', 'Me3']],
  );
}

/** (2Z,5Z)-hepta-2,5-diene with an explicit bis-allylic H on C4, plus CH3OO• (lipid autoxidation probe). */
function lipidAndPeroxyl(sc: Scene): void {
  const P: XY[] = [[-2.6, 0.5], [-1.73, 0], [-0.87, 0.5], [0, 0], [0.87, 0.5], [1.73, 0], [2.6, 0.5]];
  P.forEach((p, k) => sc.atom(`C${k + 1}`, 'C', p[0], p[1]));
  sc.bond('C1', 'C2'); sc.bond('C2', 'C3', 2); sc.bond('C3', 'C4'); sc.bond('C4', 'C5'); sc.bond('C5', 'C6', 2); sc.bond('C6', 'C7');
  sc.atom('H4', 'H', 0, -1);
  sc.bond('C4', 'H4');
  sc.atom('Cm', 'C', 2.6, -2.0);
  sc.atom('Oa', 'O', 1.6, -1.9);
  sc.atom('Ob', 'O', 1.2, -1.0, { radical: 1 });
  sc.bond('Cm', 'Oa'); sc.bond('Oa', 'Ob');
}

/** Lipid steps 1 and 2: bis-allylic H abstraction by ROO•, then pentadienyl resonance C4• → C6•. */
function pentadienylRadical(): { sc: Scene; r1: MechanismResult; r2: MechanismResult } {
  const sc = new Scene();
  lipidAndPeroxyl(sc);
  const r1 = sc.step([sc.fish(sc.A('Ob'), sc.BT('Ob', 'H4')), sc.fish(sc.B('C4', 'H4'), sc.BT('Ob', 'H4')), sc.fish(sc.B('C4', 'H4'), sc.A('C4'))]);
  const r2 = sc.step([sc.fish(sc.A('C4'), sc.B('C4', 'C5')), sc.fish(sc.B('C5', 'C6'), sc.B('C4', 'C5')), sc.fish(sc.B('C5', 'C6'), sc.A('C6'))]);
  return { sc, r1, r2 };
}

/** 2-bromobutane, stereocentre C2 at the origin, Br up, with the wedge on C2–Br or on C2–C1; HO− drawn below (backside). */
function bromobutaneAndHydroxide(sc: Scene, wedgeOn: 'Br' | 'Me'): void {
  sc.build(
    { C2: ['C', 0, 0], C1: ['C', -0.866, 0.5], C3: ['C', 0.866, 0.5], C4: ['C', 1.732, 0], Br: ['Br', 0, -1] },
    [['C2', 'C1', 1, wedgeOn === 'Me' ? 'wedge' : 'plain'], ['C2', 'C3'], ['C3', 'C4'], ['C2', 'Br', 1, wedgeOn === 'Br' ? 'wedge' : 'plain']],
  );
  sc.build({ O: ['O', 0, 2.2, -1] });
}

/** Butan-2-ol with the wedge drawn from the stereocentre C2 to O (narrow end at C2), plus a bare H+. */
function wedgedButanolAndProton(sc: Scene): void {
  sc.build(
    { C2: ['C', 0, 0], C1: ['C', -0.866, 0.5], C3: ['C', 0.866, 0.5], C4: ['C', 1.732, 0], O: ['O', 0, -1] },
    [['C2', 'C1'], ['C2', 'C3'], ['C3', 'C4'], ['C2', 'O', 1, 'wedge']],
  );
  sc.build({ H: ['H', -1.2, -2, 1] });
}

/** Me3P=CH2 (ylene form) with P at the origin and CH2 to the right. */
function ylene(sc: Scene): void {
  sc.atom('P', 'P', 0, 0);
  ([[-0.87, -0.5], [-0.87, 0.5], [0, 1]] as const).forEach(([x, y], i) => {
    sc.atom(`Me${i}`, 'C', x, y);
    sc.bond('P', `Me${i}`);
  });
  sc.atom('CH2', 'C', 1, 0);
  sc.bond('P', 'CH2', 2);
}

// Coordinates the edge probes got from layoutMol (frozen here so the geometry-dependent cases stay put).
const PROPENE: XY[] = [[0, 0.25], [0.866, -0.25], [1.732, 0.25]];
const ACETONE_AT_3: XY[] = [[3, 0.75], [3.866, 0.25], [4.732, 0.75], [3.866, -0.75]];
const ACETONE_AT_2_2: XY[] = [[2.2, 0.75], [3.066, 0.25], [3.932, 0.75], [3.066, -0.75]];
const ACETONE_ABOVE_P: XY[] = [[-0.6, -1.45], [0.266, -1.95], [1.132, -1.45], [0.266, -2.95]];
const HYDROXYPENTANAL: XY[] = [[0, -0.25], [0.866, 0.25], [1.732, -0.25], [2.598, 0.25], [3.464, -0.25], [4.33, 0.25], [5.196, -0.25]];
const ALDOLATE: XY[] = [[0.415, 0.883], [0.883, 0], [1.352, -0.883], [0, -0.469], [1.767, 0.469], [2.614, -0.062], [3.498, 0.407]];
const HEPTYL_BROMIDE: XY[] = [[3, 0.245], [3.867, -0.254], [4.732, 0.248], [5.599, -0.251], [6.464, 0.251], [7.331, -0.248], [8.196, 0.254], [9.063, -0.245]];
const NONYL_BROMIDE: XY[] = [[3, 0.249], [3.866, -0.251], [4.732, 0.25], [5.598, -0.25], [6.464, 0.25], [7.33, -0.25], [8.196, 0.25], [9.062, -0.25], [9.928, 0.251], [10.794, -0.249]];

// ═════════════════════════ textbook cases the probes found correct ═════════════════════════

describe('textbook mechanisms', () => {
  describe('substitution and elimination', () => {
    it('SN2: HO− + CH3Br → methanol + bromide', () => {
      const sc = new Scene();
      sc.build({ O: ['O', 0, 0, -1] });
      sc.build({ C: ['C', 2, 0], Br: ['Br', 3, 0] }, [['C', 'Br']]);
      const r = sc.step([sc.arrow(sc.A('O'), sc.A('C')), sc.arrow(sc.B('C', 'Br'), sc.A('Br'))]);
      expectProduct(r, 'CO.[Br-]');
      expectClean(r);
      expect(r.resonance).toBe(false);
    });

    it('SN2 with the nucleophile arrow ending between(O, C)', () => {
      const sc = new Scene();
      sc.build({ O: ['O', 0, 0, -1] });
      sc.build({ C: ['C', 2, 0], Br: ['Br', 3, 0] }, [['C', 'Br']]);
      const r = sc.step([sc.arrow(sc.A('O'), sc.BT('O', 'C')), sc.arrow(sc.B('C', 'Br'), sc.A('Br'))]);
      expectProduct(r, 'CO.[Br-]');
      expectClean(r);
    });

    it('SN1 chain: ionisation of t-BuBr, water attack, deprotonation by a second water', () => {
      const sc = new Scene();
      sc.build(
        { C: ['C', 0, 0], Me1: ['C', -1, 0], Me2: ['C', 0, -1], Me3: ['C', 0, 1], Br: ['Br', 1, 0] },
        [['C', 'Me1'], ['C', 'Me2'], ['C', 'Me3'], ['C', 'Br']],
      );
      const r1 = sc.step([sc.arrow(sc.B('C', 'Br'), sc.A('Br'))]);
      expectProduct(r1, 'C[C+](C)C.[Br-]');
      expectClean(r1);

      sc.near('Ow', 'O', 'C', 1.6, -1.6);
      const r2 = sc.step([sc.arrow(sc.A('Ow'), sc.A('C'))]);
      expectProduct(r2, 'CC(C)(C)[OH2+]');
      expectClean(r2);

      sc.near('Hox', 'H', 'Ow', 0.8, -0.6, { bondTo: 'Ow' });
      sc.near('Ob', 'O', 'Hox', 1.5, -0.3);
      const r3 = sc.step([sc.arrow(sc.A('Ob'), sc.A('Hox')), sc.arrow(sc.B('Ow', 'Hox'), sc.A('Ow'))]);
      expectProduct(r3, 'CC(C)(C)O.[OH3+]');
      expectClean(r3);
    });

    it.each(['to the Cβ–Cα bond', 'to between(Cβ, Cα)', 'with the arrow list reversed'])('E2: ethoxide + 2-bromopropane (explicit β-H), C–H arrow %s', (variant) => {
      const sc = new Scene();
      sc.build({ MeEt: ['C', -0.8, -1.5], CH2Et: ['C', 0.066, -1], O: ['O', 0.932, -1.5, -1] }, [['MeEt', 'CH2Et'], ['CH2Et', 'O']]);
      sc.build(
        { Ca: ['C', 4.6, 0.2], Br: ['Br', 4.6, 1.2], Cb: ['C', 3.734, -0.3], Cc: ['C', 5.466, -0.3], H: ['H', 3.0, -1.0] },
        [['Ca', 'Br'], ['Ca', 'Cb'], ['Ca', 'Cc'], ['Cb', 'H']],
      );
      const ids = [
        sc.arrow(sc.A('O'), sc.A('H')),
        sc.arrow(sc.B('Cb', 'H'), variant.includes('between') ? sc.BT('Cb', 'Ca') : sc.B('Cb', 'Ca')),
        sc.arrow(sc.B('Ca', 'Br'), sc.A('Br')),
      ];
      if (variant.includes('reversed')) ids.reverse();
      const r = sc.step(ids);
      expectProduct(r, 'C=CC.CCO.[Br-]');
      expectClean(r);
    });

    it('Markovnikov HBr addition: H drawn left of the CH2 (nearer the terminal carbon) gives the 2° cation', () => {
      const sc = new Scene();
      sc.smiles('p', 'C=CC', PROPENE);
      sc.build({ H: ['H', -1.3, 0.25], Br: ['Br', -2.3, 0.25] }, [['H', 'Br']]);
      const r = sc.step([sc.arrow(sc.B('p0', 'p1'), sc.A('H')), sc.arrow(sc.B('H', 'Br'), sc.A('Br'))]);
      expectProduct(r, '[Br-].C[CH+]C');
      expectClean(r);
    });
  });

  describe('acid–base and carbonyl addition', () => {
    function aceticAcidAndHydroxide(sc: Scene, typedLabel = false): void {
      sc.build(
        { Me: ['C', 0, 0], Cc: ['C', 0.866, -0.5], Oc: ['O', 0.866, -1.5], Oa: ['O', 1.732, 0], H: ['H', 2.598, -0.5] },
        [['Me', 'Cc'], ['Cc', 'Oc', 2], ['Cc', 'Oa'], ['Oa', 'H']],
      );
      if (typedLabel) sc.typed('Ohyd', 'OH-', 4.2, -0.5);
      else sc.build({ Ohyd: ['O', 4.2, -0.5, -1] });
    }

    it('acid–base: acetic acid (explicit O–H) + HO− → acetate + water', () => {
      const sc = new Scene();
      aceticAcidAndHydroxide(sc);
      const r = sc.step([sc.arrow(sc.A('Ohyd'), sc.A('H')), sc.arrow(sc.B('Oa', 'H'), sc.A('Oa'))]);
      expectProduct(r, 'CC(=O)[O-].O');
      expectClean(r);
    });

    it('acid–base with the hydroxide typed as the label "OH-"', () => {
      const sc = new Scene();
      aceticAcidAndHydroxide(sc, true);
      const r = sc.step([sc.arrow(sc.A('Ohyd'), sc.A('H')), sc.arrow(sc.B('Oa', 'H'), sc.A('Oa'))]);
      expectProduct(r, 'CC(=O)[O-].O');
      expectClean(r);
    });

    it('acid–base drawing slip: HO− arrow on the O–H bond (H happens to be the nearer end)', () => {
      const sc = new Scene();
      aceticAcidAndHydroxide(sc);
      const r = sc.step([sc.arrow(sc.A('Ohyd'), sc.B('Oa', 'H')), sc.arrow(sc.B('Oa', 'H'), sc.A('Oa'))]);
      // The engine adds an advisory "interpreted as bond formation to the nearer atom" note; only the
      // product is asserted (see 'ambiguous arrow ends' for the layout where this rule goes wrong).
      expectProduct(r, 'CC(=O)[O-].O');
    });

    it('HCl heterolysis with an explicit H: H–Cl → Cl gives H+ and Cl−', () => {
      const sc = new Scene();
      const [h, cl] = sc.smiles('hcl', '[H]Cl', { x0: 0 });
      const r = sc.step([sc.arrow(sc.B(h, cl), sc.A(cl))]);
      expectProduct(r, '[H+].[Cl-]');
      expectClean(r);
    });

    it('carbonyl protonation by H3O+, then the C=O+ → O carbocation resonance form', () => {
      const sc = new Scene();
      acetone(sc);
      sc.build({ Hh: ['H', 0, -2.4], Oh: ['O', 0, -3.4, 1] }, [['Hh', 'Oh']]);
      const r1 = sc.step([sc.arrow(sc.A('Oc'), sc.A('Hh')), sc.arrow(sc.B('Hh', 'Oh'), sc.A('Oh'))]);
      expectProduct(r1, 'CC(C)=[OH+].O');
      expectClean(r1);
      const r2 = sc.step([sc.arrow(sc.B('Cc', 'Oc'), sc.A('Oc'))]);
      expectProduct(r2, 'C[C+](C)O');
      expect(r2.resonance).toBe(true);
      expectClean(r2);
    });

    it('carbonyl protonation by a bare H+', () => {
      const sc = new Scene();
      acetone(sc);
      sc.build({ Hp: ['H', 0, -2.4, 1] });
      const r = sc.step([sc.arrow(sc.A('Oc'), sc.A('Hp'))]);
      expectProduct(r, 'CC(C)=[OH+]');
      expectClean(r);
    });

    it('cyanohydrin: NC− adds to acetone', () => {
      const sc = new Scene();
      sc.build({ N: ['N', -3.2, -0.3], Ccn: ['C', -2.2, -0.3, -1] }, [['N', 'Ccn', 3]]);
      acetone(sc, true);
      const r = sc.step([sc.arrow(sc.A('Ccn'), sc.A('Cc')), sc.arrow(sc.B('Cc', 'Oc'), sc.A('Oc'))]);
      expectProduct(r, 'CC(C)([O-])C#N');
      expectClean(r);
    });

    it('water typed as the label "OH2" attacks t-Bu+', () => {
      const sc = new Scene();
      tertButylCation(sc);
      sc.typed('W', 'OH2', 2, 0);
      const r = sc.step([sc.arrow(sc.A('W'), sc.A('Ct'))]);
      expectProduct(r, 'CC(C)(C)[OH2+]');
      expectClean(r);
    });

    it('intramolecular hemiacetal: 5-hydroxypentanal closes to the zwitterionic ring', () => {
      const sc = new Scene();
      const m = sc.smiles('m', 'OCCCCC=O', HYDROXYPENTANAL);
      const r = sc.step([sc.arrow(sc.A(m[0]), sc.A(m[5])), sc.arrow(sc.B(m[5], m[6]), sc.A(m[6]))]);
      expectProduct(r, '[O-]C1CCCC[OH+]1');
      expectClean(r);
    });
  });

  describe('multi-step polar chains', () => {
    it('saponification: addition, collapse expelling MeO−, then MeO− deprotonates the acid', () => {
      const sc = new Scene();
      methylAcetateAndHydroxide(sc);
      const r1 = sc.step([sc.arrow(sc.A('Ohyd'), sc.A('Cc')), sc.arrow(sc.B('Cc', 'Oc'), sc.A('Oc'))]);
      expectProduct(r1, 'COC(C)(O)[O-]');
      expectClean(r1);
      const r2 = sc.step([sc.arrow(sc.A('Oc'), sc.B('Cc', 'Oc')), sc.arrow(sc.B('Cc', 'Oe'), sc.A('Oe'))]);
      expectProduct(r2, 'CC(=O)O.C[O-]');
      expectClean(r2);
      sc.near('Hacid', 'H', 'Ohyd', 0, -0.9, { bondTo: 'Ohyd' });
      const r3 = sc.step([sc.arrow(sc.A('Oe'), sc.A('Hacid')), sc.arrow(sc.B('Ohyd', 'Hacid'), sc.A('Ohyd'))]);
      expectProduct(r3, 'CC(=O)[O-].CO');
      expectClean(r3);
    });

    it('Fischer esterification: protonation, MeOH addition, proton transfers, loss of water', () => {
      const sc = new Scene();
      sc.build(
        { MeAc: ['C', 0, 0], Cc: ['C', 0.866, -0.5], Oc: ['O', 0.866, -1.5], Oa: ['O', 1.732, 0] },
        [['MeAc', 'Cc'], ['Cc', 'Oc', 2], ['Cc', 'Oa']],
      );
      sc.build({ Hh: ['H', 0.866, -2.5], Oh: ['O', 0.866, -3.5, 1] }, [['Hh', 'Oh']]);
      const r1 = sc.step([sc.arrow(sc.A('Oc'), sc.A('Hh')), sc.arrow(sc.B('Hh', 'Oh'), sc.A('Oh'))]);
      expectProduct(r1, 'CC(O)=[OH+].O');
      expectClean(r1);

      sc.near('Om', 'O', 'Cc', -1.2, 1.8);
      sc.near('Cm', 'C', 'Om', -0.866, 0.5, { bondTo: 'Om' });
      const r2 = sc.step([sc.arrow(sc.A('Om'), sc.A('Cc')), sc.arrow(sc.B('Cc', 'Oc'), sc.A('Oc'))]);
      expectProduct(r2, 'C[OH+]C(C)(O)O');
      expectClean(r2);

      sc.near('Hox', 'H', 'Om', -0.3, 0.95, { bondTo: 'Om' });
      sc.near('Ow', 'O', 'Hox', -0.4, 1.4);
      const r3 = sc.step([sc.arrow(sc.A('Ow'), sc.A('Hox')), sc.arrow(sc.B('Om', 'Hox'), sc.A('Om'))]);
      expectProduct(r3, 'COC(C)(O)O.[OH3+]');
      expectClean(r3);

      sc.near('Hh2', 'H', 'Oa', 0.9, 0.5);
      sc.near('Oh2', 'O', 'Hh2', 1.0, 0, { bondTo: 'Hh2', charge: 1 });
      const r4 = sc.step([sc.arrow(sc.A('Oa'), sc.A('Hh2')), sc.arrow(sc.B('Hh2', 'Oh2'), sc.A('Oh2'))]);
      expectProduct(r4, 'COC(C)(O)[OH2+].O');
      expectClean(r4);

      const r5 = sc.step([sc.arrow(sc.A('Oc'), sc.B('Cc', 'Oc')), sc.arrow(sc.B('Cc', 'Oa'), sc.A('Oa'))]);
      expectProduct(r5, 'COC(C)=[OH+].O');
      expectClean(r5);
    });

    it('imine formation: addition, proton transfer, OH protonation, iminium, deprotonation', () => {
      const sc = new Scene();
      sc.build({ Cam: ['C', -1.866, 0.5], N: ['N', -1, 0] }, [['Cam', 'N']]);
      sc.build({ Cc: ['C', 0.8, 0.3], Oc: ['O', 0.8, -0.7], Me: ['C', 1.666, 0.8] }, [['Cc', 'Oc', 2], ['Cc', 'Me']]);
      const r1 = sc.step([sc.arrow(sc.A('N'), sc.A('Cc')), sc.arrow(sc.B('Cc', 'Oc'), sc.A('Oc'))]);
      expectProduct(r1, 'C[NH2+]C(C)[O-]');
      expectClean(r1);

      sc.near('HN', 'H', 'N', 0.2, -0.95, { bondTo: 'N' });
      const r2 = sc.step([sc.arrow(sc.A('Oc'), sc.A('HN')), sc.arrow(sc.B('N', 'HN'), sc.A('N'))]);
      expectProduct(r2, 'CNC(C)O');
      expectClean(r2);

      sc.near('Hh', 'H', 'Oc', 0.3, -0.95);
      sc.near('Oh', 'O', 'Hh', 0.2, -1.0, { bondTo: 'Hh', charge: 1 });
      const r3 = sc.step([sc.arrow(sc.A('Oc'), sc.A('Hh')), sc.arrow(sc.B('Hh', 'Oh'), sc.A('Oh'))]);
      expectProduct(r3, 'CNC(C)[OH2+].O');
      expectClean(r3);

      const r4 = sc.step([sc.arrow(sc.A('N'), sc.B('N', 'Cc')), sc.arrow(sc.B('Cc', 'Oc'), sc.A('Oc'))]);
      expectProduct(r4, 'C[NH+]=CC.O');
      expectClean(r4);

      sc.near('HN2', 'H', 'N', -0.3, -0.95, { bondTo: 'N' });
      sc.near('Ob', 'O', 'HN2', -0.4, -1.3);
      const r5 = sc.step([sc.arrow(sc.A('Ob'), sc.A('HN2')), sc.arrow(sc.B('N', 'HN2'), sc.A('N'))]);
      expectProduct(r5, 'CN=CC.[OH3+]');
      expectClean(r5);
    });

    it('enolate formation with an explicit α-H, then aldol on the product copy', () => {
      const sc = new Scene();
      sc.build(
        { Me: ['C', -0.866, 0.5], Cc: ['C', 0, 0], Oc: ['O', 0, -1], Ca: ['C', 0.866, 0.5], Ha: ['H', 1.6, -0.1] },
        [['Me', 'Cc'], ['Cc', 'Oc', 2], ['Cc', 'Ca'], ['Ca', 'Ha']],
      );
      sc.build({ Ohyd: ['O', 3.0, -0.6, -1] });
      const r1 = sc.step([
        sc.arrow(sc.A('Ohyd'), sc.A('Ha')),
        sc.arrow(sc.B('Ca', 'Ha'), sc.B('Ca', 'Cc')),
        sc.arrow(sc.B('Cc', 'Oc'), sc.A('Oc')),
      ]);
      expectProduct(r1, 'C=C(C)[O-].O');
      expectClean(r1);

      sc.near('Cald', 'C', 'Ca', 1.2, 1.4);
      sc.near('Oald', 'O', 'Cald', 0.866, -0.5);
      sc.near('Meald', 'C', 'Cald', 0.4, 0.95, { bondTo: 'Cald' });
      sc.bond('Cald', 'Oald', 2);
      const r2 = sc.step([
        sc.arrow(sc.A('Oc'), sc.B('Cc', 'Oc')),
        sc.arrow(sc.B('Cc', 'Ca'), sc.A('Cald')),
        sc.arrow(sc.B('Cald', 'Oald'), sc.A('Oald')),
      ]);
      expectProduct(r2, 'CC(=O)CC(C)[O-]');
      expectClean(r2);
    });

    it('saponification step 1 with the leaving group drawn as an "OMe" label (abbreviation as spectator)', () => {
      const sc = new Scene();
      sc.build({ MeAc: ['C', 0, 0], Cc: ['C', 0.866, -0.5], Oc: ['O', 0.866, -1.5] }, [['MeAc', 'Cc'], ['Cc', 'Oc', 2]]);
      sc.typed('OMe', 'OMe', 1.732, 0);
      sc.bond('Cc', 'OMe');
      sc.build({ Ohyd: ['O', -0.8, -2.3, -1] });
      const r = sc.step([sc.arrow(sc.A('Ohyd'), sc.A('Cc')), sc.arrow(sc.B('Cc', 'Oc'), sc.A('Oc'))]);
      expectProduct(r, 'COC(C)(O)[O-]');
      expectClean(r);
    });
  });

  describe('enolates: aldol, Michael, retro-aldol', () => {
    function aldehydeRight(sc: Scene): void {
      sc.build({ Cald: ['C', 2.6, 1.0], Oald: ['O', 2.6, 0], Meald: ['C', 3.466, 1.5] }, [['Cald', 'Oald', 2], ['Cald', 'Meald']]);
    }
    function enolate(sc: Scene, form: 'O-' | 'C-'): void {
      const o = form === 'O-';
      sc.build(
        { Ce: ['C', 0, 0], Oe: ['O', 0, -1, o ? -1 : 0], Mee: ['C', -0.866, 0.5], CH2: ['C', 0.866, 0.5, o ? 0 : -1] },
        [['Ce', 'Oe', o ? 1 : 2], ['Ce', 'Mee'], ['Ce', 'CH2', o ? 2 : 1]],
      );
    }

    it('aldol: O− enolate + acetaldehyde, three arrows', () => {
      const sc = new Scene();
      enolate(sc, 'O-');
      aldehydeRight(sc);
      const r = sc.step([
        sc.arrow(sc.A('Oe'), sc.B('Ce', 'Oe')),
        sc.arrow(sc.B('Ce', 'CH2'), sc.A('Cald')),
        sc.arrow(sc.B('Cald', 'Oald'), sc.A('Oald')),
      ]);
      expectProduct(r, 'CC(=O)CC(C)[O-]');
      expectClean(r);
    });

    it('aldol: carbanion enolate, two arrows', () => {
      const sc = new Scene();
      enolate(sc, 'C-');
      aldehydeRight(sc);
      const r = sc.step([sc.arrow(sc.A('CH2'), sc.A('Cald')), sc.arrow(sc.B('Cald', 'Oald'), sc.A('Oald'))]);
      expectProduct(r, 'CC(=O)CC(C)[O-]');
      expectClean(r);
    });

    it('aldol: O− enolate with the π arrow to between(CH2, C=O carbon)', () => {
      const sc = new Scene();
      enolate(sc, 'O-');
      aldehydeRight(sc);
      const r = sc.step([
        sc.arrow(sc.A('Oe'), sc.B('Ce', 'Oe')),
        sc.arrow(sc.B('Ce', 'CH2'), sc.BT('CH2', 'Cald')),
        sc.arrow(sc.B('Cald', 'Oald'), sc.A('Oald')),
      ]);
      expectProduct(r, 'CC(=O)CC(C)[O-]');
      expectClean(r);
    });

    it.each(['carbanion donor', 'O− donor (4 arrows)', 'carbanion donor, arrow 2 to between(Cα, C=O carbon)', 'carbanion donor, arrows reversed'])(
      'Michael addition of the acetone enolate to methyl vinyl ketone: %s',
      (variant) => {
        const sc = new Scene();
        const oDonor = variant.startsWith('O');
        enolate(sc, oDonor ? 'O-' : 'C-');
        sc.build(
          { Cb: ['C', 2.4, 1.0], Ca: ['C', 3.266, 0.5], Ck: ['C', 4.132, 1.0], Ok: ['O', 4.132, 2.0], Mek: ['C', 4.998, 0.5] },
          [['Cb', 'Ca', 2], ['Ca', 'Ck'], ['Ck', 'Ok', 2], ['Ck', 'Mek']],
        );
        const ids: number[] = [];
        if (oDonor) {
          ids.push(sc.arrow(sc.A('Oe'), sc.B('Ce', 'Oe')));
          ids.push(sc.arrow(sc.B('Ce', 'CH2'), sc.A('Cb')));
        } else ids.push(sc.arrow(sc.A('CH2'), sc.A('Cb')));
        ids.push(sc.arrow(sc.B('Cb', 'Ca'), variant.includes('between') ? sc.BT('Ca', 'Ck') : sc.B('Ca', 'Ck')));
        ids.push(sc.arrow(sc.B('Ck', 'Ok'), sc.A('Ok')));
        if (variant.includes('reversed')) ids.reverse();
        const r = sc.step(ids);
        expectProduct(r, 'CC(=O)CCC=C(C)[O-]');
        expectClean(r);
      },
    );

    it('conjugate addition with the HO− arrow on the remote C=C bond (the β carbon is the nearer end)', () => {
      const sc = new Scene();
      sc.smiles('ho', '[OH-]', [[-3, 0]]);
      const en = sc.smiles('en', 'C=CC(C)=O', [[0, 0.25], [0.866, 0.75], [1.732, 0.25], [2.598, 0.75], [1.732, -0.75]]);
      const r = sc.step([
        sc.arrow(sc.A('ho0'), sc.B(en[0], en[1])),
        sc.arrow(sc.B(en[0], en[1]), sc.B(en[1], en[2])),
        sc.arrow(sc.B(en[2], en[4]), sc.A(en[4])),
      ]);
      // advisory "interpreted as bond formation to the nearer atom" note is expected; product only
      expectProduct(r, 'CC([O-])=CCO');
    });

    it('retro-aldol fragmentation of the aldolate', () => {
      const sc = new Scene();
      const m = sc.smiles('m', 'CC([O-])(C)CC=O', ALDOLATE);
      const r = sc.step([sc.arrow(sc.A(m[2]), sc.B(m[1], m[2])), sc.arrow(sc.B(m[1], m[4]), sc.B(m[4], m[5])), sc.arrow(sc.B(m[5], m[6]), sc.A(m[6]))]);
      expectProduct(r, 'C=C[O-].CC(C)=O');
      expectClean(r);
    });
  });

  describe('resonance', () => {
    it('enolate: carbanion → oxyanion, then back again on the copy', () => {
      const sc = new Scene();
      sc.build({ CH2: ['C', 0, 0, -1], Cc: ['C', 0.866, -0.5], O: ['O', 0.866, -1.5], Me: ['C', 1.732, 0] }, [['CH2', 'Cc'], ['Cc', 'O', 2], ['Cc', 'Me']]);
      const r1 = sc.step([sc.arrow(sc.A('CH2'), sc.B('CH2', 'Cc')), sc.arrow(sc.B('Cc', 'O'), sc.A('O'))]);
      expectProduct(r1, 'C=C(C)[O-]');
      expect(r1.resonance).toBe(true);
      expectClean(r1);
      const r2 = sc.step([sc.arrow(sc.A('O'), sc.B('Cc', 'O')), sc.arrow(sc.B('CH2', 'Cc'), sc.A('CH2'))]);
      expectProduct(r2, '[CH2-]C(C)=O');
      expect(r2.resonance).toBe(true);
      expectClean(r2);
    });

    it('carboxylate: the charge moves from O2 to O1', () => {
      const sc = new Scene();
      sc.build({ Me: ['C', 0, 0], C: ['C', 0.866, -0.5], O1: ['O', 0.866, -1.5], O2: ['O', 1.732, 0, -1] }, [['Me', 'C'], ['C', 'O1', 2], ['C', 'O2']]);
      const r = sc.step([sc.arrow(sc.A('O2'), sc.B('C', 'O2')), sc.arrow(sc.B('C', 'O1'), sc.A('O1'))]);
      expectProduct(r, 'CC(=O)[O-]');
      expect(r.resonance).toBe(true);
      expectClean(r);
      expect(sc.at('O1').charge).toBe(-1);
      expect(sc.at('O2').charge).toBe(0);
    });

    it('phenoxide on a Kekulé ring: O− → ortho carbanion', () => {
      const sc = new Scene();
      hexagon(sc, 'C', 0, 0, KEKULE);
      sc.near('O', 'O', 'C1', 0, -1, { bondTo: 'C1', charge: -1 });
      const r = sc.step([sc.arrow(sc.A('O'), sc.B('C1', 'O')), sc.arrow(sc.B('C1', 'C2'), sc.A('C2'))]);
      expectProduct(r, 'O=C1C=CC=C[CH-]1');
      expect(r.resonance).toBe(true);
      expectClean(r);
    });

    it('benzene Kekulé ↔ Kekulé with three arrows', () => {
      const sc = new Scene();
      hexagon(sc, 'C', 0, 0, KEKULE);
      const r = sc.step([
        sc.arrow(sc.B('C1', 'C2'), sc.B('C2', 'C3')),
        sc.arrow(sc.B('C3', 'C4'), sc.B('C4', 'C5')),
        sc.arrow(sc.B('C5', 'C6'), sc.B('C6', 'C1')),
      ]);
      expectProduct(r, 'c1ccccc1');
      expect(r.resonance).toBe(true);
      expectClean(r);
    });

    it('DMSO: S=O → S+–O−', () => {
      const sc = new Scene();
      const d = sc.smiles('d', 'CS(C)=O', { x0: 0 });
      const r = sc.step([sc.arrow(sc.B(d[1], d[3]), sc.A(d[3]))]);
      expectProduct(r, 'C[S+](C)[O-]');
      expect(r.resonance).toBe(true);
      expectClean(r);
    });
  });

  describe('electrophilic aromatic substitution on Kekulé benzene', () => {
    it('attack on Br2 drawn upper-left (nearer C1) gives the arenium ion with Br on C1', () => {
      const { sc, r } = arenium();
      expectProduct(r, '[Br-].BrC1C=CC=C[CH+]1');
      expectClean(r);
      expect(sc.hasBond('C1', 'Bra')).toBe(true);
      expect(sc.at('C2').charge).toBe(1);
    });

    it('attack on Br2 drawn to the right (nearer C2) still gives a valid bromobenzenium ion', () => {
      // The verifier ruled this correct: same constitution, the user's explicit H just sits on the C+ carbon.
      const sc = new Scene();
      benzeneBr2(sc, KEKULE, 'right');
      const r = sc.step([sc.arrow(sc.B('C1', 'C2'), sc.A('Bra')), sc.arrow(sc.B('Bra', 'Brb'), sc.A('Brb'))]);
      expectProduct(r, '[Br-].BrC1C=CC=C[CH+]1');
      expectClean(r);
    });

    it('attack with a between(C1, Br) target puts Br on C1 whatever the layout', () => {
      const sc = new Scene();
      benzeneBr2(sc, KEKULE, 'right');
      const r = sc.step([sc.arrow(sc.B('C1', 'C2'), sc.BT('C1', 'Bra')), sc.arrow(sc.B('Bra', 'Brb'), sc.A('Brb'))]);
      expectProduct(r, '[Br-].BrC1C=CC=C[CH+]1');
      expectClean(r);
      expect(sc.hasBond('C1', 'Bra')).toBe(true);
    });

    it('attack with C1’s hydrogen left implicit keeps the H counts', () => {
      const sc = new Scene();
      hexagon(sc, 'C', 0, 0, KEKULE);
      sc.build({ Bra: ['Br', -1.2, -2.2], Brb: ['Br', -2.2, -2.8] }, [['Bra', 'Brb']]);
      const r = sc.step([sc.arrow(sc.B('C1', 'C2'), sc.BT('C1', 'Bra')), sc.arrow(sc.B('Bra', 'Brb'), sc.A('Brb'))]);
      expectProduct(r, '[Br-].BrC1C=CC=C[CH+]1');
      expectClean(r);
    });

    it('arenium resonance on the product copies: C2+ → C4+ → C6+', () => {
      const { sc } = arenium();
      const r2 = sc.step([sc.arrow(sc.B('C3', 'C4'), sc.B('C2', 'C3'))]);
      expectProduct(r2, 'BrC1C=C[CH+]C=C1');
      expect(r2.resonance).toBe(true);
      expectClean(r2);
      expect(sc.at('C4').charge).toBe(1);
      const r3 = sc.step([sc.arrow(sc.B('C5', 'C6'), sc.B('C4', 'C5'))]);
      expectProduct(r3, 'BrC1C=CC=C[CH+]1');
      expect(r3.resonance).toBe(true);
      expectClean(r3);
      expect(sc.at('C6').charge).toBe(1);
    });

    it('arenium resonance drawn with the arrow head on the C+ atom', () => {
      const sc = new Scene();
      hexagon(sc, 'C', 0, 0, [1, 1, 2, 1, 2, 1]);
      sc.at('C2').charge = 1;
      sc.build({ H1: ['H', 0, -2], Br: ['Br', -0.87, -1.5] }, [['C1', 'H1'], ['C1', 'Br']]);
      const r = sc.step([sc.arrow(sc.B('C3', 'C4'), sc.A('C2'))]);
      expectProduct(r, 'BrC1C=C[CH+]C=C1');
      expect(r.resonance).toBe(true);
      expectClean(r);
      expect(sc.at('C4').charge).toBe(1);
    });

    it.each(['bond C1–C2', 'atom C2', 'between(C1, C2)'])('re-aromatisation: Br− takes the ipso H, C–H electrons to %s', (target) => {
      const { sc } = arenium();
      const to = target.startsWith('bond') ? sc.B('C1', 'C2') : target.startsWith('atom') ? sc.A('C2') : sc.BT('C1', 'C2');
      const r = sc.step([sc.arrow(sc.A('Brb'), sc.A('H1')), sc.arrow(sc.B('C1', 'H1'), to)]);
      expectProduct(r, 'Br.Brc1ccccc1');
      expectClean(r);
    });

    it('nitration: NO2+ drawn above the ring, then water removes the ipso H (polar probe)', () => {
      const sc = new Scene();
      hexagon(sc, 'C', 0, 0, KEKULE);
      sc.build({ N: ['N', 0, -2.6, 1], Oa: ['O', -1, -2.6], Ob: ['O', 1, -2.6] }, [['N', 'Oa', 2], ['N', 'Ob', 2]]);
      const r1 = sc.step([sc.arrow(sc.B('C1', 'C2'), sc.A('N')), sc.arrow(sc.B('N', 'Oa'), sc.A('Oa'))]);
      expectProduct(r1, '[O-][N+](=O)C1C=CC=C[CH+]1');
      expectClean(r1);
      sc.near('Hipso', 'H', 'C1', -0.7, -0.7, { bondTo: 'C1' });
      sc.near('Ow', 'O', 'Hipso', -1.0, -0.6);
      const r2 = sc.step([sc.arrow(sc.A('Ow'), sc.A('Hipso')), sc.arrow(sc.B('C1', 'Hipso'), sc.B('C1', 'C2'))]);
      expectProduct(r2, '[O-][N+](=O)c1ccccc1.[OH3+]');
      expectClean(r2);
    });

    it('nitration with an explicit ipso H from the start, then deprotonation by water (aromatic probe)', () => {
      const sc = new Scene();
      hexagon(sc, 'C', 0, 0, KEKULE);
      sc.atom('H1', 'H', 0, -2);
      sc.bond('C1', 'H1');
      sc.build({ N: ['N', -1.3, -2.1, 1], Oa: ['O', -0.8, -2.97], Ob: ['O', -1.8, -1.23] }, [['N', 'Oa', 2], ['N', 'Ob', 2]]);
      const r1 = sc.step([sc.arrow(sc.B('C1', 'C2'), sc.A('N')), sc.arrow(sc.B('N', 'Oa'), sc.A('Oa'))]);
      expectProduct(r1, '[O-][N+](=O)C1C=CC=C[CH+]1');
      expectClean(r1);
      sc.atom('Ow', 'O', 1.2, -2.6);
      const r2 = sc.step([sc.arrow(sc.A('Ow'), sc.A('H1')), sc.arrow(sc.B('C1', 'H1'), sc.B('C1', 'C2'))]);
      expectProduct(r2, '[O-][N+](=O)c1ccccc1.[OH3+]');
      expectClean(r2);
    });

    it('SN2 on Kekulé benzyl bromide leaves the ring intact (control for the 1.5-bond case)', () => {
      const sc = new Scene();
      hexagon(sc, 'C', 0, 0, KEKULE);
      sc.build({ C7: ['C', 0, -2], Br: ['Br', 0.87, -2.5] }, [['C1', 'C7'], ['C7', 'Br']]);
      sc.build({ O: ['O', -1.2, -2.8, -1] });
      const r = sc.step([sc.arrow(sc.A('O'), sc.A('C7')), sc.arrow(sc.B('C7', 'Br'), sc.A('Br'))]);
      expectProduct(r, '[Br-].OCc1ccccc1');
      expectClean(r);
    });
  });

  describe('pericyclic', () => {
    function butadieneAndEthene(sc: Scene): void {
      sc.build(
        { C1: ['C', -0.866, 0.5], C2: ['C', -0.866, -0.5], C3: ['C', 0, -1], C4: ['C', 0.866, -0.5] },
        [['C1', 'C2', 2], ['C2', 'C3', 1], ['C3', 'C4', 2]],
      );
      sc.build({ C5: ['C', 1.4, 1.2], C6: ['C', 0.3, 1.9] }, [['C5', 'C6', 2]]);
    }

    it('Diels–Alder, arrow heads on a bond and between pairs', () => {
      const sc = new Scene();
      butadieneAndEthene(sc);
      const r = sc.step([
        sc.arrow(sc.B('C1', 'C2'), sc.B('C2', 'C3')),
        sc.arrow(sc.B('C3', 'C4'), sc.BT('C4', 'C5')),
        sc.arrow(sc.B('C5', 'C6'), sc.BT('C6', 'C1')),
      ]);
      expectProduct(r, 'C1=CCCCC1');
      expectClean(r);
    });

    it('Diels–Alder, arrow heads on the atoms', () => {
      const sc = new Scene();
      butadieneAndEthene(sc);
      const r = sc.step([
        sc.arrow(sc.B('C1', 'C2'), sc.B('C2', 'C3')),
        sc.arrow(sc.B('C3', 'C4'), sc.A('C5')),
        sc.arrow(sc.B('C5', 'C6'), sc.A('C1')),
      ]);
      expectProduct(r, 'C1=CCCCC1');
      expectClean(r);
    });

    it('Diels–Alder from a zig-zag (s-trans) drawing', () => {
      const sc = new Scene();
      sc.build(
        { C1: ['C', 0, 0], C2: ['C', 0.866, -0.5], C3: ['C', 1.732, 0], C4: ['C', 2.598, -0.5] },
        [['C1', 'C2', 2], ['C2', 'C3'], ['C3', 'C4', 2]],
      );
      sc.build({ C5: ['C', 2.4, 1.5], C6: ['C', 0.6, 1.5] }, [['C5', 'C6', 2]]);
      const r = sc.step([
        sc.arrow(sc.B('C1', 'C2'), sc.B('C2', 'C3')),
        sc.arrow(sc.B('C3', 'C4'), sc.BT('C4', 'C5')),
        sc.arrow(sc.B('C5', 'C6'), sc.BT('C6', 'C1')),
      ]);
      expectProduct(r, 'C1=CCCCC1');
      expectClean(r);
    });

    it('Diels–Alder with diene and dienophile drawn apart: the app tidies a regular ring', () => {
      const sc = new Scene();
      const d = sc.smiles('d', 'C=CC=C', { x0: 0, y0: 0 });
      const e = sc.smiles('e', 'C=C', { x0: 0.5, y0: 3 });
      const r = sc.step([
        sc.arrow(sc.B(d[0], d[1]), sc.BT(d[0], e[0])),
        sc.arrow(sc.B(e[0], e[1]), sc.BT(e[1], d[3])),
        sc.arrow(sc.B(d[2], d[3]), sc.B(d[1], d[2])),
      ]);
      expectProduct(r, 'C1=CCCCC1');
      expectClean(r);
      for (const len of bondLengths(r.product)) expect(len).toBeCloseTo(1, 1);
    });

    it('Cope [3,3]: 3-methylhexa-1,5-diene → hepta-1,5-diene', () => {
      const sc = new Scene();
      const pos: XY[] = [[0.866, 0.5], [0.866, -0.5], [0, -1], [-0.866, -0.5], [-0.866, 0.5], [0, 1]];
      pos.forEach(([x, y], k) => {
        const s = k === 0 || k === 5 ? 1.35 : 1;
        sc.atom(`C${k + 1}`, 'C', x * s, y * s);
      });
      sc.bond('C1', 'C2', 2); sc.bond('C2', 'C3'); sc.bond('C3', 'C4'); sc.bond('C4', 'C5'); sc.bond('C5', 'C6', 2);
      sc.atom('Me', 'C', 0, -2);
      sc.bond('C3', 'Me');
      const r = sc.step([
        sc.arrow(sc.B('C1', 'C2'), sc.BT('C1', 'C6')),
        sc.arrow(sc.B('C3', 'C4'), sc.B('C2', 'C3')),
        sc.arrow(sc.B('C5', 'C6'), sc.B('C4', 'C5')),
      ]);
      expectProduct(r, 'C=CCCC=CC');
      expectClean(r);
    });

    it('6π electrocyclisation: hexa-1,3,5-triene → cyclohexa-1,3-diene', () => {
      const sc = new Scene();
      const pos: XY[] = [[0.95, 0.75], [0.866, -0.5], [0, -1], [-0.866, -0.5], [-0.95, 0.75], [0, 1.6]];
      pos.forEach(([x, y], k) => sc.atom(`C${k + 1}`, 'C', x, y));
      sc.bond('C1', 'C2', 2); sc.bond('C2', 'C3'); sc.bond('C3', 'C4', 2); sc.bond('C4', 'C5'); sc.bond('C5', 'C6', 2);
      const r = sc.step([
        sc.arrow(sc.B('C1', 'C2'), sc.B('C2', 'C3')),
        sc.arrow(sc.B('C3', 'C4'), sc.B('C4', 'C5')),
        sc.arrow(sc.B('C5', 'C6'), sc.BT('C6', 'C1')),
      ]);
      expectProduct(r, 'C1=CC=CCC1');
      expectClean(r);
    });
  });

  describe('radicals', () => {
    // Cl2 homolysis with two fishhooks is already covered by tests/mechanism.test.ts.
    function methaneAndChlorine(sc: Scene): void {
      sc.build({ C: ['C', 0, 0], H: ['H', 1, 0] }, [['C', 'H']]);
      sc.atom('Cl', 'Cl', 2.4, 0, { radical: 1 });
    }

    it('chlorination: Cl• abstracts H from CH4 (fishhooks to between(H, Cl))', () => {
      const sc = new Scene();
      methaneAndChlorine(sc);
      const r = sc.step([sc.fish(sc.A('Cl'), sc.BT('H', 'Cl')), sc.fish(sc.B('C', 'H'), sc.BT('H', 'Cl')), sc.fish(sc.B('C', 'H'), sc.A('C'))]);
      expectProduct(r, '[CH3].Cl');
      expectClean(r);
      expect(sc.at('C').radical).toBe(1);
    });

    it('chlorination: H abstraction with the fishhook heads on the atoms', () => {
      const sc = new Scene();
      methaneAndChlorine(sc);
      const r = sc.step([sc.fish(sc.A('Cl'), sc.A('H')), sc.fish(sc.B('C', 'H'), sc.A('Cl')), sc.fish(sc.B('C', 'H'), sc.A('C'))]);
      expectProduct(r, '[CH3].Cl');
      expectClean(r);
    });

    it.each(['between targets', 'heads on atoms'])('chlorination: CH3• + Cl2 → CH3Cl + Cl• (%s)', (style) => {
      const sc = new Scene();
      sc.atom('C', 'C', 0, 0, { radical: 1 });
      sc.build({ Cla: ['Cl', 1.6, 0], Clb: ['Cl', 2.6, 0] }, [['Cla', 'Clb']]);
      const between = style.startsWith('between');
      const r = sc.step([
        sc.fish(sc.A('C'), between ? sc.BT('C', 'Cla') : sc.A('Cla')),
        sc.fish(sc.B('Cla', 'Clb'), between ? sc.BT('C', 'Cla') : sc.A('C')),
        sc.fish(sc.B('Cla', 'Clb'), sc.A('Clb')),
      ]);
      expectProduct(r, '[Cl].CCl');
      expectClean(r);
      expect(sc.at('Clb').radical).toBe(1);
    });

    it('termination: CH3• + •CH3 → ethane', () => {
      const sc = new Scene();
      sc.atom('Ca', 'C', 0, 0, { radical: 1 });
      sc.atom('Cb', 'C', 1.8, 0, { radical: 1 });
      const r = sc.step([sc.fish(sc.A('Ca'), sc.BT('Ca', 'Cb')), sc.fish(sc.A('Cb'), sc.BT('Ca', 'Cb'))]);
      expectProduct(r, 'CC');
      expectClean(r);
    });

    it('anti-Markovnikov radical addition: Br• + propene → secondary radical', () => {
      const sc = new Scene();
      sc.build({ C1: ['C', 0, 0], C2: ['C', 0.87, -0.5], C3: ['C', 1.73, 0] }, [['C1', 'C2', 2], ['C2', 'C3']]);
      sc.atom('Br', 'Br', -1.2, -0.8, { radical: 1 });
      const r = sc.step([sc.fish(sc.A('Br'), sc.BT('Br', 'C1')), sc.fish(sc.B('C1', 'C2'), sc.BT('Br', 'C1')), sc.fish(sc.B('C1', 'C2'), sc.A('C2'))]);
      expectProduct(r, 'C[CH]CBr');
      expectClean(r);
      expect(sc.at('C2').radical).toBe(1);
    });

    it('radical addition: Br• + ethene with three fishhooks', () => {
      const sc = new Scene();
      sc.atom('Br', 'Br', -1.5, 0, { radical: 1 });
      const e = sc.smiles('e', 'C=C', { x0: 0 });
      const r = sc.step([sc.fish(sc.A('Br'), sc.BT('Br', e[0])), sc.fish(sc.B(e[0], e[1]), sc.BT('Br', e[0])), sc.fish(sc.B(e[0], e[1]), sc.A(e[1]))]);
      expectProduct(r, '[CH2]CBr');
      expectClean(r);
    });

    it('lipid autoxidation: bis-allylic H abstraction by ROO•, then pentadienyl resonance', () => {
      const { sc, r1, r2 } = pentadienylRadical();
      expectProduct(r1, 'COO.CC=C[CH]C=CC');
      expectClean(r1);
      expectProduct(r2, 'CC=CC=C[CH]C');
      expect(r2.resonance).toBe(true);
      expectClean(r2);
      expect(sc.at('C6').radical).toBe(1);
    });

    it('lipid autoxidation: O2 drawn as the triplet •O–O• adds with two fishhooks', () => {
      const { sc } = pentadienylRadical();
      sc.near('O1', 'O', 'C6', 0.9, -1.3, { radical: 1 });
      sc.near('O2', 'O', 'C6', 1.9, -1.3, { radical: 1 });
      sc.bond('O1', 'O2');
      const r = sc.step([sc.fish(sc.A('C6'), sc.BT('C6', 'O1')), sc.fish(sc.A('O1'), sc.BT('C6', 'O1'))]);
      expectProduct(r, 'CC=CC=CC(C)O[O]');
      expectClean(r);
      expect(sc.at('O2').radical).toBe(1);
    });

    it('lipid autoxidation: O2 drawn as O=O adds with three fishhooks', () => {
      const { sc } = pentadienylRadical();
      sc.near('O1', 'O', 'C6', 0.9, -1.3);
      sc.near('O2', 'O', 'C6', 1.9, -1.3);
      sc.bond('O1', 'O2', 2);
      const r = sc.step([sc.fish(sc.A('C6'), sc.BT('C6', 'O1')), sc.fish(sc.B('O1', 'O2'), sc.BT('C6', 'O1')), sc.fish(sc.B('O1', 'O2'), sc.A('O2'))]);
      expectProduct(r, 'CC=CC=CC(C)O[O]');
      expectClean(r);
    });

    it('phenolic antioxidant: H donation to ROO•, then phenoxyl resonance O• → ortho C• → para C•', () => {
      const sc = new Scene();
      hexagon(sc, 'C', 0, 0, KEKULE);
      sc.build({ O: ['O', 0, -2], H: ['H', 0.87, -2.5] }, [['C1', 'O'], ['O', 'H']]);
      sc.atom('Ob', 'O', 2.1, -2.6, { radical: 1 });
      sc.build({ Oa: ['O', 3.0, -2.1], Cm: ['C', 3.9, -2.6] }, [['Ob', 'Oa'], ['Oa', 'Cm']]);
      const r1 = sc.step([sc.fish(sc.A('Ob'), sc.BT('Ob', 'H')), sc.fish(sc.B('O', 'H'), sc.BT('Ob', 'H')), sc.fish(sc.B('O', 'H'), sc.A('O'))]);
      expectProduct(r1, 'COO.[O]c1ccccc1');
      expectClean(r1);
      expect(sc.at('O').radical).toBe(1);

      const r2 = sc.step([sc.fish(sc.A('O'), sc.B('O', 'C1')), sc.fish(sc.B('C1', 'C2'), sc.B('O', 'C1')), sc.fish(sc.B('C1', 'C2'), sc.A('C2'))]);
      expectProduct(r2, 'O=C1C=CC=C[CH]1');
      expect(r2.resonance).toBe(true);
      expectClean(r2);

      const r3 = sc.step([sc.fish(sc.A('C2'), sc.B('C2', 'C3')), sc.fish(sc.B('C3', 'C4'), sc.B('C2', 'C3')), sc.fish(sc.B('C3', 'C4'), sc.A('C4'))]);
      expectProduct(r3, 'O=C1C=C[CH]C=C1');
      expect(r3.resonance).toBe(true);
      expectClean(r3);
      expect(sc.at('C4').radical).toBe(1);
    });
  });

  describe('carbenes and nitrenes', () => {
    it('α-elimination CCl3− → :CCl2 + Cl−, then cyclopropanation with that carbene', () => {
      const sc = new Scene();
      sc.build(
        { C: ['C', 0, 0, -1], Cl1: ['Cl', 1, 0], Cl2: ['Cl', -0.5, 0.87], Cl3: ['Cl', -0.5, -0.87] },
        [['C', 'Cl1'], ['C', 'Cl2'], ['C', 'Cl3']],
      );
      const r1 = sc.step([sc.arrow(sc.B('C', 'Cl1'), sc.A('Cl1'))]);
      expectProduct(r1, '[Cl-].Cl[C]Cl');
      expectClean(r1);
      expect(sc.at('C').charge).toBe(0);

      sc.build({ Ca: ['C', -2.0, -0.5], Cb: ['C', -2.0, 0.5] }, [['Ca', 'Cb', 2]]);
      const r2 = sc.step([sc.arrow(sc.A('C'), sc.A('Ca')), sc.arrow(sc.B('Ca', 'Cb'), sc.BT('Cb', 'C'))]);
      expectProduct(r2, 'ClC1(Cl)CC1');
      expectClean(r2);
    });

    it('singlet :CCl2 (radical = 2) + ethene: lone pair → Ca, π → between(Cb, C)', () => {
      const sc = new Scene();
      sc.atom('C', 'C', 0, 0, { radical: 2 });
      sc.build({ Cl1: ['Cl', 0.87, -0.5], Cl2: ['Cl', 0.87, 0.5] }, [['C', 'Cl1'], ['C', 'Cl2']]);
      sc.build({ Ca: ['C', -1.6, -0.5], Cb: ['C', -1.6, 0.5] }, [['Ca', 'Cb', 2]]);
      const r = sc.step([sc.arrow(sc.A('C'), sc.A('Ca')), sc.arrow(sc.B('Ca', 'Cb'), sc.BT('Cb', 'C'))]);
      expectProduct(r, 'ClC1(Cl)CC1');
      expectClean(r);
    });

    it('triplet CH2 (radical = 2) adds stepwise with fishhooks → 1,3-diradical', () => {
      const sc = new Scene();
      sc.atom('C', 'C', 0, 0, { radical: 2 });
      sc.build({ Ca: ['C', -1.6, -0.5], Cb: ['C', -1.6, 0.5] }, [['Ca', 'Cb', 2]]);
      const r = sc.step([sc.fish(sc.A('C'), sc.BT('C', 'Ca')), sc.fish(sc.B('Ca', 'Cb'), sc.BT('C', 'Ca')), sc.fish(sc.B('Ca', 'Cb'), sc.A('Cb'))]);
      expectProduct(r, '[CH2]C[CH2]');
      expectClean(r);
      expect(sc.at('C').radical).toBe(1);
      expect(sc.at('Cb').radical).toBe(1);
    });

    it('nitrene from methyl azide (N2 loss), then aziridination of ethene', () => {
      const sc = new Scene();
      sc.build(
        { Cm: ['C', -1, 0], N1: ['N', 0, 0, -1], N2: ['N', 1, 0, 1], N3: ['N', 2, 0] },
        [['Cm', 'N1'], ['N1', 'N2'], ['N2', 'N3', 3]],
      );
      const r1 = sc.step([sc.arrow(sc.B('N1', 'N2'), sc.A('N2'))]);
      expectProduct(r1, 'C[N].N#N');
      expectClean(r1);
      sc.build({ Ca: ['C', -0.5, -1.6], Cb: ['C', 0.5, -1.6] }, [['Ca', 'Cb', 2]]);
      const r2 = sc.step([sc.arrow(sc.A('N1'), sc.A('Ca')), sc.arrow(sc.B('Ca', 'Cb'), sc.BT('Cb', 'N1'))]);
      expectProduct(r2, 'CN1CC1');
      expectClean(r2);
    });
  });

  describe('carbocation rearrangements with an unambiguous target', () => {
    it('1,2-hydride shift drawn to between(H, C+) (polar probe geometry)', () => {
      const sc = new Scene();
      sc.build(
        { C1: ['C', -0.866, 0.5], C2: ['C', 0, 0, 1], C3: ['C', 0.866, 0.5], C4: ['C', 1.732, 0], Me3: ['C', 0.866, 1.5], H3: ['H', 0.4, 1.3] },
        [['C1', 'C2'], ['C2', 'C3'], ['C3', 'C4'], ['C3', 'Me3'], ['C3', 'H3']],
      );
      const r = sc.step([sc.arrow(sc.B('C3', 'H3'), sc.BT('H3', 'C2'))]);
      expectProduct(r, 'CC[C+](C)C');
      expectClean(r);
    });

    it('1,2-hydride shift drawn to between(H, C+) (aromatic probe geometry)', () => {
      const sc = new Scene();
      methylbutylCation(sc);
      const r = sc.step([sc.arrow(sc.B('C3', 'H3'), sc.BT('H3', 'C2'))]);
      expectProduct(r, 'CC[C+](C)C');
      expectClean(r);
    });

    it('1,2-hydride shift with the H drawn bridging (nearer C2 than C3 is) works with the atom target', () => {
      const sc = new Scene();
      methylbutylCation(sc);
      const h = sc.at('H3');
      h.x = 1.3;
      h.y = 0.35;
      const r = sc.step([sc.arrow(sc.B('C3', 'H3'), sc.A('C2'))]);
      expectProduct(r, 'CC[C+](C)C');
      expectClean(r);
    });

    it('1,2-methyl shift (neopentyl cation) drawn to between(Me, C+)', () => {
      const sc = new Scene();
      neopentylCation(sc);
      const r = sc.step([sc.arrow(sc.B('Cq', 'Me3'), sc.BT('Me3', 'Cp'))]);
      expectProduct(r, 'CC[C+](C)C');
      expectClean(r);
    });
  });

  describe('reagents: hydrides, organometallics, P, S, B and Pd', () => {
    it('hydride delivery from BH4− (explicit B–H) to acetone', () => {
      const sc = new Scene();
      acetone(sc);
      sc.build({ Hh: ['H', 0, 1.9], B: ['B', 0, 2.9, -1] }, [['Hh', 'B']]);
      const r = sc.step([sc.arrow(sc.B('B', 'Hh'), sc.A('Cc')), sc.arrow(sc.B('Cc', 'Oc'), sc.A('Oc'))]);
      expectProduct(r, 'B.CC(C)[O-]');
      expectClean(r);
    });

    it('Grignard: MeMgBr drawn below acetone with the CH3 end nearer', () => {
      const sc = new Scene();
      acetone(sc);
      sc.build({ Cg: ['C', 0, 2.0], Mg: ['Mg', 0, 3.0], Br: ['Br', 0, 4.0] }, [['Cg', 'Mg'], ['Mg', 'Br']]);
      const r = sc.step([sc.arrow(sc.B('Cg', 'Mg'), sc.A('Cc')), sc.arrow(sc.B('Cc', 'Oc'), sc.A('Oc'))]);
      expectProduct(r, 'CC(C)(C)[O-].[Mg+]Br');
      expectClean(r);
    });

    it('Grignard with a between(CH3, C=O carbon) target works even with Mg nearer', () => {
      const sc = new Scene();
      const k = sc.smiles('k', 'CC(C)=O', ACETONE_AT_3);
      sc.build({ Mg: ['Mg', 2.266, 0.25], Cg: ['C', 2.266, -0.75], Br: ['Br', 1.266, 0.25] }, [['Cg', 'Mg'], ['Mg', 'Br']]);
      const r = sc.step([sc.arrow(sc.B('Cg', 'Mg'), sc.BT('Cg', k[1])), sc.arrow(sc.B(k[1], k[3]), sc.A(k[3]))]);
      expectProduct(r, 'CC(C)(C)[O-].[Mg+]Br');
      expectClean(r);
    });

    it('Wittig: ylide zwitterion Me3P+–CH2− adds to acetone (betaine)', () => {
      const sc = new Scene();
      const y = sc.smiles('y', 'C[P+](C)(C)[CH2-]', { x0: 0 });
      const k = sc.smiles('k', 'CC(C)=O', { x0: 4 });
      const r = sc.step([sc.arrow(sc.A(y[4]), sc.A(k[1])), sc.arrow(sc.B(k[1], k[3]), sc.A(k[3]))]);
      expectProduct(r, 'CC(C)([O-])C[P+](C)(C)C');
      expectClean(r);
    });

    it('Wittig: ylene P=CH2 with the CH2 end nearer the carbonyl carbon', () => {
      const sc = new Scene();
      ylene(sc);
      const k = sc.smiles('k', 'CC(C)=O', ACETONE_AT_2_2);
      const r = sc.step([sc.arrow(sc.B('P', 'CH2'), sc.A(k[1])), sc.arrow(sc.B(k[1], k[3]), sc.A(k[3]))]);
      expectProduct(r, 'CC(C)([O-])C[P+](C)(C)C');
      expectClean(r);
    });

    it('Wittig: betaine → oxaphosphetane (O− → P+)', () => {
      const sc = new Scene();
      const b = sc.smiles('b', 'C[P+](C)(C)CC(C)(C)[O-]', { x0: 0 });
      const r = sc.step([sc.arrow(sc.A(b[8]), sc.A(b[1]))]);
      expectProduct(r, 'CC1(C)CP(C)(C)(C)O1');
      expectClean(r);
    });

    it('Wittig: oxaphosphetane → alkene + Me3P=O (10-electron P accepted)', () => {
      const sc = new Scene();
      const o = sc.smiles('o', 'CP1(C)(C)CC(C)(C)O1', { x0: 0 });
      const r = sc.step([sc.arrow(sc.B(o[1], o[4]), sc.B(o[4], o[5])), sc.arrow(sc.B(o[5], o[8]), sc.B(o[8], o[1]))]);
      expectProduct(r, 'C=C(C)C.CP(C)(C)=O');
      expectClean(r);
    });

    it('Corey–Chaykovsky: sulfonium ylide adds to acetone, then the alkoxide closes the epoxide', () => {
      const sc = new Scene();
      const y = sc.smiles('y', 'C[S+](C)[CH2-]', { x0: 0 });
      const k = sc.smiles('k', 'CC(C)=O', { x0: 3.5 });
      const r1 = sc.step([sc.arrow(sc.A(y[3]), sc.A(k[1])), sc.arrow(sc.B(k[1], k[3]), sc.A(k[3]))]);
      expectProduct(r1, 'C[S+](C)CC(C)(C)[O-]');
      expectClean(r1);
      const r2 = sc.step([sc.arrow(sc.A(k[3]), sc.A(y[3])), sc.arrow(sc.B(y[1], y[3]), sc.A(y[1]))]);
      expectProduct(r2, 'CSC.CC1(C)CO1');
      expectClean(r2);
    });

    it('sulfonyl chloride + NH3: N → S with S=O → O keeps S at 12 electrons', () => {
      const sc = new Scene();
      const s = sc.smiles('s', 'CS(=O)(=O)Cl', { x0: 0 });
      const [n] = sc.smiles('n', 'N', { x0: -3 });
      const r = sc.step([sc.arrow(sc.A(n), sc.A(s[1])), sc.arrow(sc.B(s[1], s[2]), sc.A(s[2]))]);
      expectProduct(r, 'CS([NH3+])([O-])(=O)Cl');
      expectClean(r);
    });

    it('sulfonyl chloride + NH3 without breaking S=O is flagged (14 electrons on S)', () => {
      const sc = new Scene();
      const s = sc.smiles('s', 'CS(=O)(=O)Cl', { x0: 0 });
      const [n] = sc.smiles('n', 'N', { x0: -3 });
      const r = sc.step([sc.arrow(sc.A(n), sc.A(s[1]))]);
      expectFlagged(r);
    });

    it('Cl− + PCl5 → PCl6− is accepted; Cl− + PCl6− is flagged', () => {
      const sc = new Scene();
      const p = sc.smiles('p', 'ClP(Cl)(Cl)(Cl)Cl', { x0: 0 });
      const [cl] = sc.smiles('cl', '[Cl-]', { x0: -3 });
      const r = sc.step([sc.arrow(sc.A(cl), sc.A(p[1]))]);
      expectProduct(r, 'Cl[P-](Cl)(Cl)(Cl)(Cl)Cl');
      expectClean(r);

      const sc2 = new Scene();
      const p2 = sc2.smiles('p', 'Cl[P-](Cl)(Cl)(Cl)(Cl)Cl', { x0: 0 });
      const [cl2] = sc2.smiles('cl', '[Cl-]', { x0: -3 });
      expectFlagged(sc2.step([sc2.arrow(sc2.A(cl2), sc2.A(p2[1]))]));
    });

    it('Lewis adduct NH3 + BF3', () => {
      const sc = new Scene();
      const [n] = sc.smiles('n', 'N', { x0: 0 });
      const bf = sc.smiles('bf', 'FB(F)F', { x0: 2 });
      const r = sc.step([sc.arrow(sc.A(n), sc.A(bf[1]))]);
      expectProduct(r, '[NH3+][B-](F)(F)F');
      expectClean(r);
    });

    it('NaH: the Na–H bond heterolyses toward H to give hydride', () => {
      const sc = new Scene();
      const [na, h] = sc.smiles('nah', '[Na][H]', { x0: 0 });
      const r = sc.step([sc.arrow(sc.B(na, h), sc.A(h))]);
      expectProduct(r, '[H-].[Na+]');
      expectClean(r);
    });

    it('ionising a covalently drawn Na–OMe bond gives Na+ and MeO−', () => {
      const sc = new Scene();
      const m = sc.smiles('m', '[Na]OC', { x0: 0 });
      const r = sc.step([sc.arrow(sc.B(m[0], m[1]), sc.A(m[1]))]);
      expectProduct(r, '[Na+].C[O-]');
      expectClean(r);
    });

    it('a separate Na+ counter-ion stays out of a methoxide SN2', () => {
      const sc = new Scene();
      sc.smiles('na', '[Na+]', { x0: -2 });
      const m = sc.smiles('m', 'C[O-]', { x0: 0 });
      const cb = sc.smiles('cb', 'CBr', { x0: 3 });
      const r = sc.step([sc.arrow(sc.A(m[1]), sc.A(cb[0])), sc.arrow(sc.B(cb[0], cb[1]), sc.A(cb[1]))]);
      expectProduct(r, '[Br-].COC');
      expectClean(r);
    });

    it('Pd oxidative addition of MeBr with a spectator P→Pd dative bond', () => {
      const sc = new Scene();
      sc.atom('Pd', 'Pd', 0, 0);
      const p = sc.smiles('p', 'P(C)(C)C', { x0: -2.5, y0: 0 });
      sc.bond(p[0], 'Pd', 1, 'dative');
      const cb = sc.smiles('cb', 'CBr', { x0: 1.5 });
      const r = sc.step([sc.arrow(sc.A('Pd'), sc.A(cb[0])), sc.arrow(sc.B(cb[0], cb[1]), sc.BT(cb[1], 'Pd'))]);
      expectProduct(r, 'C[Pd](Br)P(C)(C)C');
      expectClean(r);
      expect(r.product.bonds.some((b) => b.style === 'dative')).toBe(true);
    });

    it('Pd(PMe3)2 with plain Pd–P bonds: oxidative addition of MeBr', () => {
      const sc = new Scene();
      sc.atom('Pd', 'Pd', 0, 0);
      const p1 = sc.smiles('p', 'P(C)(C)C', { x0: -2.5, y0: 0 });
      sc.bond(p1[0], 'Pd');
      const p2 = sc.smiles('q', 'P(C)(C)C', { x0: -0.5, y0: 2.5 });
      sc.bond(p2[0], 'Pd');
      const cb = sc.smiles('cb', 'CBr', { x0: 1.5 });
      const r = sc.step([sc.arrow(sc.A('Pd'), sc.A(cb[0])), sc.arrow(sc.B(cb[0], cb[1]), sc.BT(cb[1], 'Pd'))]);
      // Product only: the engine also reports "Pd would have 14 valence electrons (limit 12)", which the
      // probe considered spurious for a 16-electron Pd(II) complex (not part of the verified bug list).
      expectProduct(r, 'CP(C)(C)[Pd](C)(Br)P(C)(C)C');
    });
  });

  describe('abbreviations and stereo that work', () => {
    it('Ph label as a spectator: cyanide adds to Ph–C(=O)CH3', () => {
      const sc = new Scene();
      const [cn] = sc.smiles('cn', '[C-]#N', { x0: -3 });
      const k = sc.smiles('k', 'CC(C)=O', { x0: 0 });
      sc.at(k[2]).abbrev = 'Ph';
      const r = sc.step([sc.arrow(sc.A(cn), sc.A(k[1])), sc.arrow(sc.B(k[1], k[3]), sc.A(k[3]))]);
      expectProduct(r, 'CC([O-])(C#N)c1ccccc1');
      expectClean(r);
    });

    it('PPh3 drawn with three Ph labels: SN2 on MeI', () => {
      const sc = new Scene();
      const p = sc.smiles('p', 'P(C)(C)C', { x0: 0 });
      for (const k of [1, 2, 3]) sc.at(p[k]).abbrev = 'Ph';
      const mi = sc.smiles('mi', 'CI', { x0: 3 });
      const r = sc.step([sc.arrow(sc.A(p[0]), sc.A(mi[0])), sc.arrow(sc.B(mi[0], mi[1]), sc.A(mi[1]))]);
      expectProduct(r, '[I-].C[P+](c1ccccc1)(c1ccccc1)c1ccccc1');
      expectClean(r);
    });

    it('methoxide drawn out in full (control for the "OMe" label case)', () => {
      const sc = new Scene();
      const mo = sc.smiles('mo', 'C[O-]', { x0: 0 });
      const cb = sc.smiles('cb', 'CBr', { x0: 3 });
      const r = sc.step([sc.arrow(sc.A(mo[1]), sc.A(cb[0])), sc.arrow(sc.B(cb[0], cb[1]), sc.A(cb[1]))]);
      expectProduct(r, '[Br-].COC');
      expectClean(r);
    });

    it('tetrahedral intermediate with OMe drawn out in full collapses to acid + MeO− (control)', () => {
      const sc = new Scene();
      const t = sc.smiles('t', 'CC([O-])(O)OC', { x0: 0 });
      const r = sc.step([sc.arrow(sc.A(t[2]), sc.B(t[1], t[2])), sc.arrow(sc.B(t[1], t[4]), sc.A(t[4]))]);
      expectProduct(r, 'C[O-].CC(=O)O');
      expectClean(r);
    });

    it('the wedged 2-bromobutane and butan-2-ol drawings are (R)', () => {
      const a = new Scene();
      bromobutaneAndHydroxide(a, 'Br');
      expect(cip(docToMol(a.doc).mol).centers.get(a.id('C2'))).toBe('R');
      const b = new Scene();
      bromobutaneAndHydroxide(b, 'Me');
      expect(cip(docToMol(b.doc).mol).centers.get(b.id('C2'))).toBe('R');
      const c = new Scene();
      wedgedButanolAndProton(c);
      expect(cip(docToMol(c.doc).mol).centers.get(c.id('C2'))).toBe('R');
    });

    it('SN2 with the wedge on a spectator C–CH3 and HO− drawn on the backside gives (S)', () => {
      const sc = new Scene();
      bromobutaneAndHydroxide(sc, 'Me');
      const c2 = sc.id('C2');
      const r = sc.step([sc.arrow(sc.A('O'), sc.A('C2')), sc.arrow(sc.B('C2', 'Br'), sc.A('Br'))]);
      expectProduct(r, '[Br-].CC[C@H](C)O'); // the product carries the inverted (S) centre
      expectClean(r);
      expect(cip(r.product).centers.get(c2)).toBe('S');
    });

    it('anti-periplanar E2 of (2R,3R)-2-bromo-3-methylpentane drawn with Br wedge / H hash gives the E alkene', () => {
      const sc = new Scene();
      sc.build(
        { C1: ['C', -0.866, 0.5], C2: ['C', 0, 0], C3: ['C', 0.866, 0.5], C4: ['C', 1.732, 0], C5: ['C', 2.6, 0.5], Me: ['C', 0.866, 1.5], Br: ['Br', 0, -1], H: ['H', 1.3, -0.3] },
        [['C2', 'C1'], ['C2', 'C3'], ['C3', 'C4'], ['C4', 'C5'], ['C3', 'Me'], ['C2', 'Br', 1, 'wedge'], ['C3', 'H', 1, 'hash']],
      );
      sc.build({ O: ['O', 2.5, -1.2, -1] });
      const before = cip(docToMol(sc.doc).mol).centers;
      expect([before.get(sc.id('C2')), before.get(sc.id('C3'))]).toEqual(['R', 'R']);
      const r = sc.step([sc.arrow(sc.A('O'), sc.A('H')), sc.arrow(sc.B('C3', 'H'), sc.B('C2', 'C3')), sc.arrow(sc.B('C2', 'Br'), sc.A('Br'))]);
      expectProduct(r, 'O.[Br-].CC=C(C)CC');
      expectClean(r);
      expect(cip(r.product).bonds).toEqual(['E']);
    });
  });
});

// ═════════════════════════ bugs the probes confirmed in the first engine (all fixed) ═════════════════════════

describe('ambiguous arrow ends', () => {
  // Root cause: for a bond → atom arrow whose target is outside the bond (and for a lone pair → remote
  // bond), mechanism.ts picks the bond end nearer in space (nearer(), `addBondE(nearer(srcPair, t), t)`),
  // not the chemically sensible end.

  it('1,2-hydride shift with the textbook arrow σ(C–H) → C+ (polar probe geometry)', () => {
    // today: [H+].CC=C(C)C – an E1 alkene plus a free H+, with no warning
    const sc = new Scene();
    sc.build(
      { C1: ['C', -0.866, 0.5], C2: ['C', 0, 0, 1], C3: ['C', 0.866, 0.5], C4: ['C', 1.732, 0], Me3: ['C', 0.866, 1.5], H3: ['H', 0.4, 1.3] },
      [['C1', 'C2'], ['C2', 'C3'], ['C3', 'C4'], ['C3', 'Me3'], ['C3', 'H3']],
    );
    const r = sc.step([sc.arrow(sc.B('C3', 'H3'), sc.A('C2'))]);
    expectProduct(r, 'CC[C+](C)C');
  });

  it('1,2-hydride shift with the textbook arrow, H drawn down-left of C3 (aromatic probe geometry)', () => {
    const sc = new Scene();
    methylbutylCation(sc);
    const r = sc.step([sc.arrow(sc.B('C3', 'H3'), sc.A('C2'))]);
    expectProduct(r, 'CC[C+](C)C');
  });

  it('1,2-hydride shift with the textbook arrow, H drawn up-right of C3', () => {
    const sc = new Scene();
    sc.build(
      { C1: ['C', 0, 0], C2: ['C', 0.866, -0.5, 1], C3: ['C', 1.732, 0], C4: ['C', 2.598, -0.5], C5: ['C', 1.732, 1.0], H3: ['H', 2.2, -0.75] },
      [['C1', 'C2'], ['C2', 'C3'], ['C3', 'C4'], ['C3', 'C5'], ['C3', 'H3']],
    );
    const r = sc.step([sc.arrow(sc.B('C3', 'H3'), sc.A('C2'))]);
    expectProduct(r, 'CC[C+](C)C');
  });

  it('1,2-methyl shift in the 3,3-dimethyl-2-butyl cation with the textbook arrow', () => {
    // today: [CH3+].CC=C(C)C, no warning
    const sc = new Scene();
    sc.build(
      { C1: ['C', -0.866, 0.5], C2: ['C', 0, 0, 1], C3: ['C', 0.866, 0.5], C4: ['C', 1.732, 0], Me3a: ['C', 0.866, 1.5], Me3b: ['C', 1.3, -0.4] },
      [['C1', 'C2'], ['C2', 'C3'], ['C3', 'C4'], ['C3', 'Me3a'], ['C3', 'Me3b']],
    );
    const r = sc.step([sc.arrow(sc.B('C3', 'Me3a'), sc.A('C2'))]);
    expectProduct(r, 'CC(C)[C+](C)C');
  });

  it('1,2-methyl shift in the neopentyl cation with the textbook arrow', () => {
    // today: [CH3+].C=C(C)C, no warning
    const sc = new Scene();
    neopentylCation(sc);
    const r = sc.step([sc.arrow(sc.B('Cq', 'Me3'), sc.A('Cp'))]);
    expectProduct(r, 'CC[C+](C)C');
  });

  it('aldol: aldehyde drawn above the C=C, a little nearer the internal enolate carbon', () => {
    // today: CC([O-])[C-](C)([CH2+])=O + a 10-electron warning; the O− → C–O arrow already gives the
    // internal carbon its fourth bond, so the CH2 end must attack.
    const sc = new Scene();
    sc.build(
      { CH2: ['C', 0, 0], Ce: ['C', 1, 0], Oe: ['O', 1.5, 0.866, -1], Mee: ['C', 1.5, -0.866] },
      [['Ce', 'Oe'], ['Ce', 'Mee'], ['Ce', 'CH2', 2]],
    );
    sc.build({ Cald: ['C', 0.6, -2.0], Oald: ['O', -0.266, -2.5], Meald: ['C', 1.466, -2.5] }, [['Cald', 'Oald', 2], ['Cald', 'Meald']]);
    const r = sc.step([
      sc.arrow(sc.A('Oe'), sc.B('Ce', 'Oe')),
      sc.arrow(sc.B('Ce', 'CH2'), sc.A('Cald')),
      sc.arrow(sc.B('Cald', 'Oald'), sc.A('Oald')),
    ]);
    expectProduct(r, 'CC(=O)CC(C)[O-]');
    expectClean(r);
  });

  it('saponification slip: HO− (drawn nearer the carbonyl O) aimed at the C=O bond', () => {
    // today: CO[C+](C)[O-2]O – an O–O bond; the electrophilic end is the C whose π electrons move away
    const sc = new Scene();
    methylAcetateAndHydroxide(sc);
    const r = sc.step([sc.arrow(sc.A('Ohyd'), sc.B('Cc', 'Oc')), sc.arrow(sc.B('Cc', 'Oc'), sc.A('Oc'))]);
    expectProduct(r, 'COC(C)(O)[O-]');
  });

  it('Markovnikov HBr addition: H–Br drawn above the alkene (H 1.63 from CH2, 1.60 from CH)', () => {
    // today: [Br-].CC[CH2+] – the primary cation, no warning
    const sc = new Scene();
    sc.smiles('p', 'C=CC', PROPENE);
    sc.smiles('hb', '[H]Br', [[-0.3, -1.35], [0.7, -1.35]]);
    const r = sc.step([sc.arrow(sc.B('p0', 'p1'), sc.A('hb0')), sc.arrow(sc.B('hb0', 'hb1'), sc.A('hb1'))]);
    expectProduct(r, '[Br-].C[CH+]C');
  });

  it('Markovnikov HBr addition: H drawn directly above the terminal CH2', () => {
    const sc = new Scene();
    sc.smiles('p', 'C=CC', PROPENE);
    sc.build({ H: ['H', 0, -1.25], Br: ['Br', 0, -2.25] }, [['H', 'Br']]);
    const r = sc.step([sc.arrow(sc.B('p0', 'p1'), sc.A('H')), sc.arrow(sc.B('H', 'Br'), sc.A('Br'))]);
    expectProduct(r, '[Br-].C[CH+]C');
  });

  it('Wittig: ylene P=CH2 with the carbonyl carbon drawn nearer P (1.97) than CH2 (2.08)', () => {
    // today: CC(C)([O-])P(C)(C)(C)[CH2+] – a P–C bond and a primary carbocation, no warning
    const sc = new Scene();
    ylene(sc);
    const k = sc.smiles('k', 'CC(C)=O', ACETONE_ABOVE_P);
    const r = sc.step([sc.arrow(sc.B('P', 'CH2'), sc.A(k[1])), sc.arrow(sc.B(k[1], k[3]), sc.A(k[3]))]);
    expectProduct(r, 'CC(C)([O-])C[P+](C)(C)C');
  });

  it('Grignard: Br–Mg–CH3 drawn with Mg nearer the carbonyl carbon', () => {
    // today: [CH3+].CC(C)([O-])[Mg]Br, no warning (C–Mg polarity is ignored)
    const sc = new Scene();
    const k = sc.smiles('k', 'CC(C)=O', ACETONE_AT_3);
    sc.build({ Mg: ['Mg', 2.266, 0.25], Cg: ['C', 2.266, -0.75], Br: ['Br', 1.266, 0.25] }, [['Cg', 'Mg'], ['Mg', 'Br']]);
    const r = sc.step([sc.arrow(sc.B('Cg', 'Mg'), sc.A(k[1])), sc.arrow(sc.B(k[1], k[3]), sc.A(k[3]))]);
    expectProduct(r, 'CC(C)(C)[O-].[Mg+]Br');
  });

  it('cyclopropanation with the π arrow aimed at the carbene atom', () => {
    // today: [CH2+][CH2-]=C(Cl)Cl + a 10-electron warning – the nearer alkene carbon Ca is bonded to the
    // carbene twice; the π electrons must come from the far end Cb.
    const sc = new Scene();
    sc.atom('C', 'C', 0, 0, { radical: 2 });
    sc.build({ Cl1: ['Cl', 0.87, -0.5], Cl2: ['Cl', 0.87, 0.5] }, [['C', 'Cl1'], ['C', 'Cl2']]);
    sc.build({ Ca: ['C', -1.6, -0.3], Cb: ['C', -1.6, 0.7] }, [['Ca', 'Cb', 2]]);
    const r = sc.step([sc.arrow(sc.A('C'), sc.A('Ca')), sc.arrow(sc.B('Ca', 'Cb'), sc.A('C'))]);
    expectProduct(r, 'ClC1(Cl)CC1');
  });

  it('a between(a, b) target that excludes the arrow’s source does not bond two unrelated atoms', () => {
    // Root cause: a 'between' target adds the a–b bond and ignores the source. Today the cyanide lone pair
    // is spent while a Me–O bond forms: C[C+]1[CH3-][O-2]1.[C+]#N (only octet warnings).
    const sc = new Scene();
    sc.build({ N: ['N', -3.2, -0.3], Ccn: ['C', -2.2, -0.3, -1] }, [['N', 'Ccn', 3]]);
    acetone(sc, true);
    const r = sc.step([sc.arrow(sc.A('Ccn'), sc.BT('Me2', 'Oc')), sc.arrow(sc.B('Cc', 'Oc'), sc.A('Oc'))]);
    expectFlagged(r);
    expect(sc.hasBond('Me2', 'Oc')).toBe(false);
  });
});

describe('aromatic bonds', () => {
  // Root cause: order-1.5 bonds (aromatic bond tool, keys 4 / a) are counted as 3 electrons
  // (Math.round(order * 2)); the odd-electron fallback then breaks every ring bond, adds radicals and
  // charges, and blames fishhooks nobody drew. Any arrow on a fragment containing such a ring wrecks it.

  it('nitration of benzene drawn with 1.5 bonds gives the Wheland ion', () => {
    const sc = new Scene();
    hexagon(sc, 'C', 0, 0, AROMATIC);
    sc.build({ N: ['N', 0, -2.8, 1], Oa: ['O', -1, -2.8], Ob: ['O', 1, -2.8] }, [['N', 'Oa', 2], ['N', 'Ob', 2]]);
    const r = sc.step([sc.arrow(sc.B('C1', 'C2'), sc.A('N')), sc.arrow(sc.B('N', 'Oa'), sc.A('Oa'))]);
    expectProduct(r, '[O-][N+](=O)C1C=CC=C[CH+]1');
    expectClean(r);
  });

  it('bromination of benzene drawn with 1.5 bonds gives the arenium ion', () => {
    const sc = new Scene();
    benzeneBr2(sc, AROMATIC, 'left');
    const r = sc.step([sc.arrow(sc.B('C1', 'C2'), sc.BT('C1', 'Bra')), sc.arrow(sc.B('Bra', 'Brb'), sc.A('Brb'))]);
    expectProduct(r, '[Br-].BrC1C=CC=C[CH+]1');
  });

  it('SN2 on benzyl bromide whose ring uses 1.5 bonds leaves the (untouched) ring intact', () => {
    const sc = new Scene();
    hexagon(sc, 'C', 0, 0, AROMATIC);
    sc.build({ C7: ['C', 0, -2], Br: ['Br', 0.87, -2.5] }, [['C1', 'C7'], ['C7', 'Br']]);
    sc.build({ O: ['O', -1.2, -2.8, -1] });
    const r = sc.step([sc.arrow(sc.A('O'), sc.A('C7')), sc.arrow(sc.B('C7', 'Br'), sc.A('Br'))]);
    expectProduct(r, '[Br-].OCc1ccccc1');
    expectClean(r);
  });

  it('an arrow on a 1.5-bond ring never opens the ring or invents radicals', () => {
    const sc = new Scene();
    hexagon(sc, 'C', 0, 0, AROMATIC);
    const r = sc.step([sc.arrow(sc.B('C1', 'C2'), sc.B('C2', 'C3'))]);
    expect(r.product.atoms.filter((a) => a.radical)).toEqual([]);
    for (let k = 1; k <= 6; k++) expect(sc.hasBond(`C${k}`, `C${(k % 6) + 1}`)).toBe(true);
  });

  it('protonating phenol whose ring uses 1.5 bonds (ring is a spectator)', () => {
    const sc = new Scene();
    const ph = sc.smiles('ph', 'Oc1ccccc1', { x0: 0 });
    for (const b of sc.doc.bonds.values()) {
      if (b.order === 2 || (sc.doc.atoms.get(b.a)!.el === 'C' && sc.doc.atoms.get(b.b)!.el === 'C')) b.order = 1.5;
    }
    sc.atom('Hp', 'H', -2, 0, { charge: 1 });
    const r = sc.step([sc.arrow(sc.A(ph[0]), sc.A('Hp'))]);
    expectProduct(r, '[OH2+]c1ccccc1');
    expectClean(r);
  });
});

describe('abbreviations and labels', () => {
  // Root causes: abbreviation atoms (OMe, CO2H, …) and pseudo-atoms carry no electrons and are skipped when
  // charges are recomputed (valence.ts nonBondingElectrons → 0; mechanism.ts skips a.abbrev / unknown
  // elements), and there is no charge-conservation check. Formula-style labels ('H2O', 'CN-', …) fall back
  // to R pseudo-atoms in parseAtomLabel; 'Br+' gets two implicit H.

  it('saponification collapse with the leaving group drawn as an "OMe" label expels MeO−', () => {
    // today: step 2 gives acetic acid + a neutral OMe; the total charge goes −1 → 0 with no warning
    const sc = new Scene();
    sc.build({ MeAc: ['C', 0, 0], Cc: ['C', 0.866, -0.5], Oc: ['O', 0.866, -1.5] }, [['MeAc', 'Cc'], ['Cc', 'Oc', 2]]);
    sc.typed('OMe', 'OMe', 1.732, 0);
    sc.bond('Cc', 'OMe');
    sc.build({ Ohyd: ['O', -0.8, -2.3, -1] });
    sc.step([sc.arrow(sc.A('Ohyd'), sc.A('Cc')), sc.arrow(sc.B('Cc', 'Oc'), sc.A('Oc'))]);
    const r = sc.step([sc.arrow(sc.A('Oc'), sc.B('Cc', 'Oc')), sc.arrow(sc.B('Cc', 'OMe'), sc.A('OMe'))]);
    expectProduct(r, 'CC(=O)O.C[O-]');
    expect(netCharge(r.product)).toBe(-1);
  });

  it('tetrahedral intermediate with an "OMe" label collapses to acid + MeO−', () => {
    const sc = new Scene();
    const t = sc.smiles('t', 'CC([O-])(O)C', { x0: 0 });
    sc.at(t[4]).abbrev = 'OMe';
    const r = sc.step([sc.arrow(sc.A(t[2]), sc.B(t[1], t[2])), sc.arrow(sc.B(t[1], t[4]), sc.A(t[4]))]);
    expectProduct(r, 'C[O-].CC(=O)O');
    expect(netCharge(r.product)).toBe(-1);
  });

  it('methoxide drawn as an "OMe" label with charge −1 does SN2 on CH3Br', () => {
    // today: false "OMe has no lone pair to donate", and the −1 stays on the label (total −1 → −2)
    const sc = new Scene();
    sc.atom('OMe', 'C', 0, 0, { abbrev: 'OMe', charge: -1 });
    const cb = sc.smiles('cb', 'CBr', { x0: 2 });
    const r = sc.step([sc.arrow(sc.A('OMe'), sc.A(cb[0])), sc.arrow(sc.B(cb[0], cb[1]), sc.A(cb[1]))]);
    expectProduct(r, '[Br-].COC');
    expect(netCharge(r.product)).toBe(-1);
    expectClean(r);
  });

  it('HO− pushed into a "CO2H" label atom is flagged and keeps the charge', () => {
    // today: a pentavalent carbon (after expansion) with no warning, and the −1 disappears
    const sc = new Scene();
    sc.smiles('ho', '[OH-]', { x0: -3 });
    const m = sc.smiles('m', 'CC', { x0: 0 });
    sc.at(m[1]).abbrev = 'CO2H';
    const r = sc.step([sc.arrow(sc.A('ho0'), sc.A(m[1]))]);
    expectFlagged(r);
    expect(netCharge(r.product)).toBe(-1);
  });

  it('water typed as the label "H2O" attacks t-Bu+', () => {
    // today: [*]C(C)(C)C – 'H2O' is an R pseudo-atom with no electrons
    const sc = new Scene();
    tertButylCation(sc);
    sc.typed('W', 'H2O', 2, 0);
    const r = sc.step([sc.arrow(sc.A('W'), sc.A('Ct'))]);
    expectProduct(r, 'CC(C)(C)[OH2+]');
  });

  it('cyanide typed as the label "CN-" adds to acetone', () => {
    // today: [*]C(C)(C)[O-]
    const sc = new Scene();
    acetone(sc, true);
    sc.typed('CN', 'CN-', -2.2, -0.3);
    const r = sc.step([sc.arrow(sc.A('CN'), sc.A('Cc')), sc.arrow(sc.B('Cc', 'Oc'), sc.A('Oc'))]);
    expectProduct(r, 'CC(C)([O-])C#N');
  });

  it('bromenium typed as the label "Br+" brominates benzene', () => {
    // today: the electrophile is H2Br+ (implicit H by the isoelectronic rule) → a C–BrH2 group
    const sc = new Scene();
    hexagon(sc, 'C', 0, 0, KEKULE);
    sc.atom('H1', 'H', 0, -2);
    sc.bond('C1', 'H1');
    sc.typed('Br', 'Br+', -1.2, -2.2);
    const r = sc.step([sc.arrow(sc.B('C1', 'C2'), sc.A('Br'))]);
    expectProduct(r, 'BrC1C=CC=C[CH+]1');
  });

  it('a generic electrophile typed as "E+" no longer shows "+" once bonded to the ring', () => {
    // today: the product keeps alias 'E+' next to a C+ ring carbon, which reads as a dication
    const sc = new Scene();
    hexagon(sc, 'C', 0, 0, KEKULE);
    sc.atom('H1', 'H', 0, -2);
    sc.bond('C1', 'H1');
    sc.typed('E', 'E+', -1.2, -2.2);
    sc.step([sc.arrow(sc.B('C1', 'C2'), sc.A('E'))]);
    const e = sc.at('E');
    expect(e.charge ?? 0).toBe(0);
    expect(e.alias ?? e.el).not.toContain('+');
  });
});

describe('implicit hydrogens', () => {
  // An implicit H is an arrow anchor when the arrow starts or ends on the H of an atom's label
  // ({ type: 'atom', id, h: true }): as a target it is the H atom, as a source the X–H bond.
  const H = (an: Anchor): Anchor => (an.type === 'atom' ? { ...an, h: true } : an);

  it('deprotonating acetic acid at the H of its OH label gives acetate + water', () => {
    const sc = new Scene();
    sc.smiles('ho', '[OH-]', { x0: 0 });
    const ac = sc.smiles('ac', 'CC(=O)O', { x0: 3 });
    const r = sc.step([sc.arrow(sc.A('ho0'), H(sc.A(ac[3]))), sc.arrow(H(sc.A(ac[3])), sc.A(ac[3]))]);
    expectProduct(r, 'CC(=O)[O-].O');
    expectClean(r);
    // the transferred H is folded into water's implicit hydrogens, not left as a stray atom
    expect(r.product.atoms.some((a) => a.el === 'H')).toBe(false);
  });

  it('an arrow from a base to the acid O itself (not its H) is still flagged', () => {
    const sc = new Scene();
    sc.smiles('ho', '[OH-]', { x0: 0 });
    const ac = sc.smiles('ac', 'CC(=O)O', { x0: 3 });
    const r = sc.step([sc.arrow(sc.A('ho0'), sc.A(ac[3]))]);
    expectFlagged(r);
  });

  it('chained SN1: water deprotonates the oxonium at the H of its OH2+ label', () => {
    const sc = new Scene();
    const ox = sc.smiles('ox', 'CC(C)(C)[OH2+]', { x0: 0 });
    sc.smiles('w', 'O', { x0: 6 });
    const r = sc.step([sc.arrow(sc.A('w0'), H(sc.A(ox[4]))), sc.arrow(H(sc.A(ox[4])), sc.A(ox[4]))]);
    expectProduct(r, 'CC(C)(C)O.[OH3+]');
    expectClean(r);
  });

  it('an atom with no hydrogen to give is flagged', () => {
    const sc = new Scene();
    sc.smiles('ho', '[OH-]', { x0: 0 });
    const me = sc.smiles('me', 'COC', { x0: 3 });
    const r = sc.step([sc.arrow(sc.A('ho0'), H(sc.A(me[1]))), sc.arrow(H(sc.A(me[1])), sc.A(me[1]))]);
    expectFlagged(r);
  });
});

describe('stereochemistry', () => {
  // Root cause: the engine has no stereo model. A wedge on a broken bond simply disappears (no Walden
  // inversion), and rebuilt bonds lose their a→b direction, so a wedge pointing at the reacting atom flips.

  it('SN2 on (R)-2-bromobutane with the Br on a wedge gives (S)-butan-2-ol', () => {
    const sc = new Scene();
    bromobutaneAndHydroxide(sc, 'Br');
    const c2 = sc.id('C2');
    const r = sc.step([sc.arrow(sc.A('O'), sc.A('C2')), sc.arrow(sc.B('C2', 'Br'), sc.A('Br'))]);
    expect(cip(r.product).centers.get(c2)).toBe('S');
  });

  it('protonating (R)-butan-2-ol drawn with a wedge C2 → O keeps (R)', () => {
    const sc = new Scene();
    wedgedButanolAndProton(sc);
    const c2 = sc.id('C2');
    const r = sc.step([sc.arrow(sc.A('O'), sc.A('H'))]);
    expect(cip(r.product).centers.get(c2)).toBe('R');
  });
});

describe('dative and hydrogen bonds', () => {
  // Root cause: bond electrons are Math.round(order * 2) for every bond, while the valence code counts a
  // dative bond as 0 and drops order-0 bonds, so special bond types are double counted or lost.

  it('PMe3 leaves Pd neutral when the P→Pd dative bond returns to P', () => {
    // today: [Pd].C[P-2](C)C – the donor's lone pair is counted twice
    const sc = new Scene();
    sc.atom('Pd', 'Pd', 0, 0);
    const p = sc.smiles('p', 'P(C)(C)C', { x0: -2.5, y0: 0 });
    sc.bond(p[0], 'Pd', 1, 'dative');
    const r = sc.step([sc.arrow(sc.B(p[0], 'Pd'), sc.A(p[0]))]);
    expectProduct(r, '[Pd].CP(C)C');
    expect(sc.at(p[0]).charge).toBe(0);
  });

  it('a spectator hydrogen bond survives an unrelated protonation', () => {
    const sc = new Scene();
    const w = sc.smiles('w', '[H]O', { x0: 0 });
    const k = sc.smiles('k', 'CC(C)=O', { x0: 2 });
    sc.bond(w[0], k[3], 0, 'hbond');
    sc.atom('Hp', 'H', -2, 0, { charge: 1 });
    sc.step([sc.arrow(sc.A(w[1]), sc.A('Hp'))]);
    expect(sc.hasBond(w[0], k[3])).toBe(true);
  });
});

describe('workflow and malformed input', () => {
  // A malformed step must not be presented as a valid next intermediate: some warning must exist, and it
  // must not be labelled a resonance structure (↔). Where the engine already does that, the test is a
  // none of them may be offered as a valid intermediate or as a resonance structure.

  it('SN2 with the C–Br arrow forgotten is flagged (10-electron carbon)', () => {
    const sc = new Scene();
    sc.build({ O: ['O', 0, 0, -1] });
    sc.build({ C: ['C', 2, 0], Br: ['Br', 3, 0] }, [['C', 'Br']]);
    const r = sc.step([sc.arrow(sc.A('O'), sc.A('C'))]);
    expectFlagged(r);
    expect(r.resonance).toBe(false);
  });

  it('acid–base with the O–H → O arrow forgotten is flagged (4-electron H)', () => {
    const sc = new Scene();
    sc.build(
      { Me: ['C', 0, 0], Cc: ['C', 0.866, -0.5], Oc: ['O', 0.866, -1.5], Oa: ['O', 1.732, 0], H: ['H', 2.598, -0.5] },
      [['Me', 'Cc'], ['Cc', 'Oc', 2], ['Cc', 'Oa'], ['Oa', 'H']],
    );
    sc.build({ Ohyd: ['O', 4.2, -0.5, -1] });
    const r = sc.step([sc.arrow(sc.A('Ohyd'), sc.A('H'))]);
    expectFlagged(r);
    expect(r.resonance).toBe(false);
  });

  it('SN2 whose C–Br arrow ends in empty space is flagged', () => {
    const sc = new Scene();
    const [o] = sc.smiles('o', '[OH-]', { x0: 0 });
    const cb = sc.smiles('cb', 'CBr', { x0: 3 });
    const r = sc.step([sc.arrow(sc.A(o), sc.A(cb[0])), sc.arrow(sc.B(cb[0], cb[1]), { type: 'point', x: 7, y: 1 })]);
    expectFlagged(r);
    expect(r.resonance).toBe(false);
  });

  it('an arrow from a carbocation (no lone pair) to water is flagged', () => {
    const sc = new Scene();
    const c = sc.smiles('c', 'C[C+](C)C', { x0: 0 });
    const [w] = sc.smiles('w', 'O', { x0: 4 });
    const r = sc.step([sc.arrow(sc.A(c[1]), sc.A(w))]);
    expectFlagged(r);
    expect(r.resonance).toBe(false);
  });

  it('a bond arrow aimed at an unrelated distant atom of the same chain is flagged', () => {
    const sc = new Scene();
    const h = sc.smiles('h', 'CCCCCC', { x0: 0 });
    const r = sc.step([sc.arrow(sc.B(h[0], h[1]), sc.A(h[5]))]);
    expectFlagged(r);
  });

  it('applying no arrows is flagged', () => {
    const sc = new Scene();
    sc.smiles('m', 'CBr', { x0: 0 });
    const r = sc.step([]);
    expectFlagged(r);
    expect(r.product.atoms.length).toBe(0);
  });

  it('two independent steps side by side form two groups; the right-hand one has the larger maxX', () => {
    const sc = new Scene();
    const [o] = sc.smiles('o', '[OH-]', { x0: 0 });
    const cb = sc.smiles('cb', 'CBr', { x0: 3 });
    const left = [sc.arrow(sc.A(o), sc.A(cb[0])), sc.arrow(sc.B(cb[0], cb[1]), sc.A(cb[1]))];
    const [cn] = sc.smiles('cn', '[C-]#N', { x0: 12 });
    const k = sc.smiles('k', 'CC=O', { x0: 15 });
    const right = [sc.arrow(sc.A(cn), sc.A(k[1])), sc.arrow(sc.B(k[1], k[2]), sc.A(k[2]))];
    const groups = arrowGroups(sc.doc).sort((a, b) => b.maxX - a.maxX);
    expect(groups.map((g) => [...g.arrows].sort((x, y) => x - y))).toEqual([right, left]);
  });

  it('arrows of two successive steps drawn on shared reactants are applied together', () => {
    const sc = new Scene();
    const [o] = sc.smiles('o', '[OH-]', { x0: 0 });
    const k = sc.smiles('k', 'CC=O', { x0: 3 });
    const w = sc.smiles('w', '[H]O', { x0: 7 });
    sc.arrow(sc.A(o), sc.A(k[1]));
    sc.arrow(sc.B(k[1], k[2]), sc.A(k[2]));
    sc.arrow(sc.A(k[2]), sc.A(w[0]));
    sc.arrow(sc.B(w[0], w[1]), sc.A(w[1]));
    const groups = arrowGroups(sc.doc);
    expect(groups.length).toBe(1);
    const r = sc.step(groups[0].arrows);
    expectProduct(r, '[OH-].CC(O)O');
    expectClean(r);
  });

  it('an arrow ending in empty space is not offered as a resonance structure', () => {
    // today: warning, but resonance = true (nothing changed), so the app inserts a copy with a ↔ arrow
    const sc = new Scene();
    const [o] = sc.smiles('o', '[OH-]', { x0: 0 });
    sc.smiles('cb', 'CBr', { x0: 3 });
    const r = sc.step([sc.arrow(sc.A(o), { type: 'point', x: 1.5, y: 1 })]);
    expectFlagged(r);
    expect(r.resonance).toBe(false);
  });

  it('an arrow starting in empty space is not offered as a resonance structure', () => {
    const sc = new Scene();
    const cb = sc.smiles('cb', 'CBr', { x0: 3 });
    const r = sc.step([sc.arrow({ type: 'point', x: 0, y: 0 }, sc.A(cb[0]))]);
    expectFlagged(r);
    expect(r.resonance).toBe(false);
  });

  it('a lone pair pointed back at its own atom is not offered as a resonance structure', () => {
    // today: no warning at all and resonance = true
    const sc = new Scene();
    const [o] = sc.smiles('o', '[OH-]', { x0: 0 });
    const r = sc.step([sc.arrow(sc.A(o), sc.A(o))]);
    expectFlagged(r);
    expect(r.resonance).toBe(false);
  });

  it('pushing a lone pair into a triple bond (bond order 4) is not offered as a resonance structure', () => {
    // today: "Bond order 4" warning, two electrons silently dropped (charge −1 → +1), resonance = true
    const sc = new Scene();
    const [c, n] = sc.smiles('cn', '[C-]#N', { x0: 0 });
    const r = sc.step([sc.arrow(sc.A(c), sc.B(c, n))]);
    expectFlagged(r);
    expect(r.resonance).toBe(false);
  });

  it('a π arrow to a non-adjacent bond is not offered as a resonance structure', () => {
    // not separately listed by the verifier; same "labelled resonance" defect as the two cases above
    const sc = new Scene();
    const m = sc.smiles('m', 'C=CCCC', { x0: 0 });
    const r = sc.step([sc.arrow(sc.B(m[0], m[1]), sc.B(m[3], m[4]))]);
    expectFlagged(r);
    expect(r.resonance).toBe(false);
  });

  it('a single fishhook on Cl–Cl gives the same (mirrored) result whichever atom it points at', () => {
    // today: toward the first atom → [Cl-].[Cl+]; toward the second → [Cl].[Cl] (the odd electron always
    // goes to the lower-index atom)
    const fishhookOnce = (towardFirst: boolean) => {
      const sc = new Scene();
      const [a, b] = sc.smiles('cl', 'ClCl', { x0: 0 });
      const [target, other] = towardFirst ? [a, b] : [b, a];
      const r = sc.step([sc.fish(sc.B(a, b), sc.A(target))]);
      expectFlagged(r);
      const t = sc.at(target), o = sc.at(other);
      return { target: [t.charge, t.radical ?? 0], other: [o.charge, o.radical ?? 0] };
    };
    expect(fishhookOnce(true)).toEqual(fishhookOnce(false));
  });

  // Not reproduced: which arrow group applyStep picks when nothing is selected (a two-row scheme applies
  // row 1 again instead of the row just drawn, and a second click duplicates the product). That choice is
  // made in src/app/panels/mechanism.ts with the editor's selection, so it needs the DOM/app.
});

describe('layout', () => {
  // Measured on the product as placeStep lays it out (the panel inserts exactly that).

  it('hemiacetal ring closure from a zig-zag drawing is tidied into a regular ring', () => {
    const sc = new Scene();
    const m = sc.smiles('m', 'OCCCCC=O', HYDROXYPENTANAL);
    const r = sc.step([sc.arrow(sc.A(m[0]), sc.A(m[5])), sc.arrow(sc.B(m[5], m[6]), sc.A(m[6]))]);
    for (const len of bondLengths(r.product)) expect(len).toBeCloseTo(1, 1);
    expect(closestContact(r.product, 'same fragment')).toBeGreaterThan(1.5);
  });

  it.each([
    [7, HEPTYL_BROMIDE],
    [9, NONYL_BROMIDE],
  ] as const)('SN2 on a 1-bromoalkane (C%i) with HO− drawn far left: the alcohol is not folded onto itself', (n, coords) => {
    // today: clean2D starts from an 8–10 unit long bond and leaves the O ~1.0 from a carbon it is not
    // bonded to, so the alcohol looks like a cyclic ether
    const sc = new Scene();
    const chain = sc.smiles('c', 'C'.repeat(n) + 'Br', coords);
    sc.smiles('o', '[OH-]', [[0, 0]]);
    const r = sc.step([sc.arrow(sc.A('o0'), sc.A(chain[n - 1])), sc.arrow(sc.B(chain[n - 1], chain[n]), sc.A(chain[n]))]);
    expect(closestContact(r.product, 'same fragment')).toBeGreaterThan(1.25);
  });

  it('retro-aldol fragments are moved apart instead of staying one bond length apart', () => {
    const sc = new Scene();
    const m = sc.smiles('m', 'CC([O-])(C)CC=O', ALDOLATE);
    const r = sc.step([sc.arrow(sc.A(m[2]), sc.B(m[1], m[2])), sc.arrow(sc.B(m[1], m[4]), sc.B(m[4], m[5])), sc.arrow(sc.B(m[5], m[6]), sc.A(m[6]))]);
    expect(closestContact(r.product, 'different fragments')).toBeGreaterThan(1.5);
  });
});
