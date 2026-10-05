import { describe, it, expect } from 'vitest';
import { createDoc, addAtom, addBond, insertMol, docToJSON, docFromJSON } from '../src/doc/document';
import { parseAtomLabel } from '../src/editor/ops';
import { buildScene, chargeText } from '../src/render/scene';
import { expandAbbreviations } from '../src/chem/abbreviations';
import { perceiveStereo2D } from '../src/chem/stereo2d';
import { assignCIP } from '../src/chem/cip';
import { applyArrows, arrowGroups, octetViolations, placeStep, isApplied } from '../src/doc/mechanism';
import { parseSmiles, writeSmiles, suppressHydrogens } from '../src/chem/smiles';
import { ChemDoc, Anchor } from '../src/doc/types';
import { Mol } from '../src/chem/mol';

/** Inserts a SMILES with crude coordinates (atoms on a line) – geometry is irrelevant for electron bookkeeping. */
function put(doc: ChemDoc, smiles: string, x0: number): number[] {
  const m = parseSmiles(smiles);
  m.atoms.forEach((a, i) => { a.x = x0 + i; a.y = (i % 2) * 0.5; });
  return insertMol(doc, m).atomIds;
}

function arrow(doc: ChemDoc, from: Anchor, to: Anchor, electrons: 1 | 2 = 2): number {
  const id = doc.nextId++;
  doc.curved.set(id, { id, type: 'curved', electrons, from, to, c1: { t: 0.3, h: -0.4 }, c2: { t: 0.7, h: -0.4 } });
  return id;
}

function bondId(doc: ChemDoc, a: number, b: number): number {
  for (const bd of doc.bonds.values()) if ((bd.a === a && bd.b === b) || (bd.a === b && bd.b === a)) return bd.id;
  throw new Error('no bond');
}

const smi = (m: Mol) => writeSmiles(m);

describe('arrow pushing', () => {
  it('SN2: hydroxide + bromomethane → methanol + bromide', () => {
    const doc = createDoc();
    const [o] = put(doc, '[OH-]', 0);
    const [c, br] = put(doc, 'CBr', 3);
    const a1 = arrow(doc, { type: 'atom', id: o }, { type: 'atom', id: c });
    const a2 = arrow(doc, { type: 'bond', id: bondId(doc, c, br) }, { type: 'atom', id: br });
    const r = applyArrows(doc, [a1, a2]);
    expect(r.warnings).toEqual([]);
    expect(r.resonance).toBe(false);
    expect(smi(r.product).split('.').sort()).toEqual(['CO', '[Br-]'].sort());
  });

  it('carbonyl addition gives an alkoxide', () => {
    const doc = createDoc();
    const [n] = put(doc, '[C-]#N', 0);
    const ids = put(doc, 'CC=O', 3);
    const a1 = arrow(doc, { type: 'atom', id: n }, { type: 'atom', id: ids[1] });
    const a2 = arrow(doc, { type: 'bond', id: bondId(doc, ids[1], ids[2]) }, { type: 'atom', id: ids[2] });
    const r = applyArrows(doc, [a1, a2]);
    expect(r.warnings).toEqual([]);
    expect(smi(r.product)).toBe(writeSmiles(parseSmiles('CC([O-])C#N')));
  });

  it('detects resonance (enolate)', () => {
    const doc = createDoc();
    const ids = put(doc, '[CH2-]C(C)=O', 0);
    const a1 = arrow(doc, { type: 'atom', id: ids[0] }, { type: 'bond', id: bondId(doc, ids[0], ids[1]) });
    const a2 = arrow(doc, { type: 'bond', id: bondId(doc, ids[1], ids[3]) }, { type: 'atom', id: ids[3] });
    const r = applyArrows(doc, [a1, a2]);
    expect(r.warnings).toEqual([]);
    expect(r.resonance).toBe(true);
    expect(smi(r.product)).toBe(writeSmiles(parseSmiles('C=C(C)[O-]')));
  });

  it('protonation of ammonia', () => {
    const doc = createDoc();
    const [nn] = put(doc, 'N', 0);
    const h = addAtom(doc, { el: 'H', x: 2, y: 0, charge: 1 });
    const a = arrow(doc, { type: 'atom', id: nn }, { type: 'atom', id: h.id });
    const r = applyArrows(doc, [a]);
    expect(r.warnings).toEqual([]);
    const m = r.product;
    expect(m.atoms.find((x) => x.el === 'N')!.charge).toBe(1);
    expect(m.atoms.find((x) => x.el === 'H')!.charge).toBe(0);
  });

  it('homolysis with two fishhooks gives two radicals', () => {
    const doc = createDoc();
    const [c1, c2] = put(doc, 'ClCl', 0);
    const b = bondId(doc, c1, c2);
    const a1 = arrow(doc, { type: 'bond', id: b }, { type: 'atom', id: c1 }, 1);
    const a2 = arrow(doc, { type: 'bond', id: b }, { type: 'atom', id: c2 }, 1);
    const r = applyArrows(doc, [a1, a2]);
    expect(r.warnings).toEqual([]);
    expect(r.product.bonds.length).toBe(0);
    expect(r.product.atoms.every((a) => a.radical === 1 && a.charge === 0)).toBe(true);
  });

  it('warns about octet violations', () => {
    const doc = createDoc();
    const [o] = put(doc, '[OH-]', 0);
    const ids = put(doc, 'CC=O', 3);
    // forget to push the C=O electrons onto oxygen → pentavalent carbon
    const a1 = arrow(doc, { type: 'atom', id: o }, { type: 'atom', id: ids[1] });
    const r = applyArrows(doc, [a1]);
    expect(r.warnings.some((w) => /10 valence electrons/.test(w.message))).toBe(true);
  });

  it('groups arrows by connected reactants', () => {
    const doc = createDoc();
    const [o] = put(doc, '[OH-]', 0);
    const [c, br] = put(doc, 'CBr', 3);
    arrow(doc, { type: 'atom', id: o }, { type: 'atom', id: c });
    arrow(doc, { type: 'bond', id: bondId(doc, c, br) }, { type: 'atom', id: br });
    const [o2] = put(doc, '[OH-]', 20);
    const [x] = put(doc, 'C(=O)C', 23);
    arrow(doc, { type: 'atom', id: o2 }, { type: 'atom', id: x });
    const g = arrowGroups(doc);
    expect(g.length).toBe(2);
    expect(g.map((q) => q.arrows.length).sort()).toEqual([1, 2]);
  });

  it('finds octet violations in drawings', () => {
    const doc = createDoc();
    const c = addAtom(doc, { el: 'C', x: 0, y: 0 });
    for (let k = 0; k < 5; k++) {
      const a = addAtom(doc, { el: 'C', x: Math.cos(k), y: Math.sin(k) });
      addBond(doc, c.id, a.id);
    }
    expect(octetViolations(doc).map((v) => v.atomId)).toContain(c.id);
  });

  it('marks impossible products as errors, not as valid steps', () => {
    const doc = createDoc();
    const [o] = put(doc, '[OH-]', 0);
    const ids = put(doc, 'CC=O', 3);
    const r = applyArrows(doc, [arrow(doc, { type: 'atom', id: o }, { type: 'atom', id: ids[1] })]);
    expect(r.ok).toBe(false);
    expect(r.changed).toBe(true);
    expect(r.resonance).toBe(false);
    expect(r.warnings.find((w) => /10 valence electrons/.test(w.message))?.level).toBe('error');
  });

  it('never calls a failed or empty step resonance', () => {
    // arrow into empty space: nothing changes
    const doc = createDoc();
    const [o] = put(doc, '[OH-]', 0);
    const r1 = applyArrows(doc, [arrow(doc, { type: 'atom', id: o }, { type: 'point', x: 5, y: 5 })]);
    expect(r1.ok).toBe(false);
    expect(r1.changed).toBe(false);
    expect(r1.resonance).toBe(false);
    // cyanide lone pair pushed into its own C≡N: bond order 4, σ framework unchanged
    const doc2 = createDoc();
    const [c, n] = put(doc2, '[C-]#N', 0);
    const r2 = applyArrows(doc2, [arrow(doc2, { type: 'atom', id: c }, { type: 'bond', id: bondId(doc2, c, n) })]);
    expect(r2.ok).toBe(false);
    expect(r2.resonance).toBe(false);
    // a lone pair pushed back onto its own atom changes nothing
    const doc3 = createDoc();
    const [o3] = put(doc3, '[OH-]', 0);
    const r3 = applyArrows(doc3, [arrow(doc3, { type: 'atom', id: o3 }, { type: 'atom', id: o3 })]);
    expect(r3.changed).toBe(false);
    expect(r3.resonance).toBe(false);
  });

  it('reports arrows whose anchors no longer exist', () => {
    const doc = createDoc();
    const [o] = put(doc, '[OH-]', 0);
    const [c] = put(doc, 'CBr', 3);
    const a = arrow(doc, { type: 'atom', id: o }, { type: 'atom', id: c });
    doc.atoms.delete(c);
    const r = applyArrows(doc, [a]);
    expect(r.ok).toBe(false);
    expect(r.warnings.some((w) => /no longer exists/.test(w.message))).toBe(true);
  });

  it('knows which arrow groups were already applied', () => {
    const doc = createDoc();
    const [o] = put(doc, '[OH-]', 0);
    const [c, br] = put(doc, 'CBr', 3);
    const a1 = arrow(doc, { type: 'atom', id: o }, { type: 'atom', id: c });
    const a2 = arrow(doc, { type: 'bond', id: bondId(doc, c, br) }, { type: 'atom', id: br });
    expect(arrowGroups(doc)[0].applied).toBe(false);
    const step = doc.nextId++;
    doc.arrows.set(step, { id: step, type: 'arrow', kind: 'reaction', x1: 6, y1: 0, x2: 8, y2: 0 });
    doc.curved.get(a1)!.step = step;
    expect(arrowGroups(doc)[0].applied).toBe(false); // only one of the two arrows
    doc.curved.get(a2)!.step = step;
    expect(arrowGroups(doc)[0].applied).toBe(true);
    doc.arrows.delete(step); // deleting the step's arrow makes it applicable again
    expect(arrowGroups(doc)[0].applied).toBe(false);
  });

  it('places the product to the right of the reactants', () => {
    const doc = createDoc();
    const [o] = put(doc, '[OH-]', 0);
    const [c, br] = put(doc, 'CBr', 3);
    const r = applyArrows(doc, [
      arrow(doc, { type: 'atom', id: o }, { type: 'atom', id: c }),
      arrow(doc, { type: 'bond', id: bondId(doc, c, br) }, { type: 'atom', id: br }),
    ]);
    const p = placeStep(doc, r)!;
    expect(p.arrow.kind).toBe('reaction');
    expect(p.arrow.x1).toBeGreaterThan(4);
    expect(Math.min(...p.mol.atoms.map((a) => a.x))).toBeGreaterThan(p.arrow.x2);
    expect(r.product.atoms.map((a) => a.x)).toEqual([0, 3, 4]); // the engine's product itself is not moved
  });
});

/** Atom ids by label for hand-placed drawings: spec maps label → [element, x, y, charge?]. */
function draw(doc: ChemDoc, spec: Record<string, [string, number, number, number?]>, bonds: [string, string, number?][]): Record<string, number> {
  const ids: Record<string, number> = {};
  for (const [k, [el, x, y, charge]] of Object.entries(spec)) ids[k] = addAtom(doc, { el, x, y, charge: charge ?? 0 }).id;
  for (const [a, b, o] of bonds) addBond(doc, ids[a], ids[b], o ?? 1);
  return ids;
}
const atomA = (id: number): Anchor => ({ type: 'atom', id });
const canon = (s: string) => writeSmiles(parseSmiles(s));
const productOf = (r: { product: Mol }) => writeSmiles(suppressHydrogens(r.product));

describe('ambiguous arrows are read chemically, not by distance', () => {
  it('1,2-hydride shift: the H migrates even though C3 is drawn nearer', () => {
    const doc = createDoc();
    const L = draw(doc, { C1: ['C', 0, 0], C2: ['C', 1, 0, 1], C3: ['C', 2, 0], C4: ['C', 3, 0], C5: ['C', 2, 1], H: ['H', 2, -1] },
      [['C1', 'C2'], ['C2', 'C3'], ['C3', 'C4'], ['C3', 'C5'], ['C3', 'H']]);
    const r = applyArrows(doc, [arrow(doc, { type: 'bond', id: bondId(doc, L.C3, L.H) }, atomA(L.C2))]);
    expect(r.ok).toBe(true);
    expect(productOf(r)).toBe(canon('CC[C+](C)C'));
  });

  it('1,2-methyl shift moves the methyl group', () => {
    const doc = createDoc();
    const L = draw(doc, { Cp: ['C', 0, 0, 1], Cq: ['C', 1, 0], M1: ['C', 2, 0], M2: ['C', 1, 1], M3: ['C', 1, -1] },
      [['Cp', 'Cq'], ['Cq', 'M1'], ['Cq', 'M2'], ['Cq', 'M3']]);
    const r = applyArrows(doc, [arrow(doc, { type: 'bond', id: bondId(doc, L.Cq, L.M3) }, atomA(L.Cp))]);
    expect(productOf(r)).toBe(canon('CC[C+](C)C'));
  });

  it('HBr adds Markovnikov wherever the H is drawn', () => {
    for (const hx of [0, 1, 2.2]) {
      const doc = createDoc();
      const L = draw(doc, { C1: ['C', 0, 0], C2: ['C', 1, 0], C3: ['C', 2, 0.5], H: ['H', hx, -1.3], Br: ['Br', hx, -2.3] },
        [['C1', 'C2', 2], ['C2', 'C3'], ['H', 'Br']]);
      const r = applyArrows(doc, [
        arrow(doc, { type: 'bond', id: bondId(doc, L.C1, L.C2) }, atomA(L.H)),
        arrow(doc, { type: 'bond', id: bondId(doc, L.H, L.Br) }, atomA(L.Br)),
      ]);
      expect(r.ok).toBe(true);
      expect(productOf(r).split('.').sort()).toEqual([canon('C[CH+]C'), '[Br-]'].sort());
    }
  });

  it('a Grignard attacks through carbon even with Mg drawn nearer the carbonyl', () => {
    const doc = createDoc();
    const L = draw(doc, { Me: ['C', 0, 0], Mg: ['Mg', 1, 0], Br: ['Br', 2, 0], C: ['C', 1.6, 1.2], O: ['O', 2.5, 1.7], A: ['C', 1, 2], B: ['C', 2, 0.9] },
      [['Me', 'Mg'], ['Mg', 'Br'], ['C', 'O', 2], ['C', 'A'], ['C', 'B']]);
    const r = applyArrows(doc, [
      arrow(doc, { type: 'bond', id: bondId(doc, L.Me, L.Mg) }, atomA(L.C)),
      arrow(doc, { type: 'bond', id: bondId(doc, L.C, L.O) }, atomA(L.O)),
    ]);
    expect(r.ok).toBe(true);
    expect(productOf(r).split('.').sort()).toEqual([canon('CC(C)(C)[O-]'), canon('[Mg+]Br')].sort());
  });

  it('a lone pair aimed at a C=O bond attacks the carbon', () => {
    const doc = createDoc();
    const L = draw(doc, { O1: ['O', 0, 0, -1], C: ['C', 2, 1], O: ['O', 1, 0.5], Me: ['C', 3, 1] }, [['C', 'O', 2], ['C', 'Me']]);
    const r = applyArrows(doc, [
      arrow(doc, atomA(L.O1), { type: 'bond', id: bondId(doc, L.C, L.O) }),
      arrow(doc, { type: 'bond', id: bondId(doc, L.C, L.O) }, atomA(L.O)),
    ]);
    expect(r.ok).toBe(true);
    expect(productOf(r)).toBe(canon('CC([O-])O'));
  });

  it('flags an arrow from a lone pair that points between two other atoms', () => {
    const doc = createDoc();
    const L = draw(doc, { N: ['C', 0, 0, -1], C: ['C', 2, 0], O: ['O', 3, 0], M: ['C', 2, 1] }, [['C', 'O', 2], ['C', 'M']]);
    const r = applyArrows(doc, [arrow(doc, atomA(L.N), { type: 'between', a: L.M, b: L.O })]);
    expect(r.ok).toBe(false);
    expect(r.warnings[0].message).toMatch(/between two other atoms/);
  });

  it('does not warn when both readings give the same product (benzene)', () => {
    const doc = createDoc();
    const ids = put(doc, 'C1=CC=CC=C1', 0);
    const br = put(doc, 'BrBr', 0).map((id, k) => (doc.atoms.get(id)!.y = -2 - k, id));
    const r = applyArrows(doc, [
      arrow(doc, { type: 'bond', id: bondId(doc, ids[0], ids[1]) }, atomA(br[0])),
      arrow(doc, { type: 'bond', id: bondId(doc, br[0], br[1]) }, atomA(br[1])),
    ]);
    expect(r.ok).toBe(true);
    expect(r.warnings).toEqual([]);
  });

  it('warns when only the drawing decides between two different products', () => {
    // pent-2-ene + H⁺: both ends give a secondary cation, but different ones
    const doc = createDoc();
    const L = draw(doc, { C1: ['C', 0, 0], C2: ['C', 1, 0], C3: ['C', 2, 0], C4: ['C', 3, 0], C5: ['C', 4, 0], H: ['H', 1.2, -1.5, 1] },
      [['C1', 'C2'], ['C2', 'C3', 2], ['C3', 'C4'], ['C4', 'C5']]);
    const r = applyArrows(doc, [arrow(doc, { type: 'bond', id: bondId(doc, L.C2, L.C3) }, atomA(L.H))]);
    expect(r.ok).toBe(true);
    expect(productOf(r)).toBe(canon('CC[CH+]CC')); // H drawn nearer C2 → pentan-3-yl cation
    expect(r.warnings.some((w) => w.level === 'warning' && /nearer atom was used/.test(w.message))).toBe(true);
  });
});

describe('rings drawn with delocalised (aromatic) bonds', () => {
  /** Benzene drawn with the aromatic bond tool (order 1.5), atoms on a hexagon starting at (cx, cy). */
  function aromaticBenzene(doc: ChemDoc, cx = 0, cy = 0): number[] {
    const ids = [...Array(6)].map((_, k) => addAtom(doc, { el: 'C', x: cx + Math.cos((k * Math.PI) / 3), y: cy + Math.sin((k * Math.PI) / 3) }).id);
    ids.forEach((id, k) => addBond(doc, id, ids[(k + 1) % 6], 1.5));
    return ids;
  }

  it('electrophilic attack on benzene gives the arenium ion', () => {
    const doc = createDoc();
    const ring = aromaticBenzene(doc);
    const br1 = addAtom(doc, { el: 'Br', x: 2.2, y: 0 }).id, br2 = addAtom(doc, { el: 'Br', x: 3.2, y: 0 }).id;
    addBond(doc, br1, br2);
    const r = applyArrows(doc, [
      arrow(doc, { type: 'bond', id: bondId(doc, ring[0], ring[1]) }, atomA(br1)),
      arrow(doc, { type: 'bond', id: bondId(doc, br1, br2) }, atomA(br2)),
    ]);
    expect(r.warnings).toEqual([]);
    expect(productOf(r).split('.').sort()).toEqual([canon('BrC1C=CC=C[CH+]1'), '[Br-]'].sort());
  });

  it('a spectator ring stays delocalised', () => {
    const doc = createDoc();
    const ring = aromaticBenzene(doc);
    const ch2 = addAtom(doc, { el: 'C', x: 2, y: 0 }).id, br = addAtom(doc, { el: 'Br', x: 3, y: 0 }).id;
    addBond(doc, ring[0], ch2);
    addBond(doc, ch2, br);
    const o = addAtom(doc, { el: 'O', x: 2, y: -1.5, charge: -1 }).id;
    const r = applyArrows(doc, [arrow(doc, atomA(o), atomA(ch2)), arrow(doc, { type: 'bond', id: bondId(doc, ch2, br) }, atomA(br))]);
    expect(r.warnings).toEqual([]);
    expect(productOf(r).split('.').sort()).toEqual([canon('OCc1ccccc1'), '[Br-]'].sort());
    const ringBonds = r.product.bonds.filter((b) => ring.includes(r.product.atoms[b.a].id) && ring.includes(r.product.atoms[b.b].id));
    expect(ringBonds.map((b) => b.order)).toEqual([1.5, 1.5, 1.5, 1.5, 1.5, 1.5]);
    expect(r.product.atoms.every((a) => !a.charge || a.el === 'Br')).toBe(true);
  });

  it('three resonance arrows around a delocalised ring are a valid resonance step', () => {
    const doc = createDoc();
    const ring = aromaticBenzene(doc);
    const ids = [0, 2, 4].map((k) => arrow(doc, { type: 'bond', id: bondId(doc, ring[k], ring[k + 1]) }, { type: 'bond', id: bondId(doc, ring[k + 1], ring[(k + 2) % 6]) }));
    const r = applyArrows(doc, ids);
    expect(r.ok).toBe(true);
    expect(r.resonance).toBe(true);
    expect(productOf(r)).toBe(canon('c1ccccc1'));
  });

  it('explains when a ring cannot be given alternating double bonds', () => {
    const doc = createDoc();
    const ids = [...Array(5)].map((_, k) => addAtom(doc, { el: 'C', x: Math.cos((k * 2 * Math.PI) / 5), y: Math.sin((k * 2 * Math.PI) / 5) }).id);
    ids.forEach((id, k) => addBond(doc, id, ids[(k + 1) % 5], 1.5));
    const o = addAtom(doc, { el: 'O', x: 3, y: 0, charge: -1 }).id;
    const h = addAtom(doc, { el: 'H', x: 2, y: 0, charge: 1 }).id;
    const r = applyArrows(doc, [arrow(doc, atomA(o), atomA(h)), arrow(doc, atomA(ids[0]), atomA(ids[0]))]);
    expect(r.warnings.some((w) => w.level === 'error' && /Kekulé/.test(w.message))).toBe(true);
  });
});

describe('labels, generic atoms and charge conservation', () => {
  const labelled = (doc: ChemDoc, label: string, x: number, y: number) => {
    const a = addAtom(doc, { el: 'C', x, y });
    Object.assign(a, parseAtomLabel(label));
    return a.id;
  };
  const expanded = (r: { product: Mol }) => writeSmiles(suppressHydrogens(expandAbbreviations(r.product)));

  it('a methoxide drawn as an OMe⁻ label does an SN2 and stays a label', () => {
    const doc = createDoc();
    const ome = labelled(doc, 'MeO-', 0, 0);
    const [c, br] = put(doc, 'CBr', 2);
    const r = applyArrows(doc, [arrow(doc, atomA(ome), atomA(c)), arrow(doc, { type: 'bond', id: bondId(doc, c, br) }, atomA(br))]);
    expect(r.warnings).toEqual([]);
    const lab = r.product.atoms.find((a) => a.id === ome)!;
    expect(lab).toMatchObject({ abbrev: 'OMe', charge: 0 });
    expect(expanded(r).split('.').sort()).toEqual([canon('COC'), '[Br-]'].sort());
  });

  it('an OMe leaving group leaves as methoxide', () => {
    const doc = createDoc();
    const L = draw(doc, { Me: ['C', 0, 0], C: ['C', 1, 0], O: ['O', 1, -1, -1], OH: ['O', 1, 1] }, [['Me', 'C'], ['C', 'O'], ['C', 'OH']]);
    const ome = labelled(doc, 'OMe', 2, 0);
    addBond(doc, L.C, ome);
    const r = applyArrows(doc, [arrow(doc, atomA(L.O), { type: 'bond', id: bondId(doc, L.C, L.O) }), arrow(doc, { type: 'bond', id: bondId(doc, L.C, ome) }, atomA(ome))]);
    expect(r.warnings).toEqual([]);
    expect(r.product.atoms.find((a) => a.id === ome)).toMatchObject({ abbrev: 'OMe', charge: -1 });
    expect(expanded(r).split('.').sort()).toEqual([canon('CC(=O)O'), canon('C[O-]')].sort());
  });

  it('a nucleophile attacking a label pushes the label’s own C=O onto its oxygen', () => {
    const doc = createDoc();
    const [me] = put(doc, 'C', 0);
    const co2h = labelled(doc, 'CO2H', 1, 0);
    addBond(doc, me, co2h);
    const o = addAtom(doc, { el: 'O', x: 1, y: -2, charge: -1 }).id;
    const r = applyArrows(doc, [arrow(doc, atomA(o), atomA(co2h))]);
    expect(r.ok).toBe(true);
    expect(expanded(r)).toBe(canon('CC(O)(O)[O-]'));
  });

  it('points out that arrows cannot reach inside a label', () => {
    const doc = createDoc();
    const [me] = put(doc, 'C', 0);
    const lab = labelled(doc, 'CH2OH', 1, 0);
    addBond(doc, me, lab);
    const o = addAtom(doc, { el: 'O', x: 1, y: -2, charge: -1 }).id;
    const r = applyArrows(doc, [arrow(doc, atomA(o), atomA(lab))]);
    expect(r.ok).toBe(false);
    expect(r.warnings.some((w) => /expand it/.test(w.message))).toBe(true);
  });

  it('typed reagents take part as real atoms: H2O and CN⁻', () => {
    const doc = createDoc();
    const w = labelled(doc, 'H2O', 0, 0);
    const [m1, cp, m2, m3] = put(doc, 'C[C+](C)C', 2);
    const r = applyArrows(doc, [arrow(doc, atomA(w), atomA(cp))]);
    expect(r.warnings).toEqual([]);
    expect(productOf(r)).toBe(canon('CC(C)(C)[OH2+]'));
    void [m1, m2, m3];

    const doc2 = createDoc();
    const cn = labelled(doc2, 'CN-', 0, 0);
    const ids = put(doc2, 'CC(C)=O', 2);
    const r2 = applyArrows(doc2, [arrow(doc2, atomA(cn), atomA(ids[1])), arrow(doc2, { type: 'bond', id: bondId(doc2, ids[1], ids[3]) }, atomA(ids[3]))]);
    expect(r2.warnings).toEqual([]);
    expect(expanded(r2)).toBe(canon('CC(C)([O-])C#N'));
  });

  it('a generic E⁺ electrophile bonds and becomes neutral', () => {
    const doc = createDoc();
    const ids = put(doc, 'C1=CC=CC=C1', 0);
    const e = labelled(doc, 'E+', 0, -2);
    const r = applyArrows(doc, [arrow(doc, { type: 'bond', id: bondId(doc, ids[0], ids[1]) }, atomA(e))]);
    expect(r.warnings).toEqual([]);
    expect(r.product.atoms.find((a) => a.id === e)!.charge).toBe(0);
    expect(r.product.atoms.reduce((s, a) => s + a.charge, 0)).toBe(1);
  });

  it('breaking a dative bond moves no electrons; hydrogen bonds are kept', () => {
    const doc = createDoc();
    const L = draw(doc, { P: ['P', 0, 0], A: ['C', -1, 0], B: ['C', 0, 1], C: ['C', 0, -1], Pd: ['Pd', 1.5, 0] }, [['P', 'A'], ['P', 'B'], ['P', 'C']]);
    addBond(doc, L.P, L.Pd, 1, 'dative');
    const r = applyArrows(doc, [arrow(doc, { type: 'bond', id: bondId(doc, L.P, L.Pd) }, atomA(L.P))]);
    expect(r.warnings).toEqual([]);
    expect(r.product.atoms.every((a) => !a.charge)).toBe(true);
    expect(r.product.bonds.some((b) => b.style === 'dative')).toBe(false);

    const doc2 = createDoc();
    const M = draw(doc2, { O: ['O', 0, 0], H: ['H', 1, 0, 1], W: ['O', 3, 0], C: ['C', 4, 0], O2: ['O', 5, 0] }, [['C', 'O2', 2]]);
    addBond(doc2, M.H, M.W, 0, 'hbond');
    const r2 = applyArrows(doc2, [arrow(doc2, atomA(M.O), atomA(M.H))]);
    expect(r2.product.bonds.some((b) => b.style === 'hbond')).toBe(true);
  });

  it('flags steps that would change the total charge', () => {
    const doc = createDoc();
    const [c, n] = put(doc, '[C-]#N', 0);
    const r = applyArrows(doc, [arrow(doc, atomA(c), { type: 'bond', id: bondId(doc, c, n) })]);
    expect(r.warnings.some((w) => w.level === 'error' && /total charge/.test(w.message))).toBe(true);
  });

  it('draws the charge of labels and generic atoms', () => {
    const doc = createDoc();
    labelled(doc, 'MeO-', 0, 0);
    labelled(doc, 'E+', 3, 0);
    const texts = buildScene(doc, { ink: '#000' }).prims.filter((q) => q.k === 'text').map((q) => (q as { text: string }).text);
    expect(texts).toContain(chargeText(-1));
    expect(texts).toContain(chargeText(1));
  });
});

describe('stereochemistry through a step', () => {
  const cipOf = (m: Mol, id: number) => {
    const x = m.clone();
    perceiveStereo2D(x);
    const i = x.atoms.findIndex((a) => a.id === id);
    return assignCIP(x).centers.get(i);
  };
  /** (R)-2-bromobutane with the Br on a wedge, plus hydroxide at (ox, oy). */
  function sn2(ox: number, oy: number) {
    const doc = createDoc();
    const L = draw(doc, { C2: ['C', 0, 0], C1: ['C', -0.866, 0.5], C3: ['C', 0.866, 0.5], C4: ['C', 1.732, 0], Br: ['Br', 0, -1], O: ['O', ox, oy, -1] },
      [['C2', 'C1'], ['C2', 'C3'], ['C3', 'C4']]);
    addBond(doc, L.C2, L.Br, 1, 'wedge');
    const r = applyArrows(doc, [arrow(doc, atomA(L.O), atomA(L.C2)), arrow(doc, { type: 'bond', id: bondId(doc, L.C2, L.Br) }, atomA(L.Br))]);
    return { L, r };
  }

  it('backside substitution inverts the centre however the nucleophile is drawn', () => {
    for (const [x, y] of [[0, 2.2], [0, -2.2], [-2, -1]]) {
      const { L, r } = sn2(x, y);
      expect(r.warnings).toEqual([]);
      expect(cipOf(r.product, L.C2)).toBe('S');
    }
  });

  it('keeps the colour and double-bond position of untouched bonds', () => {
    const doc = createDoc();
    // allyl alkoxide (the C=C is in the reacting molecule) + H⁺
    const L = draw(doc, { A: ['C', 0, 0], B: ['C', 1, 0], C: ['C', 2, 0.5], O: ['O', 3, 0, -1], H: ['H', 4, 0, 1] }, [['B', 'C'], ['C', 'O']]);
    const db = addBond(doc, L.A, L.B, 2);
    db.dbPos = 'left';
    db.color = '#e03131';
    const r = applyArrows(doc, [arrow(doc, atomA(L.O), atomA(L.H))]);
    const kept = r.product.bonds.find((b) => b.id === db.id)!;
    expect(kept).toMatchObject({ order: 2, dbPos: 'left', color: '#e03131' });
  });
});

describe('review fixes: bonds, stereo, labels and generic atoms', () => {
  const labelled = (doc: ChemDoc, label: string, x: number, y: number) => {
    const a = addAtom(doc, { el: 'C', x, y });
    Object.assign(a, parseAtomLabel(label));
    return a.id;
  };
  const expanded = (r: { product: Mol }) => writeSmiles(suppressHydrogens(expandAbbreviations(r.product)));
  const H = (id: number): Anchor => ({ type: 'atom', id, h: true });

  it('a proton moved along a drawn hydrogen bond leaves one covalent bond', () => {
    const doc = createDoc();
    const L = draw(doc, { O1: ['O', 0, 0, -1], H: ['H', 1.2, 0], O2: ['O', 2.2, 0] }, [['H', 'O2']]);
    addBond(doc, L.O1, L.H, 0, 'hbond');
    const r = applyArrows(doc, [arrow(doc, atomA(L.O1), { type: 'bond', id: bondId(doc, L.O1, L.H) }), arrow(doc, { type: 'bond', id: bondId(doc, L.H, L.O2) }, atomA(L.O2))]);
    expect(r.ok).toBe(true);
    expect(r.product.bonds.filter((b) => b.style === 'hbond')).toEqual([]);
    expect(productOf(r).split('.').sort()).toEqual(['O', '[OH-]'].sort());
  });

  it('breaking only a dative bond is a reaction step, not resonance', () => {
    const doc = createDoc();
    const L = draw(doc, { P: ['P', 0, 0], A: ['C', -1, 0], B: ['C', 0, 1], C: ['C', 0, -1], Pd: ['Pd', 1.5, 0] }, [['P', 'A'], ['P', 'B'], ['P', 'C']]);
    addBond(doc, L.P, L.Pd, 1, 'dative');
    const r = applyArrows(doc, [arrow(doc, { type: 'bond', id: bondId(doc, L.P, L.Pd) }, atomA(L.P))]);
    expect(r.ok).toBe(true);
    expect(r.resonance).toBe(false);
  });

  it('a stereocentre that becomes a planar cation loses its wedge', () => {
    const doc = createDoc();
    const L = draw(doc, { C2: ['C', 0, 0], C3: ['C', 0.866, 0.5], C4: ['C', 1.732, 0], Br: ['Br', 0, -1], C1: ['C', -0.866, 0.5] }, [['C2', 'C3'], ['C3', 'C4'], ['C2', 'Br']]);
    addBond(doc, L.C2, L.C1, 1, 'wedge');
    const r = applyArrows(doc, [arrow(doc, { type: 'bond', id: bondId(doc, L.C2, L.Br) }, atomA(L.Br))]);
    expect(r.product.bonds.some((b) => b.style === 'wedge' || b.style === 'hash')).toBe(false);
  });

  it('a condensed OCH3⁻ label protonated through a label H becomes methanol', () => {
    const doc = createDoc();
    const och3 = labelled(doc, 'OCH3-', 0, 0);
    const h3o = labelled(doc, 'H3O+', 3, 0);
    const r = applyArrows(doc, [arrow(doc, atomA(och3), H(h3o)), arrow(doc, H(h3o), atomA(h3o))]);
    expect(r.ok).toBe(true);
    expect(expanded(r).split('.').sort()).toEqual(['CO', 'O'].sort());
  });

  it('an OMe label whose oxygen loses its bond electrons is drawn in full, without invented hydrogens', () => {
    const doc = createDoc();
    const [c] = put(doc, 'C', 0);
    const ome = labelled(doc, 'OMe', 1, 0);
    addBond(doc, c, ome);
    const r = applyArrows(doc, [arrow(doc, { type: 'bond', id: bondId(doc, c, ome) }, atomA(c))]);
    expect(expanded(r).split('.').sort()).toEqual(['C[O+]', '[CH3-]'].sort());
  });

  it('a free PPh3 label is triphenylphosphine and alkylates as a nucleophile', () => {
    const doc = createDoc();
    const p = labelled(doc, 'PPh3', 0, 0);
    const [c, i] = put(doc, 'CI', 2);
    const r = applyArrows(doc, [arrow(doc, atomA(p), atomA(c)), arrow(doc, { type: 'bond', id: bondId(doc, c, i) }, atomA(i))]);
    expect(r.ok).toBe(true);
    expect(expanded(r).split('.').sort()).toEqual([canon('C[P+](c1ccccc1)(c1ccccc1)c1ccccc1'), '[I-]'].sort());
    expect(r.product.atoms.find((a) => a.id === p)).toMatchObject({ abbrev: 'PPh3', charge: 0 }); // shown as PPh3⁺ once bonded
  });

  it('a carboxylate written as COO⁻ counts its charge once', () => {
    const doc = createDoc();
    const L = draw(doc, { N: ['N', 0, 0], C: ['C', 1, 0], Hp: ['H', -1.5, 0, 1] }, [['N', 'C']]);
    const coo = labelled(doc, 'COO-', 2, 0);
    addBond(doc, L.C, coo);
    const r = applyArrows(doc, [arrow(doc, atomA(L.N), atomA(L.Hp))]);
    expect(r.ok).toBe(true);
    expect(expanded(r)).toBe(canon('[NH3+]CC(=O)[O-]'));
  });

  it('neutral generic nucleophiles and bases have a lone pair to give', () => {
    const doc = createDoc();
    const nu = labelled(doc, 'Nu', 0, 0);
    const [c, br] = put(doc, 'CBr', 2);
    const r = applyArrows(doc, [arrow(doc, atomA(nu), atomA(c)), arrow(doc, { type: 'bond', id: bondId(doc, c, br) }, atomA(br))]);
    expect(r.ok).toBe(true);
    expect(r.product.atoms.find((a) => a.id === nu)!.charge).toBe(1);

    const doc2 = createDoc();
    const b = labelled(doc2, 'B:', 0, 0);
    const L = draw(doc2, { Ca: ['C', 2, 0], H: ['H', 1.2, -0.6], C: ['C', 3, 0], O: ['O', 3.5, -0.87], M: ['C', 3.5, 0.87] }, [['Ca', 'H'], ['Ca', 'C'], ['C', 'O', 2], ['C', 'M']]);
    const r2 = applyArrows(doc2, [
      arrow(doc2, atomA(b), atomA(L.H)),
      arrow(doc2, { type: 'bond', id: bondId(doc2, L.Ca, L.H) }, { type: 'bond', id: bondId(doc2, L.Ca, L.C) }),
      arrow(doc2, { type: 'bond', id: bondId(doc2, L.C, L.O) }, atomA(L.O)),
    ]);
    expect(r2.ok).toBe(true);
  });

  it('nitration with an NO2⁺ label gives the nitro arenium and a nitro label', () => {
    const doc = createDoc();
    const ids = put(doc, 'C1=CC=CC=C1', 0);
    const no2 = labelled(doc, 'NO2+', 0, -2);
    const r = applyArrows(doc, [arrow(doc, { type: 'bond', id: bondId(doc, ids[0], ids[1]) }, atomA(no2))]);
    expect(r.ok).toBe(true);
    expect(expanded(r)).toBe(canon('[O-][N+](=O)C1C=CC=C[CH+]1'));
    expect(r.product.atoms.find((a) => a.id === no2)).toMatchObject({ abbrev: 'NO2' });
  });

  it('N-nitrosation with an NO⁺ label gives the N-nitrosammonium', () => {
    const doc = createDoc();
    const [m1, n] = put(doc, 'CNC', 0);
    const no = labelled(doc, 'NO+', 1, -2);
    const r = applyArrows(doc, [arrow(doc, atomA(n), atomA(no))]);
    expect(r.ok).toBe(true);
    expect(expanded(r)).toBe(canon('C[NH+](C)N=O'));
    void m1;
  });

  it('nitrite (NO2⁻ label) alkylated on nitrogen gives a nitro group', () => {
    const doc = createDoc();
    const no2 = labelled(doc, 'NO2-', 0, 0);
    const [c, br] = put(doc, 'CBr', 2);
    const r = applyArrows(doc, [arrow(doc, atomA(no2), atomA(c)), arrow(doc, { type: 'bond', id: bondId(doc, c, br) }, atomA(br))]);
    expect(r.ok).toBe(true);
    expect(expanded(r).split('.').sort()).toEqual([canon('C[N+](=O)[O-]'), '[Br-]'].sort());
  });

  it('thia-Michael with a typed MeS⁻ (methanethiolate, not mesityl)', () => {
    const doc = createDoc();
    const ids = put(doc, 'C=CC(=O)OC', 0);
    const s1 = labelled(doc, 'MeS-', 0, -2);
    const r = applyArrows(doc, [
      arrow(doc, atomA(s1), atomA(ids[0])),
      arrow(doc, { type: 'bond', id: bondId(doc, ids[0], ids[1]) }, { type: 'bond', id: bondId(doc, ids[1], ids[2]) }),
      arrow(doc, { type: 'bond', id: bondId(doc, ids[2], ids[3]) }, atomA(ids[3])),
    ]);
    expect(r.ok).toBe(true);
    expect(expanded(r)).toBe(canon('CSCC=C([O-])OC'));
  });

  it('charged rings drawn delocalised: cyclopentadienyl anion and tropylium', () => {
    const ring = (doc: ChemDoc, n: number, charge: number) => {
      const ids = [...Array(n)].map((_, k) => addAtom(doc, { el: 'C', x: Math.cos((k * 2 * Math.PI) / n), y: Math.sin((k * 2 * Math.PI) / n), charge: k === 0 ? charge : 0 }).id);
      ids.forEach((id, k) => addBond(doc, id, ids[(k + 1) % n], 1.5));
      return ids;
    };
    const doc = createDoc();
    const cp = ring(doc, 5, -1);
    const h = addAtom(doc, { el: 'H', x: 3, y: 0, charge: 1 }).id;
    const r = applyArrows(doc, [arrow(doc, atomA(cp[0]), atomA(h))]);
    expect(r.warnings).toEqual([]);
    expect(productOf(r)).toBe(canon('C1C=CC=C1'));

    const doc2 = createDoc();
    const tr = ring(doc2, 7, 1);
    const o = addAtom(doc2, { el: 'O', x: 3, y: 0, charge: -1 }).id;
    const r2 = applyArrows(doc2, [arrow(doc2, atomA(o), atomA(tr[0]))]);
    expect(r2.warnings).toEqual([]);
    expect(productOf(r2)).toBe(canon('OC1C=CC=CC=C1'));
  });

  it('a ring that cannot be given alternating bonds gives one error, not one per bond', () => {
    const doc = createDoc();
    const ids = [...Array(5)].map((_, k) => addAtom(doc, { el: 'C', x: Math.cos((k * 2 * Math.PI) / 5), y: Math.sin((k * 2 * Math.PI) / 5) }).id);
    ids.forEach((id, k) => addBond(doc, id, ids[(k + 1) % 5], 1.5));
    const o = addAtom(doc, { el: 'O', x: 3, y: 0, charge: -1 }).id;
    const hp = addAtom(doc, { el: 'H', x: 2, y: 0, charge: 1 }).id;
    const r = applyArrows(doc, [arrow(doc, atomA(o), atomA(hp)), arrow(doc, atomA(ids[0]), atomA(ids[0]))]);
    expect(r.warnings.filter((w) => /Odd number/.test(w.message))).toEqual([]);
    expect(r.warnings.filter((w) => /Kekulé/.test(w.message)).length).toBe(1);
  });
});

describe('review fixes: layout', () => {
  const closestNonBonded = (m: Mol, els: string[]) => {
    let best = Infinity;
    for (let i = 0; i < m.atoms.length; i++)
      for (let j = i + 1; j < m.atoms.length; j++) {
        if (!els.includes(m.atoms[i].el) || !els.includes(m.atoms[j].el) || m.bondBetween(i, j) >= 0) continue;
        best = Math.min(best, Math.hypot(m.atoms[i].x - m.atoms[j].x, m.atoms[i].y - m.atoms[j].y));
      }
    return best;
  };

  it('opening a protonated epoxide keeps the two oxygens apart', () => {
    const doc = createDoc();
    const L = draw(doc, { C1: ['C', 0, 0], C2: ['C', 1, 0], O: ['O', 0.5, -0.866, 1], Me: ['C', 1.87, 0.5], W: ['O', 2.2, -0.9] }, [['C1', 'C2'], ['C2', 'O'], ['C1', 'O'], ['C2', 'Me']]);
    doc.atoms.get(L.O)!.hCount = 1;
    const r = applyArrows(doc, [arrow(doc, atomA(L.W), atomA(L.C2)), arrow(doc, { type: 'bond', id: bondId(doc, L.C2, L.O) }, atomA(L.O))]);
    expect(r.ok).toBe(true);
    const p = placeStep(doc, r)!;
    expect(closestNonBonded(p.mol, ['O'])).toBeGreaterThan(1.3);
  });

  it('the step arrow starts after a wide label on the right of the reactants', () => {
    const doc = createDoc();
    const a = addAtom(doc, { el: 'C', x: -3, y: 0 });
    const b = addAtom(doc, { el: 'C', x: -2, y: 0 });
    Object.assign(a, parseAtomLabel('OtBu'));
    Object.assign(b, parseAtomLabel('OtBu'));
    addBond(doc, a.id, b.id);
    const bid = bondId(doc, a.id, b.id);
    const r = applyArrows(doc, [arrow(doc, { type: 'bond', id: bid }, atomA(a.id), 1), arrow(doc, { type: 'bond', id: bid }, atomA(b.id), 1)]);
    expect(r.ok).toBe(true);
    const p = placeStep(doc, r)!;
    expect(p.arrow.x1).toBeGreaterThan(-2 + 1.2);
  });
});


describe('review fixes: applied steps', () => {
  it('a step whose reaction arrow was deleted is forgotten on reload, so its id cannot mark it applied again', () => {
    const doc = createDoc();
    const [o, c] = put(doc, '[OH-].CBr', 0);
    arrow(doc, { type: 'atom', id: o }, { type: 'between', a: o, b: c + 1 });
    const step = doc.nextId + 50; // the deleted step's arrow had the highest id in the drawing
    for (const cv of doc.curved.values()) cv.step = step;
    expect(isApplied(doc, [...doc.curved.values()][0])).toBe(false);
    const back = docFromJSON(docToJSON(doc));
    expect([...back.curved.values()].every((cv) => cv.step === undefined)).toBe(true);
    // a new reaction arrow taking that id leaves the arrows unapplied
    back.arrows.set(step, { id: step, type: 'arrow', kind: 'reaction', x1: 0, y1: 0, x2: 1, y2: 0 } as never);
    expect(isApplied(back, [...back.curved.values()][0])).toBe(false);
  });
});
