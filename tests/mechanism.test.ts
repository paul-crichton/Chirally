import { describe, it, expect } from 'vitest';
import { createDoc, addAtom, addBond, insertMol } from '../src/doc/document';
import { applyArrows, arrowGroups, octetViolations } from '../src/doc/mechanism';
import { parseSmiles, writeSmiles } from '../src/chem/smiles';
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
});
