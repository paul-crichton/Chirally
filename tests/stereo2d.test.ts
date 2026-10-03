import { describe, it, expect } from 'vitest';
import { Mol, tetraCcwForOrder } from '../src/chem/mol';
import { parseSmiles, writeSmiles } from '../src/chem/smiles';
import { perceiveStereo2D, assignWedgesFromSpecs } from '../src/chem/stereo2d';

function lAlanineDrawing(style: 'wedge' | 'hash'): Mol {
  // centre (0,0); N lower-left, CH3 up (wedge/hash), COOH lower-right. y is DOWN.
  const m = new Mol();
  const c = m.addAtom({ el: 'C', x: 0, y: 0 });
  const n = m.addAtom({ el: 'N', x: -0.87, y: 0.5 });
  const me = m.addAtom({ el: 'C', x: 0, y: -1 });
  const cc = m.addAtom({ el: 'C', x: 0.87, y: 0.5 });
  const o1 = m.addAtom({ el: 'O', x: 1.74, y: 0 });
  const o2 = m.addAtom({ el: 'O', x: 0.87, y: 1.5 });
  m.addBond(c, n);
  m.addBond(c, me, 1, style);
  m.addBond(c, cc);
  m.addBond(cc, o1, 2);
  m.addBond(cc, o2);
  return m;
}

describe('2D stereo perception', () => {
  it('derives the same chirality as SMILES for L-alanine', () => {
    const m = lAlanineDrawing('wedge');
    const { tetra } = perceiveStereo2D(m);
    expect(tetra.length).toBe(1);
    const ref = parseSmiles('N[C@@H](C)C(=O)O').tetra[0]; // nbrs [N, H, CH3, COOH]
    // map drawing indices: N=1, CH3=2, COOH=3
    const refOrder = [1, -1, 2, 3];
    const drawnCcw = tetraCcwForOrder(tetra[0], refOrder);
    expect(drawnCcw).toBe(ref.ccw);
    expect(writeSmiles(m)).toBe(writeSmiles(parseSmiles('N[C@@H](C)C(=O)O')));
  });
  it('hash inverts the centre', () => {
    const m = lAlanineDrawing('hash');
    perceiveStereo2D(m);
    expect(writeSmiles(m)).toBe(writeSmiles(parseSmiles('N[C@H](C)C(=O)O')));
  });
  it('assigns wedges that reproduce a spec', () => {
    const m = lAlanineDrawing('wedge');
    perceiveStereo2D(m);
    const before = writeSmiles(m);
    // flip spec, re-assign wedges, perceive again
    m.tetra[0].ccw = !m.tetra[0].ccw;
    const flipped = writeSmiles(m);
    assignWedgesFromSpecs(m);
    perceiveStereo2D(m);
    expect(writeSmiles(m)).toBe(flipped);
    expect(flipped).not.toBe(before);
  });
  it('perceives E/Z from coordinates', () => {
    const m = new Mol();
    const a = m.addAtom({ el: 'C', x: 0, y: 0 });
    const b = m.addAtom({ el: 'C', x: 1, y: 0 });
    const ra = m.addAtom({ el: 'C', x: -0.5, y: -0.87 });
    const rb = m.addAtom({ el: 'C', x: 1.5, y: 0.87 });
    m.addBond(a, b, 2);
    m.addBond(a, ra);
    m.addBond(b, rb);
    perceiveStereo2D(m);
    expect(m.dbStereo.length).toBe(1);
    expect(m.dbStereo[0].cis).toBe(false);
    expect(writeSmiles(m)).toBe(writeSmiles(parseSmiles('C/C=C/C')));
  });
});
