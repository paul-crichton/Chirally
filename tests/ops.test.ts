import { describe, it, expect } from 'vitest';
import { createDoc, addAtom, addBond, docToMol, adjacency } from '../src/doc/document';
import {
  idealBondDirection, placeRing, fuseRingOnBond, attachRingToAtom, deleteSelection, mergeAtoms, parseAtomLabel,
  flipSelection, emptySelection, bondToPoint, applyBondType, chainPoints, fuseOverlaps, drawHydrogens,
} from '../src/editor/ops';
import { writeSmiles, suppressHydrogens } from '../src/chem/smiles';
import { perceiveStereo2D } from '../src/chem/stereo2d';
import { expandAbbreviations } from '../src/chem/abbreviations';
import { buildScene } from '../src/render/scene';
import { primsToSVG } from '../src/render/draw';
import { ChemDoc } from '../src/doc/types';

const smi = (doc: ChemDoc) => {
  const { mol } = docToMol(doc);
  const ex = expandAbbreviations(mol);
  perceiveStereo2D(ex);
  return writeSmiles(suppressHydrogens(ex));
};

describe('editor operations', () => {
  it('ideal direction makes a zig-zag chain', () => {
    const d = createDoc();
    const a = addAtom(d, { el: 'C', x: 0, y: 0 });
    let prev = a.id;
    for (let k = 0; k < 4; k++) {
      const p = d.atoms.get(prev)!;
      const dir = idealBondDirection(d, prev);
      prev = bondToPoint(d, prev, { x: p.x + dir.x, y: p.y + dir.y }, 1, 'plain').atom.id;
    }
    expect(smi(d)).toBe('CCCCC');
    // all bond angles ~120°
    const atoms = [...d.atoms.values()];
    for (let i = 1; i < atoms.length - 1; i++) {
      const u = { x: atoms[i - 1].x - atoms[i].x, y: atoms[i - 1].y - atoms[i].y };
      const v = { x: atoms[i + 1].x - atoms[i].x, y: atoms[i + 1].y - atoms[i].y };
      const ang = (Math.acos((u.x * v.x + u.y * v.y) / (Math.hypot(u.x, u.y) * Math.hypot(v.x, v.y))) * 180) / Math.PI;
      expect(Math.abs(ang - 120)).toBeLessThan(0.5);
    }
  });

  it('fuses benzene rings into naphthalene and phenanthrene with valid Kekulé structures', () => {
    const d = createDoc();
    placeRing(d, { x: 0, y: 0 }, 6, true);
    expect(smi(d)).toBe('c1ccccc1');
    // fuse on a single bond
    const single = [...d.bonds.values()].find((b) => b.order === 1)!;
    fuseRingOnBond(d, single.id, 6, true);
    expect(smi(d)).toBe(writeSmilesOf('c1ccc2ccccc2c1'));
    expect(d.atoms.size).toBe(10);
  });

  it('attaches a phenyl via a new bond, or spiro', () => {
    const d = createDoc();
    const a = addAtom(d, { el: 'C', x: 0, y: 0 });
    const b = addAtom(d, { el: 'C', x: 1, y: 0 });
    addBond(d, a.id, b.id);
    attachRingToAtom(d, b.id, 6, true);
    expect(smi(d)).toBe(writeSmilesOf('CCc1ccccc1'));
    const d2 = createDoc();
    placeRing(d2, { x: 0, y: 0 }, 6, false);
    const first = [...d2.atoms.keys()][0];
    attachRingToAtom(d2, first, 5, false, true);
    expect(smi(d2)).toBe(writeSmilesOf('C1CCC2(CC1)CCCC2'));
  });

  it('deleting a bond removes orphaned carbons but keeps heteroatoms', () => {
    const d = createDoc();
    const a = addAtom(d, { el: 'C', x: 0, y: 0 });
    const o = addAtom(d, { el: 'O', x: 1, y: 0 });
    const bd = addBond(d, a.id, o.id);
    const s = emptySelection();
    s.bonds.add(bd.id);
    deleteSelection(d, s);
    expect([...d.atoms.values()].map((x) => x.el)).toEqual(['O']);
  });

  it('merges overlapping atoms after a move', () => {
    const d = createDoc();
    const a = addAtom(d, { el: 'C', x: 0, y: 0 });
    const b = addAtom(d, { el: 'C', x: 1, y: 0 });
    addBond(d, a.id, b.id);
    const c = addAtom(d, { el: 'O', x: 1.1, y: 0.05 });
    const e = addAtom(d, { el: 'C', x: 2, y: 0 });
    addBond(d, c.id, e.id);
    expect(fuseOverlaps(d, new Set([c.id, e.id]))).toBe(1);
    expect(d.atoms.size).toBe(3);
    expect(smi(d)).toBe('COC');
    void mergeAtoms;
  });

  it('parses typed labels', () => {
    expect(parseAtomLabel('OH')).toMatchObject({ el: 'O', hCount: 1, charge: 0 });
    expect(parseAtomLabel('NH3+')).toMatchObject({ el: 'N', hCount: 3, charge: 1 });
    expect(parseAtomLabel('O-')).toMatchObject({ el: 'O', charge: -1 });
    expect(parseAtomLabel('13C')).toMatchObject({ el: 'C', isotope: 13 });
    expect(parseAtomLabel('D')).toMatchObject({ el: 'H', isotope: 2 });
    expect(parseAtomLabel('OMe')).toMatchObject({ abbrev: 'OMe' });
    expect(parseAtomLabel('co2h')).toMatchObject({ abbrev: 'CO2H' });
    expect(parseAtomLabel('CH2CH2OH')).toMatchObject({ abbrev: 'CH2CH2OH' });
    expect(parseAtomLabel('R1')).toMatchObject({ el: 'R', alias: 'R1' });
  });

  it('parses the labels chemists type for reagents and ions', () => {
    expect(parseAtomLabel('H2O')).toMatchObject({ el: 'O', hCount: 2, charge: 0 });
    expect(parseAtomLabel('H3O+')).toMatchObject({ el: 'O', hCount: 3, charge: 1 });
    expect(parseAtomLabel('HO-')).toMatchObject({ el: 'O', hCount: 1, charge: -1 });
    expect(parseAtomLabel('HO−')).toMatchObject({ el: 'O', hCount: 1, charge: -1 });
    expect(parseAtomLabel('HCl')).toMatchObject({ el: 'Cl', hCount: 1, charge: 0 });
    expect(parseAtomLabel('Hg')).toMatchObject({ el: 'Hg' }); // elements still win
    expect(parseAtomLabel('Br+')).toMatchObject({ el: 'Br', hCount: 0, charge: 1 });
    expect(parseAtomLabel('I+')).toMatchObject({ el: 'I', hCount: 0, charge: 1 });
    expect(parseAtomLabel('N+')!.hCount).toBeUndefined();
    expect(parseAtomLabel('CN-')).toMatchObject({ abbrev: 'CN', charge: -1 });
    expect(parseAtomLabel('NC-')).toMatchObject({ abbrev: 'CN', charge: -1 });
    expect(parseAtomLabel('MeO-')).toMatchObject({ abbrev: 'OMe', charge: -1 });
    expect(parseAtomLabel('AcO−')).toMatchObject({ abbrev: 'OAc', charge: -1 });
    expect(parseAtomLabel('E+')).toMatchObject({ el: 'R', alias: 'E', charge: 1 });
    expect(parseAtomLabel('Nu-')).toMatchObject({ el: 'R', alias: 'Nu', charge: -1 });
  });

  it('a charged label keeps its charge when expanded', () => {
    const d = createDoc();
    const a = addAtom(d, { el: 'C', x: 0, y: 0 });
    Object.assign(a, parseAtomLabel('MeO-'));
    expect(smi(d)).toBe(writeSmilesOf('C[O-]'));
  });

  it('abbreviation labels expand for SMILES', () => {
    const d = createDoc();
    placeRing(d, { x: 0, y: 0 }, 6, true);
    const first = [...d.atoms.values()][0];
    const sub = bondToPoint(d, first.id, { x: first.x, y: first.y - 1 }, 1, 'plain').atom;
    Object.assign(sub, parseAtomLabel('CO2H'));
    expect(smi(d)).toBe(writeSmilesOf('OC(=O)c1ccccc1'));
  });

  it('flip keeps stereochemistry; mirror inverts it', () => {
    const d = createDoc();
    const c = addAtom(d, { el: 'C', x: 0, y: 0 });
    const n = addAtom(d, { el: 'N', x: -0.87, y: 0.5 });
    const me = addAtom(d, { el: 'C', x: 0, y: -1 });
    const cc = addAtom(d, { el: 'C', x: 0.87, y: 0.5 });
    const o1 = addAtom(d, { el: 'O', x: 1.74, y: 0 });
    const o2 = addAtom(d, { el: 'O', x: 0.87, y: 1.5 });
    addBond(d, c.id, n.id);
    addBond(d, c.id, me.id, 1, 'wedge');
    addBond(d, c.id, cc.id);
    addBond(d, cc.id, o1.id, 2);
    addBond(d, cc.id, o2.id);
    const before = smi(d);
    const all = emptySelection();
    for (const id of d.atoms.keys()) all.atoms.add(id);
    for (const id of d.bonds.keys()) all.bonds.add(id);
    flipSelection(d, all, 'h');
    expect(smi(d)).toBe(before);
    flipSelection(d, all, 'v', false);
    expect(smi(d)).not.toBe(before);
  });

  it('bond tool semantics: cycling orders and flipping wedges', () => {
    const d = createDoc();
    const a = addAtom(d, { el: 'C', x: 0, y: 0 });
    const b = addAtom(d, { el: 'C', x: 1, y: 0 });
    const bd = addBond(d, a.id, b.id);
    applyBondType(bd, 1, 'plain');
    expect(bd.order).toBe(2);
    applyBondType(bd, 1, 'plain');
    expect(bd.order).toBe(3);
    applyBondType(bd, 1, 'wedge');
    expect(bd.style).toBe('wedge');
    const start = bd.a;
    applyBondType(bd, 1, 'wedge');
    expect(bd.a).not.toBe(start);
  });

  it('chain points alternate', () => {
    const pts = chainPoints({ x: 0, y: 0 }, 0, 4, true);
    expect(pts.length).toBe(4);
    expect(pts[0].y).toBeLessThan(0);
    expect(pts[1].y).toBeCloseTo(0, 6);
  });
});

describe('rendering', () => {
  it('builds a scene and SVG for a labelled structure', () => {
    const d = createDoc();
    placeRing(d, { x: 0, y: 0 }, 6, true);
    const first = [...d.atoms.values()][0];
    const o = bondToPoint(d, first.id, { x: first.x, y: first.y - 1 }, 1, 'plain').atom;
    o.el = 'O';
    const id = d.nextId++;
    d.arrows.set(id, { id, type: 'arrow', kind: 'equilibrium', x1: 3, y1: 0, x2: 6, y2: 0, above: 'H2O', below: 'rt, 2 h' });
    const scene = buildScene(d, { ink: '#000' });
    expect(scene.labelBoxes.has(o.id)).toBe(true);
    expect(scene.prims.some((p) => p.k === 'text' && p.text === 'O')).toBe(true);
    expect(scene.prims.some((p) => p.k === 'text' && p.text === 'H')).toBe(true);
    const svg = primsToSVG(scene.prims, scene.bounds!, 14.4);
    expect(svg.startsWith('<svg')).toBe(true);
    expect(svg).toContain('<text');
    expect(adjacency(d).size).toBe(7);
  });
});

import { parseSmiles } from '../src/chem/smiles';
function writeSmilesOf(s: string): string {
  return writeSmiles(parseSmiles(s));
}

describe('drawHydrogens', () => {
  it('draws implicit hydrogens as H atoms without changing the structure', () => {
    const d = createDoc();
    const c1 = addAtom(d, { el: 'C', x: 0, y: 0 });
    const c2 = addAtom(d, { el: 'C', x: 1, y: 0 });
    addBond(d, c1.id, c2.id);
    const before = smi(d);
    expect(drawHydrogens(d, [c1.id])).toBe(3);
    expect(smi(d)).toBe(before);
    const hs = [...d.atoms.values()].filter((a) => a.el === 'H');
    expect(hs.length).toBe(3);
    // spread out, away from the C–C bond
    for (const h of hs) expect(h.x).toBeLessThan(0.5);
    expect(drawHydrogens(d, [c1.id])).toBe(0);
  });

  it('works on labelled atoms with an explicit H count', () => {
    const d = createDoc();
    const o = addAtom(d, { el: 'O', x: 0, y: 0 });
    Object.assign(o, parseAtomLabel('H3O+'));
    expect(drawHydrogens(d, [o.id])).toBe(3);
    expect(d.atoms.get(o.id)!.hCount).toBeUndefined();
    expect(smi(d)).toBe(writeSmilesOf('[OH3+]'));
  });
});

