import { describe, it, expect } from 'vitest';
import { Mol } from '../src/chem/mol';
import { parseSmiles, writeSmiles, suppressHydrogens } from '../src/chem/smiles';
import { computeFormula } from '../src/chem/formula';
import { perceiveStereo2D } from '../src/chem/stereo2d';
import { createDoc, addAtom, addBond, insertMol, newId } from '../src/doc/document';
import type { ChemDoc, ArrowKind } from '../src/doc/types';
import {
  readMolfile, writeMolfile, readSDF, writeSDF, hasCoordinates, readRxn, writeRxn, readCDXML, writeCDXML,
  readCML, writeCML, readXYZ, writeXYZ, detectFormat, FormatError, parseXml,
} from '../src/chem/formats';

// ───────────── helpers ─────────────

const S3 = Math.sqrt(3) / 2;
const canon = (m: Mol, stereo = true) => writeSmiles(suppressHydrogens(m), { stereo });
const smi = (s: string, stereo = true) => writeSmiles(parseSmiles(s), { stereo });

/** L-alanine drawn on an exact hexagonal grid (all bonds length 1), methyl wedged. */
function alanine(style: 'wedge' | 'hash' = 'wedge'): Mol {
  const m = new Mol();
  const c = m.addAtom({ el: 'C', x: 0, y: 0 });
  const n = m.addAtom({ el: 'N', x: -S3, y: 0.5 });
  const me = m.addAtom({ el: 'C', x: 0, y: -1 });
  const cc = m.addAtom({ el: 'C', x: S3, y: 0.5 });
  const o1 = m.addAtom({ el: 'O', x: 2 * S3, y: 0 });
  const o2 = m.addAtom({ el: 'O', x: S3, y: 1.5 });
  m.addBond(c, n);
  m.addBond(c, me, 1, style);
  m.addBond(c, cc);
  m.addBond(cc, o1, 2);
  m.addBond(cc, o2);
  return m;
}

/** A small molecule exercising charges, isotopes, radicals and odd bond styles (unit bonds). */
function decorated(): Mol {
  const m = new Mol();
  const a0 = m.addAtom({ el: 'N', x: 0, y: 0, charge: 1 });
  const a1 = m.addAtom({ el: 'O', x: S3, y: 0.5, charge: -1 });
  const a2 = m.addAtom({ el: 'C', x: -S3, y: 0.5, isotope: 13 });
  const a3 = m.addAtom({ el: 'C', x: -2 * S3, y: 0, radical: 1 });
  const a4 = m.addAtom({ el: 'C', x: 0, y: -1 });
  const a5 = m.addAtom({ el: 'C', x: S3, y: -1.5, radical: 2 });
  const a6 = m.addAtom({ el: 'H', x: -S3, y: 1.5, isotope: 2 });
  const a7 = m.addAtom({ el: 'C', x: -2 * S3, y: -1 });
  const a8 = m.addAtom({ el: 'C', x: -3 * S3, y: -1.5 });
  m.addBond(a0, a1);
  m.addBond(a0, a2);
  m.addBond(a2, a3);
  m.addBond(a0, a4, 1, 'hash');
  m.addBond(a4, a5, 1, 'wavy');
  m.addBond(a2, a6);
  m.addBond(a3, a7);
  m.addBond(a7, a8, 2, 'crossed');
  return m;
}

function expectSameMol(a: Mol, b: Mol, coordTol = 1e-3) {
  expect(b.atoms.length).toBe(a.atoms.length);
  expect(b.bonds.length).toBe(a.bonds.length);
  a.atoms.forEach((x, i) => {
    const y = b.atoms[i];
    expect(y.el, `atom ${i} element`).toBe(x.el);
    expect(y.charge, `atom ${i} charge`).toBe(x.charge);
    expect(y.isotope, `atom ${i} isotope`).toBe(x.isotope);
    expect(y.radical, `atom ${i} radical`).toBe(x.radical);
    expect(Math.abs(y.x - x.x), `atom ${i} x`).toBeLessThan(coordTol);
    expect(Math.abs(y.y - x.y), `atom ${i} y`).toBeLessThan(coordTol);
  });
  a.bonds.forEach((x, i) => {
    const y = b.bonds[i];
    expect([y.a, y.b, y.order, y.style], `bond ${i}`).toEqual([x.a, x.b, x.order, x.style]);
  });
}

// ───────────── molfile ─────────────

describe('Molfile V2000', () => {
  it('round-trips atoms, bonds, orders, charges, isotopes, radicals, wedges and coordinates', () => {
    for (const m of [alanine('wedge'), alanine('hash'), decorated()]) {
      const text = writeMolfile(m);
      expect(text).toContain('V2000');
      expect(text).toContain('M  END');
      const back = readMolfile(text);
      expectSameMol(m, back);
    }
  });

  it('writes standard fixed-width fields', () => {
    const text = writeMolfile(decorated(), { title: 'decorated' });
    const lines = text.split('\n');
    expect(lines[0]).toBe('decorated');
    expect(lines[1].slice(20, 22)).toBe('2D');
    expect(lines[3]).toMatch(/^  9  8  0  0  1  0  0  0  0  0999 V2000$/);
    // coordinates are scaled to 1.5 Å with y up
    const n = lines[4], o = lines[5];
    expect(n.slice(31, 34).trim()).toBe('N');
    expect(Number(o.slice(0, 10))).toBeCloseTo(1.5 * S3, 4);
    expect(Number(o.slice(10, 20))).toBeCloseTo(-0.75, 4);
    expect(text).toMatch(/M {2}CHG {2}2 {3}1 {3}1 {3}2 {2}-1/);
    expect(text).toMatch(/M {2}ISO {2}2 {3}3 {2}13 {3}7 {3}2/);
    expect(text).toMatch(/M {2}RAD {2}2 {3}4 {3}2 {3}6 {3}3/);
    // hashed bond N→C has stereo 6, wavy 4, crossed double 3
    expect(text).toMatch(/\n {2}1 {2}5 {2}1 {2}6/);
    expect(text).toMatch(/\n {2}5 {2}6 {2}1 {2}4/);
    expect(text).toMatch(/\n {2}8 {2}9 {2}2 {2}3/);
  });

  it('preserves canonical SMILES for a range of molecules', () => {
    const list = [
      'CCO', 'c1ccccc1', 'CC(=O)Oc1ccccc1C(=O)O', 'C[N+](C)(C)C', '[O-][N+](=O)c1ccccc1', 'c1ccc2[nH]ccc2c1',
      'c1ccncc1', 'O=C1CCCCC1', 'C#N', 'C=CC=C', '[2H]C([2H])([2H])O', '[13CH4]', '[CH3]', 'CS(=O)(=O)O',
      'Cn1cnc2c1c(=O)n(C)c(=O)n2C', 'c1ccsc1', '[Na+].[Cl-]', 'OC(=O)CC(O)(CC(=O)O)C(=O)O', 'B(O)(O)c1ccccc1',
    ];
    for (const s of list) {
      const m = parseSmiles(s);
      const back = readMolfile(writeMolfile(m));
      expect(canon(back), s).toBe(smi(s));
    }
  });

  it('keeps stereo from wedges (tetrahedral) and coordinates (E/Z)', () => {
    const m = alanine();
    perceiveStereo2D(m);
    const ref = writeSmiles(m);
    expect(ref).toContain('@');
    expect(ref).toBe(smi('N[C@@H](C)C(=O)O'));
    const back = readMolfile(writeMolfile(m));
    expect(back.tetra.length).toBe(1);
    expect(writeSmiles(back)).toBe(ref);
    // E-2-butene drawn trans
    const e = new Mol();
    e.addAtom({ el: 'C', x: -S3, y: -0.5 });
    e.addAtom({ el: 'C', x: 0, y: 0 });
    e.addAtom({ el: 'C', x: S3, y: -0.5 });
    e.addAtom({ el: 'C', x: 2 * S3, y: 0 });
    e.addBond(0, 1);
    e.addBond(1, 2, 2);
    e.addBond(2, 3);
    expect(writeSmiles(readMolfile(writeMolfile(e)))).toBe(smi('C/C=C/C'));
  });

  it('derives wedges when a molecule only has stereo specs', () => {
    const m = alanine();
    perceiveStereo2D(m);
    const ref = writeSmiles(m);
    m.bonds[1].style = 'plain';
    const back = readMolfile(writeMolfile(m));
    expect(back.bonds.some((b) => b.style === 'wedge' || b.style === 'hash')).toBe(true);
    expect(writeSmiles(back)).toBe(ref);
  });

  it('handles aliases, R groups, explicit H counts and aromatic bonds', () => {
    const text = `alias test
  manual

  7  7  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    2.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    3.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    2.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    0.7500    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
    2.5981    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0  0  0  0
  2  3  4  0  0  0  0
  3  4  4  0  0  0  0
  4  5  4  0  0  0  0
  5  1  4  0  0  0  0
  1  6  1  0  0  0  0
  2  7  1  0  0  0  0
A    7
OMe
M  RGP  1   6   1
M  END
`;
    const m = readMolfile(text);
    expect(m.name).toBe('alias test');
    expect(m.atoms[5].el).toBe('R');
    expect(m.atoms[5].alias).toBe('R1');
    expect(m.atoms[6].abbrev).toBe('OMe');
    // cyclopentadienyl-like ring with 5 aromatic bonds cannot be kekulized cleanly → stays sane
    expect(m.bonds.every((b) => [1, 1.5, 2].includes(b.order))).toBe(true);
    // benzene with aromatic bond type 4 is kekulized
    const benz = readMolfile(writeMolfile(parseSmiles('c1ccccc1')).replace(/( {2}[1-6] {2}[1-6]) {2}[12]/g, '$1  4'));
    expect(benz.bonds.map((b) => b.order).sort()).toEqual([1, 1, 1, 2, 2, 2]);
  });

  it('round-trips abbreviations as superatom S-groups', () => {
    const m = alanine();
    m.atoms[4].abbrev = 'OMe';
    m.atoms[4].el = 'O';
    m.bonds[3].order = 1;
    const text = writeMolfile(m);
    expect(text).toContain('M  STY  1   1 SUP');
    expect(text).toContain('M  SMT   1 OMe');
    const back = readMolfile(text);
    expect(back.atoms.length).toBe(6);
    expect(back.atoms[4].abbrev).toBe('OMe');
    expect(Math.abs(back.atoms[4].x - m.atoms[4].x)).toBeLessThan(1e-3);
    // without expansion the label is stored as an alias atom and read back as an abbreviation
    const alias = writeMolfile(m, { expandAbbreviations: false });
    expect(alias).toMatch(/\nA {4}5\nOMe\n/);
    expect(readMolfile(alias).atoms[4].abbrev).toBe('OMe');
    // the expanded file has the full structure for other programs
    const full = readMolfile(text.replace(/M {2}S(TY|AL|BL|MT).*\n/g, ''));
    expect(computeFormula(full).formula).toBe(computeFormula(back).formula);
    expect(full.atoms.length).toBe(7);
  });

  it('leaves coordinate-less files at the origin and flags 3D files', () => {
    const m = readMolfile(writeMolfile(parseSmiles('CCO')));
    expect(hasCoordinates(m)).toBe(false);
    expect(m.atoms.every((a) => a.x === 0 && a.y === 0)).toBe(true);
    const three = `3d
  ChemWrit00000000003D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.1000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7570    0.5860    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7570    0.5860   -0.2000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END
`;
    const w = readMolfile(three);
    expect(w.props.dim).toBe('3D');
    expect(w.atoms[1].x).toBeCloseTo(0.757, 6); // not rescaled
    expect(w.atoms[1].y).toBeCloseTo(-0.586, 6); // y flipped
    expect(w.atoms[2].z).toBeCloseTo(-0.2, 6);
    const again = readMolfile(writeMolfile(w));
    expect(again.props.dim).toBe('3D');
    expectSameMol(w, again);
    expect(again.atoms[0].z).toBeCloseTo(0.1, 4);
  });
});

describe('Molfile V3000', () => {
  it('reads a V3000 CTAB with charges, isotopes, radicals, stereo and continuation lines', () => {
    const text = `acetate
  ChemWrit00000000002D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 2 C 1.299 0.75 0 0 MASS=13
M  V30 3 O 2.598 0 0 0 CHG=-1
M  V30 4 O 1.299 2.25 0 0
M  V30 5 C -1.299 0.75 0 0 RAD=2 -
M  V30 CHG=0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 2 2 4
M  V30 4 1 1 5 CFG=1
M  V30 END BOND
M  V30 END CTAB
M  END
`;
    const m = readMolfile(text);
    expect(m.name).toBe('acetate');
    expect(m.atoms.map((a) => a.el)).toEqual(['C', 'C', 'O', 'O', 'C']);
    expect(m.atoms[1].isotope).toBe(13);
    expect(m.atoms[2].charge).toBe(-1);
    expect(m.atoms[4].radical).toBe(1);
    expect(m.bonds[2].order).toBe(2);
    expect(m.bonds[3].style).toBe('wedge');
    expect(m.bonds[3].a).toBe(0);
    // median bond length 1.5 Å → 1 unit, y flipped
    expect(m.atoms[3].x).toBeCloseTo(1.299 / 1.5, 3);
    expect(m.atoms[3].y).toBeCloseTo(-1.5, 3);
  });

  it('round-trips through the V3000 writer', () => {
    for (const m of [alanine(), decorated()]) {
      const text = writeMolfile(m, { v3000: true });
      expect(text).toContain('M  V30 BEGIN CTAB');
      expectSameMol(m, readMolfile(text));
    }
    const sup = alanine();
    sup.atoms[4].abbrev = 'OMe';
    sup.bonds[3].order = 1;
    const back = readMolfile(writeMolfile(sup, { v3000: true }));
    expect(back.atoms[4].abbrev).toBe('OMe');
  });

  it('switches to V3000 automatically above 999 atoms', () => {
    const m = new Mol();
    for (let i = 0; i < 1200; i++) {
      m.addAtom({ el: 'C', x: i * S3, y: i % 2 ? 0.5 : 0 });
      if (i) m.addBond(i - 1, i);
    }
    const text = writeMolfile(m);
    expect(text).toContain('V3000');
    const back = readMolfile(text);
    expect(back.atoms.length).toBe(1200);
    expect(back.bonds.length).toBe(1199);
  });
});

// ───────────── SDF ─────────────

const ASPIRIN_SDF = `2244
  -OEChem-10032612062D

 21 21  0     0  0  0  0  0  0999 V2000
    3.7321   -0.0600    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.3301    1.4400    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5981    1.4400    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8660   -1.5600    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5981   -0.5600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.4641   -0.0600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5981   -1.5600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3301   -0.5600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.4641   -2.0600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3301   -1.5600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.4641    0.9400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8660   -0.5600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000   -0.0600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0611   -1.8700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.8671   -0.2500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.4641   -2.6800    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.8671   -1.8700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.3100    0.4769    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.4631    0.2500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.6900   -0.5969    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.3301    2.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  5  1  0  0  0  0
  1 12  1  0  0  0  0
  2 11  1  0  0  0  0
  2 21  1  0  0  0  0
  3 11  2  0  0  0  0
  4 12  2  0  0  0  0
  5  6  1  0  0  0  0
  5  7  2  0  0  0  0
  6  8  2  0  0  0  0
  6 11  1  0  0  0  0
  7  9  1  0  0  0  0
  7 14  1  0  0  0  0
  8 10  1  0  0  0  0
  8 15  1  0  0  0  0
  9 10  2  0  0  0  0
  9 16  1  0  0  0  0
 10 17  1  0  0  0  0
 12 13  1  0  0  0  0
 13 18  1  0  0  0  0
 13 19  1  0  0  0  0
 13 20  1  0  0  0  0
M  END
> <PUBCHEM_COMPOUND_CID>
2244

> <PUBCHEM_COMPOUND_CANONICALIZED>
1

> <PUBCHEM_CACTVS_COMPLEXITY>
212

> <PUBCHEM_CACTVS_HBOND_ACCEPTOR>
4

> <PUBCHEM_CACTVS_HBOND_DONOR>
1

> <PUBCHEM_CACTVS_ROTATABLE_BOND>
3

> <PUBCHEM_CACTVS_SUBSKEYS>
AAADccBwOAAAAAAAAAAAAAAAAAAAAAAAAAAwAAAAAAAAAAABAAAAGgAACAAADASAmAAyDoAABgCIAiDSCAACCAAkIAAIiAEGCMgMJzaENRqCe2Cl4BEIuYeIyCCOAAAAAAAIAAAAAAAAABAAAAAAAAAAAA==

> <PUBCHEM_IUPAC_OPENEYE_NAME>
2-acetoxybenzoic acid

> <PUBCHEM_IUPAC_CAS_NAME>
2-acetyloxybenzoic acid

> <PUBCHEM_IUPAC_NAME_MARKUP>
2-acetyloxybenzoic acid

> <PUBCHEM_IUPAC_NAME>
2-acetyloxybenzoic acid

> <PUBCHEM_IUPAC_SYSTEMATIC_NAME>
2-acetyloxybenzoic acid

> <PUBCHEM_IUPAC_TRADITIONAL_NAME>
2-acetoxybenzoic acid

> <PUBCHEM_IUPAC_INCHI>
InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)

> <PUBCHEM_IUPAC_INCHIKEY>
BSYNRYMUTXBXSQ-UHFFFAOYSA-N

> <PUBCHEM_XLOGP3>
1.2

> <PUBCHEM_EXACT_MASS>
180.04225873

> <PUBCHEM_MOLECULAR_FORMULA>
C9H8O4

> <PUBCHEM_MOLECULAR_WEIGHT>
180.16

> <PUBCHEM_SMILES>
CC(=O)OC1=CC=CC=C1C(=O)O

> <PUBCHEM_CONNECTIVITY_SMILES>
CC(=O)OC1=CC=CC=C1C(=O)O

> <PUBCHEM_CACTVS_TPSA>
63.6

> <PUBCHEM_MONOISOTOPIC_WEIGHT>
180.04225873

> <PUBCHEM_TOTAL_CHARGE>
0

> <PUBCHEM_HEAVY_ATOM_COUNT>
13

> <PUBCHEM_ATOM_DEF_STEREO_COUNT>
0

> <PUBCHEM_ATOM_UDEF_STEREO_COUNT>
0

> <PUBCHEM_BOND_DEF_STEREO_COUNT>
0

> <PUBCHEM_BOND_UDEF_STEREO_COUNT>
0

> <PUBCHEM_ISOTOPIC_ATOM_COUNT>
0

> <PUBCHEM_COMPONENT_COUNT>
1

> <PUBCHEM_CACTVS_TAUTO_COUNT>
-1

> <PUBCHEM_COORDINATE_TYPE>
1
5
255

> <PUBCHEM_BONDANNOTATIONS>
5  6  8
5  7  8
6  8  8
7  9  8
8  10  8
9  10  8

$$$$
`;

describe('SD files', () => {
  it('reads multiple records with data fields', () => {
    const a = parseSmiles('CCO');
    a.name = 'ethanol';
    a.props = { MW: '46.07', NOTE: 'line one\nline two' };
    const b = alanine();
    b.name = 'alanine';
    b.props = { ID: 'A-2' };
    const text = writeSDF([a, b]);
    expect(text.match(/\$\$\$\$/g)!.length).toBe(2);
    const mols = readSDF(text);
    expect(mols.length).toBe(2);
    expect(mols[0].name).toBe('ethanol');
    expect(mols[0].props).toEqual({ MW: '46.07', NOTE: 'line one\nline two' });
    expect(mols[1].name).toBe('alanine');
    expect(mols[1].props.ID).toBe('A-2');
    expect(canon(mols[0])).toBe(smi('CCO'));
    expectSameMol(b, mols[1]);
  });

  it('reads a PubChem record (aspirin, CID 2244)', () => {
    const mols = readSDF(ASPIRIN_SDF);
    expect(mols.length).toBe(1);
    const m = mols[0];
    expect(m.name).toBe('2244');
    expect(m.atoms.length).toBe(21);
    expect(m.bonds.length).toBe(21);
    expect(m.props.PUBCHEM_COMPOUND_CID).toBe('2244');
    expect(m.props.PUBCHEM_IUPAC_NAME).toBe('2-acetyloxybenzoic acid');
    expect(m.props.PUBCHEM_COORDINATE_TYPE).toBe('1\n5\n255');
    expect(computeFormula(m).formula).toBe('C9H8O4');
    expect(canon(m)).toBe(smi('CC(=O)Oc1ccccc1C(=O)O'));
    // heavy-atom bonds were ~1 Å: median bond length is normalised to 1
    const lens = m.bonds.map((b) => Math.hypot(m.atoms[b.a].x - m.atoms[b.b].x, m.atoms[b.a].y - m.atoms[b.b].y)).sort((p, q) => p - q);
    expect(lens[10]).toBeCloseTo(1, 6);
    // y is flipped: O1 (file y = -0.06) is above C12 (file y = -0.56), i.e. it has the smaller screen y
    expect(m.atoms[0].y).toBeLessThan(m.atoms[11].y);
    // readMolfile on an SD file returns the first record with its data
    expect(readMolfile(ASPIRIN_SDF).props.PUBCHEM_MOLECULAR_FORMULA).toBe('C9H8O4');
    // and writing it back preserves everything
    const again = readSDF(writeSDF(mols));
    expect(again[0].props.PUBCHEM_IUPAC_INCHIKEY).toBe('BSYNRYMUTXBXSQ-UHFFFAOYSA-N');
    expectSameMol(m, again[0]);
  });
});

// ───────────── RXN ─────────────

describe('RXN', () => {
  it('round-trips reactants, products and agents', () => {
    const r = {
      reactants: [parseSmiles('CCO'), parseSmiles('CC(=O)O')],
      products: [parseSmiles('CCOC(C)=O'), parseSmiles('O')],
      agents: [parseSmiles('OS(=O)(=O)O')],
    };
    const text = writeRxn(r);
    expect(text.startsWith('$RXN')).toBe(true);
    expect(text.split('\n')[4]).toBe('  2  2  1');
    const back = readRxn(text);
    expect(back.reactants.map((m) => canon(m))).toEqual(r.reactants.map((m) => writeSmiles(m)));
    expect(back.products.map((m) => canon(m))).toEqual(r.products.map((m) => writeSmiles(m)));
    expect(back.agents.map((m) => canon(m))).toEqual(r.agents.map((m) => writeSmiles(m)));
  });

  it('keeps the relative placement of components', () => {
    const a = alanine();
    const b = alanine();
    b.translate(5, 0);
    const back = readRxn(writeRxn({ reactants: [a], products: [b] }));
    expect(back.products[0].atoms[0].x - back.reactants[0].atoms[0].x).toBeCloseTo(5, 3);
    expectSameMol(a, back.reactants[0]);
  });

  it('reads V3000 reactions', () => {
    const text = `$RXN V3000

      ChemWrite100320261200

M  V30 COUNTS 1 1
M  V30 BEGIN REACTANT
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 1
M  V30 2 O 1.5 0 0 2
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 END CTAB
M  V30 END REACTANT
M  V30 BEGIN PRODUCT
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 6 0 0 1
M  V30 2 O 7.5 0 0 2
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 END BOND
M  V30 END CTAB
M  V30 END PRODUCT
M  END
`;
    const r = readRxn(text);
    expect(r.reactants.length).toBe(1);
    expect(r.products.length).toBe(1);
    expect(r.products[0].atoms[1].map).toBe(2);
    for (const m of [...r.reactants, ...r.products]) for (const a of m.atoms) delete a.map;
    expect(canon(r.reactants[0])).toBe('CO');
    expect(canon(r.products[0])).toBe('C=O');
    expect(r.products[0].atoms[0].x).toBeCloseTo(4, 6);
  });
});

// ───────────── CDXML ─────────────

/** Hexagon node positions (pt) around (cx, cy), bond length 14.4, flat top. */
function hexNodes(cx: number, cy: number): [number, number][] {
  const out: [number, number][] = [];
  for (let k = 0; k < 6; k++) {
    const a = (Math.PI / 3) * k + Math.PI / 6;
    out.push([cx + 14.4 * Math.cos(a), cy + 14.4 * Math.sin(a)]);
  }
  return out;
}

const HEX = hexNodes(100, 100);
const f = (v: number) => v.toFixed(2);
const CDXML_FIXTURE = `<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >
<CDXML
 CreationProgram="ChemDraw 20.1.1.125"
 Name="fixture.cdxml"
 BoundingBox="70 60 380 200"
 WindowPosition="0 0"
 FractionalWidths="yes"
 InterpretChemically="yes"
 LabelFont="3"
 LabelSize="10"
 LabelFace="96"
 CaptionFont="3"
 CaptionSize="10"
 HashSpacing="2.50"
 MarginWidth="1.60"
 LineWidth="0.60"
 BoldWidth="2"
 BondLength="14.40"
 BondSpacing="18"
 ChainAngle="120"
 color="0"
 bgcolor="1"
><colortable>
<color r="1" g="1" b="1"/>
<color r="0" g="0" b="0"/>
<color r="1" g="0" b="0"/>
</colortable><fonttable>
<font id="3" charset="iso-8859-1" name="Arial"/>
</fonttable><page
 id="1"
 BoundingBox="0 0 540 720"
 HeaderPosition="36"
 FooterPosition="36"
 PrintTrimMarks="yes"
 HeightPages="1"
 WidthPages="1"
><!-- benzene with a wedged methyl and a methoxy nickname -->
<fragment
 id="10"
 BoundingBox="80 70 150 130"
 Z="1"
>${HEX.map(([x, y], k) => `<n id="${11 + k}" p="${f(x)} ${f(y)}" Z="${2 + k}" AS="N"/>`).join('\n')}
<n id="17" p="${f(HEX[0][0] + 14.4)} ${f(HEX[0][1])}" Z="9" AS="N"/>
<n id="18" p="${f(HEX[3][0] - 14.4)} ${f(HEX[3][1])}" Z="10" NodeType="Nickname" NeedsClean="yes" AS="N"
><fragment id="40"><n id="41" p="70 100" Element="8" NumHydrogens="0"/><n id="42" p="60 95"/><n id="43" p="86 100" NodeType="ExternalConnectionPoint"/><b id="44" B="41" E="42"/><b id="45" B="43" E="41"/></fragment
><t p="${f(HEX[3][0] - 14.4 - 3.3)} ${f(HEX[3][1] + 3.5)}" BoundingBox="0 0 1 1" LabelJustification="Right"><s font="3" size="10" color="0" face="96">MeO</s></t></n>
<n id="19" p="${f(HEX[1][0])} ${f(HEX[1][1] + 14.4)}" Z="11" Element="7" NumHydrogens="2" Charge="1" color="4"
><t p="0 0"><s font="3" size="10" color="4" face="96">NH2+</s></t></n>
<b id="20" Z="12" B="11" E="12"/>
<b id="21" Z="13" B="12" E="13" Order="2" DoublePosition="Right"/>
<b id="22" Z="14" B="13" E="14"/>
<b id="23" Z="15" B="14" E="15" Order="2"/>
<b id="24" Z="16" B="15" E="16"/>
<b id="25" Z="17" B="16" E="11" Order="2"/>
<b id="26" Z="18" B="17" E="11" Display="WedgeEnd"/>
<b id="27" Z="19" B="14" E="18"/>
<b id="28" Z="20" B="12" E="19" Display="WedgedHashBegin"/>
<b id="29" Z="21" B="99" E="11"/>
</fragment>
<t id="30" p="190 85" BoundingBox="190 76 230 88" Z="22" Justification="Center" InterpretChemically="no"
><s font="3" size="10" color="0" face="0">H</s><s font="3" size="10" color="0" face="32">2</s><s font="3" size="10" color="0" face="0">O, </s><s font="3" size="10" color="0" face="2">heat</s></t>
<arrow
 id="31"
 BoundingBox="170 96 262 104"
 Z="23"
 FillType="None"
 ArrowheadHead="Full"
 ArrowheadType="Solid"
 HeadSize="1000"
 ArrowheadCenterSize="875"
 ArrowheadWidth="250"
 Head3D="260 100 0"
 Tail3D="172 100 0"
/>
<graphic id="32" SupersededBy="31" BoundingBox="260 100 172 100" Z="24" GraphicType="Line" ArrowType="FullHead" HeadSize="1000"/>
<t id="33" p="290 104" Z="25"><s font="3" size="12" color="3" face="96">CH3CO2H</s></t>
<graphic id="34" BoundingBox="280 140 360 180" Z="26" GraphicType="Rectangle" RectangleType="RoundEdge Dashed"/>
<arrow id="35" BoundingBox="100 140 160 190" Z="27" FillType="None" ArrowheadHead="HalfLeft" ArrowheadType="Solid" HeadSize="1000"
 AngularSize="120" Head3D="120 177.32 0" Tail3D="150 160 0" Center3D="130 160 0" MajorAxisEnd3D="150 160 0" MinorAxisEnd3D="130 180 0"/>
<graphic id="36" BoundingBox="270 100 276 106" Z="28" GraphicType="Symbol" SymbolType="Plus"/>
<graphic id="37" BoundingBox="300 60 350 60" Z="29" GraphicType="Line" ArrowType="Equilibrium"/>
<graphic id="38" BoundingBox="60 60 160 200" Z="30" GraphicType="Bracket" BracketType="SquarePair"/>
<curve id="39" Z="31" CurvePoints="200 150 210 140 230 140 240 150" ArrowheadHead="Full"/>
<scheme id="50"><step id="51" ReactionStepArrows="31" ReactionStepObjectsAboveArrow="30"/></scheme>
<unknownthing foo="bar"><n id="999"/></unknownthing>
</page></CDXML>
`;

describe('CDXML reading', () => {
  const doc = readCDXML(CDXML_FIXTURE);
  const atoms = [...doc.atoms.values()];
  const bonds = [...doc.bonds.values()];

  it('reads atoms and bonds, scaling by BondLength', () => {
    expect(doc.style.bondLengthPt).toBeCloseTo(14.4);
    expect(atoms.length).toBe(9); // 6 ring + methyl + OMe + NH2+
    expect(bonds.length).toBe(9); // dangling bond to unknown node 99 skipped
    expect(atoms[0].x).toBeCloseTo(HEX[0][0] / 14.4, 2);
    expect(atoms[0].y).toBeCloseTo(HEX[0][1] / 14.4, 2);
    const orders = bonds.slice(0, 6).map((b) => b.order);
    expect(orders).toEqual([1, 2, 1, 2, 1, 2]);
    expect(bonds[1].dbPos).toBe('right');
  });

  it('reads nickname nodes as abbreviations and heteroatom properties', () => {
    const ome = atoms[7];
    expect(ome.abbrev).toBe('OMe'); // label drawn reversed as "MeO" (right-justified)
    expect(ome.el).toBe('O');
    const n = atoms[8];
    expect(n.el).toBe('N');
    expect(n.charge).toBe(1);
    expect(n.hCount).toBe(2); // NumHydrogens kept because N+ with one bond would default to NH3+
    expect(n.color).toBe('#ff0000');
    // ring carbons carry no explicit H count
    expect(atoms.slice(0, 6).every((a) => a.hCount === undefined)).toBe(true);
  });

  it('reads wedge orientation (narrow end = bond.a)', () => {
    const wedge = bonds.find((b) => b.style === 'wedge')!;
    expect(wedge).toBeDefined();
    expect(wedge.a).toBe(atoms[0].id); // WedgeEnd: narrow end at E (ring atom 11)
    expect(wedge.b).toBe(atoms[6].id);
    const hash = bonds.find((b) => b.style === 'hash')!;
    expect(hash.a).toBe(atoms[1].id);
    expect(hash.b).toBe(atoms[8].id);
  });

  it('reads text with markup, formula text and step captions', () => {
    const texts = [...doc.texts.values()];
    const formula = texts.find((t) => t.formula)!;
    expect(formula.text).toBe('CH3CO2H');
    expect(formula.size).toBeCloseTo(1.2);
    expect(formula.x).toBeCloseTo(290 / 14.4);
    expect(texts.some((t) => t.text === '+')).toBe(true);
    // the "H2O, heat" caption is attached to the arrow through the <step>
    expect(texts.some((t) => t.text.includes('heat'))).toBe(false);
    const arrows = [...doc.arrows.values()];
    const rxn = arrows.find((a) => a.kind === 'reaction')!;
    expect(rxn.above).toBe('H_2O, *heat*');
  });

  it('reads arrows, curved arrows and shapes', () => {
    const arrows = [...doc.arrows.values()];
    expect(arrows.length).toBe(2); // the superseded <graphic> is not duplicated
    const rxn = arrows.find((a) => a.kind === 'reaction')!;
    expect(rxn.x1).toBeCloseTo(172 / 14.4);
    expect(rxn.x2).toBeCloseTo(260 / 14.4);
    expect(rxn.y1).toBeCloseTo(100 / 14.4);
    const eq = arrows.find((a) => a.kind === 'equilibrium')!;
    expect(eq.x2).toBeCloseTo(300 / 14.4); // graphic lines list the head first
    const curved = [...doc.curved.values()];
    expect(curved.length).toBe(2);
    const arc = curved.find((c) => c.electrons === 1)!;
    expect(arc.from).toEqual({ type: 'point', x: 150 / 14.4, y: 160 / 14.4 });
    if (arc.to.type !== 'point') throw new Error('expected a point anchor');
    expect(arc.to.x).toBeCloseTo(120 / 14.4);
    // the Bézier midpoint lies on the arc (radius 20 pt around the centre)
    const S = { x: 150 / 14.4, y: 160 / 14.4 }, E = { x: arc.to.x, y: arc.to.y };
    const len = Math.hypot(E.x - S.x, E.y - S.y);
    const u = { x: (E.x - S.x) / len, y: (E.y - S.y) / len };
    const cp = (q: { t: number; h: number }) => ({ x: S.x + u.x * q.t * len - u.y * q.h * len, y: S.y + u.y * q.t * len + u.x * q.h * len });
    const c1 = cp(arc.c1), c2 = cp(arc.c2);
    const mid = { x: (S.x + 3 * c1.x + 3 * c2.x + E.x) / 8, y: (S.y + 3 * c1.y + 3 * c2.y + E.y) / 8 };
    expect(Math.hypot(mid.x - 130 / 14.4, mid.y - 160 / 14.4)).toBeCloseTo(20 / 14.4, 2);
    expect(curved.find((c) => c.electrons === 2)).toBeDefined();
    const shapes = [...doc.shapes.values()];
    const rect = shapes.find((s) => s.kind === 'roundRect')!;
    expect(rect.dashed).toBe(true);
    expect(rect.x2).toBeCloseTo(360 / 14.4);
    expect(shapes.some((s) => s.kind === 'bracket')).toBe(true);
  });
});

function sampleDoc(): ChemDoc {
  const doc = createDoc();
  doc.meta.title = 'round trip';
  const m = alanine();
  m.atoms[1].charge = 1; // N+
  m.atoms[4].abbrev = 'OMe';
  m.bonds[3].order = 1;
  m.bonds[4].order = 2;
  m.bonds[4].dbPos = 'left';
  const { atomIds, bondIds } = insertMol(doc, m, 2, 3);
  const r = addAtom(doc, { el: 'R', alias: 'R1', x: 2, y: 5 });
  addBond(doc, atomIds[5], r.id, 1, 'bold');
  const iso = addAtom(doc, { el: 'C', x: 2 + 2 * S3, y: 1, isotope: 13 });
  addBond(doc, atomIds[3], iso.id);
  const add = <T extends { id: number }>(map: Map<number, T>, obj: Omit<T, 'id'>) => {
    const o = { ...obj, id: newId(doc) } as T;
    map.set(o.id, o);
    return o;
  };
  add(doc.texts, { type: 'text', x: 1, y: 8, text: 'CH_3CO_2H yields *product* **1**' });
  add(doc.texts, { type: 'text', x: 6, y: 8, text: 'Ph3P', formula: true, size: 1.5, align: 'center', color: '#0000ff' });
  const kinds: ArrowKind[] = ['reaction', 'equilibrium', 'unbalancedEq', 'retro', 'resonance', 'dashed', 'noGo', 'line'];
  kinds.forEach((kind, k) => add(doc.arrows, { type: 'arrow', kind, x1: 6, y1: 2 + k, x2: 9, y2: 2 + k, ...(kind === 'reaction' ? { above: 'H_2O', below: 'heat, 2 h' } : {}) }));
  // a curved arrow along a 90° arc from atom N to the C–C bond
  const k = (4 / 3) * Math.tan(Math.PI / 8);
  add(doc.curved, { type: 'curved', electrons: 2, from: { type: 'atom', id: atomIds[1] }, to: { type: 'bond', id: bondIds[2] }, c1: { t: k / 2, h: k / 2 }, c2: { t: 1 - k / 2, h: k / 2 } });
  add(doc.curved, { type: 'curved', electrons: 1, from: { type: 'point', x: 0, y: 0 }, to: { type: 'point', x: 1, y: 1 }, c1: { t: 0.3, h: -0.4 }, c2: { t: 0.7, h: -0.4 } });
  // exact half circles bulging to either side (the ambiguous case for arc arrows)
  add(doc.curved, { type: 'curved', electrons: 2, from: { type: 'point', x: 3, y: 9 }, to: { type: 'point', x: 5, y: 9 }, c1: { t: 0, h: 2 / 3 }, c2: { t: 1, h: 2 / 3 } });
  add(doc.curved, { type: 'curved', electrons: 2, from: { type: 'point', x: 3, y: 10 }, to: { type: 'point', x: 4, y: 11 }, c1: { t: 0, h: -2 / 3 }, c2: { t: 1, h: -2 / 3 } });
  add(doc.shapes, { type: 'shape', kind: 'rect', x1: 10, y1: 1, x2: 12, y2: 3 });
  add(doc.shapes, { type: 'shape', kind: 'ellipse', x1: 10, y1: 4, x2: 13, y2: 6, dashed: true });
  add(doc.shapes, { type: 'shape', kind: 'paren', x1: 0.5, y1: 1, x2: 4, y2: 6 });
  return doc;
}

describe('CDXML writing', () => {
  it('produces a well-formed ChemDraw document', () => {
    const xml = writeCDXML(sampleDoc());
    expect(xml).toContain('<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >');
    expect(xml).toMatch(/<CDXML[^>]*BondLength="18"/);
    expect(xml).toContain('NodeType="Nickname"');
    expect(xml).toContain('NodeType="ExternalConnectionPoint"');
    expect(xml).toContain('Display="WedgeBegin"');
    expect(xml).toContain('AngularSize=');
    const root = parseXml(xml);
    expect(root.name).toBe('CDXML');
  });

  it('round-trips a document', () => {
    const doc = sampleDoc();
    const back = readCDXML(writeCDXML(doc));
    const A = [...doc.atoms.values()], B = [...back.atoms.values()];
    expect(B.length).toBe(A.length);
    const dx = B[0].x - A[0].x, dy = B[0].y - A[0].y;
    A.forEach((a, i) => {
      const b = B[i];
      expect([b.el, b.charge, b.isotope, b.abbrev, b.alias]).toEqual([a.el, a.charge, a.isotope, a.abbrev, a.alias]);
      expect(b.x - dx).toBeCloseTo(a.x, 2);
      expect(b.y - dy).toBeCloseTo(a.y, 2);
    });
    const idx = (d: ChemDoc, id: number) => [...d.atoms.keys()].indexOf(id);
    const bondSig = (d: ChemDoc) => [...d.bonds.values()].map((b) => [idx(d, b.a), idx(d, b.b), b.order, b.style, b.dbPos ?? 'auto']);
    expect(bondSig(back)).toEqual(bondSig(doc));
    // texts
    const texts = [...back.texts.values()];
    expect(texts.map((t) => t.text).sort()).toEqual(['CH_3CO_2H yields *product* **1**', 'Ph3P']);
    const ph = texts.find((t) => t.formula)!;
    expect([ph.size, ph.align, ph.color]).toEqual([1.5, 'center', '#0000ff']);
    expect(ph.x - dx).toBeCloseTo(6, 2);
    // arrows (kinds, geometry, conditions)
    const arrows = [...back.arrows.values()];
    expect(arrows.map((a) => a.kind)).toEqual([...doc.arrows.values()].map((a) => a.kind));
    const rxn = arrows[0];
    expect([rxn.above, rxn.below]).toEqual(['H_2O', 'heat, 2 h']);
    expect(rxn.x1 - dx).toBeCloseTo(6, 2);
    expect(rxn.y2 - dy).toBeCloseTo(2, 2);
    // curved arrows: anchors become points, an exact arc survives unchanged
    const curved = [...back.curved.values()];
    expect(curved.length).toBe(4);
    const orig = [...doc.curved.values()];
    for (const k of [0, 2, 3]) {
      for (const key of ['c1', 'c2'] as const) {
        expect(curved[k][key].t, `curved ${k} ${key}.t`).toBeCloseTo(orig[k][key].t, 2);
        expect(curved[k][key].h, `curved ${k} ${key}.h`).toBeCloseTo(orig[k][key].h, 2);
      }
    }
    const c0 = orig[0];
    expect(curved[0].electrons).toBe(2);
    expect(curved[0].c1.t).toBeCloseTo(c0.c1.t, 2);
    expect(curved[0].c1.h).toBeCloseTo(c0.c1.h, 2);
    expect(curved[0].c2.t).toBeCloseTo(c0.c2.t, 2);
    expect(curved[0].c2.h).toBeCloseTo(c0.c2.h, 2);
    const nAtom = A[1];
    if (curved[0].from.type !== 'point') throw new Error('expected point anchor');
    expect(curved[0].from.x - dx).toBeCloseTo(nAtom.x, 2);
    expect(curved[1].electrons).toBe(1);
    expect(Math.sign(curved[1].c1.h)).toBe(-1);
    // shapes
    const shapes = [...back.shapes.values()];
    expect(shapes.map((s) => s.kind)).toEqual(['rect', 'ellipse', 'paren']);
    expect(shapes[1].dashed).toBe(true);
    expect(shapes[1].x2 - dx).toBeCloseTo(13, 2);
    expect(shapes[1].y1 - dy).toBeCloseTo(4, 2);
    expect(back.meta.title).toBe('round trip');
  });
});

// ───────────── CML ─────────────

describe('CML', () => {
  it('round-trips a molecule', () => {
    for (const m of [alanine(), decorated()]) {
      m.name = 'test';
      const xml = writeCML(m);
      const [back] = readCML(xml);
      expect(back.name).toBe('test');
      // CML only knows wedge/hatch bond stereo
      for (const b of m.bonds) if (b.style !== 'wedge' && b.style !== 'hash') b.style = 'plain';
      expectSameMol(m, back);
    }
    const s = parseSmiles('[NH4+].c1ccc2[nH]ccc2c1');
    expect(canon(readCML(writeCML(s))[0])).toBe(writeSmiles(s));
  });

  it('reads several molecules, array style and aromatic bonds', () => {
    const xml = `<?xml version="1.0"?>
<cml xmlns="http://www.xml-cml.org/schema" xmlns:cml="http://www.xml-cml.org/schema">
  <molecule id="m1" title="methanol">
    <atomArray>
      <atom id="a1" elementType="C" x2="0" y2="0" hydrogenCount="3"/>
      <atom id="a2" elementType="O" x2="1.5" y2="0" hydrogenCount="1"/>
    </atomArray>
    <bondArray><bond atomRefs2="a1 a2" order="1"/></bondArray>
  </molecule>
  <molecule id="m2">
    <name>benzene</name>
    <atomArray atomID="a1 a2 a3 a4 a5 a6" elementType="C C C C C C"
      x2="0 1.3 1.3 0 -1.3 -1.3" y2="1.5 0.75 -0.75 -1.5 -0.75 0.75"/>
    <bondArray atomRef1="a1 a2 a3 a4 a5 a6" atomRef2="a2 a3 a4 a5 a6 a1" order="A A A A A A"/>
  </molecule>
  <cml:molecule id="m3">
    <cml:atomArray>
      <cml:atom id="a1" elementType="N" formalCharge="1" hydrogenCount="4"/>
    </cml:atomArray>
  </cml:molecule>
</cml>`;
    const mols = readCML(xml);
    expect(mols.length).toBe(3);
    expect(mols[0].name).toBe('methanol');
    expect(canon(mols[0])).toBe('CO');
    expect(mols[1].name).toBe('benzene');
    expect(mols[1].bonds.map((b) => b.order).sort()).toEqual([1, 1, 1, 2, 2, 2]);
    expect(canon(mols[2])).toBe('[NH4+]');
  });
});

// ───────────── XYZ ─────────────

describe('XYZ', () => {
  it('perceives bonds in water and methane', () => {
    const water = readXYZ(`3
water
O    0.000000    0.000000    0.117300
H    0.000000    0.757200   -0.469200
H    0.000000   -0.757200   -0.469200
`);
    expect(water.name).toBe('water');
    expect(water.props.dim).toBe('3D');
    expect(water.bonds.length).toBe(2);
    expect(water.atoms[1].y).toBeCloseTo(-0.7572, 6); // y flipped into the y-down frame
    expect(water.atoms[1].z).toBeCloseTo(-0.4692, 6);
    expect(computeFormula(water).formula).toBe('H2O');
    const methane = readXYZ(`5
methane
6   0.0000   0.0000   0.0000
1   0.6291   0.6291   0.6291
1  -0.6291  -0.6291   0.6291
1  -0.6291   0.6291  -0.6291
1   0.6291  -0.6291  -0.6291
`);
    expect(methane.atoms.map((a) => a.el)).toEqual(['C', 'H', 'H', 'H', 'H']);
    expect(methane.bonds.length).toBe(4);
    expect(methane.bonds.every((b) => b.a === 0 || b.b === 0)).toBe(true);
  });

  it('perceives benzene with alternating bond orders', () => {
    const lines = ['12', 'benzene'];
    for (let k = 0; k < 6; k++) {
      const a = (Math.PI / 3) * k;
      lines.push(`C ${(1.39 * Math.cos(a)).toFixed(4)} ${(1.39 * Math.sin(a)).toFixed(4)} 0.0`);
      lines.push(`H ${(2.47 * Math.cos(a)).toFixed(4)} ${(2.47 * Math.sin(a)).toFixed(4)} 0.0`);
    }
    // a flat molecule: z has no spread but XYZ is always 3D
    const m = readXYZ(lines.join('\n'));
    expect(m.bonds.length).toBe(12);
    const cc = m.bonds.filter((b) => m.atoms[b.a].el === 'C' && m.atoms[b.b].el === 'C');
    expect(cc.map((b) => b.order).sort()).toEqual([1, 1, 1, 2, 2, 2]);
    expect(canon(m)).toBe('c1ccccc1');
    expect(computeFormula(m).formula).toBe('C6H6');
    // writer/reader round trip
    const again = readXYZ(writeXYZ(m));
    expect(canon(again)).toBe('c1ccccc1');
    expect(again.atoms[0].x).toBeCloseTo(m.atoms[0].x, 4);
  });

  it('assigns C=O and C≡N bonds', () => {
    const acn = readXYZ(`6
acetonitrile
C  0.000  0.000  0.000
C  1.460  0.000  0.000
N  2.616  0.000  0.000
H -0.390  1.028  0.000
H -0.390 -0.514  0.890
H -0.390 -0.514 -0.890
`);
    expect(canon(acn)).toBe(smi('CC#N'));
    const co2 = readXYZ('3\n\nO -1.16 0 0\nC 0 0 0\nO 1.16 0 0\n');
    expect(canon(co2)).toBe(smi('O=C=O'));
  });
});

// ───────────── detection ─────────────

describe('detectFormat', () => {
  it('uses the file extension', () => {
    expect(detectFormat('a.mol', '')).toBe('mol');
    expect(detectFormat('A.SDF', '')).toBe('sdf');
    expect(detectFormat('x.sd', '')).toBe('sdf');
    expect(detectFormat('r.rxn', '')).toBe('rxn');
    expect(detectFormat('d.cdxml', '')).toBe('cdxml');
    expect(detectFormat('m.cml', '')).toBe('cml');
    expect(detectFormat('w.xyz', '')).toBe('xyz');
    expect(detectFormat('s.smi', '')).toBe('smiles');
    expect(detectFormat('s.smiles', '')).toBe('smiles');
    expect(detectFormat('doc.cwj', '')).toBe('chemwrite');
    expect(detectFormat('doc.chemwrite', '')).toBe('chemwrite');
    expect(detectFormat('doc.json', '{"format":"chemwrite","version":1}')).toBe('chemwrite');
    expect(detectFormat('data.json', '{"foo":1}')).toBe('unknown');
    expect(detectFormat('multi.mol', ASPIRIN_SDF + ASPIRIN_SDF)).toBe('sdf');
  });

  it('sniffs content', () => {
    expect(detectFormat(null, writeMolfile(alanine()))).toBe('mol');
    expect(detectFormat(null, writeMolfile(alanine(), { v3000: true }))).toBe('mol');
    expect(detectFormat(null, ASPIRIN_SDF)).toBe('sdf');
    expect(detectFormat(null, writeRxn({ reactants: [alanine()], products: [alanine()] }))).toBe('rxn');
    expect(detectFormat(null, CDXML_FIXTURE)).toBe('cdxml');
    expect(detectFormat(null, writeCML(alanine()))).toBe('cml');
    expect(detectFormat(null, '<molecule><atomArray/></molecule>')).toBe('cml');
    expect(detectFormat(null, '3\nwater\nO 0 0 0\nH 0 0.75 0.5\nH 0 -0.75 0.5\n')).toBe('xyz');
    expect(detectFormat('clip.txt', 'CC(=O)Oc1ccccc1C(=O)O aspirin')).toBe('smiles');
    expect(detectFormat(null, 'CCO>>CC=O')).toBe('smiles');
    expect(detectFormat(null, '  {"format": "chemwrite", "atoms": []}')).toBe('chemwrite');
    expect(detectFormat(null, 'hello world, this is not chemistry')).toBe('unknown');
    expect(detectFormat(null, '')).toBe('unknown');
    expect(detectFormat(null, '<html><body/></html>')).toBe('unknown');
  });
});

// ───────────── edge cases ─────────────

describe('edge cases', () => {
  it('handles empty molecules and documents', () => {
    const empty = readMolfile(writeMolfile(new Mol()));
    expect(empty.atoms.length).toBe(0);
    const d = readCDXML(writeCDXML(createDoc()));
    expect(d.atoms.size + d.arrows.size + d.texts.size).toBe(0);
    expect(readCML(writeCML(parseSmiles('[He]')))[0].atoms[0].el).toBe('He');
  });

  it('reads files with CRLF line endings, a BOM and a missing header line', () => {
    const text = writeMolfile(alanine());
    const crlf = '\uFEFF' + text.replace(/\n/g, '\r\n');
    expectSameMol(alanine(), readMolfile(crlf));
    const noHeader = text.split('\n').slice(1).join('\n'); // title line lost
    expect(readMolfile(noHeader).atoms.length).toBe(6);
  });

  it('reads whitespace-separated (non fixed-width) atom and bond lines', () => {
    const text = `loose\n  prog\n\n2 1 0 0 0 0 0 0 0 0999 V2000\n0.0 0.0 0.0 C 0 0\n1.5 0.0 0.0 O 0 5\n1 2 1 0\nM  END\n`;
    const m = readMolfile(text);
    expect(m.atoms.map((a) => a.el)).toEqual(['C', 'O']);
    expect(m.atoms[1].charge).toBe(-1);
    expect(m.atoms[1].x).toBeCloseTo(1, 6);
  });
});

// ───────────── robustness ─────────────

describe('robustness', () => {
  const garbage = ['', 'garbage', 'hello\nworld\n', '\u0000\u0001\u0002', '<<<>>>', '{"a":1}', '$RXN\n', '3\nx\nO 0 0\n', '  99999  1  0  0  0  0  0  0  0  0999 V2000'];

  it('throws descriptive errors on garbage', () => {
    const readers: [string, (s: string) => unknown][] = [
      ['molfile', readMolfile], ['sdf', readSDF], ['rxn', readRxn], ['cdxml', readCDXML], ['cml', readCML], ['xyz', readXYZ],
    ];
    for (const [name, read] of readers) {
      for (const g of garbage) {
        let err: unknown = null;
        try {
          read(g);
        } catch (e) {
          err = e;
        }
        expect(err, `${name} should reject ${JSON.stringify(g)}`).toBeInstanceOf(FormatError);
        expect((err as Error).message.length).toBeGreaterThan(8);
      }
    }
    expect(() => readMolfile(CDXML_FIXTURE)).toThrow(/XML/);
    expect(() => readMolfile('$RXN\n\n\n\n  1  1\n')).toThrow(/readRxn/);
    expect(() => readCDXML('<cml><molecule/></cml>')).toThrow(/expected <CDXML>/);
    expect(() => readMolfile('x\n\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0 0.0 0.0 C\n')).toThrow(/truncated/);
  });

  it('never hangs or crashes on truncated or mangled input', () => {
    const inputs = [ASPIRIN_SDF, CDXML_FIXTURE, writeMolfile(decorated(), { v3000: true }), writeCML(alanine()), writeRxn({ reactants: [alanine()], products: [decorated()] })];
    const readers = [readMolfile, readSDF, readRxn, readCDXML, readCML, readXYZ];
    const t0 = Date.now();
    let seed = 12345;
    const rnd = () => ((seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff);
    for (const input of inputs) {
      for (let k = 0; k < 40; k++) {
        const cut = Math.floor(rnd() * input.length);
        let s = input.slice(0, cut);
        if (k % 2) {
          const p = Math.floor(rnd() * input.length);
          s = input.slice(0, p) + String.fromCharCode(32 + Math.floor(rnd() * 90)) + input.slice(p + 1);
        }
        for (const read of readers) {
          try {
            read(s);
          } catch (e) {
            expect(e).toBeInstanceOf(FormatError);
          }
        }
      }
    }
    expect(Date.now() - t0).toBeLessThan(20000);
  });
});
