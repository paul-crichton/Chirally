import { describe, it, expect } from 'vitest';
import { parseSmiles, writeSmiles, suppressHydrogens } from '../src/chem/smiles';
import { computeFormula, parseFormula } from '../src/chem/formula';
import { perceiveRings } from '../src/chem/rings';
import { perceiveAromaticity } from '../src/chem/aromaticity';
import { implicitH } from '../src/chem/valence';

const canon = (s: string) => writeSmiles(parseSmiles(s));

describe('SMILES parsing', () => {
  it('parses simple chains with implicit H', () => {
    const m = parseSmiles('CCO');
    expect(m.atoms.length).toBe(3);
    expect(implicitH(m, 0)).toBe(3);
    expect(implicitH(m, 1)).toBe(2);
    expect(implicitH(m, 2)).toBe(1);
  });
  it('kekulizes benzene', () => {
    const m = parseSmiles('c1ccccc1');
    const orders = m.bonds.map((b) => b.order).sort();
    expect(orders).toEqual([1, 1, 1, 2, 2, 2]);
    for (let i = 0; i < 6; i++) expect(implicitH(m, i)).toBe(1);
  });
  it('kekulizes pyrrole, furan, pyridine, indole, imidazole', () => {
    for (const s of ['c1cc[nH]c1', 'c1ccoc1', 'c1ccncc1', 'c1ccc2[nH]ccc2c1', 'c1c[nH]cn1', 'c1ccc2ccccc2c1', 'O=c1cccc[nH]1', 'c1ccsc1', 'C[n+]1ccccc1']) {
      const m = parseSmiles(s);
      expect(m.bonds.some((b) => b.order === 1.5), s).toBe(false);
    }
  });
  it('tolerates pyrrole written without explicit [nH]', () => {
    const m = parseSmiles('c1ccnc1');
    expect(m.bonds.some((b) => b.order === 1.5)).toBe(false);
    expect(computeFormula(m).formula).toBe('C4H5N');
  });
  it('parses charges, isotopes, brackets', () => {
    const m = parseSmiles('[NH4+].[Cl-]');
    expect(m.atoms[0].charge).toBe(1);
    expect(implicitH(m, 0)).toBe(4);
    expect(m.atoms[1].charge).toBe(-1);
    const d = parseSmiles('[2H]C([2H])([2H])O');
    expect(d.atoms[0].isotope).toBe(2);
    expect(computeFormula(d).formula).toBe('CHD3O');
  });
  it('parses ring closures with %nn and bond symbols', () => {
    const m = parseSmiles('C%10CCCCC%10');
    expect(m.bonds.length).toBe(6);
    const r = perceiveRings(m);
    expect(r.rings.length).toBe(1);
    expect(r.rings[0].length).toBe(6);
  });
  it('parses tetrahedral stereo', () => {
    const m = parseSmiles('N[C@@H](C)C(=O)O');
    expect(m.tetra.length).toBe(1);
    expect(m.tetra[0].center).toBe(1);
    expect(m.tetra[0].nbrs).toEqual([0, -1, 2, 3]);
    expect(m.tetra[0].ccw).toBe(false);
  });
  it('parses double-bond stereo', () => {
    const e = parseSmiles('F/C=C/F');
    expect(e.dbStereo.length).toBe(1);
    expect(e.dbStereo[0].cis).toBe(false);
    const z = parseSmiles('F/C=C\\F');
    expect(z.dbStereo[0].cis).toBe(true);
    const z2 = parseSmiles('C(/F)=C/F');
    expect(z2.dbStereo[0].cis).toBe(true);
  });
  it('rejects invalid input', () => {
    expect(() => parseSmiles('C1CC')).toThrow();
    expect(() => parseSmiles('C(C')).toThrow();
    expect(() => parseSmiles('[Xx]')).toThrow();
  });
});

describe('SMILES writing', () => {
  it('is canonical across input orderings', () => {
    expect(canon('OCC')).toBe(canon('CCO'));
    expect(canon('c1ccccc1O')).toBe(canon('Oc1ccccc1'));
    expect(canon('CC(=O)Oc1ccccc1C(=O)O')).toBe(canon('OC(=O)c1ccccc1OC(C)=O'));
    expect(canon('C1CCCCC1N')).toBe(canon('NC1CCCCC1'));
  });
  it('round-trips structures', () => {
    const list = [
      'CC(=O)Oc1ccccc1C(=O)O',
      'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
      'c1ccc2c(c1)[nH]c1ccccc12',
      'O=C([O-])C.[Na+]',
      'C#N',
      'C=C=C',
      'FC(F)(F)c1ccc(cc1)S(=O)(=O)N',
      'C1CC2CCC1C2',
      'OCCOCCOCCO',
      'CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]',
    ];
    for (const s of list) {
      const out = canon(s);
      expect(computeFormula(parseSmiles(out)).formula, s).toBe(computeFormula(parseSmiles(s)).formula);
      expect(canon(out), s).toBe(out);
    }
  });
  it('writes aromatic SMILES for benzene', () => {
    expect(canon('C1=CC=CC=C1')).toBe('c1ccccc1');
    expect(canon('c1cc[nH]c1')).toBe('c1cc[nH]c1');
  });
  it('preserves tetrahedral stereo through round trip', () => {
    const l = canon('N[C@@H](C)C(=O)O');
    const d = canon('N[C@H](C)C(=O)O');
    expect(l).not.toBe(d);
    expect(canon(l)).toBe(l);
    expect(canon('C[C@H](N)C(=O)O')).toBe(l);
  });
  it('preserves E/Z through round trip', () => {
    const e = canon('C/C=C/C');
    const z = canon('C/C=C\\C');
    expect(e).not.toBe(z);
    expect(canon(e)).toBe(e);
    expect(canon('C(/C)=C/C')).toBe(z);
  });
  it('suppresses explicit hydrogens', () => {
    const m = parseSmiles('[H]C([H])([H])O');
    const s = suppressHydrogens(m);
    expect(writeSmiles(s)).toBe('CO');
  });
});

describe('formula', () => {
  it('computes Hill formula and masses', () => {
    const f = computeFormula(parseSmiles('CC(=O)Oc1ccccc1C(=O)O'));
    expect(f.formula).toBe('C9H8O4');
    expect(f.mw).toBeCloseTo(180.157, 2);
    expect(f.exactMass).toBeCloseTo(180.0423, 3);
    expect(f.dbe).toBe(6);
  });
  it('handles charges and salts', () => {
    const f = computeFormula(parseSmiles('CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]'));
    expect(f.formula).toBe('C12H25NaO4S');
    expect(f.charge).toBe(0);
  });
  it('parses formula strings', () => {
    expect(parseFormula('C6H12O6').mw).toBeCloseTo(180.156, 2);
    expect(parseFormula('Ca(OH)2').formula).toBe('CaH2O2');
    expect(parseFormula('SO4^2-').charge).toBe(-2);
  });
});

describe('aromaticity', () => {
  it('perceives aromatic rings from Kekulé input', () => {
    const cases: [string, number][] = [
      ['C1=CC=CC=C1', 6], ['C1=CC=CN1', 5], ['C1=COC=C1', 5], ['C1=CC=C2C=CC=CC2=C1', 10],
      ['C1CCCCC1', 0], ['C1=CCCC=C1', 0], ['O=C1C=CC(=O)C=C1', 0], ['C1=CC=C2C=CC=C2C=C1', 10],
    ];
    for (const [s, n] of cases) {
      const m = parseSmiles(s);
      const a = perceiveAromaticity(m);
      expect(a.atoms.filter(Boolean).length, s).toBe(n);
    }
  });
});
