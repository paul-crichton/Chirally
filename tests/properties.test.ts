import { describe, it, expect } from 'vitest';
import { parseSmiles } from '../src/chem/smiles';
import { computeProperties, crippenContribs, detectFunctionalGroups } from '../src/chem/properties';

const props = (smi: string) => computeProperties(parseSmiles(smi));

// Reference values from RDKit 2026.03 (Crippen.MolLogP / MolMR, rdMolDescriptors.CalcTPSA,
// CalcNumHBD/HBA, CalcNumRotatableBonds (default = Strict), CalcNumLipinskiHBD/HBA, CalcNumRings,
// CalcNumAromaticRings, CalcFractionCSP3).
// [name, SMILES, logP, MR, TPSA, HBD, HBA, rotB, lipinskiHBD, lipinskiHBA, rings, aromRings, Fsp3]
type Row = [string, string, number, number, number, number, number, number, number, number, number, number, number];
const RDKIT: Row[] = [
  // cosmetic / skincare ingredients
  ['glycerol', 'OCC(O)CO', -1.6681, 20.1784, 60.69, 3, 3, 2, 3, 3, 0, 0, 1.0],
  ['niacinamide', 'NC(=O)c1cccnc1', 0.1805, 32.7549, 55.98, 1, 2, 1, 2, 3, 1, 1, 0.0],
  ['retinol', 'CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/CO)/C)/C', 5.5103, 93.2118, 20.23, 1, 1, 5, 1, 1, 1, 0, 0.5],
  ['salicylic acid', 'OC(=O)c1ccccc1O', 1.0904, 35.0661, 57.53, 2, 2, 1, 2, 3, 1, 1, 0.0],
  ['phenoxyethanol', 'OCCOc1ccccc1', 1.0577, 39.0228, 29.46, 1, 2, 3, 1, 2, 1, 1, 0.25],
  ['methylparaben', 'COC(=O)c1ccc(O)cc1', 1.1788, 39.4463, 46.53, 1, 3, 1, 1, 3, 1, 1, 0.125],
  ['caffeine', 'Cn1cnc2c1c(=O)n(C)c(=O)n2C', -1.0293, 51.196, 61.82, 0, 3, 0, 0, 6, 2, 2, 0.375],
  ['squalane', 'CC(C)CCCC(C)CCCC(C)CCCCC(C)CCCC(C)CCCC(C)C', 11.0844, 140.204, 0.0, 0, 0, 21, 0, 0, 0, 0, 1.0],
  ['sodium lauryl sulfate', 'CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]', 0.3881, 67.4932, 66.43, 0, 4, 12, 0, 4, 0, 0, 1.0],
  ['oxybenzone', 'COc1ccc(C(=O)c2ccccc2)c(O)c1', 2.6318, 64.5333, 46.53, 1, 3, 3, 1, 3, 2, 2, 0.0714],
  ['avobenzone', 'CC(C)(C)c1ccc(cc1)C(=O)CC(=O)c1ccc(OC)cc1', 4.4484, 91.414, 43.37, 0, 3, 5, 0, 3, 2, 2, 0.3],
  ['panthenol', 'CC(C)(CO)C(O)C(=O)NCCCO', -1.1356, 51.5861, 89.79, 4, 4, 6, 4, 5, 0, 0, 0.8889],
  ['α-tocopherol', 'Cc1c(C)c2O[C@](C)(CCC[C@H](C)CCC[C@H](C)CCCC(C)C)CCc2c(C)c1O', 8.8403, 134.3908, 29.46, 1, 2, 12, 1, 2, 2, 1, 0.7931],
  ['hyaluronan disaccharide', 'CC(=O)NC1C(O)OC(CO)C(O)C1OC1OC(C(O)C(O)C1O)C(=O)O', -5.1611, 80.7673, 215.47, 8, 11, 5, 8, 13, 2, 0, 0.8571],
  ['cetyl alcohol', 'CCCCCCCCCCCCCCCCO', 5.46, 77.3978, 20.23, 1, 1, 14, 1, 1, 0, 0, 1.0],
  ['citric acid', 'OC(=O)CC(O)(CC(=O)O)C(=O)O', -1.2485, 37.0912, 132.13, 4, 4, 5, 4, 7, 0, 0, 0.5],
  ['lactic acid', 'OC(=O)C(O)C', -0.5482, 19.3166, 57.53, 2, 2, 1, 2, 3, 0, 0, 0.6667],
  ['ascorbic acid', 'O=C1OC(C(O)=C1O)C(O)CO', -1.4074, 35.2562, 107.22, 4, 6, 2, 4, 6, 1, 0, 0.5],
  ['resveratrol', 'Oc1ccc(cc1)/C=C/c1cc(O)cc(O)c1', 2.9738, 66.8064, 60.69, 3, 3, 2, 3, 3, 2, 2, 0.0],
  ['ferulic acid', 'COc1cc(/C=C/C(=O)O)ccc1O', 1.4986, 51.3286, 66.76, 2, 3, 3, 2, 4, 1, 1, 0.1],
  ['octinoxate', 'CCCCC(CC)COC(=O)/C=C/c1ccc(OC)cc1', 4.468, 86.293, 35.53, 0, 3, 9, 0, 3, 1, 1, 0.5],
  ['octocrylene', 'CCCCC(CC)COC(=O)C(C#N)=C(c1ccccc1)c1ccccc1', 5.7717, 108.784, 50.09, 0, 3, 9, 0, 3, 2, 2, 0.3333],
  ['EDTA', 'C(C(=O)O)N(CC(=O)O)CCN(CC(=O)O)CC(=O)O', -2.0712, 63.4232, 155.68, 4, 6, 11, 4, 10, 0, 0, 0.6],
  ['triethanolamine', 'OCCN(CCO)CCO', -1.7347, 37.6974, 63.93, 3, 4, 6, 3, 4, 0, 0, 1.0],
  ['betaine', 'C[N+](C)(C)CC(=O)[O-]', -1.5575, 27.9064, 40.13, 0, 2, 2, 0, 3, 0, 0, 0.8],
  ['cetrimonium', 'C[N+](C)(C)CCCCCCCCCCCCCCCC', 6.1739, 93.2114, 0.0, 0, 0, 15, 0, 1, 0, 0, 1.0],
  ['isopropyl myristate', 'CCCCCCCCCCCCCC(=O)OC(C)C', 5.6391, 82.306, 26.3, 0, 2, 13, 0, 2, 0, 0, 0.9412],
  ['cinnamal', 'O=C/C=C/c1ccccc1', 1.8987, 41.54, 17.07, 0, 1, 2, 0, 1, 1, 1, 0.0],
  ['coumarin', 'O=c1ccc2ccccc2o1', 1.793, 42.484, 30.21, 0, 2, 0, 0, 2, 2, 2, 0.0],
  ['sorbitan monolaurate', 'CCCCCCCCCCCC(=O)OC[C@@H](O)[C@H]1OC[C@H](O)[C@H]1O', 1.932, 90.5634, 96.22, 3, 6, 13, 3, 6, 1, 0, 0.9444],
  // drugs and reference molecules
  ['benzene', 'c1ccccc1', 1.6866, 26.442, 0.0, 0, 0, 0, 0, 0, 1, 1, 0.0],
  ['ethanol', 'CCO', -0.0014, 12.7598, 20.23, 1, 1, 0, 1, 1, 0, 0, 1.0],
  ['aspirin', 'CC(=O)Oc1ccccc1C(=O)O', 1.3101, 44.7103, 63.6, 1, 3, 2, 1, 4, 1, 1, 0.1111],
  ['ibuprofen', 'CC(C)Cc1ccc(cc1)C(C)C(=O)O', 3.0732, 61.0348, 37.3, 1, 1, 4, 1, 2, 1, 1, 0.4615],
  ['phenol', 'Oc1ccccc1', 1.3922, 28.1068, 20.23, 1, 1, 0, 1, 1, 1, 1, 0.0],
  ['acetic acid', 'CC(=O)O', 0.0909, 13.3098, 37.3, 1, 1, 0, 1, 2, 0, 0, 0.5],
  ['glycine', 'NCC(=O)O', -0.9703, 16.6902, 63.32, 2, 2, 1, 3, 3, 0, 0, 0.5],
  ['paracetamol', 'CC(=O)Nc1ccc(O)cc1', 1.3506, 42.4105, 49.33, 2, 2, 1, 2, 3, 1, 1, 0.125],
  ['arginine (zwitterion)', 'NC(=[NH2+])NCCC[C@H]([NH3+])C(=O)[O-]', -5.2096, 40.0396, 131.41, 4, 2, 5, 8, 6, 0, 0, 0.6667],
  ['AMP', 'OP(O)(=O)OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1O', -1.863, 73.6551, 186.07, 5, 9, 4, 6, 12, 3, 2, 0.5],
  ['nitrobenzene', '[O-][N+](=O)c1ccccc1', 1.5948, 33.0964, 43.14, 0, 2, 1, 0, 3, 1, 1, 0.0],
  ['sulfanilamide', 'Nc1ccc(cc1)S(N)(=O)=O', -0.0838, 42.2276, 86.18, 2, 3, 1, 4, 4, 1, 1, 0.0],
  ['pyridine N-oxide', '[O-][n+]1ccccc1', 0.32, 25.402, 26.94, 0, 1, 0, 0, 2, 1, 1, 0.0],
  ['imidazole', 'c1cnc[nH]1', 0.4097, 18.5877, 28.68, 1, 1, 0, 1, 2, 1, 1, 0.0],
  ['uracil', 'O=c1cc[nH]c(=O)[nH]1', -0.9368, 27.6834, 65.72, 2, 2, 0, 2, 4, 1, 1, 0.0],
  ['indole', 'c1ccc2[nH]ccc2c1', 2.1679, 38.2987, 15.79, 1, 0, 0, 1, 1, 2, 2, 0.0],
  ['methyl parathion', 'COP(=S)(OC)Oc1ccc([N+](=O)[O-])cc1', 2.4909, 62.0274, 70.83, 0, 6, 5, 0, 6, 1, 1, 0.25],
  ['hydrocortisone acetate', 'CC(=O)OCC(=O)[C@@]1(O)CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@H]3[C@@H](O)C[C@@]21C', 2.3524, 104.6896, 100.9, 2, 6, 3, 2, 6, 4, 0, 0.7826],
  ['penicillin G', 'CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O', 0.8608, 85.8045, 86.71, 2, 4, 4, 2, 6, 3, 1, 0.4375],
  ['cetirizine', 'Clc1ccc(cc1)C(c1ccccc1)N1CCN(CC1)CCOCC(=O)O', 3.1482, 106.2048, 53.01, 1, 4, 8, 1, 5, 3, 2, 0.381],
  ['diazepam', 'CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21', 3.1538, 81.81, 32.67, 0, 2, 1, 0, 3, 3, 2, 0.125],
  ['flufenamic acid', 'OC(=O)c1ccccc1Nc1cccc(c1)C(F)(F)F', 4.1472, 68.128, 49.33, 2, 2, 3, 2, 3, 2, 2, 0.0714],
];

describe('RDKit-compatible descriptors', () => {
  for (const [name, smi, logP, mr, tpsa, hbd, hba, rb, lhbd, lhba, rings, arom, fsp3] of RDKIT) {
    it(name, () => {
      const p = props(smi);
      expect(p.logP).toBeCloseTo(logP, 3);
      expect(Math.abs(p.mr - mr)).toBeLessThan(0.005);
      expect(Math.abs(p.tpsa - tpsa)).toBeLessThan(0.006);
      expect([p.hbd, p.hba, p.rotatableBonds]).toEqual([hbd, hba, rb]);
      expect([p.lipinskiHBD, p.lipinskiHBA]).toEqual([lhbd, lhba]);
      expect([p.rings, p.aromaticRings]).toEqual([rings, arom]);
      expect(p.fractionCsp3).toBeCloseTo(fsp3, 3);
    });
  }

  it('formula, masses, heavy atoms and charge', () => {
    const p = props('CC(=O)Oc1ccccc1C(=O)O');
    expect(p.formula).toBe('C9H8O4');
    expect(p.mw).toBeCloseTo(180.159, 2);
    expect(p.exactMass).toBeCloseTo(180.0423, 3);
    expect(p.heavyAtoms).toBe(13);
    expect(props('CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]').charge).toBe(0);
    expect(props('C[N+](C)(C)CCCCCCCCCCCCCCCC').charge).toBe(1);
    expect(props('CC(C)CCCC(C)CCCC(C)CCCCC(C)CCCC(C)CCCC(C)C').formula).toBe('C30H62');
  });

  it('Lipinski and Veber rules', () => {
    expect(props('CC(=O)Oc1ccccc1C(=O)O').lipinski).toEqual({ violations: 0, pass: true });
    const sq = props('CC(C)CCCC(C)CCCC(C)CCCCC(C)CCCC(C)CCCC(C)C'); // logP 11 → 1 violation, RB 21
    expect(sq.lipinski).toEqual({ violations: 1, pass: true });
    expect(sq.veber.pass).toBe(false);
    const ha = props('CC(=O)NC1C(O)OC(CO)C(O)C1OC1OC(C(O)C(O)C1O)C(=O)O'); // TPSA 215, N+O = 13
    expect(ha.veber.pass).toBe(false);
    expect(ha.lipinski.violations).toBe(2); // HBD 8 > 5, N+O 13 > 10
  });

  it('stereocentre counts from constitution', () => {
    expect(props('CC(C)Cc1ccc(cc1)C(C)C(=O)O').stereocenters).toBe(1);
    expect(props('CC(C)C1CCC(C)CC1O').stereocenters).toBe(3); // menthol
    expect(props('OCC1OC(O)C(O)C(O)C1O').stereocenters).toBe(5); // glucopyranose
    expect(props('CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C').stereocenters).toBe(8); // cholesterol
    expect(props('OCC(O)CO').stereocenters).toBe(0);
    expect(props('OC(=O)CC(O)(CC(=O)O)C(=O)O').stereocenters).toBe(0);
    expect(props('CS(=O)c1ccc(C)cc1').stereocenters).toBe(1); // sulfoxide with lone pair
    expect(props('CS(C)=O').stereocenters).toBe(0);
    // Kekulé alternation must not split symmetry classes (diphenylmethanol has no stereocentre)
    expect(props('OC(C1=CC=CC=C1)C1=CC=CC=C1').stereocenters).toBe(0);
  });
});

describe('formulation descriptors', () => {
  it('ESOL logS follows Delaney 2004 with Crippen logP, rotatable bonds and aromatic proportion', () => {
    const p = props('CC(=O)Oc1ccccc1C(=O)O');
    // aspirin: AP = 6 aromatic / 13 heavy atoms
    const expected = 0.16 - 0.63 * p.logP - 0.0062 * p.mw + 0.066 * 2 - 0.74 * (6 / 13);
    expect(p.logS).toBeCloseTo(expected, 10);
    expect(p.logS).toBeCloseTo(-1.992, 2);
    const g = props('OCC(O)CO'); // no aromatic atoms: 0.16 − 0.63(−1.6681) − 0.0062(92.094) + 0.066·2
    expect(g.logS).toBeCloseTo(0.16 + 0.63 * 1.6681 - 0.0062 * 92.094 + 0.132, 3);
  });
  it('Potts–Guy skin permeability log Kp = −2.72 + 0.71·logP − 0.0061·MW', () => {
    const c = props('Cn1cnc2c1c(=O)n(C)c(=O)n2C'); // caffeine: −2.72 + 0.71(−1.0293) − 0.0061(194.19)
    expect(c.logKp).toBeCloseTo(-2.72 + 0.71 * -1.0293 - 0.0061 * 194.194, 3);
    expect(c.logKp).toBeCloseTo(-4.635, 2);
    const m = props('COC(=O)c1ccc(O)cc1'); // methylparaben: −2.72 + 0.71(1.1788) − 0.0061(152.15)
    expect(m.logKp).toBeCloseTo(-2.72 + 0.71 * 1.1788 - 0.0061 * 152.149, 3);
  });
  it('skin-penetration guidance: 500 Da rule and logP 1–3 window', () => {
    const nia = props('NC(=O)c1cccnc1');
    expect(nia.skinPenetration.under500Da).toBe(true);
    expect(nia.skinPenetration.logPInRange).toBe(false);
    expect(nia.skinPenetration.note).toMatch(/hydrophilic/);
    expect(props('COC(=O)c1ccc(O)cc1').skinPenetration.logPInRange).toBe(true);
    const ret = props('CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/CO)/C)/C');
    expect(ret.skinPenetration.logPInRange).toBe(false);
    expect(ret.skinPenetration.note).toMatch(/lipophilic/);
    const ha = props('CC(=O)NC1C(O)OC(CO)C(O)C1OC1OC(C(O)C(O)C1O)C(=O)O.CC(=O)NC1C(O)OC(CO)C(O)C1OC1OC(C(O)C(O)C1O)C(=O)O');
    expect(ha.skinPenetration.under500Da).toBe(false);
  });
});

describe('HLB (Davies group contributions and Griffin)', () => {
  it('sodium dodecyl sulfate: Davies 7 + 38.7 − 12×0.475 = 40.0; Griffin n/a (ionic)', () => {
    const h = props('CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]').hlb!;
    expect(h.davies).toBeCloseTo(40.0, 3);
    expect(h.griffin).toBeNull();
  });
  it('glyceryl monostearate: Davies 7 + 2.4 + 2×1.9 − 20×0.475 = 3.7; Griffin 20·75.09/358.56 = 4.19', () => {
    // Griffin hydrophilic part = glyceryl C3H5 + 2 OH (ester O stays with the acyl group), i.e.
    // Griffin's 20(1 − S/A) form; the tabulated experimental value is 3.8. (Excluding the three
    // glycerol CH groups from Davies' −0.475 count would give ≈ 5.1 instead.)
    const h = props('CCCCCCCCCCCCCCCCCC(=O)OCC(O)CO').hlb!;
    expect(h.davies).toBeCloseTo(3.7, 3);
    expect(h.griffin!).toBeCloseTo((20 * (3 * 12.0107 + 5 * 1.00794 + 2 * 17.00734)) / 358.563, 1);
    expect(h.griffin!).toBeCloseTo(4.19, 2);
  });
  it('C12E4 (Brij-30-like): Griffin 20·(4×44.05)/362.55 = 9.72; Davies 7 + 4×0.33 + 1.9 − 12×0.475 = 4.52', () => {
    const h = props('CCCCCCCCCCCCOCCOCCOCCOCCO').hlb!;
    expect(h.griffin!).toBeCloseTo(9.72, 2);
    expect(h.davies!).toBeCloseTo(4.52, 3);
    expect(h.note).toMatch(/4 EO units/);
  });
  it('cetyl alcohol: Davies 7 + 1.9 − 16×0.475 = 1.3; Griffin 20·17.01/242.45 = 1.40', () => {
    const h = props('CCCCCCCCCCCCCCCCO').hlb!;
    expect(h.davies!).toBeCloseTo(1.3, 3);
    expect(h.griffin!).toBeCloseTo(1.4, 2);
  });
  it('sorbitan monolaurate (Span 20): Davies 7 + 6.8 + 3×0.5 + 1.3 − 17×0.475 = 8.525; Griffin ≈ 8.5', () => {
    // Griffin Mh = sorbitan C6H8 + 3 OH + ring O = 147.15; M = 346.46 → 8.49 (literature 8.6)
    const h = props('CCCCCCCCCCCC(=O)OC[C@@H](O)[C@H]1OC[C@H](O)[C@H]1O').hlb!;
    expect(h.davies!).toBeCloseTo(8.525, 3);
    expect(h.griffin!).toBeCloseTo(8.49, 2);
  });
  it('soaps: sodium stearate 7 + 19.1 − 17×0.475 = 18.025; potassium stearate 20.025', () => {
    expect(props('CCCCCCCCCCCCCCCCCC(=O)[O-].[Na+]').hlb!.davies!).toBeCloseTo(18.025, 3);
    expect(props('CCCCCCCCCCCCCCCCCC(=O)[O-].[K+]').hlb!.davies!).toBeCloseTo(20.025, 3);
  });
  it('mixed EO/PO: C12–PO–EO–OH: Davies 7 + 0.33 − 0.15 + 1.9 − 12×0.475 = 3.38; Griffin 20·44.05/288.47', () => {
    const h = props('CCCCCCCCCCCCOCC(C)OCCO').hlb!;
    expect(h.davies!).toBeCloseTo(3.38, 3);
    expect(h.griffin!).toBeCloseTo((20 * 44.053) / 288.472, 2);
  });
  it('ethoxylates: steareth-20 15.3, PEG-100 stearate 18.8, polysorbate 20 ≈ 15.9 (Griffin)', () => {
    // steareth-20 C18H37(OCH2CH2)20OH: 20 × (20 × 44.053) / 1151.5 = 15.30 (literature 15.3)
    const stP = props('CCCCCCCCCCCCCCCCCCO' + 'CCO'.repeat(20));
    expect(stP.mw).toBeCloseTo(1151.5, 0);
    expect(stP.hlb!.griffin!).toBeCloseTo((20 * 20 * 44.053) / stP.mw, 2);
    expect(stP.hlb!.griffin!).toBeCloseTo(15.3, 1);
    // PEG-100 stearate: 20 × (100 × 44.053) / 4689.7 = 18.79 (literature 18.8)
    expect(props('CCCCCCCCCCCCCCCCCC(=O)O' + 'CCO'.repeat(100)).hlb!.griffin!).toBeCloseTo(18.79, 2);
    // idealised polysorbate 20 (20 EO, laurate at the end of one chain): Mh = 20 EO (881.06) + sorbitan
    // C6H8 + ring O (96.13) → 20 × 977.19 / 1227.5 = 15.92 (literature 16.7)
    const tw = props('CCCCCCCCCCCC(=O)OCCOCCOCCOCCOCCOCC(OCCOCCOCCOCCOCCO)C1OCC(OCCOCCOCCOCCOCCO)C1OCCOCCOCCOCCOCCO').hlb!;
    expect(tw.griffin!).toBeCloseTo(15.92, 2);
    expect(tw.note).toMatch(/20 EO units/);
  });
  it('sodium laureth-2 sulfate: Davies 7 + 38.7 + 2×0.33 − 12×0.475 = 40.66', () => {
    expect(props('CCCCCCCCCCCCOCCOCCOS(=O)(=O)[O-].[Na+]').hlb!.davies!).toBeCloseTo(40.66, 3);
  });
  it('glycol distearate: the glycol is not an EO unit; Griffin 20 × C2H4 / 595.0 = 0.94, Davies 7 + 2×2.4 − 36×0.475 = −5.3', () => {
    const h = props('CCCCCCCCCCCCCCCCCC(=O)OCCOC(=O)CCCCCCCCCCCCCCCCC').hlb!;
    expect(h.griffin!).toBeCloseTo((20 * (2 * 12.0107 + 4 * 1.00794)) / 594.99, 2);
    expect(h.davies!).toBeCloseTo(-5.3, 3);
    expect(h.note).not.toMatch(/EO unit/);
  });
  it('fluorosurfactant: sodium perfluorooctanoate 7 + 19.1 − 7×0.870 = 20.01', () => {
    expect(props('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)[O-].[Na+]').hlb!.davies!).toBeCloseTo(20.01, 3);
  });
  it('only amphiphiles get an HLB', () => {
    for (const s of ['OCC(O)CO', 'Cn1cnc2c1c(=O)n(C)c(=O)n2C', 'CC(C)CCCC(C)CCCC(C)CCCCC(C)CCCC(C)CCCC(C)C', 'c1ccccc1', 'CCCCO'])
      expect(props(s).hlb, s).toBeNull();
    const ipm = props('CCCCCCCCCCCCCC(=O)OC(C)C').hlb!; // oily ester: Davies 7 + 2.4 − 16×0.475 = 1.8
    expect(ipm.davies!).toBeCloseTo(1.8, 3);
    expect(ipm.griffin).toBeNull();
  });
});

describe('functional groups', () => {
  const names = (smi: string) => Object.fromEntries(detectFunctionalGroups(parseSmiles(smi)).map((g) => [g.name, g.count]));
  it('alcohols, phenols, ethers', () => {
    expect(names('CCO')).toEqual({ 'primary alcohol': 1 });
    expect(names('CC(C)O')).toEqual({ 'secondary alcohol': 1 });
    expect(names('CC(C)(C)O')).toEqual({ 'tertiary alcohol': 1 });
    expect(names('Oc1ccccc1')).toEqual({ phenol: 1, arene: 1 });
    expect(names('COc1ccccc1')).toEqual({ ether: 1, arene: 1 });
    expect(names('OCC(O)CO')).toEqual({ 'primary alcohol': 2, 'secondary alcohol': 1 });
  });
  it('carbonyl compounds and acid derivatives', () => {
    expect(names('CC=O')).toEqual({ aldehyde: 1, 'skin sensitisation alert (aldehyde)': 1 });
    expect(names('CC(C)=O')).toEqual({ ketone: 1 });
    expect(names('CC(=O)O')).toEqual({ 'carboxylic acid': 1 });
    expect(names('CC(=O)[O-].[Na+]')).toEqual({ carboxylate: 1 });
    expect(names('CCOC(=O)C')).toEqual({ ester: 1 });
    expect(names('O=C1CCCCCO1')).toEqual({ lactone: 1 });
    expect(names('CC(N)=O')).toEqual({ amide: 1 });
    expect(names('O=C1CCCN1')).toEqual({ lactam: 1 });
    expect(names('NC(N)=O')).toEqual({ urea: 1 });
    expect(names('COC(N)=O')).toEqual({ carbamate: 1 });
    expect(names('CC(=O)OC(C)=O')).toEqual({ anhydride: 1 });
    expect(names('CC(=O)Cl')).toEqual({ 'acyl halide': 1 });
  });
  it('nitrogen groups', () => {
    expect(names('CN')).toEqual({ 'primary amine': 1 });
    expect(names('CNC')).toEqual({ 'secondary amine': 1 });
    expect(names('CN(C)C')).toEqual({ 'tertiary amine': 1 });
    expect(names('C[N+](C)(C)CCCCCCCCCCCCCCCC')).toEqual({ 'quaternary ammonium': 1 });
    expect(names('Nc1ccccc1')).toEqual({ aniline: 1, arene: 1 });
    expect(names('CC#N')).toEqual({ nitrile: 1 });
    expect(names('[O-][N+](=O)c1ccccc1')).toEqual({ nitro: 1, arene: 1 });
    expect(names('CC=NC')).toEqual({ imine: 1 });
    expect(names('CC=NO')).toEqual({ oxime: 1 });
    expect(names('NN')).toEqual({ hydrazine: 1 });
    expect(names('c1ccc(cc1)N=Nc1ccccc1')).toEqual({ azo: 1, arene: 2 });
    expect(names('c1ccncc1')).toEqual({ heteroarene: 1 });
  });
  it('halogen, sulfur and phosphorus groups', () => {
    expect(names('CCCl')).toEqual({ 'alkyl halide': 1 });
    expect(names('Clc1ccccc1')).toEqual({ 'aryl halide': 1, arene: 1 });
    expect(names('CCS')).toEqual({ thiol: 1 });
    expect(names('CCSCC')).toEqual({ sulfide: 1 });
    expect(names('CS(C)=O')).toEqual({ sulfoxide: 1 });
    expect(names('CS(C)(=O)=O')).toEqual({ sulfone: 1 });
    expect(names('CS(N)(=O)=O')).toEqual({ sulfonamide: 1 });
    expect(names('OS(=O)(=O)c1ccccc1')).toEqual({ 'sulfonic acid': 1, arene: 1 });
    expect(names('CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]')).toEqual({ 'sulfate ester': 1 });
    expect(names('CCOP(=O)(OCC)OCC')).toEqual({ phosphate: 1 });
  });
  it('unsaturation, small rings, acetals, peroxides', () => {
    expect(names('C=C')).toEqual({ alkene: 1 });
    expect(names('C#C')).toEqual({ alkyne: 1 });
    expect(names('C1CO1')).toEqual({ epoxide: 1 });
    expect(names('CC(OCC)OCC')).toEqual({ 'acetal/ketal': 1 });
    expect(names('OCC1OC(O)C(O)C(O)C1O')).toEqual({ 'primary alcohol': 1, 'secondary alcohol': 3, hemiacetal: 1 });
    expect(names('COC=C')).toEqual({ ether: 1, 'enol ether': 1, alkene: 1 });
    expect(names('CC(C)(C)OO')).toEqual({ peroxide: 1 });
  });
  it('sensitisation alerts: Michael acceptors and aldehydes (cinnamal, acrylates)', () => {
    expect(names('O=C/C=C/c1ccccc1')).toEqual({
      aldehyde: 1, alkene: 1, 'α,β-unsaturated carbonyl (Michael acceptor)': 1, 'skin sensitisation alert (aldehyde)': 1, arene: 1,
    });
    expect(names('C=CC(=O)OC')['α,β-unsaturated carbonyl (Michael acceptor)']).toBe(1);
    expect(names('C=CC#N')['α,β-unsaturated carbonyl (Michael acceptor)']).toBe(1);
  });
  it('reports atom indices of each occurrence', () => {
    const g = detectFunctionalGroups(parseSmiles('CCOC(=O)C'));
    expect(g).toEqual([{ name: 'ester', count: 1, atoms: [[2, 3, 4]] }]);
    const glyc = detectFunctionalGroups(parseSmiles('OCC(O)CO')).find((x) => x.name === 'primary alcohol')!;
    expect(glyc.atoms).toEqual([[0, 1], [4, 5]]);
  });
});

describe('per-atom contributions, abbreviations and pseudo atoms', () => {
  it('crippenContribs sum to logP / MR (implicit H folded into heavy atoms)', () => {
    for (const s of ['CC(=O)Oc1ccccc1C(=O)O', 'NC(=O)c1cccnc1', 'OCCOc1ccccc1']) {
      const m = parseSmiles(s);
      const c = crippenContribs(m);
      const p = computeProperties(m);
      expect(c.logP.length).toBe(m.atoms.length);
      expect(c.logP.reduce((a, b) => a + b, 0)).toBeCloseTo(p.logP, 10);
      expect(c.mr.reduce((a, b) => a + b, 0)).toBeCloseTo(p.mr, 10);
    }
    // ethanol: CH3 (C1 0.1441 + 3×H1 0.123), CH2 (C3 −0.2035 + 2×0.123), OH (O2 −0.2893 + H2 −0.2677)
    const e = crippenContribs(parseSmiles('CCO')).logP;
    expect(e[0]).toBeCloseTo(0.1441 + 3 * 0.123, 4);
    expect(e[1]).toBeCloseTo(-0.2035 + 2 * 0.123, 4);
    expect(e[2]).toBeCloseTo(-0.2893 - 0.2677, 4);
  });
  it('expands abbreviations (OMe) and maps contributions/groups back to the label atom', () => {
    const m = parseSmiles('CC(=O)c1ccc(O)cc1');
    m.atoms[0].abbrev = 'OMe'; // → methylparaben
    const p = computeProperties(m);
    expect(p.formula).toBe('C8H8O3');
    expect(p.logP).toBeCloseTo(1.1788, 3);
    expect(p.tpsa).toBeCloseTo(46.53, 2);
    const c = crippenContribs(m);
    expect(c.logP.length).toBe(m.atoms.length);
    expect(c.logP.reduce((a, b) => a + b, 0)).toBeCloseTo(1.1788, 3);
    const ester = p.functionalGroups.find((g) => g.name === 'ester')!;
    expect(ester.atoms[0]).toEqual([0, 1, 2]); // O of OMe is the abbreviation atom 0
  });
  it('ignores pseudo atoms gracefully', () => {
    const p = props('*CCO');
    expect(p.heavyAtoms).toBe(3);
    expect(Number.isFinite(p.logP)).toBe(true);
    const r = parseSmiles('CC(=O)O');
    r.atoms[0].el = 'R';
    r.atoms[0].alias = 'R1';
    expect(() => computeProperties(r)).not.toThrow();
  });
});

describe('performance', () => {
  it('computes all properties of a ~100-atom molecule in < 20 ms', () => {
    // decapeptide-like chain, 99 heavy atoms
    const smi =
      'CC(C)C[C@H](NC(=O)[C@H](CC(C)C)NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](CO)NC(=O)[C@H](Cc1c[nH]c2ccccc12)NC(=O)[C@H](CCCCN)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@H](CCC(N)=O)NC(=O)[C@H](C)N)C(=O)N[C@@H](CCSC)C(=O)O';
    const m = parseSmiles(smi);
    expect(m.atoms.length).toBeGreaterThanOrEqual(99);
    computeProperties(m); // warm-up (pattern compilation, JIT)
    computeProperties(m);
    // best of N: robust against other test workers competing for the CPU (typical: 3–6 ms)
    let best = Infinity;
    for (let k = 0; k < 10; k++) {
      const t0 = performance.now();
      computeProperties(m);
      best = Math.min(best, performance.now() - t0);
    }
    expect(best).toBeLessThan(20);
  });
});
