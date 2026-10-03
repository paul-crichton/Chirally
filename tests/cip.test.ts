import { describe, it, expect } from 'vitest';
import { parseSmiles } from '../src/chem/smiles';
import { assignCIP, _mancudeDuplicateZ } from '../src/chem/cip';
import { Mol } from '../src/chem/mol';
import { perceiveStereo2D } from '../src/chem/stereo2d';

// Expected descriptors were checked against PubChem IUPAC names (OpenEye Lexichem) where PubChem
// provides them, and against RDKit's new CIP labeller (rdCIPLabeler, Hanson et al. 2018) for r/s,
// rule-4b and Kekulé-variant cases. Atom indices follow the SMILES atom order.

/** Centre descriptors as { atomIndex: label }. */
function centers(smi: string): Record<number, string> {
  const r = assignCIP(parseSmiles(smi));
  return Object.fromEntries([...r.centers.entries()].sort((a, b) => a[0] - b[0]));
}

/** Double-bond descriptors as { 'i-j': label } with i < j the bond's atom indices. */
function bonds(smi: string): Record<string, string> {
  const m = parseSmiles(smi);
  const r = assignCIP(m);
  const out: Record<string, string> = {};
  for (const [bi, l] of r.bonds) {
    const b = m.bonds[bi];
    out[`${Math.min(b.a, b.b)}-${Math.max(b.a, b.b)}`] = l;
  }
  return out;
}

describe('CIP: amino acids and simple centres', () => {
  it('L-alanine is S, D-alanine is R', () => {
    expect(centers('N[C@@H](C)C(=O)O')).toEqual({ 1: 'S' });
    expect(centers('N[C@H](C)C(=O)O')).toEqual({ 1: 'R' });
  });
  it('L-cysteine is R (S outranks O-bearing carboxyl), L-serine is S', () => {
    expect(centers('N[C@@H](CS)C(=O)O')).toEqual({ 1: 'R' });
    expect(centers('N[C@@H](CO)C(=O)O')).toEqual({ 1: 'S' });
  });
  it('L-threonine (2S,3R) and L-allothreonine (2S,3S)', () => {
    // PubChem CID 6288 / 99289 isomeric SMILES
    expect(centers('C[C@H]([C@@H](C(=O)O)N)O')).toEqual({ 1: 'R', 2: 'S' });
    expect(centers('C[C@@H]([C@@H](C(=O)O)N)O')).toEqual({ 1: 'S', 2: 'S' });
  });
  it('L-isoleucine (2S,3S)', () => {
    expect(centers('N[C@@H]([C@@H](C)CC)C(=O)O')).toEqual({ 1: 'S', 2: 'S' });
  });
  it('L-proline, L-malic acid, (R)-mandelic acid, Ala-Ala', () => {
    expect(centers('OC(=O)[C@@H]1CCCN1')).toEqual({ 3: 'S' });
    expect(centers('O[C@@H](CC(=O)O)C(=O)O')).toEqual({ 1: 'S' });
    expect(centers('O=C(O)[C@H](O)c1ccccc1')).toEqual({ 3: 'R' });
    expect(centers('C[C@H](N)C(=O)N[C@@H](C)C(=O)O')).toEqual({ 1: 'S', 6: 'S' });
  });
  it('(R)- and (S)-butan-2-ol', () => {
    expect(centers('C[C@@H](O)CC')).toEqual({ 1: 'R' });
    expect(centers('C[C@H](O)CC')).toEqual({ 1: 'S' });
  });
  it('(R)-glyceraldehyde', () => {
    expect(centers('O=C[C@H](O)CO')).toEqual({ 2: 'R' });
  });
  it('CHFClBr / CFClBrI', () => {
    expect(centers('F[C@@H](Cl)Br')).toEqual({ 1: 'S' });
    expect(centers('F[C@](Cl)(Br)I')).toEqual({ 1: 'S' });
  });
  it('explicit hydrogen atoms give the same result as implicit ones', () => {
    expect(centers('[H][C@](C)(O)CC')).toEqual({ 1: 'R' });
    expect(centers('C[C@@H](O)CC')).toEqual({ 1: 'R' });
  });
});

describe('CIP: drugs, terpenes and natural products', () => {
  it('(S)-ibuprofen and (S)-naproxen', () => {
    expect(centers('CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O')).toEqual({ 10: 'S' });
    expect(centers('C[C@H](C(=O)O)c1ccc2cc(OC)ccc2c1')).toEqual({ 1: 'S' });
  });
  it('(R)- and (S)-carvone', () => {
    expect(centers('CC1=CC[C@H](CC1=O)C(C)=C')).toEqual({ 4: 'R' });
    expect(centers('CC1=CC[C@@H](CC1=O)C(C)=C')).toEqual({ 4: 'S' }); // PubChem: (5S)-…
  });
  it('(S)-limonene, (S)-linalool, (S)-nicotine', () => {
    expect(centers('CC1=CC[C@H](CC1)C(C)=C')).toEqual({ 4: 'S' });
    expect(centers('CC(C)=CCC[C@](C)(O)C=C')).toEqual({ 6: 'S' });
    expect(centers('CN1CCC[C@H]1c1cccnc1')).toEqual({ 5: 'S' });
  });
  it('L-ascorbic acid: (2R)-ring, (1S)-side chain', () => {
    expect(centers('OC[C@H](O)[C@H]1OC(=O)C(O)=C1O')).toEqual({ 2: 'S', 4: 'R' });
  });
  it('(−)-menthol is (1R,2S,5R)', () => {
    // atom 9 = C1 (carbinol), atom 3 = C2 (isopropyl), atom 6 = C5 (methyl)
    expect(centers('CC(C)[C@@H]1CC[C@@H](C)C[C@H]1O')).toEqual({ 3: 'S', 6: 'R', 9: 'R' });
  });
  it('cholesterol: (3S,8S,9S,10R,13R,14S,17R) and (20R) side chain', () => {
    const c = centers('C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C');
    expect(c).toEqual({ 1: 'R', 8: 'R', 11: 'S', 12: 'R', 15: 'S', 16: 'S', 20: 'R', 23: 'S' });
  });
  it('α-D-glucopyranose is (2S,3R,4S,5S,6R)', () => {
    // atom 4 = anomeric C1 (2S), 6 = C2, 8 = C3, 10 = C4, 2 = C5
    expect(centers('OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O')).toEqual({ 2: 'R', 4: 'S', 6: 'R', 8: 'S', 10: 'S' });
  });
  it('estradiol (8R,9S,13S,14S,17S) and testosterone (8R,9S,10S,13S,14S,17S)', () => {
    expect(centers('C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc34)[C@@H]1CC[C@@H]2O')).toEqual({ 1: 'S', 4: 'S', 5: 'R', 15: 'S', 18: 'S' });
    expect(centers('C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@@]34C)[C@@H]1CC[C@@H]2O')).toEqual({ 1: 'S', 4: 'S', 5: 'R', 14: 'S', 16: 'S', 19: 'S' });
  });
  it('bicyclic terpenes: (1R,4R)-camphor and (1S,5S)-α-pinene', () => {
    expect(centers('CC1(C)[C@@H]2CC[C@@]1(C)C(=O)C2')).toEqual({ 3: 'R', 6: 'R' });
    expect(centers('CC1=CC[C@H]2C[C@@H]1C2(C)C')).toEqual({ 4: 'S', 6: 'S' });
  });
});

describe('CIP: double bonds', () => {
  it('(E)- and (Z)-but-2-ene', () => {
    expect(bonds('C/C=C/C')).toEqual({ '1-2': 'E' });
    expect(bonds('C/C=C\\C')).toEqual({ '1-2': 'Z' });
  });
  it('(E)-cinnamic acid and (2E,4E)-sorbic acid', () => {
    expect(bonds('OC(=O)/C=C/c1ccccc1')).toEqual({ '3-4': 'E' });
    expect(bonds('C/C=C/C=C/C(=O)O')).toEqual({ '1-2': 'E', '3-4': 'E' });
  });
  it('all-trans retinol has four E double bonds (ring C=C unlabelled)', () => {
    expect(bonds('CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/CO)/C)/C')).toEqual({ '9-10': 'E', '11-12': 'E', '13-14': 'E', '15-16': 'E' });
  });
  it('citral: geranial is E, neral is Z; geraniol E, nerol Z', () => {
    expect(bonds('CC(C)=CCC/C(C)=C/C=O')).toEqual({ '6-8': 'E' });
    expect(bonds('CC(C)=CCC/C(C)=C\\C=O')).toEqual({ '6-8': 'Z' });
    expect(bonds('CC(C)=CCC/C(C)=C/CO')).toEqual({ '6-8': 'E' });
    expect(bonds('CC(C)=CCC/C(C)=C\\CO')).toEqual({ '6-8': 'Z' });
  });
  it('tiglic (E) vs angelic (Z) acid: priority, not the drawn cis/trans relation, decides', () => {
    expect(bonds('OC(=O)/C(C)=C/C')).toEqual({ '3-5': 'E' });
    expect(bonds('OC(=O)/C(C)=C\\C')).toEqual({ '3-5': 'Z' });
  });
  it('oleic acid is Z; tetrasubstituted haloalkene; oximes', () => {
    expect(bonds('CCCCCCCC/C=C\\CCCCCCCC(=O)O')).toEqual({ '8-9': 'Z' });
    expect(bonds('C/C(Cl)=C(\\F)Br')).toEqual({ '1-3': 'Z' });
    expect(bonds('C/C=N/O')).toEqual({ '1-2': 'E' });
    expect(bonds('C/C=N\\O')).toEqual({ '1-2': 'Z' });
  });
  it('no E/Z when one end carries identical substituents', () => {
    expect(bonds('C=C/C=C/C')).toEqual({ '2-3': 'E' });
    expect(bonds('C/C=C(/C)C')).toEqual({});
  });
});

describe('CIP: hierarchical digraph details', () => {
  it('duplicate atoms: vinyl outranks isopropyl, (3R)-4-methylpent-1-en-3-ol', () => {
    expect(centers('CC(C)[C@@H](O)C=C')).toEqual({ 3: 'R' });
    expect(centers('CC(C)[C@H](O)C=C')).toEqual({ 3: 'S' });
  });
  it('rule 2 (isotopes): (R)-1-deuterioethanol, 13C-labelled ethanol', () => {
    expect(centers('[2H][C@@H](O)C')).toEqual({ 1: 'R' });
    expect(centers('[13CH3][C@H](O)C')).toEqual({ 1: 'R' });
    expect(centers('C[C@H](C)O')).toEqual({}); // unlabelled: not a stereocentre
  });
  it('sulfoxides use a lone pair as the lowest ligand', () => {
    expect(centers('C[S@](=O)c1ccc(C)cc1')).toEqual({ 1: 'S' }); // PubChem: (S)-methylsulfinyl
    expect(centers('C[S@@](=O)c1ccc(C)cc1')).toEqual({ 1: 'R' });
    expect(centers('CC[S@+]([O-])C')).toEqual({ 2: 'R' });
  });
  it("root multiple bonds (P=O) are treated as single: -OMe outranks =O", () => {
    expect(centers('C[P@](=O)(OC)c1ccccc1')).toEqual({ 1: 'S' });
  });
  it('quaternary ammonium centre', () => {
    expect(centers('C[N@+](CC)(CCC)c1ccccc1')).toEqual({ 1: 'S' });
  });
  it('non-stereogenic atoms carrying a stereo flag get no label', () => {
    expect(centers('C[C@H](C)O')).toEqual({});
    expect(centers('CC[C@](C)(CC)O')).toEqual({});
  });
  it('mancude rings are Kekulé-independent (averaged duplicate atomic numbers)', () => {
    const variants = ['O[C@H](c1ccccn1)c1cccnc1', 'O[C@H](C1=NC=CC=C1)C1=CN=CC=C1', 'O[C@H](C1=CC=CC=N1)C1=CC=CN=C1'];
    for (const v of variants) expect(centers(v), v).toEqual({ 1: 'S' });
    for (const v of ['C[C@H](Cl)c1ccccn1', 'C[C@H](Cl)C1=NC=CC=C1', 'C[C@H](Cl)C1=CC=CC=N1']) expect(centers(v), v).toEqual({ 1: 'S' });
    for (const v of ['Oc1ccc(cc1)[C@@H](N)C(=O)O', 'OC1=CC=C(C=C1)[C@@H](N)C(=O)O', 'OC1C=CC(=CC=1)[C@@H](N)C(=O)O'])
      expect(centers(v), v).toEqual({ 7: 'R' });
    // pyridine: C2/C6 duplicates average N and C (6.5); C3–C5 average C (6); N gets 6
    const z = _mancudeDuplicateZ(parseSmiles('C1=CC=NC=C1'));
    expect(z.map((v) => v ?? 0)).toEqual([6, 6, 6.5, 6, 6.5, 6]);
  });
});

describe('CIP: like/unlike, meso and pseudoasymmetric centres (rules 4 and 5)', () => {
  it('meso-tartaric acid is (2R,3S); L-tartaric acid (2R,3R)', () => {
    expect(centers('OC(=O)[C@@H](O)[C@@H](O)C(=O)O')).toEqual({ 3: 'S', 5: 'R' }); // PubChem CID 447315
    expect(centers('[C@@H]([C@H](C(=O)O)O)(C(=O)O)O')).toEqual({ 0: 'R', 1: 'R' }); // PubChem CID 444305
  });
  it('2,3,4-trihydroxyglutaric acids: C3 is pseudoasymmetric (r/s)', () => {
    expect(centers('OC(=O)[C@@H](O)[C@H](O)[C@@H](O)C(=O)O')).toEqual({ 3: 'S', 5: 's', 7: 'R' });
    expect(centers('OC(=O)[C@@H](O)[C@@H](O)[C@@H](O)C(=O)O')).toEqual({ 3: 'S', 5: 'r', 7: 'R' });
    // mirror image of the first meso form is the same compound: C3 stays 's'
    expect(centers('OC(=O)[C@H](O)[C@@H](O)[C@H](O)C(=O)O')).toEqual({ 3: 'R', 5: 's', 7: 'S' });
  });
  it('pentane-2,3,4-triol: C3 unlabelled in (2R,4R), pseudoasymmetric in (2R,4S)', () => {
    expect(centers('C[C@@H](O)[C@H](O)[C@@H](C)O')).toEqual({ 1: 'R', 5: 'R' });
    expect(centers('C[C@@H](O)[C@@H](O)[C@H](C)O')).toEqual({ 1: 'R', 3: 's', 5: 'S' });
    expect(centers('C[C@@H](O)[C@H](O)[C@H](C)O')).toEqual({ 1: 'R', 3: 'r', 5: 'S' });
  });
  it('1,4-disubstituted cyclohexanes get r/s via auxiliary descriptors', () => {
    expect(centers('NC[C@H]1CC[C@@H](CC1)C(=O)O')).toEqual({ 2: 'r', 5: 'r' }); // trans-tranexamic acid (1r,4r)
    expect(centers('C[C@H]1CC[C@@H](C)CC1')).toEqual({ 1: 's', 4: 's' });
    expect(centers('C[C@H]1CC[C@H](C)CC1')).toEqual({ 1: 'r', 4: 'r' });
    expect(centers('C[C@@]1(CC)CC[C@@](C)(CC)CC1')).toEqual({ 1: 'r', 6: 'r' });
  });
  it('a dienol whose centre is decided by Z > E (rule 3) is chiral (upper case)', () => {
    expect(centers('C/C=C\\[C@H](O)/C=C/C')).toEqual({ 3: 'R' });
    expect(bonds('C/C=C\\[C@H](O)/C=C/C')).toEqual({ '1-2': 'Z', '5-6': 'E' });
    expect(centers('C/C=C/[C@H](O)/C=C/C')).toEqual({}); // E,E: two identical ligands
  });
  it('inositol stereoisomers mix R/S and r/s (rule 4b like/unlike and rule 5)', () => {
    expect(centers('O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O')).toEqual({ 1: 'S', 2: 'R', 4: 'r', 6: 'S', 8: 'R', 10: 's' });
    expect(centers('O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O')).toEqual({ 1: 'S', 2: 's', 4: 'R', 6: 'S', 8: 'r', 10: 'R' });
    expect(centers('O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O')).toEqual({ 1: 'S', 2: 'R', 4: 'R', 6: 'R', 8: 'R', 10: 'S' });
  });
});

describe('CIP: input from 2D drawings and robustness', () => {
  it('works on wedge drawings after perceiveStereo2D', () => {
    // butan-2-ol drawn with C2 at the origin, OH up (wedge, towards the viewer), H implicit (away):
    // O → ethyl (lower right) → methyl (lower left) runs clockwise on screen → (R)
    const m = new Mol();
    const c2 = m.addAtom({ el: 'C', x: 0, y: 0 });
    const c1 = m.addAtom({ el: 'C', x: -0.866, y: 0.5 });
    const c3 = m.addAtom({ el: 'C', x: 0.866, y: 0.5 });
    const c4 = m.addAtom({ el: 'C', x: 1.732, y: 0 });
    const o = m.addAtom({ el: 'O', x: 0, y: -1 });
    m.addBond(c2, c1);
    m.addBond(c2, c3);
    m.addBond(c3, c4);
    m.addBond(c2, o, 1, 'wedge');
    perceiveStereo2D(m);
    expect(assignCIP(m).centers.get(c2)).toBe('R');
    m.bonds[3].style = 'hash';
    perceiveStereo2D(m);
    expect(assignCIP(m).centers.get(c2)).toBe('S');
  });
  it('never throws on malformed specs and leaves them unlabelled', () => {
    const m = parseSmiles('CC(O)CC=CC');
    m.tetra.push({ center: 1, nbrs: [0, 5, 2, -1], ccw: true }); // atom 5 is not a neighbour
    m.tetra.push({ center: 99, nbrs: [0, 1, 2, 3], ccw: true });
    m.tetra.push({ center: 1, nbrs: [0, 0, 2, -1], ccw: true });
    m.dbStereo.push({ bond: 0, a: 1, b: 2, cis: true }); // single bond
    m.dbStereo.push({ bond: 99, a: 1, b: 2, cis: true });
    let r: ReturnType<typeof assignCIP> | undefined;
    expect(() => (r = assignCIP(m))).not.toThrow();
    expect(r!.centers.size).toBe(0);
    expect(r!.bonds.size).toBe(0);
    expect(assignCIP(new Mol()).centers.size).toBe(0);
  });
  it('stays bounded on highly symmetric cages (digraph budget)', () => {
    const c60 = parseSmiles('c12c3c4c5c1c1c6c7c2c2c8c3c3c9c4c4c%10c5c5c1c1c6c6c%11c7c2c2c7c8c3c3c8c9c4c4c9c%10c5c5c1c1c6c6c%11c2c2c7c3c3c8c4c4c9c5c1c1c6c2c3c41');
    c60.tetra = c60.atoms.map((_, i) => {
      const nb = c60.neighbors(i);
      return { center: i, nbrs: [nb[0], nb[1], nb[2], -1] as [number, number, number, number], ccw: true };
    });
    const t0 = performance.now();
    const r = assignCIP(c60);
    expect(r.centers.size).toBe(0);
    expect(performance.now() - t0).toBeLessThan(5000);
  });
  it('labels a 15-centre macrolide quickly', () => {
    // erythromycin A (PubChem isomeric SMILES)
    const smi =
      'CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O';
    const t0 = performance.now();
    const r = assignCIP(parseSmiles(smi));
    expect(performance.now() - t0).toBeLessThan(500);
    expect(r.centers.size).toBe(18);
    expect([...r.centers.values()].every((l) => l === 'R' || l === 'S')).toBe(true);
  });
});
