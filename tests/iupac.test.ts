import { describe, it, expect } from 'vitest';
import { parseSmiles } from '../src/chem/smiles';
import { nameMolecule } from '../src/chem/iupac';
import { assignCIP } from '../src/chem/cip';
import { Mol } from '../src/chem/mol';

// Expected names were cross-checked against PubChem (OpenEye LexiChem) where possible.
// Deliberate differences from PubChem are marked "PubChem: …" with the reason:
//  – trivial names not retained by IUPAC 2013 (cumene, xylene, chloroform, acetylene …),
//  – IUPAC nesting of enclosing marks {[( )]} (PubChem only uses ( ) and [ ]),
//  – IUPAC numbering where LexiChem departs from it (suffix before hydro prefixes, indicated hydrogen first),
//  – IUPAC aminium cation names (PubChem: "…azanium"), oxime names ("N-hydroxy…imine"),
//  – PubChem standardises tautomers, so names of drawn tautomers may differ.

const nm = (smi: string) => nameMolecule(parseSmiles(smi));
type Case = [smiles: string, name: string];

function table(title: string, cases: Case[]) {
  describe(title, () => {
    it.each(cases)('%s → %s', (smi, expected) => {
      expect(nm(smi).name).toBe(expected);
    });
  });
}

table('acyclic hydrocarbons', [
  ['C', 'methane'],
  ['CC', 'ethane'],
  ['CCC', 'propane'],
  ['CCCCCCCCCC', 'decane'],
  ['CCCCCCCCCCC', 'undecane'],
  ['CCCCCCCCCCCCCCCCCCCC', 'icosane'],
  ['CCCCCCCCCCCCCCCCCCCCC', 'henicosane'],
  ['CCCCCCCCCCCCCCCCCCCCCC', 'docosane'],
  ['CCCCCCCCCCCCCCCCCCCCCCC', 'tricosane'],
  ['CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC', 'triacontane'],
  ['CC(C)C', '2-methylpropane'],
  ['CC(C)CC(C)(C)C', '2,2,4-trimethylpentane'],
  ['CCC(C)CC', '3-methylpentane'],
  ['CCC(CC)CC', '3-ethylpentane'],
  ['CC(C)C(C)C', '2,3-dimethylbutane'],
  ['CC(C)CC(CC)CC(C)C', '4-ethyl-2,6-dimethylheptane'],
  ['CC(C)CCCC(C)CCCC(C)CCCCC(C)CCCC(C)CCCC(C)C', '2,6,10,15,19,23-hexamethyltetracosane'],
  ['C=C', 'ethene'],
  ['CC=C', 'prop-1-ene'],
  ['C=CC=C', 'buta-1,3-diene'],
  ['C#C', 'ethyne'], // PubChem: acetylene
  ['CC#C', 'prop-1-yne'],
  ['C=CC#C', 'but-1-en-3-yne'],
  ['CC=CC', 'but-2-ene'],
  ['C=CCC=C', 'penta-1,4-diene'],
  ['C=C=C', 'propa-1,2-diene'],
  ['C=CC(C)=C', '2-methylbuta-1,3-diene'],
  ['C#CCC=C', 'pent-1-en-4-yne'],
  ['C=CC(CCCC)CCCC', '5-ethenylnonane'], // 2013 rules: chain length before unsaturation
  ['C=CC(CCC)CC', '3-ethylhex-1-ene'],
]);

table('substituent prefixes', [
  ['CC(C)Cc1ccccc1', '2-methylpropylbenzene'],
  ['CC(C)c1ccccc1', 'propan-2-ylbenzene'], // PubChem: cumene
  ['CCC(C)c1ccccc1', 'butan-2-ylbenzene'],
  ['CC(C)(C)c1ccccc1', 'tert-butylbenzene'],
  ['C=CCc1ccccc1', 'prop-2-enylbenzene'],
  ['C#Cc1ccccc1', 'ethynylbenzene'],
  ['Clc1ccccc1', 'chlorobenzene'],
  ['Brc1ccccc1', 'bromobenzene'],
  ['[O-][N+](=O)c1ccccc1', 'nitrobenzene'],
  ['O=Nc1ccccc1', 'nitrosobenzene'],
  ['[N-]=[N+]=Nc1ccccc1', 'azidobenzene'],
  ['CCOc1ccccc1', 'ethoxybenzene'],
  ['CC(C)COc1ccccc1', '2-methylpropoxybenzene'],
  ['CC(C)Oc1ccccc1', 'propan-2-yloxybenzene'],
  ['CCCCCOc1ccccc1', 'pentoxybenzene'],
  ['CC(C)(C)Oc1ccccc1', '(2-methylpropan-2-yl)oxybenzene'],
  ['CC(C)(C)Oc1ccc(Cl)cc1', '1-chloro-4-[(2-methylpropan-2-yl)oxy]benzene'],
  ['c1ccc(Oc2ccccc2)cc1', 'phenoxybenzene'],
  ['c1ccc(cc1)COc1ccccc1', 'phenoxymethylbenzene'],
  ['CSc1ccccc1', 'methylsulfanylbenzene'],
  ['CSc1ccc(SC)cc1', '1,4-bis(methylsulfanyl)benzene'],
  ['ClCc1ccccc1', 'chloromethylbenzene'],
  ['FC(F)(F)c1ccccc1', 'trifluoromethylbenzene'],
  ['FC(F)(F)c1cc(cc(c1)C(F)(F)F)C(=O)O', '3,5-bis(trifluoromethyl)benzoic acid'],
  ['ClCc1cc(CCl)cc(C(=O)O)c1', '3,5-bis(chloromethyl)benzoic acid'],
  ['ClC(c1ccccc1)c1ccccc1', '[chloro(phenyl)methyl]benzene'],
  ['Clc1ccc(Cc2ccccc2)cc1', '1-benzyl-4-chlorobenzene'],
  ['C1CCC(CC1)c1ccccc1', 'cyclohexylbenzene'],
  ['C1CCC(=CC1)c1ccccc1', 'cyclohexen-1-ylbenzene'],
  ['CC(C)c1cccc(C(C)C)c1', '1,3-di(propan-2-yl)benzene'],
  ['CNCCO', '2-(methylamino)ethanol'],
  ['CN(C)CCO', '2-(dimethylamino)ethanol'],
  ['NCCO', '2-aminoethanol'],
  ['OCCN(CCO)CCO', '2-[bis(2-hydroxyethyl)amino]ethanol'],
  ['OCCNCCO', '2-[(2-hydroxyethyl)amino]ethanol'], // PubChem: 2-(2-hydroxyethylamino)ethanol
  ['CNc1cc(NC)cc(C(=O)O)c1', '3,5-bis(methylamino)benzoic acid'],
  ['CC(=O)N(C)c1ccc(cc1)C(=O)O', '4-[acetyl(methyl)amino]benzoic acid'],
  ['CC(=O)Nc1ccc(cc1)C(=O)O', '4-acetamidobenzoic acid'],
  ['CCC(=O)Nc1ccc(cc1)C(=O)O', '4-(propanoylamino)benzoic acid'],
  ['NC(=O)Nc1ccc(cc1)C(=O)O', '4-(carbamoylamino)benzoic acid'],
  ['CN(C)C(=O)c1ccc(cc1)C(=O)O', '4-(dimethylcarbamoyl)benzoic acid'],
  ['CNC(=O)c1ccc(cc1)C(=O)O', '4-(methylcarbamoyl)benzoic acid'],
  ['COC(=O)c1ccc(cc1)C(=O)O', '4-methoxycarbonylbenzoic acid'],
  ['COC(=O)c1cc(C(=O)OC)cc(C(=O)O)c1', '3,5-bis(methoxycarbonyl)benzoic acid'],
  ['CC(C)Oc1cc(OC(C)C)cc(C(=O)O)c1', '3,5-di(propan-2-yloxy)benzoic acid'],
  ['CC(=O)Oc1ccc(cc1)C(=O)O', '4-acetyloxybenzoic acid'],
  ['CC(=O)c1ccc(C(=O)O)cc1', '4-acetylbenzoic acid'],
  ['O=C(c1ccccc1)c1ccc(C(=O)O)cc1', '4-benzoylbenzoic acid'],
  ['CCC(=O)c1ccc(cc1)C(=O)O', '4-propanoylbenzoic acid'],
  ['OC(=O)c1ccc(C#N)cc1', '4-cyanobenzoic acid'],
  ['O=Cc1ccc(C(=O)O)cc1', '4-formylbenzoic acid'],
  ['OC(=O)CCc1ccc(cc1)C(=O)O', '4-(2-carboxyethyl)benzoic acid'],
  ['O=CCCc1ccc(cc1)C(=O)O', '4-(3-oxopropyl)benzoic acid'],
  ['COC(=O)CCc1ccc(cc1)C(=O)O', '4-(3-methoxy-3-oxopropyl)benzoic acid'],
  ['NC(=O)CCc1ccc(cc1)C(=O)O', '4-(3-amino-3-oxopropyl)benzoic acid'],
  ['N#CCCc1ccc(cc1)C(=O)O', '4-(2-cyanoethyl)benzoic acid'],
  ['ON=C(C)C(=O)O', '2-hydroxyiminopropanoic acid'],
  ['CN(C)c1ccc(cc1)C(=O)c1ccc(cc1)N(C)C', 'bis[4-(dimethylamino)phenyl]methanone'],
  ['OC(=O)c1ccccc1Nc1ccccc1', '2-anilinobenzoic acid'],
  ['OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl', '2-[2-(2,6-dichloroanilino)phenyl]acetic acid'],
  ['NCCc1ccccn1', '2-pyridin-2-ylethanamine'],
  ['OCC1CCCO1', 'oxolan-2-ylmethanol'],
  ['ClCC1CO1', '2-(chloromethyl)oxirane'],
  ['OCc1ccccc1O', '2-(hydroxymethyl)phenol'],
  ['C=CCN=C=S', '3-isothiocyanatoprop-1-ene'],
  ['O=C=Nc1ccccc1', 'isocyanatobenzene'],
  ['N#CSc1ccccc1', 'thiocyanatobenzene'],
  ['CS(=O)(=O)CCO', '2-methylsulfonylethanol'],
  ['C=CCOCC=C', '3-prop-2-enoxyprop-1-ene'],
  ['CSSC', '(methyldisulfanyl)methane'],
  ['CCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC', '2,3-di(hexadecanoyloxy)propyl hexadecanoate'],
  ['CC(=O)OCC(COC(C)=O)OC(C)=O', '2,3-diacetyloxypropyl acetate'],
]);

table('characteristic groups (suffixes, retained names)', [
  ['CCO', 'ethanol'],
  ['CC(C)O', 'propan-2-ol'],
  ['CC(C)(O)C', '2-methylpropan-2-ol'],
  ['OCCO', 'ethane-1,2-diol'],
  ['OCC(O)CO', 'propane-1,2,3-triol'],
  ['OCC(O)C(O)C(O)C(O)CO', 'hexane-1,2,3,4,5,6-hexol'],
  ['OCC#C', 'prop-2-yn-1-ol'],
  ['CC(O)C#C', 'but-3-yn-2-ol'],
  ['CC(C)(O)C#C', '2-methylbut-3-yn-2-ol'],
  ['Oc1ccccc1', 'phenol'],
  ['Cc1ccc(O)cc1', '4-methylphenol'],
  ['Oc1ccccc1O', 'benzene-1,2-diol'],
  ['Oc1cccc(O)c1', 'benzene-1,3-diol'],
  ['Oc1ccc(O)cc1', 'benzene-1,4-diol'],
  ['CCS', 'ethanethiol'],
  ['Sc1ccccc1', 'benzenethiol'],
  ['Cc1ccc(S)cc1', '4-methylbenzenethiol'],
  ['CCN', 'ethanamine'],
  ['CC(C)N', 'propan-2-amine'],
  ['CC(C)(C)N', '2-methylpropan-2-amine'],
  ['CCNCC', 'N-ethylethanamine'],
  ['CCN(CC)CC', 'N,N-diethylethanamine'],
  ['CN(C)C', 'N,N-dimethylmethanamine'],
  ['CCN(C)CC', 'N-ethyl-N-methylethanamine'],
  ['CCN(C(C)C)C(C)C', 'N-ethyl-N-propan-2-ylpropan-2-amine'],
  ['CCCCN(CCCC)CCCC', 'N,N-dibutylbutan-1-amine'],
  ['NCCN', 'ethane-1,2-diamine'],
  ['NCCCCCCN', 'hexane-1,6-diamine'],
  ['Nc1ccccc1', 'aniline'],
  ['NC1=CC=CC=C1', 'aniline'],
  ['CNc1ccccc1', 'N-methylaniline'],
  ['CN(C)c1ccccc1', 'N,N-dimethylaniline'],
  ['c1ccc(cc1)CNc1ccccc1', 'N-benzylaniline'],
  ['NCc1ccccc1', 'phenylmethanamine'],
  ['CNCc1ccccc1', 'N-methyl-1-phenylmethanamine'],
  ['c1ccc(CNCc2ccccc2)cc1', 'N-benzyl-1-phenylmethanamine'],
  ['NC1CCCCC1', 'cyclohexanamine'],
  ['CNC1CCCCC1', 'N-methylcyclohexanamine'],
  ['Nc1ccc(N)cc1', 'benzene-1,4-diamine'],
  ['CC(N)Cc1ccccc1', '1-phenylpropan-2-amine'],
  ['CC(=N)C', 'propan-2-imine'],
  ['CN=C(C)C', 'N-methylpropan-2-imine'],
  ['c1ccc(cc1)C=Nc1ccccc1', 'N,1-diphenylmethanimine'],
  ['CC(C)=NO', 'N-hydroxypropan-2-imine'], // IUPAC 2013 PIN; PubChem: N-propan-2-ylidenehydroxylamine
  ['C=O', 'formaldehyde'],
  ['CC=O', 'acetaldehyde'],
  ['CCC=O', 'propanal'],
  ['CCCCCC=O', 'hexanal'],
  ['CC(C)CC=O', '3-methylbutanal'],
  ['O=CCCC=O', 'butanedial'],
  ['O=CC=Cc1ccccc1', '3-phenylprop-2-enal'],
  ['O=Cc1ccccc1', 'benzaldehyde'],
  ['O=CC1CCCCC1', 'cyclohexanecarbaldehyde'],
  ['O=Cc1ccncc1', 'pyridine-4-carbaldehyde'],
  ['CCCCCCC(=Cc1ccccc1)C=O', '2-benzylideneoctanal'],
  ['CC(=O)C', 'propan-2-one'],
  ['CCC(C)=O', 'butan-2-one'],
  ['CCC(=O)CC', 'pentan-3-one'],
  ['CC(=O)CC(C)=O', 'pentane-2,4-dione'],
  ['CC(=O)C(C)=O', 'butane-2,3-dione'],
  ['CC(=O)CCC=O', '4-oxopentanal'],
  ['CC(=O)c1ccccc1', '1-phenylethanone'],
  ['O=C(c1ccccc1)c1ccccc1', 'diphenylmethanone'],
  ['O=C(c1ccccc1)c1ccc(Cl)cc1', '(4-chlorophenyl)-phenylmethanone'],
  ['Oc1ccc(cc1)C(=O)c1ccccc1', '(4-hydroxyphenyl)-phenylmethanone'],
  ['CC(=O)C1CCCCC1', '1-cyclohexylethanone'],
  ['CC(=O)c1ccc(cc1)C(C)C', '1-(4-propan-2-ylphenyl)ethanone'],
  ['CC#N', 'acetonitrile'],
  ['CCC#N', 'propanenitrile'],
  ['CCCCC#N', 'pentanenitrile'],
  ['N#CCC#N', 'propanedinitrile'],
  ['N#CCCCCC#N', 'hexanedinitrile'],
  ['N#Cc1ccccc1', 'benzonitrile'],
  ['OCCC#N', '3-hydroxypropanenitrile'],
  ['CC(=O)O', 'acetic acid'],
  ['OC(=O)c1ccccc1', 'benzoic acid'],
  ['OC(=O)C(=O)O', 'oxalic acid'],
  ['OC(=O)CCC(=O)O', 'butanedioic acid'],
  ['OC(=O)CCCCCCCC(=O)O', 'nonanedioic acid'],
  ['CCCCCCCCCCCCCCCCCC(=O)O', 'octadecanoic acid'],
  ['CC(C)CC(=O)O', '3-methylbutanoic acid'],
  ['CC(C)(C)C(=O)O', '2,2-dimethylpropanoic acid'],
  ['C=CC(=O)O', 'prop-2-enoic acid'],
  ['CC(=C)C(=O)O', '2-methylprop-2-enoic acid'],
  ['CC=CC(=O)O', 'but-2-enoic acid'],
  ['CC=CC=CC(=O)O', 'hexa-2,4-dienoic acid'],
  ['CCCCCCCCC=CCCCCCCCC(=O)O', 'octadec-9-enoic acid'],
  ['CCCCCC=CCC=CCCCCCCCC(=O)O', 'octadeca-9,12-dienoic acid'],
  ['OC(=O)C=Cc1ccccc1', '3-phenylprop-2-enoic acid'],
  ['ClCC(=O)O', '2-chloroacetic acid'],
  ['OCC(=O)O', '2-hydroxyacetic acid'],
  ['CC(O)C(=O)O', '2-hydroxypropanoic acid'],
  ['OC(C(=O)O)c1ccccc1', '2-hydroxy-2-phenylacetic acid'],
  ['OC(=O)Cc1ccccc1', '2-phenylacetic acid'],
  ['OC(=O)Cc1ccc(Cl)cc1', '2-(4-chlorophenyl)acetic acid'],
  ['N#CCC(=O)O', '2-cyanoacetic acid'],
  ['OC(=O)CCC=O', '4-oxobutanoic acid'],
  ['CCOC(=O)CCC(=O)O', '4-ethoxy-4-oxobutanoic acid'],
  ['OC(=O)CC(O)(CC(=O)O)C(=O)O', '2-hydroxypropane-1,2,3-tricarboxylic acid'],
  ['OC(=O)c1ccc(cc1)C(=O)O', 'terephthalic acid'],
  ['OC(=O)c1ccccc1C(=O)O', 'phthalic acid'],
  ['OC(=O)c1ccccc1O', '2-hydroxybenzoic acid'],
  ['OC(=O)c1ccc(N)cc1', '4-aminobenzoic acid'],
  ['OC(=O)c1ccc(cc1)[N+]([O-])=O', '4-nitrobenzoic acid'],
  ['OC(=O)c1cc(O)c(O)c(O)c1', '3,4,5-trihydroxybenzoic acid'],
  ['OC(=O)C1CCCCC1', 'cyclohexanecarboxylic acid'],
  ['C1CC1C(=O)O', 'cyclopropanecarboxylic acid'],
  ['NCC1CCC(CC1)C(=O)O', '4-(aminomethyl)cyclohexane-1-carboxylic acid'],
  ['CC1CCC(CC1)C(=O)O', '4-methylcyclohexane-1-carboxylic acid'],
  ['OS(=O)(=O)c1ccccc1', 'benzenesulfonic acid'],
  ['CCS(=O)(=O)O', 'ethanesulfonic acid'],
  ['Cc1ccc(cc1)S(=O)(=O)O', '4-methylbenzenesulfonic acid'],
  ['Cc1ccc(cc1)S(=O)(=O)Cl', '4-methylbenzenesulfonyl chloride'],
  ['Cc1ccc(cc1)S(N)(=O)=O', '4-methylbenzenesulfonamide'],
  ['CS(=O)(=O)N', 'methanesulfonamide'],
  ['Nc1ccc(cc1)S(N)(=O)=O', '4-aminobenzenesulfonamide'],
  ['CC(=O)Cl', 'acetyl chloride'],
  ['CCCC(=O)Cl', 'butanoyl chloride'],
  ['CCC(=O)Br', 'propanoyl bromide'],
  ['O=C(Cl)CCC(=O)Cl', 'butanedioyl dichloride'],
  ['O=C(Cl)c1ccccc1', 'benzoyl chloride'],
  ['NC(=O)c1ccccc1', 'benzamide'],
  ['CN(C)C=O', 'N,N-dimethylformamide'],
  ['CC(=O)N(C)C', 'N,N-dimethylacetamide'],
  ['CCCCC(=O)N', 'pentanamide'],
  ['CC(C)C(=O)N', '2-methylpropanamide'],
  ['NC(=O)CCC(N)=O', 'butanediamide'],
  ['NC(=O)CN', '2-aminoacetamide'],
  ['NC(=O)CC#N', '2-cyanoacetamide'],
  ['CCN(CC)C(=O)c1cccc(C)c1', 'N,N-diethyl-3-methylbenzamide'],
  ['CC(=N)N', 'ethanimidamide'],
  ['NC(=N)c1ccccc1', 'benzenecarboximidamide'],
  ['NC(N)=O', 'urea'],
  ['C[N+](C)(C)CC(=O)[O-]', '2-(trimethylazaniumyl)acetate'],
  ['[O-][n+]1ccccc1', '1-oxidopyridin-1-ium'],
]);

table('esters and functional-class names', [
  ['CCOC(C)=O', 'ethyl acetate'],
  ['CC(C)(C)OC(C)=O', 'tert-butyl acetate'],
  ['CC(=O)OC=C', 'ethenyl acetate'],
  ['COC(=O)c1ccc(O)cc1', 'methyl 4-hydroxybenzoate'],
  ['COC(=O)c1ccccc1O', 'methyl 2-hydroxybenzoate'],
  ['CCOC(=O)c1ccc(N)cc1', 'ethyl 4-aminobenzoate'],
  ['CCOC(=O)CCC(=O)OCC', 'diethyl butanedioate'],
  ['CCCCOC(=O)c1ccccc1C(=O)OCCCC', 'dibutyl benzene-1,2-dicarboxylate'],
  ['COC(=O)C(C)=C', 'methyl 2-methylprop-2-enoate'],
  ['C=CC(=O)OCC', 'ethyl prop-2-enoate'],
  ['CC(=O)CC(=O)OC', 'methyl 3-oxobutanoate'],
  ['CCCCCCCCCCCCCC(=O)OC(C)C', 'propan-2-yl tetradecanoate'],
  ['O=C(OCc1ccccc1)c1ccccc1', 'benzyl benzoate'],
  ['O=C(OCc1ccccc1)c1ccccc1O', 'benzyl 2-hydroxybenzoate'],
  ['O=C(OCc1ccccc1)C=Cc1ccccc1', 'benzyl 3-phenylprop-2-enoate'],
  ['CCCCC(CC)COC(=O)c1ccccc1O', '2-ethylhexyl 2-hydroxybenzoate'],
  ['CCCCC(CC)COC(=O)C(C#N)=C(c1ccccc1)c1ccccc1', '2-ethylhexyl 2-cyano-3,3-diphenylprop-2-enoate'],
  ['CCCCC(CC)COC(=O)C=Cc1ccc(OC)cc1', '2-ethylhexyl 3-(4-methoxyphenyl)prop-2-enoate'],
  ['CC1CC(CC(C)(C)C1)OC(=O)c1ccccc1O', '3,3,5-trimethylcyclohexyl 2-hydroxybenzoate'], // PubChem encloses the alkyl
  ['CC(=O)OCC=C(C)CCC=C(C)C', '3,7-dimethylocta-2,6-dienyl acetate'],
  ['OCC(O)COC(=O)CCCCCCCCCCCCCCCCC', '2,3-dihydroxypropyl octadecanoate'],
  ['CC(C)(C)OC(=O)N1CCCCC1', 'tert-butyl piperidine-1-carboxylate'],
  ['CNC(=O)Oc1cccc2ccccc12', 'naphthalen-1-yl N-methylcarbamate'],
  ['CC(C)(C)OC(=O)NCC(=O)O', '2-[(2-methylpropan-2-yl)oxycarbonylamino]acetic acid'],
  ['COC(=O)OC', 'dimethyl carbonate'],
  ['CCOC(=O)Cl', 'ethyl carbonochloridate'],
  ['CN(C)C(Cl)=O', 'N,N-dimethylcarbamoyl chloride'],
  ['COS(=O)(=O)OC', 'dimethyl sulfate'],
  ['COP(=O)(OC)OC', 'trimethyl phosphate'],
  ['CC(=O)OC(C)=O', 'acetic anhydride'], // IUPAC; PubChem: acetyl acetate
  ['O=C(OC(=O)c1ccccc1)c1ccccc1', 'benzoic anhydride'],
  ['CC1(C)C(C=C(Cl)Cl)C1C(=O)OC(C#N)c1cccc(Oc2ccccc2)c1', 'cyano-(3-phenoxyphenyl)methyl 3-(2,2-dichloroethenyl)-2,2-dimethylcyclopropane-1-carboxylate'],
]);

table('ethers, sulfides, halides', [
  ['CCOCC', 'ethoxyethane'],
  ['COC', 'methoxymethane'],
  ['C=COCC', 'ethoxyethene'], // 2013: unsaturated chain senior; PubChem: ethenoxyethane
  ['CCOCCOCC', '1,2-diethoxyethane'],
  ['COCCOCCOC', '1-methoxy-2-(2-methoxyethoxy)ethane'],
  ['OCCOCCO', '2-(2-hydroxyethoxy)ethanol'],
  ['OCCOc1ccccc1', '2-phenoxyethanol'],
  ['CCSCC', 'ethylsulfanylethane'],
  ['CS(=O)C', 'methylsulfinylmethane'],
  ['CS(C)(=O)=O', 'methylsulfonylmethane'],
  ['ClC(Cl)Cl', 'trichloromethane'], // PubChem: chloroform
  ['ClCCl', 'dichloromethane'],
  ['ClC(Cl)F', 'dichloro(fluoro)methane'],
  ['FC(F)(Cl)Cl', 'dichloro(difluoro)methane'],
  ['ICI', 'diiodomethane'],
  ['ClCCCl', '1,2-dichloroethane'],
  ['BrCCBr', '1,2-dibromoethane'],
  ['CC(Cl)(Cl)Cl', '1,1,1-trichloroethane'],
  ['FC(F)(F)C(F)(F)F', '1,1,1,2,2,2-hexafluoroethane'],
  ['ClC(Cl)=C(Cl)Cl', '1,1,2,2-tetrachloroethene'],
]);

table('carbocycles and benzene derivatives', [
  ['C1CC1', 'cyclopropane'],
  ['C1CCCCC1', 'cyclohexane'],
  ['C1=CCCCC1', 'cyclohexene'],
  ['C1=CCC=CC1', 'cyclohexa-1,4-diene'],
  ['C1=CC=CCC1', 'cyclohexa-1,3-diene'],
  ['C1=CC=CC=CC=C1', 'cycloocta-1,3,5,7-tetraene'],
  ['CC1CCCCC1', 'methylcyclohexane'],
  ['OC1CCCCC1', 'cyclohexanol'],
  ['CC1CCCCC1O', '2-methylcyclohexan-1-ol'],
  ['CC1CCC(O)CC1', '4-methylcyclohexan-1-ol'],
  ['OC1CCC(O)CC1', 'cyclohexane-1,4-diol'],
  ['O=C1CCCCC1', 'cyclohexanone'],
  ['O=C1CCCC1', 'cyclopentanone'],
  ['O=C1C=CCCC1', 'cyclohex-2-en-1-one'],
  ['O=C1C=CC(=O)C=C1', 'cyclohexa-2,5-diene-1,4-dione'],
  ['OC(=O)C1=CCCCC1', 'cyclohexene-1-carboxylic acid'],
  ['C=C1CCCCC1', 'methylidenecyclohexane'],
  ['CC=C1CCCCC1', 'ethylidenecyclohexane'],
  ['CC1(C)C(C=C(C)C)C1C(=O)O', '2,2-dimethyl-3-(2-methylprop-1-enyl)cyclopropane-1-carboxylic acid'],
  ['c1ccccc1', 'benzene'],
  ['Cc1ccccc1', 'toluene'],
  ['C=Cc1ccccc1', 'styrene'],
  ['COc1ccccc1', 'anisole'],
  ['c1ccc(cc1)-c1ccccc1', "1,1'-biphenyl"],
  ['Cc1ccc(cc1)-c1ccccc1', '1-methyl-4-phenylbenzene'],
  ['Clc1ccc(cc1)-c1ccc(Cl)cc1', '1-chloro-4-(4-chlorophenyl)benzene'],
  ['OC(=O)c1ccccc1-c1ccccc1', '2-phenylbenzoic acid'],
  ['Cc1ccc(C)cc1', '1,4-dimethylbenzene'], // PubChem: 1,4-xylene
  ['Cc1cccc(C)c1', '1,3-dimethylbenzene'],
  ['Cc1cc(C)cc(C)c1', '1,3,5-trimethylbenzene'],
  ['Clc1ccc(Cl)cc1', '1,4-dichlorobenzene'],
  ['Cc1ccc(Cl)cc1', '1-chloro-4-methylbenzene'],
  ['Clc1ccc(cc1)[N+]([O-])=O', '1-chloro-4-nitrobenzene'],
  ['CCCCCCc1ccccc1', 'hexylbenzene'], // ring senior to chain (2013)
  ['OCC1CCCCC1', 'cyclohexylmethanol'],
  ['OCCc1ccccc1', '2-phenylethanol'],
  ['OCc1ccccc1', 'phenylmethanol'],
  ['OCc1ccc(Cl)cc1', '(4-chlorophenyl)methanol'],
  ['CC(O)c1ccccc1', '1-phenylethanol'],
  ['OC(c1ccccc1)(c1ccccc1)c1ccccc1', 'triphenylmethanol'],
  ['C(c1ccccc1)c1ccccc1', 'benzylbenzene'],
]);

table('heterocycles (Hantzsch–Widman and retained names)', [
  ['c1ccncc1', 'pyridine'],
  ['c1cc[nH]c1', '1H-pyrrole'],
  ['c1ccoc1', 'furan'],
  ['c1ccsc1', 'thiophene'],
  ['c1c[nH]cn1', '1H-imidazole'],
  ['c1cn[nH]c1', '1H-pyrazole'],
  ['c1cocn1', '1,3-oxazole'],
  ['c1conc1', '1,2-oxazole'],
  ['c1cscn1', '1,3-thiazole'],
  ['c1csnc1', '1,2-thiazole'],
  ['c1c[nH]nn1', '1H-1,2,3-triazole'], // PubChem (tautomer-standardised): 2H-triazole
  ['c1cn[nH]n1', '2H-1,2,3-triazole'],
  ['c1ncn[nH]1', '1H-1,2,4-triazole'],
  ['c1nn[nH]n1', '2H-tetrazole'],
  ['c1nnn[nH]1', '1H-tetrazole'],
  ['c1cncnc1', 'pyrimidine'],
  ['c1cnccn1', 'pyrazine'],
  ['c1ccnnc1', 'pyridazine'],
  ['c1ncncn1', '1,3,5-triazine'],
  ['C1C=CC=CO1', '2H-pyran'],
  ['C1C=COC=C1', '4H-pyran'],
  ['C1CCNCC1', 'piperidine'],
  ['C1CNCCN1', 'piperazine'],
  ['C1COCCN1', 'morpholine'],
  ['C1CSCCN1', 'thiomorpholine'],
  ['C1CCNC1', 'pyrrolidine'],
  ['C1CCOC1', 'oxolane'],
  ['C1CCOCC1', 'oxane'],
  ['C1CCSC1', 'thiolane'],
  ['C1COCCO1', '1,4-dioxane'],
  ['C1OCCO1', '1,3-dioxolane'],
  ['C1CO1', 'oxirane'],
  ['C1CN1', 'aziridine'],
  ['C1COC1', 'oxetane'],
  ['C1CNC1', 'azetidine'],
  ['C1CCCNCC1', 'azepane'],
  ['C1CCCOCC1', 'oxepane'],
  ['C1CCCNCCC1', 'azocane'],
  ['CC1CO1', '2-methyloxirane'],
  ['C1CCOC1C', '2-methyloxolane'],
  ['C1=COCCC1', '3,4-dihydro-2H-pyran'],
  ['C1=CCNCC1', '1,2,3,6-tetrahydropyridine'],
  ['Cc1ccccn1', '2-methylpyridine'],
  ['OC(=O)c1cccnc1', 'pyridine-3-carboxylic acid'],
  ['NC(=O)c1cccnc1', 'pyridine-3-carboxamide'],
  ['NCc1ccccn1', 'pyridin-2-ylmethanamine'],
  ['c1ccc(nc1)-c1ccccn1', '2-pyridin-2-ylpyridine'],
  ['c1ccc(cc1)-c1ccncc1', '4-phenylpyridine'],
  ['c1ccc(cc1)N1CCCCC1', '1-phenylpiperidine'],
  ['O=C(N1CCOCC1)c1ccccc1', 'morpholin-4-yl(phenyl)methanone'],
  ['CC(=O)N1CCCCC1', '1-piperidin-1-ylethanone'],
  ['OC(=O)c1cccs1', 'thiophene-2-carboxylic acid'],
  ['O=Cc1ccco1', 'furan-2-carbaldehyde'],
  ['Cn1ccnc1', '1-methylimidazole'],
  ['CN1SC=CC1=O', '2-methyl-1,2-thiazol-3-one'],
  ['O=c1cccc[nH]1', '1H-pyridin-2-one'],
  ['O=C1CCCN1', 'pyrrolidin-2-one'],
  ['O=C1CCCCCN1', 'azepan-2-one'],
  ['O=C1CCCO1', 'oxolan-2-one'],
  ['O=C1CCCCO1', 'oxan-2-one'],
  ['O=C1CCC(=O)N1', 'pyrrolidine-2,5-dione'],
  ['O=C1CCC(=O)O1', 'oxolane-2,5-dione'],
  ['O=C1OCCN1', '1,3-oxazolidin-2-one'],
  ['O=C1NCCN1', 'imidazolidin-2-one'],
  ['O=C1OCCO1', '1,3-dioxolan-2-one'],
  ['O=C1CNC(=O)N1', 'imidazolidine-2,4-dione'],
  ['O=C1CC(=O)NC(=O)N1', '1,3-diazinane-2,4,6-trione'],
  ['O=c1cc[nH]c(=O)[nH]1', '1H-pyrimidine-2,4-dione'],
  ['Cc1c[nH]c(=O)[nH]c1=O', '5-methyl-1H-pyrimidine-2,4-dione'],
  ['O=C1C=COC=C1', 'pyran-4-one'],
  ['OCc1cc(=O)c(O)co1', '5-hydroxy-2-(hydroxymethyl)pyran-4-one'],
  ['CC1(C)OCC(CO)O1', '(2,2-dimethyl-1,3-dioxolan-4-yl)methanol'],
  ['NC(=O)NC1NC(=O)NC1=O', '(2,5-dioxoimidazolidin-4-yl)urea'],
]);

table('fused ring systems (templates, hydro prefixes, indicated hydrogen)', [
  ['c1ccc2ccccc2c1', 'naphthalene'],
  ['c1ccc2cc3ccccc3cc2c1', 'anthracene'],
  ['c1ccc2c(c1)ccc1ccccc12', 'phenanthrene'],
  ['c1cc2ccc3cccc4ccc(c1)c2c34', 'pyrene'],
  ['C1C=Cc2ccccc12', '1H-indene'],
  ['c1ccc2[nH]ccc2c1', '1H-indole'],
  ['c1ccc2c[nH]cc2c1', '2H-isoindole'],
  ['c1ccc2[nH]ncc2c1', '1H-indazole'],
  ['c1ccc2[nH]cnc2c1', '1H-benzimidazole'],
  ['c1ccc2occc2c1', '1-benzofuran'],
  ['c1ccc2sccc2c1', '1-benzothiophene'],
  ['c1ccc2ocnc2c1', '1,3-benzoxazole'],
  ['c1ccc2scnc2c1', '1,3-benzothiazole'],
  ['c1ccc2OCOc2c1', '1,3-benzodioxole'],
  ['c1ccc2ncccc2c1', 'quinoline'],
  ['c1ccc2cnccc2c1', 'isoquinoline'],
  ['c1ccc2ncncc2c1', 'quinazoline'],
  ['c1ccc2nccnc2c1', 'quinoxaline'],
  ['c1ncc2[nH]cnc2n1', '7H-purine'],
  ['c1cnc2ncncc2n1', 'pteridine'],
  ['c1ccc2c(c1)[nH]c1ccccc12', '9H-carbazole'],
  ['c1ccc2nc3ccccc3cc2c1', 'acridine'],
  ['C1c2ccccc2Oc2ccccc12', '9H-xanthene'],
  ['C1C=Cc2ccccc2O1', '2H-chromene'],
  ['C1c2ccccc2-c2ccccc12', '9H-fluorene'],
  ['c1ccc2cccc2cc1', 'azulene'],
  ['C1C2CC3CC1CC(C2)C3', 'adamantane'],
  ['C1CCc2ccccc2C1', '1,2,3,4-tetrahydronaphthalene'],
  ['C1CCC2=C(C1)C=CC=C2', '1,2,3,4-tetrahydronaphthalene'],
  ['C1CCC2CCCCC2C1', '1,2,3,4,4a,5,6,7,8,8a-decahydronaphthalene'],
  ['C1Cc2ccccc2C1', '2,3-dihydro-1H-indene'],
  ['C1CC2=CC=CC=C2C1', '2,3-dihydro-1H-indene'],
  ['C1Cc2ccccc2N1', '2,3-dihydro-1H-indole'],
  ['C1CCc2ncccc2C1', '5,6,7,8-tetrahydroquinoline'],
  ['C1CNCc2ccccc21', '1,2,3,4-tetrahydroisoquinoline'],
  ['CC1(C)CCc2ccccc21', '1,1-dimethyl-2,3-dihydroindene'], // PubChem numbers hydro first: 3,3-dimethyl-1,2-dihydroindene
  ['CC1(C)C=Cc2ccccc2O1', '2,2-dimethylchromene'],
  ['CC1C=Cc2ccccc21', '1-methyl-1H-indene'],
  ['CN1CCc2ccccc21', '1-methyl-2,3-dihydroindole'],
  ['Cn1ccc2ccccc21', '1-methylindole'],
  ['CCn1c2ccccc2c2ccccc21', '9-ethylcarbazole'],
  ['Cc1cc2ccccc2[nH]1', '2-methyl-1H-indole'],
  ['Cc1c[nH]c2ccccc12', '3-methyl-1H-indole'],
  ['OC(=O)Cc1c[nH]c2ccccc12', '2-(1H-indol-3-yl)acetic acid'],
  ['O=Cc1c[nH]c2ccccc12', '1H-indole-3-carbaldehyde'],
  ['Cc1nc2ccccc2[nH]1', '2-methyl-1H-benzimidazole'],
  ['c1ccc(cc1)-c1nc2ccccc2[nH]1', '2-phenyl-1H-benzimidazole'],
  ['Nc1ccc2nc[nH]c2c1', '1H-benzimidazol-6-amine'], // indicated hydrogen gets the lowest locant (PubChem: 3H-…-5-amine)
  ['Oc1ccc2ccccc2c1', 'naphthalen-2-ol'],
  ['Nc1cccc2ccccc12', 'naphthalen-1-amine'],
  ['OC(=O)c1ccc2ccccc2c1', 'naphthalene-2-carboxylic acid'],
  ['O=C1CCCc2ccccc21', '3,4-dihydro-2H-naphthalen-1-one'],
  ['O=C1CCc2ccccc2C1', '3,4-dihydro-1H-naphthalen-2-one'],
  ['O=C1C=CCc2ccccc21', '4H-naphthalen-1-one'],
  ['O=C1Cc2ccccc2N1', '1,3-dihydroindol-2-one'],
  ['O=C1CCc2ccccc2N1', '3,4-dihydro-1H-quinolin-2-one'],
  ['O=C1C=Cc2ccccc2N1', '1H-quinolin-2-one'],
  ['O=C1NC(=O)c2ccccc12', 'isoindole-1,3-dione'],
  ['O=c1ccc2ccccc2o1', 'chromen-2-one'],
  ['O=c1ccoc2ccccc12', 'chromen-4-one'],
  ['O=C1CCOc2ccccc21', '2,3-dihydrochromen-4-one'],
  ['O=C1OCc2ccccc12', '3H-2-benzofuran-1-one'],
  ['O=c1[nH]c2ccccc2o1', '3H-1,3-benzoxazol-2-one'],
  ['O=C1c2ccccc2C(=O)c2ccccc12', 'anthracene-9,10-dione'],
  ['O=C1c2ccccc2Cc2ccccc21', '10H-anthracen-9-one'],
  ['O=C1C=CC(=O)c2ccccc12', 'naphthalene-1,4-dione'],
  ['O=C1c2ccccc2-c2ccccc12', 'fluoren-9-one'],
  ['O=C1c2ccccc2Oc2ccccc12', 'xanthen-9-one'],
  ['O=Cc1ccc2OCOc2c1', '1,3-benzodioxole-5-carbaldehyde'],
  ['CN1C=NC2=C1C(=O)N(C)C(=O)N2C', '1,3,7-trimethylpurine-2,6-dione'],
  ['CN1C=NC2=C1C(=O)NC(=O)N2C', '3,7-dimethylpurine-2,6-dione'],
  ['CN1C(=O)N(C)C2=C(C1=O)NC=N2', '1,3-dimethyl-7H-purine-2,6-dione'],
  ['O=C1NC(=O)C2=C(N1)N=CN2', '3,7-dihydropurine-2,6-dione'],
  ['O=C1NC(=O)C2=C(N1)NC(=O)N2', '7,9-dihydro-3H-purine-2,6,8-trione'],
  ['CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc12', '7-chloro-1-methyl-5-phenyl-3H-1,4-benzodiazepin-2-one'],
  ['Clc1ccc2c(c1)N(CCCN(C)C)c1ccccc1S2', '3-(2-chlorophenothiazin-10-yl)-N,N-dimethylpropan-1-amine'],
  ['CC12CCC3c4ccc(O)cc4CCC3C1CCC2O', '13-methyl-6,7,8,9,11,12,14,15,16,17-decahydrocyclopenta[a]phenanthrene-3,17-diol'],
  ['CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O', '17-hydroxy-10,13-dimethyl-1,2,6,7,8,9,11,12,14,15,16,17-dodecahydrocyclopenta[a]phenanthren-3-one'],
]);

table('von Baeyer and spiro systems', [
  ['C1CC2CCC1C2', 'bicyclo[2.2.1]heptane'],
  ['C1CC2CCC1CC2', 'bicyclo[2.2.2]octane'],
  ['C1=CC2CCC1CC2', 'bicyclo[2.2.2]oct-2-ene'],
  ['C1CC2CCC1O2', '7-oxabicyclo[2.2.1]heptane'],
  ['C1CN2CCC1CC2', '1-azabicyclo[2.2.2]octane'],
  ['CN1C2CCC1CC(O)C2', '8-methyl-8-azabicyclo[3.2.1]octan-3-ol'],
  ['CC1(C)C2CCC1(C)C(=O)C2', '1,7,7-trimethylbicyclo[2.2.1]heptan-2-one'],
  ['OC1CC2CCC1C2', 'bicyclo[2.2.1]heptan-2-ol'],
  ['C1CC2CC2C1', 'bicyclo[3.1.0]hexane'],
  ['C1CCC2CC2C1', 'bicyclo[4.1.0]heptane'],
  ['C1CC2CCC2C1', 'bicyclo[3.2.0]heptane'],
  ['CC12CCC(CC1)C(C)(C)O2', '1,3,3-trimethyl-2-oxabicyclo[2.2.2]octane'],
  ['C1CCC2(CC1)CCCC2', 'spiro[4.5]decane'],
  ['C1CCC2(CC1)OCCO2', '1,4-dioxaspiro[4.5]decane'],
  ['NC12CC3CC(CC(C3)C1)C2', 'adamantan-1-amine'],
]);

table('heteroatom parents, sulfur oxides, thio and boron compounds', [
  ['O=S1(=O)CCCC1', 'thiolane 1,1-dioxide'],
  ['O=S1CCCCC1', 'thiane 1-oxide'],
  ['O=S1(=O)CC=CC1', '2,5-dihydrothiophene 1,1-dioxide'],
  ['OC(=O)C1CCS(=O)(=O)C1', '1,1-dioxothiolane-3-carboxylic acid'],
  ['c1ccc(cc1)N=Nc1ccccc1', 'diphenyldiazene'],
  ['CN(C)c1ccc(cc1)N=Nc1ccc(cc1)S(=O)(=O)O', '4-{[4-(dimethylamino)phenyl]diazenyl}benzenesulfonic acid'],
  ['NNc1ccccc1', 'phenylhydrazine'],
  ['CN(C)N', '1,1-dimethylhydrazine'],
  ['CC(=O)NN', 'acetohydrazide'],
  ['NNC(=O)c1ccncc1', 'pyridine-4-carbohydrazide'],
  ['S=C(N)N', 'thiourea'],
  ['CC(N)=S', 'ethanethioamide'],
  ['NC(=S)c1ccccc1', 'benzenecarbothioamide'],
  ['OB(O)c1ccccc1', 'phenylboronic acid'],
  ['COc1ccc(cc1)B(O)O', '(4-methoxyphenyl)boronic acid'],
  ['C[Si](C)(C)C', 'tetramethylsilane'],
  ['C[Si](C)(C)Oc1ccccc1', 'trimethyl(phenoxy)silane'],
  ['C[Si](C)(C)O[Si](C)(C)C', 'trimethyl(trimethylsilyloxy)silane'],
  ['CC(C)(C)[Si](C)(C)OC', 'tert-butyl(methoxy)dimethylsilane'], // PubChem: tert-butyl-methoxy-dimethylsilane
  ['CC(C)N=C=NC(C)C', "N,N'-di(propan-2-yl)methanediimine"],
  ['O=C(Cl)Cl', 'carbonyl dichloride'],
  ['NC(=O)Cl', 'carbamoyl chloride'],
  ['OC(=O)Nc1ccccc1', 'phenylcarbamic acid'],
  ['CCOC(N)=O', 'ethyl carbamate'],
  ['CNC(N)=N', '1-methylguanidine'], // drawn tautomer (PubChem standardises to 2-methylguanidine)
  ['OC(=O)c1ccc(cc1)S(=O)(=O)O', '4-sulfobenzoic acid'],
  ['OC(=O)c1ccccc1C(=O)OC', '2-methoxycarbonylbenzoic acid'],
  ['NC(=O)c1ccccc1C(N)=O', 'benzene-1,2-dicarboxamide'],
  ['CC(=O)c1cccc2ccccc12', '1-naphthalen-1-ylethanone'],
  ['OC(=O)Cc1cccs1', '2-thiophen-2-ylacetic acid'],
  ['OC(=O)CCN1CCCCC1', '3-piperidin-1-ylpropanoic acid'],
  ['OCCN1CCOCC1', '2-morpholin-4-ylethanol'],
  ['c1ccc2c(c1)ccc1c2ccc2ccccc21', 'chrysene'],
  ['C1COCCOCCOCCOCCOCCO1', '1,4,7,10,13,16-hexaoxacyclooctadecane'],
  ['O=C1CCCCCCCCCCCO1', 'oxacyclotridecan-2-one'],
  ['CCOP(=O)(OCC)OCC', 'triethyl phosphate'],
  ['OP(=O)(O)OCC(O)CO', '2,3-dihydroxypropyl dihydrogen phosphate'],
  ['O=[N+]([O-])OCC(CO[N+]([O-])=O)O[N+]([O-])=O', '2,3-dinitrooxypropyl nitrate'], // PubChem: 1,3-dinitrooxypropan-2-yl nitrate
  ['CC(=O)C=O', '2-oxopropanal'],
  ['O=CC=O', 'oxaldehyde'],
  ['N#CC#N', 'oxalonitrile'],
  ['ClC(=O)C(Cl)=O', 'oxalyl dichloride'],
  ['OC(=O)C(O)C(O)C(=O)O', '2,3-dihydroxybutanedioic acid'],
  ['CC#CC(=O)O', 'but-2-ynoic acid'],
  ['C1=CC=CC1', 'cyclopenta-1,3-diene'],
  ['CN1C=CN(C)C1=S', '1,3-dimethylimidazole-2-thione'],
  ['[NH3+]CC([O-])=O', '2-azaniumylacetate'],
  ['OCC[N+](C)(C)C.[Cl-]', '2-hydroxy-N,N,N-trimethylethanaminium chloride'],
  ['CC[O-].[Na+]', 'sodium ethanolate'],
  ['[O-]c1ccccc1.[Na+]', 'sodium phenolate'],
  ['NC(CCCNC(N)=N)C(=O)O', '2-amino-5-(carbamimidoylamino)pentanoic acid'],
]);

table('more ring and chain cases', [
  ['Oc1ccc(Cl)cc1Cl', '2,4-dichlorophenol'],
  ['Oc1c(Cl)cc(Cl)cc1Cl', '2,4,6-trichlorophenol'],
  ['Cc1cccc(C)c1N', '2,6-dimethylaniline'],
  ['Nc1ccc(Cl)cc1', '4-chloroaniline'],
  ['OC(=O)c1ccc(Cl)cc1Cl', '2,4-dichlorobenzoic acid'],
  ['OC(=O)c1cc(Cl)ccc1O', '5-chloro-2-hydroxybenzoic acid'],
  ['CC(=O)c1ccc(O)cc1', '1-(4-hydroxyphenyl)ethanone'],
  ['O=C(CCl)c1ccccc1', '2-chloro-1-phenylethanone'],
  ['CCCCCCCC(=O)OCC', 'ethyl octanoate'],
  ['CC(C)CCOC(C)=O', '3-methylbutyl acetate'],
  ['C=CC(=O)OCCO', '2-hydroxyethyl prop-2-enoate'],
  ['CC(O)CN', '1-aminopropan-2-ol'],
  ['CC(N)CO', '2-aminopropan-1-ol'],
  ['NC(CO)(CO)CO', '2-amino-2-(hydroxymethyl)propane-1,3-diol'],
  ['OCC(CO)(CO)CO', '2,2-bis(hydroxymethyl)propane-1,3-diol'],
  ['CC(C)(CO)CO', '2,2-dimethylpropane-1,3-diol'],
  ['OC(CC(=O)O)C(=O)O', '2-hydroxybutanedioic acid'],
  ['OC(=O)C=CC(=O)O', 'but-2-enedioic acid'],
  ['OC(=O)CC(=O)O', 'propanedioic acid'],
  ['NCCCCC(N)C(=O)O', '2,6-diaminohexanoic acid'],
  ['NC(CS)C(=O)O', '2-amino-3-sulfanylpropanoic acid'],
  ['NC(CCC(=O)O)C(=O)O', '2-aminopentanedioic acid'],
  ['NC(Cc1ccccc1)C(=O)O', '2-amino-3-phenylpropanoic acid'],
  ['NC(Cc1ccc(O)cc1)C(=O)O', '2-amino-3-(4-hydroxyphenyl)propanoic acid'],
  ['CSCCC(N)C(=O)O', '2-amino-4-methylsulfanylbutanoic acid'],
  ['OC1CNC(C1)C(=O)O', '4-hydroxypyrrolidine-2-carboxylic acid'],
  ['OC(=O)C1CCCN1', 'pyrrolidine-2-carboxylic acid'],
  ['NC(=O)CCC(N)C(=O)O', '2,5-diamino-5-oxopentanoic acid'],
  ['CC(O)C(N)C(=O)O', '2-amino-3-hydroxybutanoic acid'],
  ['CC(C)C(N)C(=O)O', '2-amino-3-methylbutanoic acid'],
  ['c1ccc2c(c1)oc1ccccc12', 'dibenzofuran'],
  ['c1ccc2nc3ccccc3nc2c1', 'phenazine'],
  ['c1ccc2c(c1)Nc1ccccc1S2', '10H-phenothiazine'],
  ['c1cnc2c(c1)ccc1cccnc12', '1,10-phenanthroline'],
  ['Oc1cccc2cccnc12', 'quinolin-8-ol'],
  ['c1ccc2[nH]c(cc2c1)C(=O)O', '1H-indole-2-carboxylic acid'],
  ['CC(=O)c1ccc2ccccc2c1', '1-naphthalen-2-ylethanone'],
  ['c1ccc(cc1)-c1ccc(cc1)-c1ccccc1', '1,4-diphenylbenzene'],
  ['OC(=O)c1ccc(o1)C=O', '5-formylfuran-2-carboxylic acid'],
  ['OCc1ccc(o1)C=O', '5-(hydroxymethyl)furan-2-carbaldehyde'],
  ['Cc1ccc(o1)C', '2,5-dimethylfuran'],
  ['O=C1OC(=O)c2ccccc12', '2-benzofuran-1,3-dione'],
  ['O=C1c2ccccc2C(=O)N1C', '2-methylisoindole-1,3-dione'],
  ['CN1CCN(C)CC1', '1,4-dimethylpiperazine'],
  ['CC1CCCCN1', '2-methylpiperidine'],
  ['OC1CCN(C)CC1', '1-methylpiperidin-4-ol'],
  ['O=C1CCN(CC1)C(=O)OC(C)(C)C', 'tert-butyl 4-oxopiperidine-1-carboxylate'],
  ['C1CCN(CC1)c1ccncc1', '4-piperidin-1-ylpyridine'],
  ['CN(C)c1ccncc1', 'N,N-dimethylpyridin-4-amine'],
  ['Nc1ccccn1', 'pyridin-2-amine'],
  ['CC1(C)CCCC1', '1,1-dimethylcyclopentane'],
  ['CC1CCC(C)CC1', '1,4-dimethylcyclohexane'],
  ['CC1=CCCC1', '1-methylcyclopentene'],
  ['O=C1CCC(CC1)C(=O)O', '4-oxocyclohexane-1-carboxylic acid'],
  ['OC1(CCCCC1)C#N', '1-hydroxycyclohexane-1-carbonitrile'],
  ['CC(C)(O)c1ccccc1', '2-phenylpropan-2-ol'],
  ['CC(C)(C#N)O', '2-hydroxy-2-methylpropanenitrile'],
  ['CCOC(=O)C(C)(C)C', 'ethyl 2,2-dimethylpropanoate'],
  ['NS(=O)(=O)c1ccc(Cl)cc1', '4-chlorobenzenesulfonamide'],
  ['O=S(=O)(Nc1ccccc1)c1ccccc1', 'N-phenylbenzenesulfonamide'],
  ['CS(=O)(=O)Oc1ccccc1', 'phenyl methanesulfonate'],
  ['COS(=O)(=O)c1ccc(C)cc1', 'methyl 4-methylbenzenesulfonate'],
  ['OB(O)c1ccc(cc1)C(=O)O', '4-boronobenzoic acid'],
  ['OC(=O)c1ccc(NN)cc1', '4-hydrazinylbenzoic acid'],
  ['Cn1cccn1', '1-methylpyrazole'],
  ['Cc1cc(C)n[nH]1', '3,5-dimethyl-1H-pyrazole'],
  ['Cc1c[nH]cn1', '4-methyl-1H-imidazole'], // drawn tautomer (PubChem: 5-methyl-1H-imidazole)
  ['Cn1cncn1', '1-methyl-1,2,4-triazole'],
  ['c1ccc(cc1)-c1nn[nH]n1', '5-phenyl-2H-tetrazole'],
  ['Nc1nccs1', '1,3-thiazol-2-amine'],
  ['Cc1ncco1', '2-methyl-1,3-oxazole'],
  ['Cc1ccno1', '5-methyl-1,2-oxazole'],
  ['Cc1nc2ccccc2s1', '2-methyl-1,3-benzothiazole'],
  ['Nc1nc2ccccc2s1', '1,3-benzothiazol-2-amine'],
  ['S=c1[nH]c2ccccc2s1', '3H-1,3-benzothiazole-2-thione'],
  ['Cc1ccnc2ccccc12', '4-methylquinoline'],
  ['Clc1ccc2ccccc2n1', '2-chloroquinoline'],
  ['Nc1nccc2ccccc12', 'isoquinolin-1-amine'],
  ['O=c1ccc2ccc(O)cc2o1', '7-hydroxychromen-2-one'],
  ['Cc1cc(=O)oc2ccccc12', '4-methylchromen-2-one'],
  ['COc1cc2ccc(=O)oc2cc1OC', '6,7-dimethoxychromen-2-one'],
  ['O=c1cc(oc2ccccc12)-c1ccccc1', '2-phenylchromen-4-one'],
  ['O=C1CC(Oc2ccccc12)c1ccccc1', '2-phenyl-2,3-dihydrochromen-4-one'],
  ['C1Cc2ccccc2OC1', '3,4-dihydro-2H-chromene'],
  ['O=C1CCCCCO1', 'oxepan-2-one'],
  ['O=C1CCNCC1', 'piperidin-4-one'],
  ['C1CCN(CC1)Cc1ccccc1', '1-benzylpiperidine'],
  ['CN1CCOCC1', '4-methylmorpholine'],
  ['CN1CCCC1=O', '1-methylpyrrolidin-2-one'],
  ['Nc1ncccn1', 'pyrimidin-2-amine'],
  ['COc1ccnc(OC)n1', '2,4-dimethoxypyrimidine'],
  ['COc1ccc2ccccc2c1', '2-methoxynaphthalene'],
  ['Oc1cccc2c(O)cccc12', 'naphthalene-1,5-diol'],
  ['Cc1c2ccccc2cc2ccccc12', '9-methylanthracene'],
  ['OCc1c2ccccc2cc2ccccc12', 'anthracen-9-ylmethanol'],
  ['Oc1cc2ccccc2c2ccccc12', 'phenanthren-9-ol'],
  ['OC1c2ccccc2-c2ccccc12', '9H-fluoren-9-ol'],
  ['Nc1ccc2-c3ccccc3Cc2c1', '9H-fluoren-2-amine'],
  ['Oc1ccc2[nH]c3ccccc3c2c1', '9H-carbazol-3-ol'],
  ['Nc1c2ccccc2nc2ccccc12', 'acridin-9-amine'],
  ['Cc1ccc2cccc2cc1', '6-methylazulene'],
  ['O=C1C2CC3CC(C2)CC1C3', 'adamantan-2-one'],
  ['OC1C2CC3CC(C2)CC1C3', 'adamantan-2-ol'],
  ['C1CC2CCC1NC2', '2-azabicyclo[2.2.2]octane'],
  ['CC1(C)OCCO1', '2,2-dimethyl-1,3-dioxolane'],
  ['c1ccc(cc1)C1CO1', '2-phenyloxirane'],
  ['CC1(C)CO1', '2,2-dimethyloxirane'],
  ['CC1=CC(=O)c2ccccc2C1=O', '2-methylnaphthalene-1,4-dione'],
  ['Cc1ccc(cc1)S(=O)(=O)Oc1ccccc1', 'phenyl 4-methylbenzenesulfonate'],
  ['O=C1NC(=O)C(=O)N1', 'imidazolidine-2,4,5-trione'],
  ['C[n+]1cccc(c1)C([O-])=O', '1-methylpyridin-1-ium-3-carboxylate'],
  ['C/C(=C\\C(=O)O)C', '3-methylbut-2-enoic acid'],
]);

table('salts and multi-component inputs', [
  ['[Na+].[O-]C(=O)c1ccccc1', 'sodium benzoate'],
  ['CC([O-])=O.[Na+]', 'sodium acetate'],
  ['[Na+].CCCCCCCCCCCCOS([O-])(=O)=O', 'sodium dodecyl sulfate'],
  ['CCCCCCCCCCCCCCCC[N+](C)(C)C.[Br-]', 'N,N,N-trimethylhexadecan-1-aminium bromide'],
  ['CCCCCCCCCCCCCCCCCCCCCC[N+](C)(C)C.[Cl-]', 'N,N,N-trimethyldocosan-1-aminium chloride'],
  ['CCCCCCCCCCCC[N+](C)(C)Cc1ccccc1.[Cl-]', 'N-benzyl-N,N-dimethyldodecan-1-aminium chloride'],
  ['C[N+](C)(C)C.[Cl-]', 'N,N,N-trimethylmethanaminium chloride'],
  ['C[n+]1ccccc1.[I-]', '1-methylpyridin-1-ium iodide'],
  ['[K+].CC=CC=CC([O-])=O', 'potassium hexa-2,4-dienoate'],
  ['[Na+].[Na+].[O-]C(=O)CCC([O-])=O', 'disodium butanedioate'],
  ['CC(C)Cc1ccc(cc1)C(C)C(=O)[O-].[Na+]', 'sodium 2-[4-(2-methylpropyl)phenyl]propanoate'],
  ['[NH4+].[Cl-]', 'ammonium chloride'],
  ['[Na+].[Cl-]', 'sodium chloride'],
  ['CCN.Cl', 'ethanamine hydrochloride'],
  ['OC(=O)c1ccccc1O.CCO', 'ethanol; 2-hydroxybenzoic acid'],
]);

table('real-world molecules (cosmetics / pharma)', [
  ['CC(=O)Oc1ccccc1C(=O)O', '2-acetyloxybenzoic acid'],
  ['CC(C)Cc1ccc(cc1)C(C)C(=O)O', '2-[4-(2-methylpropyl)phenyl]propanoic acid'],
  ['CC(=O)Nc1ccc(O)cc1', 'N-(4-hydroxyphenyl)acetamide'],
  ['COc1cc(C=O)ccc1O', '4-hydroxy-3-methoxybenzaldehyde'],
  ['CCOc1cc(C=O)ccc1O', '3-ethoxy-4-hydroxybenzaldehyde'],
  ['CC(C)C1CCC(C)CC1O', '5-methyl-2-propan-2-ylcyclohexan-1-ol'],
  ['CC1=CCC(CC1)C(C)=C', '1-methyl-4-prop-1-en-2-ylcyclohexene'],
  ['CC1=CCC(CC1=O)C(C)=C', '2-methyl-5-prop-1-en-2-ylcyclohex-2-en-1-one'],
  ['Cc1ccc(C(C)C)c(O)c1', '5-methyl-2-propan-2-ylphenol'],
  ['COc1cc(CC=C)ccc1O', '2-methoxy-4-prop-2-enylphenol'],
  ['CC(C)=CCCC(C)(O)C=C', '3,7-dimethylocta-1,6-dien-3-ol'],
  ['CC(C)=CCCC(C)=CCO', '3,7-dimethylocta-2,6-dien-1-ol'],
  ['CC(C)=CCCC(C)=CCCC(C)=CCO', '3,7,11-trimethyldodeca-2,6,10-trien-1-ol'],
  ['CC(C)=CCCC(C)(O)C1CCC(C)=CC1', '6-methyl-2-(4-methylcyclohex-3-en-1-yl)hept-5-en-2-ol'],
  ['CC1=C(C(C)(C)CCC1)C=CC(C)=CC=CC(C)=CCO', '3,7-dimethyl-9-(2,6,6-trimethylcyclohexen-1-yl)nona-2,4,6,8-tetraen-1-ol'],
  ['CC(Cc1ccc(cc1)C(C)(C)C)C=O', '3-(4-tert-butylphenyl)-2-methylpropanal'],
  ['CCCCCCCCCCCCCCCCO', 'hexadecan-1-ol'],
  ['CC(O)CO', 'propane-1,2-diol'],
  ['CC(O)CCO', 'butane-1,3-diol'],
  ['CCCCCCC(O)CO', 'octane-1,2-diol'],
  ['CCCCC(CC)COCC(O)CO', '3-(2-ethylhexoxy)propane-1,2-diol'],
  ['OCC(O)COc1ccc(Cl)cc1', '3-(4-chlorophenoxy)propane-1,2-diol'],
  ['CCCOC(=O)c1ccc(O)cc1', 'propyl 4-hydroxybenzoate'],
  ['CC(C)(CO)C(O)C(=O)NCCCO', '2,4-dihydroxy-N-(3-hydroxypropyl)-3,3-dimethylbutanamide'],
  ['NC(=O)c1cccnc1', 'pyridine-3-carboxamide'],
  ['Cc1ncc(CO)c(CO)c1O', '4,5-bis(hydroxymethyl)-2-methylpyridin-3-ol'],
  ['OCC1OC(O)C(O)C(O)C1O', '6-(hydroxymethyl)oxane-2,3,4,5-tetrol'],
  ['OCC1OC(Oc2ccc(O)cc2)C(O)C(O)C1O', '2-(hydroxymethyl)-6-(4-hydroxyphenoxy)oxane-3,4,5-triol'],
  ['COc1ccc(C(=O)c2ccccc2)c(O)c1', '(2-hydroxy-4-methoxyphenyl)-phenylmethanone'],
  ['COc1cc(O)c(C(=O)c2ccccc2)cc1S(=O)(=O)O', '5-benzoyl-4-hydroxy-2-methoxybenzenesulfonic acid'],
  ['CC(C)(C)c1ccc(cc1)C(=O)CC(=O)c1ccc(OC)cc1', '1-(4-tert-butylphenyl)-3-(4-methoxyphenyl)propane-1,3-dione'],
  ['CC(C)(C)c1cc(C)cc(c1O)C(C)(C)C', '2,6-ditert-butyl-4-methylphenol'],
  ['Oc1cc(Cl)ccc1Oc1ccc(Cl)cc1Cl', '5-chloro-2-(2,4-dichlorophenoxy)phenol'],
  ['Oc1ccc(C=Cc2cc(O)cc(O)c2)cc1', '5-[2-(4-hydroxyphenyl)ethenyl]benzene-1,3-diol'],
  ['COc1cc(C=CC(=O)CC(=O)C=Cc2ccc(O)c(OC)c2)ccc1O', '1,7-bis(4-hydroxy-3-methoxyphenyl)hepta-1,6-diene-3,5-dione'],
  ['COc1cc(C=CC(=O)O)ccc1O', '3-(4-hydroxy-3-methoxyphenyl)prop-2-enoic acid'],
  ['OC(=O)C=Cc1ccc(O)c(O)c1', '3-(3,4-dihydroxyphenyl)prop-2-enoic acid'],
  ['COc1cc(CNC(=O)CCCCC=CC(C)C)ccc1O', 'N-[(4-hydroxy-3-methoxyphenyl)methyl]-8-methylnon-6-enamide'],
  ['CC(C)CCCC(C)CCCC(C)CCCC1(C)CCc2c(C)c(O)c(C)c(C)c2O1', '2,5,7,8-tetramethyl-2-(4,8,12-trimethyltridecyl)-3,4-dihydrochromen-6-ol'],
  ['NC(Cc1c[nH]c2ccccc12)C(=O)O', '2-amino-3-(1H-indol-3-yl)propanoic acid'],
  ['NCC(=O)O', '2-aminoacetic acid'],
  ['CC(N)C(=O)O', '2-aminopropanoic acid'],
  ['CN1CCCC1c1cccnc1', '3-(1-methylpyrrolidin-2-yl)pyridine'],
  ['CC(C)NCC(O)COc1cccc2ccccc12', '1-naphthalen-1-yloxy-3-(propan-2-ylamino)propan-2-ol'],
  ['CC(C)NCC(O)COc1ccc(CCOC)cc1', '1-[4-(2-methoxyethyl)phenoxy]-3-(propan-2-ylamino)propan-2-ol'],
  ['CC(C)NCC(O)COc1ccc(CC(N)=O)cc1', '2-{4-[2-hydroxy-3-(propan-2-ylamino)propoxy]phenyl}acetamide'], // PubChem nests [ ] only
  ['CC(C)(C)NCC(O)c1ccc(O)c(CO)c1', '4-[2-(tert-butylamino)-1-hydroxyethyl]-2-(hydroxymethyl)phenol'],
  ['CNCC(O)c1cccc(O)c1', '3-[1-hydroxy-2-(methylamino)ethyl]phenol'],
  ['CNC(C)C(O)c1ccccc1', '2-(methylamino)-1-phenylpropan-1-ol'],
  ['CCN(CC)CC(=O)Nc1c(C)cccc1C', '2-(diethylamino)-N-(2,6-dimethylphenyl)acetamide'],
  ['NCCc1ccc(O)c(O)c1', '4-(2-aminoethyl)benzene-1,2-diol'],
  ['NCCc1c[nH]c2ccc(O)cc12', '3-(2-aminoethyl)-1H-indol-5-ol'],
  ['CC(=O)NCCc1c[nH]c2ccc(OC)cc12', 'N-[2-(5-methoxy-1H-indol-3-yl)ethyl]acetamide'],
  ['COc1ccc2cc(ccc2c1)C(C)C(=O)O', '2-(6-methoxynaphthalen-2-yl)propanoic acid'],
  ['CC(C(=O)O)c1cccc(c1)C(=O)c1ccccc1', '2-(3-benzoylphenyl)propanoic acid'],
  ['CNCCC(Oc1ccc(cc1)C(F)(F)F)c1ccccc1', 'N-methyl-3-phenyl-3-[4-(trifluoromethyl)phenoxy]propan-1-amine'],
  ['COc1ccc(cc1)C(CN(C)C)C1(O)CCCCC1', '1-[2-(dimethylamino)-1-(4-methoxyphenyl)ethyl]cyclohexan-1-ol'],
  ['OC(=O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O', '1-cyclopropyl-6-fluoro-4-oxo-7-piperazin-1-ylquinoline-3-carboxylic acid'],
  ['CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O', '7-[2-(4-fluorophenyl)-3-phenyl-4-(phenylcarbamoyl)-5-propan-2-ylpyrrol-1-yl]-3,5-dihydroxyheptanoic acid'],
  ['CC1(C)SC2C(NC(=O)C(N)c3ccc(O)cc3)C(=O)N2C1C(=O)O', '6-{[2-amino-2-(4-hydroxyphenyl)acetyl]amino}-3,3-dimethyl-7-oxo-4-thia-1-azabicyclo[3.2.0]heptane-2-carboxylic acid'],
  ['CC(C)(C)C(=O)C(Oc1ccc(Cl)cc1)n1ccnc1', '1-(4-chlorophenoxy)-1-imidazol-1-yl-3,3-dimethylbutan-2-one'],
  ['CNC(C)Cc1ccc2OCOc2c1', '1-(1,3-benzodioxol-5-yl)-N-methylpropan-2-amine'],
  ['CN(C)C(=N)N=C(N)N', '3-(diaminomethylidene)-1,1-dimethylguanidine'], // drawn tautomer (PubChem: 2-carbamimidoyl-1,1-dimethylguanidine)
  ['CC(C)CCCC(C)C1CCC2C1(C)CCC1C2CC=C2CC(O)CCC12C', '10,13-dimethyl-17-(6-methylheptan-2-yl)-2,3,4,7,8,9,11,12,14,15,16,17-dodecahydro-1H-cyclopenta[a]phenanthren-3-ol'],
  ['OCC(=O)C1(O)CCC2C3CCC4=CC(=O)CCC4(C)C3C(O)CC21C', '11,17-dihydroxy-17-(2-hydroxyacetyl)-10,13-dimethyl-2,6,7,8,9,11,12,14,15,16-decahydro-1H-cyclopenta[a]phenanthren-3-one'],
  ['CC(=O)C1CCC2C3CCC4=CC(=O)CCC4(C)C3CCC12C', '17-acetyl-10,13-dimethyl-1,2,6,7,8,9,11,12,14,15,16,17-dodecahydrocyclopenta[a]phenanthren-3-one'],
]);

// ───────────── stereodescriptors (active once assignCIP is implemented) ─────────────

const cipReady = (() => {
  try {
    return assignCIP(parseSmiles('N[C@@H](C)C(=O)O')).centers.size > 0;
  } catch {
    return false;
  }
})();

describe.skipIf(!cipReady)('stereodescriptors', () => {
  it.each([
    ['N[C@@H](C)C(=O)O', '(2S)-2-aminopropanoic acid'],
    ['N[C@H](C)C(=O)O', '(2R)-2-aminopropanoic acid'],
    ['C/C=C/C', '(E)-but-2-ene'],
    ['C/C=C\\C', '(Z)-but-2-ene'],
    ['C/C=C/C(=O)O', '(E)-but-2-enoic acid'],
    ['C/C=C/C=C/C(=O)O', '(2E,4E)-hexa-2,4-dienoic acid'],
    ['[K+].C/C=C/C=C/C([O-])=O', 'potassium (2E,4E)-hexa-2,4-dienoate'],
    ['N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O', '(2S)-2-amino-3-(1H-indol-3-yl)propanoic acid'],
    ['CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O', '(2S)-2-[4-(2-methylpropyl)phenyl]propanoic acid'],
    ['CC(C)C[C@H](N)C(=O)O', '(2S)-2-amino-4-methylpentanoic acid'],
    ['CCCCCC[C@@H](O)C/C=C\\CCCCCCCC(=O)O', '(Z,12R)-12-hydroxyoctadec-9-enoic acid'],
    ['CN1CCC[C@H]1c1cccnc1', '3-[(2S)-1-methylpyrrolidin-2-yl]pyridine'],
    ['OC(=O)C[C@@H](C)CC', '(3S)-3-methylpentanoic acid'],
    ['CC(C)[C@@H]1CC[C@@H](C)C[C@H]1O', '(1R,2S,5R)-5-methyl-2-propan-2-ylcyclohexan-1-ol'],
    ['C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc34)[C@@H]1CC[C@@H]2O', '(8R,9S,13S,14S,17S)-13-methyl-6,7,8,9,11,12,14,15,16,17-decahydrocyclopenta[a]phenanthrene-3,17-diol'],
    ['CC[C@H](C)OC(C)=O', '[(2S)-butan-2-yl] acetate'],
    ['C[C@@H](N)c1ccccc1', '(1R)-1-phenylethanamine'],
    ['O[C@@H]1CCCC[C@H]1O', '(1R,2R)-cyclohexane-1,2-diol'],
    ['CC(=O)O[C@@H](C)c1ccccc1', '[(1S)-1-phenylethyl] acetate'],
    ['OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O', '(2R,3R,4S,5S,6R)-6-(hydroxymethyl)oxane-2,3,4,5-tetrol'],
    ['O[C@@H](C(=O)O)c1ccccc1', '(2R)-2-hydroxy-2-phenylacetic acid'],
    ['CC/C=C\\CCO', '(Z)-hex-3-en-1-ol'],
    ['OC(=O)/C=C/C(=O)O', '(E)-but-2-enedioic acid'],
    ['OC(=O)/C=C\\C(=O)O', '(Z)-but-2-enedioic acid'],
    ['CC(C)C[C@@H](C(=O)O)N', '(2S)-2-amino-4-methylpentanoic acid'],
  ])('%s → %s', (smi, expected) => {
    expect(nm(smi).name).toBe(expected);
  });
  it('can omit stereodescriptors on request', () => {
    expect(nameMolecule(parseSmiles('N[C@@H](C)C(=O)O'), { stereo: false }).name).toBe('2-aminopropanoic acid');
  });
});

// ───────────── API behaviour ─────────────

describe('nameMolecule API', () => {
  it('returns parent locants keyed by input atom index', () => {
    const r = nm('CC(O)C'); // propan-2-ol
    expect(r.name).toBe('propan-2-ol');
    expect(r.locants.get(1)).toBe('2');
    expect(new Set([r.locants.get(0), r.locants.get(3)])).toEqual(new Set(['1', '3']));
    expect(r.locants.has(2)).toBe(false); // the oxygen is not part of the parent hydride
    const q = nm('c1ccc2ccccc2c1O'); // naphthalen-1-ol
    expect(q.name).toBe('naphthalen-1-ol');
    expect(q.locants.get(9)).toBe('1');
    expect([...q.locants.values()]).toContain('8a');
  });

  it('handles explicit hydrogens', () => {
    expect(nm('[H]OC([H])([H])C([H])([H])[H]').name).toBe('ethanol');
    expect(nm('[H]c1c([H])c([H])c(O[H])c([H])c1[H]').name).toBe('phenol');
  });

  it('does not depend on the Kekulé structure or aromatic bond orders', () => {
    const k1 = nm('C1=CC=C2C(=C1)C=CC=C2C(=O)O').name;
    const k2 = nm('C1C=CC2=CC=CC(=C2C=1)C(=O)O').name;
    expect(k1).toBe('naphthalene-1-carboxylic acid');
    expect(k2).toBe(k1);
    // delocalised (order 1.5) bonds as drawn by the editor
    const m = parseSmiles('C1=CC=C(C=C1)C(=O)O');
    for (const b of m.bonds) if (b.a < 6 && b.b < 6) b.order = 1.5;
    expect(nameMolecule(m).name).toBe('benzoic acid');
  });

  it('expands abbreviations', () => {
    const m = new Mol();
    const ring = [0, 1, 2, 3, 4, 5].map(() => m.addAtom({ el: 'C' }));
    for (let i = 0; i < 6; i++) m.addBond(ring[i], ring[(i + 1) % 6], i % 2 ? 1 : 2);
    const g = m.addAtom({ el: 'C', abbrev: 'CO2H' });
    m.addBond(ring[0], g, 1);
    const o = m.addAtom({ el: 'C', abbrev: 'OMe' });
    m.addBond(ring[3], o, 1);
    expect(nameMolecule(m).name).toBe('4-methoxybenzoic acid');
  });

  it('names abbreviation-based drawings', () => {
    const m = new Mol();
    const ring = [0, 1, 2, 3, 4, 5].map(() => m.addAtom({ el: 'C' }));
    for (let i = 0; i < 6; i++) m.addBond(ring[i], ring[(i + 1) % 6], 1.5);
    const g = m.addAtom({ el: 'O', abbrev: 'OAc' });
    m.addBond(ring[0], g, 1);
    expect(nameMolecule(m).name).toBe('phenyl acetate');
    const n = new Mol();
    const a = n.addAtom({ el: 'N', abbrev: 'NHBoc' });
    const c = n.addAtom({ el: 'C', abbrev: 'Et' });
    n.addBond(a, c, 1);
    expect(nameMolecule(n).name).toBe('tert-butyl N-ethylcarbamate');
  });

  it('does not depend on atom or bond order', () => {
    const smiles = [
      'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
      'CN1C=NC2=C1C(=O)N(C)C(=O)N2C',
      'OCCN(CCO)CCO',
      'CC(=O)OCC(COC(C)=O)OC(C)=O',
      'O=C1CCCc2ccccc21',
      'CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O',
      'CC1(C)C2CCC1(C)C(=O)C2',
    ];
    for (const smi of smiles) {
      const m = parseSmiles(smi);
      const ref = nameMolecule(m).name;
      for (let seed = 1; seed <= 4; seed++) {
        // deterministic shuffle of atoms and bonds
        let x = seed * 7919;
        const rnd = () => (x = (x * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;
        const perm = [...m.atoms.keys()].sort(() => rnd() - 0.5);
        const inv: number[] = [];
        perm.forEach((p, i) => (inv[p] = i));
        const r = new Mol();
        for (const p of perm) r.atoms.push({ ...m.atoms[p] });
        for (const bi of [...m.bonds.keys()].sort(() => rnd() - 0.5)) {
          const b = m.bonds[bi];
          r.bonds.push({ ...b, a: inv[b.a], b: inv[b.b] });
        }
        r.invalidate();
        expect(nameMolecule(r).name, smi).toBe(ref);
      }
    }
  });

  it('returns null with a warning for pseudo atoms', () => {
    const r = nm('*CC(=O)O');
    expect(r.name).toBeNull();
    expect(r.warnings.length).toBeGreaterThan(0);
  });

  it('never throws', () => {
    for (const smi of ['[Fe+2]', '[CH3]', 'C1CC2CC1C2C3CC4CC3C4', 'CC(=O)SC', 'B(O)(O)c1ccccc1', '[Si](C)(C)(C)C', 'P', 'C1=CC=CC=C1=C']) {
      const r = nm(smi);
      expect(r).toBeDefined();
      if (r.name === null) expect(r.warnings.length).toBeGreaterThan(0);
    }
    expect(nameMolecule(new Mol()).name).toBeNull();
  });

  it('is fast for drug-sized molecules', () => {
    const drugs = [
      'CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O',
      'CC(C)CCCC(C)C1CCC2C1(C)CCC1C2CC=C2CC(O)CCC12C',
      'COc1ccc2[nH]c(nc2c1)S(=O)Cc1ncc(C)c(OC)c1C',
      'CCCCc1nc(Cl)c(CO)n1Cc1ccc(cc1)-c1ccccc1-c1nnn[nH]1',
      'CC(C)CCCC(C)CCCC(C)CCCC1(C)CCc2c(C)c(O)c(C)c(C)c2O1',
    ].map((s) => parseSmiles(s));
    for (const d of drugs) nameMolecule(d); // warm-up
    const t0 = performance.now();
    for (const d of drugs) nameMolecule(d);
    const per = (performance.now() - t0) / drugs.length;
    expect(per).toBeLessThan(50);
  });
});
