// Template library: ready-made structures for the template palette.
// SMILES with stereochemistry were taken from PubChem (isomeric SMILES) where available; the
// hyaluronan repeat unit and trans-tranexamic acid were assembled programmatically and checked
// against PubChem names. Coordinates are generated on demand with layoutMol and cached.
import { Mol } from './mol';
import { parseSmiles } from './smiles';
import { layoutMol } from './layout2d';

export interface TemplateItem {
  name: string;
  smiles: string;
  tags?: string[];
}

export interface TemplateGroup {
  name: string;
  items: TemplateItem[];
}

const I = (name: string, smiles: string, ...tags: string[]): TemplateItem => (tags.length ? { name, smiles, tags } : { name, smiles });

export const TEMPLATE_GROUPS: TemplateGroup[] = [
  {
    name: 'Rings',
    items: [
      I('Cyclopropane', 'C1CC1'),
      I('Cyclobutane', 'C1CCC1'),
      I('Cyclopentane', 'C1CCCC1'),
      I('Cyclohexane', 'C1CCCCC1'),
      I('Cycloheptane', 'C1CCCCCC1'),
      I('Cyclooctane', 'C1CCCCCCC1'),
      I('Benzene', 'c1ccccc1', 'aromatic'),
      I('Cyclopentadiene', 'C1C=CC=C1'),
      I('Cyclohexene', 'C1CCC=CC1'),
      I('Naphthalene', 'c1ccc2ccccc2c1', 'aromatic'),
      I('Anthracene', 'c1ccc2cc3ccccc3cc2c1', 'aromatic'),
      I('Phenanthrene', 'c1ccc2c(c1)ccc1ccccc12', 'aromatic'),
      I('Indane', 'C1Cc2ccccc2C1'),
      I('Tetralin', 'C1CCc2ccccc2C1'),
      I('Decalin', 'C1CCC2CCCCC2C1'),
      I('Norbornane', 'C1CC2CCC1C2', 'bridged'),
      I('Adamantane', 'C1C2CC3CC1CC(C2)C3', 'bridged'),
    ],
  },
  {
    name: 'Heterocycles',
    items: [
      I('Pyridine', 'c1ccncc1', 'aromatic'),
      I('Pyrrole', 'c1cc[nH]c1', 'aromatic'),
      I('Furan', 'c1ccoc1', 'aromatic'),
      I('Thiophene', 'c1ccsc1', 'aromatic'),
      I('Imidazole', 'c1c[nH]cn1', 'aromatic'),
      I('Pyrazole', 'c1cc[nH]n1', 'aromatic'),
      I('Oxazole', 'c1cocn1', 'aromatic'),
      I('Thiazole', 'c1cscn1', 'aromatic'),
      I('1,2,3-Triazole', 'c1c[nH]nn1', 'aromatic'),
      I('Tetrazole', 'c1nn[nH]n1', 'aromatic'),
      I('Pyrimidine', 'c1cncnc1', 'aromatic'),
      I('Pyrazine', 'c1cnccn1', 'aromatic'),
      I('Pyridazine', 'c1ccnnc1', 'aromatic'),
      I('1,3,5-Triazine', 'c1ncncn1', 'aromatic'),
      I('Indole', 'c1ccc2[nH]ccc2c1', 'aromatic'),
      I('Benzimidazole', 'c1ccc2[nH]cnc2c1', 'aromatic'),
      I('Benzofuran', 'c1ccc2occc2c1', 'aromatic'),
      I('Quinoline', 'c1ccc2ncccc2c1', 'aromatic'),
      I('Isoquinoline', 'c1ccc2cnccc2c1', 'aromatic'),
      I('Purine', 'c1ncc2[nH]cnc2n1', 'aromatic'),
      I('Pteridine', 'c1cnc2ncncc2n1', 'aromatic'),
      I('Piperidine', 'C1CCNCC1'),
      I('Piperazine', 'C1CNCCN1'),
      I('Morpholine', 'C1COCCN1'),
      I('Pyrrolidine', 'C1CCNC1'),
      I('Tetrahydrofuran', 'C1CCOC1'),
      I('Tetrahydropyran', 'C1CCOCC1'),
      I('1,4-Dioxane', 'C1COCCO1'),
      I('Oxirane', 'C1CO1'),
      I('Aziridine', 'C1CN1'),
      I('Coumarin', 'O=c1ccc2ccccc2o1'),
      I('Chromane', 'C1Cc2ccccc2OC1'),
    ],
  },
  {
    name: 'Amino acids',
    items: [
      I('Glycine', 'C(C(=O)O)N', 'Gly', 'G'),
      I('L-Alanine', 'C[C@@H](C(=O)O)N', 'Ala', 'A'),
      I('L-Valine', 'CC(C)[C@@H](C(=O)O)N', 'Val', 'V'),
      I('L-Leucine', 'CC(C)C[C@@H](C(=O)O)N', 'Leu', 'L'),
      I('L-Isoleucine', 'CC[C@H](C)[C@@H](C(=O)O)N', 'Ile', 'I'),
      I('L-Proline', 'C1C[C@H](NC1)C(=O)O', 'Pro', 'P'),
      I('L-Phenylalanine', 'C1=CC=C(C=C1)C[C@@H](C(=O)O)N', 'Phe', 'F'),
      I('L-Tryptophan', 'C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N', 'Trp', 'W'),
      I('L-Methionine', 'CSCC[C@@H](C(=O)O)N', 'Met', 'M'),
      I('L-Serine', 'C([C@@H](C(=O)O)N)O', 'Ser', 'S'),
      I('L-Threonine', 'C[C@H]([C@@H](C(=O)O)N)O', 'Thr', 'T'),
      I('L-Cysteine', 'C([C@@H](C(=O)O)N)S', 'Cys', 'C'),
      I('L-Tyrosine', 'C1=CC(=CC=C1C[C@@H](C(=O)O)N)O', 'Tyr', 'Y'),
      I('L-Asparagine', 'C([C@@H](C(=O)O)N)C(=O)N', 'Asn', 'N'),
      I('L-Glutamine', 'C(CC(=O)N)[C@@H](C(=O)O)N', 'Gln', 'Q'),
      I('L-Aspartic acid', 'C([C@@H](C(=O)O)N)C(=O)O', 'Asp', 'D'),
      I('L-Glutamic acid', 'C(CC(=O)O)[C@@H](C(=O)O)N', 'Glu', 'E'),
      I('L-Lysine', 'C(CCN)C[C@@H](C(=O)O)N', 'Lys', 'K'),
      I('L-Arginine', 'C(C[C@@H](C(=O)O)N)CN=C(N)N', 'Arg', 'R'),
      I('L-Histidine', 'C1=C(NC=N1)C[C@@H](C(=O)O)N', 'His', 'H'),
    ],
  },
  {
    name: 'Sugars & nucleosides',
    items: [
      I('α-D-Glucopyranose', 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O', 'sugar'),
      I('β-D-Glucopyranose', 'C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O', 'sugar'),
      I('β-D-Fructofuranose', 'C([C@@H]1[C@H]([C@@H]([C@](O1)(CO)O)O)O)O', 'sugar'),
      I('β-D-Ribofuranose', 'C([C@@H]1[C@H]([C@H]([C@@H](O1)O)O)O)O', 'sugar'),
      I('2-Deoxy-β-D-ribofuranose', 'C1[C@@H]([C@H](O[C@H]1O)CO)O', 'sugar'),
      I('Sucrose', 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O', 'sugar', 'disaccharide'),
      I('Adenosine', 'C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O)N', 'nucleoside'),
      I('Uridine', 'C1=CN(C(=O)NC1=O)[C@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O', 'nucleoside'),
    ],
  },
  {
    name: 'Nucleobases',
    items: [
      I('Adenine', 'C1=NC2=NC=NC(=C2N1)N', 'purine'),
      I('Guanine', 'C1=NC2=C(N1)C(=O)NC(=N2)N', 'purine'),
      I('Cytosine', 'C1=C(NC(=O)N=C1)N', 'pyrimidine'),
      I('Thymine', 'CC1=CNC(=O)NC1=O', 'pyrimidine'),
      I('Uracil', 'C1=CNC(=O)NC1=O', 'pyrimidine'),
    ],
  },
  {
    name: 'Cosmetic actives & excipients',
    items: [
      // vitamins & skin actives
      I('Niacinamide', 'C1=CC(=CN=C1)C(=O)N', 'vitamin B3', 'active'),
      I('Retinol', 'CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/CO)/C)/C', 'vitamin A', 'retinoid'),
      I('Retinyl palmitate', 'CCCCCCCCCCCCCCCC(=O)OC/C=C(\\C)/C=C/C=C(\\C)/C=C/C1=C(CCCC1(C)C)C', 'vitamin A ester', 'retinoid'),
      I('Tretinoin', 'CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/C(=O)O)/C)/C', 'retinoic acid', 'retinoid'),
      I('L-Ascorbic acid', 'C([C@@H]([C@@H]1C(=C(C(=O)O1)O)O)O)O', 'vitamin C', 'antioxidant'),
      I('Ascorbyl glucoside', 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)OC2=C([C@H](OC2=O)[C@H](CO)O)O)O)O)O)O', 'vitamin C derivative'),
      I('3-O-Ethyl ascorbic acid', 'CCOC1=C(C(=O)O[C@@H]1[C@H](CO)O)O', 'vitamin C derivative'),
      I('α-Tocopherol', 'CC1=C(C2=C(CC[C@@](O2)(C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C(=C1O)C)C', 'vitamin E', 'antioxidant'),
      I('Tocopheryl acetate', 'CC1=C(C(=C(C2=C1O[C@](CC2)(C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C)OC(=O)C)C', 'vitamin E ester'),
      I('D-Panthenol', 'CC(C)(CO)[C@H](C(=O)NCCCO)O', 'provitamin B5', 'humectant'),
      I('Allantoin', 'C1(C(=O)NC(=O)N1)NC(=O)N', 'soothing'),
      I('(−)-α-Bisabolol', 'CC1=CC[C@H](CC1)[C@](C)(CCC=C(C)C)O', 'soothing'),
      I('Bakuchiol', 'CC(=CCC[C@@](C)(C=C)/C=C/C1=CC=C(C=C1)O)C', 'retinol alternative'),
      I('Caffeine', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'active'),
      I('Ectoine', 'CC1=NCC[C@H](N1)C(=O)O', 'active'),
      I('Ferulic acid', 'COC1=C(C=CC(=C1)/C=C/C(=O)O)O', 'antioxidant'),
      I('Resveratrol', 'C1=CC(=CC=C1/C=C/C2=CC(=CC(=C2)O)O)O', 'antioxidant'),
      // exfoliants & brightening agents
      I('Salicylic acid', 'C1=CC=C(C(=C1)C(=O)O)O', 'BHA', 'exfoliant'),
      I('Glycolic acid', 'C(C(=O)O)O', 'AHA', 'exfoliant'),
      I('L-Lactic acid', 'C[C@@H](C(=O)O)O', 'AHA', 'exfoliant'),
      I('Mandelic acid', 'C1=CC=C(C=C1)C(C(=O)O)O', 'AHA', 'exfoliant'),
      I('Azelaic acid', 'C(CCCC(=O)O)CCCC(=O)O', 'active'),
      I('Kojic acid', 'C1=C(OC=C(C1=O)O)CO', 'brightening'),
      I('α-Arbutin', 'C1=CC(=CC=C1O)O[C@@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)O', 'brightening'),
      I('Tranexamic acid', 'NC[C@H]1CC[C@@H](CC1)C(=O)O', 'brightening', 'trans'),
      // humectants & glycols
      I('Glycerol', 'C(C(CO)O)O', 'humectant'),
      I('Urea', 'C(=O)(N)N', 'humectant'),
      I('Hyaluronic acid (repeat unit)', 'CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O[C@@H]1O[C@H](C(=O)O)[C@@H](O)[C@H](O)[C@H]1O', 'humectant', 'GlcA-β(1→3)-GlcNAc'),
      I('Propylene glycol', 'CC(CO)O', 'humectant', 'solvent'),
      I('Butylene glycol', 'CC(CCO)O', 'humectant', 'solvent'),
      I('1,2-Hexanediol', 'CCCCC(CO)O', 'humectant', 'preservative booster'),
      // emollients & lipids
      I('Squalane', 'CC(C)CCCC(C)CCCC(C)CCCCC(C)CCCC(C)CCCC(C)C', 'emollient'),
      I('Ceramide NP', 'CCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@@H]([C@@H](CCCCCCCCCCCCCC)O)O', 'barrier lipid', 'N-stearoyl phytosphingosine'),
      I('Cetyl alcohol', 'CCCCCCCCCCCCCCCCO', 'fatty alcohol', 'thickener'),
      I('Stearic acid', 'CCCCCCCCCCCCCCCCCC(=O)O', 'fatty acid'),
      I('Glyceryl monostearate', 'CCCCCCCCCCCCCCCCCC(=O)OCC(CO)O', 'emulsifier'),
      I('Isopropyl myristate', 'CCCCCCCCCCCCCC(=O)OC(C)C', 'emollient'),
      // preservatives
      I('Phenoxyethanol', 'C1=CC=C(C=C1)OCCO', 'preservative'),
      I('Ethylhexylglycerin', 'CCCCC(CC)COCC(CO)O', 'preservative booster'),
      I('Methylparaben', 'COC(=O)C1=CC=C(C=C1)O', 'preservative'),
      I('Propylparaben', 'CCCOC(=O)C1=CC=C(C=C1)O', 'preservative'),
      I('Sorbic acid', 'C/C=C/C=C/C(=O)O', 'preservative'),
      I('Potassium sorbate', 'C/C=C/C=C/C(=O)[O-].[K+]', 'preservative'),
      I('Sodium benzoate', 'C1=CC=C(C=C1)C(=O)[O-].[Na+]', 'preservative'),
      I('Benzyl alcohol', 'C1=CC=C(C=C1)CO', 'preservative'),
      I('Dehydroacetic acid', 'CC1=CC(=O)C(C(=O)O1)C(=O)C', 'preservative'),
      // UV filters
      I('Avobenzone', 'CC(C)(C)C1=CC=C(C=C1)C(=O)CC(=O)C2=CC=C(C=C2)OC', 'UV filter', 'UVA'),
      I('Octocrylene', 'CCCCC(CC)COC(=O)C(=C(C1=CC=CC=C1)C2=CC=CC=C2)C#N', 'UV filter', 'UVB'),
      I('Oxybenzone', 'COC1=CC(=C(C=C1)C(=O)C2=CC=CC=C2)O', 'UV filter', 'benzophenone-3'),
      I('Homosalate', 'CC1CC(CC(C1)(C)C)OC(=O)C2=CC=CC=C2O', 'UV filter', 'UVB'),
      I('Octisalate', 'CCCCC(CC)COC(=O)C1=CC=CC=C1O', 'UV filter', 'ethylhexyl salicylate'),
      I('Ethylhexyl methoxycinnamate', 'CCCCC(CC)COC(=O)/C=C/C1=CC=C(C=C1)OC', 'UV filter', 'octinoxate'),
      // surfactants
      I('Sodium lauryl sulfate', 'CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]', 'surfactant', 'anionic'),
      I('Sodium laureth sulfate (n=2)', 'CCCCCCCCCCCCOCCOCCOS(=O)(=O)[O-].[Na+]', 'surfactant', 'anionic'),
      I('Cocamidopropyl betaine (C12)', 'CCCCCCCCCCCC(=O)NCCC[N+](C)(C)CC(=O)[O-]', 'surfactant', 'amphoteric'),
      // fragrance & sensates
      I('(−)-Menthol', 'C[C@@H]1CC[C@H]([C@@H](C1)O)C(C)C', 'cooling', 'fragrance'),
      I('(R)-Limonene', 'CC1=CC[C@@H](CC1)C(=C)C', 'fragrance'),
      I('Linalool', 'CC(=CCCC(C)(C=C)O)C', 'fragrance'),
      I('Vanillin', 'COC1=C(C=CC(=C1)C=O)O', 'fragrance'),
      I('Coumarin', 'C1=CC=C2C(=C1)C=CC(=O)O2', 'fragrance'),
      // anti-dandruff
      I('Zinc pyrithione', 'C1=CC=[N+](C(=C1)[S-])[O-].C1=CC=[N+](C(=C1)[S-])[O-].[Zn+2]', 'anti-dandruff'),
    ],
  },
  {
    name: 'Functional groups',
    items: [
      I('Carboxylic acid', 'CC(=O)O', 'acetic acid'),
      I('Ester', 'CCOC(C)=O', 'ethyl acetate'),
      I('Amide', 'CC(=O)NC', 'N-methylacetamide'),
      I('Ketone', 'CCC(C)=O', 'butanone'),
      I('Aldehyde', 'CCC=O', 'propanal'),
      I('Nitrile', 'CCC#N', 'propionitrile'),
      I('Nitro', 'CC[N+](=O)[O-]', 'nitroethane'),
      I('Sulfonamide', 'CS(=O)(=O)NC', 'N-methylmethanesulfonamide'),
      I('Acetal', 'CC(OC)OC', '1,1-dimethoxyethane'),
      I('Imine', 'CC=NC', 'N-methylethanimine'),
      I('Oxime', 'CC(C)=NO', 'acetone oxime'),
      I('Urea', 'CNC(=O)NC', 'N,N′-dimethylurea'),
      I('Carbamate', 'CNC(=O)OC', 'methyl N-methylcarbamate'),
      I('Anhydride', 'CC(=O)OC(C)=O', 'acetic anhydride'),
      I('Alcohol', 'CCO', 'ethanol'),
      I('Ether', 'CCOCC', 'diethyl ether'),
      I('Amine', 'CCN(CC)CC', 'triethylamine'),
      I('Thiol', 'CCS', 'ethanethiol'),
      I('Acyl chloride', 'CC(=O)Cl', 'acetyl chloride'),
      I('Sulfonic acid', 'CS(=O)(=O)O', 'methanesulfonic acid'),
      I('Lactone', 'O=C1CCCO1', 'γ-butyrolactone'),
      I('Lactam', 'O=C1CCCN1', '2-pyrrolidone'),
    ],
  },
  {
    name: 'Steroids & terpenes',
    items: [
      I('Steroid nucleus (gonane)', 'C1CCC2C(C1)CCC1C2CCC2CCCC12', 'steroid'),
      I('Cholesterol', 'C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C', 'steroid'),
      I('Testosterone', 'C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C', 'steroid'),
      I('Estradiol', 'C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O', 'steroid'),
      I('β-Carotene', 'CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C2=C(CCCC2(C)C)C)\\C)\\C)/C)/C', 'carotenoid', 'terpene'),
      I('(+)-Camphor', 'C[C@@]12CC[C@@H](C1(C)C)CC2=O', 'terpene', 'bridged'),
    ],
  },
];

const cache = new Map<string, Mol>();

/** Parses and lays out a template (cached); returns a fresh copy the caller may modify. */
export function templateMol(item: TemplateItem): Mol {
  let m = cache.get(item.smiles);
  if (!m) {
    m = parseSmiles(item.smiles);
    m.name = item.name;
    layoutMol(m);
    cache.set(item.smiles, m);
  }
  const out = m.clone();
  out.name = item.name;
  return out;
}

/** All template items (flattened), e.g. for search. */
export function allTemplates(): TemplateItem[] {
  return TEMPLATE_GROUPS.flatMap((g) => g.items);
}
