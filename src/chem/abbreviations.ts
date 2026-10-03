// Contracted labels ("nicknames"/superatoms) such as OMe, CO2H, Ph, Boc.
import { Mol } from './mol';
import { parseSmiles } from './smiles';
import { element } from './elements';

export interface AbbrevDef {
  /** SMILES of the group; the FIRST atom is the attachment atom. */
  smiles: string;
  /** Human readable description */
  name: string;
  category: 'alkyl' | 'aryl' | 'functional' | 'protecting' | 'other';
}

export const ABBREVIATIONS: Record<string, AbbrevDef> = {
  Me: { smiles: 'C', name: 'methyl', category: 'alkyl' },
  Et: { smiles: 'CC', name: 'ethyl', category: 'alkyl' },
  Pr: { smiles: 'CCC', name: 'propyl', category: 'alkyl' },
  nPr: { smiles: 'CCC', name: 'n-propyl', category: 'alkyl' },
  iPr: { smiles: 'C(C)C', name: 'isopropyl', category: 'alkyl' },
  Bu: { smiles: 'CCCC', name: 'butyl', category: 'alkyl' },
  nBu: { smiles: 'CCCC', name: 'n-butyl', category: 'alkyl' },
  iBu: { smiles: 'CC(C)C', name: 'isobutyl', category: 'alkyl' },
  sBu: { smiles: 'C(C)CC', name: 'sec-butyl', category: 'alkyl' },
  tBu: { smiles: 'C(C)(C)C', name: 'tert-butyl', category: 'alkyl' },
  Pent: { smiles: 'CCCCC', name: 'pentyl', category: 'alkyl' },
  Hex: { smiles: 'CCCCCC', name: 'hexyl', category: 'alkyl' },
  Cy: { smiles: 'C1CCCCC1', name: 'cyclohexyl', category: 'alkyl' },
  Ad: { smiles: 'C12CC3CC(CC(C3)C1)C2', name: 'adamantyl', category: 'alkyl' },
  C2H5: { smiles: 'CC', name: 'ethyl', category: 'alkyl' },
  C3H7: { smiles: 'CCC', name: 'propyl', category: 'alkyl' },
  C4H9: { smiles: 'CCCC', name: 'butyl', category: 'alkyl' },
  C5H11: { smiles: 'CCCCC', name: 'pentyl', category: 'alkyl' },
  C6H13: { smiles: 'CCCCCC', name: 'hexyl', category: 'alkyl' },
  C6H11: { smiles: 'C1CCCCC1', name: 'cyclohexyl', category: 'alkyl' },
  C6H5: { smiles: 'c1ccccc1', name: 'phenyl', category: 'aryl' },
  Ph: { smiles: 'c1ccccc1', name: 'phenyl', category: 'aryl' },
  Bn: { smiles: 'Cc1ccccc1', name: 'benzyl', category: 'aryl' },
  Tol: { smiles: 'c1ccc(C)cc1', name: 'p-tolyl', category: 'aryl' },
  Mes: { smiles: 'c1c(C)cc(C)cc1C', name: 'mesityl', category: 'aryl' },
  Np: { smiles: 'c1cccc2ccccc12', name: '1-naphthyl', category: 'aryl' },
  Py: { smiles: 'c1ccccn1', name: '2-pyridyl', category: 'aryl' },
  OH: { smiles: 'O', name: 'hydroxy', category: 'functional' },
  OMe: { smiles: 'OC', name: 'methoxy', category: 'functional' },
  OEt: { smiles: 'OCC', name: 'ethoxy', category: 'functional' },
  OiPr: { smiles: 'OC(C)C', name: 'isopropoxy', category: 'functional' },
  OtBu: { smiles: 'OC(C)(C)C', name: 'tert-butoxy', category: 'functional' },
  OPh: { smiles: 'Oc1ccccc1', name: 'phenoxy', category: 'functional' },
  OBn: { smiles: 'OCc1ccccc1', name: 'benzyloxy', category: 'functional' },
  OAc: { smiles: 'OC(=O)C', name: 'acetoxy', category: 'functional' },
  OCF3: { smiles: 'OC(F)(F)F', name: 'trifluoromethoxy', category: 'functional' },
  SMe: { smiles: 'SC', name: 'methylthio', category: 'functional' },
  SPh: { smiles: 'Sc1ccccc1', name: 'phenylthio', category: 'functional' },
  NMe2: { smiles: 'N(C)C', name: 'dimethylamino', category: 'functional' },
  NEt2: { smiles: 'N(CC)CC', name: 'diethylamino', category: 'functional' },
  NHMe: { smiles: 'NC', name: 'methylamino', category: 'functional' },
  NHAc: { smiles: 'NC(=O)C', name: 'acetamido', category: 'functional' },
  NHBoc: { smiles: 'NC(=O)OC(C)(C)C', name: 'Boc-amino', category: 'protecting' },
  NHCbz: { smiles: 'NC(=O)OCc1ccccc1', name: 'Cbz-amino', category: 'protecting' },
  NHPh: { smiles: 'Nc1ccccc1', name: 'anilino', category: 'functional' },
  NO2: { smiles: '[N+](=O)[O-]', name: 'nitro', category: 'functional' },
  NO: { smiles: 'N=O', name: 'nitroso', category: 'functional' },
  N3: { smiles: 'N=[N+]=[N-]', name: 'azido', category: 'functional' },
  CN: { smiles: 'C#N', name: 'cyano', category: 'functional' },
  NC: { smiles: '[N+]#[C-]', name: 'isocyano', category: 'functional' },
  NCO: { smiles: 'N=C=O', name: 'isocyanato', category: 'functional' },
  NCS: { smiles: 'N=C=S', name: 'isothiocyanato', category: 'functional' },
  SCN: { smiles: 'SC#N', name: 'thiocyanato', category: 'functional' },
  CF3: { smiles: 'C(F)(F)F', name: 'trifluoromethyl', category: 'functional' },
  CCl3: { smiles: 'C(Cl)(Cl)Cl', name: 'trichloromethyl', category: 'functional' },
  CHO: { smiles: 'C=O', name: 'formyl', category: 'functional' },
  COOH: { smiles: 'C(=O)O', name: 'carboxy', category: 'functional' },
  CO2H: { smiles: 'C(=O)O', name: 'carboxy', category: 'functional' },
  COO: { smiles: 'C(=O)[O-]', name: 'carboxylate', category: 'functional' },
  CO2: { smiles: 'C(=O)[O-]', name: 'carboxylate', category: 'functional' },
  CO2Me: { smiles: 'C(=O)OC', name: 'methoxycarbonyl', category: 'functional' },
  COOMe: { smiles: 'C(=O)OC', name: 'methoxycarbonyl', category: 'functional' },
  CO2Et: { smiles: 'C(=O)OCC', name: 'ethoxycarbonyl', category: 'functional' },
  COOEt: { smiles: 'C(=O)OCC', name: 'ethoxycarbonyl', category: 'functional' },
  CO2tBu: { smiles: 'C(=O)OC(C)(C)C', name: 'tert-butoxycarbonyl', category: 'functional' },
  COMe: { smiles: 'C(=O)C', name: 'acetyl', category: 'functional' },
  CONH2: { smiles: 'C(=O)N', name: 'carbamoyl', category: 'functional' },
  CONMe2: { smiles: 'C(=O)N(C)C', name: 'dimethylcarbamoyl', category: 'functional' },
  COCl: { smiles: 'C(=O)Cl', name: 'chlorocarbonyl', category: 'functional' },
  SO3H: { smiles: 'S(=O)(=O)O', name: 'sulfo', category: 'functional' },
  SO3: { smiles: 'S(=O)(=O)[O-]', name: 'sulfonate', category: 'functional' },
  SO2Cl: { smiles: 'S(=O)(=O)Cl', name: 'chlorosulfonyl', category: 'functional' },
  SO2NH2: { smiles: 'S(=O)(=O)N', name: 'sulfamoyl', category: 'functional' },
  SO2Me: { smiles: 'S(=O)(=O)C', name: 'methylsulfonyl', category: 'functional' },
  OSO3: { smiles: 'OS(=O)(=O)[O-]', name: 'sulfate', category: 'functional' },
  OPO3H2: { smiles: 'OP(=O)(O)O', name: 'phosphonooxy', category: 'functional' },
  PO3H2: { smiles: 'P(=O)(O)O', name: 'phosphono', category: 'functional' },
  PPh2: { smiles: 'P(c1ccccc1)c1ccccc1', name: 'diphenylphosphino', category: 'functional' },
  PPh3: { smiles: '[P+](c1ccccc1)(c1ccccc1)c1ccccc1', name: 'triphenylphosphonio', category: 'functional' },
  Bpin: { smiles: 'B1OC(C)(C)C(C)(C)O1', name: 'pinacolatoboryl', category: 'functional' },
  B: { smiles: 'B', name: 'boryl', category: 'functional' },
  Ac: { smiles: 'C(=O)C', name: 'acetyl', category: 'protecting' },
  Bz: { smiles: 'C(=O)c1ccccc1', name: 'benzoyl', category: 'protecting' },
  Piv: { smiles: 'C(=O)C(C)(C)C', name: 'pivaloyl', category: 'protecting' },
  Boc: { smiles: 'C(=O)OC(C)(C)C', name: 'tert-butoxycarbonyl', category: 'protecting' },
  Cbz: { smiles: 'C(=O)OCc1ccccc1', name: 'benzyloxycarbonyl', category: 'protecting' },
  Z: { smiles: 'C(=O)OCc1ccccc1', name: 'benzyloxycarbonyl', category: 'protecting' },
  Fmoc: { smiles: 'C(=O)OCC1c2ccccc2-c2ccccc21', name: '9-fluorenylmethoxycarbonyl', category: 'protecting' },
  Alloc: { smiles: 'C(=O)OCC=C', name: 'allyloxycarbonyl', category: 'protecting' },
  Troc: { smiles: 'C(=O)OCC(Cl)(Cl)Cl', name: '2,2,2-trichloroethoxycarbonyl', category: 'protecting' },
  Ts: { smiles: 'S(=O)(=O)c1ccc(C)cc1', name: 'tosyl', category: 'protecting' },
  Tos: { smiles: 'S(=O)(=O)c1ccc(C)cc1', name: 'tosyl', category: 'protecting' },
  Ms: { smiles: 'S(=O)(=O)C', name: 'mesyl', category: 'protecting' },
  Tf: { smiles: 'S(=O)(=O)C(F)(F)F', name: 'triflyl', category: 'protecting' },
  Ns: { smiles: 'S(=O)(=O)c1ccc([N+](=O)[O-])cc1', name: 'nosyl', category: 'protecting' },
  OTs: { smiles: 'OS(=O)(=O)c1ccc(C)cc1', name: 'tosylate', category: 'functional' },
  OMs: { smiles: 'OS(=O)(=O)C', name: 'mesylate', category: 'functional' },
  OTf: { smiles: 'OS(=O)(=O)C(F)(F)F', name: 'triflate', category: 'functional' },
  TMS: { smiles: '[Si](C)(C)C', name: 'trimethylsilyl', category: 'protecting' },
  TES: { smiles: '[Si](CC)(CC)CC', name: 'triethylsilyl', category: 'protecting' },
  TBS: { smiles: '[Si](C)(C)C(C)(C)C', name: 'tert-butyldimethylsilyl', category: 'protecting' },
  TBDMS: { smiles: '[Si](C)(C)C(C)(C)C', name: 'tert-butyldimethylsilyl', category: 'protecting' },
  TIPS: { smiles: '[Si](C(C)C)(C(C)C)C(C)C', name: 'triisopropylsilyl', category: 'protecting' },
  TBDPS: { smiles: '[Si](c1ccccc1)(c1ccccc1)C(C)(C)C', name: 'tert-butyldiphenylsilyl', category: 'protecting' },
  OTMS: { smiles: 'O[Si](C)(C)C', name: 'trimethylsilyloxy', category: 'protecting' },
  OTBS: { smiles: 'O[Si](C)(C)C(C)(C)C', name: 'TBS ether', category: 'protecting' },
  THP: { smiles: 'C1CCCCO1', name: 'tetrahydropyranyl', category: 'protecting' },
  MOM: { smiles: 'COC', name: 'methoxymethyl', category: 'protecting' },
  OMOM: { smiles: 'OCOC', name: 'methoxymethoxy', category: 'protecting' },
  SEM: { smiles: 'COCC[Si](C)(C)C', name: '2-(trimethylsilyl)ethoxymethyl', category: 'protecting' },
  PMB: { smiles: 'Cc1ccc(OC)cc1', name: 'p-methoxybenzyl', category: 'protecting' },
  Tr: { smiles: 'C(c1ccccc1)(c1ccccc1)c1ccccc1', name: 'trityl', category: 'protecting' },
  Trt: { smiles: 'C(c1ccccc1)(c1ccccc1)c1ccccc1', name: 'trityl', category: 'protecting' },
  DMT: { smiles: 'C(c1ccccc1)(c1ccc(OC)cc1)c1ccc(OC)cc1', name: 'dimethoxytrityl', category: 'protecting' },
};

/** Case-insensitive lookup helper ("ome" → "OMe"). Exact matches win. */
export function findAbbreviation(label: string): string | null {
  if (ABBREVIATIONS[label]) return label;
  const lc = label.toLowerCase();
  for (const k of Object.keys(ABBREVIATIONS)) if (k.toLowerCase() === lc) return k;
  return null;
}

/**
 * Reverses a label for display when its bond comes from the right: OMe → MeO, CO2H → HO2C,
 * NHBoc → BocHN, CH2OH → HOH2C.
 */
export function reverseLabel(label: string): string {
  const tokens = label.match(/[a-z]*[A-Z][a-z]*\d*|\d+|[+-]/g);
  if (!tokens || tokens.join('') !== label) return label;
  // keep trailing charge at end
  const charge: string[] = [];
  while (tokens.length && /^[+-]$/.test(tokens[tokens.length - 1])) charge.unshift(tokens.pop()!);
  return tokens.reverse().join('') + charge.join('');
}

const cache = new Map<string, Mol>();
function groupMol(smiles: string): Mol {
  let m = cache.get(smiles);
  if (!m) {
    m = parseSmiles(smiles);
    cache.set(smiles, m);
  }
  return m.clone();
}

/** Parses a label such as "CH2CH2OH" or "OCH2CH3" into a linear group (first heavy atom attaches). */
export function parseCondensedLabel(label: string): Mol | null {
  const tokens = label.match(/[A-Z][a-z]?\d*/g);
  if (!tokens || tokens.join('') !== label) return null;
  const m = new Mol();
  let last = -1;
  for (const t of tokens) {
    const mm = /^([A-Z][a-z]?)(\d*)$/.exec(t)!;
    const sym = mm[1];
    const n = mm[2] ? +mm[2] : 1;
    if (!element(sym)) return null;
    if (sym === 'H') {
      if (last < 0) return null;
      m.atoms[last].hCount = (m.atoms[last].hCount ?? 0) + n;
      continue;
    }
    for (let k = 0; k < n; k++) {
      const idx = m.addAtom({ el: sym, hCount: 0 });
      if (last >= 0) m.addBond(last, idx, 1);
      last = idx;
    }
  }
  if (m.atoms.length < 2) return null;
  return m;
}

/** Returns the structure an abbreviated atom stands for (attachment = atom 0), or null. */
export function abbreviationMol(label: string): Mol | null {
  const key = findAbbreviation(label);
  if (key) return groupMol(ABBREVIATIONS[key].smiles);
  return parseCondensedLabel(label);
}

/**
 * Replaces every abbreviation atom by its full structure. The attachment atom keeps the
 * abbreviation atom's bonds. Coordinates of new atoms are laid out roughly along the outgoing
 * direction (sufficient for formula/SMILES; the editor relayouts when expanding visually).
 * Atoms of the result carry `id` of the originating abbreviation atom for new atoms.
 */
export function expandAbbreviations(input: Mol): Mol {
  if (!input.atoms.some((a) => a.abbrev)) return input;
  const mol = input.clone();
  for (let i = 0; i < mol.atoms.length; i++) {
    const a = mol.atoms[i];
    if (!a.abbrev) continue;
    const g = abbreviationMol(a.abbrev);
    if (!g) {
      // unknown label → pseudo atom
      a.el = '*';
      a.alias = a.abbrev;
      delete a.abbrev;
      continue;
    }
    // direction away from first neighbour
    const nb = mol.neighbors(i);
    let dx = 1, dy = 0;
    if (nb.length) {
      dx = a.x - mol.atoms[nb[0]].x;
      dy = a.y - mol.atoms[nb[0]].y;
      const l = Math.hypot(dx, dy) || 1;
      dx /= l;
      dy /= l;
    }
    const off = mol.atoms.length;
    // place group atoms along BFS depth with zig-zag
    const depth = new Array(g.atoms.length).fill(-1);
    depth[0] = 0;
    const q = [0];
    for (let k = 0; k < q.length; k++) for (const w of g.neighbors(q[k])) if (depth[w] < 0) { depth[w] = depth[q[k]] + 1; q.push(w); }
    let sib = 0;
    g.atoms.forEach((ga, gi) => {
      if (gi === 0) return;
      const d = depth[gi];
      const zig = (d % 2 ? 0.5 : -0.5) + (sib++ % 3) * 0.15;
      mol.atoms.push({
        ...ga,
        id: a.id,
        x: a.x + dx * d * 0.87 - dy * zig,
        y: a.y + dy * d * 0.87 + dx * zig,
      });
    });
    const att = g.atoms[0];
    const newIdx = (gi: number) => (gi === 0 ? i : off + gi - 1);
    for (const b of g.bonds) mol.bonds.push({ ...b, id: -1, a: newIdx(b.a), b: newIdx(b.b) });
    a.el = att.el;
    a.charge = att.charge;
    a.isotope = att.isotope;
    if (att.hCount !== undefined) a.hCount = att.hCount;
    else delete a.hCount;
    delete a.abbrev;
    mol.invalidate();
  }
  mol.invalidate();
  return mol;
}
