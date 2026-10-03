// Molecular property calculator: Lipinski/Veber descriptors (RDKit-compatible definitions),
// Wildman–Crippen logP/MR, Ertl TPSA, ESOL solubility, Potts–Guy skin permeability, HLB for
// amphiphiles and functional groups — with formulation (cosmetic/skincare) guidance.
import { Mol } from './mol';
import { element } from './elements';
import { implicitH } from './valence';
import { perceiveRings } from './rings';
import { perceiveAromaticity } from './aromaticity';
import { computeFormula } from './formula';
import { expandAbbreviations } from './abbreviations';
import { symmetryClasses } from './canon';
import { isPotentialStereocenter } from './stereo2d';
import { buildMatchGraph, compileSmarts, findMatches, MatchGraph } from './props/smarts';
import { crippen } from './props/crippen';
import { tpsaContribs } from './props/tpsa';
import { detectGroupsInGraph } from './props/groups';
import { computeHLB } from './props/hlb';

export interface MolProperties {
  formula: string; mw: number; exactMass: number; heavyAtoms: number; charge: number;
  /** Wildman–Crippen logP (RDKit MolLogP-compatible). */
  logP: number;
  /** Crippen molar refractivity. */
  mr: number;
  /** Ertl TPSA (N, O contributions; RDKit default). */
  tpsa: number;
  /** RDKit NumHDonors / NumHAcceptors. */
  hbd: number; hba: number;
  /** Lipinski counts: H on N/O (NH + OH bonds) and number of N + O atoms. */
  lipinskiHBD: number; lipinskiHBA: number;
  rotatableBonds: number; rings: number; aromaticRings: number; fractionCsp3: number; stereocenters: number;
  lipinski: { violations: number; pass: boolean };
  veber: { pass: boolean };
  /** ESOL (Delaney 2004): logS = 0.16 − 0.63·cLogP − 0.0062·MW + 0.066·RB − 0.74·AP (mol/L). */
  logS: number;
  /** Skin permeability (Potts & Guy 1992): log Kp (cm/h) = −2.72 + 0.71·logP − 0.0061·MW. */
  logKp: number;
  /** Only for surfactant-like/amphiphilic molecules. */
  hlb: { griffin: number | null; davies: number | null; note: string } | null;
  /** Formulation guidance: MW < 500 Da rule, logP 1–3 sweet spot. */
  skinPenetration: { under500Da: boolean; logPInRange: boolean; note: string };
  functionalGroups: { name: string; count: number; atoms: number[][] }[];
}

// RDKit definitions (rdMolDescriptors CalcNumHBD / CalcNumHBA / CalcNumRotatableBonds "Strict").
const HBD_SMARTS = '[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]';
const HBA_SMARTS =
  '[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0X2,o,s;+0])]';
const ROT_END =
  '!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])&!$([CH3])';
const ROT_SMARTS =
  `[${ROT_END}&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]=[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-,:;!@[${ROT_END}]`;

/** Expands abbreviations and returns a map from expanded atom index to input atom index. */
function expandWithMap(mol: Mol): { exp: Mol; toInput: number[] } {
  const exp = expandAbbreviations(mol);
  const toInput = exp.atoms.map((_, i) => i);
  if (exp !== mol) {
    const byId = new Map<number, number>();
    mol.atoms.forEach((a, i) => {
      if (a.abbrev) byId.set(a.id, i);
    });
    for (let i = mol.atoms.length; i < exp.atoms.length; i++) toInput[i] = byId.get(exp.atoms[i].id) ?? -1;
  }
  return { exp, toInput };
}

/** Per-atom Wildman–Crippen contributions (implicit-H contributions folded into their heavy atom). */
export function crippenContribs(mol: Mol): { logP: number[]; mr: number[] } {
  const { exp, toInput } = expandWithMap(mol);
  const c = crippen(exp);
  const logP = new Array(mol.atoms.length).fill(0);
  const mr = new Array(mol.atoms.length).fill(0);
  for (let i = 0; i < exp.atoms.length; i++) {
    const t = toInput[i];
    if (t < 0) continue;
    logP[t] += c.atomLogP[i];
    mr[t] += c.atomMR[i];
  }
  return { logP, mr };
}

function mapGroups(groups: MolProperties['functionalGroups'], toInput: number[]): MolProperties['functionalGroups'] {
  return groups.map((g) => {
    const atoms = g.atoms.map((list) => [...new Set(list.map((k) => toInput[k]).filter((k) => k >= 0))].sort((a, b) => a - b));
    return { name: g.name, count: g.count, atoms };
  });
}

/** Functional groups (atom indices refer to `mol`; atoms inside abbreviations map to the label atom). */
export function detectFunctionalGroups(mol: Mol): MolProperties['functionalGroups'] {
  const { exp, toInput } = expandWithMap(mol);
  return mapGroups(detectGroupsInGraph(buildMatchGraph(exp, 'none')), toInput);
}

function countMatches(g: MatchGraph, smarts: string): number {
  return findMatches(compileSmarts(smarts), g, true).length;
}

/**
 * Atoms that can be tetrahedral stereocentres given the constitution: potential centre elements
 * with four topologically distinct ligands (implicit H / lone pair counted as distinct ligands).
 * Aromatic bonds are normalised so the drawn Kekulé structure cannot split symmetry classes.
 */
function countStereocenters(mol: Mol, aromaticBonds: boolean[]): number {
  const m = mol.clone();
  m.atoms.forEach((a, i) => {
    if (element(a.el) && !a.abbrev) a.hCount = implicitH(mol, i);
  });
  m.bonds.forEach((b, i) => {
    if (aromaticBonds[i]) b.order = 1.5;
  });
  m.invalidate();
  const sym = symmetryClasses(m);
  let count = 0;
  for (let i = 0; i < mol.atoms.length; i++) {
    if (!isPotentialStereocenter(mol, i)) continue;
    const classes = mol.neighbors(i).map((j) => sym[j]);
    const h = implicitH(mol, i);
    for (let k = 0; k < h; k++) classes.push(-1);
    if (classes.length === 3) classes.push(-2); // lone pair (S, P, Se)
    if (classes.length === 4 && new Set(classes).size === 4) count++;
  }
  return count;
}

function skinNote(mw: number, logP: number, logKp: number): string {
  const parts: string[] = [];
  parts.push(mw < 500 ? `MW ${mw.toFixed(0)} Da < 500 Da (500-Dalton rule met)` : `MW ${mw.toFixed(0)} Da ≥ 500 Da: passive penetration of intact skin unlikely`);
  if (logP < 1) parts.push(`logP ${logP.toFixed(2)} < 1: hydrophilic, poor partitioning into the stratum corneum lipids`);
  else if (logP <= 3) parts.push(`logP ${logP.toFixed(2)} within the 1–3 optimum for dermal penetration`);
  else parts.push(`logP ${logP.toFixed(2)} > 3: lipophilic, tends to stay in the stratum corneum (reservoir) with slow transfer to the viable epidermis`);
  parts.push(`Potts–Guy log Kp ${logKp.toFixed(2)} cm/h`);
  return parts.join('; ');
}

/** Computes all properties. Abbreviations are expanded; pseudo atoms (R, *) contribute nothing. */
export function computeProperties(mol: Mol): MolProperties {
  const { exp, toInput } = expandWithMap(mol);
  const f = computeFormula(exp);
  const rings = perceiveRings(exp);
  const aro = perceiveAromaticity(exp, rings);
  const g = buildMatchGraph(exp, 'none', rings);
  const cr = crippen(exp);

  let heavyAtoms = 0, aromaticAtoms = 0, carbons = 0, csp3 = 0;
  let lipinskiHBD = 0, lipinskiHBA = 0;
  for (let i = 0; i < exp.atoms.length; i++) {
    const a = exp.atoms[i];
    const z = element(a.el)?.z ?? 0;
    if (a.abbrev || z <= 1) continue;
    heavyAtoms++;
    if (aro.atoms[i]) aromaticAtoms++;
    if (z === 7 || z === 8) {
      lipinskiHBA++;
      const k = g.fromMol[i];
      if (k >= 0) lipinskiHBD += g.hTotal[k];
    }
    if (z === 6) {
      carbons++;
      // RDKit FractionCSP3: carbon with four connections (hydrogens included)
      const k = g.fromMol[i];
      if (k >= 0 && g.degree[k] + g.hImplicit[k] === 4) csp3++;
    }
  }
  const tpsa = tpsaContribs(g).reduce((s, v) => s + v, 0);
  const hbd = countMatches(g, HBD_SMARTS);
  const hba = countMatches(g, HBA_SMARTS);
  const rotatableBonds = countMatches(g, ROT_SMARTS);
  const aromaticRings = rings.ringBonds.filter((rb) => rb.every((b) => aro.bonds[b])).length;
  const logP = cr.logP;
  const mw = f.mw;
  let violations = 0;
  if (mw > 500) violations++;
  if (logP > 5) violations++;
  if (lipinskiHBD > 5) violations++;
  if (lipinskiHBA > 10) violations++;
  const ap = heavyAtoms ? aromaticAtoms / heavyAtoms : 0;
  const logS = 0.16 - 0.63 * logP - 0.0062 * mw + 0.066 * rotatableBonds - 0.74 * ap;
  const logKp = -2.72 + 0.71 * logP - 0.0061 * mw;
  let hlb: MolProperties['hlb'] = null;
  try {
    hlb = computeHLB(g, mw);
  } catch {
    hlb = null;
  }
  return {
    formula: f.formula,
    mw,
    exactMass: f.exactMass,
    heavyAtoms,
    charge: f.charge,
    logP,
    mr: cr.mr,
    tpsa,
    hbd,
    hba,
    lipinskiHBD,
    lipinskiHBA,
    rotatableBonds,
    rings: rings.rings.length,
    aromaticRings,
    fractionCsp3: carbons ? csp3 / carbons : 0,
    stereocenters: countStereocenters(exp, aro.bonds),
    lipinski: { violations, pass: violations <= 1 },
    veber: { pass: rotatableBonds <= 10 && tpsa <= 140 },
    logS,
    logKp,
    hlb,
    skinPenetration: { under500Da: mw < 500, logPInRange: logP >= 1 && logP <= 3, note: skinNote(mw, logP, logKp) },
    functionalGroups: mapGroups(detectGroupsInGraph(g), toInput),
  };
}
