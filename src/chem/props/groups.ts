// Functional-group perception with SMARTS on an aromaticity-perceived, implicit-H match graph.
// Each definition lists one or more SMARTS; `core` selects the pattern positions reported as the
// group's atoms (context atoms such as the R in R–C(=O)–O–R' are matched but not reported).
// Matches with identical core atom sets are counted once.
import { buildMatchGraph, compileSmarts, findMatches, MatchGraph } from './smarts';
import { Mol } from '../mol';

export interface FunctionalGroupHit {
  name: string;
  count: number;
  /** One atom-index list per occurrence (indices of the Mol passed to detectGroupsInGraph). */
  atoms: number[][];
}

interface GroupDef {
  name: string;
  smarts: string[];
  /** Pattern positions to report (default: all). Applied to every SMARTS of the definition. */
  core?: number[];
}

// Reusable fragments
const NOT_ACYL = '!$([#6]=[O,S,N])';
const AMINE_N = '+0;!$(N~[!#6;!#1]);!$(N[#6]=[O,S,N]);!$(N[#6]#N)';
const ACETAL_C = '$([CX4]([OD2])[OD2])';
const HEMIACETAL_C = '$([CX4]([OD2])[OX2H1])';
const ETHER_C = `[#6;${NOT_ACYL};!${ACETAL_C};!${HEMIACETAL_C}]`;

const GROUPS: GroupDef[] = [
  { name: 'primary alcohol', smarts: ['[OX2H1][CX4;H2,H3;!$(C[OX2][#6])]'], core: [0, 1] },
  { name: 'secondary alcohol', smarts: ['[OX2H1][CX4;H1;!$(C[OX2][#6])]'], core: [0, 1] },
  { name: 'tertiary alcohol', smarts: ['[OX2H1][CX4;H0;!$(C[OX2][#6])]'], core: [0, 1] },
  { name: 'phenol', smarts: ['[OX2H1]c'], core: [0, 1] },
  { name: 'ether', smarts: [`[OD2;!r3](${ETHER_C})${ETHER_C}`], core: [0] },
  { name: 'enol ether', smarts: ['[OD2]([#6])[CX3;!$(C=[O,S,N])]=[CX3]'], core: [0, 2, 3] },
  { name: 'epoxide', smarts: ['[#6]1[OX2][#6]1'] },
  { name: 'acetal/ketal', smarts: ['[CX4;!$(C(O)(O)O)]([OX2][#6])[OX2][#6]'], core: [0, 1, 3] },
  { name: 'hemiacetal', smarts: ['[CX4]([OX2H1])[OX2][#6]'], core: [0, 1, 2] },
  { name: 'peroxide', smarts: ['[OX2][OX2]'] },
  { name: 'aldehyde', smarts: ['[CX3H1](=O)[#6]', '[CX3H2]=O'], core: [0, 1] },
  { name: 'ketone', smarts: ['[#6][CX3](=O)[#6]'], core: [1, 2] },
  { name: 'carboxylic acid', smarts: ['[CX3](=O)[OX2H1]'] },
  { name: 'carboxylate', smarts: ['[CX3](=O)[OX1-]', '[CX3](=O)[OX2][#3,#11,#19,#37,#55]'], core: [0, 1, 2] },
  { name: 'ester', smarts: ['[CX3;$([CH1]),$(C[#6])](=O)-!@[OX2H0][#6;!$(C=[O,S])]'], core: [0, 1, 2] },
  { name: 'lactone', smarts: ['[CX3;$([CH1]),$(C[#6])](=O)-@[OX2H0][#6;!$(C=[O,S])]'], core: [0, 1, 2] },
  { name: 'anhydride', smarts: ['[CX3](=O)[OX2][CX3]=O'] },
  { name: 'acyl halide', smarts: ['[CX3](=O)[F,Cl,Br,I]'] },
  { name: 'amide', smarts: ['[CX3;$([CH1]),$(C[#6])](=O)-!@[NX3]'] },
  { name: 'lactam', smarts: ['[CX3;$([CH1]),$(C[#6])](=O)-@[NX3]'] },
  { name: 'urea', smarts: ['[NX3][CX3](=[OX1])[NX3]'] },
  { name: 'carbamate', smarts: ['[NX3][CX3](=[OX1])[OX2]'] },
  { name: 'primary amine', smarts: [`[NX3;H2;${AMINE_N};!$(Na)][#6]`], core: [0] },
  { name: 'secondary amine', smarts: [`[NX3;H1;${AMINE_N};!$(Na)]([#6])[#6]`], core: [0] },
  { name: 'tertiary amine', smarts: [`[NX3;H0;${AMINE_N};!$(Na)]([#6])([#6])[#6]`], core: [0] },
  { name: 'aniline', smarts: [`[NX3;${AMINE_N}]c`], core: [0] },
  { name: 'ammonium (protonated amine)', smarts: ['[NX4+;H1,H2,H3;!$(N~[!#6;!#1])]'] },
  { name: 'quaternary ammonium', smarts: ['[NX4+;H0]([#6])([#6])([#6])[#6]'], core: [0] },
  { name: 'nitrile', smarts: ['[NX1]#[CX2]'] },
  { name: 'nitro', smarts: ['[NX3+](=O)[OX1-]', '[NX3](=O)=O'] },
  { name: 'imine', smarts: ['[CX3;!$(C-[#7,#8,#16])]=[NX2;!$(N-[#7,#8])]'] },
  { name: 'oxime', smarts: ['[CX3]=[NX2][OX2]'] },
  { name: 'hydrazone', smarts: ['[CX3]=[NX2][NX3]'] },
  { name: 'hydrazine', smarts: ['[NX3][NX3]'] },
  { name: 'azo', smarts: ['[#6][NX2]=[NX2][#6]'], core: [1, 2] },
  { name: 'alkyl halide', smarts: ['[F,Cl,Br,I][CX4]'], core: [0] },
  { name: 'aryl halide', smarts: ['[F,Cl,Br,I]c'], core: [0] },
  { name: 'thiol', smarts: ['[SX2H1][#6]'], core: [0] },
  { name: 'sulfide', smarts: ['[SX2H0]([#6])[#6]'], core: [0] },
  { name: 'disulfide', smarts: ['[SX2][SX2]'] },
  { name: 'sulfoxide', smarts: ['[SX3](=O)([#6])[#6]', '[SX3+]([OX1-])([#6])[#6]'], core: [0, 1] },
  { name: 'sulfone', smarts: ['[SX4](=O)(=O)([#6])[#6]'], core: [0, 1, 2] },
  { name: 'sulfonamide', smarts: ['[SX4](=O)(=O)[NX3]'] },
  { name: 'sulfonic acid', smarts: ['[SX4](=O)(=O)([#6])[OX2H1]'], core: [0, 1, 2, 4] },
  { name: 'sulfonate', smarts: ['[SX4](=O)(=O)([#6])[OX1-,$([OX2][#6]),$([OX2][#3,#11,#19])]'], core: [0, 1, 2, 4] },
  { name: 'sulfate ester', smarts: ['[SX4](=O)(=O)([OX2][#6])[OX2H1,OX1-,$([OX2][#6]),$([OX2][#3,#11,#19])]'], core: [0, 1, 2, 3, 5] },
  { name: 'phosphate', smarts: ['[PX4](=[OX1])([OX2,OX1-])([OX2,OX1-])[OX2,OX1-]'] },
  { name: 'alkene', smarts: ['[CX3]=[CX3]'] },
  { name: 'alkyne', smarts: ['[CX2]#[CX2]'] },
  {
    name: 'α,β-unsaturated carbonyl (Michael acceptor)',
    smarts: ['[CX3]=[CX3]-[CX3]=[OX1]', '[CX3]=[CX3]-[CX2]#[NX1]', '[CX3]=[CX3]-[SX4](=O)=O', '[CX3]=[CX3]-[NX3+](=O)[OX1-]'],
    core: [0, 1, 2, 3],
  },
  { name: 'skin sensitisation alert (aldehyde)', smarts: ['[CX3H1](=O)[#6]', '[CX3H2]=O'], core: [0, 1] },
];

/**
 * Detects functional groups on a prepared match graph (hMode 'none'). Atom indices are graph
 * source indices (Mol atom indices). Aromatic rings are reported as 'arene' / 'heteroarene'.
 */
export function detectGroupsInGraph(g: MatchGraph): FunctionalGroupHit[] {
  const out: FunctionalGroupHit[] = [];
  for (const def of GROUPS) {
    const seen = new Set<string>();
    const atoms: number[][] = [];
    for (const s of def.smarts) {
      const pat = compileSmarts(s);
      for (const m of findMatches(pat, g, true)) {
        const core = (def.core ? def.core.filter((k) => k < m.length).map((k) => m[k]) : m).map((k) => g.src[k]);
        const sorted = [...new Set(core)].sort((a, b) => a - b);
        const key = sorted.join(',');
        if (seen.has(key)) continue;
        seen.add(key);
        atoms.push(sorted);
      }
    }
    if (atoms.length) out.push({ name: def.name, count: atoms.length, atoms });
  }
  // Aromatic rings (SSSR rings whose atoms are all aromatic in the match graph)
  const arene: number[][] = [];
  const hetero: number[][] = [];
  for (const ring of g.rings.rings) {
    const gi = ring.map((a) => g.fromMol[a]);
    if (gi.some((k) => k < 0 || !g.arom[k])) continue;
    const list = [...ring].sort((a, b) => a - b);
    if (gi.every((k) => g.z[k] === 6)) arene.push(list);
    else hetero.push(list);
  }
  if (arene.length) out.push({ name: 'arene', count: arene.length, atoms: arene });
  if (hetero.length) out.push({ name: 'heteroarene', count: hetero.length, atoms: hetero });
  return out;
}

/** Convenience wrapper for a Mol whose abbreviations are already expanded. */
export function detectGroups(mol: Mol): FunctionalGroupHit[] {
  return detectGroupsInGraph(buildMatchGraph(mol, 'none'));
}
