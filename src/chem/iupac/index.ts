// Local IUPAC name generation.
//
// Public API: nameMolecule(mol, { stereo }) → { name, warnings, locants }. Never throws.
//
// Conventions (details in namer.ts / ringnames.ts):
//  • IUPAC 2013 substitutive nomenclature. Where PubChem (OpenEye LexiChem) uses a form that is also
//    acceptable IUPAC, that form is preferred so names can be cross-checked against PubChem:
//    "2-acetyloxybenzoic acid", "5-methyl-2-propan-2-ylcyclohexan-1-ol", "1-phenylethanone" (ketone
//    locant omitted on two-carbon parents), "prop-2-enyl", "benzyl", "tert-butyl", "pentoxy",
//    a hyphen after a closing bracket ("(4-hydroxyphenyl)-phenylmethanone"), "cyclohexene".
//  • Deliberate differences from PubChem: trivial names not retained by IUPAC 2013 are not used
//    (propan-2-ylbenzene not cumene, 1,4-dimethylbenzene not xylene, ethyne, trichloromethane);
//    enclosing marks nest as ( ) [ ] { } (PubChem only uses ( ) and [ ]); cations are "-aminium"
//    (PubChem: "azanium"); oximes are "N-hydroxy…imine"; anhydrides are "acetic anhydride";
//    numbering always follows the IUPAC rule order (indicated hydrogen, suffixes, hydro/ene, prefixes).
//  • Hydro prefixes and indicated hydrogen follow PubChem: ring positions that are saturated only
//    because they carry =O, a non-hydrogen substituent, or cannot form a ring double bond are implied
//    and not cited ("1,3,7-trimethylpurine-2,6-dione", "chromen-2-one", "1-methylindole"); added
//    hydrogen is written as indicated hydrogen in front ("1H-pyridin-2-one",
//    "3,4-dihydro-2H-naphthalen-1-one"). Hydro prefixes are cited next to the parent, after the
//    alphabetised substituent prefixes ("6-methoxy-1,2,3,4-tetrahydronaphthalene").
//  • Stereodescriptors: "(2S)-…", "(2E,4E)-…"; a single E/Z descriptor of one parent is cited without
//    locant and first ("(E)-but-2-ene", "(Z,12R)-…"). Descriptors of substituents are placed inside
//    the substituent's enclosing marks ("3-[(2S)-1-methylpyrrolidin-2-yl]pyridine").
//  • Salts: cations then anions, alphabetical, with multiplying prefixes ("sodium benzoate",
//    "disodium butanedioate", "N,N,N-trimethylhexadecan-1-aminium bromide"). Water / hydrogen halides
//    accompanying a single species become "hydrate" / "hydrochloride". Other multi-component inputs
//    are joined with "; " in alphabetical order ("ethanol; 2-hydroxybenzoic acid").
//  • Unsupported structures (name = null, reason in warnings): pseudo atoms, radicals, charged carbon,
//    metals/organometallics, polycyclic von Baeyer systems (> 2 rings) and fused ring systems without a
//    template, thioesters/thioacids, phosphines and other uncommon heteroatom functions.
import { Mol } from '../mol';
import { expandAbbreviations } from '../abbreviations';
import { assignCIP } from '../cip';
import { element } from '../elements';
import { buildGraph, NGraph } from './graph';
import { Namer, NamingError, CIPInfo } from './namer';
import { multiplier, multiplierComplex } from './numerals';
import { startsAmbiguous, enclose, alphaKey } from './assemble';

export interface IupacResult {
  /** null when the structure is outside what the namer supports */
  name: string | null;
  /** e.g. "stereodescriptors in substituents omitted" */
  warnings: string[];
  /** Locants of the parent hydride: input atom index → locant label ("1", "4a", "N"…) */
  locants: Map<number, string>;
}

interface ComponentName {
  name: string;
  charge: number;
  key: string;
  formula: string;
  locants: Map<number, string>; // graph index → label
  graphIdx: number[];
}

/** Names a structure. Never throws. */
export function nameMolecule(mol: Mol, opts: { stereo?: boolean } = {}): IupacResult {
  const warnings: string[] = [];
  try {
    return nameInner(mol, opts, warnings);
  } catch (e) {
    const msg = e instanceof NamingError ? e.message : 'internal error: ' + (e instanceof Error ? e.message : String(e));
    return { name: null, warnings: [...warnings, msg], locants: new Map() };
  }
}

function nameInner(mol: Mol, opts: { stereo?: boolean }, warnings: string[]): IupacResult {
  if (!mol.atoms.length) return { name: null, warnings: ['empty structure'], locants: new Map() };
  const em = expandAbbreviations(mol);
  if (em.atoms.some((a) => a.el === 'R' || a.el === '*' || a.alias)) {
    return { name: null, warnings: ['structure contains pseudo atoms (R/*) – no IUPAC name'], locants: new Map() };
  }
  const bg = buildGraph(em);
  warnings.push(...bg.warnings);
  if (!bg.graph) return { name: null, warnings: [...warnings, bg.error ?? 'unsupported structure'], locants: new Map() };
  const G = bg.graph;

  // CIP descriptors (graph indices)
  let cipAll: CIPInfo | null = null;
  const wantStereo = opts.stereo !== false && (em.tetra.length > 0 || em.dbStereo.length > 0);
  if (wantStereo) {
    try {
      const res = assignCIP(em);
      const inv = new Map<number, number>();
      G.src.forEach((s, i) => inv.set(s, i));
      cipAll = { centers: new Map(), bonds: new Map() };
      for (const [a, d] of res.centers) {
        const gi = inv.get(a);
        if (gi !== undefined) cipAll.centers.set(gi, d);
      }
      for (const [b, d] of res.bonds) {
        const bd = em.bonds[b];
        if (!bd) continue;
        const i = inv.get(bd.a), j = inv.get(bd.b);
        if (i === undefined || j === undefined) continue;
        cipAll.bonds.set(Math.min(i, j) + ',' + Math.max(i, j), d);
      }
      if (!cipAll.centers.size && !cipAll.bonds.size) cipAll = null;
    } catch {
      warnings.push('stereodescriptors could not be assigned');
      cipAll = null;
    }
  }

  const comps = G.mol.components();
  const names: ComponentName[] = [];
  let stereoMissing = false;
  for (const comp of comps) {
    const { graph: sg, idx } = subGraph(G, comp);
    const formula = formulaKey(sg);
    const charge = sg.charge.reduce((a, b) => a + b, 0);
    const inorg = inorganicName(sg, formula, charge);
    if (inorg) {
      names.push({ name: inorg, charge, key: formula, formula, locants: new Map(), graphIdx: idx });
      continue;
    }
    // map CIP to the component
    let cip: CIPInfo | null = null;
    if (cipAll) {
      cip = { centers: new Map(), bonds: new Map() };
      const local = new Map<number, number>();
      idx.forEach((gi, li) => local.set(gi, li));
      for (const [a, d] of cipAll.centers) if (local.has(a)) cip.centers.set(local.get(a)!, d);
      for (const [k, d] of cipAll.bonds) {
        const [i, j] = k.split(',').map(Number);
        if (local.has(i) && local.has(j)) {
          const li = local.get(i)!, lj = local.get(j)!;
          cip.bonds.set(Math.min(li, lj) + ',' + Math.max(li, lj), d);
        }
      }
    }
    const namer = new Namer(sg, cip);
    const r = namer.nameAll();
    warnings.push(...namer.warnings);
    if (cip) {
      const total = cip.centers.size + cip.bonds.size;
      if (countDescriptors(r.name) < total) stereoMissing = true;
    }
    names.push({ name: r.name, charge, key: formula + '|' + r.name, formula, locants: r.locants, graphIdx: idx });
  }
  if (stereoMissing) warnings.push('some stereodescriptors could not be placed in the name and were omitted');

  const name = combine(names);
  const locants = new Map<number, string>();
  for (const c of names) {
    for (const [li, lab] of c.locants) {
      const src = G.src[c.graphIdx[li]];
      if (src !== undefined && src < mol.atoms.length) locants.set(src, lab);
    }
  }
  return { name, warnings: [...new Set(warnings)], locants };
}

/** Number of CIP descriptors cited in a name ("(2S,3R)-…", "[(E)-…]"). */
function countDescriptors(name: string): number {
  let n = 0;
  for (const m of name.matchAll(/\(([^()]*)\)-/g)) {
    const parts = m[1].split(',');
    if (parts.every((p) => /^(\d+[a-z]?'*)?[RSrsEZ]$/.test(p))) n += parts.length;
  }
  return n;
}

/** Component subgraph with index map (local → whole-graph index). */
function subGraph(G: NGraph, comp: number[]): { graph: NGraph; idx: number[] } {
  if (comp.length === G.n) return { graph: G, idx: comp.map((_, i) => i) };
  const { mol, map } = G.mol.subset(comp);
  const idx: number[] = new Array(mol.atoms.length);
  map.forEach((nw, old) => {
    if (nw >= 0) idx[nw] = old;
  });
  const srcBond: number[] = [];
  const set = new Set(comp);
  G.mol.bonds.forEach((b, bi) => {
    if (set.has(b.a) && set.has(b.b)) srcBond.push(G.srcBond[bi]);
  });
  const src = idx.map((i) => G.src[i]);
  return { graph: new NGraph(mol, src, srcBond), idx };
}

function formulaKey(g: NGraph): string {
  const counts = new Map<string, number>();
  let h = 0;
  for (let i = 0; i < g.n; i++) {
    counts.set(g.el[i], (counts.get(g.el[i]) ?? 0) + 1);
    h += g.h[i];
  }
  if (h) counts.set('H', (counts.get('H') ?? 0) + h);
  const charge = g.charge.reduce((a, b) => a + b, 0);
  return [...counts.entries()].sort().map(([e, n]) => e + n).join('') + (charge ? (charge > 0 ? '+' : '') + charge : '');
}

const INORGANIC: Record<string, string> = {
  'H2O1': 'water', 'H3N1': 'ammonia', 'H4N1+1': 'ammonium', 'H1O1-1': 'hydroxide', 'Cl1H1': 'hydrogen chloride',
  'Br1H1': 'hydrogen bromide', 'H1I1': 'hydrogen iodide', 'F1H1': 'hydrogen fluoride', 'H2S1': 'hydrogen sulfide',
  'H2O2': 'hydrogen peroxide', 'H2O4S1': 'sulfuric acid', 'O4S1-2': 'sulfate', 'H1O4S1-1': 'hydrogen sulfate',
  'H1N1O3': 'nitric acid', 'N1O3-1': 'nitrate', 'N1O2-1': 'nitrite', 'H3O4P1': 'phosphoric acid', 'O4P1-3': 'phosphate',
  'H1O4P1-2': 'hydrogen phosphate', 'H2O4P1-1': 'dihydrogen phosphate', 'C1O3-2': 'carbonate', 'C1H1O3-1': 'hydrogen carbonate',
  'C1O2': 'carbon dioxide', 'C1Cl2O1': 'carbonyl dichloride', 'B1H3O3': 'boric acid', 'H3P1': 'phosphane',
  'Cl2O1S1': 'thionyl dichloride', 'Cl2O2S1': 'sulfuryl dichloride', 'Cl3O1P1': 'phosphoryl trichloride', 'Cl3P1': 'trichlorophosphane', 'C1O1': 'carbon monoxide', 'C1N1-1': 'cyanide', 'Cl1O4-1': 'perchlorate', 'B1F4-1': 'tetrafluoroborate',
  'F6P1-1': 'hexafluorophosphate', 'O3S1-2': 'sulfite', 'H1O3S1-1': 'hydrogen sulfite', 'Cl1O3-1': 'chlorate', 'Cl1O1-1': 'hypochlorite',
  'O1-2': 'oxide', 'S1-2': 'sulfide', 'H1-1': 'hydride', 'C1H1N1': 'formonitrile', 'H1+1': 'hydron', 'C1H4': 'methane',
  'C1S2': 'carbon disulfide', 'O3': 'ozone', 'O2': 'dioxygen', 'N2': 'dinitrogen', 'H2': 'dihydrogen', 'Cl2': 'chlorine', 'Br2': 'bromine', 'I2': 'iodine',
};

const SIMPLE_METALS = new Set(['Li', 'Na', 'K', 'Rb', 'Cs', 'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Al', 'Zn', 'Ag', 'Ga', 'Cd']);

function inorganicName(g: NGraph, formula: string, charge: number): string | null {
  if (INORGANIC[formula] && !(formula === 'C1H4')) return INORGANIC[formula];
  if (g.n === 1) {
    const el = g.el[0];
    const e = element(el);
    if (!e) return null;
    if (charge < 0 && g.h[0] === 0) {
      const an: Record<string, string> = { F: 'fluoride', Cl: 'chloride', Br: 'bromide', I: 'iodide' };
      if (an[el] && charge === -1) return an[el];
    }
    if (charge > 0 && g.h[0] === 0 && el !== 'C' && el !== 'N' && el !== 'O' && el !== 'S' && el !== 'P') {
      const nm = e.name.toLowerCase();
      return SIMPLE_METALS.has(el) ? nm : `${nm}(${charge}+)`;
    }
  }
  return null;
}

/** Combines component names (salts, hydrates, mixtures). */
function combine(names: ComponentName[]): string {
  if (names.length === 1) return names[0].name;
  // group identical components
  const groups: { c: ComponentName; n: number }[] = [];
  for (const c of names) {
    const g = groups.find((x) => x.c.key === c.key);
    if (g) g.n++;
    else groups.push({ c, n: 1 });
  }
  if (groups.length === 1 && groups[0].c.charge === 0) return groups[0].c.name;
  const mult = (name: string, n: number) => {
    if (n <= 1) return name;
    if (startsAmbiguous(name) || /\s/.test(name) || /[()[\]]/.test(name)) return multiplierComplex(n) + enclose(name);
    return multiplier(n) + name;
  };
  const sortAlpha = (a: { c: ComponentName }, b: { c: ComponentName }) => (alphaKey(a.c.name) < alphaKey(b.c.name) ? -1 : 1);
  const cations = groups.filter((x) => x.c.charge > 0).sort(sortAlpha);
  const anions = groups.filter((x) => x.c.charge < 0).sort(sortAlpha);
  const neutral = groups.filter((x) => x.c.charge === 0);
  const addends: string[] = [];
  const main: string[] = [];
  const isSalt = cations.length > 0 && anions.length > 0;
  const nonAddend = neutral.filter((x) => !['water', 'hydrogen chloride', 'hydrogen bromide', 'hydrogen iodide'].includes(x.c.name));
  const speciesCount = (isSalt ? 1 : cations.length + anions.length) + nonAddend.length;
  for (const x of neutral) {
    if (speciesCount === 1 && x.c.name === 'water') addends.push(multiplier(x.n) + 'hydrate');
    else if (speciesCount === 1 && x.c.name === 'hydrogen chloride') addends.push(multiplier(x.n) + 'hydrochloride');
    else if (speciesCount === 1 && x.c.name === 'hydrogen bromide') addends.push(multiplier(x.n) + 'hydrobromide');
    else if (speciesCount === 1 && x.c.name === 'hydrogen iodide') addends.push(multiplier(x.n) + 'hydroiodide');
    else main.push(x.n > 1 ? x.c.name + ' (' + x.n + ')' : x.c.name);
  }
  let saltName = '';
  if (isSalt) {
    saltName = [...cations.map((x) => mult(x.c.name, x.n)), ...anions.map((x) => mult(x.c.name, x.n))].join(' ');
  } else {
    for (const x of [...cations, ...anions]) main.unshift(x.n > 1 ? x.c.name + ' (' + x.n + ')' : x.c.name);
  }
  main.sort((a, b) => (alphaKey(a) < alphaKey(b) ? -1 : alphaKey(a) > alphaKey(b) ? 1 : a < b ? -1 : 1));
  const parts = [saltName, ...main].filter(Boolean);
  let out = parts.join('; ');
  if (addends.length) out += ' ' + addends.join(' ');
  return out;
}
