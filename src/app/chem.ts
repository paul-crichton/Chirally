// App-level chemistry facade: glue between the document model and the chemistry modules.
import { Mol } from '../chem/mol';
import { parseSmiles, writeSmiles, suppressHydrogens, parseReactionSmiles } from '../chem/smiles';
import { layoutMol } from '../chem/layout2d';
import { clean2D } from '../chem/clean2d';
import { perceiveStereo2D, assignWedgesFromSpecs, explicitFusionHydrogens } from '../chem/stereo2d';
import { expandAbbreviations, abbreviationMol, attachmentCharge } from '../chem/abbreviations';
import { computeFormula, FormulaInfo } from '../chem/formula';
import { computeProperties, MolProperties } from '../chem/properties';
import { nameMolecule } from '../chem/iupac';
import { assignCIP } from '../chem/cip';
import { ChemDoc } from '../doc/types';
import { docToMol, adjacency, fragments, insertMol, createDoc, fragmentOf } from '../doc/document';
import * as F from '../chem/formats';
import { buildScene } from '../render/scene';

/** Parses SMILES (or reaction SMILES) and computes 2D coordinates. */
export function smilesToMol(s: string): Mol {
  const m = parseSmiles(s);
  layoutMol(m);
  explicitFusionHydrogens(m);
  return m;
}

/** Heuristic: does this text look like a SMILES string? */
export function looksLikeSmiles(s: string): boolean {
  const t = s.trim();
  if (!t || /\s/.test(t) || t.length > 2000) return false;
  if (!/^[A-Za-z0-9@+\-\[\]\(\)=#$:/\\%.*]+$/.test(t)) return false;
  // pure words such as "aspirin" are names, not SMILES
  if (/^[a-z]{4,}$/i.test(t) && !/^(C|O|N|S|P|F|Cl|Br|I|c|n|o|s)+$/.test(t)) return false;
  try {
    parseSmiles(t);
    return true;
  } catch {
    return false;
  }
}

/** Mol ready for chemistry: abbreviations expanded, stereo perceived from wedges, explicit H suppressed. */
export function chemMol(doc: ChemDoc, atomIds?: Iterable<number>): Mol {
  const { mol } = docToMol(doc, atomIds);
  const ex = expandAbbreviations(mol);
  perceiveStereo2D(ex);
  return suppressHydrogens(ex);
}

export function canonicalSmiles(mol: Mol): string {
  return writeSmiles(mol);
}

export interface FragmentAnalysis {
  atomIds: number[];
  smiles: string;
  formula: FormulaInfo;
  name: string | null;
  nameWarnings: string[];
  props: MolProperties | null;
  hasPseudo: boolean;
  /** document atom id for each atom index of the analysed molecule (used to map highlights back) */
  idMap: number[];
}

export interface Analysis {
  fragments: FragmentAnalysis[];
  /** Combined formula of everything analysed. */
  total: FormulaInfo | null;
  smiles: string;
  /** atom id → CIP label; bond id → E/Z */
  stereo: Map<number, string>;
  bondStereo: Map<number, string>;
  /** atom id → IUPAC locant of the parent structure */
  locants: Map<number, string>;
}

const nameCache = new Map<string, { name: string | null; warnings: string[]; locants: Map<number, string> }>();

/** Full analysis of the given atoms (grouped by connected fragment). */
export function analyze(doc: ChemDoc, atomIds: number[]): Analysis {
  const ids = new Set(atomIds);
  const adj = adjacency(doc);
  const frags = fragments(doc, adj).filter((f) => f.some((id) => ids.has(id)));
  const out: FragmentAnalysis[] = [];
  const stereo = new Map<number, string>();
  const bondStereo = new Map<number, string>();
  const locants = new Map<number, string>();
  for (const f of frags) {
    const { mol } = docToMol(doc, f);
    const ex = expandAbbreviations(mol);
    perceiveStereo2D(ex);
    // CIP labels on the drawn atoms (indices of expanded mol keep the original atoms first)
    try {
      const cip = assignCIP(ex);
      for (const [i, lab] of cip.centers) {
        const a = ex.atoms[i];
        if (i < mol.atoms.length) stereo.set(a.id, lab);
      }
      for (const [bi, lab] of cip.bonds) {
        const b = ex.bonds[bi];
        if (b && b.id > 0) bondStereo.set(b.id, lab);
      }
    } catch {
      /* CIP is optional */
    }
    const sup = suppressHydrogens(ex);
    let smiles = '';
    try {
      smiles = writeSmiles(sup);
    } catch {
      smiles = '';
    }
    const formula = computeFormula(mol);
    let name: string | null = null;
    let nameWarnings: string[] = [];
    if (!formula.hasPseudo) {
      const cached = nameCache.get(smiles);
      if (cached) {
        name = cached.name;
        nameWarnings = cached.warnings;
        for (const [i, l] of cached.locants) if (sup.atoms[i]) locants.set(sup.atoms[i].id, l);
      } else {
        try {
          const r = nameMolecule(sup, { stereo: true });
          name = r.name;
          nameWarnings = r.warnings;
          nameCache.set(smiles, { name, warnings: nameWarnings, locants: r.locants });
          for (const [i, l] of r.locants) if (sup.atoms[i]) locants.set(sup.atoms[i].id, l);
        } catch (e) {
          nameWarnings = [String(e)];
        }
      }
    }
    let props: MolProperties | null = null;
    if (!formula.hasPseudo && sup.atoms.length) {
      try {
        props = computeProperties(sup);
      } catch {
        props = null;
      }
    }
    out.push({ atomIds: f, smiles, formula, name, nameWarnings, props, hasPseudo: formula.hasPseudo, idMap: sup.atoms.map((a) => a.id) });
  }
  let total: FormulaInfo | null = null;
  if (out.length) {
    const { mol } = docToMol(doc, frags.flat());
    total = computeFormula(mol);
  }
  return { fragments: out, total, smiles: out.map((f) => f.smiles).filter(Boolean).join('.'), stereo, bondStereo, locants };
}

/** Cleans (regularises) the given atoms of the document in place. */
export function cleanAtoms(doc: ChemDoc, atomIds: number[]): void {
  const adj = adjacency(doc);
  const ids = new Set(atomIds);
  // clean whole fragments touched by the selection, but only move the selected atoms
  const frags = fragments(doc, adj).filter((f) => f.some((id) => ids.has(id)));
  for (const f of frags) {
    const { mol, index } = docToMol(doc, f);
    const moving = new Set<number>();
    for (const id of f) if (ids.has(id)) moving.add(index.get(id)!);
    const all = moving.size === f.length;
    clean2D(mol, all ? {} : { atoms: moving });
    mol.atoms.forEach((a, i) => {
      const da = doc.atoms.get(a.id)!;
      da.x = a.x;
      da.y = a.y;
      void i;
    });
    // wedges may have been re-assigned to keep stereo
    for (const b of mol.bonds) {
      const db = doc.bonds.get(b.id);
      if (!db) continue;
      db.style = b.style;
      const A = mol.atoms[b.a].id, B = mol.atoms[b.b].id;
      db.a = A;
      db.b = B;
    }
  }
}

/** Expands an abbreviation atom into real atoms with a fresh local layout. Returns false if not possible. */
export function expandLabel(doc: ChemDoc, atomId: number): boolean {
  const a = doc.atoms.get(atomId);
  if (!a?.abbrev) return false;
  const g = abbreviationMol(a.abbrev);
  if (!g) return false;
  const adj = adjacency(doc);
  // the whole fragment takes part (fixed) so the expanded group avoids the rest of the drawing
  const frag = new Set<number>(fragmentOf(doc, atomId, adj));
  // Build a local mol: neighbours (fixed) + group atoms (to lay out)
  const { mol, index } = docToMol(doc, frag);
  const ai = index.get(atomId)!;
  const off = mol.atoms.length;
  for (let k = 1; k < g.atoms.length; k++) mol.atoms.push({ ...g.atoms[k], id: -k, x: 0, y: 0 });
  for (const b of g.bonds) mol.bonds.push({ ...b, id: -1, a: b.a === 0 ? ai : off + b.a - 1, b: b.b === 0 ? ai : off + b.b - 1 });
  const att = g.atoms[0];
  // the label's own charge sits on the attachment atom (a free onium label is the neutral parent)
  const charge = attachmentCharge(g, a.charge, [...doc.bonds.values()].some((b) => b.a === atomId || b.b === atomId), a.abbrev);
  mol.atoms[ai] = { ...mol.atoms[ai], el: att.el, charge, hCount: att.hCount, abbrev: undefined };
  mol.invalidate();
  const fixed = new Set<number>();
  for (let i = 0; i < off; i++) fixed.add(i);
  layoutMol(mol, { fixed });
  // write back
  a.el = att.el;
  a.charge = charge;
  if (att.hCount !== undefined) a.hCount = att.hCount;
  else delete a.hCount;
  delete a.abbrev;
  const newIds = new Map<number, number>();
  newIds.set(ai, atomId);
  for (let i = off; i < mol.atoms.length; i++) {
    const m = mol.atoms[i];
    const id = doc.nextId++;
    doc.atoms.set(id, { ...m, id });
    newIds.set(i, id);
  }
  for (const b of mol.bonds) {
    if (b.id !== -1) continue;
    const id = doc.nextId++;
    doc.bonds.set(id, { ...b, id, a: newIds.get(b.a)!, b: newIds.get(b.b)! });
  }
  return true;
}

// ───────────── import ─────────────

export interface ImportResult {
  doc?: ChemDoc;
  mols?: Mol[];
  kind: string;
}

/** Converts arbitrary chemical text (or file contents) to molecules or a document. */
export function importText(text: string, fileName: string | null = null, opts: { suppressH?: boolean } = {}): ImportResult {
  const r = importTextRaw(text, fileName);
  if (opts.suppressH && r.mols) r.mols = r.mols.map(stripHydrogens);
  return r;
}

/** Removes ordinary explicit hydrogens (as found in PubChem SDF files) while keeping stereo wedges. */
export function stripHydrogens(m: Mol): Mol {
  if (!m.atoms.some((a) => a.el === 'H')) return m;
  if (!m.tetra.length && !m.dbStereo.length) perceiveStereo2D(m);
  const s = suppressHydrogens(m);
  if (s.tetra.length) {
    assignWedgesFromSpecs(s);
    explicitFusionHydrogens(s);
  }
  s.name = m.name;
  s.props = m.props;
  return s;
}

function importTextRaw(text: string, fileName: string | null): ImportResult {
  const fmt = F.detectFormat(fileName, text);
  switch (fmt) {
    case 'chirally':
      return { kind: 'chirally' };
    case 'cdxml':
      return { doc: F.readCDXML(text), kind: 'CDXML' };
    case 'mol': {
      const m = F.readMolfile(text);
      return { mols: [ensureLayout(m)], kind: 'MOL' };
    }
    case 'sdf':
      return { mols: F.readSDF(text).map(ensureLayout), kind: 'SDF' };
    case 'rxn': {
      const r = F.readRxn(text);
      return { doc: reactionDoc(r.reactants.map(ensureLayout), r.products.map(ensureLayout), r.agents.map(ensureLayout)), kind: 'RXN' };
    }
    case 'cml':
      return { mols: F.readCML(text).map(ensureLayout), kind: 'CML' };
    case 'xyz': {
      const m = F.readXYZ(text);
      return { mols: [to2D(m)], kind: 'XYZ' };
    }
    case 'smiles':
    default: {
      const lines = text.split(/\r?\n/).map((l) => l.trim()).filter(Boolean);
      if (lines.length === 1 && lines[0].split(/\s+/)[0].includes('>')) {
        const r = parseReactionSmiles(lines[0].split(/\s+/)[0]);
        const lay = (ms: Mol[]) => ms.map((m) => (layoutMol(m), m));
        return { doc: reactionDoc(lay(r.reactants), lay(r.products), lay(r.agents)), kind: 'reaction SMILES' };
      }
      const mols: Mol[] = [];
      for (const l of lines) {
        const [smi, ...rest] = l.split(/\s+/);
        const m = smilesToMol(smi);
        if (rest.length) m.name = rest.join(' ');
        mols.push(m);
      }
      if (!mols.length) throw new Error('Nothing to import');
      return { mols, kind: 'SMILES' };
    }
  }
}

/** Lays out molecules that came without coordinates; 3D molecules are flattened. */
export function ensureLayout(m: Mol): Mol {
  const has2D = m.atoms.length < 2 || m.atoms.some((a) => Math.abs(a.x) > 1e-4 || Math.abs(a.y) > 1e-4);
  if (m.props['dim'] === '3D') return to2D(m);
  if (!has2D) {
    perceiveStereo2D(m);
    layoutMol(m);
  }
  return m;
}

/** Converts a 3D structure to a clean 2D depiction preserving stereochemistry. */
export function to2D(m3: Mol): Mol {
  const m = suppressHydrogens(m3);
  // stereo from 3D coordinates is assigned by the formats module when available; fall back to none
  for (const a of m.atoms) delete a.z;
  layoutMol(m);
  return m;
}

/** Builds a reaction scheme document: reactants + … → products. */
export function reactionDoc(reactants: Mol[], products: Mol[], agents: Mol[] = []): ChemDoc {
  const doc = createDoc();
  let x = 0;
  const place = (m: Mol) => {
    // use rendered bounds so atom labels (OH, NH2…) are accounted for
    const tmp = createDoc();
    insertMol(tmp, m);
    const sb = buildScene(tmp, { ink: '#000', showErrors: false }).bounds;
    const bb = m.bbox();
    const left = sb ? sb.x1 : bb.minX, right = sb ? sb.x2 : bb.maxX;
    insertMol(doc, m, x - left, -(bb.minY + bb.maxY) / 2);
    x += right - left;
  };
  const plus = () => {
    x += 0.7;
    const id = doc.nextId++;
    doc.texts.set(id, { id, type: 'text', x, y: doc.style.fontSize * 1.4 * 0.36, text: '+', size: 1.4, align: 'center' });
    x += 0.7;
  };
  reactants.forEach((m, i) => {
    if (i) plus();
    place(m);
  });
  x += 0.6;
  const aid = doc.nextId++;
  const above = agents.map((m) => {
    try {
      return computeFormula(m).formula;
    } catch {
      return '';
    }
  }).filter(Boolean).join(', ');
  doc.arrows.set(aid, { id: aid, type: 'arrow', kind: 'reaction', x1: x, y1: 0, x2: x + 3, y2: 0, above: above || undefined });
  x += 3.6;
  products.forEach((m, i) => {
    if (i) plus();
    place(m);
  });
  return doc;
}

export { writeSmiles, F as formats };
