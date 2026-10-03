// MDL molfile (V2000 + V3000) and SD file reader/writer.
//
// Coordinates: molfiles store Å with y UP. On import y is flipped and 2D structures are scaled so
// the median bond length is 1 model unit (3D structures keep Å and their z; props.dim = '3D').
// Wedge semantics match ours: the narrow end of a wedge is the bond's FIRST atom (= bond.a).
import { Atom, Mol, BondStyle } from '../mol';
import { BY_SYMBOL, IMPLICIT_H_ELEMENTS, element } from '../elements';
import { bondOrderSum } from '../valence';
import { expandAbbreviations } from '../abbreviations';
import { assignWedgesFromSpecs } from '../stereo2d';
import {
  FormatError, guard, splitLines, normalizeImportedCoords, kekulizeFlagged, dropDefaultHCounts, applyAlias,
  contractSuperatoms, perceiveStereo, hasCoordinates, defaultImplicitH, SupGroup, normalizeElement,
  DEFAULT_BOND_LENGTH_ANGSTROM, pad, fixed,
} from './common';

export { hasCoordinates };

export interface MolfileWriteOptions {
  /** Force the V3000 format (automatic above 999 atoms or bonds). */
  v3000?: boolean;
  /** Header line 1 (default: mol.name). */
  title?: string;
  /** Length (Å) of one model unit for 2D structures (default 1.5). */
  bondLengthAngstrom?: number;
  /**
   * Expand abbreviation atoms into their full structure (default true). Expanded groups are
   * written as superatom (SUP) S-groups so other programs can show them contracted. With false,
   * abbreviations are written as alias (pseudo) atoms.
   */
  expandAbbreviations?: boolean;
}

// ───────────────────────────── reading ─────────────────────────────

/** Intermediate result of CTAB parsing (raw Å coordinates, y up). */
export interface RawCtab {
  mol: Mol;
  aromatic: boolean[];
  /** Valence field per atom (V2000 vvv / V3000 VAL=); 0 = unspecified, 15/-1 = zero. */
  valence: number[];
  /** Pending aliases (atom index → text), applied after kekulisation. */
  aliases: Map<number, string>;
  sups: SupGroup[];
  /** Keys of bonded atom pairs (duplicate detection without rebuilding adjacency). */
  bondKeys: Set<number>;
}

const isV30 = (line: string) => line.startsWith('M  V30');

/** Reads the first molecule of a molfile or SD file (V2000 or V3000). Data fields go to mol.props. */
export function readMolfile(text: string): Mol {
  return guard('Molfile', () => {
    const lines = splitLines(text);
    sniffWrongFormat(lines, 'Molfile');
    let end = lines.findIndex((l) => l.startsWith('$$$$'));
    if (end < 0) end = lines.length;
    return parseRecord(lines.slice(0, end), 'Molfile');
  });
}

/** Reads all records of an SD file. Records that fail to parse are skipped (unless all fail). */
export function readSDF(text: string): Mol[] {
  return guard('SDF', () => {
    const lines = splitLines(text);
    sniffWrongFormat(lines, 'SDF');
    const mols: Mol[] = [];
    const errors: string[] = [];
    let start = 0;
    let rec = 0;
    while (start < lines.length) {
      let end = start;
      while (end < lines.length && !lines[end].startsWith('$$$$')) end++;
      const recLines = lines.slice(start, end);
      if (recLines.some((l) => l.trim())) {
        rec++;
        try {
          mols.push(parseRecord(recLines, 'SDF'));
        } catch (e) {
          errors.push(`record ${rec}: ${(e as Error).message.replace(/^(SDF|Molfile): /, '')}`);
        }
      }
      start = end + 1;
    }
    if (!mols.length) throw new FormatError('SDF', errors.length ? errors[0] : 'no records found');
    return mols;
  });
}

function sniffWrongFormat(lines: string[], fmt: string): void {
  const head = lines.slice(0, 5).join('\n');
  if (lines[0]?.startsWith('$RXN')) throw new FormatError(fmt, 'this is a reaction (RXN) file — use readRxn');
  if (/<\?xml|<CDXML|<cml|<molecule/i.test(head)) throw new FormatError(fmt, 'this looks like an XML file (CDXML/CML), not a molfile');
  if (head.trimStart().startsWith('{')) throw new FormatError(fmt, 'this looks like a JSON file, not a molfile');
}

/** One SD record: CTAB + data items. */
function parseRecord(lines: string[], fmt: string): Mol {
  const { raw, end } = parseMolBlock(lines, 0, fmt);
  normalizeImportedCoords([raw.mol]);
  const mol = finishCtab(raw);
  readDataItems(lines, end, mol);
  return mol;
}

function readDataItems(lines: string[], from: number, mol: Mol): void {
  for (let p = from; p < lines.length; p++) {
    const line = lines[p];
    if (!line.startsWith('>')) continue;
    const m = /<([^>]*)>/.exec(line);
    const name = (m ? m[1] : line.slice(1)).trim();
    const vals: string[] = [];
    while (p + 1 < lines.length && lines[p + 1].trim() !== '' && !/^>\s.*<.*>/.test(lines[p + 1])) vals.push(lines[++p]);
    if (name) mol.props[name] = vals.join('\n');
  }
}

/** Index of the counts line of a molfile starting at `from` (standard: from + 3), or -1. */
function findCountsLine(lines: string[], from: number): number {
  const std = lines[from + 3];
  if (std !== undefined && isCountsLine(std)) return from + 3;
  // tolerate missing/extra header lines when the version tag is present
  for (let k = from; k < Math.min(lines.length, from + 8); k++) if (/V[23]000/.test(lines[k]) && isCountsLine(lines[k])) return k;
  return -1;
}

function isCountsLine(line: string): boolean {
  if (/^\s*\d+\s+\d+.*V3000/.test(line)) return true;
  const a = line.slice(0, 3), b = line.slice(3, 6);
  if (/^ *\d+$/.test(a) && /^ *\d+$/.test(b)) return true;
  return /^\s*\d+\s+\d+(\s|$)/.test(line);
}

const int = (s: string | undefined): number => {
  if (s === undefined) return NaN;
  const t = s.trim();
  return /^[+-]?\d+$/.test(t) ? parseInt(t, 10) : NaN;
};
const num = (s: string | undefined): number => {
  if (s === undefined) return NaN;
  const t = s.trim();
  return t === '' ? NaN : Number(t);
};

/** Parses header + CTAB (V2000 or V3000) beginning at `from`. Returns the line index after M  END. */
export function parseMolBlock(lines: string[], from: number, fmt = 'Molfile'): { raw: RawCtab; end: number } {
  const ci = findCountsLine(lines, from);
  if (ci < 0) {
    if (lines.slice(from).join('').trim() === '') throw new FormatError(fmt, 'empty input');
    throw new FormatError(fmt, 'no valid counts line found (expected a 3-line header followed by "aaabbb… V2000/V3000")');
  }
  const title = ci - 3 >= from ? lines[ci - 3].trim() : '';
  const counts = lines[ci];
  const res = /V3000/.test(counts) ? parseV3000(lines, ci, fmt) : parseV2000(lines, ci, fmt);
  res.raw.mol.name = title;
  return res;
}

// ───── symbols & codes ─────

/** Atom fields for a molfile atom symbol (elements, D/T, R groups, query/pseudo atoms). */
function atomFromSymbol(sym: string): Partial<Atom> & { el: string } {
  if (sym === 'D') return { el: 'H', isotope: 2 };
  if (sym === 'T') return { el: 'H', isotope: 3 };
  if (BY_SYMBOL.has(sym)) return { el: sym };
  if (sym === 'R#' || sym === 'R') return { el: 'R' };
  if (/^R\d+$/.test(sym)) return { el: 'R', alias: sym };
  const norm = normalizeElement(sym);
  if (norm && sym.length <= 2 && sym !== 'LP') return { el: norm };
  return sym === '*' ? { el: '*' } : { el: '*', alias: sym };
}

/** Applies a molfile bond type (shared by V2000/V3000). */
function bondFromType(type: number): { order: number; style: BondStyle; aromatic: boolean } {
  switch (type) {
    case 1: return { order: 1, style: 'plain', aromatic: false };
    case 2: return { order: 2, style: 'plain', aromatic: false };
    case 3: return { order: 3, style: 'plain', aromatic: false };
    case 4: return { order: 1.5, style: 'plain', aromatic: true };
    case 6: return { order: 1.5, style: 'plain', aromatic: true }; // single or aromatic (query)
    case 7: return { order: 1.5, style: 'plain', aromatic: true }; // double or aromatic (query)
    case 9: return { order: 1, style: 'dative', aromatic: false };
    case 10: return { order: 0, style: 'hbond', aromatic: false };
    case 0: return { order: 0, style: 'plain', aromatic: false };
    default: return { order: 1, style: 'plain', aromatic: false }; // 5 (S/D), 8 (any)
  }
}

function newRaw(): RawCtab {
  return { mol: new Mol(), aromatic: [], valence: [], aliases: new Map(), sups: [], bondKeys: new Set() };
}

function addRawBond(raw: RawCtab, a: number, b: number, type: number, fmt: string, lineNo: number): number {
  const n = raw.mol.atoms.length;
  if (!(a >= 0 && a < n && b >= 0 && b < n)) throw new FormatError(fmt, `bond on line ${lineNo + 1} references a missing atom`);
  const key = a < b ? a * 2097152 + b : b * 2097152 + a;
  if (a === b || raw.bondKeys.has(key)) return -1;
  raw.bondKeys.add(key);
  const t = bondFromType(type);
  const bi = raw.mol.addBond(a, b, t.order, t.style);
  raw.aromatic[bi] = t.aromatic;
  return bi;
}

// ───── V2000 ─────

function parseV2000(lines: string[], ci: number, fmt: string): { raw: RawCtab; end: number } {
  const counts = lines[ci];
  let na = int(counts.slice(0, 3));
  let nb = int(counts.slice(3, 6));
  if (!Number.isFinite(na) || !Number.isFinite(nb)) {
    const t = counts.trim().split(/\s+/);
    na = int(t[0]);
    nb = int(t[1]);
  }
  if (!Number.isFinite(na) || !Number.isFinite(nb) || na < 0 || nb < 0) throw new FormatError(fmt, `invalid counts line "${counts}"`);
  if (ci + na + nb > lines.length - 1) throw new FormatError(fmt, `file truncated: counts line announces ${na} atoms and ${nb} bonds`);
  const raw = newRaw();
  const mol = raw.mol;
  const massDiffIso = new Set<number>();
  let p = ci + 1;

  for (let k = 0; k < na; k++, p++) {
    const line = lines[p];
    if (line === undefined) throw new FormatError(fmt, 'atom block truncated');
    let x = num(line.slice(0, 10)), y = num(line.slice(10, 20)), z = num(line.slice(20, 30));
    let sym = line.slice(31, 34).trim();
    let md = int(line.slice(34, 36)), cc = int(line.slice(36, 39));
    let vv = int(line.slice(48, 51)), mm = int(line.slice(60, 63));
    if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z) || !/^[A-Za-z*#][A-Za-z0-9#]*$/.test(sym)) {
      // not fixed-width: fall back to whitespace-separated fields
      const t = line.trim().split(/\s+/);
      x = num(t[0]); y = num(t[1]); z = num(t[2]); sym = t[3] ?? '';
      md = int(t[4]); cc = int(t[5]); vv = int(t[9]); mm = int(t[13]);
      if (!Number.isFinite(x) || !Number.isFinite(y) || !sym) throw new FormatError(fmt, `cannot parse atom line ${p + 1}: "${line}"`);
      if (!Number.isFinite(z)) z = 0;
    }
    const fields = atomFromSymbol(sym);
    const i = mol.addAtom({ ...fields, x, y, z });
    const a = mol.atoms[i];
    if (md && Number.isFinite(md) && element(a.el) && a.isotope === undefined) {
      a.isotope = Math.round(element(a.el)!.mass) + md;
      massDiffIso.add(i);
    }
    if (cc >= 1 && cc <= 7) {
      if (cc === 4) a.radical = 1;
      else a.charge = 4 - cc;
    }
    raw.valence[i] = Number.isFinite(vv) ? vv : 0;
    if (mm > 0) a.map = mm;
  }

  for (let k = 0; k < nb; k++, p++) {
    const line = lines[p];
    if (line === undefined) throw new FormatError(fmt, 'bond block truncated');
    let a = int(line.slice(0, 3)), b = int(line.slice(3, 6)), type = int(line.slice(6, 9)), st = int(line.slice(9, 12));
    if (!Number.isFinite(a) || !Number.isFinite(b) || !Number.isFinite(type)) {
      const t = line.trim().split(/\s+/);
      a = int(t[0]); b = int(t[1]); type = int(t[2]); st = int(t[3]);
      if (!Number.isFinite(a) || !Number.isFinite(b) || !Number.isFinite(type)) throw new FormatError(fmt, `cannot parse bond line ${p + 1}: "${line}"`);
    }
    const bi = addRawBond(raw, a - 1, b - 1, type, fmt, p);
    if (bi < 0) continue;
    const bd = mol.bonds[bi];
    if (bd.order === 1 && bd.style === 'plain') {
      if (st === 1) bd.style = 'wedge';
      else if (st === 6) bd.style = 'hash';
      else if (st === 4) bd.style = 'wavy';
    } else if (bd.order === 2 && st === 3) bd.style = 'crossed';
  }

  // properties block
  const chg: [number, number][] = [];
  const rad: [number, number][] = [];
  const iso: [number, number][] = [];
  const sg = new Map<number, { type: string; atoms: number[]; label: string }>();
  const sgroup = (k: number) => {
    let g = sg.get(k);
    if (!g) sg.set(k, (g = { type: '', atoms: [], label: '' }));
    return g;
  };
  const pairs = (line: string): [number, number][] => {
    const t = line.slice(6).trim().split(/\s+/).map((s) => int(s));
    const out: [number, number][] = [];
    for (let k = 1; k + 1 < t.length && out.length < t[0]; k += 2) out.push([t[k], t[k + 1]]);
    return out;
  };
  for (; p < lines.length; p++) {
    const line = lines[p];
    if (line.startsWith('M  END')) {
      p++;
      break;
    }
    if (line.startsWith('$$$$') || line.startsWith('$MOL') || line.startsWith('>')) break;
    const tag = line.slice(0, 6);
    if (tag === 'M  CHG') chg.push(...pairs(line));
    else if (tag === 'M  RAD') rad.push(...pairs(line));
    else if (tag === 'M  ISO') iso.push(...pairs(line));
    else if (tag === 'M  RGP') {
      for (const [ai, r] of pairs(line)) {
        const a = mol.atoms[ai - 1];
        if (a && a.el === 'R') a.alias = 'R' + r;
      }
    } else if (tag === 'M  STY') {
      const t = line.slice(6).trim().split(/\s+/);
      for (let k = 1; k + 1 < t.length; k += 2) sgroup(int(t[k])).type = t[k + 1];
    } else if (tag === 'M  SAL') {
      const t = line.slice(6).trim().split(/\s+/).map((s) => int(s));
      const g = sgroup(t[0]);
      for (let k = 2; k < t.length && k < 2 + t[1]; k++) g.atoms.push(t[k] - 1);
    } else if (tag === 'M  SMT') {
      const m = /^M {2}SMT\s+(\d+)\s?(.*)$/.exec(line);
      if (m) sgroup(+m[1]).label = m[2].trim();
    } else if (line.startsWith('A  ')) {
      const ai = int(line.slice(3, 6)) || int(line.slice(3).trim().split(/\s+/)[0]);
      const text = lines[p + 1];
      p++;
      if (ai >= 1 && ai <= mol.atoms.length && text !== undefined && text.trim()) raw.aliases.set(ai - 1, text.trim());
    } else if (line.startsWith('G  ')) {
      p++; // obsolete group abbreviation: skip its text line
    } else if (line.startsWith('S  SKP')) {
      p += Math.max(0, int(line.slice(6)) || 0);
    }
  }

  // M  CHG / M  RAD supersede the atom block's charge and radical values; M  ISO its mass differences
  if (chg.length || rad.length) for (const a of mol.atoms) {
    a.charge = 0;
    delete a.radical;
  }
  for (const [ai, v] of chg) if (mol.atoms[ai - 1]) mol.atoms[ai - 1].charge = v;
  for (const [ai, v] of rad) if (mol.atoms[ai - 1]) setRadical(mol.atoms[ai - 1], v);
  if (iso.length) for (const i of massDiffIso) delete mol.atoms[i].isotope;
  for (const [ai, v] of iso) if (mol.atoms[ai - 1] && v > 0) mol.atoms[ai - 1].isotope = v;
  for (const g of sg.values()) if (g.type === 'SUP' && g.label) raw.sups.push({ atoms: g.atoms, label: g.label });
  return { raw, end: p };
}

/** Molfile radical code (1 singlet, 2 doublet, 3 triplet) → unpaired electrons. */
function setRadical(a: Atom, code: number): void {
  if (code === 2) a.radical = 1;
  else if (code === 1 || code === 3) a.radical = 2;
  else delete a.radical;
}

// ───── V3000 ─────

/** Collects logical "M  V30" lines (continuations joined) from `from` until M  END. */
export function collectV30(lines: string[], from: number): { v30: string[]; aliases: [number, string][]; end: number } {
  const v30: string[] = [];
  const aliases: [number, string][] = [];
  let buf = '';
  let p = from;
  for (; p < lines.length; p++) {
    const line = lines[p];
    if (line.startsWith('M  END')) {
      p++;
      break;
    }
    if (line.startsWith('$$$$') || line.startsWith('$MOL') || line.startsWith('>')) break;
    if (line.startsWith('A  ')) {
      const ai = int(line.slice(3).trim().split(/\s+/)[0]);
      if (ai > 0 && lines[p + 1] !== undefined) aliases.push([ai, lines[p + 1].trim()]);
      p++;
      continue;
    }
    if (!isV30(line)) continue;
    const content = line.slice(6).replace(/^ /, '');
    if (content.endsWith('-')) {
      buf += content.slice(0, -1);
      continue;
    }
    v30.push((buf + content).trim());
    buf = '';
  }
  if (buf) v30.push(buf.trim());
  return { v30, aliases, end: p };
}

/** Splits a V3000 line into tokens, keeping "quoted strings" and (parenthesised lists) together. */
function tokenizeV30(s: string): string[] {
  const out: string[] = [];
  let i = 0;
  const n = s.length;
  while (i < n) {
    while (i < n && s[i] === ' ') i++;
    if (i >= n) break;
    let tok = '';
    while (i < n && s[i] !== ' ') {
      if (s[i] === '"') {
        let j = i + 1;
        while (j < n && !(s[j] === '"' && s[j + 1] !== '"')) j += s[j] === '"' ? 2 : 1;
        tok += s.slice(i, Math.min(j + 1, n));
        i = j + 1;
      } else if (s[i] === '(' || s[i] === '[') {
        const close = s[i] === '(' ? ')' : ']';
        const j = s.indexOf(close, i);
        const e = j < 0 ? n : j + 1;
        tok += s.slice(i, e);
        i = e;
      } else tok += s[i++];
    }
    out.push(tok);
  }
  return out;
}

const unquote = (s: string) => (s.length >= 2 && s[0] === '"' && s[s.length - 1] === '"' ? s.slice(1, -1).replace(/""/g, '"') : s);
const listOf = (s: string): number[] => {
  const t = s.replace(/^\(|\)$/g, '').trim().split(/\s+/).map((x) => int(x));
  return t.slice(1, 1 + (t[0] || 0)).filter((x) => Number.isFinite(x));
};

function parseKV(tokens: string[]): Map<string, string> {
  const kv = new Map<string, string>();
  for (const t of tokens) {
    const k = t.indexOf('=');
    if (k > 0) kv.set(t.slice(0, k).toUpperCase(), t.slice(k + 1));
  }
  return kv;
}

function parseV3000(lines: string[], ci: number, fmt: string): { raw: RawCtab; end: number } {
  const { v30, aliases, end } = collectV30(lines, ci + 1);
  const k = v30.findIndex((l) => l.startsWith('BEGIN CTAB'));
  if (k < 0) throw new FormatError(fmt, 'V3000 file without "BEGIN CTAB"');
  const { raw } = parseV3000Ctab(v30, k, fmt);
  for (const [ai, text] of aliases) if (ai <= raw.mol.atoms.length && text) raw.aliases.set(ai - 1, text);
  return { raw, end };
}

/** Parses one V3000 CTAB from logical lines; v30[i] must be "BEGIN CTAB". Returns index after END CTAB. */
export function parseV3000Ctab(v30: string[], i: number, fmt: string): { raw: RawCtab; next: number } {
  const raw = newRaw();
  const mol = raw.mol;
  const idx = new Map<number, number>(); // file atom index → our index
  const bidx = new Map<number, number>();
  let section = '';
  let p = i + 1;
  for (; p < v30.length; p++) {
    const line = v30[p];
    if (line.startsWith('END CTAB')) {
      p++;
      break;
    }
    if (line.startsWith('BEGIN ')) {
      section = line.slice(6).trim().split(/\s+/)[0];
      continue;
    }
    if (line.startsWith('END ')) {
      section = '';
      continue;
    }
    const t = tokenizeV30(line);
    if (section === 'ATOM') {
      if (t[1] === 'NOT' && t.length > 2) t.splice(1, 2, 'NOT ' + t[2]); // negated atom list
      if (t.length < 5) throw new FormatError(fmt, `malformed V3000 atom line "${line}"`);
      const fileIdx = int(t[0]);
      const type = unquote(t[1]);
      const x = num(t[2]), y = num(t[3]), z = num(t[4]);
      if (!Number.isFinite(fileIdx) || !Number.isFinite(x) || !Number.isFinite(y)) throw new FormatError(fmt, `malformed V3000 atom line "${line}"`);
      const fields = /^\[|^NOT\b/i.test(type) ? { el: '*', alias: type } : atomFromSymbol(type);
      const ai = mol.addAtom({ ...fields, x, y, z: Number.isFinite(z) ? z : 0 });
      idx.set(fileIdx, ai);
      const a = mol.atoms[ai];
      const map = int(t[5]);
      if (map > 0) a.map = map;
      const kv = parseKV(t.slice(6));
      if (kv.has('CHG')) a.charge = int(kv.get('CHG')) || 0;
      if (kv.has('RAD')) setRadical(a, int(kv.get('RAD')));
      if (kv.has('MASS')) {
        const m = Math.round(num(kv.get('MASS')));
        if (m > 0) a.isotope = m;
      }
      const val = int(kv.get('VAL'));
      raw.valence[ai] = Number.isFinite(val) ? (val === -1 ? 15 : val) : 0;
      if (kv.has('RGROUPS') && a.el === 'R') {
        const r = listOf(kv.get('RGROUPS')!);
        if (r.length) a.alias = 'R' + r[0];
      }
    } else if (section === 'BOND') {
      if (t.length < 4) throw new FormatError(fmt, `malformed V3000 bond line "${line}"`);
      const a = idx.get(int(t[2])), b = idx.get(int(t[3]));
      if (a === undefined || b === undefined) throw new FormatError(fmt, `V3000 bond "${line}" references a missing atom`);
      const bi = addRawBond(raw, a, b, int(t[1]), fmt, p);
      if (bi < 0) continue;
      bidx.set(int(t[0]), bi);
      const cfg = int(parseKV(t.slice(4)).get('CFG'));
      const bd = mol.bonds[bi];
      if (bd.order === 1 && bd.style === 'plain') {
        if (cfg === 1) bd.style = 'wedge';
        else if (cfg === 3) bd.style = 'hash';
        else if (cfg === 2) bd.style = 'wavy';
      } else if (bd.order === 2 && cfg === 2) bd.style = 'crossed';
    } else if (section === 'SGROUP') {
      if (t[1] !== 'SUP') continue;
      const kv = parseKV(t.slice(3));
      const atoms = listOf(kv.get('ATOMS') ?? '').map((x) => idx.get(x) ?? -1);
      const label = unquote(kv.get('LABEL') ?? '');
      if (label) raw.sups.push({ atoms, label });
    }
  }
  return { raw, next: p };
}

/**
 * Post-processing after coordinates were normalised: kekulise aromatic bonds, apply the valence
 * field and aliases, contract known superatoms, perceive stereo from wedges/geometry.
 */
export function finishCtab(raw: RawCtab): Mol {
  const mol = raw.mol;
  mol.invalidate();
  kekulizeFlagged(mol, raw.aromatic);
  // valence field → explicit hydrogen count
  raw.valence.forEach((v, i) => {
    if (!v) return;
    const a = mol.atoms[i];
    if (a.el === 'R' || a.el === '*') return;
    a.hCount = v === 15 ? 0 : Math.max(0, v - Math.round(bondOrderSum(mol, i)));
  });
  dropDefaultHCounts(mol);
  for (const [i, text] of raw.aliases) applyAlias(mol.atoms[i], text);
  if (raw.sups.length) contractSuperatoms(mol, raw.sups);
  mol.invalidate();
  perceiveStereo(mol);
  return mol;
}

// ───────────────────────────── writing ─────────────────────────────

interface Prepared {
  mol: Mol;
  /** Superatom groups (atom indices) with labels. */
  sups: { atoms: number[]; label: string }[];
  is3D: boolean;
  /** Writable coordinates (Å, y up). */
  xyz: [number, number, number][];
  /** Bonds that are written (zero-order bonds skipped unless hbond in V3000). */
  chiral: boolean;
}

/** Expands abbreviations, assigns wedges if needed, converts coordinates. */
function prepare(input: Mol, opts: MolfileWriteOptions): Prepared {
  let mol = input;
  const sups: Prepared['sups'] = [];
  if (input.atoms.some((a) => a.abbrev)) {
    if (opts.expandAbbreviations !== false) {
      const tmp = input.clone();
      const labels = new Map<number, string>();
      tmp.atoms.forEach((a, i) => {
        a.id = i + 1;
        if (a.abbrev) labels.set(i, a.abbrev);
      });
      const n0 = tmp.atoms.length;
      mol = expandAbbreviations(tmp);
      for (const [i, label] of labels) {
        if (mol.atoms[i].el === '*' && mol.atoms[i].alias === label) continue; // unknown label → alias atom
        const atoms = [i];
        for (let j = n0; j < mol.atoms.length; j++) if (mol.atoms[j].id === i + 1) atoms.push(j);
        sups.push({ atoms, label });
      }
    } else {
      mol = input.clone();
      for (const a of mol.atoms) {
        if (!a.abbrev) continue;
        a.alias = a.abbrev;
        a.el = 'R';
        delete a.abbrev;
      }
    }
  }
  mol.invalidate();
  const threeD = mol.props.dim === '3D' || mol.atoms.some((a) => a.z !== undefined && Math.abs(a.z) > 1e-9);
  // Stereo given only as specs (e.g. from SMILES) but drawn without wedges: derive wedges.
  if (!threeD && mol.tetra.length && hasCoordinates(mol) && !mol.bonds.some((b) => b.style === 'wedge' || b.style === 'hash')) {
    mol = mol === input ? input.clone() : mol;
    try {
      assignWedgesFromSpecs(mol);
    } catch {
      /* best effort */
    }
  }
  const L = opts.bondLengthAngstrom ?? DEFAULT_BOND_LENGTH_ANGSTROM;
  const xyz: [number, number, number][] = mol.atoms.map((a) =>
    threeD ? [a.x, -a.y, a.z ?? 0] : [a.x * L, -a.y * L, 0],
  );
  const chiral = mol.bonds.some((b) => b.style === 'wedge' || b.style === 'hash');
  return { mol, sups, is3D: threeD, xyz, chiral };
}

/** Molfile bond type for a bond, or 0 to skip it. */
function bondType(order: number, style: BondStyle, v3000: boolean): number {
  if (style === 'hbond' || order === 0) return v3000 && style === 'hbond' ? 10 : 0;
  if (style === 'dative') return 9;
  if (order === 1.5) return 4;
  if (order === 2 || order === 3) return order;
  return 1;
}

/** Query-atom symbols that molfiles can store directly. */
const QUERY_SYMBOLS = new Set(['A', 'Q', 'L', 'LP', 'X', 'M', 'AH', 'QH', 'XH', 'MH']);

function atomSymbol(a: Atom): string {
  if (a.el === 'R') return a.alias && /^R\d+$/.test(a.alias) ? 'R#' : 'R';
  if (a.el === '*') return a.alias && QUERY_SYMBOLS.has(a.alias) ? a.alias : '*';
  return BY_SYMBOL.has(a.el) ? a.el : '*';
}

/** Valence-field value expressing a non-default explicit H count (0 = none). */
function valenceField(mol: Mol, i: number): number {
  const a = mol.atoms[i];
  if (a.hCount === undefined || a.el === 'R' || a.el === '*') return 0;
  if (!IMPLICIT_H_ELEMENTS.has(a.el) && a.hCount === 0) return 0;
  if (defaultImplicitH(mol, i) === a.hCount) return 0;
  const v = Math.round(bondOrderSum(mol, i)) + a.hCount;
  return v === 0 ? 15 : v > 14 ? 0 : v;
}

function headerLines(title: string, threeD: boolean): string[] {
  const d = new Date();
  const p2 = (v: number) => String(v).padStart(2, '0');
  const stamp = p2(d.getMonth() + 1) + p2(d.getDate()) + p2(d.getFullYear() % 100) + p2(d.getHours()) + p2(d.getMinutes());
  return [title.replace(/[\r\n]+/g, ' ').slice(0, 80), '  ChemWrit' + stamp + (threeD ? '3D' : '2D'), ''];
}

/** Writes a molfile (V2000 by default; V3000 when requested or for > 999 atoms/bonds). */
export function writeMolfile(mol: Mol, opts: MolfileWriteOptions = {}): string {
  const prep = prepare(mol, opts);
  const big = prep.mol.atoms.length > 999 || prep.mol.bonds.length > 999;
  const title = opts.title ?? mol.name ?? '';
  const head = headerLines(title, prep.is3D);
  if (opts.v3000 || big) {
    return [...head, '  0  0  0     0  0            999 V3000', ...v3000CtabLines(prep), 'M  END', ''].join('\n');
  }
  return [...head, ...v2000Body(prep), ''].join('\n');
}

function v2000Body(prep: Prepared): string[] {
  const { mol, xyz } = prep;
  const out: string[] = [];
  const bonds: number[] = [];
  mol.bonds.forEach((b, bi) => bondType(b.order, b.style, false) && bonds.push(bi));
  const bondNo = new Map<number, number>();
  bonds.forEach((bi, k) => bondNo.set(bi, k + 1));
  out.push(pad(mol.atoms.length, 3) + pad(bonds.length, 3) + '  0  0' + pad(prep.chiral ? 1 : 0, 3) + '  0  0  0  0  0999 V2000');
  const chg: [number, number][] = [], rad: [number, number][] = [], iso: [number, number][] = [], rgp: [number, number][] = [];
  const aliases: [number, string][] = [];
  mol.atoms.forEach((a, i) => {
    const sym = atomSymbol(a);
    const c = a.charge >= -3 && a.charge <= 3 && a.charge !== 0 ? 4 - a.charge : 0;
    const [x, y, z] = xyz[i];
    out.push(
      fixed(x, 4).padStart(10) + fixed(y, 4).padStart(10) + fixed(z, 4).padStart(10) + ' ' + sym.padEnd(3) + ' 0' + pad(c, 3) +
        '  0  0  0' + pad(valenceField(mol, i), 3) + '  0  0  0' + pad(a.map ?? 0, 3) + '  0  0',
    );
    if (a.charge) chg.push([i + 1, a.charge]);
    if (a.radical) rad.push([i + 1, a.radical === 1 ? 2 : 3]);
    if (a.isotope) iso.push([i + 1, a.isotope]);
    if (sym === 'R#') rgp.push([i + 1, +a.alias!.slice(1)]);
    else if ((a.el === 'R' || a.el === '*') && a.alias && sym !== a.alias) aliases.push([i + 1, a.alias]);
  });
  for (const bi of bonds) {
    const b = mol.bonds[bi];
    const type = bondType(b.order, b.style, false);
    let st = 0;
    if (type === 1) st = b.style === 'wedge' ? 1 : b.style === 'hash' ? 6 : b.style === 'wavy' ? 4 : 0;
    else if (type === 2 && b.style === 'crossed') st = 3;
    out.push(pad(b.a + 1, 3) + pad(b.b + 1, 3) + pad(type, 3) + pad(st, 3) + '  0  0  0');
  }
  const prop = (tag: string, list: [number, number][]) => {
    for (let k = 0; k < list.length; k += 8) {
      const chunk = list.slice(k, k + 8);
      out.push(`M  ${tag}` + pad(chunk.length, 3) + chunk.map(([i, v]) => ' ' + pad(i, 3) + ' ' + pad(v, 3)).join(''));
    }
  };
  for (const [i, text] of aliases) out.push('A  ' + pad(i, 3), text);
  prop('CHG', chg);
  prop('RAD', rad);
  prop('ISO', iso);
  prop('RGP', rgp);
  if (prep.sups.length) {
    for (let k = 0; k < prep.sups.length; k += 8) {
      const chunk = prep.sups.slice(k, k + 8);
      out.push('M  STY' + pad(chunk.length, 3) + chunk.map((_, j) => ' ' + pad(k + j + 1, 3) + ' SUP').join(''));
    }
    prep.sups.forEach((g, k) => {
      const set = new Set(g.atoms);
      const xb = bonds.filter((bi) => set.has(mol.bonds[bi].a) !== set.has(mol.bonds[bi].b)).map((bi) => bondNo.get(bi)!);
      for (let j = 0; j < g.atoms.length; j += 15) {
        const chunk = g.atoms.slice(j, j + 15);
        out.push('M  SAL ' + pad(k + 1, 3) + pad(chunk.length, 3) + chunk.map((i) => ' ' + pad(i + 1, 3)).join(''));
      }
      for (let j = 0; j < xb.length; j += 15) {
        const chunk = xb.slice(j, j + 15);
        out.push('M  SBL ' + pad(k + 1, 3) + pad(chunk.length, 3) + chunk.map((b) => ' ' + pad(b, 3)).join(''));
      }
      out.push('M  SMT ' + pad(k + 1, 3) + ' ' + g.label);
    });
  }
  out.push('M  END');
  return out;
}

/** Wraps a V3000 line at 80 columns using '-' continuation. */
function v30(content: string): string[] {
  const out: string[] = [];
  let rest = content;
  while (rest.length > 72) {
    out.push('M  V30 ' + rest.slice(0, 72) + '-');
    rest = rest.slice(72);
  }
  out.push('M  V30 ' + rest);
  return out;
}

const quoteV30 = (s: string) => (/[\s"()=]/.test(s) || s === '' ? '"' + s.replace(/"/g, '""') + '"' : s);

/** "BEGIN CTAB" … "END CTAB" lines (V3000) for a prepared molecule. */
function v3000CtabLines(prep: Prepared): string[] {
  const { mol, xyz } = prep;
  const out: string[] = [];
  const bonds: number[] = [];
  mol.bonds.forEach((b, bi) => bondType(b.order, b.style, true) && bonds.push(bi));
  const bondNo = new Map<number, number>();
  bonds.forEach((bi, k) => bondNo.set(bi, k + 1));
  out.push('M  V30 BEGIN CTAB');
  out.push(...v30(`COUNTS ${mol.atoms.length} ${bonds.length} ${prep.sups.length} 0 ${prep.chiral ? 1 : 0}`));
  out.push('M  V30 BEGIN ATOM');
  mol.atoms.forEach((a, i) => {
    let sym = atomSymbol(a);
    if ((a.el === 'R' || a.el === '*') && a.alias && sym !== 'R#' && sym !== a.alias) sym = quoteV30(a.alias);
    const [x, y, z] = xyz[i];
    let s = `${i + 1} ${sym} ${fixed(x, 4)} ${fixed(y, 4)} ${fixed(z, 4)} ${a.map ?? 0}`;
    if (a.charge) s += ` CHG=${a.charge}`;
    if (a.radical) s += ` RAD=${a.radical === 1 ? 2 : 3}`;
    if (a.isotope) s += ` MASS=${a.isotope}`;
    const v = valenceField(mol, i);
    if (v) s += ` VAL=${v === 15 ? -1 : v}`;
    if (sym === 'R#') s += ` RGROUPS=(1 ${+a.alias!.slice(1)})`;
    out.push(...v30(s));
  });
  out.push('M  V30 END ATOM');
  if (bonds.length) {
    out.push('M  V30 BEGIN BOND');
    for (const bi of bonds) {
      const b = mol.bonds[bi];
      const type = bondType(b.order, b.style, true);
      let cfg = 0;
      if (type === 1) cfg = b.style === 'wedge' ? 1 : b.style === 'hash' ? 3 : b.style === 'wavy' ? 2 : 0;
      else if (type === 2 && b.style === 'crossed') cfg = 2;
      out.push(...v30(`${bondNo.get(bi)} ${type} ${b.a + 1} ${b.b + 1}${cfg ? ` CFG=${cfg}` : ''}`));
    }
    out.push('M  V30 END BOND');
  }
  if (prep.sups.length) {
    out.push('M  V30 BEGIN SGROUP');
    prep.sups.forEach((g, k) => {
      const set = new Set(g.atoms);
      const xb = bonds.filter((bi) => set.has(mol.bonds[bi].a) !== set.has(mol.bonds[bi].b)).map((bi) => bondNo.get(bi)!);
      let s = `${k + 1} SUP 0 ATOMS=(${g.atoms.length} ${g.atoms.map((i) => i + 1).join(' ')})`;
      if (xb.length) s += ` XBONDS=(${xb.length} ${xb.join(' ')})`;
      s += ` LABEL=${quoteV30(g.label)}`;
      out.push(...v30(s));
    });
    out.push('M  V30 END SGROUP');
  }
  out.push('M  V30 END CTAB');
  return out;
}

/** CTAB lines of a molecule for embedding in V3000 RXN files. */
export function v3000CtabFor(mol: Mol, opts: MolfileWriteOptions = {}): string[] {
  return v3000CtabLines(prepare(mol, opts));
}

/** Writes an SD file; mol.name becomes the title, mol.props the data fields. */
export function writeSDF(mols: Mol[]): string {
  let out = '';
  for (const m of mols) {
    out += writeMolfile(m, { title: m.name });
    for (const [k, v] of Object.entries(m.props)) {
      if (k === 'dim') continue; // internal flag
      const value = String(v).replace(/\r\n?/g, '\n').replace(/\n\s*\n/g, '\n').replace(/\n+$/, '');
      out += `> <${k}>\n${value}\n\n`;
    }
    out += '$$$$\n';
  }
  return out;
}
