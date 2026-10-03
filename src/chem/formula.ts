import { element, ISOTOPE_MASS } from './elements';
import { Mol } from './mol';
import { implicitH } from './valence';
import { expandAbbreviations } from './abbreviations';

export interface FormulaInfo {
  /** element symbol (isotopes as e.g. "[13C]" or "D") → count */
  counts: Map<string, number>;
  charge: number;
  /** Hill-order formula, plain text, e.g. "C9H8O4" */
  formula: string;
  /** Hill formula as HTML with <sub>/<sup> */
  html: string;
  /** Average molecular weight (g/mol) */
  mw: number;
  /** Monoisotopic exact mass */
  exactMass: number;
  /** m/z of the (charged) species = exactMass / |charge| (electron mass corrected) */
  mz: number;
  /** Mass percentages by element */
  composition: { el: string; pct: number }[];
  /** Degrees of unsaturation (rings + π bonds) for neutral C/H/N/O/X/S compounds, null if undefined */
  dbe: number | null;
  /** true if the structure contains pseudo atoms (R, *, unknown labels) */
  hasPseudo: boolean;
}

const ELECTRON_MASS = 0.00054857990946;

function isoKey(el: string, iso?: number): string {
  if (!iso) return el;
  if (el === 'H' && iso === 2) return 'D';
  if (el === 'H' && iso === 3) return 'T';
  return `[${iso}${el}]`;
}

function massOf(key: string, avg: boolean): number {
  if (key === 'D') return ISOTOPE_MASS['H2'];
  if (key === 'T') return ISOTOPE_MASS['H3'];
  const m = /^\[(\d+)([A-Za-z]+)\]$/.exec(key);
  if (m) {
    const exact = ISOTOPE_MASS[m[2] + m[1]];
    return exact ?? +m[1];
  }
  const e = element(key);
  if (!e) return 0;
  return avg ? e.mass : e.mono;
}

/** Computes formula, masses and composition for all atoms of `mol` (abbreviations expanded). */
export function computeFormula(input: Mol): FormulaInfo {
  const mol = expandAbbreviations(input);
  const counts = new Map<string, number>();
  let charge = 0;
  let hasPseudo = false;
  const add = (k: string, n: number) => counts.set(k, (counts.get(k) ?? 0) + n);
  for (let i = 0; i < mol.atoms.length; i++) {
    const a = mol.atoms[i];
    if (!element(a.el)) {
      hasPseudo = true;
      continue;
    }
    add(isoKey(a.el, a.isotope), 1);
    const h = implicitH(mol, i);
    if (h) add('H', h);
    charge += a.charge;
  }
  return summarize(counts, charge, hasPseudo, mol);
}

/** Builds FormulaInfo from element counts (e.g. parsed from a formula string). */
export function summarize(counts: Map<string, number>, charge: number, hasPseudo = false, mol?: Mol): FormulaInfo {
  const keys = [...counts.keys()];
  const hasC = keys.some((k) => k === 'C' || /\dC\]$/.test(k));
  const order = (k: string): [number, string] => {
    const base = k.replace(/^\[\d+/, '').replace(/\]$/, '');
    const isC = base === 'C';
    const isH = base === 'H' || k === 'D' || k === 'T';
    if (hasC && isC) return [0, k];
    if (hasC && isH) return [1, k === 'H' ? 'H' : 'H' + k];
    return [2, base + k];
  };
  keys.sort((a, b) => {
    const oa = order(a), ob = order(b);
    if (oa[0] !== ob[0]) return oa[0] - ob[0];
    return oa[1] < ob[1] ? -1 : oa[1] > ob[1] ? 1 : 0;
  });
  let formula = '';
  let html = '';
  let mw = 0;
  let exact = 0;
  for (const k of keys) {
    const n = counts.get(k)!;
    formula += k + (n > 1 ? n : '');
    const kh = k.replace(/^\[(\d+)([A-Za-z]+)\]$/, '<sup>$1</sup>$2');
    html += kh + (n > 1 ? `<sub>${n}</sub>` : '');
    mw += massOf(k, true) * n;
    exact += massOf(k, false) * n;
  }
  if (charge) {
    const c = (Math.abs(charge) > 1 ? Math.abs(charge) : '') + (charge > 0 ? '+' : '−');
    formula += c.replace('−', '-');
    html += `<sup>${c}</sup>`;
  }
  const composition = keys.map((k) => ({ el: k, pct: mw ? (100 * massOf(k, true) * counts.get(k)!) / mw : 0 }));
  const exactCharged = exact - charge * ELECTRON_MASS;
  const mz = charge ? exactCharged / Math.abs(charge) : exact;
  // DBE = C - H/2 - X/2 + N/2 + 1 (for C, H, N, O, S, halogens, Si, P(3))
  let dbe: number | null = null;
  const tally = (els: string[]) => els.reduce((s, e) => s + (counts.get(e) ?? 0), 0);
  const allowed = new Set(['C', 'H', 'D', 'T', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'I', 'Si', 'P', 'Se']);
  if (keys.every((k) => allowed.has(k.replace(/^\[\d+/, '').replace(/\]$/, '')))) {
    const C = tally(['C']) + tally(['Si']) + keys.filter((k) => /\dC\]$/.test(k)).reduce((s, k) => s + counts.get(k)!, 0);
    const H = tally(['H', 'D', 'T']);
    const X = tally(['F', 'Cl', 'Br', 'I']);
    const N = tally(['N', 'P']);
    dbe = C - (H + X) / 2 + N / 2 + 1;
    if (charge) dbe += 0; // charged species: report the neutral formula value
  }
  void mol;
  return { counts, charge, formula, html, mw, exactMass: exactCharged, mz, composition, dbe, hasPseudo };
}

/** Parses a molecular formula string such as "C6H12O6", "CH3COO-", "C2H5OH", "Ca(OH)2", "[13C]H4". */
export function parseFormula(text: string): FormulaInfo {
  const s = text.replace(/\s+/g, '');
  let i = 0;
  const parseGroup = (): Map<string, number> => {
    const m = new Map<string, number>();
    while (i < s.length) {
      const c = s[i];
      if (c === '(') {
        i++;
        const inner = parseGroup();
        if (s[i] !== ')') throw new Error('Unbalanced parenthesis');
        i++;
        const n = readNum();
        for (const [k, v] of inner) m.set(k, (m.get(k) ?? 0) + v * n);
      } else if (c === ')') {
        return m;
      } else if (c === '[') {
        const end = s.indexOf(']', i);
        const tok = s.slice(i, end + 1);
        i = end + 1;
        const n = readNum();
        m.set(tok, (m.get(tok) ?? 0) + n);
      } else if (/[A-Z]/.test(c)) {
        let sym = c;
        if (/[a-z]/.test(s[i + 1] ?? '')) sym += s[i + 1];
        if (!element(sym) && sym.length === 2 && element(c)) sym = c;
        i += sym.length;
        if (!element(sym) && sym !== 'D' && sym !== 'T') throw new Error(`Unknown element ${sym}`);
        const n = readNum();
        m.set(sym, (m.get(sym) ?? 0) + n);
      } else if (c === '+' || c === '-' || c === '·' || c === '.' || c === '^') {
        return m;
      } else {
        throw new Error(`Unexpected '${c}'`);
      }
    }
    return m;
  };
  const readNum = (): number => {
    const m = /^\d+/.exec(s.slice(i));
    if (!m) return 1;
    i += m[0].length;
    return +m[0];
  };
  const counts = parseGroup();
  let charge = 0;
  const rest = s.slice(i).replace(/^\^/, '');
  const cm = /^(\d*)([+-])$|^([+-]+)$|^([+-])(\d+)$/.exec(rest);
  if (cm) {
    if (cm[2]) charge = (cm[2] === '+' ? 1 : -1) * (cm[1] ? +cm[1] : 1);
    else if (cm[3]) charge = (cm[3][0] === '+' ? 1 : -1) * cm[3].length;
    else charge = (cm[4] === '+' ? 1 : -1) * +cm[5];
  }
  return summarize(counts, charge);
}
