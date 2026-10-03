// Name assembly helpers: alphanumerical ordering, enclosing marks, multiplying prefixes, elision.
import { multiplier, multiplierComplex } from './numerals';

const MULT_START = /^(di|tri|tetra|penta|hexa|hepta|octa|nona|deca|undeca|dodeca|icosa|bis|tris)/;

/** A substituent prefix as produced by the namer. */
export interface SubName {
  /** Prefix text without outer enclosing marks, e.g. "methyl", "2-methylpropyl", "(2-methylpropan-2-yl)oxy". */
  text: string;
  /** Contains substituent prefixes of its own → "bis/tris" when multiplied. */
  compound: boolean;
  /** Must be enclosed when cited after a locant (amino-like prefixes such as "methylamino"). */
  enclose: boolean;
  /** "…yl" with free valence at C1 of a chain or phenyl: may contract with "oxy" (methoxy, phenoxy). */
  oxyContract?: boolean;
  /** Unsubstituted phenyl (for "anilino"/"benzyl"). */
  phenyl?: boolean;
  /** Ring-substituted phenyl ("4-chlorophenyl") – for "anilino". */
  phenylLike?: boolean;
  /** Locant of the free valence in the substituent's own numbering (skeletal substituents). */
  fvLoc?: number;
}

/** Sort key for alphanumerical ordering of prefixes (P-14.5): italic parts, locants and marks ignored. */
export function alphaKey(text: string): string {
  return text
    .replace(/\b(tert|sec|cis|trans)-/g, '')
    .replace(/\d+[a-z]?H(?=[-,)\]}])/g, '')
    .replace(/(^|[^A-Za-z])[A-Z]'*(?=[,-])/g, '$1')
    .replace(/[^A-Za-z]/g, '')
    .toLowerCase();
}

/** Maximum nesting depth of enclosing marks in a string. */
function bracketDepth(s: string): number {
  let d = 0, max = 0;
  for (const c of s) {
    if (c === '(' || c === '[' || c === '{') max = Math.max(max, ++d);
    else if (c === ')' || c === ']' || c === '}') d--;
  }
  return max;
}

/** Encloses text in the next level of enclosing marks: ( ) → [ ] → { } → ( ) … */
export function enclose(s: string): string {
  const d = bracketDepth(s);
  const k = d % 3;
  return (['(', '[', '{'][k]) + s + ([')', ']', '}'][k]);
}

/** True if a prefix needs marks to be separated from a preceding locant. */
export function startsAmbiguous(s: string): boolean {
  return /^[\d([{]/.test(s) || /^[A-Z]'*[,-]/.test(s);
}

export function startsWithVowel(s: string): boolean {
  return /^[aeiouy]/.test(s);
}

/** Drops the final "e" of a parent name before a suffix beginning with a vowel or "y". */
export function elide(parent: string, suffix: string): string {
  let first = suffix;
  if (first.startsWith('-')) {
    const k = first.indexOf('-', 1);
    if (k > 0) first = first.slice(k + 1);
  }
  if (parent.endsWith('e') && startsWithVowel(first)) return parent.slice(0, -1) + suffix;
  return parent + suffix;
}

export interface PrefixItem {
  sub: SubName;
  /** Locant label, '' when locants are omitted. */
  locant: string;
  /** Sort value of the locant (for ordering locants inside a group). */
  value: number;
}

/**
 * Formats the detachable prefixes: groups identical prefixes, adds multiplying prefixes and
 * enclosing marks, sorts alphanumerically and joins them.
 * `locantless` – the prefixes carry no locants by nature (methane-type parents, substituents of
 * amino/carbamoyl/silyl groups): consecutive simple prefixes are separated by enclosing marks
 * ("dichloro(fluoro)methane", "ethyl(methyl)amino").
 */
export function formatPrefixes(items: PrefixItem[], locantless: boolean): string {
  const groups = new Map<string, PrefixItem[]>();
  for (const it of items) {
    const k = it.sub.text;
    if (!groups.has(k)) groups.set(k, []);
    groups.get(k)!.push(it);
  }
  const entries = [...groups.values()].map((list) => {
    list.sort((a, b) => a.value - b.value || (a.locant < b.locant ? -1 : a.locant > b.locant ? 1 : 0));
    const sub = list[0].sub;
    const t = sub.text;
    const n = list.length;
    const locs = list.map((x) => x.locant).filter((x) => x !== '');
    const hasLoc = locs.length > 0;
    const digitStart = /^\d/.test(t) || /^[A-Z]'*[,-]/.test(t);
    const bracketStart = /^[([{]/.test(t);
    let body: string;
    if (n > 1 && sub.compound) body = multiplierComplex(n) + enclose(t);
    // simple prefixes are enclosed after a multiplier if they contain locants/marks or themselves begin
    // with a multiplier-like syllable ("di(propan-2-yl)", "di(hexadecanoyloxy)", P-16.5.1)
    else if (n > 1) body = multiplier(n) + (/[\d([{]/.test(t) || sub.enclose || MULT_START.test(t) ? enclose(t) : t);
    else if (hasLoc) body = digitStart || bracketStart || sub.enclose ? enclose(t) : t;
    else if (locantless) body = digitStart || sub.enclose ? enclose(t) : t;
    else body = /[([{]/.test(t.slice(1)) && !bracketStart ? enclose(t) : t;
    return { key: alphaKey(t), text: t, locs, body };
  });
  entries.sort((a, b) => (a.key < b.key ? -1 : a.key > b.key ? 1 : a.text < b.text ? -1 : a.text > b.text ? 1 : 0));
  let out = '';
  let sepParen = false; // the last piece was a simple prefix enclosed only to separate it
  for (const e of entries) {
    const piece = (e.locs.length ? e.locs.join(',') + '-' : '') + e.body;
    const wasSep = sepParen;
    sepParen = false;
    if (!out) {
      out = piece;
      continue;
    }
    if (e.locs.length) out += '-' + piece;
    else if (wasSep && !/^[([{]/.test(piece)) out += piece;
    else if (/[)\]}]$/.test(out)) out += '-' + piece;
    else if (locantless && !/^[([{]/.test(piece)) {
      out += enclose(piece);
      sepParen = true;
    } else if (locantless) out += '-' + piece;
    else out += piece;
  }
  return out;
}

/** "2,3-" style locant list. */
export function locList(locs: string[]): string {
  return locs.join(',');
}
