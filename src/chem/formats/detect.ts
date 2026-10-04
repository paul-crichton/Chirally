// File-format detection by extension and content sniffing.
import { parseSmiles, parseReactionSmiles } from '../smiles';

export type FormatId = 'chirally' | 'mol' | 'sdf' | 'rxn' | 'cdxml' | 'cml' | 'xyz' | 'smiles' | 'unknown';

const BY_EXTENSION: Record<string, FormatId> = {
  mol: 'mol', mdl: 'mol', sdf: 'sdf', sd: 'sdf', rxn: 'rxn', cdxml: 'cdxml', cml: 'cml', xyz: 'xyz',
  smi: 'smiles', smiles: 'smiles', chirally: 'chirally', json: 'chirally',
};

/**
 * Guesses the format of a file. A known extension wins (except that a ".mol" file holding several
 * records is reported as SDF, and a ".json" file that is not a Chirally document as unknown);
 * otherwise the content is sniffed.
 */
export function detectFormat(fileName: string | null, text: string): FormatId {
  const ext = fileName ? /\.([A-Za-z0-9]+)$/.exec(fileName.trim())?.[1]?.toLowerCase() : undefined;
  const byExt = ext ? BY_EXTENSION[ext] : undefined;
  const sniffed = sniff(text ?? '');
  if (byExt) {
    if (byExt === 'mol' && (sniffed === 'sdf' || sniffed === 'rxn')) return sniffed;
    if (byExt === 'chirally' && ext === 'json' && sniffed !== 'chirally') return 'unknown';
    return byExt;
  }
  return sniffed;
}

function sniff(raw: string): FormatId {
  const text = raw.charCodeAt(0) === 0xfeff ? raw.slice(1) : raw;
  const head = text.slice(0, 4096);
  const trimmed = head.trimStart();
  if (!trimmed) return 'unknown';
  if (trimmed.startsWith('{')) return /"format"\s*:\s*"chirally"/.test(head) ? 'chirally' : 'unknown';
  if (trimmed.startsWith('$RXN')) return 'rxn';
  if (trimmed.startsWith('<')) {
    if (/<CDXML[\s>]/.test(head)) return 'cdxml';
    if (/<(?:\w+:)?(?:cml|molecule)[\s>]/.test(head)) return 'cml';
    return 'unknown';
  }
  const lines = head.replace(/\r\n?/g, '\n').split('\n');
  // molfile: counts line (normally the 4th) carrying V2000/V3000, or a fixed-width counts line
  const counts = lines.findIndex((l, i) => i <= 6 && /^\s*\d+\s*\d+.*V[23]000\s*$/.test(l));
  const fixedCounts = lines.length > 3 && /^[ \d]{2}\d[ \d]{2}\d/.test(lines[3]) && /^\s*\d+\s+\d+/.test(lines[3]);
  if (counts >= 0 || fixedCounts) return /^\$\$\$\$/m.test(text) ? 'sdf' : 'mol';
  if (/^\$\$\$\$/m.test(text) && /M {2}END/.test(text)) return 'sdf';
  // XYZ: atom count, comment, then "El x y z"
  if (/^\s*\d+\s*$/.test(lines[0]) && lines.length > 2 && /^\s*([A-Za-z]{1,3}|\d{1,3})\s+[-+]?\d*\.?\d+([eE][-+]?\d+)?\s+[-+]?\d*\.?\d+([eE][-+]?\d+)?\s+[-+]?\d*\.?\d+/.test(lines[2])) {
    return 'xyz';
  }
  // SMILES: first token of the first line (a name may follow after whitespace)
  const first = lines[0].trim().split(/\s+/)[0];
  if (first && first.length < 5000 && /^[A-Za-z0-9@+\-[\]()=#$:/\\%.*>~]+$/.test(first) && /[A-Za-z*]/.test(first)) {
    try {
      if (first.includes('>')) parseReactionSmiles(first);
      else parseSmiles(first);
      return 'smiles';
    } catch {
      return 'unknown';
    }
  }
  return 'unknown';
}
