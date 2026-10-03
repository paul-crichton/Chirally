// Text measurement and light-weight rich text (sub/superscripts) shared by canvas and SVG output.

export interface TextRun {
  text: string;
  sub?: boolean;
  sup?: boolean;
  bold?: boolean;
  italic?: boolean;
}

let measureCtx: CanvasRenderingContext2D | OffscreenCanvasRenderingContext2D | null | undefined;
const widthCache = new Map<string, number>();

// Helvetica advance widths (per 1000 em) for ASCII 32..126, used when no canvas is available (tests / SSR).
const HELV = [
  278, 278, 355, 556, 556, 889, 667, 191, 333, 333, 389, 584, 278, 333, 278, 278, 556, 556, 556, 556, 556, 556, 556, 556,
  556, 556, 278, 278, 584, 584, 584, 556, 1015, 667, 667, 722, 722, 667, 611, 778, 722, 278, 500, 667, 556, 833, 722,
  778, 667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611, 278, 278, 278, 469, 556, 333, 556, 556, 500, 556, 556,
  278, 556, 556, 222, 222, 500, 222, 833, 556, 556, 556, 556, 333, 500, 278, 556, 500, 722, 500, 500, 500, 334, 260,
  334, 584,
];

function getCtx() {
  if (measureCtx !== undefined) return measureCtx;
  try {
    if (typeof OffscreenCanvas !== 'undefined') measureCtx = new OffscreenCanvas(8, 8).getContext('2d');
    else if (typeof document !== 'undefined') measureCtx = document.createElement('canvas').getContext('2d');
    else measureCtx = null;
  } catch {
    measureCtx = null;
  }
  return measureCtx;
}

export function fontString(size: number, family: string, bold = false, italic = false): string {
  return `${italic ? 'italic ' : ''}${bold ? 'bold ' : ''}${size}px ${family}`;
}

/** Width of `text` at font size `size` (any unit; result in the same unit). */
export function measureText(text: string, size: number, family: string, bold = false, italic = false): number {
  const key = `${bold ? 'b' : ''}${italic ? 'i' : ''}|${family}|${text}`;
  let w = widthCache.get(key);
  if (w === undefined) {
    const ctx = getCtx();
    if (ctx) {
      ctx.font = fontString(100, family, bold, italic);
      w = ctx.measureText(text).width / 100;
    } else {
      w = 0;
      for (const ch of text) {
        const c = ch.charCodeAt(0);
        w += (c >= 32 && c <= 126 ? HELV[c - 32] : 600) / 1000;
      }
      if (bold) w *= 1.06;
    }
    widthCache.set(key, w);
  }
  return w * size;
}

/**
 * Parses light markup: `_x` / `_{xy}` subscript, `^x` / `^{xy}` superscript, `**bold**`, `*italic*`.
 * Backslash escapes the next character.
 */
export function parseMarkup(s: string): TextRun[] {
  const runs: TextRun[] = [];
  let bold = false;
  let italic = false;
  let buf = '';
  const flush = () => {
    if (buf) runs.push({ text: buf, bold: bold || undefined, italic: italic || undefined });
    buf = '';
  };
  for (let i = 0; i < s.length; i++) {
    const c = s[i];
    if (c === '\\' && i + 1 < s.length) {
      buf += s[++i];
    } else if (c === '*' && s[i + 1] === '*') {
      flush();
      bold = !bold;
      i++;
    } else if (c === '*') {
      flush();
      italic = !italic;
    } else if ((c === '_' || c === '^') && i + 1 < s.length) {
      flush();
      let t: string;
      if (s[i + 1] === '{') {
        const end = s.indexOf('}', i + 2);
        t = end < 0 ? s.slice(i + 2) : s.slice(i + 2, end);
        i = end < 0 ? s.length : end;
      } else {
        // a run of digits / +/- counts as one script token
        const m = /^[0-9]+[+\-−]?|^[+\-−]|^./.exec(s.slice(i + 1))!;
        t = m[0];
        i += t.length;
      }
      runs.push({ text: t, sub: c === '_' || undefined, sup: c === '^' || undefined, bold: bold || undefined, italic: italic || undefined });
    } else {
      buf += c;
    }
  }
  flush();
  return runs;
}

/** Formats a chemical formula: digits after an element/bracket become subscripts, trailing charge superscript. */
export function formulaRuns(s: string): TextRun[] {
  const runs: TextRun[] = [];
  const m = /^(.*?)(\d*[+\-−])?$/.exec(s)!;
  let body = m[1];
  const charge = m[2];
  if (charge && /\d$/.test(body) === false && body.length === 0) body = s;
  let buf = '';
  let prevLetter = false;
  for (let i = 0; i < body.length; i++) {
    const c = body[i];
    if (/\d/.test(c) && prevLetter) {
      if (buf) runs.push({ text: buf });
      buf = '';
      let d = c;
      while (i + 1 < body.length && /\d/.test(body[i + 1])) d += body[++i];
      runs.push({ text: d, sub: true });
      prevLetter = true;
      continue;
    }
    buf += c;
    prevLetter = /[A-Za-z)\]]/.test(c);
  }
  if (buf) runs.push({ text: buf });
  if (charge && body !== s) runs.push({ text: charge.replace('-', '−'), sup: true });
  return runs;
}

/** Splits text into lines and runs (markup or formula mode). */
export function richLines(text: string, formula = false): TextRun[][] {
  return text.split('\n').map((line) => (formula ? formulaRuns(line) : parseMarkup(line)));
}

export const SCRIPT_SCALE = 0.7;
export const SUB_SHIFT = 0.28; // fraction of font size, downward
export const SUP_SHIFT = 0.42; // upward

export function runsWidth(runs: TextRun[], size: number, family: string): number {
  let w = 0;
  for (const r of runs) {
    const s = r.sub || r.sup ? size * SCRIPT_SCALE : size;
    w += measureText(r.text, s, family, r.bold, r.italic);
  }
  return w;
}
