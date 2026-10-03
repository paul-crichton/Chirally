// Minimal, dependency-free XML parser and writer helpers (no DOMParser, so it runs in Node too).
//
// Supports elements, attributes ('…' or "…"), self-closing tags, text, CDATA, comments,
// processing instructions, DOCTYPE (incl. internal subset) and the predefined/numeric entities.
// The parser is linear-time (index scanning, no backtracking regexes) and lenient about
// mismatched closing tags so that slightly broken files from other programs still load.
import { FormatError } from './common';

export interface XmlElement {
  name: string;
  attrs: Record<string, string>;
  children: XmlNode[];
}

export type XmlNode = XmlElement | string;

const NAMED_ENTITIES: Record<string, string> = { lt: '<', gt: '>', amp: '&', quot: '"', apos: "'", nbsp: ' ' };

/** Replaces XML entities (&amp; &#65; &#x41; …) in raw text. Unknown entities are kept verbatim. */
export function decodeEntities(s: string): string {
  if (s.indexOf('&') < 0) return s;
  return s.replace(/&(#x[0-9a-fA-F]+|#[0-9]+|[A-Za-z][A-Za-z0-9]*);/g, (m, ent: string) => {
    if (ent[0] === '#') {
      const code = ent[1] === 'x' || ent[1] === 'X' ? parseInt(ent.slice(2), 16) : parseInt(ent.slice(1), 10);
      if (!Number.isFinite(code) || code < 0 || code > 0x10ffff) return m;
      try {
        return String.fromCodePoint(code);
      } catch {
        return m;
      }
    }
    return NAMED_ENTITIES[ent] ?? m;
  });
}

/** Escapes text for use in XML content or attribute values. */
export function escapeXml(s: string): string {
  return s
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/\n/g, '&#10;')
    .replace(/\r/g, '&#13;');
}

/** Escapes text content (newlines kept literally). */
export function escapeText(s: string): string {
  return s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}

const isNameChar = (c: number) =>
  (c >= 48 && c <= 57) || // 0-9
  (c >= 65 && c <= 90) || // A-Z
  (c >= 97 && c <= 122) || // a-z
  c === 95 || c === 58 || c === 45 || c === 46 || // _ : - .
  c > 127;

const isSpace = (c: number) => c === 32 || c === 9 || c === 10 || c === 13;

/**
 * Parses an XML document and returns its root element. Throws FormatError on input that does not
 * contain a well-formed root element.
 */
export function parseXml(text: string, formatName = 'XML'): XmlElement {
  const s = text.charCodeAt(0) === 0xfeff ? text.slice(1) : text;
  const n = s.length;
  const fail = (msg: string, pos: number): never => {
    const line = s.slice(0, pos).split('\n').length;
    throw new FormatError(formatName, `${msg} (line ${line})`);
  };
  // Synthetic document node collects top-level elements.
  const doc: XmlElement = { name: '#document', attrs: {}, children: [] };
  const stack: XmlElement[] = [doc];
  let i = 0;

  while (i < n) {
    const lt = s.indexOf('<', i);
    const textEnd = lt < 0 ? n : lt;
    if (textEnd > i) {
      const raw = s.slice(i, textEnd);
      if (stack.length > 1) stack[stack.length - 1].children.push(decodeEntities(raw));
      else if (raw.trim()) fail('Text outside of the root element', i);
    }
    if (lt < 0) break;
    i = lt;
    if (s.startsWith('<!--', i)) {
      const end = s.indexOf('-->', i + 4);
      if (end < 0) fail('Unterminated comment', i);
      i = end + 3;
    } else if (s.startsWith('<![CDATA[', i)) {
      const end = s.indexOf(']]>', i + 9);
      if (end < 0) fail('Unterminated CDATA section', i);
      if (stack.length > 1) stack[stack.length - 1].children.push(s.slice(i + 9, end));
      i = end + 3;
    } else if (s.startsWith('<?', i)) {
      const end = s.indexOf('?>', i + 2);
      if (end < 0) fail('Unterminated processing instruction', i);
      i = end + 2;
    } else if (s.startsWith('<!', i)) {
      // DOCTYPE or other declaration; may contain an internal subset in [...] and quoted strings
      let j = i + 2;
      let depth = 0;
      let quote = 0;
      for (; j < n; j++) {
        const c = s.charCodeAt(j);
        if (quote) {
          if (c === quote) quote = 0;
        } else if (c === 34 || c === 39) quote = c;
        else if (c === 91) depth++; // [
        else if (c === 93) depth--; // ]
        else if (c === 62 && depth <= 0) break; // >
      }
      if (j >= n) fail('Unterminated declaration', i);
      i = j + 1;
    } else if (s.charCodeAt(i + 1) === 47) {
      // closing tag </name>
      const end = s.indexOf('>', i + 2);
      if (end < 0) fail('Unterminated closing tag', i);
      const name = s.slice(i + 2, end).trim();
      // pop to the matching element; ignore stray closing tags
      let k = stack.length - 1;
      while (k > 0 && stack[k].name !== name) k--;
      if (k > 0) stack.length = k;
      i = end + 1;
    } else {
      // opening tag
      let j = i + 1;
      while (j < n && isNameChar(s.charCodeAt(j))) j++;
      const name = s.slice(i + 1, j);
      if (!name) fail("Invalid character after '<'", i);
      const el: XmlElement = { name, attrs: {}, children: [] };
      let selfClosing = false;
      for (;;) {
        while (j < n && isSpace(s.charCodeAt(j))) j++;
        if (j >= n) fail(`Unterminated tag <${name}>`, i);
        const c = s.charCodeAt(j);
        if (c === 62) { // >
          j++;
          break;
        }
        if (c === 47) { // />
          if (s.charCodeAt(j + 1) !== 62) fail(`Malformed tag <${name}>`, j);
          selfClosing = true;
          j += 2;
          break;
        }
        const aStart = j;
        while (j < n && isNameChar(s.charCodeAt(j))) j++;
        const aName = s.slice(aStart, j);
        if (!aName) fail(`Malformed attribute in <${name}>`, j);
        while (j < n && isSpace(s.charCodeAt(j))) j++;
        if (s.charCodeAt(j) !== 61) {
          // attribute without value (HTML style) – tolerate
          el.attrs[aName] = '';
          continue;
        }
        j++;
        while (j < n && isSpace(s.charCodeAt(j))) j++;
        const q = s.charCodeAt(j);
        if (q === 34 || q === 39) {
          const end = s.indexOf(q === 34 ? '"' : "'", j + 1);
          if (end < 0) fail(`Unterminated attribute value in <${name}>`, j);
          el.attrs[aName] = decodeEntities(s.slice(j + 1, end));
          j = end + 1;
        } else {
          // unquoted value – tolerate
          const vStart = j;
          while (j < n && !isSpace(s.charCodeAt(j)) && s.charCodeAt(j) !== 62) j++;
          el.attrs[aName] = decodeEntities(s.slice(vStart, j));
        }
      }
      stack[stack.length - 1].children.push(el);
      if (!selfClosing) stack.push(el);
      i = j;
    }
  }
  const root = doc.children.find((c): c is XmlElement => typeof c !== 'string');
  if (!root) throw new FormatError(formatName, 'No XML root element found');
  return root;
}

/** Element name without namespace prefix ('cml:atom' → 'atom'). */
export function localName(name: string): string {
  const k = name.indexOf(':');
  return k < 0 ? name : name.slice(k + 1);
}

/** Child elements (optionally filtered by local name, case-sensitive). */
export function childElements(el: XmlElement, name?: string): XmlElement[] {
  const out: XmlElement[] = [];
  for (const c of el.children) if (typeof c !== 'string' && (!name || localName(c.name) === name)) out.push(c);
  return out;
}

export function firstChild(el: XmlElement, name: string): XmlElement | undefined {
  for (const c of el.children) if (typeof c !== 'string' && localName(c.name) === name) return c;
  return undefined;
}

/** All descendant elements with the given local name (depth-first, document order). */
export function descendants(el: XmlElement, name: string): XmlElement[] {
  const out: XmlElement[] = [];
  // iterative pre-order traversal (children pushed in reverse so they pop in document order)
  const stack: XmlElement[] = childElements(el).reverse();
  while (stack.length) {
    const e = stack.pop()!;
    if (localName(e.name) === name) out.push(e);
    const kids = childElements(e);
    for (let k = kids.length - 1; k >= 0; k--) stack.push(kids[k]);
  }
  return out;
}

/** Concatenated text of an element and its descendants. */
export function textContent(el: XmlElement): string {
  let s = '';
  for (const c of el.children) s += typeof c === 'string' ? c : textContent(c);
  return s;
}

/** Attribute string builder: skips undefined/null values. */
export function attrs(list: Record<string, string | number | undefined | null | false>): string {
  let out = '';
  for (const [k, v] of Object.entries(list)) {
    if (v === undefined || v === null || v === false) continue;
    out += ` ${k}="${escapeXml(String(v))}"`;
  }
  return out;
}
