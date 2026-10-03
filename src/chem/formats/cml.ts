// Chemical Markup Language (CML): <molecule><atomArray>/<bondArray>, several molecules per file.
// Supports the attribute-per-atom style, the old array-attribute style (atomID="a1 a2" …) and
// CML1 <string builtin="…"> children. Coordinates are Å with y up (as in molfiles).
import { Atom, Mol } from '../mol';
import { BY_SYMBOL } from '../elements';
import { implicitH, explicitHNeighbors } from '../valence';
import { expandAbbreviations } from '../abbreviations';
import { parseXml, XmlElement, childElements, firstChild, descendants, localName, textContent, attrs, escapeText } from './xml';
import {
  FormatError, guard, normalizeImportedCoords, kekulizeFlagged, dropDefaultHCounts, perceiveStereo, normalizeElement,
  DEFAULT_BOND_LENGTH_ANGSTROM, fixed, applyAlias,
} from './common';

/** Reads every molecule in a CML document. */
export function readCML(xml: string): Mol[] {
  return guard('CML', () => {
    const root = parseXml(xml, 'CML');
    const all = localName(root.name) === 'molecule' ? [root, ...descendants(root, 'molecule')] : descendants(root, 'molecule');
    // keep "leaf" molecules that carry atoms (container molecules group sub-molecules)
    const mols = all.filter((m) => firstChild(m, 'atomArray') || childElements(m, 'atom').length);
    if (!mols.length) throw new FormatError('CML', 'no <molecule> with an <atomArray> found');
    return mols.map(readMolecule);
  });
}

/** Atom/bond properties may be attributes or CML1 builtin children. */
function prop(el: XmlElement, name: string): string | undefined {
  if (el.attrs[name] !== undefined) return el.attrs[name];
  for (const c of childElements(el)) if (c.attrs.builtin === name) return textContent(c).trim();
  return undefined;
}

interface AtomRec {
  id: string;
  el?: string;
  x2?: string; y2?: string; x3?: string; y3?: string; z3?: string;
  charge?: string; hcount?: string; isotope?: string; spin?: string; label?: string;
}

function readMolecule(m: XmlElement): Mol {
  const mol = new Mol();
  mol.name = m.attrs.title ?? (firstChild(m, 'name') ? textContent(firstChild(m, 'name')!).trim() : '') ?? '';
  const recs: AtomRec[] = [];
  const atomArray = firstChild(m, 'atomArray');
  const atomEls = atomArray ? childElements(atomArray, 'atom') : childElements(m, 'atom');
  if (atomEls.length) {
    for (const a of atomEls) {
      const lab = firstChild(a, 'label');
      recs.push({
        id: a.attrs.id ?? `a${recs.length + 1}`,
        el: prop(a, 'elementType'),
        x2: prop(a, 'x2'), y2: prop(a, 'y2'), x3: prop(a, 'x3'), y3: prop(a, 'y3'), z3: prop(a, 'z3'),
        charge: prop(a, 'formalCharge'), hcount: prop(a, 'hydrogenCount'),
        isotope: prop(a, 'isotopeNumber') ?? prop(a, 'isotope'), spin: prop(a, 'spinMultiplicity'),
        label: lab?.attrs.value ?? (lab ? textContent(lab).trim() : undefined),
      });
    }
  } else if (atomArray) {
    // array style: <atomArray atomID="a1 a2" elementType="C O" x2="…" …/>
    const split = (k: string) => (atomArray.attrs[k] ?? '').trim().split(/\s+/).filter(Boolean);
    const ids = split('atomID');
    const cols: Record<string, string[]> = {};
    for (const k of ['elementType', 'x2', 'y2', 'x3', 'y3', 'z3', 'formalCharge', 'hydrogenCount', 'isotopeNumber']) cols[k] = split(k);
    ids.forEach((id, i) =>
      recs.push({
        id, el: cols.elementType[i], x2: cols.x2[i], y2: cols.y2[i], x3: cols.x3[i], y3: cols.y3[i], z3: cols.z3[i],
        charge: cols.formalCharge[i], hcount: cols.hydrogenCount[i], isotope: cols.isotopeNumber[i],
      }),
    );
  }
  const index = new Map<string, number>();
  const hcounts: (number | undefined)[] = [];
  for (const r of recs) {
    const sym = r.el ?? 'C';
    const el = sym === 'D' || sym === 'T' ? 'H' : normalizeElement(sym);
    const has3 = r.x3 !== undefined && r.y3 !== undefined;
    const x = Number(has3 ? r.x3 : r.x2 ?? 0) || 0;
    const y = Number(has3 ? r.y3 : r.y2 ?? 0) || 0;
    const z = has3 ? Number(r.z3 ?? 0) || 0 : 0;
    const atom: Partial<Atom> & { el: string } = { el: el ?? (sym === 'R' ? 'R' : '*'), x, y, z, charge: parseInt(r.charge ?? '0', 10) || 0 };
    if (sym === 'D') atom.isotope = 2;
    if (sym === 'T') atom.isotope = 3;
    if (!el && sym !== 'R' && sym !== 'Du' && sym !== '*') atom.alias = sym;
    const iso = parseInt(r.isotope ?? '', 10);
    if (iso > 0) atom.isotope = iso;
    const spin = parseInt(r.spin ?? '', 10);
    if (spin === 2) atom.radical = 1;
    else if (spin === 3) atom.radical = 2;
    const i = mol.addAtom(atom);
    if (r.label && (atom.el === 'R' || atom.el === '*')) applyAlias(mol.atoms[i], r.label);
    index.set(r.id, i);
    const h = parseInt(r.hcount ?? '', 10);
    hcounts[i] = Number.isFinite(h) ? h : undefined;
  }
  // bonds
  const aromatic: boolean[] = [];
  const bonded = new Set<string>();
  const bondArray = firstChild(m, 'bondArray');
  const bondEls = bondArray ? childElements(bondArray, 'bond') : childElements(m, 'bond');
  const addBond = (refs: string[], order: string | undefined, stereo: string | undefined, stereoRefs?: string[]) => {
    let a = index.get(refs[0]), b = index.get(refs[1]);
    if (a === undefined || b === undefined || a === b) return;
    const key = a < b ? `${a},${b}` : `${b},${a}`;
    if (bonded.has(key)) return;
    bonded.add(key);
    const o = (order ?? '1').toUpperCase();
    let ord = 1;
    let aro = false;
    if (o === '2' || o === 'D') ord = 2;
    else if (o === '3' || o === 'T') ord = 3;
    else if (o === 'A' || o === '1.5') { ord = 1.5; aro = true; }
    else if (o === '0') ord = 0;
    let style: 'plain' | 'wedge' | 'hash' = 'plain';
    const st = (stereo ?? '').trim().toUpperCase();
    if (ord === 1 && (st === 'W' || st === 'H')) {
      style = st === 'W' ? 'wedge' : 'hash';
      // optional atomRefs2 on bondStereo gives the narrow end first
      const sa = stereoRefs ? index.get(stereoRefs[0]) : undefined;
      if (sa === b) [a, b] = [b, a];
    }
    const bi = mol.addBond(a, b, ord, style);
    aromatic[bi] = aro;
  };
  if (bondEls.length) {
    for (const b of bondEls) {
      const refs = (prop(b, 'atomRefs2') ?? prop(b, 'atomRefs') ?? '').trim().split(/\s+/);
      const bs = firstChild(b, 'bondStereo');
      const sRefs = bs?.attrs.atomRefs2?.trim().split(/\s+/);
      addBond(refs, prop(b, 'order'), bs ? textContent(bs) : prop(b, 'stereo'), sRefs);
    }
  } else if (bondArray) {
    const split = (k: string) => (bondArray.attrs[k] ?? '').trim().split(/\s+/).filter(Boolean);
    const r1 = split('atomRef1'), r2 = split('atomRef2'), ord = split('order');
    r1.forEach((r, i) => addBond([r, r2[i]], ord[i], undefined));
  }
  normalizeImportedCoords([mol]);
  // hydrogenCount is the total H count (implicit + explicit H atoms); known before kekulisation
  hcounts.forEach((h, i) => {
    if (h === undefined || mol.atoms[i].el === 'R' || mol.atoms[i].el === '*') return;
    mol.atoms[i].hCount = Math.max(0, h - explicitHNeighbors(mol, i));
  });
  kekulizeFlagged(mol, aromatic);
  dropDefaultHCounts(mol);
  perceiveStereo(mol);
  return mol;
}

/** Writes one molecule as a CML document (abbreviations expanded). */
export function writeCML(input: Mol): string {
  const mol = input.atoms.some((a) => a.abbrev) ? expandAbbreviations(input) : input;
  mol.invalidate();
  const threeD = mol.props.dim === '3D' || mol.atoms.some((a) => a.z);
  const L = DEFAULT_BOND_LENGTH_ANGSTROM;
  const out: string[] = ['<?xml version="1.0" encoding="UTF-8"?>', '<cml xmlns="http://www.xml-cml.org/schema">'];
  out.push(`  <molecule${attrs({ id: 'm1', title: mol.name || undefined })}>`);
  if (mol.atoms.length) {
    out.push('    <atomArray>');
    mol.atoms.forEach((a, i) => {
      const pseudo = !BY_SYMBOL.has(a.el);
      const coords = threeD
        ? { x3: fixed(a.x, 4), y3: fixed(-a.y, 4), z3: fixed(a.z ?? 0, 4) }
        : { x2: fixed(a.x * L, 4), y2: fixed(-a.y * L, 4) };
      const h = pseudo ? undefined : implicitH(mol, i) + explicitHNeighbors(mol, i);
      const at = attrs({
        id: `a${i + 1}`,
        elementType: pseudo ? 'R' : a.el,
        ...coords,
        formalCharge: a.charge || undefined,
        hydrogenCount: h,
        isotopeNumber: a.isotope,
        spinMultiplicity: a.radical ? a.radical + 1 : undefined,
      });
      if (pseudo && a.alias) out.push(`      <atom${at}><label${attrs({ value: a.alias })}/></atom>`);
      else out.push(`      <atom${at}/>`);
    });
    out.push('    </atomArray>');
  }
  if (mol.bonds.length) {
    out.push('    <bondArray>');
    mol.bonds.forEach((b, i) => {
      const order = b.order === 1.5 ? 'A' : b.order === 0 ? '0' : String(Math.min(3, Math.max(1, Math.round(b.order))));
      const at = attrs({ id: `b${i + 1}`, atomRefs2: `a${b.a + 1} a${b.b + 1}`, order });
      if (b.style === 'wedge' || b.style === 'hash') out.push(`      <bond${at}><bondStereo>${escapeText(b.style === 'wedge' ? 'W' : 'H')}</bondStereo></bond>`);
      else out.push(`      <bond${at}/>`);
    });
    out.push('    </bondArray>');
  }
  out.push('  </molecule>', '</cml>', '');
  return out.join('\n');
}
