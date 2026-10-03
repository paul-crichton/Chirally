// MDL reaction files (RXN V2000 and V3000).
// All components share one coordinate frame: they are scaled/flipped together on import so their
// relative placement is preserved.
import { Mol } from '../mol';
import { FormatError, guard, splitLines, normalizeImportedCoords } from './common';
import { parseMolBlock, collectV30, parseV3000Ctab, finishCtab, writeMolfile, v3000CtabFor, RawCtab } from './molfile';

export interface Reaction {
  reactants: Mol[];
  products: Mol[];
  agents: Mol[];
}

const int = (s: string | undefined) => (s !== undefined && /^\s*\d+\s*$/.test(s) ? parseInt(s, 10) : NaN);

/** Reads an RXN file (V2000 or V3000). */
export function readRxn(text: string): Reaction {
  return guard('RXN', () => {
    const lines = splitLines(text);
    let start = lines.findIndex((l) => l.trim() !== '');
    if (start < 0 || !lines[start].startsWith('$RXN')) throw new FormatError('RXN', 'missing "$RXN" header line');
    const groups = /V3000/.test(lines[start]) ? readV3000(lines, start) : readV2000(lines, start);
    const all = [...groups.reactants, ...groups.products, ...groups.agents];
    normalizeImportedCoords(all.map((r) => r.mol));
    const fin = (list: RawCtab[]) => list.map((r) => finishCtab(r));
    return { reactants: fin(groups.reactants), products: fin(groups.products), agents: fin(groups.agents) };
  });
}

interface RawGroups {
  reactants: RawCtab[];
  products: RawCtab[];
  agents: RawCtab[];
}

function readV2000(lines: string[], start: number): RawGroups {
  const counts = lines[start + 4];
  if (counts === undefined) throw new FormatError('RXN', 'file truncated before the counts line');
  let r = int(counts.slice(0, 3)), p = int(counts.slice(3, 6)), a = int(counts.slice(6, 9));
  if (!Number.isFinite(r) || !Number.isFinite(p)) {
    const t = counts.trim().split(/\s+/);
    r = int(t[0]);
    p = int(t[1]);
    a = int(t[2]);
  }
  if (!Number.isFinite(r) || !Number.isFinite(p)) throw new FormatError('RXN', `invalid counts line "${counts}"`);
  if (!Number.isFinite(a)) a = 0;
  const mols: RawCtab[] = [];
  for (let k = start + 5; k < lines.length; k++) {
    if (!lines[k].startsWith('$MOL')) continue;
    const { raw, end } = parseMolBlock(lines, k + 1, 'RXN');
    mols.push(raw);
    k = end - 1;
  }
  if (mols.length < r + p) throw new FormatError('RXN', `expected ${r + p + a} molecules, found ${mols.length}`);
  // V2000 order: reactants, products, then agents
  return { reactants: mols.slice(0, r), products: mols.slice(r, r + p), agents: mols.slice(r + p, r + p + a) };
}

function readV3000(lines: string[], start: number): RawGroups {
  const { v30 } = collectV30(lines, start + 1);
  const out: RawGroups = { reactants: [], products: [], agents: [] };
  let section: RawCtab[] | null = null;
  for (let k = 0; k < v30.length; k++) {
    const l = v30[k];
    if (l.startsWith('BEGIN REACTANT')) section = out.reactants;
    else if (l.startsWith('BEGIN PRODUCT')) section = out.products;
    else if (l.startsWith('BEGIN AGENT')) section = out.agents;
    else if (/^END (REACTANT|PRODUCT|AGENT)/.test(l)) section = null;
    else if (l.startsWith('BEGIN CTAB')) {
      const { raw, next } = parseV3000Ctab(v30, k, 'RXN');
      (section ?? out.agents).push(raw);
      k = next - 1;
    }
  }
  if (!out.reactants.length && !out.products.length) throw new FormatError('RXN', 'no reactant or product CTAB blocks found');
  return out;
}

/** Writes an RXN file (V2000; V3000 when any component exceeds 999 atoms/bonds). */
export function writeRxn(r: { reactants: Mol[]; products: Mol[]; agents?: Mol[] }): string {
  const agents = r.agents ?? [];
  const all = [...r.reactants, ...r.products, ...agents];
  const big = all.some((m) => m.atoms.length > 999 || m.bonds.length > 999);
  const d = new Date();
  const p2 = (v: number) => String(v).padStart(2, '0');
  const stamp = p2(d.getMonth() + 1) + p2(d.getDate()) + d.getFullYear() + p2(d.getHours()) + p2(d.getMinutes());
  if (big) {
    const out = ['$RXN V3000', '', '      ChemWrite' + stamp, '', `M  V30 COUNTS ${r.reactants.length} ${r.products.length}${agents.length ? ' ' + agents.length : ''}`];
    const block = (name: string, list: Mol[]) => {
      if (!list.length) return;
      out.push(`M  V30 BEGIN ${name}`);
      for (const m of list) out.push(...v3000CtabFor(m));
      out.push(`M  V30 END ${name}`);
    };
    block('REACTANT', r.reactants);
    block('PRODUCT', r.products);
    block('AGENT', agents);
    out.push('M  END', '');
    return out.join('\n');
  }
  const pad3 = (v: number) => String(v).padStart(3);
  let s = '$RXN\n\n      ChemWrite' + stamp + '\n\n';
  s += pad3(r.reactants.length) + pad3(r.products.length) + (agents.length ? pad3(agents.length) : '') + '\n';
  for (const m of all) s += '$MOL\n' + writeMolfile(m, { title: m.name });
  return s;
}
