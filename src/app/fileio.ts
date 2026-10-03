// Open / save / export / share-link helpers.
import { ChemDoc } from '../doc/types';
import { docToJSON, docFromJSON, docToMol, createDoc, insertMol } from '../doc/document';
import { buildScene } from '../render/scene';
import { primsToSVG, primsToCanvas } from '../render/draw';
import { downloadBlob, toast } from './dom';
import { formats, chemMol, canonicalSmiles } from './chem';
import { expandAbbreviations } from '../chem/abbreviations';
import { perceiveStereo2D } from '../chem/stereo2d';
import { Mol } from '../chem/mol';

export const OPEN_ACCEPT = '.cwj,.chemwrite,.json,.mol,.sdf,.sd,.rxn,.cdxml,.cml,.xyz,.smi,.smiles,.txt';

export function pickFile(accept = OPEN_ACCEPT): Promise<File | null> {
  return new Promise((resolve) => {
    const input = document.createElement('input');
    input.type = 'file';
    input.accept = accept;
    input.onchange = () => resolve(input.files?.[0] ?? null);
    input.addEventListener('cancel', () => resolve(null));
    input.click();
  });
}

// ───────────── native save with File System Access API when available ─────────────

let fileHandle: any = null;

export function resetFileHandle(): void {
  fileHandle = null;
}

export function setFileHandle(h: any): void {
  fileHandle = h;
}

export async function saveNative(doc: ChemDoc, saveAs = false): Promise<string | null> {
  const json = docToJSON(doc, true);
  const name = safeName(doc.meta.title || 'structure') + '.cwj';
  const w = window as any;
  if (w.showSaveFilePicker) {
    try {
      if (!fileHandle || saveAs) {
        fileHandle = await w.showSaveFilePicker({
          suggestedName: name,
          types: [{ description: 'ChemWrite document', accept: { 'application/json': ['.cwj'] } }],
        });
      }
      const ws = await fileHandle.createWritable();
      await ws.write(json);
      await ws.close();
      return fileHandle.name ?? name;
    } catch (e) {
      if ((e as Error).name === 'AbortError') return null;
      // fall through to download
    }
  }
  downloadBlob(new Blob([json], { type: 'application/json' }), name);
  return name;
}

export function safeName(s: string): string {
  return s.replace(/[^\w\-. ]+/g, '_').trim().slice(0, 80) || 'structure';
}

// ───────────── export ─────────────

export type ExportFormat = 'svg' | 'png' | 'mol' | 'sdf' | 'smiles' | 'cdxml' | 'rxn' | 'cwj' | 'cml';

export interface ExportOptions {
  /** Restrict to these atom/object ids (selection); null = everything */
  subDoc?: ChemDoc | null;
  pngScale?: number;
  transparent?: boolean;
}

/** Renders a document as an SVG string (always light/black-on-white palette). */
export function docSVG(doc: ChemDoc, transparent = true): string {
  const scene = buildScene(doc, { ink: '#000000', showErrors: false });
  const b = scene.bounds ?? { x1: 0, y1: 0, x2: 1, y2: 1 };
  return primsToSVG(scene.prims, b, doc.style.bondLengthPt, { padding: 0.4, background: transparent ? null : '#ffffff', title: doc.meta.title });
}

export function docPNG(doc: ChemDoc, scale = 2, transparent = false): Promise<Blob> {
  const scene = buildScene(doc, { ink: '#000000', showErrors: false });
  const b = scene.bounds ?? { x1: 0, y1: 0, x2: 1, y2: 1 };
  // 96 dpi screen reference: 1pt = 1.333 px
  const px = doc.style.bondLengthPt * 1.3333 * scale;
  const canvas = primsToCanvas(scene.prims, b, px, { padding: 0.4, background: transparent ? null : '#ffffff', maxSize: 8000 });
  return new Promise((resolve, reject) => canvas.toBlob((bl) => (bl ? resolve(bl) : reject(new Error('PNG export failed'))), 'image/png'));
}

/** Mols of all fragments in a document (abbreviations expanded, stereo perceived from wedges). */
export function docMols(doc: ChemDoc): Mol[] {
  const { mol } = docToMol(doc);
  return mol.components().map((c) => {
    const sub = mol.subset(c).mol;
    const ex = expandAbbreviations(sub);
    perceiveStereo2D(ex);
    return ex;
  });
}

export async function exportDoc(doc: ChemDoc, fmt: ExportFormat, opts: ExportOptions = {}): Promise<void> {
  const d = opts.subDoc ?? doc;
  const base = safeName(doc.meta.title || 'structure');
  switch (fmt) {
    case 'svg':
      downloadBlob(new Blob([docSVG(d, opts.transparent !== false)], { type: 'image/svg+xml' }), base + '.svg');
      break;
    case 'png':
      downloadBlob(await docPNG(d, opts.pngScale ?? 3, !!opts.transparent), base + '.png');
      break;
    case 'mol': {
      const { mol } = docToMol(d);
      downloadBlob(new Blob([formats.writeMolfile(mol, { title: doc.meta.title })], { type: 'chemical/x-mdl-molfile' }), base + '.mol');
      break;
    }
    case 'sdf':
      downloadBlob(new Blob([formats.writeSDF(docMols(d))], { type: 'chemical/x-mdl-sdfile' }), base + '.sdf');
      break;
    case 'cml': {
      const { mol } = docToMol(d);
      downloadBlob(new Blob([formats.writeCML(mol)], { type: 'chemical/x-cml' }), base + '.cml');
      break;
    }
    case 'smiles': {
      const smi = docMols(d).map((m) => canonicalSmiles(m)).join('\n');
      downloadBlob(new Blob([smi + '\n'], { type: 'chemical/x-daylight-smiles' }), base + '.smi');
      break;
    }
    case 'cdxml':
      downloadBlob(new Blob([formats.writeCDXML(d)], { type: 'chemical/x-cdxml' }), base + '.cdxml');
      break;
    case 'rxn': {
      const r = splitReaction(d);
      if (!r) {
        toast('Draw a reaction arrow with reactants on the left and products on the right', 'error');
        return;
      }
      downloadBlob(new Blob([formats.writeRxn(r)], { type: 'chemical/x-mdl-rxnfile' }), base + '.rxn');
      break;
    }
    case 'cwj':
      downloadBlob(new Blob([docToJSON(d, true)], { type: 'application/json' }), base + '.cwj');
      break;
  }
}

/** Splits a document into reactants / products around its first reaction arrow. */
export function splitReaction(doc: ChemDoc): { reactants: Mol[]; products: Mol[]; agents: Mol[] } | null {
  const arrow = [...doc.arrows.values()].find((a) => a.kind !== 'line' && a.kind !== 'resonance' && a.kind !== 'retro');
  if (!arrow) return null;
  const ux = arrow.x2 - arrow.x1, uy = arrow.y2 - arrow.y1;
  const L = Math.hypot(ux, uy) || 1;
  const mols = docMols(doc);
  const reactants: Mol[] = [], products: Mol[] = [], agents: Mol[] = [];
  for (const m of mols) {
    const c = m.atoms.reduce((s, a) => ({ x: s.x + a.x / m.atoms.length, y: s.y + a.y / m.atoms.length }), { x: 0, y: 0 });
    const t = ((c.x - arrow.x1) * ux + (c.y - arrow.y1) * uy) / L;
    if (t < 0) reactants.push(m);
    else if (t > L) products.push(m);
    else agents.push(m);
  }
  if (!reactants.length && !products.length) return null;
  return { reactants, products, agents };
}

/** Reaction SMILES of a drawn scheme (reactants>agents>products). */
export function reactionSmiles(doc: ChemDoc): string | null {
  const r = splitReaction(doc);
  if (!r) return null;
  const s = (ms: Mol[]) => ms.map((m) => canonicalSmiles(m)).join('.');
  return `${s(r.reactants)}>${s(r.agents)}>${s(r.products)}`;
}

// ───────────── share links ─────────────

async function deflate(text: string): Promise<Uint8Array> {
  const cs = new (window as any).CompressionStream('deflate-raw');
  const stream = new Blob([text]).stream().pipeThrough(cs);
  return new Uint8Array(await new Response(stream).arrayBuffer());
}

async function inflate(bytes: Uint8Array): Promise<string> {
  const ds = new (window as any).DecompressionStream('deflate-raw');
  const stream = new Blob([bytes as BlobPart]).stream().pipeThrough(ds);
  return new Response(stream).text();
}

const b64url = (bytes: Uint8Array) => {
  let s = '';
  for (let i = 0; i < bytes.length; i++) s += String.fromCharCode(bytes[i]);
  return btoa(s).replace(/\+/g, '-').replace(/\//g, '_').replace(/=+$/, '');
};
const fromB64url = (s: string) => {
  const bin = atob(s.replace(/-/g, '+').replace(/_/g, '/'));
  const out = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
  return out;
};

export async function shareLink(doc: ChemDoc): Promise<string> {
  const json = docToJSON(doc);
  let payload: string;
  if ((window as any).CompressionStream) payload = 'z' + b64url(await deflate(json));
  else payload = 'j' + b64url(new TextEncoder().encode(json));
  const url = new URL(location.href);
  url.hash = 'doc=' + payload;
  return url.toString();
}

export async function docFromHash(hash: string): Promise<ChemDoc | null> {
  const m = /doc=([A-Za-z0-9_\-]+)/.exec(hash);
  if (!m) {
    const s = /smiles=([^&]+)/.exec(hash);
    if (s) {
      const { smilesToMol } = await import('./chem');
      const doc = createDoc();
      insertMol(doc, smilesToMol(decodeURIComponent(s[1])));
      return doc;
    }
    return null;
  }
  const p = m[1];
  const bytes = fromB64url(p.slice(1));
  const json = p[0] === 'z' ? await inflate(bytes) : new TextDecoder().decode(bytes);
  return docFromJSON(json);
}

export { chemMol };
