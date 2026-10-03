// PubChem PUG REST / PUG View client (CORS-enabled, no key required).
// Docs: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest

const BASE = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug';
const VIEW = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view';
const AUTO = 'https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound';

export interface CompoundSummary {
  cid: number;
  title?: string;
  iupacName?: string;
  formula?: string;
  mw?: number;
  exactMass?: number;
  smiles?: string;
  connectivitySmiles?: string;
  inchi?: string;
  inchiKey?: string;
  xlogp?: number;
  tpsa?: number;
  hbd?: number;
  hba?: number;
  charge?: number;
  complexity?: number;
}

export interface Hazards {
  signal?: string;
  pictograms: { code: string; url: string; name: string }[];
  statements: string[];
  precautionary?: string;
  source?: string;
}

const PICTO_NAMES: Record<string, string> = {
  GHS01: 'Explosive', GHS02: 'Flammable', GHS03: 'Oxidizer', GHS04: 'Compressed gas', GHS05: 'Corrosive',
  GHS06: 'Acute toxicity', GHS07: 'Irritant / harmful', GHS08: 'Health hazard', GHS09: 'Environmental hazard',
};

const PROPS = 'Title,IUPACName,MolecularFormula,MolecularWeight,ExactMass,SMILES,ConnectivitySMILES,InChI,InChIKey,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,Charge,Complexity';

const cache = new Map<string, Promise<any>>();

async function getJSON(url: string, opts: RequestInit = {}): Promise<any> {
  const key = url + (opts.body ? '|' + String(opts.body) : '');
  if (!opts.method || opts.method === 'GET') {
    const c = cache.get(key);
    if (c) return c;
  }
  const p = (async () => {
    let lastErr: unknown;
    for (let attempt = 0; attempt < 3; attempt++) {
      const res = await fetch(url, opts).catch((e) => {
        lastErr = e;
        return null;
      });
      if (!res) {
        await new Promise((r) => setTimeout(r, 400 * (attempt + 1)));
        continue;
      }
      if (res.status === 404) throw new NotFoundError();
      if (res.status === 503 || res.status === 429) {
        await new Promise((r) => setTimeout(r, 800 * (attempt + 1)));
        continue;
      }
      if (!res.ok) {
        let msg = `PubChem error ${res.status}`;
        try {
          const j = await res.json();
          msg = j?.Fault?.Message ? `PubChem: ${j.Fault.Message}${j.Fault.Details ? ' — ' + j.Fault.Details.join(' ') : ''}` : msg;
        } catch {
          /* ignore */
        }
        throw new Error(msg);
      }
      const ct = res.headers.get('content-type') ?? '';
      return ct.includes('json') ? res.json() : res.text();
    }
    throw lastErr instanceof Error ? lastErr : new Error('PubChem is not reachable (offline?)');
  })();
  if (!opts.method || opts.method === 'GET') {
    cache.set(key, p);
    p.catch(() => cache.delete(key));
  }
  return p;
}

export class NotFoundError extends Error {
  constructor() {
    super('No matching compound found in PubChem');
  }
}

function toSummary(p: any): CompoundSummary {
  return {
    cid: p.CID,
    title: p.Title,
    iupacName: p.IUPACName,
    formula: p.MolecularFormula,
    mw: p.MolecularWeight !== undefined ? +p.MolecularWeight : undefined,
    exactMass: p.ExactMass !== undefined ? +p.ExactMass : undefined,
    smiles: p.SMILES ?? p.IsomericSMILES ?? p.CanonicalSMILES,
    connectivitySmiles: p.ConnectivitySMILES ?? p.CanonicalSMILES,
    inchi: p.InChI,
    inchiKey: p.InChIKey,
    xlogp: p.XLogP,
    tpsa: p.TPSA,
    hbd: p.HBondDonorCount,
    hba: p.HBondAcceptorCount,
    charge: p.Charge,
    complexity: p.Complexity,
  };
}

/** Name suggestions for a partial query. */
export async function autocomplete(q: string, limit = 8): Promise<string[]> {
  if (q.trim().length < 2) return [];
  const j = await getJSON(`${AUTO}/${encodeURIComponent(q.trim())}/json?limit=${limit}`);
  return j?.dictionary_terms?.compound ?? [];
}

/** Search by name, CAS number, InChIKey, CID or formula-like text. Returns up to `limit` compounds. */
export async function search(query: string, limit = 10): Promise<CompoundSummary[]> {
  const q = query.trim();
  if (!q) return [];
  let cids: number[] = [];
  if (/^\d+$/.test(q)) cids = [+q];
  else if (/^[A-Z]{14}-[A-Z]{10}-[A-Z]$/.test(q)) {
    const j = await getJSON(`${BASE}/compound/inchikey/${q}/cids/JSON`);
    cids = j.IdentifierList?.CID ?? [];
  } else {
    // name (also matches CAS numbers and many synonyms)
    try {
      const j = await getJSON(`${BASE}/compound/name/${encodeURIComponent(q)}/cids/JSON`);
      cids = j.IdentifierList?.CID ?? [];
    } catch (e) {
      if (!(e instanceof NotFoundError)) throw e;
    }
    if (!cids.length && /^[A-Z][A-Za-z0-9()]*$/.test(q)) {
      // try molecular formula (asynchronous fast search)
      try {
        const j = await getJSON(`${BASE}/compound/fastformula/${encodeURIComponent(q)}/cids/JSON?MaxRecords=${limit}`);
        cids = j.IdentifierList?.CID ?? [];
      } catch {
        /* ignore */
      }
    }
  }
  if (!cids.length) return [];
  return properties(cids.slice(0, limit));
}

export async function properties(cids: number[]): Promise<CompoundSummary[]> {
  if (!cids.length) return [];
  const j = await getJSON(`${BASE}/compound/cid/${cids.join(',')}/property/${PROPS}/JSON`);
  return (j.PropertyTable?.Properties ?? []).map(toSummary);
}

/** Exact-structure lookup by SMILES (POST so that any SMILES is accepted). */
export async function identifyBySmiles(smiles: string): Promise<CompoundSummary | null> {
  try {
    const j = await getJSON(`${BASE}/compound/smiles/cids/JSON`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/x-www-form-urlencoded' },
      body: 'smiles=' + encodeURIComponent(smiles),
    });
    const cid = j.IdentifierList?.CID?.[0];
    if (!cid) return null;
    const [s] = await properties([cid]);
    return s ?? null;
  } catch (e) {
    if (e instanceof NotFoundError) return null;
    throw e;
  }
}

/** 2D structure (SDF) of a compound. */
export async function sdf2d(cid: number): Promise<string> {
  return getJSON(`${BASE}/compound/cid/${cid}/SDF?record_type=2d`);
}

/** 3D conformer (SDF) of a compound, if PubChem has one. */
export async function sdf3d(cid: number): Promise<string> {
  return getJSON(`${BASE}/compound/cid/${cid}/SDF?record_type=3d`);
}

export async function synonyms(cid: number): Promise<string[]> {
  try {
    const j = await getJSON(`${BASE}/compound/cid/${cid}/synonyms/JSON`);
    return j.InformationList?.Information?.[0]?.Synonym ?? [];
  } catch {
    return [];
  }
}

/** CAS registry numbers among the synonyms (validated with the CAS check digit). */
export function casNumbers(syns: string[]): string[] {
  const out: string[] = [];
  for (const s of syns) {
    const m = /^(\d{2,7})-(\d{2})-(\d)$/.exec(s);
    if (!m) continue;
    const digits = (m[1] + m[2]).split('').reverse();
    const sum = digits.reduce((acc, d, i) => acc + +d * (i + 1), 0);
    if (sum % 10 === +m[3] && !out.includes(s)) out.push(s);
  }
  return out;
}

/** GHS hazard classification (aggregated by PubChem from ECHA C&L and other sources). */
export async function hazards(cid: number): Promise<Hazards | null> {
  let j: any;
  try {
    j = await getJSON(`${VIEW}/data/compound/${cid}/JSON?heading=GHS+Classification`);
  } catch (e) {
    if (e instanceof NotFoundError) return null;
    throw e;
  }
  const out: Hazards = { pictograms: [], statements: [] };
  const walk = (sec: any) => {
    for (const s of sec.Section ?? []) walk(s);
    for (const inf of sec.Information ?? []) {
      const name: string = inf.Name ?? '';
      const swm = inf.Value?.StringWithMarkup ?? [];
      if (name.startsWith('Pictogram')) {
        for (const x of swm)
          for (const m of x.Markup ?? []) {
            const url: string = m.URL ?? '';
            const code = /GHS\d\d/.exec(url)?.[0];
            if (code && !out.pictograms.some((p) => p.code === code)) out.pictograms.push({ code, url, name: PICTO_NAMES[code] ?? code });
          }
      } else if (name === 'Signal' && !out.signal) {
        out.signal = swm[0]?.String;
      } else if (name === 'GHS Hazard Statements' && !out.statements.length) {
        out.statements = swm.map((x: any) => x.String).filter(Boolean);
      } else if (name === 'Precautionary Statement Codes' && !out.precautionary) {
        out.precautionary = swm[0]?.String;
      }
    }
  };
  walk(j.Record ?? {});
  out.source = 'PubChem (aggregated GHS notifications)';
  if (!out.pictograms.length && !out.statements.length && !out.signal) return null;
  return out;
}

/** Structurally similar compounds (2D Tanimoto ≥ threshold). */
export async function similar(smiles: string, threshold = 90, max = 12): Promise<CompoundSummary[]> {
  try {
    const j = await getJSON(`${BASE}/compound/fastsimilarity_2d/smiles/cids/JSON?Threshold=${threshold}&MaxRecords=${max}`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/x-www-form-urlencoded' },
      body: 'smiles=' + encodeURIComponent(smiles),
    });
    const cids: number[] = j.IdentifierList?.CID ?? [];
    return properties(cids.slice(0, max));
  } catch (e) {
    if (e instanceof NotFoundError) return [];
    throw e;
  }
}

export const imageUrl = (cid: number, size = 200) => `${BASE}/compound/cid/${cid}/PNG?image_size=${size}x${size}`;
export const pageUrl = (cid: number) => `https://pubchem.ncbi.nlm.nih.gov/compound/${cid}`;
