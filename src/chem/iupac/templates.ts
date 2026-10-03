// Ring-system templates with fixed IUPAC numbering.
//
// Each template is described by its rings, written as cycles of locant labels, plus the
// heteroatom positions. The skeleton graph (atoms + ring bonds) is matched against a molecule's
// ring system by graph isomorphism; every isomorphism yields one admissible numbering.
// Hydrogenation state is NOT part of the template: names are built for the mancude parent and
// hydro prefixes / indicated hydrogen are derived from the actual structure.

export interface TemplateDef {
  /** Parent name of the mancude ring system (without indicated hydrogen). */
  name: string;
  /** Rings as comma-separated locant cycles. */
  rings: string[];
  /** Heteroatoms: locant → element symbol. */
  het?: Record<string, string>;
  /** Saturated template (no hydro prefixes, unsaturation not supported), e.g. adamantane. */
  saturated?: boolean;
}

const NAPH = ['1,2,3,4,4a,8a', '4a,5,6,7,8,8a'];
const INDENE = ['1,2,3,3a,7a', '3a,4,5,6,7,7a'];
const ANTH = ['1,2,3,4,4a,9a', '4a,10,10a,8a,9,9a', '10a,5,6,7,8,8a'];
const PHEN = ['1,2,3,4,4a,10a', '4a,4b,8a,9,10,10a', '4b,5,6,7,8,8a'];
const FLUO = ['1,2,3,4,4a,9a', '4a,4b,8a,9,9a', '4b,5,6,7,8,8a'];
const PHENAZ = ['1,2,3,4,4a,10a', '4a,5,5a,9a,10,10a', '5a,6,7,8,9,9a'];
const PHENTRL = ['1,2,3,4,4a,10b', '4a,5,6,6a,10a,10b', '6a,7,8,9,10,10a'];
const DBF = ['1,2,3,4,4a,9b', '4a,5,5a,9a,9b', '5a,6,7,8,9,9a'];
const INDOLIZ = ['1,2,3,4,8a', '4,5,6,7,8,8a'];

export const TEMPLATES: TemplateDef[] = [
  // ── 6-6 ──
  { name: 'naphthalene', rings: NAPH },
  { name: 'quinoline', rings: NAPH, het: { '1': 'N' } },
  { name: 'isoquinoline', rings: NAPH, het: { '2': 'N' } },
  { name: 'quinazoline', rings: NAPH, het: { '1': 'N', '3': 'N' } },
  { name: 'quinoxaline', rings: NAPH, het: { '1': 'N', '4': 'N' } },
  { name: 'cinnoline', rings: NAPH, het: { '1': 'N', '2': 'N' } },
  { name: 'phthalazine', rings: NAPH, het: { '2': 'N', '3': 'N' } },
  { name: '1,5-naphthyridine', rings: NAPH, het: { '1': 'N', '5': 'N' } },
  { name: '1,6-naphthyridine', rings: NAPH, het: { '1': 'N', '6': 'N' } },
  { name: '1,7-naphthyridine', rings: NAPH, het: { '1': 'N', '7': 'N' } },
  { name: '1,8-naphthyridine', rings: NAPH, het: { '1': 'N', '8': 'N' } },
  { name: '2,6-naphthyridine', rings: NAPH, het: { '2': 'N', '6': 'N' } },
  { name: '2,7-naphthyridine', rings: NAPH, het: { '2': 'N', '7': 'N' } },
  { name: 'pteridine', rings: NAPH, het: { '1': 'N', '3': 'N', '5': 'N', '8': 'N' } },
  { name: 'chromene', rings: NAPH, het: { '1': 'O' } },
  { name: 'isochromene', rings: NAPH, het: { '2': 'O' } },
  { name: 'thiochromene', rings: NAPH, het: { '1': 'S' } },
  { name: '1,4-benzodioxine', rings: NAPH, het: { '1': 'O', '4': 'O' } },
  { name: '1,4-benzoxazine', rings: NAPH, het: { '1': 'O', '4': 'N' } },
  { name: '1,4-benzothiazine', rings: NAPH, het: { '1': 'S', '4': 'N' } },
  { name: 'quinolizine', rings: ['1,2,3,4,5,9a', '5,6,7,8,9,9a'], het: { '5': 'N' } },
  // ── 6-5 ──
  { name: 'indene', rings: INDENE },
  { name: 'indole', rings: INDENE, het: { '1': 'N' } },
  { name: 'isoindole', rings: INDENE, het: { '2': 'N' } },
  { name: 'indazole', rings: INDENE, het: { '1': 'N', '2': 'N' } },
  { name: 'benzimidazole', rings: INDENE, het: { '1': 'N', '3': 'N' } },
  { name: 'benzotriazole', rings: INDENE, het: { '1': 'N', '2': 'N', '3': 'N' } },
  { name: '1-benzofuran', rings: INDENE, het: { '1': 'O' } },
  { name: '2-benzofuran', rings: INDENE, het: { '2': 'O' } },
  { name: '1-benzothiophene', rings: INDENE, het: { '1': 'S' } },
  { name: '2-benzothiophene', rings: INDENE, het: { '2': 'S' } },
  { name: '1,3-benzoxazole', rings: INDENE, het: { '1': 'O', '3': 'N' } },
  { name: '1,3-benzothiazole', rings: INDENE, het: { '1': 'S', '3': 'N' } },
  { name: '1,2-benzoxazole', rings: INDENE, het: { '1': 'O', '2': 'N' } },
  { name: '1,2-benzothiazole', rings: INDENE, het: { '1': 'S', '2': 'N' } },
  { name: '1,3-benzodioxole', rings: INDENE, het: { '1': 'O', '3': 'O' } },
  { name: '2,1,3-benzothiadiazole', rings: INDENE, het: { '1': 'N', '2': 'S', '3': 'N' } },
  { name: 'pyrrolo[2,3-b]pyridine', rings: INDENE, het: { '1': 'N', '7': 'N' } },
  { name: 'thieno[3,2-c]pyridine', rings: INDENE, het: { '1': 'S', '5': 'N' } },
  { name: 'thieno[2,3-c]pyridine', rings: INDENE, het: { '1': 'S', '6': 'N' } },
  { name: 'pyrazolo[3,4-d]pyrimidine', rings: INDENE, het: { '1': 'N', '2': 'N', '5': 'N', '7': 'N' } },
  { name: 'indolizine', rings: INDOLIZ, het: { '4': 'N' } },
  { name: 'imidazo[1,2-a]pyridine', rings: INDOLIZ, het: { '1': 'N', '4': 'N' } },
  { name: 'imidazo[1,2-a]pyrimidine', rings: INDOLIZ, het: { '1': 'N', '4': 'N', '8': 'N' } },
  { name: '[1,2,4]triazolo[4,3-a]pyridine', rings: INDOLIZ, het: { '1': 'N', '2': 'N', '4': 'N' } },
  { name: 'purine', rings: ['1,2,3,4,5,6', '4,5,7,8,9'], het: { '1': 'N', '3': 'N', '7': 'N', '9': 'N' } },
  // ── 5-5, 7-5, 7-6, 7-6 hetero ──
  { name: 'pentalene', rings: ['1,2,3,3a,6a', '3a,4,5,6,6a'] },
  { name: 'azulene', rings: ['1,2,3,3a,8a', '3a,4,5,6,7,8,8a'] },
  { name: 'benzo[7]annulene', rings: ['1,2,3,4,4a,9a', '4a,5,6,7,8,9,9a'] },
  { name: '1,4-benzodiazepine', rings: ['1,2,3,4,5,5a,9a', '5a,6,7,8,9,9a'], het: { '1': 'N', '4': 'N' } },
  { name: '1,5-benzodiazepine', rings: ['1,2,3,4,5,5a,9a', '5a,6,7,8,9,9a'], het: { '1': 'N', '5': 'N' } },
  // ── tricyclic ──
  { name: 'anthracene', rings: ANTH },
  { name: 'acridine', rings: ANTH, het: { '10': 'N' } },
  { name: 'xanthene', rings: ANTH, het: { '10': 'O' } },
  { name: 'thioxanthene', rings: ANTH, het: { '10': 'S' } },
  { name: 'phenanthrene', rings: PHEN },
  { name: 'fluorene', rings: FLUO },
  { name: 'carbazole', rings: FLUO, het: { '9': 'N' } },
  { name: 'dibenzofuran', rings: DBF, het: { '5': 'O' } },
  { name: 'dibenzothiophene', rings: DBF, het: { '5': 'S' } },
  { name: 'phenazine', rings: PHENAZ, het: { '5': 'N', '10': 'N' } },
  { name: 'phenothiazine', rings: PHENAZ, het: { '5': 'S', '10': 'N' } },
  { name: 'phenoxazine', rings: PHENAZ, het: { '5': 'O', '10': 'N' } },
  { name: 'phenoxathiine', rings: PHENAZ, het: { '5': 'O', '10': 'S' } },
  { name: 'thianthrene', rings: PHENAZ, het: { '5': 'S', '10': 'S' } },
  { name: 'dibenzo[b,e][1,4]dioxine', rings: PHENAZ, het: { '5': 'O', '10': 'O' } },
  { name: '1,10-phenanthroline', rings: PHENTRL, het: { '1': 'N', '10': 'N' } },
  { name: 'phenanthridine', rings: PHENTRL, het: { '5': 'N' } },
  { name: 'pyrido[3,4-b]indole', rings: ['1,2,3,4,4a,9a', '4a,4b,8a,9,9a', '4b,5,6,7,8,8a'], het: { '2': 'N', '9': 'N' } },
  // ── tetracyclic (steroid skeleton; PubChem uses steroid numbering for this system) ──
  { name: 'cyclopenta[a]phenanthrene', rings: ['1,2,3,4,5,10', '5,6,7,8,9,10', '8,9,11,12,13,14', '13,14,15,16,17'] },
  { name: 'chrysene', rings: ['1,2,3,4,4a,12a', '4a,4b,10b,11,12,12a', '4b,5,6,6a,10a,10b', '6a,7,8,9,10,10a'] },
  { name: 'tetracene', rings: ['1,2,3,4,4a,12a', '4a,5,5a,11a,12,12a', '5a,6,6a,10a,11,11a', '6a,7,8,9,10,10a'] },
  { name: 'triphenylene', rings: ['1,2,3,4,4a,12b', '4a,4b,8a,8b,12a,12b', '4b,5,6,7,8,8a', '8b,9,10,11,12,12a'] },
  // ── pyrene ──
  { name: 'pyrene', rings: ['1,2,3,3a,10b,10a', '3a,4,5,5a,10c,10b', '5a,6,7,8,8a,10c', '8a,9,10,10a,10b,10c'] },
  // ── saturated polycycles ──
  { name: 'adamantane', rings: ['1,2,3,4,5,6,7,8', '1,9,5,4,3,2', '3,10,7,6,5,4'], saturated: true },
];

export interface CompiledTemplate {
  def: TemplateDef;
  labels: string[];
  els: string[];
  adj: number[][];
  nBonds: number;
}

let compiled: CompiledTemplate[] | null = null;

/** Locant label → sortable number ("4a" → 4.01, "10b" → 10.02). */
export function locantValue(label: string): number {
  const m = /^(\d+)([a-z]?)/.exec(label);
  if (!m) return 1e6;
  return +m[1] + (m[2] ? (m[2].charCodeAt(0) - 96) / 100 : 0);
}

export function compiledTemplates(): CompiledTemplate[] {
  if (compiled) return compiled;
  compiled = TEMPLATES.map((def) => {
    const labels: string[] = [];
    const idx = new Map<string, number>();
    const edges = new Set<string>();
    const adj: number[][] = [];
    const get = (l: string) => {
      let i = idx.get(l);
      if (i === undefined) {
        i = labels.length;
        idx.set(l, i);
        labels.push(l);
        adj.push([]);
      }
      return i;
    };
    for (const r of def.rings) {
      const cyc = r.split(',').map((s) => get(s.trim()));
      for (let k = 0; k < cyc.length; k++) {
        const a = cyc[k], b = cyc[(k + 1) % cyc.length];
        const key = a < b ? a + '-' + b : b + '-' + a;
        if (edges.has(key)) continue;
        edges.add(key);
        adj[a].push(b);
        adj[b].push(a);
      }
    }
    const els = labels.map((l) => def.het?.[l] ?? 'C');
    return { def, labels, els, adj, nBonds: edges.size };
  });
  return compiled;
}
