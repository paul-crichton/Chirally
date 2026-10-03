// Web worker for 3D structure generation / force-field optimisation (keeps the UI responsive).
import { MolData, toData, fromData } from './moldata';
import { generate3D, optimizeGeometry, conformerSearch } from '../chem/3d';

export type WorkerRequest =
  | { id: number; type: 'generate'; mol: MolData; conformers: number; seed: number }
  | { id: number; type: 'optimize'; mol: MolData }
  | { id: number; type: 'conformers'; mol: MolData; n: number; seed: number };

const ctx = self as unknown as { postMessage: (m: unknown) => void; onmessage: ((e: MessageEvent<WorkerRequest>) => void) | null };

ctx.onmessage = (e) => {
  const req = e.data;
  try {
    if (req.type === 'generate') {
      const r = generate3D(fromData(req.mol), { conformers: req.conformers, seed: req.seed });
      ctx.postMessage({ id: req.id, ok: true, mol: toData(r.mol), result: r.result });
    } else if (req.type === 'optimize') {
      const m = fromData(req.mol);
      const result = optimizeGeometry(m, { maxIter: 2000 });
      ctx.postMessage({ id: req.id, ok: true, mol: toData(m), result });
    } else if (req.type === 'conformers') {
      const m = fromData(req.mol);
      const r = conformerSearch(m, { n: req.n, seed: req.seed, maxTimeMs: 8000 });
      const result = optimizeGeometry(r.best, { maxIter: 500 });
      ctx.postMessage({ id: req.id, ok: true, mol: toData(r.best), result, energies: r.energies });
    }
  } catch (err) {
    ctx.postMessage({ id: req.id, ok: false, error: String((err as Error)?.message ?? err) });
  }
};
