// Public API of the 3D module: structure generation, UFF optimisation, conformer search, XYZ export.
import { Mol } from '../mol';
import { embed3D, embedCoordinates, EmbedOptions } from './embed';
import { UFF, UFFOptions, TermEnergies, emptyTerms, coordsOf, setCoords } from './uff';
import { lbfgs } from './minimize';
import { conformerSearch, ConformerOptions } from './conformers';
import { stereoReport, StereoReport } from './stereo3d';

export { embed3D, embedCoordinates, conformerSearch, stereoReport, UFF, coordsOf, setCoords };
export type { EmbedOptions, ConformerOptions, StereoReport, UFFOptions, TermEnergies };
export { addExplicitHydrogens } from './hydrogens';
export { typeUFF } from './uffTyping';
export { dihedral, angle, tetraSatisfied, dbSatisfied } from './stereo3d';

export interface OptResult {
  /** Total UFF energy (kcal/mol). */
  energy: number;
  terms: Record<'stretch' | 'bend' | 'torsion' | 'oop' | 'vdw' | 'elec', number>;
  converged: boolean;
  iterations: number;
  /** RMS gradient component at the end (kcal/mol/Å). */
  rmsGrad: number;
}

export interface OptimizeOptions extends UFFOptions {
  /** Maximum L-BFGS iterations (default 2000). */
  maxIter?: number;
  /** Convergence threshold on the RMS gradient (default 1e-3 kcal/mol/Å). */
  gradTol?: number;
  onProgress?: (it: number, e: number) => void;
}

/** Minimises the UFF energy of an H-explicit 3D molecule in place (atoms[i].x/y/z). */
export function optimizeGeometry(mol3d: Mol, opts: OptimizeOptions = {}): OptResult {
  const ff = new UFF(mol3d, opts);
  const x = coordsOf(mol3d);
  const r = lbfgs((xx, g) => ff.energy(xx, g), x, {
    maxIter: opts.maxIter ?? 2000,
    gradTol: opts.gradTol ?? 1e-3,
    maxStep: 0.3,
    onProgress: opts.onProgress,
  });
  setCoords(mol3d, x);
  const terms = emptyTerms();
  const energy = ff.energy(x, null, terms);
  return { energy, terms, converged: r.converged, iterations: r.iterations, rmsGrad: r.rmsGrad };
}

/** UFF energy (and per-term breakdown) of the current coordinates, without moving atoms. */
export function uffEnergy(mol3d: Mol, opts: UFFOptions = {}): { energy: number; terms: TermEnergies } {
  const ff = new UFF(mol3d, opts);
  const terms = emptyTerms();
  const energy = ff.energy(coordsOf(mol3d), null, terms);
  return { energy, terms };
}

/** XYZ file text: atom count, comment line, then "El x y z" (Å) per atom. */
export function toXYZ(mol3d: Mol, comment?: string): string {
  const lines = [String(mol3d.atoms.length), (comment ?? mol3d.name ?? '').replace(/[\r\n]+/g, ' ')];
  for (const a of mol3d.atoms) {
    const el = a.el === 'R' || a.el === '*' ? 'X' : a.el;
    const f = (v: number | undefined) => (v ?? 0).toFixed(6).padStart(12);
    lines.push(`${el.padEnd(2)} ${f(a.x)} ${f(a.y)} ${f(a.z)}`);
  }
  return lines.join('\n') + '\n';
}

/**
 * Convenience: embed + optimise (+ optional quick conformer search when `conformers` > 1).
 * Stereo specs are verified on the final structure; on failure a different seed is tried.
 */
export function generate3D(mol: Mol, opts: { conformers?: number; seed?: number } = {}): { mol: Mol; result: OptResult } {
  const seed = opts.seed ?? 1;
  let fallback: { mol: Mol; result: OptResult } | null = null;
  for (let attempt = 0; attempt < 3; attempt++) {
    let m = embed3D(mol, { seed: seed + attempt * 7919 });
    let result = optimizeGeometry(m);
    if ((opts.conformers ?? 0) > 1) {
      const cs = conformerSearch(m, { n: opts.conformers, seed: seed + attempt });
      m = cs.best;
      result = optimizeGeometry(m);
    }
    const out = { mol: m, result };
    if (stereoReport(m).allOK) return out;
    fallback ??= out;
  }
  return fallback!;
}
