// 3D structure generation, UFF geometry optimisation, conformer search and interactive viewer.
import type { App } from '../app';
import { h, svgEl, toast, fmt, downloadBlob } from '../dom';
import { ICONS } from '../icons';
import { Mol } from '../../chem/mol';
import { Viewer3D, Style3D } from '../../viewer3d/viewer';
import { toXYZ, generate3D, optimizeGeometry, conformerSearch } from '../../chem/3d';
import { MolData, toData, fromData } from '../../workers/moldata';
import { writeMolfile } from '../../chem/formats';

interface WorkerReply {
  id: number;
  ok: boolean;
  mol?: MolData;
  result?: { energy: number; terms: Record<string, number>; converged: boolean; iterations: number; rmsGrad: number };
  energies?: number[];
  error?: string;
}

export class Viewer3DPanel {
  el: HTMLElement;
  private viewer: Viewer3D | null = null;
  private host: HTMLElement;
  private info: HTMLElement;
  private mol3d: Mol | null = null;
  private worker: Worker | null = null;
  private pending = new Map<number, { resolve: (r: WorkerReply) => void; fallback: () => WorkerReply; timer: ReturnType<typeof setTimeout> }>();
  private reqId = 0;
  private busy = false;
  private source = '';

  constructor(private app: App) {
    this.host = h('div', { class: 'viewer-host', tabindex: 0, 'aria-label': '3D viewer — drag to rotate, scroll or pinch to zoom' });
    this.info = h('div', { class: 'viewer-info muted small' }, 'Generate a 3D model of the selection (or the whole drawing). The geometry is optimised with the UFF force field in your browser.');
    const styleSel = h('select', { class: 'input input-small', 'aria-label': 'Display style' },
      h('option', { value: 'ballstick' }, 'Ball & stick'),
      h('option', { value: 'sticks' }, 'Sticks'),
      h('option', { value: 'spacefill' }, 'Space-filling'),
      h('option', { value: 'wireframe' }, 'Wireframe'),
    ) as HTMLSelectElement;
    styleSel.addEventListener('change', () => this.viewer?.setStyle(styleSel.value as Style3D));
    const check = (label: string, on: boolean, fn: (v: boolean) => void) => {
      const cb = h('input', { type: 'checkbox', checked: on }) as HTMLInputElement;
      cb.addEventListener('change', () => fn(cb.checked));
      return h('label', { class: 'check' }, cb, label);
    };
    this.el = h(
      'section',
      { class: 'panel panel-3d', 'aria-label': '3D model' },
      h('div', { class: 'panel-body' },
        h('div', { class: 'row-actions wrap' },
          h('button', { class: 'btn btn-primary', onclick: () => this.generate() }, svgEl(ICONS.cube), 'Generate 3D'),
          h('button', { class: 'btn', onclick: () => this.optimize(), title: 'Further minimise the current geometry' }, 'Optimize'),
          h('button', { class: 'btn', onclick: () => this.conformers(), title: 'Random-perturbation conformer search; keeps the lowest-energy conformer' }, 'Conformer search'),
        ),
        this.host,
        h('div', { class: 'row-actions wrap' },
          styleSel,
          check('Labels', false, (v) => this.viewer?.setShowLabels(v)),
          check('Hydrogens', true, (v) => this.viewer?.setShowHydrogens(v)),
          check('Spin', false, (v) => this.viewer?.setAutoRotate(v)),
          check('Measure', false, (v) => this.viewer?.setMeasureMode(v)),
        ),
        this.info,
        h('div', { class: 'row-actions wrap' },
          h('button', { class: 'btn btn-small', onclick: () => this.exportXYZ() }, 'XYZ'),
          h('button', { class: 'btn btn-small', onclick: () => this.exportMol() }, 'MOL (3D)'),
          h('button', { class: 'btn btn-small', onclick: () => this.exportPNG() }, 'PNG'),
          h('button', { class: 'btn btn-small', onclick: () => this.viewer?.resetView() }, 'Reset view'),
        ),
      ),
    );
  }

  /** Lazily creates the viewer when the panel becomes visible. */
  activate(): void {
    if (!this.viewer) {
      this.viewer = new Viewer3D(this.host);
      this.viewer.onMeasure = (t) => {
        this.info.textContent = t;
      };
      if (this.mol3d) this.viewer.setMolecule(this.mol3d);
    }
  }

  private getWorker(): Worker | null {
    if (this.worker) return this.worker;
    try {
      this.worker = new Worker(new URL('../../workers/chem3d.worker.ts', import.meta.url), { type: 'module' });
      this.worker.onmessage = (e: MessageEvent<WorkerReply>) => {
        const p = this.pending.get(e.data.id);
        if (p) {
          this.pending.delete(e.data.id);
          clearTimeout(p.timer);
          p.resolve(e.data);
        }
      };
      this.worker.onerror = (e) => {
        // e.g. the worker chunk 404s in a tab opened before a redeploy, or offline without it cached:
        // answer the waiting requests on the main thread now instead of after the 60 s timeout
        e.preventDefault();
        this.worker?.terminate();
        this.worker = null;
        for (const p of this.pending.values()) {
          clearTimeout(p.timer);
          p.resolve(p.fallback());
        }
        this.pending.clear();
      };
    } catch {
      this.worker = null;
    }
    return this.worker;
  }

  private run(req: Record<string, unknown>, fallback: () => WorkerReply): Promise<WorkerReply> {
    const w = this.getWorker();
    if (!w) return Promise.resolve(fallback());
    const id = ++this.reqId;
    return new Promise((resolve) => {
      const timer = setTimeout(() => {
        // worker hung: fall back to main thread
        this.pending.delete(id);
        resolve(fallback());
      }, 60000);
      this.pending.set(id, { resolve, fallback, timer });
      w.postMessage({ ...req, id });
    });
  }

  async generate(): Promise<void> {
    if (this.busy) return;
    const mol = this.app.chemMolForSelection();
    if (!mol || !mol.atoms.length) return toast('Draw a structure first', 'error');
    if (mol.atoms.some((a) => a.el === 'R' || a.el === '*')) return toast('3D needs a fully specified structure (no R groups)', 'error');
    this.activate();
    this.setBusy(true, 'Embedding and optimising (UFF)…');
    const t0 = performance.now();
    const r = await this.run({ type: 'generate', mol: toData(mol), conformers: mol.atoms.length < 40 ? 8 : 3, seed: 7 }, () => {
      const g = generate3D(mol, { conformers: 3, seed: 7 });
      return { id: 0, ok: true, mol: toData(g.mol), result: g.result as WorkerReply['result'] };
    });
    this.setBusy(false);
    if (!r.ok || !r.mol) return this.fail(r.error);
    this.source = this.app.currentSmiles();
    this.show(fromData(r.mol), r.result, performance.now() - t0);
  }

  async optimize(): Promise<void> {
    if (!this.mol3d || this.busy) return this.generate();
    const m = this.mol3d;
    this.setBusy(true, 'Optimising…');
    const t0 = performance.now();
    const r = await this.run({ type: 'optimize', mol: toData(m) }, () => {
      const c = fromData(toData(m));
      const result = optimizeGeometry(c, { maxIter: 2000 });
      return { id: 0, ok: true, mol: toData(c), result: result as WorkerReply['result'] };
    });
    this.setBusy(false);
    if (!r.ok || !r.mol) return this.fail(r.error);
    this.show(fromData(r.mol), r.result, performance.now() - t0);
  }

  async conformers(): Promise<void> {
    if (!this.mol3d) await this.generate();
    if (!this.mol3d || this.busy) return;
    const m = this.mol3d;
    this.setBusy(true, 'Searching conformers…');
    const t0 = performance.now();
    const r = await this.run({ type: 'conformers', mol: toData(m), n: 30, seed: Date.now() % 100000 }, () => {
      const c = conformerSearch(fromData(toData(m)), { n: 10, seed: 3, maxTimeMs: 4000 });
      const result = optimizeGeometry(c.best, { maxIter: 500 });
      return { id: 0, ok: true, mol: toData(c.best), result: result as WorkerReply['result'], energies: c.energies };
    });
    this.setBusy(false);
    if (!r.ok || !r.mol) return this.fail(r.error);
    this.show(fromData(r.mol), r.result, performance.now() - t0, r.energies);
  }

  private show(m: Mol, result: WorkerReply['result'] | undefined, ms: number, energies?: number[]): void {
    this.mol3d = m;
    this.viewer?.setMolecule(m);
    if (result) {
      const t = result.terms;
      this.info.innerHTML = '';
      this.info.append(
        h('div', null, h('b', null, `E(UFF) = ${fmt(result.energy, 2)} kcal/mol`), ` · ${result.converged ? 'converged' : 'not fully converged'} in ${result.iterations} steps · ${fmt(ms / 1000, 2)} s`),
        h('div', null, `stretch ${fmt(t.stretch, 1)} · bend ${fmt(t.bend, 1)} · torsion ${fmt(t.torsion, 1)} · oop ${fmt(t.oop, 1)} · vdW ${fmt(t.vdw, 1)}`),
        energies && energies.length ? h('div', null, `${energies.length} conformers; lowest ${fmt(Math.min(...energies), 2)}, highest ${fmt(Math.max(...energies), 2)} kcal/mol`) : '',
        h('div', null, 'Drag to rotate · scroll/pinch to zoom · Measure: click 2–4 atoms for distance/angle/dihedral'),
      );
    }
  }

  private setBusy(on: boolean, msg = ''): void {
    this.busy = on;
    this.el.classList.toggle('busy', on);
    if (on) this.info.textContent = msg;
  }

  private fail(err?: string): void {
    toast('3D generation failed' + (err ? ': ' + err : ''), 'error');
    this.info.textContent = err ?? 'Failed';
  }

  private exportXYZ(): void {
    if (!this.mol3d) return toast('Generate a 3D model first', 'error');
    downloadBlob(new Blob([toXYZ(this.mol3d, this.source || 'ChemWrite 3D model')], { type: 'chemical/x-xyz' }), 'structure3d.xyz');
  }

  private exportMol(): void {
    if (!this.mol3d) return toast('Generate a 3D model first', 'error');
    const m = this.mol3d;
    m.props = { ...m.props, dim: '3D' };
    downloadBlob(new Blob([writeMolfile(m, { title: 'ChemWrite 3D (UFF)' })], { type: 'chemical/x-mdl-molfile' }), 'structure3d.mol');
  }

  private exportPNG(): void {
    if (!this.viewer || !this.mol3d) return toast('Generate a 3D model first', 'error');
    const url = this.viewer.toPNG(2);
    const a = h('a', { href: url, download: 'structure3d.png' });
    document.body.appendChild(a);
    a.click();
    a.remove();
  }
}
