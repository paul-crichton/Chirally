// Mechanism workbench: apply electron-pushing arrows to generate intermediates, check electron counts.
import type { App } from '../app';
import { h, clear, svgEl, toast } from '../dom';
import { ICONS } from '../icons';
import { applyArrows, arrowGroups, isApplied, octetViolations, placeStep, MechanismResult, MechanismWarning, StepPlacement } from '../../doc/mechanism';
import { createDoc, insertMol, docBounds } from '../../doc/document';
import { ChemDoc } from '../../doc/types';

/** A step whose product is not a valid structure: shown as a ghost until the user inserts or discards it. */
interface PendingStep {
  result: MechanismResult;
  placement: StepPlacement;
  arrowIds: number[];
  /** Document revision the preview was computed for; any edit invalidates it. */
  rev: number;
}

export class MechanismPanel {
  el: HTMLElement;
  private status: HTMLElement;
  private warnings: HTMLElement;
  private pendingBox: HTMLElement;
  private pending: PendingStep | null = null;

  constructor(private app: App) {
    this.status = h('div', { class: 'mech-status' });
    this.warnings = h('div', { class: 'mech-warnings' });
    this.pendingBox = h('div', { class: 'mech-pending', hidden: true });
    const ed = app.editor;
    ed.on('change', () => {
      if (this.pending && this.pending.rev !== ed.rev) this.clearPending();
    });
    this.el = h(
      'section',
      { class: 'panel', 'aria-label': 'Mechanism' },
      h('div', { class: 'panel-body' },
        h('p', { class: 'small' },
          'Draw curved arrows from a ', h('b', null, 'lone pair (atom)'), ' or ', h('b', null, 'bond'), ' to an ', h('b', null, 'atom'), ', a ', h('b', null, 'bond'),
          ' or the ', h('b', null, 'space between two atoms'), '. Chirally does the electron bookkeeping and draws the next intermediate for you.'),
        h('div', { class: 'row-actions' },
          h('button', { class: 'btn', onclick: () => { ed.settings.curved = 2; ed.setTool('curved'); } }, svgEl(ICONS.curved2), 'Electron pair'),
          h('button', { class: 'btn', onclick: () => { ed.settings.curved = 1; ed.setTool('curved'); } }, svgEl(ICONS.curved1), 'Single electron'),
        ),
        h('button', { class: 'btn btn-primary btn-block', onclick: () => this.applyStep() }, svgEl(ICONS.mech), 'Apply arrows → next intermediate'),
        this.status,
        this.pendingBox,
        this.warnings,
        h('div', { class: 'row-actions wrap' },
          h('button', { class: 'btn btn-small', onclick: () => this.checkElectrons() }, svgEl(ICONS.check), 'Check electron counts'),
          h('button', { class: 'btn btn-small', onclick: () => this.toggleLonePairs() }, svgEl(ICONS.lonepair), 'Toggle all lone pairs'),
          h('button', { class: 'btn btn-small', onclick: () => { ed.settings.shape = 'tsBracket'; ed.setTool('shape'); } }, svgEl(ICONS.tsBracket), 'TS brackets'),
        ),
        h('details', { class: 'sub' }, h('summary', null, 'How arrows are interpreted'),
          h('ul', { class: 'small' },
            h('li', null, 'Lone pair → atom: forms a new bond (nucleophile attacks electrophile).'),
            h('li', null, 'Lone pair → adjacent bond: makes it a π bond.'),
            h('li', null, 'Bond → one of its atoms: heterolysis, the electrons stay on that atom.'),
            h('li', null, 'Bond → another bond / atom: shifts the π/σ electrons (resonance, additions).'),
            h('li', null, 'Fishhooks move one electron each: two fishhooks from one bond = homolysis (radicals).'),
            h('li', null, 'Proton transfers: start or end an arrow on the H of a label (OH, NH₃⁺, OH₂⁺) — as a target it is that H, as a source its X–H bond. For C–H bonds, right-click the atom → Draw hydrogens as atoms.'),
            h('li', null, 'If an arrow could mean two different bonds, Chirally keeps octets intact first, then follows the chemistry (1,2-shifts, hydrogen transfer, bond polarity, Markovnikov and Michael selectivity); it tells you when only the drawing decided.'),
            h('li', null, 'Charges, radicals and lone pairs of the product are recomputed from electron counts; octet violations are flagged.'),
            h('li', null, 'Rings drawn with delocalised (aromatic) bonds are given alternating double bonds for the step; rings the arrows don’t touch stay delocalised.'),
            h('li', null, 'Labels such as OMe, CN⁻ or H₂O take part as the atoms they stand for; E⁺ and Nu⁻ act as one-bond atoms. The total charge must not change.'),
            h('li', null, 'Stereocentres keep their wedges when their bonds are untouched; a substitution at a stereocentre (one bond formed, one broken) inverts it, as in SN2.'),
            h('li', null, 'If the result breaks the octet rule or electrons go missing, it is previewed in red instead of being added; you can still insert it.'),
            h('li', null, 'Apply picks the most recently drawn arrows that have not been applied yet; select arrows to apply a particular step (again).'),
          ),
        ),
      ),
    );
  }

  update(): void {
    const groups = arrowGroups(this.app.editor.doc);
    const n = groups.reduce((s, g) => s + g.arrows.length, 0);
    const done = groups.filter((g) => g.applied).length;
    this.status.textContent = n
      ? `${n} curved arrow${n > 1 ? 's' : ''} in ${groups.length} step${groups.length > 1 ? 's' : ''}${done ? ` · ${done === groups.length ? 'all' : done} applied` : ''}`
      : 'No curved arrows yet.';
  }

  /**
   * Applies `ids`, or else the selected curved arrows, or else the most recently drawn group of arrows that
   * has not been applied yet. Valid steps are added to the drawing; invalid ones are only previewed.
   */
  applyStep(ids?: number[]): void {
    const ed = this.app.editor;
    const doc = ed.doc;
    let arrowIds = (ids ?? [...ed.sel.objects]).filter((id) => doc.curved.has(id));
    if (arrowIds.length) {
      if (arrowIds.some((id) => isApplied(doc, doc.curved.get(id)!)) && !confirm('These arrows were already applied. Apply them again?')) return;
    } else {
      const groups = arrowGroups(doc);
      if (!groups.length) return toast('Draw electron-pushing arrows first (curved-arrow tool, Shift+A)', 'error');
      const open = groups.filter((g) => !g.applied);
      if (!open.length) return toast('All curved arrows have been applied. Select arrows to apply them again.', 'info', 4000);
      open.sort((a, b) => Math.max(...b.arrows) - Math.max(...a.arrows) || b.maxX - a.maxX);
      arrowIds = open[0].arrows;
    }
    this.clearPending();
    const r = applyArrows(doc, arrowIds);
    this.showWarnings(r.warnings);
    const placement = r.changed ? placeStep(doc, r) : null;
    if (!placement) {
      toast(r.ok ? 'These arrows don’t change anything' : 'No intermediate generated — see the Mechanism panel', r.ok ? 'info' : 'error', 4000);
      return;
    }
    if (r.ok) {
      this.insertStep(r, placement, arrowIds);
      return;
    }
    this.pending = { result: r, placement, arrowIds, rev: ed.rev };
    this.showPending();
    toast('Not a valid intermediate — previewed in red, not added', 'error', 4000);
  }

  /** Applies the arrow group that contains the given curved arrow (context menu). */
  applyGroupOf(curvedId: number): void {
    const g = arrowGroups(this.app.editor.doc).find((x) => x.arrows.includes(curvedId));
    if (g) this.applyStep(g.arrows);
  }

  private insertStep(r: MechanismResult, placement: StepPlacement, arrowIds: number[]): void {
    const ed = this.app.editor;
    const doc = ed.doc;
    ed.begin();
    const aid = doc.nextId++;
    doc.arrows.set(aid, { id: aid, type: 'arrow', ...placement.arrow });
    this.addProduct(doc, placement);
    // remember the step so Apply moves on to the next one
    for (const id of arrowIds) {
      const c = doc.curved.get(id);
      if (c) c.step = aid;
    }
    ed.commit(r.resonance ? 'Resonance structure' : 'Mechanism step');
    ed.fitToContent();
    const summary = r.summary.length ? r.summary.join('; ') : 'no change';
    toast(`${r.resonance ? 'Resonance structure' : 'Intermediate'} generated: ${summary}`, r.warnings.length ? 'info' : 'success', 4000);
  }

  /** Adds the placed product to `target`, keeping lone-pair display from the reactant atoms. Returns the new atom ids. */
  private addProduct(target: ChemDoc, placement: StepPlacement): number[] {
    const src = this.app.editor.doc;
    const mol = placement.mol.clone();
    for (const a of mol.atoms) delete (a as { lonePairs?: boolean }).lonePairs;
    const { atomIds } = insertMol(target, mol);
    placement.mol.atoms.forEach((a, i) => {
      if (src.atoms.get(a.id)?.lonePairs) target.atoms.get(atomIds[i])!.lonePairs = true;
    });
    return atomIds;
  }

  private showPending(): void {
    const p = this.pending!;
    const ed = this.app.editor;
    const ghost = createDoc({ ...ed.doc.style });
    const aid = ghost.nextId++;
    ghost.arrows.set(aid, { id: aid, type: 'arrow', ...p.placement.arrow });
    const atomIds = this.addProduct(ghost, p.placement);
    const bad = new Set(p.result.warnings.filter((w) => w.level === 'error').flatMap((w) => w.atomIds ?? []));
    const halo = new Map<number, string>();
    p.placement.mol.atoms.forEach((a, i) => {
      if (bad.has(a.id)) halo.set(atomIds[i], 'rgba(224,49,49,0.35)');
    });
    ed.setGhost(ghost, { halo, label: 'Preview — not a valid intermediate' });
    // make sure the preview is in view
    const all = docBounds(ed.doc), gb = docBounds(ghost);
    if (all && gb) ed.fitToContent({ x1: Math.min(all.minX, gb.minX), y1: Math.min(all.minY, gb.minY) - 0.8, x2: Math.max(all.maxX, gb.maxX), y2: Math.max(all.maxY, gb.maxY) });
    clear(this.pendingBox);
    this.pendingBox.append(
      h('div', { class: 'small' }, h('b', null, 'Not added: '), 'the product breaks the rules listed below. Fix the arrows and apply again, or insert it anyway.'),
      h('div', { class: 'row-actions' },
        h('button', { class: 'btn btn-small', onclick: () => this.insertPending() }, 'Insert anyway'),
        h('button', { class: 'btn btn-small', onclick: () => this.clearPending() }, 'Discard'),
      ),
    );
    this.pendingBox.hidden = false;
  }

  private insertPending(): void {
    const p = this.pending;
    if (!p || p.rev !== this.app.editor.rev) {
      this.clearPending();
      toast('The drawing changed — apply the arrows again', 'info');
      return;
    }
    this.clearPending();
    this.insertStep(p.result, p.placement, p.arrowIds);
  }

  private clearPending(): void {
    if (!this.pending && this.pendingBox.hidden) return;
    this.pending = null;
    this.app.editor.setGhost(null);
    clear(this.pendingBox);
    this.pendingBox.hidden = true;
  }

  private showWarnings(ws: MechanismWarning[]): void {
    clear(this.warnings);
    const halo = new Map<number, string>();
    for (const w of ws) {
      const isError = w.level === 'error';
      this.warnings.appendChild(h('div', { class: `${isError ? 'err' : 'warn'} small` }, (isError ? '✖ ' : '⚠ ') + w.message));
      for (const id of w.atomIds ?? []) halo.set(id, isError ? 'rgba(224,49,49,0.25)' : 'rgba(232,89,12,0.22)');
    }
    this.app.setWarningHalos(halo.size ? halo : null);
  }

  checkElectrons(): void {
    const v = octetViolations(this.app.editor.doc);
    clear(this.warnings);
    if (!v.length) {
      this.warnings.appendChild(h('div', { class: 'ok small' }, '✓ No octet/duet violations found.'));
      this.app.setWarningHalos(null);
      return;
    }
    const halo = new Map<number, string>();
    for (const x of v) {
      const a = this.app.editor.doc.atoms.get(x.atomId)!;
      this.warnings.appendChild(h('div', { class: 'warn small' }, `⚠ ${a.el} has ${x.count} valence electrons (limit ${x.limit})`));
      halo.set(x.atomId, 'rgba(224,49,49,0.25)');
    }
    this.app.setWarningHalos(halo);
  }

  toggleLonePairs(): void {
    const ed = this.app.editor;
    ed.doc.style.showLonePairs = !ed.doc.style.showLonePairs;
    ed.touch();
    this.app.stylePanel?.refresh();
  }
}
