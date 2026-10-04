// Mechanism workbench: apply electron-pushing arrows to generate intermediates, check electron counts.
import type { App } from '../app';
import { h, clear, svgEl, toast } from '../dom';
import { ICONS } from '../icons';
import { applyArrows, arrowGroups, octetViolations, MechanismWarning } from '../../doc/mechanism';
import { insertMol, docBounds } from '../../doc/document';
import { clean2D } from '../../chem/clean2d';

export class MechanismPanel {
  el: HTMLElement;
  private status: HTMLElement;
  private warnings: HTMLElement;

  constructor(private app: App) {
    this.status = h('div', { class: 'mech-status' });
    this.warnings = h('div', { class: 'mech-warnings' });
    const ed = app.editor;
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
            h('li', null, 'Charges, radicals and lone pairs of the product are recomputed from electron counts; octet violations are flagged.'),
          ),
        ),
      ),
    );
  }

  update(): void {
    const groups = arrowGroups(this.app.editor.doc);
    const n = groups.reduce((s, g) => s + g.arrows.length, 0);
    this.status.textContent = n ? `${n} curved arrow${n > 1 ? 's' : ''} in ${groups.length} step${groups.length > 1 ? 's' : ''}` : 'No curved arrows yet.';
  }

  /** Applies the selected curved arrows, or the right-most group of arrows. */
  applyStep(): void {
    const ed = this.app.editor;
    const doc = ed.doc;
    let arrowIds = [...ed.sel.objects].filter((id) => doc.curved.has(id));
    if (!arrowIds.length) {
      const groups = arrowGroups(doc);
      if (!groups.length) return toast('Draw electron-pushing arrows first (curved-arrow tool, Shift+A)', 'error');
      groups.sort((a, b) => b.maxX - a.maxX);
      arrowIds = groups[0].arrows;
    }
    const r = applyArrows(doc, arrowIds);
    this.showWarnings(r.warnings);
    if (!r.product.atoms.length) return;
    // place the product to the right of everything involved
    const ids = new Set(r.reactantAtomIds);
    const b = docBounds(doc, ids);
    const all = docBounds(doc);
    if (!b || !all) return;
    const width = b.maxX - b.minX;
    const arrowGap = 3.4;
    // start right of the reactants, but after anything else drawn in the same horizontal band
    let bandMax = b.maxX;
    const inBand = (y1: number, y2: number) => y2 >= b.minY - 1.5 && y1 <= b.maxY + 1.5;
    for (const a of doc.atoms.values()) if (!ids.has(a.id) && a.x > b.minX && inBand(a.y, a.y)) bandMax = Math.max(bandMax, a.x);
    for (const o of doc.arrows.values()) if (Math.max(o.x1, o.x2) > b.minX && inBand(Math.min(o.y1, o.y2), Math.max(o.y1, o.y2))) bandMax = Math.max(bandMax, o.x1, o.x2);
    for (const t of doc.texts.values()) if (t.x > b.maxX && inBand(t.y, t.y)) bandMax = Math.max(bandMax, t.x + 1);
    const startX = bandMax + 0.8;
    void all;
    // tidy bonds that were created between separate fragments
    const formedLong = r.product.bonds.some((bd) => {
      const A = r.product.atoms[bd.a], B = r.product.atoms[bd.b];
      return Math.hypot(A.x - B.x, A.y - B.y) > 1.6;
    });
    if (formedLong) {
      try {
        clean2D(r.product);
      } catch {
        /* keep raw geometry */
      }
    }
    const pb = r.product.bbox();
    const midY = (b.minY + b.maxY) / 2;
    ed.begin();
    const aid = doc.nextId++;
    doc.arrows.set(aid, { id: aid, type: 'arrow', kind: r.resonance ? 'resonance' : 'reaction', x1: startX, y1: midY, x2: startX + arrowGap - 0.8, y2: midY });
    const dx = startX + arrowGap - pb.minX + 0.2;
    const dy = midY - (pb.minY + pb.maxY) / 2;
    // product atoms are new copies
    for (const a of r.product.atoms) delete (a as { lonePairs?: boolean }).lonePairs;
    const { atomIds } = insertMol(doc, r.product, dx, dy);
    // keep lone-pair display for atoms that showed them in the reactant
    r.product.atoms.forEach((a, i) => {
      const src = doc.atoms.get(a.id);
      if (src?.lonePairs) doc.atoms.get(atomIds[i])!.lonePairs = true;
    });
    ed.commit(r.resonance ? 'Resonance structure' : 'Mechanism step');
    void width;
    ed.fitToContent();
    const summary = r.summary.length ? r.summary.join('; ') : 'no change';
    toast(`${r.resonance ? 'Resonance structure' : 'Intermediate'} generated: ${summary}`, r.warnings.length ? 'info' : 'success', 4000);
  }

  private showWarnings(ws: MechanismWarning[]): void {
    clear(this.warnings);
    const halo = new Map<number, string>();
    for (const w of ws) {
      this.warnings.appendChild(h('div', { class: 'warn small' }, '⚠ ' + w.message));
      for (const id of w.atomIds ?? []) halo.set(id, 'rgba(224,49,49,0.25)');
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
