// Live structure analysis: name, formula, masses, properties, formulation insights, reaction balance.
import type { App } from '../app';
import { h, clear, fmt, copyText, svgEl } from '../dom';
import { ICONS } from '../icons';
import { analyze, Analysis, FragmentAnalysis } from '../chem';
import { summarize, FormulaInfo } from '../../chem/formula';
import { docToMol } from '../../doc/document';
import { computeFormula } from '../../chem/formula';

export class AnalysisPanel {
  el: HTMLElement;
  private body: HTMLElement;
  last: Analysis | null = null;

  constructor(private app: App) {
    this.body = h('div', { class: 'panel-body' });
    this.el = h('section', { class: 'panel', 'aria-label': 'Analysis' }, this.body);
  }

  update(): void {
    const ed = this.app.editor;
    const doc = ed.doc;
    clear(this.body);
    if (!doc.atoms.size) {
      this.last = null;
      this.body.appendChild(this.emptyState());
      this.app.onAnalysis(null);
      return;
    }
    const ids = ed.analysisAtomIds();
    let a: Analysis;
    try {
      a = analyze(doc, ids);
    } catch (e) {
      this.body.appendChild(h('p', { class: 'muted' }, 'Analysis failed: ' + String(e)));
      return;
    }
    this.last = a;
    this.app.onAnalysis(a);
    const scope = ed.hasSelection() ? 'Selection' : 'Drawing';
    const head = h('div', { class: 'scope' }, h('span', { class: 'chip' }, scope), a.fragments.length > 1 ? h('span', { class: 'muted' }, `${a.fragments.length} fragments`) : null);
    this.body.appendChild(head);
    if (a.fragments.length > 1 && a.total) this.body.appendChild(this.totalCard(a.total));
    const frags = a.fragments.slice(0, 6);
    frags.forEach((f) => this.body.appendChild(this.fragmentCard(f, a.fragments.length > 1)));
    if (a.fragments.length > frags.length) this.body.appendChild(h('p', { class: 'muted' }, `…and ${a.fragments.length - frags.length} more fragments`));
    const rx = this.reactionCard();
    if (rx) this.body.appendChild(rx);
  }

  private emptyState(): HTMLElement {
    return h(
      'div',
      { class: 'empty-state' },
      svgEl(ICONS.flask),
      h('h3', null, 'Draw or import a structure'),
      h('p', null, 'Name, formula, mass, properties and formulation insights appear here as you draw.'),
      h('div', { class: 'empty-actions' },
        h('button', { class: 'btn btn-primary', onclick: () => this.app.openPalette('') }, svgEl(ICONS.search), 'Search PubChem / insert SMILES'),
        h('button', { class: 'btn', onclick: () => this.app.loadExample() }, svgEl(ICONS.sparkle), 'Load an example'),
      ),
    );
  }

  private row(label: string, value: Node | string, copy?: string, title?: string): HTMLElement {
    return h(
      'div',
      { class: 'kv', title },
      h('span', { class: 'k' }, label),
      h('span', { class: 'v' }, value),
      copy ? h('button', { class: 'icon-btn tiny', title: 'Copy', 'aria-label': `Copy ${label}`, onclick: () => copyText(copy, label) }, svgEl(ICONS.copy)) : null,
    );
  }

  private formulaHTML(f: FormulaInfo): HTMLElement {
    return h('span', { class: 'formula', html: f.html });
  }

  private totalCard(t: FormulaInfo): HTMLElement {
    return h('div', { class: 'card' }, h('h4', null, 'All fragments'), this.row('Formula', this.formulaHTML(t), t.formula), this.row('MW', `${fmt(t.mw, 2)} g/mol`), this.row('Exact mass', fmt(t.exactMass, 4)));
  }

  private fragmentCard(f: FragmentAnalysis, multi: boolean): HTMLElement {
    const card = h('div', { class: 'card' });
    card.addEventListener('mouseenter', () => multi && this.app.highlightAtoms(f.atomIds, 'rgba(250,176,5,0.28)'));
    card.addEventListener('mouseleave', () => multi && this.app.highlightAtoms(null));
    // name
    const nameBox = h('div', { class: 'name-box' });
    if (f.name) {
      nameBox.appendChild(h('div', { class: 'iupac', title: 'IUPAC name (generated locally)' }, f.name));
      nameBox.appendChild(
        h('div', { class: 'name-actions' },
          h('button', { class: 'btn btn-small', onclick: () => copyText(f.name!, 'Name') }, svgEl(ICONS.copy), 'Copy'),
          h('button', { class: 'btn btn-small', onclick: () => this.app.identify(f.smiles) }, svgEl(ICONS.pubchem), 'PubChem'),
          h('button', { class: 'btn btn-small', title: 'Show IUPAC locants on the structure', onclick: () => this.app.toggleLocants() }, '1,2,3'),
        ),
      );
    } else if (f.hasPseudo) {
      nameBox.appendChild(h('div', { class: 'muted' }, 'Generic structure (contains R/X groups) — no name.'));
    } else {
      nameBox.appendChild(h('div', { class: 'muted' }, 'No local IUPAC name for this structure.'));
      nameBox.appendChild(h('button', { class: 'btn btn-small', onclick: () => this.app.identify(f.smiles) }, svgEl(ICONS.pubchem), 'Look up name in PubChem'));
    }
    if (f.nameWarnings.length) nameBox.appendChild(h('div', { class: 'muted small' }, f.nameWarnings.join('; ')));
    card.appendChild(nameBox);
    const fm = f.formula;
    card.appendChild(this.row('Formula', this.formulaHTML(fm), fm.formula));
    card.appendChild(this.row('Mol. weight', `${fmt(fm.mw, 2)} g/mol`, fmt(fm.mw, 3)));
    card.appendChild(this.row('Exact mass', fmt(fm.exactMass, 4), fmt(fm.exactMass, 5)));
    if (fm.charge) card.appendChild(this.row('m/z', fmt(fm.mz, 4)));
    if (f.smiles) card.appendChild(this.row('SMILES', h('code', { class: 'smiles' }, f.smiles), f.smiles));
    const p = f.props;
    if (p) {
      const det = h('details', { class: 'sub', open: true }, h('summary', null, 'Properties'));
      det.appendChild(this.row('cLogP', fmt(p.logP, 2), undefined, 'Wildman–Crippen octanol/water partition coefficient'));
      det.appendChild(this.row('TPSA', `${fmt(p.tpsa, 1)} Å²`, undefined, 'Topological polar surface area (Ertl)'));
      det.appendChild(this.row('H-bond donors / acceptors', `${p.hbd} / ${p.hba}`));
      det.appendChild(this.row('Rotatable bonds', String(p.rotatableBonds)));
      det.appendChild(this.row('Rings (aromatic)', `${p.rings} (${p.aromaticRings})`));
      det.appendChild(this.row('Stereocentres', String(p.stereocenters)));
      det.appendChild(this.row('Fraction sp³ C', fmt(p.fractionCsp3, 2)));
      if (fm.dbe !== null) det.appendChild(this.row('Degrees of unsaturation', fmt(fm.dbe, fm.dbe % 1 ? 1 : 0)));
      det.appendChild(this.row('logS (ESOL)', `${fmt(p.logS, 2)} (${solubilityClass(p.logS)})`, undefined, 'Estimated aqueous solubility, log mol/L (Delaney 2004)'));
      det.appendChild(this.row('Molar refractivity', fmt(p.mr, 1)));
      det.appendChild(this.row('Lipinski / Veber', `${p.lipinski.pass ? '✓' : '✗'} (${p.lipinski.violations} violation${p.lipinski.violations === 1 ? '' : 's'}) / ${p.veber.pass ? '✓' : '✗'}`));
      card.appendChild(det);
      // formulation insights
      const fo = h('details', { class: 'sub', open: true }, h('summary', null, 'Formulation insights'));
      const sp = p.skinPenetration;
      fo.appendChild(this.row('MW < 500 Da', sp.under500Da ? '✓ yes' : '✗ no', undefined, 'The "500 Dalton rule" for skin penetration'));
      fo.appendChild(this.row('logP 1–3', sp.logPInRange ? '✓ in range' : '✗ outside', undefined, 'Sweet spot for passive dermal penetration'));
      fo.appendChild(this.row('log Kp (Potts–Guy)', `${fmt(p.logKp, 2)} cm/h`, undefined, 'Estimated skin permeability coefficient'));
      if (sp.note) fo.appendChild(h('p', { class: 'muted small' }, sp.note));
      if (p.hlb) {
        const parts: string[] = [];
        if (p.hlb.davies !== null) parts.push(`Davies ${fmt(p.hlb.davies, 1)}`);
        if (p.hlb.griffin !== null) parts.push(`Griffin ${fmt(p.hlb.griffin, 1)}`);
        fo.appendChild(this.row('HLB', parts.join(' · ') || '—', undefined, 'Hydrophilic–lipophilic balance'));
        const v = p.hlb.griffin ?? p.hlb.davies;
        if (v !== null) fo.appendChild(h('p', { class: 'muted small' }, hlbUse(v)));
        if (p.hlb.note) fo.appendChild(h('p', { class: 'muted small' }, p.hlb.note));
      }
      card.appendChild(fo);
      if (p.functionalGroups.length) {
        const fg = h('details', { class: 'sub', open: true }, h('summary', null, 'Functional groups'));
        const chips = h('div', { class: 'chips' });
        for (const g of p.functionalGroups) {
          const chip = h('button', { class: 'chip chip-btn' + (/alert|Michael|sensiti/i.test(g.name) ? ' chip-warn' : ''), title: 'Hover to highlight on the structure' }, g.count > 1 ? `${g.name} ×${g.count}` : g.name);
          // atoms are indices into the analysed (H-suppressed, expanded) molecule → map back through ids where possible
          chip.addEventListener('mouseenter', () => this.app.highlightFunctionalGroup(f, g.atoms));
          chip.addEventListener('mouseleave', () => this.app.highlightAtoms(null));
          chips.appendChild(chip);
        }
        fg.appendChild(chips);
        card.appendChild(fg);
      }
      const ea = h('details', { class: 'sub' }, h('summary', null, 'Elemental analysis'));
      ea.appendChild(h('div', { class: 'ea' }, fm.composition.map((c) => h('span', null, `${c.el} ${fmt(c.pct, 2)}%`))));
      card.appendChild(ea);
    }
    return card;
  }

  /** Mass balance and atom economy for drawn reactions (left vs right of the first reaction arrow). */
  private reactionCard(): HTMLElement | null {
    const doc = this.app.editor.doc;
    const arrow = [...doc.arrows.values()].find((a) => a.kind === 'reaction' || a.kind === 'equilibrium' || a.kind === 'unbalancedEq' || a.kind === 'dashed');
    if (!arrow) return null;
    const ux = arrow.x2 - arrow.x1, uy = arrow.y2 - arrow.y1;
    const L = Math.hypot(ux, uy) || 1;
    const side = (x: number, y: number) => ((x - arrow.x1) * ux + (y - arrow.y1) * uy) / L;
    const left: number[] = [], right: number[] = [];
    for (const a of doc.atoms.values()) {
      const s = side(a.x, a.y);
      if (s < 0) left.push(a.id);
      else if (s > L) right.push(a.id);
    }
    if (!left.length || !right.length) return null;
    const fl = computeFormula(docToMol(doc, left).mol);
    const fr = computeFormula(docToMol(doc, right).mol);
    const els = new Set([...fl.counts.keys(), ...fr.counts.keys()]);
    const diff = new Map<string, number>();
    for (const e of els) {
      const d = (fr.counts.get(e) ?? 0) - (fl.counts.get(e) ?? 0);
      if (d) diff.set(e, d);
    }
    const card = h('div', { class: 'card' }, h('h4', null, 'Reaction'));
    card.appendChild(this.row('Reactants', this.formulaHTML(fl)));
    card.appendChild(this.row('Products', this.formulaHTML(fr)));
    if (!diff.size && fl.charge === fr.charge) card.appendChild(h('p', { class: 'ok' }, '✓ Balanced (atoms and charge)'));
    else {
      const lost = new Map<string, number>(), gained = new Map<string, number>();
      for (const [e, d] of diff) (d < 0 ? lost : gained).set(e, Math.abs(d));
      const txt = (m: Map<string, number>) => summarize(m, 0).formula;
      card.appendChild(
        h('p', { class: 'warn' },
          '⚠ Not balanced. ',
          lost.size ? `Missing on product side: ${txt(lost)}. ` : '',
          gained.size ? `Extra on product side: ${txt(gained)}. ` : '',
          fl.charge !== fr.charge ? `Charge ${fl.charge} → ${fr.charge}.` : '',
        ),
      );
    }
    if (fl.mw > 0) card.appendChild(this.row('Atom economy', `${fmt(Math.min(100, (100 * fr.mw) / fl.mw), 1)} %`, undefined, 'Mass of products ÷ mass of reactants (all drawn products counted as desired)'));
    return card;
  }
}

function solubilityClass(logS: number): string {
  if (logS > 0) return 'highly soluble';
  if (logS > -2) return 'soluble';
  if (logS > -4) return 'moderately soluble';
  if (logS > -6) return 'poorly soluble';
  return 'insoluble';
}

function hlbUse(v: number): string {
  if (v < 3) return 'HLB < 3: antifoam / very lipophilic';
  if (v < 7) return 'HLB 3–6: W/O (water-in-oil) emulsifier';
  if (v < 9) return 'HLB 7–9: wetting / spreading agent';
  if (v < 13) return 'HLB 8–16: O/W (oil-in-water) emulsifier';
  if (v < 16) return 'HLB 13–16: detergent / O/W emulsifier';
  return 'HLB 15–18+: solubiliser / hydrotrope';
}
