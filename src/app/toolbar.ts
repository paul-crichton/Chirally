// Left tool palette with ChemDraw-style flyout groups.
import type { Editor } from '../editor/editor';
import { h, svgEl } from './dom';
import { ICONS, elementIcon } from './icons';
import type { BondStyle } from '../chem/mol';
import type { ArrowKind, ShapeKind } from '../doc/types';

export interface ToolVariant {
  id: string;
  label: string;
  icon: string;
  key?: string;
  apply: (ed: Editor) => void;
  isActive: (ed: Editor) => boolean;
}

export interface ToolGroup {
  id: string;
  label: string;
  variants: ToolVariant[];
  /** index of the variant shown on the button */
  current: number;
  /** extra content for the flyout (e.g. periodic table button) */
  extra?: (close: () => void) => HTMLElement;
  columns?: number;
}

const bondVariant = (id: string, label: string, order: number, style: BondStyle, key?: string): ToolVariant => ({
  id, label, icon: ICONS[id], key,
  apply: (ed) => {
    ed.settings.bond = { order, style };
    ed.setTool('bond');
  },
  isActive: (ed) => ed.toolId === 'bond' && ed.settings.bond.order === order && ed.settings.bond.style === style,
});

const ringVariant = (id: string, label: string, size: number, aromatic: boolean, key?: string): ToolVariant => ({
  id, label, icon: ICONS[id], key,
  apply: (ed) => {
    ed.settings.ring = { size, aromatic };
    ed.setTool('ring');
  },
  isActive: (ed) => ed.toolId === 'ring' && !('chair' in ed.settings.ring) && ed.settings.ring.size === size && ed.settings.ring.aromatic === aromatic,
});

const arrowVariant = (id: string, label: string, kind: ArrowKind, key?: string): ToolVariant => ({
  id, label, icon: ICONS[id], key,
  apply: (ed) => {
    ed.settings.arrow = kind;
    ed.setTool('arrow');
  },
  isActive: (ed) => ed.toolId === 'arrow' && ed.settings.arrow === kind,
});

const shapeVariant = (id: string, label: string, kind: ShapeKind): ToolVariant => ({
  id, label, icon: ICONS[id],
  apply: (ed) => {
    ed.settings.shape = kind;
    ed.setTool('shape');
  },
  isActive: (ed) => ed.toolId === 'shape' && ed.settings.shape === kind,
});

export const QUICK_ELEMENTS = ['C', 'H', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'B', 'Si', 'Na', 'K', 'Li', 'Mg'];
export const QUICK_LABELS = ['Me', 'Et', 'iPr', 'tBu', 'Ph', 'Bn', 'OH', 'OMe', 'NH2', 'NO2', 'CN', 'CF3', 'CO2H', 'CO2Me', 'CHO', 'OAc', 'Ac', 'Boc', 'Ts', 'TBS', 'R', 'X'];

const elementVariant = (el: string): ToolVariant => ({
  id: 'el-' + el, label: `Atom: ${el}`, icon: elementIcon(el),
  apply: (ed) => {
    ed.settings.atom = { el };
    ed.setTool('atom');
  },
  isActive: (ed) => ed.toolId === 'atom' && !ed.settings.atom.label && ed.settings.atom.el === el,
});

const labelVariant = (label: string): ToolVariant => ({
  id: 'lab-' + label, label: `Label: ${label}`, icon: elementIcon(label),
  apply: (ed) => {
    ed.settings.atom = { el: 'C', label };
    ed.setTool('atom');
  },
  isActive: (ed) => ed.toolId === 'atom' && ed.settings.atom.label === label,
});

export function buildToolGroups(openPeriodic: () => void, openLibrary: () => void): ToolGroup[] {
  return [
    {
      id: 'select', label: 'Selection', current: 0,
      variants: [
        { id: 'select', label: 'Select / move (rectangle)', icon: ICONS.select, key: 'V', apply: (ed) => ed.setTool('select'), isActive: (ed) => ed.toolId === 'select' },
        { id: 'lasso', label: 'Lasso select', icon: ICONS.lasso, apply: (ed) => ed.setTool('lasso'), isActive: (ed) => ed.toolId === 'lasso' },
        { id: 'pan', label: 'Pan (or hold Space)', icon: ICONS.pan, apply: (ed) => ed.setTool('pan'), isActive: (ed) => ed.toolId === 'pan' },
      ],
    },
    {
      id: 'erase', label: 'Eraser', current: 0,
      variants: [{ id: 'erase', label: 'Eraser', icon: ICONS.erase, key: 'E', apply: (ed) => ed.setTool('erase'), isActive: (ed) => ed.toolId === 'erase' }],
    },
    {
      id: 'bond', label: 'Bonds', current: 0, columns: 4,
      variants: [
        bondVariant('single', 'Single bond', 1, 'plain', '1'),
        bondVariant('double', 'Double bond', 2, 'plain', '2'),
        bondVariant('triple', 'Triple bond', 3, 'plain', '3'),
        bondVariant('aromatic', 'Aromatic / delocalised bond', 1.5, 'plain'),
        bondVariant('wedge', 'Wedged bond (toward viewer)', 1, 'wedge', 'W'),
        bondVariant('hash', 'Hashed wedge (away from viewer)', 1, 'hash', 'Q'),
        bondVariant('hollow', 'Hollow wedge', 1, 'hollow'),
        bondVariant('wavy', 'Wavy bond (unknown stereo)', 1, 'wavy'),
        bondVariant('bold', 'Bold bond', 1, 'bold'),
        bondVariant('dashed', 'Dashed / partial bond', 1, 'dashed'),
        bondVariant('dative', 'Dative (coordinate) bond', 1, 'dative'),
        bondVariant('hbond', 'Hydrogen bond', 0, 'hbond'),
        bondVariant('crossed', 'Crossed double bond (E/Z unknown)', 2, 'crossed'),
      ],
    },
    {
      id: 'chain', label: 'Chain', current: 0,
      variants: [{ id: 'chain', label: 'Carbon chain', icon: ICONS.chain, key: 'G', apply: (ed) => ed.setTool('chain'), isActive: (ed) => ed.toolId === 'chain' }],
    },
    {
      id: 'ring', label: 'Rings', current: 0, columns: 3,
      variants: [
        ringVariant('benzene', 'Benzene', 6, true, 'R'),
        ringVariant('ring3', 'Cyclopropane', 3, false),
        ringVariant('ring4', 'Cyclobutane', 4, false),
        ringVariant('ring5', 'Cyclopentane', 5, false),
        ringVariant('ring6', 'Cyclohexane', 6, false),
        ringVariant('ring7', 'Cycloheptane', 7, false),
        ringVariant('ring8', 'Cyclooctane', 8, false),
        ringVariant('cyclopentadiene', 'Cyclopentadiene', 5, true),
        {
          id: 'chair', label: 'Cyclohexane chair', icon: ICONS.chair,
          apply: (ed) => {
            ed.settings.ring = { chair: true };
            ed.setTool('ring');
          },
          isActive: (ed) => ed.toolId === 'ring' && 'chair' in ed.settings.ring,
        },
      ],
    },
    {
      id: 'atom', label: 'Atoms & labels', current: 2, columns: 6,
      variants: [...QUICK_ELEMENTS.map(elementVariant), ...QUICK_LABELS.map(labelVariant)],
      extra: (close) =>
        h('button', { class: 'btn btn-small flyout-wide', onclick: () => { close(); openPeriodic(); } }, svgEl(ICONS.periodic), ' Periodic table…'),
    },
    {
      id: 'charge', label: 'Charges & electrons', current: 0,
      variants: [
        { id: 'chargePlus', label: 'Positive charge (+)', icon: ICONS.chargePlus, key: '+', apply: (ed) => { ed.settings.charge = 1; ed.setTool('charge'); }, isActive: (ed) => ed.toolId === 'charge' && ed.settings.charge === 1 },
        { id: 'chargeMinus', label: 'Negative charge (−)', icon: ICONS.chargeMinus, key: '−', apply: (ed) => { ed.settings.charge = -1; ed.setTool('charge'); }, isActive: (ed) => ed.toolId === 'charge' && ed.settings.charge === -1 },
        { id: 'radical', label: 'Radical / carbene', icon: ICONS.radical, apply: (ed) => ed.setTool('radical'), isActive: (ed) => ed.toolId === 'radical' },
        { id: 'lonepair', label: 'Lone pairs', icon: ICONS.lonepair, apply: (ed) => ed.setTool('lonepair'), isActive: (ed) => ed.toolId === 'lonepair' },
      ],
    },
    {
      id: 'curved', label: 'Electron-pushing arrows', current: 0,
      variants: [
        { id: 'curved2', label: 'Curved arrow — electron pair', icon: ICONS.curved2, key: '⇧A', apply: (ed) => { ed.settings.curved = 2; ed.setTool('curved'); }, isActive: (ed) => ed.toolId === 'curved' && ed.settings.curved === 2 },
        { id: 'curved1', label: 'Fishhook — single electron', icon: ICONS.curved1, apply: (ed) => { ed.settings.curved = 1; ed.setTool('curved'); }, isActive: (ed) => ed.toolId === 'curved' && ed.settings.curved === 1 },
      ],
    },
    {
      id: 'arrow', label: 'Reaction arrows', current: 0, columns: 4,
      variants: [
        arrowVariant('reaction', 'Reaction arrow', 'reaction', 'A'),
        arrowVariant('equilibrium', 'Equilibrium', 'equilibrium'),
        arrowVariant('unbalancedEq', 'Unbalanced equilibrium', 'unbalancedEq'),
        arrowVariant('retro', 'Retrosynthetic arrow', 'retro'),
        arrowVariant('resonance', 'Resonance arrow', 'resonance'),
        arrowVariant('dashedArrow', 'Dashed arrow (multi-step / proposed)', 'dashed'),
        arrowVariant('noGo', 'No-go (failed) reaction', 'noGo'),
        arrowVariant('line', 'Line', 'line'),
      ],
    },
    {
      id: 'text', label: 'Text', current: 0,
      variants: [
        { id: 'text', label: 'Text / atom labels', icon: ICONS.text, key: 'T', apply: (ed) => ed.setTool('text'), isActive: (ed) => ed.toolId === 'text' },
        { id: 'plus', label: 'Plus sign', icon: ICONS.plus, apply: (ed) => ed.setTool('plus'), isActive: (ed) => ed.toolId === 'plus' },
      ],
    },
    {
      id: 'shape', label: 'Brackets, shapes & orbitals', current: 0, columns: 3,
      variants: [
        shapeVariant('bracket', 'Brackets [ ]', 'bracket'),
        shapeVariant('tsBracket', 'Transition state [ ]‡', 'tsBracket'),
        shapeVariant('paren', 'Parentheses ( )', 'paren'),
        shapeVariant('brace', 'Braces { }', 'brace'),
        shapeVariant('rect', 'Rectangle', 'rect'),
        shapeVariant('roundRect', 'Rounded rectangle', 'roundRect'),
        shapeVariant('ellipse', 'Ellipse', 'ellipse'),
        shapeVariant('orbitalP', 'p orbital', 'orbitalP'),
        shapeVariant('orbitalS', 's orbital', 'orbitalS'),
      ],
    },
    {
      id: 'template', label: 'Templates library', current: 0,
      variants: [{ id: 'template', label: 'Template library (amino acids, sugars, heterocycles, cosmetic actives…)', icon: ICONS.library, apply: () => openLibrary(), isActive: (ed) => ed.toolId === 'template' }],
    },
  ];
}

export class Toolbar {
  el: HTMLElement;
  private buttons = new Map<string, HTMLButtonElement>();
  private flyout: HTMLElement | null = null;

  constructor(private ed: Editor, private groups: ToolGroup[]) {
    this.el = h('nav', { class: 'toolbar', 'aria-label': 'Drawing tools' });
    for (const g of groups) {
      const btn = h('button', { class: 'tool-btn', 'data-group': g.id }) as HTMLButtonElement;
      btn.addEventListener('click', (e) => {
        if ((e.target as HTMLElement).closest('.tool-more')) {
          this.openFlyout(g, btn);
          return;
        }
        if (g.id === 'template') {
          g.variants[0].apply(ed);
          return;
        }
        const v = g.variants[g.current];
        // a second click on an active multi-variant group opens the flyout (fast access on touch)
        if (v.isActive(ed) && g.variants.length > 1) this.openFlyout(g, btn);
        else v.apply(ed);
      });
      btn.addEventListener('contextmenu', (e) => {
        e.preventDefault();
        if (g.variants.length > 1) this.openFlyout(g, btn);
      });
      let pressTimer: number | null = null;
      btn.addEventListener('pointerdown', (e) => {
        if (g.variants.length < 2 || e.pointerType === 'mouse') return;
        pressTimer = window.setTimeout(() => this.openFlyout(g, btn), 450);
      });
      const cancel = () => {
        if (pressTimer) clearTimeout(pressTimer);
        pressTimer = null;
      };
      btn.addEventListener('pointerup', cancel);
      btn.addEventListener('pointerleave', cancel);
      this.buttons.set(g.id, btn);
      this.el.appendChild(btn);
    }
    ed.on('tool', () => this.refresh());
    document.addEventListener('pointerdown', (e) => {
      if (this.flyout && !this.flyout.contains(e.target as Node) && !(e.target as HTMLElement).closest('.tool-btn')) this.closeFlyout();
    });
    this.refresh();
  }

  refresh(): void {
    for (const g of this.groups) {
      const btn = this.buttons.get(g.id)!;
      const idx = g.variants.findIndex((v) => v.isActive(this.ed));
      if (idx >= 0) g.current = idx;
      const v = g.variants[g.current];
      btn.innerHTML = v.icon + (g.variants.length > 1 ? '<span class="tool-more" aria-hidden="true"></span>' : '');
      btn.classList.toggle('active', idx >= 0);
      const tip = `${v.label}${v.key ? ` (${v.key})` : ''}${g.variants.length > 1 ? ' — right-click or click again for more' : ''}`;
      btn.title = tip;
      btn.setAttribute('aria-label', v.label);
      btn.setAttribute('aria-pressed', String(idx >= 0));
    }
  }

  private openFlyout(g: ToolGroup, anchor: HTMLElement): void {
    this.closeFlyout();
    const close = () => this.closeFlyout();
    const fly = h('div', { class: 'flyout', role: 'menu', 'aria-label': g.label });
    fly.appendChild(h('div', { class: 'flyout-title' }, g.label));
    const grid = h('div', { class: 'flyout-grid', style: `grid-template-columns: repeat(${g.columns ?? Math.min(4, g.variants.length)}, 40px)` });
    g.variants.forEach((v, i) => {
      const b = h('button', { class: 'tool-btn' + (v.isActive(this.ed) ? ' active' : ''), title: v.label + (v.key ? ` (${v.key})` : ''), 'aria-label': v.label, html: v.icon });
      b.addEventListener('click', () => {
        g.current = i;
        v.apply(this.ed);
        this.refresh();
        close();
      });
      grid.appendChild(b);
    });
    fly.appendChild(grid);
    if (g.extra) fly.appendChild(g.extra(close));
    document.body.appendChild(fly);
    const r = anchor.getBoundingClientRect();
    const fr = fly.getBoundingClientRect();
    const vertical = getComputedStyle(this.el).flexDirection === 'column';
    let x = vertical ? r.right + 6 : r.left;
    let y = vertical ? r.top : r.top - fr.height - 6;
    x = Math.max(4, Math.min(x, window.innerWidth - fr.width - 4));
    y = Math.max(4, Math.min(y, window.innerHeight - fr.height - 4));
    fly.style.left = x + 'px';
    fly.style.top = y + 'px';
    this.flyout = fly;
  }

  closeFlyout(): void {
    this.flyout?.remove();
    this.flyout = null;
  }
}
