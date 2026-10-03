// Template library: rings, heterocycles, amino acids, sugars, cosmetic actives… plus user templates.
import type { App } from '../app';
import { h, clear, toast, svgEl } from '../dom';
import { ICONS } from '../icons';
import { TEMPLATE_GROUPS, templateMol, TemplateItem } from '../../chem/templates';
import { Mol } from '../../chem/mol';
import { createDoc, insertMol, docToMol } from '../../doc/document';
import { buildScene } from '../../render/scene';
import { primsToSVG } from '../../render/draw';
import { parseSmiles } from '../../chem/smiles';
import { layoutMol } from '../../chem/layout2d';

const USER_KEY = 'chemwrite:userTemplates';

interface UserTemplate {
  name: string;
  /** serialized atoms/bonds (Mol with coordinates) */
  atoms: Mol['atoms'];
  bonds: Mol['bonds'];
}

const thumbCache = new Map<string, string>();

export function molThumbnail(mol: Mol, key?: string): string {
  if (key && thumbCache.has(key)) return thumbCache.get(key)!;
  const doc = createDoc();
  insertMol(doc, mol);
  const scene = buildScene(doc, { ink: '#222', showErrors: false });
  const b = scene.bounds ?? { x1: 0, y1: 0, x2: 1, y2: 1 };
  const svg = primsToSVG(scene.prims, b, 12, { padding: 0.35 });
  const url = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svg);
  if (key) thumbCache.set(key, url);
  return url;
}

export class LibraryPanel {
  el: HTMLElement;
  private body: HTMLElement;
  private filter: HTMLInputElement;
  private built = false;

  constructor(private app: App) {
    this.filter = h('input', { type: 'search', class: 'input', placeholder: 'Filter templates (e.g. indole, retinol, alanine)…', 'aria-label': 'Filter templates' }) as HTMLInputElement;
    this.filter.addEventListener('input', () => this.applyFilter());
    this.body = h('div', { class: 'library' });
    this.el = h(
      'section',
      { class: 'panel', 'aria-label': 'Template library' },
      h('div', { class: 'panel-body' },
        this.filter,
        h('div', { class: 'row-actions' }, h('button', { class: 'btn btn-small', onclick: () => this.saveSelection() }, svgEl(ICONS.save), 'Save selection as template')),
        h('p', { class: 'muted small' }, 'Click a template, then click the canvas to place it (or click an atom to attach it). Drag to rotate while placing.'),
        this.body,
      ),
    );
  }

  /** Builds thumbnails lazily the first time the panel is shown. */
  activate(): void {
    if (this.built) return;
    this.built = true;
    this.render();
  }

  private render(): void {
    clear(this.body);
    const user = this.loadUser();
    if (user.length) {
      this.body.appendChild(this.groupEl('My templates', user.map((u) => ({ name: u.name, mol: userMol(u), user: true }))));
    }
    for (const g of TEMPLATE_GROUPS) {
      this.body.appendChild(this.groupEl(g.name, g.items.map((it) => ({ name: it.name, item: it }))));
    }
  }

  private groupEl(title: string, items: { name: string; item?: TemplateItem; mol?: Mol; user?: boolean }[]): HTMLElement {
    const grid = h('div', { class: 'tpl-grid' });
    const det = h('details', { class: 'tpl-group', open: title === 'Rings' || title === 'My templates' || title.startsWith('Cosmetic') }, h('summary', null, `${title} (${items.length})`), grid);
    const fill = () => {
      if (grid.childElementCount) return;
      for (const it of items) {
        let mol: Mol;
        try {
          mol = it.mol ?? templateMol(it.item!);
        } catch {
          continue;
        }
        const btn = h('button', { class: 'tpl', title: it.name, 'data-name': it.name.toLowerCase() + ' ' + (it.item?.tags ?? []).join(' ').toLowerCase() },
          h('img', { src: molThumbnail(mol, it.user ? undefined : 'tpl:' + it.name), alt: '', loading: 'lazy' }),
          h('span', null, it.name),
        );
        btn.addEventListener('click', () => {
          this.app.editor.settings.template = { name: it.name, mol };
          this.app.editor.setTool('template');
          this.app.editor.setStatus(`Template “${it.name}”: click to place, click an atom to attach, drag to rotate`);
          this.app.collapseSidebarOnMobile();
        });
        if (it.user) {
          btn.addEventListener('contextmenu', (e) => {
            e.preventDefault();
            if (confirm(`Delete template “${it.name}”?`)) {
              this.saveUser(this.loadUser().filter((u) => u.name !== it.name));
              this.render();
            }
          });
        }
        grid.appendChild(btn);
      }
    };
    if ((det as HTMLDetailsElement).open) fill();
    det.addEventListener('toggle', () => (det as HTMLDetailsElement).open && fill());
    return det;
  }

  private applyFilter(): void {
    const q = this.filter.value.trim().toLowerCase();
    for (const det of this.body.querySelectorAll('details')) {
      const d = det as HTMLDetailsElement;
      if (q) d.open = true;
      let any = false;
      for (const b of d.querySelectorAll('.tpl')) {
        const show = !q || (b as HTMLElement).dataset.name!.includes(q);
        (b as HTMLElement).style.display = show ? '' : 'none';
        any ||= show;
      }
      d.style.display = any || !q ? '' : 'none';
    }
  }

  private loadUser(): UserTemplate[] {
    try {
      return JSON.parse(localStorage.getItem(USER_KEY) ?? '[]');
    } catch {
      return [];
    }
  }

  private saveUser(list: UserTemplate[]): void {
    try {
      localStorage.setItem(USER_KEY, JSON.stringify(list));
    } catch {
      toast('Could not save templates (storage full or blocked)', 'error');
    }
  }

  private saveSelection(): void {
    const ed = this.app.editor;
    if (!ed.hasSelection()) return toast('Select a structure first', 'error');
    const ids = ed.analysisAtomIds();
    const { mol } = docToMol(ed.doc, ids);
    const name = prompt('Template name:', 'My template');
    if (!name) return;
    const list = this.loadUser().filter((u) => u.name !== name);
    list.unshift({ name, atoms: mol.atoms, bonds: mol.bonds });
    this.saveUser(list.slice(0, 60));
    this.render();
    toast(`Saved template “${name}”`, 'success');
  }
}

function userMol(u: UserTemplate): Mol {
  const m = new Mol();
  m.atoms = u.atoms.map((a) => ({ ...a }));
  m.bonds = u.bonds.map((b) => ({ ...b }));
  return m;
}

/** Quick helper used by the command palette: template by name → Mol. */
export function findTemplate(name: string): Mol | null {
  const q = name.trim().toLowerCase();
  for (const g of TEMPLATE_GROUPS) for (const it of g.items) if (it.name.toLowerCase() === q) return templateMol(it);
  return null;
}

export function allTemplates(): { group: string; item: TemplateItem }[] {
  const out: { group: string; item: TemplateItem }[] = [];
  for (const g of TEMPLATE_GROUPS) for (const it of g.items) out.push({ group: g.name, item: it });
  return out;
}

export { parseSmiles, layoutMol };
