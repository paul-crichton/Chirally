// Command palette (Ctrl+K): commands, templates, SMILES insertion and PubChem name → structure.
import type { App } from './app';
import { h, clear, debounce } from './dom';
import { looksLikeSmiles } from './chem';
import { allTemplates } from './panels/library';
import * as PC from '../services/pubchem';

export interface Command {
  id: string;
  title: string;
  keys?: string;
  section?: string;
  run: () => void;
}

interface Item {
  title: string;
  subtitle?: string;
  keys?: string;
  run: () => void;
  score: number;
}

function fuzzy(q: string, s: string): number {
  if (!q) return 1;
  const ql = q.toLowerCase(), sl = s.toLowerCase();
  if (sl === ql) return 100;
  if (sl.startsWith(ql)) return 50 - sl.length * 0.01;
  const idx = sl.indexOf(ql);
  if (idx >= 0) return 30 - idx * 0.1;
  // subsequence
  let j = 0, score = 0;
  for (let i = 0; i < sl.length && j < ql.length; i++) if (sl[i] === ql[j]) { j++; score++; }
  return j === ql.length ? 5 + score * 0.1 : 0;
}

export class CommandPalette {
  private overlay: HTMLElement | null = null;
  private input!: HTMLInputElement;
  private list!: HTMLElement;
  private items: Item[] = [];
  private active = 0;
  private remote: Item[] = [];
  private seq = 0;

  constructor(private app: App, private commands: () => Command[]) {}

  open(initial = ''): void {
    if (this.overlay) return;
    this.input = h('input', { class: 'palette-input', placeholder: 'Type a command, a compound name (PubChem), a SMILES string, or a template…', 'aria-label': 'Command palette', spellcheck: false, autocomplete: 'off' }) as HTMLInputElement;
    this.list = h('div', { class: 'palette-list', role: 'listbox' });
    const box = h('div', { class: 'palette', role: 'dialog', 'aria-label': 'Command palette' }, this.input, this.list, h('div', { class: 'palette-foot muted small' }, '↑↓ navigate · Enter run · Esc close'));
    this.overlay = h('div', { class: 'modal-overlay palette-overlay' }, box);
    this.overlay.addEventListener('pointerdown', (e) => {
      if (e.target === this.overlay) this.close();
    });
    document.body.appendChild(this.overlay);
    this.input.value = initial;
    const remoteSearch = debounce(() => this.searchRemote(), 280);
    this.input.addEventListener('input', () => {
      this.remote = [];
      this.refresh();
      remoteSearch();
    });
    this.input.addEventListener('keydown', (e) => {
      e.stopPropagation();
      if (e.key === 'Escape') this.close();
      else if (e.key === 'ArrowDown') {
        this.active = Math.min(this.items.length - 1, this.active + 1);
        this.paint();
        e.preventDefault();
      } else if (e.key === 'ArrowUp') {
        this.active = Math.max(0, this.active - 1);
        this.paint();
        e.preventDefault();
      } else if (e.key === 'Enter') {
        e.preventDefault();
        this.runActive();
      }
    });
    this.refresh();
    if (initial) remoteSearch();
    this.input.focus();
  }

  close(): void {
    this.overlay?.remove();
    this.overlay = null;
    this.app.editor.canvas.focus({ preventScroll: true });
  }

  private runActive(): void {
    const it = this.items[this.active];
    if (!it) return;
    this.close();
    it.run();
  }

  private refresh(): void {
    const q = this.input.value.trim();
    const items: Item[] = [];
    if (q && looksLikeSmiles(q)) {
      items.push({ title: `Insert SMILES: ${q}`, subtitle: 'Parse and lay out this structure', score: 200, run: () => this.app.insertFromText(q, null) });
    }
    if (q) {
      items.push({
        title: `Insert “${q}” from PubChem`,
        subtitle: 'Name, CAS number, InChIKey or CID → structure',
        score: looksLikeSmiles(q) ? 150 : 180,
        run: () => this.app.insertFromPubChem(q),
      });
    }
    for (const c of this.commands()) {
      const s = fuzzy(q, c.title);
      if (s > 0) items.push({ title: c.title, subtitle: c.section, keys: c.keys, score: s + (q ? 0 : 10), run: c.run });
    }
    if (q.length >= 2) {
      for (const t of allTemplates()) {
        const s = Math.max(fuzzy(q, t.item.name), ...(t.item.tags ?? []).map((x) => fuzzy(q, x) * 0.8));
        if (s >= 20) items.push({ title: `Template: ${t.item.name}`, subtitle: t.group, score: s, run: () => this.app.insertTemplateByName(t.item.name) });
      }
    }
    items.push(...this.remote);
    items.sort((a, b) => b.score - a.score);
    this.items = items.slice(0, 40);
    this.active = 0;
    this.paint();
  }

  private async searchRemote(): Promise<void> {
    const q = this.input.value.trim();
    const my = ++this.seq;
    if (q.length < 3 || looksLikeSmiles(q)) return;
    try {
      const terms = await PC.autocomplete(q, 6);
      if (my !== this.seq || !this.overlay) return;
      this.remote = terms
        .filter((t) => t.toLowerCase() !== q.toLowerCase())
        .map((t, i) => ({ title: `PubChem: ${t}`, subtitle: 'Insert structure', score: 60 - i, run: () => this.app.insertFromPubChem(t) }));
      this.refresh();
    } catch {
      /* offline */
    }
  }

  private paint(): void {
    clear(this.list);
    this.items.forEach((it, i) => {
      const row = h('div', { class: 'palette-item' + (i === this.active ? ' active' : ''), role: 'option', 'aria-selected': String(i === this.active) },
        h('div', null, h('div', { class: 'palette-title' }, it.title), it.subtitle ? h('div', { class: 'muted small' }, it.subtitle) : null),
        it.keys ? h('kbd', null, it.keys) : null,
      );
      row.addEventListener('pointerenter', () => {
        this.active = i;
        for (const [k, el] of [...this.list.children].entries()) el.classList.toggle('active', k === i);
      });
      row.addEventListener('click', () => {
        this.active = i;
        this.runActive();
      });
      this.list.appendChild(row);
    });
    const act = this.list.children[this.active] as HTMLElement | undefined;
    act?.scrollIntoView({ block: 'nearest' });
  }
}
