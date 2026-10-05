// Application shell: layout, panels, commands, shortcuts, clipboard, persistence.
import { Editor, LIGHT_THEME, DARK_THEME } from '../editor/editor';
import { h, svgEl, toast, debounce, copyText } from './dom';
import { ICONS } from './icons';
import { Toolbar, buildToolGroups, QUICK_ELEMENTS } from './toolbar';
import { AnalysisPanel } from './panels/analysis';
import { PubChemPanel } from './panels/pubchem';
import { Viewer3DPanel } from './panels/viewer3d';
import { MechanismPanel } from './panels/mechanism';
import { StylePanel } from './panels/style';
import { LibraryPanel, findTemplate, molThumbnail } from './panels/library';
import { CommandPalette, Command } from './palette';
import { showMenu, MenuItem } from './contextmenu';
import { periodicTableDialog, exportDialog, shortcutsDialog, shareDialog, versionsDialog, welcomeDialog, copyImage } from './dialogs';
import { pickFile, saveNative, resetFileHandle, docFromHash, exportDoc, safeName } from './fileio';
import { Analysis, FragmentAnalysis, importText, chemMol, cleanAtoms, expandLabel, smilesToMol, canonicalSmiles, reactionDoc, formats } from './chem';
import { createDoc, docFromJSON, docToJSON, docToMol, insertMol, isDocEmpty, docBounds, cloneDoc, serializeDoc, deserializeDoc } from '../doc/document';
import { ChemDoc, STYLE_PRESETS, DocStyle } from '../doc/types';
import { Mol } from '../chem/mol';
import { Hit } from '../editor/types';
import {
  flipSelection, rotateSelection, drawHydrogens, selectedAtomIds, emptySelection, setElement, applyBondType, deleteSelection, scaleSelection,
  attachRingToAtom, fuseRingOnBond,
} from '../editor/ops';
import { editTextObject, editArrowCaption } from '../editor/tools';
import * as PC from '../services/pubchem';
import { layoutMol } from '../chem/layout2d';
import { parseSmiles } from '../chem/smiles';

const AUTOSAVE_KEY = 'chirally:autosave';
const VERSIONS_KEY = 'chirally:versions';
const PREFS_KEY = 'chirally:prefs';
const SEEN_KEY = 'chirally:welcomed';

type TabId = 'analysis' | 'pubchem' | 'mechanism' | '3d' | 'library' | 'style';

interface Prefs {
  theme?: 'light' | 'dark' | 'auto';
  style?: Partial<DocStyle>;
  grid?: { show: boolean; snap: boolean };
  wheelZooms?: boolean;
  sidebar?: boolean;
  tab?: TabId;
}

export class App {
  editor!: Editor;
  root: HTMLElement;
  analysisPanel!: AnalysisPanel;
  pubchemPanel!: PubChemPanel;
  viewerPanel!: Viewer3DPanel;
  mechanismPanel!: MechanismPanel;
  stylePanel!: StylePanel;
  libraryPanel!: LibraryPanel;
  palette!: CommandPalette;
  private toolbar!: Toolbar;
  private sidebar!: HTMLElement;
  private tabs = new Map<TabId, { btn: HTMLButtonElement; panel: HTMLElement }>();
  private currentTab: TabId = 'analysis';
  private statusEl!: HTMLElement;
  private zoomEl!: HTMLElement;
  private selBar!: HTMLElement;
  private titleEl!: HTMLInputElement;
  private undoBtn!: HTMLButtonElement;
  private redoBtn!: HTMLButtonElement;
  private prefs: Prefs = {};
  private highlight: Map<number, string> | null = null;
  private warningHalos: Map<number, string> | null = null;
  private showLocants = false;
  private lastAnalysis: Analysis | null = null;
  private themeMode: 'light' | 'dark' | 'auto' = 'auto';
  private scheduleAnalysis: ReturnType<typeof debounce>;
  private scheduleAutosave: ReturnType<typeof debounce>;
  private lastVersionTime = 0;
  dirty = false;

  constructor(root: HTMLElement) {
    this.root = root;
    this.loadPrefs();
    this.scheduleAnalysis = debounce(() => this.runAnalysis(), 220);
    this.scheduleAutosave = debounce(() => this.autosave(), 600);
    this.buildLayout();
    this.bindGlobal();
    this.applyTheme();
    this.restore();
  }

  // ───────────── layout ─────────────

  private buildLayout(): void {
    const stage = h('main', { class: 'stage', 'aria-label': 'Canvas' });
    const canvasHost = h('div', { class: 'canvas-host' });
    stage.appendChild(canvasHost);
    this.editor = new Editor(canvasHost);
    const ed = this.editor;

    // top bar
    this.titleEl = h('input', { class: 'doc-title', value: 'Untitled', 'aria-label': 'Document title', spellcheck: false }) as HTMLInputElement;
    this.titleEl.addEventListener('change', () => {
      ed.doc.meta.title = this.titleEl.value.trim() || 'Untitled';
      this.scheduleAutosave();
    });
    this.titleEl.addEventListener('keydown', (e) => {
      if (e.key === 'Enter') (e.target as HTMLInputElement).blur();
    });
    const btn = (icon: string, label: string, run: () => void, cls = 'icon-btn') => {
      const b = h('button', { class: cls, title: label, 'aria-label': label, onclick: run, html: ICONS[icon] }) as HTMLButtonElement;
      return b;
    };
    this.undoBtn = btn('undo', 'Undo (Ctrl+Z)', () => ed.undo());
    this.redoBtn = btn('redo', 'Redo (Ctrl+Shift+Z)', () => ed.redo());
    const fileMenuBtn = h('button', { class: 'btn btn-ghost', 'aria-haspopup': 'menu' }, 'File');
    fileMenuBtn.addEventListener('click', () => {
      const r = fileMenuBtn.getBoundingClientRect();
      showMenu(r.left, r.bottom + 4, this.fileMenu());
    });
    const editMenuBtn = h('button', { class: 'btn btn-ghost', 'aria-haspopup': 'menu' }, 'Edit');
    editMenuBtn.addEventListener('click', () => {
      const r = editMenuBtn.getBoundingClientRect();
      showMenu(r.left, r.bottom + 4, this.editMenu());
    });
    const viewMenuBtn = h('button', { class: 'btn btn-ghost hide-sm', 'aria-haspopup': 'menu' }, 'View');
    viewMenuBtn.addEventListener('click', () => {
      const r = viewMenuBtn.getBoundingClientRect();
      showMenu(r.left, r.bottom + 4, this.viewMenu());
    });
    const search = h('button', { class: 'search-trigger', onclick: () => this.openPalette(''), 'aria-label': 'Search PubChem or run a command (Ctrl+K)' },
      svgEl(ICONS.search), h('span', { class: 'hide-sm' }, 'Search PubChem, SMILES or commands'), h('kbd', { class: 'hide-sm' }, 'Ctrl K'));
    const themeBtn = btn(this.isDark() ? 'sun' : 'moon', 'Toggle dark mode', () => this.toggleTheme());
    themeBtn.classList.add('theme-btn');
    const topbar = h('header', { class: 'topbar' },
      h('div', { class: 'brand', title: 'Chirally' }, h('span', { class: 'logo', html: ICONS.benzene }), h('span', { class: 'brand-name hide-sm' }, 'Chirally')),
      fileMenuBtn, editMenuBtn, viewMenuBtn,
      h('div', { class: 'sep hide-sm' }),
      this.undoBtn, this.redoBtn,
      btn('clean', 'Clean structure (Ctrl+Shift+K)', () => this.clean(), 'icon-btn hide-xs'),
      h('div', { class: 'title-wrap hide-sm' }, this.titleEl),
      h('div', { class: 'spacer' }),
      search,
      btn('save', 'Save (Ctrl+S)', () => this.save(), 'icon-btn hide-xs'),
      btn('export', 'Export (Ctrl+E)', () => exportDialog(this), 'icon-btn'),
      btn('share', 'Share link', () => shareDialog(this), 'icon-btn hide-xs'),
      themeBtn,
      btn('help', 'Keyboard shortcuts & help (?)', () => shortcutsDialog(), 'icon-btn hide-xs'),
      btn('panel', 'Toggle side panel', () => this.toggleSidebar(), 'icon-btn panel-toggle'),
    );

    // toolbar
    this.toolbar = new Toolbar(ed, buildToolGroups(() => this.openPeriodic(), () => this.showTab('library', true)));

    // sidebar
    this.analysisPanel = new AnalysisPanel(this);
    this.pubchemPanel = new PubChemPanel(this);
    this.viewerPanel = new Viewer3DPanel(this);
    this.mechanismPanel = new MechanismPanel(this);
    this.stylePanel = new StylePanel(this);
    this.libraryPanel = new LibraryPanel(this);
    const tabDefs: [TabId, string, string, HTMLElement][] = [
      ['analysis', 'Analysis', 'info', this.analysisPanel.el],
      ['pubchem', 'PubChem', 'pubchem', this.pubchemPanel.el],
      ['mechanism', 'Mechanism', 'mech', this.mechanismPanel.el],
      ['3d', '3D', 'cube', this.viewerPanel.el],
      ['library', 'Library', 'library', this.libraryPanel.el],
      ['style', 'Style', 'style', this.stylePanel.el],
    ];
    const tabBar = h('div', { class: 'tabs', role: 'tablist' });
    const panels = h('div', { class: 'tab-panels' });
    for (const [id, label, icon, panel] of tabDefs) {
      const b = h('button', { class: 'tab', role: 'tab', 'aria-selected': 'false', title: label }, svgEl(ICONS[icon]), h('span', null, label)) as HTMLButtonElement;
      b.addEventListener('click', () => this.showTab(id));
      tabBar.appendChild(b);
      panel.setAttribute('role', 'tabpanel');
      panel.hidden = true;
      panels.appendChild(panel);
      this.tabs.set(id, { btn: b, panel });
    }
    this.sidebar = h('aside', { class: 'sidebar', 'aria-label': 'Panels' },
      h('div', { class: 'sheet-handle', onclick: () => this.toggleSidebar() }),
      tabBar, panels);

    // stage overlays
    this.zoomEl = h('span', { class: 'zoom-label' }, '100%');
    const zoomBox = h('div', { class: 'zoom-box' },
      btn('zoomOut', 'Zoom out', () => ed.zoomBy(0.8)),
      h('button', { class: 'btn btn-ghost zoom-pct', title: 'Reset to 100%', onclick: () => ed.setZoomPercent(100) }, this.zoomEl),
      btn('zoomIn', 'Zoom in', () => ed.zoomBy(1.25)),
      btn('fit', 'Fit to window (Ctrl+0)', () => ed.fitToContent()),
    );
    this.selBar = h('div', { class: 'sel-bar', hidden: true });
    stage.append(zoomBox, this.selBar);
    this.statusEl = h('div', { class: 'statusbar', role: 'status' });

    this.root.append(topbar, this.toolbar.el, stage, this.sidebar, this.statusEl);

    // editor events
    ed.on('change', () => {
      this.dirty = true;
      this.scheduleAnalysis();
      this.updateSelBar();
    });
    ed.on('commit', () => {
      this.scheduleAutosave();
      this.updateUndo();
    });
    ed.on('selection', () => {
      this.scheduleAnalysis();
      this.updateSelBar();
    });
    ed.on('view', () => {
      this.zoomEl.textContent = ed.zoomPercent + '%';
      this.updateSelBar();
    });
    ed.on('status', () => (this.statusEl.textContent = ed.status));
    ed.on('tool', () => this.updateSelBar());
    ed.contextMenuHandler = (sx, sy, hit) => this.contextMenu(sx, sy, hit);
    ed.emptyHint = ed.isTouchDevice
      ? ['Tap to draw a bond', 'Pinch to zoom · long-press for options', 'Search PubChem with the 🔍 button']
      : ['Click to draw a bond — drag to set its angle', 'Hover an atom and type O, N, Cl… to change it', 'Ctrl+K: search PubChem or paste a SMILES · drop a file to open it'];
    this.palette = new CommandPalette(this, () => this.commands());
    this.showTab(this.prefs.tab ?? 'analysis');
    if (this.prefs.sidebar === false || window.innerWidth < 760) this.root.classList.add('sidebar-hidden');
    this.updateUndo();
  }

  showTab(id: TabId, ensureOpen = false): void {
    this.currentTab = id;
    for (const [k, t] of this.tabs) {
      const on = k === id;
      t.btn.classList.toggle('active', on);
      t.btn.setAttribute('aria-selected', String(on));
      t.panel.hidden = !on;
    }
    if (ensureOpen && this.root.classList.contains('sidebar-hidden')) this.toggleSidebar();
    if (id === '3d') this.viewerPanel.activate();
    if (id === 'library') this.libraryPanel.activate();
    if (id === 'mechanism') this.mechanismPanel.update();
    if (id === 'style') this.stylePanel.refresh();
    if (id === 'analysis') this.scheduleAnalysis();
    this.prefs.tab = id;
    this.persistPrefs();
  }

  toggleSidebar(): void {
    this.root.classList.toggle('sidebar-hidden');
    const open = !this.root.classList.contains('sidebar-hidden');
    this.prefs.sidebar = open;
    this.persistPrefs();
    if (open) this.scheduleAnalysis();
    // on phones the panel is a bottom sheet: keep the drawing visible above it
    if (window.innerWidth < 760) {
      const sheet = this.sidebar.getBoundingClientRect().height || window.innerHeight * 0.62;
      this.editor.view.oy += (open ? -1 : 1) * sheet * 0.45;
      this.editor.requestRender();
    }
    setTimeout(() => this.editor.requestRender(), 250);
  }

  /** Height (px) of the canvas hidden behind the panel when it is a bottom sheet (phones); 0 otherwise. */
  canvasCoveredBottom(): number {
    if (window.innerWidth >= 760 || this.root.classList.contains('sidebar-hidden')) return 0;
    const c = this.editor.canvas.getBoundingClientRect(), s = this.sidebar.getBoundingClientRect();
    return s.top < c.bottom && s.left <= c.left + 8 ? c.bottom - s.top : 0;
  }

  collapseSidebarOnMobile(): void {
    if (window.innerWidth < 760) this.root.classList.add('sidebar-hidden');
  }

  private updateUndo(): void {
    const hs = this.editor.history;
    this.undoBtn.disabled = !hs.canUndo;
    this.redoBtn.disabled = !hs.canRedo;
    this.undoBtn.title = hs.undoLabel ? `Undo ${hs.undoLabel} (Ctrl+Z)` : 'Undo (Ctrl+Z)';
    this.redoBtn.title = hs.redoLabel ? `Redo ${hs.redoLabel} (Ctrl+Shift+Z)` : 'Redo (Ctrl+Shift+Z)';
  }

  // ───────────── floating selection bar ─────────────

  private updateSelBar(): void {
    const ed = this.editor;
    const box = ed.selectionScreenBox();
    const show = !!box && (ed.toolId === 'select' || ed.toolId === 'lasso') && !ed.isDragging && !ed.isEditingInline;
    if (!show) {
      this.selBar.hidden = true;
      return;
    }
    const atoms = selectedAtomIds(ed.doc, ed.sel);
    const singleAtom = ed.sel.atoms.size === 1 && ed.sel.bonds.size === 0 && ed.sel.objects.size === 0 ? [...ed.sel.atoms][0] : null;
    const singleBond = ed.sel.bonds.size === 1 && ed.sel.atoms.size <= 2 && ed.sel.objects.size === 0 ? [...ed.sel.bonds][0] : null;
    const key = `${singleAtom}|${singleBond}|${atoms.size}|${ed.sel.objects.size}`;
    if (this.selBar.dataset.key !== key) {
      this.selBar.dataset.key = key;
      this.selBar.innerHTML = '';
      const b = (icon: string, label: string, run: () => void) => h('button', { class: 'icon-btn', title: label, 'aria-label': label, html: ICONS[icon], onclick: run });
      if (singleAtom !== null) {
        for (const el of ['C', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'P', 'H']) {
          this.selBar.appendChild(h('button', { class: 'el-btn', title: `Change to ${el}`, onclick: () => ed.mutate('Change element', (d) => setElement(d, singleAtom, el)) }, el));
        }
        this.selBar.appendChild(h('button', { class: 'el-btn', title: 'Charge +', onclick: () => ed.mutate('Charge', (d) => (d.atoms.get(singleAtom)!.charge += 1)) }, '+'));
        this.selBar.appendChild(h('button', { class: 'el-btn', title: 'Charge −', onclick: () => ed.mutate('Charge', (d) => (d.atoms.get(singleAtom)!.charge -= 1)) }, '−'));
        this.selBar.appendChild(b('text', 'Edit label', () => ed.editAtomLabel(singleAtom)));
        if (ed.doc.atoms.get(singleAtom)?.abbrev) this.selBar.appendChild(b('expand', 'Expand label', () => this.expandSelectedLabels()));
      } else if (singleBond !== null) {
        const set = (o: number, s: 'plain' | 'wedge' | 'hash' | 'wavy' | 'bold') => () => ed.mutate('Change bond', (d) => applyBondType(d.bonds.get(singleBond)!, o, s));
        this.selBar.append(b('single', 'Single', set(1, 'plain')), b('double', 'Double', set(2, 'plain')), b('triple', 'Triple', set(3, 'plain')), b('wedge', 'Wedge', set(1, 'wedge')), b('hash', 'Hash', set(1, 'hash')), b('wavy', 'Wavy', set(1, 'wavy')));
      }
      if (atoms.size > 1 || ed.sel.objects.size) {
        this.selBar.append(
          b('clean', 'Clean structure', () => this.clean()),
          b('flipH', 'Flip horizontal', () => ed.mutate('Flip', (d) => flipSelection(d, ed.sel, 'h'))),
          b('flipV', 'Flip vertical', () => ed.mutate('Flip', (d) => flipSelection(d, ed.sel, 'v'))),
          b('rotate', 'Rotate 30°', () => this.rotateSel(Math.PI / 6)),
        );
      }
      this.selBar.append(
        b('duplicate', 'Duplicate (Ctrl+D)', () => this.duplicate()),
        b('copy', 'Copy (Ctrl+C)', () => this.copy()),
        atoms.size > 1 ? b('cube', 'View in 3D', () => { this.showTab('3d', true); this.viewerPanel.generate(); }) : null as unknown as HTMLElement,
        b('trash', 'Delete', () => ed.deleteSelected()),
      );
      [...this.selBar.childNodes].forEach((n) => { if (!n) this.selBar.removeChild(n as Node); });
    }
    this.selBar.hidden = false;
    const stage = this.selBar.parentElement!;
    const bw = this.selBar.offsetWidth || 300;
    let x = (box!.x1 + box!.x2) / 2 - bw / 2;
    // keep clear of the rotate handle drawn just above the selection box
    let y = box!.y1 - (ed.rotateHandle() ? 82 : 54);
    if (y < 6) y = box!.y2 + 10;
    x = Math.max(6, Math.min(x, stage.clientWidth - bw - 6));
    y = Math.max(6, Math.min(y, stage.clientHeight - 50));
    this.selBar.style.transform = `translate(${x}px, ${y}px)`;
  }

  private rotateSel(ang: number): void {
    const ed = this.editor;
    const c = ed.selectionCenter();
    if (c) ed.mutate('Rotate', (d) => rotateSelection(d, ed.sel, c, ang));
  }

  // ───────────── analysis & annotations ─────────────

  private runAnalysis(): void {
    const sidebarOpen = !this.root.classList.contains('sidebar-hidden');
    const needed = (this.currentTab === 'analysis' && sidebarOpen) || this.editor.doc.style.showStereoLabels || this.showLocants;
    if (needed) this.analysisPanel.update();
    if (this.currentTab === 'mechanism') this.mechanismPanel.update();
    this.updateTitle();
  }

  onAnalysis(a: Analysis | null): void {
    this.lastAnalysis = a;
    this.pushAnnotations();
  }

  private pushAnnotations(): void {
    const a = this.lastAnalysis;
    const halo = new Map<number, string>();
    if (this.warningHalos) for (const [k, v] of this.warningHalos) halo.set(k, v);
    if (this.highlight) for (const [k, v] of this.highlight) halo.set(k, v);
    this.editor.setAnnotations({
      stereo: a?.stereo,
      bondLabels: a?.bondStereo,
      notes: this.showLocants ? a?.locants : undefined,
      halo: halo.size ? halo : undefined,
    });
  }

  highlightAtoms(ids: number[] | null, color = 'rgba(250,176,5,0.35)'): void {
    this.highlight = ids ? new Map(ids.map((id) => [id, color])) : null;
    this.pushAnnotations();
  }

  highlightFunctionalGroup(f: FragmentAnalysis, groups: number[][]): void {
    const ids = new Set<number>();
    for (const g of groups) for (const i of g) if (f.idMap[i] !== undefined) ids.add(f.idMap[i]);
    this.highlightAtoms([...ids], 'rgba(18,184,134,0.32)');
  }

  setWarningHalos(m: Map<number, string> | null): void {
    this.warningHalos = m;
    this.pushAnnotations();
  }

  toggleLocants(): void {
    this.showLocants = !this.showLocants;
    this.pushAnnotations();
    if (this.showLocants && !this.lastAnalysis?.locants.size) toast('No locants available for this structure', 'info');
  }

  currentSmiles(selectionOnly = false): string {
    const ids = selectionOnly || this.editor.hasSelection() ? this.editor.analysisAtomIds() : [...this.editor.doc.atoms.keys()];
    if (!ids.length) return '';
    try {
      return canonicalSmiles(chemMol(this.editor.doc, ids));
    } catch {
      return '';
    }
  }

  /** Mol of the selection (or everything) for 3D: abbreviations expanded, stereo perceived. */
  chemMolForSelection(): Mol | null {
    const ids = this.editor.analysisAtomIds();
    if (!ids.length) return null;
    return chemMol(this.editor.doc, ids);
  }

  /** A standalone document containing just the selection. */
  selectionDoc(): ChemDoc {
    const ed = this.editor;
    if (!ed.hasSelection()) return ed.doc;
    return deserializeDoc(JSON.parse(ed.selectionToJSON()));
  }

  identify(smiles: string): void {
    this.showTab('pubchem', true);
    this.pubchemPanel.identify(smiles);
  }

  openPalette(q = ''): void {
    this.palette.open(q);
  }

  openPeriodic(): void {
    periodicTableDialog((el) => {
      const ed = this.editor;
      const atoms = [...selectedAtomIds(ed.doc, ed.sel)];
      if (atoms.length && ed.toolId === 'select') {
        ed.mutate('Change element', (d) => atoms.forEach((id) => setElement(d, id, el)));
      } else {
        ed.settings.atom = { el };
        ed.setTool('atom');
      }
    });
  }

  // ───────────── insertion / import ─────────────

  /** Inserts structures from any supported text (SMILES, MOL, SDF, CDXML, …). */
  insertFromText(text: string, fileName: string | null, name?: string, opts: { suppressH?: boolean } = {}): void {
    try {
      const r = importText(text, fileName, opts);
      if (r.kind === 'chirally') {
        this.mergeDoc(docFromJSON(text));
        return;
      }
      if (r.doc) {
        this.mergeDoc(r.doc);
        toast(`Imported ${r.kind}`, 'success');
        return;
      }
      const mols = r.mols ?? [];
      if (!mols.length) throw new Error('No structures found');
      // lay several molecules out in rows (grid) so large imports stay readable
      const combined = new Mol();
      const boxes = mols.map((m) => m.bbox());
      const area = boxes.reduce((s, b) => s + (b.maxX - b.minX + 2) * (b.maxY - b.minY + 2), 0);
      const rowWidth = Math.max(14, Math.sqrt(area) * 1.6);
      let x = 0, y = 0, rowH = 0;
      mols.forEach((m, i) => {
        const bb = boxes[i];
        const w = bb.maxX - bb.minX, hgt = bb.maxY - bb.minY;
        if (x > 0 && x + w > rowWidth) {
          x = 0;
          y += rowH + 2.2;
          rowH = 0;
        }
        m.translate(x - bb.minX, y - bb.minY);
        x += w + 2;
        rowH = Math.max(rowH, hgt);
        combined.append(m);
      });
      this.editor.insertMolecule(combined, undefined, name ? `Insert ${name}` : 'Insert structure');
      this.ensureVisible();
      if (name && isDocEmptyExcept(this.editor.doc, combined.atoms.length) && this.editor.doc.meta.title === 'Untitled') {
        this.editor.doc.meta.title = name;
        this.updateTitle();
      }
    } catch (e) {
      toast((e as Error).message || 'Could not read that structure', 'error', 4000);
    }
  }

  async insertFromPubChem(q: string): Promise<void> {
    toast(`Searching PubChem for “${q}”…`);
    try {
      const list = await PC.search(q, 1);
      if (!list.length) return toast(`Nothing found in PubChem for “${q}”`, 'error');
      await this.pubchemPanel.insert(list[0]);
    } catch (e) {
      toast((e as Error).message, 'error');
    }
  }

  insertTemplateByName(name: string): void {
    const m = findTemplate(name);
    if (!m) return;
    this.editor.insertMolecule(m.clone(), undefined, `Insert ${name}`);
    this.ensureVisible();
  }

  /** Adds another document's content to the right of the current drawing. */
  mergeDoc(other: ChemDoc): void {
    const ed = this.editor;
    if (isDocEmpty(ed.doc)) {
      const keepStyle = ed.doc.style;
      ed.setDoc(other);
      if (!other.style) ed.doc.style = keepStyle;
      this.updateTitle();
      setTimeout(() => ed.fitToContent(), 0);
      return;
    }
    const json = JSON.stringify(serializeDoc(other));
    const b = docBounds(other);
    const cur = docBounds(ed.doc)!;
    const w = b ? b.maxX - b.minX : 0;
    ed.pasteDocJSON(json, { x: cur.maxX + 2 + w / 2, y: (cur.minY + cur.maxY) / 2 });
    this.ensureVisible();
  }

  private ensureVisible(): void {
    const ed = this.editor;
    const box = ed.selectionScreenBox();
    if (!box) return;
    const W = ed.canvas.clientWidth, H = ed.canvas.clientHeight;
    if (box.x1 < 0 || box.y1 < 0 || box.x2 > W || box.y2 > H) ed.fitToContent();
  }

  async openFile(file?: File | null): Promise<void> {
    const f = file ?? (await pickFile());
    if (!f) return;
    const text = await f.text();
    const isNative = /\.(chirally|json)$/i.test(f.name) || text.trimStart().startsWith('{');
    if (isNative) {
      try {
        const doc = docFromJSON(text);
        if (!isDocEmpty(this.editor.doc) && !confirm('Replace the current drawing? (Cancel adds it to the drawing instead)')) {
          this.mergeDoc(doc);
          return;
        }
        resetFileHandle();
        this.editor.setDoc(doc);
        this.updateTitle();
        setTimeout(() => this.editor.fitToContent(), 0);
        toast(`Opened ${f.name}`, 'success');
      } catch (e) {
        toast('Could not open file: ' + (e as Error).message, 'error');
      }
      return;
    }
    this.insertFromText(text, f.name, f.name.replace(/\.[^.]+$/, ''));
  }

  // ───────────── commands ─────────────

  newDoc(): void {
    if (!isDocEmpty(this.editor.doc) && !confirm('Start a new drawing? The current one stays in Version history.')) return;
    this.recordVersion(true);
    resetFileHandle();
    const style = { ...this.editor.doc.style };
    const d = createDoc(style);
    this.editor.setDoc(d);
    this.editor.view = { scale: 42, ox: this.editor.canvas.clientWidth / 2, oy: this.editor.canvas.clientHeight / 2 };
    this.editor.requestRender();
    this.updateTitle();
  }

  async save(saveAs = false): Promise<void> {
    const name = await saveNative(this.editor.doc, saveAs);
    if (name) {
      toast(`Saved ${name}`, 'success');
      this.dirty = false;
    }
  }

  clean(): void {
    const ed = this.editor;
    const ids = ed.hasSelection() ? [...selectedAtomIds(ed.doc, ed.sel)] : [...ed.doc.atoms.keys()];
    if (!ids.length) return;
    try {
      ed.mutate('Clean structure', (d) => cleanAtoms(d, ids));
    } catch (e) {
      ed.cancelChange();
      toast('Clean failed: ' + (e as Error).message, 'error');
    }
  }

  expandSelectedLabels(): void {
    const ed = this.editor;
    const ids = [...selectedAtomIds(ed.doc, ed.sel)].filter((id) => ed.doc.atoms.get(id)?.abbrev);
    const targets = ids.length ? ids : ed.hover?.kind === 'atom' && ed.doc.atoms.get(ed.hover.id)?.abbrev ? [ed.hover.id] : [];
    if (!targets.length) return toast('Select an abbreviation label (e.g. OMe, Boc) to expand', 'info');
    ed.mutate('Expand label', (d) => targets.forEach((id) => expandLabel(d, id)));
  }

  /** Draws implicit H as H atoms on the given atoms, the selection, or the hovered atom (for curved arrows on C–H). */
  drawHydrogenAtoms(ids?: number[]): void {
    const ed = this.editor;
    const sel = [...selectedAtomIds(ed.doc, ed.sel)];
    const targets = ids ?? (sel.length ? sel : ed.hover?.kind === 'atom' ? [ed.hover.id] : []);
    if (!targets.length) return toast('Select the atoms whose hydrogens you want to draw', 'info');
    let n = 0;
    ed.mutate('Draw hydrogens', (d) => (n = drawHydrogens(d, targets)));
    if (!n) toast('Those atoms have no implicit hydrogens', 'info');
  }

  duplicate(): void {
    const ed = this.editor;
    if (!ed.hasSelection()) return;
    const json = ed.selectionToJSON();
    const c = ed.selectionCenter()!;
    const b = ed.selectionScreenBox()!;
    const w = (b.x2 - b.x1) / ed.view.scale;
    ed.pasteDocJSON(json, { x: c.x + w + 0.5, y: c.y });
  }

  copy(cut = false): void {
    const ed = this.editor;
    if (!ed.doc.atoms.size && !ed.doc.arrows.size && !ed.doc.texts.size) return;
    const json = ed.selectionToJSON();
    const smi = this.currentSmiles();
    ed.setInternalClipboard(json, smi);
    const CI = (window as any).ClipboardItem;
    const plain = smi || ' ';
    if (navigator.clipboard && CI) {
      const items: Record<string, Blob> = { 'text/plain': new Blob([plain], { type: 'text/plain' }) };
      navigator.clipboard.write([new CI(items)]).catch(() => navigator.clipboard.writeText(plain).catch(() => undefined));
    } else navigator.clipboard?.writeText(plain).catch(() => undefined);
    if (cut) ed.deleteSelected();
    toast(cut ? 'Cut' : smi ? 'Copied (paste into Chirally, or as SMILES elsewhere)' : 'Copied', 'success', 1600);
  }

  async paste(text?: string): Promise<void> {
    const ed = this.editor;
    const internal = ed.getInternalClipboard();
    let t = text;
    if (t === undefined) {
      try {
        t = await navigator.clipboard.readText();
      } catch {
        t = undefined;
      }
    }
    const at = ed.mouseInside ? ed.mouse : undefined;
    if (internal && (t === undefined || t.trim() === internal.text.trim() || !t.trim())) {
      ed.pasteDocJSON(internal.json, at);
      return;
    }
    if (t && t.trim()) this.insertFromText(t, null);
  }

  loadExample(): void {
    const ed = this.editor;
    try {
      const doc = createDoc(ed.doc.style);
      // 1) aspirin synthesis scheme
      const sal = smilesToMol('OC(=O)c1ccccc1O');
      const anh = smilesToMol('CC(=O)OC(C)=O');
      const asp = smilesToMol('CC(=O)Oc1ccccc1C(=O)O');
      const scheme = reactionDoc([sal, anh], [asp]);
      const arrow = [...scheme.arrows.values()][0];
      arrow.above = 'H_2SO_4 (cat.)';
      arrow.below = '85 °C, 15 min';
      const sj = serializeDoc(scheme);
      const s2 = deserializeDoc(sj);
      for (const [k, v] of s2.atoms) doc.atoms.set(k, v);
      for (const [k, v] of s2.bonds) doc.bonds.set(k, v);
      for (const [k, v] of s2.arrows) doc.arrows.set(k, v);
      for (const [k, v] of s2.texts) doc.texts.set(k, v);
      doc.nextId = s2.nextId;
      const title = doc.nextId++;
      doc.texts.set(title, { id: title, type: 'text', x: -0.5, y: -3.2, text: '**Synthesis of aspirin** (try: select a product → Analysis tab)', size: 0.9 });
      // 2) SN2 mechanism with curved arrows (apply them from the Mechanism tab)
      const oh = parseSmiles('[OH-]');
      oh.atoms[0].x = 0;
      oh.atoms[0].y = 0;
      oh.atoms[0].lonePairs = true;
      const ohIds = insertMol(doc, oh, 0, 5.5).atomIds;
      const ebr = smilesToMol('CCBr');
      const bb = ebr.bbox();
      const ebIds = insertMol(doc, ebr, 2.2 - bb.minX, 5.5 - (bb.minY + bb.maxY) / 2).atomIds;
      // find C attached to Br
      const brIdx = ebr.atoms.findIndex((a) => a.el === 'Br');
      const cIdx = ebr.neighbors(brIdx)[0];
      const cId = ebIds[cIdx], brId = ebIds[brIdx];
      doc.atoms.get(brId)!.lonePairs = true;
      let cbBond = 0;
      for (const b of doc.bonds.values()) if ((b.a === cId && b.b === brId) || (b.a === brId && b.b === cId)) cbBond = b.id;
      const a1 = doc.nextId++;
      doc.curved.set(a1, { id: a1, type: 'curved', electrons: 2, from: { type: 'atom', id: ohIds[0] }, to: { type: 'atom', id: cId }, c1: { t: 0.25, h: -0.55 }, c2: { t: 0.75, h: -0.55 } });
      const a2 = doc.nextId++;
      doc.curved.set(a2, { id: a2, type: 'curved', electrons: 2, from: { type: 'bond', id: cbBond }, to: { type: 'atom', id: brId }, c1: { t: 0.3, h: 0.55 }, c2: { t: 0.95, h: 0.55 } });
      const t2 = doc.nextId++;
      doc.texts.set(t2, { id: t2, type: 'text', x: -0.5, y: 3.6, text: '**S_N2 mechanism** — open the Mechanism tab and press “Apply arrows”', size: 0.9 });
      doc.meta.title = 'Chirally examples';
      if (!isDocEmpty(ed.doc)) this.recordVersion(true);
      ed.setDoc(doc);
      this.updateTitle();
      setTimeout(() => ed.fitToContent(), 0);
    } catch (e) {
      toast('Could not build the example: ' + (e as Error).message, 'error');
    }
  }

  commands(): Command[] {
    const ed = this.editor;
    const C = (id: string, title: string, run: () => void, keys?: string, section = 'Command'): Command => ({ id, title, run, keys, section });
    return [
      C('new', 'New drawing', () => this.newDoc(), undefined, 'File'),
      C('open', 'Open file… (Chirally, MOL, SDF, CDXML, RXN, CML, XYZ, SMILES)', () => this.openFile(), 'Ctrl O', 'File'),
      C('save', 'Save', () => this.save(), 'Ctrl S', 'File'),
      C('saveas', 'Save as…', () => this.save(true), 'Ctrl Shift S', 'File'),
      C('export', 'Export…', () => exportDialog(this), 'Ctrl E', 'File'),
      C('export-svg', 'Export SVG', () => exportDoc(ed.doc, 'svg'), undefined, 'File'),
      C('export-png', 'Export PNG', () => exportDoc(ed.doc, 'png'), undefined, 'File'),
      C('export-cdxml', 'Export ChemDraw CDXML', () => exportDoc(ed.doc, 'cdxml'), undefined, 'File'),
      C('export-mol', 'Export MOL file', () => exportDoc(ed.doc, 'mol'), undefined, 'File'),
      C('copy-image', 'Copy as image', () => copyImage(this), undefined, 'Edit'),
      C('copy-smiles', 'Copy SMILES', () => copyText(this.currentSmiles(), 'SMILES'), undefined, 'Edit'),
      C('copy-mol', 'Copy MOL block', () => copyText(formats.writeMolfile(docToMol(ed.doc, ed.analysisAtomIds()).mol), 'MOL block'), undefined, 'Edit'),
      C('share', 'Share link', () => shareDialog(this), undefined, 'File'),
      C('versions', 'Version history', () => versionsDialog(this), undefined, 'File'),
      C('undo', 'Undo', () => ed.undo(), 'Ctrl Z', 'Edit'),
      C('redo', 'Redo', () => ed.redo(), 'Ctrl Shift Z', 'Edit'),
      C('selectall', 'Select all', () => ed.selectAll(), 'Ctrl A', 'Edit'),
      C('duplicate', 'Duplicate selection', () => this.duplicate(), 'Ctrl D', 'Edit'),
      C('delete', 'Delete selection', () => ed.deleteSelected(), 'Del', 'Edit'),
      C('clean', 'Clean structure', () => this.clean(), 'Ctrl Shift K', 'Structure'),
      C('expand', 'Expand abbreviation label', () => this.expandSelectedLabels(), undefined, 'Structure'),
      C('draw-h', 'Draw hydrogens as atoms (for arrows on C–H bonds)', () => this.drawHydrogenAtoms(), undefined, 'Structure'),
      C('fliph', 'Flip horizontal', () => ed.mutate('Flip', (d) => flipSelection(d, ed.sel, 'h')), undefined, 'Structure'),
      C('flipv', 'Flip vertical', () => ed.mutate('Flip', (d) => flipSelection(d, ed.sel, 'v')), undefined, 'Structure'),
      C('mirror', 'Mirror image (enantiomer)', () => ed.mutate('Mirror', (d) => flipSelection(d, ed.sel.atoms.size ? ed.sel : allSel(ed.doc), 'h', false)), undefined, 'Structure'),
      C('rot90', 'Rotate 90° clockwise', () => this.rotateSel(Math.PI / 2), undefined, 'Structure'),
      C('scaleup', 'Scale selection up 10%', () => ed.mutate('Scale', (d) => scaleSelection(d, ed.sel, 1.1)), undefined, 'Structure'),
      C('scaledown', 'Scale selection down 10%', () => ed.mutate('Scale', (d) => scaleSelection(d, ed.sel, 1 / 1.1)), undefined, 'Structure'),
      C('3d', 'Generate 3D model', () => { this.showTab('3d', true); this.viewerPanel.generate(); }, undefined, 'Structure'),
      C('identify', 'Identify structure in PubChem', () => this.identify(this.currentSmiles()), undefined, 'PubChem'),
      C('similar', 'Find similar compounds in PubChem', () => { this.showTab('pubchem', true); this.pubchemPanel.similarCurrent(); }, undefined, 'PubChem'),
      C('mech-apply', 'Apply electron-pushing arrows → next intermediate', () => this.mechanismPanel.applyStep(), undefined, 'Mechanism'),
      C('mech-check', 'Check electron counts (octet rule)', () => { this.showTab('mechanism', true); this.mechanismPanel.checkElectrons(); }, undefined, 'Mechanism'),
      C('lonepairs', 'Toggle all lone pairs', () => this.mechanismPanel.toggleLonePairs(), undefined, 'View'),
      C('stereo', 'Toggle R/S and E/Z labels', () => { ed.doc.style.showStereoLabels = !ed.doc.style.showStereoLabels; ed.touch(); this.scheduleAnalysis(); }, undefined, 'View'),
      C('locants', 'Toggle IUPAC locant numbering', () => this.toggleLocants(), undefined, 'View'),
      C('periodic', 'Periodic table', () => this.openPeriodic(), undefined, 'Tools'),
      C('library', 'Template library', () => this.showTab('library', true), undefined, 'Tools'),
      C('fit', 'Zoom to fit', () => ed.fitToContent(), 'Ctrl 0', 'View'),
      C('zoom100', 'Zoom 100%', () => ed.setZoomPercent(100), 'Ctrl 1', 'View'),
      C('grid', 'Toggle grid', () => { ed.grid.show = !ed.grid.show; ed.requestRender(); this.persistPrefs(); }, undefined, 'View'),
      C('theme', 'Toggle dark mode', () => this.toggleTheme(), undefined, 'View'),
      C('panel', 'Toggle side panel', () => this.toggleSidebar(), undefined, 'View'),
      C('shortcuts', 'Keyboard shortcuts', () => shortcutsDialog(), '?', 'Help'),
      C('welcome', 'Welcome / quick tour', () => welcomeDialog(this), undefined, 'Help'),
      C('example', 'Load example drawing', () => this.loadExample(), undefined, 'Help'),
      ...Object.keys(STYLE_PRESETS).map((k) => C('style-' + k, `Apply style: ${k}`, () => { ed.mutate('Style', (d) => Object.assign(d.style, STYLE_PRESETS[k])); this.stylePanel.refresh(); }, undefined, 'Style')),
    ];
  }

  private fileMenu(): MenuItem[] {
    return [
      { label: 'New', run: () => this.newDoc() },
      { label: 'Open…', keys: 'Ctrl O', run: () => this.openFile() },
      { label: 'Save', keys: 'Ctrl S', run: () => this.save() },
      { label: 'Save as…', keys: 'Ctrl ⇧ S', run: () => this.save(true) },
      { separator: true, label: '' },
      { label: 'Export…', keys: 'Ctrl E', run: () => exportDialog(this) },
      { label: 'Copy as image', run: () => copyImage(this) },
      { label: 'Share link…', run: () => shareDialog(this) },
      { separator: true, label: '' },
      { label: 'Version history…', run: () => versionsDialog(this) },
      { label: 'Load example', run: () => this.loadExample() },
    ];
  }

  private editMenu(): MenuItem[] {
    const ed = this.editor;
    return [
      { label: 'Undo', keys: 'Ctrl Z', run: () => ed.undo(), disabled: !ed.history.canUndo },
      { label: 'Redo', keys: 'Ctrl ⇧ Z', run: () => ed.redo(), disabled: !ed.history.canRedo },
      { separator: true, label: '' },
      { label: 'Cut', keys: 'Ctrl X', run: () => this.copy(true) },
      { label: 'Copy', keys: 'Ctrl C', run: () => this.copy() },
      { label: 'Paste', keys: 'Ctrl V', run: () => this.paste() },
      { label: 'Duplicate', keys: 'Ctrl D', run: () => this.duplicate() },
      { label: 'Select all', keys: 'Ctrl A', run: () => ed.selectAll() },
      { separator: true, label: '' },
      { label: 'Clean structure', keys: 'Ctrl ⇧ K', run: () => this.clean() },
      { label: 'Expand label', run: () => this.expandSelectedLabels() },
      { label: 'Draw hydrogens as atoms', run: () => this.drawHydrogenAtoms() },
      { label: 'Flip horizontal', run: () => ed.mutate('Flip', (d) => flipSelection(d, ed.sel.atoms.size ? ed.sel : allSel(d), 'h')) },
      { label: 'Flip vertical', run: () => ed.mutate('Flip', (d) => flipSelection(d, ed.sel.atoms.size ? ed.sel : allSel(d), 'v')) },
      { label: 'Command palette…', keys: 'Ctrl K', run: () => this.openPalette() },
    ];
  }

  private viewMenu(): MenuItem[] {
    const ed = this.editor;
    const st = ed.doc.style;
    const tog = (label: string, on: boolean, run: () => void): MenuItem => ({ label: (on ? '✓ ' : '   ') + label, run });
    return [
      { label: 'Zoom to fit', keys: 'Ctrl 0', run: () => ed.fitToContent() },
      { label: 'Zoom 100%', keys: 'Ctrl 1', run: () => ed.setZoomPercent(100) },
      { separator: true, label: '' },
      tog('R/S & E/Z labels', st.showStereoLabels, () => { st.showStereoLabels = !st.showStereoLabels; ed.touch(); this.scheduleAnalysis(); }),
      tog('Lone pairs', st.showLonePairs, () => this.mechanismPanel.toggleLonePairs()),
      tog('IUPAC locants', this.showLocants, () => this.toggleLocants()),
      tog('Colour atoms', st.colorAtoms, () => { st.colorAtoms = !st.colorAtoms; ed.touch(); }),
      tog('Grid', ed.grid.show, () => { ed.grid.show = !ed.grid.show; ed.requestRender(); this.persistPrefs(); }),
      tog('Dark mode', this.isDark(), () => this.toggleTheme()),
      tog('Side panel', !this.root.classList.contains('sidebar-hidden'), () => this.toggleSidebar()),
    ];
  }

  private contextMenu(sx: number, sy: number, hit: Hit): void {
    const ed = this.editor;
    const r = ed.canvas.getBoundingClientRect();
    const x = r.left + sx, y = r.top + sy;
    const items: MenuItem[] = [];
    if (hit?.kind === 'atom') {
      const id = hit.id;
      const a = ed.doc.atoms.get(id)!;
      items.push(
        { label: 'Edit label…', keys: 'Enter', run: () => ed.editAtomLabel(id) },
        { label: 'Element', submenu: [...QUICK_ELEMENTS.slice(0, 12).map((el) => ({ label: el, run: () => ed.mutate('Change element', (d) => setElement(d, id, el)) })), { label: 'Periodic table…', run: () => { ed.setSelection({ atoms: new Set([id]), bonds: new Set(), objects: new Set() }); ed.setTool('select'); this.openPeriodic(); } }] },
        { label: 'Charge', submenu: [-2, -1, 0, 1, 2].map((c) => ({ label: c > 0 ? `+${c}` : String(c), run: () => ed.mutate('Charge', (d) => (d.atoms.get(id)!.charge = c)) })) },
        { label: a.radical ? 'Remove radical' : 'Add radical', run: () => ed.mutate('Radical', (d) => { const t = d.atoms.get(id)!; if (t.radical) delete t.radical; else t.radical = 1; }) },
        { label: a.lonePairs ? 'Hide lone pairs' : 'Show lone pairs', run: () => ed.mutate('Lone pairs', (d) => { const t = d.atoms.get(id)!; t.lonePairs = !t.lonePairs; if (!t.lonePairs) delete t.lonePairs; }) },
        { label: 'Isotope…', run: () => { const v = prompt('Mass number (empty = natural):', a.isotope ? String(a.isotope) : ''); if (v === null) return; ed.mutate('Isotope', (d) => { const t = d.atoms.get(id)!; if (v.trim()) t.isotope = +v; else delete t.isotope; }); } },
        { label: 'Attach phenyl', keys: 'R', run: () => ed.mutate('Add phenyl', (d) => attachRingToAtom(d, id, 6, true)) },
      );
      if (a.abbrev) items.push({ label: `Expand “${a.abbrev}”`, run: () => ed.mutate('Expand label', (d) => expandLabel(d, id)) });
      else items.push({ label: 'Draw hydrogens as atoms', run: () => this.drawHydrogenAtoms([id]) });
      items.push({ label: 'Select fragment', run: () => ed.selectFragments([id]) }, { separator: true, label: '' }, { label: 'Delete atom', keys: 'Del', run: () => ed.mutate('Delete', (d) => deleteSelection(d, { atoms: new Set([id]), bonds: new Set(), objects: new Set() })) });
    } else if (hit?.kind === 'bond') {
      const id = hit.id;
      const set = (o: number, s: Parameters<typeof applyBondType>[2]) => () => ed.mutate('Change bond', (d) => { const b = d.bonds.get(id)!; b.order = o; b.style = s; });
      items.push(
        { label: 'Single', keys: '1', run: set(1, 'plain') },
        { label: 'Double', keys: '2', run: set(2, 'plain') },
        { label: 'Triple', keys: '3', run: set(3, 'plain') },
        { label: 'Stereo', submenu: [
          { label: 'Wedge', run: set(1, 'wedge') }, { label: 'Hash', run: set(1, 'hash') }, { label: 'Wavy', run: set(1, 'wavy') },
          { label: 'Bold', run: set(1, 'bold') }, { label: 'Dashed', run: set(1, 'dashed') }, { label: 'Hollow wedge', run: set(1, 'hollow') },
          { label: 'Crossed double (E/Z unknown)', run: set(2, 'crossed') },
        ] },
        { label: 'Flip direction', keys: 'F', run: () => ed.mutate('Flip bond', (d) => { const b = d.bonds.get(id)!; const t = b.a; b.a = b.b; b.b = t; }) },
        { label: 'Double-bond position', submenu: (['auto', 'left', 'right', 'center'] as const).map((p) => ({ label: p, run: () => ed.mutate('Bond position', (d) => (d.bonds.get(id)!.dbPos = p)) })) },
        { label: 'Fuse benzene', keys: 'R', run: () => ed.mutate('Fuse ring', (d) => fuseRingOnBond(d, id, 6, true, ed.mouse)) },
        { label: 'Select fragment', run: () => ed.selectFragments([ed.doc.bonds.get(id)!.a]) },
        { separator: true, label: '' },
        { label: 'Delete bond', keys: 'Del', run: () => ed.mutate('Delete', (d) => deleteSelection(d, { atoms: new Set(), bonds: new Set([id]), objects: new Set() })) },
      );
    } else if (hit && (hit.kind === 'arrow' || hit.kind === 'text' || hit.kind === 'shape' || hit.kind === 'curved')) {
      const id = hit.id;
      if (hit.kind === 'arrow') {
        const a = ed.doc.arrows.get(id)!;
        items.push(
          { label: 'Edit reagents (above)…', run: () => editArrowCaption(ed, id, { x: (a.x1 + a.x2) / 2, y: Math.min(a.y1, a.y2) - 1 }) },
          { label: 'Edit conditions (below)…', run: () => editArrowCaption(ed, id, { x: (a.x1 + a.x2) / 2, y: Math.max(a.y1, a.y2) + 1 }) },
          { label: 'Arrow type', submenu: (['reaction', 'equilibrium', 'unbalancedEq', 'retro', 'resonance', 'dashed', 'noGo', 'line'] as const).map((k) => ({ label: k, run: () => ed.mutate('Arrow type', (d) => (d.arrows.get(id)!.kind = k)) })) },
        );
      }
      if (hit.kind === 'text') items.push({ label: 'Edit text…', run: () => editTextObject(ed, id) }, { label: 'Toggle formula formatting', run: () => ed.mutate('Formula style', (d) => { const t = d.texts.get(id)!; t.formula = !t.formula; }) });
      if (hit.kind === 'curved') items.push(
        { label: 'Flip curvature', run: () => ed.mutate('Flip arrow', (d) => { const c = d.curved.get(id)!; c.c1.h = -c.c1.h; c.c2.h = -c.c2.h; }) },
        { label: 'Toggle pair / single electron', run: () => ed.mutate('Arrow electrons', (d) => { const c = d.curved.get(id)!; c.electrons = c.electrons === 2 ? 1 : 2; }) },
        { label: 'Apply this arrow group', run: () => this.mechanismPanel.applyGroupOf(id) },
      );
      items.push({ separator: true, label: '' }, { label: 'Delete', keys: 'Del', run: () => ed.mutate('Delete', (d) => deleteSelection(d, { atoms: new Set(), bonds: new Set(), objects: new Set([id]) })) });
    } else {
      const hasSel = ed.hasSelection();
      items.push(
        { label: 'Paste', keys: 'Ctrl V', run: () => this.paste() },
        { label: 'Select all', keys: 'Ctrl A', run: () => ed.selectAll() },
        { label: hasSel ? 'Clean selection' : 'Clean all', keys: 'Ctrl ⇧ K', run: () => this.clean() },
        { label: 'Insert from PubChem / SMILES…', keys: 'Ctrl K', run: () => this.openPalette() },
        { label: 'Zoom to fit', keys: 'Ctrl 0', run: () => ed.fitToContent() },
      );
      if (hasSel) items.push({ separator: true, label: '' }, { label: 'Copy', run: () => this.copy() }, { label: 'Copy as image', run: () => copyImage(this, true) }, { label: 'Generate 3D', run: () => { this.showTab('3d', true); this.viewerPanel.generate(); } }, { label: 'Delete selection', run: () => ed.deleteSelected() });
    }
    showMenu(x, y, items);
  }

  // ───────────── global input ─────────────

  private bindGlobal(): void {
    const ed = this.editor;
    window.addEventListener('keydown', (e) => {
      if (!ed.isCanvasContext(e)) return;
      const mod = e.ctrlKey || e.metaKey;
      const k = e.key.toLowerCase();
      if (mod && k === 'k' && !e.shiftKey) { e.preventDefault(); this.openPalette(); return; }
      if (mod && k === 'k' && e.shiftKey) { e.preventDefault(); this.clean(); return; }
      if (mod && k === 's') { e.preventDefault(); this.save(e.shiftKey); return; }
      if (mod && k === 'o') { e.preventDefault(); this.openFile(); return; }
      if (mod && k === 'e') { e.preventDefault(); exportDialog(this); return; }
      if (mod && k === 'd') { e.preventDefault(); this.duplicate(); return; }
      if (mod && k === 'f') { e.preventDefault(); this.showTab('pubchem', true); this.pubchemPanel.focus(); return; }
      if (!mod && e.key === '?') { e.preventDefault(); shortcutsDialog(); return; }
      if (mod && (k === 'c' || k === 'x' || k === 'v')) return; // handled by clipboard events
      if (ed.handleKey(e)) e.preventDefault();
    });
    document.addEventListener('copy', (e) => {
      if (!ed.isCanvasContext(e as unknown as KeyboardEvent) || window.getSelection()?.toString()) return;
      e.preventDefault();
      this.copyToEvent(e as ClipboardEvent);
    });
    document.addEventListener('cut', (e) => {
      if (!ed.isCanvasContext(e as unknown as KeyboardEvent) || window.getSelection()?.toString()) return;
      e.preventDefault();
      this.copyToEvent(e as ClipboardEvent);
      ed.deleteSelected();
    });
    document.addEventListener('paste', (e) => {
      if (!ed.isCanvasContext(e as unknown as KeyboardEvent)) return;
      e.preventDefault();
      const ce = e as ClipboardEvent;
      const t = ce.clipboardData?.getData('text/plain') ?? '';
      const files = ce.clipboardData?.files;
      if (files && files.length && !t) {
        this.openFile(files[0]);
        return;
      }
      this.paste(t);
    });
    // drag & drop files
    const stage = ed.container;
    stage.addEventListener('dragover', (e) => {
      e.preventDefault();
      stage.classList.add('drop');
    });
    stage.addEventListener('dragleave', () => stage.classList.remove('drop'));
    stage.addEventListener('drop', (e) => {
      e.preventDefault();
      stage.classList.remove('drop');
      const f = e.dataTransfer?.files?.[0];
      if (f) this.openFile(f);
      else {
        const t = e.dataTransfer?.getData('text/plain');
        if (t) this.insertFromText(t, null);
      }
    });
    window.addEventListener('beforeunload', () => this.autosave());
    window.addEventListener('hashchange', () => this.loadHash());
    matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => this.applyTheme());
    window.addEventListener('resize', () => this.updateSelBar());
  }

  private copyToEvent(e: ClipboardEvent): void {
    const ed = this.editor;
    const json = ed.selectionToJSON();
    const smi = this.currentSmiles();
    ed.setInternalClipboard(json, smi);
    e.clipboardData?.setData('text/plain', smi);
    e.clipboardData?.setData('application/x-chirally+json', json);
    toast('Copied — pastes as a structure here or as SMILES elsewhere', 'success', 1600);
  }

  // ───────────── theme & prefs ─────────────

  isDark(): boolean {
    if (this.themeMode === 'auto') return matchMedia('(prefers-color-scheme: dark)').matches;
    return this.themeMode === 'dark';
  }

  toggleTheme(): void {
    this.themeMode = this.isDark() ? 'light' : 'dark';
    this.prefs.theme = this.themeMode;
    this.persistPrefs();
    this.applyTheme();
  }

  private applyTheme(): void {
    const dark = this.isDark();
    document.documentElement.dataset.theme = dark ? 'dark' : 'light';
    this.editor.theme = dark ? DARK_THEME : LIGHT_THEME;
    this.editor.requestRender();
    const tb = this.root.querySelector('.theme-btn');
    if (tb) tb.innerHTML = dark ? ICONS.sun : ICONS.moon;
  }

  private loadPrefs(): void {
    try {
      this.prefs = JSON.parse(localStorage.getItem(PREFS_KEY) ?? '{}');
    } catch {
      this.prefs = {};
    }
    this.themeMode = this.prefs.theme ?? 'auto';
  }

  persistPrefs(): void {
    if (!this.editor) return;
    this.prefs.grid = { show: this.editor.grid.show, snap: this.editor.grid.snap };
    this.prefs.wheelZooms = this.editor.wheelZooms;
    this.prefs.style = { ...this.editor.doc.style };
    try {
      localStorage.setItem(PREFS_KEY, JSON.stringify(this.prefs));
    } catch {
      /* storage blocked */
    }
  }

  // ───────────── persistence ─────────────

  private autosave(): void {
    try {
      const doc = this.editor.doc;
      if (isDocEmpty(doc)) {
        localStorage.removeItem(AUTOSAVE_KEY);
        return;
      }
      localStorage.setItem(AUTOSAVE_KEY, docToJSON(doc));
      this.recordVersion(false);
    } catch {
      /* storage full or blocked */
    }
  }

  /** Keeps a rolling history of snapshots (every ≥ 2 minutes of work, or on demand). */
  private recordVersion(force: boolean): void {
    const now = Date.now();
    if (!force && now - this.lastVersionTime < 120000) return;
    const doc = this.editor.doc;
    if (isDocEmpty(doc)) return;
    this.lastVersionTime = now;
    try {
      const list = this.listVersions();
      list.unshift({ time: now, title: doc.meta.title, atoms: doc.atoms.size, json: docToJSON(doc), thumb: thumbOf(doc) });
      localStorage.setItem(VERSIONS_KEY, JSON.stringify(list.slice(0, 25)));
    } catch {
      /* ignore */
    }
  }

  listVersions(): { time: number; title: string; atoms: number; json: string; thumb: string }[] {
    try {
      return JSON.parse(localStorage.getItem(VERSIONS_KEY) ?? '[]');
    } catch {
      return [];
    }
  }

  restoreVersion(time: number): void {
    const v = this.listVersions().find((x) => x.time === time);
    if (!v) return;
    this.recordVersion(true);
    this.editor.begin();
    const d = docFromJSON(v.json);
    this.editor.doc = d;
    this.editor.commit('Restore version');
    this.updateTitle();
    setTimeout(() => this.editor.fitToContent(), 0);
    toast('Version restored (undo with Ctrl+Z)', 'success');
  }

  private async restore(): Promise<void> {
    const ed = this.editor;
    if (this.prefs.grid) Object.assign(ed.grid, this.prefs.grid);
    if (this.prefs.wheelZooms !== undefined) ed.wheelZooms = this.prefs.wheelZooms;
    if (this.prefs.style) ed.doc.style = { ...STYLE_PRESETS.Chirally, ...this.prefs.style };
    if (await this.loadHash()) return;
    try {
      const saved = localStorage.getItem(AUTOSAVE_KEY);
      if (saved) {
        const d = docFromJSON(saved);
        ed.setDoc(d);
        this.updateTitle();
        setTimeout(() => ed.fitToContent(), 30);
        ed.setStatus('Restored your last drawing (autosaved in this browser)');
      }
    } catch {
      /* corrupted autosave */
    }
    ed.setTool('bond');
    if (!localStorage.getItem(SEEN_KEY)) {
      try {
        localStorage.setItem(SEEN_KEY, '1');
      } catch {
        /* ignore */
      }
      if (!navigator.webdriver) setTimeout(() => welcomeDialog(this), 300);
    }
    this.stylePanel.refresh();
  }

  private async loadHash(): Promise<boolean> {
    if (!location.hash) return false;
    try {
      const d = await docFromHash(location.hash);
      if (!d) return false;
      if (!isDocEmpty(this.editor.doc)) this.recordVersion(true);
      this.editor.setDoc(d);
      this.updateTitle();
      setTimeout(() => this.editor.fitToContent(), 30);
      history.replaceState(null, '', location.pathname + location.search);
      toast('Opened shared drawing', 'success');
      return true;
    } catch (e) {
      toast('Could not open the shared link: ' + (e as Error).message, 'error');
      return false;
    }
  }

  private updateTitle(): void {
    const t = this.editor.doc.meta.title || 'Untitled';
    if (this.titleEl && document.activeElement !== this.titleEl) this.titleEl.value = t;
    document.title = `${t} — Chirally`;
  }
}

function allSel(doc: ChemDoc) {
  const s = emptySelection();
  for (const id of doc.atoms.keys()) s.atoms.add(id);
  for (const id of doc.bonds.keys()) s.bonds.add(id);
  return s;
}

function isDocEmptyExcept(doc: ChemDoc, n: number): boolean {
  return doc.atoms.size === n && !doc.arrows.size && !doc.texts.size;
}

function thumbOf(doc: ChemDoc): string {
  try {
    const { mol } = docToMol(doc);
    return molThumbnail(mol);
  } catch {
    return '';
  }
}

export { safeName, cloneDoc, layoutMol };
