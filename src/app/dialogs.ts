// Modal dialogs: periodic table, export, keyboard shortcuts, share link, version history, welcome.
import type { App } from './app';
import { h, modal, toast, copyText, svgEl } from './dom';
import { ICONS } from './icons';
import { PERIODIC_LAYOUT, element, labelColor } from '../chem/elements';
import { ExportFormat, exportDoc, docPNG, docSVG, shareLink, reactionSmiles } from './fileio';
import { ChemDoc } from '../doc/types';

export function periodicTableDialog(onPick: (el: string) => void): void {
  const grid = h('div', { class: 'ptable', role: 'grid' });
  const info = h('div', { class: 'pt-info muted small' }, 'Hover an element for details');
  let close: () => void = () => undefined;
  for (const row of PERIODIC_LAYOUT) {
    for (const sym of row) {
      if (!sym) {
        grid.appendChild(h('span', { class: 'pt-gap' }));
        continue;
      }
      const e = element(sym)!;
      const cls = e.group === 0 ? 'f' : e.group <= 2 ? 's' : e.group <= 12 ? 'd' : 'p';
      const b = h('button', { class: `pt-cell pt-${cls}`, style: `--el:${labelColor(sym)}`, title: `${e.name} (${e.z})` },
        h('span', { class: 'pt-z' }, String(e.z)), h('span', { class: 'pt-sym' }, sym));
      b.addEventListener('mouseenter', () => {
        info.textContent = `${e.z} ${e.name} — ${e.mass} g/mol · EN ${e.en || '—'} · covalent r ${e.covRadius} Å`;
      });
      b.addEventListener('click', () => {
        onPick(sym);
        close();
      });
      grid.appendChild(b);
    }
  }
  close = modal('Periodic table', h('div', null, grid, info), { wide: true });
}

export function exportDialog(app: App): void {
  const ed = app.editor;
  const hasSel = ed.hasSelection();
  let fmt: ExportFormat = 'png';
  const formats: { id: ExportFormat; label: string; desc: string }[] = [
    { id: 'png', label: 'PNG image', desc: 'Raster image for slides, documents, ELNs' },
    { id: 'svg', label: 'SVG vector', desc: 'Publication-quality vector graphics' },
    { id: 'cdxml', label: 'ChemDraw CDXML', desc: 'Open in ChemDraw / ChemOffice' },
    { id: 'mol', label: 'MDL Molfile', desc: 'V2000/V3000 connection table' },
    { id: 'sdf', label: 'SD file', desc: 'One record per fragment' },
    { id: 'smiles', label: 'SMILES', desc: 'Canonical, with stereo' },
    { id: 'rxn', label: 'RXN file', desc: 'Reaction (left/right of arrow)' },
    { id: 'cml', label: 'CML', desc: 'Chemical Markup Language' },
    { id: 'chirally', label: 'Chirally (.chirally)', desc: 'Native, lossless' },
  ];
  const list = h('div', { class: 'fmt-list', role: 'radiogroup' });
  const scale = h('select', { class: 'input input-small', 'aria-label': 'PNG scale' }, [1, 2, 3, 4, 6].map((s) => h('option', { value: s, selected: s === 3 }, `${s}× (${Math.round(96 * s)} dpi)`))) as HTMLSelectElement;
  const transparent = h('input', { type: 'checkbox' }) as HTMLInputElement;
  const onlySel = h('input', { type: 'checkbox', checked: hasSel, disabled: !hasSel }) as HTMLInputElement;
  const preview = h('div', { class: 'export-preview' });
  const refreshPreview = () => {
    const d = onlySel.checked ? app.selectionDoc() : ed.doc;
    preview.innerHTML = '';
    const img = h('img', { alt: 'Preview', src: 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(docSVG(d, true)) });
    preview.appendChild(img);
  };
  formats.forEach((f) => {
    const r = h('input', { type: 'radio', name: 'fmt', value: f.id, checked: f.id === fmt }) as HTMLInputElement;
    r.addEventListener('change', () => (fmt = f.id));
    list.appendChild(h('label', { class: 'fmt' }, r, h('span', null, h('b', null, f.label), h('span', { class: 'muted small' }, f.desc))));
  });
  onlySel.addEventListener('change', refreshPreview);
  const body = h('div', { class: 'export-dlg' },
    h('div', null, list),
    h('div', null,
      preview,
      h('label', { class: 'check' }, onlySel, 'Selection only'),
      h('label', { class: 'check' }, transparent, 'Transparent background (PNG/SVG)'),
      h('label', { class: 'field' }, h('span', null, 'PNG resolution'), scale),
      h('div', { class: 'row-actions wrap' },
        h('button', { class: 'btn btn-small', onclick: () => copyImage(app, onlySel.checked) }, svgEl(ICONS.copy), 'Copy image'),
        h('button', { class: 'btn btn-small', onclick: () => copyText(app.currentSmiles(onlySel.checked), 'SMILES') }, 'Copy SMILES'),
        h('button', { class: 'btn btn-small', onclick: () => { const r = reactionSmiles(ed.doc); if (r) copyText(r, 'Reaction SMILES'); else toast('No reaction arrow found', 'error'); } }, 'Copy reaction SMILES'),
      ),
    ),
  );
  refreshPreview();
  modal('Export', body, {
    wide: true,
    actions: [
      { label: 'Cancel', onClick: () => undefined },
      {
        label: 'Export',
        primary: true,
        onClick: () => {
          exportDoc(ed.doc, fmt, { subDoc: onlySel.checked ? app.selectionDoc() : null, pngScale: +scale.value, transparent: transparent.checked }).catch((e) => toast(String(e), 'error'));
        },
      },
    ],
  });
}

export async function copyImage(app: App, selectionOnly = false): Promise<void> {
  const d: ChemDoc = selectionOnly ? app.selectionDoc() : app.editor.hasSelection() ? app.selectionDoc() : app.editor.doc;
  try {
    const png = await docPNG(d, 3, false);
    const svg = docSVG(d, true);
    const CI = (window as any).ClipboardItem;
    const item = new CI({ 'image/png': png, 'text/plain': new Blob([app.currentSmiles()], { type: 'text/plain' }) });
    await navigator.clipboard.write([item]);
    void svg;
    toast('Image copied — paste into Word, PowerPoint or an ELN', 'success');
  } catch {
    toast('Your browser blocked image copying — use Export → PNG instead', 'error');
  }
}

export const SHORTCUTS: [string, string][] = [
  ['Hover atom + C N O S P F H I', 'Change element (L = Cl, B = Br, Shift+S = Si)'],
  ['Hover atom + 1 / 2 / 3', 'Sprout single / double / triple bond'],
  ['Hover atom + + / − / 0', 'Charge up / down / reset'],
  ['Hover atom + . / ,', 'Radical / show lone pairs'],
  ['Hover atom + M E T A X', 'Me, Et, tBu, Ac, Ph labels'],
  ['Hover atom + R', 'Attach a phenyl ring'],
  ['Hover atom + Enter (or double-click)', 'Type a label (OMe, CO2H, NHBoc…)'],
  ['Hover bond + 1 2 3 4', 'Single / double / triple / aromatic'],
  ['Hover bond + W / H / Y / B / D', 'Wedge / hash / wavy / bold / dashed (W again flips)'],
  ['Hover bond + R or 5–8', 'Fuse benzene / n-membered ring'],
  ['Hover + Delete', 'Delete atom or bond'],
  ['V · E · G · R · T · A · Shift+A', 'Select · eraser · chain · ring · text · arrow · curved arrow'],
  ['1 · 2 · 3 · W · Q', 'Bond tools: single · double · triple · wedge · hash'],
  ['C N O S …(nothing hovered)', 'Atom tool with that element'],
  ['Space + drag / middle drag', 'Pan'],
  ['Wheel / pinch / Ctrl + = −', 'Zoom'],
  ['Ctrl+0 / Ctrl+1', 'Fit to window / 100%'],
  ['Ctrl+Z / Ctrl+Shift+Z', 'Undo / redo'],
  ['Ctrl+C / X / V / D', 'Copy / cut / paste / duplicate'],
  ['Ctrl+A', 'Select all'],
  ['Ctrl+Shift+K', 'Clean structure'],
  ['Ctrl+K', 'Command palette (PubChem search, SMILES, commands)'],
  ['Ctrl+S / Ctrl+O / Ctrl+E', 'Save / open / export'],
  ['Arrow keys (+Shift)', 'Nudge selection'],
  ['Alt while dragging', 'Disable angle/length snapping'],
  ['Shift while rotating', '15° steps'],
  ['Esc', 'Cancel · deselect · select tool'],
];

export function shortcutsDialog(): void {
  const t = h('table', { class: 'shortcuts' }, SHORTCUTS.map(([k, d]) => h('tr', null, h('td', null, h('kbd', null, k)), h('td', null, d))));
  modal('Keyboard shortcuts', t, { wide: true });
}

export async function shareDialog(app: App): Promise<void> {
  try {
    const url = await shareLink(app.editor.doc);
    const input = h('input', { class: 'input', value: url, readonly: true, 'aria-label': 'Share link' }) as HTMLInputElement;
    const note = url.length > 8000 ? h('p', { class: 'warn small' }, 'This drawing is large — very long links may not work in every app. Consider saving a .chirally file instead.') : null;
    modal('Share link', h('div', null, h('p', { class: 'small' }, 'Anyone with this link opens a copy of the drawing. Nothing is uploaded — the drawing is encoded in the link itself.'), input, note), {
      actions: [{ label: 'Copy link', primary: true, onClick: () => { copyText(url, 'Link'); } }],
    });
    input.select();
  } catch (e) {
    toast('Could not create link: ' + String(e), 'error');
  }
}

export function versionsDialog(app: App): void {
  const versions = app.listVersions();
  if (!versions.length) return toast('No saved versions yet — versions are recorded automatically as you work', 'info');
  const list = h('div', { class: 'versions' });
  let close: () => void = () => undefined;
  for (const v of versions) {
    list.appendChild(
      h('div', { class: 'version' },
        h('img', { src: v.thumb, alt: '' }),
        h('div', null, h('b', null, new Date(v.time).toLocaleString()), h('div', { class: 'muted small' }, `${v.atoms} atoms · ${v.title}`)),
        h('button', { class: 'btn btn-small', onclick: () => { app.restoreVersion(v.time); close(); } }, 'Restore'),
      ),
    );
  }
  close = modal('Version history', list, { wide: true });
}

export function welcomeDialog(app: App): void {
  const body = h('div', { class: 'welcome' },
    h('p', null, 'Chirally is a structure editor in your browser: draw molecules and reaction mechanisms, look compounds up in PubChem, and get IUPAC names, properties and 3D models as you draw.'),
    h('ul', null,
      h('li', null, h('b', null, 'Draw: '), 'click to add bonds, drag to draw at 15° steps, hover an atom and type a letter (O, N, Cl…) to change it.'),
      h('li', null, h('b', null, 'Mechanisms: '), 'draw curved arrows (Shift+A), then ', h('i', null, 'Apply arrows'), ' in the Mechanism tab. Chirally works out the next intermediate for you.'),
      h('li', null, h('b', null, 'PubChem: '), 'press Ctrl+K and type a name or CAS number (e.g. “niacinamide”, “50-78-2”).'),
      h('li', null, h('b', null, 'Formulation: '), 'the Analysis tab shows cLogP, TPSA, skin-permeability estimates and HLB for surfactants.'),
      h('li', null, h('b', null, 'Touch: '), 'tap to draw, pinch to zoom, two fingers to pan, long-press for the context menu.'),
    ),
  );
  modal('Welcome to Chirally', body, {
    actions: [
      { label: 'Keyboard shortcuts', onClick: () => { shortcutsDialog(); } },
      { label: 'Load an example', onClick: () => { app.loadExample(); } },
      { label: 'Start drawing', primary: true, onClick: () => undefined },
    ],
  });
}
