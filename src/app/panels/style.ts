// Document style & display options.
import type { App } from '../app';
import { h, clear } from '../dom';
import { STYLE_PRESETS, DocStyle } from '../../doc/types';

export class StylePanel {
  el: HTMLElement;
  private body: HTMLElement;

  constructor(private app: App) {
    this.body = h('div', { class: 'panel-body' });
    this.el = h('section', { class: 'panel', 'aria-label': 'Style' }, this.body);
    this.refresh();
  }

  refresh(): void {
    const ed = this.app.editor;
    const st = ed.doc.style;
    clear(this.body);
    const setStyle = (patch: Partial<DocStyle>, label = 'Change style') => {
      ed.begin();
      Object.assign(ed.doc.style, patch);
      ed.commit(label);
      this.app.persistPrefs();
    };
    const preset = h('select', { class: 'input', 'aria-label': 'Style preset' }, Object.keys(STYLE_PRESETS).map((k) => h('option', { value: k, selected: st.name === k }, k))) as HTMLSelectElement;
    preset.addEventListener('change', () => {
      setStyle({ ...STYLE_PRESETS[preset.value] }, 'Apply style preset');
      this.refresh();
    });
    const toggle = (label: string, key: keyof DocStyle, hint?: string) => {
      const cb = h('input', { type: 'checkbox', checked: !!st[key] }) as HTMLInputElement;
      cb.addEventListener('change', () => setStyle({ [key]: cb.checked } as Partial<DocStyle>));
      return h('label', { class: 'check', title: hint }, cb, label);
    };
    const slider = (label: string, key: keyof DocStyle, min: number, max: number, step: number, fmtv: (v: number) => string) => {
      const out = h('span', { class: 'muted small' }, fmtv(st[key] as number));
      const r = h('input', { type: 'range', min, max, step, value: st[key] as number, 'aria-label': label }) as HTMLInputElement;
      r.addEventListener('input', () => {
        (ed.doc.style as any)[key] = +r.value;
        out.textContent = fmtv(+r.value);
        ed.touch();
      });
      r.addEventListener('change', () => this.app.persistPrefs());
      return h('label', { class: 'slider' }, h('span', null, label), r, out);
    };
    const carbons = h('select', { class: 'input input-small', 'aria-label': 'Carbon labels' },
      h('option', { value: 'none', selected: st.showCarbons === 'none' }, 'Skeletal (hidden)'),
      h('option', { value: 'terminal', selected: st.showCarbons === 'terminal' }, 'Terminal CH₃'),
      h('option', { value: 'all', selected: st.showCarbons === 'all' }, 'All carbons'),
    ) as HTMLSelectElement;
    carbons.addEventListener('change', () => setStyle({ showCarbons: carbons.value as DocStyle['showCarbons'] }));

    this.body.append(
      h('h4', null, 'Drawing style'),
      h('label', { class: 'field' }, h('span', null, 'Preset'), preset),
      slider('Line width', 'lineWidth', 0.02, 0.12, 0.005, (v) => (v * st.bondLengthPt).toFixed(2) + ' pt'),
      slider('Bold width', 'boldWidth', 0.08, 0.3, 0.01, (v) => (v * st.bondLengthPt).toFixed(1) + ' pt'),
      slider('Double-bond spacing', 'bondSpacing', 0.1, 0.3, 0.01, (v) => Math.round(v * 100) + '%'),
      slider('Label size', 'fontSize', 0.35, 0.9, 0.01, (v) => (v * st.bondLengthPt).toFixed(1) + ' pt'),
      slider('Bond length (export)', 'bondLengthPt', 10, 40, 0.2, (v) => v.toFixed(1) + ' pt'),
      h('h4', null, 'Display'),
      h('label', { class: 'field' }, h('span', null, 'Carbon labels'), carbons),
      toggle('Colour atoms by element', 'colorAtoms'),
      toggle('Show implicit hydrogens on labels', 'showImplicitH'),
      toggle('Show R/S and E/Z labels', 'showStereoLabels', 'CIP descriptors computed live'),
      toggle('Show all lone pairs', 'showLonePairs'),
      toggle('Aromatic rings as circles', 'aromaticCircles'),
      toggle('Flag valence errors', 'showValenceErrors'),
      h('h4', null, 'Canvas'),
      this.gridControls(),
      this.wheelControl(),
    );
  }

  private gridControls(): HTMLElement {
    const ed = this.app.editor;
    const show = h('input', { type: 'checkbox', checked: ed.grid.show }) as HTMLInputElement;
    const snap = h('input', { type: 'checkbox', checked: ed.grid.snap }) as HTMLInputElement;
    show.addEventListener('change', () => {
      ed.grid.show = show.checked;
      ed.requestRender();
      this.app.persistPrefs();
    });
    snap.addEventListener('change', () => {
      ed.grid.snap = snap.checked;
      this.app.persistPrefs();
    });
    return h('div', null, h('label', { class: 'check' }, show, 'Show grid'), h('label', { class: 'check' }, snap, 'Snap new atoms & objects to grid'));
  }

  private wheelControl(): HTMLElement {
    const ed = this.app.editor;
    const sel = h('select', { class: 'input input-small', 'aria-label': 'Mouse wheel' },
      h('option', { value: 'zoom', selected: ed.wheelZooms }, 'Mouse wheel zooms'),
      h('option', { value: 'pan', selected: !ed.wheelZooms }, 'Mouse wheel pans (Ctrl+wheel zooms)'),
    ) as HTMLSelectElement;
    sel.addEventListener('change', () => {
      ed.wheelZooms = sel.value === 'zoom';
      this.app.persistPrefs();
    });
    return h('label', { class: 'field' }, h('span', null, 'Wheel'), sel);
  }
}
