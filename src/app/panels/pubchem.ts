// PubChem search, structure import, identification of drawn structures, CAS numbers and GHS hazards.
import type { App } from '../app';
import { h, clear, fmt, copyText, svgEl, debounce, toast } from '../dom';
import { ICONS } from '../icons';
import * as PC from '../../services/pubchem';

export class PubChemPanel {
  el: HTMLElement;
  private input: HTMLInputElement;
  private suggest: HTMLElement;
  private results: HTMLElement;
  private detail: HTMLElement;
  private seq = 0;

  constructor(private app: App) {
    this.input = h('input', {
      type: 'search', class: 'input', placeholder: 'Name, CAS number, InChIKey or CID…', 'aria-label': 'Search PubChem',
      autocomplete: 'off', spellcheck: false,
    }) as HTMLInputElement;
    this.suggest = h('div', { class: 'suggest', role: 'listbox' });
    this.results = h('div', { class: 'results' });
    this.detail = h('div', { class: 'detail' });
    const doSuggest = debounce(() => this.loadSuggestions(), 220);
    this.input.addEventListener('input', () => doSuggest());
    this.input.addEventListener('keydown', (e) => {
      if (e.key === 'Enter') {
        e.preventDefault();
        this.suggest.innerHTML = '';
        this.search(this.input.value);
      }
    });
    this.el = h(
      'section',
      { class: 'panel', 'aria-label': 'PubChem' },
      h('div', { class: 'panel-body' },
        h('div', { class: 'search-row' }, this.input, h('button', { class: 'btn btn-primary', onclick: () => this.search(this.input.value) }, 'Search')),
        this.suggest,
        h('div', { class: 'row-actions' },
          h('button', { class: 'btn btn-small', onclick: () => this.identifyCurrent() }, svgEl(ICONS.sparkle), 'Identify drawn structure'),
          h('button', { class: 'btn btn-small', onclick: () => this.similarCurrent() }, svgEl(ICONS.search), 'Find similar'),
        ),
        this.results,
        this.detail,
        h('p', { class: 'muted small' }, 'Data from PubChem (NCBI). Hazard data are aggregated GHS notifications — always consult the supplier SDS.'),
      ),
    );
  }

  focus(): void {
    this.input.focus();
    this.input.select();
  }

  private async loadSuggestions(): Promise<void> {
    const q = this.input.value.trim();
    const my = ++this.seq;
    if (q.length < 2) {
      this.suggest.innerHTML = '';
      return;
    }
    try {
      const terms = await PC.autocomplete(q);
      if (my !== this.seq) return;
      clear(this.suggest);
      for (const t of terms) {
        const b = h('button', { class: 'suggest-item', role: 'option' }, t);
        b.addEventListener('click', () => {
          this.input.value = t;
          this.suggest.innerHTML = '';
          this.search(t);
        });
        this.suggest.appendChild(b);
      }
    } catch {
      /* suggestions are best-effort */
    }
  }

  async search(q: string): Promise<void> {
    q = q.trim();
    if (!q) return;
    this.input.value = q;
    clear(this.detail);
    this.results.innerHTML = '<div class="loading">Searching PubChem…</div>';
    try {
      const list = await PC.search(q, 12);
      this.showResults(list, `No PubChem results for “${q}”.`);
      if (list.length === 1) this.showDetail(list[0]);
    } catch (e) {
      this.results.innerHTML = '';
      this.results.appendChild(h('p', { class: 'warn' }, (e as Error).message));
    }
  }

  private showResults(list: PC.CompoundSummary[], emptyMsg: string): void {
    clear(this.results);
    if (!list.length) {
      this.results.appendChild(h('p', { class: 'muted' }, emptyMsg));
      return;
    }
    for (const c of list) {
      const item = h(
        'div',
        { class: 'result', tabindex: 0 },
        h('img', { src: PC.imageUrl(c.cid, 120), alt: c.title ?? `CID ${c.cid}`, loading: 'lazy', width: 64, height: 64 }),
        h('div', { class: 'result-text' },
          h('div', { class: 'result-title' }, c.title ?? `CID ${c.cid}`),
          h('div', { class: 'muted small', html: `${formulaHtml(c.formula ?? '')} · ${c.mw ? fmt(c.mw, 2) : '?'} g/mol · CID ${c.cid}` }),
        ),
        h('button', { class: 'btn btn-small', title: 'Insert into drawing', onclick: (e: Event) => { e.stopPropagation(); this.insert(c); } }, 'Insert'),
      );
      item.addEventListener('click', () => this.showDetail(c));
      item.addEventListener('keydown', (e) => {
        if ((e as KeyboardEvent).key === 'Enter') this.showDetail(c);
      });
      this.results.appendChild(item);
    }
  }

  async insert(c: PC.CompoundSummary): Promise<void> {
    try {
      const sdf = await PC.sdf2d(c.cid);
      this.app.insertFromText(sdf, `${c.title ?? 'CID' + c.cid}.sdf`, c.title ?? undefined, { suppressH: true });
      toast(`Inserted ${c.title ?? 'CID ' + c.cid}`, 'success');
    } catch (e) {
      // fall back to SMILES
      if (c.smiles) {
        this.app.insertFromText(c.smiles, null, c.title);
      } else toast((e as Error).message, 'error');
    }
  }

  async showDetail(c: PC.CompoundSummary): Promise<void> {
    clear(this.detail);
    const card = h('div', { class: 'card' });
    this.detail.appendChild(card);
    card.appendChild(
      h('div', { class: 'detail-head' },
        h('img', { src: PC.imageUrl(c.cid, 240), alt: '', width: 120, height: 120 }),
        h('div', null,
          h('h4', null, c.title ?? `CID ${c.cid}`),
          h('div', { class: 'row-actions' },
            h('button', { class: 'btn btn-primary btn-small', onclick: () => this.insert(c) }, 'Insert'),
            h('a', { class: 'btn btn-small', href: PC.pageUrl(c.cid), target: '_blank', rel: 'noopener' }, 'Open in PubChem ↗'),
          ),
        ),
      ),
    );
    const kv = (k: string, v: string | undefined, copy = true) =>
      v ? h('div', { class: 'kv' }, h('span', { class: 'k' }, k), h('span', { class: 'v' }, v), copy ? h('button', { class: 'icon-btn tiny', 'aria-label': 'Copy ' + k, onclick: () => copyText(v, k) }, svgEl(ICONS.copy)) : null) : null;
    card.appendChild(h('div', null,
      kv('IUPAC name', c.iupacName),
      kv('Formula', c.formula),
      kv('Mol. weight', c.mw ? `${fmt(c.mw, 2)} g/mol` : undefined),
      kv('XLogP3', c.xlogp !== undefined ? String(c.xlogp) : undefined, false),
      kv('TPSA', c.tpsa !== undefined ? `${c.tpsa} Å²` : undefined, false),
      kv('SMILES', c.smiles),
      kv('InChIKey', c.inchiKey),
      kv('CID', String(c.cid)),
    ));
    const casBox = h('div', { class: 'kv' }, h('span', { class: 'k' }, 'CAS RN'), h('span', { class: 'v muted' }, 'loading…'));
    card.appendChild(casBox);
    const synBox = h('details', { class: 'sub' }, h('summary', null, 'Synonyms'));
    card.appendChild(synBox);
    const hazBox = h('div', { class: 'hazards' }, h('div', { class: 'muted small' }, 'Loading GHS hazards…'));
    card.appendChild(hazBox);
    PC.synonyms(c.cid).then((syns) => {
      const cas = PC.casNumbers(syns);
      casBox.querySelector('.v')!.textContent = cas.length ? cas.slice(0, 3).join(', ') : 'not listed';
      casBox.querySelector('.v')!.classList.toggle('muted', !cas.length);
      if (cas.length) casBox.appendChild(h('button', { class: 'icon-btn tiny', 'aria-label': 'Copy CAS', onclick: () => copyText(cas[0], 'CAS number') }, svgEl(ICONS.copy)));
      synBox.appendChild(h('div', { class: 'syn-list' }, syns.slice(0, 40).map((s) => h('span', { class: 'chip' }, s))));
    });
    PC.hazards(c.cid)
      .then((hz) => {
        clear(hazBox);
        if (!hz) {
          hazBox.appendChild(h('p', { class: 'muted small' }, 'No GHS classification reported in PubChem.'));
          return;
        }
        hazBox.appendChild(h('h5', null, 'GHS hazards ', hz.signal ? h('span', { class: 'signal signal-' + hz.signal.toLowerCase() }, hz.signal) : null));
        if (hz.pictograms.length)
          hazBox.appendChild(h('div', { class: 'pictos' }, hz.pictograms.map((p) => h('img', { src: p.url, alt: p.name, title: `${p.code}: ${p.name}`, width: 44, height: 44 }))));
        if (hz.statements.length) hazBox.appendChild(h('ul', { class: 'hstatements' }, hz.statements.slice(0, 12).map((s) => h('li', null, s))));
        if (hz.precautionary) hazBox.appendChild(h('details', { class: 'sub' }, h('summary', null, 'Precautionary codes'), h('p', { class: 'small' }, hz.precautionary)));
      })
      .catch(() => {
        clear(hazBox);
        hazBox.appendChild(h('p', { class: 'muted small' }, 'Hazard data unavailable.'));
      });
  }

  async identify(smiles: string): Promise<void> {
    if (!smiles) {
      toast('Draw a structure first', 'error');
      return;
    }
    clear(this.detail);
    this.results.innerHTML = '<div class="loading">Looking up the drawn structure…</div>';
    try {
      const c = await PC.identifyBySmiles(smiles);
      clear(this.results);
      if (!c) {
        this.results.appendChild(h('p', { class: 'muted' }, 'This exact structure is not in PubChem — it may be novel. Try “Find similar”.'));
        return;
      }
      this.results.appendChild(h('p', { class: 'ok' }, `✓ Found: ${c.title ?? 'CID ' + c.cid}`));
      this.showDetail(c);
    } catch (e) {
      this.results.innerHTML = '';
      this.results.appendChild(h('p', { class: 'warn' }, (e as Error).message));
    }
  }

  identifyCurrent(): void {
    this.identify(this.app.currentSmiles());
  }

  async similarCurrent(): Promise<void> {
    const smi = this.app.currentSmiles();
    if (!smi) return toast('Draw a structure first', 'error');
    clear(this.detail);
    this.results.innerHTML = '<div class="loading">Searching similar compounds (Tanimoto ≥ 90%)…</div>';
    try {
      const list = await PC.similar(smi, 90, 12);
      this.showResults(list, 'No similar compounds found.');
    } catch (e) {
      this.results.innerHTML = '';
      this.results.appendChild(h('p', { class: 'warn' }, (e as Error).message));
    }
  }
}

function formulaHtml(f: string): string {
  return f.replace(/(\d+)/g, '<sub>$1</sub>').replace(/([+-])$/, '<sup>$1</sup>');
}
