// Tiny DOM helpers (no framework).

type Attrs = Record<string, unknown> & { class?: string; style?: string | Partial<CSSStyleDeclaration> };
type Child = Node | string | number | null | undefined | false | Child[];

export function h<K extends keyof HTMLElementTagNameMap>(tag: K, attrs: Attrs | null = null, ...children: Child[]): HTMLElementTagNameMap[K] {
  const el = document.createElement(tag);
  if (attrs) {
    for (const [k, v] of Object.entries(attrs)) {
      if (v === undefined || v === null || v === false) continue;
      if (k === 'class') el.className = String(v);
      else if (k === 'style') {
        if (typeof v === 'string') el.setAttribute('style', v);
        else Object.assign(el.style, v);
      } else if (k.startsWith('on') && typeof v === 'function') {
        el.addEventListener(k.slice(2).toLowerCase(), v as EventListener);
      } else if (k === 'html') {
        el.innerHTML = String(v);
      } else if (k in el && typeof v !== 'string') {
        (el as any)[k] = v;
      } else {
        el.setAttribute(k, v === true ? '' : String(v));
      }
    }
  }
  append(el, children);
  return el;
}

function append(el: Node, children: Child[]): void {
  for (const c of children) {
    if (c === null || c === undefined || c === false) continue;
    if (Array.isArray(c)) append(el, c);
    else el.appendChild(typeof c === 'object' ? c : document.createTextNode(String(c)));
  }
}

/** Creates an element from an SVG/HTML string. */
export function svgEl(markup: string): HTMLElement {
  const t = document.createElement('template');
  t.innerHTML = markup.trim();
  return t.content.firstElementChild as HTMLElement;
}

export function clear(el: Element): void {
  while (el.firstChild) el.removeChild(el.firstChild);
}

let toastHost: HTMLElement | null = null;

/** Small transient notification. */
export function toast(msg: string, kind: 'info' | 'error' | 'success' = 'info', ms = 2600): void {
  if (!toastHost) {
    toastHost = h('div', { class: 'toasts', role: 'status', 'aria-live': 'polite' });
    document.body.appendChild(toastHost);
  }
  const t = h('div', { class: `toast toast-${kind}` }, msg);
  toastHost.appendChild(t);
  requestAnimationFrame(() => t.classList.add('show'));
  setTimeout(() => {
    t.classList.remove('show');
    setTimeout(() => t.remove(), 300);
  }, ms);
}

export function escapeHtml(s: string): string {
  return s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
}

export async function copyText(text: string, what = 'Copied'): Promise<void> {
  try {
    await navigator.clipboard.writeText(text);
    toast(`${what} to clipboard`, 'success');
  } catch {
    // fallback
    const ta = h('textarea', { style: 'position:fixed;left:-9999px' }) as HTMLTextAreaElement;
    ta.value = text;
    document.body.appendChild(ta);
    ta.select();
    try {
      document.execCommand('copy');
      toast(`${what} to clipboard`, 'success');
    } catch {
      toast('Copy failed', 'error');
    }
    ta.remove();
  }
}

export function debounce<T extends (...a: any[]) => void>(fn: T, ms: number): T & { flush(): void } {
  let t: number | null = null;
  let lastArgs: any[] = [];
  const d = ((...a: any[]) => {
    lastArgs = a;
    if (t !== null) clearTimeout(t);
    t = window.setTimeout(() => {
      t = null;
      fn(...lastArgs);
    }, ms);
  }) as T & { flush(): void };
  d.flush = () => {
    if (t !== null) {
      clearTimeout(t);
      t = null;
      fn(...lastArgs);
    }
  };
  return d;
}

export function fmt(v: number, digits = 2): string {
  if (!isFinite(v)) return '—';
  return v.toFixed(digits);
}

/** Generic modal dialog. Returns a close function. */
export function modal(title: string, body: Node, opts: { wide?: boolean; onClose?: () => void; actions?: { label: string; primary?: boolean; onClick: () => void | boolean }[] } = {}): () => void {
  const overlay = h('div', { class: 'modal-overlay' });
  const close = () => {
    overlay.remove();
    document.removeEventListener('keydown', onKey, true);
    opts.onClose?.();
  };
  const onKey = (e: KeyboardEvent) => {
    if (e.key === 'Escape') {
      e.stopPropagation();
      close();
    }
  };
  document.addEventListener('keydown', onKey, true);
  const dlg = h(
    'div',
    { class: 'modal' + (opts.wide ? ' modal-wide' : ''), role: 'dialog', 'aria-modal': 'true', 'aria-label': title },
    h('div', { class: 'modal-head' }, h('h2', null, title), h('button', { class: 'icon-btn', 'aria-label': 'Close', onclick: close, html: '&times;' })),
    h('div', { class: 'modal-body' }, body),
    opts.actions
      ? h(
          'div',
          { class: 'modal-actions' },
          opts.actions.map((a) =>
            h('button', { class: a.primary ? 'btn btn-primary' : 'btn', onclick: () => { if (a.onClick() !== false) close(); } }, a.label),
          ),
        )
      : null,
  );
  overlay.addEventListener('pointerdown', (e) => {
    if (e.target === overlay) close();
  });
  overlay.appendChild(dlg);
  document.body.appendChild(overlay);
  const focusable = dlg.querySelector('input, select, textarea, button.btn-primary') as HTMLElement | null;
  focusable?.focus();
  return close;
}

export function downloadBlob(blob: Blob, filename: string): void {
  const url = URL.createObjectURL(blob);
  const a = h('a', { href: url, download: filename });
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 5000);
}
