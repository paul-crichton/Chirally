// Right-click / long-press context menu.
import { h } from './dom';

export interface MenuItem {
  label: string;
  keys?: string;
  run?: () => void;
  submenu?: MenuItem[];
  disabled?: boolean;
  separator?: boolean;
}

let current: HTMLElement | null = null;

export function closeMenu(): void {
  current?.remove();
  current = null;
}

export function showMenu(x: number, y: number, items: MenuItem[]): void {
  closeMenu();
  const menu = buildMenu(items);
  document.body.appendChild(menu);
  const r = menu.getBoundingClientRect();
  menu.style.left = Math.max(4, Math.min(x, window.innerWidth - r.width - 4)) + 'px';
  menu.style.top = Math.max(4, Math.min(y, window.innerHeight - r.height - 4)) + 'px';
  current = menu;
  const off = (e: Event) => {
    if (current && !current.contains(e.target as Node)) {
      closeMenu();
      document.removeEventListener('pointerdown', off, true);
      document.removeEventListener('keydown', esc, true);
    }
  };
  const esc = (e: KeyboardEvent) => {
    if (e.key === 'Escape') {
      closeMenu();
      document.removeEventListener('keydown', esc, true);
    }
  };
  setTimeout(() => {
    document.addEventListener('pointerdown', off, true);
    document.addEventListener('keydown', esc, true);
  }, 0);
}

function buildMenu(items: MenuItem[]): HTMLElement {
  const menu = h('div', { class: 'ctx-menu', role: 'menu' });
  for (const it of items) {
    if (it.separator) {
      menu.appendChild(h('div', { class: 'ctx-sep', role: 'separator' }));
      continue;
    }
    const row = h('button', { class: 'ctx-item', role: 'menuitem', disabled: it.disabled }, h('span', null, it.label), it.keys ? h('kbd', null, it.keys) : null, it.submenu ? h('span', { class: 'ctx-arrow' }, '›') : null);
    if (it.submenu) {
      const sub = buildMenu(it.submenu);
      sub.classList.add('ctx-sub');
      row.appendChild(sub);
      row.addEventListener('click', (e) => {
        e.stopPropagation();
        row.classList.toggle('open');
      });
    } else {
      row.addEventListener('click', () => {
        closeMenu();
        it.run?.();
      });
    }
    menu.appendChild(row);
  }
  return menu;
}
