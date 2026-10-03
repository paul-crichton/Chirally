// Keyboard handling: ChemDraw/Ketcher-style hover hotkeys plus global tool shortcuts.
import type { Editor } from './editor';
import { applyBondType, idealBondDirection, bondToPoint, attachRingToAtom, fuseRingOnBond, selectedAtomIds, emptySelection, setElement } from './ops';
import { add } from '../render/geom';
import { BondStyle } from '../chem/mol';

/** Element hotkeys (hover an atom or select atoms, then press). */
export const ELEMENT_KEYS: Record<string, string> = {
  c: 'C', n: 'N', o: 'O', s: 'S', p: 'P', f: 'F', h: 'H', i: 'I', l: 'Cl', b: 'Br', d: 'D',
  C: 'Cl', B: 'B', S: 'Si', N: 'Na', K: 'K', L: 'Li', M: 'Mg', Z: 'Zn',
};

/** Abbreviation hotkeys on hovered atoms. */
export const LABEL_KEYS: Record<string, string> = { m: 'Me', e: 'Et', t: 'tBu', a: 'Ac', x: 'Ph', P: 'Ph', A: 'OAc', E: 'CO2Et', O: 'OMe', T: 'TMS' };

const BOND_KEYS: Record<string, { order: number; style: BondStyle }> = {
  '1': { order: 1, style: 'plain' },
  '2': { order: 2, style: 'plain' },
  '3': { order: 3, style: 'plain' },
  '4': { order: 1.5, style: 'plain' },
  w: { order: 1, style: 'wedge' },
  h: { order: 1, style: 'hash' },
  y: { order: 1, style: 'wavy' },
  b: { order: 1, style: 'bold' },
  d: { order: 1, style: 'dashed' },
  a: { order: 1.5, style: 'plain' },
  x: { order: 2, style: 'crossed' },
  k: { order: 1, style: 'dative' },
};

export function handleEditorKey(ed: Editor, ev: KeyboardEvent): boolean {
  if (!ed.isCanvasContext(ev)) return false;
  const k = ev.key;
  const mod = ev.ctrlKey || ev.metaKey;

  if (mod) {
    switch (k.toLowerCase()) {
      case 'z':
        if (ev.shiftKey) ed.redo();
        else ed.undo();
        return true;
      case 'y':
        ed.redo();
        return true;
      case 'a':
        ed.selectAll();
        return true;
      case '0':
        ed.fitToContent();
        return true;
      case '=':
      case '+':
        ed.zoomBy(1.25);
        return true;
      case '-':
        ed.zoomBy(0.8);
        return true;
      case '1':
        ed.setZoomPercent(100);
        return true;
    }
    return false;
  }
  if (ev.altKey) return false;

  // arrow-key nudging
  if (k.startsWith('Arrow') && ed.hasSelection()) {
    const s = ev.shiftKey ? 0.5 : 0.05;
    const d = { ArrowLeft: [-s, 0], ArrowRight: [s, 0], ArrowUp: [0, -s], ArrowDown: [0, s] }[k] as [number, number];
    ed.nudgeSelection(d[0], d[1]);
    return true;
  }

  const h = ed.hover;
  // ───── hovered atom ─────
  if (h?.kind === 'atom' && ed.doc.atoms.has(h.id)) {
    const id = h.id;
    if (k === 'Delete' || k === 'Backspace') {
      const s = emptySelection();
      s.atoms.add(id);
      ed.setSelection(s);
      ed.deleteSelected();
      ed.hover = null;
      return true;
    }
    if (k === 'Enter' || k === 'F2') {
      ed.editAtomLabel(id);
      return true;
    }
    if (k === '1' || k === '2' || k === '3') {
      const order = +k;
      ed.mutate('Add bond', (d) => {
        const a = d.atoms.get(id)!;
        bondToPoint(d, id, add(a, idealBondDirection(d, id, order)), order, 'plain');
      });
      return true;
    }
    if (k === '+' || k === '=') {
      ed.mutate('Charge', (d) => (d.atoms.get(id)!.charge += 1));
      return true;
    }
    if (k === '-' || k === '_') {
      ed.mutate('Charge', (d) => (d.atoms.get(id)!.charge -= 1));
      return true;
    }
    if (k === '0') {
      ed.mutate('Reset charge', (d) => (d.atoms.get(id)!.charge = 0));
      return true;
    }
    if (k === '.') {
      ed.mutate('Radical', (d) => {
        const a = d.atoms.get(id)!;
        a.radical = ((a.radical ?? 0) + 1) % 3;
        if (!a.radical) delete a.radical;
      });
      return true;
    }
    if (k === ',') {
      ed.mutate('Lone pairs', (d) => {
        const a = d.atoms.get(id)!;
        a.lonePairs = !a.lonePairs;
        if (!a.lonePairs) delete a.lonePairs;
      });
      return true;
    }
    if (k === 'r') {
      ed.mutate('Add phenyl', (d) => attachRingToAtom(d, id, 6, true));
      return true;
    }
    if (ELEMENT_KEYS[k]) {
      const el = ELEMENT_KEYS[k];
      ed.mutate('Change element', (d) => {
        if (el === 'D') {
          setElement(d, id, 'H');
          d.atoms.get(id)!.isotope = 2;
        } else setElement(d, id, el);
      });
      return true;
    }
    if (LABEL_KEYS[k]) {
      ed.mutate('Set label', () => ed.applyAtomLabel(id, LABEL_KEYS[k]));
      return true;
    }
    if (/^[A-Za-z]$/.test(k)) {
      // any other letter starts label editing with that letter
      ed.editAtomLabel(id, k.toUpperCase());
      return true;
    }
  }

  // ───── hovered bond ─────
  if (h?.kind === 'bond' && ed.doc.bonds.has(h.id)) {
    const id = h.id;
    if (k === 'Delete' || k === 'Backspace') {
      const s = emptySelection();
      s.bonds.add(id);
      ed.setSelection(s);
      ed.deleteSelected();
      ed.hover = null;
      return true;
    }
    if (BOND_KEYS[k]) {
      const t = BOND_KEYS[k];
      ed.mutate('Change bond', (d) => applyBondType(d.bonds.get(id)!, t.order, t.style));
      return true;
    }
    if (k === 'f') {
      ed.mutate('Flip bond', (d) => {
        const b = d.bonds.get(id)!;
        const t = b.a;
        b.a = b.b;
        b.b = t;
      });
      return true;
    }
    if (k === 'r' || k === '5' || k === '6' || k === '7' || k === '8') {
      const n = k === 'r' ? 6 : +k;
      ed.mutate('Fuse ring', (d) => fuseRingOnBond(d, id, n, k === 'r', ed.mouse));
      return true;
    }
    if (k === ' ') return false;
  }

  // ───── selection-wide element change ─────
  const selAtoms = selectedAtomIds(ed.doc, ed.sel);
  if (selAtoms.size && ELEMENT_KEYS[k] && ed.toolId === 'select') {
    const el = ELEMENT_KEYS[k];
    ed.mutate('Change element', (d) => {
      for (const id of selAtoms) setElement(d, id, el === 'D' ? 'H' : el);
    });
    return true;
  }

  // ───── global ─────
  if (k === 'Escape') {
    if (ed.isDragging) {
      ed.tool.cancel?.();
      ed.cancelChange();
    } else if (ed.hasSelection()) ed.clearSelection();
    else ed.setTool('select');
    return true;
  }
  if (k === 'Delete' || k === 'Backspace') {
    if (ed.hasSelection()) {
      ed.deleteSelected();
      return true;
    }
    return false;
  }
  const toolKeys: Record<string, () => void> = {
    v: () => ed.setTool('select'),
    e: () => ed.setTool('erase'),
    '1': () => setBond(ed, 1, 'plain'),
    '2': () => setBond(ed, 2, 'plain'),
    '3': () => setBond(ed, 3, 'plain'),
    w: () => setBond(ed, 1, 'wedge'),
    q: () => setBond(ed, 1, 'hash'),
    g: () => ed.setTool('chain'),
    r: () => {
      ed.settings.ring = { size: 6, aromatic: true };
      ed.setTool('ring');
    },
    t: () => ed.setTool('text'),
    a: () => ed.setTool('arrow'),
    A: () => ed.setTool('curved'),
    '+': () => {
      ed.settings.charge = 1;
      ed.setTool('charge');
    },
    '-': () => {
      ed.settings.charge = -1;
      ed.setTool('charge');
    },
  };
  if (toolKeys[k]) {
    toolKeys[k]();
    return true;
  }
  if (ELEMENT_KEYS[k] && k !== 'd') {
    ed.settings.atom = { el: ELEMENT_KEYS[k] };
    ed.setTool('atom');
    return true;
  }
  return false;
}

function setBond(ed: Editor, order: number, style: BondStyle): void {
  ed.settings.bond = { order, style };
  ed.setTool('bond');
}
