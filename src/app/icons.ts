// Inline SVG icons (24×24, stroke = currentColor).

const svg = (body: string, extra = '') =>
  `<svg viewBox="0 0 24 24" width="22" height="22" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" ${extra}>${body}</svg>`;

const poly = (n: number, r = 8, cx = 12, cy = 12.5, start = -Math.PI / 2) => {
  const pts: string[] = [];
  for (let k = 0; k < n; k++) {
    const a = start + (2 * Math.PI * k) / n;
    pts.push(`${(cx + r * Math.cos(a)).toFixed(2)},${(cy + r * Math.sin(a)).toFixed(2)}`);
  }
  return pts;
};

export const ICONS: Record<string, string> = {
  select: svg('<path d="M5 3l12 8-5.5 1.2L9 18z" fill="currentColor" stroke="none"/><path d="M12 13l4 6"/>'),
  lasso: svg('<path d="M6 15c-3-3 0-9 6-9s9 3 7 6-7 4-10 3" stroke-dasharray="3 2"/><path d="M9 15l-2 5"/>'),
  erase: svg('<path d="M4 15l8-8 7 7-5 5H8z"/><path d="M9 11l6 6"/><path d="M14 19h6"/>'),
  pan: svg('<path d="M8 12V6a1.5 1.5 0 013 0v5V4.5a1.5 1.5 0 013 0V11V6a1.5 1.5 0 013 0v8c0 4-2.5 6-6 6s-5-2-7-5l-1.5-3a1.5 1.5 0 012.6-1.5L8 13"/>'),
  single: svg('<path d="M5 18L19 6"/>'),
  double: svg('<path d="M4 16L16 5"/><path d="M8 19L20 8"/>'),
  triple: svg('<path d="M3.5 14.5L14 4.5"/><path d="M6 18L18 7"/><path d="M9.5 20.5L20 10.5"/>'),
  aromatic: svg('<path d="M4 16L16 5"/><path d="M8 19L20 8" stroke-dasharray="2.5 2"/>'),
  wedge: svg('<path d="M5 18L19 4.5l1.5 2.8z" fill="currentColor"/>'),
  hash: svg('<path d="M6.5 16.5l.6.7M9 13.5l1.4 1.5M11.5 10.5l2 2.2M14 7.5l2.8 3M16.5 4.7l3.5 3.7"/>'),
  hollow: svg('<path d="M5 18L19 4.5l1.5 2.8z"/>'),
  wavy: svg('<path d="M4 18c1.5-3 3-1 4.5-3.5S10 11 12 10s2-3.5 4-4.5 2.5-.5 4-2"/>'),
  bold: svg('<path d="M5 18L19 6" stroke-width="4.5"/>'),
  dashed: svg('<path d="M5 18L19 6" stroke-dasharray="2.5 2.5"/>'),
  dative: svg('<path d="M4 18L17 7"/><path d="M13 6.5l5.5-1-1.2 5.5z" fill="currentColor"/>'),
  hbond: svg('<path d="M5 18L19 6" stroke-dasharray="0.5 3" stroke-width="2.4"/>'),
  crossed: svg('<path d="M4 15L20 9"/><path d="M4 9L20 15"/>'),
  chain: svg('<path d="M3 15l4.5-6 4.5 6 4.5-6 4.5 6"/>'),
  ring3: svg(`<polygon points="${poly(3, 8, 12, 14).join(' ')}"/>`),
  ring4: svg(`<polygon points="${poly(4, 8.5, 12, 12, -Math.PI / 4).join(' ')}"/>`),
  ring5: svg(`<polygon points="${poly(5).join(' ')}"/>`),
  ring6: svg(`<polygon points="${poly(6).join(' ')}"/>`),
  ring7: svg(`<polygon points="${poly(7).join(' ')}"/>`),
  ring8: svg(`<polygon points="${poly(8, 8.5).join(' ')}"/>`),
  benzene: svg(`<polygon points="${poly(6).join(' ')}"/><path d="M12 7.5L16.3 10M16.3 15L12 17.5M7.7 15V10"/>`.replace('M12 7.5L16.3 10M16.3 15L12 17.5M7.7 15V10', inner6())),
  chair: svg('<path d="M2.5 14l4-4.5 6 2 3.5-3.5 5.5 1.5-4 4.5-6-2-3.5 3.5z"/>'),
  cyclopentadiene: svg(`<polygon points="${poly(5).join(' ')}"/>`),
  template: svg('<path d="M4 5h6v14H4zM10 5h4v14h-4zM14 6l4-1 2.5 13.5-4 1z"/>'),
  atom: svg('<text x="12" y="17" text-anchor="middle" font-size="14" font-weight="700" fill="currentColor" stroke="none" font-family="system-ui,sans-serif">N</text>'),
  periodic: svg('<rect x="3" y="4" width="4" height="4"/><rect x="17" y="4" width="4" height="4"/><rect x="3" y="10" width="18" height="4"/><rect x="6" y="16" width="12" height="3"/>'),
  chargePlus: svg('<circle cx="12" cy="12" r="8"/><path d="M12 8v8M8 12h8"/>'),
  chargeMinus: svg('<circle cx="12" cy="12" r="8"/><path d="M8 12h8"/>'),
  radical: svg('<text x="10" y="18" text-anchor="middle" font-size="13" font-weight="700" fill="currentColor" stroke="none" font-family="system-ui,sans-serif">C</text><circle cx="18" cy="7" r="1.9" fill="currentColor" stroke="none"/>'),
  lonepair: svg('<text x="12" y="18" text-anchor="middle" font-size="13" font-weight="700" fill="currentColor" stroke="none" font-family="system-ui,sans-serif">O</text><circle cx="9.5" cy="4.5" r="1.4" fill="currentColor" stroke="none"/><circle cx="14.5" cy="4.5" r="1.4" fill="currentColor" stroke="none"/>'),
  text: svg('<path d="M5 6V4h14v2M12 4v16M9 20h6"/>'),
  plus: svg('<path d="M12 5v14M5 12h14" stroke-width="2.2"/>'),
  reaction: svg('<path d="M3 12h15"/><path d="M15 8l5 4-5 4z" fill="currentColor"/>'),
  equilibrium: svg('<path d="M3 10h16l-4-3.5"/><path d="M21 14H5l4 3.5"/>'),
  unbalancedEq: svg('<path d="M3 10h16l-4-3.5"/><path d="M17 14H7l4 3.5"/>'),
  retro: svg('<path d="M3 9.5h14M3 14.5h14"/><path d="M15 5.5l5.5 6.5-5.5 6.5"/>'),
  resonance: svg('<path d="M6 12h12"/><path d="M8 8l-5 4 5 4z M16 8l5 4-5 4z" fill="currentColor"/>'),
  dashedArrow: svg('<path d="M3 12h14" stroke-dasharray="3 2.5"/><path d="M15 8l5 4-5 4z" fill="currentColor"/>'),
  noGo: svg('<path d="M3 12h15"/><path d="M15 8l5 4-5 4z" fill="currentColor"/><path d="M8 8l4 8M12 8l-4 8"/>'),
  line: svg('<path d="M4 12h16"/>'),
  curved2: svg('<path d="M4 17C5 6 15 4 18 10"/><path d="M14.5 9.5l4.5 2 .5-5z" fill="currentColor"/>'),
  curved1: svg('<path d="M4 17C5 6 15 4 18 10"/><path d="M14.5 9.5l4.5 2 .3-3z" fill="currentColor"/>'),
  bracket: svg('<path d="M7 4H4v16h3M17 4h3v16h-3"/>'),
  tsBracket: svg('<path d="M6 5H3.5v15H6M15 5h2.5v15H15"/><text x="21" y="9" text-anchor="middle" font-size="7" fill="currentColor" stroke="none">‡</text>'),
  paren: svg('<path d="M7 4c-3 4-3 12 0 16M17 4c3 4 3 12 0 16"/>'),
  brace: svg('<path d="M8 4c-2 0-2 1-2 3v3c0 1-1 2-2 2 1 0 2 1 2 2v3c0 2 0 3 2 3M16 4c2 0 2 1 2 3v3c0 1 1 2 2 2-1 0-2 1-2 2v3c0 2 0 3-2 3"/>'),
  rect: svg('<rect x="4" y="6" width="16" height="12"/>'),
  roundRect: svg('<rect x="4" y="6" width="16" height="12" rx="3"/>'),
  ellipse: svg('<ellipse cx="12" cy="12" rx="8.5" ry="6"/>'),
  orbitalP: svg('<path d="M12 12C8 9 9 3 12 3s4 6 0 9z" fill="currentColor" fill-opacity="0.35"/><path d="M12 12c-4 3-3 9 0 9s4-6 0-9z"/>'),
  orbitalS: svg('<circle cx="12" cy="12" r="7"/>'),
  undo: svg('<path d="M9 14L4 9l5-5"/><path d="M4 9h10a6 6 0 010 12h-3"/>'),
  redo: svg('<path d="M15 14l5-5-5-5"/><path d="M20 9H10a6 6 0 000 12h3"/>'),
  open: svg('<path d="M3 7a2 2 0 012-2h4l2 2h8a2 2 0 012 2v8a2 2 0 01-2 2H5a2 2 0 01-2-2z"/>'),
  save: svg('<path d="M5 3h11l3 3v13a2 2 0 01-2 2H5a2 2 0 01-2-2V5a2 2 0 012-2z"/><path d="M7 3v5h8V3M7 21v-7h10v7"/>'),
  newDoc: svg('<path d="M14 3H6a2 2 0 00-2 2v14a2 2 0 002 2h12a2 2 0 002-2V9z"/><path d="M14 3v6h6M12 12v6M9 15h6"/>'),
  export: svg('<path d="M12 3v12M7 8l5-5 5 5"/><path d="M5 15v4a2 2 0 002 2h10a2 2 0 002-2v-4"/>'),
  share: svg('<circle cx="18" cy="5" r="3"/><circle cx="6" cy="12" r="3"/><circle cx="18" cy="19" r="3"/><path d="M8.6 13.5l6.8 4M15.4 6.5l-6.8 4"/>'),
  clean: svg('<path d="M12 3l1.8 4.5L18 9l-4.2 1.5L12 15l-1.8-4.5L6 9l4.2-1.5z"/><path d="M18 15l.9 2.1L21 18l-2.1.9L18 21l-.9-2.1L15 18l2.1-.9z"/>'),
  fit: svg('<path d="M4 9V4h5M20 9V4h-5M4 15v5h5M20 15v5h-5"/>'),
  zoomIn: svg('<circle cx="11" cy="11" r="7"/><path d="M16 16l5 5M11 8v6M8 11h6"/>'),
  zoomOut: svg('<circle cx="11" cy="11" r="7"/><path d="M16 16l5 5M8 11h6"/>'),
  search: svg('<circle cx="11" cy="11" r="7"/><path d="M16 16l5 5"/>'),
  sun: svg('<circle cx="12" cy="12" r="4"/><path d="M12 2v2M12 20v2M4.9 4.9l1.4 1.4M17.7 17.7l1.4 1.4M2 12h2M20 12h2M4.9 19.1l1.4-1.4M17.7 6.3l1.4-1.4"/>'),
  moon: svg('<path d="M20 14.5A8 8 0 019.5 4a8 8 0 1010.5 10.5z"/>'),
  help: svg('<circle cx="12" cy="12" r="9"/><path d="M9.5 9a2.5 2.5 0 015 .5c0 1.5-2.5 2-2.5 3.5M12 17h.01"/>'),
  panel: svg('<rect x="3" y="4" width="18" height="16" rx="2"/><path d="M15 4v16"/>'),
  menu: svg('<path d="M4 6h16M4 12h16M4 18h16"/>'),
  flipH: svg('<path d="M12 3v18" stroke-dasharray="2 2"/><path d="M9 7L4 12l5 5zM15 7l5 5-5 5z"/>'),
  flipV: svg('<path d="M3 12h18" stroke-dasharray="2 2"/><path d="M7 9l5-5 5 5zM7 15l5 5 5-5z"/>'),
  rotate: svg('<path d="M20 11a8 8 0 10-2.3 5.7"/><path d="M20 4v7h-7"/>'),
  copy: svg('<rect x="8" y="8" width="12" height="12" rx="2"/><path d="M16 8V6a2 2 0 00-2-2H6a2 2 0 00-2 2v8a2 2 0 002 2h2"/>'),
  duplicate: svg('<rect x="8" y="8" width="12" height="12" rx="2"/><path d="M14 11v6M11 14h6M16 8V6a2 2 0 00-2-2H6a2 2 0 00-2 2v8a2 2 0 002 2h2"/>'),
  trash: svg('<path d="M4 7h16M10 11v6M14 11v6M5 7l1 13h12l1-13M9 7V4h6v3"/>'),
  cube: svg('<path d="M12 2l9 5v10l-9 5-9-5V7z"/><path d="M12 22V12M21 7l-9 5-9-5"/>'),
  flask: svg('<path d="M9 3h6M10 3v6l-5.5 9.5A2 2 0 006.2 21h11.6a2 2 0 001.7-2.5L14 9V3"/><path d="M7.5 15h9"/>'),
  info: svg('<circle cx="12" cy="12" r="9"/><path d="M12 11v6M12 7.5h.01"/>'),
  mech: svg('<path d="M5 18C6 8 14 5 18 9"/><path d="M15 8.5l4 1.6.6-4.4z" fill="currentColor"/><circle cx="5" cy="19.5" r="1.2" fill="currentColor"/>'),
  style: svg('<circle cx="12" cy="12" r="9"/><circle cx="8" cy="10" r="1.3" fill="currentColor"/><circle cx="12" cy="7" r="1.3" fill="currentColor"/><circle cx="16" cy="10" r="1.3" fill="currentColor"/><path d="M12 21a3 3 0 010-6h2a3 3 0 003-3"/>'),
  library: svg('<path d="M4 19V5a2 2 0 012-2h13v16H6a2 2 0 00-2 2zm0 0a2 2 0 002 2h13"/>'),
  pubchem: svg('<circle cx="10" cy="10" r="6"/><path d="M14.5 14.5L20 20"/><path d="M7.5 10h5M10 7.5v5"/>'),
  expand: svg('<path d="M6 9l6 6 6-6"/>'),
  close: svg('<path d="M6 6l12 12M18 6L6 18"/>'),
  grid: svg('<path d="M4 4h16v16H4zM4 12h16M12 4v16"/>'),
  keyboard: svg('<rect x="2" y="6" width="20" height="12" rx="2"/><path d="M6 10h.01M10 10h.01M14 10h.01M18 10h.01M7 14h10"/>'),
  link: svg('<path d="M10 14a4 4 0 005.7 0l3-3a4 4 0 00-5.7-5.7l-1 1"/><path d="M14 10a4 4 0 00-5.7 0l-3 3a4 4 0 005.7 5.7l1-1"/>'),
  sparkle: svg('<path d="M12 3l2 6 6 2-6 2-2 6-2-6-6-2 6-2z"/>'),
  check: svg('<path d="M5 12l5 5L20 7"/>'),
  warn: svg('<path d="M12 3l10 18H2z"/><path d="M12 10v5M12 18h.01"/>'),
};

function inner6(): string {
  // alternating inner double-bond lines of a hexagon (vertex at top)
  const outer = poly(6).map((p) => p.split(',').map(Number));
  const cx = 12, cy = 12.5;
  const inset = (p: number[]) => [cx + (p[0] - cx) * 0.7, cy + (p[1] - cy) * 0.7];
  let d = '';
  for (const k of [0, 2, 4]) {
    const a = inset(outer[k]), b = inset(outer[(k + 1) % 6]);
    d += `M${a[0].toFixed(2)} ${a[1].toFixed(2)}L${b[0].toFixed(2)} ${b[1].toFixed(2)}`;
  }
  return d;
}

/** Icon for an element symbol, drawn as text (colour applied by CSS). */
export function elementIcon(sym: string): string {
  const size = sym.length > 2 ? 9 : sym.length > 1 ? 11 : 14;
  return svg(`<text x="12" y="${12 + size * 0.36}" text-anchor="middle" font-size="${size}" font-weight="700" fill="currentColor" stroke="none" font-family="system-ui,sans-serif">${sym}</text>`);
}
