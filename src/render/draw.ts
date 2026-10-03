// Primitive renderers: Canvas 2D (interactive view, PNG) and SVG (vector export).
import { Prim, PathSeg } from './scene';
import { fontString } from './text';
import { Box } from './geom';

export interface ViewTransform {
  scale: number; // px per model unit
  ox: number;
  oy: number;
}

function tracePath(ctx: CanvasRenderingContext2D, d: PathSeg[], T: ViewTransform): void {
  const X = (x: number) => x * T.scale + T.ox;
  const Y = (y: number) => y * T.scale + T.oy;
  ctx.beginPath();
  for (const s of d) {
    switch (s[0]) {
      case 'M': ctx.moveTo(X(s[1]), Y(s[2])); break;
      case 'L': ctx.lineTo(X(s[1]), Y(s[2])); break;
      case 'Q': ctx.quadraticCurveTo(X(s[1]), Y(s[2]), X(s[3]), Y(s[4])); break;
      case 'C': ctx.bezierCurveTo(X(s[1]), Y(s[2]), X(s[3]), Y(s[4]), X(s[5]), Y(s[6])); break;
      case 'Z': ctx.closePath(); break;
    }
  }
}

/** Draws primitives on a canvas context. `minLine` clamps hairlines when zoomed out. */
export function drawPrims(ctx: CanvasRenderingContext2D, prims: Prim[], T: ViewTransform, minLine = 0.6): void {
  const X = (x: number) => x * T.scale + T.ox;
  const Y = (y: number) => y * T.scale + T.oy;
  const W = (w: number) => Math.max(minLine, w * T.scale);
  ctx.lineJoin = 'round';
  for (const p of prims) {
    switch (p.k) {
      case 'line': {
        ctx.beginPath();
        ctx.moveTo(X(p.x1), Y(p.y1));
        ctx.lineTo(X(p.x2), Y(p.y2));
        ctx.strokeStyle = p.color;
        ctx.lineWidth = W(p.w);
        ctx.lineCap = p.cap ?? 'round';
        ctx.setLineDash(p.dash ? p.dash.map((v) => v * T.scale) : []);
        ctx.stroke();
        break;
      }
      case 'poly': {
        ctx.beginPath();
        for (let i = 0; i < p.pts.length; i += 2) {
          if (i === 0) ctx.moveTo(X(p.pts[0]), Y(p.pts[1]));
          else ctx.lineTo(X(p.pts[i]), Y(p.pts[i + 1]));
        }
        if (p.closed) ctx.closePath();
        if (p.fill) {
          ctx.fillStyle = p.fill;
          ctx.fill();
        }
        if (p.stroke) {
          ctx.strokeStyle = p.stroke;
          ctx.lineWidth = W(p.w ?? 0.05);
          ctx.lineCap = 'round';
          ctx.setLineDash(p.dash ? p.dash.map((v) => v * T.scale) : []);
          ctx.stroke();
        }
        break;
      }
      case 'path': {
        tracePath(ctx, p.d, T);
        if (p.fill) {
          ctx.fillStyle = p.fill;
          ctx.fill();
        }
        if (p.stroke) {
          ctx.strokeStyle = p.stroke;
          ctx.lineWidth = W(p.w ?? 0.05);
          ctx.lineCap = p.cap ?? 'round';
          ctx.setLineDash(p.dash ? p.dash.map((v) => v * T.scale) : []);
          ctx.stroke();
        }
        break;
      }
      case 'circle': {
        ctx.beginPath();
        ctx.arc(X(p.x), Y(p.y), Math.max(0.5, p.r * T.scale), 0, Math.PI * 2);
        if (p.fill) {
          ctx.fillStyle = p.fill;
          ctx.fill();
        }
        if (p.stroke) {
          ctx.strokeStyle = p.stroke;
          ctx.lineWidth = W(p.w ?? 0.05);
          ctx.setLineDash(p.dash ? p.dash.map((v) => v * T.scale) : []);
          ctx.stroke();
        }
        break;
      }
      case 'text': {
        ctx.font = fontString(p.size * T.scale, p.family, p.bold, p.italic);
        ctx.fillStyle = p.color;
        ctx.textBaseline = 'alphabetic';
        ctx.textAlign = 'left';
        ctx.fillText(p.text, X(p.x), Y(p.y));
        break;
      }
    }
  }
  ctx.setLineDash([]);
}

const esc = (s: string) => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
const f = (v: number) => (Math.round(v * 100) / 100).toString();

/**
 * Serializes primitives to a standalone SVG. `unit` = output units (pt) per model unit.
 * `background` null = transparent.
 */
export function primsToSVG(prims: Prim[], bounds: Box, unit: number, opts: { padding?: number; background?: string | null; title?: string } = {}): string {
  const pad = opts.padding ?? 0.5;
  const x0 = bounds.x1 - pad, y0 = bounds.y1 - pad;
  const Wm = bounds.x2 - bounds.x1 + 2 * pad, Hm = bounds.y2 - bounds.y1 + 2 * pad;
  const X = (x: number) => f((x - x0) * unit);
  const Y = (y: number) => f((y - y0) * unit);
  const S = (v: number) => f(v * unit);
  const out: string[] = [];
  out.push(`<svg xmlns="http://www.w3.org/2000/svg" width="${f(Wm * unit)}pt" height="${f(Hm * unit)}pt" viewBox="0 0 ${f(Wm * unit)} ${f(Hm * unit)}">`);
  if (opts.title) out.push(`<title>${esc(opts.title)}</title>`);
  if (opts.background) out.push(`<rect width="100%" height="100%" fill="${opts.background}"/>`);
  const dashAttr = (d?: number[]) => (d ? ` stroke-dasharray="${d.map((v) => S(v)).join(' ')}"` : '');
  for (const p of prims) {
    switch (p.k) {
      case 'line':
        out.push(`<line x1="${X(p.x1)}" y1="${Y(p.y1)}" x2="${X(p.x2)}" y2="${Y(p.y2)}" stroke="${p.color}" stroke-width="${S(p.w)}" stroke-linecap="${p.cap ?? 'round'}"${dashAttr(p.dash)}/>`);
        break;
      case 'poly': {
        const pts: string[] = [];
        for (let i = 0; i < p.pts.length; i += 2) pts.push(`${X(p.pts[i])},${Y(p.pts[i + 1])}`);
        const tag = p.closed ? 'polygon' : 'polyline';
        out.push(`<${tag} points="${pts.join(' ')}" fill="${p.fill ?? 'none'}"${p.stroke ? ` stroke="${p.stroke}" stroke-width="${S(p.w ?? 0.05)}" stroke-linejoin="round" stroke-linecap="round"` : ''}${dashAttr(p.dash)}/>`);
        break;
      }
      case 'path': {
        const d = p.d
          .map((s) => {
            switch (s[0]) {
              case 'M': return `M${X(s[1])} ${Y(s[2])}`;
              case 'L': return `L${X(s[1])} ${Y(s[2])}`;
              case 'Q': return `Q${X(s[1])} ${Y(s[2])} ${X(s[3])} ${Y(s[4])}`;
              case 'C': return `C${X(s[1])} ${Y(s[2])} ${X(s[3])} ${Y(s[4])} ${X(s[5])} ${Y(s[6])}`;
              case 'Z': return 'Z';
            }
            return '';
          })
          .join(' ');
        out.push(`<path d="${d}" fill="${p.fill ?? 'none'}"${p.stroke ? ` stroke="${p.stroke}" stroke-width="${S(p.w ?? 0.05)}" stroke-linecap="${p.cap ?? 'round'}" stroke-linejoin="round"` : ''}${dashAttr(p.dash)}/>`);
        break;
      }
      case 'circle':
        out.push(`<circle cx="${X(p.x)}" cy="${Y(p.y)}" r="${S(p.r)}" fill="${p.fill ?? 'none'}"${p.stroke ? ` stroke="${p.stroke}" stroke-width="${S(p.w ?? 0.05)}"` : ''}${dashAttr(p.dash)}/>`);
        break;
      case 'text':
        out.push(`<text x="${X(p.x)}" y="${Y(p.y)}" font-family="${esc(p.family)}" font-size="${S(p.size)}" fill="${p.color}"${p.bold ? ' font-weight="bold"' : ''}${p.italic ? ' font-style="italic"' : ''}>${esc(p.text)}</text>`);
        break;
    }
  }
  out.push('</svg>');
  return out.join('\n');
}

/** Renders primitives into a new canvas (for PNG export / thumbnails). */
export function primsToCanvas(prims: Prim[], bounds: Box, pxPerUnit: number, opts: { padding?: number; background?: string | null; maxSize?: number } = {}): HTMLCanvasElement {
  const pad = opts.padding ?? 0.5;
  let scale = pxPerUnit;
  const Wm = bounds.x2 - bounds.x1 + 2 * pad, Hm = bounds.y2 - bounds.y1 + 2 * pad;
  if (opts.maxSize) scale = Math.min(scale, opts.maxSize / Math.max(Wm, Hm));
  const canvas = document.createElement('canvas');
  canvas.width = Math.max(1, Math.ceil(Wm * scale));
  canvas.height = Math.max(1, Math.ceil(Hm * scale));
  const ctx = canvas.getContext('2d')!;
  if (opts.background) {
    ctx.fillStyle = opts.background;
    ctx.fillRect(0, 0, canvas.width, canvas.height);
  }
  drawPrims(ctx, prims, { scale, ox: (-bounds.x1 + pad) * scale, oy: (-bounds.y1 + pad) * scale }, 0.4);
  return canvas;
}
