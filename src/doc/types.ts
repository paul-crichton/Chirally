// Document model: everything that can live on a Chirally canvas.
// Coordinates are in model units (1 = standard bond length), y DOWN.
import type { Atom, Bond } from '../chem/mol';

/** Atom on the canvas. `id` is unique within the document. */
export type DocAtom = Atom;
/** Bond on the canvas. `a`/`b` are atom IDs (not indices). */
export type DocBond = Bond;

export type ArrowKind =
  | 'reaction'      // →
  | 'equilibrium'   // ⇌ (two half-headed arrows)
  | 'unbalancedEq'  // ⇌ with longer forward arrow
  | 'retro'         // ⇒ retrosynthetic (open double-shaft arrow)
  | 'resonance'     // ↔
  | 'dashed'        // ⇢ hypothetical / multistep
  | 'noGo'          // → crossed out
  | 'line';         // plain line, no heads

export interface ArrowObj {
  id: number;
  type: 'arrow';
  kind: ArrowKind;
  x1: number;
  y1: number;
  x2: number;
  y2: number;
  /** Reagents/conditions written above and below the arrow (multi-line allowed with \n). */
  above?: string;
  below?: string;
  color?: string;
}

/** Where a curved (electron-pushing) arrow starts or ends. */
export type Anchor =
  | { type: 'atom'; id: number }                       // lone pair / atom centre
  | { type: 'bond'; id: number }                       // bonding pair (bond midpoint)
  | { type: 'between'; a: number; b: number }          // new bond between two atoms (target only)
  | { type: 'point'; x: number; y: number };           // free point

export interface CurvedArrowObj {
  id: number;
  type: 'curved';
  /** 2 = full head (electron pair), 1 = fishhook (single electron). */
  electrons: 1 | 2;
  from: Anchor;
  to: Anchor;
  /**
   * Control-point offsets relative to the straight line between the endpoints, expressed in a frame
   * where u = unit vector from start to end and v = u rotated +90°: cp = start + u*t*len + v*h*len.
   * Two control points for a cubic Bézier.
   */
  c1: { t: number; h: number };
  c2: { t: number; h: number };
  color?: string;
}

export interface TextObj {
  id: number;
  type: 'text';
  x: number;
  y: number;
  /** Text with light markup: _x / _{..} subscript, ^x / ^{..} superscript, **bold**, *italic*. */
  text: string;
  /** Font size relative to the document label size (1 = same as atom labels). */
  size?: number;
  align?: 'left' | 'center' | 'right';
  /** Auto-format as a chemical formula (digits subscript, trailing charges superscript). */
  formula?: boolean;
  bold?: boolean;
  italic?: boolean;
  color?: string;
}

export type ShapeKind =
  | 'rect'
  | 'roundRect'
  | 'ellipse'
  | 'line'
  | 'bracket'      // [ ] pair around the box
  | 'paren'        // ( ) pair
  | 'brace'        // { } pair
  | 'tsBracket'    // [ ]‡ transition-state brackets
  | 'orbitalP'     // p orbital (two lobes) centred on box
  | 'orbitalS';    // s orbital (circle)

export interface ShapeObj {
  id: number;
  type: 'shape';
  kind: ShapeKind;
  x1: number;
  y1: number;
  x2: number;
  y2: number;
  /** Optional label, e.g. '‡', 'n', a charge */
  label?: string;
  dashed?: boolean;
  color?: string;
  fill?: string;
  /** Rotation in radians around the centre (orbitals) */
  angle?: number;
}

export interface DocStyle {
  name: string;
  /** Physical bond length for export, in points (ACS 1996: 14.4pt). */
  bondLengthPt: number;
  /** All following sizes are fractions of the bond length. */
  lineWidth: number;
  boldWidth: number;
  /** Distance between double-bond lines. */
  bondSpacing: number;
  hashSpacing: number;
  /** Gap left around atom labels. */
  labelMargin: number;
  fontSize: number;
  fontFamily: string;
  /** Wedge wide-end width */
  wedgeWidth: number;
  colorAtoms: boolean;
  /** Show carbon labels: 'none' (skeletal), 'terminal' (CH3 at chain ends), 'all'. */
  showCarbons: 'none' | 'terminal' | 'all';
  showImplicitH: boolean;
  showStereoLabels: boolean;
  showLonePairs: boolean;
  showAtomNumbers: boolean;
  /** Draw aromatic rings with an inner circle instead of Kekulé double bonds. */
  aromaticCircles: boolean;
  /** Mark atoms with valence errors in red. */
  showValenceErrors: boolean;
}

export interface ChemDoc {
  version: 1;
  atoms: Map<number, DocAtom>;
  bonds: Map<number, DocBond>;
  arrows: Map<number, ArrowObj>;
  curved: Map<number, CurvedArrowObj>;
  texts: Map<number, TextObj>;
  shapes: Map<number, ShapeObj>;
  nextId: number;
  style: DocStyle;
  meta: { title: string; created?: string; modified?: string; notes?: string };
}

export type DocObject = ArrowObj | CurvedArrowObj | TextObj | ShapeObj;

export const STYLE_PRESETS: Record<string, DocStyle> = {
  Chirally: {
    name: 'Chirally',
    bondLengthPt: 18,
    lineWidth: 0.055,
    boldWidth: 0.16,
    bondSpacing: 0.18,
    hashSpacing: 0.13,
    labelMargin: 0.1,
    fontSize: 0.56,
    fontFamily: 'Arial, Helvetica, sans-serif',
    wedgeWidth: 0.2,
    colorAtoms: true,
    showCarbons: 'none',
    showImplicitH: true,
    showStereoLabels: false,
    showLonePairs: false,
    showAtomNumbers: false,
    aromaticCircles: false,
    showValenceErrors: true,
  },
  'ACS 1996': {
    name: 'ACS 1996',
    bondLengthPt: 14.4,
    lineWidth: 0.6 / 14.4,
    boldWidth: 2 / 14.4,
    bondSpacing: 0.18,
    hashSpacing: 2.5 / 14.4,
    labelMargin: 1.6 / 14.4,
    fontSize: 10 / 14.4,
    fontFamily: 'Arial, Helvetica, sans-serif',
    wedgeWidth: 0.2,
    colorAtoms: false,
    showCarbons: 'none',
    showImplicitH: true,
    showStereoLabels: false,
    showLonePairs: false,
    showAtomNumbers: false,
    aromaticCircles: false,
    showValenceErrors: true,
  },
  RSC: {
    name: 'RSC',
    bondLengthPt: 14.4,
    lineWidth: 0.6 / 14.4,
    boldWidth: 2 / 14.4,
    bondSpacing: 0.2,
    hashSpacing: 2.7 / 14.4,
    labelMargin: 1.6 / 14.4,
    fontSize: 7 / 14.4 * 1.3,
    fontFamily: 'Helvetica, Arial, sans-serif',
    wedgeWidth: 0.22,
    colorAtoms: false,
    showCarbons: 'none',
    showImplicitH: true,
    showStereoLabels: false,
    showLonePairs: false,
    showAtomNumbers: false,
    aromaticCircles: false,
    showValenceErrors: true,
  },
  Presentation: {
    name: 'Presentation',
    bondLengthPt: 30,
    lineWidth: 0.075,
    boldWidth: 0.2,
    bondSpacing: 0.2,
    hashSpacing: 0.15,
    labelMargin: 0.1,
    fontSize: 0.6,
    fontFamily: 'Arial, Helvetica, sans-serif',
    wedgeWidth: 0.24,
    colorAtoms: true,
    showCarbons: 'none',
    showImplicitH: true,
    showStereoLabels: true,
    showLonePairs: false,
    showAtomNumbers: false,
    aromaticCircles: false,
    showValenceErrors: true,
  },
};
