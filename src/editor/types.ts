import type { BondStyle } from '../chem/mol';
import type { ArrowKind, ShapeKind } from '../doc/types';
import type { Pt } from '../render/geom';
import type { Mol } from '../chem/mol';

export type ToolId =
  | 'select'
  | 'lasso'
  | 'erase'
  | 'bond'
  | 'chain'
  | 'ring'
  | 'template'
  | 'atom'
  | 'charge'
  | 'radical'
  | 'lonepair'
  | 'text'
  | 'arrow'
  | 'curved'
  | 'shape'
  | 'plus'
  | 'pan';

export type RingKind = { size: number; aromatic: boolean } | { chair: true };

export interface ToolSettings {
  bond: { order: number; style: BondStyle };
  ring: RingKind;
  template: { name: string; mol: Mol } | null;
  atom: { el: string; label?: string };
  charge: 1 | -1;
  arrow: ArrowKind;
  curved: 1 | 2;
  shape: ShapeKind;
}

export type Hit =
  | { kind: 'atom'; id: number }
  | { kind: 'bond'; id: number }
  | { kind: 'arrow'; id: number; part: 'start' | 'end' | 'body' }
  | { kind: 'curved'; id: number; part: 'c1' | 'c2' | 'body' }
  | { kind: 'text'; id: number }
  | { kind: 'shape'; id: number; part: 'body' | 'p1' | 'p2' }
  | { kind: 'rotate' }
  | null;

export interface PEvent {
  /** model coordinates */
  x: number;
  y: number;
  /** screen (CSS px, relative to canvas) */
  sx: number;
  sy: number;
  shift: boolean;
  alt: boolean;
  mod: boolean; // ctrl or cmd
  button: number;
  touch: boolean;
  hit: Hit;
  /** pointer down position (model) during drags */
  start?: Pt;
}

export interface Tool {
  id: ToolId;
  cursor?: string;
  down?(e: PEvent): void;
  drag?(e: PEvent): void;
  up?(e: PEvent): void;
  hover?(e: PEvent): void;
  dblclick?(e: PEvent): void;
  cancel?(): void;
  /** draw previews in screen space */
  overlay?(ctx: CanvasRenderingContext2D): void;
  /** status-bar hint */
  hint?(): string;
}
