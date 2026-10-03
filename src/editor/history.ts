import { ChemDoc } from '../doc/types';
import { cloneDoc } from '../doc/document';

interface Entry {
  doc: ChemDoc;
  label: string;
}

/** Snapshot-based undo/redo. Documents are small, so whole-document snapshots are simple and robust. */
export class History {
  private undoStack: Entry[] = [];
  private redoStack: Entry[] = [];
  private pending: ChemDoc | null = null;
  limit = 300;

  /** Remember the state before a change starts (idempotent until commit/cancel). */
  begin(doc: ChemDoc): void {
    if (!this.pending) this.pending = cloneDoc(doc);
  }

  get inProgress(): boolean {
    return this.pending !== null;
  }

  /** Records the change started with begin(). */
  commit(label: string): void {
    if (!this.pending) return;
    this.undoStack.push({ doc: this.pending, label });
    if (this.undoStack.length > this.limit) this.undoStack.shift();
    this.redoStack = [];
    this.pending = null;
  }

  /** Abandons a change; returns the pre-change document so the caller can restore it. */
  cancel(): ChemDoc | null {
    const p = this.pending;
    this.pending = null;
    return p;
  }

  undo(current: ChemDoc): ChemDoc | null {
    const e = this.undoStack.pop();
    if (!e) return null;
    this.redoStack.push({ doc: cloneDoc(current), label: e.label });
    return e.doc;
  }

  redo(current: ChemDoc): ChemDoc | null {
    const e = this.redoStack.pop();
    if (!e) return null;
    this.undoStack.push({ doc: cloneDoc(current), label: e.label });
    return e.doc;
  }

  get canUndo(): boolean {
    return this.undoStack.length > 0;
  }
  get canRedo(): boolean {
    return this.redoStack.length > 0;
  }
  get undoLabel(): string | null {
    return this.undoStack[this.undoStack.length - 1]?.label ?? null;
  }
  get redoLabel(): string | null {
    return this.redoStack[this.redoStack.length - 1]?.label ?? null;
  }
  clear(): void {
    this.undoStack = [];
    this.redoStack = [];
    this.pending = null;
  }
}
