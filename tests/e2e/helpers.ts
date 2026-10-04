import { Page, expect } from '@playwright/test';

/** Opens the app with a clean slate (no autosave, no welcome dialog). */
export async function openApp(page: Page): Promise<void> {
  await page.addInitScript(() => {
    localStorage.clear();
    localStorage.setItem('chirally:welcomed', '1');
  });
  await page.goto('/');
  await expect(page.locator('canvas.chirally-canvas')).toBeVisible();
}

/** Screen position of a model coordinate. */
export async function screenOf(page: Page, x: number, y: number): Promise<{ x: number; y: number }> {
  return page.evaluate(([mx, my]) => {
    const app = (window as any).chirally;
    const p = app.editor.toScreen({ x: mx, y: my });
    const r = app.editor.canvas.getBoundingClientRect();
    return { x: r.left + p.x, y: r.top + p.y };
  }, [x, y]);
}

export async function atomCount(page: Page): Promise<number> {
  return page.evaluate(() => (window as any).chirally.editor.doc.atoms.size);
}

export async function smiles(page: Page): Promise<string> {
  return page.evaluate(() => (window as any).chirally.currentSmiles());
}

export async function atomScreen(page: Page, index: number): Promise<{ x: number; y: number }> {
  return page.evaluate((i) => {
    const app = (window as any).chirally;
    const a = [...app.editor.doc.atoms.values()][i];
    const p = app.editor.toScreen(a);
    const r = app.editor.canvas.getBoundingClientRect();
    return { x: r.left + p.x, y: r.top + p.y };
  }, index);
}
