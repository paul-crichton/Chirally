import { test, expect, Page } from '@playwright/test';
import { openApp, atomCount } from './helpers';

test.beforeEach(async ({ page }) => {
  await openApp(page);
  await page.locator('.tab', { hasText: 'Mechanism' }).click();
});

/** Hydroxide and acetaldehyde with only the O⁻ → C arrow drawn (the C=O → O arrow is forgotten). */
async function drawIncompleteAddition(page: Page): Promise<void> {
  await page.evaluate(() => {
    const ed = (window as any).chirally.editor;
    ed.mutate('setup', (d: any) => {
      const add = (el: string, x: number, y: number, charge = 0) => {
        const id = d.nextId++;
        d.atoms.set(id, { id, el, x, y, charge });
        return id;
      };
      const bond = (a: number, b: number, order = 1) => {
        const id = d.nextId++;
        d.bonds.set(id, { id, a, b, order, style: 'plain' });
      };
      const o = add('O', 0, 0, -1);
      const c1 = add('C', 3, 0.5), c2 = add('C', 3.87, 0), o2 = add('O', 4.73, 0.5);
      bond(c1, c2);
      bond(c2, o2, 2);
      const id = d.nextId++;
      d.curved.set(id, { id, type: 'curved', electrons: 2, from: { type: 'atom', id: o }, to: { type: 'atom', id: c2 }, c1: { t: 0.25, h: -0.5 }, c2: { t: 0.75, h: -0.5 } });
    });
    ed.fitToContent();
  });
}

const hasGhost = (page: Page) => page.evaluate(() => (window as any).chirally.editor.hasGhost);
const reactionArrows = (page: Page) => page.evaluate(() => (window as any).chirally.editor.doc.arrows.size);

test('an impossible intermediate is previewed, not added', async ({ page }) => {
  await drawIncompleteAddition(page);
  const before = await atomCount(page);
  await page.getByRole('button', { name: /Apply arrows/ }).click();
  expect(await atomCount(page)).toBe(before);
  expect(await reactionArrows(page)).toBe(0);
  expect(await hasGhost(page)).toBe(true);
  await expect(page.locator('.mech-pending')).toBeVisible();
  await expect(page.locator('.mech-warnings .err')).toContainText(['10 valence electrons']);
  await page.screenshot({ path: 'test-results/mechanism-ghost.png' });

  // insert anyway → added as a normal reaction step and remembered as applied
  await page.getByRole('button', { name: 'Insert anyway' }).click();
  expect(await atomCount(page)).toBe(before + 4);
  expect(await hasGhost(page)).toBe(false);
  await expect(page.locator('.mech-pending')).toBeHidden();
  const kind = await page.evaluate(() => [...(window as any).chirally.editor.doc.arrows.values()][0].kind);
  expect(kind).toBe('reaction');
  await page.getByRole('button', { name: /Apply arrows/ }).click();
  await expect(page.locator('.toast', { hasText: 'All curved arrows have been applied' })).toBeVisible();
  expect(await atomCount(page)).toBe(before + 4);
});

test('discarding or editing the drawing clears the preview', async ({ page }) => {
  await drawIncompleteAddition(page);
  await page.getByRole('button', { name: /Apply arrows/ }).click();
  await page.getByRole('button', { name: 'Discard' }).click();
  expect(await hasGhost(page)).toBe(false);
  await expect(page.locator('.mech-pending')).toBeHidden();

  await page.getByRole('button', { name: /Apply arrows/ }).click();
  expect(await hasGhost(page)).toBe(true);
  await page.evaluate(() => (window as any).chirally.editor.mutate('edit', (d: any) => d.atoms.set(d.nextId, { id: d.nextId++, el: 'N', x: 0, y: 4, charge: 0 })));
  expect(await hasGhost(page)).toBe(false);
  await expect(page.locator('.mech-pending')).toBeHidden();
});

test('applied steps are remembered; undo makes them applicable again', async ({ page }) => {
  await page.evaluate(() => (window as any).chirally.loadExample());
  const before = await atomCount(page);
  await page.getByRole('button', { name: /Apply arrows/ }).click();
  const after = await atomCount(page);
  expect(after).toBeGreaterThan(before);
  await expect(page.locator('.mech-status')).toContainText('all applied');

  await page.getByRole('button', { name: /Apply arrows/ }).click();
  expect(await atomCount(page)).toBe(after);
  await expect(page.locator('.toast', { hasText: 'All curved arrows have been applied' })).toBeVisible();

  await page.keyboard.press('Control+z');
  expect(await atomCount(page)).toBe(before);
  await page.getByRole('button', { name: /Apply arrows/ }).click();
  expect(await atomCount(page)).toBe(after);
});

test('arrows can start and end on the H of an OH label', async ({ page }) => {
  await page.evaluate(() => {
    const app = (window as any).chirally;
    const ed = app.editor;
    ed.mutate('setup', (d: any) => {
      const add = (el: string, x: number, y: number, charge = 0) => {
        const id = d.nextId++;
        d.atoms.set(id, { id, el, x, y, charge });
        return id;
      };
      const bond = (a: number, b: number, order = 1) => {
        const id = d.nextId++;
        d.bonds.set(id, { id, a, b, order, style: 'plain' });
      };
      add('O', -1, 0, -1); // hydroxide, left
      const c1 = add('C', 2, 0), c2 = add('C', 2.87, 0.5), o1 = add('O', 2.87, 1.5), o2 = add('O', 3.73, 0);
      bond(c1, c2);
      bond(c2, o1, 2);
      bond(c2, o2); // acid O–H with its H drawn on the right of the label
    });
    ed.fitToContent();
    ed.settings.curved = 2;
    ed.setTool('curved');
  });
  const pts = await page.evaluate(() => {
    const ed = (window as any).chirally.editor;
    const atoms = [...ed.doc.atoms.values()];
    const ho = atoms[0], acidO = atoms[4];
    const box = ed.getScene().labelBoxes.get(acidO.id).find((b: any) => b.role === 'h');
    const r = ed.canvas.getBoundingClientRect();
    const scr = (p: any) => { const s = ed.toScreen(p); return { x: r.left + s.x, y: r.top + s.y }; };
    return { ho: scr(ho), o: scr(acidO), h: scr({ x: (box.x1 + box.x2) / 2, y: (box.y1 + box.y2) / 2 }) };
  });
  const drag = async (a: { x: number; y: number }, b: { x: number; y: number }) => {
    await page.mouse.move(a.x, a.y);
    await page.mouse.down();
    await page.mouse.move((a.x + b.x) / 2, a.y - 30, { steps: 4 });
    await page.mouse.move(b.x, b.y, { steps: 4 });
    await page.mouse.up();
  };
  await drag(pts.ho, pts.h); // base → the H
  await drag(pts.h, pts.o); // the O–H bond → O
  const anchors = await page.evaluate(() => [...(window as any).chirally.editor.doc.curved.values()].map((c: any) => [c.from, c.to]));
  expect(anchors[0][1]).toMatchObject({ type: 'atom', h: true });
  expect(anchors[1][0]).toMatchObject({ type: 'atom', h: true });
  expect(anchors[1][1].h).toBeFalsy();
  await page.getByRole('button', { name: /Apply arrows/ }).click();
  await expect(page.locator('.mech-warnings .err')).toHaveCount(0);
  const smi = await page.evaluate(() => (window as any).chirally.currentSmiles());
  expect(smi.split('.').sort()).toEqual(expect.arrayContaining(['CC([O-])=O', 'O']));
  await page.screenshot({ path: 'test-results/mechanism-h-anchor.png' });
});

