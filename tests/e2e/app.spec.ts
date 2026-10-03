import { test, expect } from '@playwright/test';
import { openApp, atomCount, smiles, atomScreen } from './helpers';

test.beforeEach(async ({ page }) => {
  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e)));
  (page as any)._errors = errors;
  await openApp(page);
});

test.afterEach(async ({ page }) => {
  expect((page as any)._errors, 'uncaught page errors').toEqual([]);
});

test('draws a bond with a click and changes an atom with a hotkey', async ({ page }) => {
  const canvas = page.locator('canvas.cw-canvas');
  const box = (await canvas.boundingBox())!;
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
  expect(await atomCount(page)).toBe(2);
  // sprout a third atom by clicking the second atom
  const a1 = await atomScreen(page, 1);
  await page.mouse.click(a1.x, a1.y);
  expect(await atomCount(page)).toBe(3);
  // hover last atom and press "o"
  const a2 = await atomScreen(page, 2);
  await page.mouse.move(a2.x, a2.y);
  await page.keyboard.press('o');
  expect(await smiles(page)).toBe('CCO');
  await expect(page.locator('.iupac')).toHaveText('ethanol');
  await page.keyboard.press('Control+z');
  expect(await smiles(page)).toBe('CCC');
  await page.keyboard.press('Control+Shift+z');
  expect(await smiles(page)).toBe('CCO');
});

test('ring tool places benzene and analysis shows formula', async ({ page }) => {
  await page.keyboard.press('r');
  const box = (await page.locator('canvas.cw-canvas').boundingBox())!;
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
  expect(await atomCount(page)).toBe(6);
  await expect(page.locator('.iupac')).toHaveText('benzene');
  await expect(page.locator('.formula').first()).toContainText('C6H6');
});

test('command palette inserts SMILES and names it', async ({ page }) => {
  await page.keyboard.press('Control+k');
  await page.locator('.palette-input').fill('CC(=O)Oc1ccccc1C(=O)O');
  await page.keyboard.press('Enter');
  expect(await atomCount(page)).toBe(13);
  await expect(page.locator('.iupac')).toHaveText('2-acetyloxybenzoic acid');
  await expect(page.locator('.formula').first()).toContainText('C9H8O4');
});

test('drag-drawing snaps bond angles and joins atoms', async ({ page }) => {
  const box = (await page.locator('canvas.cw-canvas').boundingBox())!;
  const cx = box.x + box.width / 2, cy = box.y + box.height / 2;
  await page.mouse.move(cx, cy);
  await page.mouse.down();
  await page.mouse.move(cx + 30, cy + 3, { steps: 4 });
  await page.mouse.move(cx + 60, cy + 4, { steps: 4 });
  await page.mouse.up();
  const angle = await page.evaluate(() => {
    const atoms = [...(window as any).chemwrite.editor.doc.atoms.values()];
    return (Math.atan2(atoms[1].y - atoms[0].y, atoms[1].x - atoms[0].x) * 180) / Math.PI;
  });
  expect(Math.abs(angle % 15)).toBeLessThan(0.01);
});

test('mechanism example: apply arrows generates the SN2 product', async ({ page }) => {
  await page.evaluate(() => (window as any).chemwrite.loadExample());
  await page.locator('.tab', { hasText: 'Mechanism' }).click();
  const before = await atomCount(page);
  await page.getByRole('button', { name: /Apply arrows/ }).click();
  const after = await atomCount(page);
  expect(after).toBeGreaterThan(before);
  const product = await page.evaluate(() => {
    const app = (window as any).chemwrite;
    const doc = app.editor.doc;
    const arrows = [...doc.arrows.values()];
    const last = arrows[arrows.length - 1];
    const ids = [...doc.atoms.values()].filter((a: any) => a.x > Math.max(last.x1, last.x2)).map((a: any) => a.id);
    return ids.length;
  });
  expect(product).toBeGreaterThan(0);
  await expect(page.locator('.mech-warnings .warn')).toHaveCount(0);
});

test('save/load round trip via JSON', async ({ page }) => {
  await page.evaluate(() => (window as any).chemwrite.insertFromText('OC(=O)c1ccccc1O', null));
  const json = await page.evaluate(() => {
    const m = (window as any).chemwrite;
    return JSON.stringify((window as any).chemwrite.editor.doc.atoms.size) + '|' + m.currentSmiles();
  });
  expect(json).toContain('|');
  const smi = await smiles(page);
  await page.evaluate(() => {
    const app = (window as any).chemwrite;
    const data = app.editor.selectionToJSON();
    app.newDoc = app.newDoc.bind(app);
    app.editor.selectAll();
    app.editor.deleteSelected();
    app.editor.pasteDocJSON(data);
  });
  expect(await smiles(page)).toBe(smi);
});

test('export SVG downloads a file', async ({ page }) => {
  await page.evaluate(() => (window as any).chemwrite.insertFromText('c1ccccc1O', null));
  await page.keyboard.press('Control+e');
  await page.locator('.fmt', { hasText: 'SVG vector' }).click();
  const [download] = await Promise.all([page.waitForEvent('download'), page.getByRole('button', { name: 'Export', exact: true }).click()]);
  expect(download.suggestedFilename()).toMatch(/\.svg$/);
});

test('3D model generation reports a UFF energy', async ({ page }) => {
  await page.evaluate(() => (window as any).chemwrite.insertFromText('CCO', null));
  await page.locator('.tab', { hasText: '3D' }).click();
  await page.getByRole('button', { name: 'Generate 3D' }).click();
  await expect(page.locator('.viewer-info')).toContainText('kcal/mol', { timeout: 30000 });
});

test('screenshot of example drawing', async ({ page }) => {
  await page.evaluate(() => (window as any).chemwrite.loadExample());
  await page.waitForTimeout(400);
  await page.screenshot({ path: 'test-results/example-desktop.png' });
});
