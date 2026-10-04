import { test, expect } from '@playwright/test';
import { openApp, atomCount } from './helpers';

test('mobile layout: bottom toolbar, tap to draw, panel sheet', async ({ page }) => {
  await openApp(page);
  const toolbar = page.locator('.toolbar');
  const tb = (await toolbar.boundingBox())!;
  const vp = page.viewportSize()!;
  expect(tb.y + tb.height).toBeGreaterThan(vp.height - 120);
  const box = (await page.locator('canvas.chirally-canvas').boundingBox())!;
  await page.touchscreen.tap(box.x + box.width / 2, box.y + box.height / 2);
  expect(await atomCount(page)).toBe(2);
  await page.locator('.panel-toggle').tap();
  await expect(page.locator('.sidebar')).toBeInViewport();
  await page.screenshot({ path: 'test-results/mobile.png' });
});
