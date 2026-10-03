import { defineConfig } from 'vite';

// GitHub Pages sites are public even when the repository is private, so the Pages
// build (PAGES=1) ships without source maps.
const pages = process.env.PAGES === '1';

export default defineConfig({
  base: './',
  build: {
    target: 'es2022',
    sourcemap: !pages,
    chunkSizeWarningLimit: 900,
  },
  test: {
    globals: true,
    environment: 'node',
    include: ['tests/**/*.test.ts', 'src/**/*.test.ts'],
  },
});
