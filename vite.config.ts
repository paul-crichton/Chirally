import { defineConfig, type Plugin } from 'vite';
import { readFileSync, writeFileSync, existsSync } from 'node:fs';
import { join } from 'node:path';

// GitHub Pages sites are public even when the repository is private, so the Pages
// build (PAGES=1) ships without source maps.
const pages = process.env.PAGES === '1';

/** Stamps public/sw.js in the build output with an id derived from the emitted (content-hashed) files. */
function serviceWorkerBuildId(): Plugin {
  let id = 'dev';
  return {
    name: 'sw-build-id',
    apply: 'build',
    generateBundle(_opts, bundle) {
      const names = Object.keys(bundle).sort().join('|');
      let h = 0x811c9dc5;
      for (let i = 0; i < names.length; i++) h = Math.imul(h ^ names.charCodeAt(i), 0x01000193) >>> 0;
      id = h.toString(36);
    },
    writeBundle(opts) {
      const file = join(opts.dir ?? 'dist', 'sw.js');
      if (existsSync(file)) writeFileSync(file, readFileSync(file, 'utf8').replace('__BUILD_ID__', id));
    },
  };
}

export default defineConfig({
  base: './',
  plugins: [serviceWorkerBuildId()],
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
