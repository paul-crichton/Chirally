// ChemWrite service worker: offline support for the app shell (same-origin GET requests only).
// The cache name is stamped with a build id at build time (see vite.config.ts), so each deployment
// gets a fresh cache and the previous one is removed on activation.
const PREFIX = 'chemwrite-';
const CACHE = PREFIX + '__BUILD_ID__';

self.addEventListener('install', () => self.skipWaiting());
self.addEventListener('activate', (e) => {
  // Only remove ChemWrite's own old caches: Cache Storage is shared by every site on the origin
  // (e.g. all GitHub Pages project sites of one user live on <user>.github.io).
  e.waitUntil(
    caches
      .keys()
      .then((keys) => Promise.all(keys.filter((k) => k.startsWith(PREFIX) && k !== CACHE).map((k) => caches.delete(k))))
      .then(() => self.clients.claim()),
  );
});

self.addEventListener('fetch', (e) => {
  const req = e.request;
  const url = new URL(req.url);
  if (req.method !== 'GET' || url.origin !== self.location.origin) return; // PubChem etc. always go to the network
  if (req.mode === 'navigate') {
    // network first so new versions arrive, cache fallback when offline
    e.respondWith(
      fetch(req)
        .then((res) => {
          if (res.ok && res.type === 'basic') {
            const copy = res.clone();
            caches.open(CACHE).then((c) => c.put(req, copy));
          }
          return res;
        })
        .catch(() => caches.match(req).then((r) => r || caches.match('./'))),
    );
    return;
  }
  // hashed build assets: cache first
  e.respondWith(
    caches.match(req).then(
      (hit) =>
        hit ||
        fetch(req).then((res) => {
          if (res.ok && res.type === 'basic') {
            const copy = res.clone();
            caches.open(CACHE).then((c) => c.put(req, copy));
          }
          return res;
        }),
    ),
  );
});
