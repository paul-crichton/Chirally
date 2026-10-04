// Chirally service worker: offline support for the app shell (same-origin GET requests only).
// The cache name and the precache list are stamped at build time (see vite.config.ts).
const PREFIX = 'chirally-';
const CACHE = PREFIX + '__BUILD_ID__';
const PRECACHE = /*__PRECACHE__*/ [];
const SHELL = new URL('./', self.location).href; // the app's index.html, e.g. https://<user>.github.io/Chirally/
const ASSETS = new URL('./assets/', self.location).pathname; // content-hashed build output

self.addEventListener('install', (e) => {
  // Precache this build so the app works offline right after the first visit and right after an update.
  e.waitUntil(
    caches
      .open(CACHE)
      .then((c) => c.addAll(PRECACHE.map((u) => new Request(u, { cache: 'reload' }))))
      .then(() => self.skipWaiting()),
  );
});
self.addEventListener('activate', (e) => {
  // Only remove Chirally's own old caches: Cache Storage is shared by every site on the origin.
  e.waitUntil(
    caches
      .keys()
      .then((keys) => Promise.all(keys.filter((k) => k.startsWith(PREFIX) && k !== CACHE).map((k) => caches.delete(k))))
      .then(() => self.clients.claim()),
  );
});

const put = (e, key, res) => {
  if (res.ok && res.type === 'basic') {
    const copy = res.clone();
    e.waitUntil(caches.open(CACHE).then((c) => c.put(key, copy)));
  }
  return res;
};
const fromCache = (key) => caches.open(CACHE).then((c) => c.match(key));

self.addEventListener('fetch', (e) => {
  const req = e.request;
  const url = new URL(req.url);
  if (req.method !== 'GET' || url.origin !== self.location.origin) return; // PubChem etc. always go to the network
  if (req.mode === 'navigate') {
    // network first so new versions arrive; one cache entry for the shell (no ?query / #doc= keys)
    e.respondWith(
      // no-cache: revalidate (cheap 304) instead of trusting GitHub Pages' max-age=600, which can hand back
      // an old index.html whose hashed assets were deleted by the deploy and by the new cache
      fetch(req, { cache: 'no-cache' })
        .then((res) => put(e, SHELL, res))
        .catch(() => fromCache(SHELL).then((r) => r || Response.error())),
    );
    return;
  }
  if (url.pathname.startsWith(ASSETS)) {
    // content-hashed build assets never change: cache first
    e.respondWith(fromCache(req).then((hit) => hit || fetch(req).then((res) => put(e, req, res))));
    return;
  }
  // everything else (manifest, icon): network first, cache when offline
  e.respondWith(
    fetch(req)
      .then((res) => put(e, req, res))
      .catch(() => fromCache(req).then((r) => r || Response.error())),
  );
});
