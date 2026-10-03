import './styles.css';
import { App } from './app/app';

const root = document.getElementById('app')!;
const app = new App(root);
// expose for debugging and end-to-end tests
(window as unknown as { chemwrite: App }).chemwrite = app;

// offline support in production builds
if (import.meta.env.PROD && 'serviceWorker' in navigator) {
  window.addEventListener('load', () => navigator.serviceWorker.register('./sw.js').catch(() => undefined));
}
