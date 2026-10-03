import './styles.css';
import { App } from './app/app';

const root = document.getElementById('app')!;
const app = new App(root);
// expose for debugging and end-to-end tests
(window as unknown as { chemwrite: App }).chemwrite = app;
