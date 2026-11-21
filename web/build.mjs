// Bundles src/app.js (with three.js) into one self-contained HTML file:
//   ../docs/index.html   (GitHub Pages; also opens directly from disk)
import { build } from 'esbuild';
import { readFileSync, writeFileSync, mkdirSync } from 'node:fs';

const res = await build({
  entryPoints: ['src/app.js'], bundle: true, minify: true, format: 'iife',
  target: ['es2020'], write: false, legalComments: 'eof',
});
const js = res.outputFiles[0].text.replace(/<\/script/gi, '<\\/script');
const b64 = (f) => readFileSync(f).toString('base64');
const html = readFileSync('src/index.html', 'utf8')
  .replace('/*FAVICON*/', () => 'data:image/png;base64,' + b64('../docs/logo/lambda-favicon-64.png'))
  .replace('/*LOGO*/', () => 'data:image/png;base64,' + b64('../docs/logo/lambda-logo-light-300w@2x.png'))
  .replace('/*BUNDLE*/', () => js);
const banner = '<!-- LAMBDA web (built by web/build.mjs from web/src; three.js, MIT licence) -->\n';
mkdirSync('../docs', { recursive: true });
writeFileSync('../docs/index.html', banner + html);
console.log('wrote ../docs/index.html', (html.length / 1024).toFixed(0), 'kB');
