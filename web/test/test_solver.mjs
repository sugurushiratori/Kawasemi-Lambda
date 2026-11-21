// Compares the browser solver (src/solver.js) with the exact roots used for the
// C++ tests (../tests/reference/chandrasekhar_roots.csv): relative errors of
// lambda, Gamma and omega below 1e-12 at Nr = 48, as for the C++ solver.
// Also checks the eigenfunction (kinematic condition, div u = 0) and Legendre functions.
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { solveMode, radial, legendre } from '../src/solver.js';

const here = fileURLToPath(new URL('.', import.meta.url));
const csv = readFileSync(here + '../../tests/reference/chandrasekhar_roots.csv', 'utf8');
let n = 0, fail = 0, tmax = 0;
for (const line of csv.split('\n')) {
  if (!line || line[0] === '#' || line[0] === 'W') continue;
  const [w, l, re, im] = line.split(',');
  const exRe = Number(re), exIm = Number(im);
  const t0 = performance.now();
  const m = solveMode(Number(l), Math.sqrt(Number(w)), 48);
  const dt = performance.now() - t0; tmax = Math.max(tmax, dt);
  const z = m.lambda;
  const eL = Math.hypot(z.re - exRe, z.im - exIm) / Math.hypot(exRe, exIm);
  const eG = Math.abs(z.re / exRe - 1);
  const eW = exIm !== 0 ? Math.abs(z.im / exIm - 1) : Math.abs(z.im) / Math.hypot(exRe, exIm);
  const ok = eL < 1e-12 && eG < 1e-12 && eW < 1e-12;
  console.log(`We ${w.padEnd(10)} l ${l}  lambda ${z.re.toExponential(15)} ${z.im.toExponential(15)}  dl ${eL.toExponential(1)} dG ${eG.toExponential(1)} dw ${eW.toExponential(1)}  bwd ${m.backwardError.toExponential(1)}  ${dt.toFixed(0)} ms  ${ok ? 'ok' : 'FAIL'}`);
  n++; if (!ok) fail++;

  // kinematic condition u(1) = lambda d (d = 1)
  const a = radial(m, 1.0);
  const ek = Math.hypot(a.u[0] - z.re, a.u[1] - z.im) / Math.hypot(z.re, z.im);
  if (!(ek < 1e-11)) { console.log('  kinematic FAIL', ek); fail++; }
}
// div u = 0 for the poloidal field: u' + 2u/r - L v/r = 0 at an interior point
{
  const m = solveMode(3, 0.0316, 48), r = 0.6, h = 1e-5, L = 12;
  const a = radial(m, r), ap = radial(m, r + h), am = radial(m, r - h);
  const du = [(ap.u[0] - am.u[0]) / (2 * h), (ap.u[1] - am.u[1]) / (2 * h)];
  const div = Math.hypot(du[0] + 2 * a.u[0] / r - L * a.v[0] / r, du[1] + 2 * a.u[1] / r - L * a.v[1] / r);
  const e = div / Math.hypot(m.lambda.re, m.lambda.im);
  console.log(`continuity residual ${e.toExponential(1)} ${e < 1e-6 ? 'ok' : 'FAIL'}`); if (!(e < 1e-6)) fail++;
}
// Legendre closed forms
{
  let e = 0;
  for (const th of [0, 0.3, 1.1, 2.0, Math.PI]) {
    const mu = Math.cos(th), s = Math.sin(th);
    e = Math.max(e, Math.abs(legendre(2, 0, th).P - (3 * mu * mu - 1) / 2));
    e = Math.max(e, Math.abs(legendre(2, 1, th).P + 3 * mu * s));
    e = Math.max(e, Math.abs(legendre(4, 3, th).P + 105 * mu * s * s * s));
  }
  console.log(`Legendre ${e.toExponential(1)} ${e < 1e-12 ? 'ok' : 'FAIL'}`); if (!(e < 1e-12)) fail++;
}
console.log(`${n - fail} of ${n} roots passed; slowest ${tmax.toFixed(0)} ms`);
process.exit(fail ? 1 : 0);
