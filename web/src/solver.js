// LAMBDA in the browser: the same one-dimensional Chebyshev collocation as the
// C++ solver (src/solver/Solver.cpp), in plain JavaScript.
//
// Units: length R, time R^2/nu, pressure rho nu^2/R^2. lambda = -Gamma + i omega.
// Unknowns: Chebyshev coefficients of u(r), v(r), p(r) (Nr each) and the surface
// displacement d; u_r = u Y, u_tan = v r grad Y, p = p Y, surface r = 1 + d Y.

// ---------------------------------------------------------------- radial basis
export function points(Nr) {
  const r = new Float64Array(Nr);
  for (let j = 0; j < Nr; j++) {
    let s = j / (Nr - 1);
    if (j === 0) s = 1 / (2 * (Nr - 1));
    r[j] = Math.cos((Math.PI / 2) * (1 - s));
  }
  r[Nr - 1] = 1;
  return r;
}

// T_{2n+parity}(x) and its first two derivatives, n = 0..Nr-1
export function basis(Nr, parity, x) {
  const K = 2 * Nr;
  const t = new Float64Array(K), t1 = new Float64Array(K), t2 = new Float64Array(K);
  t[0] = 1; t[1] = x; t1[1] = 1;
  for (let k = 1; k + 1 < K; k++) {
    t[k + 1] = 2 * x * t[k] - t[k - 1];
    t1[k + 1] = 2 * t[k] + 2 * x * t1[k] - t1[k - 1];
    t2[k + 1] = 4 * t1[k] + 2 * x * t2[k] - t2[k - 1];
  }
  const T = new Float64Array(Nr), T1 = new Float64Array(Nr), T2 = new Float64Array(Nr);
  for (let n = 0; n < Nr; n++) {
    T[n] = t[2 * n + parity]; T1[n] = t1[2 * n + parity]; T2[n] = t2[2 * n + parity];
  }
  return { T, T1, T2 };
}

// ---------------------------------------------------------------- pencil A x = lambda B x
export function pencil(l, Oh, Nr) {
  const r = points(Nr);
  const L = l * (l + 1), We = Oh * Oh;
  const n = 3 * Nr + 1;
  const iu = 0, iv = Nr, ip = 2 * Nr, id = 3 * Nr;
  const A = new Float64Array(n * n), B = new Float64Array(n * n);
  const set = (M, i, j, v) => { M[i * n + j] = v; };
  for (let j = 0; j < Nr; j++) {
    const x = r[j], surface = j === Nr - 1;
    const V = basis(Nr, (l + 1) % 2, x), S = basis(Nr, l % 2, x);
    const rc = j, rr = Nr + j, rt = 2 * Nr + j;
    for (let k = 0; k < Nr; k++) {
      const u0 = V.T[k], u1 = V.T1[k], u2 = V.T2[k], p0 = S.T[k], p1 = S.T1[k];
      set(A, rc, iu + k, u1 + (2 * u0) / x);
      set(A, rc, iv + k, (-L * u0) / x);
      if (!surface) {
        set(A, rr, iu + k, (-L / (x * x)) * u0); set(B, rr, iu + k, u0);
        set(A, rr, iv + k, (L / x) * (u1 + u0 / x));
        set(A, rr, ip + k, -p1);
        set(A, rt, iu + k, -u1 / x);
        set(A, rt, iv + k, u2 + (2 * u1) / x); set(B, rt, iv + k, u0);
        set(A, rt, ip + k, -p0 / x);
      } else {
        set(A, rr, ip + k, p0);
        set(A, rr, iu + k, -2 * u1);
        set(A, rt, iv + k, u1 - u0);
        set(A, rt, iu + k, u0);
        set(A, id, iu + k, u0);
      }
    }
    if (surface) set(A, rr, id, (-(l - 1) * (l + 2)) / We);
  }
  set(B, id, id, 1);
  for (let i = 0; i < n; i++) {
    let s = 0;
    for (let j = 0; j < n; j++) s = Math.max(s, Math.abs(A[i * n + j]), Math.abs(B[i * n + j]));
    if (s > 0) for (let j = 0; j < n; j++) { A[i * n + j] /= s; B[i * n + j] /= s; }
  }
  return { A, B, n, Nr, l };
}

// ---------------------------------------------------------------- complex dense LU
function luFactor(P, sr, si) {
  // M = A - sigma B (complex), LU with partial pivoting, stored in place
  const { A, B, n } = P;
  const re = new Float64Array(n * n), im = new Float64Array(n * n);
  for (let k = 0; k < n * n; k++) { re[k] = A[k] - sr * B[k]; im[k] = -si * B[k]; }
  const piv = new Int32Array(n);
  for (let k = 0; k < n; k++) {
    let p = k, best = -1;
    for (let i = k; i < n; i++) {
      const a = Math.hypot(re[i * n + k], im[i * n + k]);
      if (a > best) { best = a; p = i; }
    }
    piv[k] = p;
    if (p !== k) {
      for (let j = 0; j < n; j++) {
        let t = re[k * n + j]; re[k * n + j] = re[p * n + j]; re[p * n + j] = t;
        t = im[k * n + j]; im[k * n + j] = im[p * n + j]; im[p * n + j] = t;
      }
    }
    let dr = re[k * n + k], di = im[k * n + k];
    if (dr === 0 && di === 0) { dr = 1e-300; re[k * n + k] = dr; }
    const den = dr * dr + di * di;
    for (let i = k + 1; i < n; i++) {
      const ar = re[i * n + k], ai = im[i * n + k];
      const mr = (ar * dr + ai * di) / den, mi = (ai * dr - ar * di) / den;
      re[i * n + k] = mr; im[i * n + k] = mi;
      if (mr === 0 && mi === 0) continue;
      for (let j = k + 1; j < n; j++) {
        const br = re[k * n + j], bi = im[k * n + j];
        re[i * n + j] -= mr * br - mi * bi;
        im[i * n + j] -= mr * bi + mi * br;
      }
    }
  }
  return { re, im, piv, n };
}

function luSolve(F, xr, xi) {
  const { re, im, piv, n } = F;
  for (let k = 0; k < n; k++) {
    const p = piv[k];
    if (p !== k) { let t = xr[k]; xr[k] = xr[p]; xr[p] = t; t = xi[k]; xi[k] = xi[p]; xi[p] = t; }
  }
  for (let i = 0; i < n; i++) {
    let sr = xr[i], si = xi[i];
    for (let j = 0; j < i; j++) {
      const ar = re[i * n + j], ai = im[i * n + j];
      sr -= ar * xr[j] - ai * xi[j]; si -= ar * xi[j] + ai * xr[j];
    }
    xr[i] = sr; xi[i] = si;
  }
  for (let i = n - 1; i >= 0; i--) {
    let sr = xr[i], si = xi[i];
    for (let j = i + 1; j < n; j++) {
      const ar = re[i * n + j], ai = im[i * n + j];
      sr -= ar * xr[j] - ai * xi[j]; si -= ar * xi[j] + ai * xr[j];
    }
    const dr = re[i * n + i], di = im[i * n + i], den = dr * dr + di * di;
    xr[i] = (sr * dr + si * di) / den; xi[i] = (si * dr - sr * di) / den;
  }
}

function matvec(M, n, xr, xi, yr, yi) {
  for (let i = 0; i < n; i++) {
    let sr = 0, si = 0;
    for (let j = 0; j < n; j++) { const m = M[i * n + j]; sr += m * xr[j]; si += m * xi[j]; }
    yr[i] = sr; yi[i] = si;
  }
}

// least-squares eigenvalue (Bx)^H A x / ||Bx||^2 and relative backward error
function rayleigh(P, xr, xi) {
  const { A, B, n } = P;
  const ar = new Float64Array(n), ai = new Float64Array(n), br = new Float64Array(n), bi = new Float64Array(n);
  matvec(A, n, xr, xi, ar, ai); matvec(B, n, xr, xi, br, bi);
  let nr = 0, ni = 0, d = 0, na = 0;
  for (let i = 0; i < n; i++) {
    nr += br[i] * ar[i] + bi[i] * ai[i]; ni += br[i] * ai[i] - bi[i] * ar[i];
    d += br[i] * br[i] + bi[i] * bi[i]; na += ar[i] * ar[i] + ai[i] * ai[i];
  }
  const lr = nr / d, li = ni / d;
  let res = 0;
  for (let i = 0; i < n; i++) {
    const rr = ar[i] - (lr * br[i] - li * bi[i]), ri = ai[i] - (lr * bi[i] + li * br[i]);
    res += rr * rr + ri * ri;
  }
  const eta = Math.sqrt(res) / (Math.sqrt(na) + Math.hypot(lr, li) * Math.sqrt(d));
  return { re: lr, im: li, eta };
}

function normalize(xr, xi) {
  let s = 0;
  for (let i = 0; i < xr.length; i++) s += xr[i] * xr[i] + xi[i] * xi[i];
  s = Math.sqrt(s);
  for (let i = 0; i < xr.length; i++) { xr[i] /= s; xi[i] /= s; }
}

// Eigenpair nearest to the shift (sr, si): inverse iteration with the shift fixed
// (converges to the nearest eigenvalue), then with the shift updated by the
// Rayleigh quotient (fast final convergence).
export function nearestEigenpair(P, sr, si, opts = {}) {
  const { n, B } = P;
  const maxFixed = opts.maxFixed ?? 60, refine = opts.refine ?? true;
  let xr = new Float64Array(n), xi = new Float64Array(n);
  for (let i = 0; i < n; i++) { xr[i] = 1 + 0.5 * Math.sin(1.3 * i + 0.7); xi[i] = 0.3 * Math.cos(0.9 * i); }
  normalize(xr, xi);
  const tr = new Float64Array(n), ti = new Float64Array(n);
  let F = luFactor(P, sr * (1 + 1e-13), si * (1 + 1e-13));
  let lam = { re: sr, im: si, eta: Infinity }, prev = null, it = 0;
  for (; it < maxFixed; it++) {
    matvec(B, n, xr, xi, tr, ti); luSolve(F, tr, ti);
    xr.set(tr); xi.set(ti); normalize(xr, xi);
    lam = rayleigh(P, xr, xi);
    if (prev && Math.hypot(lam.re - prev.re, lam.im - prev.im) < 1e-9 * Math.hypot(lam.re, lam.im)) break;
    prev = lam;
  }
  let best = { lam, xr: Float64Array.from(xr), xi: Float64Array.from(xi) };
  if (refine) {
    for (let k = 0; k < 8; k++) {
      F = luFactor(P, lam.re * (1 + 1e-14), lam.im * (1 + 1e-14));
      matvec(B, n, xr, xi, tr, ti); luSolve(F, tr, ti);
      xr.set(tr); xi.set(ti); normalize(xr, xi);
      const l2 = rayleigh(P, xr, xi);
      const stalled = l2.eta >= 0.5 * best.lam.eta && k >= 1;
      if (l2.eta < best.lam.eta) best = { lam: l2, xr: Float64Array.from(xr), xi: Float64Array.from(xi) };
      lam = l2;
      if (best.lam.eta < 1e-15 || stalled) break;
    }
  }
  return { ...best, iterations: it };
}

export function lambGuess(l, Oh) {
  return { re: -(l - 1) * (2 * l + 1), im: Math.sqrt(l * (l - 1) * (l + 2)) / Oh };
}

const isReal = (z) => Math.abs(z.im) <= 1e-10 * Math.hypot(z.re, z.im);

// The shape-oscillation eigenmode of degree l (same root selection as the C++ solver):
// Lamb's value as shift for Oh <= 0.02; otherwise the root is followed in log(Oh) from
// Oh = 0.01. After the pair has become real (aperiodic), the slower real root is followed.
export function solveMode(l, Oh, Nr = 48, guess = null) {
  const kFrom = 0.01;
  let shift;
  if (guess) {
    shift = guess;
  } else if (Oh <= 2 * kFrom) {
    shift = lambGuess(l, Oh);
  } else {
    const steps = Math.ceil(20 * Math.log10(Oh / kFrom));
    let o = kFrom;
    let z = nearestEigenpair(pencil(l, o, Nr), lambGuess(l, o).re, lambGuess(l, o).im, { refine: false }).lam;
    let zprev = null;
    for (let s = 1; s <= steps; s++) {
      const on = kFrom * Math.pow(Oh / kFrom, s / steps);
      let pr;
      if (isReal(z) || Math.abs(z.im) < 1e-3 * Math.abs(z.re)) {
        // aperiodic: extrapolate along the real axis, biased to the slower root
        const dz = zprev ? z.re - zprev.re : 0;
        pr = { re: z.re + dz + 0.02 * Math.abs(z.re), im: 0 };
      } else {
        pr = { re: z.re, im: (z.im * o) / on };
      }
      const wasComplex = !(isReal(z) || Math.abs(z.im) < 1e-3 * Math.abs(z.re));
      const Pn = pencil(l, on, Nr);
      zprev = z;
      z = nearestEigenpair(Pn, pr.re, pr.im, { refine: false }).lam;
      if (isReal(z)) {
        z = { re: z.re, im: 0, eta: z.eta };
        if (wasComplex) {
          // the pair has just split into two real roots placed about the old real part:
          // find the partner and keep the slower (less damped) one
          const c = zprev.re;
          const w = nearestEigenpair(Pn, 2 * c - z.re, 0, { refine: false }).lam;
          if (isReal(w) && w.re > z.re) z = { re: w.re, im: 0, eta: w.eta };
        }
      }
      o = on;
    }
    shift = z;
  }
  const P = pencil(l, Oh, Nr);
  let r = nearestEigenpair(P, shift.re, shift.im);
  if (!guess && isReal(r.lam)) {
    // aperiodic: return the least-damped real root (as the C++ solver); scan shifts
    // between this root and zero for real roots further to the right
    for (let k = 1; k <= 12; k++) {
      const sr = r.lam.re * (1 - k / 13);
      const w = nearestEigenpair(P, sr, 0, { refine: false });
      if (isReal(w.lam) && w.lam.re > r.lam.re * (1 - 1e-9) + 1e-12 && Math.abs(w.lam.re - r.lam.re) > 1e-8 * Math.abs(r.lam.re))
        r = nearestEigenpair(P, w.lam.re, 0);
    }
  }
  let { lam, xr, xi } = r;
  if (lam.im < 0) { lam = { ...lam, im: -lam.im }; xi = xi.map((v) => -v); }
  if (isReal(lam)) lam = { ...lam, im: 0 };
  // normalise the eigenvector to d = 1
  const n = P.n, dr = xr[n - 1], di = xi[n - 1], dd = dr * dr + di * di;
  const cr = new Float64Array(n), ci = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    cr[i] = (xr[i] * dr + xi[i] * di) / dd; ci[i] = (xi[i] * dr - xr[i] * di) / dd;
  }
  return { l, Oh, Nr, lambda: lam, backwardError: lam.eta, coeffs: { re: cr, im: ci } };
}

// ---------------------------------------------------------------- eigenfunctions
// radial amplitudes u, v, p, and g = ((r v)' - u)/r (vorticity factor), complex
export function radial(mode, r) {
  const { l, Nr, coeffs } = mode;
  const V = basis(Nr, (l + 1) % 2, r), S = basis(Nr, l % 2, r);
  let ur = 0, ui = 0, vr = 0, vi = 0, v1r = 0, v1i = 0, pr = 0, pi = 0;
  for (let k = 0; k < Nr; k++) {
    ur += V.T[k] * coeffs.re[k]; ui += V.T[k] * coeffs.im[k];
    vr += V.T[k] * coeffs.re[Nr + k]; vi += V.T[k] * coeffs.im[Nr + k];
    v1r += V.T1[k] * coeffs.re[Nr + k]; v1i += V.T1[k] * coeffs.im[Nr + k];
    pr += S.T[k] * coeffs.re[2 * Nr + k]; pi += S.T[k] * coeffs.im[2 * Nr + k];
  }
  const rr = Math.max(r, 1e-12);
  // (r v)' - u = v + r v' - u
  const gr = (vr + r * v1r - ur) / rr, gi = (vi + r * v1i - ui) / rr;
  return { u: [ur, ui], v: [vr, vi], p: [pr, pi], g: [gr, gi] };
}

// associated Legendre P_l^m(cos theta) (Condon-Shortley), dP/dtheta, P/sin(theta)
export function legendre(l, m, theta) {
  const mu = Math.cos(theta), s = Math.sin(theta);
  const K = m + 1;
  const D = [];
  for (let k = 0; k <= K; k++) {
    const row = new Float64Array(l + 2);
    let pm1 = 0, p0 = k === 0 ? 1 : 0;
    row[0] = p0;
    for (let n = 0; n < l + 1; n++) {
      const dk1 = k > 0 ? D[k - 1][n] : 0;
      const p1 = ((2 * n + 1) * (mu * p0 + k * dk1) - n * pm1) / (n + 1);
      row[n + 1] = p1; pm1 = p0; p0 = p1;
    }
    D.push(row);
  }
  const sign = m % 2 ? -1 : 1;
  if (m === 0) return { P: D[0][l], dP: -s * D[1][l], Ps: 0 };
  const sm1 = Math.pow(s, m - 1);
  return {
    P: sign * sm1 * s * D[m][l],
    Ps: sign * sm1 * D[m][l],
    dP: sign * (m * sm1 * mu * D[m][l] - sm1 * s * s * D[m + 1][l]),
  };
}

export function legendreNorm(l, m) {
  let pmax = 0;
  for (let i = 0; i <= 2000; i++) pmax = Math.max(pmax, Math.abs(legendre(l, m, (Math.PI * i) / 2000).P));
  return 1 / pmax;
}

// Lamb and Rayleigh values in viscous units
export function lambRayleigh(l, Oh) {
  return { gamma: (l - 1) * (2 * l + 1), omega: Math.sqrt(l * (l - 1) * (l + 2)) / Oh };
}
