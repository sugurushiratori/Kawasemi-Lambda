// LAMBDA web: interactive eigenmodes of a viscous drop.
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { solveMode, radial, legendre, legendreNorm } from './solver.js';
import { MATERIALS } from './materials.js';
import { TEXT } from './i18n.js';

// ------------------------------------------------------------------ state
const state = {
  lang: (navigator.language || '').startsWith('ja') ? 'ja' : 'en',
  material: 'Ni',
  rho: 7890, mu: 4.7, sigma: 1.739, R: 1.0,
  useOh: false, Oh: 0.01,
  l: 2, m: 0, Nr: 48,
  amp: 0.12, playing: true, speed: 0.5, damping: false,
  surfColor: 'disp', opacity: 0.55, grid: true,
  vectors: 'section', vscale: 3.0,
  section: 'speed',
};
let mode = null;          // result of solveMode
let fields = null;        // precomputed amplitudes for drawing
let ynorm = 1;
const t = () => TEXT[state.lang];

// ------------------------------------------------------------------ three.js scene
const canvas = document.getElementById('canvas');
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
const scene = new THREE.Scene();
scene.background = new THREE.Color().setRGB(0.32, 0.34, 0.43, THREE.SRGBColorSpace);   // ParaView's default background
const camera = new THREE.PerspectiveCamera(35, 1, 0.01, 100);
camera.position.set(3.2, 1.6, 4.2);
const controls = new OrbitControls(camera, canvas);
controls.enableDamping = true;
// ParaView-like lighting: a headlight that follows the camera, plus some ambient light
scene.add(new THREE.AmbientLight(0xffffff, 0.9));
const headlight = new THREE.DirectionalLight(0xffffff, 2.2);
scene.add(headlight); scene.add(headlight.target);

// orientation axes in the lower-left corner, as in ParaView (X red, Y yellow, Z green; Z = polar axis)
const axScene = new THREE.Scene();
const axCam = new THREE.PerspectiveCamera(30, 1, 0.1, 20);
axScene.add(new THREE.AmbientLight(0xffffff, 1.2));
const axLight = new THREE.DirectionalLight(0xffffff, 1.5); axScene.add(axLight); axScene.add(axLight.target);
function axisArrow(dir, color, label) {
  const g = new THREE.Group();
  const mat = new THREE.MeshPhongMaterial({ color });
  const shaft = new THREE.Mesh(new THREE.CylinderGeometry(0.035, 0.035, 0.75, 12), mat); shaft.position.y = 0.375;
  const head = new THREE.Mesh(new THREE.ConeGeometry(0.09, 0.25, 16), mat); head.position.y = 0.875;
  g.add(shaft, head);
  g.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir);
  const cv = document.createElement('canvas'); cv.width = cv.height = 64;
  const cx = cv.getContext('2d'); cx.fillStyle = '#fff'; cx.font = 'bold 44px sans-serif'; cx.textAlign = 'center'; cx.textBaseline = 'middle'; cx.fillText(label, 32, 34);
  const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: new THREE.CanvasTexture(cv), depthTest: false }));
  sp.scale.set(0.3, 0.3, 1); sp.position.copy(dir.clone().multiplyScalar(1.18));
  axScene.add(g, sp);
}
// three.js y is the polar axis (Z); x stays X; three.js z is -Y for a right-handed (X, Y, Z)
axisArrow(new THREE.Vector3(1, 0, 0), 0xe63b2e, 'X');
axisArrow(new THREE.Vector3(0, 0, -1), 0xf2d22e, 'Y');
axisArrow(new THREE.Vector3(0, 1, 0), 0x3cc24a, 'Z');

// polar axis: +y. point(r, theta, phi) = (r sin th cos ph, r cos th, r sin th sin ph)
const P3 = (r, th, ph) => [r * Math.sin(th) * Math.cos(ph), r * Math.cos(th), r * Math.sin(th) * Math.sin(ph)];
function toCart(ur, ut, uf, th, ph) {   // complex spherical components -> complex Cartesian
  const st = Math.sin(th), ct = Math.cos(th), sp = Math.sin(ph), cp = Math.cos(ph);
  const er = [st * cp, ct, st * sp], et = [ct * cp, -st, ct * sp], ef = [-sp, 0, cp];
  const out = [];
  for (let k = 0; k < 3; k++)
    out.push([ur[0] * er[k] + ut[0] * et[k] + uf[0] * ef[k], ur[1] * er[k] + ut[1] * et[k] + uf[1] * ef[k]]);
  return out;
}
const cmul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
const cabs = (a) => Math.hypot(a[0], a[1]);

// colormaps
// ParaView's default 'Cool to Warm' colour map (Moreland), sampled
const CTW = [[0.2298, 0.2987, 0.7537], [0.3537, 0.4729, 0.8674], [0.4845, 0.6283, 0.9536], [0.6201, 0.7509, 0.9935],
  [0.7468, 0.8318, 0.9590], [0.8392, 0.8508, 0.8803], [0.9020, 0.8064, 0.7505], [0.9576, 0.7266, 0.6170],
  [0.9562, 0.5976, 0.4767], [0.9006, 0.4446, 0.3418], [0.7987, 0.2763, 0.2219], [0.7057, 0.0156, 0.1502]];
function coolwarm(x) { // x in [0, 1]
  x = Math.max(0, Math.min(1, x)) * (CTW.length - 1);
  const i = Math.min(CTW.length - 2, Math.floor(x)), f = x - i;
  return CTW[i].map((c, k) => c + f * (CTW[i + 1][k] - c));
}
// vertex colours are interpreted as linear: convert the (sRGB) colour map
const toLinear = (c) => (c <= 0.04045 ? c / 12.92 : Math.pow((c + 0.055) / 1.055, 2.4));
const LUT = Array.from({ length: 256 }, (_, i) => coolwarm(i / 255).map(toLinear));
const lut = (x) => LUT[Math.max(0, Math.min(255, Math.round(255 * x)))];
const diverging = (x) => lut(0.5 + 0.5 * x);   // x in [-1, 1]
const sequential = (x) => lut(x);              // x in [0, 1]
const CTW_CSS = 'linear-gradient(to top,' + CTW.map((c, k) => `rgb(${c.map((v) => Math.round(255 * v)).join(',')}) ${(100 * k / (CTW.length - 1)).toFixed(1)}%`).join(',') + ')';

// surface mesh
const surfGeom = new THREE.SphereGeometry(1, 160, 120);
const nSurf = surfGeom.attributes.position.count;
const surfBase = Float32Array.from(surfGeom.attributes.position.array);
surfGeom.setAttribute('color', new THREE.BufferAttribute(new Float32Array(nSurf * 3), 3));
const surfMat = new THREE.MeshPhongMaterial({ vertexColors: true, shininess: 30, specular: 0x333333, transparent: true, side: THREE.DoubleSide, depthWrite: false });
const surface = new THREE.Mesh(surfGeom, surfMat);
scene.add(surface);
const surfAngles = new Float32Array(nSurf * 2);
for (let i = 0; i < nSurf; i++) {
  const x = surfBase[3 * i], y = surfBase[3 * i + 1], z = surfBase[3 * i + 2];
  surfAngles[2 * i] = Math.acos(Math.max(-1, Math.min(1, y)));
  surfAngles[2 * i + 1] = Math.atan2(z, x);
}

// grid of latitude/longitude lines on the deformed surface, with hidden lines removed
const G_NT = 72, G_NP = 144, G_STEP = 6;              // 2.5 deg sampling, lines every 15 deg
const gridNodes = (G_NT + 1) * (G_NP + 1);
const gridAng = [];
for (let i = 0; i <= G_NT; i++) for (let j = 0; j <= G_NP; j++) gridAng.push([(Math.PI * i) / G_NT, (2 * Math.PI * j) / G_NP]);
const gridSegs = [];
for (let i = G_STEP; i < G_NT; i += G_STEP) for (let j = 0; j < G_NP; j++) gridSegs.push([i * (G_NP + 1) + j, i * (G_NP + 1) + j + 1]);
for (let j = 0; j < G_NP; j += G_STEP) for (let i = 0; i < G_NT; i++) gridSegs.push([i * (G_NP + 1) + j, (i + 1) * (G_NP + 1) + j]);
const gridGeom = new THREE.BufferGeometry();
gridGeom.setAttribute('position', new THREE.BufferAttribute(new Float32Array(gridSegs.length * 6), 3));
const gridLines = new THREE.LineSegments(gridGeom, new THREE.LineBasicMaterial({ color: 0x000000, transparent: true, opacity: 0.5 }));
gridLines.renderOrder = 2;
scene.add(gridLines);
const gridPos = new Float32Array(gridNodes * 3), gridNor = new Float32Array(gridNodes * 3);

// section disc in the plane z = 0 (phi = 0 for x > 0, phi = pi for x < 0)
const SEC_NR = 56, SEC_NA = 241;
const secR = Array.from({ length: SEC_NR }, (_, i) => Math.sin((Math.PI / 2) * (i / (SEC_NR - 1))));
const secGeom = new THREE.BufferGeometry();
{
  const pos = new Float32Array(SEC_NR * SEC_NA * 3), idx = [];
  for (let j = 0; j < SEC_NA; j++) for (let i = 0; i < SEC_NR; i++) {
    const a = (2 * Math.PI * j) / (SEC_NA - 1), k = j * SEC_NR + i;
    pos[3 * k] = secR[i] * Math.sin(a); pos[3 * k + 1] = secR[i] * Math.cos(a);
    if (i < SEC_NR - 1 && j < SEC_NA - 1) {
      const p = k, q = k + 1, r = k + SEC_NR, s = k + SEC_NR + 1;
      idx.push(p, r, q, q, r, s);
    }
  }
  secGeom.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  secGeom.setAttribute('color', new THREE.BufferAttribute(new Float32Array(pos.length), 3));
  secGeom.setIndex(idx);
}
const secMesh = new THREE.Mesh(secGeom, new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide }));
scene.add(secMesh);

// arrows (instanced)
const MAX_ARROWS = 700;
function arrowGeometry() {
  const shaft = new THREE.CylinderGeometry(0.018, 0.018, 0.72, 8); shaft.translate(0, 0.36, 0);
  const head = new THREE.ConeGeometry(0.05, 0.28, 12); head.translate(0, 0.86, 0);
  const g = new THREE.BufferGeometry();
  const merge = (geoms) => {
    let n = 0; for (const x of geoms) n += x.attributes.position.count;
    const pos = new Float32Array(n * 3), nor = new Float32Array(n * 3); const idx = []; let off = 0;
    for (const x of geoms) {
      pos.set(x.attributes.position.array, off * 3); nor.set(x.attributes.normal.array, off * 3);
      for (const i of x.index.array) idx.push(i + off); off += x.attributes.position.count;
    }
    g.setAttribute('position', new THREE.BufferAttribute(pos, 3)); g.setAttribute('normal', new THREE.BufferAttribute(nor, 3)); g.setIndex(idx);
  };
  merge([shaft, head]); return g;
}
const arrows = new THREE.InstancedMesh(arrowGeometry(), new THREE.MeshPhongMaterial({ color: 0xf2f2f2, shininess: 20 }), MAX_ARROWS);
arrows.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
scene.add(arrows);

// polar axis line
{
  const g = new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(0, -1.45, 0), new THREE.Vector3(0, 1.45, 0)]);
  const axis = new THREE.Line(g, new THREE.LineDashedMaterial({ color: 0xc8ccd6, dashSize: 0.05, gapSize: 0.04 }));
  axis.computeLineDistances();
  scene.add(axis);
}

// ------------------------------------------------------------------ field precomputation
function angular(th, ph) { // Y, dY/dtheta, (1/sin) dY/dphi as complex numbers
  const L = legendre(state.l, state.m, th);
  const e = [ynorm * Math.cos(state.m * ph), ynorm * Math.sin(state.m * ph)];
  return { Y: [L.P * e[0], L.P * e[1]], dY: [L.dP * e[0], L.dP * e[1]], sY: cmul([0, state.m * L.Ps], e) };
}

function pointFields(r, th, ph, rad, ang) {
  const a = rad || radial(mode, r), A = ang || angular(th, ph);
  const ur = cmul(a.u, A.Y), ut = cmul(a.v, A.dY), uf = cmul(a.v, A.sY);
  const vel = toCart(ur, ut, uf, th, ph);
  const p = cmul(a.p, A.Y);
  const wt = cmul([-a.g[0], -a.g[1]], A.sY), wf = cmul(a.g, A.dY);
  const vort = toCart([0, 0], wt, wf, th, ph);
  return { vel, p, vort, Y: A.Y };
}

function fibonacci(n) {
  const out = [], ga = Math.PI * (3 - Math.sqrt(5));
  for (let i = 0; i < n; i++) {
    const y = 1 - (2 * (i + 0.5)) / n;
    out.push([Math.acos(y), (i * ga) % (2 * Math.PI)]);
  }
  return out;
}

function precompute() {
  ynorm = legendreNorm(state.l, state.m);
  // surface
  const rad1 = radial(mode, 1);
  const surf = [];
  let pmaxS = 0;
  for (let i = 0; i < nSurf; i++) {
    const th = surfAngles[2 * i], ph = surfAngles[2 * i + 1];
    const A = angular(th, ph), p = cmul(rad1.p, A.Y);
    surf.push({ Y: A.Y, p });
    pmaxS = Math.max(pmaxS, cabs(p));
  }
  // section
  const radCache = secR.map((r) => radial(mode, Math.max(r, 1e-6)));
  const sec = []; let vmax = 0, pmax = 0, wmax = 0;
  for (let j = 0; j < SEC_NA; j++) {
    const a = (2 * Math.PI * j) / (SEC_NA - 1);
    const th = a <= Math.PI ? a : 2 * Math.PI - a, ph = a <= Math.PI ? 0 : Math.PI;
    const A = angular(th, ph);
    for (let i = 0; i < SEC_NR; i++) {
      const f = pointFields(secR[i], th, ph, radCache[i], A);
      f.th = th; f.ph = ph; f.r = secR[i];
      sec.push(f);
      vmax = Math.max(vmax, Math.hypot(...f.vel.map(cabs)));
      pmax = Math.max(pmax, cabs(f.p)); wmax = Math.max(wmax, Math.hypot(...f.vort.map(cabs)));
    }
  }
  // arrow sample points
  const sets = { section: [], volume: [], surface: [] };
  for (const r of [0.22, 0.38, 0.52, 0.65, 0.77, 0.87, 0.95]) {
    const n = Math.max(8, Math.round(40 * r));
    for (let k = 0; k < n; k++) {
      const a = (2 * Math.PI * (k + 0.5)) / n;
      const th = a <= Math.PI ? a : 2 * Math.PI - a, ph = a <= Math.PI ? 0 : Math.PI;
      sets.section.push({ r, th, ph });
    }
  }
  for (const [r, n] of [[0.35, 40], [0.62, 110], [0.88, 220]]) for (const [th, ph] of fibonacci(n)) sets.volume.push({ r, th, ph });
  for (const [th, ph] of fibonacci(320)) sets.surface.push({ r: 1, th, ph });
  for (const key of Object.keys(sets)) sets[key] = sets[key].slice(0, MAX_ARROWS).map((s) => ({ ...s, ...pointFields(s.r, s.th, s.ph), pos: P3(1, s.th, s.ph) }));
  const gridY = gridAng.map(([th, ph]) => angular(th, ph).Y);
  fields = { surf, pmaxS, sec, vmax, pmax, wmax, sets, gridY, lamAbs: Math.hypot(mode.lambda.re, mode.lambda.im) };
}

// ------------------------------------------------------------------ animation
let phase = 0, decayArg = 0, last = performance.now();
const re = (q, c, s) => q[0] * c - q[1] * s;  // Re[q e^{i tau}]
const tmpM = new THREE.Matrix4(), tmpQ = new THREE.Quaternion(), up = new THREE.Vector3(0, 1, 0), dirV = new THREE.Vector3();
const tmpP = new THREE.Vector3(), tmpS = new THREE.Vector3();

function frame(now) {
  const dt = Math.min(0.1, (now - last) / 1000); last = now;
  controls.update();
  if (mode && fields) {
    const aperiodic = mode.lambda.im === 0;
    const zeta = aperiodic ? 0 : -mode.lambda.re / mode.lambda.im;   // Gamma / omega
    if (state.playing) {
      const dphi = 2 * Math.PI * state.speed * dt;
      if (aperiodic) decayArg += dphi / (2 * Math.PI) * 1.5;
      else { phase += dphi; if (state.damping) decayArg += zeta * dphi; }
      if (Math.exp(-decayArg) < 0.03) { decayArg = 0; phase = 0; }
    }
    const fac = (state.damping || aperiodic ? Math.exp(-decayArg) : 1);
    const c = Math.cos(phase), s = Math.sin(phase);
    const eps = state.amp * fac;

    // surface
    const pos = surfGeom.attributes.position.array, col = surfGeom.attributes.color.array;
    for (let i = 0; i < nSurf; i++) {
      const f = fields.surf[i], disp = re(f.Y, c, s);
      const k = 1 + eps * disp;
      pos[3 * i] = surfBase[3 * i] * k; pos[3 * i + 1] = surfBase[3 * i + 1] * k; pos[3 * i + 2] = surfBase[3 * i + 2] * k;
      let rgb;
      if (state.surfColor === 'disp') rgb = diverging(disp * (fac > 0 ? 1 : 0));
      else if (state.surfColor === 'pres') rgb = diverging(re(f.p, c, s) / (fields.pmaxS || 1));
      else rgb = [1, 1, 1];
      col[3 * i] = rgb[0]; col[3 * i + 1] = rgb[1]; col[3 * i + 2] = rgb[2];
    }
    surfGeom.attributes.position.needsUpdate = true; surfGeom.attributes.color.needsUpdate = true;
    surfGeom.computeVertexNormals();
    surfMat.opacity = state.opacity; surfMat.depthWrite = state.opacity > 0.98; surface.visible = state.opacity > 0.01;

    // surface grid: deformed lattice, normals from neighbouring nodes, back-facing segments hidden
    gridLines.visible = state.grid;
    if (state.grid) {
      for (let k = 0; k < gridNodes; k++) {
        const [th, ph] = gridAng[k], rr = 1.003 * (1 + eps * re(fields.gridY[k], c, s));
        const p = P3(rr, th, ph);
        gridPos[3 * k] = p[0]; gridPos[3 * k + 1] = p[1]; gridPos[3 * k + 2] = p[2];
      }
      const W = G_NP + 1;
      for (let i = 0; i <= G_NT; i++) for (let j = 0; j <= G_NP; j++) {
        const k = i * W + j;
        if (i === 0 || i === G_NT) { gridNor[3 * k] = gridPos[3 * k]; gridNor[3 * k + 1] = gridPos[3 * k + 1]; gridNor[3 * k + 2] = gridPos[3 * k + 2]; continue; }
        const a = (i + 1) * W + j, b = (i - 1) * W + j;
        const jp = j === G_NP ? 1 : j + 1, jm = j === 0 ? G_NP - 1 : j - 1;
        const cc = i * W + jp, d = i * W + jm;
        const tx = gridPos[3 * a] - gridPos[3 * b], ty = gridPos[3 * a + 1] - gridPos[3 * b + 1], tz = gridPos[3 * a + 2] - gridPos[3 * b + 2];
        const ux = gridPos[3 * cc] - gridPos[3 * d], uy = gridPos[3 * cc + 1] - gridPos[3 * d + 1], uz = gridPos[3 * cc + 2] - gridPos[3 * d + 2];
        // outward normal = d/dphi x d/dtheta
        gridNor[3 * k] = uy * tz - uz * ty; gridNor[3 * k + 1] = uz * tx - ux * tz; gridNor[3 * k + 2] = ux * ty - uy * tx;
      }
      const gp = gridGeom.attributes.position.array, cam = camera.position;
      for (let q = 0; q < gridSegs.length; q++) {
        const [a, b] = gridSegs[q];
        const mx = 0.5 * (gridPos[3 * a] + gridPos[3 * b]), my = 0.5 * (gridPos[3 * a + 1] + gridPos[3 * b + 1]), mz = 0.5 * (gridPos[3 * a + 2] + gridPos[3 * b + 2]);
        const nx = gridNor[3 * a] + gridNor[3 * b], ny = gridNor[3 * a + 1] + gridNor[3 * b + 1], nz = gridNor[3 * a + 2] + gridNor[3 * b + 2];
        const visible = nx * (cam.x - mx) + ny * (cam.y - my) + nz * (cam.z - mz) > 0;
        const o = 6 * q;
        if (visible) {
          gp[o] = gridPos[3 * a]; gp[o + 1] = gridPos[3 * a + 1]; gp[o + 2] = gridPos[3 * a + 2];
          gp[o + 3] = gridPos[3 * b]; gp[o + 4] = gridPos[3 * b + 1]; gp[o + 5] = gridPos[3 * b + 2];
        } else {
          gp[o] = gp[o + 3] = mx; gp[o + 1] = gp[o + 4] = my; gp[o + 2] = gp[o + 5] = mz;
        }
      }
      gridGeom.attributes.position.needsUpdate = true;
      gridGeom.computeBoundingSphere();
    }

    // section
    secMesh.visible = state.section !== 'none';
    if (secMesh.visible) {
      const sp = secGeom.attributes.position.array, sc = secGeom.attributes.color.array;
      for (let k = 0; k < fields.sec.length; k++) {
        const f = fields.sec[k], kk = 1 + eps * re(f.Y, c, s);
        const p = P3(f.r * kk, f.th, f.ph);
        sp[3 * k] = p[0]; sp[3 * k + 1] = p[1]; sp[3 * k + 2] = 0;
        let rgb;
        if (state.section === 'pres') rgb = diverging(re(f.p, c, s) / fields.pmax);
        else if (state.section === 'speed') rgb = sequential(Math.hypot(...f.vel.map((q) => re(q, c, s))) / fields.vmax);
        else rgb = sequential(Math.hypot(...f.vort.map((q) => re(q, c, s))) / fields.wmax);
        sc[3 * k] = rgb[0]; sc[3 * k + 1] = rgb[1]; sc[3 * k + 2] = rgb[2];
      }
      secGeom.attributes.position.needsUpdate = true; secGeom.attributes.color.needsUpdate = true;
    }

    // arrows
    const set = state.vectors === 'none' ? [] : fields.sets[state.vectors];
    arrows.count = set.length;
    for (let k = 0; k < set.length; k++) {
      const a = set[k], kk = 1 + eps * re(a.Y, c, s);
      const v = a.vel.map((q) => re(q, c, s));
      const len = (state.vscale * eps * Math.hypot(...v)) / fields.lamAbs;
      const p = P3(a.r * kk, a.th, a.ph);
      if (len < 1e-4) { tmpM.makeScale(0, 0, 0); }
      else {
        dirV.set(v[0], v[1], v[2]).normalize();
        tmpQ.setFromUnitVectors(up, dirV);
        const w = Math.min(1, Math.max(0.3, len * 3));   // thinner arrows when short
        tmpM.compose(tmpP.set(p[0], p[1], p[2]), tmpQ, tmpS.set(w, len, w));
      }
      arrows.setMatrixAt(k, tmpM);
    }
    arrows.instanceMatrix.needsUpdate = true;
  }
  headlight.position.copy(camera.position); headlight.target.position.copy(controls.target);
  renderer.setScissorTest(false);
  renderer.render(scene, camera);
  // orientation axes
  const size = Math.round(Math.min(150, 0.22 * Math.min(renderer.domElement.clientWidth, renderer.domElement.clientHeight)));
  const dir = camera.position.clone().sub(controls.target).normalize();
  axCam.position.copy(dir.multiplyScalar(5.6)); axCam.up.copy(camera.up); axCam.lookAt(0, 0, 0);
  axLight.position.copy(axCam.position);
  renderer.autoClear = false; renderer.clearDepth();
  renderer.setViewport(8, 8, size, size); renderer.setScissor(8, 8, size, size); renderer.setScissorTest(true);
  renderer.render(axScene, axCam);
  renderer.setScissorTest(false); renderer.autoClear = true;
  const v = document.getElementById('view'); renderer.setViewport(0, 0, v.clientWidth, v.clientHeight);
  requestAnimationFrame(frame);
}

function resize() {
  const v = document.getElementById('view');
  const w = v.clientWidth, h = v.clientHeight;
  renderer.setSize(w, h, false); camera.aspect = w / h; camera.updateProjectionMatrix();
}
window.addEventListener('resize', resize);

// ------------------------------------------------------------------ physics summary
function oh() {
  if (state.useOh) return state.Oh;
  return (state.mu * 1e-3) / Math.sqrt(state.rho * state.sigma * state.R * 1e-3);
}
const fmt = (x, d = 4) => (Math.abs(x) >= 1e4 || (Math.abs(x) < 1e-3 && x !== 0)) ? x.toExponential(d - 1) : x.toPrecision(d);
const pct = (x) => (x >= 0 ? '+' : '') + (100 * x).toFixed(Math.abs(x) < 0.001 ? 4 : 2) + ' %';

function results() {
  const T = t(), l = state.l, lam = mode.lambda, Oh = mode.Oh;
  const A = (l - 1) * (2 * l + 1), B = l * (l - 1) * (l + 2);
  const rows = [];
  rows.push(['Oh', fmt(Oh, 4)]);
  const aperiodic = lam.im === 0;
  if (!state.useOh) {
    const R = state.R * 1e-3, rho = state.rho, mu = state.mu * 1e-3, sigma = state.sigma;
    const tv = (rho * R * R) / mu;
    const G = -lam.re / tv, f = lam.im / (2 * Math.PI * tv);
    const GL = (A * mu) / (rho * R * R), wR = Math.sqrt((B * sigma) / (rho * R ** 3)), fR = wR / (2 * Math.PI);
    rows.push([T.Gamma, fmt(G) + ' s⁻¹']);
    rows.push([T.freq, aperiodic ? '—' : fmt(f) + ' Hz']);
    rows.push([T.tau, fmt(1000 / G) + ' ms']);
    rows.push([T.lamb, fmt(GL) + ' s⁻¹ <span class="note">(' + T.dev + ' ' + pct(G / GL - 1) + ')</span>']);
    rows.push([T.rayleigh, fmt(fR) + ' Hz <span class="note">(' + T.dev + ' ' + pct(f / fR - 1) + ')</span>']);
    rows.push([T.nuLamb, T.err + ' ' + pct(G / GL - 1)]);
    if (!aperiodic) rows.push([T.sigRay, T.err + ' ' + pct((lam.im / (wR * tv)) ** 2 - 1)]);
  } else {
    const wR = Math.sqrt(B) / Oh;
    rows.push([T.Gamma + ' R²/ν', fmt(-lam.re, 6)]);
    rows.push(['ω R²/ν', aperiodic ? '—' : fmt(lam.im, 6)]);
    rows.push([T.lamb + ' R²/ν', fmt(A, 4) + ' <span class="note">(' + T.dev + ' ' + pct(-lam.re / A - 1) + ')</span>']);
    if (!aperiodic) rows.push(['ω_R R²/ν', fmt(wR, 6) + ' <span class="note">(' + T.dev + ' ' + pct(lam.im / wR - 1) + ')</span>']);
  }
  let html = '<table class="res">' + rows.map(([k, v]) => `<tr><td>${k}</td><td class="v">${v}</td></tr>`).join('') + '</table>';
  if (aperiodic) html += `<div class="warn">${T.aperiodic}</div>`;
  html += `<div class="note">λ R²/ν = ${lam.re.toPrecision(12)} ${lam.im >= 0 ? '+' : '−'} ${Math.abs(lam.im).toPrecision(12)} i</div>`;
  return html;
}

// ------------------------------------------------------------------ panel
const panel = document.getElementById('panel');
const toolbar = document.getElementById('toolbar');
function opt(v, label, cur) { return `<option value="${v}"${String(v) === String(cur) ? ' selected' : ''}>${label}</option>`; }
function render() {
  const T = t();
  document.documentElement.lang = state.lang; document.title = T.title;
  const mats = MATERIALS.map((m) => opt(m.key, `${m.key} (${m.T} K)`, state.material)).join('') + opt('custom', T.custom, state.material);
  const ls = [2, 3, 4, 5, 6, 7, 8].map((l) => opt(l, l, state.l)).join('');
  const ms = Array.from({ length: state.l + 1 }, (_, m) => opt(m, m, state.m)).join('');
  panel.innerHTML = `
  <div class="lang"><button data-lang="en" class="${state.lang === 'en' ? 'on' : ''}">EN</button> <button data-lang="ja" class="${state.lang === 'ja' ? 'on' : ''}">日本語</button></div>
  <h1>${T.title}</h1><div class="sub">${T.subtitle}</div>
  <h2>${T.drop}</h2>
  <label class="check"><input type="checkbox" id="useOh" ${state.useOh ? 'checked' : ''}>${T.useOh}</label>
  <div id="physBox" style="display:${state.useOh ? 'none' : 'block'}">
    <label>${T.material}</label><select id="material">${mats}</select>
    <div class="row"><div><label>${T.rho}</label><input type="number" id="rho" value="${state.rho}" step="any"></div>
      <div><label>${T.mu}</label><input type="number" id="mu" value="${state.mu}" step="any"></div></div>
    <div class="row"><div><label>${T.sigma}</label><input type="number" id="sigma" value="${state.sigma}" step="any"></div>
      <div><label>${T.R}</label><input type="number" id="R" value="${state.R}" step="any"></div></div>
  </div>
  <div id="ohBox" style="display:${state.useOh ? 'block' : 'none'}"><label>${T.Oh}</label><input type="number" id="Oh" value="${state.Oh}" step="any"></div>
  <h2>${T.mode}</h2>
  <div class="row"><div><label>${T.l}</label><select id="l">${ls}</select></div><div><label>${T.m}</label><select id="m">${ms}</select></div></div>
  <div class="note">${T.mNote}</div>
  <h2>${T.results}</h2><div id="results">${mode ? results() : T.computing}</div>
  <h2>${T.advanced}</h2>
  <label>${T.Nr}</label><select id="Nr">${[32, 40, 48, 64, 80].map((n) => opt(n, n, state.Nr)).join('')}</select>
  <footer>${T.note}<br><br>LAMBDA · <a href="https://github.com/sugurushiratori/Kawasemi-Lambda">github.com/sugurushiratori/Kawasemi-Lambda</a></footer>`;
  toolbar.innerHTML = `
  <div class="grp"><button class="play" id="play" title="${state.playing ? T.pause : T.play}">${state.playing ? '❚❚' : '▶'}</button>
    <label>${T.speed}<input type="range" id="speed" min="0.05" max="2" step="0.05" value="${state.speed}"></label>
    <label class="check"><input type="checkbox" id="damping" ${state.damping ? 'checked' : ''}>${T.damping}</label></div>
  <div class="grp"><label>${T.amp} <b id="ampV">${state.amp.toFixed(2)}</b> R<input type="range" id="amp" min="0" max="0.35" step="0.01" value="${state.amp}"></label></div>
  <div class="grp"><label>${T.surfColor}<select id="surfColor">${opt('disp', T.sDisp, state.surfColor)}${opt('pres', T.sPres, state.surfColor)}${opt('none', T.sNone, state.surfColor)}</select></label>
    <label>${T.opacity}<input type="range" id="opacity" min="0" max="1" step="0.05" value="${state.opacity}"></label>
    <label class="check"><input type="checkbox" id="grid" ${state.grid ? 'checked' : ''}>${T.grid}</label></div>
  <div class="grp"><label>${T.vectors}<select id="vectors">${opt('none', T.vNone, state.vectors)}${opt('section', T.vSection, state.vectors)}${opt('volume', T.vVolume, state.vectors)}${opt('surface', T.vSurface, state.vectors)}</select></label>
    <label>${T.vscale}<input type="range" id="vscale" min="0.5" max="10" step="0.5" value="${state.vscale}"></label></div>
  <div class="grp"><label>${T.section}<select id="section">${opt('none', T.secNone, state.section)}${opt('speed', T.secSpeed, state.section)}${opt('pres', T.secPres, state.section)}${opt('vort', T.secVort, state.section)}</select></label>
    <button class="ghost" id="resetView">${T.resetView}</button></div>`;
  bind();
  overlay();
}

function scalarBar(title, lo, mid, hi) {
  return `<div class="sbar"><div class="sbt">${title}</div><div class="sbb"><div class="sbg" style="background:${CTW_CSS}"></div>
    <div class="sbl"><span>${hi}</span><span>${mid}</span><span>${lo}</span></div></div></div>`;
}
function overlay() {
  const T = t();
  const bars = [];
  if (state.surfColor === 'disp') bars.push(scalarBar(T.barDisp, '−1', '0', '1'));
  if (state.surfColor === 'pres') bars.push(scalarBar(T.barPresS, '−1', '0', '1'));
  if (state.section === 'pres') bars.push(scalarBar(T.barPres, '−1', '0', '1'));
  if (state.section === 'speed') bars.push(scalarBar(T.barSpeed, '0', '0.5', '1'));
  if (state.section === 'vort') bars.push(scalarBar(T.barVort, '0', '0.5', '1'));
  document.getElementById('bars').innerHTML = bars.join('');
  document.getElementById('overlay').innerHTML = T.mouse;
}

let timer = null;
function scheduleSolve() {
  clearTimeout(timer);
  document.getElementById('status').textContent = t().computing;
  timer = setTimeout(solve, 120);
}
function solve() {
  const Oh = oh();
  const ok = state.useOh ? Oh > 0 : state.rho > 0 && state.mu > 0 && state.sigma > 0 && state.R > 0;
  if (!ok || !isFinite(Oh)) { document.getElementById('results').innerHTML = `<div class="warn">${t().invalid}</div>`; document.getElementById('status').textContent = ''; return; }
  const t0 = performance.now();
  mode = solveMode(state.l, Oh, state.Nr);
  precompute();
  phase = 0; decayArg = 0;
  document.getElementById('results').innerHTML = results();
  document.getElementById('status').textContent = `ℓ = ${state.l}, m = ${state.m} · ${(performance.now() - t0).toFixed(0)} ms`;
}

function bind() {
  const $ = (id) => document.getElementById(id);
  panel.querySelectorAll('[data-lang]').forEach((b) => b.onclick = () => { state.lang = b.dataset.lang; render(); });
  $('useOh').onchange = (e) => { state.useOh = e.target.checked; if (state.useOh) state.Oh = Number(oh().toPrecision(4)); render(); scheduleSolve(); };
  $('material').onchange = (e) => {
    state.material = e.target.value;
    const m = MATERIALS.find((x) => x.key === state.material);
    if (m) { state.rho = m.rho; state.mu = m.mu; state.sigma = m.sigma; }
    render(); scheduleSolve();
  };
  for (const k of ['rho', 'mu', 'sigma', 'R']) $(k).oninput = (e) => { state[k] = Number(e.target.value); state.material = k === 'R' ? state.material : 'custom'; $('material').value = state.material; scheduleSolve(); };
  $('Oh').oninput = (e) => { state.Oh = Number(e.target.value); scheduleSolve(); };
  $('l').onchange = (e) => { state.l = Number(e.target.value); state.m = Math.min(state.m, state.l); render(); scheduleSolve(); };
  $('m').onchange = (e) => { state.m = Number(e.target.value); if (mode) { precompute(); document.getElementById('status').textContent = `ℓ = ${state.l}, m = ${state.m}`; } };
  $('Nr').onchange = (e) => { state.Nr = Number(e.target.value); scheduleSolve(); };
  $('amp').oninput = (e) => { state.amp = Number(e.target.value); $('ampV').textContent = state.amp.toFixed(2); };
  $('play').onclick = () => { state.playing = !state.playing; $('play').textContent = state.playing ? '❚❚' : '▶'; $('play').title = state.playing ? t().pause : t().play; };
  $('grid').onchange = (e) => { state.grid = e.target.checked; };
  $('resetView').onclick = () => { camera.position.set(3.2, 1.6, 4.2); controls.target.set(0, 0, 0); controls.update(); };
  $('damping').onchange = (e) => { state.damping = e.target.checked; decayArg = 0; };
  $('speed').oninput = (e) => { state.speed = Number(e.target.value); };
  $('surfColor').onchange = (e) => { state.surfColor = e.target.value; overlay(); };
  $('opacity').oninput = (e) => { state.opacity = Number(e.target.value); };
  $('vectors').onchange = (e) => { state.vectors = e.target.value; };
  $('vscale').oninput = (e) => { state.vscale = Number(e.target.value); };
  $('section').onchange = (e) => { state.section = e.target.value; overlay(); };
}

render();
resize();
solve();
requestAnimationFrame(frame);
