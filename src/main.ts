import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { bindUI, readParams } from "./ui";
import { estimateMaxCurvature, samplePolyline } from "./lissajous";
import { downloadSTL, exportSTLBinary } from "./stl";
import { runSelfTests } from "./selftest";
import MeshWorker from "./worker?worker";
import type { WorkerRequest, WorkerResponse } from "./worker";
import type { Params, ValidationResult } from "./types";

const canvas = document.getElementById("viewport") as HTMLCanvasElement;
const warningsEl = document.getElementById("warnings") as HTMLElement;
const statsEl = document.getElementById("stats-body") as HTMLElement;
const buildBtn = document.getElementById("build") as HTMLButtonElement;
const exportBtn = document.getElementById("export") as HTMLButtonElement;

const renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
renderer.setClearColor(0x0b0d10, 1);

const scene = new THREE.Scene();

const camera = new THREE.PerspectiveCamera(45, 1, 0.01, 100);
camera.position.set(2.4, 2.0, 2.6);

const controls = new OrbitControls(camera, canvas);
controls.enableDamping = true;
controls.dampingFactor = 0.08;

// Ambient pur pour garantir un plancher de luminosité — évite les faces
// complètement noires quand les normales pointent "ailleurs" (non-manifold/
// non-orientable : on préfère voir la matière, même mal orientée).
scene.add(new THREE.AmbientLight(0xffffff, 0.8));

// HemisphereLight : sky/ground gris clairs, les deux côtés restent visibles.
const hemi = new THREE.HemisphereLight(0xffffff, 0x707580, 1.0);
scene.add(hemi);

const keyLight = new THREE.DirectionalLight(0xffffff, 1.5);
keyLight.position.set(3, 4, 2);
scene.add(keyLight);

const fillLight = new THREE.DirectionalLight(0xffffff, 0.6);
fillLight.position.set(-2, -1, -3);
scene.add(fillLight);

const grid = new THREE.GridHelper(4, 20, 0x2a2f36, 0x1a1d21);
scene.add(grid);
const axes = new THREE.AxesHelper(1.2);
scene.add(axes);

const curveGeom = new THREE.BufferGeometry();
curveGeom.setAttribute("position", new THREE.BufferAttribute(new Float32Array(0), 3));
const curveLine = new THREE.Line(curveGeom, new THREE.LineBasicMaterial({ color: 0x4ab0ff }));
scene.add(curveLine);

const meshGeom = new THREE.BufferGeometry();
// flatShading: les normales sont dérivées par triangle directement dans le
// fragment shader (dFdx/dFdy des positions). Immunise contre les vertices
// non-orientables qui produisent des normales NaN en smooth shading.
const meshMat = new THREE.MeshStandardMaterial({
  color: 0xc8cdd3,
  metalness: 0.1,
  roughness: 0.6,
  flatShading: true,
  side: THREE.DoubleSide,
});
const meshObj = new THREE.Mesh(meshGeom, meshMat);
meshObj.visible = false;
scene.add(meshObj);

function sanitizeNormals(geom: THREE.BufferGeometry): void {
  const attr = geom.getAttribute("normal") as THREE.BufferAttribute | undefined;
  if (!attr) return;
  const arr = attr.array as Float32Array;
  for (let i = 0; i < arr.length; i += 3) {
    const x = arr[i], y = arr[i + 1], z = arr[i + 2];
    const bad = !Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)
      || (x * x + y * y + z * z) < 1e-10;
    if (bad) {
      arr[i] = 0; arr[i + 1] = 1; arr[i + 2] = 0;
    }
  }
  attr.needsUpdate = true;
}

function setWarnings(messages: { text: string; level: "warn" | "err" }[]): void {
  warningsEl.innerHTML = "";
  for (const m of messages) {
    const div = document.createElement("div");
    div.className = m.level === "err" ? "warn err" : "warn";
    div.textContent = m.text;
    warningsEl.appendChild(div);
  }
}

// h_cible = min(R, max(k, R/3)) / 3, plancher par 2/résolution
function computeHTarget(params: Params): number {
  const floor = Math.max(params.k, params.R / 3);
  const raw = Math.min(params.R, floor) / 3;
  const resCap = 2 / params.resolution;
  return Math.max(raw, resCap);
}

function updateCurve(params: Params): void {
  const pts = samplePolyline(params, computeHTarget(params));
  const attr = curveGeom.getAttribute("position") as THREE.BufferAttribute;
  if (attr.array.length !== pts.length) {
    curveGeom.setAttribute("position", new THREE.BufferAttribute(pts, 3));
  } else {
    (attr.array as Float32Array).set(pts);
    attr.needsUpdate = true;
  }
  curveGeom.computeBoundingSphere();
}

function updateReachWarnings(params: Params): void {
  const kmax = estimateMaxCurvature(params);
  const reachLimit = 0.8 / Math.max(kmax, 1e-6);
  const warnings: { text: string; level: "warn" | "err" }[] = [];
  if (params.R > reachLimit) {
    warnings.push({
      text: `R=${params.R.toFixed(3)} > 0.8/κ_max=${reachLimit.toFixed(3)} — tube se recouvre hors croisements.`,
      level: "err",
    });
  } else if (params.R > 0.6 / kmax) {
    warnings.push({
      text: `R approche 0.8/κ_max (=${reachLimit.toFixed(3)}).`,
      level: "warn",
    });
  }
  setWarnings(warnings);
}

// ── Worker, debounce, requêtes id-taggées ──────────────────────────────
const worker = new MeshWorker();
let nextReqId = 1;
let inFlightId = 0;
let pendingParams: Params | null = null;
let debounceTimer: number | null = null;
let lastMesh: { positions: Float32Array; indices: Uint32Array; validation: ValidationResult } | null = null;

function ok(b: boolean): string {
  return b ? "✓" : "✗";
}

function renderStats(params: Params, r: WorkerResponse | null): void {
  const kmax = estimateMaxCurvature(params);
  const header = [
    `κ_max  = ${kmax.toFixed(3)}`,
    `h      = ${computeHTarget(params).toFixed(4)}`,
    r ? `grid   = ${r.gridN}³` : `grid   = ${params.resolution}³ (pending)`,
  ];
  if (!r) {
    statsEl.textContent = [...header, "", "computing…"].join("\n");
    return;
  }
  const v = r.validation;
  statsEl.textContent = [
    ...header,
    `V = ${v.V}    E = ${v.E}    F = ${v.F}`,
    `χ = V−E+F = ${v.euler}`,
    `watertight ${ok(v.watertight)}  manifold ${ok(v.manifold)}  orientable ${ok(v.orientable)}`,
    `degenerate: ${v.degenerateCount}`,
    ``,
    `sample  : ${r.timings.sample.toFixed(1)} ms`,
    `sdf     : ${r.timings.sdf.toFixed(1)} ms`,
    `dc      : ${r.timings.dc.toFixed(1)} ms`,
    `validate: ${r.timings.validate.toFixed(1)} ms`,
    `total   : ${r.timings.total.toFixed(1)} ms`,
  ].join("\n");
}

worker.onerror = (ev) => {
  // eslint-disable-next-line no-console
  console.error("[worker]", ev.message ?? ev);
  inFlightId = 0;
  if (pendingParams) {
    const next = pendingParams;
    pendingParams = null;
    dispatchBuild(next);
  }
};

worker.onmessage = (ev: MessageEvent<WorkerResponse>) => {
  const r = ev.data;
  if (r.id !== inFlightId) return; // réponse obsolète
  inFlightId = 0;

  meshGeom.setAttribute("position", new THREE.BufferAttribute(r.positions, 3));
  meshGeom.setIndex(new THREE.BufferAttribute(r.indices, 1));
  meshGeom.computeVertexNormals();
  // Sanitize : si le mesh est non-orientable, computeVertexNormals peut produire
  // des normales NaN ou nulles (faces voisines qui s'annulent). Remplace par un
  // fallback arbitraire pour éviter les pixels noirs dans le shader.
  sanitizeNormals(meshGeom);
  meshGeom.computeBoundingSphere();
  meshObj.visible = true;

  lastMesh = { positions: r.positions, indices: r.indices, validation: r.validation };
  renderStats(readParams(), r);

  // Export toujours autorisé (même mesh invalide). Tooltip indicatif seulement.
  const printable = r.validation.watertight && r.validation.manifold;
  exportBtn.title = printable
    ? "Mesh watertight + manifold, imprimable."
    : "Mesh non-manifold/non-watertight : exportable mais probablement à réparer (Blender 3D Print Toolbox, Meshmixer) avant impression.";

  if (pendingParams) {
    const next = pendingParams;
    pendingParams = null;
    dispatchBuild(next);
  }
};

function dispatchBuild(params: Params): void {
  if (inFlightId !== 0) {
    pendingParams = params;
    return;
  }
  const id = nextReqId++;
  inFlightId = id;
  renderStats(params, null);
  const req: WorkerRequest = { id, params, hTarget: computeHTarget(params) };
  worker.postMessage(req);
}

function scheduleBuild(params: Params, delay: number): void {
  if (debounceTimer !== null) clearTimeout(debounceTimer);
  debounceTimer = window.setTimeout(() => {
    debounceTimer = null;
    dispatchBuild(params);
  }, delay);
}

// ── Wiring ──────────────────────────────────────────────────────────────
function resize(): void {
  const w = canvas.clientWidth;
  const h = canvas.clientHeight;
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
}
window.addEventListener("resize", resize);
resize();

bindUI((params) => {
  updateCurve(params);
  updateReachWarnings(params);
  renderStats(params, null);
  meshObj.visible = false;
  scheduleBuild(params, 200);
});

buildBtn.addEventListener("click", () => {
  if (debounceTimer !== null) { clearTimeout(debounceTimer); debounceTimer = null; }
  dispatchBuild(readParams());
});

exportBtn.addEventListener("click", () => {
  if (!lastMesh) return;
  const params = readParams();
  const buffer = exportSTLBinary(lastMesh.positions, lastMesh.indices, params.scaleMm);
  const fmt = (n: number): string => n.toFixed(2).replace(".", "_");
  const name = `lissajous_p${fmt(params.p)}_q${fmt(params.q)}_r${fmt(params.r)}_R${fmt(params.R)}_${params.resolution}_mode${params.mode}.stl`;
  downloadSTL(buffer, name);
});

const initialParams = readParams();
updateCurve(initialParams);
updateReachWarnings(initialParams);
renderStats(initialParams, null);
dispatchBuild(initialParams);

if (import.meta.env.DEV) {
  runSelfTests();
}

function frame(): void {
  controls.update();
  renderer.render(scene, camera);
  requestAnimationFrame(frame);
}
frame();
