import {
  curvature,
  estimateMaxCurvature,
  gamma,
  gammaDoublePrime,
  gammaPrime,
} from "./lissajous";
import { buildSDF } from "./sdf";
import { dualContour } from "./dualcontouring";
import { validateMesh } from "./validate";
import type { Grid, Params } from "./types";

function finiteDiff1(
  t: number,
  params: Params,
  axis: "x" | "y" | "z",
  h = 1e-5,
): number {
  return (gamma(t + h, params)[axis] - gamma(t - h, params)[axis]) / (2 * h);
}

function finiteDiff2(
  t: number,
  params: Params,
  axis: "x" | "y" | "z",
  h = 1e-4,
): number {
  const fp = gamma(t + h, params)[axis];
  const f0 = gamma(t, params)[axis];
  const fm = gamma(t - h, params)[axis];
  return (fp - 2 * f0 + fm) / (h * h);
}

function voxelCenter(grid: Grid, i: number, j: number, k: number): [number, number, number] {
  return [
    grid.originX + (i + 0.5) * grid.h,
    grid.originY + (j + 0.5) * grid.h,
    grid.originZ + (k + 0.5) * grid.h,
  ];
}

function fieldAt(field: Float32Array, grid: Grid, i: number, j: number, k: number): number {
  return field[i + j * grid.nx + k * grid.nx * grid.ny];
}

function distPointSegment(
  px: number, py: number, pz: number,
  ax: number, ay: number, az: number,
  bx: number, by: number, bz: number,
): number {
  const vx = bx - ax, vy = by - ay, vz = bz - az;
  const wx = px - ax, wy = py - ay, wz = pz - az;
  const vv = vx * vx + vy * vy + vz * vz;
  let t = vv < 1e-20 ? 0 : (wx * vx + wy * vy + wz * vz) / vv;
  if (t < 0) t = 0;
  else if (t > 1) t = 1;
  const ex = wx - t * vx, ey = wy - t * vy, ez = wz - t * vz;
  return Math.sqrt(ex * ex + ey * ey + ez * ez);
}

function nearestVoxel(grid: Grid, x: number, y: number, z: number): [number, number, number] {
  const i = Math.max(0, Math.min(grid.nx - 1, Math.round((x - grid.originX) / grid.h - 0.5)));
  const j = Math.max(0, Math.min(grid.ny - 1, Math.round((y - grid.originY) / grid.h - 0.5)));
  const k = Math.max(0, Math.min(grid.nz - 1, Math.round((z - grid.originZ) / grid.h - 0.5)));
  return [i, j, k];
}

function lissajousTests(): Record<string, number> {
  const params: Params = {
    p: 5, q: 4, r: 3,
    phix: 0.3, phiy: 0.7, phiz: 1.1,
    L: 2,
    R: 0.15, k: 0,
    mode: "A",
    resolution: 128,
    scaleMm: 80,
  };

  let maxErr1 = 0;
  let maxErr2 = 0;
  for (let i = 0; i < 50; i++) {
    const t = (i / 50) * Math.PI * 2;
    const g1 = gammaPrime(t, params);
    const g2 = gammaDoublePrime(t, params);
    maxErr1 = Math.max(
      maxErr1,
      Math.abs(g1.x - finiteDiff1(t, params, "x")),
      Math.abs(g1.y - finiteDiff1(t, params, "y")),
      Math.abs(g1.z - finiteDiff1(t, params, "z")),
    );
    maxErr2 = Math.max(
      maxErr2,
      Math.abs(g2.x - finiteDiff2(t, params, "x")),
      Math.abs(g2.y - finiteDiff2(t, params, "y")),
      Math.abs(g2.z - finiteDiff2(t, params, "z")),
    );
  }

  return {
    gammaPrimeErrMax: maxErr1,
    gammaDoublePrimeErrMax: maxErr2,
    curvatureAtT: curvature(0.42, params),
    kmaxCoarse: estimateMaxCurvature(params, 500),
    kmaxDense: estimateMaxCurvature(params, 20000),
  };
}

function sdfSingleSegmentTest(): {
  maxErrInBand: number;
  sampled: number;
  farVoxelsReported: number;
} {
  // Segment [(-0.5,0,0), (+0.5,0,0)], R=0.1, mode A.
  const ax = -0.5, ay = 0, az = 0, bx = 0.5, by = 0, bz = 0;
  const polyline = new Float32Array([ax, ay, az, bx, by, bz]);
  const R = 0.1;
  const grid: Grid = {
    nx: 96, ny: 48, nz: 48,
    h: 0.03,
    originX: -1.44, originY: -0.72, originZ: -0.72,
  };
  const field = buildSDF(polyline, R, 0, grid);

  // Pour chaque voxel dans la bande (d ≤ R+3h), comparer field[idx] à d-R calculé.
  const band = R + 3 * grid.h;
  let maxErr = 0;
  let sampled = 0;
  let farCount = 0;
  for (let k = 0; k < grid.nz; k++) {
    for (let j = 0; j < grid.ny; j++) {
      for (let i = 0; i < grid.nx; i++) {
        const [cx, cy, cz] = voxelCenter(grid, i, j, k);
        const dExact = distPointSegment(cx, cy, cz, ax, ay, az, bx, by, bz);
        const expected = dExact - R;
        const got = fieldAt(field, grid, i, j, k);
        if (dExact <= band) {
          sampled++;
          maxErr = Math.max(maxErr, Math.abs(got - expected));
        } else if (got < 1e9) {
          farCount++;
        }
      }
    }
  }
  return { maxErrInBand: maxErr, sampled, farVoxelsReported: farCount };
}

function sdfSphereTest(): {
  maxErrInBand: number;
  sampled: number;
} {
  // Polyligne dégénérée (point). SDF = |p| - R.
  const polyline = new Float32Array([0, 0, 0, 0, 0, 0]);
  const R = 0.3;
  const grid: Grid = {
    nx: 64, ny: 64, nz: 64,
    h: 0.04,
    originX: -1.28, originY: -1.28, originZ: -1.28,
  };
  const field = buildSDF(polyline, R, 0, grid);

  const band = R + 3 * grid.h;
  let maxErr = 0;
  let sampled = 0;
  for (let k = 0; k < grid.nz; k++) {
    for (let j = 0; j < grid.ny; j++) {
      for (let i = 0; i < grid.nx; i++) {
        const [cx, cy, cz] = voxelCenter(grid, i, j, k);
        const dExact = Math.sqrt(cx * cx + cy * cy + cz * cz);
        if (dExact <= band) {
          sampled++;
          maxErr = Math.max(maxErr, Math.abs(fieldAt(field, grid, i, j, k) - (dExact - R)));
        }
      }
    }
  }
  return { maxErrInBand: maxErr, sampled };
}

function sdfSmoothModeTest(): Record<string, number> {
  // Segment X seul. fB ≡ fA partout (un seul segment ⇒ smin = min).
  const R = 0.1;
  const k = 0.05;
  const polylineX = new Float32Array([-0.5, 0, 0, 0.5, 0, 0]);
  const grid: Grid = {
    nx: 64, ny: 64, nz: 64,
    h: 0.03,
    originX: -0.96, originY: -0.96, originZ: -0.96,
  };
  const fieldA = buildSDF(polylineX, R, 0, grid);
  const fieldB = buildSDF(polylineX, R, k, grid);

  let maxDiff = 0;
  for (let idx = 0; idx < fieldA.length; idx++) {
    if (fieldA[idx] < 1e9) maxDiff = Math.max(maxDiff, Math.abs(fieldA[idx] - fieldB[idx]));
  }

  // Deux segments croisés (mode A) : comparer f_A au min des deux SDF individuelles.
  const polylineY = new Float32Array([0, -0.5, 0, 0, 0.5, 0]);
  const combined = new Float32Array([-0.5, 0, 0, 0.5, 0, 0, 0, -0.5, 0, 0, 0.5, 0]);
  const fX = buildSDF(polylineX, R, 0, grid);
  const fY = buildSDF(polylineY, R, 0, grid);
  const fAB = buildSDF(combined, R, 0, grid);
  // Note : combined contient un segment parasite [(+0.5,0,0)→(0,-0.5,0)] qu'on
  // ne peut pas éliminer sans une API multi-polyligne. On sonde uniquement les
  // voxels proches de l'origine, où l'AABB du segment parasite est également
  // présent et donc pollue potentiellement le test. On accepte fAB ≤ min(fX,fY).
  const [i0, j0, k0] = nearestVoxel(grid, 0, 0, 0);
  const fCenterX = fieldAt(fX, grid, i0, j0, k0);
  const fCenterY = fieldAt(fY, grid, i0, j0, k0);
  const fCenterAB = fieldAt(fAB, grid, i0, j0, k0);

  // Mode B au centre : smin doit creuser par rapport à min. Deux segments
  // équidistants ⇒ smin ≈ d_min - k·log(2) - R (+ contribution parasite).
  const fCenterABsmooth = fieldAt(buildSDF(combined, R, k, grid), grid, i0, j0, k0);

  return {
    singleSegment_AvsB_maxDiff: maxDiff,   // doit être ≈ 0
    fCenterX,
    fCenterY,
    fCenterAB,                             // ≤ min(fCenterX, fCenterY)
    unionDeltaFromMin: fCenterAB - Math.min(fCenterX, fCenterY), // ≤ 0 attendu
    fCenterABsmooth,                       // < fCenterAB attendu
    smoothVsSharpAtCenter: fCenterABsmooth - fCenterAB, // < 0 attendu
  };
}

function dcSphereTest(): Record<string, number> {
  // SDF analytique d'une sphère : f(p) = |p| - R
  const R = 0.4;
  const grid: Grid = {
    nx: 48, ny: 48, nz: 48,
    h: 0.04,
    originX: -0.96, originY: -0.96, originZ: -0.96,
  };
  const field = new Float32Array(grid.nx * grid.ny * grid.nz);
  for (let k = 0; k < grid.nz; k++) {
    const z = grid.originZ + (k + 0.5) * grid.h;
    for (let j = 0; j < grid.ny; j++) {
      const y = grid.originY + (j + 0.5) * grid.h;
      for (let i = 0; i < grid.nx; i++) {
        const x = grid.originX + (i + 0.5) * grid.h;
        field[i + j * grid.nx + k * grid.nx * grid.ny] = Math.sqrt(x * x + y * y + z * z) - R;
      }
    }
  }

  const mesh = dualContour(field, grid);
  const V = mesh.positions.length / 3;
  const F = mesh.indices.length / 3;
  // Watertight ⇒ E = 3F/2 ⇒ χ = V - F/2.
  const eulerIfClosed = V - F / 2;

  // Distance max des sommets à la sphère de rayon R (doit être < h pour un mesh correct).
  let maxRadiusErr = 0;
  for (let i = 0; i < V; i++) {
    const x = mesh.positions[i * 3];
    const y = mesh.positions[i * 3 + 1];
    const z = mesh.positions[i * 3 + 2];
    const r = Math.sqrt(x * x + y * y + z * z);
    maxRadiusErr = Math.max(maxRadiusErr, Math.abs(r - R));
  }

  return {
    V, F,
    eulerIfClosed,
    maxRadiusErr,
    hGrid: grid.h,
  };
}

function dcCubeTest(): Record<string, number> {
  // SDF d'un cube [-a, a]³ : f(p) = max(|p|_∞) - a (approx pour l'intérieur).
  // Mode A-like : arêtes vives, QEF doit les préserver.
  const a = 0.3;
  const grid: Grid = {
    nx: 32, ny: 32, nz: 32,
    h: 0.05,
    originX: -0.8, originY: -0.8, originZ: -0.8,
  };
  const field = new Float32Array(grid.nx * grid.ny * grid.nz);
  for (let k = 0; k < grid.nz; k++) {
    const z = grid.originZ + (k + 0.5) * grid.h;
    for (let j = 0; j < grid.ny; j++) {
      const y = grid.originY + (j + 0.5) * grid.h;
      for (let i = 0; i < grid.nx; i++) {
        const x = grid.originX + (i + 0.5) * grid.h;
        // SDF du cube (Inigo Quilez)
        const dx = Math.abs(x) - a;
        const dy = Math.abs(y) - a;
        const dz = Math.abs(z) - a;
        const outside = Math.hypot(Math.max(dx, 0), Math.max(dy, 0), Math.max(dz, 0));
        const inside = Math.min(Math.max(dx, Math.max(dy, dz)), 0);
        field[i + j * grid.nx + k * grid.nx * grid.ny] = outside + inside;
      }
    }
  }
  const mesh = dualContour(field, grid);
  const V = mesh.positions.length / 3;
  const F = mesh.indices.length / 3;
  return { V, F, eulerIfClosed: V - F / 2 };
}

export function runSelfTests(): void {
  // eslint-disable-next-line no-console
  console.group("[selftest] lissajous");
  // eslint-disable-next-line no-console
  console.table(lissajousTests());
  // eslint-disable-next-line no-console
  console.groupEnd();

  // eslint-disable-next-line no-console
  console.group("[selftest] sdf single segment (maxErrInBand doit être ≈ 0)");
  // eslint-disable-next-line no-console
  console.table(sdfSingleSegmentTest());
  // eslint-disable-next-line no-console
  console.groupEnd();

  // eslint-disable-next-line no-console
  console.group("[selftest] sdf sphere (maxErrInBand doit être ≈ 0)");
  // eslint-disable-next-line no-console
  console.table(sdfSphereTest());
  // eslint-disable-next-line no-console
  console.groupEnd();

  // eslint-disable-next-line no-console
  console.group("[selftest] sdf smooth mode");
  // eslint-disable-next-line no-console
  console.table(sdfSmoothModeTest());
  // eslint-disable-next-line no-console
  console.groupEnd();

  // eslint-disable-next-line no-console
  console.group("[selftest] dc sphere (attendu : eulerIfClosed = 2, maxRadiusErr ≲ h)");
  // eslint-disable-next-line no-console
  console.table(dcSphereTest());
  // eslint-disable-next-line no-console
  console.groupEnd();

  // eslint-disable-next-line no-console
  console.group("[selftest] dc cube (attendu : eulerIfClosed = 2, arêtes vives)");
  // eslint-disable-next-line no-console
  console.table(dcCubeTest());
  // eslint-disable-next-line no-console
  console.groupEnd();

  // eslint-disable-next-line no-console
  console.group("[selftest] validate sphere (attendu : watertight/manifold/orientable=true, χ=2)");
  // eslint-disable-next-line no-console
  console.table(validateSphereTest());
  // eslint-disable-next-line no-console
  console.groupEnd();

  // eslint-disable-next-line no-console
  console.group("[selftest] validate broken (attendu : watertight=false)");
  // eslint-disable-next-line no-console
  console.table(validateBrokenTest());
  // eslint-disable-next-line no-console
  console.groupEnd();
}

function validateSphereTest(): Record<string, number | boolean> {
  const R = 0.4;
  const grid: Grid = {
    nx: 48, ny: 48, nz: 48,
    h: 0.04,
    originX: -0.96, originY: -0.96, originZ: -0.96,
  };
  const field = new Float32Array(grid.nx * grid.ny * grid.nz);
  for (let k = 0; k < grid.nz; k++) {
    const z = grid.originZ + (k + 0.5) * grid.h;
    for (let j = 0; j < grid.ny; j++) {
      const y = grid.originY + (j + 0.5) * grid.h;
      for (let i = 0; i < grid.nx; i++) {
        const x = grid.originX + (i + 0.5) * grid.h;
        field[i + j * grid.nx + k * grid.nx * grid.ny] = Math.sqrt(x * x + y * y + z * z) - R;
      }
    }
  }
  const mesh = dualContour(field, grid);
  const v = validateMesh(mesh.positions, mesh.indices);
  return {
    V: v.V, E: v.E, F: v.F,
    euler: v.euler,
    watertight: v.watertight,
    manifold: v.manifold,
    orientable: v.orientable,
    degenerateCount: v.degenerateCount,
  };
}

function validateBrokenTest(): Record<string, number | boolean> {
  // Construit un tétraèdre, retire une face ⇒ mesh ouvert (non-watertight).
  // 4 sommets, 3 triangles restants (au lieu de 4).
  const positions = new Float32Array([
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    0, 0, 1,
  ]);
  const indices = new Uint32Array([
    0, 2, 1,
    0, 1, 3,
    0, 3, 2,
    // triangle (1,2,3) retiré
  ]);
  const v = validateMesh(positions, indices);
  return {
    V: v.V, E: v.E, F: v.F,
    euler: v.euler,
    watertight: v.watertight,
    orientable: v.orientable,
  };
}
