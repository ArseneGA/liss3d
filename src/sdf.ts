import { curveHalfExtent } from "./lissajous";
import type { Grid, Params, SDFResult } from "./types";

const FAR = 1e10;

/**
 * Grille cubique centrée sur l'origine. Demi-côté = max(A, B, C) + R + 3h.
 * Résolution n = params.resolution (cap utilisateur).
 */
export function buildGrid(params: Params): Grid {
  const n = params.resolution;
  const half = Math.max(curveHalfExtent(params), 1e-3);
  const hEstimate = (2 * half) / n;
  const padding = params.R + 3 * hEstimate;
  const span = 2 * half + 2 * padding;
  const h = span / n;
  const origin = -(half + padding);
  return {
    nx: n, ny: n, nz: n,
    h,
    originX: origin, originY: origin, originZ: origin,
  };
}

/**
 * Rasterise chaque segment dans son AABB (étendu de R+3h) et pousse la distance
 * point-segment dans `dmin`. Retourne le nombre de voxels touchés (debug).
 */
function rasterizeDmin(
  polyline: Float32Array,
  R: number,
  grid: Grid,
  dmin: Float32Array,
): void {
  const { nx, ny, nz, h, originX, originY, originZ } = grid;
  const nSegments = polyline.length / 3 - 1;
  const expand = R + 3 * h;
  const nxny = nx * ny;

  for (let s = 0; s < nSegments; s++) {
    const ax = polyline[s * 3];
    const ay = polyline[s * 3 + 1];
    const az = polyline[s * 3 + 2];
    const bx = polyline[(s + 1) * 3];
    const by = polyline[(s + 1) * 3 + 1];
    const bz = polyline[(s + 1) * 3 + 2];

    const minX = Math.min(ax, bx) - expand;
    const maxX = Math.max(ax, bx) + expand;
    const minY = Math.min(ay, by) - expand;
    const maxY = Math.max(ay, by) + expand;
    const minZ = Math.min(az, bz) - expand;
    const maxZ = Math.max(az, bz) + expand;

    const i0 = Math.max(0, Math.floor((minX - originX) / h));
    const i1 = Math.min(nx - 1, Math.floor((maxX - originX) / h));
    const j0 = Math.max(0, Math.floor((minY - originY) / h));
    const j1 = Math.min(ny - 1, Math.floor((maxY - originY) / h));
    const k0 = Math.max(0, Math.floor((minZ - originZ) / h));
    const k1 = Math.min(nz - 1, Math.floor((maxZ - originZ) / h));

    if (i0 > i1 || j0 > j1 || k0 > k1) continue;

    const vx = bx - ax;
    const vy = by - ay;
    const vz = bz - az;
    const vDotV = vx * vx + vy * vy + vz * vz;
    const invVdV = vDotV < 1e-20 ? 0 : 1 / vDotV;

    for (let kv = k0; kv <= k1; kv++) {
      const pz = originZ + (kv + 0.5) * h;
      const slab = kv * nxny;
      const wzBase = pz - az;
      for (let jv = j0; jv <= j1; jv++) {
        const py = originY + (jv + 0.5) * h;
        const row = slab + jv * nx;
        const wyBase = py - ay;
        for (let iv = i0; iv <= i1; iv++) {
          const px = originX + (iv + 0.5) * h;
          const wx = px - ax;
          let t = (wx * vx + wyBase * vy + wzBase * vz) * invVdV;
          if (t < 0) t = 0;
          else if (t > 1) t = 1;
          const ex = wx - t * vx;
          const ey = wyBase - t * vy;
          const ez = wzBase - t * vz;
          const d = Math.sqrt(ex * ex + ey * ey + ez * ez);
          const idx = iv + row;
          if (d < dmin[idx]) dmin[idx] = d;
        }
      }
    }
  }
}

/**
 * Pass 2 mode B : accumule S = Σ exp(−(d_i − d_min)/k) dans S[].
 * Pruning : segments avec d − d_min ≥ k·20 contribuent < 2e-9, ignorés.
 */
function rasterizeSmoothSum(
  polyline: Float32Array,
  R: number,
  k: number,
  grid: Grid,
  dmin: Float32Array,
  S: Float32Array,
): void {
  const { nx, ny, nz, h, originX, originY, originZ } = grid;
  const nSegments = polyline.length / 3 - 1;
  const expand = R + 3 * h;
  const threshold = k * 20;
  const invK = 1 / k;
  const nxny = nx * ny;

  for (let s = 0; s < nSegments; s++) {
    const ax = polyline[s * 3];
    const ay = polyline[s * 3 + 1];
    const az = polyline[s * 3 + 2];
    const bx = polyline[(s + 1) * 3];
    const by = polyline[(s + 1) * 3 + 1];
    const bz = polyline[(s + 1) * 3 + 2];

    const minX = Math.min(ax, bx) - expand;
    const maxX = Math.max(ax, bx) + expand;
    const minY = Math.min(ay, by) - expand;
    const maxY = Math.max(ay, by) + expand;
    const minZ = Math.min(az, bz) - expand;
    const maxZ = Math.max(az, bz) + expand;

    const i0 = Math.max(0, Math.floor((minX - originX) / h));
    const i1 = Math.min(nx - 1, Math.floor((maxX - originX) / h));
    const j0 = Math.max(0, Math.floor((minY - originY) / h));
    const j1 = Math.min(ny - 1, Math.floor((maxY - originY) / h));
    const k0 = Math.max(0, Math.floor((minZ - originZ) / h));
    const k1 = Math.min(nz - 1, Math.floor((maxZ - originZ) / h));

    if (i0 > i1 || j0 > j1 || k0 > k1) continue;

    const vx = bx - ax;
    const vy = by - ay;
    const vz = bz - az;
    const vDotV = vx * vx + vy * vy + vz * vz;
    const invVdV = vDotV < 1e-20 ? 0 : 1 / vDotV;

    for (let kv = k0; kv <= k1; kv++) {
      const pz = originZ + (kv + 0.5) * h;
      const slab = kv * nxny;
      const wzBase = pz - az;
      for (let jv = j0; jv <= j1; jv++) {
        const py = originY + (jv + 0.5) * h;
        const row = slab + jv * nx;
        const wyBase = py - ay;
        for (let iv = i0; iv <= i1; iv++) {
          const px = originX + (iv + 0.5) * h;
          const wx = px - ax;
          let t = (wx * vx + wyBase * vy + wzBase * vz) * invVdV;
          if (t < 0) t = 0;
          else if (t > 1) t = 1;
          const ex = wx - t * vx;
          const ey = wyBase - t * vy;
          const ez = wzBase - t * vz;
          const d = Math.sqrt(ex * ex + ey * ey + ez * ez);
          const idx = iv + row;
          const delta = d - dmin[idx];
          if (delta < threshold) {
            S[idx] += Math.exp(-delta * invK);
          }
        }
      }
    }
  }
}

/**
 * f_k(p) = smin_k({d(p, seg_i)}) − R
 *
 * Mode A (k=0) : smin = min, arêtes vives.
 * Mode B (k>0) : smin = −k·log(Σ exp(−d_i/k)), log-sum-exp stabilisé par
 * soustraction de d_min avant exp.
 *
 * Layout : idx = i + j·nx + k·nx·ny.
 */
export function buildSDF(
  polyline: Float32Array,
  R: number,
  k: number,
  grid: Grid,
): Float32Array {
  const N = grid.nx * grid.ny * grid.nz;
  const dmin = new Float32Array(N);
  dmin.fill(FAR);

  rasterizeDmin(polyline, R, grid, dmin);

  const field = new Float32Array(N);
  if (k <= 0) {
    for (let i = 0; i < N; i++) field[i] = dmin[i] - R;
    return field;
  }

  const S = new Float32Array(N);
  rasterizeSmoothSum(polyline, R, k, grid, dmin, S);

  for (let i = 0; i < N; i++) {
    field[i] = S[i] > 0 ? dmin[i] - k * Math.log(S[i]) - R : dmin[i] - R;
  }
  return field;
}

export function buildSDFResult(
  polyline: Float32Array,
  params: Params,
  grid?: Grid,
): SDFResult {
  const g = grid ?? buildGrid(params);
  const field = buildSDF(polyline, params.R, params.k, g);
  return { field, grid: g };
}
