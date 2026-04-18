import { samplePolyline } from "./lissajous";
import { buildGrid, buildSDF } from "./sdf";
import { dualContour } from "./dualcontouring";
import { validateMesh } from "./validate";
import type { Params, ValidationResult } from "./types";

export interface WorkerRequest {
  id: number;
  params: Params;
  hTarget: number;
}

export interface WorkerResponse {
  id: number;
  positions: Float32Array;
  indices: Uint32Array;
  validation: ValidationResult;
  timings: {
    sample: number;
    sdf: number;
    dc: number;
    validate: number;
    total: number;
  };
  gridN: number;
}

const ctx = self as unknown as {
  onmessage: ((ev: MessageEvent<WorkerRequest>) => void) | null;
  postMessage: (msg: WorkerResponse, transfer: Transferable[]) => void;
};

ctx.onmessage = (ev) => {
  const { id, params, hTarget } = ev.data;

  const t0 = performance.now();
  const polyline = samplePolyline(params, hTarget);
  const t1 = performance.now();
  const grid = buildGrid(params);
  const field = buildSDF(polyline, params.R, params.k, grid);
  const t2 = performance.now();
  const mesh = dualContour(field, grid);
  const t3 = performance.now();
  const validation = validateMesh(mesh.positions, mesh.indices);
  const t4 = performance.now();

  const response: WorkerResponse = {
    id,
    positions: mesh.positions,
    indices: mesh.indices,
    validation,
    timings: {
      sample: t1 - t0,
      sdf: t2 - t1,
      dc: t3 - t2,
      validate: t4 - t3,
      total: t4 - t0,
    },
    gridN: grid.nx,
  };

  ctx.postMessage(response, [response.positions.buffer, response.indices.buffer]);
};
