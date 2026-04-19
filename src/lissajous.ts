import type { Params, Vec3 } from "./types";

function tMax(params: Params): number {
  return params.L * Math.PI;
}

// γ(t) = (A·sin(p·t + φx·π), B·sin(q·t + φy·π), C·sin(r·t + φz·π))
export function gamma(t: number, params: Params): Vec3 {
  return {
    x: params.A * Math.sin(params.p * t + params.phix * Math.PI),
    y: params.B * Math.sin(params.q * t + params.phiy * Math.PI),
    z: params.C * Math.sin(params.r * t + params.phiz * Math.PI),
  };
}

// γ'(t) = (A·p·cos(…), B·q·cos(…), C·r·cos(…))
export function gammaPrime(t: number, params: Params): Vec3 {
  return {
    x: params.A * params.p * Math.cos(params.p * t + params.phix * Math.PI),
    y: params.B * params.q * Math.cos(params.q * t + params.phiy * Math.PI),
    z: params.C * params.r * Math.cos(params.r * t + params.phiz * Math.PI),
  };
}

// γ''(t) = (−A·p²·sin(…), −B·q²·sin(…), −C·r²·sin(…))
export function gammaDoublePrime(t: number, params: Params): Vec3 {
  return {
    x: -params.A * params.p * params.p * Math.sin(params.p * t + params.phix * Math.PI),
    y: -params.B * params.q * params.q * Math.sin(params.q * t + params.phiy * Math.PI),
    z: -params.C * params.r * params.r * Math.sin(params.r * t + params.phiz * Math.PI),
  };
}

// κ(t) = |γ' × γ''| / |γ'|³
export function curvature(t: number, params: Params): number {
  const d1 = gammaPrime(t, params);
  const d2 = gammaDoublePrime(t, params);
  const cx = d1.y * d2.z - d1.z * d2.y;
  const cy = d1.z * d2.x - d1.x * d2.z;
  const cz = d1.x * d2.y - d1.y * d2.x;
  const num = Math.hypot(cx, cy, cz);
  const denom = Math.hypot(d1.x, d1.y, d1.z);
  const d3 = denom * denom * denom;
  return d3 < 1e-12 ? 0 : num / d3;
}

export function estimateMaxCurvature(params: Params, nSamples = 2000): number {
  const tm = tMax(params);
  let kmax = 0;
  for (let i = 0; i < nSamples; i++) {
    const t = (i / nSamples) * tm;
    const k = curvature(t, params);
    if (k > kmax) kmax = k;
  }
  return kmax;
}

// Borne sup de |γ'(t)| : √((A·p)² + (B·q)² + (C·r)²).
export function maxSpeed(params: Params): number {
  return Math.hypot(params.A * params.p, params.B * params.q, params.C * params.r);
}

/**
 * Demi-extent maximal de la courbe en monde (max(A, B, C)). Utilisé pour
 * dimensionner la grille SDF.
 */
export function curveHalfExtent(params: Params): number {
  return Math.max(params.A, params.B, params.C);
}

/**
 * Polyligne sur t ∈ [0, L·π]. Si params.N > 0, on prend N points (N−1 segments).
 * Sinon, N est calculé automatiquement pour garantir erreur Hausdorff corde
 * < hTarget/4 sur la pire courbure.
 */
export function samplePolyline(params: Params, hTarget: number): Float32Array {
  const tm = tMax(params);
  let n: number; // nombre de segments (points = n + 1)
  if (params.N > 0) {
    n = Math.max(1, Math.floor(params.N - 1));
  } else {
    const kmax = Math.max(estimateMaxCurvature(params), 1e-6);
    const vmax = Math.max(maxSpeed(params), 1e-6);
    const chordMax = Math.sqrt((2 * hTarget) / kmax);
    n = Math.max(64, Math.ceil((tm * vmax) / chordMax));
  }
  const out = new Float32Array((n + 1) * 3);
  for (let i = 0; i <= n; i++) {
    const t = (i / n) * tm;
    const g = gamma(t, params);
    out[i * 3] = g.x;
    out[i * 3 + 1] = g.y;
    out[i * 3 + 2] = g.z;
  }
  return out;
}

export function samplePolylineUniform(params: Params, n: number): Float32Array {
  const tm = tMax(params);
  const out = new Float32Array((n + 1) * 3);
  for (let i = 0; i <= n; i++) {
    const t = (i / n) * tm;
    const g = gamma(t, params);
    out[i * 3] = g.x;
    out[i * 3 + 1] = g.y;
    out[i * 3 + 2] = g.z;
  }
  return out;
}
