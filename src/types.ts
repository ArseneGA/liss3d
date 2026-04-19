export type Mode = "A" | "B";

export interface Params {
  p: number; // fréquence en x
  q: number; // fréquence en y
  r: number; // fréquence en z
  A: number; // amplitude en x
  B: number; // amplitude en y
  C: number; // amplitude en z
  phix: number;
  phiy: number;
  phiz: number;
  L: number; // longueur de la courbe : t ∈ [0, L·π], L entier
  R: number; // rayon du tube
  k: number;
  mode: Mode;
  resolution: number;
  scaleMm: number;
}

export interface Grid {
  nx: number;
  ny: number;
  nz: number;
  h: number;
  originX: number;
  originY: number;
  originZ: number;
}

export interface Mesh {
  positions: Float32Array;
  indices: Uint32Array;
}

export interface SDFResult {
  field: Float32Array;
  grid: Grid;
}

export interface ValidationResult {
  V: number;
  E: number;
  F: number;
  euler: number;
  watertight: boolean;
  manifold: boolean;
  orientable: boolean;
  degenerateCount: number;
  selfIntersects?: boolean;
}

export interface Vec3 {
  x: number;
  y: number;
  z: number;
}
