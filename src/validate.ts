import type { ValidationResult } from "./types";

// Encode (u, v) → numérique, u ∈ [0, 2^22). Produits ≤ 2^44 < 2^53 (Number safe).
const ENC = 1 << 22;

/**
 * Validation topologique d'un mesh triangulé.
 *
 * - V, E, F, χ = V − E + F
 * - watertight : chaque arête non-orientée apparaît exactement 2 fois
 * - orientable : chaque arête orientée (u→v) apparaît exactement 1 fois
 * - manifold : le link de chaque sommet est un cycle unique
 * - degenerate : triangles d'aire < 1e-10
 */
export function validateMesh(
  positions: Float32Array,
  indices: Uint32Array,
): ValidationResult {
  const V = positions.length / 3;
  const F = indices.length / 3;

  if (V >= ENC) {
    throw new Error(`validateMesh: V=${V} ≥ 2^22, encoding overflow`);
  }

  // Arêtes non-orientées (count) et orientées (count)
  const undirected = new Map<number, number>();
  const directed = new Map<number, number>();

  // CSR : sommet → triangles incidents
  const degree = new Int32Array(V);
  for (let f = 0; f < F; f++) {
    degree[indices[f * 3]]++;
    degree[indices[f * 3 + 1]]++;
    degree[indices[f * 3 + 2]]++;
  }
  const offset = new Int32Array(V + 1);
  for (let i = 0; i < V; i++) offset[i + 1] = offset[i] + degree[i];
  const incident = new Int32Array(offset[V]);
  const cursor = new Int32Array(V);
  for (let f = 0; f < F; f++) {
    for (let c = 0; c < 3; c++) {
      const vi = indices[f * 3 + c];
      incident[offset[vi] + cursor[vi]] = f;
      cursor[vi]++;
    }
  }

  let degenerateCount = 0;
  for (let f = 0; f < F; f++) {
    const i0 = indices[f * 3];
    const i1 = indices[f * 3 + 1];
    const i2 = indices[f * 3 + 2];

    // aire
    const ax = positions[i0 * 3], ay = positions[i0 * 3 + 1], az = positions[i0 * 3 + 2];
    const bx = positions[i1 * 3], by = positions[i1 * 3 + 1], bz = positions[i1 * 3 + 2];
    const cx = positions[i2 * 3], cy = positions[i2 * 3 + 1], cz = positions[i2 * 3 + 2];
    const ux = bx - ax, uy = by - ay, uz = bz - az;
    const vvx = cx - ax, vvy = cy - ay, vvz = cz - az;
    const crx = uy * vvz - uz * vvy;
    const cry = uz * vvx - ux * vvz;
    const crz = ux * vvy - uy * vvx;
    const area2 = crx * crx + cry * cry + crz * crz;
    if (0.25 * area2 < 1e-20) degenerateCount++;

    // 3 arêtes (i0→i1), (i1→i2), (i2→i0)
    for (let e = 0; e < 3; e++) {
      const u = indices[f * 3 + e];
      const v = indices[f * 3 + ((e + 1) % 3)];
      const dKey = u * ENC + v;
      directed.set(dKey, (directed.get(dKey) ?? 0) + 1);
      const lo = u < v ? u : v;
      const hi = u < v ? v : u;
      const uKey = lo * ENC + hi;
      undirected.set(uKey, (undirected.get(uKey) ?? 0) + 1);
    }
  }

  const E = undirected.size;
  const euler = V - E + F;

  let watertight = true;
  for (const c of undirected.values()) {
    if (c !== 2) { watertight = false; break; }
  }

  let orientable = true;
  for (const c of directed.values()) {
    if (c !== 1) { orientable = false; break; }
  }

  // Manifold : link de chaque sommet = cycle unique.
  // Pour chaque sommet v, les triangles incidents définissent des arêtes
  // opposées (a,b). Ces arêtes doivent former un cycle unique.
  let manifold = true;
  for (let vi = 0; vi < V && manifold; vi++) {
    const start = offset[vi];
    const end = offset[vi + 1];
    if (start === end) continue;

    const linkAdj = new Map<number, number[]>();
    let bad = false;
    for (let k = start; k < end; k++) {
      const f = incident[k];
      const i0 = indices[f * 3];
      const i1 = indices[f * 3 + 1];
      const i2 = indices[f * 3 + 2];
      let a: number, b: number;
      if (i0 === vi) { a = i1; b = i2; }
      else if (i1 === vi) { a = i2; b = i0; }
      else { a = i0; b = i1; }

      let la = linkAdj.get(a);
      if (!la) { la = []; linkAdj.set(a, la); }
      la.push(b);
      let lb = linkAdj.get(b);
      if (!lb) { lb = []; linkAdj.set(b, lb); }
      lb.push(a);
    }
    // degré 2 partout
    for (const list of linkAdj.values()) {
      if (list.length !== 2) { bad = true; break; }
    }
    if (bad) { manifold = false; break; }

    // cycle unique : parcourir depuis un point, compter les pas
    let first = -1;
    for (const k of linkAdj.keys()) { first = k; break; }
    if (first < 0) continue;
    let prev = -1;
    let curr = first;
    let steps = 0;
    const n = linkAdj.size;
    while (steps < n) {
      const neigh = linkAdj.get(curr)!;
      const next = neigh[0] !== prev ? neigh[0] : neigh[1];
      prev = curr;
      curr = next;
      steps++;
      if (curr === first) break;
    }
    if (curr !== first || steps !== n) { manifold = false; break; }
  }

  return { V, E, F, euler, watertight, manifold, orientable, degenerateCount };
}
