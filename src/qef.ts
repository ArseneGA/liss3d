/**
 * Diagonalisation Jacobi 3×3 symétrique.
 * A = U · diag(d) · Uᵀ, U rotation (colonnes = vecteurs propres).
 *
 * Convention : `u` stocké row-major ⇒ U[i][j] = u[i·3 + j].
 */
export function jacobi3Sym(
  a00: number, a11: number, a22: number,
  a01: number, a02: number, a12: number,
): {
  d: [number, number, number];
  u: [number, number, number, number, number, number, number, number, number];
} {
  let U00 = 1, U01 = 0, U02 = 0;
  let U10 = 0, U11 = 1, U12 = 0;
  let U20 = 0, U21 = 0, U22 = 1;
  let d0 = a00, d1 = a11, d2 = a22;
  let o01 = a01, o02 = a02, o12 = a12;

  for (let sweep = 0; sweep < 8; sweep++) {
    const magOff = Math.abs(o01) + Math.abs(o02) + Math.abs(o12);
    if (magOff < 1e-14) break;

    // rotation (0,1)
    if (Math.abs(o01) > 1e-20) {
      const theta = (d1 - d0) / (2 * o01);
      const sgn = theta >= 0 ? 1 : -1;
      const t = sgn / (Math.abs(theta) + Math.sqrt(1 + theta * theta));
      const c = 1 / Math.sqrt(1 + t * t);
      const s = t * c;
      const nd0 = d0 - t * o01;
      const nd1 = d1 + t * o01;
      const no02 = c * o02 - s * o12;
      const no12 = s * o02 + c * o12;
      d0 = nd0; d1 = nd1; o02 = no02; o12 = no12; o01 = 0;
      // colonnes 0 et 1 de U
      const a = c * U00 - s * U01, b = s * U00 + c * U01;
      const aa = c * U10 - s * U11, bb = s * U10 + c * U11;
      const aaa = c * U20 - s * U21, bbb = s * U20 + c * U21;
      U00 = a; U01 = b; U10 = aa; U11 = bb; U20 = aaa; U21 = bbb;
    }

    // rotation (0,2)
    if (Math.abs(o02) > 1e-20) {
      const theta = (d2 - d0) / (2 * o02);
      const sgn = theta >= 0 ? 1 : -1;
      const t = sgn / (Math.abs(theta) + Math.sqrt(1 + theta * theta));
      const c = 1 / Math.sqrt(1 + t * t);
      const s = t * c;
      const nd0 = d0 - t * o02;
      const nd2 = d2 + t * o02;
      const no01 = c * o01 - s * o12;
      const no12 = s * o01 + c * o12;
      d0 = nd0; d2 = nd2; o01 = no01; o12 = no12; o02 = 0;
      const a = c * U00 - s * U02, b = s * U00 + c * U02;
      const aa = c * U10 - s * U12, bb = s * U10 + c * U12;
      const aaa = c * U20 - s * U22, bbb = s * U20 + c * U22;
      U00 = a; U02 = b; U10 = aa; U12 = bb; U20 = aaa; U22 = bbb;
    }

    // rotation (1,2)
    if (Math.abs(o12) > 1e-20) {
      const theta = (d2 - d1) / (2 * o12);
      const sgn = theta >= 0 ? 1 : -1;
      const t = sgn / (Math.abs(theta) + Math.sqrt(1 + theta * theta));
      const c = 1 / Math.sqrt(1 + t * t);
      const s = t * c;
      const nd1 = d1 - t * o12;
      const nd2 = d2 + t * o12;
      const no01 = c * o01 - s * o02;
      const no02 = s * o01 + c * o02;
      d1 = nd1; d2 = nd2; o01 = no01; o02 = no02; o12 = 0;
      const a = c * U01 - s * U02, b = s * U01 + c * U02;
      const aa = c * U11 - s * U12, bb = s * U11 + c * U12;
      const aaa = c * U21 - s * U22, bbb = s * U21 + c * U22;
      U01 = a; U02 = b; U11 = aa; U12 = bb; U21 = aaa; U22 = bbb;
    }
  }

  return {
    d: [d0, d1, d2],
    u: [U00, U01, U02, U10, U11, U12, U20, U21, U22],
  };
}

/**
 * QEF : v* = argmin Σ (n_e · (v − p_e))²
 *
 * A = Σ n n^T    (3×3 symétrique PSD, possiblement rang 1 ou 2)
 * b = Σ (n·p) n  (vecteur 3)
 *
 * Solution biaisée vers `cBias` (centre cellule) pour stabilité dans les modes
 * dégénérés : v = cBias + pinv(A) · (b − A · cBias),
 * avec pinv via SVD tronqué (σ_i < σ_max · 0.1 ⇒ 0).
 */
export function solveQEF(
  points: Float32Array,
  normals: Float32Array,
  count: number,
  cBiasX: number,
  cBiasY: number,
  cBiasZ: number,
): [number, number, number] {
  let A00 = 0, A11 = 0, A22 = 0, A01 = 0, A02 = 0, A12 = 0;
  let bx = 0, by = 0, bz = 0;
  for (let i = 0; i < count; i++) {
    const nx = normals[i * 3];
    const ny = normals[i * 3 + 1];
    const nz = normals[i * 3 + 2];
    const px = points[i * 3];
    const py = points[i * 3 + 1];
    const pz = points[i * 3 + 2];
    A00 += nx * nx; A11 += ny * ny; A22 += nz * nz;
    A01 += nx * ny; A02 += nx * nz; A12 += ny * nz;
    const dp = nx * px + ny * py + nz * pz;
    bx += dp * nx; by += dp * ny; bz += dp * nz;
  }

  // rhs = b − A · cBias
  const r0 = bx - (A00 * cBiasX + A01 * cBiasY + A02 * cBiasZ);
  const r1 = by - (A01 * cBiasX + A11 * cBiasY + A12 * cBiasZ);
  const r2 = bz - (A02 * cBiasX + A12 * cBiasY + A22 * cBiasZ);

  const { d, u } = jacobi3Sym(A00, A11, A22, A01, A02, A12);
  const sigMax = Math.max(Math.abs(d[0]), Math.abs(d[1]), Math.abs(d[2]));
  const thresh = sigMax * 0.1;
  const invd0 = Math.abs(d[0]) > thresh ? 1 / d[0] : 0;
  const invd1 = Math.abs(d[1]) > thresh ? 1 / d[1] : 0;
  const invd2 = Math.abs(d[2]) > thresh ? 1 / d[2] : 0;

  // y = Uᵀ · rhs
  const y0 = u[0] * r0 + u[3] * r1 + u[6] * r2;
  const y1 = u[1] * r0 + u[4] * r1 + u[7] * r2;
  const y2 = u[2] * r0 + u[5] * r1 + u[8] * r2;
  // z = diag(invd) · y
  const z0 = invd0 * y0, z1 = invd1 * y1, z2 = invd2 * y2;
  // dv = U · z
  const dvx = u[0] * z0 + u[1] * z1 + u[2] * z2;
  const dvy = u[3] * z0 + u[4] * z1 + u[5] * z2;
  const dvz = u[6] * z0 + u[7] * z1 + u[8] * z2;
  return [cBiasX + dvx, cBiasY + dvy, cBiasZ + dvz];
}
