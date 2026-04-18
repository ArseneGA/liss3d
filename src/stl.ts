/**
 * STL binaire :
 *   80 octets header + uint32 nTriangles
 *   puis, par triangle :
 *     float32[3] normale + 3 × float32[3] sommet + uint16 attribute
 *
 * Le mesh est recentré et mis à l'échelle pour que sa dimension max vaille
 * `scaleMm`. Normales recalculées par produit vectoriel (b−a) × (c−a).
 */
export function exportSTLBinary(
  positions: Float32Array,
  indices: Uint32Array,
  scaleMm: number,
): ArrayBuffer {
  const F = indices.length / 3;
  const buffer = new ArrayBuffer(80 + 4 + F * 50);
  const view = new DataView(buffer);

  // bbox pour recentrage + scale uniforme
  let minX = Infinity, minY = Infinity, minZ = Infinity;
  let maxX = -Infinity, maxY = -Infinity, maxZ = -Infinity;
  for (let i = 0; i < positions.length; i += 3) {
    const x = positions[i], y = positions[i + 1], z = positions[i + 2];
    if (x < minX) minX = x; if (x > maxX) maxX = x;
    if (y < minY) minY = y; if (y > maxY) maxY = y;
    if (z < minZ) minZ = z; if (z > maxZ) maxZ = z;
  }
  const cx = (minX + maxX) / 2;
  const cy = (minY + maxY) / 2;
  const cz = (minZ + maxZ) / 2;
  const size = Math.max(maxX - minX, maxY - minY, maxZ - minZ);
  const scale = size > 0 ? scaleMm / size : 1;

  view.setUint32(80, F, true);
  let off = 84;

  for (let f = 0; f < F; f++) {
    const i0 = indices[f * 3];
    const i1 = indices[f * 3 + 1];
    const i2 = indices[f * 3 + 2];
    const ax = (positions[i0 * 3] - cx) * scale;
    const ay = (positions[i0 * 3 + 1] - cy) * scale;
    const az = (positions[i0 * 3 + 2] - cz) * scale;
    const bx = (positions[i1 * 3] - cx) * scale;
    const by = (positions[i1 * 3 + 1] - cy) * scale;
    const bz = (positions[i1 * 3 + 2] - cz) * scale;
    const qx = (positions[i2 * 3] - cx) * scale;
    const qy = (positions[i2 * 3 + 1] - cy) * scale;
    const qz = (positions[i2 * 3 + 2] - cz) * scale;

    const ux = bx - ax, uy = by - ay, uz = bz - az;
    const vx = qx - ax, vy = qy - ay, vz = qz - az;
    let nx = uy * vz - uz * vy;
    let ny = uz * vx - ux * vz;
    let nz = ux * vy - uy * vx;
    const nl = Math.sqrt(nx * nx + ny * ny + nz * nz);
    if (nl > 0) { nx /= nl; ny /= nl; nz /= nl; }

    view.setFloat32(off, nx, true); off += 4;
    view.setFloat32(off, ny, true); off += 4;
    view.setFloat32(off, nz, true); off += 4;
    view.setFloat32(off, ax, true); off += 4;
    view.setFloat32(off, ay, true); off += 4;
    view.setFloat32(off, az, true); off += 4;
    view.setFloat32(off, bx, true); off += 4;
    view.setFloat32(off, by, true); off += 4;
    view.setFloat32(off, bz, true); off += 4;
    view.setFloat32(off, qx, true); off += 4;
    view.setFloat32(off, qy, true); off += 4;
    view.setFloat32(off, qz, true); off += 4;
    view.setUint16(off, 0, true); off += 2;
  }

  return buffer;
}

export function downloadSTL(buffer: ArrayBuffer, filename: string): void {
  const blob = new Blob([buffer], { type: "application/octet-stream" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}
