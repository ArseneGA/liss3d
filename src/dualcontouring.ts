import { solveQEF } from "./qef";
import type { Grid, Mesh } from "./types";

const MAX_EDGES_PER_CELL = 12;

/**
 * 12 arêtes d'une cellule : paires (corner_a, corner_b) parmi les 8 coins.
 *
 * Numérotation des coins dans l'ordre (i,j,k), (i+1,j,k), (i+1,j+1,k),
 * (i,j+1,k), (i,j,k+1), (i+1,j,k+1), (i+1,j+1,k+1), (i,j+1,k+1).
 */
const CELL_EDGES: ReadonlyArray<readonly [number, number]> = [
  [0, 1], [1, 2], [3, 2], [0, 3],
  [4, 5], [5, 6], [7, 6], [4, 7],
  [0, 4], [1, 5], [2, 6], [3, 7],
];

const CORNER_OFFSETS: ReadonlyArray<readonly [number, number, number]> = [
  [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
  [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
];

/**
 * Dual Contouring : extrait un mesh de l'isosurface f = 0 à partir d'un champ
 * scalaire échantillonné sur une grille régulière.
 *
 * Référence : Ju, Losasso, Schaefer, Warren, "Dual Contouring of Hermite Data",
 * SIGGRAPH 2002.
 */
export function dualContour(field: Float32Array, grid: Grid): Mesh {
  const { nx, ny, nz, h, originX, originY, originZ } = grid;
  const nxny = nx * ny;
  const cnx = nx - 1, cny = ny - 1, cnz = nz - 1;
  const cellCount = cnx * cny * cnz;

  const cellVertexIdx = new Int32Array(cellCount);
  cellVertexIdx.fill(-1);

  const positions: number[] = [];
  const indices: number[] = [];

  // Buffers réutilisés entre cellules pour les arêtes traversées.
  const edgePts = new Float32Array(MAX_EDGES_PER_CELL * 3);
  const edgeNormals = new Float32Array(MAX_EDGES_PER_CELL * 3);

  // ∇f au point de grille (i,j,k) via différences centrées (clamp aux bords).
  const grad = (i: number, j: number, k: number, out: [number, number, number]): void => {
    const ip = Math.min(nx - 1, i + 1);
    const im = Math.max(0, i - 1);
    const jp = Math.min(ny - 1, j + 1);
    const jm = Math.max(0, j - 1);
    const kp = Math.min(nz - 1, k + 1);
    const km = Math.max(0, k - 1);
    const invDx = 1 / ((ip - im) * h);
    const invDy = 1 / ((jp - jm) * h);
    const invDz = 1 / ((kp - km) * h);
    out[0] = (field[ip + j * nx + k * nxny] - field[im + j * nx + k * nxny]) * invDx;
    out[1] = (field[i + jp * nx + k * nxny] - field[i + jm * nx + k * nxny]) * invDy;
    out[2] = (field[i + j * nx + kp * nxny] - field[i + j * nx + km * nxny]) * invDz;
  };

  const cornerIdx = (ci: number, cj: number, ck: number, corner: number): number => {
    const o = CORNER_OFFSETS[corner];
    return (ci + o[0]) + (cj + o[1]) * nx + (ck + o[2]) * nxny;
  };

  const gradA: [number, number, number] = [0, 0, 0];
  const gradB: [number, number, number] = [0, 0, 0];
  const v = new Float32Array(8);

  // --- Pass 1 : pour chaque cellule avec changement de signe, calculer le sommet QEF.
  for (let ck = 0; ck < cnz; ck++) {
    for (let cj = 0; cj < cny; cj++) {
      for (let ci = 0; ci < cnx; ci++) {
        // Lire les 8 valeurs, construire le masque de signes.
        let mask = 0;
        for (let c = 0; c < 8; c++) {
          const val = field[cornerIdx(ci, cj, ck, c)];
          v[c] = val;
          if (val < 0) mask |= 1 << c;
        }
        if (mask === 0 || mask === 255) continue;

        // Positions des 8 coins en world-space.
        const cox = originX + ci * h;
        const coy = originY + cj * h;
        const coz = originZ + ck * h;

        let edgeCount = 0;
        for (let e = 0; e < 12; e++) {
          const a = CELL_EDGES[e][0];
          const b = CELL_EDGES[e][1];
          const va = v[a];
          const vb = v[b];
          if ((va < 0) === (vb < 0)) continue;
          const t = va / (va - vb);
          const oa = CORNER_OFFSETS[a];
          const ob = CORNER_OFFSETS[b];
          const pax = cox + oa[0] * h;
          const pay = coy + oa[1] * h;
          const paz = coz + oa[2] * h;
          const pbx = cox + ob[0] * h;
          const pby = coy + ob[1] * h;
          const pbz = coz + ob[2] * h;
          const epx = pax + t * (pbx - pax);
          const epy = pay + t * (pby - pay);
          const epz = paz + t * (pbz - paz);

          grad(ci + oa[0], cj + oa[1], ck + oa[2], gradA);
          grad(ci + ob[0], cj + ob[1], ck + ob[2], gradB);
          let nx_ = gradA[0] + t * (gradB[0] - gradA[0]);
          let ny_ = gradA[1] + t * (gradB[1] - gradA[1]);
          let nz_ = gradA[2] + t * (gradB[2] - gradA[2]);
          const nlen = Math.hypot(nx_, ny_, nz_);
          if (nlen > 1e-12) {
            nx_ /= nlen; ny_ /= nlen; nz_ /= nlen;
          } else {
            // fallback : normale alignée sur l'arête
            nx_ = ob[0] - oa[0];
            ny_ = ob[1] - oa[1];
            nz_ = ob[2] - oa[2];
            const nl = Math.hypot(nx_, ny_, nz_);
            nx_ /= nl; ny_ /= nl; nz_ /= nl;
          }

          edgePts[edgeCount * 3] = epx;
          edgePts[edgeCount * 3 + 1] = epy;
          edgePts[edgeCount * 3 + 2] = epz;
          edgeNormals[edgeCount * 3] = nx_;
          edgeNormals[edgeCount * 3 + 1] = ny_;
          edgeNormals[edgeCount * 3 + 2] = nz_;
          edgeCount++;
        }

        if (edgeCount === 0) continue;

        // cBias = moyenne des points d'intersection.
        let cx = 0, cy = 0, cz = 0;
        for (let i = 0; i < edgeCount; i++) {
          cx += edgePts[i * 3];
          cy += edgePts[i * 3 + 1];
          cz += edgePts[i * 3 + 2];
        }
        cx /= edgeCount; cy /= edgeCount; cz /= edgeCount;

        const [vx, vy, vz] = solveQEF(edgePts, edgeNormals, edgeCount, cx, cy, cz);

        // Clamp dans la cellule étendue d'un demi-voxel.
        const minX = cox - 0.5 * h, maxX = cox + 1.5 * h;
        const minY = coy - 0.5 * h, maxY = coy + 1.5 * h;
        const minZ = coz - 0.5 * h, maxZ = coz + 1.5 * h;
        const fx = vx < minX ? minX : vx > maxX ? maxX : vx;
        const fy = vy < minY ? minY : vy > maxY ? maxY : vy;
        const fz = vz < minZ ? minZ : vz > maxZ ? maxZ : vz;

        const vertIdx = positions.length / 3;
        positions.push(fx, fy, fz);
        cellVertexIdx[ci + cj * cnx + ck * cnx * cny] = vertIdx;
      }
    }
  }

  // --- Pass 2 : pour chaque arête de grille intérieure avec changement de signe,
  // émettre un quad connectant les sommets des 4 cellules adjacentes.
  const cellIdx = (ci: number, cj: number, ck: number): number =>
    ci + cj * cnx + ck * cnx * cny;

  const emitQuad = (
    a: number, b: number, c: number, d: number,
    flip: boolean,
  ): void => {
    // a, b, c, d sont dans l'ordre CCW vu depuis l'extérieur (gradient +)
    const q = flip ? [a, d, c, b] : [a, b, c, d];
    // Triangulation par diagonale la plus courte
    const pa = q[0] * 3, pb = q[1] * 3, pc = q[2] * 3, pd = q[3] * 3;
    const dAC = sqDist(positions, pa, pc);
    const dBD = sqDist(positions, pb, pd);
    if (dAC < dBD) {
      indices.push(q[0], q[1], q[2], q[0], q[2], q[3]);
    } else {
      indices.push(q[1], q[2], q[3], q[1], q[3], q[0]);
    }
  };

  // Arêtes +X : grille (i, j, k) → (i+1, j, k), partagées par 4 cellules
  for (let k = 1; k < nz - 1; k++) {
    for (let j = 1; j < ny - 1; j++) {
      for (let i = 0; i < nx - 1; i++) {
        const f0 = field[i + j * nx + k * nxny];
        const f1 = field[(i + 1) + j * nx + k * nxny];
        if ((f0 < 0) === (f1 < 0)) continue;
        const v00 = cellVertexIdx[cellIdx(i, j - 1, k - 1)];
        const v10 = cellVertexIdx[cellIdx(i, j, k - 1)];
        const v11 = cellVertexIdx[cellIdx(i, j, k)];
        const v01 = cellVertexIdx[cellIdx(i, j - 1, k)];
        if (v00 < 0 || v10 < 0 || v11 < 0 || v01 < 0) continue;
        // CCW vu depuis +X : (+Y,−Z) → (+Y,+Z) → (−Y,+Z) → (−Y,−Z)
        //                  = cell(i,j,k−1), cell(i,j,k), cell(i,j−1,k), cell(i,j−1,k−1)
        emitQuad(v10, v11, v01, v00, f0 >= 0);
      }
    }
  }
  // Arêtes +Y
  for (let k = 1; k < nz - 1; k++) {
    for (let j = 0; j < ny - 1; j++) {
      for (let i = 1; i < nx - 1; i++) {
        const f0 = field[i + j * nx + k * nxny];
        const f1 = field[i + (j + 1) * nx + k * nxny];
        if ((f0 < 0) === (f1 < 0)) continue;
        const v00 = cellVertexIdx[cellIdx(i - 1, j, k - 1)];
        const v10 = cellVertexIdx[cellIdx(i, j, k - 1)];
        const v11 = cellVertexIdx[cellIdx(i, j, k)];
        const v01 = cellVertexIdx[cellIdx(i - 1, j, k)];
        if (v00 < 0 || v10 < 0 || v11 < 0 || v01 < 0) continue;
        // CCW vu depuis +Y : (−X,+Z) → (+X,+Z) → (+X,−Z) → (−X,−Z)
        //                  = cell(i−1,j,k), cell(i,j,k), cell(i,j,k−1), cell(i−1,j,k−1)
        emitQuad(v01, v11, v10, v00, f0 >= 0);
      }
    }
  }
  // Arêtes +Z
  for (let k = 0; k < nz - 1; k++) {
    for (let j = 1; j < ny - 1; j++) {
      for (let i = 1; i < nx - 1; i++) {
        const f0 = field[i + j * nx + k * nxny];
        const f1 = field[i + j * nx + (k + 1) * nxny];
        if ((f0 < 0) === (f1 < 0)) continue;
        const v00 = cellVertexIdx[cellIdx(i - 1, j - 1, k)];
        const v10 = cellVertexIdx[cellIdx(i, j - 1, k)];
        const v11 = cellVertexIdx[cellIdx(i, j, k)];
        const v01 = cellVertexIdx[cellIdx(i - 1, j, k)];
        if (v00 < 0 || v10 < 0 || v11 < 0 || v01 < 0) continue;
        // CCW vu depuis +Z : (−X,−Y) → (+X,−Y) → (+X,+Y) → (−X,+Y)
        //                  = cell(i−1,j−1,k), cell(i,j−1,k), cell(i,j,k), cell(i−1,j,k)
        emitQuad(v00, v10, v11, v01, f0 >= 0);
      }
    }
  }

  return {
    positions: new Float32Array(positions),
    indices: new Uint32Array(indices),
  };
}

function sqDist(positions: number[], ia: number, ib: number): number {
  const dx = positions[ia] - positions[ib];
  const dy = positions[ia + 1] - positions[ib + 1];
  const dz = positions[ia + 2] - positions[ib + 2];
  return dx * dx + dy * dy + dz * dz;
}
