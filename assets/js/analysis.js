// Numerical helpers: LOSVD moments, Gauss-Hermite fits, statistics, spatial
// grids, interpolation, contours and MGE surface-brightness models.

const DELAUNATOR_URL = "https://cdn.jsdelivr.net/npm/delaunator@5.0.1/+esm";

// ---------------------------------------------------------------- geometry

export function rotatePoint(x, y, cos, sin) {
  return [x * cos - y * sin, x * sin + y * cos];
}

export function rotationTrig(deg) {
  const t = (Number(deg) || 0) * Math.PI / 180;
  return [Math.cos(t), Math.sin(t)];
}

// Median nearest-neighbour distance (sampled for large inputs).
export function medianSpacing(xs, ys) {
  const idx = [];
  for (let i = 0; i < xs.length; i += 1) {
    if (Number.isFinite(xs[i]) && Number.isFinite(ys[i])) idx.push(i);
  }
  if (idx.length < 2) return 0.1;
  const step = Math.max(1, Math.floor(idx.length / 400));
  const d = [];
  for (let a = 0; a < idx.length; a += step) {
    const i = idx[a];
    let best = Infinity;
    for (const j of idx) {
      if (j === i) continue;
      const d2 = (xs[i] - xs[j]) ** 2 + (ys[i] - ys[j]) ** 2;
      if (d2 > 0 && d2 < best) best = d2;
    }
    if (Number.isFinite(best)) d.push(Math.sqrt(best));
  }
  d.sort((p, q) => p - q);
  return d.length ? d[Math.floor(d.length / 2)] : 0.1;
}

// For every bin of `a`, the index of the nearest bin of `b` within `tol` (or -1).
export function matchBins(a, b, tol) {
  const out = new Int32Array(a.x.length).fill(-1);
  if (a === b) {
    for (let i = 0; i < out.length; i += 1) out[i] = i;
    return out;
  }
  const tol2 = tol * tol;
  for (let i = 0; i < a.x.length; i += 1) {
    const xi = a.x[i];
    const yi = a.y[i];
    if (!Number.isFinite(xi) || !Number.isFinite(yi)) continue;
    // Fast path: identical binning schemes share indices.
    if (i < b.x.length && (b.x[i] - xi) ** 2 + (b.y[i] - yi) ** 2 <= tol2) {
      out[i] = i;
      continue;
    }
    let best = -1;
    let bestD2 = tol2;
    for (let j = 0; j < b.x.length; j += 1) {
      const d2 = (b.x[j] - xi) ** 2 + (b.y[j] - yi) ** 2;
      if (d2 <= bestD2) {
        bestD2 = d2;
        best = j;
      }
    }
    out[i] = best;
  }
  return out;
}

// ---------------------------------------------------------------- LOSVDs

// Sort a velocity axis ascending; returns the permutation and sorted axis.
export function sortVelocityAxis(xvel) {
  const order = Array.from(xvel, (_, i) => i).sort((i, j) => xvel[i] - xvel[j]);
  return {order, xvel: Float64Array.from(order, (i) => xvel[i])};
}

export function trapz(y, x, offset = 0, n = x.length) {
  let s = 0;
  for (let i = 1; i < n; i += 1) {
    const a = y[offset + i - 1];
    const b = y[offset + i];
    if (Number.isFinite(a) && Number.isFinite(b)) s += 0.5 * (a + b) * (x[i] - x[i - 1]);
  }
  return s;
}

// Mean, dispersion, standardized skewness and excess kurtosis of a
// tabulated LOSVD (same definition as the original dashboard builder).
export function losvdMoments(xvel, losvd, offset = 0) {
  const n = xvel.length;
  let area = 0;
  let m1 = 0;
  let count = 0;
  for (let i = 0; i < n; i += 1) {
    const y = losvd[offset + i];
    if (!Number.isFinite(y) || !Number.isFinite(xvel[i])) continue;
    area += y;
    m1 += xvel[i] * y;
    count += 1;
  }
  const nan = {vel: NaN, sigma: NaN, skew: NaN, kurt: NaN};
  if (count < 3 || !(area > 0)) return nan;
  const mean = m1 / area;
  let m2 = 0;
  for (let i = 0; i < n; i += 1) {
    const y = losvd[offset + i];
    if (Number.isFinite(y)) m2 += (xvel[i] - mean) ** 2 * y;
  }
  const variance = m2 / area;
  if (!(variance > 0)) return {...nan, vel: mean};
  const sigma = Math.sqrt(variance);
  let m3 = 0;
  let m4 = 0;
  for (let i = 0; i < n; i += 1) {
    const y = losvd[offset + i];
    if (!Number.isFinite(y)) continue;
    const w = (xvel[i] - mean) / sigma;
    m3 += w * w * w * y;
    m4 += w * w * w * w * y;
  }
  return {vel: mean, sigma, skew: m3 / area, kurt: m4 / area - 3};
}

const SQRT2PI = Math.sqrt(2 * Math.PI);
const SQRT2 = Math.SQRT2;
const SQRT6 = Math.sqrt(6);
const SQRT24 = Math.sqrt(24);

// Gauss-Hermite series (van der Marel & Franx 1993), truncated at h4.
export function ghProfile(v, amp, vel, sigma, h3, h4) {
  const w = (v - vel) / sigma;
  const w2 = w * w;
  const H3 = (2 * SQRT2 * w2 * w - 3 * SQRT2 * w) / SQRT6;
  const H4 = (4 * w2 * w2 - 12 * w2 + 3) / SQRT24;
  return amp / (SQRT2PI * sigma) * Math.exp(-0.5 * w2) * (1 + h3 * H3 + h4 * H4);
}

function solveLinear(A, b) {
  const n = b.length;
  const M = A.map((row, i) => [...row, b[i]]);
  for (let c = 0; c < n; c += 1) {
    let piv = c;
    for (let r = c + 1; r < n; r += 1) if (Math.abs(M[r][c]) > Math.abs(M[piv][c])) piv = r;
    if (Math.abs(M[piv][c]) < 1e-300) return null;
    [M[c], M[piv]] = [M[piv], M[c]];
    for (let r = 0; r < n; r += 1) {
      if (r === c) continue;
      const f = M[r][c] / M[c][c];
      for (let k = c; k <= n; k += 1) M[r][k] -= f * M[c][k];
    }
  }
  return M.map((row, i) => row[n] / row[i]);
}

// Levenberg-Marquardt least-squares fit of a Gauss-Hermite profile to a
// tabulated LOSVD. Returns {amp, vel, sigma, h3, h4, ok}.
export function fitGaussHermite(xvel, losvd, offset = 0) {
  const n = xvel.length;
  const xs = [];
  const ys = [];
  for (let i = 0; i < n; i += 1) {
    const y = losvd[offset + i];
    if (Number.isFinite(y)) {
      xs.push(xvel[i]);
      ys.push(y);
    }
  }
  const bad = {amp: NaN, vel: NaN, sigma: NaN, h3: NaN, h4: NaN, ok: false};
  if (xs.length < 6) return bad;
  const m = losvdMoments(Float64Array.from(xs), Float64Array.from(ys));
  if (!Number.isFinite(m.sigma)) return bad;
  const area = trapz(ys, xs);
  let p = [area > 0 ? area : 1, m.vel, m.sigma, 0, 0];

  const residuals = (q) => xs.map((x, i) => ghProfile(x, q[0], q[1], q[2], q[3], q[4]) - ys[i]);
  const chi2Of = (r) => r.reduce((s, v) => s + v * v, 0);
  let r = residuals(p);
  let chi2 = chi2Of(r);
  let lambda = 1e-3;
  const dv = Math.abs(xs[1] - xs[0]) || 1;
  const steps = [1e-6, 1e-4 * dv, 1e-4 * dv, 1e-6, 1e-6];

  for (let iter = 0; iter < 200; iter += 1) {
    const J = xs.map(() => new Array(5).fill(0));
    for (let k = 0; k < 5; k += 1) {
      const h = steps[k] * Math.max(1, k === 0 ? Math.abs(p[0]) : 1);
      const qp = p.slice();
      qp[k] += h;
      const rp = residuals(qp);
      for (let i = 0; i < xs.length; i += 1) J[i][k] = (rp[i] - r[i]) / h;
    }
    const JTJ = Array.from({length: 5}, () => new Array(5).fill(0));
    const JTr = new Array(5).fill(0);
    for (let i = 0; i < xs.length; i += 1) {
      for (let a = 0; a < 5; a += 1) {
        JTr[a] += J[i][a] * r[i];
        for (let b = 0; b < 5; b += 1) JTJ[a][b] += J[i][a] * J[i][b];
      }
    }
    let improved = false;
    for (let tries = 0; tries < 12; tries += 1) {
      const A = JTJ.map((row, a) => row.map((v, b) => (a === b ? v * (1 + lambda) + 1e-30 : v)));
      const delta = solveLinear(A, JTr.map((v) => -v));
      if (!delta) {
        lambda *= 10;
        continue;
      }
      const q = p.map((v, k) => v + delta[k]);
      if (!(q[2] > 0)) {
        lambda *= 10;
        continue;
      }
      const rq = residuals(q);
      const c2 = chi2Of(rq);
      if (c2 < chi2) {
        const rel = (chi2 - c2) / Math.max(chi2, 1e-300);
        p = q;
        r = rq;
        chi2 = c2;
        lambda = Math.max(lambda / 10, 1e-12);
        improved = true;
        if (rel < 1e-12) iter = Infinity;
        break;
      }
      lambda *= 10;
    }
    if (!improved) break;
  }
  return {amp: p[0], vel: p[1], sigma: Math.abs(p[2]), h3: p[3], h4: p[4], ok: true};
}

// ---------------------------------------------------------------- statistics

export function finiteSorted(values, mask = null) {
  const out = [];
  for (let i = 0; i < values.length; i += 1) {
    if (mask && mask[i]) continue;
    if (Number.isFinite(values[i])) out.push(values[i]);
  }
  return Float64Array.from(out).sort();
}

export function percentileSorted(sorted, p) {
  if (!sorted.length) return NaN;
  const pos = (sorted.length - 1) * p / 100;
  const lo = Math.floor(pos);
  const hi = Math.ceil(pos);
  return sorted[lo] + (sorted[hi] - sorted[lo]) * (pos - lo);
}

export function limitsFor(values, mode = "robust") {
  const s = finiteSorted(values);
  if (!s.length) return [0, 1];
  let lo;
  let hi;
  if (mode === "minmax") {
    lo = s[0];
    hi = s[s.length - 1];
  } else if (mode === "symmetric") {
    const abs = Float64Array.from(s, Math.abs).sort();
    const a = percentileSorted(abs, 98) || abs[abs.length - 1];
    lo = -a;
    hi = a;
  } else {
    lo = percentileSorted(s, 2);
    hi = percentileSorted(s, 98);
  }
  if (!(hi > lo)) {
    lo = s[0];
    hi = s[s.length - 1];
  }
  if (!(hi > lo)) {
    const c = Number.isFinite(lo) ? lo : 0;
    const d = Math.abs(c) * 0.05 || 1;
    lo = c - d;
    hi = c + d;
  }
  return [lo, hi];
}

// ---------------------------------------------------------------- grids

function minPositiveStep(values) {
  const u = Array.from(new Set(Array.from(values, (v) => Math.round(v * 1e5)))).sort((a, b) => a - b);
  let best = Infinity;
  for (let i = 1; i < u.length; i += 1) {
    const d = u[i] - u[i - 1];
    if (d > 0 && d < best) best = d;
  }
  return Number.isFinite(best) ? best / 1e5 : NaN;
}

// Places scattered spaxel centres on a regular grid. Returns null when the
// spaxels are not on a regular lattice.
export function buildSpaxelGrid(xs, ys, binIds, flux) {
  const n = xs.length;
  if (!n) return null;
  const dx = minPositiveStep(xs);
  const dy = minPositiveStep(ys);
  if (!(dx > 0) || !(dy > 0)) return null;
  let xmin = Infinity;
  let xmax = -Infinity;
  let ymin = Infinity;
  let ymax = -Infinity;
  for (let i = 0; i < n; i += 1) {
    xmin = Math.min(xmin, xs[i]);
    xmax = Math.max(xmax, xs[i]);
    ymin = Math.min(ymin, ys[i]);
    ymax = Math.max(ymax, ys[i]);
  }
  const nx = Math.round((xmax - xmin) / dx) + 1;
  const ny = Math.round((ymax - ymin) / dy) + 1;
  if (nx * ny > 6e6) return null;
  const bin = new Int32Array(nx * ny).fill(-1);
  const fl = new Float64Array(nx * ny).fill(NaN);
  for (let i = 0; i < n; i += 1) {
    const ix = Math.round((xs[i] - xmin) / dx);
    const iy = Math.round((ys[i] - ymin) / dy);
    const k = iy * nx + ix;
    bin[k] = binIds ? Number(binIds[i]) : -1;
    if (flux) fl[k] = flux[i];
  }
  return {nx, ny, x0: xmin, y0: ymin, dx, dy, bin, flux: flux ? fl : null};
}

export function gridCorners(grid) {
  const x0 = grid.x0 - grid.dx / 2;
  const y0 = grid.y0 - grid.dy / 2;
  const x1 = grid.x0 + (grid.nx - 0.5) * grid.dx;
  const y1 = grid.y0 + (grid.ny - 0.5) * grid.dy;
  return [[x0, y0], [x1, y0], [x1, y1], [x0, y1]];
}

let delaunatorPromise = null;
function loadDelaunator() {
  if (!delaunatorPromise) delaunatorPromise = import(DELAUNATOR_URL).then((m) => m.default);
  return delaunatorPromise;
}

// Linear (Delaunay barycentric) interpolation of bin values onto a grid.
export async function interpolateBins(xs, ys, values, spacing) {
  const Delaunator = await loadDelaunator();
  const px = [];
  const py = [];
  const pv = [];
  for (let i = 0; i < xs.length; i += 1) {
    if (Number.isFinite(xs[i]) && Number.isFinite(ys[i]) && Number.isFinite(values[i])) {
      px.push(xs[i]);
      py.push(ys[i]);
      pv.push(values[i]);
    }
  }
  if (px.length < 3) return null;
  const coords = new Float64Array(px.length * 2);
  px.forEach((v, i) => {
    coords[2 * i] = v;
    coords[2 * i + 1] = py[i];
  });
  const tri = new Delaunator(coords).triangles;
  const xmin = Math.min(...px);
  const xmax = Math.max(...px);
  const ymin = Math.min(...py);
  const ymax = Math.max(...py);
  const span = Math.max(xmax - xmin, ymax - ymin);
  const cell = Math.max(span / 420, Math.min(spacing / 3, span / 120));
  const nx = Math.max(2, Math.ceil((xmax - xmin) / cell) + 1);
  const ny = Math.max(2, Math.ceil((ymax - ymin) / cell) + 1);
  const out = new Float64Array(nx * ny).fill(NaN);
  const maxEdge2 = (4 * spacing) ** 2;
  for (let t = 0; t < tri.length; t += 3) {
    const a = tri[t];
    const b = tri[t + 1];
    const c = tri[t + 2];
    const ax = px[a];
    const ay = py[a];
    const bx = px[b];
    const by = py[b];
    const cx = px[c];
    const cy = py[c];
    if ((ax - bx) ** 2 + (ay - by) ** 2 > maxEdge2
      || (bx - cx) ** 2 + (by - cy) ** 2 > maxEdge2
      || (cx - ax) ** 2 + (cy - ay) ** 2 > maxEdge2) continue;
    const det = (by - cy) * (ax - cx) + (cx - bx) * (ay - cy);
    if (Math.abs(det) < 1e-18) continue;
    const i0 = Math.max(0, Math.floor((Math.min(ax, bx, cx) - xmin) / cell));
    const i1 = Math.min(nx - 1, Math.ceil((Math.max(ax, bx, cx) - xmin) / cell));
    const j0 = Math.max(0, Math.floor((Math.min(ay, by, cy) - ymin) / cell));
    const j1 = Math.min(ny - 1, Math.ceil((Math.max(ay, by, cy) - ymin) / cell));
    for (let j = j0; j <= j1; j += 1) {
      const y = ymin + j * cell;
      for (let i = i0; i <= i1; i += 1) {
        const x = xmin + i * cell;
        const l1 = ((by - cy) * (x - cx) + (cx - bx) * (y - cy)) / det;
        const l2 = ((cy - ay) * (x - cx) + (ax - cx) * (y - cy)) / det;
        const l3 = 1 - l1 - l2;
        if (l1 < -1e-9 || l2 < -1e-9 || l3 < -1e-9) continue;
        out[j * nx + i] = l1 * pv[a] + l2 * pv[b] + l3 * pv[c];
      }
    }
  }
  return {nx, ny, x0: xmin, y0: ymin, dx: cell, dy: cell, values: out};
}

// ---------------------------------------------------------------- contours

// Marching squares. Returns, per level, a flat array of segment endpoints
// [x1, y1, x2, y2, ...] in grid world coordinates.
export function marchingSquares(grid, values, levels) {
  const {nx, ny, x0, y0, dx, dy} = grid;
  const out = [];
  for (const level of levels) {
    const seg = [];
    for (let j = 0; j < ny - 1; j += 1) {
      for (let i = 0; i < nx - 1; i += 1) {
        const v0 = values[j * nx + i];
        const v1 = values[j * nx + i + 1];
        const v2 = values[(j + 1) * nx + i + 1];
        const v3 = values[(j + 1) * nx + i];
        if (!Number.isFinite(v0) || !Number.isFinite(v1) || !Number.isFinite(v2) || !Number.isFinite(v3)) continue;
        let code = 0;
        if (v0 > level) code |= 1;
        if (v1 > level) code |= 2;
        if (v2 > level) code |= 4;
        if (v3 > level) code |= 8;
        if (code === 0 || code === 15) continue;
        const xa = x0 + i * dx;
        const ya = y0 + j * dy;
        const lerp = (a, b) => (level - a) / (b - a);
        const e = [
          () => [xa + lerp(v0, v1) * dx, ya],
          () => [xa + dx, ya + lerp(v1, v2) * dy],
          () => [xa + lerp(v3, v2) * dx, ya + dy],
          () => [xa, ya + lerp(v0, v3) * dy]
        ];
        const push = (p, q) => {
          const a = e[p]();
          const b = e[q]();
          seg.push(a[0], a[1], b[0], b[1]);
        };
        const centre = (v0 + v1 + v2 + v3) / 4 > level;
        switch (code) {
          case 1: case 14: push(0, 3); break;
          case 2: case 13: push(0, 1); break;
          case 3: case 12: push(1, 3); break;
          case 4: case 11: push(1, 2); break;
          case 6: case 9: push(0, 2); break;
          case 7: case 8: push(2, 3); break;
          case 5:
            if (centre) { push(0, 1); push(2, 3); } else { push(0, 3); push(1, 2); }
            break;
          case 10:
            if (centre) { push(0, 3); push(1, 2); } else { push(0, 1); push(2, 3); }
            break;
          default: break;
        }
      }
    }
    if (seg.length) out.push({level, segments: Float32Array.from(seg)});
  }
  return out;
}

// Isophote levels spaced by `stepMag` magnitudes below the peak.
export function magnitudeLevels(values, stepMag = 0.5, maxLevels = 12) {
  let peak = -Infinity;
  let floor = Infinity;
  for (const v of values) {
    if (Number.isFinite(v) && v > 0) {
      if (v > peak) peak = v;
      if (v < floor) floor = v;
    }
  }
  if (!Number.isFinite(peak)) return [];
  const levels = [];
  for (let k = 1; k <= maxLevels; k += 1) {
    const level = peak * 10 ** (-0.4 * stepMag * k);
    if (level <= floor) break;
    levels.push(level);
  }
  return levels;
}

// ---------------------------------------------------------------- MGE

// Surface brightness of a set of MGE Gaussians [{sigma, q, lum}] on a grid
// in the display frame, with the major axis at `paDeg` from +x.
export function mgeModelGrid(gaussians, extent, paDeg = 0, n = 220) {
  const [xmin, xmax, ymin, ymax] = extent;
  const dx = (xmax - xmin) / (n - 1);
  const dy = (ymax - ymin) / (n - 1);
  const [c, s] = rotationTrig(paDeg);
  const values = new Float64Array(n * n);
  for (let j = 0; j < n; j += 1) {
    const y = ymin + j * dy;
    for (let i = 0; i < n; i += 1) {
      const x = xmin + i * dx;
      const xp = x * c + y * s;
      const yp = -x * s + y * c;
      let sum = 0;
      for (const g of gaussians) {
        const amp = g.lum / (2 * Math.PI * g.sigma * g.sigma * g.q);
        sum += amp * Math.exp(-0.5 * ((xp / g.sigma) ** 2 + (yp / (g.sigma * g.q)) ** 2));
      }
      values[j * n + i] = sum;
    }
  }
  return {grid: {nx: n, ny: n, x0: xmin, y0: ymin, dx, dy}, values};
}
