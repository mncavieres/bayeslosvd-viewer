// Converts every supported input (BAYES-LOSVD HDF5 summaries, viewer bundles,
// legacy dashboard payloads, CSV tables, MGE tables) into one result model:
//
// {
//   uid, label, source, kind, fitType, fitMode, regularization, regularized,
//   attrs, scalars, nbins,
//   x, y, snr, binFlux           Float64Array(nbins), raw (unrotated) arcsec
//   hasOutput                    Uint8Array(nbins)
//   xvel, nvel, rows, losvd      sorted velocity axis; Float64Array(nbins*rows*nvel), unit area
//   bands                        true when percentile rows differ (sampling runs)
//   quantities                   Map key -> Float64Array(nbins)
//   paramErrors                  Map key -> {lo, hi} (16/84th percentile rows)
//   grid                         spaxel lattice {nx, ny, x0, y0, dx, dy, bin, flux} or null
//   staticOverlays               [{name, frame:"raw", contours}] from legacy dashboards
//   spectra                      {fileId} when spectra can be fetched from the worker
//   wcs                          {cd: [CD1_1, CD1_2, CD2_1, CD2_2] (deg/pix, rows East, North), source} or null
// }

import {buildSpaxelGrid, losvdMoments, medianSpacing, rotatePoint, rotationTrig, sortVelocityAxis, trapz} from "./analysis.js";
import {apply, skyMatrix} from "./frames.js";

export const BUNDLE_FORMAT = "bayeslosvd-viewer-bundle";
export const PERCENTILES = [0.1, 15.9, 50.0, 84.1, 99.9];
const MEDIAN_ROW = 2;

let uidCounter = 0;
function nextUid() {
  uidCounter += 1;
  return `r${uidCounter}`;
}

function stripExtension(name) {
  return String(name || "result").replace(/\.(h5|hdf5|hdf|he5|json|csv|html?)$/i, "");
}

function defaultLabel(fileName) {
  return stripExtension(fileName).replace(/_results$/i, "");
}

export function isRegularized(...parts) {
  return parts.some((p) => /regulari[sz]/i.test(String(p || "")));
}

function nanArray(n) {
  return new Float64Array(n).fill(NaN);
}

function toFloatArray(values, n = values ? values.length : 0) {
  const out = nanArray(n);
  if (!values) return out;
  for (let i = 0; i < n; i += 1) {
    const v = values[i];
    out[i] = v === null || v === undefined || v === "" ? NaN : Number(v);
  }
  return out;
}

// Sort the LOSVD velocity axis, normalize each bin to unit area (using the
// median row) and derive LOSVD moments.
function finalizeLosvd(result, xvelRaw, losvdRaw, rows) {
  const nbins = result.nbins;
  const nvel = xvelRaw ? xvelRaw.length : 0;
  result.nvel = nvel;
  result.rows = rows;
  result.bands = false;
  ["vel", "sigma", "skew", "kurt"].forEach((k) => result.quantities.set(k, nanArray(nbins)));
  if (!nvel || !losvdRaw) {
    result.xvel = new Float64Array(0);
    result.losvd = null;
    return;
  }
  const {order, xvel} = sortVelocityAxis(xvelRaw);
  result.xvel = xvel;
  const losvd = nanArray(nbins * rows * nvel);
  const medRow = rows > MEDIAN_ROW ? MEDIAN_ROW : 0;
  const vel = result.quantities.get("vel");
  const sigma = result.quantities.get("sigma");
  const skew = result.quantities.get("skew");
  const kurt = result.quantities.get("kurt");
  for (let b = 0; b < nbins; b += 1) {
    if (!result.hasOutput[b]) continue;
    for (let r = 0; r < rows; r += 1) {
      for (let v = 0; v < nvel; v += 1) {
        losvd[(b * rows + r) * nvel + v] = losvdRaw[(b * rows + r) * nvel + order[v]];
      }
    }
    const medOffset = (b * rows + medRow) * nvel;
    const norm = trapz(losvd, xvel, medOffset, nvel);
    if (Number.isFinite(norm) && norm !== 0) {
      for (let k = b * rows * nvel; k < (b + 1) * rows * nvel; k += 1) losvd[k] /= norm;
    }
    if (!result.bands && rows > 1) {
      for (let r = 0; r < rows && !result.bands; r += 1) {
        if (r === medRow) continue;
        for (let v = 0; v < nvel; v += 1) {
          const a = losvd[(b * rows + r) * nvel + v];
          const m = losvd[medOffset + v];
          if (Number.isFinite(a) && Number.isFinite(m) && Math.abs(a - m) > 1e-9 * (Math.abs(m) + 1e-12)) {
            result.bands = true;
            break;
          }
        }
      }
    }
    const mom = losvdMoments(xvel, losvd, medOffset);
    vel[b] = mom.vel;
    sigma[b] = mom.sigma;
    skew[b] = mom.skew;
    kurt[b] = mom.kurt;
  }
  result.losvd = losvd;
  result.medianRow = medRow;
}

function baseResult(fields) {
  const nbins = fields.x.length;
  const result = {
    uid: nextUid(),
    label: fields.label,
    source: fields.source,
    kind: fields.kind,
    fitType: fields.fitType || "",
    fitMode: fields.fitMode || "",
    regularization: fields.regularization || "",
    attrs: fields.attrs || {},
    scalars: fields.scalars || {},
    nbins,
    x: fields.x,
    y: fields.y,
    snr: fields.snr || nanArray(nbins),
    binFlux: fields.binFlux || nanArray(nbins),
    hasOutput: fields.hasOutput || new Uint8Array(nbins).fill(1),
    quantities: new Map(),
    paramErrors: new Map(),
    grid: fields.grid || null,
    staticOverlays: fields.staticOverlays || [],
    spectra: fields.spectra || null,
    notes: fields.notes || [],
    wcs: fields.wcs || null
  };
  result.regularized = isRegularized(result.regularization, result.fitType, result.label, result.attrs.model_code);
  result.r = Float64Array.from(result.x, (xv, i) => Math.hypot(xv, result.y[i]));
  result.spacing = medianSpacing(result.x, result.y);
  return result;
}

function addCommonQuantities(result) {
  if (Array.from(result.snr).some(Number.isFinite)) result.quantities.set("snr", result.snr);
  if (Array.from(result.binFlux).some((v) => Number.isFinite(v) && v > 0)) {
    result.quantities.set("flux", Float64Array.from(result.binFlux, (v) => (v > 0 ? Math.log10(v) : NaN)));
  } else if (result.grid?.flux && !result.grid.imageOnly) {
    // derive per-bin flux from the spaxel lattice
    const sum = new Float64Array(result.nbins);
    const g = result.grid;
    for (let k = 0; k < g.bin.length; k += 1) {
      const b = g.bin[k];
      if (b >= 0 && b < result.nbins && Number.isFinite(g.flux[k])) sum[b] += g.flux[k];
    }
    result.binFlux = sum;
    result.quantities.set("flux", Float64Array.from(sum, (v) => (v > 0 ? Math.log10(v) : NaN)));
  }
}

// ---------------------------------------------------------------- HDF5

// Optional sky orientation stored as root attributes WCS_CD1_1 ... WCS_CD2_2
// (the spatial CD matrix of the cube, deg/pixel; see tools/add_wcs.py).
function wcsFromAttrs(attrs, fileName) {
  const get = (k) => {
    const key = Object.keys(attrs).find((a) => a.toUpperCase() === k);
    return key === undefined ? NaN : Number(attrs[key]);
  };
  const cd = ["WCS_CD1_1", "WCS_CD1_2", "WCS_CD2_1", "WCS_CD2_2"].map(get);
  if (!cd.every(Number.isFinite) || cd[0] * cd[3] - cd[1] * cd[2] === 0) return null;
  const src = Object.keys(attrs).find((a) => a.toUpperCase() === "WCS_SOURCE");
  return {cd, source: src ? `${attrs[src]} (via ${fileName})` : `attributes of ${fileName}`};
}

const PARAM_QUANTITIES = {
  vel: "fit_vel", sigma: "fit_sigma", h1: "fit_h1", h2: "fit_h2", h3: "fit_h3", h4: "fit_h4", h5: "fit_h5", h6: "fit_h6"
};

export function resultFromHdf5Summary(summary, {fileId, fileName, label}) {
  const s = summary;
  const nbins = s.nbins;
  let grid = null;
  if (s.spaxelX && s.spaxelY && s.spaxelX.length === s.spaxelY.length) {
    grid = buildSpaxelGrid(s.spaxelX, s.spaxelY, s.spaxelBin, s.spaxelFlux);
  }
  const attrs = s.attrs || {};
  const result = baseResult({
    label: label || defaultLabel(fileName),
    source: fileName,
    kind: "hdf5",
    fitType: attrs.fit_type,
    fitMode: attrs.fit_mode,
    regularization: attrs.regularization,
    attrs,
    scalars: s.scalars,
    x: s.xbin,
    y: s.ybin,
    snr: s.binSnr && s.binSnr.length === nbins ? s.binSnr : null,
    binFlux: s.binFlux && s.binFlux.length === nbins ? s.binFlux : null,
    hasOutput: s.hasOutput,
    grid,
    spectra: s.hasSpectra ? {fileId} : null,
    wcs: wcsFromAttrs(attrs, fileName)
  });
  finalizeLosvd(result, s.xvel, s.losvd, s.rows);

  // Model parameters stored per bin (e.g. vel/sigma/h3/h4 for GH fits).
  Object.entries(s.params || {}).forEach(([key, flat]) => {
    const qkey = PARAM_QUANTITIES[key];
    if (!qkey) return;
    const rows = s.paramRows[key] || 1;
    const medRow = rows >= 5 ? MEDIAN_ROW : 0;
    const median = nanArray(nbins);
    const lo = nanArray(nbins);
    const hi = nanArray(nbins);
    for (let b = 0; b < nbins; b += 1) {
      median[b] = flat[b * 5 + medRow];
      if (rows >= 5) {
        lo[b] = flat[b * 5 + 1];
        hi[b] = flat[b * 5 + 3];
      }
    }
    // Skip parameters that are not physically meaningful maps (e.g. a "sigma"
    // hyper-parameter of a regularized model is ~0.05, not km/s).
    if (key === "sigma" && !s.params.vel) return;
    if (!median.some(Number.isFinite)) return;
    result.quantities.set(qkey, median);
    if (rows >= 5 && lo.some((v, i) => v !== hi[i])) result.paramErrors.set(qkey, {lo, hi});
  });
  addCommonQuantities(result);
  const done = result.hasOutput.reduce((a, b) => a + b, 0);
  if (done < nbins) result.notes.push(`${done} of ${nbins} bins have output (partial run)`);
  return result;
}

// ---------------------------------------------------------------- legacy dashboards

// The original build_kinematics_dashboard.py stores bin coordinates rotated
// by DATA.rotationDeg; undo that so all results share the raw frame.
export function resultsFromLegacyPayload(data, fileName) {
  const rotation = Number(data.rotationDeg) || 0;
  const [c, s] = rotationTrig(-rotation);
  const unrotate = (x, y) => rotatePoint(x, y, c, s);
  const results = [];
  let imageGrid = null;
  if (data.image && Array.isArray(data.image.values) && data.image.shape) {
    const [ny, nx] = data.image.shape;
    const pix = Number(data.image.pixscaleArcsec) || 0.1;
    const center = data.image.center || {row: 0, col: 0};
    // Image pixels live in the raw frame: x = (col - c0) * pix, y = (row - r0) * pix
    imageGrid = {
      nx,
      ny,
      x0: -center.col * pix,
      y0: -center.row * pix,
      dx: pix,
      dy: pix,
      bin: new Int32Array(nx * ny).fill(-1),
      flux: toFloatArray(data.image.values, nx * ny),
      imageOnly: true
    };
  }
  const staticOverlays = [];
  if (data.mge && Array.isArray(data.mge.contours) && data.mge.contours.length) {
    staticOverlays.push({
      name: "MGE (from dashboard)",
      kind: "mge",
      contours: data.mge.contours.map((contour) => ({
        level: contour.level,
        polylines: contour.segments.map((seg) => seg.map(([x, y]) => unrotate(x, y)))
      }))
    });
  }

  (data.results || []).forEach((r, idx) => {
    const bins = r.bins || [];
    const n = bins.length;
    const x = nanArray(n);
    const y = nanArray(n);
    const snr = nanArray(n);
    bins.forEach((bin, i) => {
      const [xr, yr] = unrotate(Number(bin.x), Number(bin.y));
      x[i] = xr;
      y[i] = yr;
      snr[i] = bin.snr === null || bin.snr === undefined ? NaN : Number(bin.snr);
    });
    const hasOutput = new Uint8Array(n);
    const xvel = r.xvel && r.xvel.length ? Float64Array.from(r.xvel) : null;
    let losvd = null;
    if (xvel && r.losvds) {
      losvd = nanArray(n * xvel.length);
      Object.entries(r.losvds).forEach(([key, arr]) => {
        const b = Number(key);
        if (!Number.isInteger(b) || b < 0 || b >= n || !arr) return;
        hasOutput[b] = 1;
        for (let v = 0; v < xvel.length; v += 1) losvd[b * xvel.length + v] = arr[v] === null ? NaN : arr[v];
      });
    } else {
      hasOutput.fill(1);
    }
    const result = baseResult({
      label: r.label || `${defaultLabel(fileName)} ${idx + 1}`,
      source: r.source || fileName,
      kind: "legacy",
      fitType: r.fitType,
      regularization: r.regularization,
      x,
      y,
      snr,
      hasOutput,
      grid: imageGrid,
      staticOverlays,
      notes: [`imported from legacy dashboard ${fileName}`]
    });
    finalizeLosvd(result, xvel, losvd, 1);
    if (!losvd && r.moments) {
      // No LOSVDs: take the moment maps as stored (legacy h3/h4 are skew/kurt).
      const map = {vel: "vel", sigma: "sigma", h3: "skew", h4: "kurt"};
      Object.entries(map).forEach(([from, to]) => {
        if (r.moments[from]) result.quantities.set(to, toFloatArray(r.moments[from], n));
      });
    }
    addCommonQuantities(result);
    results.push(result);
  });
  return {results, settings: {rotation}};
}

export function extractLegacyPayloadFromHtml(text) {
  const marker = "const DATA = ";
  const start = text.indexOf(marker);
  if (start < 0) return null;
  const from = start + marker.length;
  const end = text.indexOf(";\n", from);
  if (end < 0) return null;
  return JSON.parse(text.slice(from, end));
}

// ---------------------------------------------------------------- CSV

function parseCsvLine(line, delimiter) {
  const values = [];
  let current = "";
  let quoted = false;
  for (let i = 0; i < line.length; i += 1) {
    const ch = line[i];
    if (ch === '"') {
      if (quoted && line[i + 1] === '"') {
        current += '"';
        i += 1;
      } else {
        quoted = !quoted;
      }
    } else if (ch === delimiter && !quoted) {
      values.push(current.trim());
      current = "";
    } else {
      current += ch;
    }
  }
  values.push(current.trim());
  return values;
}

export function parseTable(text) {
  const lines = text.split(/\r?\n/).filter((l) => l.trim() && !l.trim().startsWith("#"));
  if (!lines.length) throw new Error("The table is empty");
  const delimiter = lines[0].includes(",") ? "," : (lines[0].includes("\t") ? "\t" : null);
  const split = (line) => (delimiter ? parseCsvLine(line, delimiter) : line.trim().split(/\s+/));
  const headers = split(lines[0]);
  const rows = lines.slice(1).map(split);
  return {headers, rows};
}

const norm = (name) => String(name || "").trim().toLowerCase().replace(/[^a-z0-9]+/g, "");

function findColumn(headers, aliases) {
  const normalized = headers.map(norm);
  for (const alias of aliases) {
    const idx = normalized.indexOf(norm(alias));
    if (idx >= 0) return idx;
  }
  return -1;
}

export function looksLikeMgeTable(headers) {
  return findColumn(headers, ["sigma_arcsec", "sigmaarcsec", "sigma"]) >= 0
    && findColumn(headers, ["q_obs", "qobs", "q"]) >= 0
    && findColumn(headers, ["x", "xrot", "x_arcsec"]) < 0;
}

export function mgeFromTable(text, fileName) {
  const {headers, rows} = parseTable(text);
  const iSigma = findColumn(headers, ["sigma_arcsec", "sigmaarcsec", "sigma"]);
  const iQ = findColumn(headers, ["q_obs", "qobs", "q"]);
  const iL = findColumn(headers, ["luminosity_lsun", "luminosity", "lum", "total_flux", "total_counts", "counts", "flux", "l"]);
  const iSurf = findColumn(headers, ["surf", "surface_density", "sigma0", "i0", "peak"]);
  if (iSigma < 0 || iQ < 0) throw new Error("MGE table needs sigma and q columns");
  const gaussians = [];
  rows.forEach((row) => {
    const sigma = Number(row[iSigma]);
    const q = Number(row[iQ]);
    let lum = iL >= 0 ? Number(row[iL]) : NaN;
    if (!Number.isFinite(lum) && iSurf >= 0) lum = 2 * Math.PI * sigma * sigma * q * Number(row[iSurf]);
    if (!Number.isFinite(lum)) lum = 1;
    if (sigma > 0 && q > 0 && lum > 0) gaussians.push({sigma, q, lum});
  });
  if (!gaussians.length) throw new Error("MGE table has no valid Gaussians");
  return {source: fileName, gaussians};
}

// Generic bin table: x/y columns plus any number of numeric columns.
export function resultFromTable(text, fileName) {
  const {headers, rows} = parseTable(text);
  const ix = findColumn(headers, ["x", "xbin", "x_arcsec", "xarcsec", "xrot", "x_rot", "xrotarcsec"]);
  const iy = findColumn(headers, ["y", "ybin", "y_arcsec", "yarcsec", "yrot", "y_rot", "yrotarcsec"]);
  if (ix < 0 || iy < 0) throw new Error("The table needs x and y columns (arcsec)");
  const aliases = {
    vel: ["vel", "velocity", "v", "meanvelocity", "v_kms"],
    sigma: ["sigma", "dispersion", "sig", "sigma_kms"],
    h3: ["h3"],
    h4: ["h4"],
    skew: ["skew", "skewness"],
    kurt: ["kurt", "kurtosis"],
    snr: ["snr", "s/n", "sn", "signaltonoise"]
  };
  const used = new Set([ix, iy]);
  const columns = [];
  Object.entries(aliases).forEach(([key, names]) => {
    const i = findColumn(headers, names);
    if (i >= 0 && !used.has(i)) {
      used.add(i);
      columns.push([key === "h3" || key === "h4" ? `table_${key}` : key, i]);
    }
  });
  const iBin = findColumn(headers, ["bin", "binid", "bin_id", "id"]);
  if (iBin >= 0) used.add(iBin);
  headers.forEach((h, i) => {
    if (!used.has(i)) columns.push([`col:${h}`, i]);
  });
  const keep = rows.filter((row) => Number.isFinite(Number(row[ix])) && Number.isFinite(Number(row[iy])));
  if (!keep.length) throw new Error("The table has no rows with finite x and y");
  const x = Float64Array.from(keep, (row) => Number(row[ix]));
  const y = Float64Array.from(keep, (row) => Number(row[iy]));
  const result = baseResult({label: defaultLabel(fileName), source: fileName, kind: "table", fitType: "table", x, y});
  finalizeLosvd(result, null, null, 1);
  result.quantities.clear();
  columns.forEach(([key, i]) => {
    const values = Float64Array.from(keep, (row) => (row[i] === "" ? NaN : Number(row[i])));
    if (!values.some(Number.isFinite)) return;
    if (key === "snr") result.snr = values;
    result.quantities.set(key, values);
  });
  return result;
}

// ---------------------------------------------------------------- bundles

function round(values, digits = 6) {
  if (!values) return null;
  return Array.from(values, (v) => (Number.isFinite(v) ? Number(v.toPrecision(digits)) : null));
}

export function makeBundle(results, {settings, mge}) {
  return {
    format: BUNDLE_FORMAT,
    version: 1,
    created: new Date().toISOString(),
    settings,
    mge: mge ? {source: mge.source, gaussians: mge.gaussians.map((g) => [g.sigma, g.q, g.lum])} : null,
    results: results.map((r) => {
      // Drop duplicated percentile rows (MAP runs store five identical rows).
      const rows = r.bands ? r.rows : 1;
      let losvd = null;
      if (r.losvd) {
        const out = new Float64Array(r.nbins * rows * r.nvel);
        for (let b = 0; b < r.nbins; b += 1) {
          for (let k = 0; k < rows; k += 1) {
            const src = (b * r.rows + (rows === 1 ? r.medianRow : k)) * r.nvel;
            out.set(r.losvd.subarray(src, src + r.nvel), (b * rows + k) * r.nvel);
          }
        }
        losvd = round(out, 7);
      }
      const quantities = {};
      r.quantities.forEach((values, key) => {
        if (["vel", "sigma", "skew", "kurt", "flux", "snr"].includes(key) && r.losvd) return;
        if (key.startsWith("gh_")) return;
        quantities[key] = round(values);
      });
      return {
        label: r.label,
        source: r.source,
        fitType: r.fitType,
        fitMode: r.fitMode,
        regularization: r.regularization,
        attrs: r.attrs,
        scalars: r.scalars,
        x: round(r.x),
        y: round(r.y),
        snr: round(r.snr),
        binFlux: round(r.binFlux),
        hasOutput: Array.from(r.hasOutput),
        xvel: round(r.xvel),
        rows,
        losvd,
        quantities,
        paramErrors: Object.fromEntries(Array.from(r.paramErrors, ([k, v]) => [k, {lo: round(v.lo), hi: round(v.hi)}])),
        wcs: r.wcs,
        grid: r.grid ? {
          nx: r.grid.nx,
          ny: r.grid.ny,
          x0: r.grid.x0,
          y0: r.grid.y0,
          dx: r.grid.dx,
          dy: r.grid.dy,
          bin: Array.from(r.grid.bin),
          flux: round(r.grid.flux, 5),
          imageOnly: Boolean(r.grid.imageOnly)
        } : null,
        notes: r.notes
      };
    })
  };
}

export function resultsFromBundle(bundle, fileName) {
  if (bundle.version !== 1) throw new Error(`Unsupported bundle version ${bundle.version}`);
  const results = (bundle.results || []).map((r) => {
    const n = r.x.length;
    const grid = r.grid ? {
      ...r.grid,
      bin: Int32Array.from(r.grid.bin),
      flux: r.grid.flux ? toFloatArray(r.grid.flux) : null
    } : null;
    const result = baseResult({
      label: r.label,
      source: r.source || fileName,
      kind: "bundle",
      fitType: r.fitType,
      fitMode: r.fitMode,
      regularization: r.regularization,
      attrs: r.attrs,
      scalars: r.scalars,
      x: toFloatArray(r.x),
      y: toFloatArray(r.y),
      snr: toFloatArray(r.snr, n),
      binFlux: toFloatArray(r.binFlux, n),
      hasOutput: Uint8Array.from(r.hasOutput || new Array(n).fill(1)),
      grid,
      notes: [...(r.notes || []), `loaded from bundle ${fileName}`],
      wcs: r.wcs && Array.isArray(r.wcs.cd) ? r.wcs : null
    });
    // Bundles store an already sorted axis; finalizeLosvd re-sorts (no-op) and renormalizes.
    finalizeLosvd(result, r.xvel ? toFloatArray(r.xvel) : null, r.losvd ? toFloatArray(r.losvd) : null, r.rows || 1);
    Object.entries(r.quantities || {}).forEach(([key, values]) => result.quantities.set(key, toFloatArray(values, n)));
    Object.entries(r.paramErrors || {}).forEach(([key, v]) => {
      result.paramErrors.set(key, {lo: toFloatArray(v.lo, n), hi: toFloatArray(v.hi, n)});
    });
    addCommonQuantities(result);
    return result;
  });
  const mge = bundle.mge ? {
    source: bundle.mge.source,
    gaussians: bundle.mge.gaussians.map(([sigma, q, lum]) => ({sigma, q, lum}))
  } : null;
  return {results, settings: bundle.settings || {}, mge};
}

// ---------------------------------------------------------------- exports

// displayMatrix: raw (x, y) -> displayed frame (2x2, row-major).
export function resultToCsv(result, displayMatrix, quantityKeys) {
  const M = displayMatrix || [1, 0, 0, 1];
  const W = result.wcs ? skyMatrix(result.wcs.cd) : null;
  const sky = W ? ["dRA_arcsec", "dDec_arcsec"] : [];
  const header = ["bin", "x_arcsec", "y_arcsec", "x_display", "y_display", ...sky, "has_output", ...quantityKeys];
  const lines = [header.join(",")];
  for (let b = 0; b < result.nbins; b += 1) {
    const [xd, yd] = apply(M, result.x[b], result.y[b]);
    const skyCols = W ? apply(W, result.x[b], result.y[b]).map((v) => v.toFixed(6)) : [];
    const vals = quantityKeys.map((k) => {
      const v = result.quantities.get(k)?.[b];
      return Number.isFinite(v) ? Number(v.toPrecision(8)) : "";
    });
    lines.push([b, result.x[b], result.y[b], xd.toFixed(6), yd.toFixed(6), ...skyCols, result.hasOutput[b], ...vals].join(","));
  }
  return `${lines.join("\n")}\n`;
}
