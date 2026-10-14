import {
  fitGaussHermite, interpolateBins, limitsFor, magnitudeLevels, marchingSquares, matchBins, mgeModelGrid,
  rotatePoint, rotationTrig, gridCorners
} from "./analysis.js";
import {COLORMAP_NAMES, getLut, lutIndex} from "./colormaps.js";
import {
  BUNDLE_FORMAT, PERCENTILES, extractLegacyPayloadFromHtml, looksLikeMgeTable, makeBundle, mgeFromTable, parseTable,
  resultFromHdf5Summary, resultFromTable, resultToCsv, resultsFromBundle, resultsFromLegacyPayload
} from "./data.js";
import {Hdf5Client} from "./hdf5-client.js";
import {drawAxes, drawColorbar, formatValue, setupCanvas, theme} from "./plot.js";

const REPO_URL = "https://github.com/mncavieres/bayeslosvd-viewer";
const DEMO_FILES = [
  "demo/demo_disk_sampling_results.hdf5",
  "demo/demo_disk_regularized_results.hdf5"
];
const DEMO_MGE = "demo/demo_mge.csv";
const DEMO_ROTATION = -25; // aligns the synthetic disc's major axis with +x

const BIN_COLORS = ["#0f766e", "#b6465f", "#2f5aa8", "#c27418", "#6f4aa7", "#54733c", "#8f3b32", "#2f747f"];
const RESULT_COLORS = ["#2563eb", "#dc2626", "#059669", "#9333ea", "#ea580c", "#0891b2", "#be185d", "#65a30d"];
const MAX_SELECTED = 8;

// ------------------------------------------------------------------ quantities

const QUANTITIES = {
  vel: {label: "V", name: "Mean velocity", unit: "km/s", cmap: "RdBu_r", limits: "symmetric", group: "LOSVD moments",
    desc: "First moment of the median LOSVD"},
  sigma: {label: "σ", name: "Velocity dispersion", unit: "km/s", cmap: "magma", limits: "robust", group: "LOSVD moments",
    desc: "Square root of the second central moment of the median LOSVD"},
  skew: {label: "Skew", name: "LOSVD skewness", unit: "", cmap: "coolwarm", limits: "symmetric", group: "LOSVD moments",
    desc: "Standardized third moment (≈ 4√3·h3 for small h3)"},
  kurt: {label: "Kurt", name: "LOSVD excess kurtosis", unit: "", cmap: "coolwarm", limits: "symmetric", group: "LOSVD moments",
    desc: "Standardized fourth moment minus 3 (≈ 8√6·h4 for small h4)"},
  gh_vel: {label: "V", name: "Gauss–Hermite V", unit: "km/s", cmap: "RdBu_r", limits: "symmetric", group: "Gauss–Hermite fit",
    desc: "Least-squares Gauss–Hermite fit (to h4) of the median LOSVD"},
  gh_sigma: {label: "σ", name: "Gauss–Hermite σ", unit: "km/s", cmap: "magma", limits: "robust", group: "Gauss–Hermite fit",
    desc: "Least-squares Gauss–Hermite fit (to h4) of the median LOSVD"},
  gh_h3: {label: "h3", name: "Gauss–Hermite h3", unit: "", cmap: "coolwarm", limits: "symmetric", group: "Gauss–Hermite fit",
    desc: "Least-squares Gauss–Hermite fit (to h4) of the median LOSVD"},
  gh_h4: {label: "h4", name: "Gauss–Hermite h4", unit: "", cmap: "coolwarm", limits: "symmetric", group: "Gauss–Hermite fit",
    desc: "Least-squares Gauss–Hermite fit (to h4) of the median LOSVD"},
  fit_vel: {label: "V", name: "Model V", unit: "km/s", cmap: "RdBu_r", limits: "symmetric", group: "Model parameters",
    desc: "Parameter stored in out/<bin>/vel (parametric fits)"},
  fit_sigma: {label: "σ", name: "Model σ", unit: "km/s", cmap: "magma", limits: "robust", group: "Model parameters",
    desc: "Parameter stored in out/<bin>/sigma (parametric fits)"},
  fit_h1: {label: "h1", name: "Model h1", unit: "", cmap: "coolwarm", limits: "symmetric", group: "Model parameters",
    desc: "Parameter stored in out/<bin>/h1 (GHfree fits)"},
  fit_h2: {label: "h2", name: "Model h2", unit: "", cmap: "coolwarm", limits: "symmetric", group: "Model parameters",
    desc: "Parameter stored in out/<bin>/h2 (GHfree fits)"},
  fit_h3: {label: "h3", name: "Model h3", unit: "", cmap: "coolwarm", limits: "symmetric", group: "Model parameters",
    desc: "Parameter stored in out/<bin>/h3"},
  fit_h4: {label: "h4", name: "Model h4", unit: "", cmap: "coolwarm", limits: "symmetric", group: "Model parameters",
    desc: "Parameter stored in out/<bin>/h4"},
  fit_h5: {label: "h5", name: "Model h5", unit: "", cmap: "coolwarm", limits: "symmetric", group: "Model parameters", desc: ""},
  fit_h6: {label: "h6", name: "Model h6", unit: "", cmap: "coolwarm", limits: "symmetric", group: "Model parameters", desc: ""},
  table_h3: {label: "h3", name: "h3 (table)", unit: "", cmap: "coolwarm", limits: "symmetric", group: "Table", desc: ""},
  table_h4: {label: "h4", name: "h4 (table)", unit: "", cmap: "coolwarm", limits: "symmetric", group: "Table", desc: ""},
  snr: {label: "S/N", name: "Bin S/N", unit: "", cmap: "viridis", limits: "robust", group: "Data", desc: "in/bin_snr"},
  flux: {label: "Flux", name: "log₁₀ flux", unit: "dex", cmap: "gray", limits: "robust", group: "Data",
    desc: "Spaxel flux (in/flux) or bin flux (in/bin_flux), log scale"}
};
const QUANTITY_ORDER = Object.keys(QUANTITIES);

function quantityInfo(key) {
  if (QUANTITIES[key]) return QUANTITIES[key];
  const name = key.replace(/^col:/, "");
  return {label: name, name, unit: "", cmap: "viridis", limits: "robust", group: "Table", desc: `Column "${name}"`};
}

// ------------------------------------------------------------------ state

const app = {results: [], client: null, mge: null};
const state = {
  activeUid: null,
  compareUid: "",
  quantity: "vel",
  view: "spaxels",
  cmap: {},
  reverse: false,
  limitsMode: null,
  vmin: null,
  vmax: null,
  rotation: 0,
  maskSigma: false,
  maskRadius: 1.0,
  minSnr: NaN,
  showIsophotes: false,
  showCenters: false,
  showMge: false,
  mgePa: 0,
  selected: [],
  active: null,
  zoom: null,
  drag: null,
  losvdAll: false,
  specRange: null,
  specDrag: null,
  mapView: null,
  specView: null
};

const $ = (id) => document.getElementById(id);
const el = {};
[
  "openBtn", "fileInput", "exportBtn", "exportMenu", "themeBtn", "repoLink", "statusBar", "statusText", "progressFill",
  "welcome", "dropzone", "welcomeOpenBtn", "demoBtn", "workspace", "resultList", "compareSelect", "compareHint",
  "quantityButtons", "viewButtons", "cmapSelect", "cmapReverse", "vminInput", "vmaxInput", "limitButtons",
  "rotationInput", "maskSigma", "maskRadius", "minSnr", "maskHint", "showIsophotes", "showCenters", "showMge",
  "mgeLoadBtn", "mgeInput", "mgePa", "mgeStatus", "mapTitle", "mapSubtitle", "mapStatus", "mapCanvas", "mapTooltip",
  "colorbarCanvas", "zoomInBtn", "zoomOutBtn", "zoomResetBtn", "binInput", "addBinBtn", "clearBinsBtn", "binChips",
  "metrics", "losvdCanvas", "losvdLegend", "losvdAll", "losvdHint", "spectrumCanvas", "spectrumStatus", "specResetBtn",
  "resultInfoBody", "dropOverlay", "addMoreBtn"
].forEach((id) => {
  el[id] = $(id);
});

// ------------------------------------------------------------------ status

let statusTimer = null;
function setStatus(text, {progress = null, error = false, sticky = false} = {}) {
  clearTimeout(statusTimer);
  el.statusBar.hidden = !text;
  el.statusBar.classList.toggle("error", error);
  el.statusText.textContent = text || "";
  el.progressFill.parentElement.hidden = progress === null;
  if (progress !== null) el.progressFill.style.width = `${Math.round(Math.max(0, Math.min(1, progress)) * 100)}%`;
  if (text && !sticky) statusTimer = setTimeout(() => setStatus(""), error ? 12000 : 3500);
}

// ------------------------------------------------------------------ loading

function client() {
  if (!app.client) app.client = new Hdf5Client();
  return app.client;
}

function extensionOf(name) {
  const m = /\.([a-z0-9]+)$/i.exec(name || "");
  return m ? m[1].toLowerCase() : "";
}

async function loadHdf5(source, name) {
  const onProgress = (msg) => {
    if (msg.stage === "download") {
      setStatus(`Downloading ${name} — ${(msg.done / 1e6).toFixed(1)} of ${(msg.total / 1e6).toFixed(1)} MB`,
        {progress: msg.done / msg.total, sticky: true});
    } else {
      setStatus(`Reading ${name} — bin ${msg.done} of ${msg.total}`, {progress: msg.done / msg.total, sticky: true});
    }
  };
  setStatus(`Opening ${name}…`, {progress: 0, sticky: true});
  const {id, summary} = await client().open(source, onProgress);
  return {results: [resultFromHdf5Summary(summary, {fileId: id, fileName: name})], settings: null, mge: null};
}

function parseJsonInput(text, name) {
  const parsed = JSON.parse(text);
  if (parsed && parsed.format === BUNDLE_FORMAT) return resultsFromBundle(parsed, name);
  if (parsed && Array.isArray(parsed.results)) return resultsFromLegacyPayload(parsed, name);
  if (parsed && Array.isArray(parsed.bins)) return resultsFromLegacyPayload({results: [parsed]}, name);
  throw new Error("Unrecognized JSON: expected a viewer bundle or a legacy dashboard payload");
}

async function readSource(source) {
  const name = source.name;
  const ext = extensionOf(name);
  if (["h5", "hdf5", "hdf", "he5"].includes(ext)) return loadHdf5(source.file ? {file: source.file} : {url: source.url}, name);
  const text = source.file ? await source.file.text() : await fetch(source.url).then((r) => {
    if (!r.ok) throw new Error(`HTTP ${r.status} for ${source.url}`);
    return r.text();
  });
  if (ext === "json") return parseJsonInput(text, name);
  if (ext === "html" || ext === "htm") {
    const payload = extractLegacyPayloadFromHtml(text);
    if (!payload) throw new Error("No dashboard data (const DATA = …) found in this HTML file");
    return resultsFromLegacyPayload(payload, name);
  }
  const {headers} = parseTable(text);
  if (looksLikeMgeTable(headers)) return {results: [], settings: null, mge: mgeFromTable(text, name)};
  return {results: [resultFromTable(text, name)], settings: null, mge: null};
}

async function loadSources(sources) {
  let loaded = 0;
  let firstNew = null;
  for (const source of sources) {
    try {
      const {results, settings, mge} = await readSource(source);
      if (mge) setMge(mge);
      if (results.length) {
        firstNew = firstNew || results[0];
        addResults(results, settings);
        loaded += results.length;
      }
      setStatus(results.length ? `Loaded ${source.name}` : `Loaded MGE table ${source.name}`);
    } catch (error) {
      console.error(error);
      setStatus(`Could not load ${source.name}: ${error.message}`, {error: true});
      break;
    }
  }
  // Show the first file of a batch rather than the last one.
  if (firstNew && state.activeUid !== firstNew.uid) setActive(firstNew.uid);
  if (loaded) setStatus(`Loaded ${loaded} result${loaded === 1 ? "" : "s"}`);
}

function loadFiles(fileList) {
  const files = Array.from(fileList || []);
  if (!files.length) return;
  loadSources(files.map((file) => ({file, name: file.name})));
}

function addResults(results, settings) {
  const first = app.results.length === 0;
  results.forEach((r, i) => {
    r.color = RESULT_COLORS[(app.results.length + i) % RESULT_COLORS.length];
  });
  app.results.push(...results);
  if (first && settings) {
    if (Number.isFinite(settings.rotation)) state.rotation = settings.rotation;
    if (Number.isFinite(settings.mgePa)) state.mgePa = settings.mgePa;
  }
  el.welcome.hidden = true;
  el.workspace.hidden = false;
  setActive(results[0].uid, {keepSelection: !first});
  if (first && state.selected.length === 0) selectDefaultBin();
}

function removeResult(uid) {
  const idx = app.results.findIndex((r) => r.uid === uid);
  if (idx < 0) return;
  const [r] = app.results.splice(idx, 1);
  if (r.spectra && app.client) app.client.close(r.spectra.fileId).catch(() => {});
  app.results.forEach((o) => o._matchCache?.delete(uid));
  if (state.compareUid === uid) state.compareUid = "";
  if (!app.results.length) {
    state.activeUid = null;
    el.workspace.hidden = true;
    el.welcome.hidden = false;
    return;
  }
  if (state.activeUid === uid) setActive(app.results[Math.min(idx, app.results.length - 1)].uid, {keepSelection: true});
  else renderAll();
}

function setMge(mge) {
  app.mge = mge;
  app.mgeCache = null;
  state.showMge = true;
  if (!el.workspace.hidden) renderAll();
}

// ------------------------------------------------------------------ result helpers

function activeResult() {
  return app.results.find((r) => r.uid === state.activeUid) || null;
}

function compareResult() {
  if (!state.compareUid || state.compareUid === state.activeUid) return null;
  return app.results.find((r) => r.uid === state.compareUid) || null;
}

function matchFor(a, b) {
  if (a === b) return null;
  a._matchCache = a._matchCache || new Map();
  if (!a._matchCache.has(b.uid)) {
    a._matchCache.set(b.uid, matchBins(a, b, 0.5 * Math.min(a.spacing, b.spacing)));
  }
  return a._matchCache.get(b.uid);
}

function ensureGh(r) {
  if (r._gh || !r.losvd) return;
  const n = r.nbins;
  const out = {gh_vel: new Float64Array(n).fill(NaN), gh_sigma: new Float64Array(n).fill(NaN),
    gh_h3: new Float64Array(n).fill(NaN), gh_h4: new Float64Array(n).fill(NaN)};
  for (let b = 0; b < n; b += 1) {
    if (!r.hasOutput[b]) continue;
    const fit = fitGaussHermite(r.xvel, r.losvd, (b * r.rows + r.medianRow) * r.nvel);
    out.gh_vel[b] = fit.vel;
    out.gh_sigma[b] = fit.sigma;
    out.gh_h3[b] = fit.h3;
    out.gh_h4[b] = fit.h4;
  }
  Object.entries(out).forEach(([k, v]) => r.quantities.set(k, v));
  r._gh = true;
}

function availableQuantities(r) {
  const keys = new Set(r.quantities.keys());
  if (r.losvd) ["gh_vel", "gh_sigma", "gh_h3", "gh_h4"].forEach((k) => keys.add(k));
  if (r.grid?.flux) keys.add("flux");
  const ordered = QUANTITY_ORDER.filter((k) => keys.has(k));
  const extra = Array.from(keys).filter((k) => !QUANTITIES[k]);
  return [...ordered, ...extra];
}

function quantityValues(r, key) {
  if (key.startsWith("gh_")) ensureGh(r);
  return r.quantities.get(key) || null;
}

function maskKey() {
  return `${state.maskSigma}|${state.maskRadius}|${state.minSnr}`;
}

// 1 = hidden. Bins without output are always hidden.
function binMask(r) {
  const key = maskKey();
  if (r._mask && r._maskKey === key) return r._mask;
  const m = new Uint8Array(r.nbins);
  for (let b = 0; b < r.nbins; b += 1) m[b] = r.hasOutput[b] ? 0 : 1;
  r._maskInfo = null;
  if (state.maskSigma) {
    const sigma = r.quantities.get("sigma") || r.quantities.get("fit_sigma");
    if (sigma) {
      let cmax = -Infinity;
      for (let b = 0; b < r.nbins; b += 1) {
        if (r.r[b] <= state.maskRadius && Number.isFinite(sigma[b]) && !m[b]) cmax = Math.max(cmax, sigma[b]);
      }
      let count = 0;
      if (Number.isFinite(cmax)) {
        for (let b = 0; b < r.nbins; b += 1) {
          if (!m[b] && sigma[b] > cmax) {
            m[b] = 1;
            count += 1;
          }
        }
      }
      r._maskInfo = {count, cmax};
    }
  }
  if (Number.isFinite(state.minSnr)) {
    for (let b = 0; b < r.nbins; b += 1) if (!(r.snr[b] >= state.minSnr)) m[b] = 1;
  }
  r._mask = m;
  r._maskKey = key;
  return m;
}

// Per-bin values for the current quantity, with masking and differencing applied.
function binValues(r) {
  const key = state.quantity;
  const base = quantityValues(r, key);
  if (!base) return null;
  const mask = binMask(r);
  const ref = compareResult();
  const out = new Float64Array(r.nbins).fill(NaN);
  let refValues = null;
  let match = null;
  if (ref) {
    refValues = quantityValues(ref, key);
    match = refValues ? matchFor(r, ref) : null;
  }
  for (let b = 0; b < r.nbins; b += 1) {
    if (mask[b]) continue;
    if (refValues) {
      const j = match ? match[b] : b;
      out[b] = j >= 0 ? base[b] - refValues[j] : NaN;
    } else {
      out[b] = base[b];
    }
  }
  return out;
}

function isDiffMode() {
  const ref = compareResult();
  return Boolean(ref && quantityValues(ref, state.quantity));
}

function spaxelGridUsable(r) {
  return Boolean(r.grid && !r.grid.imageOnly);
}

// Decide what the map draws: a lattice image or bin markers.
function mapLayer(r) {
  const key = state.quantity;
  const values = binValues(r);
  const diff = isDiffMode();
  if (key === "flux" && r.grid?.flux && state.view !== "bins" && !diff) {
    const g = r.grid;
    const mask = binMask(r);
    const gv = new Float64Array(g.flux.length).fill(NaN);
    for (let k = 0; k < gv.length; k += 1) {
      const b = g.bin[k];
      if (b >= 0 && mask[b]) continue;
      const f = g.flux[k];
      gv[k] = f > 0 ? Math.log10(f) : NaN;
    }
    return {type: "grid", grid: g, gridValues: gv, values, smooth: false, note: "spaxel flux"};
  }
  if (!values) return {type: "none", values: null};
  if (state.view === "spaxels" && spaxelGridUsable(r)) {
    const g = r.grid;
    const gv = new Float64Array(g.bin.length).fill(NaN);
    for (let k = 0; k < gv.length; k += 1) {
      const b = g.bin[k];
      if (b >= 0 && b < values.length) gv[k] = values[b];
    }
    return {type: "grid", grid: g, gridValues: gv, values, smooth: false, footprint: true};
  }
  if (state.view === "interp") {
    const cacheKey = `${key}|${maskKey()}|${state.compareUid}`;
    r._interp = r._interp || new Map();
    const cached = r._interp.get(cacheKey);
    if (cached && cached !== "pending") {
      return {type: "grid", grid: cached, gridValues: cached.values, values, smooth: true};
    }
    if (!cached) {
      r._interp.set(cacheKey, "pending");
      interpolateBins(r.x, r.y, values, r.spacing).then((grid) => {
        if (grid) r._interp.set(cacheKey, grid);
        else r._interp.delete(cacheKey);
        renderMap();
      }).catch((error) => {
        r._interp.delete(cacheKey);
        setStatus(`Interpolation failed: ${error.message}`, {error: true});
      });
    }
    return {type: "bins", values, note: "interpolating…"};
  }
  const note = state.view === "spaxels" ? "no spaxel lattice in this result; showing bins" : "";
  return {type: "bins", values, note};
}

// ------------------------------------------------------------------ limits & colormaps

function currentCmap() {
  const info = quantityInfo(state.quantity);
  if (isDiffMode()) return state.cmap[`diff:${state.quantity}`] || "RdBu_r";
  return state.cmap[state.quantity] || info.cmap;
}

function defaultLimitsMode() {
  return isDiffMode() ? "symmetric" : quantityInfo(state.quantity).limits;
}

function resetLimitsMode() {
  state.limitsMode = defaultLimitsMode();
}

function resolveLimits(values) {
  if (state.limitsMode === "manual" && Number.isFinite(state.vmin) && Number.isFinite(state.vmax) && state.vmax > state.vmin) {
    return [state.vmin, state.vmax];
  }
  const mode = state.limitsMode === "manual" ? defaultLimitsMode() : (state.limitsMode || defaultLimitsMode());
  const lims = values ? limitsFor(values, mode) : [0, 1];
  state.vmin = lims[0];
  state.vmax = lims[1];
  return lims;
}

// ------------------------------------------------------------------ map geometry

function displayExtent(r) {
  const [c, s] = rotationTrig(state.rotation);
  let xmin = Infinity;
  let xmax = -Infinity;
  let ymin = Infinity;
  let ymax = -Infinity;
  const add = (x, y) => {
    const [xd, yd] = rotatePoint(x, y, c, s);
    if (!Number.isFinite(xd) || !Number.isFinite(yd)) return;
    xmin = Math.min(xmin, xd);
    xmax = Math.max(xmax, xd);
    ymin = Math.min(ymin, yd);
    ymax = Math.max(ymax, yd);
  };
  for (let b = 0; b < r.nbins; b += 1) add(r.x[b], r.y[b]);
  if (r.grid && (spaxelGridUsable(r) || state.quantity === "flux")) gridCorners(r.grid).forEach(([x, y]) => add(x, y));
  if (!Number.isFinite(xmin)) return {xmin: -1, xmax: 1, ymin: -1, ymax: 1};
  const pad = 0.04 * Math.max(xmax - xmin, ymax - ymin, 0.2);
  return {xmin: xmin - pad, xmax: xmax + pad, ymin: ymin - pad, ymax: ymax + pad};
}

function computeView(width, height, bounds) {
  const pad = {left: 58, right: 12, top: 12, bottom: 48};
  const pw = Math.max(50, width - pad.left - pad.right);
  const ph = Math.max(50, height - pad.top - pad.bottom);
  const bw = bounds.xmax - bounds.xmin;
  const bh = bounds.ymax - bounds.ymin;
  const k = Math.min(pw / bw, ph / bh);
  const cx = (bounds.xmin + bounds.xmax) / 2;
  const cy = (bounds.ymin + bounds.ymax) / 2;
  const [cos, sin] = rotationTrig(state.rotation);
  return {
    left: pad.left, top: pad.top, width: pw, height: ph, k, cos, sin,
    xmin: cx - pw / k / 2, xmax: cx + pw / k / 2, ymin: cy - ph / k / 2, ymax: cy + ph / k / 2
  };
}

function rawToScreen(v, x, y) {
  const xr = x * v.cos - y * v.sin;
  const yr = x * v.sin + y * v.cos;
  return [v.left + (xr - v.xmin) * v.k, v.top + (v.ymax - yr) * v.k];
}

function displayToScreen(v, xd, yd) {
  return [v.left + (xd - v.xmin) * v.k, v.top + (v.ymax - yd) * v.k];
}

function screenToDisplay(v, sx, sy) {
  return [v.xmin + (sx - v.left) / v.k, v.ymax - (sy - v.top) / v.k];
}

function screenToRaw(v, sx, sy) {
  const [xd, yd] = screenToDisplay(v, sx, sy);
  return [xd * v.cos + yd * v.sin, -xd * v.sin + yd * v.cos];
}

// ------------------------------------------------------------------ map drawing

const offscreen = document.createElement("canvas");

// "#rrggbb", "#rgb" or "rgb(r, g, b)" -> [r, g, b]
function parseColor(color) {
  const c = String(color).trim();
  if (c.startsWith("#")) {
    const h = c.length === 4 ? c.slice(1).split("").map((d) => d + d).join("") : c.slice(1, 7);
    return [0, 2, 4].map((i) => parseInt(h.slice(i, i + 2), 16));
  }
  const m = c.match(/\d+(\.\d+)?/g);
  return m ? m.slice(0, 3).map(Number) : [128, 128, 128];
}

function drawGridImage(ctx, v, dpr, grid, gridValues, lut, vmin, vmax, smooth, footprintColor) {
  const {nx, ny} = grid;
  offscreen.width = nx;
  offscreen.height = ny;
  const octx = offscreen.getContext("2d");
  const img = octx.createImageData(nx, ny);
  const data = img.data;
  const fp = footprintColor ? parseColor(footprintColor) : null;
  for (let k = 0; k < gridValues.length; k += 1) {
    const val = gridValues[k];
    const o = k * 4;
    if (Number.isFinite(val)) {
      const i = lutIndex(val, vmin, vmax) * 3;
      data[o] = lut[i];
      data[o + 1] = lut[i + 1];
      data[o + 2] = lut[i + 2];
      data[o + 3] = 255;
    } else if (fp && grid.bin && grid.bin[k] >= 0) {
      data[o] = fp[0];
      data[o + 1] = fp[1];
      data[o + 2] = fp[2];
      data[o + 3] = 255;
    }
  }
  octx.putImageData(img, 0, 0);
  ctx.save();
  const k = v.k;
  ctx.setTransform(
    dpr * k * v.cos, dpr * -k * v.sin, dpr * -k * v.sin, dpr * -k * v.cos,
    dpr * (v.left - k * v.xmin), dpr * (v.top + k * v.ymax)
  );
  ctx.transform(grid.dx, 0, 0, grid.dy, grid.x0 - grid.dx / 2, grid.y0 - grid.dy / 2);
  ctx.imageSmoothingEnabled = smooth;
  ctx.drawImage(offscreen, 0, 0);
  ctx.restore();
}

function drawBinMarkers(ctx, v, r, values, lut, vmin, vmax, t) {
  const radius = Math.max(2.2, Math.min(10, 0.48 * r.spacing * v.k));
  ctx.save();
  ctx.lineWidth = radius > 4 ? 0.8 : 0.5;
  ctx.strokeStyle = t.plotBg;
  for (let b = 0; b < r.nbins; b += 1) {
    const val = values[b];
    if (!Number.isFinite(val)) continue;
    const [px, py] = rawToScreen(v, r.x[b], r.y[b]);
    const i = lutIndex(val, vmin, vmax) * 3;
    ctx.fillStyle = `rgb(${lut[i]},${lut[i + 1]},${lut[i + 2]})`;
    ctx.beginPath();
    ctx.arc(px, py, radius, 0, Math.PI * 2);
    ctx.fill();
    if (radius > 3) ctx.stroke();
  }
  ctx.restore();
}

function drawSegmentsRaw(ctx, v, contours, color, width, dash) {
  ctx.save();
  ctx.strokeStyle = color;
  ctx.lineWidth = width;
  ctx.setLineDash(dash);
  ctx.beginPath();
  contours.forEach(({segments}) => {
    for (let i = 0; i < segments.length; i += 4) {
      const [ax, ay] = rawToScreen(v, segments[i], segments[i + 1]);
      const [bx, by] = rawToScreen(v, segments[i + 2], segments[i + 3]);
      ctx.moveTo(ax, ay);
      ctx.lineTo(bx, by);
    }
  });
  ctx.stroke();
  ctx.restore();
}

function isophotes(r) {
  if (r._isophotes !== undefined) return r._isophotes;
  r._isophotes = null;
  if (r.grid?.flux) {
    const levels = magnitudeLevels(r.grid.flux, 0.5, 14);
    r._isophotes = marchingSquares(r.grid, r.grid.flux, levels);
  }
  return r._isophotes;
}

function mgeContours(r) {
  if (!app.mge) return null;
  const e = displayExtent(r);
  const half = Math.max(Math.abs(e.xmin), Math.abs(e.xmax), Math.abs(e.ymin), Math.abs(e.ymax)) * 1.05;
  const key = `${state.mgePa}|${half.toFixed(3)}`;
  if (app.mgeCache?.key === key) return app.mgeCache.contours;
  const {grid, values} = mgeModelGrid(app.mge.gaussians, [-half, half, -half, half], state.mgePa, 260);
  const contours = marchingSquares(grid, values, magnitudeLevels(values, 0.5, 14));
  app.mgeCache = {key, contours};
  return contours;
}

function drawOverlays(ctx, v, r, t) {
  let any = false;
  if (state.showIsophotes) {
    const iso = isophotes(r);
    if (iso?.length) {
      drawSegmentsRaw(ctx, v, iso, "#f59e0b", 1.3, [5, 3]);
      any = true;
    }
  }
  if (state.showMge) {
    const mge = mgeContours(r);
    if (mge?.length) {
      ctx.save();
      ctx.strokeStyle = t.ink;
      ctx.lineWidth = 1.3;
      ctx.beginPath();
      mge.forEach(({segments}) => {
        for (let i = 0; i < segments.length; i += 4) {
          const [ax, ay] = displayToScreen(v, segments[i], segments[i + 1]);
          const [bx, by] = displayToScreen(v, segments[i + 2], segments[i + 3]);
          ctx.moveTo(ax, ay);
          ctx.lineTo(bx, by);
        }
      });
      ctx.stroke();
      ctx.restore();
      any = true;
    } else {
      const legacy = r.staticOverlays.find((o) => o.kind === "mge");
      if (legacy) {
        ctx.save();
        ctx.strokeStyle = t.ink;
        ctx.lineWidth = 1.3;
        legacy.contours.forEach((c) => c.polylines.forEach((line) => {
          ctx.beginPath();
          line.forEach(([x, y], i) => {
            const [px, py] = rawToScreen(v, x, y);
            if (i === 0) ctx.moveTo(px, py);
            else ctx.lineTo(px, py);
          });
          ctx.stroke();
        }));
        ctx.restore();
        any = true;
      }
    }
  }
  if (state.showCenters) {
    ctx.save();
    ctx.fillStyle = t.ink;
    ctx.globalAlpha = 0.45;
    const mask = binMask(r);
    for (let b = 0; b < r.nbins; b += 1) {
      if (mask[b]) continue;
      const [px, py] = rawToScreen(v, r.x[b], r.y[b]);
      ctx.fillRect(px - 1, py - 1, 2, 2);
    }
    ctx.restore();
  }
  if (any) {
    const [px, py] = rawToScreen(v, 0, 0);
    ctx.save();
    ctx.strokeStyle = t.ink;
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(px - 7, py);
    ctx.lineTo(px + 7, py);
    ctx.moveTo(px, py - 7);
    ctx.lineTo(px, py + 7);
    ctx.stroke();
    ctx.restore();
  }
}

function drawSelectionRings(ctx, v, r) {
  state.selected.forEach((b, idx) => {
    if (b >= r.nbins) return;
    const [px, py] = rawToScreen(v, r.x[b], r.y[b]);
    const activeBin = b === state.active;
    ctx.save();
    ctx.lineWidth = activeBin ? 4.5 : 3.5;
    ctx.strokeStyle = "rgba(255,255,255,0.85)";
    ctx.beginPath();
    ctx.arc(px, py, activeBin ? 9 : 7, 0, Math.PI * 2);
    ctx.stroke();
    ctx.lineWidth = activeBin ? 2.6 : 1.8;
    ctx.strokeStyle = BIN_COLORS[idx % BIN_COLORS.length];
    ctx.stroke();
    ctx.restore();
  });
}

function renderMap() {
  const r = activeResult();
  if (!r) return;
  const t = theme();
  const {ctx, width, height, dpr} = setupCanvas(el.mapCanvas);
  ctx.clearRect(0, 0, width, height);
  const base = displayExtent(r);
  const v = computeView(width, height, state.zoom || base);
  state.mapView = v;
  ctx.fillStyle = t.plotBg;
  ctx.fillRect(v.left, v.top, v.width, v.height);

  const info = quantityInfo(state.quantity);
  const diff = isDiffMode();
  const layer = mapLayer(r);
  // Limits follow the per-bin values so they do not jump between display modes;
  // the spaxel flux image is the exception (it has no per-bin equivalent).
  const limitSource = layer.note === "spaxel flux" ? layer.gridValues : layer.values;
  const [vmin, vmax] = resolveLimits(limitSource);
  syncLimitInputs();
  const lut = getLut(currentCmap(), state.reverse);

  ctx.save();
  ctx.beginPath();
  ctx.rect(v.left, v.top, v.width, v.height);
  ctx.clip();
  if (layer.type === "grid") {
    drawGridImage(ctx, v, dpr, layer.grid, layer.gridValues, lut, vmin, vmax, layer.smooth, layer.footprint ? t.empty : null);
    // the transform above is reset by restore(); re-apply the base dpr scale for vector drawing
  } else if (layer.type === "bins") {
    drawBinMarkers(ctx, v, r, layer.values, lut, vmin, vmax, t);
  }
  drawOverlays(ctx, v, r, t);
  drawSelectionRings(ctx, v, r);
  ctx.restore();

  if (layer.type === "none") {
    ctx.fillStyle = t.muted;
    ctx.font = `14px ${t.font}`;
    ctx.textAlign = "center";
    ctx.fillText(`${info.name} is not available for this result`, v.left + v.width / 2, v.top + v.height / 2);
  }
  const rotated = Math.abs(state.rotation) > 1e-9;
  drawAxes(ctx, {left: v.left, top: v.top, width: v.width, height: v.height, xmin: v.xmin, xmax: v.xmax, ymin: v.ymin, ymax: v.ymax}, t, {
    xlabel: rotated ? "x′ (arcsec)" : "x (arcsec)",
    ylabel: rotated ? "y′ (arcsec)" : "y (arcsec)"
  });
  const unit = info.unit ? (diff ? `Δ ${info.unit}` : info.unit) : (diff ? `Δ ${info.label}` : info.label);
  drawColorbar(el.colorbarCanvas, lut, vmin, vmax, unit, t, {top: v.top, bottom: height - v.top - v.height});

  // titles & status
  const ref = compareResult();
  el.mapTitle.textContent = diff ? `${info.name}: difference` : info.name;
  el.mapSubtitle.textContent = diff ? `${r.label} − ${ref.label}` : r.label;
  const parts = [];
  if (layer.note) parts.push(layer.note);
  if (ref && !diff) parts.push(`${ref.label} has no ${info.name}; showing ${r.label} only`);
  if (diff) {
    const match = matchFor(r, ref);
    const matched = match ? match.reduce((n, j) => n + (j >= 0 ? 1 : 0), 0) : r.nbins;
    parts.push(`${matched}/${r.nbins} bins matched`);
  }
  if (state.maskSigma && r._maskInfo) parts.push(`${r._maskInfo.count} σ-outlier bins hidden (σ > ${formatValue(r._maskInfo.cmax, "km/s")})`);
  if (Number.isFinite(state.minSnr)) parts.push(`S/N ≥ ${state.minSnr}`);
  if (state.zoom) parts.push(`zoom ${formatValue(v.xmax - v.xmin)}″ × ${formatValue(v.ymax - v.ymin)}″`);
  if (state.showMge && !app.mge && !r.staticOverlays.length) parts.push("load an MGE table to draw MGE contours");
  el.mapStatus.textContent = parts.join(" · ");
}

// ------------------------------------------------------------------ hit testing & map interaction

function binAtScreen(sx, sy) {
  const r = activeResult();
  const v = state.mapView;
  if (!r || !v) return null;
  if (sx < v.left || sx > v.left + v.width || sy < v.top || sy > v.top + v.height) return null;
  const [x, y] = screenToRaw(v, sx, sy);
  const mask = binMask(r);
  if (spaxelGridUsable(r) && state.view === "spaxels") {
    const g = r.grid;
    const ix = Math.round((x - g.x0) / g.dx);
    const iy = Math.round((y - g.y0) / g.dy);
    if (ix >= 0 && iy >= 0 && ix < g.nx && iy < g.ny) {
      const b = g.bin[iy * g.nx + ix];
      if (b >= 0 && b < r.nbins) return mask[b] ? null : b;
    }
  }
  let best = null;
  let bestD2 = (Math.max(2 * r.spacing, 10 / v.k)) ** 2;
  for (let b = 0; b < r.nbins; b += 1) {
    if (mask[b]) continue;
    const d2 = (r.x[b] - x) ** 2 + (r.y[b] - y) ** 2;
    if (d2 < bestD2) {
      bestD2 = d2;
      best = b;
    }
  }
  return best;
}

function canvasPoint(event, canvas) {
  const rect = canvas.getBoundingClientRect();
  return [event.clientX - rect.left, event.clientY - rect.top];
}

function zoomMap(factor, anchor) {
  const r = activeResult();
  const v = state.mapView;
  if (!r || !v) return;
  const cur = {xmin: v.xmin, xmax: v.xmax, ymin: v.ymin, ymax: v.ymax};
  const [ax, ay] = anchor || [(cur.xmin + cur.xmax) / 2, (cur.ymin + cur.ymax) / 2];
  const next = {
    xmin: ax - (ax - cur.xmin) * factor,
    xmax: ax + (cur.xmax - ax) * factor,
    ymin: ay - (ay - cur.ymin) * factor,
    ymax: ay + (cur.ymax - ay) * factor
  };
  const base = displayExtent(r);
  const baseSpan = Math.max(base.xmax - base.xmin, base.ymax - base.ymin);
  const span = Math.max(next.xmax - next.xmin, next.ymax - next.ymin);
  if (span >= baseSpan * 1.02) state.zoom = null;
  else if (span > baseSpan / 400) state.zoom = next;
  renderMap();
}

function updateTooltip(event) {
  const [sx, sy] = canvasPoint(event, el.mapCanvas);
  const b = binAtScreen(sx, sy);
  const r = activeResult();
  if (b === null || !r) {
    el.mapTooltip.hidden = true;
    return;
  }
  const info = quantityInfo(state.quantity);
  const values = binValues(r);
  const val = values ? values[b] : NaN;
  const prefix = isDiffMode() ? "Δ" : "";
  let text = `Bin ${b} · ${prefix}${info.label} ${formatValue(val, info.unit)}`;
  if (state.quantity === "flux" && r.grid?.flux && state.view !== "bins") {
    text = `Bin ${b} · log flux ${formatValue(val)}`;
  }
  el.mapTooltip.textContent = text;
  el.mapTooltip.hidden = false;
  const wrap = el.mapCanvas.parentElement.getBoundingClientRect();
  const tipW = el.mapTooltip.offsetWidth;
  const left = Math.min(sx + 14, wrap.width - tipW - 4);
  el.mapTooltip.style.left = `${Math.max(4, left)}px`;
  el.mapTooltip.style.top = `${Math.max(4, sy - 30)}px`;
}

function bindMapInteraction() {
  const c = el.mapCanvas;
  c.addEventListener("wheel", (event) => {
    const v = state.mapView;
    if (!v) return;
    const [sx, sy] = canvasPoint(event, c);
    if (sx < v.left || sx > v.left + v.width || sy < v.top || sy > v.top + v.height) return;
    event.preventDefault();
    zoomMap(event.deltaY < 0 ? 0.8 : 1.25, screenToDisplay(v, sx, sy));
  }, {passive: false});
  c.addEventListener("pointerdown", (event) => {
    const [sx, sy] = canvasPoint(event, c);
    const v = state.mapView;
    if (!v) return;
    state.drag = {id: event.pointerId, sx, sy, view: {...v}, moved: false};
    c.setPointerCapture(event.pointerId);
  });
  c.addEventListener("pointermove", (event) => {
    if (!state.drag || state.drag.id !== event.pointerId) {
      updateTooltip(event);
      return;
    }
    const [sx, sy] = canvasPoint(event, c);
    const dx = sx - state.drag.sx;
    const dy = sy - state.drag.sy;
    if (!state.drag.moved && Math.hypot(dx, dy) < 4) return;
    state.drag.moved = true;
    c.classList.add("panning");
    el.mapTooltip.hidden = true;
    const v0 = state.drag.view;
    state.zoom = {
      xmin: v0.xmin - dx / v0.k, xmax: v0.xmax - dx / v0.k,
      ymin: v0.ymin + dy / v0.k, ymax: v0.ymax + dy / v0.k
    };
    renderMap();
  });
  const end = (event) => {
    if (!state.drag || state.drag.id !== event.pointerId) return;
    const moved = state.drag.moved;
    state.drag = null;
    c.classList.remove("panning");
    if (!moved && event.type === "pointerup") {
      const [sx, sy] = canvasPoint(event, c);
      const b = binAtScreen(sx, sy);
      if (b !== null) selectBin(b);
    }
  };
  c.addEventListener("pointerup", end);
  c.addEventListener("pointercancel", end);
  c.addEventListener("pointerleave", () => {
    el.mapTooltip.hidden = true;
  });
  c.addEventListener("dblclick", () => {
    state.zoom = null;
    renderMap();
  });
}

// ------------------------------------------------------------------ selection

function selectBin(b, {replace = false} = {}) {
  const r = activeResult();
  if (!r || !Number.isInteger(b) || b < 0 || b >= r.nbins) return;
  if (replace) state.selected = [b];
  else if (!state.selected.includes(b)) {
    state.selected.push(b);
    if (state.selected.length > MAX_SELECTED) state.selected.shift();
  }
  state.active = b;
  el.binInput.value = String(b);
  state.specRange = null;
  renderAll();
}

function selectDefaultBin() {
  const r = activeResult();
  if (!r) return;
  // the bin closest to the centre that has output
  let best = null;
  let bestR = Infinity;
  for (let b = 0; b < r.nbins; b += 1) {
    if (r.hasOutput[b] && r.r[b] < bestR) {
      bestR = r.r[b];
      best = b;
    }
  }
  if (best !== null) selectBin(best, {replace: true});
}

function removeSelected(b) {
  state.selected = state.selected.filter((x) => x !== b);
  if (state.active === b) state.active = state.selected.length ? state.selected[state.selected.length - 1] : null;
  renderAll();
}

function setActive(uid, {keepSelection = true} = {}) {
  const prev = activeResult();
  state.activeUid = uid;
  const r = activeResult();
  if (!r) return;
  if (keepSelection && prev && prev !== r) {
    const match = matchBins(prev, r, 0.5 * Math.min(prev.spacing, r.spacing));
    const mapped = state.selected.map((b) => match[b]).filter((b) => b >= 0);
    state.active = state.active !== null && match[state.active] >= 0 ? match[state.active] : (mapped[mapped.length - 1] ?? null);
    state.selected = Array.from(new Set(mapped));
  } else if (!keepSelection) {
    state.selected = [];
    state.active = null;
  }
  if (!availableQuantities(r).includes(state.quantity)) state.quantity = availableQuantities(r)[0] || "vel";
  if (state.compareUid === uid) state.compareUid = "";
  state.zoom = null;
  state.specRange = null;
  resetLimitsMode();
  renderAll();
}

// ------------------------------------------------------------------ side panels

function binMetrics(r, b) {
  const rows = [];
  const add = (name, value, title = "") => rows.push({name, value, title});
  add("Bin", String(b));
  add("Radius", formatValue(r.r[b], "″"));
  add("S/N", formatValue(r.snr[b]));
  if (r.losvd) {
    const q = (k) => r.quantities.get(k)?.[b];
    add("V", formatValue(q("vel"), "km/s"), QUANTITIES.vel.desc);
    add("σ", formatValue(q("sigma"), "km/s"), QUANTITIES.sigma.desc);
    add("Skew", formatValue(q("skew")), QUANTITIES.skew.desc);
    add("Kurt", formatValue(q("kurt")), QUANTITIES.kurt.desc);
    if (r.hasOutput[b]) {
      const gh = fitGaussHermite(r.xvel, r.losvd, (b * r.rows + r.medianRow) * r.nvel);
      add("GH V", formatValue(gh.vel, "km/s"), QUANTITIES.gh_vel.desc);
      add("GH σ", formatValue(gh.sigma, "km/s"), QUANTITIES.gh_vel.desc);
      add("GH h3", formatValue(gh.h3), QUANTITIES.gh_vel.desc);
      add("GH h4", formatValue(gh.h4), QUANTITIES.gh_vel.desc);
    }
  }
  const moments = ["vel", "sigma", "skew", "kurt"];
  r.quantities.forEach((values, key) => {
    const extra = key.startsWith("fit_") || key.startsWith("col:") || key.startsWith("table_");
    if (!extra && !(moments.includes(key) && !r.losvd)) return;
    const info = quantityInfo(key);
    let text = formatValue(values[b], info.unit);
    const err = r.paramErrors.get(key);
    if (err && Number.isFinite(err.lo[b]) && Number.isFinite(err.hi[b])) {
      text = `${formatValue(values[b])} (+${formatValue(err.hi[b] - values[b])}/−${formatValue(values[b] - err.lo[b])})${info.unit ? ` ${info.unit}` : ""}`;
    }
    add(info.group === "Model parameters" ? `Fit ${info.label}` : info.name, text, info.desc);
  });
  return rows;
}

function renderBinPanel() {
  const r = activeResult();
  el.binChips.innerHTML = "";
  state.selected.forEach((b, idx) => {
    const chip = document.createElement("span");
    chip.className = `chip${b === state.active ? " active" : ""}`;
    chip.style.setProperty("--chip-color", BIN_COLORS[idx % BIN_COLORS.length]);
    const label = document.createElement("button");
    label.type = "button";
    label.className = "chip-label";
    label.textContent = `Bin ${b}`;
    label.addEventListener("click", () => {
      state.active = b;
      state.specRange = null;
      renderAll();
    });
    const close = document.createElement("button");
    close.type = "button";
    close.className = "chip-close";
    close.setAttribute("aria-label", `Remove bin ${b}`);
    close.textContent = "×";
    close.addEventListener("click", () => removeSelected(b));
    chip.append(label, close);
    el.binChips.appendChild(chip);
  });
  el.metrics.innerHTML = "";
  if (!r || state.active === null || state.active >= r.nbins) {
    el.metrics.innerHTML = '<p class="hint">Click the map or enter a bin number.</p>';
    return;
  }
  if (!r.hasOutput[state.active]) {
    el.metrics.innerHTML = `<p class="hint">Bin ${state.active} has no output in this result.</p>`;
    return;
  }
  binMetrics(r, state.active).forEach(({name, value, title}) => {
    const div = document.createElement("div");
    div.className = "metric";
    if (title) div.title = title;
    const n = document.createElement("div");
    n.className = "metric-name";
    n.textContent = name;
    const val = document.createElement("div");
    val.className = "metric-value";
    val.textContent = value;
    div.append(n, val);
    el.metrics.appendChild(div);
  });
}

function messageOnCanvas(canvas, text) {
  const t = theme();
  const {ctx, width, height} = setupCanvas(canvas);
  ctx.clearRect(0, 0, width, height);
  ctx.fillStyle = t.muted;
  ctx.font = `13px ${t.font}`;
  ctx.textAlign = "center";
  ctx.textBaseline = "middle";
  ctx.fillText(text, width / 2, height / 2);
}

function losvdSeries() {
  const r = activeResult();
  if (!r) return [];
  if (state.losvdAll && state.active !== null) {
    const series = [];
    app.results.forEach((o) => {
      if (!o.losvd) return;
      const b = o === r ? state.active : matchBins({x: [r.x[state.active]], y: [r.y[state.active]]}, o, 0.5 * Math.min(r.spacing, o.spacing))[0];
      if (b >= 0 && o.hasOutput[b]) series.push({r: o, b, color: o.color, label: `${o.label} (bin ${b})`, bold: o === r});
    });
    return series;
  }
  if (!r.losvd) return [];
  return state.selected
    .map((b, idx) => ({r, b, color: BIN_COLORS[idx % BIN_COLORS.length], label: `Bin ${b}`, bold: b === state.active}))
    .filter((s) => s.b < r.nbins && r.hasOutput[s.b]);
}

function renderLosvd() {
  const r = activeResult();
  el.losvdLegend.innerHTML = "";
  if (!r) return;
  if (!r.losvd && !state.losvdAll) {
    messageOnCanvas(el.losvdCanvas, "This result has no LOSVDs");
    return;
  }
  const series = losvdSeries();
  if (!series.length) {
    messageOnCanvas(el.losvdCanvas, "No bins selected");
    return;
  }
  const t = theme();
  const {ctx, width, height} = setupCanvas(el.losvdCanvas);
  ctx.clearRect(0, 0, width, height);
  let xmin = Infinity;
  let xmax = -Infinity;
  let ymax = 0;
  series.forEach(({r: o, b}) => {
    xmin = Math.min(xmin, o.xvel[0]);
    xmax = Math.max(xmax, o.xvel[o.nvel - 1]);
    const top = o.bands ? o.rows - 1 : o.medianRow;
    for (let v = 0; v < o.nvel; v += 1) {
      const y = o.losvd[(b * o.rows + top) * o.nvel + v];
      if (Number.isFinite(y)) ymax = Math.max(ymax, y);
    }
  });
  if (!(ymax > 0)) ymax = 1;
  const plot = {left: 58, top: 10, width: width - 70, height: height - 58, xmin, xmax, ymin: 0, ymax: ymax * 1.08};
  ctx.fillStyle = t.plotBg;
  ctx.fillRect(plot.left, plot.top, plot.width, plot.height);
  const xs = (val) => plot.left + (val - plot.xmin) / (plot.xmax - plot.xmin) * plot.width;
  const ys = (val) => plot.top + (plot.ymax - val) / (plot.ymax - plot.ymin) * plot.height;
  ctx.save();
  ctx.beginPath();
  ctx.rect(plot.left, plot.top, plot.width, plot.height);
  ctx.clip();
  ctx.strokeStyle = t.grid;
  ctx.setLineDash([3, 3]);
  ctx.beginPath();
  ctx.moveTo(xs(0), plot.top);
  ctx.lineTo(xs(0), plot.top + plot.height);
  ctx.stroke();
  ctx.setLineDash([]);
  const band = (o, b, lo, hi, alpha, color) => {
    ctx.fillStyle = color;
    ctx.globalAlpha = alpha;
    ctx.beginPath();
    for (let v = 0; v < o.nvel; v += 1) {
      const y = o.losvd[(b * o.rows + hi) * o.nvel + v];
      if (v === 0) ctx.moveTo(xs(o.xvel[v]), ys(y));
      else ctx.lineTo(xs(o.xvel[v]), ys(y));
    }
    for (let v = o.nvel - 1; v >= 0; v -= 1) ctx.lineTo(xs(o.xvel[v]), ys(o.losvd[(b * o.rows + lo) * o.nvel + v]));
    ctx.closePath();
    ctx.fill();
    ctx.globalAlpha = 1;
  };
  series.forEach(({r: o, b, color, bold}) => {
    if (o.bands && o.rows >= 5 && (bold || series.length === 1)) {
      band(o, b, 0, 4, 0.12, color);
      band(o, b, 1, 3, 0.25, color);
    }
  });
  series.forEach(({r: o, b, color, bold}) => {
    ctx.strokeStyle = color;
    ctx.lineWidth = bold ? 2.4 : 1.5;
    ctx.beginPath();
    let started = false;
    for (let v = 0; v < o.nvel; v += 1) {
      const y = o.losvd[(b * o.rows + o.medianRow) * o.nvel + v];
      if (!Number.isFinite(y)) continue;
      if (!started) ctx.moveTo(xs(o.xvel[v]), ys(y));
      else ctx.lineTo(xs(o.xvel[v]), ys(y));
      started = true;
    }
    ctx.stroke();
  });
  ctx.restore();
  drawAxes(ctx, plot, t, {xlabel: "velocity (km/s)", ylabel: "LOSVD", yTicks: 4});
  series.forEach(({color, label}) => {
    const item = document.createElement("span");
    item.className = "legend-item";
    item.style.setProperty("--swatch", color);
    item.textContent = label;
    el.losvdLegend.appendChild(item);
  });
  const anyBands = series.some(({r: o}) => o.bands);
  el.losvdHint.textContent = anyBands
    ? `Median LOSVD; shaded: ${PERCENTILES[1]}–${PERCENTILES[3]} and ${PERCENTILES[0]}–${PERCENTILES[4]} percentiles`
    : "Median LOSVD (MAP run: no credible bands stored)";
}

// ------------------------------------------------------------------ spectrum

function requestSpectrum(r, b) {
  if (!r.spectra) return null;
  r._spectra = r._spectra || new Map();
  if (r._spectra.has(b)) return r._spectra.get(b);
  r._spectra.set(b, "loading");
  client().spectrum(r.spectra.fileId, b).then((data) => {
    r._spectra.set(b, data);
    if (activeResult() === r && state.active === b) renderSpectrum();
  }).catch((error) => {
    r._spectra.set(b, {error: error.message});
    if (activeResult() === r && state.active === b) renderSpectrum();
  });
  return "loading";
}

function renderSpectrum() {
  const r = activeResult();
  el.spectrumStatus.textContent = "";
  state.specView = null;
  if (!r) return;
  if (!r.spectra) {
    messageOnCanvas(el.spectrumCanvas, "Spectral fits are available for results loaded from HDF5 files");
    return;
  }
  const b = state.active;
  if (b === null || b >= r.nbins) {
    messageOnCanvas(el.spectrumCanvas, "Select a bin to see its spectral fit");
    return;
  }
  if (!r.hasOutput[b]) {
    messageOnCanvas(el.spectrumCanvas, `Bin ${b} has no fit in this result`);
    return;
  }
  const data = requestSpectrum(r, b);
  if (data === "loading") {
    messageOnCanvas(el.spectrumCanvas, `Loading spectrum of bin ${b}…`);
    return;
  }
  if (data.error) {
    messageOnCanvas(el.spectrumCanvas, `Could not read spectrum: ${data.error}`);
    return;
  }
  const {wave, obs, err, bestfit, bestfitRows, mask} = data;
  if (!wave || !obs) {
    messageOnCanvas(el.spectrumCanvas, "No observed spectrum stored for this bin");
    return;
  }
  const n = wave.length;
  const fitRows = bestfitRows || 1;
  const medRow = fitRows >= 5 ? 2 : 0;
  const fit = (row, i) => (bestfit ? bestfit[row * n + i] : NaN);
  const fitted = new Uint8Array(n);
  if (mask && mask.length) mask.forEach((i) => {
    if (i >= 0 && i < n) fitted[i] = 1;
  });
  else fitted.fill(1);

  const full = [wave[0], wave[n - 1]];
  const [x0, x1] = state.specRange || full;
  let ylo = Infinity;
  let yhi = -Infinity;
  let rlo = 0;
  let rhi = 0;
  let chi2 = 0;
  let nchi = 0;
  for (let i = 0; i < n; i += 1) {
    const f = fit(medRow, i);
    if (fitted[i] && Number.isFinite(obs[i]) && Number.isFinite(f) && err && err[i] > 0) {
      chi2 += ((obs[i] - f) / err[i]) ** 2;
      nchi += 1;
    }
    if (wave[i] < x0 || wave[i] > x1) continue;
    [obs[i], f].forEach((val) => {
      if (Number.isFinite(val)) {
        ylo = Math.min(ylo, val);
        yhi = Math.max(yhi, val);
      }
    });
    if (fitted[i] && Number.isFinite(obs[i]) && Number.isFinite(f)) {
      rlo = Math.min(rlo, obs[i] - f);
      rhi = Math.max(rhi, obs[i] - f);
    }
  }
  if (!Number.isFinite(ylo)) {
    ylo = 0;
    yhi = 1;
  }
  const ypad = 0.06 * (yhi - ylo || 1);
  const rabs = Math.max(Math.abs(rlo), Math.abs(rhi)) * 1.15 || 0.01;

  const t = theme();
  const {ctx, width, height} = setupCanvas(el.spectrumCanvas);
  ctx.clearRect(0, 0, width, height);
  const left = 58;
  const pw = width - left - 14;
  const mainH = Math.round((height - 58) * 0.7);
  const resH = height - 58 - mainH - 8;
  const top = {left, top: 10, width: pw, height: mainH, xmin: x0, xmax: x1, ymin: ylo - ypad, ymax: yhi + ypad};
  const bottom = {left, top: 10 + mainH + 8, width: pw, height: resH, xmin: x0, xmax: x1, ymin: -rabs, ymax: rabs};
  state.specView = {left, width: pw, x0, x1, full};
  const xs = (w) => left + (w - x0) / (x1 - x0) * pw;

  [top, bottom].forEach((p) => {
    ctx.fillStyle = t.plotBg;
    ctx.fillRect(p.left, p.top, p.width, p.height);
  });
  // shade pixels excluded from the fit
  const shade = (p) => {
    ctx.save();
    ctx.fillStyle = t.empty;
    let i = 0;
    while (i < n) {
      if (fitted[i]) {
        i += 1;
        continue;
      }
      const start = i;
      while (i < n && !fitted[i]) i += 1;
      const wa = start > 0 ? 0.5 * (wave[start - 1] + wave[start]) : wave[start];
      const wb = i < n ? 0.5 * (wave[i - 1] + wave[i]) : wave[n - 1];
      const a = Math.max(p.left, xs(wa));
      const c = Math.min(p.left + p.width, xs(wb));
      if (c > a) ctx.fillRect(a, p.top, c - a, p.height);
    }
    ctx.restore();
  };
  shade(top);
  shade(bottom);

  const line = (p, getY, color, widthPx, onlyFitted = false) => {
    const ys = (val) => p.top + (p.ymax - val) / (p.ymax - p.ymin) * p.height;
    ctx.save();
    ctx.beginPath();
    ctx.rect(p.left, p.top, p.width, p.height);
    ctx.clip();
    ctx.strokeStyle = color;
    ctx.lineWidth = widthPx;
    ctx.beginPath();
    let pen = false;
    for (let i = 0; i < n; i += 1) {
      const val = getY(i);
      if (!Number.isFinite(val) || (onlyFitted && !fitted[i])) {
        pen = false;
        continue;
      }
      if (!pen) ctx.moveTo(xs(wave[i]), ys(val));
      else ctx.lineTo(xs(wave[i]), ys(val));
      pen = true;
    }
    ctx.stroke();
    ctx.restore();
  };
  if (fitRows >= 5) {
    // credible band of the best-fit model
    const ys = (val) => top.top + (top.ymax - val) / (top.ymax - top.ymin) * top.height;
    ctx.save();
    ctx.beginPath();
    ctx.rect(top.left, top.top, top.width, top.height);
    ctx.clip();
    ctx.fillStyle = t.accent;
    ctx.globalAlpha = 0.25;
    ctx.beginPath();
    for (let i = 0; i < n; i += 1) {
      const val = fit(3, i);
      if (i === 0) ctx.moveTo(xs(wave[i]), ys(val));
      else ctx.lineTo(xs(wave[i]), ys(val));
    }
    for (let i = n - 1; i >= 0; i -= 1) ctx.lineTo(xs(wave[i]), ys(fit(1, i)));
    ctx.closePath();
    ctx.fill();
    ctx.restore();
  }
  line(top, (i) => obs[i], t.ink, 1);
  line(top, (i) => fit(medRow, i), t.accent, 1.5);
  // residuals with ±1σ error band
  if (err) {
    const ys = (val) => bottom.top + (bottom.ymax - val) / (bottom.ymax - bottom.ymin) * bottom.height;
    ctx.save();
    ctx.beginPath();
    ctx.rect(bottom.left, bottom.top, bottom.width, bottom.height);
    ctx.clip();
    ctx.fillStyle = t.grid;
    ctx.beginPath();
    for (let i = 0; i < n; i += 1) ctx.lineTo(xs(wave[i]), ys(Number.isFinite(err[i]) ? err[i] : 0));
    for (let i = n - 1; i >= 0; i -= 1) ctx.lineTo(xs(wave[i]), ys(Number.isFinite(err[i]) ? -err[i] : 0));
    ctx.closePath();
    ctx.fill();
    ctx.restore();
  }
  line(bottom, () => 0, t.axis, 1);
  line(bottom, (i) => obs[i] - fit(medRow, i), t.accent, 1, true);
  drawAxes(ctx, top, t, {ylabel: "flux", xTicks: 8, xFormat: () => ""});
  drawAxes(ctx, bottom, t, {xlabel: "wavelength (Å)", ylabel: "resid.", xTicks: 8, yTicks: 2});
  el.spectrumStatus.textContent = `Bin ${b}${nchi ? ` · χ²/N = ${(chi2 / nchi).toFixed(2)} over ${nchi} fitted pixels` : ""}`;
}

function bindSpectrumInteraction() {
  const c = el.spectrumCanvas;
  const toWave = (sx) => {
    const v = state.specView;
    return v.x0 + (sx - v.left) / v.width * (v.x1 - v.x0);
  };
  const setRange = (lo, hi) => {
    const v = state.specView;
    const [f0, f1] = v.full;
    const span = Math.min(hi - lo, f1 - f0);
    if (span >= (f1 - f0) * 0.999) {
      state.specRange = null;
    } else {
      let a = lo;
      let b = lo + span;
      if (a < f0) {
        a = f0;
        b = f0 + span;
      }
      if (b > f1) {
        b = f1;
        a = f1 - span;
      }
      state.specRange = [a, b];
    }
    renderSpectrum();
  };
  c.addEventListener("wheel", (event) => {
    const v = state.specView;
    if (!v) return;
    const [sx] = canvasPoint(event, c);
    if (sx < v.left || sx > v.left + v.width) return;
    event.preventDefault();
    const w = toWave(sx);
    const f = event.deltaY < 0 ? 0.8 : 1.25;
    setRange(w - (w - v.x0) * f, w + (v.x1 - w) * f);
  }, {passive: false});
  c.addEventListener("pointerdown", (event) => {
    if (!state.specView) return;
    const [sx] = canvasPoint(event, c);
    state.specDrag = {id: event.pointerId, sx, x0: state.specView.x0, x1: state.specView.x1};
    c.setPointerCapture(event.pointerId);
  });
  c.addEventListener("pointermove", (event) => {
    const d = state.specDrag;
    if (!d || d.id !== event.pointerId || !state.specView) return;
    const [sx] = canvasPoint(event, c);
    const dw = (sx - d.sx) / state.specView.width * (d.x1 - d.x0);
    setRange(d.x0 - dw, d.x1 - dw);
  });
  const end = () => {
    state.specDrag = null;
  };
  c.addEventListener("pointerup", end);
  c.addEventListener("pointercancel", end);
  c.addEventListener("dblclick", () => {
    state.specRange = null;
    renderSpectrum();
  });
}

// ------------------------------------------------------------------ sidebar rendering

function renderResultList() {
  el.resultList.innerHTML = "";
  app.results.forEach((r) => {
    const li = document.createElement("li");
    li.className = `result-item${r.uid === state.activeUid ? " active" : ""}`;
    li.style.setProperty("--result-color", r.color);
    const main = document.createElement("button");
    main.type = "button";
    main.className = "result-main";
    main.title = r.source;
    const name = document.createElement("span");
    name.className = "result-name";
    name.textContent = r.label;
    const meta = document.createElement("span");
    meta.className = "result-meta";
    const bits = [r.fitType || r.kind, `${r.nbins} bins`];
    if (r.regularized) bits.push("regularized");
    if (r.bands) bits.push("credible bands");
    meta.textContent = bits.join(" · ");
    main.append(name, meta);
    main.addEventListener("click", () => setActive(r.uid));
    main.addEventListener("dblclick", () => {
      const label = window.prompt("Rename result", r.label);
      if (label && label.trim()) {
        r.label = label.trim();
        renderAll();
      }
    });
    const remove = document.createElement("button");
    remove.type = "button";
    remove.className = "icon-btn";
    remove.title = "Remove result";
    remove.setAttribute("aria-label", `Remove ${r.label}`);
    remove.textContent = "×";
    remove.addEventListener("click", () => removeResult(r.uid));
    li.append(main, remove);
    el.resultList.appendChild(li);
  });

  el.compareSelect.innerHTML = "";
  const none = document.createElement("option");
  none.value = "";
  none.textContent = "None";
  el.compareSelect.appendChild(none);
  app.results.forEach((r) => {
    if (r.uid === state.activeUid) return;
    const o = document.createElement("option");
    o.value = r.uid;
    o.textContent = r.label;
    el.compareSelect.appendChild(o);
  });
  el.compareSelect.value = state.compareUid;
  el.compareSelect.disabled = app.results.length < 2;
  el.compareHint.textContent = app.results.length < 2
    ? "Load a second result to show difference maps."
    : (state.compareUid ? "Map shows active − reference; bins matched by position." : "");
}

function renderQuantityButtons() {
  const r = activeResult();
  el.quantityButtons.innerHTML = "";
  if (!r) return;
  const keys = availableQuantities(r);
  let group = null;
  let row = null;
  keys.forEach((key) => {
    const info = quantityInfo(key);
    if (info.group !== group) {
      group = info.group;
      const g = document.createElement("div");
      g.className = "q-group-label";
      g.textContent = group;
      el.quantityButtons.appendChild(g);
      row = document.createElement("div");
      row.className = "segmented wrap";
      el.quantityButtons.appendChild(row);
    }
    const btn = document.createElement("button");
    btn.type = "button";
    btn.textContent = info.label;
    btn.title = `${info.name}${info.desc ? ` — ${info.desc}` : ""}`;
    btn.classList.toggle("active", key === state.quantity);
    btn.setAttribute("aria-pressed", String(key === state.quantity));
    btn.addEventListener("click", () => {
      state.quantity = key;
      resetLimitsMode();
      if (key.startsWith("gh_") && !r._gh) {
        setStatus("Fitting Gauss–Hermite profiles…", {sticky: true});
        setTimeout(() => {
          ensureGh(r);
          setStatus("");
          renderAll();
        }, 20);
        return;
      }
      renderAll();
    });
    row.appendChild(btn);
  });
}

function syncControls() {
  el.viewButtons.querySelectorAll("button").forEach((b) => {
    b.classList.toggle("active", b.dataset.view === state.view);
    b.setAttribute("aria-pressed", String(b.dataset.view === state.view));
  });
  el.cmapSelect.value = currentCmap();
  el.cmapReverse.checked = state.reverse;
  el.limitButtons.querySelectorAll("button").forEach((b) => b.classList.toggle("active", b.dataset.mode === state.limitsMode));
  el.rotationInput.value = String(state.rotation);
  el.maskSigma.checked = state.maskSigma;
  el.maskRadius.value = String(state.maskRadius);
  el.showIsophotes.checked = state.showIsophotes;
  el.showCenters.checked = state.showCenters;
  const r = activeResult();
  const hasMge = Boolean(app.mge) || Boolean(r?.staticOverlays.some((o) => o.kind === "mge"));
  el.showMge.disabled = !hasMge;
  el.showMge.checked = state.showMge && hasMge;
  el.showIsophotes.disabled = !r?.grid?.flux;
  el.mgePa.value = String(state.mgePa);
  el.mgeStatus.textContent = app.mge
    ? `${app.mge.source}: ${app.mge.gaussians.length} Gaussians`
    : (hasMge ? "Using MGE contours stored in the dashboard file." : "Columns: sigma_arcsec, q_obs, luminosity (or total_flux).");
  el.losvdAll.checked = state.losvdAll;
  el.losvdAll.disabled = app.results.filter((o) => o.losvd).length < 2;
}

function syncLimitInputs() {
  if (document.activeElement !== el.vminInput) el.vminInput.value = Number.isFinite(state.vmin) ? Number(state.vmin.toPrecision(5)) : "";
  if (document.activeElement !== el.vmaxInput) el.vmaxInput.value = Number.isFinite(state.vmax) ? Number(state.vmax.toPrecision(5)) : "";
  el.limitButtons.querySelectorAll("button").forEach((b) => b.classList.toggle("active", b.dataset.mode === state.limitsMode));
}

function renderInfo() {
  const r = activeResult();
  el.resultInfoBody.innerHTML = "";
  if (!r) return;
  const rows = [
    ["Source", r.source],
    ["Fit type", r.fitType],
    ["Fit mode", r.fitMode],
    ["Regularization", r.regularization || (r.regularized ? "yes" : "none")],
    ["Bins", `${r.nbins} (${r.hasOutput.reduce((a, b) => a + b, 0)} with output)`],
    ["Velocity axis", r.nvel ? `${r.nvel} channels, ${formatValue(r.xvel[0])} to ${formatValue(r.xvel[r.nvel - 1])} km/s` : "–"],
    ["LOSVD rows", r.losvd ? (r.bands ? `${r.rows} percentiles (${PERCENTILES.join(", ")})` : "median only (MAP)") : "–"],
    ["Spaxels", r.grid ? `${r.grid.nx} × ${r.grid.ny} lattice, ${formatValue(r.grid.dx)}″ pixels` : "–"]
  ];
  Object.entries(r.scalars || {}).forEach(([k, v]) => rows.push([`in/${k}`, String(v)]));
  Object.entries(r.attrs || {}).forEach(([k, v]) => rows.push([`attr ${k}`, typeof v === "object" ? JSON.stringify(v) : String(v)]));
  r.notes.forEach((note) => rows.push(["Note", note]));
  const dl = document.createElement("dl");
  rows.forEach(([k, v]) => {
    if (v === undefined || v === null || v === "") return;
    const dt = document.createElement("dt");
    dt.textContent = k;
    const dd = document.createElement("dd");
    dd.textContent = v;
    dl.append(dt, dd);
  });
  el.resultInfoBody.appendChild(dl);
}

function renderAll() {
  if (!activeResult()) return;
  renderResultList();
  renderQuantityButtons();
  syncControls();
  renderMap();
  renderBinPanel();
  renderLosvd();
  renderSpectrum();
  renderInfo();
}

// ------------------------------------------------------------------ export

function download(blob, name) {
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = name;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 2000);
}

function safeName(text) {
  return String(text).replace(/[^a-z0-9._-]+/gi, "_").replace(/^_+|_+$/g, "") || "result";
}

function exportPng() {
  const r = activeResult();
  const map = el.mapCanvas;
  const bar = el.colorbarCanvas;
  const out = document.createElement("canvas");
  out.width = map.width + bar.width;
  out.height = Math.max(map.height, bar.height);
  const ctx = out.getContext("2d");
  ctx.fillStyle = getComputedStyle(document.body).backgroundColor;
  ctx.fillRect(0, 0, out.width, out.height);
  ctx.drawImage(map, 0, 0);
  ctx.drawImage(bar, map.width, 0);
  out.toBlob((blob) => download(blob, `${safeName(r.label)}_${safeName(state.quantity)}.png`), "image/png");
}

function exportCsv() {
  const r = activeResult();
  const keys = availableQuantities(r).filter((k) => r.quantities.has(k) || k.startsWith("gh_"));
  keys.forEach((k) => quantityValues(r, k));
  const text = resultToCsv(r, state.rotation, keys.filter((k) => r.quantities.has(k)));
  download(new Blob([text], {type: "text/csv"}), `${safeName(r.label)}_maps.csv`);
}

function exportBundle() {
  const bundle = makeBundle(app.results, {settings: {rotation: state.rotation, mgePa: state.mgePa}, mge: app.mge});
  const text = JSON.stringify(bundle);
  download(new Blob([text], {type: "application/json"}), "bayeslosvd_viewer_session.json");
  setStatus(`Saved session bundle (${(text.length / 1e6).toFixed(1)} MB; spectra are not included)`);
}

// ------------------------------------------------------------------ theme

const THEMES = ["auto", "light", "dark"];
function applyTheme(mode) {
  if (mode === "auto") document.documentElement.removeAttribute("data-theme");
  else document.documentElement.setAttribute("data-theme", mode);
  el.themeBtn.title = `Theme: ${mode}`;
  el.themeBtn.setAttribute("aria-label", `Theme: ${mode}`);
  el.themeBtn.dataset.mode = mode;
  try {
    localStorage.setItem("bayeslosvd-viewer-theme", mode);
  } catch (error) {
    // storage unavailable
  }
}

// ------------------------------------------------------------------ wiring

function bindControls() {
  el.openBtn.addEventListener("click", () => el.fileInput.click());
  el.welcomeOpenBtn.addEventListener("click", () => el.fileInput.click());
  el.addMoreBtn.addEventListener("click", () => el.fileInput.click());
  el.fileInput.addEventListener("change", () => {
    loadFiles(el.fileInput.files);
    el.fileInput.value = "";
  });
  el.demoBtn.addEventListener("click", loadDemo);

  // drag & drop anywhere
  let dragDepth = 0;
  window.addEventListener("dragenter", (event) => {
    if (!event.dataTransfer?.types?.includes("Files")) return;
    dragDepth += 1;
    el.dropOverlay.hidden = false;
  });
  window.addEventListener("dragleave", () => {
    dragDepth = Math.max(0, dragDepth - 1);
    if (!dragDepth) el.dropOverlay.hidden = true;
  });
  window.addEventListener("dragover", (event) => event.preventDefault());
  window.addEventListener("drop", (event) => {
    event.preventDefault();
    dragDepth = 0;
    el.dropOverlay.hidden = true;
    loadFiles(event.dataTransfer?.files);
  });

  // export menu
  el.exportBtn.addEventListener("click", (event) => {
    event.stopPropagation();
    const list = el.exportMenu.querySelector(".menu-list");
    list.hidden = !list.hidden;
    el.exportBtn.setAttribute("aria-expanded", String(!list.hidden));
  });
  document.addEventListener("click", () => {
    el.exportMenu.querySelector(".menu-list").hidden = true;
    el.exportBtn.setAttribute("aria-expanded", "false");
  });
  el.exportMenu.querySelectorAll("[data-export]").forEach((b) => b.addEventListener("click", () => {
    if (!activeResult()) {
      setStatus("Load a result first", {error: true});
      return;
    }
    ({png: exportPng, csv: exportCsv, bundle: exportBundle})[b.dataset.export]();
  }));

  el.themeBtn.addEventListener("click", () => {
    const next = THEMES[(THEMES.indexOf(el.themeBtn.dataset.mode || "auto") + 1) % THEMES.length];
    applyTheme(next);
    renderAll();
  });
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", renderAll);

  el.compareSelect.addEventListener("change", () => {
    state.compareUid = el.compareSelect.value;
    resetLimitsMode();
    renderAll();
  });
  el.viewButtons.querySelectorAll("button").forEach((b) => b.addEventListener("click", () => {
    state.view = b.dataset.view;
    renderAll();
  }));
  COLORMAP_NAMES.forEach((name) => {
    const o = document.createElement("option");
    o.value = name;
    o.textContent = name;
    el.cmapSelect.appendChild(o);
  });
  el.cmapSelect.addEventListener("change", () => {
    state.cmap[isDiffMode() ? `diff:${state.quantity}` : state.quantity] = el.cmapSelect.value;
    renderMap();
  });
  el.cmapReverse.addEventListener("change", () => {
    state.reverse = el.cmapReverse.checked;
    renderMap();
  });
  el.limitButtons.querySelectorAll("button").forEach((b) => b.addEventListener("click", () => {
    state.limitsMode = b.dataset.mode;
    renderMap();
  }));
  const manual = () => {
    const lo = Number(el.vminInput.value);
    const hi = Number(el.vmaxInput.value);
    if (el.vminInput.value !== "" && el.vmaxInput.value !== "" && Number.isFinite(lo) && Number.isFinite(hi) && hi > lo) {
      state.limitsMode = "manual";
      state.vmin = lo;
      state.vmax = hi;
      renderMap();
    }
  };
  el.vminInput.addEventListener("change", manual);
  el.vmaxInput.addEventListener("change", manual);
  el.rotationInput.addEventListener("change", () => {
    const val = Number(el.rotationInput.value);
    state.rotation = Number.isFinite(val) ? val : 0;
    state.zoom = null;
    renderAll();
  });
  el.maskSigma.addEventListener("change", () => {
    state.maskSigma = el.maskSigma.checked;
    resetLimitsMode();
    renderAll();
  });
  el.maskRadius.addEventListener("change", () => {
    const val = Number(el.maskRadius.value);
    state.maskRadius = val > 0 ? val : 1;
    renderAll();
  });
  el.minSnr.addEventListener("change", () => {
    state.minSnr = el.minSnr.value === "" ? NaN : Number(el.minSnr.value);
    resetLimitsMode();
    renderAll();
  });
  el.showIsophotes.addEventListener("change", () => {
    state.showIsophotes = el.showIsophotes.checked;
    renderMap();
  });
  el.showCenters.addEventListener("change", () => {
    state.showCenters = el.showCenters.checked;
    renderMap();
  });
  el.showMge.addEventListener("change", () => {
    state.showMge = el.showMge.checked;
    renderMap();
  });
  el.mgeLoadBtn.addEventListener("click", () => el.mgeInput.click());
  el.mgeInput.addEventListener("change", async () => {
    const file = el.mgeInput.files?.[0];
    el.mgeInput.value = "";
    if (!file) return;
    try {
      setMge(mgeFromTable(await file.text(), file.name));
      setStatus(`Loaded MGE table ${file.name}`);
    } catch (error) {
      setStatus(`Could not read MGE table: ${error.message}`, {error: true});
    }
  });
  el.mgePa.addEventListener("change", () => {
    const val = Number(el.mgePa.value);
    state.mgePa = Number.isFinite(val) ? val : 0;
    renderMap();
  });
  el.zoomInBtn.addEventListener("click", () => zoomMap(0.6));
  el.zoomOutBtn.addEventListener("click", () => zoomMap(1.6));
  el.zoomResetBtn.addEventListener("click", () => {
    state.zoom = null;
    renderMap();
  });
  const addFromInput = () => {
    const b = Number(el.binInput.value);
    if (Number.isInteger(b)) selectBin(b);
  };
  el.addBinBtn.addEventListener("click", addFromInput);
  el.binInput.addEventListener("keydown", (event) => {
    if (event.key === "Enter") addFromInput();
  });
  el.clearBinsBtn.addEventListener("click", () => {
    state.selected = [];
    state.active = null;
    renderAll();
  });
  el.losvdAll.addEventListener("change", () => {
    state.losvdAll = el.losvdAll.checked;
    renderLosvd();
  });
  el.specResetBtn.addEventListener("click", () => {
    state.specRange = null;
    renderSpectrum();
  });
  document.addEventListener("keydown", (event) => {
    if (event.key === "Escape" && !el.workspace.hidden && !/INPUT|SELECT|TEXTAREA/.test(document.activeElement?.tagName || "")) {
      state.selected = [];
      state.active = null;
      renderAll();
    }
  });

  let resizeFrame = null;
  new ResizeObserver(() => {
    cancelAnimationFrame(resizeFrame);
    resizeFrame = requestAnimationFrame(() => {
      if (!el.workspace.hidden) {
        renderMap();
        renderLosvd();
        renderSpectrum();
      }
    });
  }).observe(el.workspace);
}

async function loadDemo() {
  const base = new URL(".", window.location.href);
  await loadSources(DEMO_FILES.map((path) => ({url: new URL(path, base).href, name: path.split("/").pop()})));
  try {
    const text = await fetch(new URL(DEMO_MGE, base)).then((r) => (r.ok ? r.text() : Promise.reject(new Error(r.status))));
    setMge(mgeFromTable(text, DEMO_MGE.split("/").pop()));
  } catch (error) {
    // optional
  }
  state.rotation = DEMO_ROTATION;
  renderAll();
}

// URL parameters: ?src=<url>&src=<url>&rotation=-18&mge=<url>&demo=1
async function loadFromQuery() {
  const params = new URLSearchParams(window.location.search);
  if (params.get("demo")) {
    await loadDemo();
  }
  const sources = params.getAll("src").map((url) => ({url, name: decodeURIComponent(url.split("?")[0].split("/").pop() || "remote")}));
  if (sources.length) await loadSources(sources);
  const mge = params.get("mge");
  if (mge) {
    try {
      const text = await fetch(mge).then((r) => r.text());
      setMge(mgeFromTable(text, mge.split("/").pop()));
    } catch (error) {
      setStatus(`Could not load MGE table: ${error.message}`, {error: true});
    }
  }
  if (params.has("rotation")) {
    const rot = Number(params.get("rotation"));
    if (Number.isFinite(rot)) {
      state.rotation = rot;
      renderAll();
    }
  }
}

function init() {
  let stored = "auto";
  try {
    stored = localStorage.getItem("bayeslosvd-viewer-theme") || "auto";
  } catch (error) {
    // storage unavailable
  }
  applyTheme(THEMES.includes(stored) ? stored : "auto");
  if (REPO_URL) el.repoLink.href = REPO_URL;
  else el.repoLink.hidden = true;
  bindControls();
  bindMapInteraction();
  bindSpectrumInteraction();
  loadFromQuery();
}

init();

// Exposed for debugging and automated checks.
window.bayeslosvdViewer = {app, state, loadSources, renderAll};
