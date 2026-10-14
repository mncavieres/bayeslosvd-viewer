/* Classic Web Worker that reads BAYES-LOSVD *_results.hdf5 files with h5wasm.
 *
 * Local files are mounted with Emscripten's WORKERFS, so only the parts of the
 * file that are actually read get pulled into memory (results files can be
 * hundreds of MB). Remote files (URLs) are downloaded into memory.
 */
/* global h5wasm */

const H5WASM_URL = "https://cdn.jsdelivr.net/npm/h5wasm@0.10.3/dist/iife/h5wasm.js";

let loadError = null;
try {
  self.importScripts(H5WASM_URL);
} catch (error) {
  loadError = `Could not load the HDF5 reader (h5wasm) from ${H5WASM_URL}: ${error.message}`;
}

const files = new Map();
let counter = 0;

function post(type, requestId, payload = {}, transfer = []) {
  self.postMessage({type, requestId, ...payload}, transfer);
}

function toNumber(v) {
  return typeof v === "bigint" ? Number(v) : v;
}

function toFloat64(value) {
  if (value === null || value === undefined) return null;
  if (value instanceof Float64Array) return value;
  if (typeof value === "number" || typeof value === "bigint") return Float64Array.of(toNumber(value));
  return Float64Array.from(value, toNumber);
}

function jsonSafe(value) {
  if (value === null || value === undefined) return null;
  if (typeof value === "bigint") return Number(value);
  if (typeof value === "number" || typeof value === "string" || typeof value === "boolean") return value;
  if (ArrayBuffer.isView(value) || Array.isArray(value)) {
    const arr = Array.from(value, jsonSafe);
    return arr.length === 1 ? arr[0] : arr;
  }
  return String(value);
}

function dataset(h5, path) {
  try {
    const node = h5.get(path);
    return node && node.type === "Dataset" ? node : null;
  } catch (error) {
    return null;
  }
}

function read(h5, path) {
  const ds = dataset(h5, path);
  return ds ? ds.value : null;
}

function readScalar(h5, path) {
  const v = read(h5, path);
  if (v === null) return null;
  if (typeof v === "number" || typeof v === "bigint") return toNumber(v);
  return v.length ? toNumber(v[0]) : null;
}

function readAttrs(node) {
  const out = {};
  try {
    const attrs = node.attrs || {};
    for (const [key, attr] of Object.entries(attrs)) {
      try {
        out[key] = jsonSafe(attr.value);
      } catch (error) {
        out[key] = "(unreadable)";
      }
    }
  } catch (error) {
    // ignore
  }
  return out;
}

async function fetchWithProgress(url, requestId) {
  const response = await fetch(url);
  if (!response.ok) throw new Error(`HTTP ${response.status} while fetching ${url}`);
  const total = Number(response.headers.get("content-length")) || 0;
  if (!response.body || !total) return new Uint8Array(await response.arrayBuffer());
  const reader = response.body.getReader();
  const out = new Uint8Array(total);
  let received = 0;
  let lastPost = 0;
  for (;;) {
    const {done, value} = await reader.read();
    if (done) break;
    if (received + value.length > out.length) {
      // content-length was wrong (e.g. compressed transfer); fall back to concatenation
      const grown = new Uint8Array(Math.max(out.length * 2, received + value.length));
      grown.set(out.subarray(0, received));
      return finishGrowing(reader, grown, value, received, requestId, total);
    }
    out.set(value, received);
    received += value.length;
    if (received - lastPost > 2e6) {
      lastPost = received;
      post("progress", requestId, {stage: "download", done: received, total});
    }
  }
  return out.subarray(0, received);
}

async function finishGrowing(reader, buffer, firstChunk, received, requestId, total) {
  let buf = buffer;
  buf.set(firstChunk, received);
  let size = received + firstChunk.length;
  for (;;) {
    const {done, value} = await reader.read();
    if (done) break;
    if (size + value.length > buf.length) {
      const grown = new Uint8Array(Math.max(buf.length * 2, size + value.length));
      grown.set(buf.subarray(0, size));
      buf = grown;
    }
    buf.set(value, size);
    size += value.length;
    post("progress", requestId, {stage: "download", done: size, total});
  }
  return buf.subarray(0, size);
}

async function openFile(msg) {
  if (loadError) throw new Error(loadError);
  const {FS} = await h5wasm.ready;
  counter += 1;
  const id = `f${counter}`;
  let path;
  let mountDir = null;
  if (msg.file) {
    mountDir = `/mnt_${id}`;
    FS.mkdir(mountDir);
    FS.mount(FS.filesystems.WORKERFS, {files: [msg.file]}, mountDir);
    path = `${mountDir}/${msg.file.name}`;
  } else {
    const bytes = msg.buffer ? new Uint8Array(msg.buffer) : await fetchWithProgress(msg.url, msg.requestId);
    path = `/mem_${id}.h5`;
    FS.writeFile(path, bytes);
  }

  let h5;
  try {
    h5 = new h5wasm.File(path, "r");
  } catch (error) {
    cleanup({path, mountDir});
    throw new Error(`Not a readable HDF5 file (${error.message || error})`);
  }
  const entry = {h5, path, mountDir, cache: {}};
  files.set(id, entry);
  try {
    return {id, summary: summarize(h5, msg.requestId)};
  } catch (error) {
    closeFile(id);
    throw error;
  }
}

function summarize(h5, requestId) {
  const xbin = read(h5, "in/xbin");
  const ybin = read(h5, "in/ybin");
  if (!xbin || !ybin) {
    throw new Error("This HDF5 file has no in/xbin and in/ybin datasets; it does not look like a BAYES-LOSVD results file.");
  }
  const nbins = xbin.length;
  const xvelRaw = read(h5, "in/xvel");
  const nvel = xvelRaw ? xvelRaw.length : 0;

  let outKeys = [];
  try {
    outKeys = h5.get("out").keys();
  } catch (error) {
    outKeys = [];
  }
  const binKeys = outKeys.map(Number).filter((k) => Number.isInteger(k) && k >= 0 && k < nbins).sort((a, b) => a - b);

  // Discover the per-bin datasets and the number of percentile rows.
  let rows = 1;
  const paramKeys = [];
  if (binKeys.length) {
    const first = h5.get(`out/${binKeys[0]}`);
    const keys = first.keys();
    const losvdDs = dataset(h5, `out/${binKeys[0]}/losvd`);
    if (losvdDs && losvdDs.shape.length === 2) rows = losvdDs.shape[0];
    for (const key of keys) {
      const ds = dataset(h5, `out/${binKeys[0]}/${key}`);
      if (ds && (ds.shape.length === 0 || (ds.shape.length === 1 && ds.shape[0] <= 9))) paramKeys.push(key);
    }
  }
  const paramRows = {};
  const params = {};
  paramKeys.forEach((key) => {
    params[key] = new Float64Array(nbins * 5).fill(NaN);
    paramRows[key] = 1;
  });

  const losvd = new Float64Array(nbins * rows * nvel).fill(NaN);
  const hasOutput = new Uint8Array(nbins);
  let lastPost = 0;
  binKeys.forEach((bin, n) => {
    const raw = read(h5, `out/${bin}/losvd`);
    if (raw && nvel) {
      hasOutput[bin] = 1;
      const nr = Math.min(rows, Math.floor(raw.length / nvel));
      for (let r = 0; r < nr; r += 1) {
        for (let v = 0; v < nvel; v += 1) {
          losvd[(bin * rows + r) * nvel + v] = toNumber(raw[r * nvel + v]);
        }
      }
    }
    paramKeys.forEach((key) => {
      const val = read(h5, `out/${bin}/${key}`);
      if (val === null) return;
      const arr = typeof val === "number" || typeof val === "bigint" ? [val] : val;
      const nr = Math.min(5, arr.length);
      paramRows[key] = Math.max(paramRows[key], nr);
      for (let r = 0; r < nr; r += 1) params[key][bin * 5 + r] = toNumber(arr[r]);
    });
    if (n - lastPost >= 100) {
      lastPost = n;
      post("progress", requestId, {stage: "bins", done: n, total: binKeys.length});
    }
  });

  const hasSpectra = Boolean(dataset(h5, "in/spec_obs")) && binKeys.length > 0
    && Boolean(dataset(h5, `out/${binKeys[0]}/bestfit`));

  const scalars = {};
  ["lmin", "lmax", "velscale", "snr", "porder", "ntemp", "npix", "npix_obs", "nmask", "nvel", "nbins", "nspec", "pca_target", "pca_variance"]
    .forEach((key) => {
      const v = readScalar(h5, `in/${key}`);
      if (v !== null && v !== undefined) scalars[key] = v;
    });

  const summary = {
    attrs: readAttrs(h5),
    scalars,
    nbins,
    nvel,
    rows,
    xbin: toFloat64(xbin),
    ybin: toFloat64(ybin),
    binSnr: toFloat64(read(h5, "in/bin_snr")),
    binFlux: toFloat64(read(h5, "in/bin_flux")),
    xvel: toFloat64(xvelRaw),
    spaxelX: toFloat64(read(h5, "in/x")),
    spaxelY: toFloat64(read(h5, "in/y")),
    spaxelBin: toFloat64(read(h5, "in/binID")),
    spaxelFlux: toFloat64(read(h5, "in/flux")),
    losvd,
    hasOutput,
    params,
    paramRows,
    hasSpectra
  };
  return summary;
}

function spectrum(msg) {
  const entry = files.get(msg.id);
  if (!entry) throw new Error("File is no longer open");
  const {h5, cache} = entry;
  if (!cache.waveObs) {
    const w = toFloat64(read(h5, "in/wave_obs"));
    // BAYES-LOSVD stores ln(lambda); convert to linear wavelength.
    const logScale = w && w.length && Math.max(...w) < 25;
    cache.waveObs = w ? (logScale ? w.map(Math.exp) : w) : null;
    cache.mask = toFloat64(read(h5, "in/mask"));
    cache.specDs = dataset(h5, "in/spec_obs");
    cache.errDs = dataset(h5, "in/sigma_obs");
  }
  const bin = msg.bin;
  const column = (ds) => {
    if (!ds) return null;
    if (ds.shape.length !== 2) return null;
    // spec_obs is stored as (npix_obs, nbins)
    const byColumn = ds.shape[1] > bin && ds.shape[0] === (cache.waveObs ? cache.waveObs.length : ds.shape[0]);
    return toFloat64(byColumn ? ds.slice([[], [bin, bin + 1]]) : ds.slice([[bin, bin + 1], []]));
  };
  const obs = column(cache.specDs);
  const err = column(cache.errDs);
  const bestfitDs = dataset(h5, `out/${bin}/bestfit`);
  const bestfit = bestfitDs ? toFloat64(bestfitDs.value) : null;
  const bestfitRows = bestfitDs ? (bestfitDs.shape.length === 2 ? bestfitDs.shape[0] : 1) : 0;
  const wave = cache.waveObs ? Float64Array.from(cache.waveObs) : null;
  const mask = cache.mask ? Float64Array.from(cache.mask) : null;
  const transfer = [obs, err, bestfit, wave, mask].filter(Boolean).map((a) => a.buffer);
  return {payload: {bin, wave, obs, err, bestfit, bestfitRows, mask}, transfer};
}

function cleanup({path, mountDir}) {
  try {
    const {FS} = h5wasm;
    if (mountDir) {
      FS.unmount(mountDir);
      FS.rmdir(mountDir);
    } else if (path) {
      FS.unlink(path);
    }
  } catch (error) {
    // ignore
  }
}

function closeFile(id) {
  const entry = files.get(id);
  if (!entry) return;
  try {
    entry.h5.close();
  } catch (error) {
    // ignore
  }
  cleanup(entry);
  files.delete(id);
}

self.onmessage = async (event) => {
  const msg = event.data;
  try {
    if (msg.type === "open") {
      const {id, summary} = await openFile(msg);
      const transfer = Object.values(summary)
        .filter((v) => ArrayBuffer.isView(v))
        .map((v) => v.buffer);
      Object.values(summary.params).forEach((v) => transfer.push(v.buffer));
      post("result", msg.requestId, {payload: {id, summary}}, transfer);
    } else if (msg.type === "spectrum") {
      const {payload, transfer} = spectrum(msg);
      post("result", msg.requestId, {payload}, transfer);
    } else if (msg.type === "close") {
      closeFile(msg.id);
      post("result", msg.requestId, {payload: true});
    }
  } catch (error) {
    post("error", msg.requestId, {message: error && error.message ? error.message : String(error)});
  }
};
