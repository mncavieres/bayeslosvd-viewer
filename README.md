# BAYES-LOSVD Viewer

An interactive, browser-based viewer for [BAYES-LOSVD](https://github.com/jfalconbarroso/BAYES-LOSVD) results.
Open one or more `*_results.hdf5` files and explore kinematic maps, line-of-sight velocity distributions
(LOSVDs), and spectral fits, with no installation, no server and no upload. Files are read locally in
your browser with [h5wasm](https://github.com/usnistgov/h5wasm).

**Live app:** https://mncavieres.github.io/bayeslosvd-viewer/ (append `?demo=1` to load synthetic demo data)

## Features

- **Reads BAYES-LOSVD output directly.** Drop `*_results.hdf5` files on the page. Local files are mounted
  lazily inside a Web Worker, so multi-hundred-MB MUSE results open in a few seconds.
- **Maps.** Mean velocity, dispersion, skewness and excess kurtosis of each bin's median LOSVD;
  Gauss–Hermite V, σ, h3, h4 fitted to the LOSVDs; model parameters stored in the file (e.g. `vel`,
  `sigma`, `h3`, `h4` for GH fits); bin S/N; and log flux.
- **Three display modes.** Spaxels (each spaxel coloured by its bin, using `in/x`, `in/y`, `in/binID`),
  bin centres, or Delaunay-interpolated maps.
- **LOSVD panel.** Up to eight bins at once, with 68% and 99.8% credible bands when the run stored
  percentiles (sampling runs). Optionally overlays the same sky position from every loaded result.
- **Spectral fit panel.** Observed spectrum, median best fit (with band), residuals, pixels excluded
  from the fit shaded, and χ²/N. Scroll to zoom, drag to pan.
- **Comparison.** Difference maps between any two results (bins matched by position).
- **Display controls.** Rotation of the map frame, colormaps (incl. *sauron*), robust / min–max /
  symmetric / manual colour limits, σ-outlier and S/N masks, flux isophotes, MGE contour overlays.
- **Export.** PNG of the map, CSV of all maps for a result, and a compact JSON *session bundle*
  (≈1–2 MB even for large MUSE runs) that anyone can reopen in the viewer.
- **Backwards compatible.** Opens the self-contained HTML dashboards made by
  `build_kinematics_dashboard.py` and their JSON payloads.
- Light and dark themes; works on tablets and phones.

## Using it

### Online

Open the live app and drag files onto the page, or use **Open files**. Nothing leaves your machine.

### Locally

The app is plain HTML/CSS/JavaScript with no build step, but it must be served over HTTP (browsers do
not run module Web Workers from `file://`):

```bash
git clone https://github.com/mncavieres/bayeslosvd-viewer.git
cd bayeslosvd-viewer
python3 -m http.server 8000
# open http://localhost:8000
```

### URL parameters

| Parameter | Example | Effect |
|---|---|---|
| `src` (repeatable) | `?src=https://example.org/run_results.hdf5` | Load files from URLs (HDF5, JSON bundle, HTML dashboard, CSV). The host must allow CORS. |
| `mge` | `?mge=https://example.org/mge.csv` | Load an MGE table. |
| `rotation` | `?rotation=-18` | Rotate the map frame (degrees, counter-clockwise). |
| `demo` | `?demo=1` | Load the synthetic demo data. |

`src` is handy for sharing: put a session bundle (`.json`) in a public place (a GitHub repo, Zenodo)
and send a link such as `https://mncavieres.github.io/bayeslosvd-viewer/?src=<raw-url-of-bundle>`.

## Supported inputs

| File | What is read |
|---|---|
| `.hdf5`, `.h5` | BAYES-LOSVD results (see below). |
| `.json` | Session bundles exported by this viewer, or the `DATA` payload of a legacy dashboard. |
| `.html` | Legacy dashboards from `build_kinematics_dashboard.py` (the embedded `const DATA = …`). |
| `.csv`, `.txt` | Bin tables: `x`, `y` (arcsec) plus any numeric columns, each becoming a map (`vel`, `sigma`, `h3`, `h4`, `snr` are recognised). Or MGE tables: `sigma_arcsec`, `q_obs` and `luminosity_Lsun` / `total_flux` (or `surf`). |

### Expected HDF5 layout

This is the layout written by `bayes_losvd_run.py`. Only `in/xbin`, `in/ybin`, `in/xvel` and
`out/<bin>/losvd` are required; everything else is used when present.

```
in/xbin, in/ybin          bin centres [arcsec]
in/xvel                   LOSVD velocity axis [km/s] (any order)
in/bin_snr, in/bin_flux   per-bin S/N and flux
in/x, in/y, in/binID      spaxel centres and bin membership  -> spaxel maps
in/flux                   spaxel flux                        -> flux map, isophotes
in/wave_obs               ln(wavelength) of the fitted spectra
in/spec_obs, in/sigma_obs observed spectra and errors, shape (npix_obs, nbins)
in/mask                   indices of fitted pixels
out/<bin>/losvd           (5, nvel) percentiles 0.1, 15.9, 50, 84.1, 99.9 (MAP runs: 5 identical rows)
out/<bin>/bestfit         (5, npix_obs) best-fit model percentiles
out/<bin>/vel, sigma, h1…h6   model parameters (parametric fits), shown as "Model parameters"
attributes fit_type, fit_mode, regularization   shown in the result list and details
```

Partial runs (only some bins in `out/`) are fine; missing bins are shown as empty.

## Definitions

- **LOSVD moments.** Each bin's median LOSVD is normalized to unit area. V is its mean, σ its standard
  deviation, *skew* the standardized third moment and *kurt* the standardized fourth moment minus 3.
  These are not Gauss–Hermite coefficients (for small values skew ≈ 4√3·h3 and kurt ≈ 8√6·h4). The
  definitions match `build_kinematics_dashboard.py`.
- **Gauss–Hermite fit.** A least-squares fit of a Gauss–Hermite series truncated at h4
  (van der Marel & Franx 1993) to the median LOSVD, computed in the browser. It is a quick
  characterisation, not BAYES-LOSVD's own Bayesian GH fit (`bayes_losvd_ghfit.py`).
- **σ outliers.** Bins whose σ exceeds the largest σ found within a chosen radius of the centre
  (default 1″), as in the original dashboards.
- **Rotation.** Bin and spaxel coordinates are rotated counter-clockwise by the given angle
  (x′ = x cos θ − y sin θ, y′ = x sin θ + y cos θ). Use it to align the major axis with x′.
- **MGE overlay.** Contours of the MGE surface brightness, centred on (0, 0), with the major axis at the
  given PA measured counter-clockwise from +x′ in the displayed (rotated) frame.
- **Isophotes and MGE contours** are drawn every 0.5 mag below the peak.

## Session bundles

**Export → Session bundle** writes every loaded result (bins, LOSVDs, model parameters, spaxel
lattice, flux, metadata), the MGE table and the display rotation to one JSON file. Duplicated MAP
percentile rows are stored once. Spectra are not included; share the HDF5 files if the spectral fits
are needed.

## Development

```
index.html              page shell
assets/css/app.css      styles (light/dark themes)
assets/js/main.js       state, UI, map/LOSVD/spectrum rendering, interaction
assets/js/data.js       input formats -> common result model; bundles; CSV export
assets/js/analysis.js   moments, Gauss–Hermite fits, statistics, grids, contours, MGE
assets/js/hdf5-worker.js  Web Worker that reads HDF5 with h5wasm (WORKERFS for local files)
assets/js/colormaps.js  colormap lookup tables
assets/js/plot.js       canvas axes, ticks and colorbars
demo/                   synthetic BAYES-LOSVD files (tools/make_demo_data.py)
tests/                  browser test page + reference values (tools/make_reference.py)
```

External code is loaded from jsDelivr at pinned versions: `h5wasm@0.10.3` (HDF5 reader, loaded by the
worker) and `delaunator@5.0.1` (triangulation for interpolated maps).

**Tests.** Serve the repository and open `/tests/`. The page loads the demo files through the real
worker and compares moments and Gauss–Hermite fits against NumPy/SciPy values in
`tests/reference.json`. To regenerate the demo data and the reference values:

```bash
python tools/make_demo_data.py      # needs numpy, h5py
python tools/make_reference.py      # needs numpy, h5py, scipy
```

**Deployment.** `.github/workflows/pages.yml` publishes the repository to GitHub Pages on every push
to `main`. In the repository settings, set *Pages → Source* to *GitHub Actions*.

## Citation

If you use BAYES-LOSVD, please cite Falcón-Barroso & Martig 2021, A&A, 646, A31.
