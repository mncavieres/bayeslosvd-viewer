# BAYES-LOSVD Viewer

An interactive, browser-based viewer for [BAYES-LOSVD](https://github.com/jfalconbarroso/BAYES-LOSVD) results.
Open one or more `*_results.hdf5` files and explore kinematic maps, line-of-sight velocity distributions
(LOSVDs), and spectral fits, with no installation, no server and no upload. Files are read locally in
your browser with [h5wasm](https://github.com/usnistgov/h5wasm).

**Live app:** https://mncavieres.github.io/bayeslosvd-viewer/ ·
**[Open the Sombrero demo](https://mncavieres.github.io/bayeslosvd-viewer/?demo=1)**

The demo opens in the major-axis orientation (PA 90°); switch to *Sky* or *Data* in the sidebar. It loads BAYES-LOSVD fits to the JWST/NIRSpec G235H IFU cube of the Sombrero galaxy (M104), with
the AGN contribution subtracted: a MAP fit without regularization (`GaussianMAP`) and one with
random-walk regularization (`RWregularizedMAP`), plus the galaxy's MGE model. Use *Compare with* to see
how regularization changes the maps.

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
- **Sky orientation.** Switch between the data frame (pixel axes, optionally rotated), the sky
  (North up, East left, ΔRA/ΔDec axes, compass) and a major-axis frame (galaxy PA horizontal),
  using the WCS of the cube. Results from different instruments share one sky frame.
- **Display controls.** Colormaps (incl. *sauron*), robust / min–max /
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
| `orientation` | `?orientation=sky` | `data`, `sky` (North up, East left) or `major` (major axis horizontal). |
| `pa` | `?pa=90` | Galaxy major-axis PA (degrees East of North), used by `orientation=major`. |
| `rotation` | `?rotation=-18` | Rotation of the data orientation (degrees, counter-clockwise). |
| `demo` | `?demo=1` | Load the Sombrero NIRSpec demo. |

`src` is handy for sharing: put a session bundle (`.json`) in a public place (a GitHub repo, Zenodo)
and send a link such as `https://mncavieres.github.io/bayeslosvd-viewer/?src=<raw-url-of-bundle>`.

## Supported inputs

| File | What is read |
|---|---|
| `.hdf5`, `.h5` | BAYES-LOSVD results (see below). |
| `.json` | Session bundles exported by this viewer, or the `DATA` payload of a legacy dashboard. |
| `.html` | Legacy dashboards from `build_kinematics_dashboard.py` (the embedded `const DATA = …`). |
| `.fits` | The cube used for the fit (or any image on the same pixel grid). Only the header is read, for the sky orientation (see below). |
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
attributes WCS_CD1_1 … WCS_CD2_2, WCS_SOURCE      optional sky orientation (see below)
```

Partial runs (only some bins in `out/`) are fine; missing bins are shown as empty.

## Orientation

BAYES-LOSVD's instrument readers define `x = (column − c0)·pixel` and `y = (row − r0)·pixel`, offsets
along the cube's pixel axes from the brightest spaxel. The cube's spatial WCS matrix therefore maps
(x, y) onto the sky, and the viewer offers three orientations:

| Orientation | Axes | Needs |
|---|---|---|
| **Data** | x, y as stored, optionally rotated by *Data rotation* | — |
| **Sky** | ΔRA (increasing to the left, East) and ΔDec (North up) | WCS |
| **Major axis** | the galaxy major axis along +x′, keeping North in the upper half (tick *Flip* for the other end) | WCS and galaxy PA |

The WCS can come from

- **the FITS file of the cube**, dropped on the page or loaded with *Load WCS from FITS* (only the
  header is read). It is applied to every loaded result whose spaxel size matches the FITS pixel scale;
- **attributes in the results file**, written once with

  ```bash
  python tools/add_wcs.py run_results.hdf5 cube.fits    # needs h5py, astropy
  ```

  or while slimming with `tools/slim_results.py … --wcs-fits cube.fits`.

In the Sky and Major-axis orientations, results from different instruments are drawn in the same frame
and matched by sky position (difference maps, *Same bin in all results*). Each result is still centred on
its own brightest spaxel, so this assumes those coincide on the sky. The MGE overlay follows the galaxy
major axis; its *Angle* is measured from that axis.

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
- **Data rotation.** In the Data orientation, coordinates are rotated counter-clockwise by the given
  angle (x′ = x cos θ − y sin θ, y′ = x sin θ + y cos θ).
- **MGE overlay.** Contours of the MGE surface brightness centred on (0, 0). Their major axis lies along
  the galaxy major axis plus *Angle* when the galaxy PA and WCS are known, otherwise at *Angle* from +x.
- **Isophotes and MGE contours** are drawn every 0.5 mag below the peak.

## Sharing results files

BAYES-LOSVD results contain arrays the viewer never reads (templates, convolved spectra) and, for MAP
runs, five identical percentile rows per bin. `tools/slim_results.py` writes a compact copy (typically
6–7× smaller) that is easy to host or send, and drops attributes containing local file paths:

```bash
python tools/slim_results.py run_results.hdf5 run_slim_results.hdf5   # needs numpy, h5py
```

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
assets/js/frames.js     FITS header WCS reader and orientation matrices
assets/js/analysis.js   moments, Gauss–Hermite fits, statistics, grids, contours, MGE
assets/js/hdf5-worker.js  Web Worker that reads HDF5 with h5wasm (WORKERFS for local files)
assets/js/colormaps.js  colormap lookup tables
assets/js/plot.js       canvas axes, ticks and colorbars
demo/                   Sombrero NIRSpec demo (slimmed results + MGE table)
tests/                  browser test page, synthetic data (tests/data) and reference values
```

External code is loaded from jsDelivr at pinned versions: `h5wasm@0.10.3` (HDF5 reader, loaded by the
worker) and `delaunator@5.0.1` (triangulation for interpolated maps).

**Tests.** Serve the repository and open `/tests/`. The page loads the synthetic files in `tests/data/` through the real
worker and compares moments and Gauss–Hermite fits against NumPy/SciPy values in
`tests/reference.json`. To regenerate the synthetic data and the reference values:

```bash
python tools/make_demo_data.py      # needs numpy, h5py
python tools/make_reference.py      # needs numpy, h5py, scipy
```

**Deployment.** `.github/workflows/pages.yml` publishes the repository to GitHub Pages on every push
to `main`. It replaces `__VERSION__` in `index.html` with the commit hash, so browsers fetch a
consistent set of modules after each update instead of mixing cached ones. In the repository settings, set *Pages → Source* to *GitHub Actions*.

## Citation

If you use BAYES-LOSVD, please cite Falcón-Barroso & Martig 2021, A&A, 646, A31.
