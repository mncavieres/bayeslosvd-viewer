#!/usr/bin/env python3
"""Write a compact copy of a BAYES-LOSVD results file for sharing or hosting.

Keeps only what the viewer reads and shrinks it:

* in/: bin and spaxel geometry, flux, S/N, velocity axis, observed spectra,
  errors, fit mask and scalar settings (templates and preprocessing arrays
  are dropped);
* out/<bin>/: losvd, bestfit and model parameters (vel, sigma, h1..h6);
  conv_spec, poly, spec, weights and coefs are dropped;
* percentile rows that are all identical (MAP runs) are stored once;
* spectra and best fits are stored as float32, everything gzip-compressed;
* attributes holding local file paths are removed;
* optionally, the sky orientation of the cube is stored (--wcs-fits).

Usage:
    python tools/slim_results.py INPUT_results.hdf5 OUTPUT_results.hdf5 [--wcs-fits CUBE.fits] [--keep-attr-paths]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np

IN_ARRAYS = ["xbin", "ybin", "bin_snr", "bin_flux", "xvel", "x", "y", "binID", "flux", "wave_obs", "mask"]
IN_SPECTRA = ["spec_obs", "sigma_obs"]
IN_SCALARS = ["lmin", "lmax", "velscale", "snr", "porder", "ntemp", "npix", "npix_obs", "nmask", "nvel", "nbins",
              "nspec", "pca_target", "pca_variance"]
OUT_KEEP = ["losvd", "bestfit", "vel", "sigma", "h1", "h2", "h3", "h4", "h5", "h6"]
FLOAT32 = {"spec_obs", "sigma_obs", "bestfit"}
GZIP = {"compression": "gzip", "compression_opts": 6, "shuffle": True}


def dedupe_rows(arr):
    """Collapse identical percentile rows (axis 0) to a single row."""
    if arr.ndim >= 1 and arr.shape[0] > 1 and np.all(arr == arr[:1]):
        return arr[:1]
    return arr


def looks_like_path(value):
    text = value.decode() if isinstance(value, bytes) else str(value)
    return "/" in text and (text.startswith("/") or ":\\" in text)


def slim(src: Path, dst: Path, keep_attr_paths=False, wcs_fits: Path | None = None):
    with h5py.File(src, "r") as fin, h5py.File(dst, "w") as fout:
        for key, value in fin.attrs.items():
            if keep_attr_paths or not looks_like_path(value):
                fout.attrs[key] = value
        fout.attrs["slimmed_from"] = src.name
        if wcs_fits is not None:
            from add_wcs import celestial_cd, write_wcs_attrs  # needs astropy

            write_wcs_attrs(fout, celestial_cd(wcs_fits), wcs_fits.name)
        gin = fout.create_group("in")
        for key in IN_ARRAYS + IN_SPECTRA + IN_SCALARS:
            if f"in/{key}" not in fin:
                continue
            data = fin[f"in/{key}"][()]
            if key in FLOAT32:
                data = np.asarray(data, dtype=np.float32)
            if np.ndim(data) == 0:
                gin.create_dataset(key, data=data)
            else:
                gin.create_dataset(key, data=data, **GZIP)
        gout = fout.create_group("out")
        for bin_key in fin["out"].keys():
            src_bin = fin[f"out/{bin_key}"]
            dst_bin = gout.create_group(bin_key)
            for key in OUT_KEEP:
                if key not in src_bin:
                    continue
                data = dedupe_rows(np.asarray(src_bin[key][()]))
                if key in FLOAT32:
                    data = data.astype(np.float32)
                if data.size > 16:
                    dst_bin.create_dataset(key, data=data, **GZIP)
                else:
                    dst_bin.create_dataset(key, data=data)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--keep-attr-paths", action="store_true", help="keep attributes that contain file paths")
    parser.add_argument("--wcs-fits", type=Path, default=None, help="cube whose WCS orients the maps on the sky")
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    slim(args.input, args.output, keep_attr_paths=args.keep_attr_paths, wcs_fits=args.wcs_fits)
    before = args.input.stat().st_size / 1e6
    after = args.output.stat().st_size / 1e6
    print(f"{args.output}: {before:.1f} MB -> {after:.1f} MB")


if __name__ == "__main__":
    main()
