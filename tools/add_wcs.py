#!/usr/bin/env python3
"""Store the sky orientation of a cube in a BAYES-LOSVD results file.

BAYES-LOSVD's (x, y) are offsets along the cube's pixel axes, so the spatial
CD matrix of the cube maps them to sky offsets. This writes it as root
attributes, which the viewer uses for its Sky and Major-axis orientations:

    WCS_CD1_1, WCS_CD1_2, WCS_CD2_1, WCS_CD2_2   deg/pixel; rows RA (East), Dec (North)
    WCS_SOURCE                                   name of the FITS file

Usage:
    python tools/add_wcs.py RESULTS.hdf5 CUBE.fits [--hdu N]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


def celestial_cd(fits_path: Path, hdu=None):
    """Spatial CD matrix (deg/pixel, rows East/North) of the first HDU with a celestial WCS."""
    with fits.open(fits_path) as hdul:
        candidates = [hdul[hdu]] if hdu is not None else list(hdul)
        for h in candidates:
            if "CTYPE1" not in h.header and "CTYPE2" not in h.header:
                continue
            wcs = WCS(h.header).celestial
            if wcs.naxis == 2 and wcs.wcs.lng >= 0:
                cd = np.asarray(wcs.pixel_scale_matrix, dtype=float)
                if wcs.wcs.lng == 1:  # Dec listed before RA: put East first
                    cd = cd[::-1]
                return cd
    raise ValueError(f"no celestial WCS found in {fits_path}")


def write_wcs_attrs(h5file: h5py.File, cd, source_name: str):
    for i in range(2):
        for j in range(2):
            h5file.attrs[f"WCS_CD{i + 1}_{j + 1}"] = float(cd[i, j])
    h5file.attrs["WCS_SOURCE"] = source_name


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("results", type=Path)
    parser.add_argument("fits", type=Path)
    parser.add_argument("--hdu", type=int, default=None, help="HDU index holding the WCS (default: first with RA/Dec)")
    args = parser.parse_args()
    cd = celestial_cd(args.fits, args.hdu)
    with h5py.File(args.results, "a") as f:
        write_wcs_attrs(f, cd, args.fits.name)
    scale = np.sqrt(abs(np.linalg.det(cd))) * 3600
    pa_y = np.degrees(np.arctan2(cd[0, 1], cd[1, 1])) % 360
    print(f"{args.results}: pixel {scale:.4f} arcsec, data +y axis at PA {pa_y:.2f} deg")


if __name__ == "__main__":
    main()
