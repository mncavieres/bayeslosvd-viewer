#!/usr/bin/env python3
"""Compute reference values from the demo files for tests/index.html.

Moments use the same definitions as the viewer (and the original
build_kinematics_dashboard.py); Gauss-Hermite parameters come from
scipy.optimize.curve_fit.

Usage:  python tools/make_reference.py
"""

import json
from pathlib import Path

import h5py
import numpy as np
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parent.parent
BINS = [0, 5, 40, 87, 120, 173]


def gauss_hermite(v, amp, vel, sigma, h3, h4):
    w = (v - vel) / sigma
    H3 = (2 * np.sqrt(2) * w**3 - 3 * np.sqrt(2) * w) / np.sqrt(6)
    H4 = (4 * w**4 - 12 * w**2 + 3) / np.sqrt(24)
    return amp / (np.sqrt(2 * np.pi) * sigma) * np.exp(-0.5 * w**2) * (1 + h3 * H3 + h4 * H4)


def moments(x, y):
    area = y.sum()
    mean = (x * y).sum() / area
    sigma = np.sqrt(((x - mean) ** 2 * y).sum() / area)
    w = (x - mean) / sigma
    return mean, sigma, (w**3 * y).sum() / area, (w**4 * y).sum() / area - 3


def main():
    out = {}
    for name in ["demo_disk_sampling_results.hdf5", "demo_disk_regularized_results.hdf5"]:
        with h5py.File(ROOT / "demo" / name, "r") as f:
            xvel = f["in/xvel"][:]
            order = np.argsort(xvel)
            x = xvel[order]
            rows = []
            for b in BINS:
                losvd = f[f"out/{b}/losvd"][2][order]
                losvd = losvd / np.trapz(losvd, x)
                vel, sigma, skew, kurt = moments(x, losvd)
                popt, _ = curve_fit(gauss_hermite, x, losvd, p0=[np.trapz(losvd, x), vel, sigma, 0, 0], maxfev=20000)
                rows.append({
                    "bin": b, "vel": vel, "sigma": sigma, "skew": skew, "kurt": kurt,
                    "gh": {"vel": popt[1], "sigma": abs(popt[2]), "h3": popt[3], "h4": popt[4]},
                })
            out[name] = {"nbins": int(f["in/nbins"][()]), "bins": rows}
    path = ROOT / "tests" / "reference.json"
    path.write_text(json.dumps(out, indent=1) + "\n")
    print(path)


if __name__ == "__main__":
    main()
