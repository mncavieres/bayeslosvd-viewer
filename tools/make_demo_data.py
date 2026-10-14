#!/usr/bin/env python3
"""Write small synthetic BAYES-LOSVD results files for the viewer demo.

The files follow the layout written by BAYES-LOSVD's ``bayes_losvd_run.py``:

    in/   binning, spaxel and spectral inputs (xbin, ybin, x, y, binID, flux,
          bin_snr, xvel, wave_obs [ln lambda], spec_obs, sigma_obs, mask, ...)
    out/<bin>/  losvd (5 x nvel percentiles), bestfit (5 x npix_obs), ...

The galaxy is a rotating, inclined disc with a central dispersion peak whose
major axis lies 25 degrees counter-clockwise from +x. Two files are produced:

* ``demo_disk_sampling_results.hdf5``: "sampling" run with credible bands.
* ``demo_disk_regularized_results.hdf5``: MAP run (five identical rows),
  random-walk regularized, slightly smoother LOSVDs.

Usage:  python tools/make_demo_data.py [output_dir]
"""

from __future__ import annotations

import sys
from pathlib import Path

import h5py
import numpy as np

C_KMS = 299792.458
PA_DEG = 25.0
Q_DISK = 0.65
PERCENTILE_Z = np.array([-3.09, -1.0, 0.0, 1.0, 3.09])  # 0.1, 15.9, 50, 84.1, 99.9 percent
MGE = [  # sigma [arcsec], q, total luminosity
    (0.15, 0.90, 4.0),
    (0.45, 0.80, 6.0),
    (1.10, 0.68, 9.0),
    (2.60, 0.62, 12.0),
]


def major_minor(x, y):
    t = np.deg2rad(PA_DEG)
    return x * np.cos(t) + y * np.sin(t), -x * np.sin(t) + y * np.cos(t)


def surface_brightness(x, y):
    xp, yp = major_minor(x, y)
    total = np.zeros_like(x, dtype=float)
    for sigma, q, lum in MGE:
        total += lum / (2 * np.pi * sigma**2 * q) * np.exp(-0.5 * ((xp / sigma) ** 2 + (yp / (sigma * q)) ** 2))
    return total


def kinematics(x, y):
    xp, yp = major_minor(x, y)
    radius = np.hypot(xp, yp / Q_DISK)
    cos_phi = np.divide(xp, radius, out=np.zeros_like(xp), where=radius > 0)
    vel = 230.0 * np.tanh(radius / 0.7) * cos_phi
    sigma = 75.0 + 150.0 * np.exp(-radius / 0.9)
    h3 = -0.09 * vel / 230.0
    h4 = 0.05 * np.exp(-radius / 1.5) - 0.01
    return vel, sigma, h3, h4


def gauss_hermite(v, vel, sigma, h3, h4):
    w = (v - vel) / sigma
    H3 = (2 * np.sqrt(2) * w**3 - 3 * np.sqrt(2) * w) / np.sqrt(6)
    H4 = (4 * w**4 - 12 * w**2 + 3) / np.sqrt(24)
    return np.exp(-0.5 * w**2) / (np.sqrt(2 * np.pi) * sigma) * (1 + h3 * H3 + h4 * H4)


def block_binning(x, y):
    """Bins of growing size with radius (a crude stand-in for Voronoi bins)."""
    r = np.hypot(x, y)
    step = np.min(np.diff(np.unique(np.round(x, 6))))
    ix = np.round((x - x.min()) / step).astype(int)
    iy = np.round((y - y.min()) / step).astype(int)
    block = np.select([r < 0.9, r < 2.0, r < 3.0], [1, 2, 3], default=4)
    keys = [f"{b}:{i // b}:{j // b}" for b, i, j in zip(block, ix, iy)]
    order = {}
    bin_id = np.empty(len(x), dtype=np.int64)
    for k, key in enumerate(keys):
        bin_id[k] = order.setdefault(key, len(order))
    return bin_id


def template(wave):
    """Continuum with the Ca II triplet and a few weaker absorption lines."""
    lines = [(8498.0, 0.38, 1.4), (8542.1, 0.55, 1.8), (8662.1, 0.48, 1.7),
             (8514.1, 0.10, 0.9), (8582.3, 0.08, 0.9), (8688.6, 0.12, 1.0), (8736.0, 0.07, 0.8)]
    spec = np.ones_like(wave)
    for centre, depth, width in lines:
        spec -= depth * np.exp(-0.5 * ((wave - centre) / width) ** 2)
    return spec


def write_file(path, *, sampling, rng):
    pix = 0.25
    grid = np.arange(-15, 16) * pix
    xx, yy = np.meshgrid(grid, grid)
    x = xx.ravel()
    y = yy.ravel()
    flux = surface_brightness(x, y)
    bin_id = block_binning(x, y)
    nbins = int(bin_id.max() + 1)

    bin_flux = np.bincount(bin_id, weights=flux, minlength=nbins)
    xbin = np.bincount(bin_id, weights=flux * x, minlength=nbins) / bin_flux
    ybin = np.bincount(bin_id, weights=flux * y, minlength=nbins) / bin_flux
    bin_snr = 25.0 * np.sqrt(bin_flux / np.median(bin_flux))

    vel, sigma, h3, h4 = kinematics(xbin, ybin)
    # A handful of outer bins with spurious high dispersion (for the sigma mask demo).
    outer = np.argsort(np.hypot(xbin, ybin))[-4:]
    sigma[outer] = 420.0

    xvel = np.arange(1200.0, -1200.1, -50.0)  # BAYES-LOSVD stores a descending axis
    nvel = len(xvel)
    dv = 50.0

    ln_wave = np.linspace(np.log(8450.0), np.log(8760.0), 420)
    wave = np.exp(ln_wave)
    npix = len(wave)
    keep = ~(((wave > 8600) & (wave < 8612)) | (wave < 8458) | (wave > 8752))
    mask = np.flatnonzero(keep).astype(np.int64)
    continuum_slope = (wave - wave.mean()) / np.ptp(wave)

    spec_obs = np.empty((npix, nbins))
    sigma_obs = np.empty((npix, nbins))
    losvd_out = np.empty((nbins, 5, nvel))
    bestfit_out = np.empty((nbins, 5, npix))

    for b in range(nbins):
        true = np.clip(gauss_hermite(xvel, vel[b], sigma[b], h3[b], h4[b]), 0, None)
        true /= true.sum() * dv
        if sampling:
            noise_level = 0.12 * true.max() * 30.0 / bin_snr[b]
            # Noise follows the profile so the wings stay near zero (flat noise
            # would bias the LOSVD second moment at low S/N).
            envelope = np.sqrt(true / true.max())
            median = np.clip(true + rng.normal(0, 0.3 * noise_level, nvel) * envelope, 0, None)
            width = noise_level * (0.05 + 0.95 * envelope)
            rows = np.clip(median[None, :] + PERCENTILE_Z[:, None] * width[None, :], 0, None)
            rows[2] = median
        else:
            kernel = np.array([0.25, 0.5, 0.25])
            smooth = np.convolve(true, kernel, mode="same")
            rows = np.repeat(smooth[None, :], 5, axis=0)
        losvd_out[b] = rows

        def model_for(lv):
            weights = lv / max(lv.sum(), 1e-30)
            spec = np.zeros(npix)
            for w, v in zip(weights, xvel):
                if w > 0:
                    spec += w * template(wave * np.exp(-v / C_KMS))
            return spec * (1.0 + 0.04 * continuum_slope)

        clean = model_for(true)
        noise = 1.0 / bin_snr[b]
        spec_obs[:, b] = clean + rng.normal(0, noise, npix)
        spec_obs[~keep, b] += np.where(rng.random((~keep).sum()) < 0.2, rng.normal(0, 4 * noise, (~keep).sum()), 0)
        sigma_obs[:, b] = noise
        fit = model_for(rows[2])
        if sampling:
            spread = 0.25 * noise
            bestfit_out[b] = fit[None, :] + PERCENTILE_Z[:, None] * spread
        else:
            bestfit_out[b] = np.repeat(fit[None, :], 5, axis=0)

    with h5py.File(path, "w") as f:
        f.attrs["fit_type"] = "SP" if sampling else "RWregularized"
        f.attrs["fit_mode"] = "NUTS sampling (synthetic)" if sampling else "MAP optimization (synthetic)"
        if not sampling:
            f.attrs["regularization"] = "RW"
        f.attrs["description"] = "Synthetic demo data for bayeslosvd-viewer; not a real galaxy."
        g = f.create_group("in")
        comp = {"compression": "gzip", "compression_opts": 6}
        g.create_dataset("x", data=x)
        g.create_dataset("y", data=y)
        g.create_dataset("binID", data=bin_id)
        g.create_dataset("flux", data=flux.astype(np.float32))
        g.create_dataset("xbin", data=xbin)
        g.create_dataset("ybin", data=ybin)
        g.create_dataset("bin_flux", data=bin_flux)
        g.create_dataset("bin_snr", data=bin_snr)
        g.create_dataset("xvel", data=xvel)
        g.create_dataset("wave_obs", data=ln_wave)
        g.create_dataset("wave", data=wave)
        g.create_dataset("mask", data=mask)
        g.create_dataset("spec_obs", data=spec_obs.astype(np.float32), **comp)
        g.create_dataset("sigma_obs", data=sigma_obs.astype(np.float32), **comp)
        for key, value in {
            "nbins": nbins, "nvel": nvel, "npix_obs": npix, "npix": npix, "nmask": len(mask),
            "velscale": float(np.diff(ln_wave).mean() * C_KMS), "lmin": 8450.0, "lmax": 8760.0,
            "porder": 1, "ntemp": 1, "snr": 25.0,
        }.items():
            g.create_dataset(key, data=value)
        out = f.create_group("out")
        for b in range(nbins):
            ob = out.create_group(str(b))
            ob.create_dataset("losvd", data=losvd_out[b], **comp)
            ob.create_dataset("bestfit", data=bestfit_out[b].astype(np.float32), **comp)
    return nbins


def write_mge(path):
    lines = ["sigma_arcsec,q_obs,luminosity_Lsun"]
    lines += [f"{s},{q},{lum}" for s, q, lum in MGE]
    path.write_text("\n".join(lines) + "\n")


def main():
    out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent.parent / "demo"
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20261014)
    n1 = write_file(out_dir / "demo_disk_sampling_results.hdf5", sampling=True, rng=rng)
    n2 = write_file(out_dir / "demo_disk_regularized_results.hdf5", sampling=False, rng=rng)
    write_mge(out_dir / "demo_mge.csv")
    print(f"wrote {n1} + {n2} bins to {out_dir}")


if __name__ == "__main__":
    main()
