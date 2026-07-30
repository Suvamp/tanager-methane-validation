"""
tanager_diagnostics.py
======================

Reusable CH4 detection + validation chain for Planet Tanager-1 hyperspectral
scenes.

The point of this module is that detection alone is not evidence. A matched
filter will happily flag dark soil, shadow, or carbonate mineralogy as a
"plume". This module bundles detection with three independent checks that a
false positive has to survive:

    1. Sign convention      - methane absorbs, so the target is negative and
                              radiance-weighted (t = -k * mu). An unweighted
                              positive target preferentially selects dark pixels.
    2. Brightness matching  - compare flagged pixels against background pixels of
                              similar continuum brightness, not the scene mean.
                              Removes the albedo confound.
    3. HITRAN shape test    - correlate the observed absorption depth against a
                              real CH4 spectrum convolved to Tanager resolution.
                              Catches surface minerals that absorb in the same
                              broad region (calcite ~2340 nm is the common one).

Usage
-----
    from tanager_diagnostics import run_validation

    result = run_validation("../data/raw/<scene>.h5")
    print(result.verdict)

Author: Suvam Patel, Oregon State University
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import h5py
import numpy as np

# numpy 2.0 renamed trapz -> trapezoid; support both
_trap = getattr(np, "trapezoid", None) or np.trapz

# ── HDF5 paths inside a Tanager Basic Radiance file ────────────────────────
_TOA = "HDFEOS/SWATHS/HYP/Data Fields/toa_radiance"
_LAT = "HDFEOS/SWATHS/HYP/Geolocation Fields/Latitude"
_LON = "HDFEOS/SWATHS/HYP/Geolocation Fields/Longitude"

FILL = -9999.0

# CH4 absorption window used for detection
WL_MIN, WL_MAX = 2050.0, 2450.0
# Sub-window where CH4 absorption is concentrated
BAND_LO, BAND_HI = 2200.0, 2400.0
# Tanager spectral resolution (FWHM, nm) - from toa_radiance 'fwhm' attribute
FWHM_NM = 5.5


# ══════════════════════════════════════════════════════════════════════════
# Loading
# ══════════════════════════════════════════════════════════════════════════

def load_scene(path):
    """
    Read a Tanager Basic Radiance HDF5 file.

    Returns (rad, wl, lat, lon) where rad is (bands, lines, samples) with
    fill values converted to NaN. Note bands-first ordering.
    """
    path = Path(path)
    with h5py.File(path, "r") as f:
        toa = f[_TOA]
        rad = toa[:].astype(np.float32)
        wl = np.asarray(toa.attrs["wavelengths"], dtype=np.float64)
        lat = f[_LAT][:]
        lon = f[_LON][:]

    rad[rad == FILL] = np.nan
    return rad, wl, lat, lon


def swir_matrix(rad, wl, wl_min=WL_MIN, wl_max=WL_MAX):
    """
    Subset to the SWIR window and reshape to (n_pixels, n_bands).

    Returns (X, wl_swir, valid, shape) where valid marks pixels with no NaNs
    and shape is the original (lines, samples) for reshaping results back.
    """
    sel = (wl >= wl_min) & (wl <= wl_max)
    wl_swir = wl[sel]
    sub = rad[sel, :, :]
    n_bands, n_lines, n_samples = sub.shape

    X = sub.reshape(n_bands, -1).T
    valid = ~np.any(np.isnan(X), axis=1)
    return X, wl_swir, valid, (n_lines, n_samples)


# ══════════════════════════════════════════════════════════════════════════
# Target construction
# ══════════════════════════════════════════════════════════════════════════

def hitran_ch4_k(wl_swir, cache_dir="../data/reference/hitran",
                 p_atm=1.0, T_K=288.0, fwhm_nm=FWHM_NM):
    """
    HITRAN-derived CH4 absorption coefficient spectrum, convolved to Tanager's
    instrument response and resampled onto wl_swir. Normalized to max 1.

    Requires `pip install hitran-api` and internet on first call; the line data
    is cached in cache_dir afterward.

    Falls back to None if hapi is unavailable, so callers can degrade to the
    Gaussian placeholder rather than crashing.
    """
    try:
        from hapi import db_begin, fetch, absorptionCoefficient_Voigt
    except ImportError:
        return None

    Path(cache_dir).mkdir(parents=True, exist_ok=True)
    db_begin(str(cache_dir))

    # Wavenumber range covering wl_swir, with margin for the convolution wings
    nu_lo = 1e7 / (wl_swir.max() + 10 * fwhm_nm)
    nu_hi = 1e7 / (wl_swir.min() - 10 * fwhm_nm)

    if not (Path(cache_dir) / "CH4.data").exists():
        fetch("CH4", 6, 1, nu_lo, nu_hi)

    nu, coef = absorptionCoefficient_Voigt(
        SourceTables="CH4",
        Environment={"p": p_atm, "T": T_K},
        WavenumberStep=0.01,
        HITRAN_units=False,
    )

    wl_h = 1e7 / np.asarray(nu)
    order = np.argsort(wl_h)
    wl_h, k_h = wl_h[order], np.asarray(coef)[order]

    return _convolve_to_bands(wl_h, k_h, wl_swir, fwhm_nm)


def _convolve_to_bands(wl_fine, k_fine, wl_bands, fwhm_nm=FWHM_NM):
    """Gaussian-convolve a fine spectrum onto instrument band centers."""
    sigma = fwhm_nm / 2.355
    W = np.exp(-0.5 * ((wl_fine[None, :] - wl_bands[:, None]) / sigma) ** 2)
    num = _trap(W * k_fine[None, :], wl_fine, axis=1)
    den = _trap(W, wl_fine, axis=1)
    k = num / den
    return k / k.max()


def gaussian_ch4_k(wl_swir):
    """
    Crude 3-Gaussian stand-in for the CH4 absorption spectrum.

    PLACEHOLDER ONLY. Retained so the pipeline runs without hapi, and so the
    paper can show what a naive target produces. Do not use for any claim -
    the band centers are approximate and the shape has no line structure,
    which is precisely what lets surface mineralogy pass as gas.
    """
    k = np.zeros(len(wl_swir))
    k += 1.0 * np.exp(-0.5 * ((wl_swir - 2304) / 40) ** 2)
    k += 0.4 * np.exp(-0.5 * ((wl_swir - 2204) / 25) ** 2)
    k += 0.3 * np.exp(-0.5 * ((wl_swir - 2380) / 20) ** 2)
    return k / k.max()


def build_target(k, mu):
    """
    Matched filter target: t = -k * mu, unit-normalized.

    Negative because methane absorbs. Multiplied by the background mean because
    the radiance perturbation from a plume scales with the brightness of the
    surface underneath it - an unweighted target biases toward dark pixels.
    """
    t = -k * mu
    return t / (np.linalg.norm(t) + 1e-12)


# ══════════════════════════════════════════════════════════════════════════
# Detection
# ══════════════════════════════════════════════════════════════════════════

def background_stats(X_valid, reg_frac=1e-4):
    """Background mean and regularized inverse covariance."""
    mu = X_valid.mean(axis=0)
    cov = np.cov((X_valid - mu).T)
    n = cov.shape[0]
    cov = cov + (reg_frac * np.trace(cov) / n) * np.eye(n)
    return mu, np.linalg.inv(cov)


def matched_filter(X, mu, cov_inv, t):
    """MF(x) = t' S^-1 (x - mu) / (t' S^-1 t). Positive => absorption."""
    Si_t = cov_inv @ t
    return ((X - mu) @ Si_t) / (t @ Si_t + 1e-12)


def threshold_mask(mf_map, n_sigma=2.0, smooth_sigma=1.5):
    """Threshold the MF map at mean + n_sigma*std, with light smoothing."""
    from scipy.ndimage import gaussian_filter

    m, s = np.nanmean(mf_map), np.nanstd(mf_map)
    thr = m + n_sigma * s
    smooth = gaussian_filter(np.nan_to_num(mf_map, nan=0.0), sigma=smooth_sigma)
    return smooth > thr, thr


# ══════════════════════════════════════════════════════════════════════════
# Validation
# ══════════════════════════════════════════════════════════════════════════

def brightness_matched_ratio(X, valid, plume_flat, wl_swir, k_neighbors=25):
    """
    Spectral ratio of flagged pixels to brightness-matched background pixels.

    Each flagged pixel is paired with the k nearest background pixels by
    continuum brightness (measured outside the CH4 bands). This removes the
    albedo confound: without it, any dark patch produces a depressed ratio at
    every wavelength.

    Returns (ratio, outside_mean, inside_mean).
    """
    from scipy.spatial import cKDTree

    cont = (wl_swir < BAND_LO - 50) | (wl_swir > BAND_HI + 20)
    bright = np.nanmean(X[:, cont], axis=1)

    p_idx = np.where(plume_flat & valid)[0]
    b_idx = np.where((~plume_flat) & valid)[0]

    if len(p_idx) == 0 or len(b_idx) < k_neighbors:
        return None, np.nan, np.nan

    tree = cKDTree(bright[b_idx][:, None])
    _, nn = tree.query(bright[p_idx][:, None], k=k_neighbors)
    matched = b_idx[np.ravel(nn)]

    ratio = np.nanmean(X[p_idx], axis=0) / np.nanmean(X[matched], axis=0)

    inside = (wl_swir >= BAND_LO) & (wl_swir <= BAND_HI)
    return ratio, float(np.nanmean(ratio[~inside])), float(np.nanmean(ratio[inside]))


def hitran_shape_test(ratio, k_ref, wl_swir):
    """
    Compare observed absorption depth against a reference CH4 spectrum.

    Returns (correlation, observed_peak_nm, reference_peak_nm).

    Correlation alone is not sufficient - two curves that both rise across
    2200-2350 nm will correlate even if their peaks are 40 nm apart. Always
    check the peak positions too.
    """
    depth = 1.0 - ratio
    r = float(np.corrcoef(depth, k_ref)[0, 1])
    return r, float(wl_swir[np.argmax(depth)]), float(wl_swir[np.argmax(k_ref)])


# ══════════════════════════════════════════════════════════════════════════
# Orchestration
# ══════════════════════════════════════════════════════════════════════════

@dataclass
class ValidationResult:
    scene: str
    n_plume_pixels: int = 0
    threshold: float = np.nan
    ratio_outside: float = np.nan
    ratio_inside: float = np.nan
    differential: float = np.nan
    hitran_r: float = np.nan
    observed_peak_nm: float = np.nan
    hitran_peak_nm: float = np.nan
    peak_offset_nm: float = np.nan
    used_hitran_target: bool = False
    verdict: str = ""
    arrays: dict = field(default_factory=dict, repr=False)

    def report(self):
        L = [
            "=" * 62,
            f"VALIDATION — {self.scene}",
            "=" * 62,
            f"  Flagged pixels           : {self.n_plume_pixels:,}",
            f"  MF threshold             : {self.threshold:.4f}",
            "",
            "  Brightness-matched control",
            f"    ratio outside CH4 bands: {self.ratio_outside:.3f}   (want ~1.00)",
            f"    ratio inside  CH4 bands: {self.ratio_inside:.3f}",
            f"    differential           : {self.differential:.4f}",
            "",
            "  HITRAN shape test",
            f"    correlation            : {self.hitran_r:.3f}",
            f"    observed peak          : {self.observed_peak_nm:.0f} nm",
            f"    HITRAN peak            : {self.hitran_peak_nm:.0f} nm",
            f"    offset                 : {self.peak_offset_nm:.0f} nm   (want <15)",
            "",
            f"  VERDICT: {self.verdict}",
            "=" * 62,
        ]
        return "\n".join(L)


def _verdict(res):
    """
    Decision rule, stated before looking at any particular scene.

    Thresholds are deliberately conservative. The failure mode being guarded
    against is a broad surface absorption feature passing as gas, which is why
    peak position is weighted as heavily as correlation.
    """
    if res.n_plume_pixels < 5:
        return "NO DETECTION - too few pixels above threshold"

    albedo_ok = abs(res.ratio_outside - 1.0) < 0.03
    if not albedo_ok:
        return ("REJECT - albedo confound; flagged pixels differ from background "
                "outside the CH4 bands, so brightness is driving selection")

    if res.differential < 0.005:
        return "REJECT - no meaningful absorption inside the CH4 bands"

    peak_ok = res.peak_offset_nm < 15
    corr_ok = res.hitran_r > 0.80

    if corr_ok and peak_ok:
        return "ACCEPT - consistent with CH4 in both shape and peak position"
    if peak_ok and res.hitran_r > 0.60:
        return "MARGINAL - peak aligns but shape correlation is weak; inspect manually"
    if corr_ok and not peak_ok:
        return (f"REJECT - correlation is driven by shared slope, not shared features; "
                f"peak is {res.peak_offset_nm:.0f} nm from HITRAN")
    return ("REJECT - absorption feature does not match CH4; check surface "
            "mineralogy (calcite ~2340 nm is common in caliche soils)")


def run_validation(scene_path, n_sigma=2.0, use_hitran=True,
                   hitran_cache="../data/reference/hitran", verbose=True):
    """
    Full detection + validation chain on one scene.

    Returns a ValidationResult. The `arrays` field carries mf_map, plume_mask,
    ratio, k_ref, and wl_swir for plotting.
    """
    scene_path = Path(scene_path)
    rad, wl, lat, lon = load_scene(scene_path)
    X, wl_swir, valid, shape = swir_matrix(rad, wl)

    mu, cov_inv = background_stats(X[valid])

    k_ref = hitran_ch4_k(wl_swir, cache_dir=hitran_cache) if use_hitran else None
    used_hitran = k_ref is not None
    if k_ref is None:
        k_ref = gaussian_ch4_k(wl_swir)

    t = build_target(k_ref, mu)

    mf_flat = np.full(X.shape[0], np.nan)
    mf_flat[valid] = matched_filter(X[valid], mu, cov_inv, t)
    mf_map = mf_flat.reshape(shape)

    plume_mask, thr = threshold_mask(mf_map, n_sigma=n_sigma)
    plume_flat = plume_mask.ravel()

    res = ValidationResult(
        scene=scene_path.stem,
        n_plume_pixels=int(plume_mask.sum()),
        threshold=float(thr),
        used_hitran_target=used_hitran,
    )

    ratio, out_m, in_m = brightness_matched_ratio(X, valid, plume_flat, wl_swir)
    if ratio is not None:
        res.ratio_outside = out_m
        res.ratio_inside = in_m
        res.differential = out_m - in_m
        r, obs_pk, ref_pk = hitran_shape_test(ratio, k_ref, wl_swir)
        res.hitran_r = r
        res.observed_peak_nm = obs_pk
        res.hitran_peak_nm = ref_pk
        res.peak_offset_nm = abs(obs_pk - ref_pk)

    res.verdict = _verdict(res)
    res.arrays = dict(mf_map=mf_map, plume_mask=plume_mask, ratio=ratio,
                      k_ref=k_ref, wl_swir=wl_swir, lat=lat, lon=lon)

    if verbose:
        print(res.report())
        if not used_hitran:
            print("  WARNING: hapi unavailable - used Gaussian placeholder target.")
            print("           Verdict is not trustworthy. pip install hitran-api")
    return res


def plot_validation(res, out_path=None):
    """Two-panel figure: brightness-matched ratio, and observed vs HITRAN shape."""
    import matplotlib.pyplot as plt

    a = res.arrays
    if a.get("ratio") is None:
        print("No ratio to plot.")
        return

    wl_s, ratio, k_ref = a["wl_swir"], a["ratio"], a["k_ref"]
    fig, ax = plt.subplots(1, 2, figsize=(14, 4.5))

    ax[0].plot(wl_s, ratio, lw=2)
    ax[0].axhline(1.0, color="k", ls="--", lw=0.8)
    ax[0].axvspan(BAND_LO, BAND_HI, alpha=0.15, color="orange",
                  label="CH4 absorption")
    ax[0].set_xlabel("Wavelength (nm)")
    ax[0].set_ylabel("Plume / brightness-matched background")
    ax[0].set_title(f"Brightness-matched ratio\n"
                    f"outside {res.ratio_outside:.3f} | inside {res.ratio_inside:.3f}")
    ax[0].legend()
    ax[0].grid(alpha=0.3)

    depth = 1.0 - ratio
    ax[1].plot(wl_s, depth / np.nanmax(depth), lw=2, label="observed")
    ax[1].plot(wl_s, k_ref, lw=2, label="HITRAN CH4"
               if res.used_hitran_target else "Gaussian placeholder")
    ax[1].axvline(res.observed_peak_nm, color="C0", ls=":", lw=1.2)
    ax[1].axvline(res.hitran_peak_nm, color="C1", ls=":", lw=1.2)
    ax[1].set_xlabel("Wavelength (nm)")
    ax[1].set_ylabel("Normalized absorption")
    ax[1].set_title(f"Shape test — r = {res.hitran_r:.3f}, "
                    f"peak offset {res.peak_offset_nm:.0f} nm")
    ax[1].legend()
    ax[1].grid(alpha=0.3)

    fig.suptitle(f"{res.scene}   —   {res.verdict}", fontsize=11, y=1.02)
    fig.tight_layout()
    if out_path:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {out_path}")
    plt.show()