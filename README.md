# Detecting Methane Point Sources in Tanager-1 Hyperspectral Imagery

A matched-filter CH₄ detection pipeline for Planet Tanager-1 open data, with a
validation chain for separating gas absorption from surface mineralogy.

**Suvam Patel, Oregon State University** - entry in progress for the Planet
Tanager Open Data Competition (deadline 2026-08-31).

---

## Status: work in progress, seeking review

**The detection method is not yet validated. No headline result in this repo
should be cited.** This is stated up front because the repository is public for
mentor feedback, and several earlier results were wrong in ways that took a
while to find.

| Component | Status |
|---|---|
| STAC query, scene loading, HDF5 handling | Working, verified |
| Matched filter sign convention (`t = −k⊙μ`) | Fixed and verified |
| HITRAN target construction | Verified stable across p, T, FWHM |
| Column-wise (push-broom) background estimation | Implemented, bias fixed |
| Brightness-matched control | Working |
| Shape test against HITRAN | Working; thresholds unsourced |
| IME quantification arithmetic | Verified against outputs |
| EPA facility attribution | **Broken** — dead endpoint returns empty |
| **Positive control** | **Never established** — see Open Questions |

### Results known to be invalid

- The four Permian clusters totalling 29.5 kg CH₄/hr, and the later sign-fixed
  29-cluster / 214.9 kg/hr version. The validation chain rejects both: broad
  absorption peaking near 2340 nm against HITRAN's 2370 nm, no line structure,
  consistent with calcite in caliche soils. Carbon Mapper's own QC documentation
  lists soil false positives as a known artifact requiring analyst review.
- The "sub-threshold sources invisible to GHGRP" argument. It rested on an EPA
  FRS endpoint returning a 302 with an empty body, which the parser read as
  "no facilities found". That was a silent failure, not a finding.
- The conversion factor `CF = 8.0e-6 kg/m² per MF unit` is arbitrary, exceeds
  the range quoted in the project's own documentation, and was misattributed.

---

## What is actually interesting here

Three findings that survive scrutiny, and are the likely basis of the submission:

1. **A naive matched filter on Tanager selects terrain, not gas.** Scene-wide
   background statistics fill the top of the score distribution with push-broom
   detector striping; per-column statistics remove it, leaving field boundaries
   and drainage patterns. Neither is methane.

2. **Correlation against a CH₄ spectrum is not sufficient evidence.** A calcite
   absorption feature reaches r ≈ 0.75 against HITRAN CH₄ over 2050–2450 nm,
   because both curves rise across 2200–2350 nm. Peak position and line
   structure discriminate where correlation does not.

3. **STAC `notes` can describe the collect rather than the frame.** Scene
   `20241121_183741_33_4001` is annotated "Plume A is a controlled release", but
   the release lies ~0.7 km outside the frame footprint. Planet's own
   `ql_ch4_json` for the frame is an empty FeatureCollection. Nothing in the item
   metadata flags the discrepancy. Notebook 06 documents the check.

---

## Layout

```
notebooks/
  01_stac_query.ipynb          catalog query, scene download
  02_scene_inspection.ipynb    HDF5 loading, RGB quicklook
  03_matched_filter.ipynb      CH4 detection
  04_ime_quantification.ipynb  ERA5 wind, IME emission rates
  05_facility_attribution.ipynb EPA GHGRP join    [endpoint broken]
  06_validation.ipynb          method characterization, control-scene audit
  07_positive_control.ipynb    [not yet written]
src/
  tanager_diagnostics.py       detection + validation module
data/raw/                      Tanager HDF5 (gitignored, ~0.6 GB per scene)
data/reference/                ERA5 GRIB, HITRAN cache (gitignored, regenerable)
data/processed/                masks, CSVs
outputs/figures/
```

## Reproducing

```bash
conda create -n tanager python=3.11 && conda activate tanager
pip install -r requirements.txt
```

Requires credentials for ERA5 (free Copernicus account, `~/.cdsapirc`, and the
dataset licence must be accepted on the CDS site). HITRAN line data downloads on
first use and caches. Carbon Mapper's catalog API needs no auth.

Run notebooks from `notebooks/`; paths are relative (`../data`, `../outputs`).
Raw scenes are gitignored — notebook 01 downloads them.

**Known environment issue:** on python.org builds, `pystac` uses `urllib`, which
skips certifi and raises `SSLCertVerificationError`. Set `SSL_CERT_FILE` to
`certifi.where()`. Conda builds are unaffected.

---

## Open questions - feedback especially welcome here

1. **Positive control.** Nothing has yet shown the chain would *accept* a real
   plume. Next step is scene `20241004_081921_28_4001` (Iraq, 12 cataloged
   in-frame plumes): measure what scene percentile known plumes reach, and
   calibrate thresholds on held-out plumes rather than the evaluation set.

2. **Decision thresholds are unsourced.** `_verdict()` uses r > 0.80, peak offset
   < 15 nm, albedo tolerance 0.03, differential > 0.005. There is no literature
   standard - published practice is manual analyst review. The correlation is
   continuum-inflated over the full window; computing it over ~2150–2450 nm would
   help. The peak metric is a single argmax and is noise-fragile. Bootstrapped
   null distributions would be more defensible than fiat.

3. **Two targets, one pipeline.** Notebooks 03–05 still use a Gaussian
   placeholder while 06 uses the HITRAN target from the module. 03–05 should be
   ported onto `tanager_diagnostics`.

4. **Framing.** A validated-negative result plus a reusable validation chain may
   be a stronger submission than a detection claim. The competition accepts a
   script or methodology adaptation, not only a case study. Views welcome.

## Data and licensing

Tanager-1 imagery: Planet Labs PBC, CC BY 4.0, via the Tanager Open Data STAC
catalog. ERA5: Copernicus Climate Change Service. Facility data: EPA Envirofacts.
Plume catalog: Carbon Mapper.