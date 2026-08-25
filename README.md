# Detecting Methane Point Sources in Tanager-1 Hyperspectral Imagery

A matched-filter CH₄ detection pipeline for Planet Tanager-1 open data, with a
validation chain for separating gas absorption from surface mineralogy.

**Suvam Patel, Oregon State University** - entry in progress for the Planet
Tanager Open Data Competition (deadline 2026-08-31).

---

## Headline finding

**Detection works. Automated spectral confirmation does not.**

On a scene with 12 independently cataloged in-frame plumes, the column-wise
matched filter locates all 12 to within 3-19 m, at a median scene percentile of
99.94. The same plumes cannot be told apart from background by correlation
against a HITRAN CH₄ spectrum: against a size- and seeding-matched null, confirmed
plumes score a **median r = 0.635 versus the null's 0.723**, and only 1 of 10
exceeds the null's 95th percentile.

The gap between those two sentences is the result. A matched filter can put a
real plume in the top 0.1% of a scene while the spectrum extracted at that
location carries no usable CH₄ line structure - so ranking cannot be confirmed
by shape, and the confirmation step this repository was built around does not
work at Tanager's SNR and 30-43 m GSD.

This is a negative result and it is not being softened. It is also, as far as the
repository can tell, the correct one: it comes from the positive control, not
from the scenes the method rejected.

---

## Status

| Component | Status |
|---|---|
| STAC query, scene loading, HDF5 handling | Working, verified |
| Matched filter sign convention (`t = −k⊙μ`) | Fixed and verified |
| HITRAN target construction | Verified stable across p, T, FWHM |
| Column-wise (push-broom) background estimation | Working; normalization bias fixed |
| Ledoit-Wolf per-column shrinkage | Working, replaces fixed 0.15 |
| Brightness-matched control | Working |
| **Positive control (detection)** | **Established - 12/12 plumes found, median 99.94th pct** |
| **Positive control (spectral confirmation)** | **Established as failing — plumes do not separate from a matched null** |
| Shape test against HITRAN | Implemented; **shown insufficient**, see above |
| IME quantification arithmetic | Verified against outputs |
| EPA facility attribution | **Broken** — notebook 05 still calls a dead endpoint |
| Notebooks 03-05 | **Superseded** — see Pipeline provenance |

---

## The positive control

Notebook `07_positive_control.ipynb`, scene `20241004_081921_28_4001` (Iraq,
Energy & Mining, 2024-10-04). Ground truth is Planet's **per-frame**
`ql_ch4_json` - 12 plume features, frame-specific by construction - cross-checked
against the Carbon Mapper catalog filtered to the scene date. (Unfiltered, that
bbox returns 53 plumes spanning multiple dates; 14 are same-day. The date filter
is not optional at producing infrastructure, where the same wellhead leaks
repeatedly.)

### Detection: passes

| Measure | Result |
|---|---|
| Plumes located in frame | 12 / 12 |
| Geolocation offset | 3-19 m (sub-pixel to half-pixel at 43.6 m GSD) |
| Median scene percentile | 99.94 |
| Above 99th percentile | 9 / 12 |
| Above 99.9th percentile | 7 / 12 |
| MF score range | 0.0069 - 0.2354 |

### Spectral confirmation: fails

The null is 200 background regions grown from above-99th-percentile seeds by the
same region-growing rule used on the plumes, excluded to >15 px from any known
plume - so it is matched to the plumes in both region size (median 12 px vs 16 px)
and in how the seed was selected. This matters: an earlier null drawn from
*single random pixels* gives a median r of 0.355, against which plumes look
separable. They are not. Matching the null to the plumes' size and seeding is
what removes the apparent separation.

| Measure | Confirmed plumes | Matched null |
|---|---|---|
| Median lag-optimized r vs HITRAN | 0.635 | 0.723 |
| 95th percentile | - | 0.796 |
| 99th percentile | - | 0.851 |
| Above null 95th pct | 1 / 10 | - |

Background clusters tightly at r = 0.72-0.80. Plumes spread from 0.41 to 0.85.
The plume distribution sits *below* the background distribution and is wider than
it - the correlation statistic is measuring continuum shape shared by everything
in the scene, and the plumes' real absorption adds enough spectral noise to pull
them off it.

Section 7 puts a price on it. Thresholds fitted to admit the known plumes
(`r ≥ −0.032`, peak offset `≤ 64.5 nm`, differential `≥ −0.0041`) pass 5 of 6
held-out plumes - and **68 of 200 random background patches, 34%**. There is no
setting of these thresholds that keeps the plumes and rejects the scene.

### Detection tracks column density, not emission rate

Two cataloged sources of essentially identical strength - 1150.3 and 1151.9 kg/hr
- score 0.0221 and 0.1648, a 7.5× spread. The filter responds to the CH₄ column
density in the pixel, which depends on wind speed, plume age, and how the plume
falls across a 43.6 m grid. Reported emission rate is a poor predictor of
detectability, and MF score should not be read as a proxy for source strength.

---

## Pipeline provenance - read before citing any number

**Notebooks 03-05 are superseded and do not reproduce the paper's results.** They
predate the module and were never ported. Specifically:

- **03** builds its target with `build_ch4_absorption_coeffs()`, a placeholder sum
  of Gaussians. Its own comment says "replace with HITRAN". Notebooks 06 and 07
  use `hitran_ch4_k()` from `src/tanager_diagnostics.py`. These are different
  detectors and their outputs are not comparable.
- **04** applies `CF = 8.0e-6 kg/m² per MF unit`. This factor is arbitrary,
  exceeds the range quoted elsewhere in the project's own documentation, and was
  misattributed to a source that contains no such value. Every emission rate
  downstream of it inherits that.
- **05** calls `ofmpub.epa.gov/frs_public2/...` at cell 4, which returns a 302
  with an empty body. The parser reads that as "no facilities found". A working
  replacement (`data.epa.gov/efservice/...`) is already used at cell 10 of the
  same notebook, so the fix is local.

The results in this README come from **06 and 07 only**. Treat 03-05 as a record
of how the analysis started, not as a pipeline that produces current numbers.

### Results known to be invalid

- The four Permian clusters totalling 29.5 kg CH₄/hr, and the later sign-fixed
  29-cluster / 214.9 kg/hr version. The validation chain rejects both: broad
  absorption peaking near 2340 nm against HITRAN's 2370 nm, no line structure,
  consistent with calcite in caliche soils. Carbon Mapper's own QC documentation
  lists soil false positives as a known artifact requiring analyst review.
  Note that after notebook 07, this rejection is on weaker ground than it looked:
  the shape test that produced it also rejects 11 of 12 genuine plumes.
- The "sub-threshold sources invisible to GHGRP" argument. It rested on the dead
  EPA endpoint above. That was a silent failure, not a finding.

---

## What is actually interesting here

1. **A naive matched filter on Tanager selects terrain, not gas.** Scene-wide
   background statistics fill the top of the score distribution with push-broom
   detector striping; per-column statistics remove it, leaving field boundaries
   and drainage patterns. Neither is methane.

2. **Correlation against a CH₄ spectrum is not sufficient evidence - in either
   direction.** A calcite feature reaches r ≈ 0.75 against HITRAN CH₄ because both
   curves rise across 2200-2350 nm. Notebook 07 closes the loop: genuine plumes
   reach a *lower* median r than matched background. The statistic does not
   separate the classes, so neither its high values nor its low ones carry
   information.

3. **A null has to be matched to the thing it is a null for.** Single random
   pixels give median r = 0.355; regions matched to the plumes in size and seeding
   give 0.723. The first null would have supported a positive claim. The
   difference between them is the entire result.

4. **STAC `notes` can describe the collect rather than the frame.** Scene
   `20241121_183741_33_4001` is annotated "Plume A is a controlled release", but
   the release lies ~0.7 km outside the frame footprint, and Planet's own
   `ql_ch4_json` for the frame is an empty FeatureCollection. Nothing in the item
   metadata flags the discrepancy. Notebook 06 Section 6 documents the check.

---

## Layout

```
notebooks/
  01_stac_query.ipynb          catalog query, scene download
  02_scene_inspection.ipynb    HDF5 loading, RGB quicklook
  03_matched_filter.ipynb      CH4 detection          [superseded: Gaussian target]
  04_ime_quantification.ipynb  ERA5 wind, IME rates   [superseded: arbitrary CF]
  05_facility_attribution.ipynb EPA GHGRP join        [superseded: dead endpoint]
  06_validation.ipynb          method characterization, Casa Grande out-of-frame audit
  07_positive_control.ipynb    positive control on the Iraq scene  <- current result
src/
  tanager_diagnostics.py       detection + validation module
data/raw/                      Tanager HDF5 (gitignored, ~0.6 GB per scene)
data/reference/                ERA5 GRIB, HITRAN cache (gitignored, regenerable)
data/processed/                masks, CSVs, cached MF maps
outputs/figures/
```

`src/tanager_diagnostics.py` provides scene I/O (`load_scene`, `swir_matrix`),
target construction (`hitran_ch4_k`, `build_target`), detection
(`matched_filter`, `matched_filter_columnwise`, `ledoit_wolf_cov`,
`score_at_coords`), and validation (`brightness_matched_ratio`,
`hitran_shape_test`, `_verdict`, `run_validation`).

**Known gap:** the three functions carrying notebook 07's central result —
`grow()` (region growing from a seed), `shape_test_windowed()` (windowed,
lag-optimized correlation), and `mask_radius()` — exist only as notebook cells,
and `grow()` is defined twice. They are not importable, not tested, and not
reusable. They belong in `tanager_diagnostics.py`; see Open questions.

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

Notebook 07 caches its matched-filter map to
`data/processed/mf_columnwise_20241004_081921_28_4001.npy` (3.2 MB, committed), so
a rerun reproduces the analysis without recomputing the filter. It does still need
the raw scene for the spectral sections.

**Known environment issue:** on python.org builds, `pystac` uses `urllib`, which
skips certifi and raises `SSLCertVerificationError`. Set `SSL_CERT_FILE` to
`certifi.where()`. Conda builds are unaffected.

---

## Open questions - feedback especially welcome here

1. **What should replace the shape test?** This is no longer "are the thresholds
   right". Notebook 07 shows the statistic itself does not separate genuine plumes
   from matched background, so retuning `r`, peak offset, and differential cannot
   fix it - thresholds loose enough to keep the plumes admit 34% of random
   background. The open question is whether *any* per-region spectral test is
   viable at this SNR, or whether confirmation has to come from spatial structure
   (plume morphology, wind alignment, downwind persistence) or from repeat
   acquisition instead. Published practice is manual analyst review, which is
   consistent with there being no such statistic.

2. **Does the null result hold on other scenes?** It currently rests on one scene,
   one sensor geometry (43.6 m GSD, 29.9° off-nadir), and one surface type. The
   candidate replications are `20250226_065447_19_4001` (8 in-frame plumes) and
   `20250423_134021_00` / `134026_31` (Brazil, 1-2). If plumes separate from a
   matched null on a nadir, low-haze, dark-surface scene, the result is a
   statement about acquisition conditions rather than about Tanager.

3. **Is the result an artifact of the lag search?** `shape_test_windowed()` takes
   the maximum correlation over lags of ±8 bands (±40 nm). That is a generous
   search, and it is applied to the null as well as to the plumes - which is
   exactly what lets background reach a median of 0.723. Rerunning the comparison
   with `max_lag=0` would show whether the lag freedom is inflating the null more
   than it helps the plumes. If the separation reappears at zero lag, the
   conclusion narrows to "the lag-optimized statistic fails" rather than "the
   correlation statistic fails".

4. **Move the notebook-only analysis functions into the module.** `grow`,
   `shape_test_windowed`, and `mask_radius` produce the central result and live in
   notebook cells, one of them duplicated. Nothing else can import them and
   nothing tests them.

5. **Two targets, one pipeline.** Notebooks 03-05 still use the Gaussian
   placeholder while 06-07 use the HITRAN target from the module. Porting 03-05
   onto `tanager_diagnostics` would also fix the dead EPA endpoint and remove the
   unsupported conversion factor.

6. **Framing.** A characterized negative result - detection validated, automated
   confirmation shown to fail, with the matched-null methodology that demonstrates
   it - may be a stronger submission than a detection claim. The competition
   accepts a script or methodology adaptation, not only a case study. Views
   welcome.

## Data and licensing

Tanager-1 imagery: Planet Labs PBC, CC BY 4.0, via the Tanager Open Data STAC
catalog. ERA5: Copernicus Climate Change Service. Facility data: EPA Envirofacts.
Plume catalog: Carbon Mapper.
