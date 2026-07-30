# Tanager Methane Attribution — Environment Setup & STAC Query
## Permian Basin Case Study

---

## 1. Environment Setup

### Prerequisites
- Python 3.10+
- conda or venv recommended

### Create & Activate Environment

```bash
# Using conda (recommended)
conda create -n tanager python=3.11
conda activate tanager

# OR using venv
python -m venv tanager_env
source tanager_env/bin/activate  # Linux/Mac
tanager_env\Scripts\activate     # Windows
```

### Install Dependencies

```bash
pip install \
  pystac-client \
  pystac \
  rasterio \
  xarray \
  h5py \
  numpy \
  pandas \
  geopandas \
  shapely \
  matplotlib \
  spectral \
  requests \
  jupyter
```

**What each package does:**

| Package | Purpose |
|---|---|
| `pystac-client` | Query the Tanager Open Data STAC catalog |
| `rasterio` | Read/write raster data |
| `xarray` + `h5py` | Load Tanager HDF5 data cubes |
| `geopandas` | Load facility shapefiles (EPA GHGRP, HIFLD) |
| `spectral` | Matched filter & hyperspectral analysis |
| `numpy` / `pandas` | Data manipulation |
| `matplotlib` | Visualization |

---

## 2. STAC Catalog Query

### 2a. Explore the Catalog Structure

```python
import pystac_client
import json

# Open the Tanager Open Data STAC catalog (static catalog)
catalog = pystac_client.Client.open(
    "https://storage.googleapis.com/open-cogs/planet-stac/catalog.json"
)

print(f"Catalog title: {catalog.title}")
print(f"Catalog description: {catalog.description}")

# List all available collections
print("\nAvailable collections:")
for collection in catalog.get_collections():
    print(f"  - {collection.id}: {collection.description}")
```

### 2b. Search for Permian Basin Scenes

```python
import pystac_client
import geopandas as gpd
from shapely.geometry import box

# Permian Basin bounding box (West TX / SE New Mexico)
# [min_lon, min_lat, max_lon, max_lat]
PERMIAN_BBOX = [-104.5, 30.5, -100.5, 33.5]

catalog = pystac_client.Client.open(
    "https://storage.googleapis.com/open-cogs/planet-stac/catalog.json"
)

# Search across all collections in the bbox
search = catalog.search(
    bbox=PERMIAN_BBOX,
    max_items=100,
)

items = list(search.items())
print(f"\nFound {len(items)} scene(s) in Permian Basin bbox\n")

for item in items:
    print(f"ID:          {item.id}")
    print(f"Date:        {item.datetime}")
    print(f"BBox:        {item.bbox}")
    print(f"Collection:  {item.collection_id}")
    print(f"Assets:      {list(item.assets.keys())}")
    print("-" * 60)
```

### 2c. Save Scene Inventory to CSV

```python
import pandas as pd

records = []
for item in items:
    records.append({
        "id": item.id,
        "datetime": item.datetime,
        "bbox_west": item.bbox[0],
        "bbox_south": item.bbox[1],
        "bbox_east": item.bbox[2],
        "bbox_north": item.bbox[3],
        "collection": item.collection_id,
        "assets": ", ".join(item.assets.keys()),
    })

df = pd.DataFrame(records)
df.to_csv("permian_tanager_scenes.csv", index=False)
print(df.to_string())
```

### 2d. Download a Scene's HDF5 Asset

```python
import requests
import os

def download_asset(item, asset_key="basic_radiance_hdf5", out_dir="./data/raw"):
    """Download a Tanager scene asset locally."""
    os.makedirs(out_dir, exist_ok=True)

    asset = item.assets.get(asset_key)
    if asset is None:
        print(f"Asset '{asset_key}' not found. Available: {list(item.assets.keys())}")
        return None

    url = asset.href
    filename = os.path.join(out_dir, f"{item.id}_{asset_key}.h5")

    if os.path.exists(filename):
        print(f"Already downloaded: {filename}")
        return filename

    print(f"Downloading {item.id} → {filename}")
    response = requests.get(url, stream=True)
    response.raise_for_status()

    with open(filename, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)

    print(f"  ✓ Done ({os.path.getsize(filename) / 1e6:.1f} MB)")
    return filename

# Download the first Permian scene found
if items:
    download_asset(items[0])
```

---

## 3. Load & Inspect a Tanager Scene

```python
import h5py
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

def load_tanager_scene(filepath):
    """Load a Tanager HDF5 file into an xarray Dataset."""
    with h5py.File(filepath, "r") as f:
        print("HDF5 structure:")
        f.visit(lambda name: print(f"  {name}"))

        # Core radiance cube: shape (lines, samples, bands)
        radiance = f["radiance"][:]           # adjust key if needed
        wavelengths = f["wavelengths"][:]     # band center wavelengths in nm
        lat = f["latitude"][:]
        lon = f["longitude"][:]

    ds = xr.Dataset(
        {"radiance": (["line", "sample", "band"], radiance)},
        coords={
            "wavelength_nm": ("band", wavelengths),
            "lat": (["line", "sample"], lat),
            "lon": (["line", "sample"], lon),
        }
    )
    return ds


def plot_rgb_quicklook(ds, out_path="quicklook_rgb.png"):
    """Plot a false-color RGB quicklook using SWIR/NIR/Red bands."""
    wl = ds["wavelength_nm"].values
    rad = ds["radiance"].values

    # Approximate band indices for RGB quicklook
    r_idx = np.argmin(np.abs(wl - 650))    # Red ~650 nm
    g_idx = np.argmin(np.abs(wl - 550))    # Green ~550 nm
    b_idx = np.argmin(np.abs(wl - 450))    # Blue ~450 nm

    rgb = np.stack([
        rad[:, :, r_idx],
        rad[:, :, g_idx],
        rad[:, :, b_idx]
    ], axis=-1)

    # Normalize
    rgb = (rgb - rgb.min()) / (rgb.max() - rgb.min() + 1e-9)
    rgb = np.clip(rgb, 0, 1)

    plt.figure(figsize=(10, 8))
    plt.imshow(rgb)
    plt.title("Tanager RGB Quicklook")
    plt.axis("off")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.show()
    print(f"Saved: {out_path}")
```

---

## 4. Recommended Project Directory Structure

```
tanager_methane/
├── data/
│   ├── raw/              # Downloaded .h5 Tanager scenes
│   ├── processed/        # Orthorectified, clipped outputs
│   └── reference/        # EPA GHGRP, HIFLD infrastructure shapefiles
├── notebooks/
│   ├── 01_stac_query.ipynb
│   ├── 02_scene_inspection.ipynb
│   ├── 03_matched_filter.ipynb
│   ├── 04_ime_quantification.ipynb
│   └── 05_facility_attribution.ipynb
├── scripts/
│   ├── stac_query.py
│   ├── matched_filter.py
│   ├── ime_quantification.py
│   └── attribution.py
├── outputs/
│   ├── plume_maps/
│   └── figures/
├── requirements.txt
└── README.md
```

---

## 5. Reference Data Sources to Download Next

| Dataset | Source | Use |
|---|---|---|
| EPA GHGRP facilities | https://ghgdata.epa.gov/ghgp/main.do | Ground-truth reported emissions |
| HIFLD O&G infrastructure | https://hifld-geoplatform.opendata.arcgis.com | Wellheads, compressors, pipelines |
| ERA5 wind reanalysis | https://cds.climate.copernicus.eu | Wind speed/direction for IME |
| Carbon Mapper plumes | https://data.carbonmapper.org | Pre-identified Tanager plume locations |

---

## Next Step: Step 3 — Methane Detection with Matched Filter

Once you have scenes loaded, the next notebook (`03_matched_filter.ipynb`) will apply a
matched filter to the SWIR bands (~2300 nm) to detect CH₄ absorption features.
