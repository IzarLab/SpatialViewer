# SpatialViewer

An R package providing an interactive Shiny application for exploring spatial transcriptomics data stored in Seurat objects.

---

## Installation

```r
install.packages("devtools")
devtools::install_github("IzarLab/SpatialViewer")
```

> **Note on Seurat:** Seurat requires system libraries (HDF5, libcurl) that may need separate installation. See the [official Seurat install guide](https://satijalab.org/seurat/articles/install) if installation fails.

---

## Quick Start

```r
library(SpatialViewer)
launch_spatial_viewer()
```

The app opens in your browser with a file picker. Select your Seurat RDS file, optionally a polygon file, and click **Load dataset**. To switch datasets, use the **Load Dataset** panel in the sidebar and click **Reload** — the app restarts automatically with the new file.

You can also pass a path directly to skip the picker:

```r
launch_spatial_viewer("path/to/your_seurat_object.rds")
```

---

## App Features

### Spatial View tab

- **Color by**: Choose any metadata column (categorical or continuous) or search by gene name to color cells by log-normalized expression.
- **Coordinate system**: Toggle between physical slide coordinates (X/Y) and a dimensionality reduction embedding (e.g. UMAP).
- **Display mode**: Show cells as centroids (points) or as polygon boundaries (when a polygon file is provided).
- **Cell type grouping**: Pin selected cell types to a fixed color while others are colored by the active variable.
- **Color scale lock**: Fix the color scale to a specific range for cross-variable comparisons.
- **Multi-gene mode**: Select multiple genes and color cells by "any expressed" (binary).

### Summary Statistics tab

Distribution plots for any metadata column or gene, with an optional grouping variable:

| Reference \ Secondary | None | Categorical | Continuous |
|---|---|---|---|
| Categorical | Bar chart | Stacked bar (% and counts) | Violin + Box |
| Continuous | Violin + Box | Violin + Box | — |

---

## Data Requirements

The input must be a Seurat object (≥ v5.0) saved as an `.rds` file with:

- **Spatial coordinates**: Two numeric columns in `obj@meta.data`. Defaults are `x_centroid` and `y_centroid`; override with `x_col` and `y_col` parameters.
- **Counts assay**: A counts matrix accessible via `obj[["RNA"]]$counts` (or the assay set via the `assay` parameter).
- **Dimensionality reduction** (optional): A reduction slot (e.g. `umap`) in `obj@reductions`. The app auto-detects a fallback if the named slot is absent.

**Optional polygon file**: A CSV or gzipped CSV (`.csv.gz`) with columns:

| Column | Description |
|---|---|
| `cell` | Cell identifier matching a metadata column in the Seurat object |
| `x_global_px` | X coordinate of polygon vertex |
| `y_global_px` | Y coordinate of polygon vertex |

---

## Function Reference

```r
launch_spatial_viewer(
  seurat_path      = NULL,   # if NULL, opens with interactive file picker
  exclude_vars     = character(0),
  exclude_patterns = character(0),
  x_col            = "x_centroid",
  y_col            = "y_centroid",
  reduction        = "umap",
  assay            = "RNA",
  polygon_path     = NULL,
  cell_id_col      = NULL,
  celltype_col     = NULL
)
```

---

## Configuration

For repeated use, copy the bundled `config.R` template, set your paths, and source it:

```r
file.copy(system.file("config.R", package = "SpatialViewer"), "~/my_config.R")
# edit ~/my_config.R, then:
source("~/my_config.R")
launch_spatial_viewer(
  seurat_path      = DATA_FILE,
  exclude_vars     = EXCLUDE_VARS,
  exclude_patterns = EXCLUDE_PATTERNS,
  x_col            = X_COL,
  y_col            = Y_COL,
  reduction        = REDUCTION,
  polygon_path     = POLYGON_FILE,
  cell_id_col      = CELL_ID_COL
)
```

Alternatively, run via `Rscript` from a terminal (useful for scripting or remote sessions):

```bash
Rscript $(Rscript -e "cat(system.file('run_app.R', package='SpatialViewer'))") \
  /path/to/data.rds /path/to/polygons.csv.gz
```

---

## Troubleshooting

**Seurat installation fails**

Seurat requires system libraries that may not be present. On Linux/macOS:

```bash
# Ubuntu/Debian
sudo apt-get install libhdf5-dev libcurl4-openssl-dev

# macOS (Homebrew)
brew install hdf5
```

Then retry `install.packages("Seurat")`. See the [Seurat install guide](https://satijalab.org/seurat/articles/install) for full details.

**"Column not found" or blank plot**

The coordinate columns (`x_col`, `y_col`) must exist in `obj@meta.data`. Check available columns:

```r
obj <- readRDS("your_object.rds")
colnames(obj@meta.data)
```

Pass the correct column names to `launch_spatial_viewer()`.

**Seurat v4 object**

The package requires Seurat ≥ 5.0. If you have a v4 object, update it first:

```r
obj <- UpdateSeuratObject(obj)
saveRDS(obj, "updated_object.rds")
```

**App is slow or crashes with large datasets**

For objects with more than ~500,000 cells or very large counts matrices:

- Subset to cells or features of interest before saving.
- Omit the polygon file — polygon rendering is the most memory-intensive feature.
- Run on a machine with at least 32 GB RAM for large datasets.
- On Windows: `memory.limit(size = 64000)` before launching.
