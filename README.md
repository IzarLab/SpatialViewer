# SpatialViewer

An R package providing an interactive Shiny application for exploring spatial transcriptomics data stored in Seurat objects.

---

## Installation

**Requirements:** R ≥ 4.4.0

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

The app opens in your browser with a file picker. Select your Seurat RDS file, optionally a polygon file, and click **Load dataset**. Expand **Advanced settings** to change coordinate columns, assay, reduction, or polygon cell ID column before loading. To switch datasets, use the **Load Dataset** panel in the sidebar and click **Reload** — the app restarts with the new file.

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
- **Layer coloring**: Assign a different continuous variable (or fixed color) to each cell type group simultaneously — e.g. color T cells by activation score and tumor cells by IFN-γ response in the same view. Each layer supports its own palette and optional scale lock. Works in both centroid and polygon display modes.
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

- **Spatial coordinates**: Two numeric columns in `obj@meta.data`. Defaults are `x_slide_mm` and `y_slide_mm`; change in **Advanced settings** or via the `x_col`/`y_col` parameters.
- **Counts assay**: A counts matrix accessible via `obj[["RNA"]]$counts` (or the assay set in **Advanced settings** / via the `assay` parameter).
- **Dimensionality reduction** (optional): A reduction slot (e.g. `umap`) in `obj@reductions`. The app auto-detects a fallback if the named slot is absent.

**Optional polygon file**: A CSV or gzipped CSV (`.csv.gz`) with columns:

| Column | Description |
|---|---|
| `cell` | Cell identifier matching a metadata column in the Seurat object |
| `x_global_px` | X coordinate of polygon vertex |
| `y_global_px` | Y coordinate of polygon vertex |

Set the matching metadata column in **Advanced settings** (Cell ID column) or via `cell_id_col`.

---

## Function Reference

```r
launch_spatial_viewer(
  seurat_path      = NULL,  # if NULL, opens with interactive file picker
  exclude_vars     = NULL,  # character vector of column names to hide (programmatic only)
  exclude_patterns = NULL,  # substrings to match and hide; also in Advanced settings
  x_col            = NULL,  # default "x_slide_mm"; also in Advanced settings
  y_col            = NULL,  # default "y_slide_mm"; also in Advanced settings
  reduction        = NULL,  # default "umap";       also in Advanced settings
  assay            = NULL,  # default "RNA";        also in Advanced settings
  continuous_pals  = NULL,  # additional continuous palette names to add to the dropdown
  polygon_path     = NULL,  # path to polygon CSV or .csv.gz
  cell_id_col      = NULL,  # metadata column matching polygon 'cell'; also in Advanced settings
  celltype_col     = NULL,  # pre-selected cell type grouping column
  config_path      = NULL   # path to a config.R file; auto-detects "config.R" in working directory
)
```

Explicit arguments always take priority over config file values, which in turn take priority over built-in defaults.

---

## Configuration

For repeated use with a specific dataset, copy the bundled `config.R` template into your working directory and edit it:

```r
file.copy(system.file("config.R", package = "SpatialViewer"), "config.R")
```

Then launch normally — the file is detected automatically:

```r
launch_spatial_viewer()
```

To use a config file stored elsewhere:

```r
launch_spatial_viewer(config_path = "~/datasets/lung/config.R")
```

Explicit arguments always override config values, so you can mix both:

```r
launch_spatial_viewer(x_col = "x_global_mm")  # everything else from config.R
```

---

## Troubleshooting

**R version too old**

The `Matrix` package (a dependency of Seurat) requires R ≥ 4.4.0. If you see an error like `package 'Matrix' requires R version 4.4.0 or higher`, update R from [cran.r-project.org](https://cran.r-project.org/) before installing.

---

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

Set the correct names in **Advanced settings** before loading, or pass them as arguments:

```r
launch_spatial_viewer(x_col = "x_slide_mm", y_col = "y_slide_mm")
```

**Seurat v4 object**

The package requires Seurat ≥ 5.0. If you have a v4 object, update it first:

```r
obj <- UpdateSeuratObject(obj)
saveRDS(obj, "updated_object.rds")
```

**App is slow or crashes with large datasets (it runs locally and performance depends on the machine)**

For objects with more than ~500,000 cells or very large counts matrices:

- Subset to cells or features of interest before saving.
- Omit the polygon file — polygon rendering is the most memory-intensive feature.
- Run on a machine with at least 32 GB RAM for large datasets.
- On Windows: `memory.limit(size = 64000)` before launching.
