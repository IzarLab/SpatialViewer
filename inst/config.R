# SpatialViewer configuration
# Place this file in your working directory — launch_spatial_viewer() will
# detect and load it automatically. Or pass config_path explicitly:
#   launch_spatial_viewer(config_path = "path/to/this/config.R")
# Any value set here can be overridden by the Advanced settings panel in the
# app or by passing an argument directly to launch_spatial_viewer().

DATA_FILE <- "/path/to/your_seurat_object.rds"

# Metadata columns to hide from the display dropdowns.
# All other obj@meta.data columns appear automatically.
EXCLUDE_VARS <- c(
  # "fov", "Area", "AspectRatio"
)

# Also hide any column whose name contains one of these substrings.
EXCLUDE_PATTERNS <- character(0)

# Spatial coordinate columns in obj@meta.data
X_COL <- "x_slide_mm"
Y_COL <- "y_slide_mm"

# Dimensionality reduction slot to show (auto-detects a fallback if absent)
REDUCTION <- "umap"

# Seurat assay name to read counts from
ASSAY <- "RNA"

# Polygon overlay — set to NULL to disable.
# CSV or .csv.gz with columns: cell, x_global_px, y_global_px
POLYGON_FILE <- NULL  # e.g., "/path/to/polygons.csv.gz"
CELL_ID_COL  <- NULL  # obj@meta.data column matching the 'cell' column above

# Pre-selected cell type grouping column (NULL = user picks in the app)
CELLTYPE_COL <- NULL
