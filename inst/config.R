# SpatialViewer configuration
# Copy this file locally, edit the values below, then source() it before
# calling launch_spatial_viewer(). See README for details.

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

# Polygon overlay — set to NULL to disable.
# CSV or .csv.gz with columns: cell, x_global_px, y_global_px
POLYGON_FILE <- NULL  # e.g., "/path/to/polygons.csv.gz"
CELL_ID_COL  <- NULL  # obj@meta.data column matching the 'cell' column above

# Pre-selected cell type grouping column (NULL = user picks in the app)
CELLTYPE_COL <- NULL
