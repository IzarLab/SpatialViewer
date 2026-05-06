options(shiny.launch.browser = TRUE)

library(SpatialViewer)

config_path <- system.file("config.R", package = "SpatialViewer")
if (nzchar(config_path)) {
  source(config_path)
} else {
  source("config.R")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[1])) {
  DATA_FILE <- args[1]
  message("Loading dataset from CLI argument: ", DATA_FILE)
}
if (length(args) >= 2 && nzchar(args[2])) {
  POLYGON_FILE <- args[2]
  message("Loading polygon file from CLI argument: ", POLYGON_FILE)
}

# In-app reload overrides everything (user explicitly chose a new file)
reload_path_file <- file.path(tempdir(), ".spatial_viewer_path")
if (file.exists(reload_path_file)) {
  DATA_FILE <- readLines(reload_path_file, n = 1)
  file.remove(reload_path_file)
  message("Loading dataset from reload request: ", DATA_FILE)
}

reload_poly_file <- file.path(tempdir(), ".spatial_viewer_polypath")
if (file.exists(reload_poly_file)) {
  poly_path <- readLines(reload_poly_file, n = 1)
  file.remove(reload_poly_file)
  if (nzchar(poly_path)) {
    POLYGON_FILE <- poly_path
    message("Loading polygon file from reload request: ", POLYGON_FILE)
  } else {
    POLYGON_FILE <- NULL
    message("Polygon file cleared via reload request")
  }
}

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
