#' Launch the Spatial Transcriptomics Viewer
#'
#' Opens an interactive Shiny application for exploring a Seurat spatial
#' transcriptomics object. Supports coloring by metadata or gene expression,
#' spatial and dimensionality reduction views, polygon overlays, and a
#' Summary Statistics tab.
#'
#' @param seurat_path Character or `NULL`. Path to a Seurat object saved as an
#'   `.rds` file. If `NULL` (the default), the app opens with a file picker so
#'   the dataset can be selected interactively.
#' @param exclude_vars Character vector. Metadata column names to omit from
#'   display dropdowns. Default: `character(0)`.
#' @param exclude_patterns Character vector. Metadata columns whose names contain
#'   any of these substrings are also excluded. Default: `character(0)`.
#' @param x_col Character. Metadata column for spatial X coordinates.
#'   Default: `"x_centroid"`.
#' @param y_col Character. Metadata column for spatial Y coordinates.
#'   Default: `"y_centroid"`.
#' @param reduction Character. Dimensionality reduction slot name (e.g. `"umap"`).
#'   Auto-detects a fallback if the named slot is absent. Default: `"umap"`.
#' @param assay Character. Seurat assay name to read counts from.
#'   Default: `"RNA"`.
#' @param continuous_pals Character vector of additional continuous palette names
#'   to add to the UI dropdown. Default: `NULL`.
#' @param polygon_path Character or `NULL`. Path to a CSV or gzipped CSV of
#'   cell boundary polygon vertices (columns: `cell`, `x_global_px`,
#'   `y_global_px`). Default: `NULL`.
#' @param cell_id_col Character or `NULL`. Seurat metadata column matching the
#'   `cell` column in the polygon file. Default: `NULL`.
#' @param celltype_col Character or `NULL`. Metadata column pre-selected for
#'   cell type grouping in the UI. Default: `NULL`.
#'
#' @return Called for its side effect: launches a blocking Shiny application.
#'
#' @import data.table
#'
#' @examples
#' \dontrun{
#' # Minimal usage
#' launch_spatial_viewer("path/to/seurat_object.rds")
#'
#' # Custom coordinates + polygon overlay
#' launch_spatial_viewer(
#'   seurat_path  = "path/to/seurat_object.rds",
#'   x_col        = "x_slide_mm",
#'   y_col        = "y_slide_mm",
#'   reduction    = "scpearson_umap_nobatch",
#'   polygon_path = "path/to/polygons.csv.gz",
#'   cell_id_col  = "cell_ID"
#' )
#' }
#'
#' @export
launch_spatial_viewer <- function(seurat_path      = NULL,
                                  exclude_vars     = NULL,
                                  exclude_patterns = NULL,
                                  x_col            = NULL,
                                  y_col            = NULL,
                                  reduction        = NULL,
                                  assay            = NULL,
                                  continuous_pals  = NULL,
                                  polygon_path     = NULL,
                                  cell_id_col      = NULL,
                                  celltype_col     = NULL,
                                  config_path      = if (file.exists("config.R")) "config.R" else NULL) {
  library(shiny)
  library(plotly)
  library(RColorBrewer)
  library(viridis)
  library(colourpicker)
  library(ggplot2)
  library(Seurat)
  library(data.table)
  library(shinyFiles)

  # Load config file and apply as defaults for any parameters not explicitly provided
  if (!is.null(config_path) && file.exists(config_path)) {
    cfg <- new.env(parent = emptyenv())
    source(config_path, local = cfg)
    if (is.null(seurat_path)      && exists("DATA_FILE",        envir = cfg)) seurat_path      <- cfg$DATA_FILE
    if (is.null(exclude_vars)     && exists("EXCLUDE_VARS",     envir = cfg)) exclude_vars     <- cfg$EXCLUDE_VARS
    if (is.null(exclude_patterns) && exists("EXCLUDE_PATTERNS", envir = cfg)) exclude_patterns <- cfg$EXCLUDE_PATTERNS
    if (is.null(x_col)            && exists("X_COL",            envir = cfg)) x_col            <- cfg$X_COL
    if (is.null(y_col)            && exists("Y_COL",            envir = cfg)) y_col            <- cfg$Y_COL
    if (is.null(reduction)        && exists("REDUCTION",        envir = cfg)) reduction        <- cfg$REDUCTION
    if (is.null(assay)            && exists("ASSAY",            envir = cfg)) assay            <- cfg$ASSAY
    if (is.null(polygon_path)     && exists("POLYGON_FILE",     envir = cfg)) polygon_path     <- cfg$POLYGON_FILE
    if (is.null(cell_id_col)      && exists("CELL_ID_COL",      envir = cfg)) cell_id_col      <- cfg$CELL_ID_COL
    if (is.null(celltype_col)     && exists("CELLTYPE_COL",     envir = cfg)) celltype_col     <- cfg$CELLTYPE_COL
  }

  # Hardcoded fallbacks (lowest priority — config and explicit args override these)
  if (is.null(exclude_vars))     exclude_vars     <- character(0)
  if (is.null(exclude_patterns)) exclude_patterns <- character(0)
  if (is.null(x_col))            x_col            <- "x_slide_mm"
  if (is.null(y_col))            y_col            <- "y_slide_mm"
  if (is.null(reduction))        reduction        <- "umap"
  if (is.null(assay))            assay            <- "RNA"
  if (is.null(cell_id_col))      cell_id_col      <- "cell_ID_new"


  repeat {
    .sv_path <- file.path(tempdir(), ".spatial_viewer_path")
    if (file.exists(.sv_path)) {
      seurat_path <- readLines(.sv_path, n = 1); file.remove(.sv_path)
    }
    .sv_poly <- file.path(tempdir(), ".spatial_viewer_polypath")
    if (file.exists(.sv_poly)) {
      pp <- readLines(.sv_poly, n = 1); file.remove(.sv_poly)
      polygon_path <- if (nzchar(pp)) pp else NULL
    }
    .sv_settings <- file.path(tempdir(), ".spatial_viewer_settings.rds")
    if (file.exists(.sv_settings)) {
      .cfg <- readRDS(.sv_settings); file.remove(.sv_settings)
      x_col            <- .cfg$x_col
      y_col            <- .cfg$y_col
      reduction        <- .cfg$reduction
      assay            <- .cfg$assay
      cell_id_col      <- .cfg$cell_id_col
      exclude_patterns <- .cfg$exclude_patterns
    }

    if (is.null(seurat_path)) {
      volumes_l <- c(Home = path.expand("~"), getVolumes()())
      land_ui <- fluidPage(
        tags$head(tags$style(HTML("
          .adv-section > summary {
            cursor: pointer; color: #888; font-size: 12px;
            list-style: none; padding: 2px 0;
          }
          .adv-section > summary::-webkit-details-marker { display: none; }
          .adv-section > summary:hover { color: #337ab7; }
          .adv-section .form-group { margin-bottom: 6px; }
          .adv-section label { font-size: 12px; color: #555; margin-bottom: 2px; }
          .adv-section input[type=’text’] { font-size: 12px; height: 28px; padding: 2px 6px; }
        "))),
        titlePanel("Spatial Transcriptomics Viewer"),
        sidebarLayout(
          sidebarPanel(
            tags$p(tags$strong("Select a dataset to get started.")),
            shinyFilesButton("file_pick", "Browse for Seurat RDS...",
                             title = "Select a Seurat RDS file", multiple = FALSE),
            tags$br(), tags$br(),
            verbatimTextOutput("current_file_path", placeholder = TRUE),
            tags$hr(),
            tags$p(style = "color:#666; font-size:12px;", "Polygon file (optional):"),
            shinyFilesButton("poly_file_pick", "Browse for polygon CSV...",
                             title = "Select a polygon CSV file", multiple = FALSE),
            tags$br(), tags$br(),
            verbatimTextOutput("current_poly_path", placeholder = TRUE),
            actionButton("clear_poly_btn", "Clear", icon = icon("times")),
            tags$hr(),
            tags$details(
              class = "adv-section",
              tags$summary("▶ Advanced settings"),
              tags$div(
                style = "padding: 8px 0 4px 0;",
                textInput("x_col_input",    "X coordinate column:",       value = x_col),
                textInput("y_col_input",    "Y coordinate column:",       value = y_col),
                textInput("reduction_input","Dimensionality reduction:",  value = reduction),
                textInput("assay_input",    "Seurat assay:",              value = assay),
                textInput("cell_id_col_input", "Cell ID column (for polygons):",
                          value = if (!is.null(cell_id_col)) cell_id_col else ""),
                textInput("exclude_patterns_input", "Hide columns matching (comma-separated):",
                          value = paste(exclude_patterns, collapse = ", "))
              )
            ),
            tags$br(),
            actionButton("reload_btn", "Load dataset", icon = icon("play"),
                         style = "color: white; background-color: #337ab7;"),
            width = 3
          ),
          mainPanel(
            tags$div(
              style = "padding: 80px 40px; text-align: center; color: #aaa;",
              tags$h3("No dataset loaded"),
              tags$p("Select a Seurat RDS file and click ‘Load dataset’.")
            ),
            width = 9
          )
        )
      )
      land_server <- function(input, output, session) {
        shinyFileChoose(input, "file_pick", roots = volumes_l,
                        filetypes = c("RDS", "rds"))
        shinyFileChoose(input, "poly_file_pick", roots = volumes_l,
                        filetypes = c("csv", "CSV", "gz"))
        sel_path <- reactiveVal("")
        sel_poly <- reactiveVal("")
        observeEvent(input$file_pick, {
          fi <- parseFilePaths(volumes_l, input$file_pick)
          if (nrow(fi) > 0) sel_path(as.character(fi$datapath))
        })
        observeEvent(input$poly_file_pick, {
          fi <- parseFilePaths(volumes_l, input$poly_file_pick)
          if (nrow(fi) > 0) sel_poly(as.character(fi$datapath))
        })
        observeEvent(input$clear_poly_btn, sel_poly(""))
        output$current_file_path <- renderText(sel_path())
        output$current_poly_path <- renderText({
          p <- sel_poly(); if (nzchar(p)) p else "(none)"
        })
        observeEvent(input$reload_btn, {
          p <- sel_path()
          if (nzchar(p) && file.exists(p)) {
            writeLines(p, file.path(tempdir(), ".spatial_viewer_path"))
            writeLines(sel_poly(), file.path(tempdir(), ".spatial_viewer_polypath"))
            raw_patterns <- trimws(strsplit(input$exclude_patterns_input, ",")[[1]])
            saveRDS(
              list(
                x_col            = input$x_col_input,
                y_col            = input$y_col_input,
                reduction        = input$reduction_input,
                assay            = input$assay_input,
                cell_id_col      = if (nzchar(input$cell_id_col_input)) input$cell_id_col_input else NULL,
                exclude_patterns = raw_patterns[nzchar(raw_patterns)]
              ),
              file.path(tempdir(), ".spatial_viewer_settings.rds")
            )
            stopApp()
          }
        })
      }
      runApp(shinyApp(land_ui, land_server))
      next
    }

  # Load Seurat object and extract metadata
  obj <- readRDS(seurat_path)
  df <- obj@meta.data

  # Derive display variables: all metadata columns minus exclusions
  all_meta <- colnames(df)
  all_meta <- setdiff(all_meta, c(x_col, y_col))
  if (length(exclude_patterns) > 0) {
    pattern_hits <- unique(unlist(lapply(exclude_patterns, function(p)
      all_meta[grepl(p, all_meta, fixed = TRUE)])))
    all_meta <- setdiff(all_meta, pattern_hits)
  }
  metadata_vars <- setdiff(all_meta, exclude_vars)

  # Add log10(nCount_RNA)
  if ("nCount_RNA" %in% metadata_vars) {
    df$log10_nCount_RNA <- log10(df$nCount_RNA + 1)
    metadata_vars <- c(metadata_vars, "log10_nCount_RNA")
  }

  # Extract dimensionality reduction embeddings (auto-detect if configured one is missing)
  available_reductions <- Reductions(obj)
  if (!reduction %in% available_reductions && length(available_reductions) > 0) {
    # Prefer reductions with "umap" in the name, then "tsne", then first available
    umap_hits <- grep("umap", available_reductions, ignore.case = TRUE, value = TRUE)
    tsne_hits <- grep("tsne", available_reductions, ignore.case = TRUE, value = TRUE)
    fallback <- if (length(umap_hits) > 0) umap_hits[1]
                else if (length(tsne_hits) > 0) tsne_hits[1]
                else available_reductions[1]
    message("Note: reduction '", reduction, "' not found. Using '", fallback, "' instead.")
    message("Available reductions: ", paste(available_reductions, collapse = ", "))
    reduction <- fallback
  }

  has_reduction <- tryCatch({
    emb <- as.data.frame(Embeddings(obj, reduction = reduction))
    colnames(emb)[1:2] <- c("Dim1", "Dim2")
    df$Dim1 <- emb$Dim1
    df$Dim2 <- emb$Dim2
    rm(emb)
    TRUE
  }, error = function(e) {
    message("Note: no usable reduction found. Dimensionality reduction view will be unavailable.")
    FALSE
  })

  # Extract and log-normalize counts matrix for gene expression visualization
  counts_raw <- Matrix::t(obj[[assay]]$counts)
  scale_factor <- mean(df$nCount_RNA)
  scale_row <- scale_factor / df$nCount_RNA
  counts_norm <- counts_raw
  if (!inherits(counts_norm, "dgCMatrix")) {
    counts_norm <- Matrix::Matrix(counts_norm, sparse = TRUE)
  }
  counts_norm@x <- counts_norm@x * scale_row[counts_norm@i + 1L]
  counts_norm@x <- log1p(counts_norm@x)
  gene_names <- sort(colnames(counts_norm))
  rm(counts_raw)

  rm(obj)

  # Load polygon data if provided
  has_polygons <- FALSE
  poly_dt <- NULL
  if (!is.null(polygon_path) && file.exists(polygon_path) && !is.null(cell_id_col)) {
    if (!cell_id_col %in% colnames(df)) {
      message("Warning: cell_id_col '", cell_id_col, "' not found in Seurat metadata. ",
              "Polygons disabled. Available columns: ",
              paste(colnames(df), collapse = ", "))
    } else {
      message("Loading polygon data from: ", polygon_path)
      poly_dt <- fread(polygon_path)
      valid_cells <- df[[cell_id_col]]
      poly_dt <- poly_dt[poly_dt[["cell"]] %in% valid_cells, ]
      if (nrow(poly_dt) > 0) {
        has_polygons <- TRUE
        setkey(poly_dt, "cell")
        poly_centroids <- poly_dt[, list(cx = mean(x_global_px), cy = mean(y_global_px)), by = "cell"]
        message("Polygon data loaded: ", length(unique(poly_dt[["cell"]])),
                " cells, ", nrow(poly_dt), " vertices")
      } else {
        message("Warning: No polygon cells matched Seurat cell IDs. Polygons disabled.")
        poly_dt <- NULL
      }
    }
  }

  # Validate required coordinate columns
  required_cols <- c(x_col, y_col)
  missing_required <- setdiff(required_cols, colnames(df))
  if (length(missing_required) > 0) {
    stop("Required coordinate columns not found: ",
         paste(missing_required, collapse = ", "))
  }

  # Drop any metadata_vars not present in the data, with a note
  missing_meta <- setdiff(metadata_vars, colnames(df))
  if (length(missing_meta) > 0) {
    message("Note: the following metadata columns were not found and will ",
            "be skipped: ", paste(missing_meta, collapse = ", "))
    metadata_vars <- setdiff(metadata_vars, missing_meta)
  }

  # Subset to only needed columns
  keep_cols <- c(x_col, y_col,
                 if (has_reduction) c("Dim1", "Dim2"),
                 if (has_polygons && !is.null(cell_id_col)) cell_id_col,
                 metadata_vars)
  keep_cols <- unique(keep_cols)
  df <- df[, keep_cols, drop = FALSE]

  # ---- Palette definitions ----

  # Master list of all supported continuous palettes
  all_palette_names <- c("Viridis", "Inferno", "Plasma", "Magma", "Cividis",
                          "Blues", "Greens", "Reds", "Purples", "Oranges",
                          "YlOrRd", "YlGnBu", "Hot", "Electric", "Jet", "Rainbow")

  # Subset that Plotly supports natively (pass name as string)
  plotly_native_scales <- c("Viridis", "Inferno", "Plasma", "Magma", "Cividis",
                             "Blues", "Greens", "Reds", "YlOrRd", "YlGnBu",
                             "Hot", "Electric", "Jet", "Rainbow")

  # Default choices shown in the UI

  cont_choices <- c("Inferno", "Viridis", "Plasma", "Magma",
                    "Blues", "Greens", "Reds", "Purples",
                    "YlOrRd", "Hot")
  if (!is.null(continuous_pals)) {
    matched <- sapply(continuous_pals, function(p) {
      idx <- match(tolower(p), tolower(all_palette_names))
      if (!is.na(idx)) all_palette_names[idx] else NA_character_
    })
    cont_choices <- unique(c(cont_choices, matched[!is.na(matched)]))
  }

  # Helper: generate n categorical colors from Set3 with ramp extension
  make_cat_colors <- function(n_levels) {
    base_pal <- brewer.pal(min(max(n_levels, 3), 12), "Set3")
    if (n_levels > 12) colorRampPalette(base_pal)(n_levels) else base_pal[seq_len(n_levels)]
  }

  # Helper: return a color-generating function for a continuous palette
  get_cont_pal_func <- function(pal_name) {
    switch(pal_name,
      "Viridis"  = viridis::viridis,
      "Inferno"  = viridis::inferno,
      "Plasma"   = viridis::plasma,
      "Magma"    = viridis::magma,
      "Cividis"  = viridis::cividis,
      "Blues"    = colorRampPalette(c("#f7fbff","#c6dbef","#6baed6","#2171b5","#084594")),
      "Greens"   = colorRampPalette(c("#f7fcf5","#c7e9c0","#74c476","#238b45","#005a32")),
      "Reds"     = colorRampPalette(c("#fff5f0","#fee0d2","#fc9272","#de2d26","#a50f15")),
      "Purples"  = colorRampPalette(c("#fcfbfd","#dadaeb","#9e9ac8","#6a51a3","#3f007d")),
      "Oranges"  = colorRampPalette(c("#fff5eb","#fdd0a2","#fdae6b","#e6550d","#8c2d04")),
      "YlOrRd"   = colorRampPalette(c("#ffffb2","#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026")),
      "YlGnBu"   = colorRampPalette(c("#ffffd9","#c7e9b4","#41b6c4","#225ea8","#081d58")),
      "Hot"      = colorRampPalette(c("#000000","#8b0000","#ff0000","#ff8c00","#ffff00","#ffffff")),
      "Electric" = colorRampPalette(c("#000000","#1a0099","#0000ff","#00ffff","#ffff00")),
      "Jet"      = colorRampPalette(c("#00007f","#0000ff","#00ffff","#7fff7f","#ffff00","#ff0000","#7f0000")),
      "Rainbow"  = colorRampPalette(c("#ff0000","#ff7f00","#ffff00","#00ff00","#0000ff","#8b00ff")),
      viridis::viridis  # fallback
    )
  }

  # Helper: return a Plotly colorscale — string for built-in names, [[value,color]] list for custom
  get_plotly_colorscale <- function(pal_name) {
    if (pal_name %in% plotly_native_scales) return(pal_name)
    fn <- get_cont_pal_func(pal_name)
    n <- 10
    cols <- fn(n)
    vals <- seq(0, 1, length.out = n)
    lapply(seq_along(vals), function(i) list(vals[i], cols[i]))
  }

  # UI
  ui <- fluidPage(
    tags$head(tags$style(HTML("
      .file-picker-section > summary {
        cursor: pointer;
        padding: 4px 6px;
        border-radius: 4px;
        font-weight: 600;
        color: #337ab7;
        list-style: none;
      }
      .file-picker-section > summary::-webkit-details-marker { display: none; }
      .file-picker-section > summary:hover {
        background-color: #f0f0f0;
      }
      .file-picker-section > summary .fa-caret-right {
        transition: transform 0.2s;
        margin-right: 4px;
      }
      .file-picker-section[open] > summary .fa-caret-right {
        transform: rotate(90deg);
      }
      /* #selected_gene .selectize-input .item { display: none; } */
      /* Compact grouping rows: kill default margins inside flex cells */
      .grp-row .form-group,
      .grp-row .checkbox { margin: 0; }
      .grp-row .checkbox label { padding-left: 20px; min-height: 0; }
      .grp-row .colourpicker-input-container { width: 40px; }
    "))),
    titlePanel("Spatial Transcriptomics Viewer"),
    sidebarLayout(
      sidebarPanel(
        # Tab 1 controls (Spatial View)
        conditionalPanel(
          condition = "input.main_tabs == 'Spatial View'",
          uiOutput("coord_type_ui"),
          uiOutput("display_mode_ui"),
          selectizeInput("color_var", "Color by:",
            choices = list(
              "Metadata" = setNames(metadata_vars, metadata_vars),
              "Gene Expression" = c("Search gene..." = "GENE_SELECTOR")
            )
          ),
          conditionalPanel(
            condition = "input.color_var == 'GENE_SELECTOR'",
            selectizeInput("selected_gene", "Search gene:",
              choices = NULL, multiple = TRUE,
              options = list(placeholder = "Type gene name...", maxOptions = 50, plugins = list("remove_button"))
            ),
            # uiOutput("mg_gene_list_ui")  # Disabled: remove_button plugin is sufficient
          ),
          uiOutput("point_size_ui"),
          uiOutput("polygon_controls_ui"),
          uiOutput("cont_palette_ui"),
          uiOutput("shuffle_ui"),
          uiOutput("highlight_ui"),
          uiOutput("color_pickers_ui"),
          uiOutput("scale_lock_ui"),
          uiOutput("scale_range_ui"),
          uiOutput("celltype_grouping_ui")
        ),
        # Tab 2 controls (Summary Statistics)
        conditionalPanel(
          condition = "input.main_tabs == 'Summary Statistics'",
          selectizeInput("ref_var", "Reference variable:",
            choices = list(
              "Metadata" = setNames(metadata_vars, metadata_vars),
              "Gene Expression" = c("Search gene..." = "GENE_SELECTOR_REF")
            )
          ),
          conditionalPanel(
            condition = "input.ref_var == 'GENE_SELECTOR_REF'",
            selectizeInput("selected_gene_ref", "Search gene:",
              choices = NULL,
              options = list(placeholder = "Type gene name...", maxOptions = 50)
            )
          ),
          uiOutput("secondary_var_ui"),
          uiOutput("secondary_gene_ui")
        ),
        # Load Dataset section (collapsed at bottom of sidebar)
        tags$hr(),
        tags$details(
          class = "file-picker-section",
          tags$summary(icon("caret-right"), tags$strong("Load Dataset")),
          tags$div(
            style = "padding: 8px 0 4px 0;",
            # Seurat file picker
            tags$details(
              class = "file-picker-section",
              tags$summary(icon("caret-right"), "Seurat RDS File"),
              tags$div(
                style = "padding: 8px 0 4px 8px;",
                shinyFilesButton("file_pick", "Browse...",
                                 title = "Select a Seurat RDS file",
                                 multiple = FALSE),
                tags$br(), tags$br(),
                verbatimTextOutput("current_file_path", placeholder = TRUE)
              )
            ),
            # Polygon file picker
            tags$details(
              class = "file-picker-section",
              tags$summary(icon("caret-right"), "Polygon File (optional)"),
              tags$div(
                style = "padding: 8px 0 4px 8px;",
                shinyFilesButton("poly_file_pick", "Browse...",
                                 title = "Select a polygon CSV file",
                                 multiple = FALSE),
                tags$br(), tags$br(),
                verbatimTextOutput("current_poly_path", placeholder = TRUE),
                actionButton("clear_poly_btn", "Clear", icon = icon("times"))
              )
            ),
            # Shared reload button
            tags$div(
              style = "margin-top: 8px;",
              actionButton("reload_btn", "Reload", icon = icon("refresh"))
            )
          )
        ),
        width = 3
      ),
      mainPanel(
        h4(textOutput("cell_count_text")),
        tabsetPanel(
          id = "main_tabs",
          tabPanel("Spatial View",
            plotlyOutput("spatial_plot", height = "800px")
          ),
          tabPanel("Summary Statistics",
            uiOutput("summary_plots_ui")
          )
        ),
        width = 9
      )
    )
  )

  # Server
  server <- function(input, output, session) {

    # Track zoom state to preserve across variable changes
    zoom_state <- reactiveVal(NULL)
    # Track dragmode (zoom/pan) so re-renders don't reset the toolbar selection
    dragmode_state <- reactiveVal(NULL)

    # Server-side gene search population
    updateSelectizeInput(session, "selected_gene", choices = gene_names, server = TRUE)
    updateSelectizeInput(session, "selected_gene_ref", choices = gene_names, server = TRUE)
    updateSelectizeInput(session, "selected_gene_sec", choices = gene_names, server = TRUE)

    # Multi-gene active tracking
    mg_active <- reactiveVal(character(0))

    # Sync mg_active when genes are added/removed via selectize
    observeEvent(input$selected_gene, {
      current_active <- mg_active()
      selected <- input$selected_gene
      if (is.null(selected)) selected <- character(0)
      new_genes <- setdiff(selected, current_active)
      updated <- intersect(current_active, selected)
      mg_active(c(updated, new_genes))
    }, ignoreNULL = FALSE)

    # Remove gene via × button (JS sets input$mg_remove_gene to "geneName|timestamp")
    observeEvent(input$mg_remove_gene, {
      gene <- sub("\\|.*", "", input$mg_remove_gene)
      current <- input$selected_gene
      updated <- setdiff(current, gene)
      updateSelectizeInput(session, "selected_gene", selected = updated)
      mg_active(setdiff(mg_active(), gene))
    })

    # Toggle gene via checkbox (JS sets input$mg_toggle_gene to "geneName|TRUE/FALSE|timestamp")
    observeEvent(input$mg_toggle_gene, {
      parts <- strsplit(input$mg_toggle_gene, "\\|")[[1]]
      gene <- parts[1]
      checked <- parts[2] == "true"
      active <- mg_active()
      if (checked && !gene %in% active) {
        mg_active(c(active, gene))
      } else if (!checked && gene %in% active) {
        mg_active(setdiff(active, gene))
      }
    })

    # Render gene list with checkboxes and × buttons
    output$mg_gene_list_ui <- renderUI({
      genes <- input$selected_gene
      if (is.null(genes) || length(genes) == 0) return(NULL)
      active <- mg_active()
      tags$div(style = "max-height: 200px; overflow-y: auto; margin-top: 4px;",
        lapply(genes, function(g) {
          is_active <- g %in% active
          safe_g <- gsub("'", "\\\\'", g)
          tags$div(
            style = paste0("display:flex; align-items:center; padding:1px 0;",
                           if (!is_active) " opacity:0.4;" else ""),
            tags$input(type = "checkbox", checked = if (is_active) NA else NULL,
              style = "margin-right: 6px;",
              onchange = sprintf("Shiny.setInputValue('mg_toggle_gene', '%s|' + this.checked + '|' + Date.now())", safe_g)
            ),
            tags$span(style = "flex:1; font-family:monospace; font-size:12px;", g),
            tags$a(href = "#", style = "color:#cc0000; text-decoration:none; font-weight:bold; padding:0 4px;",
              onclick = sprintf("Shiny.setInputValue('mg_remove_gene', '%s|' + Date.now()); return false;", safe_g),
              HTML("&times;")
            )
          )
        }),
        if (length(active) >= 2) {
          tags$p(style = "color:#666; font-size:11px; margin-top:4px;",
            paste0("Any expressed (", length(active), " genes)"))
        }
      )
    })

    # File pickers — browse for .RDS and polygon files
    volumes <- c(Home = path.expand("~"), getVolumes()())
    shinyFileChoose(input, "file_pick", roots = volumes,
                    filetypes = c("RDS", "rds"))
    shinyFileChoose(input, "poly_file_pick", roots = volumes,
                    filetypes = c("csv", "CSV", "gz"))

    selected_path <- reactiveVal(seurat_path)
    selected_poly_path <- reactiveVal(if (!is.null(polygon_path)) polygon_path else "")

    observeEvent(input$file_pick, {
      file_info <- parseFilePaths(volumes, input$file_pick)
      if (nrow(file_info) > 0) {
        selected_path(as.character(file_info$datapath))
      }
    })

    observeEvent(input$poly_file_pick, {
      file_info <- parseFilePaths(volumes, input$poly_file_pick)
      if (nrow(file_info) > 0) {
        selected_poly_path(as.character(file_info$datapath))
      }
    })

    observeEvent(input$clear_poly_btn, {
      selected_poly_path("")
    })

    output$current_file_path <- renderText({
      selected_path()
    })

    output$current_poly_path <- renderText({
      p <- selected_poly_path()
      if (nzchar(p)) p else "(none)"
    })

    observeEvent(input$reload_btn, {
      new_path <- selected_path()
      if (nzchar(new_path) && file.exists(new_path)) {
        reload_file <- file.path(tempdir(), ".spatial_viewer_path")
        writeLines(new_path, reload_file)
        # Write polygon path (empty string means no polygons)
        poly_reload_file <- file.path(tempdir(), ".spatial_viewer_polypath")
        writeLines(selected_poly_path(), poly_reload_file)
        stopApp()
      }
    })

    # Coordinate system selector — hide dimred option if reduction not available
    output$coord_type_ui <- renderUI({
      if (has_reduction) {
        radioButtons("coord_type", "Coordinate System:",
                     choices = c("Spatial (X/Y)" = "spatial",
                                 "Dimensionality Reduction" = "dimred"),
                     selected = "spatial")
      } else {
        radioButtons("coord_type", "Coordinate System:",
                     choices = c("Spatial (X/Y)" = "spatial"),
                     selected = "spatial")
      }
    })

    # Display mode selector (centroids vs polygons) — only when polygons available and spatial mode
    output$display_mode_ui <- renderUI({
      req(input$coord_type == "spatial")
      if (has_polygons) {
        radioButtons("display_mode", "Display mode:",
                     choices = c("Centroids" = "centroids", "Polygons" = "polygons"),
                     selected = "centroids")
      }
    })

    # Shuffle button — hidden in polygon mode
    output$shuffle_ui <- renderUI({
      if (!in_polygon_mode()) {
        actionButton("shuffle", "Shuffle cell order")
      }
    })

    # Point size slider — hidden in polygon mode
    output$point_size_ui <- renderUI({
      if (!in_polygon_mode()) {
        sliderInput("point_size", "Point size:",
                    min = 0.5, max = 7, value = 2.5, step = 0.5)
      }
    })

    # Polygon-specific controls (border color, thickness, render threshold)
    output$polygon_controls_ui <- renderUI({
      if (in_polygon_mode()) {
        tagList(
          tags$hr(),
          tags$strong("Polygon settings:"),
          colourInput("poly_border_color", "Border color:", value = "#000000",
                      showColour = "both"),
          sliderInput("poly_border_width", "Border width:",
                      min = 0, max = 2, value = 0.3, step = 0.1),
          sliderInput("poly_threshold", "Max cells to render:",
                      min = 5000, max = 250000, value = 30000, step = 5000),
          uiOutput("poly_status_text"),
          tags$hr()
        )
      }
    })

    # Helper: is polygon mode active?
    in_polygon_mode <- reactive({
      has_polygons &&
        !is.null(input$coord_type) && input$coord_type == "spatial" &&
        !is.null(input$display_mode) && input$display_mode == "polygons"
    })

    var_name <- reactive(input$color_var)

    showing_gene <- reactive({
      input$color_var == "GENE_SELECTOR" &&
        length(mg_active()) == 1
    })

    showing_multi_gene <- reactive({
      input$color_var == "GENE_SELECTOR" &&
        length(mg_active()) >= 2
    })

    var_is_numeric <- reactive({
      if (showing_gene()) return(TRUE)
      if (showing_multi_gene()) return(FALSE)
      is.numeric(df[[var_name()]])
    })

    # Coordinate system toggle
    plot_coords <- reactive({
      req(input$coord_type)
      if (input$coord_type == "spatial") {
        list(x = x_col, y = y_col)
      } else {
        list(x = "Dim1", y = "Dim2")
      }
    })

    row_order <- reactiveVal(seq_len(nrow(df)))
    observeEvent(input$shuffle, {
      row_order(sample(nrow(df)))
    })

    # Reset zoom and dragmode when switching coordinate systems or display mode
    observeEvent(input$coord_type, {
      zoom_state(NULL)
      dragmode_state(NULL)
    })
    observeEvent(input$display_mode, {
      zoom_state(NULL)
      dragmode_state(NULL)
    })

    # Capture zoom events to preserve zoom state across variable changes
    observeEvent(event_data("plotly_relayout", source = "spatial"), {
      ed <- event_data("plotly_relayout", source = "spatial")
      if (!is.null(ed)) {
        # Preserve dragmode (zoom/pan) across re-renders
        if (!is.null(ed[["dragmode"]])) {
          dragmode_state(ed[["dragmode"]])
        }
        if (any(grepl("range\\[", names(ed)))) {
          zs <- zoom_state()
          # If ranges match current zoom_state, this is a "Reset axes" click
          # (plotly restores the initial layout which has our explicit ranges)
          if (!is.null(zs) &&
              isTRUE(all.equal(ed[["xaxis.range[0]"]],
                               zs[["xaxis.range[0]"]])) &&
              isTRUE(all.equal(ed[["xaxis.range[1]"]],
                               zs[["xaxis.range[1]"]])) &&
              isTRUE(all.equal(ed[["yaxis.range[0]"]],
                               zs[["yaxis.range[0]"]])) &&
              isTRUE(all.equal(ed[["yaxis.range[1]"]],
                               zs[["yaxis.range[1]"]]))) {
            zoom_state(NULL)
          } else {
            zoom_state(ed)
          }
        } else if (isTRUE(ed[["xaxis.autorange"]]) ||
                   isTRUE(ed[["yaxis.autorange"]])) {
          # User double-clicked to reset
          zoom_state(NULL)
        }
      }
    })

    output$cont_palette_ui <- renderUI({
      if (var_is_numeric()) {
        selectInput("cont_palette", "Color palette:", choices = cont_choices)
      }
    })

    output$highlight_ui <- renderUI({
      if (!var_is_numeric() && showing_multi_gene()) {
        # "Any expressed" binary mode
        selectInput("highlight_val", "Highlight value:",
                    choices = c("None", "Expressed", "Not expressed"))
      } else if (!var_is_numeric() && !showing_multi_gene()) {
        lvls <- sort(unique(as.character(df[[var_name()]])))
        selectInput("highlight_val", "Highlight value:",
                    choices = c("None", lvls))
      }
    })

    output$color_pickers_ui <- renderUI({
      if (!var_is_numeric() && !showing_gene() && !showing_multi_gene()) {
        lvls <- sort(unique(as.character(df[[var_name()]])))
        pal <- make_cat_colors(length(lvls))
        tagList(
          tags$hr(),
          tags$strong("Custom colors:"),
          lapply(seq_along(lvls), function(i) {
            colourInput(
              inputId = paste0("col_", i),
              label = lvls[i],
              value = pal[i],
              showColour = "both"
            )
          })
        )
      }
    })

    # Helper: build categorical palette with custom color overrides
    get_cat_palette <- function() {
      if (showing_multi_gene()) {
        return(c("Expressed" = "#E41A1C", "Not expressed" = "#D9D9D9"))
      }
      lvls <- sort(unique(as.character(df[[var_name()]])))
      pal <- make_cat_colors(length(lvls))
      names(pal) <- lvls
      # Override with user-picked colors
      for (i in seq_along(lvls)) {
        val <- input[[paste0("col_", i)]]
        if (!is.null(val)) pal[lvls[i]] <- val
      }
      pal
    }

    # --- Feature: Fixed color scale ---

    global_val_range <- reactive({
      if (!var_is_numeric()) return(c(0, 1))
      vals <- if (showing_gene()) {
        as.numeric(counts_norm[, mg_active()[1]])
      } else {
        df[[var_name()]]
      }
      range(vals, na.rm = TRUE)
    })

    output$scale_lock_ui <- renderUI({
      if (var_is_numeric()) {
        checkboxInput("scale_lock", "Fix color scale", value = FALSE)
      }
    })

    output$scale_range_ui <- renderUI({
      if (!isTRUE(input$scale_lock) || !var_is_numeric()) return(NULL)
      grange <- global_val_range()
      tagList(
        fluidRow(
          column(6, numericInput("scale_min", "Min:", value = round(grange[1], 4), step = 0.001)),
          column(6, numericInput("scale_max", "Max:", value = round(grange[2], 4), step = 0.001))
        ),
        actionButton("scale_reset", "Reset to data range",
                     style = "font-size:11px; padding:2px 6px; margin-bottom:6px;")
      )
    })

    observeEvent(input$scale_reset, {
      grange <- global_val_range()
      updateNumericInput(session, "scale_min", value = round(grange[1], 4))
      updateNumericInput(session, "scale_max", value = round(grange[2], 4))
    })

    observeEvent(input$color_var, {
      updateCheckboxInput(session, "scale_lock", value = FALSE)
    }, ignoreInit = TRUE)

    # --- Feature: Cell type grouping ---

    cat_grouping_cols <- Filter(function(v) {
      x <- df[[v]]
      (is.character(x) || is.factor(x)) && length(unique(x)) <= 50
    }, metadata_vars)

    output$celltype_grouping_ui <- renderUI({
      tagList(
        tags$hr(),
        tags$strong("Cell type grouping"),
        checkboxInput("group_mode", "Enable", value = FALSE),
        conditionalPanel(
          condition = "input.group_mode == true",
          if (length(cat_grouping_cols) > 0) {
            selectInput("group_col", "Cell type column:", choices = cat_grouping_cols,
                        selected = if (!is.null(celltype_col) && celltype_col %in% cat_grouping_cols)
                                     celltype_col else cat_grouping_cols[1])
          } else {
            tags$p(style = "color:#888; font-size:11px;",
                   "No categorical columns with \u226450 levels found.")
          },
          fluidRow(
            style = "margin-bottom:4px;",
            column(6, actionButton("group_set_all_fixed", "All \u2192 Fixed",
                                   style = "font-size:10px; padding:2px 5px; width:100%;")),
            column(6, actionButton("group_set_all_normal", "All \u2192 Normal",
                                   style = "font-size:10px; padding:2px 5px; width:100%;"))
          ),
          tags$div(
            style = "display:flex; align-items:center; gap:6px; margin-bottom:4px;",
            tags$span(style = "font-size:11px; color:#555;", "Background color:"),
            colourInput("group_other_color", NULL, value = "#D9D9D9", showColour = "both")
          ),
          tags$p(style = "color:#666; font-size:11px; margin:2px 0 4px 0;",
                 "Normal = use current color variable. Fixed = render as a solid color."),
          uiOutput("group_types_ui")
        )
      )
    })

    # Stable, sorted list of cell-type levels for the selected grouping column.
    # Both the UI and cell_group_map() index into this list, so per-type inputs
    # use integer indices (grp_normal_1, grp_color_1, ...) — no name sanitization,
    # no collisions.
    group_levels <- reactive({
      if (!isTRUE(input$group_mode)) return(character(0))
      col <- input$group_col
      if (is.null(col) || !nzchar(col) || !col %in% colnames(df)) return(character(0))
      sort(unique(as.character(df[[col]])))
    })

    # Track the previous background color so propagation only overwrites per-type
    # swatches that were still at the old default — user-set overrides (e.g. red
    # for tumor) survive a background color change.
    prev_bg <- reactiveVal("#D9D9D9")

    output$group_types_ui <- renderUI({
      req(isTRUE(input$group_mode))
      lvls <- group_levels()
      if (length(lvls) == 0) return(NULL)
      default_bg <- isolate({
        v <- input$group_other_color
        if (!is.null(v) && nzchar(v)) v else "#D9D9D9"
      })
      tags$div(
        style = "max-height:300px; overflow-y:auto; overflow-x:auto; margin-top:4px; border:1px solid #ddd; border-radius:3px; padding:6px;",
        tags$div(
          style = "min-width:260px; display:flex; align-items:center; gap:8px; font-size:10px; color:#888; padding:0 0 4px 0; border-bottom:1px solid #eee; margin-bottom:4px;",
          tags$span(style = "flex:1;", "Cell type"),
          tags$span(style = "width:50px; text-align:center;", "Normal"),
          tags$span(style = "width:40px; text-align:center;", "Color")
        ),
        lapply(seq_along(lvls), function(i) {
          v <- lvls[i]
          tags$div(
            class = "grp-row",
            style = "min-width:260px; display:flex; align-items:center; gap:8px; padding:2px 0;",
            tags$span(
              style = "flex:1; font-size:11px; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; min-width:0;",
              title = v, v
            ),
            tags$div(
              style = "width:50px; display:flex; justify-content:center;",
              checkboxInput(paste0("grp_normal_", i), NULL, value = FALSE)
            ),
            tags$div(
              style = "width:90px;",
              colourInput(paste0("grp_color_", i), NULL, value = default_bg,
                          showColour = "both")
            )
          )
        })
      )
    })

    cell_group_map <- reactive({
      if (!isTRUE(input$group_mode)) return(NULL)
      lvls <- group_levels()
      if (length(lvls) == 0) return(NULL)
      col <- input$group_col
      vals <- as.character(df[[col]])
      bg_color <- if (!is.null(input$group_other_color) && nzchar(input$group_other_color))
                    input$group_other_color else "#D9D9D9"
      result <- rep("primary", nrow(df))
      any_configured <- FALSE
      for (i in seq_along(lvls)) {
        v <- lvls[i]
        idx <- which(vals == v)
        normal_flag <- input[[paste0("grp_normal_", i)]]
        color_val   <- input[[paste0("grp_color_", i)]]
        if (!is.null(normal_flag) || !is.null(color_val)) any_configured <- TRUE
        if (!isTRUE(normal_flag)) {
          result[idx] <- if (!is.null(color_val) && nzchar(color_val)) color_val else bg_color
        }
      }
      # UI not yet wired — don't apply grouping on first paint
      if (!any_configured) return(NULL)
      result
    })

    observeEvent(input$group_set_all_fixed, {
      lvls <- group_levels()
      if (length(lvls) == 0) return()
      for (i in seq_along(lvls)) {
        updateCheckboxInput(session, paste0("grp_normal_", i), value = FALSE)
      }
    })

    observeEvent(input$group_set_all_normal, {
      lvls <- group_levels()
      if (length(lvls) == 0) return()
      for (i in seq_along(lvls)) {
        updateCheckboxInput(session, paste0("grp_normal_", i), value = TRUE)
      }
    })

    # Propagate shared background color only to swatches still at the old
    # default — per-type overrides (e.g. tumor = red) are preserved.
    observeEvent(input$group_other_color, {
      if (!isTRUE(input$group_mode)) return()
      new_color <- input$group_other_color
      if (is.null(new_color) || !nzchar(new_color)) return()
      lvls <- group_levels()
      old_color <- prev_bg()
      for (i in seq_along(lvls)) {
        current <- isolate(input[[paste0("grp_color_", i)]])
        if (is.null(current) || !nzchar(current) ||
            tolower(current) == tolower(old_color)) {
          updateColourInput(session, paste0("grp_color_", i), value = new_color)
        }
      }
      prev_bg(new_color)
    }, ignoreInit = TRUE)

    # Generalized palette builder for any categorical variable (used by Tab 2)
    build_palette_for <- function(var_name_str) {
      lvls <- sort(unique(as.character(df[[var_name_str]])))
      pal <- make_cat_colors(length(lvls))
      names(pal) <- lvls
      # If this variable matches Tab 1's color_var, apply user overrides
      if (!is.null(input$color_var) && var_name_str == input$color_var && !is.numeric(df[[var_name_str]])) {
        for (i in seq_along(lvls)) {
          val <- input[[paste0("col_", i)]]
          if (!is.null(val)) pal[lvls[i]] <- val
        }
      }
      pal
    }

    # ---- Tab 2: Summary Statistics ----

    # Cell count (displayed on both tabs)
    output$cell_count_text <- renderText({
      paste0("Total cells: ", format(nrow(df), big.mark = ","))
    })

    # Tab 2 reactive helpers — gene detection
    ref_showing_gene <- reactive({
      input$ref_var == "GENE_SELECTOR_REF" &&
        !is.null(input$selected_gene_ref) &&
        input$selected_gene_ref != ""
    })

    sec_showing_gene <- reactive({
      !is.null(input$sec_var) &&
        input$sec_var == "GENE_SELECTOR_SEC" &&
        !is.null(input$selected_gene_sec) &&
        input$selected_gene_sec != ""
    })

    ref_is_numeric <- reactive({
      if (ref_showing_gene()) return(TRUE)
      req(input$ref_var)
      input$ref_var != "GENE_SELECTOR_REF" && is.numeric(df[[input$ref_var]])
    })

    sec_is_numeric <- reactive({
      if (sec_showing_gene()) return(TRUE)
      req(input$sec_var)
      if (input$sec_var == "None") return(NULL)
      input$sec_var != "GENE_SELECTOR_SEC" && is.numeric(df[[input$sec_var]])
    })

    # Tab 2 value and label helpers
    ref_values <- reactive({
      if (ref_showing_gene()) {
        as.numeric(counts_norm[, input$selected_gene_ref])
      } else {
        req(input$ref_var != "GENE_SELECTOR_REF")
        df[[input$ref_var]]
      }
    })

    ref_label <- reactive({
      if (ref_showing_gene()) {
        paste0(input$selected_gene_ref, " (log-normalized)")
      } else {
        input$ref_var
      }
    })

    sec_values <- reactive({
      if (sec_showing_gene()) {
        as.numeric(counts_norm[, input$selected_gene_sec])
      } else if (!is.null(input$sec_var) && input$sec_var != "None" && input$sec_var != "GENE_SELECTOR_SEC") {
        df[[input$sec_var]]
      } else {
        NULL
      }
    })

    sec_label <- reactive({
      if (sec_showing_gene()) {
        paste0(input$selected_gene_sec, " (log-normalized)")
      } else {
        input$sec_var
      }
    })

    plot_combo <- reactive({
      req(input$ref_var)
      # Wait for gene selection if gene selector is chosen
      if (input$ref_var == "GENE_SELECTOR_REF") req(ref_showing_gene())
      sec <- if (is.null(input$sec_var) || input$sec_var == "None") {
        "none"
      } else {
        if (input$sec_var == "GENE_SELECTOR_SEC") req(sec_showing_gene())
        if (sec_is_numeric()) "continuous" else "categorical"
      }
      ref <- if (ref_is_numeric()) "continuous" else "categorical"
      paste(ref, sec, sep = "_")
    })

    # Secondary variable dropdown
    output$secondary_var_ui <- renderUI({
      req(input$ref_var)
      choices <- list(
        "None" = "None",
        "Metadata" = setNames(metadata_vars, metadata_vars),
        "Gene Expression" = c("Search gene..." = "GENE_SELECTOR_SEC")
      )
      selectizeInput("sec_var", "Secondary variable:", choices = choices, selected = "None")
    })

    # Secondary gene search (conditional)
    output$secondary_gene_ui <- renderUI({
      req(input$sec_var == "GENE_SELECTOR_SEC")
      selectizeInput("selected_gene_sec", "Search gene:",
        choices = NULL,
        options = list(placeholder = "Type gene name...", maxOptions = 50)
      )
    })

    # Re-populate secondary gene search when it appears
    observeEvent(input$sec_var, {
      if (!is.null(input$sec_var) && input$sec_var == "GENE_SELECTOR_SEC") {
        updateSelectizeInput(session, "selected_gene_sec", choices = gene_names, server = TRUE)
      }
    })

    # Dynamic plot container
    output$summary_plots_ui <- renderUI({
      combo <- plot_combo()

      switch(combo,
        "categorical_none" = tagList(
          plotlyOutput("summary_bar", height = "500px")
        ),
        "categorical_continuous" = tagList(
          plotlyOutput("summary_violin_cat_cont", height = "500px"),
          plotlyOutput("summary_box_cat_cont", height = "500px")
        ),
        "categorical_categorical" = tagList(
          plotlyOutput("summary_prop_bar", height = "500px"),
          plotlyOutput("summary_abs_bar", height = "500px")
        ),
        "continuous_none" = fluidRow(
          column(6, plotlyOutput("summary_violin_cont", height = "400px")),
          column(6, plotlyOutput("summary_box_cont", height = "400px"))
        ),
        "continuous_continuous" = tagList(
          tags$div(
            style = "padding: 40px; text-align: center; color: #888;",
            tags$h4("No summary plots available for two continuous variables."),
            tags$p("Select a categorical variable as either Reference or Secondary to see distribution plots.")
          )
        ),
        "continuous_categorical" = tagList(
          plotlyOutput("summary_violin_cont_cat", height = "500px"),
          plotlyOutput("summary_box_cont_cat", height = "500px")
        )
      )
    })

    # ---- Tab 2 Plot Renderers ----

    # Categorical ref, no secondary: Bar chart of counts
    output$summary_bar <- renderPlotly({
      req(input$ref_var, !ref_is_numeric())

      rv <- ref_values()
      rl <- ref_label()
      counts <- as.data.frame(table(rv), stringsAsFactors = FALSE)
      colnames(counts) <- c("Category", "Count")
      counts <- counts[order(-counts$Count), ]
      counts$Category <- factor(counts$Category, levels = counts$Category)
      total <- sum(counts$Count)
      counts$Label <- paste0(round(100 * counts$Count / total, 1), "%")

      pal <- build_palette_for(input$ref_var)

      p <- ggplot(counts, aes(x = Category, y = Count, fill = Category)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = Label), vjust = -0.3, size = 3) +
        scale_fill_manual(values = pal) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none") +
        labs(title = paste("Cell counts by", rl), x = rl, y = "Count")

      ggplotly(p, tooltip = c("x", "y")) %>% layout(xaxis = list(tickangle = -45))
    })

    # Categorical ref + Continuous secondary: Violin
    output$summary_violin_cat_cont <- renderPlotly({
      req(input$ref_var, input$sec_var, !ref_is_numeric(), sec_is_numeric())

      rl <- ref_label()
      sl <- sec_label()
      plot_data <- data.frame(ref = ref_values(), sec = sec_values())
      pal <- build_palette_for(input$ref_var)

      p <- ggplot(plot_data, aes(x = ref, y = sec, fill = ref)) +
        geom_violin(trim = FALSE, width = 0.9) +
        stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "black") +
        scale_fill_manual(values = pal) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none") +
        labs(title = paste(sl, "distribution by", rl), x = rl, y = sl)

      ggplotly(p)
    })

    # Categorical ref + Continuous secondary: Box
    output$summary_box_cat_cont <- renderPlotly({
      req(input$ref_var, input$sec_var, !ref_is_numeric(), sec_is_numeric())

      rl <- ref_label()
      sl <- sec_label()
      plot_data <- data.frame(ref = ref_values(), sec = sec_values())
      pal <- build_palette_for(input$ref_var)

      p <- ggplot(plot_data, aes(x = ref, y = sec, fill = ref)) +
        geom_boxplot(fatten = NULL, outlier.size = 0.5, width = 0.8) +
        stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.8) +
        scale_fill_manual(values = pal) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none") +
        labs(title = paste(sl, "distribution by", rl), x = rl, y = sl)

      ggplotly(p)
    })

    # Categorical ref + Categorical secondary: 100% stacked bar (proportional)
    output$summary_prop_bar <- renderPlotly({
      req(input$ref_var, input$sec_var, !ref_is_numeric(), !is.null(sec_is_numeric()), !sec_is_numeric())

      rl <- ref_label()
      sl <- sec_label()
      rv <- ref_values()
      sv <- sec_values()

      ct <- as.data.frame(table(rv, sv), stringsAsFactors = FALSE)
      colnames(ct) <- c("Reference", "Secondary", "Count")

      ct <- do.call(rbind, lapply(split(ct, ct$Reference), function(x) {
        x$Pct <- 100 * x$Count / sum(x$Count)
        x$Label <- ifelse(x$Pct >= 3, paste0(round(x$Pct, 1), "%"), "")
        x
      }))

      sec_pal <- build_palette_for(input$sec_var)

      p <- ggplot(ct, aes(x = Reference, y = Pct, fill = Secondary)) +
        geom_bar(stat = "identity", position = "fill") +
        geom_text(aes(label = Label), position = position_fill(vjust = 0.5), size = 2.5) +
        scale_fill_manual(values = sec_pal) +
        scale_y_continuous(labels = scales::percent_format()) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = paste("Composition of", sl, "by", rl),
             x = rl, y = "Percentage", fill = sl)

      ggplotly(p)
    })

    # Categorical ref + Categorical secondary: Absolute stacked bar (counts)
    output$summary_abs_bar <- renderPlotly({
      req(input$ref_var, input$sec_var, !ref_is_numeric(), !is.null(sec_is_numeric()), !sec_is_numeric())

      rl <- ref_label()
      sl <- sec_label()
      rv <- ref_values()
      sv <- sec_values()

      ct <- as.data.frame(table(rv, sv), stringsAsFactors = FALSE)
      colnames(ct) <- c("Reference", "Secondary", "Count")
      sec_pal <- build_palette_for(input$sec_var)

      p <- ggplot(ct, aes(x = Reference, y = Count, fill = Secondary)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = sec_pal) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = paste("Cell counts of", sl, "by", rl),
             x = rl, y = "Count", fill = sl)

      ggplotly(p)
    })

    # Continuous ref, no secondary: Violin
    output$summary_violin_cont <- renderPlotly({
      req(input$ref_var, ref_is_numeric())

      rl <- ref_label()
      plot_data <- data.frame(ref = ref_values())

      p <- ggplot(plot_data, aes(x = "", y = ref)) +
        geom_violin(fill = "#4ECDC4", trim = FALSE) +
        stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "black") +
        theme_minimal() +
        labs(title = paste("Distribution of", rl), x = "", y = rl)

      ggplotly(p)
    })

    # Continuous ref, no secondary: Box
    output$summary_box_cont <- renderPlotly({
      req(input$ref_var, ref_is_numeric())

      rl <- ref_label()
      plot_data <- data.frame(ref = ref_values())

      p <- ggplot(plot_data, aes(x = "", y = ref)) +
        geom_boxplot(fill = "#4ECDC4", fatten = NULL, outlier.size = 0.5) +
        stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.8) +
        theme_minimal() +
        labs(title = paste("Distribution of", rl), x = "", y = rl)

      ggplotly(p)
    })

    # Continuous ref + Categorical secondary: Violin
    output$summary_violin_cont_cat <- renderPlotly({
      req(input$ref_var, input$sec_var, ref_is_numeric(), !is.null(sec_is_numeric()), !sec_is_numeric())

      rl <- ref_label()
      sl <- sec_label()
      plot_data <- data.frame(ref = ref_values(), sec = sec_values())
      sec_pal <- build_palette_for(input$sec_var)

      p <- ggplot(plot_data, aes(x = sec, y = ref, fill = sec)) +
        geom_violin(trim = FALSE, width = 0.9) +
        stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "black") +
        scale_fill_manual(values = sec_pal) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none") +
        labs(title = paste(rl, "distribution by", sl), x = sl, y = rl)

      ggplotly(p)
    })

    # Continuous ref + Categorical secondary: Box
    output$summary_box_cont_cat <- renderPlotly({
      req(input$ref_var, input$sec_var, ref_is_numeric(), !is.null(sec_is_numeric()), !sec_is_numeric())

      rl <- ref_label()
      sl <- sec_label()
      plot_data <- data.frame(ref = ref_values(), sec = sec_values())
      sec_pal <- build_palette_for(input$sec_var)

      p <- ggplot(plot_data, aes(x = sec, y = ref, fill = sec)) +
        geom_boxplot(fatten = NULL, outlier.size = 0.5, width = 0.8) +
        stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.8) +
        scale_fill_manual(values = sec_pal) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none") +
        labs(title = paste(rl, "distribution by", sl), x = sl, y = rl)

      ggplotly(p)
    })

    # Helper: build NA-separated polygon vertex vectors for a set of cells
    build_polygon_traces <- function(cell_ids) {
      if (length(cell_ids) == 0) return(list(x = numeric(0), y = numeric(0)))
      sub <- poly_dt[poly_dt[["cell"]] %in% cell_ids, ]
      result <- sub[, list(
        x = c(x_global_px, x_global_px[1L], NA_real_),
        y = c(y_global_px, y_global_px[1L], NA_real_)
      ), by = "cell"]
      list(x = result$x, y = result$y)
    }

    # Main plot
    output$spatial_plot <- renderPlotly({
      coords <- plot_coords()
      vname <- var_name()
      ord <- row_order()
      sz <- isolate(input$point_size)
      polygon_mode <- in_polygon_mode()

      plot_df <- df[ord, , drop = FALSE]

      # Get plot values: gene expression, multi-gene, or metadata
      if (showing_gene()) {
        plot_vals <- as.numeric(counts_norm[ord, mg_active()[1]])
        color_title <- paste0(mg_active()[1], "\n(log-normalized)")
      } else if (showing_multi_gene()) {
        active_genes <- mg_active()
        sub_mat <- counts_norm[ord, active_genes, drop = FALSE]
        expressed <- Matrix::rowSums(sub_mat) > 0
        plot_vals <- factor(ifelse(expressed, "Expressed", "Not expressed"),
                            levels = c("Expressed", "Not expressed"))
        color_title <- paste0("Any expressed\n(", length(active_genes), " genes)")
      } else {
        req(vname != "GENE_SELECTOR")
        plot_vals <- plot_df[[vname]]
        color_title <- vname
      }

      if (polygon_mode) {
        # ---- POLYGON RENDERING ----
        cell_ids_ordered <- df[[cell_id_col]][ord]
        border_col <- if (!is.null(input$poly_border_color)) input$poly_border_color else "#000000"
        border_w <- if (!is.null(input$poly_border_width)) input$poly_border_width else 0.3
        threshold <- if (!is.null(input$poly_threshold)) input$poly_threshold else 5000

        # Determine which cells are in the current zoom window
        zs <- zoom_state()
        if (!is.null(zs) &&
            !is.null(zs[["xaxis.range[0]"]]) && !is.null(zs[["xaxis.range[1]"]]) &&
            !is.null(zs[["yaxis.range[0]"]]) && !is.null(zs[["yaxis.range[1]"]])) {
          xmin <- zs[["xaxis.range[0]"]]
          xmax <- zs[["xaxis.range[1]"]]
          ymin <- zs[["yaxis.range[0]"]]
          ymax <- zs[["yaxis.range[1]"]]
          visible_cells <- poly_centroids$cell[
            poly_centroids$cx >= xmin & poly_centroids$cx <= xmax &
            poly_centroids$cy >= ymin & poly_centroids$cy <= ymax]
        } else {
          visible_cells <- poly_centroids$cell
        }
        n_visible <- length(visible_cells)

        # Update status text
        output$poly_status_text <- renderUI({
          if (n_visible > threshold) {
            tags$p(style = "color: #cc6600; font-weight: bold;",
              paste0("Cells in view: ", format(n_visible, big.mark = ","),
                     " (above threshold of ", format(threshold, big.mark = ","),
                     "). Zoom in to render polygons."))
          } else {
            tags$p(style = "color: #228B22;",
              paste0("Rendering ", format(n_visible, big.mark = ","), " cell polygons."))
          }
        })

        if (n_visible > threshold) {
          # Too many cells — show centroid scatter as placeholder with message
          # Use polygon centroids for positioning
          vis_idx <- match(visible_cells, df[[cell_id_col]])
          vis_idx <- vis_idx[!is.na(vis_idx)]
          vis_centroids <- poly_centroids[poly_centroids[["cell"]] %in% visible_cells, ]

          p <- plot_ly(source = "spatial") %>%
            add_trace(
              x = vis_centroids$cx, y = vis_centroids$cy,
              type = "scattergl", mode = "markers",
              marker = list(size = 2, color = "#999999", opacity = 0.4),
              showlegend = FALSE, hoverinfo = "none"
            ) %>%
            add_annotations(
              text = paste0("Zoom in to render polygons\n(",
                            format(n_visible, big.mark = ","),
                            " cells > threshold of ",
                            format(threshold, big.mark = ","), ")"),
              xref = "paper", yref = "paper", x = 0.5, y = 0.5,
              showarrow = FALSE,
              font = list(size = 16, color = "#cc6600")
            )
        } else {
          # Filter ordered data to only visible cells
          vis_mask <- cell_ids_ordered %in% visible_cells
          vis_cell_ids <- cell_ids_ordered[vis_mask]
          vis_plot_vals <- plot_vals[vis_mask]

          # Cell type grouping: split visible cells into primary vs fixed-color
          group_map_r <- cell_group_map()
          if (!is.null(group_map_r)) {
            vis_group <- group_map_r[ord][vis_mask]
            vis_primary_mask <- vis_group == "primary"
          } else {
            vis_group <- NULL
            vis_primary_mask <- rep(TRUE, length(vis_cell_ids))
          }
          primary_cell_ids    <- vis_cell_ids[vis_primary_mask]
          primary_plot_vals_p <- vis_plot_vals[vis_primary_mask]

          p <- plot_ly(source = "spatial")

          # Render fixed-color (Other/Solid) polygons first (background layer)
          if (!is.null(vis_group) && any(!vis_primary_mask)) {
            for (hex_col in sort(unique(vis_group[!vis_primary_mask]))) {
              cids <- vis_cell_ids[vis_group == hex_col]
              trace_data <- build_polygon_traces(cids)
              if (length(trace_data$x) > 0) {
                p <- p %>% add_trace(
                  x = trace_data$x, y = trace_data$y,
                  type = "scatter", mode = "lines",
                  fill = "toself",
                  fillcolor = adjustcolor(hex_col, alpha.f = 0.85),
                  line = list(color = border_col, width = border_w),
                  showlegend = FALSE, hoverinfo = "none"
                )
              }
            }
          }

          if (length(primary_cell_ids) > 0) {
            if (var_is_numeric()) {
              # Continuous variable: bin into color groups
              pal_name <- if (!is.null(input$cont_palette)) input$cont_palette else "Inferno"
              n_bins <- 100

              # Scale range: fixed or auto-from-visible
              use_fixed_scale <- isTRUE(input$scale_lock) &&
                !is.null(input$scale_min) && !is.null(input$scale_max) &&
                !is.na(input$scale_min) && !is.na(input$scale_max) &&
                input$scale_min < input$scale_max
              if (use_fixed_scale) {
                val_range <- c(input$scale_min, input$scale_max)
                vals_for_bins <- pmin(pmax(primary_plot_vals_p, val_range[1]), val_range[2])
              } else {
                val_range <- range(primary_plot_vals_p, na.rm = TRUE)
                vals_for_bins <- primary_plot_vals_p
              }

              if (val_range[1] == val_range[2]) {
                bins <- rep(1L, length(vals_for_bins))
              } else {
                bins <- as.integer(cut(vals_for_bins, breaks = n_bins, include.lowest = TRUE))
              }
              bins[is.na(bins)] <- 1L

              pal_func <- get_cont_pal_func(pal_name)
              bin_colors <- pal_func(n_bins)

              for (b in sort(unique(bins))) {
                idx <- which(bins == b)
                cids <- primary_cell_ids[idx]
                trace_data <- build_polygon_traces(cids)
                if (length(trace_data$x) > 0) {
                  p <- p %>% add_trace(
                    x = trace_data$x, y = trace_data$y,
                    type = "scatter", mode = "lines",
                    fill = "toself",
                    fillcolor = adjustcolor(bin_colors[b], alpha.f = 0.85),
                    line = list(color = border_col, width = border_w),
                    showlegend = FALSE, hoverinfo = "none"
                  )
                }
              }
              # Invisible trace for colorbar — use explicit colorscale so the bar
              # always matches the polygon fills regardless of Plotly's native
              # name support (e.g. "Inferno" string falls back to gray-to-red).
              n_bar <- 20
              bar_hex <- pal_func(n_bar)
              colorscale_bar <- lapply(seq_len(n_bar),
                function(i) list((i - 1) / (n_bar - 1), bar_hex[i]))
              p <- p %>% add_trace(
                x = c(0, 0), y = c(0, 0),
                type = "scattergl", mode = "markers",
                marker = list(
                  size = 0.001, opacity = 0,
                  color = val_range,
                  colorscale = colorscale_bar,
                  colorbar = list(title = color_title),
                  showscale = TRUE
                ),
                showlegend = FALSE, hoverinfo = "none"
              )
            } else {
              # Categorical variable
              pal <- get_cat_palette()
              present_lvls <- intersect(names(pal), unique(as.character(primary_plot_vals_p)))
              if (length(present_lvls) > 0) pal <- pal[present_lvls]
              fvals <- factor(primary_plot_vals_p, levels = present_lvls)
              highlight <- if (!is.null(input$highlight_val)) input$highlight_val else "None"

              if (highlight != "None" && highlight %in% levels(fvals)) {
                for (lvl in levels(fvals)) {
                  idx <- which(fvals == lvl)
                  cids <- primary_cell_ids[idx]
                  trace_data <- build_polygon_traces(cids)
                  if (length(trace_data$x) > 0) {
                    is_highlight <- (lvl == highlight)
                    alpha <- if (is_highlight) 0.85 else 0.15
                    fill_col <- adjustcolor(pal[lvl], alpha.f = alpha)
                    line_col <- if (is_highlight) border_col else adjustcolor(border_col, alpha.f = 0.15)
                    p <- p %>% add_trace(
                      x = trace_data$x, y = trace_data$y,
                      type = "scatter", mode = "lines",
                      fill = "toself", fillcolor = fill_col,
                      line = list(color = line_col, width = border_w),
                      name = if (is_highlight) lvl else "Other",
                      showlegend = is_highlight, hoverinfo = "none"
                    )
                  }
                }
              } else {
                for (lvl in levels(fvals)) {
                  idx <- which(fvals == lvl)
                  cids <- primary_cell_ids[idx]
                  trace_data <- build_polygon_traces(cids)
                  if (length(trace_data$x) > 0) {
                    p <- p %>% add_trace(
                      x = trace_data$x, y = trace_data$y,
                      type = "scatter", mode = "lines",
                      fill = "toself",
                      fillcolor = adjustcolor(pal[lvl], alpha.f = 0.85),
                      line = list(color = border_col, width = border_w),
                      name = lvl, hoverinfo = "none"
                    )
                  }
                }
              }
            }
          }
        }

        # Polygon mode axis config (px coordinates, 1:1 aspect)
        xaxis_config <- list(title = "x_global_px", zeroline = FALSE)
        yaxis_config <- list(title = "y_global_px", zeroline = FALSE,
                             scaleanchor = "x", scaleratio = 1)

      } else {
        # ---- CENTROID RENDERING ----
        group_map_r <- cell_group_map()

        if (!is.null(group_map_r)) {
          plot_group <- group_map_r[ord]
          primary_mask_c <- plot_group == "primary"
        } else {
          plot_group <- NULL
          primary_mask_c <- rep(TRUE, nrow(plot_df))
        }

        p <- plot_ly(source = "spatial")

        # Render fixed-color (Other/Solid) cells first (background layer)
        if (!is.null(plot_group) && any(!primary_mask_c)) {
          for (hex_col in sort(unique(plot_group[!primary_mask_c]))) {
            idx <- which(plot_group == hex_col)
            p <- p %>% add_trace(
              x = plot_df[[coords$x]][idx],
              y = plot_df[[coords$y]][idx],
              type = "scattergl", mode = "markers",
              marker = list(size = sz, color = hex_col),
              showlegend = FALSE, hoverinfo = "none"
            )
          }
        }

        primary_plot_df    <- plot_df[primary_mask_c, , drop = FALSE]
        primary_plot_vals_c <- plot_vals[primary_mask_c]

        if (nrow(primary_plot_df) > 0) {
          if (var_is_numeric()) {
            pal_name <- if (!is.null(input$cont_palette)) input$cont_palette else "Inferno"
            n_bar <- 20
            cent_bar_hex <- get_cont_pal_func(pal_name)(n_bar)
            cent_colorscale <- lapply(seq_len(n_bar),
              function(i) list((i - 1) / (n_bar - 1), cent_bar_hex[i]))
            marker_list <- list(
              size = sz,
              color = primary_plot_vals_c,
              colorscale = cent_colorscale,
              colorbar = list(title = color_title)
            )
            if (isTRUE(input$scale_lock) &&
                !is.null(input$scale_min) && !is.null(input$scale_max) &&
                !is.na(input$scale_min) && !is.na(input$scale_max)) {
              marker_list$cmin <- input$scale_min
              marker_list$cmax <- input$scale_max
            }
            p <- p %>% add_trace(
              x = primary_plot_df[[coords$x]],
              y = primary_plot_df[[coords$y]],
              type = "scattergl", mode = "markers",
              marker = marker_list,
              hoverinfo = "none"
            )
          } else {
            pal <- get_cat_palette()
            present_lvls <- intersect(names(pal), unique(as.character(primary_plot_vals_c)))
            if (length(present_lvls) > 0) pal <- pal[present_lvls]
            fvals <- factor(primary_plot_vals_c, levels = present_lvls)
            highlight <- if (!is.null(input$highlight_val)) input$highlight_val else "None"

            if (highlight != "None" && highlight %in% levels(fvals)) {
              bg_idx <- which(fvals != highlight)
              fg_idx <- which(fvals == highlight)
              p <- p %>%
                add_trace(
                  x = primary_plot_df[[coords$x]][bg_idx],
                  y = primary_plot_df[[coords$y]][bg_idx],
                  type = "scattergl", mode = "markers",
                  marker = list(size = sz,
                                color = as.character(pal[as.character(fvals[bg_idx])]),
                                opacity = 0.2),
                  name = "Other", hoverinfo = "none"
                ) %>%
                add_trace(
                  x = primary_plot_df[[coords$x]][fg_idx],
                  y = primary_plot_df[[coords$y]][fg_idx],
                  type = "scattergl", mode = "markers",
                  marker = list(size = sz, color = pal[highlight], opacity = 1.0),
                  name = highlight, hoverinfo = "none"
                )
            } else {
              point_colors <- as.character(pal[as.character(fvals)])
              p <- p %>% add_trace(
                x = primary_plot_df[[coords$x]],
                y = primary_plot_df[[coords$y]],
                type = "scattergl", mode = "markers",
                marker = list(size = sz, color = point_colors),
                hoverinfo = "none", showlegend = FALSE
              )
              for (lvl in levels(fvals)) {
                p <- p %>% add_trace(
                  x = primary_plot_df[[coords$x]][match(lvl, as.character(fvals))],
                  y = primary_plot_df[[coords$y]][match(lvl, as.character(fvals))],
                  type = "scattergl", mode = "markers",
                  marker = list(size = sz, color = pal[lvl]),
                  name = lvl, hoverinfo = "none"
                )
              }
            }
          }
        }

        # Centroid mode axis config
        yaxis_config <- list(title = coords$y, zeroline = FALSE)
        if (input$coord_type == "spatial") {
          yaxis_config$scaleanchor <- "x"
          yaxis_config$scaleratio <- 1
        }
        xaxis_config <- list(title = coords$x, zeroline = FALSE)
      }

      # Re-apply stored zoom state if available
      zs <- zoom_state()
      if (!is.null(zs)) {
        if (!is.null(zs[["xaxis.range[0]"]]) && !is.null(zs[["xaxis.range[1]"]])) {
          xaxis_config$range <- c(zs[["xaxis.range[0]"]], zs[["xaxis.range[1]"]])
        }
        if (!is.null(zs[["yaxis.range[0]"]]) && !is.null(zs[["yaxis.range[1]"]])) {
          yaxis_config$range <- c(zs[["yaxis.range[0]"]], zs[["yaxis.range[1]"]])
        }
      }

      layout_args <- list(
        xaxis = xaxis_config,
        yaxis = yaxis_config,
        legend = list(itemsizing = "constant")
      )
      dm <- dragmode_state()
      if (!is.null(dm)) layout_args$dragmode <- dm

      do.call(layout, c(list(p), layout_args)) %>%
        event_register("plotly_relayout")
    })

    # Update marker size via plotly proxy to preserve zoom state
    observeEvent(input$point_size, {
      if (!in_polygon_mode()) {
        plotlyProxy("spatial_plot", session) %>%
          plotlyProxyInvoke("restyle", list(marker.size = input$point_size))
      }
    }, ignoreInit = TRUE)
  }

  runApp(shinyApp(ui, server))

  if (!file.exists(file.path(tempdir(), ".spatial_viewer_path"))) break
  seurat_path <- NULL
  } # end repeat
  invisible(NULL)
}
