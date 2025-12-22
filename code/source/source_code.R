suppressPackageStartupMessages({
  # Core and data structures
  library(magrittr)
  library(job)
  library(tidyverse)
  # library(data.table)
  # library(matrixStats)
  
  # Visualization
  library(ggplot2)
  library(ggfortify)
  library(ggpubr)
  library(ggsci)
  library(RColorBrewer)
  library(viridis)
  library(pals)
  # library(scales)
  library(pheatmap)
  library(patchwork)
  # 
  # Dimensionality reduction and clustering
  library(Rtsne)
  library(umap)
  # 
  # Statistics and bioinformatics
  library(broom)
  # library(limma)
  # library(sva)
  # library(preprocessCore)
  # library(BiocParallel)
  # library(parallel)
  # library(qcrlscR)
  # library(imputeLCMD)
  # library(ruv)
})

# # Set conflict resolution preferences
# conflicted::conflicts_prefer(
#   base::intersect, base::setdiff, base::union,
#   magrittr::set_names,
#   dplyr::select, dplyr::filter, dplyr::mutate, dplyr::rename,
#   dplyr::summarise, dplyr::summarize, dplyr::arrange, dplyr::desc, dplyr::slice,
#   dplyr::count, dplyr::group_by, dplyr::ungroup, dplyr::left_join,
#   dplyr::inner_join, dplyr::full_join, dplyr::semi_join, dplyr::anti_join,
#   dplyr::distinct, dplyr::relocate, dplyr::rename_with, dplyr::pull,
#   tidyr::pivot_longer, tidyr::pivot_wider
# )


# Universal variables -------
ggscipal_funs <- grep("^pal_", ls("package:ggsci"), value = TRUE)
all_colors <- lapply(ggscipal_funs, function(fname) {
  fun <- get(fname, envir = asNamespace("ggsci"))
  args <- names(formals(fun))
  
  # default for n=10
  n <- 10
  
  # if selecting palette, transport default
  if ("palette" %in% args) {
    tryCatch(fun(palette = "default")(n),
             error = function(e) NULL)
  } else {
    tryCatch(fun()(n), error = function(e) NULL)
  }
}) %>% unlist() %>% .[!is.na(.)] %>% unique()
set.seed(0); all_colors_shuffle <- sample(all_colors)

mycolors <- unique(c(ggsci::pal_npg()(10), ggsci::pal_d3()(10), ggsci::pal_cosmic()(10)))
pca_colors <- mycolors[1:22] %>% setNames(c('instrument:anatomical_classification', 'trans:anatomical_classification', 'batch:anatomical_classification', 'sample_type:anatomical_classification', 'new:anatomical_classification', 'anatomical_classification', 'batch:sample_type', 'trans:sample_type', 'instrument:sample_type', 'instrument:batch', 'instrument:trans', 'sample_type:new', 'batch:new', 'trans:batch', 'trans:new', 'batch', 'instrument:new', 'trans', 'sample_type', 'instrument', 'new', 'resid'))

default_shapes <- c(16, 15, 17, 18, 7, 8, 3, 4)

pool_N <- c('K20200623yuel_TPHP_DIA_pool1_Slot2-48_1_795', 'M20200623yuel_TPHP_DIA_pool_Slot2-48_1_640', 'K20200624yuel_TPHP_DIA_pool_Slot2-48_1_840', 'M20200624yuel_TPHP_DIA_pool_Slot2-48_1_655', 'M20200625yuel_TPHP_DIA_pool_Slot2-48_1_670', 'M20200625yuel_TPHP_DIA_pool_Slot2-48_1_686', 'K20200626yuel_TPHP_DIA_pool_Slot2-48_1_887', 'K20200627yuel_TPHP_DIA_pool_Slot2-48_1_905', 'M20200627yuel_TPHP_DIA_pool_Slot2-48_1_718', 'M20200628yuel_TPHP_DIA_pool_Slot2-48_1_733', 'K20200630yuel_TPHP_DIA_pool_Slot2-48_1_942', 'M20200630yuel_TPHP_DIA_pool_Slot2-48_1_764', 'K20200701yuel_TPHP_DIA_pool_Slot2-48_1_965', 'M20200701yuel_TPHP_DIA_pool_Slot2-48_1_783', 'K20200702yuel_TPHP_DIA_pool_Slot2-48_1_988', 'M20200702yuel_TPHP_DIA_pool_Slot2-48_1_808', 'K20200704yuel_TPHP_DIA_pool_Slot2-48_1_1004', 'K20200704yuel_TPHP_DIA_pool_Slot2-48_1_996', 'K20200705yuel_TPHP_DIA_pool_Slot2-48_1_1018', 'K20200706yuel_TPHP_DIA_pool_Slot2-48_1_1048', 'M20200706yuel_TPHP_DIA_pool_Slot2-48_1_845', 'M20200704yuel_TPHP_DIA_pool_Slot2-48_1_875', 'K20200706yuel_TPHP_DIA_pool_Slot2-48_1_1069', 'M20200713yuel_TPHP_DIA_pool_Slot2-48_1_984', 'K20200705yuel_TPHP_DIA_pool_Slot2-48_1_1195', 'K20200713yuel_TPHP_DIA_pool_Slot2-48_1_1186', 'K20200716yuel_TPHP_DIA_pool_Slot2-48_1_1213', 'M20200715yuel_TPHP_DIA_pool_Slot2-48_1_1023', 'K20200716yuel_TPHP_DIA_pool_Slot2-48_1_1228', 'K20200716yuel_TPHP_DIA_pool_Slot2-48_1_1249', 'M20200717yuel_TPHP_DIA_pool_Slot2-48_1_1048', 'K20200718yuel_TPHP_DIA_pool_Slot2-48_1_1261', 'M20200718yuel_TPHP_DIA_pool_Slot2-48_1_1072', 'K20200716yuel_TPHP_DIA_pool_Slot2-48_1_1276', 'M20200719yuel_TPHP_DIA_pool_Slot2-48_1_1088', 'K20200718yuel_TPHP_DIA_pool_Slot2-48_1_1293', 'M20200720yuel_TPHP_DIA_pool_Slot2-48_1_1111', 'M20200721yuel_TPHP_DIA_pool_Slot2-48_1_1127', 'K20200718yuel_TPHP_DIA_pool_Slot2-48_1_1341', 'M20200723yuel_TPHP_DIA__pool_Slot2-48_1_1164', 'M20200723yuel_TPHP_DIA_pool_Slot2-48_1_1153', 'M20200723yuel_TPHP_DIA__pool_Slot2-48_1_1169', 'M20200723yuel_TPHP_DIA_pool_Slot2-48_1_1184', 'M20200820yuel_TPHP_DIA_pool_Slot1-54_1_4192', 'M20210103yuel_TPHP_DIA_pool_Slot1-54_1_4236', 'M20210104yuel_TPHP_DIA_pool_Slot1-54_1_4267', 'M20210115yuel_TPHP_DIA_pool_Slot1-54_1_4430', 'M20210116yuel_TPHP_DIA_pool_Slot1-54_1_4448', 'M20210118yuel_TPHP_DIA_pool_Slot1-54_1_4497', 'N20210528yuel_TPHP_DIA_pool_Slot1-7_1_5731', 'N20210529yuel_TPHP_DIA_pool_Slot1-7_1_5758', 'N20210608yuel_TPHP_DIA_pool_Slot1-7_1_5918', 'N20210609yuel_TPHP_DIA_pool_Slot1-7_1_5945', 'N20210610yuel_TPHP_DIA_pool_Slot1-7_1_5961', 'N20210611yuel_TPHP_DIA_pool_Slot1-7_1_5984', 'N20210617yuel_TPHP_DIA_pool_Slot1-7_1_6087')
pancancer_low_pearson <- c('N20250708yuel_TPHP_RCA_90mDIA_b1_76_Slot2-27_1_34640')

# Project-specific variables -----
# prot-info
# dfprot <- rio::import('../0_process_DIA-NN_data/output/protein_info_from_reports.csv')
dfprot <- rio::import('../0_process_DIA-NN_data/output/protein_info_from_reports_V1009.csv')
libprot <- rio::import('../../TPL/libs/20220616_fraglib/protein.tsv')

# names to abbr.
df_abbr <- read.delim('../0_process_DIA-NN_data/input/labels/sample_types_abbr_20250909.txt', stringsAsFactors = F, check.names = F, na.strings = '')
nm2abbr <- df_abbr$'Abbr' %>% setNames(df_abbr$'Entire')

# organ/cancer colors
df_color <- rio::import('../0_process_DIA-NN_data/input/labels/PUH_tissue_colorset_20251009.xlsx')
organ_color <- str_c('#', df_color$color) %>%
  setNames(df_color$class_abbr) %>% .[!is.na(names(.))]

dfc_color <- rio::import('../0_process_DIA-NN_data/input/labels/PUH_cancer_colorset_20251009.xlsx')
cancer_color <- str_c('#', dfc_color$color) %>% setNames(dfc_color$cancer_abbr)
cancer_color['OC'] <- '#9D7DDB'

cluster_color <- c(P1 = '#CA1D1F', P2 = '#F4AA63', P3 = '#006835', P4 = '#7E2F8E', P5 = '#2E7CAF')
sample_type_color <- c(T = '#CF3F54', NT = '#99CFD1', N = '#03A2B3', `F` = '#B09C85')
sample_color <- c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#F39800') %>% 
  setNames(c('T', 'NT', 'N', 'F', 'Brain')) # for tSNE colors
tphp_color_list <- list(
  cancer <- cancer_color,
  ana_class = organ_color,
  sample = sample_color,
  pancancer_cluster = cluster_color
)
# saveRDS(tphp_color_list, 'tphp_color_list.rds')

.Dataset_palette <- function(n, shuffle_seed = 0L, include = NULL, exclude = NULL) {
  # Guardrails
  stopifnot(is.numeric(n), length(n) == 1, n >= 1)
  if (!requireNamespace("ggsci", quietly = TRUE)) {
    stop("Package 'ggsci' is required. Install it via install.packages('ggsci').")
  }
  
  # Discover all ggsci palette factories (pal_*)
  pal_funs <- grep("^pal_", ls("package:ggsci"), value = TRUE)
  
  # Optional include/exclude filters by function name (e.g., "pal_npg", "pal_jco")
  if (!is.null(include)) pal_funs <- intersect(pal_funs, include)
  if (!is.null(exclude)) pal_funs <- setdiff(pal_funs, exclude)
  
  # Collect colors from each pal_* (default n=10, palette="default" when available)
  collect_one <- function(fname) {
    fun <- get(fname, envir = asNamespace("ggsci"))
    args <- names(formals(fun))
    n_default <- 10L
    res <- tryCatch({
      if ("palette" %in% args) fun(palette = "default")(n_default) else fun()(n_default)
    }, error = function(e) NULL)
    res
  }
  
  all_cols <- unlist(lapply(pal_funs, collect_one), use.names = FALSE)
  all_cols <- unique(all_cols[!is.na(all_cols)])
  
  # Fallback if something odd happened
  if (!length(all_cols)) {
    if (requireNamespace("viridisLite", quietly = TRUE)) return(viridisLite::viridis(n))
    # Okabe–Ito in base R (R >= 4.0): colorblind-friendly fallback
    return(grDevices::hcl.colors(n, palette = "Okabe-Ito"))
  }
  
  # Deterministic shuffle without disturbing the user's RNG state
  has_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_rng <- if (has_rng) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (has_rng) {
      assign(".Random.seed", old_rng, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  
  set.seed(as.integer(shuffle_seed))
  all_cols <- sample(all_cols)
  
  # Ensure length n (repeat if n exceeds the aggregated pool)
  if (n > length(all_cols)) {
    all_cols <- rep(all_cols, length.out = n)
  } else {
    all_cols <- all_cols[seq_len(n)]
  }
  
  unname(all_cols)
}


# General functions ----------
.is_log2_like <- function(x, probs = c(0.05, 0.5, 0.95), upper = 40, spread = 20) {
  # Flatten lists; coerce non-numeric with warnings suppressed (caller is informed via logic below)
  if (!is.atomic(x)) x <- unlist(x, use.names = FALSE)
  if (!is.numeric(x)) x <- suppressWarnings(as.numeric(x))
  
  x <- x[is.finite(x)]
  if (!length(x)) return(FALSE)
  
  q <- stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  # Heuristic: median <= upper and (q95 - q05) <= spread
  (q[3] <= upper) && ((q[3] - q[1]) <= spread)
}

.to_log2 <- function(X, offset = 1) {
  transform_vec <- function(v) {
    if (!is.numeric(v)) v <- suppressWarnings(as.numeric(v))
    v_fin <- v[is.finite(v)]
    if (!length(v_fin)) {
      message("All values are non-finite or missing; skipping transform for this vector.")
      return(v)
    }
    if (.is_log2_like(v_fin)) {
      message("Already in logarithmic scale (heuristic); skipped transform.")
      return(v)
    }
    if (any(v < -offset, na.rm = TRUE)) {
      warning("Values < -offset detected; log2(x + offset) will produce NaN for those entries.")
    }
    log2(v + offset)
  }
  
  # Vector
  if (is.null(dim(X))) {
    return(transform_vec(X))
  }
  
  # Matrix: decide on whole-matrix heuristic, then transform elementwise
  if (is.matrix(X)) {
    if (.is_log2_like(as.vector(X))) {
      message("Already in logarithmic scale (heuristic, matrix); skipped transform.")
      return(X)
    }
    if (any(X < -offset, na.rm = TRUE)) {
      warning("Values < -offset detected in matrix; log2(x + offset) will produce NaN for those entries.")
    }
    storage.mode(X) <- "double"
    return(log2(X + offset))
  }
  
  # Data frame: transform numeric-coercible columns independently
  if (is.data.frame(X)) {
    out <- X
    # Identify columns that can be meaningfully coerced to numeric (ignoring all-NA coercions)
    can_num <- vapply(X, function(col) {
      if (is.numeric(col)) return(TRUE)
      tmp <- suppressWarnings(as.numeric(col))
      any(is.finite(tmp), na.rm = TRUE)
    }, logical(1))
    
    for (nm in names(X)[can_num]) {
      out[[nm]] <- transform_vec(X[[nm]])
    }
    return(out)
  }
  
  # Fallback: try vector transform
  transform_vec(X)
}

# # Detect & convert to log2 if needed (robust heuristic)
# .to_log2 <- function(X) {
#   X <- as.matrix(X)
#   v <- as.numeric(X[is.finite(X)])
#   if (!length(v)) stop("No finite values in input matrix.")
#   q10 <- suppressWarnings(stats::quantile(v, 0.10, names = FALSE))
#   q90 <- suppressWarnings(stats::quantile(v, 0.90, names = FALSE))
#   is_log_like <- is.finite(q10) && is.finite(q90) && (q90 < 40) && ((q90 - q10) < 30)
#   if (is_log_like) return(X)
#   Xlog <- X
#   Xlog[!is.finite(Xlog) | Xlog <= 0] <- NA_real_
#   Xlog <- log2(Xlog)
#   Xlog
# }

get_outliers <- function(vec, coef = 1.5){
  # outliers based on Q1 and Q3
  stats <- quantile(vec, na.rm = T)
  iqr <- diff(stats[c(2, 4)])
  ret <- c(stats[2] - coef * iqr, stats[4] + coef * iqr)
  return(ret)
}

cv <- function(x, na.rm = T) {
  sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)
}

hex_col <- function(ls){
  # transform color name from 0-255 decimal-format to hexadecimal-format
  for(e in ls){
    if((FALSE %in% unique(e %in% 0:255)) | (length(e) != 3)){
      stop('wrong decimal-format (0-255)')
    }
  }
  hex <- lapply(ls, function(vec){
    stringr::str_c('#', stringr::str_c(as.hexmode(vec), collapse = ''))
  })
  return(unlist(hex))
  
}

removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
removeRowsNa <- function(x, threshold = 0.5) {
  # stopifnot(is.data.frame(x))
  keep <- rowMeans(is.na(x)) <= threshold
  x[keep, , drop = FALSE]
}
removeColsNa <- function(x, threshold = 0.5) {
  # stopifnot(is.data.frame(x))
  keep <- colMeans(is.na(x)) <= threshold
  x[, keep, drop = FALSE]
}

get_random_colors <- function(n, from_color_palette = brewer.pal(8, "Set2")){
  # get random colors
  my_colors <- sample(from_color_palette, 1)
  i <- 1
  while(length(my_colors) < n){
    color <- sample(from_color_palette, 1)
    if(color != my_colors[i]){
      my_colors[i+1] <- color
      i <- i + 1
    }
  }
  return(my_colors)
}

export_plots_separated_legends <- function(plots,
                                      filename = "plots_with_legends.pdf",
                                      width = 6,
                                      height = 6,
                                      legends_mode = c("per_plot", "all_last")) {
  # plots: list of ggplot objects
  # filename: output pdf file
  # width/height: size in inches
  # legends_mode:
  #   "per_plot"  -> each plot (no legend) followed by its legend
  #   "all_last"  -> all plots (no legend), legends collected at end page(s)
  
  # library(ggplot2)
  # library(cowplot)
  # library(gridExtra)
  
  legends_mode <- match.arg(legends_mode)
  
  graphics.off()
  pdf(filename, width = width, height = height)
  
  if (legends_mode == "per_plot") {
    for (i in seq_along(plots)) {
      # Plot without legend
      p_no_legend <- plots[[i]] + theme(legend.position = "none")
      print(p_no_legend)
      
      # Extract legend
      legend <- cowplot::get_legend(plots[[i]] + theme(legend.position = "right"))
      gridExtra::grid.arrange(legend)  # print legend on new page
    }
  } else if (legends_mode == "all_last") {
    # First print all plots without legends
    for (i in seq_along(plots)) {
      p_no_legend <- plots[[i]] + theme(legend.position = "none")
      print(p_no_legend)
    }
    # Collect all legends
    legends <- lapply(plots, function(p) cowplot::get_legend(p + theme(legend.position = "right")))
    # Arrange all legends together on final page
    gridExtra::grid.arrange(grobs = legends, ncol = 1)
  }
  
  dev.off()
  message(sprintf("Saved to %s", filename))
}

make_color_shape_maps <- function(levels_vec,
                                  base_colors = all_colors,
                                  n_color_max = 100,
                                  shapes = default_shapes) {
  K <- length(levels_vec)
  
  # Ensure enough colors by repeating if needed
  if (length(base_colors) < n_color_max) {
    base_colors <- rep(base_colors, length.out = n_color_max)
  }
  
  # Color assignment: cycle through up to n_color_max
  color_idx <- ((seq_len(K) - 1) %% n_color_max) + 1
  color_map <- base_colors[color_idx]
  names(color_map) <- levels_vec
  
  # Shape assignment: group every n_color_max into one shape
  if (K > n_color_max) {
    n_groups <- ceiling(K / n_color_max)
    if (length(shapes) < n_groups) {
      shapes <- rep(shapes, length.out = n_groups)
    }
    shape_idx <- ((seq_len(K) - 1) %/% n_color_max) + 1
    shape_map <- shapes[shape_idx]
  } else {
    # All same shape if under threshold
    shape_map <- rep(16, K)
  }
  names(shape_map) <- levels_vec
  
  list(colors = color_map, shapes = shape_map)
}

merge_protein_matrices_dt <- function(matrix_list, sources = NULL) {
  if (length(matrix_list) == 0) return(NULL)
  if (is.null(sources)) sources <- paste0("S", seq_along(matrix_list))
  
  # Convert each matrix to long DT: protein, sample(colname), value
  dts <- Map(function(m, s) {
    if (is.null(rownames(m))) stop("Each matrix must have rownames (protein IDs).")
    DT <- as.data.table(m, keep.rownames = "protein")
    setDT(DT)
    # Melt to long, keeping sample names
    long <- melt(DT, id.vars = "protein", variable.name = "sample", value.name = "value")
    long[, sample := paste0(s, "_", as.character(sample))]
    long
  }, matrix_list, sources)
  
  all_long <- rbindlist(dts, use.names = TRUE, fill = TRUE)
  
  # Cast back to wide: one column per sample, full join on protein
  wide <- dcast(all_long, protein ~ sample, value.var = "value")
  # Use protein as rownames
  rn <- wide$protein
  wide[, protein := NULL]
  mat <- as.matrix(wide)
  rownames(mat) <- rn
  storage.mode(mat) <- "double"
  mat
}

read_excel_allsheets <- function(filename, tibble = FALSE) {
  # Read all sheets from excel table(s) --
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
read_excel_allsheets_from_list <- function(flist, tibble = FALSE) {
  lapply(flist, read_excel_allsheets) %>% do.call(c, .) %>% return
}

#' Write one data.frame or a list of data.frames to XLSX with font control
#'
#' @param x            A data.frame or a named/unnamed list of data.frames.
#' @param path         Output .xlsx path.
#' @param latin_font   Latin font family name (e.g., "Arial").
#' @param font_size    Numeric font size for body cells.
#' @param header_bold  Logical; bold column headers.
#' @param header_bg    Header background fill color (e.g., "#F2F2F2" or NULL for none).
#' @param wrap_text    Logical; wrap text in body cells.
#' @param as_table     Logical; write each data.frame as Excel Table (openxlsx only).
#' @param freeze_row   Integer; freeze panes below this row (e.g., 2 to freeze header).
#' @param auto_width   Logical; auto-fit column widths after writing.
#' @param engine       "openxlsx" (default) or "openxlsx2".
#' @param sheet_names  Optional character vector of sheet names; recycled if needed.
#' @return             (Invisible) path to the written .xlsx file.
#' @details
#' - This function sets the workbook "base font" and applies explicit cell styles.
#' - CJK (Chinese/Japanese/Korean) glyphs are typically displayed via Excel fallback
#'   when the specified Latin font has no CJK glyphs. The actual fallback font is
#'   OS/Excel dependent and not controlled by this function.
#' - For `openxlsx2`, the function sets the theme's latin base font; per-cell styling
#'   uses `wb_set_cell_style`.
.write_excel <- function(
    x, path,
    latin_font  = "Arial",
    font_size   = 11,
    header_bold = TRUE,
    header_bg   = NULL,
    wrap_text   = FALSE,
    as_table    = FALSE,
    freeze_row  = 2L,
    auto_width  = TRUE,
    engine      = c("openxlsx", "openxlsx2"),
    sheet_names = NULL
) {
  stopifnot(is.character(path), length(path) == 1L)
  engine <- match.arg(engine)
  
  # Normalize input to a named list of data.frames
  if (inherits(x, "data.frame")) {
    x_list <- list(x)
    if (is.null(sheet_names)) sheet_names <- "Sheet1"
  } else if (is.list(x) && all(vapply(x, function(d) inherits(d, "data.frame"), logical(1)))) {
    x_list <- x
    if (is.null(sheet_names)) {
      sheet_names <- names(x_list)
      if (is.null(sheet_names) || any(nzchar(sheet_names) == FALSE)) {
        sheet_names <- paste0("Sheet", seq_along(x_list))
      }
    }
  } else {
    stop("`x` must be a data.frame or a list of data.frames.")
  }
  # Ensure unique, Excel-safe sheet names (max 31 chars, no []:*?/\\)
  sanitize_name <- function(s) {
    s <- gsub("[\\[\\]:\\*\\?/\\\\]", "_", s)
    s <- substr(s, 1L, 31L)
    s[nchar(s) == 0] <- "Sheet"
    s
  }
  sheet_names <- sanitize_name(make.unique(sheet_names))
  
  if (engine == "openxlsx") {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' is required for engine='openxlsx'. Please install it.")
    }
    wb <- openxlsx::createWorkbook()
    # Set workbook base (default) font
    openxlsx::modifyBaseFont(wb,
                             fontName  = latin_font,
                             fontSize  = font_size,
                             fontColour = "black"
    )
    
    # Define header and body styles
    headerStyle <- openxlsx::createStyle(
      fontName = latin_font,
      fontSize = font_size,
      textDecoration = if (header_bold) "bold" else NULL,
      fgFill = header_bg,
      halign = "center",
      valign = "center",
      wrapText = TRUE
    )
    bodyStyle <- openxlsx::createStyle(
      fontName = latin_font,
      fontSize = font_size,
      wrapText = wrap_text,
      valign = "center"
    )
    
    # Write each data.frame to its own worksheet
    for (i in seq_along(x_list)) {
      df <- x_list[[i]]
      sh <- sheet_names[[i]]
      openxlsx::addWorksheet(wb, sh)
      
      if (isTRUE(as_table)) {
        # writeData with asTable applies a built-in table style; headerStyle not used by tables
        openxlsx::writeDataTable(wb, sh, df, tableStyle = "TableStyleMedium2", withFilter = TRUE)
      } else {
        openxlsx::writeData(wb, sh, df, headerStyle = headerStyle)
        # Apply body style to entire data region (excluding header)
        n_rows <- nrow(df)
        n_cols <- ncol(df)
        if (n_rows > 0 && n_cols > 0) {
          openxlsx::addStyle(
            wb, sh, style = bodyStyle,
            rows = 2:(n_rows + 1L), cols = 1:n_cols,
            gridExpand = TRUE, stack = TRUE
          )
        }
      }
      
      # Freeze panes below header if requested
      if (!is.null(freeze_row) && is.finite(freeze_row) && freeze_row >= 2L) {
        openxlsx::freezePane(wb, sh, firstRow = TRUE)
      }
      
      # Auto-fit column widths
      if (isTRUE(auto_width)) {
        openxlsx::setColWidths(wb, sh, cols = 1:max(1L, ncol(df)), widths = "auto")
      }
    }
    
    openxlsx::saveWorkbook(wb, file = path, overwrite = TRUE)
    return(invisible(path))
  }
  
  # ---- openxlsx2 branch ----
  if (!requireNamespace("openxlsx2", quietly = TRUE)) {
    stop("Package 'openxlsx2' is required for engine='openxlsx2'. Please install it.")
  }
  wb <- openxlsx2::wb_workbook()
  
  # Set base font for latin major/minor (theme)
  openxlsx2::wb_set_base_font(
    wb,
    font_name = latin_font,
    font_size = font_size
  )
  
  # Predefine header/body styles for openxlsx2
  # (wb_set_cell_style works with 'openxlsx2::create_style' or style lists)
  hdr_style <- openxlsx2::create_style(
    font = openxlsx2::wb_font(
      name = latin_font,
      b = isTRUE(header_bold),
      sz = font_size
    ),
    alignment = openxlsx2::wb_alignment(horizontal = "center", vertical = "center", wrap_text = TRUE),
    fill = if (!is.null(header_bg)) openxlsx2::wb_fill(type = "solid", fgColor = header_bg) else NULL
  )
  body_style <- openxlsx2::create_style(
    font = openxlsx2::wb_font(name = latin_font, sz = font_size),
    alignment = openxlsx2::wb_alignment(vertical = "center", wrap_text = isTRUE(wrap_text))
  )
  
  for (i in seq_along(x_list)) {
    df <- x_list[[i]]
    sh <- sheet_names[[i]]
    wb$add_worksheet(sh)
    
    # Write data (as table or plain)
    if (isTRUE(as_table)) {
      wb$add_data_table(sh, x = df, withFilter = TRUE, table_style = "TableStyleMedium2")
    } else {
      wb$add_data(sh, x = df, startRow = 1, startCol = 1, colNames = TRUE, rowNames = FALSE)
      # Apply header style (row 1)
      if (ncol(df) > 0) {
        openxlsx2::wb_set_cell_style(
          wb, sh,
          dims = openxlsx2::wb_dims(1, 1, 1, ncol(df)),
          style = hdr_style
        )
      }
      # Apply body style
      if (nrow(df) > 0 && ncol(df) > 0) {
        openxlsx2::wb_set_cell_style(
          wb, sh,
          dims = openxlsx2::wb_dims(2, 1, nrow(df) + 1L, ncol(df)),
          style = body_style
        )
      }
    }
    
    # Freeze header row
    if (!is.null(freeze_row) && is.finite(freeze_row) && freeze_row >= 2L) {
      wb$freeze_pane(sh, firstRow = TRUE)
    }
    
    # Auto-fit widths (heuristic by measuring string lengths)
    if (isTRUE(auto_width)) {
      # Estimate widths using max(nchar) per column (cap at, say, 60)
      est_width <- function(v, header) {
        max_len <- max(nchar(as.character(v %||% "")), na.rm = TRUE)
        max( min(max_len + 2L, 60L), nchar(header) + 2L )
      }
      `%||%` <- function(a, b) if (is.null(a)) b else a
      widths <- Map(est_width, as.list(df), colnames(df))
      wb$set_col_widths(sh, cols = seq_along(df), widths = unlist(widths))
    }
  }
  
  openxlsx2::wb_save(wb, file = path, overwrite = TRUE)
  invisible(path)
}

# --- Main function: PCA, t-SNE, UMAP with adaptive color+shape ---
beca.DR <- function(X,
                    meta_df,
                    id_col,
                    var_col,
                    date_col = NULL,
                    seed = NA,
                    n_color_max = 100,
                    shape_set = default_shapes) {
  # X: features (rows) × samples (columns) matrix
  # meta_df: data.frame with sample metadata (must contain id_col)
  # var_col: character vector of categorical columns to color by
  # date_col: optional date column to map with continuous scale
  # seed: random seed for t-SNE/UMAP
  # n_color_max: max distinct colors before using shapes
  # shape_set: vector of plotting symbols
  
  
  # Transpose so rows = samples
  t_expr <- t(X)
  
  # ---- PCA ----
  message("Computing PCA...")
  # if (!is.na(seed)) set.seed(seed)  # for reproducibility
  pca_res <- prcomp(t_expr, center = TRUE, scale. = TRUE)
  pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2]) %>% 
    rownames_to_column(id_col) %>% 
    inner_join(meta_df)
  
  # ---- t-SNE ----
  message("Computing t-SNE...")
  if (!is.na(seed)) set.seed(seed)  # for reproducibility
  tsne_res <- Rtsne(t_expr, dims = 2, perplexity = 30, verbose = FALSE)
  tsne_df <- data.frame(Dim1 = tsne_res$Y[,1], Dim2 = tsne_res$Y[,2], row.names = rownames(t_expr)) %>% 
    rownames_to_column(id_col) %>% 
    inner_join(meta_df)
  
  # ---- UMAP ----
  message("Computing UMAP...")
  # if (!is.na(seed)) set.seed(seed)  # for reproducibility
  umap_res <- umap(t_expr, config = umap.defaults, preserve.seed = T)
  umap_df <- data.frame(UM1 = umap_res$layout[,1], UM2 = umap_res$layout[,2]) %>% 
    rownames_to_column(id_col) %>% 
    inner_join(meta_df)
  
  # ---- Plotting helpers ----
  draw_categorical <- function(df, x, y, group_var, title_prefix) {
    levels_vec <- levels(as.factor(df[[group_var]]))
    maps <- make_color_shape_maps(levels_vec,
                                  base_colors = all_colors,
                                  n_color_max = n_color_max,
                                  shapes = c(16, 15, 17, 18, 7, 8))
    ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
      geom_point(aes(color = .data[[group_var]],
                     shape = .data[[group_var]]),
                 size = 2, alpha = 0.8) +
      scale_color_manual(values = maps$colors) +
      scale_shape_manual(values = maps$shapes) +
      labs(title = paste0(title_prefix, " grouped by ", group_var)) +
      theme_bw() +
      theme(text = element_text(size = 10))
  }
  
  draw_continuous <- function(df, x, y, date_var, title_prefix) {
    date_range <- range(df[[date_var]], na.rm = TRUE)
    breaks <- date_breaks("2 month")(date_range)
    ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
      geom_point(aes(color = .data[[date_var]]), size = 2, alpha = 0.8) +
      scale_color_viridis_c(option = "G",
                            begin  = 0.05,
                            end    = 0.95,
                            name   = date_var,
                            breaks = as.numeric(breaks),
                            labels = date_format("%Y-%m-%d")(breaks)) +
      labs(title = paste0(title_prefix, " by ", date_var)) +
      theme_bw() +
      theme(text = element_text(size = 10))
  }
  
  # ---- Generate PCA plots ----
  pca_var_explained <- summary(pca_res)$importance[2, 1:2] * 100
  xlabel_pca <- sprintf("PC1 (%.1f%%)", pca_var_explained[1])
  ylabel_pca <- sprintf("PC2 (%.1f%%)", pca_var_explained[2])
  
  plots_pca <- map(var_col, ~
                     draw_categorical(pca_df, "PC1", "PC2", .x, "PCA") +
                     labs(x = xlabel_pca, y = ylabel_pca)
  ) %>% set_names(var_col)
  
  if (!is.null(date_col)) {
    plots_pca[[date_col]] <- draw_continuous(pca_df, "PC1", "PC2", date_col, "PCA") +
      labs(x = xlabel_pca, y = ylabel_pca)
  }
  
  # ---- Generate t-SNE plots ----
  plots_tsne <- map(var_col, ~
                      draw_categorical(tsne_df, "Dim1", "Dim2", .x, "t-SNE")
  ) %>% set_names(var_col)
  
  if (!is.null(date_col)) {
    plots_tsne[[date_col]] <- draw_continuous(tsne_df, "Dim1", "Dim2", date_col, "t-SNE")
  }
  
  # ---- Generate UMAP plots ----
  plots_umap <- map(var_col, ~
                      draw_categorical(umap_df, "UM1", "UM2", .x, "UMAP")
  ) %>% set_names(var_col)
  
  if (!is.null(date_col)) {
    plots_umap[[date_col]] <- draw_continuous(umap_df, "UM1", "UM2", date_col, "UMAP")
  }
  
  # Return list of results and plots
  list(
    pca_res    = pca_res,    pca_df    = pca_df,    plots_pca    = plots_pca,
    tsne_res   = tsne_res,   tsne_df   = tsne_df,   plots_tsne   = plots_tsne,
    umap_res   = umap_res,   umap_df   = umap_df,   plots_umap   = plots_umap
  )
  
  # how to run ---
  # res_dr <- beca.DR(
  #   X,
  #   meta_df,
  #   id_col     = "file_id",
  #   var_col    = c("instrument", "trans", "batch", "sample_type"),
  #   date_col   = "DateTime",
  #   seed       = 10,
  #   n_color_max = 100,
  #   shape_set   = c(16, 15, 17, 18, 7, 8)
  # )
  # export_plots_separated_legends(res_dr$plots_pca, "pca_plots.pdf", width = 6, height = 6)
  # export_plots_separated_legends(res_dr$plots_tsne, "tsne_plots.pdf", width = 6, height = 6)
  # export_plots_separated_legends(res_dr$plots_umap, "umap_plots.pdf", width = 6, height = 6)
}

beca.DR_old <- function(X, meta_df, id_col, var_col, date_col = NULL, seed = NA){
  # Input:
  # X: @rows ~ features, @cols ~ files
  # ##
  # meta_df <- metadata %>%
  #   mutate(yearmonth = str_c(year, month),
  #          DateTime = with_tz(ymd_hms(DateTime, tz = 'Asia/Shanghai'), 'UTC')) %>% 
  #   select(file_id, instrument, trans, batch, year, new, sample_type, DateTime)
  # id_col <- 'file_id'
  # date_col <- 'DateTime'
  # var_col <- colnames(meta_df) %>% setdiff(c(id_col, date_col))
  # seed <- 10
  # 
  # res_dr_ls <- list()
  # for(i in seq_along(names(pm_list))){
  #   cat(i, '...\n')
  #   nm <- names(pm_list)[i]
  #   X <- pm_list[[nm]]
  #   
  #   res_dr_ls[[i]] <- beca.DR(X, meta_df, id_col, var_col, date_col, seed)
  # }
  # names(res_dr_ls) <- names(pm_list)
  #
  # Transpose to prepare data for PCA/t-SNE/UMAP (samples as rows, features as columns)
  t_expr <- t(X) # each row = one sample, each column = a feature
  
  # PCA computation
  print('PCA computation...')
  pca_res <- prcomp(t_expr, center = TRUE, scale. = FALSE)
  pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2]) %>% 
    rownames_to_column(id_col) %>% 
    inner_join(meta_df)
  
  # Plot PCA
  plots_pca <- lapply(seq_along(var_col), function(i){
    ggplot(pca_df) +
      aes(x = PC1, y = PC2) +
      geom_point(aes_string(color = var_col[i]), size = 2, alpha = 0.8) +
      labs(title = str_c("PCA: Samples colored by ", var_col[i]),
           x = str_c("PC1 (", round(summary(pca_res)$importance[2,1]*100,1), "%)"),
           y = str_c("PC2 (", round(summary(pca_res)$importance[2,2]*100,1), "%)")) +
      scale_color_manual(values = mycolors) +
      theme_bw() +
      theme(text = element_text(size = 10))
  }) %>% setNames(var_col)
  if(!is.null(date_col)){
    plots_pca$DateTime <- ggplot(pca_df) +
      aes(x = PC1, y = PC2) +
      geom_point(aes_string(color = date_col), size = 2, alpha = 0.8) +
      labs(title = str_c("PCA: Samples colored by ", date_col),
           x = str_c("PC1 (", round(summary(pca_res)$importance[2,1]*100,1), "%)"),
           y = str_c("PC2 (", round(summary(pca_res)$importance[2,2]*100,1), "%)")) +
      scale_color_viridis_c(
        option = 'G', begin = 0.05, end = 0.95,
        name   = "Sample date",
        breaks  = as.numeric(date_breaks("2 month")(range(pca_df$DateTime))),
        labels  = date_format("%Y-%m-%d")(date_breaks("2 month")(range(pca_df$DateTime)))
      ) +
      theme_bw() +
      theme(text = element_text(size = 10))
  }
  # ggpubr::ggarrange(plotlist = plots_pca)
  
  # t-SNE computation (perplexity adjusted for dataset size)
  print('t-SNE computation...')
  if(is.integer(seed)) set.seed(seed)  # for reproducibility
  tsne_res <- Rtsne(t_expr, dims = 2, perplexity = 30, verbose = FALSE)
  tsne_df <- data.frame(Dim1 = tsne_res$Y[,1], Dim2 = tsne_res$Y[,2], row.names = rownames(t_expr)) %>% 
    rownames_to_column(id_col) %>% 
    inner_join(meta_df)
  # Plot t-SNE
  plots_tsne <- lapply(seq_along(var_col), function(i){
    ggplot(tsne_df) +
      aes(x = Dim1, y = Dim2) +
      geom_point(aes_string(color = var_col[i]), size = 2, alpha = 0.8) +
      labs(title = str_c("t-SNE: Samples colored by ", date_col), x = "t-SNE Dim1", y = "t-SNE Dim2") +
      scale_color_manual(values = mycolors) +
      theme_bw() +
      theme(text = element_text(size = 10))
  }) %>% setNames(var_col)
  if(!is.null(date_col)){
    plots_tsne$DateTime <- ggplot(tsne_df) +
      aes(x = Dim1, y = Dim2) +
      geom_point(aes_string(color = date_col), size = 2, alpha = 0.8) +
      labs(title = str_c("t-SNE: Samples colored by ", date_col), x = "t-SNE Dim1", y = "t-SNE Dim2") +
      scale_color_viridis_c(
        option = 'G', begin = 0.05, end = 0.95,
        name   = "Sample date",
        breaks  = as.numeric(date_breaks("2 month")(range(pca_df$DateTime))),
        labels  = date_format("%Y-%m-%d")(date_breaks("2 month")(range(pca_df$DateTime)))
      ) +
      theme_bw() +
      theme(text = element_text(size = 10))
  }
  # ggpubr::ggarrange(plotlist = plots_tsne)
  
  # UMAP computation
  print('UMAP computation...')
  if(is.integer(seed)) set.seed(seed)  # for reproducibility
  umap_res <- umap(t_expr)             # default 2D UMAP
  umap_df <- data.frame(UM1 = umap_res$layout[,1], UM2 = umap_res$layout[,2]) %>% 
    rownames_to_column(id_col) %>% 
    inner_join(meta_df)
  # Plot UMAP
  plots_umap <- lapply(seq_along(var_col), function(i){
    ggplot(umap_df) +
      aes(x = UM1, y = UM2) +
      geom_point(aes_string(color = var_col[i]), size = 2, alpha = 0.8) +
      labs(title = str_c("UMAP: Samples colored by ", date_col), x = "UMAP1", y = "UMAP2") +
      scale_color_manual(values = mycolors) +
      theme_bw() +
      theme(text = element_text(size = 10))
  }) %>% setNames(var_col)
  if(!is.null(date_col)){
    plots_umap$DateTime <- ggplot(umap_df) +
      aes(x = UM1, y = UM2) +
      geom_point(aes_string(color = date_col), size = 2, alpha = 0.8) +
      labs(title = str_c("UMAP: Samples colored by ", date_col), x = "UMAP1", y = "UMAP2") +
      scale_color_viridis_c(
        option = 'G', begin = 0.05, end = 0.95,
        name   = "Sample date",
        breaks  = as.numeric(date_breaks("2 month")(range(umap_df$DateTime))),
        labels  = date_format("%Y-%m-%d")(date_breaks("2 month")(range(umap_df$DateTime)))
      ) +
      theme_bw() +
      theme(text = element_text(size = 10))
  }
  # ggpubr::ggarrange(plotlist = plots_umap)
  
  ret <- list(pca_res = pca_res, pca_df = pca_df, plots_pca = plots_pca,
              tsne_res = tsne_res, tsne_df = tsne_df, plots_tsne = plots_tsne,
              umap_res = umap_res, umap_df = umap_df, plots_umap = plots_umap)
  return(ret)
}

#' ge.plot: run dimension reduction once or over a parameter grid and write multiple PDFs
#' Dependencies: Rtsne, RColorBrewer
#'
#' @param data   numeric matrix/data.frame, features x samples (rows=features, cols=samples)
#' @param type   factor/character length = ncol(data), class per sample
#' @param title  base title for output files (will be sanitized)
#' @param label  optional character vector of point labels (length = ncol(data))
#' @param seed   integer base seed; for grid mode uses seed + i - 1
#' @param grid   data.frame of parameters. Supported column names (any subset):
#'               perplexity, eta, exaggeration, exaggeration_factor, theta, max_iter, initial_dims
#'               (exaggeration / exaggeration_factor 都可；若两者皆缺，则用默认)
#' @param log1p  logical, whether to apply log1p to data before scaling
#' @param center logical, center features when scaling (row-wise)
#' @param scale_features logical, z-score scale each feature row-wise
#' @param n_pcs  integer or NULL. If set, do PCA on samples and feed top PCs to t-SNE (pca=FALSE).
#' @param outdir output directory (created if missing)
#' @return list(embeddings=<list of matrices>, params=<data.frame>, call=<match.call()>)
#'         In single-run mode, also writes one PDF like before.
ge.plot <- function(
    data, type, title = "", label = NULL, seed = 2023,
    grid = NULL,
    log1p = FALSE, center = TRUE, scale_features = TRUE,
    n_pcs = NULL, outdir = "."
) {
  # ---- checks ----
  stopifnot(ncol(data) == length(type))
  if (!is.null(label)) stopifnot(length(label) == ncol(data))
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- sanitize title -> file prefix ----
  base_title <- if (nzchar(title)) title else "TSNE"
  base_title <- gsub("[^A-Za-z0-9._-]+", "_", base_title)
  
  # ---- preprocessing: rows=features, cols=samples -> samples x features ----
  df <- as.matrix(data)
  if (log1p) df <- log1p(df)
  df[is.na(df)] <- 0
  
  if (scale_features) {
    # scale each feature across samples: result has same dim as df
    df <- t(apply(df, 1, function(x) as.numeric(scale(x, center = center, scale = TRUE))))
  }
  colnames(df) <- type
  X_sf <- t(df)  # samples x features (what Rtsne expects)
  
  # Optional PCA (on samples)
  if (!is.null(n_pcs)) {
    n_pcs <- max(1L, min(n_pcs, min(nrow(X_sf), ncol(X_sf))))
    pc <- stats::prcomp(X_sf, center = FALSE, scale. = FALSE)
    X_in <- pc$x[, seq_len(n_pcs), drop = FALSE]  # samples x n_pcs
    use_pca <- FALSE
    init_dims <- ncol(X_in)
  } else {
    X_in <- X_sf
    use_pca <- TRUE     # let Rtsne do its internal PCA
    init_dims <- min(50L, ncol(X_in))  # Rtsne default
  }
  
  # ---- color mapping for classes ----
  uniq_types <- unique(type)
  k <- length(uniq_types)
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    pal <- grDevices::rainbow(k)
  } else {
    max_set <- min(max(3, k), 12)
    pal <- RColorBrewer::brewer.pal(max_set, "Set3")
    if (k > max_set) {
      # extend palette deterministically
      pal <- grDevices::colorRampPalette(pal)(k)
    } else {
      pal <- pal[seq_len(k)]
    }
  }
  col_map <- stats::setNames(pal, uniq_types)
  col_vec <- unname(col_map[type])
  
  # ---- default grid if not provided (backward compatibility) ----
  if (is.null(grid)) {
    grid <- data.frame(
      perplexity = max(5, floor((nrow(X_in) - 1) / 3)),
      eta = 200, theta = 0.5, max_iter = 1000,
      exaggeration = 12, stringsAsFactors = FALSE
    )
  } else {
    grid <- as.data.frame(grid, stringsAsFactors = FALSE)
  }
  
  # helper: pick column if present else default
  pick <- function(df, i, nm, default) {
    if (nm %in% names(df)) df[i, nm] else default
  }
  
  # results containers
  Y_list <- vector("list", nrow(grid))
  pdf_paths <- character(nrow(grid))
  seeds <- integer(nrow(grid))
  
  n <- nrow(X_in)  # samples
  safe_perp_upper <- max(2, floor((n - 1) / 3))  # conservative upper bound
  
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, , drop = FALSE]
    perp <- min(as.numeric(pick(grid, i, "perplexity", safe_perp_upper)), safe_perp_upper)
    eta  <- as.numeric(pick(grid, i, "eta", 200))
    theta <- as.numeric(pick(grid, i, "theta", 0.5))
    max_iter <- as.integer(pick(grid, i, "max_iter", 1000))
    # accept either 'exaggeration' or 'exaggeration_factor'
    exag <- if ("exaggeration" %in% names(grid)) {
      as.numeric(grid[i, "exaggeration"])
    } else if ("exaggeration_factor" %in% names(grid)) {
      as.numeric(grid[i, "exaggeration_factor"])
    } else 12
    
    # allow overriding internal PCA initial_dims from grid (only used when use_pca=TRUE)
    init_dims_i <- if ("initial_dims" %in% names(grid)) {
      as.integer(grid[i, "initial_dims"])
    } else init_dims
    init_dims_i <- max(2L, min(init_dims_i, ncol(X_in)))
    
    seed_i <- seed + i - 1L
    set.seed(seed_i)
    
    ts <- Rtsne::Rtsne(
      X_in,
      dims                = 2L,
      perplexity          = perp,
      theta               = theta,
      max_iter            = max_iter,
      eta                 = eta,
      exaggeration_factor = exag,
      check_duplicates    = FALSE,
      pca                 = use_pca,
      pca_center          = FALSE,
      pca_scale           = FALSE,
      initial_dims        = init_dims_i,
      verbose             = FALSE
    )
    
    # ---- plotting ----
    # file name with params
    file_tag <- sprintf(
      "perp%s_eta%s_exag%s_theta%s_iter%s",
      format(perp, trim = TRUE, scientific = FALSE),
      format(eta, trim = TRUE, scientific = FALSE),
      format(exag, trim = TRUE, scientific = FALSE),
      format(theta, trim = TRUE, scientific = FALSE),
      max_iter
    )
    pdf_path <- file.path(outdir, sprintf("%s_TSNE_%s.pdf", base_title, file_tag))
    grDevices::pdf(pdf_path, height = 7, width = 10)
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
    
    if (is.null(label)) {
      graphics::par(mai = c(0.4, 0.5, 0.2, 2.4))
      graphics::plot(
        ts$Y,
        col = col_vec, main = "t-SNE", pch = 16, cex = 0.6,
        xlab = "t-SNE 1", ylab = "t-SNE 2",
        cex.axis = 1.2, cex.lab = 1.2
      )
      graphics::legend(
        "topright",
        legend = names(col_map),
        fill = unname(col_map), bty = "n", cex = 0.9
      )
    } else {
      graphics::plot(
        ts$Y, col = col_vec, main = "t-SNE",
        pch = 20, cex = 1.6, xlab = "t-SNE 1", ylab = "t-SNE 2",
        cex.axis = 1.2, cex.lab = 1.2
      )
      graphics::text(ts$Y, pos = 1, labels = label, col = "grey30", cex = 0.8)
      graphics::legend("topright", legend = names(col_map), fill = unname(col_map),
                       lty = 1, lwd = 1, bty = "n")
    }
    grDevices::dev.off()
    
    # collect
    Y_list[[i]] <- ts$Y
    pdf_paths[i] <- pdf_path
    seeds[i] <- seed_i
  }
  
  # assemble params+outputs
  # normalize column names (use 'exaggeration' in output)
  if ("exaggeration_factor" %in% names(grid) && !("exaggeration" %in% names(grid))) {
    grid$exaggeration <- grid$exaggeration_factor
    grid$exaggeration_factor <- NULL
  }
  out_df <- grid
  out_df$seed <- seeds
  out_df$pdf <- pdf_paths
  
  return(invisible(list(
    embeddings = Y_list,   # list, same order as rows of grid
    params     = out_df,   # data.frame with pdf paths
    call       = match.call()
  )))
}




# Project-specific functions -------------
get_abbr <- function(df, sample_type, df_abbr = NULL){
  # get abbreviation of sample types, return input name which not matched
  # df_abbr$Entrie should be total lowercase !!!!
  if(is.null(df_abbr)){
    df_abbr <- read.delim('//172.16.13.114/share/members/jiangwenhao/TPHP/input/sample_types_abbr_20220722.txt', stringsAsFactors = F, check.names = F, na.strings = '')
    # df_abbr[is.na(df_abbr)] <- 'NA'
  }
  h <- hash::hash(keys = df_abbr$'Entire',
                  values = df_abbr$'Abbr')
  df[, sample_type] <- unlist(lapply(df[, sample_type], function(e){
    if(is.null(h[[tolower(e)]])){
      return(e)
    }else{
      return(h[[tolower(e)]])
    }
  }))
  return(df)
}

my_plot <- function(df){
  # Scatter plot of protein abundance
  df %<>%
    add_column(label = stringr::str_c(df$organ, ' (', df$cancer_type, ')'), .before = 1) %>%
    dplyr::mutate(cancer_type = factor(cancer_type, levels = c('N', 'Adj', 'C'))) %>%
    arrange(organ, cancer_type) %>%
    dplyr::mutate(cancer_type = as.character(cancer_type))
  
  #remove organ with all subtypes >= 90% NA ratio
  na_ratio <- function(x) sum(is.na(x)) / length(x)
  organ_selected <- df %>% dplyr::group_by(organ, cancer_type) %>%
    summarise_at(vars(2), list(na_ratio = 'na_ratio')) %>%
    ungroup() %>%
    filter(na_ratio < 0.9) %>%
    pull(organ) %>%
    unique()
  
  df %<>% filter(organ %in% organ_selected)
  
  rlt <- list()
  #df[is.na(df)] <- 6.759341
  df_fillrate <- df
  my_uni <- colnames(df)[ncol(df)]
  my_title <- my_uniprot2prot(my_uni)
  colnames(df_fillrate)[ncol(df_fillrate)] <- colnames(df)[ncol(df)] <- c('intensity')
  
  
  df <- na.omit(df)
  df_fillrate <- df_fillrate[df_fillrate$organ %in% unique(df$organ), c(2, 3, 4)]
  
  # color setting
  df$intensity <- log10(2 ^ df$intensity) # log2 -> log10
  # pseudo_ht_color <- c('#F01025', '#EAD41A', '#1E72BA')
  pseudo_ht_color <- c('#F01025', '#EAD41A', '#1E72BA')
  names(pseudo_ht_color) <- c('C', 'Adj', 'N')
  my_fills <- pseudo_ht_color[df$cancer_type]
  names(my_fills) <- df$cancer_type
  
  # text plotting
  my_texts_pos <- 1:length(unique(df$organ)) - 0.5
  my_vlines <- 1:length(unique(df$organ))
  x <- c()
  for(i in 1:length(unique(df$organ))){
    df_tmp <- df[df$organ == unique(df$organ)[i],]
    x_tmp <- seq(from = i-1, to = i, length.out = nrow(df_tmp)+2)[2:(nrow(df_tmp)+1)]
    x <- append(x, x_tmp)
  }
  df$x <- x
  
  # barplot preparation
  df_bar <- lapply(unique(df_fillrate$organ), function(e){
    df_tmp <- df_fillrate[df_fillrate$organ == e, c('cancer_type', 'intensity')]
    df_ttmp <- lapply(unique(df_tmp$cancer_type), function(ee){
      vec_ttmp <- df_tmp$intensity[df_tmp$cancer_type == ee]
      return(sum(!is.na(vec_ttmp)) / length(vec_ttmp))
    }) %>% as.data.frame(stringsAsFactor = F)
    colnames(df_ttmp) <- unique(df_tmp$cancer_type)
    rownames(df_ttmp) <- e
    return(df_ttmp)
  }) %>% do.call(plyr::rbind.fill, .)
  df_bar$organ <- unique(df_fillrate$organ)
  df_bar %<>% reshape2::melt()
  colnames(df_bar)[2:3] <- c('cancer_type', 'fillingvalue')
  df_bar$cancer_type <- factor(df_bar$cancer_type, levels = c('C', 'Adj', 'N'), ordered = T)
  rlt$table <- df_bar %<>% dplyr::arrange(organ, cancer_type)
  df_bar %<>% na.omit
  df_bar$label <- stringr::str_c(df_bar$organ, df_bar$cancer_type, sep = '_')
  pseudo_ht_bar <- seq(from = 0, to = length(unique(df_bar$organ)), by = 1/4) %>% .[. %% 1 != 0]
  names(pseudo_ht_bar) <- stringr::str_c(lapply(unique(df_bar$organ), function(e) rep(e, 3)) %>% unlist, c('N', 'Adj', 'C'), sep = '_')
  
  df_bar$x <- pseudo_ht_bar[df_bar$label]
  df_bar <- df_bar[df_bar$fillingvalue != 0, ]
  nmax <- max(df$intensity); nmin <- min(df$intensity)
  #df_bar$fillingvalue_lt <- my_linear_trans(c(df_bar$fillingvalue, 0.5), nmax, nmin) %>% .[1:(length(.) - 1)]
  df_bar$fillingvalue_lt <- df_bar$fillingvalue * (nmax - nmin) + nmin
  # double coordinates
  fillingvalue_x <- max(df$x) + 0.5
  fillingvalue_scale <- seq(nmin, nmax, length.out = 6)
  fillingvalue_str <- stringr::str_c(seq(0, 100, length.out = 6), '%', sep = '')
  fillingtitle <- '1-missing rate'
  fillingtitle_pos <- mean(c(nmin, nmax))
  #hline0.5_lt <- my_linear_trans(c(df_bar$fillingvalue, 0.5), nmax, nmin) %>% .[length(.)]
  
  # figure
  tmp <- stringr::str_c("+annotate('text', x = ", my_texts_pos, ", y = max(df$intensity) * 1.02, label = '", unique(sort(df$organ)), "', color = '#000000', size = 4)", collapse = ' ') # labels of sample type
  p <- ggplot()+
    geom_bar(data = df_bar, mapping = aes(x = x, y = fillingvalue_lt, fill = cancer_type), fill = pseudo_ht_color[df_bar$cancer_type], stat = 'identity', alpha = 0.2, width = 0.25)+
    geom_bar(data = df_bar, mapping = aes(x = x, y = nmin, fill = cancer_type), fill = '#FFFFFF', stat = 'identity', width = 0.25)+# block non-missing rate lower than 0
    geom_point(data = df, mapping = aes(x = x, y = intensity), fill = my_fills, color = '#000000', shape = 21, size = 5)+
    geom_vline(xintercept = my_vlines, color = '#000000', linetype = 'dashed', size = 0.5)+
    #geom_hline(yintercept = hline0.5_lt, color = pseudo_ht_color['C'], linetype = 'dashed', size = 0.5, alpha = 0.5)+
    #annotate('text', x = 0.2, y = hline0.5_lt * 1.01, label = 'Missing value < 50%', color = pseudo_ht_color['C'], size = 4, alpha = 0.5, hjust = 0)+
    annotate('text', x = fillingvalue_x * 1.06, y = fillingtitle_pos, label = fillingtitle, size = 4.5, hjust = 0.5, vjust = 0.5, angle = 270)+# mock y-axis title
    annotate('text', x = fillingvalue_x * 1.01, y = fillingvalue_scale, label = fillingvalue_str, size = 4, hjust = 0, vjust = 0.5)+# mock y-axis text
    labs(
      y = 'log10 Intensity', title = my_title
    )+
    geom_vline(xintercept = c(placeholder_x = fillingvalue_x * 1.1), color = '#FFFFFF')+# placeholder for x-axis
    coord_cartesian(ylim = c(nmin, nmax * 1.02))+
    scale_x_continuous(expand = c(0, 0))+
    theme(panel.grid = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 12,color = 'black'),
          axis.line = element_line(color = 'black'),
          axis.line.x = element_blank(),
          plot.subtitle = element_text(size = 30, hjust = 0, color = 'black')
    )+
    theme(legend.text = element_text(size = 12, color = 'black'), legend.position = 'top',
          legend.title = element_text(size = 15, color = 'black'))+
    theme(axis.title.x = element_blank(),
          #axis.text.x = element_text(size = 12, vjust = 0.5 ,color = 'black', angle = 90),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
    )+
    theme(axis.title.y = element_text(size = 12, hjust = 0.5, color = 'black', angle = 90),
          axis.text.y = element_text(size = 10,color = 'black', angle = 0),
          axis.ticks.y = element_blank()
    )
  rlt$p <- eval(parse(text = stringr::str_c('p', tmp)))
  return(rlt)
}

info_of_pm <- function(pm){
  info <- data.frame(Filename = colnames(pm))
  info$instrument <- as.character(stringr::str_extract(colnames(pm), '^[A-Z]+'))
  info$date <- as.character(stringr::str_extract(colnames(pm), '[0-9]+')) %>% factor(., levels = unique(sort(.)), ordered = T)
  info$year <- stringr::str_sub(info$date, end = 4)
  info$month <- stringr::str_sub(info$date, end = -3)
  
  info$batch <- NA
  info$batch[grep("M202006|K202006|K202007|M202007|K202008|M202008|N202008", info$Filename)] <- 'b1'
  info$batch[grep("N202011", info$Filename)] <- 'b2'
  info$batch[grep("M202101", info$Filename)] <- 'b3'
  info$batch[grep("N202104", info$Filename)] <- 'b4'
  info$batch[grep("N202105|N202106", info$Filename)] <- 'b5'
  info$batch[is.na(info$batch)] <- 'b6'
  
  info$trans <- NA
  info$trans[grep("b1|b2", info$batch)] <- "T1"
  info$trans[grep("b3|b4", info$batch)] <- "T2"
  info$trans[grep("b5", info$batch)] <- "T3"
  info$trans[grep("b6", info$batch)] <- "T4"
  info$trans[grep("^N202203", info$Filename)] <- "T5"
  
  info$new <- ifelse(grepl('^2022', info$month), 'new', 'old')
  return(info)
}

create_dataframe_from_list <- function(named_list, method = 'extend') {
  method_list <- c('extend', 'trim')
  if (!(method %in% method_list)){
    stop('method should be extend or trim')
  }
  if (method == 'extend'){
    max_length <- max(sapply(named_list, length))
    new_list <- lapply(named_list, function(v){
      c(v, rep(NA, max_length - length(v)))
    })
  } else if (method == 'trim'){
    min_length <- min(sapply(named_list, length))
    new_list <- lapply(named_list, function(v) head(v, min_length))
  }
  
  result_df <- as.data.frame(new_list)
  return(result_df)
}


easyp2p <- function(pep_data, roundn = NULL){
  # Summarise peptide matrix to protein --
  # @peptide matrix (row: peptides; column: files)
  #   - peptide sequence as the first column
  #   - protein id as the second column
  #   - data should be in raw scale
  
  # Step 1
  t1 <- proc.time()
  cat("1) Data preparation: \n"); print(proc.time() - t1)
  
  pep_data[is.infinite(as.matrix(pep_data))] <- NA # assign infinity as NA
  pep_data[is.na(pep_data)] <- NA
  pep_data <- pep_data[complete.cases(pep_data[, 1]), ] # remove NA peptides
  pep_data <- pep_data[!grepl("^1/CON", pep_data[, 2], fixed = F), ] # remove contaminants
  pep_data[pep_data == 0] <- NA # assign zeros as NA
  
  
  # Step 2
  cat("2) log2 transformation: \n"); print(proc.time() - t1)
  pep_data_log2 <- log2(pep_data[, 3:ncol(pep_data), drop = F]) %>% tibble::add_column(., prot = pep_data$prot, .before = 1)
  rownames(pep_data_log2) <- pep_data[, 1]
  
  
  # Step 3
  cat("3) Quantile normalization: \n"); print(proc.time() - t1)
  pep_data_log2_qn <- preprocessCore::normalize.quantiles(as.matrix(pep_data_log2))
  colnames(pep_data_log2_qn) <- colnames(pep_data_log2)
  rownames(pep_data_log2_qn) <- rownames(pep_data_log2)
  
  data_tech_rep <- cbind(pep_data[, 1:2], pep_data_log2_qn)
  #is.null(tech_rep_f)
  #is.null(batchf)
  
  rm(list = ls()[grep('^pep_data', ls())])
  
  
  # Step 4
  cat("4) Arrangement: \n"); print(proc.time() - t1)
  data <- data_tech_rep
  colnames(data)[1:2] <- c("tg", "prot")
  
  n <- ncol(data)
  pep2 <- apply(data[, -c(1, 2), drop = F], 1, function(x) { # log2 then mean of all files
    NAs <- length(which(is.na(x)))
    meanexpr1 <- sum(as.numeric(x), na.rm = TRUE) / (n - NAs)
    meanexpr2 <- sum(as.numeric(x), na.rm = TRUE) / n
    d <- c(NAs, meanexpr1, meanexpr2)
    return(d)
  })
  pep2 <- t(pep2)
  colnames(pep2) <- c("NAs", "meanexpr1", "meanexpr2")
  pep_expr <- cbind(data[, 1], pep2, data[, c(-1)])
  
  #order by pg ,#NA,intesity
  pep_order <- pep_expr[order(pep_expr[, 5], pep_expr[, 2], -pep_expr[, 3]), ]
  colnames(pep_order)[1] <- "tg"
  pep_order2 <- pep_order[, c(-2, -3, -4)]
  
  rm(list = ls()[grep('^data', ls())])
  rm(list = c('pep2', 'pep_expr', 'pep_order'))
  gc() # Garbage Collection
  
  # Step 5
  cat("5) Top 3 peptides selection: \n"); print(proc.time() - t1)
  pep_order2_dup1 <- pep_order2[duplicated(pep_order2$prot), ]
  pep_order2_dup2 <- pep_order2_dup1[duplicated(pep_order2_dup1$prot), ]
  pep_order2_dup3 <- pep_order2_dup2[duplicated(pep_order2_dup2$prot), ]
  pep_order2_top3 <- pep_order2[!(rownames(pep_order2) %in% rownames(pep_order2_dup3)), ]
  
  rm(list = ls()[grep('^pep_order2_dup', ls())])
  
  
  # Step 6
  cat("6) Protein matrix generation: \n"); print(proc.time() - t1)
  pep_order2_top3 <- pep_order2_top3[c("prot", "tg", colnames(pep_order2_top3)[3:ncol(pep_order2_top3)])]
  pep_order2_top3[pep_order2_top3 == 0] <- NA
  
  lr_top3 <- "top3"
  if(lr_top3 == "top3"){ # mean of top 3
    top3_mean <- plyr::ddply(pep_order2_top3,
                             .variables = "prot",
                             .fun = function(df_sub){
                               mean_ls <- colMeans(df_sub[, -c(1, 2), drop = F], na.rm = T)
                               if(!is.null(roundn)){
                                 mean_ls <- round(mean_ls, roundn)
                               }
                               return(mean_ls)
                             })
    df_prot <- top3_mean
    #readr::write_csv(top3_mean, 'prot_matrix_top3.csv', na = '')
    
  }else{ # LR
    prot_matrix <- pep2prot(pep_order2_top3)
    prot_matrix <- prot_matrix[, -2]
    
    df_prot <- prot_matrix
    #readr::write_csv(prot_matrix, 'prot_matrix_lr.csv', na = '')
  }
  
  return(df_prot)
  
}


prot2gene <- function(uniprot_ids, dataset = "hsapiens_gene_ensembl"){
  library(biomaRt)
  mart <- useMart("ensembl", dataset = dataset)
  uniprot_ids <- rownames(expr)
  mapping <- getBM(attributes = c("uniprotswissprot", "hgnc_symbol"),
                   filters = "uniprotswissprot",
                   values = uniprot_ids,
                   mart = mart) %>% setNames(c('protein', 'gene'))
  mapping %>%
    group_by(protein) %>%
    summarise(
      gene = str_c(unique(na.omit(gene)), collapse = "/"),
      .groups = "drop"
    ) %>% 
    as.data.frame() %>% 
    magrittr::set_rownames(.$protein) %>% 
    .[uniprot_ids, ]
}


my_uniprot2prot <- function(vec_uni){
  # crawl protein name from uniprot.org
  vec_prot <- lapply(vec_uni, function(uniprotid){
    # my_html <- rvest::read_html(stringr::str_glue('https://www.uniprot.org/uniprotkb/{uniprotid}/entry'))
    my_html <- rvest::read_html(stringr::str_glue('https://rest.uniprot.org/genecentric/{uniprotid}')) # 2022-08-08 update; get json
    my_json <- my_html %>%
      rvest::html_text() %>%
      jsonlite::fromJSON()
    
    prot <- my_json$canonicalProtein$proteinName
    return(prot)
  }) %>% unlist
  return(stringr::str_c(vec_uni, vec_prot, sep = '_'))
}

ulhen_class<-function(mat_median, fct = 3){
  # version 20251021 faster
  j=1
  ## precompute global minima used repeatedly
  .min_all <- min(mat_median[,-1, drop = FALSE], na.rm = TRUE)
  # .min_all <- min(mat_median[,-c(1, ncol(mat_median)), drop = FALSE], na.rm = TRUE)
  
  ## clas: one data.frame per assayed column, as before
  clas <- plyr::llply(seq_len(ncol(mat_median) - 1L), function(.k) {
    message('Analysis of protein No.', .k)
    num <- mat_median[[.k + 1L]]
    
    ## build t1 (no factors)
    t1 <- data.frame(tissue_type = mat_median[[1L]],
                     abundance   = num,
                     classification = NA_character_,
                     stringsAsFactors = FALSE)
    
    ## sort once
    o  <- order(t1$abundance, decreasing = TRUE, na.last = TRUE)
    t1 <- t1[o, , drop = FALSE]
    ab <- t1$abundance
    n  <- length(ab)
    
    ## -------- class 1 (not detected) + class 2 (tissue enriched), vectorized --------
    ## class 1
    t1$classification[ab <= .min_all] <- "not detected"
    
    ## class 2 (respect original NA behavior: if any NA exists in column, none get enriched)
    if (!anyNA(ab)) {
      ## since ab is sorted decreasing, “max of others” is ab[2] for i=1; otherwise ab[1]
      ## (O(1) per position instead of recomputing max(ab[-i]))
      max_others <- rep.int(ab[1L], n)
      if (n >= 2L) max_others[1L] <- ab[2L]
      t1$classification[ab >= (fct * max_others)] <- "tissue enriched"
    }
    
    c12 <- t1[!is.na(t1$classification), , drop = FALSE]
    
    ## -------- class 3 (group enriched), vectorized with cumulative sums --------
    max7 <- t1[is.na(t1$classification), , drop = FALSE]
    if (nrow(max7) > 7L) {
      cs    <- cumsum(max7$abundance)           # running sum of top k
      upto  <- seq_len(min(6L, nrow(max7) - 2L))# n in 1:6 where nextp exists
      if (length(upto)) {
        cond <- (cs[upto + 1L] / (upto + 1L)) >= (fct * max7$abundance[upto + 2L])
        if (any(cond)) {
          kstar <- max(which(cond)) + 1L        # largest n+1 satisfying condition
          max7$classification[seq_len(kstar)] <- "group enriched"
        }
      }
    } else if (nrow(max7) >= 2L && nrow(max7) <= 7L) {
      max7$classification <- "group enriched"
    }
    
    ## -------- class 4 (expressed in all tissues) --------
    if (!any(t1$abundance <= .min_all, na.rm = TRUE)) {
      max7$classification[is.na(max7$classification)] <- "Expressed in all tissues"
    }
    c34 <- stats::na.omit(max7)
    
    ## -------- class 5 (tissue enhanced) + class 6 (mixed) --------
    max7_na <- max7[is.na(max7$classification), , drop = FALSE]
    if (nrow(max7_na) >= 1L) {
      thr <- fct * mean(t1$abundance, na.rm = TRUE)
      max7_na$classification <- ifelse(max7_na$abundance >= thr,
                                       "tissue enhanced", "mixed")
    }
    
    com_max <- do.call(rbind, list(c12, c34, max7_na))
    
    ## keep your side-effect counter and progress reporter
    j <<- j + 1L
    colnames(com_max)[3L] <- colnames(mat_median)[j]
    # if (j %% 1000L == 0L) {
    #   cat(paste("Processed", j - 1L, "/", ncol(mat_median) - 1L), "\n")
    # }
    com_max[, c(1L, 3L), drop = FALSE]
  })
  
  print("Finished classification and start merging classification matrix")
  
  clas.all<-clas%>% reduce_right(full_join, by = "tissue_type")
  return(clas.all)
}
# ulhen_class<-function(mat_median, fct = 3){
#   # version 20251017
#   j=1
#   clas<-apply(mat_median[,-1], 2, function(num){
#     t1<-as.data.frame(matrix(ncol=3,nrow = nrow(mat_median)))
#     colnames(t1)<-c("tissue_type","abundance","classification")
#     t1$tissue_type<-mat_median[[1]]
#     t1$abundance<-num
#     t1<-t1[order(t1$abundance,decreasing = T),]
#     # class 1(not detected) and 2(tissue enriched)
#     for (i in 1:nrow(mat_median)){
#       if (t1$abundance[i]<=min(mat_median[,-1],na.rm = T)) {t1[i,3]<-"not detected"
#       } else if(t1$abundance[i]>=(fct*max(t1$abundance[-i]))) {t1[i,3]<-"tissue enriched"
#       } else {t1[i,3]<-NA}
#     }
#     c12<-t1[!is.na(t1$classification),]
#     
#     # class 3(group enriched)
#     max7<-t1[is.na(t1$classification),]
#     if (nrow(max7)>7){
#       high=max7$abundance[1]
#       for (n in 1:6){
#         high=high+max7$abundance[n+1]
#         nextp=max7$abundance[n+2]
#         # fill=max7$classification[1:(n+1)]
#         if (high/(n+1) >= fct*nextp)  {max7$classification[1:(n+1)]<-"group enriched"}
#       } 
#     } else if (nrow(max7)>=2 & nrow(max7)<=7) {max7$classification<-"group enriched"}
#     
#     #class 4 (expressed in all)
#     if (sum(as.logical(which(t1$abundance<=min(mat_median[,-c(1,ncol(mat_median))],na.rm = T))))==0){
#       max7$classification[is.na(max7$classification)]<-"Expressed in all tissues"}
#     c34<-na.omit(max7)
#     
#     #class 5 (tissue enriched) and class 6(mixed)
#     
#     max7_na<-max7[is.na(max7$classification),]
#     if (nrow(max7_na)>=1){
#       for (i in 1:nrow(max7_na)){
#         if(max7_na$abundance[i]>= (fct*mean(t1$abundance))) {max7_na$classification[i]<-"tissue enhanced"
#         } else {max7_na$classification[i]<-"mixed"}
#       }
#     }
#     com_max<-do.call(rbind,list(c12,c34,max7_na))
#     j <<- j+1
#     colnames(com_max)[3]<-colnames(mat_median)[j]
#     if (j %% 1000 ==0){
#       cat(paste("Processed",j-1,"/",ncol(mat_median)-1,sep = " "),"\n")}
#     
#     com_max[,c(1,3)]
#     
#     # colnames(class[[j]])[2]<-colnames(mat_median)[j+1]
#   })
#   
#   print("Finished classification and start merging classification matrix")
#   
#   clas.all<-clas%>% reduce_right(full_join, by = "tissue_type")
#   return(clas.all)
# }

ulhen_class_v2 <- function(mat_median, fct = 3){
  # version 20251017
  # Keep a counter consistent with original semantics
  j <- 1
  
  # Map over all expression columns (excluding the first tissue_type column)
  clas <- purrr::map(names(mat_median)[-1], function(col_nm){
    # Build working table: tissue_type + target column as abundance; initialize classification
    t1 <- tibble::tibble(
      tissue_type    = mat_median[[1]],
      abundance      = mat_median[[col_nm]],
      classification = NA_character_
    ) %>%
      dplyr::arrange(dplyr::desc(abundance))
    
    # Side-effect print (kept from original)
    print(t1)
    
    # Global minimum across all numeric value cells except the first id column
    global_min <- suppressWarnings(
      min(as.matrix(mat_median[, -1, drop = FALSE]), na.rm = TRUE)
    )
    
    # For each row, compute the maximum abundance among all *other* rows
    # (needed for the "tissue enriched" rule). Keep it vectorized with map_dbl.
    other_max <- purrr::map_dbl(seq_len(nrow(t1)), function(i){
      max(t1$abundance[-i], na.rm = TRUE)
    })
    
    # Class 1 (not detected) and Class 2 (tissue enriched)
    t1 <- t1 %>%
      dplyr::mutate(
        classification = dplyr::case_when(
          abundance <= global_min ~ "not detected",
          abundance >= fct * other_max ~ "tissue enriched",
          TRUE ~ classification
        )
      )
    
    c12 <- t1 %>% dplyr::filter(!is.na(classification))
    
    # Candidates for further classification (still NA)
    max7 <- t1 %>% dplyr::filter(is.na(classification))
    
    # Class 3 (group enriched)
    if (nrow(max7) > 7){
      # Iteratively evaluate whether the mean of the top k abundances (k in 1..7)
      # is >= fct * the (k+1)-th abundance; if so, mark top k as group enriched.
      high <- max7$abundance[1]
      for (n in 1:6){
        high  <- high + max7$abundance[n + 1]
        nextp <- max7$abundance[n + 2]
        if ( (high / (n + 1)) >= (fct * nextp) ){
          max7$classification[1:(n + 1)] <- "group enriched"
        }
      }
    } else if (nrow(max7) >= 2 && nrow(max7) <= 7){
      max7$classification <- "group enriched"
    }
    
    # Class 4 (Expressed in all tissues)
    # Condition preserved from original: check if no abundance is <= the minimum of all
    # matrix entries excluding the first (id) column and the last column.
    min_excl_last <- suppressWarnings(
      min(as.matrix(mat_median[, -c(1, ncol(mat_median)), drop = FALSE]), na.rm = TRUE)
    )
    if ( sum(as.logical(which(t1$abundance <= min_excl_last))) == 0 ){
      max7$classification[is.na(max7$classification)] <- "Expressed in all tissues"
    }
    
    c34 <- max7 %>% tidyr::drop_na(classification)
    
    # Class 5 (tissue enhanced) and Class 6 (mixed)
    max7_na <- max7 %>% dplyr::filter(is.na(classification))
    if (nrow(max7_na) >= 1){
      thr <- fct * mean(t1$abundance, na.rm = TRUE)
      max7_na <- max7_na %>%
        dplyr::mutate(
          classification = dplyr::if_else(abundance >= thr, "tissue enhanced", "mixed")
        )
    }
    
    # Combine all assignments and rename the classification column to the current feature name
    com_max <- dplyr::bind_rows(c12, c34, max7_na)
    
    j <<- j + 1
    colnames(com_max)[3] <- colnames(mat_median)[j]
    
    if (j %% 1000 == 0){
      cat(paste("Processed", j - 1, "/", ncol(mat_median) - 1, sep = " "), "\n")
    }
    
    # Return two-column frame: tissue_type + classification-for-this-feature
    com_max[, c(1, 3)]
  })
  
  # Informative message (kept from original)
  print("Finished classification and start merging classification matrix")
  
  # Merge all per-feature classification columns by tissue_type using full joins
  clas.all <- purrr::reduce_right(clas, dplyr::full_join, by = "tissue_type")
  
  return(clas.all)
}


ulhen_class_v3 <- function(mat_median, mat_missing, fct = 3) {
  # version 20251030 - based on missing rate (>0.5) for "not detected"
  j <- 1
  
  ## 检查输入一致性
  if (!all(dim(mat_median) == dim(mat_missing))) {
    stop("mat_median and mat_missing must have the same dimensions.")
  }
  if (!all(mat_median[[1]] == mat_missing[[1]])) {
    stop("The first column (tissue_type) of both matrices must match.")
  }
  
  ## 分类
  clas <- plyr::llply(seq_len(ncol(mat_median) - 1L), function(.k) {
    message('Analysis of protein No.', .k)
    num <- mat_median[[.k + 1L]]
    miss_rate <- mat_missing[[.k + 1L]]
    
    ## build t1
    t1 <- data.frame(
      tissue_type = mat_median[[1L]],
      abundance   = num,
      missing_rate = miss_rate,
      classification = NA_character_,
      stringsAsFactors = FALSE
    )
    
    ## sort by abundance
    o  <- order(t1$abundance, decreasing = TRUE, na.last = TRUE)
    t1 <- t1[o, , drop = FALSE]
    ab <- t1$abundance
    n  <- length(ab)
    
    ## -------- class 1 (not detected): based on missing rate > 0.5 --------
    t1$classification[t1$missing_rate > 0.5] <- "not detected"
    
    ## -------- class 2 (tissue enriched) --------
    if (!anyNA(ab)) {
      max_others <- rep.int(ab[1L], n)
      if (n >= 2L) max_others[1L] <- ab[2L]
      t1$classification[ab >= (fct * max_others)] <- "tissue enriched"
    }
    
    c12 <- t1[!is.na(t1$classification), , drop = FALSE]
    
    ## -------- class 3 (group enriched) --------
    max7 <- t1[is.na(t1$classification), , drop = FALSE]
    if (nrow(max7) > 7L) {
      cs    <- cumsum(max7$abundance)
      upto  <- seq_len(min(6L, nrow(max7) - 2L))
      if (length(upto)) {
        cond <- (cs[upto + 1L] / (upto + 1L)) >= (fct * max7$abundance[upto + 2L])
        if (any(cond)) {
          kstar <- max(which(cond)) + 1L
          max7$classification[seq_len(kstar)] <- "group enriched"
        }
      }
    } else if (nrow(max7) >= 2L && nrow(max7) <= 7L) {
      max7$classification <- "group enriched"
    }
    
    ## -------- class 4 (expressed in all tissues) --------
    # 若该蛋白在所有组织中缺失率均 <= 0.5，则认为“expressed in all tissues”
    if (!any(t1$missing_rate > 0.5, na.rm = TRUE)) {
      max7$classification[is.na(max7$classification)] <- "Expressed in all tissues"
    }
    c34 <- stats::na.omit(max7)
    
    ## -------- class 5 (tissue enhanced) + class 6 (mixed) --------
    max7_na <- max7[is.na(max7$classification), , drop = FALSE]
    if (nrow(max7_na) >= 1L) {
      thr <- fct * mean(t1$abundance, na.rm = TRUE)
      max7_na$classification <- ifelse(
        max7_na$abundance >= thr, "tissue enhanced", "mixed"
      )
    }
    
    com_max <- do.call(rbind, list(c12, c34, max7_na))
    
    j <<- j + 1L
    colnames(com_max)[4L] <- colnames(mat_median)[j]
    com_max[, c("tissue_type", colnames(com_max)[4L]), drop = FALSE]
  })
  
  print("Finished classification and start merging classification matrix")
  
  clas.all <- clas %>% reduce_right(full_join, by = "tissue_type")
  return(clas.all)
}

plot_CA_boxviolin <- function(df, cancer_abbr, na_cutoff = 0.5, refer = 'local', anno_txt = NULL, drawviolin = F){
  df %<>%
    add_column(label = stringr::str_c(df$cancer_abbr, ' (', df$sample_type, ')'), .before = 1) %>%
    dplyr::mutate(sample_type = factor(sample_type, levels = c('NT', 'T'))) %>%
    arrange(cancer_abbr, sample_type) %>%
    dplyr::mutate(sample_type = as.character(sample_type))
  
  #remove cancer_abbr with all subtypes > na_cutoff NA ratio
  na_ratio <- function(x) sum(is.na(x)) / length(x)
  cancer_abbr_selected <- df %>% dplyr::group_by(cancer_abbr, sample_type) %>%
    summarise_at(vars(2), list(na_ratio = 'na_ratio')) %>%
    ungroup() %>%
    filter(na_ratio <= na_cutoff) %>%
    pull(cancer_abbr) %>%
    unique()
  
  df %<>% filter(cancer_abbr %in% cancer_abbr_selected)
  
  rlt <- list()
  my_uni <- colnames(df) %>% tail(1)
  set.seed(2023)
  if(!exists('na_filling')){
    na_filling <- min(df[, ncol(df)], na.rm = T)
  }
  df[is.na(df %>% pull(all_of(my_uni))), my_uni] <- na_filling+log10(rnorm(sum(is.na(df %>% select(all_of(my_uni)))),mean=1,sd=0.05))
  df_fillrate <- df
  if(refer == 'uniprot'){
    my_title <- my_uniprot2prot(my_uni)
  } else if(refer == 'local'){
    my_title <- tryCatch(clusterProfiler::bitr(my_uni, fromType = "UNIPROT",toType = "SYMBOL",OrgDb = "org.Hs.eg.db",drop = T) %>% unite(label) %>% pull(),
                         error = function(e) my_uni)
  }
  
  if(!is.null(cancer_abbr)){
    my_title <- str_c(my_title, cancer_abbr, sep = '_')
  }
  colnames(df_fillrate)[ncol(df_fillrate)] <- colnames(df)[ncol(df)] <- c('intensity')
  
  
  # df <- na.omit(df)
  
  
  library(ggplot2)
  df_lines_wide <- df[df$cancer_abbr == cancer_abbr, ] %>%
    pivot_wider(id_cols = c('patient_ID'), names_from = 'sample_type', values_from = 'intensity',
                values_fn = mean) %>%
    mutate(Type = ifelse(`NT` < `T`, 'Up', 'Down'))
  
  df_lines <- df_lines_wide %>%
    pivot_longer(c('NT', 'T'), names_to = 'sample_type', values_to = 'intensity') %>%
    left_join(df %>% select(-sample_type, -intensity, -label) %>% distinct(), by = 'patient_ID') %>%
    mutate(label = str_glue("{cancer_abbr} ({sample_type})"))
  
  
  p <- ggplot(df) +
    aes(x = label, y = intensity) +
    # geom_violin(width = 1, color = '#888888') +
    geom_boxplot(outlier.shape = NA, width = 0.3, color = '#888888') +
    stat_summary(fun = mean, geom = 'point', size = 2, color = '#666666')+ # #BF2626
    geom_jitter(aes(color = sample_type), size=1.5, alpha=0.1) +
    geom_line(data = df_lines %>% filter(Type == 'Up'),
              aes(x = label, y = intensity, group = patient_ID), color = '#CF3F54', alpha = 0.3) +
    geom_line(data = df_lines %>% filter(Type == 'Down'),
              aes(x = label, y = intensity, group = patient_ID), color = '#99CFD1', alpha = 0.3) +
    scale_fill_manual(values = c(T = '#CF3F54', NT = '#99CFD1', N = '#03A2B3')) +
    scale_color_manual(values = c(T = '#CF3F54', NT = '#99CFD1', N = '#03A2B3')) +
    geom_vline(xintercept = str_which(unique(df$label), '\\(T\\)') + 0.5, color = '#000000', linetype = 'dashed') +
    annotate('text', x = seq(1.5, by = 2, length.out = length(unique(df$cancer_abbr))), y = max(df$intensity, na.rm = T) * 1.1, label = unique(df$cancer_abbr), size = 4, hjust = 0.5, vjust = 1) +
    labs(
      y = 'log10 Intensity', title = my_title
    )+
    theme(panel.grid = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 12,color = 'black'),
          axis.line = element_line(color = 'black'),
          axis.line.x = element_blank(),
          plot.subtitle = element_text(size = 30, hjust = 0, color = 'black')
    )+
    theme(legend.text = element_text(size = 12, color = 'black'), legend.position = 'none',
          legend.title = element_text(size = 15, color = 'black'))+
    theme(axis.title.x = element_blank(),
          #axis.text.x = element_text(size = 12, vjust = 0.5 ,color = 'black', angle = 90),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
    )+
    theme(axis.title.y = element_text(size = 12, hjust = 0.5, color = 'black', angle = 90),
          axis.text.y = element_text(size = 10,color = 'black', angle = 0),
          axis.ticks.y = element_blank()
    )
  
  if(!is.null(anno_txt)){
    # p <- p + annotate('text', x = (4 * which(unique(df$cancer_abbr) == cancer_abbr) - 1) / 2, y = max(df$intensity, na.rm = T) * 1.03, label = anno_txt, size = 3, hjust = 0.5, vjust = 0.5)
    p <- p + ggpubr::geom_signif(comparisons = list(unique(df$label[df$cancer_abbr == cancer_abbr])),
                                 annotations = anno_txt, size = 0.5, textsize = 2.5,
                                 y_position = max(df$intensity[df$cancer_abbr == cancer_abbr], na.rm = T) * 1,# tip_length = 0.04,
                                 hjust = 0.5, vjust = 0)
  }
  if(drawviolin){
    p <- p + geom_violin(width = 1, color = '#888888')
  }
  rlt$p <- p
  return(rlt)
}

plot_CA_boxviolin_multi <- function(df, cancer_abbrs = NULL, na_cutoff = 0.5, 
                                    refer = 'local', anno_txt = NULL, 
                                    drawviolin = FALSE, protein_id = NULL){
  
  # If cancer_abbrs not specified, use all available
  if(is.null(cancer_abbrs)){
    cancer_abbrs <- unique(df$cancer_abbr)
  }
  
  # Filter for specified cancer types
  df <- df %>% filter(cancer_abbr %in% cancer_abbrs)
  
  df %<>%
    add_column(label = stringr::str_c(df$cancer_abbr, ' (', df$sample_type, ')'), 
               .before = 1) %>%
    dplyr::mutate(sample_type = factor(sample_type, levels = c('NT', 'T'))) %>%
    arrange(cancer_abbr, sample_type) %>%
    dplyr::mutate(sample_type = as.character(sample_type))
  
  # Remove cancer_abbr with all subtypes > na_cutoff NA ratio
  na_ratio <- function(x) sum(is.na(x)) / length(x)
  cancer_abbr_selected <- df %>% 
    dplyr::group_by(cancer_abbr, sample_type) %>%
    summarise_at(vars(2), list(na_ratio = 'na_ratio')) %>%
    ungroup() %>%
    filter(na_ratio <= na_cutoff) %>%
    pull(cancer_abbr) %>%
    unique()
  
  df %<>% filter(cancer_abbr %in% cancer_abbr_selected)
  
  rlt <- list()
  my_uni <- colnames(df) %>% tail(1)
  
  # Use provided protein_id if available
  if(!is.null(protein_id)){
    my_uni <- protein_id
  }
  
  set.seed(2023)
  if(!exists('na_filling')){
    na_filling <- min(df[, ncol(df)], na.rm = T)
  }
  df[is.na(df %>% pull(all_of(my_uni))), my_uni] <- 
    na_filling + log10(rnorm(sum(is.na(df %>% select(all_of(my_uni)))), 
                             mean=1, sd=0.05))
  
  df_fillrate <- df
  
  # Get protein name for title
  if(refer == 'uniprot'){
    my_title <- my_uniprot2prot(my_uni)
  } else if(refer == 'local'){
    my_title <- tryCatch(
      clusterProfiler::bitr(my_uni, fromType = "UNIPROT", toType = "SYMBOL",
                            OrgDb = "org.Hs.eg.db", drop = T) %>% 
        unite(label) %>% pull(),
      error = function(e) my_uni
    )
  }
  
  colnames(df_fillrate)[ncol(df_fillrate)] <- colnames(df)[ncol(df)] <- c('intensity')
  
  # Process lines for all cancer types
  df_lines_list <- list()
  for(ca in cancer_abbr_selected){
    df_ca <- df[df$cancer_abbr == ca, ]
    
    df_lines_wide <- df_ca %>%
      pivot_wider(id_cols = c('patient_ID'), names_from = 'sample_type', 
                  values_from = 'intensity', values_fn = mean) %>%
      mutate(Type = ifelse(`NT` < `T`, 'Up', 'Down'))
    
    df_lines <- df_lines_wide %>%
      pivot_longer(c('NT', 'T'), names_to = 'sample_type', values_to = 'intensity') %>%
      left_join(df_ca %>% select(-sample_type, -intensity, -label) %>% distinct(), 
                by = 'patient_ID') %>%
      mutate(label = str_glue("{cancer_abbr} ({sample_type})"))
    
    df_lines_list[[ca]] <- df_lines
  }
  
  df_lines_all <- bind_rows(df_lines_list)
  
  # Create plot
  p <- ggplot(df) +
    aes(x = label, y = intensity) +
    geom_boxplot(outlier.shape = NA, width = 0.3, color = '#888888') +
    stat_summary(fun = mean, geom = 'point', size = 2, color = '#666666') +
    geom_jitter(aes(color = sample_type), size=1.5, alpha=0.1) +
    geom_line(data = df_lines_all %>% filter(Type == 'Up'),
              aes(x = label, y = intensity, group = patient_ID), 
              color = '#CF3F54', alpha = 0.3) +
    geom_line(data = df_lines_all %>% filter(Type == 'Down'),
              aes(x = label, y = intensity, group = patient_ID), 
              color = '#99CFD1', alpha = 0.3) +
    scale_fill_manual(values = c(T = '#CF3F54', NT = '#99CFD1', N = '#03A2B3')) +
    scale_color_manual(values = c(T = '#CF3F54', NT = '#99CFD1', N = '#03A2B3')) +
    geom_vline(xintercept = str_which(unique(df$label), '\\(T\\)') + 0.5, 
               color = '#000000', linetype = 'dashed') +
    annotate('text', x = seq(1.5, by = 2, length.out = length(cancer_abbr_selected)), 
             y = max(df$intensity, na.rm = T) * 1.2, 
             label = cancer_abbr_selected, size = 4, hjust = 0.5, vjust = 1) +
    labs(y = 'log10 Intensity', title = my_title) +
    theme(panel.grid = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 12, color = 'black'),
          axis.line = element_line(color = 'black'),
          axis.line.x = element_blank(),
          plot.subtitle = element_text(size = 30, hjust = 0, color = 'black')) +
    theme(legend.text = element_text(size = 12, color = 'black'), 
          legend.position = 'none',
          legend.title = element_text(size = 15, color = 'black')) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    theme(axis.title.y = element_text(size = 12, hjust = 0.5, color = 'black', angle = 90),
          axis.text.y = element_text(size = 10, color = 'black', angle = 0),
          axis.ticks.y = element_blank())
  
  # Add annotations if provided (modified for multiple cancer types)
  if(!is.null(anno_txt)){
    if(is.list(anno_txt)){
      for(i in seq_along(cancer_abbr_selected)){
        if(i <= length(anno_txt)){
          ca <- cancer_abbr_selected[i]
          p <- p + ggpubr::geom_signif(
            comparisons = list(unique(df$label[df$cancer_abbr == ca])),
            annotations = anno_txt[[i]], size = 0.5, textsize = 2.5,
            # y_position = max(df$intensity[df$cancer_abbr == ca], na.rm = T) * 1.1,
            y_position = max(df$intensity, na.rm = T) * 1,
            hjust = 0.5, vjust = 0
          )
        }
      }
    }
  }
  
  if(drawviolin){
    p <- p + geom_violin(width = 1, color = '#888888')
  }
  
  rlt$p <- p
  return(rlt)
}

plot_CA_boxviolin_new <- function(df, cancer_abbr, pepseq = NULL, na_cutoff = 0.5, refer = 'local', anno_txt = NULL, drawviolin = F, drawline = F){
  # for paired boxplot
  df %<>%
    add_column(label = stringr::str_c(df$cancer_abbr, ' (', df$sample_type, ')'), .before = 1) %>%
    dplyr::mutate(sample_type = factor(sample_type, levels = c('NT', 'T'))) %>%
    arrange(cancer_abbr, sample_type) %>%
    dplyr::mutate(sample_type = as.character(sample_type))
  
  #remove cancer_abbr with all subtypes > na_cutoff NA ratio
  na_ratio <- function(x) sum(is.na(x)) / length(x)
  cancer_abbr_selected <- df %>% dplyr::group_by(cancer_abbr, sample_type) %>%
    summarise_at(vars(2), list(na_ratio = 'na_ratio')) %>%
    ungroup() %>%
    filter(na_ratio <= na_cutoff) %>%
    pull(cancer_abbr) %>%
    unique()
  
  df %<>% filter(cancer_abbr %in% cancer_abbr_selected)
  
  rlt <- list()
  my_uni <- colnames(df) %>% tail(1)
  set.seed(2023)
  if(!exists('na_filling')){
    na_filling <- min(df[, ncol(df)], na.rm = T)
  }
  df[is.na(df %>% pull(all_of(my_uni))), my_uni] <- na_filling+log10(rnorm(sum(is.na(df %>% select(all_of(my_uni)))),mean=1,sd=0.05))
  df_fillrate <- df
  if(refer == 'uniprot'){
    my_title <- my_uniprot2prot(my_uni)
  } else if(refer == 'local'){
    my_title <- clusterProfiler::bitr(my_uni, fromType = "UNIPROT",toType = "SYMBOL",OrgDb = "org.Hs.eg.db",drop = T) %>% unite(label) %>% pull()
  }
  if(!is.null(pepseq)){
    my_title <- str_c(my_title, pepseq, sep = '_')
  }
  # 
  # if(!is.null(cancer_abbr)){
  #   my_title <- str_c(my_title, cancer_abbr, sep = '_')
  # }
  colnames(df_fillrate)[ncol(df_fillrate)] <- colnames(df)[ncol(df)] <- c('intensity')
  
  
  # df <- na.omit(df)
  
  
  library(ggplot2)
  df_lines_wide <- df[df$cancer_abbr == cancer_abbr, ] %>%
    pivot_wider(id_cols = c('patient_ID'), names_from = 'sample_type', values_from = 'intensity',
                values_fn = mean) %>%
    mutate(Type = ifelse(`NT` < `T`, 'Up', 'Down'))
  
  df_lines <- df_lines_wide %>%
    pivot_longer(c('NT', 'T'), names_to = 'sample_type', values_to = 'intensity') %>%
    left_join(df %>% select(-sample_type, -intensity, -label) %>% distinct(), by = 'patient_ID') %>%
    mutate(label = str_glue("{cancer_abbr} ({sample_type})"))
  
  
  p <- ggplot(df) +
    aes(x = label, y = intensity) +
    # geom_violin(width = 1, color = '#888888') +
    # stat_summary(fun = mean, geom = 'point', size = 2, color = '#000000')+ # #BF2626
    # geom_jitter(aes(color = sample_type), size=1.5, alpha=0.1) +
    geom_jitter(color = '#000000', size=0.3, alpha=0.5) +
    # scale_fill_manual(values = c(T = '#CF3F54', NT = '#99CFD1', N = '#03A2B3')) +
    # scale_color_manual(values = c(T = '#CF3F54', NT = '#99CFD1', N = '#03A2B3')) +
    geom_vline(xintercept = str_which(unique(df$label), '\\(T\\)') + 0.5, color = '#000000', linetype = 'dashed') +
    annotate('text', x = seq(1.5, by = 2, length.out = length(unique(df$cancer_abbr))), y = max(df$intensity, na.rm = T) * 1.1, label = unique(df$cancer_abbr), size = 4, hjust = 0.5, vjust = 1) +
    labs(
      y = 'log10 Intensity', title = my_title
    )+
    theme(panel.grid = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 12,color = 'black'),
          axis.line = element_line(color = 'black'),
          axis.line.x = element_blank(),
          plot.subtitle = element_text(size = 30, hjust = 0, color = 'black')
    )+
    theme(legend.text = element_text(size = 12, color = 'black'), legend.position = 'none',
          legend.title = element_text(size = 15, color = 'black'))+
    theme(axis.title.x = element_blank(),
          #axis.text.x = element_text(size = 12, vjust = 0.5 ,color = 'black', angle = 90),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
    )+
    theme(axis.title.y = element_text(size = 12, hjust = 0.5, color = 'black', angle = 90),
          axis.text.y = element_text(size = 10,color = 'black', angle = 0),
          axis.ticks.y = element_blank()
    )
  
  if(!is.null(anno_txt)){
    # p <- p + annotate('text', x = (4 * which(unique(df$cancer_abbr) == cancer_abbr) - 1) / 2, y = max(df$intensity, na.rm = T) * 1.03, label = anno_txt, size = 3, hjust = 0.5, vjust = 0.5)
    p <- p + ggpubr::geom_signif(comparisons = list(unique(df$label[df$cancer_abbr == cancer_abbr])),
                                 annotations = anno_txt, size = 0.5, textsize = 2.5,
                                 y_position = max(df$intensity[df$cancer_abbr == cancer_abbr], na.rm = T) * 1,# tip_length = 0.04,
                                 hjust = 0.5, vjust = 0)
  }
  if(drawviolin){
    p <- p + geom_violin(width = 1, color = '#888888')
  }
  if(drawline){
    p <- p +
      geom_line(data = df_lines %>% filter(Type == 'Up'),
                aes(x = label, y = intensity, group = patient_ID), color = '#E41A1C', alpha = 0.2) +
      geom_line(data = df_lines %>% filter(Type == 'Down'),
                aes(x = label, y = intensity, group = patient_ID), color = '#377EB8', alpha = 0.2)
  }
  p <- p + geom_boxplot(aes(color = sample_type), outlier.shape = NA, width = 0.3, fill = 'white') +
    scale_color_manual(values = c(T = '#E41A1C', NT = '#377EB8'))
  rlt$p <- p
  return(rlt)
}


# NAguideR methods -------------------
nafunctions<-function(x,method="zero"){
  # cite: Shisheng Wang, Wenxue Li, Liqiang Hu, Jingqiu Cheng, Hao Yang, Yansheng Liu, NAguideR: performing and prioritizing missing value imputations for consistent bottom-up proteomic analyses, Nucleic Acids Research, gkaa498, https://doi.org/10.1093/nar/gkaa498.
  # Code from GitHub: Commit 15ec862 wangshisheng authored on Aug 20, 2021
  df<-df1<-as.data.frame(x)
  method<-tolower(method)
  if(method=="zero"){
    df[is.na(df)]<-0
  }
  else if(method=="minimum"){
    df[is.na(df)]<-min(df1,na.rm = TRUE)
  }
  else if(method=="colmedian"){
    library(e1071)
    df<-impute(df1,what ="median")
  }
  else if(method=="rowmedian"){
    library(e1071)
    dfx<-impute(t(df1),what ="median")
    df<-t(dfx)
  }
  else if(method=="knnmethod"){
    library(impute)
    data_zero1<-impute.knn(as.matrix(df1),k = 10, rowmax = 1, colmax = 1)#rowmax = 0.9, colmax = 0.9
    df<-data_zero1$data
  }
  else if(method=="seqknn"){
    library(SeqKnn)
    df <- SeqKNN(df1,k = 10)
  }
  else if(method=="bpca"){
    library(pcaMethods)
    data_zero1<-pcaMethods::pca(as.matrix(df1), nPcs = ncol(df1)-1, method = "bpca", maxSteps =100)
    df<-completeObs(data_zero1)
  }
  else if(method=="svdmethod"){
    library(pcaMethods)
    data_zero1<-pcaMethods::pca(as.matrix(df1), nPcs = ncol(df1)-1, method = "svdImpute")
    df<-completeObs(data_zero1)
  }
  else if(method=="lls"){
    library(pcaMethods)
    data_zero1<-llsImpute(t(df1), k = 10)
    df<-t(completeObs(data_zero1))
  }
  else if(method=="mle"){
    library(norm)
    xxm<-as.matrix(df1)
    ss <- norm::prelim.norm(xxm)
    thx <- norm::em.norm(ss)
    norm::rngseed(123)
    df <- norm::imp.norm(ss, thx, xxm)
  }
  else if(method=="qrilc"){
    library(imputeLCMD)
    xxm<-t(df1)
    data_zero1 <- imputeLCMD::impute.QRILC(xxm, tune.sigma = 1)[[1]]
    df<-t(data_zero1)
  }
  else if(method=="mindet"){
    library(imputeLCMD)
    xxm<-as.matrix(df1)
    df <- imputeLCMD::impute.MinDet(xxm, q = 0.01)
  }
  else if(method=="minprob"){
    library(imputeLCMD)
    xxm<-as.matrix(df1)
    df <- imputeLCMD::impute.MinProb(xxm, q = 0.01, tune.sigma = 1)
  }
  else if(method=="irm"){
    library(VIM)
    df <- irmi(df1, trace = TRUE,imp_var=FALSE)
    rownames(df)<-rownames(df1)
  }
  else if(method=="impseq"){
    library(rrcovNA)
    df <- impSeq(df1)
  }
  else if(method=="impseqrob"){
    library(rrcovNA)
    data_zero1 <- impSeqRob(df1, alpha=0.9)
    df<-data_zero1$x
  }
  else if(method=="mice-norm"){
    library(mice)
    minum<-5
    datareadmi<-mice(df1,m=minum,seed = 1234, method ="norm")
    newdatareadmi<-0
    for (i in 1:minum) {
      newdatareadmi<-complete(datareadmi,action = i)+newdatareadmi
    }
    df<-newdatareadmi/minum
    rownames(df)<-rownames(df1)
  }
  else if(method=="mice-cart"){
    library(mice)
    minum<-5
    datareadmi<-mice(df1,m=minum,seed = 1234, method ="cart")
    newdatareadmi<-0
    for (i in 1:minum) {
      newdatareadmi<-complete(datareadmi,action = i)+newdatareadmi
    }
    df<-newdatareadmi/minum
    rownames(df)<-rownames(df1)
  }
  else if(method=="trknn"){
    # source('Trunc_KNN/Imput_funcs.r') # only call the minimal requirements
    {
      ##################################################################################
      #### MLE for the Truncated Normal
      #### Creating a Function that Returns the Log Likelihood, Gradient and
      #### Hessian Functions
      ##################################################################################
      
      ## data = numeric vector
      ## t    = truncation limits
      mklhood <- function(data, t, ...) {
        
        data <- na.omit(data)
        n <- length(data)
        t <- sort(t)
        
        psi<-function(y, mu, sigma){
          exp(-(y-mu)^2/(2*sigma^2))/(sigma*sqrt(2*pi))
        }
        
        psi.mu<-function(y,mu,sigma){
          exp(-(y-mu)^2/(2*sigma^2)) * ((y-mu)/(sigma^3*sqrt(2*pi)))
        }
        
        psi.sigma<-function(y,mu,sigma){
          exp(-(y-mu)^2/(2*sigma^2)) *
            (((y-mu)^2)/(sigma^4*sqrt(2*pi)) - 1/(sigma^2*sqrt(2*pi)))
        }
        
        psi2.mu<-function(y,mu,sigma){
          exp(-(y - mu)^2/(2*sigma^2)) *
            (((y - mu)^2)/(sigma^5*sqrt(2*pi))-1/(sigma^3*sqrt(2*pi)))
        }
        
        psi2.sigma<-function(y,mu,sigma){
          exp(-(y-mu)^2/(2*sigma^2)) *
            ((2)/(sigma^3*sqrt(2*pi)) - (5*(y-mu))/(sigma^5*sqrt(2*pi)) +
               ((y-mu)^4)/(sigma^7*sqrt(2*pi)))
        }
        
        psi12.musig<-function(y,mu,sigma){
          exp(-(y-mu)^2/(2*sigma^2)) *
            (((y-mu)^3)/(sigma^6*sqrt(2*pi)) - (3*(y-mu))/(sigma^4*sqrt(2*pi)))
        }
        
        ll.tnorm2<-function(p){
          out <- (-n*log(pnorm(t[2],p[1],p[2])-pnorm(t[1],p[1],p[2]))) -
            (n*log(sqrt(2*pi*p[2]^2))) - (sum((data-p[1])^2)/(2*p[2]^2))
          -1*out
        }
        
        grad.tnorm<-function(p){
          g1 <- (-n*(integrate(psi.mu,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value) /
                   (pnorm(max(t),p[1],p[2])-pnorm(min(t),p[1],p[2]))) - ((n*p[1]-sum(data))/p[2]^2)
          g2 <- (-n*(integrate(psi.sigma,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value) /
                   (pnorm(max(t),p[1],p[2])-pnorm(min(t),p[1],p[2]))) - ((n)/(p[2])) + ((sum((data-p[1])^2))/(p[2]^3))
          out <- c(g1,g2)
          return(out)
        }
        
        hessian.tnorm<-function(p){
          
          h1<- -n*(integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value *
                     integrate(psi2.mu,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value -
                     integrate(psi.mu,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value^2) /
            (integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value^2) -
            n/(p[2]^2)
          
          h3<- -n*(integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value *
                     integrate(psi12.musig,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value -
                     integrate(psi.mu,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value *
                     integrate(psi.sigma,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value) /
            (integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value^2) +
            (2*(n*p[1]-sum(data)))/(p[2]^3)
          
          h2<- -n*(integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value *
                     integrate(psi2.sigma,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value -
                     integrate(psi.sigma,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value^2) /
            (integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value^2) +
            (n)/(p[2]^2)-(3*sum((data-p[1])^2))/(p[2]^4)
          
          H<-matrix(0,nrow=2,ncol=2)
          H[1,1]<-h1
          H[2,2]<-h2
          H[1,2]<-H[2,1]<-h3
          return(H)
        }
        
        
        return(list(ll.tnorm2 = ll.tnorm2, grad.tnorm = grad.tnorm, hessian.tnorm = hessian.tnorm))
      }
      ##################################################################################
      ###### Newton Raphson Function
      ###### This takes in the Objects Returned from mklhood Function above
      ##################################################################################
      
      NewtonRaphsonLike <- function(lhood, p, tol = 1e-07, maxit = 100) {
        
        cscore <- lhood$grad.tnorm(p)
        if(sum(abs(cscore)) < tol)
          return(list(estimate = p, value = lhood$ll.tnorm2(p), iter = 0))
        cur <- p
        for(i in 1:maxit) {
          inverseHess <- solve(lhood$hessian.tnorm(cur))
          cscore <- lhood$grad.tnorm(cur)
          new <- cur - cscore %*% inverseHess
          if (new[2] <= 0) stop("Sigma < 0")
          cscore <- lhood$grad.tnorm(new)
          
          if(((abs(lhood$ll.tnorm2(cur)- lhood$ll.tnorm2(new))/(lhood$ll.tnorm2(cur))) < tol))
            return(list(estimate = new, value= lhood$ll.tnorm2(new), iter = i))
          cur <- new
        }
        
        return(list(estimate = new, value= lhood$ll.tnorm2(new), iter = i))
      }
      
      ##################################################################################
      ###### Based on the MLE Functions (mklhood) and NewtonRaphson Function
      ###### (NewtonRaphsonLike), This function estimates the MEAN and SD from the
      ###### Truncated using Newton Raphson.
      ##################################################################################
      
      ## missingdata = matrix where rows = features, columns = samples
      ## perc = if %MVs > perc then just sample mean / SD
      ## iter = # iterations in NR algorithm
      
      EstimatesComputation <- function(missingdata, perc, iter=50) {
        
        ## 2 column matrix where column 1 = means, column 2 = SD
        ParamEstim <- matrix(NA, nrow = nrow(missingdata), ncol = 2)
        nsamp <- ncol(missingdata)
        
        ## sample means / SDs
        ParamEstim[,1] <- rowMeans(missingdata, na.rm = TRUE)
        ParamEstim[,2] <- apply(missingdata, 1, function(x) sd(x, na.rm = TRUE))
        
        ## Case 1: missing % > perc => use sample mean / SD
        na.sum <- apply(missingdata, 1, function(x) sum(is.na(x)))
        idx1 <- which(na.sum/nsamp >= perc)
        
        ## Case 2: sample mean > 3 SD away from LOD => use sample mean / SD
        lod <- min(missingdata, na.rm=TRUE) ## why use the min of whole data set??????
        idx2 <- which(ParamEstim[,1] > 3*ParamEstim[,2] + lod)
        
        ## Case 3: for all others, use NR method to obtain truncated mean / SD estimate
        idx.nr <- setdiff(1:nrow(missingdata), c(idx1, idx2))
        ## t = limits of integration (LOD and upper)
        upplim <- max(missingdata, na.rm=TRUE) + 2*max(ParamEstim[,2])
        for (i in idx.nr) {
          Likelihood <- mklhood(missingdata[i,], t=c(lod, upplim))
          res <- tryCatch(NewtonRaphsonLike(Likelihood, p = ParamEstim[i,]),
                          error = function(e) 1000)
          
          if (length(res) == 1) {
            next
          } else if (res$iter >= iter) {
            next
          } else {
            ParamEstim[i,] <- as.numeric(res$estimate)
          }
        }
        return(ParamEstim)
      }
      
      
      
      ####################################################################
      #### This Function imputes the data BASED on KNN-EUCLIDEAN
      ####################################################################
      
      ## data = data set to be imputed, where rows = features, columns = samples
      ## k    = number of neighbors for imputing values
      ## rm.na, rm.nan, rm.inf = whether NA, NaN, and Inf values should be imputed
      
      
      KNNEuc <- function (data, k, rm.na = TRUE, rm.nan = TRUE, rm.inf = TRUE) {
        
        nr <- dim(data)[1]
        
        imp.knn <- data
        imp.knn[is.finite(data) == FALSE] <- NA
        t.data<-t(data)
        
        mv.ind <- which(is.na(imp.knn), arr.ind = TRUE)
        arrays <- unique(mv.ind[, 2])
        array.ind <- match(arrays, mv.ind[, 2])
        nfeatures <- 1:nr
        
        for (i in 1:length(arrays)) {
          set <- array.ind[i]:min((array.ind[(i + 1)] - 1), dim(mv.ind)[1], na.rm = TRUE)
          cand.features <- nfeatures[-unique(mv.ind[set, 1])]
          cand.vectors <- t.data[,cand.features]
          exp.num <- arrays[i]
          
          for (j in set) {
            feature.num <- mv.ind[j, 1]
            tar.vector <- data[feature.num,]
            
            dist <- sqrt(colMeans((tar.vector-cand.vectors)^2, na.rm = TRUE))
            dist[is.nan(dist) | is.na(dist)] <- Inf
            dist[dist==0] <- ifelse(is.finite(min(dist[dist>0])), min(dist[dist>0])/2, 1)
            
            if (sum(is.finite(dist)) < k) {
              stop(message = "Fewer than K finite distances found")
            }
            k.features.ind <- order(dist)[1:k]
            k.features <- cand.features[k.features.ind]
            wghts <- 1/dist[k.features.ind]/sum(1/dist[k.features.ind])
            imp.knn[feature.num, exp.num] <- wghts %*% data[k.features, exp.num]
          }
        }
        
        if (!rm.na) {
          imp.knn[is.na(data) == TRUE & is.nan(data) == FALSE] <- NA
        }
        if (!rm.inf) {
          index <- is.finite(data) == FALSE & is.na(data) == FALSE &
            is.nan(data) == FALSE
          imp.knn[index] <- data[index]
        }
        if (!rm.nan) {
          imp.knn[is.nan(data) == TRUE] <- NaN
        }
        return(imp.knn)
      }
      
      ####################################################################
      #### This Function imputes the data based on KNN-CORRELATION or
      #### KNN-TRUNCATION. The Parameter Estimates based on the Truncated
      #### Normal from EstimateComputation function is run on this function
      ####################################################################
      
      imputeKNN <- function (data, k , distance = "correlation",
                             rm.na = TRUE, rm.nan = TRUE, rm.inf = TRUE, perc=1,...) {
        
        if (!(is.matrix(data))) {
          stop(message = paste(deparse(substitute(data)),
                               " is not a matrix.", sep = ""))
        }
        
        distance <- match.arg(distance, c("correlation","truncation"))
        
        nr <- dim(data)[1]
        if (k < 1 | k > nr) {
          stop(message = "k should be between 1 and the number of rows")
        }
        
        if (distance=="correlation"){
          genemeans<-rowMeans(data,na.rm=TRUE)
          genesd<-apply(data, 1, function(x) sd(x, na.rm = TRUE))
          data<-(data-genemeans)/genesd
        }
        
        if (distance=="truncation"){
          
          ParamMat <- EstimatesComputation(data, perc = perc)
          
          genemeans<-ParamMat[,1]
          genesd<-ParamMat[,2]
          data<-(data-genemeans)/genesd
        }
        
        imp.knn <- data
        imp.knn[is.finite(data) == FALSE] <- NA
        t.data<-t(data)
        
        mv.ind <- which(is.na(imp.knn), arr.ind = TRUE)
        arrays <- unique(mv.ind[, 2])
        array.ind <- match(arrays, mv.ind[, 2])
        ngenes <- 1:nr
        
        for (i in 1:length(arrays)) {
          set <- array.ind[i]:min((array.ind[(i + 1)] - 1), dim(mv.ind)[1],
                                  na.rm = TRUE)
          cand.genes <- ngenes[-unique(mv.ind[set, 1])]
          cand.vectors <- t.data[,cand.genes]
          exp.num<- arrays[i]
          for (j in set) {
            
            gene.num <- mv.ind[j, 1]
            tar.vector <- data[gene.num,]
            
            r <- (cor(cand.vectors,tar.vector, use = "pairwise.complete.obs"))
            dist <- switch(distance,
                           correlation = (1 - abs(r)),
                           truncation = (1 - abs(r)))
            dist[is.nan(dist) | is.na(dist)] <- Inf
            dist[dist==0]<-ifelse(is.finite(min(dist[dist>0])), min(dist[dist>0])/2, 1)
            dist[abs(r) == 1] <- Inf
            
            if (sum(is.finite(dist)) < k) {
              stop(message = "Fewer than K finite distances found")
            }
            k.genes.ind <- order(dist)[1:k]
            k.genes <- cand.genes[k.genes.ind]
            
            wghts <- (1/dist[k.genes.ind]/sum(1/dist[k.genes.ind])) * sign(r[k.genes.ind])
            imp.knn[gene.num, exp.num] <- wghts %*% data[k.genes, exp.num]
          }
        }
        
        if (distance=="correlation") {
          imp.knn <- (imp.knn * genesd) + genemeans
        }
        
        if(distance=="truncation") {
          imp.knn <- (imp.knn * genesd) + genemeans
        }
        
        if (!rm.na) {
          imp.knn[is.na(data) == TRUE & is.nan(data) == FALSE] <- NA
        }
        if (!rm.inf) {
          index <- is.finite(data) == FALSE & is.na(data) == FALSE &
            is.nan(data) == FALSE
          imp.knn[index] <- data[index]
        }
        if (!rm.nan) {
          imp.knn[is.nan(data) == TRUE] <- NaN
        }
        return(imp.knn)
      }
    }
    sim_trKNN_wrapper <- function(data) {
      result <- data %>% as.matrix %>% t %>% imputeKNN(., k=10, distance='truncation', perc=0) %>% t
      return(result)
    }
    df1x <- sim_trKNN_wrapper(t(df1))
    df<-as.data.frame(t(df1x))
  }
  else if(method=="rf"){
    library(missForest)
    data_zero1 <- missForest(t(df1), maxiter =10,ntree = input$rfntrees,mtry=floor(nrow(df1)^(1/3)),verbose = TRUE)
    df<-t(data_zero1$ximp)
  }
  else if(method=="pi"){
    width <- input$piwidth
    downshift <- input$pidownshift
    for(i in 1:ncol(df1)){
      temp <- df1[[i]]
      if(sum(is.na(temp))>0){
        temp.sd <- width * sd(temp[!is.na(temp)], na.rm = TRUE)
        temp.mean <- mean(temp[!is.na(temp)], na.rm = TRUE) - downshift * sd(temp[!is.na(temp)], na.rm = TRUE)
        n.missing <- sum(is.na(temp))
        temp[is.na(temp)] <- rnorm(n.missing, mean = temp.mean, sd = temp.sd)
        df[[i]]<-temp
      }
    }
    df
  }
  else if(method=="grr"){
    library(DreamAI)
    df<-impute.RegImpute(data=as.matrix(df1), fillmethod = "row_mean", maxiter_RegImpute = 10,conv_nrmse = 1e-03)
  }
  else if(method=="gms"){
    library(GMSimpute)
    df<-GMS.Lasso(df1,nfolds=3,log.scale=FALSE,TS.Lasso=TRUE)
  }
  else{
    stop("Unspported methods so far~~")
  }
  df<-as.data.frame(df)
  df
}



#' Impute missing values in a data frame using various methods
#'
#' @param x A data frame or matrix containing missing values.
#' @param method Character string specifying the imputation method.
#'   One of: "zero", "minimum", "colmedian", "rowmedian", "knn",
#'   "seqknn", "bpca", "svd", "lls", "mle", "qrilc", "mindet", "minprob",
#'   "irm", "impseq", "impseqrob", "mice_norm", "mice_cart", "trknn",
#'   "rf", "pi", "grr", "gms".
#' @return A data frame of the same dimension as `x`, with NAs imputed.
#' @references
#' Shisheng Wang et al., *NAguideR: performing and prioritizing missing value imputations*,
#' Nucleic Acids Research, 2021; doi:10.1093/nar/gkaa498
#' @examples
#' df_filled <- nafunctions(my_data, method = "knn")
#' @export
naimpute.methods <- c("zero", "minimum", "colmedian", "rowmedian",
                      "knn", "seqknn", "bpca", "svd", "lls",
                      "mle", "qrilc", "mindet", "minprob",
                      "irm", "impseq", "impseqrob",
                      "mice_norm", "mice_cart", "trknn",
                      "rf", "pi", "grr", "gms")
naimpute <- function(x,
                     method = c("zero", "minimum", "colmedian", "rowmedian",
                                "knn", "seqknn", "bpca", "svd", "lls",
                                "mle", "qrilc", "mindet", "minprob",
                                "irm", "impseq", "impseqrob",
                                "mice_norm", "mice_cart", "trknn",
                                "rf", "pi", "grr", "gms")) {
  method <- match.arg(tolower(method), method)
  df_orig <- as.data.frame(x)
  df_work <- df_orig
  
  imputed <- switch(
    method,
    
    # ---- Simple Replacements ----
    zero = {
      df_work[is.na(df_work)] <- 0
      df_work
    },
    minimum = {
      df_work[is.na(df_work)] <- min(df_orig, na.rm = TRUE)
      df_work
    },
    
    # ---- Median Imputation ----
    colmedian = {
      if (!requireNamespace("e1071", quietly = TRUE)) {
        stop("Package 'e1071' required for colmedian", call. = FALSE)
      }
      e1071::impute(df_orig, what = "median")
    },
    rowmedian = {
      if (!requireNamespace("e1071", quietly = TRUE)) {
        stop("Package 'e1071' required for rowmedian", call. = FALSE)
      }
      t(e1071::impute(t(df_orig), what = "median"))
    },
    
    # ---- KNN-based Methods ----
    knn = {
      if (!requireNamespace("impute", quietly = TRUE)) {
        stop("Package 'impute' required for knn", call. = FALSE)
      }
      impute::impute.knn(as.matrix(df_orig), k = 10)$data
    },
    seqknn = {
      if (!requireNamespace("SeqKnn", quietly = TRUE)) {
        stop("Package 'SeqKnn' required for seqknn", call. = FALSE)
      }
      SeqKnn::SeqKNN(df_orig, k = 10)
    },
    
    # ---- PCA / SVD Methods ----
    bpca = {
      if (!requireNamespace("pcaMethods", quietly = TRUE)) {
        stop("Package 'pcaMethods' required for bpca", call. = FALSE)
      }
      pca <- pcaMethods::pca(as.matrix(df_orig),
                             nPcs = ncol(df_orig) - 1,
                             method = "bpca", maxSteps = 100)
      pcaMethods::completeObs(pca)
    },
    svd = {
      if (!requireNamespace("pcaMethods", quietly = TRUE)) {
        stop("Package 'pcaMethods' required for svd", call. = FALSE)
      }
      pca <- pcaMethods::pca(as.matrix(df_orig),
                             nPcs = ncol(df_orig) - 1,
                             method = "svdImpute")
      pcaMethods::completeObs(pca)
    },
    lls = {
      if (!requireNamespace("pcaMethods", quietly = TRUE)) {
        stop("Package 'pcaMethods' required for lls", call. = FALSE)
      }
      # Local least squares on transposed data
      pca <- pcaMethods::llsImpute(t(df_orig), k = 10)
      t(pcaMethods::completeObs(pca))
    },
    
    # ---- Model-based Methods ----
    mle = {
      if (!requireNamespace("norm", quietly = TRUE)) {
        stop("Package 'norm' required for mle", call. = FALSE)
      }
      mat <- as.matrix(df_orig)
      pre <- norm::prelim.norm(mat)
      em  <- norm::em.norm(pre)
      norm::rngseed(123)
      norm::imp.norm(pre, em, mat)
    },
    qrilc = {
      if (!requireNamespace("imputeLCMD", quietly = TRUE)) {
        stop("Package 'imputeLCMD' required for qrilc", call. = FALSE)
      }
      tmp   <- t(df_orig)
      filled <- imputeLCMD::impute.QRILC(tmp, tune.sigma = 1)[[1]]
      t(filled)
    },
    mindet = {
      if (!requireNamespace("imputeLCMD", quietly = TRUE)) {
        stop("Package 'imputeLCMD' required for mindet", call. = FALSE)
      }
      imputeLCMD::impute.MinDet(as.matrix(df_orig), q = 0.01)
    },
    minprob = {
      if (!requireNamespace("imputeLCMD", quietly = TRUE)) {
        stop("Package 'imputeLCMD' required for minprob", call. = FALSE)
      }
      imputeLCMD::impute.MinProb(as.matrix(df_orig), q = 0.01, tune.sigma = 1)
    },
    
    # ---- Iterative / Sequential Methods ----
    irm = {
      if (!requireNamespace("VIM", quietly = TRUE)) {
        stop("Package 'VIM' required for irm", call. = FALSE)
      }
      out <- VIM::irmi(df_orig, trace = TRUE, imp_var = FALSE)
      `rownames<-`(out, rownames(df_orig))
    },
    impseq = {
      if (!requireNamespace("rrcovNA", quietly = TRUE)) {
        stop("Package 'rrcovNA' required for impseq", call. = FALSE)
      }
      rrcovNA::impSeq(df_orig)
    },
    impseqrob = {
      if (!requireNamespace("rrcovNA", quietly = TRUE)) {
        stop("Package 'rrcovNA' required for impseqrob", call. = FALSE)
      }
      out <- rrcovNA::impSeqRob(df_orig, alpha = 0.9)
      out$x
    },
    
    # ---- MICE Methods ----
    mice_norm = ,
    mice_cart = {
      if (!requireNamespace("mice", quietly = TRUE)) {
        stop("Package 'mice' required for mice methods", call. = FALSE)
      }
      m <- 5
      meth <- if (method == "mice_norm") "norm" else "cart"
      imp <- mice::mice(df_orig, m = m, seed = 1234, method = meth)
      avg <- Reduce(`+`, lapply(1:m, function(i) mice::complete(imp, action = i)))
      avg / m
    },
    
    # ---- Specialty Methods ----
    trknn = {
      # Assumes that 'Trunc_KNN/Imput_funcs.r' defines imputeKNN()
      # source('Trunc_KNN/Imput_funcs.r') # only call the minimal requirements
      {
        ##################################################################################
        #### MLE for the Truncated Normal
        #### Creating a Function that Returns the Log Likelihood, Gradient and
        #### Hessian Functions
        ##################################################################################
        
        ## data = numeric vector
        ## t    = truncation limits
        mklhood <- function(data, t, ...) {
          
          data <- na.omit(data)
          n <- length(data)
          t <- sort(t)
          
          psi<-function(y, mu, sigma){
            exp(-(y-mu)^2/(2*sigma^2))/(sigma*sqrt(2*pi))
          }
          
          psi.mu<-function(y,mu,sigma){
            exp(-(y-mu)^2/(2*sigma^2)) * ((y-mu)/(sigma^3*sqrt(2*pi)))
          }
          
          psi.sigma<-function(y,mu,sigma){
            exp(-(y-mu)^2/(2*sigma^2)) *
              (((y-mu)^2)/(sigma^4*sqrt(2*pi)) - 1/(sigma^2*sqrt(2*pi)))
          }
          
          psi2.mu<-function(y,mu,sigma){
            exp(-(y - mu)^2/(2*sigma^2)) *
              (((y - mu)^2)/(sigma^5*sqrt(2*pi))-1/(sigma^3*sqrt(2*pi)))
          }
          
          psi2.sigma<-function(y,mu,sigma){
            exp(-(y-mu)^2/(2*sigma^2)) *
              ((2)/(sigma^3*sqrt(2*pi)) - (5*(y-mu))/(sigma^5*sqrt(2*pi)) +
                 ((y-mu)^4)/(sigma^7*sqrt(2*pi)))
          }
          
          psi12.musig<-function(y,mu,sigma){
            exp(-(y-mu)^2/(2*sigma^2)) *
              (((y-mu)^3)/(sigma^6*sqrt(2*pi)) - (3*(y-mu))/(sigma^4*sqrt(2*pi)))
          }
          
          ll.tnorm2<-function(p){
            out <- (-n*log(pnorm(t[2],p[1],p[2])-pnorm(t[1],p[1],p[2]))) -
              (n*log(sqrt(2*pi*p[2]^2))) - (sum((data-p[1])^2)/(2*p[2]^2))
            -1*out
          }
          
          grad.tnorm<-function(p){
            g1 <- (-n*(integrate(psi.mu,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value) /
                     (pnorm(max(t),p[1],p[2])-pnorm(min(t),p[1],p[2]))) - ((n*p[1]-sum(data))/p[2]^2)
            g2 <- (-n*(integrate(psi.sigma,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value) /
                     (pnorm(max(t),p[1],p[2])-pnorm(min(t),p[1],p[2]))) - ((n)/(p[2])) + ((sum((data-p[1])^2))/(p[2]^3))
            out <- c(g1,g2)
            return(out)
          }
          
          hessian.tnorm<-function(p){
            
            h1<- -n*(integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value *
                       integrate(psi2.mu,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value -
                       integrate(psi.mu,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value^2) /
              (integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value^2) -
              n/(p[2]^2)
            
            h3<- -n*(integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value *
                       integrate(psi12.musig,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value -
                       integrate(psi.mu,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value *
                       integrate(psi.sigma,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value) /
              (integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value^2) +
              (2*(n*p[1]-sum(data)))/(p[2]^3)
            
            h2<- -n*(integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value *
                       integrate(psi2.sigma,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value -
                       integrate(psi.sigma,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value^2) /
              (integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value^2) +
              (n)/(p[2]^2)-(3*sum((data-p[1])^2))/(p[2]^4)
            
            H<-matrix(0,nrow=2,ncol=2)
            H[1,1]<-h1
            H[2,2]<-h2
            H[1,2]<-H[2,1]<-h3
            return(H)
          }
          
          
          return(list(ll.tnorm2 = ll.tnorm2, grad.tnorm = grad.tnorm, hessian.tnorm = hessian.tnorm))
        }
        ##################################################################################
        ###### Newton Raphson Function
        ###### This takes in the Objects Returned from mklhood Function above
        ##################################################################################
        
        NewtonRaphsonLike <- function(lhood, p, tol = 1e-07, maxit = 100) {
          
          cscore <- lhood$grad.tnorm(p)
          if(sum(abs(cscore)) < tol)
            return(list(estimate = p, value = lhood$ll.tnorm2(p), iter = 0))
          cur <- p
          for(i in 1:maxit) {
            inverseHess <- solve(lhood$hessian.tnorm(cur))
            cscore <- lhood$grad.tnorm(cur)
            new <- cur - cscore %*% inverseHess
            if (new[2] <= 0) stop("Sigma < 0")
            cscore <- lhood$grad.tnorm(new)
            
            if(((abs(lhood$ll.tnorm2(cur)- lhood$ll.tnorm2(new))/(lhood$ll.tnorm2(cur))) < tol))
              return(list(estimate = new, value= lhood$ll.tnorm2(new), iter = i))
            cur <- new
          }
          
          return(list(estimate = new, value= lhood$ll.tnorm2(new), iter = i))
        }
        
        ##################################################################################
        ###### Based on the MLE Functions (mklhood) and NewtonRaphson Function
        ###### (NewtonRaphsonLike), This function estimates the MEAN and SD from the
        ###### Truncated using Newton Raphson.
        ##################################################################################
        
        ## missingdata = matrix where rows = features, columns = samples
        ## perc = if %MVs > perc then just sample mean / SD
        ## iter = # iterations in NR algorithm
        
        EstimatesComputation <- function(missingdata, perc, iter=50) {
          
          ## 2 column matrix where column 1 = means, column 2 = SD
          ParamEstim <- matrix(NA, nrow = nrow(missingdata), ncol = 2)
          nsamp <- ncol(missingdata)
          
          ## sample means / SDs
          ParamEstim[,1] <- rowMeans(missingdata, na.rm = TRUE)
          ParamEstim[,2] <- apply(missingdata, 1, function(x) sd(x, na.rm = TRUE))
          
          ## Case 1: missing % > perc => use sample mean / SD
          na.sum <- apply(missingdata, 1, function(x) sum(is.na(x)))
          idx1 <- which(na.sum/nsamp >= perc)
          
          ## Case 2: sample mean > 3 SD away from LOD => use sample mean / SD
          lod <- min(missingdata, na.rm=TRUE) ## why use the min of whole data set??????
          idx2 <- which(ParamEstim[,1] > 3*ParamEstim[,2] + lod)
          
          ## Case 3: for all others, use NR method to obtain truncated mean / SD estimate
          idx.nr <- setdiff(1:nrow(missingdata), c(idx1, idx2))
          ## t = limits of integration (LOD and upper)
          upplim <- max(missingdata, na.rm=TRUE) + 2*max(ParamEstim[,2])
          for (i in idx.nr) {
            Likelihood <- mklhood(missingdata[i,], t=c(lod, upplim))
            res <- tryCatch(NewtonRaphsonLike(Likelihood, p = ParamEstim[i,]),
                            error = function(e) 1000)
            
            if (length(res) == 1) {
              next
            } else if (res$iter >= iter) {
              next
            } else {
              ParamEstim[i,] <- as.numeric(res$estimate)
            }
          }
          return(ParamEstim)
        }
        
        
        
        ####################################################################
        #### This Function imputes the data BASED on KNN-EUCLIDEAN
        ####################################################################
        
        ## data = data set to be imputed, where rows = features, columns = samples
        ## k    = number of neighbors for imputing values
        ## rm.na, rm.nan, rm.inf = whether NA, NaN, and Inf values should be imputed
        
        
        KNNEuc <- function (data, k, rm.na = TRUE, rm.nan = TRUE, rm.inf = TRUE) {
          
          nr <- dim(data)[1]
          
          imp.knn <- data
          imp.knn[is.finite(data) == FALSE] <- NA
          t.data<-t(data)
          
          mv.ind <- which(is.na(imp.knn), arr.ind = TRUE)
          arrays <- unique(mv.ind[, 2])
          array.ind <- match(arrays, mv.ind[, 2])
          nfeatures <- 1:nr
          
          for (i in 1:length(arrays)) {
            set <- array.ind[i]:min((array.ind[(i + 1)] - 1), dim(mv.ind)[1], na.rm = TRUE)
            cand.features <- nfeatures[-unique(mv.ind[set, 1])]
            cand.vectors <- t.data[,cand.features]
            exp.num <- arrays[i]
            
            for (j in set) {
              feature.num <- mv.ind[j, 1]
              tar.vector <- data[feature.num,]
              
              dist <- sqrt(colMeans((tar.vector-cand.vectors)^2, na.rm = TRUE))
              dist[is.nan(dist) | is.na(dist)] <- Inf
              dist[dist==0] <- ifelse(is.finite(min(dist[dist>0])), min(dist[dist>0])/2, 1)
              
              if (sum(is.finite(dist)) < k) {
                stop(message = "Fewer than K finite distances found")
              }
              k.features.ind <- order(dist)[1:k]
              k.features <- cand.features[k.features.ind]
              wghts <- 1/dist[k.features.ind]/sum(1/dist[k.features.ind])
              imp.knn[feature.num, exp.num] <- wghts %*% data[k.features, exp.num]
            }
          }
          
          if (!rm.na) {
            imp.knn[is.na(data) == TRUE & is.nan(data) == FALSE] <- NA
          }
          if (!rm.inf) {
            index <- is.finite(data) == FALSE & is.na(data) == FALSE &
              is.nan(data) == FALSE
            imp.knn[index] <- data[index]
          }
          if (!rm.nan) {
            imp.knn[is.nan(data) == TRUE] <- NaN
          }
          return(imp.knn)
        }
        
        ####################################################################
        #### This Function imputes the data based on KNN-CORRELATION or
        #### KNN-TRUNCATION. The Parameter Estimates based on the Truncated
        #### Normal from EstimateComputation function is run on this function
        ####################################################################
        
        imputeKNN <- function (data, k , distance = "correlation",
                               rm.na = TRUE, rm.nan = TRUE, rm.inf = TRUE, perc=1,...) {
          
          if (!(is.matrix(data))) {
            stop(message = paste(deparse(substitute(data)),
                                 " is not a matrix.", sep = ""))
          }
          
          distance <- match.arg(distance, c("correlation","truncation"))
          
          nr <- dim(data)[1]
          if (k < 1 | k > nr) {
            stop(message = "k should be between 1 and the number of rows")
          }
          
          if (distance=="correlation"){
            genemeans<-rowMeans(data,na.rm=TRUE)
            genesd<-apply(data, 1, function(x) sd(x, na.rm = TRUE))
            data<-(data-genemeans)/genesd
          }
          
          if (distance=="truncation"){
            
            ParamMat <- EstimatesComputation(data, perc = perc)
            
            genemeans<-ParamMat[,1]
            genesd<-ParamMat[,2]
            data<-(data-genemeans)/genesd
          }
          
          imp.knn <- data
          imp.knn[is.finite(data) == FALSE] <- NA
          t.data<-t(data)
          
          mv.ind <- which(is.na(imp.knn), arr.ind = TRUE)
          arrays <- unique(mv.ind[, 2])
          array.ind <- match(arrays, mv.ind[, 2])
          ngenes <- 1:nr
          
          for (i in 1:length(arrays)) {
            set <- array.ind[i]:min((array.ind[(i + 1)] - 1), dim(mv.ind)[1],
                                    na.rm = TRUE)
            cand.genes <- ngenes[-unique(mv.ind[set, 1])]
            cand.vectors <- t.data[,cand.genes]
            exp.num<- arrays[i]
            for (j in set) {
              
              gene.num <- mv.ind[j, 1]
              tar.vector <- data[gene.num,]
              
              r <- (cor(cand.vectors,tar.vector, use = "pairwise.complete.obs"))
              dist <- switch(distance,
                             correlation = (1 - abs(r)),
                             truncation = (1 - abs(r)))
              dist[is.nan(dist) | is.na(dist)] <- Inf
              dist[dist==0]<-ifelse(is.finite(min(dist[dist>0])), min(dist[dist>0])/2, 1)
              dist[abs(r) == 1] <- Inf
              
              if (sum(is.finite(dist)) < k) {
                stop(message = "Fewer than K finite distances found")
              }
              k.genes.ind <- order(dist)[1:k]
              k.genes <- cand.genes[k.genes.ind]
              
              wghts <- (1/dist[k.genes.ind]/sum(1/dist[k.genes.ind])) * sign(r[k.genes.ind])
              imp.knn[gene.num, exp.num] <- wghts %*% data[k.genes, exp.num]
            }
          }
          
          if (distance=="correlation") {
            imp.knn <- (imp.knn * genesd) + genemeans
          }
          
          if(distance=="truncation") {
            imp.knn <- (imp.knn * genesd) + genemeans
          }
          
          if (!rm.na) {
            imp.knn[is.na(data) == TRUE & is.nan(data) == FALSE] <- NA
          }
          if (!rm.inf) {
            index <- is.finite(data) == FALSE & is.na(data) == FALSE &
              is.nan(data) == FALSE
            imp.knn[index] <- data[index]
          }
          if (!rm.nan) {
            imp.knn[is.nan(data) == TRUE] <- NaN
          }
          return(imp.knn)
        }
      }
      result <- df_orig %>%
        as.matrix() %>%
        t() %>%
        imputeKNN(k = 10, distance = "truncation", perc = 0) %>%
        t()
      as.data.frame(result)
    },
    rf = {
      if (!requireNamespace("missForest", quietly = TRUE)) {
        stop("Package 'missForest' required for rf", call. = FALSE)
      }
      out <- missForest::missForest(t(df_orig), maxiter = 10,
                                    mtry = floor(nrow(df_orig)^(1/3)),
                                    verbose = TRUE)
      t(out$ximp)
    },
    pi = {
      # Paramters `width` and `downshift` should be passed in via `...` or environment
      width     <- getOption("nafuncs.pi.width", default = 1)
      downshift <- getOption("nafuncs.pi.downshift", default = 1)
      df_work[] <- lapply(df_orig, function(col) {
        if (anyNA(col)) {
          non_na   <- col[!is.na(col)]
          sd_val   <- width * sd(non_na, na.rm = TRUE)
          mean_val <- mean(non_na, na.rm = TRUE) - downshift * sd(non_na, na.rm = TRUE)
          col[is.na(col)] <- rnorm(sum(is.na(col)), mean = mean_val, sd = sd_val)
        }
        col
      })
      df_work
    },
    grr = {
      if (!requireNamespace("DreamAI", quietly = TRUE)) {
        stop("Package 'DreamAI' required for grr", call. = FALSE)
      }
      if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Package 'glmnet' required for grr", call. = FALSE)
      }
      if (!requireNamespace("import", quietly=TRUE)) install.packages("import")
      import::into(.into = "DreamAI", .from = "glmnet", cv.glmnet) # I really don’t want to attach glmnet, just inject the function into DreamAI’s namespace
      DreamAI::impute.RegImpute(as.matrix(df_orig),
                                fillmethod = "row_mean",
                                maxiter_RegImpute = 10,
                                conv_nrmse = 1e-03)
    },
    gms = {
      if (!requireNamespace("GMSimpute", quietly = TRUE)) {
        stop("Package 'GMSimpute' required for gms", call. = FALSE)
      }
      GMSimpute::GMS.Lasso(df_orig, nfolds = 3, log.scale = FALSE, TS.Lasso = TRUE)
    },
    
    stop(sprintf("Unsupported method '%s'", method), call. = FALSE)
  )
  
  as.data.frame(imputed)
}

