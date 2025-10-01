#### Tcrit chlorophyll fluorescence — Batch (base_dir version) ####
# Raw files:    scripts/Thermo_Water/data_raw/    YYYYMMDD_runN_data.TXT
# Label files:  scripts/Thermo_Water/data_labels/ YYYYMMDD_runN_labels.csv  (Grid.ID, Sample)
# Outputs:      scripts/Thermo_Water/data_processed/<YYYYMMDD_runN_results>/ (PNGs + per-run CSV)
# Master CSV:   scripts/Thermo_Water/data_processed/_tcrit_master_results.csv
# Last updated: 2025-10-01

suppressPackageStartupMessages({
  library(segmented)
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(janitor)
  library(here)
  library(tools)
})

# -------------------------------- #
#            PATHS/IO              #
# -------------------------------- #
# Make sure your project root is detected correctly by {here}.
# If needed, uncomment the next line once in your project to lock the root:
# here::i_am("scripts/Thermo_Water/tcrit_batch.R")

base_dir      <- here::here("scripts", "Thermo_Water")
raw_dir       <- file.path(base_dir, "data_raw")
labels_dir    <- file.path(base_dir, "data_labels")
processed_dir <- file.path(base_dir, "data_processed")

dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(labels_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------- #
#        USER PARAMETERS           #
# -------------------------------- #
xlow           <- 31    # °C lower bound for analysis window
xhigh          <- 85    # °C upper bound for analysis window
maxthreshold   <- 0.9   # scaled fluorescence threshold (0–1)
plots_per_page <- 16

# -------------------------------- #
#         HELPER FUNCTIONS         #
# -------------------------------- #

rescale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

# Collapse each 5-row block by mean
condense_5 <- function(df){
  stopifnot(nrow(df) %% 5 == 0)
  block_id <- rep(seq_len(nrow(df)/5), each = 5)
  as.data.frame(lapply(df, function(col) tapply(col, block_id, mean, na.rm = TRUE)))
}

# Expected 36 positions (row-major): A1..A8, B1..B8, C1..C8, D1..D8, E1..E4
grid_ids_36 <- function(){
  c(paste0("A", 1:8),
    paste0("B", 1:8),
    paste0("C", 1:8),
    paste0("D", 1:8),
    paste0("E", 1:8),
    paste0("F", 1:8))
}

# Parse "YYYYMMDD_runN_data.TXT" (case-insensitive handled at list.files)
parse_data_fname <- function(fname){
  m <- stringr::str_match(fname, "^(\\d{8})_run(\\d+)_data\\.(TXT|txt)$")
  if (any(is.na(m))) return(NULL)
  list(date = m[2], run = m[3])
}

# Rename fluorscale columns using labels (Grid.ID -> Sample)
rename_from_labels <- function(fluorscale, labels_df){
  labels_df <- labels_df %>% janitor::clean_names()
  if (!all(c("grid_id", "sample") %in% names(labels_df))){
    stop("Label file must have columns 'Grid.ID' and 'Sample' (any case).")
  }
  n_areas <- ncol(fluorscale) - 1
  if (n_areas != 36) stop(sprintf("Expected 36 areas, found %d.", n_areas))
  
  grids <- grid_ids_36()
  sample_names <- tibble(grid_id = grids) %>%
    left_join(labels_df %>% select(grid_id, sample), by = "grid_id") %>%
    pull(sample)
  
  if (anyNA(sample_names)){
    miss <- grids[is.na(sample_names)]
    stop(sprintf("Missing Sample names for Grid.ID(s): %s", paste(miss, collapse = ", ")))
  }
  
  names(fluorscale)[-1] <- make.unique(as.character(sample_names))
  fluorscale
}

# Core per-run processing
process_one_run <- function(raw_path, labels_path, processed_dir){
  raw_file <- basename(raw_path)
  p <- parse_data_fname(raw_file)
  if (is.null(p)) stop("Unrecognized data filename: ", raw_file)
  
  date_chr <- p$date
  run_num  <- p$run
  
  # Read raw kinetics: skip first 2 metadata lines, tab-separated, keep headers
  rawFC <- suppressMessages(read.table(raw_path, skip = 2, header = TRUE, sep = "\t", check.names = FALSE))
  if (ncol(rawFC) < 2) stop("Raw file has fewer than 2 columns after header.")
  
  # Drop timestamp, rename first to 'temp'
  rawFC <- rawFC[-1]
  names(rawFC)[1] <- "temp"
  
  # Condense 5x and subset temperature window
  fluor <- condense_5(rawFC) %>% dplyr::filter(temp > xlow, temp < xhigh)
  
  # Scale 0–1 for all channels except temp
  fluorscale <- fluor
  for (i in 2:ncol(fluorscale)) fluorscale[[i]] <- rescale01(fluorscale[[i]])
  
  # Read labels and rename columns
  labels_df  <- suppressMessages(readr::read_csv(labels_path, show_col_types = FALSE))
  fluorscale <- rename_from_labels(fluorscale, labels_df)
  
  # Output directory for this run
  run_tag     <- paste0(date_chr, "_run", run_num, "_results")
  run_out_dir <- file.path(processed_dir, run_tag)
  dir.create(run_out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ---------- Page PNGs: scaled fluorescence ----------
  graph_index <- 0L
  page_index  <- 1L
  open_dev    <- FALSE
  
  for (i in 2:ncol(fluorscale)){
    if ((graph_index %% plots_per_page) == 0){
      if (open_dev) dev.off()
      png(filename = file.path(run_out_dir, sprintf("fluor_page_%02d.png", page_index)),
          width = 1600, height = 1600, res = 150)
      par(mfrow = c(4,4), mar = c(3,4,2,1))
      open_dev  <- TRUE
      page_index <- page_index + 1L
    }
    plot(fluorscale$temp, fluorscale[[i]],
         type = "l", xlim = c(xlow, xhigh), ylim = c(0,1),
         ylab = names(fluorscale)[i], xlab = "Temperature (°C)")
    abline(v = xlow, lty = 3, col = "blue")
    graph_index <- graph_index + 1L
  }
  if (open_dev) dev.off()
  
  # ---------- Segmented regression per channel ----------
  res <- tibble::tibble(
    Date  = date_chr,
    Run   = paste0("run", run_num),
    LeafID = names(fluorscale)[-1],
    Tcrit = NA_real_,
    StdErr = NA_real_,
    T50 = NA_real_,
    xlow = xlow,
    xhigh = xhigh,
    maxthreshold = maxthreshold
  )
  
  for (i in 2:ncol(fluorscale)){
    nm <- names(fluorscale)[i]
    # first temp where channel exceeds threshold
    idx_thr <- which(fluorscale[[i]] > maxthreshold)
    if (!length(idx_thr)) next
    tempatthreshold <- fluorscale$temp[idx_thr[1]]
    
    fluor_sub <- fluorscale |>
      dplyr::filter(temp > xlow, temp < tempatthreshold) |>
      dplyr::filter(.data[[nm]] < maxthreshold)
    
    if (nrow(fluor_sub) < 10) next
    
    # T50
    t50 <- fluor_sub$temp[ which.min(abs(fluor_sub[[nm]] - 0.5)) ]
    
    # Fit segmented
    model1 <- try(lm(fluor_sub[[nm]] ~ temp, data = fluor_sub), silent = TRUE)
    if (inherits(model1, "try-error")) next
    seg_model1 <- try(segmented::segmented(model1, seg.Z = ~ temp, data = fluor_sub), silent = TRUE)
    if (inherits(seg_model1, "try-error")) next
    
    psi <- seg_model1$psi
    tcrit1    <- if (is.matrix(psi) && nrow(psi) >= 1) round(psi[1, "Est."], 3) else NA_real_
    tcriterr1 <- if (is.matrix(psi) && "St.Err" %in% colnames(psi)) round(psi[1, "St.Err"], 3) else NA_real_
    
    # Fitted lines for plot
    fitted_val1 <- fitted(seg_model1)
    breakmodel1 <- data.frame(Temperature = fluor_sub$temp, fluor_sub = fitted_val1)
    
    # Save PNG per channel
    png(filename = file.path(run_out_dir, paste0(gsub("[^A-Za-z0-9_\\-]+", "_", nm), "_seg.png")),
        width = 900, height = 600)
    par(mar = c(4,4,2,1))
    plot(fluorscale[[i]] ~ fluorscale$temp,
         ylab = nm, xlab = "Temperature (°C)",
         ylim = c(0,1), xlim = c(xlow, xhigh), type = "l")
    abline(v = xlow, lty = 3, col = "blue")
    abline(v = tempatthreshold, lty = 3, col = "blue")
    lines(breakmodel1$Temperature, breakmodel1$fluor_sub, col = "red", lwd = 2)
    points(t50,
           fluor_sub[[nm]][ which.min(abs(fluor_sub[[nm]] - 0.5)) ],
           pch = 22, bg = "green", cex = 1.2)
    legend("top",
           legend = c(
             sprintf("Tcrit = %s °C", ifelse(is.na(tcrit1), "NA", format(tcrit1, nsmall = 3))),
             sprintf("T50 = %s °C",   ifelse(is.na(t50),    "NA", format(round(t50,3), nsmall = 3)))
           ),
           bty = "n")
    dev.off()
    
    # Fill results
    res$Tcrit[res$LeafID == nm]  <- tcrit1
    res$StdErr[res$LeafID == nm] <- tcriterr1
    res$T50[res$LeafID == nm]    <- round(t50, 3)
  }
  
  # Write per-run CSV
  per_run_csv <- file.path(run_out_dir, paste0(tools::file_path_sans_ext(raw_file), "_results.csv"))
  readr::write_csv(res, per_run_csv)
  
  invisible(res)
}

# -------------------------------- #
#             DRIVER               #
# -------------------------------- #

# Find all raw data files under base_dir/data_raw/ (case-insensitive)
raw_files <- list.files(
  raw_dir,
  pattern = "(?i)^[0-9]{8}_run[0-9]+_data\\.txt$",
  full.names = TRUE
)

if (!length(raw_files)) {
  message("No raw data files found in ", raw_dir,
          ". Expecting 'YYYYMMDD_runN_data.TXT'.")
}

# Expected label path for each raw file (base_dir/data_labels/)
label_for_raw <- function(raw_path){
  b <- basename(raw_path)
  p <- parse_data_fname(b)
  if (is.null(p)) return(NA_character_)
  file.path(labels_dir, paste0(p$date, "_run", p$run, "_labels.csv"))
}

pairs <- tibble::tibble(
  raw    = raw_files,
  labels = purrr::map_chr(raw_files, label_for_raw),
  have_labels = file.exists(labels)
)

if (any(!pairs$have_labels)){
  warning("Missing label files for:\n  - ",
          paste(basename(pairs$raw[!pairs$have_labels]), collapse = "\n  - "),
          "\nLooked in: ", labels_dir)
}

to_process <- dplyr::filter(pairs, have_labels)

all_results <- vector("list", length = 0)

for (k in seq_len(nrow(to_process))){
  raw_path    <- to_process$raw[k]
  labels_path <- to_process$labels[k]
  cat("\n=== Processing:", basename(raw_path), "with", basename(labels_path), "===\n")
  res <- try(process_one_run(raw_path, labels_path, processed_dir), silent = TRUE)
  if (inherits(res, "try-error")){
    message("  -> ERROR on ", basename(raw_path), ": ", as.character(res))
  } else {
    all_results[[length(all_results) + 1]] <- res
  }
}

# Master CSV in base_dir/data_processed/
if (length(all_results)){
  master <- dplyr::bind_rows(all_results) %>% dplyr::arrange(Date, Run, LeafID)
  readr::write_csv(master, file.path(processed_dir, "_tcrit_master_results.csv"))
  cat("\nWrote master results to ", file.path(processed_dir, "_tcrit_master_results.csv"), "\n", sep = "")
} else {
  cat("\nNo runs completed successfully—check errors above.\n")
}
