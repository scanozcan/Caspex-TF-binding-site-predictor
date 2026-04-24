# =============================================================================
# caspex_analysis_hocomoco.R
# CasPEX Binding Zone Predictor — HOCOMOCO backend
#
# Same spatial model, same deconvolution, same plots as caspex_analysis.R.
# The only difference is the motif database: TF position weight matrices
# are read from a locally cached HOCOMOCO MEME bundle instead of fetched
# from JASPAR REST. This avoids per-TF network latency and works offline
# after the first download.
#
# Usage:
#   source("caspex_analysis_hocomoco.R")        # also sources caspex_analysis.R
#
#   result <- run_caspex_hocomoco(
#     gene             = "ATP7B",
#     grnas            = c(R1 = "GGG...", R2 = "...", ...),
#     data_files       = c(R1 = "Region1.txt", ...),
#     out_dir          = "caspex_output_hocomoco",
#     hocomoco_version = "v12",       # "v12" (default) or "v11"
#     hocomoco_species = "human"      # "human" (default) or "mouse"
#   )
#
# The first call downloads the HOCOMOCO MEME bundle (~few MB) to a local
# cache directory (see tools::R_user_dir("caspex", "cache")). Subsequent
# calls are fast. To force a redownload, call
# `download_hocomoco_bundle(version, species, force = TRUE)`.
#
# Differences vs the JASPAR version at the interface level:
#   - Motif source: HOCOMOCO v12 (CORE) or v11 (full)
#   - One "best-rated" PWM per gene is kept (HOCOMOCO rates matrices A > B > C > D)
#   - `result$motif_source` now records the active database
#   - Plot titles saying "JASPAR" are rewritten to "HOCOMOCO vXX" post-hoc
# =============================================================================

# --- 1. Source the base pipeline so we inherit every non-motif function ----
.caspex_base <- "caspex_analysis.R"
if (!exists("run_caspex", mode = "function")) {
  if (!file.exists(.caspex_base))
    stop("caspex_analysis_hocomoco.R expects '", .caspex_base,
         "' in the same directory. Missing file.")
  source(.caspex_base, chdir = FALSE)
}

# Sanity check the handful of helpers we hook into
for (.needed in c("run_caspex", "scan_sequence", "score_pwm_positions"))
  if (!exists(.needed, mode = "function"))
    stop("caspex_analysis.R did not define `", .needed,
         "` — cannot continue.")
rm(.needed, .caspex_base)

# --- 2. Package dependencies (base already loads httr/jsonlite/ggplot2…) ---
suppressPackageStartupMessages({
  if (!requireNamespace("utils", quietly = TRUE))
    stop("base R 'utils' package required")
})

# =============================================================================
# 3. HOCOMOCO URL table and local cache
# =============================================================================

#' Resolve the MEME bundle URL(s) for a HOCOMOCO version × species.
#' Returns a character vector — we try each URL in order until one works,
#' since HOCOMOCO has moved paths between releases.
#' @keywords internal
.hocomoco_url <- function(version, species) {
  stopifnot(species %in% c("human", "mouse"))
  stopifnot(version %in% c("v11", "v12"))
  if (version == "v12") {
    if (species == "human") {
      return(c(
        # current (April 2026) location — under formatted_motifs/
        paste0("https://hocomoco12.autosome.org/final_bundle/",
               "hocomoco12/H12CORE/formatted_motifs/",
               "H12CORE_meme_format.meme"),
        # legacy location (pre-2024) kept as fallback
        paste0("https://hocomoco12.autosome.org/final_bundle/",
               "hocomoco12/H12CORE/H12CORE_meme_format.meme")))
    } else {
      # HOCOMOCO v12 CORE is human-focused; INVIVO bundle contains mouse PWMs.
      return(c(
        paste0("https://hocomoco12.autosome.org/final_bundle/",
               "hocomoco12/H12INVIVO/formatted_motifs/",
               "H12INVIVO_meme_format.meme"),
        paste0("https://hocomoco12.autosome.org/final_bundle/",
               "hocomoco12/H12CORE_mouse/formatted_motifs/",
               "H12CORE_mouse_meme_format.meme"),
        paste0("https://hocomoco12.autosome.org/final_bundle/",
               "hocomoco12/H12CORE_mouse/",
               "H12CORE_mouse_meme_format.meme")))
    }
  }
  if (version == "v11") {
    return(switch(species,
      "human" = paste0("https://hocomoco11.autosome.org/final_bundle/",
                       "hocomoco11/full/HUMAN/mono/",
                       "HOCOMOCOv11_full_HUMAN_mono_meme_format.meme"),
      "mouse" = paste0("https://hocomoco11.autosome.org/final_bundle/",
                       "hocomoco11/full/MOUSE/mono/",
                       "HOCOMOCOv11_full_MOUSE_mono_meme_format.meme")))
  }
}

#' Local cache path for a HOCOMOCO bundle. Uses tools::R_user_dir() so the
#' cache survives across sessions and is OS-appropriate.
#' @keywords internal
.hocomoco_cache_path <- function(version, species) {
  base <- tools::R_user_dir("caspex", which = "cache")
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  file.path(base, sprintf("hocomoco_%s_%s.meme", version, species))
}

#' Download a HOCOMOCO MEME bundle to the local cache (once).
#'
#' @param version   "v12" (default) or "v11"
#' @param species   "human" (default) or "mouse"
#' @param force     Redownload even if the cache file already exists
#' @param from_file If the automatic download is blocked (corporate firewall,
#'                  offline VM, etc.), pass the path to a manually-downloaded
#'                  MEME file and it will be copied into the cache slot.
#' @return The path to the cached bundle (invisibly)
download_hocomoco_bundle <- function(version = "v12", species = "human",
                                      force = FALSE, from_file = NULL) {
  dest <- .hocomoco_cache_path(version, species)
  if (file.exists(dest) && !force) return(invisible(dest))

  # --- Manual-file shortcut ------------------------------------------------
  if (!is.null(from_file)) {
    if (!file.exists(from_file))
      stop("from_file does not exist: ", from_file)
    if (file.size(from_file) < 1000)
      stop("from_file looks too small to be a MEME bundle: ", from_file)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    ok <- file.copy(from_file, dest, overwrite = TRUE)
    if (!ok) stop("Failed to copy ", from_file, " -> ", dest)
    message("Cached MEME bundle from local file:\n  ", from_file,
            "\n  -> ", dest, " (",
            format(file.size(dest), big.mark = ","), " bytes)")
    return(invisible(dest))
  }

  # --- Network download: try each candidate URL in order ------------------
  urls <- .hocomoco_url(version, species)
  last_err <- NULL
  for (url in urls) {
    message("Downloading HOCOMOCO ", version, " (", species, ")\n  ", url)
    tmp <- tempfile(fileext = ".meme")
    ok <- tryCatch({
      # libcurl follows redirects and handles HTTPS cleanly on all platforms
      utils::download.file(url, tmp, mode = "wb", quiet = FALSE,
                           method = "libcurl")
      TRUE
    }, error = function(e) { last_err <<- e; FALSE },
       warning = function(w) { last_err <<- w; FALSE })
    # Be strict: a MEME bundle is always > 100 KB
    if (ok && file.exists(tmp) && file.size(tmp) > 1e5) {
      dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
      file.rename(tmp, dest)
      message("  cached at ", dest, " (",
              format(file.size(dest), big.mark = ","), " bytes)")
      return(invisible(dest))
    }
    try(unlink(tmp), silent = TRUE)
    message("  ... that URL did not return a valid MEME bundle, trying next.")
  }

  # --- All URLs failed — give the user an actionable message --------------
  stop(sprintf(
"Could not download HOCOMOCO %s %s from any known URL.

Tried:
%s

Workaround — download the MEME bundle manually in a browser and point
the pipeline at the local file:

  download_hocomoco_bundle(
    version  = %s,
    species  = %s,
    from_file = \"/path/to/your/downloaded.meme\"
  )

Typical download page: https://hocomoco12.autosome.org/downloads_v12
Last error was:
  %s
",
    version, species,
    paste0("  - ", urls, collapse = "\n"),
    dQuote(version, q = FALSE), dQuote(species, q = FALSE),
    if (!is.null(last_err)) conditionMessage(last_err) else "(unknown)"),
    call. = FALSE)
}

#' Clear the in-memory HOCOMOCO cache (does NOT delete the on-disk file).
clear_hocomoco_memory_cache <- function() {
  rm(list = ls(.hocomoco_env), envir = .hocomoco_env)
  invisible(TRUE)
}

# =============================================================================
# 4. MEME parser (subset — sufficient for HOCOMOCO format)
# =============================================================================

#' Parse a MEME-format motif file into a named list of PPMs.
#'
#' Return value: list keyed by motif ID, each element a list with
#'   $id      : full motif identifier (e.g. "SOX2.H12CORE.0.P.B")
#'   $alt     : the optional alt-name on the MOTIF line
#'   $gene    : extracted gene symbol (uppercase)
#'   $ppm     : 4 × w numeric matrix of base probabilities (rows A,C,G,T)
#'   $w       : motif width
parse_meme_file <- function(path) {
  lines <- readLines(path, warn = FALSE)
  motifs <- list()
  n <- length(lines)
  i <- 1
  while (i <= n) {
    if (grepl("^MOTIF\\s", lines[i])) {
      toks <- strsplit(trimws(sub("^MOTIF\\s+", "", lines[i])), "\\s+")[[1]]
      motif_id <- toks[1]
      alt_name <- if (length(toks) >= 2) toks[2] else motif_id
      # Scan for letter-probability matrix header
      j <- i + 1
      while (j <= n && !grepl("^letter-probability matrix", lines[j])) j <- j + 1
      if (j > n) { i <- i + 1; next }
      wmatch <- regmatches(lines[j], regexpr("w=\\s*[0-9]+", lines[j]))
      if (length(wmatch) == 0) { i <- j + 1; next }
      w <- as.integer(sub("w=\\s*", "", wmatch))
      if (is.na(w) || w <= 0) { i <- j + 1; next }
      ppm <- matrix(NA_real_, 4, w)
      ok <- TRUE
      for (k in seq_len(w)) {
        if (j + k > n) { ok <- FALSE; break }
        vals <- suppressWarnings(scan(text = lines[j + k],
                                       what = numeric(),
                                       quiet = TRUE, n = 4))
        if (length(vals) != 4) { ok <- FALSE; break }
        ppm[, k] <- vals
      }
      if (ok) {
        motifs[[motif_id]] <- list(id    = motif_id,
                                   alt   = alt_name,
                                   gene  = .extract_gene_symbol(motif_id),
                                   ppm   = ppm,
                                   w     = w)
      }
      i <- j + w + 1
    } else {
      i <- i + 1
    }
  }
  motifs
}

#' Extract a gene symbol from a HOCOMOCO motif identifier.
#' Handles both v11 (`SOX2_HUMAN.H11MO.0.A`) and v12 (`SOX2.H12CORE.0.P.B`)
#' naming conventions.
#' @keywords internal
.extract_gene_symbol <- function(motif_id) {
  g <- strsplit(motif_id, "\\.")[[1]][1]
  g <- sub("_(HUMAN|MOUSE|RAT)$", "", g, ignore.case = TRUE)
  toupper(g)
}

#' HOCOMOCO quality rating from the trailing letter of the motif id.
#' Lower is better (A = 1, B = 2, …). Missing or non-letter → 99 (worst).
#' @keywords internal
.extract_quality <- function(motif_id) {
  parts <- strsplit(motif_id, "\\.")[[1]]
  last  <- parts[length(parts)]
  if (nchar(last) == 1 && last %in% LETTERS)
    return(match(last, LETTERS))
  99
}

# =============================================================================
# 5. Cached loader (on-disk bundle → in-memory, best-rated PWM per gene)
# =============================================================================

.hocomoco_env <- new.env(parent = emptyenv())

#' Load HOCOMOCO PWMs into memory, indexed by gene symbol.
#' Downloads the bundle on first use. Within one session, subsequent calls
#' are free (cached in `.hocomoco_env`).
load_hocomoco_pwms <- function(version = "v12", species = "human") {
  key <- paste(version, species, sep = "_")
  if (!is.null(.hocomoco_env[[key]])) return(.hocomoco_env[[key]])

  bundle <- .hocomoco_cache_path(version, species)
  if (!file.exists(bundle)) download_hocomoco_bundle(version, species)

  message("Parsing HOCOMOCO ", version, " (", species, ") bundle...")
  motifs <- parse_meme_file(bundle)
  if (length(motifs) == 0)
    stop("No motifs parsed from ", bundle,
         " — file may be corrupt. Try `download_hocomoco_bundle(..., force=TRUE)`")

  # Collapse to one PWM per gene, preferring the lowest (= best) quality code.
  by_gene <- list()
  for (nm in names(motifs)) {
    g <- motifs[[nm]]$gene
    q <- .extract_quality(nm)
    if (is.null(by_gene[[g]]) || q < by_gene[[g]]$q)
      by_gene[[g]] <- list(motif = motifs[[nm]], q = q)
  }
  result <- lapply(by_gene, `[[`, "motif")
  message("  -> ", length(result), " unique TFs from ",
          length(motifs), " total matrices")
  .hocomoco_env[[key]] <- result
  result
}

# =============================================================================
# 6. fetch_hocomoco_pwm — drop-in replacement for fetch_jaspar_pwm
# =============================================================================

#' Fetch a HOCOMOCO PWM for a TF symbol.
#'
#' Return value mirrors fetch_jaspar_pwm():
#'   list(id, name, pwm = 4 × L log-odds matrix, len)
fetch_hocomoco_pwm <- function(tf_name, version = "v12", species = "human") {
  pwms <- load_hocomoco_pwms(version, species)
  key <- toupper(tf_name)
  m <- pwms[[key]]
  if (is.null(m)) return(NULL)

  # PPM -> log-odds PWM with uniform background (same convention as JASPAR code)
  ppm <- pmax(m$ppm, 1e-6)
  pwm <- log2(ppm / 0.25)
  list(id = m$id, name = tf_name, pwm = pwm, len = m$w)
}

# =============================================================================
# 7. run_motif_scan_hocomoco — same signature & return shape as run_motif_scan
# =============================================================================

#' HOCOMOCO drop-in replacement for run_motif_scan().
#'
#' Returns a named list keyed by TF, each element
#'   list(pwm = <fetch_hocomoco_pwm output>, hits = <TSS-relative positions>,
#'        n_hits = integer)
run_motif_scan_hocomoco <- function(tf_names, promoter_info,
                                     threshold_frac = 0.80,
                                     version = "v12",
                                     species = "human") {
  message("\nRunning HOCOMOCO ", version, " motif scan for ",
          length(tf_names), " TFs...")
  out <- list()
  for (tf in tf_names) {
    message("  ", tf, " ...", appendLF = FALSE)
    pwm <- tryCatch(fetch_hocomoco_pwm(tf, version, species),
                    error = function(e) {
                      message(" error: ", conditionMessage(e)); NULL
                    })
    if (is.null(pwm)) { message(" no motif found"); next }
    hits <- scan_sequence(promoter_info, pwm, threshold_frac)
    message(" [", pwm$id, "] ", length(hits), " hit(s)")
    out[[tf]] <- list(pwm = pwm, hits = hits, n_hits = length(hits))
  }
  out
}

#' Build a closure with the (version, species) baked in, signature-compatible
#' with the base run_motif_scan().
#' @keywords internal
.make_hocomoco_scan <- function(version = "v12", species = "human") {
  function(tf_names, promoter_info, threshold_frac = 0.80) {
    run_motif_scan_hocomoco(tf_names, promoter_info, threshold_frac,
                             version = version, species = species)
  }
}

# =============================================================================
# 8. Post-hoc title fixer (swap "JASPAR" → "HOCOMOCO" in plot titles)
# =============================================================================

.rewrite_jaspar_in_plots <- function(plots, tag) {
  fix_one <- function(p) {
    if (inherits(p, "ggplot")) {
      if (!is.null(p$labels$title))
        p$labels$title <- gsub("JASPAR", tag, p$labels$title, fixed = TRUE)
      if (!is.null(p$labels$subtitle))
        p$labels$subtitle <- gsub("JASPAR", tag, p$labels$subtitle,
                                   fixed = TRUE)
    }
    p
  }
  walker <- function(x) {
    if (inherits(x, "ggplot")) return(fix_one(x))
    if (is.list(x)) return(lapply(x, walker))
    x
  }
  walker(plots)
}

# =============================================================================
# 9. run_caspex_hocomoco — top-level entry point
# =============================================================================

#' Run the full CasPEX pipeline using HOCOMOCO as the motif source.
#'
#' All base arguments are forwarded to run_caspex(). The motif-scanning
#' step is transparently swapped from JASPAR to HOCOMOCO. Two additional
#' arguments control the HOCOMOCO source:
#'
#' @param ...                Any argument accepted by run_caspex()
#' @param hocomoco_version   "v12" (default) or "v11"
#' @param hocomoco_species   "human" (default) or "mouse"
#' @param hocomoco_from_file Optional path to a manually-downloaded MEME
#'                           bundle. Use this if automatic download is
#'                           blocked (corporate firewall, offline VM).
#'                           The file is copied into the cache slot for
#'                           `(version, species)` before the pipeline runs.
#' @return The same invisible list as run_caspex(), plus:
#'           $motif_source   e.g. "HOCOMOCO v12 (human)"
run_caspex_hocomoco <- function(...,
                                 hocomoco_version = "v12",
                                 hocomoco_species = "human",
                                 hocomoco_from_file = NULL) {
  if (!exists("run_caspex", envir = globalenv(), mode = "function"))
    stop("run_caspex() not found in globalenv — did caspex_analysis.R load?")

  # If a manual MEME file was provided, seed the cache from it first.
  if (!is.null(hocomoco_from_file)) {
    download_hocomoco_bundle(hocomoco_version, hocomoco_species,
                             force = TRUE, from_file = hocomoco_from_file)
    clear_hocomoco_memory_cache()   # drop any stale in-memory PWM dict
  }

  # Ensure the bundle is available *before* we start, so we fail fast if
  # the download can't be fetched rather than mid-run.
  invisible(load_hocomoco_pwms(hocomoco_version, hocomoco_species))

  # Swap run_motif_scan in globalenv so that the downstream reference inside
  # run_caspex resolves to the HOCOMOCO version. Restore on exit.
  ge <- globalenv()
  orig <- ge$run_motif_scan
  on.exit(assign("run_motif_scan", orig, envir = ge), add = TRUE)
  assign("run_motif_scan",
         .make_hocomoco_scan(hocomoco_version, hocomoco_species),
         envir = ge)

  # Run
  tag <- paste0("HOCOMOCO ", hocomoco_version, " (", hocomoco_species, ")")
  message("\n========== CasPEX run with ", tag, " ==========\n")

  # --- Capture and eagerly force `...` so we fail fast with a clear error
  #     if the caller accidentally included a literal `...` placeholder
  #     (a common copy-paste trap, e.g. `grnas = c(R1 = "...", ...)`).
  dots <- list(...)
  for (nm in names(dots)) {
    if (is.null(nm) || !nzchar(nm)) next
    ok <- tryCatch({ force(dots[[nm]]); TRUE },
      error = function(e) {
        msg <- conditionMessage(e)
        if (grepl("'...' used in an incorrect context", msg, fixed = TRUE)) {
          stop("Argument `", nm, "` contains a literal `...` placeholder.\n",
               "R treats `...` as a special symbol, so e.g. ",
               "`c(R1 = \"Region1.txt\", ...)` fails at the top level.\n",
               "Replace any literal `...` in your call with actual values.\n",
               "Original R error: ", msg, call. = FALSE)
        }
        stop("Could not evaluate argument `", nm, "`: ", msg, call. = FALSE)
      })
  }

  res <- do.call(run_caspex, dots)

  # Rewrite plot titles that mention "JASPAR" so the output PDFs are
  # truthful about the motif source.
  if (!is.null(res$plots))
    res$plots <- .rewrite_jaspar_in_plots(res$plots, tag)

  # Record which motif database was used
  res$motif_source    <- tag
  res$hocomoco_version <- hocomoco_version
  res$hocomoco_species <- hocomoco_species

  message("\n(Motif source: ", tag, ")")
  invisible(res)
}

# =============================================================================
# 10. Convenience helpers
# =============================================================================

#' List all TFs available in the loaded HOCOMOCO bundle (by gene symbol).
list_hocomoco_tfs <- function(version = "v12", species = "human") {
  sort(names(load_hocomoco_pwms(version, species)))
}

#' Check how many TFs in a given set are present in HOCOMOCO.
check_hocomoco_coverage <- function(tf_names,
                                     version = "v12", species = "human") {
  tfs <- toupper(tf_names)
  available <- names(load_hocomoco_pwms(version, species))
  found <- tfs %in% available
  message(sum(found), " / ", length(tfs),
          " TFs found in HOCOMOCO ", version, " (", species, ")")
  missing <- tfs[!found]
  if (length(missing) > 0)
    message("  not found: ", paste(head(missing, 20), collapse = ", "),
            if (length(missing) > 20)
              paste0("  (+", length(missing) - 20, " more)") else "")
  invisible(data.frame(tf = tfs, in_hocomoco = found,
                        stringsAsFactors = FALSE))
}
