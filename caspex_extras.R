# =============================================================================
# caspex_extras.R
# Supplementary result plots & diagnostics for the CasPEX pipeline.
# Source AFTER caspex_analysis.R:
#
#   source("caspex_analysis.R")
#   source("caspex_extras.R")
#   result <- run_caspex(gene = "ATP7B", grnas = ..., data_files = ...)
#   extras <- run_caspex_extras(result, out_dir = "caspex_output/extras")
#
# Every top-level function takes the invisible list returned by run_caspex()
# and returns either a ggplot object, a named list of ggplots, or (for
# expensive computations) a results object that its matching plot_* function
# can render.
#
# Contents
#   A. Biological interpretation
#       plot_tf_one_pager, plot_tf_family_enrichment,
#       plot_event_density, plot_composite_vs_specificity
#   B. Statistical robustness
#       run_permutation_null, plot_permutation_null,
#       run_sigma_sensitivity, plot_sigma_sensitivity,
#       run_event_jackknife, plot_event_jackknife,
#       plot_nnls_residual
#   C. QC / data sanity
#       plot_volcano_per_region, plot_region_correlation,
#       plot_pval_histograms, plot_motif_vs_nnls
#   D. Convenience wrapper
#       run_caspex_extras
# =============================================================================

suppressPackageStartupMessages({
  # Self-sufficient: attach everything the extras functions need, so this
  # file works even if caspex_analysis.R's session was only partially
  # restored (e.g. after saveRDS(result) + restart + readRDS(result)).
  .need_pkgs <- c("ggplot2", "dplyr", "tidyr", "scales", "patchwork",
                  "httr", "jsonlite")
  .missing <- .need_pkgs[
    !vapply(.need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(.missing) > 0)
    stop("caspex_extras.R needs package(s): ",
         paste(.missing, collapse = ", "),
         ". Install with install.packages(c(",
         paste0("'", .missing, "'", collapse = ", "), ")).")
  for (.p in .need_pkgs) library(.p, character.only = TRUE)
  rm(.need_pkgs, .missing, .p)
})

# Verify caspex_analysis.R has been sourced in this session — we need its
# helpers (plot_binding_deconvolution, build_caspex_signal, predict_binding_events,
# compute_region_weight, score_pwm_positions, theme_caspex, COLS).
.caspex_required <- c("plot_binding_deconvolution", "build_caspex_signal",
                      "predict_binding_events", "compute_region_weight",
                      "score_pwm_positions", "theme_caspex", "COLS")
.caspex_missing <- .caspex_required[!vapply(.caspex_required, exists,
                                             logical(1), mode = "any")]
if (length(.caspex_missing) > 0)
  stop("caspex_extras.R requires caspex_analysis.R to be sourced first. ",
       "Missing: ", paste(.caspex_missing, collapse = ", "),
       ". Run: source('caspex_analysis.R'); source('caspex_extras.R')")
rm(.caspex_required, .caspex_missing)

# Null-coalescing operator, used throughout below.
`%||%` <- function(a, b) if (is.null(a)) b else a

# =============================================================================
# A.1  Per-TF one-pager dashboard
# =============================================================================

#' One page summarising a single TF: per-region bars, spatial footprint,
#' binding-event track, and a compact stats block.
#'
#' @param result   The invisible list returned by run_caspex()
#' @param tf_name  TF symbol to summarise
#' @return A patchwork ggplot assembly
plot_tf_one_pager <- function(result, tf_name) {
  if (!tf_name %in% result$long_data$protein)
    stop(tf_name, " not found in result$long_data")

  ld       <- result$long_data
  pos_map  <- result$pos_map
  motif_hits <- if (!is.null(result$motif_results[[tf_name]]))
    result$motif_results[[tf_name]]$hits else integer(0)

  # Top panel: deconvolution view (full signal + events + motif ticks).
  # Honour whichever binding-path mode produced `result` — otherwise the
  # one-pager would silently draw smoothed-NNLS events even when the main
  # run used coverage-aware scoring, and would disagree with result$binding_events.
  p_top <- plot_binding_deconvolution(
    tf_name, ld, pos_map, motif_hits,
    weight_mode      = result$weight_mode      %||% "mod_t",
    coverage_correct = result$coverage_correct %||% FALSE,
    cov_floor        = result$cov_floor        %||% 0.05,
    kernel_sigma     = result$kernel_sigma     %||% 250
  ) + labs(title = NULL, subtitle = NULL)

  # Middle panel: per-region bar of logFC with p-value stars
  df <- ld[ld$protein == tf_name & ld$region %in% names(pos_map), ]
  df$pos    <- as.numeric(pos_map[df$region])
  df <- df[order(df$pos), ]
  df$region <- factor(df$region, levels = df$region)
  df$stars  <- cut(df$pval,
                   breaks = c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf),
                   labels = c("****", "***", "**", "*", "ns"))
  lfc_range <- range(c(0, df$lfc), na.rm = TRUE)
  p_bars <- ggplot(df, aes(x = region, y = lfc)) +
    geom_col(aes(fill = pval <= 0.05), alpha = 0.8, width = 0.7) +
    geom_text(aes(y = ifelse(lfc >= 0, lfc + diff(lfc_range) * 0.02,
                                       lfc - diff(lfc_range) * 0.02),
                  label = stars,
                  vjust = ifelse(lfc >= 0, 0, 1)),
              size = 3.2, fontface = "bold") +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    scale_fill_manual(values = c(`TRUE` = COLS$high, `FALSE` = COLS$guide),
                      labels = c(`TRUE` = "p \u2264 0.05",
                                 `FALSE` = "p > 0.05"),
                      name = NULL) +
    labs(x = NULL, y = "logFC") +
    theme_caspex() +
    theme(legend.position = "bottom")

  # Bottom panel: stats block as text
  sp_row <- result$spatial_df[result$spatial_df$protein == tf_name, , drop = FALSE]
  if (nrow(sp_row) == 0) sp_row <- data.frame(centroid = NA, spread = NA,
                                              composite = NA, n_regions = NA)
  ev <- if (!is.null(result$binding_events))
    result$binding_events[result$binding_events$tf == tf_name, , drop = FALSE] else
      data.frame()
  pwm_id <- if (!is.null(result$motif_results[[tf_name]]))
    result$motif_results[[tf_name]]$pwm$id else "n/a"
  mode_lbl <- if (isTRUE(result$coverage_correct))
    paste0("coverage-aware (cov_floor=",
           result$cov_floor %||% 0.05, ")")
  else "default NNLS"
  stats_txt <- sprintf(
    "%s | centroid %s bp | spread %s bp | composite %s | n_regions %s\nmatrix %s | PWM hits %d | events %d | mode: %s",
    tf_name,
    format(sp_row$centroid[1], nsmall = 1),
    format(sp_row$spread[1], nsmall = 1),
    format(sp_row$composite[1], nsmall = 3),
    sp_row$n_regions[1],
    pwm_id, length(motif_hits), nrow(ev), mode_lbl
  )
  p_stats <- ggplot() +
    annotate("text", x = 0, y = 0.5, label = stats_txt,
             hjust = 0, vjust = 0.5, size = 3.4,
             family = "mono", color = COLS$neutral) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_void()

  (p_top / p_bars / p_stats) +
    plot_layout(heights = c(6, 3, 0.8)) +
    plot_annotation(
      title    = paste0(tf_name, " — CasPEX summary"),
      subtitle = paste0("Region enrichment · spatial footprint · binding-event calls"),
      theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                       plot.subtitle = element_text(color = "grey50"))
    )
}

# =============================================================================
# A.2  TF-family enrichment among the motif-scanned TFs
# =============================================================================

#' Barplot of TF-family representation among the top-ranked TFs.
#'
#' Queries JASPAR for each TF's "family" annotation (cached across calls
#' via a list in the returned plot's environment). If JASPAR family info
#' is unavailable, the bar is tagged "unknown".
#'
#' @param result The result object from run_caspex()
#' @param top_n  Number of top-composite TFs to include (default 30)
#' @return ggplot
plot_tf_family_enrichment <- function(result, top_n = 30) {
  spatial <- head(result$spatial_df, top_n)
  tfs <- as.character(spatial$protein)

  fam <- vapply(tfs, function(tf) {
    tryCatch({
      u <- paste0("https://jaspar.elixir.no/api/v1/matrix/?tf_name=", tf,
                  "&collection=CORE&format=json&page_size=1")
      r <- httr::GET(u, httr::timeout(10))
      if (httr::status_code(r) != 200) return("unknown")
      js <- httr::content(r, "parsed", simplifyVector = FALSE)
      if (length(js$results) == 0) return("unknown")
      f <- js$results[[1]]$family
      if (is.null(f) || length(f) == 0) return("unknown")
      paste(unlist(f), collapse = "/")
    }, error = function(e) "unknown")
  }, character(1))

  df <- data.frame(tf = tfs, family = fam, stringsAsFactors = FALSE)
  fam_tab <- as.data.frame(table(df$family))
  names(fam_tab) <- c("family", "n")
  fam_tab <- fam_tab[order(fam_tab$n, decreasing = TRUE), ]
  fam_tab$family <- factor(fam_tab$family, levels = fam_tab$family)

  ggplot(fam_tab, aes(x = family, y = n)) +
    geom_col(fill = COLS$guide, alpha = 0.85) +
    geom_text(aes(label = n), vjust = -0.3, size = 3, fontface = "bold",
              color = COLS$neutral) +
    labs(x = "JASPAR TF family",
         y = paste0("# TFs in top-", top_n),
         title = paste0("TF family representation among top-", top_n,
                        " composite TFs"),
         subtitle = paste0(nrow(df), " TFs queried | ",
                           sum(df$family == "unknown"),
                           " with no family annotation")) +
    theme_caspex() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# =============================================================================
# A.3  TSS-relative event density
# =============================================================================

#' Histogram of predicted binding-event positions across all TFs.
#'
#' @param result   Result object
#' @param binwidth Histogram binwidth in bp (default 50)
#' @param weighted Weight each event by its NNLS weight (default TRUE)
plot_event_density <- function(result, binwidth = 50, weighted = TRUE) {
  ev <- result$binding_events
  if (is.null(ev) || nrow(ev) == 0)
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "No binding events called"))
  ev$w <- if (weighted) ev$weight else 1

  ggplot(ev, aes(x = position, weight = w)) +
    geom_histogram(binwidth = binwidth, fill = COLS$guide,
                   color = "white", alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = COLS$tss, linewidth = 0.7) +
    annotate("text", x = 0, y = 0, label = " TSS",
             color = COLS$tss, hjust = 0, vjust = -0.5, size = 3.2,
             fontface = "bold") +
    scale_x_continuous(labels = scales::comma) +
    labs(x = "Position (bp, TSS-relative)",
         y = if (weighted) "Summed event weight" else "# events",
         title = "TSS-relative distribution of predicted binding events",
         subtitle = sprintf("%d events across %d TFs | binwidth %d bp | %sweighted",
                            nrow(ev), length(unique(ev$tf)), binwidth,
                            if (weighted) "" else "un")) +
    theme_caspex()
}

# =============================================================================
# A.4  Composite vs specificity scatter
# =============================================================================

#' 2-D scatter of composite enrichment vs per-TF specificity score.
#'
#' Specificity is defined per TF as max over regions of
#'   lfc_R − mean_{R' != R}(lfc_R').
plot_composite_vs_specificity <- function(result, top_label = 15) {
  ld <- result$long_data[result$long_data$protein %in% result$spatial_df$protein &
                            result$long_data$region %in% names(result$pos_map), ]

  spec <- vapply(unique(ld$protein), function(p) {
    sub <- ld[ld$protein == p, ]
    if (nrow(sub) < 2) return(NA_real_)
    mean_by_r <- tapply(sub$lfc, sub$region, mean, na.rm = TRUE)
    rs <- names(mean_by_r)
    vals <- vapply(rs, function(r) mean_by_r[r] - mean(mean_by_r[rs != r]),
                   numeric(1))
    max(vals, na.rm = TRUE)
  }, numeric(1))

  df <- data.frame(
    protein   = names(spec),
    specificity = as.numeric(spec),
    composite = result$spatial_df$composite[
      match(names(spec), result$spatial_df$protein)],
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$composite) & !is.na(df$specificity), ]
  df$score <- df$composite * pmax(df$specificity, 0)   # for labelling rank
  df <- df[order(df$score, decreasing = TRUE), ]
  df$label <- ifelse(seq_len(nrow(df)) <= top_label, df$protein, "")

  ggplot(df, aes(x = composite, y = specificity)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey60", linewidth = 0.3) +
    geom_point(aes(size = composite), color = COLS$guide,
               alpha = 0.55, shape = 16) +
    geom_text(aes(label = label), vjust = -0.9,
              size = 3, fontface = "bold", color = COLS$neutral) +
    scale_size_continuous(range = c(1.2, 7), guide = "none") +
    labs(x = "Composite enrichment (global strength)",
         y = "Max per-region specificity (lfc_R \u2212 mean lfc_other)",
         title = "Broadly strong vs region-specific TFs",
         subtitle = paste0("Top-right: strong + region-focal | ",
                           "Top-left: focal but weaker | ",
                           "Bottom-right: diffuse")) +
    theme_caspex()
}

# =============================================================================
# B.1  Permutation null for the composite score
# =============================================================================

#' Permutation null: shuffle region labels B times, recompute the spatial
#' model, and build an empirical null distribution of composite scores.
#'
#' Returns a list with the observed composite per TF and the per-TF
#' permutation p-value (= fraction of null scores >= observed).
run_permutation_null <- function(result, n_perm = 500, seed = 42,
                                  weight_mode = NULL) {
  set.seed(seed)
  weight_mode <- weight_mode %||% result$weight_mode %||% "mod_t"
  ld <- result$long_data[result$long_data$isTF, ]
  ld <- ld[ld$region %in% names(result$pos_map), ]
  pos <- result$pos_map
  proteins <- unique(ld$protein)

  observed <- setNames(result$spatial_df$composite,
                        as.character(result$spatial_df$protein))

  message("Running ", n_perm, " permutations across ", length(proteins), " TFs...")
  nulls <- matrix(NA_real_, nrow = length(proteins), ncol = n_perm,
                  dimnames = list(proteins, NULL))
  regs <- names(pos)

  for (b in seq_len(n_perm)) {
    if (b %% max(1, floor(n_perm / 10)) == 0)
      message("  permutation ", b, " / ", n_perm)
    perm_pos <- setNames(sample(pos), regs)
    per_tf <- vapply(proteins, function(p) {
      sub <- ld[ld$protein == p, ]
      sub$pos <- perm_pos[as.character(sub$region)]
      sub <- sub[!is.na(sub$pos) & !is.na(sub$lfc), ]
      if (nrow(sub) < 2) return(NA_real_)
      # Mirror compute_spatial(): floor lfc BEFORE weighting, then floor w
      sub$lfc <- pmax(sub$lfc, 0)
      w <- pmax(compute_region_weight(sub, mode = weight_mode), 0)
      if (sum(w) < 1e-8) return(0)
      mean(w) * log1p(nrow(sub))
    }, numeric(1))
    nulls[, b] <- per_tf
  }

  pvals <- vapply(proteins, function(p) {
    obs <- observed[[p]]
    if (is.null(obs) || is.na(obs)) return(NA_real_)
    v <- nulls[p, ]
    v <- v[!is.na(v)]
    if (length(v) == 0) return(NA_real_)
    (sum(v >= obs) + 1) / (length(v) + 1)
  }, numeric(1))

  out <- data.frame(
    protein   = names(pvals),
    composite = observed[names(pvals)],
    perm_p    = pvals,
    null_mean = rowMeans(nulls, na.rm = TRUE),
    null_sd   = apply(nulls, 1, sd, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  out$perm_fdr <- p.adjust(out$perm_p, method = "BH")
  out <- out[order(out$perm_p), ]
  list(summary = out, nulls = nulls, n_perm = n_perm, weight_mode = weight_mode)
}

#' Plot the permutation null: per-TF observed composite overlaid on its
#' null distribution. Shows top_n most-significant TFs.
plot_permutation_null <- function(perm_result, top_n = 20) {
  s <- head(perm_result$summary, top_n)
  s$protein <- factor(s$protein, levels = rev(s$protein))
  s$sig <- s$perm_fdr <= 0.05

  ggplot(s, aes(y = protein)) +
    geom_errorbarh(aes(xmin = null_mean - null_sd,
                       xmax = null_mean + null_sd),
                   color = "grey70", height = 0) +
    geom_point(aes(x = null_mean), color = "grey50", size = 2, shape = 4) +
    geom_point(aes(x = composite, fill = sig),
               shape = 21, size = 3.5, color = "black", stroke = 0.4) +
    scale_fill_manual(values = c(`TRUE` = COLS$high, `FALSE` = COLS$guide),
                      labels = c(`TRUE` = "FDR \u2264 0.05",
                                 `FALSE` = "FDR > 0.05"),
                      name = NULL) +
    geom_text(aes(x = composite,
                  label = sprintf("p=%.3g, q=%.3g", perm_p, perm_fdr)),
              hjust = -0.1, size = 2.8, color = COLS$neutral) +
    labs(x = "Composite score (observed vs null mean \u00b1 1 sd)",
         y = NULL,
         title = paste0("Permutation null (B=", perm_result$n_perm,
                        ") for composite score"),
         subtitle = paste0("Top ", top_n, " TFs by permutation p-value | ",
                           "weight mode: ", perm_result$weight_mode)) +
    theme_caspex() +
    theme(legend.position = "bottom",
          panel.grid.major.y = element_blank())
}

# =============================================================================
# B.2  Kernel sigma sensitivity
# =============================================================================

#' For each TF, rerun the deconvolution at multiple kernel widths and
#' capture the top event position per sigma. Stable centroids = robust.
#'
#' Auto-detects coverage-aware mode from `result$coverage_correct` and uses
#' the matching decoder at each σ, so the sensitivity grid reflects the
#' same mode that produced `result$binding_events`.
run_sigma_sensitivity <- function(result,
                                   sigmas = c(100, 200, 250, 300, 500),
                                   tfs    = NULL,
                                   weight_mode      = NULL,
                                   coverage_correct = NULL,
                                   cov_floor        = NULL,
                                   edge_guard_frac  = NULL) {
  weight_mode      <- weight_mode      %||% result$weight_mode      %||% "mod_t"
  coverage_correct <- coverage_correct %||% result$coverage_correct %||% FALSE
  cov_floor        <- cov_floor        %||% result$cov_floor        %||% 0.05
  edge_guard_frac  <- edge_guard_frac  %||% result$edge_guard_frac  %||% 0.15
  if (is.null(tfs)) {
    # use top-15 TFs in motif_results
    tfs <- intersect(head(as.character(result$spatial_df$protein), 15),
                     names(result$motif_results))
  }
  if (length(tfs) == 0) return(invisible(NULL))

  message("Sigma sensitivity across ", length(sigmas),
          " kernels \u00d7 ", length(tfs), " TFs",
          if (coverage_correct) "  [coverage-aware]" else "  [default NNLS]",
          "...")
  rows <- list()
  for (tf in tfs) {
    hits <- result$motif_results[[tf]]$hits
    if (is.null(hits)) hits <- integer(0)
    for (sg in sigmas) {
      ev <- if (coverage_correct) {
        predict_binding_events_coverage_aware(
          tf, result$long_data, result$pos_map, hits,
          kernel_sigma = sg, weight_mode = weight_mode,
          cov_floor = cov_floor,
          edge_guard_frac = edge_guard_frac)
      } else {
        predict_binding_events(tf, result$long_data, result$pos_map,
                               hits, kernel_sigma = sg,
                               weight_mode = weight_mode)
      }
      if (nrow(ev) == 0) next
      ev <- ev[order(ev$weight, decreasing = TRUE), ]
      top_n <- min(3, nrow(ev))
      rows[[length(rows) + 1]] <- data.frame(
        tf = tf, sigma = sg,
        rank = seq_len(top_n),
        position = ev$position[seq_len(top_n)],
        weight = ev$weight[seq_len(top_n)],
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) return(invisible(NULL))
  do.call(rbind, rows)
}

plot_sigma_sensitivity <- function(sigma_result) {
  if (is.null(sigma_result) || nrow(sigma_result) == 0)
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "No sigma-sensitivity data"))

  ggplot(sigma_result, aes(x = sigma, y = position,
                            group = interaction(tf, rank),
                            color = factor(rank))) +
    geom_line(alpha = 0.5) +
    geom_point(aes(size = weight), alpha = 0.85) +
    facet_wrap(~ tf, ncol = 3, scales = "free_y") +
    scale_color_manual(values = c(`1` = COLS$high, `2` = COLS$mid,
                                   `3` = COLS$guide),
                       name = "event rank") +
    scale_size_continuous(range = c(1.5, 5), guide = "none") +
    labs(x = "Kernel \u03c3 (bp)",
         y = "Event position (bp, TSS-relative)",
         title = "Kernel \u03c3 sensitivity of deconvolved binding events",
         subtitle = "Flat lines = robust call | large drift = artefact of smoothing width") +
    theme_caspex() +
    theme(legend.position = "bottom")
}

# =============================================================================
# B.3  Event jackknife (leave-one-region-out)
# =============================================================================

#' Drop each region in turn, rerun the deconvolution, and count how often
#' each event position survives within `tol_bp` bp.
#'
#' Auto-detects coverage-aware mode from `result$coverage_correct` so the
#' leave-one-out re-fits use the same binding path as the original call.
run_event_jackknife <- function(result, tfs = NULL, tol_bp = 100,
                                 weight_mode      = NULL,
                                 coverage_correct = NULL,
                                 cov_floor        = NULL,
                                 edge_guard_frac  = NULL) {
  weight_mode      <- weight_mode      %||% result$weight_mode      %||% "mod_t"
  coverage_correct <- coverage_correct %||% result$coverage_correct %||% FALSE
  cov_floor        <- cov_floor        %||% result$cov_floor        %||% 0.05
  edge_guard_frac  <- edge_guard_frac  %||% result$edge_guard_frac  %||% 0.15
  if (is.null(tfs))
    tfs <- unique(result$binding_events$tf)
  if (length(tfs) == 0) return(invisible(NULL))

  rs <- names(result$pos_map)
  message("Jackknife across ", length(rs), " regions \u00d7 ",
          length(tfs), " TFs",
          if (coverage_correct) "  [coverage-aware]" else "  [default NNLS]",
          "...")

  rows <- list()
  for (tf in tfs) {
    obs <- result$binding_events[result$binding_events$tf == tf, ]
    hits <- if (!is.null(result$motif_results[[tf]]))
      result$motif_results[[tf]]$hits else integer(0)

    surv <- integer(nrow(obs))
    for (r in rs) {
      pos_sub <- result$pos_map[setdiff(names(result$pos_map), r)]
      ev <- if (coverage_correct) {
        predict_binding_events_coverage_aware(
          tf, result$long_data, pos_sub, hits,
          weight_mode = weight_mode, cov_floor = cov_floor,
          edge_guard_frac = edge_guard_frac)
      } else {
        predict_binding_events(tf, result$long_data, pos_sub,
                               hits, weight_mode = weight_mode)
      }
      if (nrow(ev) == 0) next
      for (j in seq_len(nrow(obs))) {
        if (any(abs(ev$position - obs$position[j]) <= tol_bp))
          surv[j] <- surv[j] + 1L
      }
    }
    rows[[length(rows) + 1]] <- data.frame(
      tf = tf,
      position = obs$position,
      weight   = obs$weight,
      n_surv   = surv,
      n_drops  = length(rs),
      frac     = surv / length(rs),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(invisible(NULL))
  do.call(rbind, rows)
}

plot_event_jackknife <- function(jk_result, top_n = 40) {
  if (is.null(jk_result) || nrow(jk_result) == 0)
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "No jackknife data"))
  jk <- jk_result[order(jk_result$frac, jk_result$weight,
                         decreasing = TRUE), ]
  jk <- head(jk, top_n)
  jk$label <- sprintf("%s @ %+.0f bp", jk$tf, jk$position)
  jk$label <- factor(jk$label, levels = rev(jk$label))

  ggplot(jk, aes(x = frac, y = label, fill = frac)) +
    geom_col(alpha = 0.9) +
    geom_text(aes(label = sprintf("%d/%d", n_surv, n_drops)),
              hjust = -0.1, size = 2.6, color = COLS$neutral) +
    scale_x_continuous(limits = c(0, 1.15),
                       breaks = seq(0, 1, 0.25),
                       labels = scales::percent) +
    scale_fill_gradient(low = COLS$low, high = COLS$high, guide = "none") +
    labs(x = "Fraction of jackknife replicates in which the event survived",
         y = NULL,
         title = "Event robustness to leave-one-region-out",
         subtitle = paste0("Top-", top_n,
                           " events by survival fraction | tolerance: \u00b1100 bp")) +
    theme_caspex() +
    theme(panel.grid.major.y = element_blank())
}

# =============================================================================
# B.4  NNLS residual plot for a TF
# =============================================================================

#' Overlay the observed signal s(x), the NNLS reconstruction X beta, and
#' the residual for a single TF. A large residual = JASPAR motifs for
#' this TF cannot fully account for the observed CasPEX signal.
plot_nnls_residual <- function(result, tf_name, kernel_sigma = 250,
                                weight_mode = NULL) {
  weight_mode <- weight_mode %||% result$weight_mode %||% "mod_t"
  if (!tf_name %in% result$long_data$protein)
    stop(tf_name, " not in long_data")

  x_grid <- seq(-2500, 500, by = 5)
  sig <- build_caspex_signal(tf_name, result$long_data, result$pos_map,
                              x_grid, kernel_sigma, weight_mode)
  hits <- if (!is.null(result$motif_results[[tf_name]]))
    result$motif_results[[tf_name]]$hits else integer(0)
  hits <- hits[!is.na(hits) & hits >= min(x_grid) & hits <= max(x_grid)]

  reconstruction <- numeric(length(x_grid))
  if (length(hits) > 0 && requireNamespace("nnls", quietly = TRUE)) {
    X <- vapply(hits,
                function(m) exp(-0.5 * ((x_grid - m) / kernel_sigma)^2),
                numeric(length(x_grid)))
    if (length(hits) == 1) X <- matrix(X, ncol = 1)
    fit <- nnls::nnls(X, sig$y)
    reconstruction <- as.numeric(X %*% fit$x)
  }
  residual <- sig$y - reconstruction

  df <- rbind(
    data.frame(x = x_grid, y = sig$y,         lab = "Observed s(x)"),
    data.frame(x = x_grid, y = reconstruction, lab = "NNLS reconstruction X\u03b2"),
    data.frame(x = x_grid, y = residual,      lab = "Residual")
  )
  df$lab <- factor(df$lab,
                   levels = c("Observed s(x)",
                              "NNLS reconstruction X\u03b2",
                              "Residual"))
  ss_obs <- sum(sig$y^2)
  ss_res <- sum(residual^2)
  r2 <- if (ss_obs > 0) 1 - ss_res / ss_obs else NA_real_

  ggplot(df, aes(x = x, y = y)) +
    geom_area(aes(fill = lab), alpha = 0.4, color = NA) +
    geom_line(aes(color = lab), linewidth = 0.5) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = COLS$tss, linewidth = 0.6) +
    facet_wrap(~ lab, ncol = 1, scales = "free_y") +
    scale_color_manual(values = c("Observed s(x)"                = COLS$guide,
                                   "NNLS reconstruction X\u03b2" = COLS$low,
                                   "Residual"                    = COLS$high),
                       guide = "none") +
    scale_fill_manual(values = c("Observed s(x)"                = COLS$guide,
                                  "NNLS reconstruction X\u03b2" = COLS$low,
                                  "Residual"                    = COLS$high),
                      guide = "none") +
    scale_x_continuous(labels = scales::comma) +
    labs(x = "Position (bp, TSS-relative)", y = "Signal (a.u.)",
         title = paste0(tf_name, " \u2014 NNLS fit quality"),
         subtitle = sprintf("%d motif candidates | fit R\u00b2 = %.3f | \u03c3 = %d bp",
                            length(hits), r2, kernel_sigma)) +
    theme_caspex() +
    theme(strip.text = element_text(face = "bold"))
}

# =============================================================================
# C.1  Volcano plots per region
# =============================================================================

#' One volcano per region, TFs highlighted.
plot_volcano_per_region <- function(result, pval_thresh = 0.05,
                                     lfc_thresh = 1, label_top = 8) {
  ld <- result$long_data
  ld$logp <- -log10(pmax(ld$pval, 1e-16))
  ld$sig  <- ld$pval <= pval_thresh & abs(ld$lfc) >= lfc_thresh
  ld$kind <- ifelse(!ld$sig, "ns",
                    ifelse(ld$isTF, "TF (sig)", "non-TF (sig)"))
  ld$kind <- factor(ld$kind, levels = c("ns", "non-TF (sig)", "TF (sig)"))

  # Label top-N TFs per region by -log10(p)
  ld$label <- ""
  for (r in unique(ld$region)) {
    sub <- ld[ld$region == r & ld$isTF & ld$sig, ]
    sub <- sub[order(sub$logp, decreasing = TRUE), ]
    top <- head(sub$protein, label_top)
    ld$label[ld$region == r & ld$protein %in% top & ld$isTF & ld$sig] <-
      as.character(ld$protein[ld$region == r & ld$protein %in% top &
                                ld$isTF & ld$sig])
  }

  ggplot(ld, aes(x = lfc, y = logp)) +
    geom_point(aes(color = kind, size = kind, alpha = kind)) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh),
               linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_hline(yintercept = -log10(pval_thresh),
               linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_text(aes(label = label), vjust = -0.8, size = 2.6,
              color = COLS$neutral, fontface = "bold",
              check_overlap = TRUE) +
    facet_wrap(~ region, scales = "free") +
    scale_color_manual(values = c(ns = "grey80",
                                   `non-TF (sig)` = COLS$mid,
                                   `TF (sig)`     = COLS$high),
                       name = NULL) +
    scale_size_manual(values = c(ns = 0.6, `non-TF (sig)` = 1.1,
                                  `TF (sig)` = 1.6),
                      guide = "none") +
    scale_alpha_manual(values = c(ns = 0.3, `non-TF (sig)` = 0.75,
                                   `TF (sig)` = 0.95),
                      guide = "none") +
    labs(x = "logFC", y = "-log10(p)",
         title = "Per-region volcano plots",
         subtitle = sprintf("dashed lines: |logFC| \u2265 %g and p \u2264 %g",
                            lfc_thresh, pval_thresh)) +
    theme_caspex() +
    theme(legend.position = "bottom")
}

# =============================================================================
# C.2  Region correlation heatmap
# =============================================================================

#' Pearson correlation heatmap of region logFC profiles.
plot_region_correlation <- function(result, method = "pearson") {
  ld <- result$long_data
  wide <- reshape(ld[, c("protein", "region", "lfc")],
                  idvar = "protein", timevar = "region",
                  direction = "wide")
  mat <- as.matrix(wide[, -1])
  colnames(mat) <- sub("^lfc\\.", "", colnames(mat))
  mat <- mat[complete.cases(mat), , drop = FALSE]
  cor_mat <- cor(mat, method = method)

  # To data.frame
  df <- as.data.frame(as.table(cor_mat))
  names(df) <- c("R1", "R2", "r")

  ggplot(df, aes(x = R1, y = R2, fill = r)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", r)),
              size = 3, color = "grey20") +
    scale_fill_gradient2(low = COLS$low, mid = "white", high = COLS$high,
                         midpoint = 0, limits = c(-1, 1),
                         name = method) +
    labs(x = NULL, y = NULL,
         title = "Inter-region logFC correlation",
         subtitle = paste0(nrow(mat),
                           " proteins with complete observations | ",
                           method)) +
    coord_fixed() +
    theme_caspex() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
}

# =============================================================================
# C.3  P-value distribution per region
# =============================================================================

#' Faceted p-value histogram — should be roughly uniform with a peak near 0.
plot_pval_histograms <- function(result, binwidth = 0.02) {
  ld <- result$long_data
  ggplot(ld, aes(x = pval)) +
    geom_histogram(binwidth = binwidth, fill = COLS$guide, color = "white",
                   alpha = 0.85) +
    geom_vline(xintercept = 0.05, linetype = "dashed",
               color = COLS$high, linewidth = 0.4) +
    facet_wrap(~ region) +
    labs(x = "p-value", y = "# proteins",
         title = "Per-region p-value distribution (QC)",
         subtitle = "Expected: uniform over (0.05, 1], with excess near 0") +
    theme_caspex()
}

# =============================================================================
# C.4  Motif-strength vs NNLS-weight
# =============================================================================

#' Scatter of PWM log-odds score against NNLS beta per motif hit,
#' faceted by TF. Reinforces: bubble size is NOT PWM score.
plot_motif_vs_nnls <- function(result, tfs = NULL,
                                kernel_sigma = 250,
                                weight_mode = NULL) {
  weight_mode <- weight_mode %||% result$weight_mode %||% "mod_t"
  if (is.null(tfs))
    tfs <- intersect(unique(result$binding_events$tf),
                      names(result$motif_results))
  if (length(tfs) == 0)
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "No motif/event data"))

  x_grid <- seq(-2500, 500, by = 5)
  rows <- list()
  for (tf in tfs) {
    mr <- result$motif_results[[tf]]
    if (is.null(mr)) next
    hits <- mr$hits
    if (length(hits) == 0) next

    # Recompute PWM score per hit on the promoter sequence
    pwm <- mr$pwm$pwm; L <- mr$pwm$len
    if (is.null(pwm) || is.null(L)) next
    seq_chars <- strsplit(result$promoter_info$seq, "")[[1]]
    rev_map   <- c(A = "T", T = "A", G = "C", C = "G", N = "N")
    rev_chars <- rev(rev_map[seq_chars])
    tss_i <- result$promoter_info$tss_offset + 1
    # fwd scores indexed by 1-based window start
    fwd_scores <- score_pwm_positions(seq_chars, pwm)
    rev_scores <- score_pwm_positions(rev_chars, pwm)
    n <- length(seq_chars)

    fwd_pos_to_tssrel <- function(i) i - tss_i
    rev_pos_to_tssrel <- function(i) n - (i + L - 2) - tss_i

    # For each hit position (TSS-relative), get the best PWM score across strands
    pwm_score <- vapply(hits, function(h) {
      fi <- which(abs(fwd_pos_to_tssrel(seq_along(fwd_scores)) - h) <= 2)
      ri <- which(abs(rev_pos_to_tssrel(seq_along(rev_scores)) - h) <= 2)
      cands <- c(fwd_scores[fi], rev_scores[ri])
      if (length(cands) == 0) return(NA_real_)
      max(cands)
    }, numeric(1))

    # NNLS β per hit
    if (!requireNamespace("nnls", quietly = TRUE)) next
    sig <- build_caspex_signal(tf, result$long_data, result$pos_map,
                                x_grid, kernel_sigma, weight_mode)
    valid <- hits >= min(x_grid) & hits <= max(x_grid)
    hits_v <- hits[valid]; pwm_score_v <- pwm_score[valid]
    if (length(hits_v) == 0) next
    X <- vapply(hits_v,
                function(m) exp(-0.5 * ((x_grid - m) / kernel_sigma)^2),
                numeric(length(x_grid)))
    if (length(hits_v) == 1) X <- matrix(X, ncol = 1)
    fit <- nnls::nnls(X, sig$y)

    rows[[length(rows) + 1]] <- data.frame(
      tf = tf, hit_pos = hits_v,
      pwm_score = pwm_score_v,
      beta = fit$x,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0)
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "No motif/beta pairs"))
  df <- do.call(rbind, rows)
  df$survived <- df$beta > 0

  ggplot(df, aes(x = pwm_score, y = beta)) +
    geom_point(aes(color = survived), size = 2, alpha = 0.8) +
    facet_wrap(~ tf, scales = "free") +
    scale_color_manual(values = c(`TRUE` = COLS$high, `FALSE` = "grey70"),
                       labels = c(`TRUE` = "\u03b2 > 0 (kept)",
                                  `FALSE` = "\u03b2 = 0"),
                       name = NULL) +
    labs(x = "JASPAR PWM log-odds score",
         y = "NNLS coefficient \u03b2",
         title = "Motif match strength vs NNLS weight",
         subtitle = "Same PWM score, different \u03b2 \u2192 NNLS is reflecting the signal, not the PWM") +
    theme_caspex() +
    theme(legend.position = "bottom")
}

# =============================================================================
# A.5  TF motif co-occurrence matrix
# =============================================================================

#' Heat-map of TF-pair co-occurrence of predicted binding events.
#'
#' For every pair of TFs with called events, count how many events of one
#' TF fall within `tol_bp` of an event of the other TF. The resulting
#' symmetric matrix (proportion-scaled) is a candidate-cofactor screen:
#' pairs lighting up have overlapping binding footprints on this promoter.
#'
#' Mode-agnostic: works identically for default and coverage-aware events
#' because it only consumes the `tf` / `position` columns of
#' `result$binding_events`.
plot_tf_cooccurrence <- function(result, tol_bp = 100, min_events = 1,
                                  top_tfs = 30) {
  ev <- result$binding_events
  if (is.null(ev) || nrow(ev) == 0)
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "No events"))
  # Restrict to TFs with >= min_events and (if too many) the top-N by
  # total weight so the heat-map stays readable.
  tot_w <- tapply(ev$weight, ev$tf, sum)
  n_ev  <- table(ev$tf)
  keep  <- names(n_ev)[n_ev >= min_events]
  tot_w <- tot_w[keep]
  keep  <- names(head(sort(tot_w, decreasing = TRUE), top_tfs))
  ev    <- ev[ev$tf %in% keep, ]
  tfs   <- sort(unique(as.character(ev$tf)))
  n     <- length(tfs)
  if (n < 2)
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "Need \u2265 2 TFs with events"))

  M <- matrix(0L, n, n, dimnames = list(tfs, tfs))
  for (i in seq_len(n)) {
    pi_ <- ev$position[ev$tf == tfs[i]]
    for (j in seq_len(n)) {
      if (i == j) next
      pj_ <- ev$position[ev$tf == tfs[j]]
      # Count i-events whose nearest j-event is within tol_bp
      M[i, j] <- sum(vapply(pi_, function(p)
        any(abs(pj_ - p) <= tol_bp), logical(1)))
    }
  }
  # Symmetric proportion of i-events co-located with any j-event
  denom <- pmax(as.integer(n_ev[tfs]), 1L)
  frac  <- sweep(M, 1, denom, FUN = "/")
  # Symmetrize by max so the heat-map reads consistently in both directions
  frac_sym <- pmax(frac, t(frac))
  diag(frac_sym) <- NA

  df <- as.data.frame(as.table(frac_sym))
  names(df) <- c("tf_i", "tf_j", "frac")

  ggplot(df, aes(x = tf_i, y = tf_j, fill = frac)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient(low = "white", high = COLS$high, na.value = "grey95",
                        labels = scales::percent, name = "co-loc",
                        limits = c(0, 1)) +
    labs(x = NULL, y = NULL,
         title = "TF co-occurrence of binding events",
         subtitle = paste0("% of events within \u00b1", tol_bp,
                           " bp of another TF's event | top-", top_tfs,
                           " TFs by summed event weight")) +
    coord_fixed() +
    theme_caspex() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid  = element_blank())
}

# =============================================================================
# A.6  Ranked event table (confidence-scored)
# =============================================================================

#' Rank every called binding event by a composite confidence score and
#' emit a barplot + CSV of the top-N events.
#'
#' Score:   conf = z_beta + z_surv
#'   z_beta = per-TF-normalized β   (how strong is this event for this TF)
#'   z_surv = jackknife survival fraction (B.3)           [default 1 if absent]
#'
#' The jackknife result is accepted as an optional second argument. If
#' provided, it will be merged by (tf, position) within ±tol_bp so events
#' inherit their survival fraction; otherwise survival defaults to 1 and
#' the ranking reduces to per-TF-normalized β.
rank_binding_events <- function(result, jk_result = NULL, tol_bp = 100,
                                 top_n = 50) {
  ev <- result$binding_events
  if (is.null(ev) || nrow(ev) == 0)
    return(invisible(NULL))
  out <- ev
  # Per-TF max-normalized β: places a TF's strongest event at 1
  tf_max <- tapply(out$weight, out$tf, max, na.rm = TRUE)
  out$beta_norm <- as.numeric(out$weight / pmax(tf_max[as.character(out$tf)], 1e-9))
  # Attach jackknife survival if we have it
  if (!is.null(jk_result) && nrow(jk_result) > 0) {
    out$surv <- vapply(seq_len(nrow(out)), function(i) {
      m <- jk_result$tf == out$tf[i] &
           abs(jk_result$position - out$position[i]) <= tol_bp
      if (!any(m)) return(NA_real_)
      mean(jk_result$frac[m], na.rm = TRUE)
    }, numeric(1))
  } else {
    out$surv <- NA_real_
  }
  # Final confidence. Treat missing survival as 1.0 so events without
  # jackknife support are not penalised relative to ones that also lack it.
  surv_fill <- ifelse(is.na(out$surv), 1, out$surv)
  out$confidence <- out$beta_norm + surv_fill
  out <- out[order(out$confidence, decreasing = TRUE), , drop = FALSE]
  top <- head(out, top_n)
  top$label <- sprintf("%s @ %+.0f", top$tf, top$position)
  top$label <- factor(top$label, levels = rev(top$label))

  p <- ggplot(top, aes(x = confidence, y = label, fill = beta_norm)) +
    geom_col(alpha = 0.9) +
    geom_text(aes(label = sprintf("\u03b2* %.2f | surv %s",
                                   beta_norm,
                                   ifelse(is.na(surv), "\u2014",
                                          sprintf("%.2f", surv)))),
              hjust = -0.05, size = 2.6, color = COLS$neutral) +
    scale_fill_gradient(low = COLS$guide, high = COLS$high, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(x = "Confidence = per-TF-normalized \u03b2 + jackknife survival",
         y = NULL,
         title = paste0("Top-", top_n, " ranked binding events"),
         subtitle = paste0("Higher = stronger + more reproducible | ",
                           "mode: ",
                           if (isTRUE(result$coverage_correct))
                             "coverage-aware" else "default NNLS")) +
    theme_caspex() +
    theme(panel.grid.major.y = element_blank())
  list(ranked = out, plot = p)
}

# =============================================================================
# D.1  Coverage-rescue audit scatter (coverage-aware mode only)
# =============================================================================

#' β vs local_coverage scatter coloured by distance-to-nearest-gRNA.
#'
#' Each point is one called event from a coverage-aware run. Rescued calls
#' — the ones that depend on dividing by a small `C(x)` to survive — land
#' in the top-left (high β, low local_coverage) and are coloured by how
#' far they sit from the nearest gRNA cut site. This is the single most
#' direct sanity check for the coverage correction: it makes "this call
#' only exists because cov_floor clamped the denominator" visible in one
#' view instead of requiring a cross-reference against C(x).
plot_coverage_rescue_scatter <- function(result, top_label = 15) {
  if (!isTRUE(result$coverage_correct))
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "D1 is coverage-aware-only"))
  ev <- result$binding_events
  need <- c("local_coverage", "distance_to_nearest_grna")
  if (is.null(ev) || nrow(ev) == 0 || !all(need %in% names(ev)))
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "No coverage-aware event table"))

  ev$rescued <- ev$local_coverage <= (result$cov_floor %||% 0.05) * 1.5
  # rank by β descending to pick top labels
  ord <- order(ev$weight, decreasing = TRUE)
  ev$label <- ""
  ev$label[ord[seq_len(min(top_label, nrow(ev)))]] <-
    sprintf("%s @ %+.0f",
            ev$tf[ord[seq_len(min(top_label, nrow(ev)))]],
            ev$position[ord[seq_len(min(top_label, nrow(ev)))]])

  ggplot(ev, aes(x = local_coverage, y = weight)) +
    geom_vline(xintercept = (result$cov_floor %||% 0.05),
               linetype = "dashed", color = COLS$high, linewidth = 0.3) +
    annotate("text",
             x = (result$cov_floor %||% 0.05),
             y = Inf,
             label = " cov_floor",
             color = COLS$high, hjust = 0, vjust = 1.4, size = 3) +
    geom_point(aes(color = distance_to_nearest_grna,
                   size  = weight,
                   shape = rescued),
               alpha = 0.85, stroke = 0.3) +
    geom_text(aes(label = label), vjust = -0.9,
              size = 2.6, fontface = "bold", color = COLS$neutral,
              check_overlap = TRUE) +
    scale_color_gradient(low = COLS$guide, high = COLS$high,
                         name = "dist to gRNA (bp)") +
    scale_size_continuous(range = c(1.2, 6), guide = "none") +
    scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16),
                       labels = c(`TRUE` = "near floor (rescued)",
                                  `FALSE` = "above floor"),
                       name = NULL) +
    labs(x = "local coverage C(event_pos)",
         y = "event weight \u03b2",
         title = "Coverage-rescue audit",
         subtitle = paste0("cov_floor = ", result$cov_floor %||% 0.05,
                           " | triangles sit near the floor and were ",
                           "amplified by the s/C correction")) +
    theme_caspex()
}

# =============================================================================
# D.2  cov_floor sensitivity sweep (coverage-aware mode only)
# =============================================================================

#' Re-score events at multiple cov_floor values and track event stability.
#'
#' Analog of B.2 (sigma sensitivity), but sweeping the coverage-floor
#' parameter that controls how aggressively distal / gap binders are
#' rescued. Lower floor = more distal rescues but more tail noise.
#' Robust calls drift little across the sweep; calls that appear only at
#' cov_floor ≤ 0.02 (say) should be treated as floor-sensitive.
run_covfloor_sensitivity <- function(result,
                                      floors = c(0.02, 0.05, 0.10, 0.20),
                                      tfs    = NULL,
                                      weight_mode = NULL,
                                      edge_guard_frac = NULL) {
  if (!isTRUE(result$coverage_correct)) {
    message("cov_floor sensitivity skipped: result was not coverage-aware.")
    return(invisible(NULL))
  }
  weight_mode     <- weight_mode     %||% result$weight_mode     %||% "mod_t"
  edge_guard_frac <- edge_guard_frac %||% result$edge_guard_frac %||% 0.15
  if (is.null(tfs))
    tfs <- intersect(unique(result$binding_events$tf),
                      names(result$motif_results))
  if (length(tfs) == 0) return(invisible(NULL))

  message("cov_floor sensitivity: ", length(floors), " floors \u00d7 ",
          length(tfs), " TFs...")
  rows <- list()
  for (tf in tfs) {
    hits <- result$motif_results[[tf]]$hits
    if (is.null(hits)) hits <- integer(0)
    for (fl in floors) {
      ev <- predict_binding_events_coverage_aware(
        tf, result$long_data, result$pos_map, hits,
        weight_mode = weight_mode, cov_floor = fl,
        edge_guard_frac = edge_guard_frac)
      if (nrow(ev) == 0) next
      ev <- ev[order(ev$weight, decreasing = TRUE), , drop = FALSE]
      top <- head(ev, 3)
      rows[[length(rows) + 1]] <- data.frame(
        tf = tf, cov_floor = fl,
        rank = seq_len(nrow(top)),
        position = top$position,
        weight   = top$weight,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) return(invisible(NULL))
  do.call(rbind, rows)
}

plot_covfloor_sensitivity <- function(cf_result) {
  if (is.null(cf_result) || nrow(cf_result) == 0)
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "No cov_floor-sensitivity data"))
  ggplot(cf_result, aes(x = cov_floor, y = position,
                         group = interaction(tf, rank),
                         color = factor(rank))) +
    geom_line(alpha = 0.5) +
    geom_point(aes(size = weight), alpha = 0.85) +
    facet_wrap(~ tf, ncol = 3, scales = "free_y") +
    scale_color_manual(values = c(`1` = COLS$high, `2` = COLS$mid,
                                   `3` = COLS$guide),
                       name = "event rank") +
    scale_size_continuous(range = c(1.5, 5), guide = "none") +
    scale_x_log10() +
    labs(x = "cov_floor (log scale)",
         y = "event position (bp, TSS-relative)",
         title = "cov_floor sensitivity of coverage-aware events",
         subtitle = "Flat lines = robust call | drift to the left = floor-sensitive rescue") +
    theme_caspex() +
    theme(legend.position = "bottom")
}

# =============================================================================
# D.3  s(x) / C(x) / β(x) per-TF stack (coverage-aware mode only)
# =============================================================================

#' Three-panel per-TF diagnostic: signal s(x), coverage C(x), ratio β(x).
#'
#' Analog of B.4 (NNLS residual) for coverage-aware mode. Stacking the
#' three panels with called events overlaid on the β panel shows exactly
#' what the s/C correction is doing at each event: where labeling
#' opportunity was scarce versus where β is simply tracking strong
#' enrichment.
plot_coverage_stack <- function(result, tf_name, kernel_sigma = NULL,
                                 weight_mode = NULL) {
  if (!isTRUE(result$coverage_correct))
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "D3 is coverage-aware-only"))
  weight_mode  <- weight_mode  %||% result$weight_mode  %||% "mod_t"
  kernel_sigma <- kernel_sigma %||% result$kernel_sigma %||% 250

  if (!tf_name %in% result$long_data$protein)
    stop(tf_name, " not in long_data")

  x_grid <- seq(-2500, 500, by = 5)
  sig <- build_caspex_signal(tf_name, result$long_data, result$pos_map,
                              x_grid, kernel_sigma, weight_mode)
  # C(x) using the same helper the pipeline uses internally
  cov <- if (exists("compute_coverage")) {
    compute_coverage(result$pos_map, x_grid, kernel_sigma)
  } else {
    # Fallback: recompute inline
    Reduce("+", lapply(result$pos_map,
      function(p) exp(-0.5 * ((x_grid - p) / kernel_sigma)^2)))
  }
  floor_val <- (result$cov_floor %||% 0.05) * max(cov)
  beta_curve <- sig$y / pmax(cov, floor_val)

  ev <- result$binding_events[result$binding_events$tf == tf_name, ,
                              drop = FALSE]

  df <- rbind(
    data.frame(x = x_grid, y = sig$y,       panel = "s(x) — signal"),
    data.frame(x = x_grid, y = cov,         panel = "C(x) — coverage"),
    data.frame(x = x_grid, y = beta_curve,  panel = "\u03b2(x) = s/C")
  )
  df$panel <- factor(df$panel,
                     levels = c("s(x) — signal",
                                "C(x) — coverage",
                                "\u03b2(x) = s/C"))
  ev_df <- if (nrow(ev) > 0)
    data.frame(x = ev$position, weight = ev$weight,
               panel = factor("\u03b2(x) = s/C",
                              levels = levels(df$panel)))
  else NULL

  p <- ggplot(df, aes(x = x, y = y)) +
    geom_area(aes(fill = panel), alpha = 0.35, color = NA) +
    geom_line(aes(color = panel), linewidth = 0.5) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = COLS$tss, linewidth = 0.6) +
    facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    scale_color_manual(values = c("s(x) — signal"     = COLS$guide,
                                   "C(x) — coverage"  = COLS$low,
                                   "\u03b2(x) = s/C"  = COLS$high),
                       guide = "none") +
    scale_fill_manual(values = c("s(x) — signal"      = COLS$guide,
                                  "C(x) — coverage"   = COLS$low,
                                  "\u03b2(x) = s/C"   = COLS$high),
                      guide = "none") +
    scale_x_continuous(labels = scales::comma) +
    labs(x = "position (bp, TSS-relative)", y = NULL,
         title = paste0(tf_name,
                        " \u2014 coverage-aware decomposition"),
         subtitle = sprintf("\u03c3 = %d bp | cov_floor = %g | %d event(s)",
                             kernel_sigma, result$cov_floor %||% 0.05,
                             nrow(ev))) +
    theme_caspex() +
    theme(strip.text = element_text(face = "bold"))
  if (!is.null(ev_df))
    p <- p + geom_point(data = ev_df,
                        aes(x = x, y = 0, size = weight),
                        inherit.aes = FALSE,
                        shape = 21, fill = COLS$high, color = "black",
                        stroke = 0.4, alpha = 0.9) +
      scale_size_continuous(range = c(2, 6), guide = "none")
  p
}

# =============================================================================
# E.  Mode comparison (default vs coverage-aware)
# =============================================================================

#' Pair-up a default-mode and a coverage-aware result object, match events
#' by (tf, rounded position), and emit a per-TF stacked count barplot plus
#' a position-jitter plot and a CSV of mode-only calls.
#'
#' Intended to be called from a driver that has already produced two
#' `run_caspex()` results — see `1-try-compare.R`.
run_caspex_extras_compare <- function(res_default, res_coverage,
                                      out_dir = "caspex_output_compare/extras",
                                      bin_bp  = 50) {
  if (missing(res_default) || missing(res_coverage))
    stop("run_caspex_extras_compare() needs both result objects. ",
         "Pass the output of run_caspex(..., coverage_correct=FALSE) ",
         "and run_caspex(..., coverage_correct=TRUE).")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir_abs <- normalizePath(out_dir, mustWork = FALSE)
  message("\n=== Mode-comparison extras: writing to ", out_dir_abs, " ===")

  ev_d <- res_default$binding_events
  ev_c <- res_coverage$binding_events
  bin_key <- function(df) paste0(df$tf, "|",
                                 round(df$position / bin_bp) * bin_bp)
  if (nrow(ev_d) > 0) ev_d$key <- bin_key(ev_d) else ev_d$key <- character(0)
  if (nrow(ev_c) > 0) ev_c$key <- bin_key(ev_c) else ev_c$key <- character(0)

  # Per-TF count table across the three provenance classes
  tfs <- sort(unique(c(as.character(ev_d$tf), as.character(ev_c$tf))))
  cnt <- data.frame(tf = tfs,
                     default_only  = 0L,
                     coverage_only = 0L,
                     both          = 0L,
                     stringsAsFactors = FALSE)
  for (i in seq_along(tfs)) {
    tf <- tfs[i]
    kd <- ev_d$key[ev_d$tf == tf]
    kc <- ev_c$key[ev_c$tf == tf]
    cnt$both[i]          <- length(intersect(kd, kc))
    cnt$default_only[i]  <- length(setdiff(kd, kc))
    cnt$coverage_only[i] <- length(setdiff(kc, kd))
  }
  cnt$total <- cnt$default_only + cnt$coverage_only + cnt$both
  cnt <- cnt[order(cnt$total, decreasing = TRUE), , drop = FALSE]
  write.csv(cnt, file.path(out_dir, "E1_per_tf_mode_counts.csv"),
            row.names = FALSE)

  # Bar plot: top-30 TFs by total event count, stacked by provenance
  top <- head(cnt, 30)
  long <- data.frame(
    tf     = rep(top$tf, 3),
    status = factor(rep(c("default_only","coverage_only","both"),
                         each = nrow(top)),
                    levels = c("default_only","coverage_only","both")),
    n      = c(top$default_only, top$coverage_only, top$both)
  )
  long$tf <- factor(long$tf, levels = rev(top$tf))

  p_bar <- ggplot(long, aes(x = n, y = tf, fill = status)) +
    geom_col(alpha = 0.9) +
    scale_fill_manual(values = c(default_only  = COLS$guide,
                                  coverage_only = COLS$high,
                                  both          = COLS$low),
                      name = NULL) +
    labs(x = "# events",
         y = NULL,
         title = "Per-TF event counts by mode",
         subtitle = paste0("Top-30 TFs | matching tolerance \u00b1",
                           bin_bp / 2, " bp")) +
    theme_caspex() +
    theme(legend.position = "bottom",
          panel.grid.major.y = element_blank())
  ggsave(file.path(out_dir, "E1_per_tf_mode_counts.pdf"),
         p_bar, width = 9, height = 8)

  # Jitter plot: position vs mode for the top-30 TFs
  jd <- rbind(
    if (nrow(ev_d) > 0)
      data.frame(tf = ev_d$tf, position = ev_d$position,
                 weight = ev_d$weight, mode = "default",
                 stringsAsFactors = FALSE)
    else data.frame(tf=character(0), position=numeric(0),
                    weight=numeric(0), mode=character(0)),
    if (nrow(ev_c) > 0)
      data.frame(tf = ev_c$tf, position = ev_c$position,
                 weight = ev_c$weight, mode = "coverage",
                 stringsAsFactors = FALSE)
    else data.frame(tf=character(0), position=numeric(0),
                    weight=numeric(0), mode=character(0))
  )
  jd <- jd[jd$tf %in% top$tf, ]
  jd$tf <- factor(jd$tf, levels = rev(top$tf))
  p_jit <- ggplot(jd, aes(x = position, y = tf, color = mode,
                           size = weight)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = COLS$tss, linewidth = 0.4) +
    geom_point(alpha = 0.75,
               position = position_nudge(
                 y = ifelse(jd$mode == "default", -0.18, 0.18))) +
    scale_color_manual(values = c(default  = COLS$guide,
                                   coverage = COLS$high),
                      name = NULL) +
    scale_size_continuous(range = c(1, 5), guide = "none") +
    labs(x = "position (bp, TSS-relative)", y = NULL,
         title = "Event positions by mode",
         subtitle = "Jittered vertically for readability; size \u221d \u03b2") +
    theme_caspex() +
    theme(legend.position = "bottom",
          panel.grid.major.y = element_blank())
  ggsave(file.path(out_dir, "E2_event_positions_by_mode.pdf"),
         p_jit, width = 10, height = 8)

  # Tables of mode-only calls (the distal rescues are in coverage_only)
  if (nrow(ev_d) > 0) {
    write.csv(ev_d[!(ev_d$key %in% ev_c$key), ,
                    drop = FALSE][, setdiff(names(ev_d), "key")],
              file.path(out_dir, "E3_default_only_events.csv"),
              row.names = FALSE)
  }
  if (nrow(ev_c) > 0) {
    write.csv(ev_c[!(ev_c$key %in% ev_d$key), ,
                    drop = FALSE][, setdiff(names(ev_c), "key")],
              file.path(out_dir, "E3_coverage_only_events.csv"),
              row.names = FALSE)
  }

  message("  wrote per-TF counts, bar + jitter plots, and mode-only CSVs.")
  invisible(list(counts = cnt, bar = p_bar, jitter = p_jit,
                 out_dir = out_dir_abs))
}

# =============================================================================
# F.  Convenience wrapper
# =============================================================================

#' Produce every extra plot as a separate PDF in `out_dir`.
#'
#' Auto-detects whether `result` came from a default-mode or coverage-aware
#' `run_caspex()` call (via `result$coverage_correct`) and routes steps
#' accordingly:
#'
#'   * Default mode (coverage_correct = FALSE):
#'       A1-A6, B1-B4, C1-C4      (B4 NNLS-residual and C4 motif-vs-\u03b2
#'                                  are NNLS-specific)
#'   * Coverage-aware (coverage_correct = TRUE):
#'       A1-A6, B1-B3, C1-C3, D1-D3
#'       (B4/C4 skipped: they inspect NNLS internals that don't exist in
#'        coverage-aware mode; D1-D3 are coverage-aware diagnostics)
#'
#' @param result        From run_caspex()
#' @param out_dir       Output directory (created if missing)
#' @param n_perm        Number of permutations for the null distribution
#' @param sigmas        Kernel widths for the sensitivity grid
#' @param cov_floors    cov_floor sweep used by D2 (coverage mode only)
#' @param one_pager_tfs TFs for which to emit a one-pager (default: top-10
#'                      composite TFs \u222a every motif-scanned TF)
#' @param skip          Character vector of step names to skip. Any of:
#'                      "one_pager","family","event_density","scatter",
#'                      "cooccurrence","ranked_events",
#'                      "permutation","sigma","jackknife","residual",
#'                      "volcano","correlation","pvalhist","motif_vs_beta",
#'                      "cov_rescue","cov_floor_sweep","cov_stack"
#' @return An invisible list of all generated objects
run_caspex_extras <- function(result,
                              out_dir        = "caspex_output/extras",
                              n_perm         = 500,
                              sigmas         = c(100, 200, 250, 300, 500),
                              cov_floors     = c(0.02, 0.05, 0.10, 0.20),
                              one_pager_tfs  = NULL,
                              skip           = character(0)) {
  if (missing(result) || is.null(result))
    stop("run_caspex_extras() needs the object returned by run_caspex(). ",
         "Call run_caspex_extras(result) — do not pass NULL.")
  for (f in c("long_data", "spatial_df", "pos_map"))
    if (is.null(result[[f]]))
      stop("result$", f, " is NULL. Did run_caspex() finish successfully?")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir_abs <- normalizePath(out_dir, mustWork = FALSE)
  is_cov <- isTRUE(result$coverage_correct)
  mode_lbl <- if (is_cov)
    paste0("coverage-aware (cov_floor=",
           result$cov_floor %||% 0.05, ")")
  else "default NNLS"
  message("\n=== CasPEX extras: writing to ", out_dir_abs, " ===")
  message("    binding mode: ", mode_lbl)

  # Auto-skip steps that don't apply to the detected mode. Users can
  # override by passing `skip` themselves — this only ADDS to the skip list.
  auto_skip <- if (is_cov) {
    # NNLS-internals diagnostics don't apply to coverage-aware runs
    c("residual", "motif_vs_beta")
  } else {
    # D.x coverage-aware diagnostics don't apply to default-mode runs
    c("cov_rescue", "cov_floor_sweep", "cov_stack")
  }
  skip <- unique(c(skip, auto_skip))

  out         <- list()
  step_status <- list()                   # name -> "ok" | "skipped" | "failed: ..."
  do_step     <- function(name) !(name %in% skip)

  # Wrap each step so one failure can't kill the rest. `expr` must write
  # its own file(s); on error we record the reason and continue.
  safely <- function(name, expr) {
    if (!do_step(name)) { step_status[[name]] <<- "skipped"; return(invisible(NULL)) }
    tryCatch({
      force(expr)
      step_status[[name]] <<- "ok"
    }, error = function(e) {
      step_status[[name]] <<- paste0("failed: ", conditionMessage(e))
      message("  !! step '", name, "' failed: ", conditionMessage(e))
      try(grDevices::dev.off(), silent = TRUE)  # close any half-open PDF
    })
  }

  # --- A.1 per-TF one-pagers (multi-page PDF)
  safely("one_pager", {
    if (is.null(one_pager_tfs)) {
      # Default: top-10 composite TFs ∪ every TF that was motif-scanned.
      # `motif_results` already contains the full union of the three
      # selection buckets (common + shared-focal + region-specific), so
      # this automatically produces a one-pager for every TF that shows
      # up on any deck in the main pipeline.
      top  <- head(as.character(result$spatial_df$protein), 10)
      scan <- if (!is.null(result$motif_results))
        names(result$motif_results) else character(0)
      one_pager_tfs <- unique(c(top, scan))
    }
    one_pager_tfs <- intersect(one_pager_tfs, result$long_data$protein)
    message("  [A1] one-pagers for ", length(one_pager_tfs), " TFs")
    pdf(file.path(out_dir, "A1_tf_one_pagers.pdf"), width = 11, height = 9)
    for (tf in one_pager_tfs) {
      p <- tryCatch(plot_tf_one_pager(result, tf),
                    error = function(e) {
                      message("    skipping ", tf, ": ", conditionMessage(e))
                      NULL
                    })
      if (!is.null(p)) print(p)
    }
    dev.off()
    out$one_pager_tfs <<- one_pager_tfs
  })

  # --- A.2 TF-family enrichment (may take time; queries JASPAR)
  safely("family", {
    message("  [A2] TF-family enrichment (JASPAR queries)")
    p <- plot_tf_family_enrichment(result)
    ggsave(file.path(out_dir, "A2_tf_family.pdf"), p, width = 10, height = 5.5)
    out$family <<- p
  })

  # --- A.3 event density
  safely("event_density", {
    message("  [A3] event density")
    p <- plot_event_density(result)
    ggsave(file.path(out_dir, "A3_event_density.pdf"), p, width = 10, height = 4)
    out$event_density <<- p
  })

  # --- A.4 composite vs specificity
  safely("scatter", {
    message("  [A4] composite vs specificity")
    p <- plot_composite_vs_specificity(result)
    ggsave(file.path(out_dir, "A4_composite_vs_specificity.pdf"),
           p, width = 8.5, height = 6.5)
    out$scatter <<- p
  })

  # --- A.5 TF-pair co-occurrence heatmap (mode-agnostic)
  safely("cooccurrence", {
    message("  [A5] TF-pair co-occurrence")
    p <- plot_tf_cooccurrence(result)
    ggsave(file.path(out_dir, "A5_tf_cooccurrence.pdf"),
           p, width = 9, height = 8)
    out$cooccurrence <<- p
  })

  # --- B.1 permutation null
  safely("permutation", {
    message("  [B1] permutation null (B=", n_perm, ")")
    perm <- run_permutation_null(result, n_perm = n_perm)
    p <- plot_permutation_null(perm)
    ggsave(file.path(out_dir, "B1_permutation_null.pdf"), p, width = 9, height = 7)
    write.csv(perm$summary,
              file.path(out_dir, "B1_permutation_null_summary.csv"),
              row.names = FALSE)
    out$permutation <<- perm
  })

  # --- B.2 sigma sensitivity
  safely("sigma", {
    message("  [B2] sigma sensitivity")
    sg <- run_sigma_sensitivity(result, sigmas = sigmas)
    if (!is.null(sg)) {
      p <- plot_sigma_sensitivity(sg)
      ggsave(file.path(out_dir, "B2_sigma_sensitivity.pdf"),
             p, width = 12, height = 8)
      out$sigma <<- sg
    }
  })

  # --- B.3 event jackknife
  safely("jackknife", {
    message("  [B3] event jackknife")
    jk <- run_event_jackknife(result)
    if (!is.null(jk)) {
      p <- plot_event_jackknife(jk)
      ggsave(file.path(out_dir, "B3_event_jackknife.pdf"),
             p, width = 8, height = 10)
      write.csv(jk, file.path(out_dir, "B3_event_jackknife.csv"),
                row.names = FALSE)
      out$jackknife <<- jk
    }
  })

  # --- A.6 ranked events — consumes B.3 if present. Placed AFTER jackknife
  # so the confidence score can include the survival fraction; if jackknife
  # failed or was skipped, the score falls back to per-TF-normalized beta.
  safely("ranked_events", {
    message("  [A6] ranked event table")
    rk <- rank_binding_events(result, jk_result = out$jackknife)
    if (!is.null(rk)) {
      ggsave(file.path(out_dir, "A6_ranked_events.pdf"),
             rk$plot, width = 9, height = 10)
      write.csv(rk$ranked, file.path(out_dir, "A6_ranked_events.csv"),
                row.names = FALSE)
      out$ranked_events <<- rk
    }
  })

  # --- B.4 NNLS residual — one page per TF with motif hits
  safely("residual", {
    tfs_r <- intersect(unique(result$binding_events$tf),
                        names(result$motif_results))
    if (length(tfs_r) > 0) {
      message("  [B4] NNLS residual for ", length(tfs_r), " TFs")
      pdf(file.path(out_dir, "B4_nnls_residual.pdf"), width = 9, height = 7)
      for (tf in tfs_r) {
        p <- tryCatch(plot_nnls_residual(result, tf),
                      error = function(e) NULL)
        if (!is.null(p)) print(p)
      }
      dev.off()
      out$residual_tfs <<- tfs_r
    } else {
      message("  [B4] skipped: no TFs with both events and motif hits")
    }
  })

  # --- C.1 volcano per region
  safely("volcano", {
    message("  [C1] volcano per region")
    p <- plot_volcano_per_region(result)
    ggsave(file.path(out_dir, "C1_volcano_per_region.pdf"),
           p, width = 12, height = 8)
    out$volcano <<- p
  })

  # --- C.2 region correlation
  safely("correlation", {
    message("  [C2] region correlation")
    p <- plot_region_correlation(result)
    ggsave(file.path(out_dir, "C2_region_correlation.pdf"),
           p, width = 7, height = 6)
    out$correlation <<- p
  })

  # --- C.3 p-value histograms
  safely("pvalhist", {
    message("  [C3] p-value histograms")
    p <- plot_pval_histograms(result)
    ggsave(file.path(out_dir, "C3_pval_histograms.pdf"),
           p, width = 11, height = 7)
    out$pvalhist <<- p
  })

  # --- C.4 motif score vs NNLS beta
  safely("motif_vs_beta", {
    message("  [C4] motif score vs NNLS beta")
    p <- plot_motif_vs_nnls(result)
    ggsave(file.path(out_dir, "C4_motif_vs_nnls_beta.pdf"),
           p, width = 12, height = 9)
    out$motif_vs_beta <<- p
  })

  # --- D.1 coverage-rescue audit scatter (auto-skipped for default mode)
  safely("cov_rescue", {
    message("  [D1] coverage-rescue audit scatter")
    p <- plot_coverage_rescue_scatter(result)
    ggsave(file.path(out_dir, "D1_coverage_rescue.pdf"),
           p, width = 9, height = 7)
    out$cov_rescue <<- p
  })

  # --- D.2 cov_floor sensitivity sweep (auto-skipped for default mode)
  safely("cov_floor_sweep", {
    message("  [D2] cov_floor sensitivity (floors: ",
            paste(cov_floors, collapse = ", "), ")")
    cf <- run_covfloor_sensitivity(result, floors = cov_floors)
    if (!is.null(cf)) {
      p <- plot_covfloor_sensitivity(cf)
      ggsave(file.path(out_dir, "D2_covfloor_sensitivity.pdf"),
             p, width = 12, height = 8)
      write.csv(cf, file.path(out_dir, "D2_covfloor_sensitivity.csv"),
                row.names = FALSE)
      out$cov_floor_sweep <<- cf
    }
  })

  # --- D.3 s(x) / C(x) / beta(x) per-TF stack (auto-skipped for default mode)
  safely("cov_stack", {
    tfs_s <- intersect(unique(result$binding_events$tf),
                        names(result$motif_results))
    if (length(tfs_s) > 0) {
      message("  [D3] s/C/beta stack for ", length(tfs_s), " TFs")
      pdf(file.path(out_dir, "D3_coverage_stack.pdf"), width = 9, height = 9)
      for (tf in tfs_s) {
        p <- tryCatch(plot_coverage_stack(result, tf),
                      error = function(e) {
                        message("    skipping ", tf, ": ",
                                conditionMessage(e))
                        NULL
                      })
        if (!is.null(p)) print(p)
      }
      dev.off()
      out$cov_stack_tfs <<- tfs_s
    } else {
      message("  [D3] skipped: no TFs with both events and motif hits")
    }
  })

  # Final tally — what worked, what didn't, what was skipped
  message("\n--- Extras summary -----------------------------------------------")
  for (nm in names(step_status))
    message(sprintf("  %-14s : %s", nm, step_status[[nm]]))
  files_written <- list.files(out_dir, full.names = FALSE)
  message("  files written  : ", length(files_written))
  if (length(files_written) > 0)
    message("    ", paste(files_written, collapse = ", "))
  else
    message("    (none — check the errors above)")
  message("------------------------------------------------------------------\n")

  out$step_status    <- step_status
  out$files_written  <- files_written
  out$out_dir        <- out_dir_abs
  invisible(out)
}
