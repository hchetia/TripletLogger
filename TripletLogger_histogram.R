#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# TripletLogger_histogram.R
#
# Summary-histogram companion for TripletLogger.
#
# Two ways to use it:
#   1. Standalone, on existing TripletLogger output:
#        Rscript TripletLogger_histogram.R \
#          --freq    sample_CAGTRIPLETLogger.v1.0.long_NumReadsPerRepeat.csv \
#          --metrics sample_CAGTRIPLETLogger.v1.0.long_RepeatMetrics.csv
#
#   2. Sourced by TripletLogger, which calls plot_repeat_histogram() directly
#      when run with histo = "ON".
#
# The CLI block at the bottom only runs when this file is executed directly
# (sys.nframe() == 0); when source()'d it merely defines plot_repeat_histogram().
# ══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(ggplot2))

# ── Plotting function ─────────────────────────────────────────────────────────
#' Plot a repeat-length read-count histogram with optional allele-peak markers.
#'
#' @param freq_df      data.frame with integer columns `Repeat_Length` and
#'                     `NumReads` (the *_NumReadsPerRepeat.csv contents).
#' @param allele1_peak,allele2_peak  repeat lengths of the called allele peaks;
#'                     pass NA to omit a marker.
#' @param freq_min,freq_max  x-axis (repeat-length) limits. NA = auto from the
#'                     populated range, padded, and widened to include peaks.
#' @param reads_min,reads_max  y-axis (read-count) limits. Default 0-100. Bars
#'                     above reads_max are clipped, not dropped; peak markers are
#'                     drawn at their true positions regardless.
#' @param format      "png" or "svg".
#' @param out_file    output path (extension should match `format`).
#' @param sample_id   title text.
#' @param subtitle    optional subtitle (e.g. "CAG | long").
#' @return out_file (invisibly).
plot_repeat_histogram <- function(freq_df,
                                  allele1_peak = NA_real_,
                                  allele2_peak = NA_real_,
                                  freq_min     = 0,
                                  freq_max     = NA_real_,
                                  reads_min    = 0,
                                  reads_max    = 1000,
                                  format       = c("png", "svg"),
                                  out_file,
                                  sample_id    = "",
                                  subtitle     = NULL,
                                  width        = 8,
                                  height       = 5,
                                  dpi          = 150) {
 
 format <- match.arg(format)
 
 if (!all(c("Repeat_Length", "NumReads") %in% names(freq_df))) {
  stop("freq_df must contain columns 'Repeat_Length' and 'NumReads'.")
 }
 freq_df <- data.frame(
  Repeat_Length = as.numeric(freq_df$Repeat_Length),
  NumReads      = as.numeric(freq_df$NumReads)
 )
 freq_df <- freq_df[is.finite(freq_df$Repeat_Length) &
                     is.finite(freq_df$NumReads), , drop = FALSE]
 if (nrow(freq_df) == 0L) stop("No usable rows in freq_df.")
 
 peaks <- c(allele1_peak, allele2_peak)
 peaks <- suppressWarnings(as.numeric(peaks))
 
 # ── x-limits ────────────────────────────────────────────────────────────────
 populated <- freq_df$Repeat_Length[freq_df$NumReads > 0]
 if (length(populated) == 0L) populated <- freq_df$Repeat_Length
 if (is.na(freq_min)) {
  freq_min <- min(c(populated, peaks[is.finite(peaks)]), na.rm = TRUE) - 2
 }
 if (is.na(freq_max)) {
  freq_max <- max(c(populated, peaks[is.finite(peaks)]), na.rm = TRUE) + 2
 }
 freq_min <- max(0, floor(freq_min))
 freq_max <- ceiling(freq_max)
 
 # x-axis breaks every 5 units across the visible range
 x_breaks <- seq(floor(freq_min / 5) * 5, ceiling(freq_max / 5) * 5, by = 5)
 
 # ── base bar layer ────────────────────────────────────────────────────────────
 p <- ggplot(freq_df, aes(x = .data$Repeat_Length, y = .data$NumReads)) +
  geom_col(width = 0.9, fill = "#3B6FB6", colour = NA)
 
 # ── allele-peak markers ───────────────────────────────────────────────────────
 peak_cols <- c("#D1495B", "#2E7D32")   # allele 1, allele 2
 peak_lab  <- c("Est.Allele1", "Est.Allele2")
 label_y   <- reads_min + 0.95 * (reads_max - reads_min)
 for (k in seq_along(peaks)) {
  pk <- peaks[k]
  if (!is.finite(pk)) next
  p <- p +
   geom_vline(xintercept = pk, colour = peak_cols[k],
              linetype = "dashed", linewidth = 0.7) +
   annotate("text", x = pk, y = label_y,
            label = sprintf("%s: %g", peak_lab[k], pk),
            colour = peak_cols[k], hjust = -0.05, vjust = 1,
            size = 2.6, fontface = "bold")
 }
 
 # ── scales / theme ────────────────────────────────────────────────────────────
 p <- p +
  scale_x_continuous(breaks = x_breaks) +
  coord_cartesian(xlim = c(freq_min, freq_max),
                  ylim = c(reads_min, reads_max), expand = FALSE) +
  labs(
   title    = if (nzchar(sample_id)) sample_id else NULL,
   subtitle = subtitle,
   x        = "Repeat length (triplets)",
   y        = "Number of reads"
  ) +
  theme_bw(base_size = 12) +
  theme(
   panel.grid.minor = element_blank(),
   plot.title       = element_text(face = "bold"),
   plot.subtitle    = element_text(colour = "grey30")
  )
 
 # ── write via base graphics device (no svglite/ragg dependency) ──────────────
 dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
 if (format == "png") {
  grDevices::png(out_file, width = width * dpi, height = height * dpi,
                 res = dpi)
 } else {
  grDevices::svg(out_file, width = width, height = height)
 }
 on.exit(grDevices::dev.off(), add = TRUE)
 print(p)
 
 invisible(out_file)
}

# ══════════════════════════════════════════════════════════════════════════════
# CLI entry point — runs ONLY when executed directly (not when source()'d).
# ══════════════════════════════════════════════════════════════════════════════
if (sys.nframe() == 0L) {
 
 suppressPackageStartupMessages(library(optparse))
 
 opt_list <- list(
  make_option("--freq", type = "character", default = NULL,
              help = "Path to *_NumReadsPerRepeat.csv (required)"),
  make_option("--metrics", type = "character", default = NULL,
              help = "Path to *_RepeatMetrics.csv (optional; draws allele peaks)"),
  make_option("--out", type = "character", default = NULL,
              help = "Output image path [default: derived from --freq]"),
  make_option("--freqMin", type = "double", default = 0,
              help = "Min repeat length on x-axis [default: 0]"),
  make_option("--freqMax", type = "double", default = NA,
              help = "Max repeat length on x-axis [default: auto]"),
  make_option("--readsMin", type = "double", default = 0,
              help = "Min read count on y-axis [default: 0]"),
  make_option("--readsMax", type = "double", default = 1000,
              help = "Max read count on y-axis [default: 1000]"),
  make_option("--format", type = "character", default = "png",
              help = "Output format: png or svg [default: png]"),
  make_option("--sampleId", type = "character", default = NULL,
              help = "Title for the plot [default: derived from --freq]")
 )
 
 opt <- parse_args(OptionParser(option_list = opt_list))
 
 if (is.null(opt$freq)) stop("--freq is required.", call. = FALSE)
 if (!file.exists(opt$freq)) stop(sprintf("Not found: %s", opt$freq), call. = FALSE)
 fmt <- tolower(opt$format)
 if (!fmt %in% c("png", "svg")) stop("--format must be png or svg.", call. = FALSE)
 
 freq_df <- utils::read.csv(opt$freq, stringsAsFactors = FALSE)
 
 a1 <- NA_real_; a2 <- NA_real_
 if (!is.null(opt$metrics) && file.exists(opt$metrics)) {
  m <- utils::read.csv(opt$metrics, stringsAsFactors = FALSE)
  if ("Estimated_Allele_1_Peak" %in% names(m)) a1 <- suppressWarnings(as.numeric(m$Estimated_Allele_1_Peak[1]))
  if ("Estimated_Allele_2_Peak" %in% names(m)) a2 <- suppressWarnings(as.numeric(m$Estimated_Allele_2_Peak[1]))
 }
 
 out_file <- opt$out
 if (is.null(out_file)) {
  out_file <- sub("_NumReadsPerRepeat\\.csv$", paste0("_Histogram.", fmt), opt$freq)
  if (identical(out_file, opt$freq)) out_file <- paste0(opt$freq, "_Histogram.", fmt)
 }
 
 sample_id <- opt$sampleId
 if (is.null(sample_id)) sample_id <- sub("_NumReadsPerRepeat\\.csv$", "", basename(opt$freq))
 
 plot_repeat_histogram(
  freq_df      = freq_df,
  allele1_peak = a1,
  allele2_peak = a2,
  freq_min     = opt$freqMin,
  freq_max     = opt$freqMax,
  reads_min    = opt$readsMin,
  reads_max    = opt$readsMax,
  format       = fmt,
  out_file     = out_file,
  sample_id    = sample_id
 )
 
 message(sprintf("Histogram written to: %s", out_file))
}