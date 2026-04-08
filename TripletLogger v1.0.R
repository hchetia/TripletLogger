#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# TRIPLET LOGGER v1.0 — R implementation
# Alignment-free estimation of CAG or CTG repeats from targeted amplicon short-read  (MiSeq) or long-read (Nanopore) sequencing data.
#
# Miseq-appropriate defaults:
#   - Raw frequency range (default 1–120 repeats) 
#   - Threshold metrics: config-specific 'base threshold' (CAG: >=35; CTG: >=50),
#   - plus >110 (all relative to the 'base threshold')
# 
# ONT-appropriate defaults:
#   - Reduced chunk size (20,000) for memory safety with long reads
#   - Extended raw frequency range (default 1–1000) for large expansions
#   - Threshold metrics: config-specific base (CAG: >=35; CTG: >=50),
#     plus >110, >150 (all relative to the base threshold)
#
# Counting method:
#   Only pure target triplets (CAG or CTG) within each accepted regex match are counted.  Non-target triplets (sequencing errors, flank elements) are tolerated for match continuity but are NOT included in the reported repeat length.
#   Type A (total non-target) and Type B (consecutive non-target) error
#   thresholds gate whether a match is accepted.
#
# Triplet-specific logic:
#   CAG mode — regex: CAG(CAG|CA[ATCGN]|C[ATCGN]G|[ATCGN]AG)+
#              counting: pure CAG triplets only (errors excluded)
#              Type B: reject entire match if consecutive errors > threshold
#              flank tract: CCG-like  (trim absorbed CCG repeats)
#              flank exclusion: ...CAA-CAG at end of tract
#              default max_type_b: 1
#
#   CTG mode — regex: CTG(CTG|CT[ATCGN]|C[ATCGN]G|[ATCGN]TG)+
#              counting: total triplet count (including substitutions)
#              Type B: truncate tract at the max_type_b-th consecutive-
#                      error pair (matching original CTG Logger behaviour)
#              flank tract: CGG-like  (trim absorbed CGG repeats)
#              default max_type_b: 2
#
# Length-adaptive Type A threshold (Fix G):
#   effective_max_type_a = max(max_type_a, floor(n_triplets * type_a_rate))
#
# Dependencies: ShortRead, Biostrings, dplyr, optparse.
# ═════════════════════════════════════════════════════════════════════

VERSION <- "1.0"

# ── Load required packages ──
suppressPackageStartupMessages({
 library(ShortRead)
 library(Biostrings)
 library(dplyr)
 library(optparse)
})

# ══════════════════════════════════════════════════════════════════════════════
# TRIPLET-TYPE CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

#' Return a configuration list for the requested triplet type.
#'
#' Each config contains:
#'   target      – the pure triplet to count (e.g. "CAG")
#'   pattern     – N-aware regex for the repeat tract
#'   flank_ref   – reference triplet for the downstream flank tract
#'                 (CCG for CAG, CGG for CTG)
#'   interrupt   – interruption triplets that signal the end of the pure
#'                 repeat tract (CAA/CAN/NAA for CAG; TTG/NTG/TNG for CTG)
#'   label       – human-readable label used in messages and filenames

get_triplet_config <- function(triplet_type = c("CAG", "CTG")) {
 triplet_type <- match.arg(toupper(triplet_type), c("CAG", "CTG"))
 
 if (triplet_type == "CAG") {
  list(
   target           = "CAG",
   pattern          = "CAG(CAG|CA[ATCGN]|C[ATCGN]G|[ATCGN]AG)+",
   flank_ref        = c("C", "C", "G"),
   interrupt        = c("CAA", "CAN", "NAA"),
   label            = "CAG",
   count_all        = FALSE,   # report pure-target triplets only
   truncate_type_b  = FALSE,   # reject entire match on Type B excess
   default_max_type_b = 1L,
   subtract_n       = 0L       # no library-prep adjustment for CAG
  )
 } else {
  list(
   target           = "CTG",
   pattern          = "CTG(CTG|CT[ATCGN]|C[ATCGN]G|[ATCGN]TG)+",
   flank_ref        = c("C", "G", "G"),
   interrupt        = c("TTG", "NTG", "TNG"),
   label            = "CTG",
   count_all        = TRUE,    # report total triplet count (incl. errors)
   truncate_type_b  = TRUE,    # truncate tract at Nth consecutive-error pair
   default_max_type_b = 2L,
   subtract_n       = 6L       # library preparation strategy-dependent adjustment for CTG
  )
 }
}

# ══════════════════════════════════════════════════════════════════════════════
# CORE FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

# ── Reverse complement (both strands searched) ───────────────────────────────

rc <- function(seq) {
 chartr("ACGTacgtN", "TGCAtgcaN",
        paste(rev(strsplit(seq, "")[[1]]), collapse = ""))
}

# ── Flank-tract detector (generalised CCG-like / CGG-like) ───────────────────
#
# Returns TRUE if `triplet` is a 0- or 1-substitution/N variant of the
# flank reference BUT is NOT the target triplet itself.

is_flank_like <- function(triplet, config) {
 if (triplet == config$target) return(FALSE)
 chars      <- strsplit(triplet, "")[[1]]
 mismatches <- sum(chars != config$flank_ref & chars != "N")
 mismatches <= 1L
}

# ── Triplet extraction (N-aware, adaptive correction) ────────────────────────

#' Extract the longest valid repeat match from a single sequence.
#' For CAG: reported length = count of pure target triplets only.
#' For CTG: reported length = total triplet count (incl. errors) after
#'          truncation at the max_type_b-th consecutive-error pair.
#' Returns a named list: length (integer or NA), hits_end, rejected_triplets.

extract_triplet_from_read <- function(seq,
                                      config,
                                      max_type_a  = 5L,
                                      max_type_b  = NULL,
                                      type_a_rate = 0.05) {
 
 if (is.null(max_type_b)) max_type_b <- config$default_max_type_b
 
 matches  <- gregexpr(config$pattern, seq, perl = TRUE)
 raw_hits <- regmatches(seq, matches)[[1]]
 
 if (length(raw_hits) == 0L) {
  return(list(length = NA_integer_, hits_end = FALSE,
              rejected_triplets = NA_integer_))
 }
 
 starts   <- as.integer(matches[[1]])
 read_len <- nchar(seq)
 
 best_length            <- NA_integer_
 best_hits_end          <- FALSE
 longest_rejected_count <- NA_integer_
 
 for (idx in seq_along(raw_hits)) {
  match_str  <- raw_hits[idx]
  usable_len <- (nchar(match_str) %/% 3L) * 3L
  if (usable_len < 3L) next
  match_str <- substr(match_str, 1L, usable_len)
  
  triplets <- vapply(
   seq(1L, usable_len, by = 3L),
   function(i) substr(match_str, i, i + 2L),
   character(1)
  )
  
  # ── Trim absorbed flank tract (Fix F generalised) ─────────────────
  consec_flank <- 0L
  trim_at      <- length(triplets) + 1L
  for (ti in seq_along(triplets)) {
   if (is_flank_like(triplets[ti], config)) {
    consec_flank <- consec_flank + 1L
    if (consec_flank >= 2L) {
     trim_at <- ti - consec_flank + 1L
     break
    }
   } else {
    consec_flank <- 0L
   }
  }
  if (trim_at <= length(triplets)) {
   triplets   <- triplets[seq_len(trim_at - 1L)]
   usable_len <- length(triplets) * 3L
  }
  if (length(triplets) < 1L) next
  
  # ── Type B handling: truncate (CTG) vs reject (CAG) ─────────────
  if (isTRUE(config$truncate_type_b)) {
   # CTG mode: scan for consecutive-error pairs and truncate at
   # the max_type_b-th occurrence (matching old CTG Logger behaviour).
   mut      <- vapply(triplets, function(t) t != config$target, logical(1))
   bcount   <- 0L
   keep_len <- length(triplets)
   
   if (length(triplets) >= 2L) {
    for (ti in 2:length(triplets)) {
     if (mut[ti] && mut[ti - 1L]) {
      bcount <- bcount + 1L
      if (bcount >= max_type_b) {
       keep_len <- ti - 1L
       break
      }
     }
    }
   }
   triplets   <- triplets[seq_len(keep_len)]
   usable_len <- length(triplets) * 3L
   if (length(triplets) < 1L) next
   
   # Type A check on the (possibly truncated) tract
   type_a               <- sum(triplets != config$target)
   n_triplets           <- length(triplets)
   effective_max_type_a <- max(max_type_a,
                               as.integer(floor(n_triplets * type_a_rate)))
   
   if (type_a > effective_max_type_a) {
    if (is.na(longest_rejected_count) || n_triplets > longest_rejected_count) {
     longest_rejected_count <- n_triplets
    }
    next
   }
   
   # CTG: report total triplet count (including errors)
   report_count <- n_triplets
   
  } else {
   # CAG mode (original behaviour): reject entire match on excess errors
   
   # ── Error counting ──────────────────────────────────────────────
   type_a   <- 0L
   type_b   <- 0L
   prev_err <- FALSE
   
   for (t in triplets) {
    is_err <- (t != config$target)
    if (is_err) {
     type_a <- type_a + 1L
     if (prev_err) type_b <- type_b + 1L
    }
    prev_err <- is_err
   }
   
   # ── Fix G: Length-adaptive Type A threshold ─────────────────────
   n_triplets           <- length(triplets)
   effective_max_type_a <- max(max_type_a,
                               as.integer(floor(n_triplets * type_a_rate)))
   
   if (type_a > effective_max_type_a || type_b > max_type_b) {
    if (is.na(longest_rejected_count) || n_triplets > longest_rejected_count) {
     longest_rejected_count <- n_triplets
    }
    next
   }
   
   # ── Pure target count ───────────────────────────────────────────
   pure_count <- sum(triplets == config$target)
   
   # ── Flank target exclusion ──────────────────────────────────────
   #   CAG: penultimate in {CAA,CAN,NAA} + last==CAG  →  subtract 1
   #   Also checks for flank-like last + interrupt ante-penultimate.
   n_tr <- length(triplets)
   if (n_tr >= 2L) {
    last   <- triplets[n_tr]
    penult <- triplets[n_tr - 1L]
    is_interrupt <- penult %in% config$interrupt
    
    if (is_interrupt && last == config$target) {
     pure_count <- pure_count - 1L
    } else if (is_flank_like(last, config) && n_tr >= 3L) {
     ante   <- triplets[n_tr - 2L]
     middle <- triplets[n_tr - 1L]
     if (ante %in% config$interrupt && middle == config$target) {
      pure_count <- pure_count - 1L
     }
    }
   }
   
   report_count <- pure_count
  }
  
  # ── Does the match reach the end of the read? ─────────────────────
  match_end_pos <- starts[idx] + usable_len - 1L
  hits_read_end <- (read_len - match_end_pos) <= 3L
  
  if (is.na(best_length) || report_count > best_length) {
   best_length   <- report_count
   best_hits_end <- hits_read_end
  }
 }
 
 rejected_out <- if (!is.na(best_length)) NA_integer_ else longest_rejected_count
 
 list(length = best_length, hits_end = best_hits_end,
      rejected_triplets = rejected_out)
}

#' Search both orientations, return best result.

extract_triplet_both_strands <- function(seq,
                                         config,
                                         max_type_a  = 5L,
                                         max_type_b  = NULL,
                                         type_a_rate = 0.05) {
 
 if (is.null(max_type_b)) max_type_b <- config$default_max_type_b
 
 fwd <- extract_triplet_from_read(seq,     config, max_type_a, max_type_b,
                                  type_a_rate)
 rev <- extract_triplet_from_read(rc(seq), config, max_type_a, max_type_b,
                                  type_a_rate)
 
 fwd_len <- ifelse(is.na(fwd$length), 0L, fwd$length)
 rev_len <- ifelse(is.na(rev$length), 0L, rev$length)
 
 if (fwd_len >= rev_len) {
  return(list(
   length            = if (fwd_len > 0L) fwd_len else NA_integer_,
   hits_end          = fwd$hits_end,
   rejected_triplets = if (fwd_len > 0L) NA_integer_ else fwd$rejected_triplets
  ))
 } else {
  return(list(
   length            = if (rev_len > 0L) rev_len else NA_integer_,
   hits_end          = rev$hits_end,
   rejected_triplets = if (rev_len > 0L) NA_integer_ else rev$rejected_triplets
  ))
 }
}

# ── Allele calling (kernel density peak detection) ────────────────────────────

call_alleles <- function(repeat_lengths, frequencies, bw = 1.5) {
 
 total <- sum(frequencies)
 if (length(repeat_lengths) < 3 || total < 20) {
  return(data.frame(Allele = NA, Peak_Repeat = NA, Proportion = NA,
                    Call = "insufficient_data", stringsAsFactors = FALSE))
 }
 
 obs <- rep(repeat_lengths, times = frequencies)
 d   <- density(obs, bw = bw, n = 1024)
 
 dy    <- diff(d$y)
 peaks <- which(dy[-length(dy)] > 0 & dy[-1] < 0) + 1L
 peak_x <- d$x[peaks]
 peak_y <- d$y[peaks]
 
 keep   <- peak_y >= 0.01 * max(peak_y)
 peak_x <- peak_x[keep]
 peak_y <- peak_y[keep]
 
 if (length(peak_x) == 0L) {
  return(data.frame(Allele = NA, Peak_Repeat = NA, Proportion = NA,
                    Call = "no_peaks", stringsAsFactors = FALSE))
 }
 
 ord    <- order(peak_y, decreasing = TRUE)
 peak_x <- peak_x[ord][seq_len(min(2, length(ord)))]
 # Snap each KDE peak to the nearest observed modal repeat length
 peak_x <- vapply(peak_x, function(px) {
  in_window <- abs(repeat_lengths - px) <= bw
  if (any(in_window)) {
   repeat_lengths[in_window][which.max(frequencies[in_window])]
  } else {
   round(px)
  }
 }, numeric(1))
 peak_x <- sort(unique(as.integer(peak_x)))
 
 assignments <- vapply(obs, function(o) peak_x[which.min(abs(o - peak_x))],
                       numeric(1))
 props <- table(assignments) / length(obs)
 
 call <- if (length(peak_x) == 1L) "unimodal" else "multimodal"
 
 data.frame(
  Allele      = seq_along(peak_x),
  Peak_Repeat = peak_x,
  Proportion  = as.numeric(props),
  Call        = call,
  stringsAsFactors = FALSE
 )
}

# ── Streaming FASTQ reader (ShortRead backend) ───────────────────────────────

fastq_streamer <- function(filepath, chunk_size = 20000L, q_threshold = 0) {
 
 sr_stream <- FastqStreamer(filepath, n = chunk_size)
 
 list(
  yield = function() {
   chunk <- yield(sr_stream)
   if (length(chunk) == 0L) return(NULL)
   
   n_total <- length(chunk)
   
   qual_matrix <- as(quality(chunk), "matrix")
   mean_quals  <- rowMeans(qual_matrix, na.rm = TRUE)
   keep        <- mean_quals >= q_threshold
   n_passed    <- sum(keep)
   
   if (n_passed == 0L) {
    return(list(
     data     = data.frame(id = character(0), seq = character(0),
                           qual = character(0), stringsAsFactors = FALSE),
     n_total  = n_total,
     n_passed = 0L
    ))
   }
   
   chunk <- chunk[keep]
   seqs  <- as.character(sread(chunk))
   ids   <- as.character(id(chunk))
   quals <- as.character(chunk@quality@quality)
   
   list(
    data     = data.frame(id = ids, seq = seqs, qual = quals,
                          stringsAsFactors = FALSE),
    n_total  = n_total,
    n_passed = n_passed
   )
  },
  close = function() close(sr_stream)
 )
}

# ── Main processor ────────────────────────────────────────────────────────────

process_fastq_file <- function(fastq_path,
                               config,
                               sample_id         = NULL,
                               q_threshold       = 20,
                               min_repeat_length = 1L,
                               max_type_a        = 5L,
                               max_type_b        = NULL,
                               type_a_rate       = 0.05,
                               chunk_size        = 20000L,
                               allele_bw         = 1.5,
                               freq_range_min    = 1L,
                               freq_range_max    = 1000L,
                               read_type         = "long",
                               output_dir        = ".") {
 
 if (is.null(max_type_b)) max_type_b <- config$default_max_type_b
 
 triplet_label <- config$label
 read_type     <- tolower(read_type)
 
 # ── Short-read override: cap freq_range_max at 120 ──────────────────────
 if (read_type == "short") {
  freq_range_max <- min(freq_range_max, 120L)
 }
 
 # ── Sample ID ────────────────────────────────────────────────────────────
 if (is.null(sample_id)) {
  sample_id <- sub("\\.(fastq|fq)(\\.gz)?$", "", basename(fastq_path))
  sample_id <- sub("\\.gz$", "", sample_id)
 }
 
 # ── Streaming reader ─────────────────────────────────────────────────────
 stream <- fastq_streamer(fastq_path, chunk_size, q_threshold)
 on.exit(stream$close())
 
 total_reads     <- 0L
 reads_passing   <- 0L
 all_lengths     <- integer(0)
 
 n_bridged_reads        <- 0L
 censored_reads         <- 0L
 threshold_rejected     <- 0L
 rejected_triplet_sizes <- integer(0)
 
 chunk_num <- 0L
 repeat {
  result_chunk <- stream$yield()
  if (is.null(result_chunk)) break
  
  chunk_num     <- chunk_num + 1L
  n             <- result_chunk$n_total
  n_passed      <- result_chunk$n_passed
  total_reads   <- total_reads + n
  reads_passing <- reads_passing + n_passed
  chunk         <- result_chunk$data
  
  if (nrow(chunk) == 0L) {
   message(sprintf("  chunk %d: %d reads, 0 passed Q filter", chunk_num, n))
   next
  }
  
  for (i in seq_len(nrow(chunk))) {
   seq_i <- chunk$seq[i]
   
   result <- extract_triplet_both_strands(seq_i, config,
                                          max_type_a, max_type_b,
                                          type_a_rate)
   
   if (!is.na(result$length)) {
    all_lengths <- c(all_lengths, result$length)
    
    if (grepl("N", seq_i, fixed = TRUE)) {
     n_bridged_reads <- n_bridged_reads + 1L
    }
    if (isTRUE(result$hits_end)) {
     censored_reads <- censored_reads + 1L
    }
   } else if (!is.na(result$rejected_triplets)) {
    threshold_rejected     <- threshold_rejected + 1L
    rejected_triplet_sizes <- c(rejected_triplet_sizes,
                                result$rejected_triplets)
   }
  }
  
  message(sprintf("  chunk %d: %d reads processed (%d passed Q filter)",
                  chunk_num, n, n_passed))
 }
 
 reads_failed <- total_reads - reads_passing
 
 message(sprintf(
  "\n[%s] Total reads: %d | Passed Q>=%d: %d | Removed: %d",
  sample_id, total_reads, q_threshold, reads_passing, reads_failed))
 
 if (reads_passing == 0L) {
  warning(sprintf("No reads passed quality filter for %s.", sample_id))
  return(invisible(NULL))
 }
 
 reads_with_match <- length(all_lengths)
 message(sprintf(
  "[%s] Reads with %s match: %d / %d (%.2f%%)",
  sample_id, triplet_label, reads_with_match, reads_passing,
  100 * reads_with_match / reads_passing))
 
 message(sprintf(
  "[%s] N-bridged: %d | End-of-read censored: %d",
  sample_id, n_bridged_reads, censored_reads))
 
 # ── Threshold rejection diagnostics ──────────────────────────────────────
 message(sprintf(
  "[%s] Threshold-rejected reads: %d (type_a_rate=%.2f, max_type_a=%d, max_type_b=%d)",
  sample_id, threshold_rejected, type_a_rate, max_type_a, max_type_b))
 
 if (threshold_rejected > 0L) {
  message(sprintf(
   "[%s]   Rejected match sizes (triplets): min=%d, median=%d, max=%d",
   sample_id,
   min(rejected_triplet_sizes),
   as.integer(median(rejected_triplet_sizes)),
   max(rejected_triplet_sizes)))
 }
 
 if (reads_with_match == 0L) {
  warning(sprintf("No %s structures detected in %s.", triplet_label, sample_id))
  return(invisible(NULL))
 }
 
 # ── Build frequency table ────────────────────────────────────────────────
 all_lengths <- all_lengths[all_lengths > 0L]
 
 # ── Library preparation strategy-dependent adjustment (CTG: subtract_n = 6) ──
 subtract_n         <- config$subtract_n
 reads_dropped_adj  <- 0L
 if (subtract_n > 0L) {
  all_lengths       <- all_lengths - subtract_n
  n_before          <- length(all_lengths)
  all_lengths       <- all_lengths[all_lengths >= 1L]
  reads_dropped_adj <- n_before - length(all_lengths)
  message(sprintf("[%s] Applied library-prep strategy-dependent adjustment: subtract_n = %d",
                  sample_id, subtract_n))
 }
 
 freq_tbl    <- as.data.frame(table(all_lengths), stringsAsFactors = FALSE)
 names(freq_tbl) <- c("Repeat_Length", "NumReads")
 freq_tbl$Repeat_Length <- as.integer(freq_tbl$Repeat_Length)
 full_freq <- freq_tbl
 
 # ── Raw frequency table (continuous range) ───────────────────────
 # Spans freq_range_min:freq_range_max, extended if observed data falls
 # outside the default range.  Every integer size gets a row.
 # Floor is clamped to 1 so that post-subtraction negatives/zero are never reported.
 effective_min <- max(1L, min(freq_range_min,
                              if (nrow(freq_tbl) > 0L) min(freq_tbl$Repeat_Length) else freq_range_min))
 effective_max <- max(freq_range_max,
                      if (nrow(freq_tbl) > 0L) max(freq_tbl$Repeat_Length) else freq_range_max)
 
 raw_freq <- data.frame(Repeat_Length = seq(effective_min, effective_max),
                        stringsAsFactors = FALSE)
 raw_freq <- merge(raw_freq, full_freq, by = "Repeat_Length", all.x = TRUE)
 raw_freq$NumReads[is.na(raw_freq$NumReads)] <- 0L
 raw_freq <- raw_freq[order(raw_freq$Repeat_Length), ]
 
 # ── Apply user-configurable floor ────────────────────────────────────────
 freq_tbl <- freq_tbl[freq_tbl$Repeat_Length >= min_repeat_length, ]
 
 if (nrow(freq_tbl) == 0L) {
  warning(sprintf("No repeats >= %d after correction in %s.",
                  min_repeat_length, sample_id))
  return(invisible(NULL))
 }
 
 # ── Mode ─────────────────────────────────────────────────────────────────
 mode_repeat    <- freq_tbl$Repeat_Length[which.max(freq_tbl$NumReads)]
 longest_repeat <- max(freq_tbl$Repeat_Length)
 
 # ── Threshold metrics (config-specific base threshold) ───────────────────
 # CAG: base threshold >=35  (>34), labelled "atleast.35"
 # CTG: base threshold >=50  (>49), labelled "atleast.50"
 if (triplet_label == "CAG") {
  base_threshold     <- 35L
  base_threshold_lab <- "atleast.35"
 } else {
  base_threshold     <- 50L
  base_threshold_lab <- "atleast.50"
 }
 
 reads_over_base <- sum(freq_tbl$NumReads[freq_tbl$Repeat_Length >= base_threshold])
 reads_over_110  <- sum(freq_tbl$NumReads[freq_tbl$Repeat_Length >  110])
 reads_over_150  <- sum(freq_tbl$NumReads[freq_tbl$Repeat_Length >  150])
 
 # Percentages relative to reads >= base_threshold
 pct_over_110 <- if (reads_over_base > 0) {
  (reads_over_110 / reads_over_base) * 100
 } else {
  NA_real_
 }
 
 pct_over_150 <- if (reads_over_base > 0) {
  (reads_over_150 / reads_over_base) * 100
 } else {
  NA_real_
 }
 
 # ── Allele calling ──────────────────────────────────────────────────────
 allele_df <- call_alleles(freq_tbl$Repeat_Length,
                           freq_tbl$NumReads,
                           bw = allele_bw)
 
 # ── Write outputs ────────────────────────────────────────────────────────
 dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
 base_name <- sub("\\.(fastq|fq)(\\.gz)?$", "", basename(fastq_path))
 base_name <- sub("\\.gz$", "", base_name)
 
 tool_name   <- paste0(triplet_label, "_TRIPLETLOGGER")
 sample_text <- paste0(tool_name, " v", VERSION, " : ", sample_id)
 
 # ── Build metric column names using config-specific label ────────────────
 col_reads_base <- paste0("Reads_", base_threshold_lab)
 col_pct_110    <- paste0("Pct_Over_110_of_", base_threshold_lab)
 col_pct_150    <- paste0("Pct_Over_150_of_", base_threshold_lab)
 
 # Metrics CSV
 metrics_df <- data.frame(
  Sample_ID                  = sample_text,
  Triplet_Type               = triplet_label,
  Total_Reads                = total_reads,
  Reads_Passed_QFilter       = reads_passing,
  Reads_Failed_QFilter       = reads_failed,
  Reads_With_Match           = reads_with_match,
  Q_Threshold                = q_threshold,
  N_Bridged_Reads            = n_bridged_reads,
  EndOfRead_Censored_Reads   = censored_reads,
  Threshold_Rejected_Reads   = threshold_rejected,
  Type_A_Rate                = type_a_rate,
  Max_Type_A_Floor           = max_type_a,
  Max_Type_B                 = max_type_b,
  Min_Repeat_Length          = min_repeat_length,
  Mode_Repeat_Size           = mode_repeat,
  Longest_Repeat             = longest_repeat,
  stringsAsFactors           = FALSE
 )
 
 # ── Add config- and readtype-specific threshold columns ─────────────────
 # Reads dropped by library-prep subtraction (0 for CAG)
 metrics_df[["Reads_Dropped_By_Adjustment"]] <- reads_dropped_adj
 
 # Base threshold reads (always reported)
 metrics_df[[col_reads_base]] <- reads_over_base
 
 # Reads >110 (always reported)
 metrics_df[["Reads_Over_110"]] <- reads_over_110
 
 # Reads >150: omitted for short reads
 if (read_type != "short") {
  metrics_df[["Reads_Over_150"]] <- reads_over_150
 }
 
 # Pct >110 of base (always reported)
 metrics_df[[col_pct_110]] <- ifelse(is.na(pct_over_110), NA,
                                     round(pct_over_110, 2))
 
 # Pct >150 of base: omitted for short reads
 if (read_type != "short") {
  metrics_df[[col_pct_150]] <- ifelse(is.na(pct_over_150), NA,
                                      round(pct_over_150, 2))
 }
 
 # Allele call columns
 metrics_df$Estimated_Allele_Call         <- allele_df$Call[1]
 metrics_df$Estimated_Allele_1_Peak       <- allele_df$Peak_Repeat[1]
 metrics_df$Estimated_Allele_1_Proportion <- allele_df$Proportion[1]
 metrics_df$Estimated_Allele_2_Peak       <- ifelse(nrow(allele_df) >= 2,
                                                    allele_df$Peak_Repeat[2], NA)
 metrics_df$Estimated_Allele_2_Proportion <- ifelse(nrow(allele_df) >= 2,
                                                    allele_df$Proportion[2], NA)
 
 # ── File tag includes read_type suffix ───────────────────────────────────
 file_tag <- paste0(triplet_label, "TRIPLETLogger.v", VERSION, ".", read_type)
 
 metrics_csv    <- file.path(output_dir,
                             paste0(base_name, "_", file_tag, "_RepeatMetrics.csv"))
 raw_freq_csv   <- file.path(output_dir,
                             paste0(base_name, "_", file_tag, "_NumReadsPerRepeat.csv"))
 
 write.csv(metrics_df,  metrics_csv,   row.names = FALSE)
 write.csv(raw_freq,    raw_freq_csv,  row.names = FALSE)
 
 # ── Write rejected match sizes (always written, even if 0 rejected) ─────
 rej_csv <- file.path(output_dir,
                      paste0(base_name, "_", file_tag,
                             "_ThresholdRejected.csv"))
 if (threshold_rejected > 0L) {
  rej_tbl <- as.data.frame(table(rejected_triplet_sizes), stringsAsFactors = FALSE)
  names(rej_tbl) <- c("Rejected_Match_Triplets", "NumReads")
  rej_tbl$Rejected_Match_Triplets <- as.integer(rej_tbl$Rejected_Match_Triplets)
  rej_tbl <- rej_tbl[order(rej_tbl$Rejected_Match_Triplets), ]
 } else {
  rej_tbl <- data.frame(Rejected_Match_Triplets = integer(0),
                        NumReads = integer(0),
                        stringsAsFactors = FALSE)
 }
 write.csv(rej_tbl, rej_csv, row.names = FALSE)
 message(sprintf("[%s] Rejected match sizes written to: %s",
                 sample_id, rej_csv))
 
 # ── Summary messages ─────────────────────────────────────────────────────
 message(sprintf(
  "\n[%s] Predicted Allele call : Peaks %s",
  sample_id,
  paste(allele_df$Peak_Repeat, collapse = ", ")))
 
 if (reads_over_base > 0) {
  if (read_type == "short") {
   message(sprintf(
    "[%s] Reads %s: %d | >110: %d (%.2f%% of %s)",
    sample_id, base_threshold_lab, reads_over_base,
    reads_over_110, pct_over_110, base_threshold_lab))
  } else {
   message(sprintf(
    "[%s] Reads %s: %d | >110: %d (%.2f%% of %s) | >150: %d (%.2f%% of %s)",
    sample_id, base_threshold_lab, reads_over_base,
    reads_over_110, pct_over_110, base_threshold_lab,
    reads_over_150, pct_over_150, base_threshold_lab))
  }
 } else {
  message(sprintf("[%s] Reads %s: 0 | >110: 0 | >150: 0",
                  sample_id, base_threshold_lab))
 }
 
 message(sprintf("[%s] Raw frequency range: %d–%d",
                 sample_id, effective_min, effective_max))
 message(sprintf("[%s] Read type: %s", sample_id, read_type))
 message(sprintf("[%s] Output: %s/", sample_id, output_dir))
 
 invisible(list(
  sample_id              = sample_id,
  triplet_type           = triplet_label,
  read_type              = read_type,
  metrics                = metrics_df,
  raw_freq               = raw_freq,
  n_bridged_reads        = n_bridged_reads,
  censored_reads         = censored_reads,
  threshold_rejected     = threshold_rejected,
  rejected_triplet_sizes = rejected_triplet_sizes
 ))
}

# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

option_list <- list(
 make_option(c("-i", "--input"), type = "character", default = NULL,
             help = "Input FASTQ file (or .gz) [required]",
             metavar = "FILE"),
 make_option(c("-o", "--output"), type = "character", default = ".",
             help = "Output directory [default: current dir]",
             metavar = "DIR"),
 make_option(c("-t", "--tripletType"), type = "character", default = "CAG",
             help = "Triplet repeat type: CAG or CTG [default: CAG]",
             metavar = "TYPE"),
 make_option(c("-r", "--readType"), type = "character", default = "long",
             help = "Read type: long or short [default: long]",
             metavar = "TYPE"),
 make_option(c("-q", "--qThreshold"), type = "integer", default = 20L,
             help = "Minimum mean base-quality to keep a read [default: 20]",
             metavar = "INT"),
 make_option(c("-m", "--minRepeatLength"), type = "integer", default = 1L,
             help = "Minimum repeat length to report [default: 1]",
             metavar = "INT"),
 make_option(c("--maxTypeA"), type = "integer", default = 5L,
             help = "Fixed floor for max Type A errors [default: 5]",
             metavar = "INT"),
 make_option(c("--maxTypeB"), type = "integer", default = NULL,
             help = "Max consecutive non-target (Type B) errors [default: 1 for CAG, 2 for CTG]",
             metavar = "INT"),
 make_option(c("--typeArate"), type = "double", default = 0.05,
             help = "Per-triplet error tolerance for adaptive Type A [default: 0.05]",
             metavar = "FLOAT"),
 make_option(c("--chunkSize"), type = "integer", default = 20000L,
             help = "FASTQ streaming chunk size [default: 20000]",
             metavar = "INT"),
 make_option(c("--alleleBW"), type = "double", default = 1.5,
             help = "Kernel density bandwidth for allele calling [default: 1.5]",
             metavar = "FLOAT"),
 make_option(c("--freqRangeMin"), type = "integer", default = 1L,
             help = "Minimum repeat length for raw frequency table [default: 1]",
             metavar = "INT"),
 make_option(c("--freqRangeMax"), type = "integer", default = 1000L,
             help = "Maximum repeat length for raw frequency table [default: 1000; short reads capped at 120]",
             metavar = "INT"),
 make_option(c("-s", "--sampleId"), type = "character", default = NULL,
             help = "Sample ID override [default: derived from filename]",
             metavar = "STRING")
)

parser <- OptionParser(
 usage       = "usage: %prog -i <input.fastq.gz> [options]",
 option_list = option_list,
 description = paste0(
  "Triplet Logger (Long Read) v", VERSION,
  " — Alignment-free estimation of CAG/CTG repeats from long-read (ONT) amplicon data.")
)

# Only parse when run from the command line (not when source()'d)
if (!interactive()) {
 opts <- parse_args(parser)
 
 if (is.null(opts$input)) {
  print_help(parser)
  stop("Error: --input (-i) is required.", call. = FALSE)
 }
 
 # Validate triplet type
 triplet_type <- toupper(opts$tripletType)
 if (!triplet_type %in% c("CAG", "CTG")) {
  stop("Error: --tripletType must be CAG or CTG.", call. = FALSE)
 }
 
 # Validate read type
 read_type <- tolower(opts$readType)
 if (!read_type %in% c("long", "short")) {
  stop("Error: --readType must be long or short.", call. = FALSE)
 }
 
 config <- get_triplet_config(triplet_type)
 
 # Resolve max_type_b: use config default if not explicitly provided
 max_type_b <- if (is.null(opts$maxTypeB)) config$default_max_type_b else opts$maxTypeB
 
 message(sprintf("\n%s TRIPLETLOGGER v%s  [triplet: %s | read type: %s]",
                 config$label, VERSION, config$label, read_type))
 message(sprintf("Regex: %s", config$pattern))
 if (isTRUE(config$count_all)) {
  message("Counting mode: total triplets (including substitutions)")
 } else {
  message("Counting mode: pure target triplets only")
 }
 if (isTRUE(config$truncate_type_b)) {
  message(sprintf("Type B mode: truncate at %d-th consecutive-error pair\n",
                  max_type_b))
 } else {
  message(sprintf("Type B mode: reject match if consecutive errors > %d\n",
                  max_type_b))
 }
 
 library(tictoc)
 tic()
 
 process_fastq_file(
  fastq_path        = opts$input,
  config            = config,
  sample_id         = opts$sampleId,
  q_threshold       = opts$qThreshold,
  min_repeat_length = opts$minRepeatLength,
  max_type_a        = opts$maxTypeA,
  max_type_b        = max_type_b,
  type_a_rate       = opts$typeArate,
  chunk_size        = opts$chunkSize,
  allele_bw         = opts$alleleBW,
  freq_range_min    = opts$freqRangeMin,
  freq_range_max    = opts$freqRangeMax,
  read_type         = read_type,
  output_dir        = opts$output
 )
 
 toc()
}