# TripletLogger v1.0
An R-based analytical tool for alignment-free estimation of CAG and CTG trinucleotide repeats from targeted amplicon MiSeq data.

## Overview
TripletLogger scans FASTQ reads for trinucleotide repeat tracts using a fuzzy regular expression matching approach. It bypasses the need to align, sort and index reads to a reference genome for genotyping. The tool retains the longest valid repeat match per read and generates a frequency distribution of allele sizes; allele calls based on that distribution; and a set of metrics, namely, Repeat Metrics, that provides a crucial quantification "Percent Reads over 110". This metric is a unique quantification of the percentage of long repeats that lies in the edge of the miseq assay limit and cannot be quantified using alignment-based approaches like BWA MEM.

TripletLogger currently handles both **CAG** and **CTG** repeats. The `--tripletType` flag selects the mode, which controls the regex, counting method, error handling, and flank-trimming behavior.

---

## Requirements

**R version:** ≥ 4.0

**R packages:**

| Package | Source | Purpose |
|---|---|---|
| ShortRead | Bioconductor | Streaming FASTQ reader |
| Biostrings | Bioconductor | Sequence utilities |
| dplyr | CRAN | Data manipulation |
| optparse | CRAN | Command-line argument parsing |
| tictoc | CRAN | Runtime reporting |

### Installation

```r
# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ShortRead", "Biostrings"))

# CRAN packages
install.packages(c("dplyr", "optparse", "tictoc"))
```

---

## Quick start

### Step 1 — Prepare your input

TripletLogger accepts a single FASTQ file per run (`.fastq`, `.fq`, or gzip-compressed `.fastq.gz` / `.fq.gz`). These should be reads from a targeted amplicon sequencing experiment containing CAG or CTG repeat tracts.

### Step 2 — Run the script

**CAG repeats (default):**

```bash
Rscript TripletLogger_V2_1d.R \
  -i sample.fastq.gz \
  -o results/
```

**CTG repeats:**

```bash
Rscript TripletLogger_V2_1d.R \
  -i sample.fastq.gz \
  -t CTG \
  -o results/
```

### Step 3 — Examine the output

Output CSVs are written to the directory specified by `-o`. See the [Output files](#output-files) section below.

---

## Command-line options

| Flag | Long form | Type | Default | Description |
|---|---|---|---|---|
| `-i` | `--input` | string | *(required)* | Path to input FASTQ file |
| `-o` | `--output` | string | `.` | Output directory |
| `-t` | `--tripletType` | string | `CAG` | Repeat type: `CAG` or `CTG` |
| `-q` | `--qThreshold` | integer | `30` | Minimum mean Phred quality per read |
| `-m` | `--minRepeatLength` | integer | `10` | Minimum repeat length to include in filtered output |
| `-s` | `--sampleId` | string | *(from filename)* | Override the sample ID label |
| | `--maxTypeA` | integer | `5` | Fixed floor for maximum allowed Type A errors |
| | `--maxTypeB` | integer | `1` (CAG) / `2` (CTG) | Type B error threshold (see [Error handling](#error-handling)) |
| | `--typeArate` | float | `0.05` | Per-triplet error tolerance for adaptive Type A scaling |
| | `--chunkSize` | integer | `50000` | Number of reads per streaming chunk |
| | `--alleleBW` | float | `1.5` | Kernel density bandwidth for allele peak detection |
| | `--freqRangeMin` | integer | `10` | Lower bound of the raw frequency table range |
| | `--freqRangeMax` | integer | `120` | Upper bound of the raw frequency table range |

---

## How it works

### Step-by-step pipeline

1. **Streaming read input** — Reads are loaded in chunks (default 50,000) to limit memory usage.

2. **Quality filtering** — Reads whose mean Phred quality falls below `--qThreshold` are discarded.

3. **Regex matching** — Each read is searched on both the forward and reverse-complement strand using a mode-specific N-aware regex:
   - CAG: `CAG(CAG|CA[ATCGN]|C[ATCGN]G|[ATCGN]AG)+`
   - CTG: `CTG(CTG|CT[ATCGN]|C[ATCGN]G|[ATCGN]TG)+`

4. **Flank-tract trimming** — Runs of ≥2 consecutive flank-like triplets (CCG-like for CAG; CGG-like for CTG) are trimmed from the end of the match to avoid counting absorbed flanking sequence.

5. **Error filtering** — Matches are evaluated against Type A and Type B error thresholds (see below). The behaviour differs by mode.

6. **Repeat counting** — The repeat length for the read is determined from the best passing match.

7. **Allele calling** — A kernel density estimate is fitted to the repeat-length distribution. The top one or two peaks are reported, and the sample is classified as homozygous or heterozygous.

8. **Output** — Frequency tables, metrics, allele calls, and diagnostics are written as CSVs.

### Error handling

Sequencing errors within the repeat tract produce triplets that are not the exact target motif. TripletLogger classifies these into two types:

- **Type A** — Total count of non-target triplets in a match.
- **Type B** — Instances where two non-target triplets occur consecutively.

The Type A threshold scales with tract length to accommodate the higher absolute error count expected in long expanded alleles:

```
effective_max_type_a = max(maxTypeA, floor(n_triplets × typeArate))
```

**CAG mode:** If a match exceeds either the Type A or Type B threshold, the entire match is rejected. The default Type B threshold is 1.

**CTG mode:** The tract is truncated at the position of the *N*-th consecutive-error pair (where *N* = `--maxTypeB`, default 2). The Type A check is then applied to the truncated tract. This preserves the counting behaviour of the original CTG Logger.

### Counting method

| | CAG mode | CTG mode |
|---|---|---|
| **What is counted** | Pure target (`CAG`) triplets only | All triplets (including substitutions) |
| **Flank exclusion** | Subtracts 1 if CAA–CAG pattern found at tract end | Not applied |
| **Example:** `CTG-CTG-CTN-CTG` | — | Reports **4** |
| **Example:** `CAG-CAG-CAN-CAG` | Reports **3** | — |

---

## Output files

All output files are named with the pattern:

```
<sample>_<TYPE>Logger.v2.1d_<suffix>.csv
```

where `<TYPE>` is `CAG` or `CTG`.

| File suffix | Description |
|---|---|
| `_RepeatMetrics.csv` | Summary metrics for the sample: read counts, quality stats, mode, allele calls, threshold parameters |
| `_FrequencyOfRepeats.csv` | Repeat-length frequency table filtered by `--minRepeatLength` (observed lengths only) |
| `_FullDistribution.csv` | Unfiltered frequency table (observed lengths only, no minimum applied) |
| `_RawFrequency.csv` | Zero-filled frequency table spanning a continuous integer range (default 10–120); every size has a row even if frequency is 0 |
| `_AlleleCalls.csv` | Allele peak positions, proportions, and homozygous/heterozygous classification |
| `_ThresholdRejected.csv` | *(conditional)* Triplet counts of matches that were rejected by error thresholds; written only if rejections occurred |

### Raw frequency table

The `_RawFrequency.csv` file provides a continuous range of repeat sizes with zero-filling. This is useful for downstream plotting and cross-sample comparisons where a consistent row structure is needed.

The default range is 10–120. If observed data falls outside this range, the table extends automatically. The bounds can be changed with `--freqRangeMin` and `--freqRangeMax`.

---

## Usage examples

### Basic CAG analysis

```bash
Rscript TripletLogger_V2_1d.R -i HD_sample.fastq.gz -o cag_results/
```

### CTG analysis with custom quality threshold

```bash
Rscript TripletLogger_V2_1d.R \
  -i DM1_sample.fastq.gz \
  -t CTG \
  -q 20 \
  -o ctg_results/
```

### Custom raw frequency range

```bash
Rscript TripletLogger_V2_1d.R \
  -i sample.fastq.gz \
  -t CTG \
  --freqRangeMin 5 \
  --freqRangeMax 200 \
  -o results/
```

### Relaxed error thresholds for long expansions

```bash
Rscript TripletLogger_V2_1d.R \
  -i expanded_sample.fastq.gz \
  --typeArate 0.08 \
  --maxTypeA 8 \
  -o results/
```

### Batch processing in R

```r
source("TripletLogger_V2_1d.R")

config <- get_triplet_config("CTG")

fastq_files <- list.files(path = "data/", pattern = "\\.fastq\\.gz$",
                           full.names = TRUE)

results <- list()
for (f in fastq_files) {
  results[[f]] <- process_fastq_file(
    fastq_path        = f,
    config            = config,
    q_threshold       = 30,
    min_repeat_length = 10L,
    freq_range_min    = 10L,
    freq_range_max    = 120L,
    output_dir        = "batch_results/"
  )
}
```

---

## Interpreting the metrics

The `_RepeatMetrics.csv` file contains a single-row summary per sample. Key fields:

| Field | Meaning |
|---|---|
| `Total_Reads` | Total reads in the FASTQ file |
| `Reads_Passed_QFilter` | Reads passing the mean quality threshold |
| `Reads_With_Match` | Reads containing at least one valid repeat match |
| `Mode_Repeat_Size` | Most frequent repeat length (after filtering) |
| `Allele_Call` | `homozygous`, `heterozygous`, `insufficient_data`, or `no_peaks` |
| `Allele_1_Peak` / `Allele_2_Peak` | Repeat length at each allele peak |
| `Allele_1_Proportion` / `Allele_2_Proportion` | Fraction of reads assigned to each allele |
| `Threshold_Rejected_Reads` | Reads with a regex match that failed error thresholds |
| `N_Bridged_Reads` | Reads containing N bases where a repeat was still detected |
| `EndOfRead_Censored_Reads` | Reads where the repeat tract extends to the end of the read (length may be underestimated) |

---

## Troubleshooting

**No reads pass quality filter** — Lower `-q` (e.g. `-q 20`) or check your sequencing run quality.

**No repeat structures detected** — Verify the correct `--tripletType` for your amplicon. Check that primer/adapter trimming has not removed the repeat tract.

**High threshold-rejected count** — The error thresholds may be too strict for your data. Try increasing `--typeArate` (e.g. `0.08`) or `--maxTypeA`. Review the `_ThresholdRejected.csv` file to see the sizes of rejected matches.

**Many end-of-read censored reads** — Read length may be too short to span the full repeat tract in expanded alleles. The reported lengths for these reads are likely underestimates.

---

## License

See [LICENSE](LICENSE) for details.
