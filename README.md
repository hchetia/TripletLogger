# TripletLogger v1.0

**Alignment-free estimation of CAG / CTG repeat lengths from targeted amplicon sequencing data (Illumina MiSeq short reads or Oxford Nanopore long reads).**

Github page- https://github.com/hchetia/TripletLogger

TripletLogger scans FASTQ reads directly with an N-aware regular expression, counts pure target triplets within each accepted match, and reports per-read repeat lengths, a raw repeat-length frequency table, summary metrics, and an estimated allele call (kernel-density peak picking). No alignment to a reference is required.

---

## Features

- Supports both **CAG** and **CTG** repeats with triplet-specific logic (regex, counting rule, flank trimming, Type B handling).
- Works on **short reads** (MiSeq) and **long reads** (ONT) via a single `--readType` switch with sensible defaults for each.
- Tolerates sequencing errors through configurable **Type A** (total non-target) and **Type B** (consecutive non-target) thresholds, with a length-adaptive Type A floor.
- Trims absorbed downstream flank tracts (CCG-like for CAG, CGG-like for CTG).
- Streams FASTQ in chunks — memory-safe for large ONT runs.
- Produces a **summary histogram** of the repeat-length distribution with the estimated allele peak(s) marked (on by default; PNG or SVG).
- Outputs per-sample files: repeat metrics, reads-per-repeat-length, rejected-match sizes, and the histogram image.

## Requirements

- R ≥ 4.4 (the bundled `renv.lock` pins R 4.4.2)
- Bioconductor: `ShortRead`, `Biostrings`
- CRAN: `optparse`, `tictoc`, `ggplot2`

Manual install:

```r
install.packages(c("optparse", "tictoc", "ggplot2"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ShortRead", "Biostrings"))
```

Alternatively, restore the exact pinned environment with renv (see below).

## Reproducible environment (renv)

The repository ships an [`renv`](https://rstudio.github.io/renv/) lockfile that pins the exact CRAN and Bioconductor package versions (and the R version) used for the analysis, so the environment can be reproduced precisely.

From the repository root:

```r
install.packages("renv")   # only if renv is not already available
renv::restore()            # installs the pinned packages into a project-local library
```

`renv::restore()` reads `renv.lock` and installs the recorded packages into a project-local library. The `.Rprofile` in the repository root activates this library automatically when R is started from the project directory, so no further setup is required. To confirm the environment matches the lockfile:

```r
renv::status()             # should report the project is in a consistent state
```

If you prefer not to use renv, the manual install in [Requirements](#requirements) reproduces the same package set, though without version pinning.

## Usage

```bash
Rscript TripletLogger_v1_0.R -i <input.fastq.gz> [options]
```

### Required

| Flag | Description |
|---|---|
| `-i`, `--input` | Input FASTQ file (plain or gzipped) |

### Common options

| Flag | Default | Description |
|---|---|---|
| `-o`, `--output` | `.` | Output directory |
| `-t`, `--tripletType` | `CAG` | Repeat type: `CAG` or `CTG` |
| `-r`, `--readType` | `long` | `long` (ONT) or `short` (MiSeq) |
| `-q`, `--qThreshold` | `20` | Minimum mean base quality per read |
| `-s`, `--sampleId` | filename | Override sample ID |

### Advanced options

| Flag | Default | Description |
|---|---|---|
| `-m`, `--minRepeatLength` | `1` | Minimum repeat length to report |
| `--maxTypeA` | `5` | Fixed floor for Type A errors (total non-target triplets) |
| `--maxTypeB` | `1` (CAG) / `2` (CTG) | Max consecutive non-target triplets |
| `--typeArate` | `0.05` | Per-triplet error tolerance for adaptive Type A |
| `--chunkSize` | `20000` | FASTQ streaming chunk size |
| `--alleleBW` | `1.5` | Kernel density bandwidth for allele calling |
| `--freqRangeMin` | `1` | Minimum repeat length for frequency table |
| `--freqRangeMax` | `1000` | Maximum repeat length (short reads capped at 120) |

### Histogram options

A summary histogram is generated **by default** alongside the CSV outputs. It plots read count against repeat length, marks the estimated allele peak(s), uses x-axis tick labels every 5 repeat units, and reports the triplet type, read type, and percentage of reads over 110 repeats in the subtitle. The y-axis is capped at `--histoReadsMax`; bars above the cap are clipped, but the allele-peak markers are always drawn at their true repeat length.

| Flag | Default | Description |
|---|---|---|
| `--histo` | `ON` | Generate the summary histogram (`ON` or `OFF`) |
| `--histoFreqMin` | `0` | Histogram x-axis minimum repeat length |
| `--histoFreqMax` | auto | Histogram x-axis maximum repeat length (auto-fits the data and always includes the called peaks) |
| `--histoReadsMin` | `0` | Histogram y-axis minimum read count |
| `--histoReadsMax` | `1000` | Histogram y-axis maximum read count |
| `--histoFormat` | `png` | Histogram image format (`png` or `svg`) |

Disable the histogram with `--histo OFF`. To regenerate a histogram from existing TripletLogger output without re-running the full pipeline, the plotting module can be run on its own:

```bash
Rscript TripletLogger_histogram.R \
  --freq    <sample>_<TYPE>_TRIPLETLogger.v1.0.<readtype>_NumReadsPerRepeat.csv \
  --metrics <sample>_<TYPE>_TRIPLETLogger.v1.0.<readtype>_RepeatMetrics.csv
```

`--metrics` is optional (without it, the histogram is drawn without allele-peak markers). The standalone module accepts the same axis and format controls via `--freqMin`, `--freqMax`, `--readsMin`, `--readsMax`, and `--format`, plus `--out` and `--sampleId`. `TripletLogger_histogram.R` must be kept in the same directory as `TripletLogger_v1_0.R`, which sources it automatically.

## Examples

CAG repeats from a Nanopore run:

```bash
Rscript TripletLogger_v1_0.R -i sample01.fastq.gz -o results/ -t CAG -r long
```

CTG repeats from a MiSeq run:

```bash
Rscript TripletLogger_v1_0.R -i sample02.fastq.gz -o results/ -t CTG -r short
```

CAG run with the histogram exported as SVG and the y-axis rescaled:

```bash
Rscript TripletLogger_v1_0.R -i sample03.fastq.gz -o results/ -t CAG -r long \
  --histoFormat svg --histoReadsMax 500
```

## Output

For each input FASTQ, TripletLogger writes three CSVs and one .png file to the output directory:

- `<sample>_<TYPE>_TRIPLETLogger.v1.0.<readtype>_RepeatMetrics.csv` — summary metrics and estimated allele call(s).
- `<sample>_<TYPE>_TRIPLETLogger.v1.0.<readtype>_NumReadsPerRepeat.csv` — raw repeat-length frequency distribution.
- `<sample>_<TYPE>_TRIPLETLogger.v1.0.<readtype>_ThresholdRejected.csv` — repeat sizes of matches rejected by error thresholds.
- `<sample>_<TYPE>_TRIPLETLogger.v1.0.<readtype>_Histogram.png` — histogram of the repeats and their respective read counts (`.svg` instead of `.png` when `--histoFormat svg`).

## Version
v1.0

## License
MIT