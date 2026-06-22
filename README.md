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
 
R packages get updated over time, and a script that runs perfectly today can behave differently — or stop working — once one of the packages it relies on changes. To prevent that, this repository includes a file called `renv.lock`. Think of it as a recipe card that records the *exact* version of every R package (and the version of R itself) that TripletLogger was built and tested with.
 
The [renv](https://rstudio.github.io/renv/) package can read that recipe and install the matching versions on your computer, so the tool works the same way for you as it does for us. You only need to set this up **once**.
 
**Step 1 — Get the code and find its folder.**
Download or clone this repository. The folder it ends up in is your "project folder" — it's the one containing `renv.lock` and the `TripletLogger_v1_0.R` script.
 
**Step 2 — Open R *inside* that folder.**
This part matters: renv only works when R is started from the project folder. The easiest way is to double-click the `TripletLogger.Rproj` file, which opens RStudio already pointed at the right place. (If you use a terminal instead, `cd` into the folder first, then start R.)
 
**Step 3 — Install renv (skip if you already have it).**
 
```r
install.packages("renv")
```
 
**Step 4 — Install the recorded packages.**
 
```r
renv::restore()
```
 
This downloads the correct version of each package into a private library that belongs to this project only — it won't change or interfere with the packages you use for other R work. The first run can take a few minutes; if it asks you to confirm, type `y` and press Enter.
 
That's the whole setup. From now on, every time you open R from this folder, renv automatically loads the right packages for you (a small hidden `.Rprofile` file in the folder handles this), so you can just run TripletLogger as shown in [Usage](#usage) below.
 
To double-check it worked, run:
 
```r
renv::status()
```
 
If it says the project is "in a consistent state," everything matches the recipe and you're good to go.
 
**Don't want to use renv?** That's fine — you can skip all of the above and simply install the packages yourself with the commands in [Requirements](#requirements). The tool will run exactly the same way; you just won't have the versions locked down, so a future package update could, in principle, change its behavior.

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