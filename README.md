# TripletLogger v1.0

**Alignment-free estimation of CAG / CTG repeat lengths from targeted amplicon sequencing data (Illumina MiSeq short reads or Oxford Nanopore long reads).**

TripletLogger scans FASTQ reads directly with an N-aware regular expression, counts pure target triplets within each accepted match, and reports per-read repeat lengths, a raw repeat-length frequency table, summary metrics, and an estimated allele call (kernel-density peak picking). No alignment to a reference is required.

---

## Features

- Supports both **CAG** and **CTG** repeats with triplet-specific logic (regex, counting rule, flank trimming, Type B handling).
- Works on **short reads** (MiSeq) and **long reads** (ONT) via a single `--readType` switch with sensible defaults for each.
- Tolerates sequencing errors through configurable **Type A** (total non-target) and **Type B** (consecutive non-target) thresholds, with a length-adaptive Type A floor.
- Trims absorbed downstream flank tracts (CCG-like for CAG, CGG-like for CTG).
- Streams FASTQ in chunks — memory-safe for large ONT runs.
- Outputs per-sample CSVs: repeat metrics, reads-per-repeat-length, and rejected-match sizes.

## Requirements

- R ≥ 4.0
- Bioconductor: `ShortRead`, `Biostrings`
- CRAN: `dplyr`, `optparse`, `tictoc`

Install:

```r
install.packages(c("dplyr", "optparse", "tictoc"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ShortRead", "Biostrings"))
```

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

## Examples

CAG repeats from a Nanopore run:

```bash
Rscript TripletLogger_v1_0.R -i sample01.fastq.gz -o results/ -t CAG -r long
```

CTG repeats from a MiSeq run:

```bash
Rscript TripletLogger_v1_0.R -i sample02.fastq.gz -o results/ -t CTG -r short
```

## Output

For each input FASTQ, TripletLogger writes three CSVs to the output directory:

- `<sample>_<TYPE>TRIPLETLogger.v1.0.<readtype>_RepeatMetrics.csv` — summary metrics and estimated allele call(s).
- `<sample>_<TYPE>TRIPLETLogger.v1.0.<readtype>_NumReadsPerRepeat.csv` — raw repeat-length frequency distribution.
- `<sample>_<TYPE>TRIPLETLogger.v1.0.<readtype>_ThresholdRejected.csv` — sizes of matches rejected by error thresholds.


## Version
v1.0

## License
MIT