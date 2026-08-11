# TripletLogger v1.0

**Alignment-free estimation of CAG / CTG repeat lengths from targeted amplicon sequencing data (Illumina MiSeq short reads or Oxford Nanopore long reads).**

Github page- https://github.com/hchetia/TripletLogger

TripletLogger scans FASTQ reads directly with an N-aware regular expression, counts pure target triplets within each accepted match, and reports per-read repeat lengths, a raw repeat-length frequency table, summary metrics, and an estimated allele call (kernel-density peak picking). No alignment to a reference is required.

---
## If you find the tool useful, please consider citing us:
Chetia, H., Kus, L., Sipos, E., Heintz, N. and Mätlik, K. (2026), Expanded ATXN3 CAG Repeat is Stable in Human Purkinje Cells. Mov Disord. https://doi.org/10.1002/mds.70465

---
## Features

- Supports both **CAG** and **CTG** repeats with triplet-specific logic (regex, counting rule, flank trimming, Type B handling).
- Successfully tested on both **CAG** and **CTG** repeats from human _HTT_ exon1 and _ATXN3_ exon10 targeted sequencing data. Any other loci remain untested at the moment.
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

 <details><summary><b>If you are a Macintosh user running into compilation errors, click here</b></summary>
If you are setting up this project on macOS and running `renv::restore()` after opening the `.Rproj` file, you may encounter a compilation error related to missing C++ headers. This guide walks through the expected error and how to fix it.

---

### Expected Error

After running `renv::restore()`, you may see an error like this:

```
Error: Error installing package 'Rcpp':
================================
* installing *source* package 'Rcpp' ...
** libs
using C++ compiler: 'Apple clang version 21.0.0 (clang-2100.1.1.101)'
using SDK: 'MacOSX26.5.sdk'
clang++ -arch x86_64 -std=gnu++17 ...
In file included from api.cpp:26:
In file included from ../inst/include/Rcpp.h:27:
In file included from ../inst/include/RcppCommon.h:30:
In file included from ../inst/include/Rcpp/r/headers.h:67:
../inst/include/Rcpp/platform/compiler.h:37:10: fatal error: 'cmath' file not found
   37 | #include <cmath>
      |          ^~~~~~~
1 error generated.
make: *** [api.o] Error 1
ERROR: compilation failed for package 'Rcpp'
```

### What Is Most Likely Causing This

This is a **macOS SDK / Xcode Command Line Tools (CLT) misconfiguration** — not a problem with R or renv. The C++ standard library header `<cmath>` lives inside the macOS SDK, and clang cannot resolve the path to it. This commonly happens after a major macOS upgrade where the CLTs are not fully re-initialised.

---

### Fix

Work through the steps below in order, stopping once `renv::restore()` succeeds.

#### Step 1 — Reinstall Command Line Tools

Open **Terminal** and run:

```bash
sudo rm -rf /Library/Developer/CommandLineTools
xcode-select --install
```

A GUI prompt will appear — follow it to complete the installation (this takes a few minutes). Afterwards, confirm it worked:

```bash
xcode-select -p        # should return a path to the CLTs
xcrun --show-sdk-path  # should return a valid SDK path, e.g. .../MacOSX26.sdk
```

#### Step 2 — Verify C++ Headers Are Reachable

Run this quick smoke test in Terminal:

```bash
echo '#include <cmath>
int main(){ return 0; }' | clang++ -x c++ -
```

If this compiles silently with no errors, the system is fixed. If it still fails, continue to Step 3.

#### Step 3 — Set `SDKROOT` Explicitly for R

R's build environment sometimes does not inherit the correct `SDKROOT` after a major macOS upgrade. Set it permanently in `~/.R/Makevars`:

```bash
mkdir -p ~/.R
echo "SDKROOT=$(xcrun --show-sdk-path)" >> ~/.R/Makevars
```

Verify the file looks correct:

```bash
cat ~/.R/Makevars
# Expected output (path may differ by macOS version):
# SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX26.sdk
```

#### Step 4 — Retry the Restore

Back in RStudio, restart your R session first:

- **Session → Restart R** (or `Cmd + Shift + F10`)

Then run:

```r
renv::restore()
```

This should now complete successfully.

---

### Note for Apple Silicon (M-series) Users

If you are on an M1/M2/M3/M4 Mac, check which build of R you have installed. Compiling with `-arch x86_64` means you are either on Intel or running the x86_64 R build under Rosetta 2. Installing the **native ARM64 build of R** from [mac.r-project.org](https://mac.r-project.org) is recommended — it avoids a broader class of Rosetta-related compilation issues and runs significantly faster. This is not the cause of the `cmath` error, but is worth addressing to prevent other issues down the line.
 </details>

---

<details><summary><b>Reproducing the renv if you are a Windows OS user </b></summary>
 
These steps rebuild the exact set of R packages TripletLogger needs, using the
`renv.lock` file included in this repository. Work through them **in order**.
After each step there's a **✓ Check** describing what success looks like — don't
move on until you see it. If a step fails, jump to **Troubleshooting** at the
bottom and find the matching symptom.

You'll need permission to install software on the PC.

---

### 1. Install the prerequisites

Install these three, in this order:

1. **R 4.4.x** (not the newest release). The lockfile targets Bioconductor 3.20, which pairs with R 4.4.
 - Go to <https://cran.r-project.org/bin/windows/base/>. If the top link is a 4.4.x version, use it. If it's newer (4.5+), click **Previous releases** and
pick the latest **4.4.x**. Run the installer with all defaults.
2. **Rtools44** — the compiler toolkit some packages need.
- Download and run the installer from
<https://cran.r-project.org/bin/windows/Rtools/rtools44/>. Accept all defaults.
3. **RStudio** — the window you'll type commands into.
- Download and run from <https://posit.co/download/rstudio-desktop/>.
✓ **Check:** the Start menu lists **R x64 4.4.x** and **RStudio**.

---

### 2. Get the project and open it in RStudio

1. On this repository's GitHub page, click the green **Code** button →
**Download ZIP**.
2. Right-click the ZIP → **Extract All**. Put it somewhere simple, e.g.
`C:\Users\<you>\Documents\`.
3. In RStudio: **File → Open Project…** and open the project folder. If the repo
contains a `.Rproj` file, open that. If it doesn't, instead run this in the Console (use **your** path, with forward slashes):

```r
setwd("C:/Users/<you>/Documents/TripletLogger-main/TripletLogger-main")
```
✓ **Check:** when the project loads, the Console prints a line mentioning **renv** (e.g. *"Project ... loaded. [renv ...]"*). If it offers to install
renv, allow it.

**This matters:** every command below must be run with this project open. renv keeps packages in a *project-local* library, so commands like `library(...)` only see the project's packages when the project is active.

---
 
 ### 3. Recreate the environment
 
 Run one command in the Console:
 
 ```r
renv::restore()
```

If it asks *"Do you want to proceed?"*, type `y` and press Enter. This downloads
and installs everything in `renv.lock` and can take several minutes. Let it finish.

✓ **Check:** it ends **without** a red `Error:` line. If you got an error, go to
**Troubleshooting**.

**Golden rule:** never fix a missing package with `install.packages(...)` while this project is open. Always use `renv::restore()` or `renv::install(...)`. Mixing the two is the most common way this setup breaks.

---
 
 ### 4. Confirm it worked
 
 With the project still open, run:
 
 ```r
renv::status()
library(ShortRead)
```

✓ **Check:**
 - `renv::status()` reports the project is **consistent / in sync** (no long list of mismatches).
- `library(ShortRead)` loads with **no error** (startup messages are normal).

If `library(ShortRead)` says *"there is no package called 'ShortRead'"*, you are
almost certainly not in the project — see Troubleshooting **C**.

---
 
 ### 5. Run TripletLogger (from a terminal, not the R Console)
 
 The program is launched from a Windows terminal. Typing `Rscript ...` into the R
Console gives an `unexpected symbol` error — that's the wrong place.

1. In RStudio, click the **Terminal** tab (next to **Console**), or open **Command Prompt** from the Start menu.
2. Move into the project folder and run it (adjust the input file name):

```bash
cd C:\Users\<you>\Documents\TripletLogger-main\TripletLogger-main
Rscript TripletLogger_v1_0.R -i your_reads.fastq -o results/ -t CAG -r short
```

If Windows says **`Rscript` is not recognized**, use the full path: 
For example, if your R version is 4.4.2, do
```bash
"C:\Program Files\R\R-4.4.2\bin\Rscript.exe" TripletLogger_v1_0.R -i your_reads.fastq -o results/ -t CAG -r short
```

✓ **Check:** a `results/` folder appears with output files in it.

---

### Troubleshooting

Find your symptom below.

#### A. `restore()` stops with a download error on a Bioconductor package
*Example: it fails on `BiocParallel`, and packages that depend on it
(`Rsamtools`, `GenomicAlignments`, `ShortRead`) then report "dependency failed".*

**Cause:** the lockfile pins a Bioconductor package to a patch version that the
Bioconductor 3.20 server no longer hosts, so the download 404s and `restore()`
stops.

**Fix:** install that one package through renv (which fetches the version
currently available), update the lockfile to match, then finish the restore.
Replace `BiocParallel` with whatever package name appeared in the error:

```r
renv::install("bioc::BiocParallel")   # the package named in the error
renv::snapshot()                       # type y when prompted — updates renv.lock
renv::restore()                        # installs the remaining packages
```

Repeat for any other Bioconductor package that fails the same way. After this,
`renv::status()` should report the project in sync.

This deliberately updates `renv.lock` to the versions that actually installed.
That's the right trade-off for getting it running; it means your environment may differ by a patch version from the original lockfile.

#### B. Errors mentioning *"compilation failed"* or *"package not available"* during install
**Cause:** a package had to be built from source and the compiler toolkit
(Rtools44) is missing or not detected.

**Fix:**
 1. Re-run the **Rtools44** installer from Step 1, accept defaults, then fully
close and reopen RStudio (so it picks up the toolchain).
2. Re-run `renv::restore()`.
3. (Optional check) if you have the `pkgbuild` package, you can confirm the
toolchain is visible to R with:
 ```r
pkgbuild::has_build_tools(debug = TRUE)
```
A result of `TRUE` means the compiler is set up correctly. (Don't rely on
   searching for `make` on the PATH — R can find Rtools without it, so that test
   can report a false problem.)

#### C. `library(...)` says *"there is no package called ..."* even though install succeeded
**Cause:** you're running in a plain R session, not inside the TripletLogger project, so R is looking in the wrong (system) library.
**Fix:** reopen the project (**File → Open Project…**, or rerun the `setwd(...)` from Step 2) and confirm the Console shows the **renv** activation line, then try again. Everything must be run with the project active.

#### D. `Rscript TripletLogger_v1_0.R ...` gives *"unexpected symbol"*
**Cause:** you typed a terminal command at the R Console `>` prompt.
**Fix:** run it in the **Terminal** tab or Command Prompt instead — see Step 5.

#### E. A warning that renv X.Y.Z was loaded but the project records a different version
*Example: "renv 1.2.3 was loaded ... project is configured to use renv 1.1.6".*
**Cause:** harmless version drift between your installed renv and the one recorded in the lockfile.
**Fix:** safe to ignore. To silence it, record the version you actually have:
```r
renv::record(paste0("renv@", as.character(packageVersion("renv"))))
```

### Three rules to remember
1. **Match the R version** to the project (R 4.4.x here) — don't grab the newest.
2. **Inside the renv project, use only `renv::restore()` / `renv::install()`** —
   never `install.packages()`.
3. **Run the tool from a terminal, not the R Console.**

---
</details>


## Don't want to use renv?
That's fine — you can skip all of the above and simply install the packages yourself with the commands in [Requirements](## Requirements). The tool will run exactly the same way; you just won't have the versions locked down, so a future package update could, in principle, change its behavior.

## Running Triplet Logger interactively in RStudio
 
The `Rscript` command above is the normal way to run TripletLogger. If you would rather work inside an interactive R session — the R console or RStudio — you must **load the script's functions first** with `source()`. Calling a function such as `process_fastq_file()` directly, without sourcing the script, produces:
 
```
Error: could not find function "process_fastq_file"
```
 
This happens because the functions only exist once the script is sourced. The command-line block at the bottom of the script is intentionally skipped in an interactive session (it is guarded by `if (!interactive())`), so sourcing the script does **not** trigger the "`--input` is required" error — it simply makes every function available to you.
 
**Step 1 — Set the working directory to the project folder** (the folder containing `TripletLogger_v1_0.R` and `TripletLogger_histogram.R`). This matters: in an interactive session the script sources the histogram module from the working directory, so the histogram is only generated when R is started in — or pointed at — the project folder.
 
```r
setwd("/path/to/TripletLogger")
```
 
**Step 2 — Load the functions.**
 
```r
source("TripletLogger_v1_0.R")
```
 
**Step 3 — Build the triplet configuration.** On the command line the `-t`/`--tripletType` flag does this for you. In the console you build it yourself and pass it to `process_fastq_file()` as the `config` argument, which has no default:
 
```r
config <- get_triplet_config("CAG")   # or "CTG"
```
 
**Step 4 — Run the tool.** Every other argument has the same default as the corresponding command-line flag, so you only set the ones you want to change:
 
```r
process_fastq_file(
  fastq_path = "sample01.fastq.gz",
  config     = config,
  read_type  = "short",      # "long" (ONT, the default) or "short" (MiSeq)
  output_dir = "results"
)
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
