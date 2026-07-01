# Replication of the Horvath (2013) pan-tissue epigenetic clock on public blood methylation data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21090241.svg)](https://doi.org/10.5281/zenodo.21090241)

Status: first reproducible replication pass completed.

This repository is an open-science replication study. The goal is to
independently run the Horvath (2013) pan-tissue epigenetic clock on a public
whole-blood methylation dataset and report whether the published result holds
under this pipeline.

## Scientific question

Can the Horvath (2013) pan-tissue DNA methylation clock replicate on public
whole-blood methylation data from GSE40279 / Hannum 2013 with accuracy close to
the originally reported blood performance?

Note: this study applies the Horvath (2013) 353-CpG pan-tissue clock to the
GSE40279 dataset. GSE40279 is the dataset from Hannum et al. (2013); we are NOT
using the separate 71-CpG Hannum clock here. Horvath's clock was validated on
blood, which is why it is appropriate for this whole-blood dataset.

## Dataset

Primary dataset:

- GEO accession: `GSE40279`
- Study: Hannum et al. 2013
- Tissue: whole blood
- Platform: Illumina HumanMethylation450
- Approximate sample count: 656
- Age range: approximately 19 to 101 years

Data handling:

- Raw GEO downloads and derived cache files are stored under gitignored
  `data/raw/`.
- `make data` downloads/parses the public GEO data when needed and writes a
  derived cache for repeat runs.
- The signal-intensity file is not required for the beta-value clock path.
- Sample metadata is parsed from the family SOFT file; age parsing fails loudly
  if a clean age field cannot be found.
- The analysis uses a memory-conscious 353-CpG beta subset matching the Horvath
  coefficient list.

Fallback datasets to document if needed:

- `GSE55763`
- `GSE69914`

## Replication criterion, declared before results

This project treats the replication as successful if the computed median
absolute error is within approximately 1 year of the published reference
performance for blood and the correlation remains high.

Reference target to compare against:

- Horvath 2013 blood performance: approximately MAE 3.6 years, Pearson r above
  0.9.

## Current headline result

Initial hand-implementation result on the local GSE40279 Horvath-353 beta
subset:

- Median absolute error: 3.856 years
- Pearson r: 0.9183
- RMSE: 6.322 years
- n samples: 656

Under the pre-declared criterion, this initial hand implementation replicates
the reference result.

Lightweight independent cross-check against Biolearn's bundled Horvathv1
coefficient data:

- n compared: 656
- mean absolute difference vs our predictions: 1.083e-17 years
- max absolute difference vs our predictions: 3.553e-15 years
- Pearson r between prediction vectors: 1.0

## Reproduce

Create or update the conda environment:

```bash
conda env create -f environment.yml
conda activate epi-clock-replication
```

Prepare data, run analysis, and run the independent cross-check:

```bash
make all
```

Individual targets:

```bash
make data       # download/parse/cache data if needed
make analyze    # run hand implementation from cache
make crosscheck # run lightweight Biolearn coefficient-data cross-check
make run        # analyze + crosscheck
make clean      # remove generated results and derived cache
```

The first `make all` run can be slow because it parses the large GEO SOFT file.
Repeat runs use the derived cache and are faster.

Docker target:

```bash
docker build -t epigenetic-clock-replication .
docker run --rm epigenetic-clock-replication
```

Note: the Dockerfile is provided for convenience but has not yet been
build-verified on this machine (local Docker/WSL not set up). Treat the Docker
path as untested.

## Repository layout

```text
epigenetic-clock-replication/
  README.md
  LICENSE
  environment.yml
  requirements.txt
  Dockerfile
  Makefile
  data/
    raw/                         # gitignored public GEO downloads
    coefficients/                # committed clock coefficients go here
  src/
    download_data.py
    horvath_clock.py
    analyze.py
    crosscheck.py
  results/                       # generated outputs
```

## Coefficient source

Working coefficient source:

- `data/coefficients/horvath2013_coefficients.csv`
- copied from Biolearn 0.9.1 package data file `biolearn/data/Horvath1.csv`
- Biolearn license reported by package metadata: `new BSD`
- Biolearn model definition: `Horvathv1`, year 2013, source URL:
  `https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115`
- Original Horvath 2013 article and supplementary coefficient files are
  attributed to Horvath 2013. The article is distributed under the Creative
  Commons Attribution License 2.0, requiring proper citation of the original
  work.

The file contains:

- 353 CpG IDs
- one coefficient per CpG
- one intercept

Transform convention:

- Biolearn stores 353 CpG coefficients without an intercept row.
- Its model definition applies `anti_trafo(sum + 0.696)`.
- This repository represents `+0.696` as an explicit `intercept` row.
- Training-set means are not provided in this coefficient source, so missing CpG
  handling must be logged in the final analysis.

## Current outputs

Generated by `python src/analyze.py`:

- `results/metrics.csv`
- `results/REPLICATION_NOTE.md`
- `results/predicted_vs_chronological_age.svg`

Generated by `python src/crosscheck.py`:

- `results/crosscheck_biolearn_lightweight.csv`

Cache verification:

- `CACHE_REPORT.md`
- cache equivalence passed with exact metadata comparison and exact beta matrix
  comparison (`atol=0.0`, `rtol=0.0`)
- cold `make all`: 296.111 seconds
- warm `make all`: 80.578 seconds
- warm `make data`: 4.738 seconds, prints `cache present, skipping parse`

Cross-check note:

- Importing `biolearn.model` requires `torch`.
- Instead of installing the heavy dependency immediately, this repository reads
  Biolearn's bundled `Horvath1.csv` directly and applies the Biolearn model
  convention `anti_trafo(sum + 0.696)` in independent code.
- This verifies the coefficient/intercept/transform/orientation path while
  avoiding a large dependency install.

Plot/runtime note:

- Some native numerical/plotting stacks can be brittle across local platforms.
  The implementation uses explicit multiply+sum, closed-form simple-regression
  slope, and dependency-light SVG output.
- This is an implementation workaround, not a change in the Horvath clock
  formula.

## Implementation notes

The clock implementation must preserve two common failure points:

1. Missing CpGs must be logged and imputed transparently.
2. Horvath's transformed-age scale must be inverted correctly.

The age transformation uses `adult.age = 20`:

```text
F(age) = log(age + 1) - log(adult.age + 1), if age <= adult.age
F(age) = (age - adult.age) / (adult.age + 1), if age > adult.age
```

The linear predictor is on the transformed scale and must be inverse-transformed
back to years.

## Limitations and what would falsify this

This replication can fail or partially replicate for legitimate reasons:

- coefficient-source mismatch
- incorrect or unavailable sample age metadata
- No BMIQ / adjusted-BMIQ normalization is applied; beta values are used as
  provided by GSE40279. The original Horvath pipeline used BMIQ normalization,
  so this is a methodological difference from the original.
- Only the 353-CpG Horvath subset is loaded, not the full methylation matrix;
  missing-CpG handling is logged at runtime.
- missing 353-clock CpGs
- batch effects
- tissue/platform differences
- age-transform implementation mistakes

Falsification conditions:

- age metadata cannot be parsed robustly
- the coefficient source cannot be validated
- the hand implementation disagrees materially with an independent
  implementation after using the same data and coefficients
- the computed metrics fall clearly outside the declared replication criterion

## Citations

- Horvath, S. (2013). DNA methylation age of human tissues and cell types.
  Genome Biology, 14(10), R115.
  https://doi.org/10.1186/gb-2013-14-10-r115
- Hannum, G., Guinney, J., Zhao, L., et al. (2013). Genome-wide methylation
  profiles reveal quantitative views of human aging rates. Molecular Cell,
  49(2), 359-367. https://doi.org/10.1016/j.molcel.2012.10.016
- Gene Expression Omnibus accession GSE40279 (data source for this
  replication).
- Original Horvath clock coefficients are (c) 2013 Steve Horvath, released
  under CC BY 2.0; obtained here via the Biolearn 0.9.1 bundled data file
  biolearn/data/Horvath1.csv (Biolearn license: new BSD).
