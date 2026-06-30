# Cache layer verification report

Date: 2026-06-29 / 2026-06-30 local session

## Scope

Added a derived cache for the GSE40279 Horvath 353-CpG replication pipeline.

Cache files:

- `data/raw/samples_metadata.parquet`
- `data/raw/betas_subset.parquet`

These files are derived cache, not source truth. They stay under gitignored
`data/raw/`. The raw GEO files and Horvath 353-CpG coefficient list remain the
source of truth.

`data/processed/` is only a working copy for inspection/compatibility.

## Cache equivalence proof

Command run:

```bash
python tests/test_cache_equivalence.py
```

Result:

- passed
- elapsed: 259.891 seconds
- metadata comparison: exact, including dtype
- beta matrix comparison: exact with `atol=0.0`, `rtol=0.0`

The test rebuilt data from raw GEO files in memory and compared it with the
cache-loaded data:

- samples shape: 656 rows
- beta matrix shape: 353 x 656

## Cold and warm make all

Cold run command:

```bash
make clean
make all
```

Cold run result:

- passed
- elapsed: 296.111 seconds
- cache rebuilt from raw GEO files

Warm run command:

```bash
make all
```

Warm run result:

- passed
- elapsed: 80.578 seconds
- `make data` warm path prints: `cache present, skipping parse`
- warm `make data` alone elapsed: 4.738 seconds

## Final numbers unchanged

`results/metrics.csv` after cold and warm runs:

- MAE: 3.856329584054624
- Pearson r: 0.9183434720662662
- Pearson p-value: 1.8569795187660721e-265
- RMSE: 6.322123159246127
- slope: 0.7965811960029956
- intercept: 10.693559505245958
- n_samples: 656

`results/crosscheck_biolearn_lightweight.csv`:

- compared samples: 656
- mean absolute difference: 1.0831444142684454e-17 years
- max absolute difference: 3.552713678800501e-15 years
- Pearson r between prediction vectors: 1.0

## Caveats

- Warm `make all` is faster than cold, but can still take meaningful time
  because it starts multiple Python processes and imports pandas/scipy. The
  cache removes the expensive raw parse/extraction step; it does not eliminate
  Python import/runtime overhead.
- The metadata parser skips embedded SOFT sample-table lines once
  `!sample_table_begin` is reached. Cache equivalence passed after this change,
  so the optimization did not alter the derived data.
- Full Biolearn model API cross-check is still not used because importing
  `biolearn.model` requires `torch`. The current cross-check is the lightweight
  Biolearn coefficient-data path.
