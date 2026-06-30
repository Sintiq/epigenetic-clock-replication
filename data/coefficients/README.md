# Horvath 2013 coefficients

File:

- `horvath2013_coefficients.csv`

Status:

- Added and structurally validated.

Contents:

- 353 CpG coefficients
- one explicit intercept row
- columns: `cpg`, `coefficient`, `mean`

Source used for this working copy:

- Biolearn 0.9.1 package, file `biolearn/data/Horvath1.csv`
- Biolearn package license reported by `pip show biolearn`: `new BSD`
- Biolearn model definition: `Horvathv1`, year 2013, source
  `https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115`
- Original Horvath 2013 article and supplementary coefficient files are
  attributed to Horvath 2013. The article is distributed under the Creative
  Commons Attribution License 2.0, requiring proper citation of the original
  work.

Important transform convention:

- Biolearn stores 353 CpG coefficients in `Horvath1.csv` without an intercept row.
- Its `Horvathv1` model applies `anti_trafo(sum + 0.696)`.
- This repository writes that `+0.696` as an explicit `intercept` row so the
  hand implementation can compute:

```text
linear_predictor = intercept + sum(beta_i * coefficient_i)
DNAmAge = inverse_horvath_age_transform(linear_predictor)
```

The `mean` column is currently empty because the Biolearn coefficient file does
not provide per-CpG training means. Missing CpGs must therefore be logged and
handled transparently by the analysis.
