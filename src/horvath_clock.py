"""Hand implementation of the Horvath 2013 pan-tissue clock.

The clock predicts transformed age. The public output must be inverse
transformed back to age in years.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd


LOGGER = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_COEFFICIENTS_PATH = PROJECT_ROOT / "data" / "coefficients" / "horvath2013_coefficients.csv"
ADULT_AGE = 20.0


def inverse_horvath_age_transform(linear_predictor: np.ndarray | pd.Series) -> np.ndarray:
    """Invert Horvath's transformed-age scale.

    Horvath 2013 uses adult.age = 20 and:

        F(age) = log(age + 1) - log(adult.age + 1), if age <= adult.age
        F(age) = (age - adult.age) / (adult.age + 1), if age > adult.age

    The model's linear predictor is on the F(age) scale. The inverse is:

        age = exp(L + log(adult.age + 1)) - 1, if L <= 0
        age = L * (adult.age + 1) + adult.age, if L > 0

    Args:
        linear_predictor: model output on Horvath's transformed-age scale.

    Returns:
        Predicted DNAmAge in years.
    """
    values = np.asarray(linear_predictor, dtype=float)
    young = np.exp(values + np.log(ADULT_AGE + 1.0)) - 1.0
    adult = values * (ADULT_AGE + 1.0) + ADULT_AGE
    return np.where(values <= 0.0, young, adult)


def load_coefficients(path: Path = DEFAULT_COEFFICIENTS_PATH) -> tuple[float, pd.DataFrame]:
    """Load Horvath coefficients.

    Expected CSV columns:
        - cpg
        - coefficient
        - optional mean

    Intercept can be provided as a row with cpg == "intercept" or "(Intercept)".
    """
    if not path.exists():
        raise FileNotFoundError(
            f"Missing coefficient file: {path}. "
            "Add the validated 353-CpG Horvath 2013 coefficient table before analysis."
        )

    coeffs = pd.read_csv(path)
    required = {"cpg", "coefficient"}
    missing = required - set(coeffs.columns)
    if missing:
        raise ValueError(f"Coefficient file missing required columns: {sorted(missing)}")

    cpg_norm = coeffs["cpg"].astype(str).str.strip().str.lower()
    intercept_mask = cpg_norm.isin({"intercept", "(intercept)"})
    if intercept_mask.sum() != 1:
        raise ValueError("Coefficient file must contain exactly one intercept row")

    intercept = float(coeffs.loc[intercept_mask, "coefficient"].iloc[0])
    cpg_coeffs = coeffs.loc[~intercept_mask].copy()
    if len(cpg_coeffs) != 353:
        raise ValueError(f"Expected 353 CpG coefficients, found {len(cpg_coeffs)}")
    if cpg_coeffs["cpg"].duplicated().any():
        duplicates = cpg_coeffs.loc[cpg_coeffs["cpg"].duplicated(), "cpg"].tolist()
        raise ValueError(f"Duplicate CpG IDs in coefficient file: {duplicates[:10]}")

    return intercept, cpg_coeffs


def predict_age(
    betas_matrix: pd.DataFrame,
    coefficients_path: Path = DEFAULT_COEFFICIENTS_PATH,
) -> pd.Series:
    """Predict DNAmAge for samples from a CpG x sample beta matrix.

    Missing CpGs are imputed with a training-set mean if the coefficient table
    provides a `mean` column. If unavailable, this scaffold uses the global mean
    of observed beta values and logs the count. That fallback must be reported in
    the final replication note if used.
    """
    intercept, coeffs = load_coefficients(coefficients_path)

    if betas_matrix.index.has_duplicates:
        raise ValueError("Beta matrix CpG index contains duplicates")

    cpg_ids = coeffs["cpg"].astype(str).tolist()
    present = [cpg for cpg in cpg_ids if cpg in betas_matrix.index]
    missing = [cpg for cpg in cpg_ids if cpg not in betas_matrix.index]

    if missing:
        LOGGER.warning("Missing Horvath CpGs: %d", len(missing))

    aligned = pd.DataFrame(index=cpg_ids, columns=betas_matrix.columns, dtype=float)
    if present:
        aligned.loc[present] = betas_matrix.loc[present].astype(float)

    if missing:
        if "mean" in coeffs.columns and coeffs["mean"].notna().all():
            mean_map = coeffs.set_index("cpg")["mean"].astype(float)
            for cpg in missing:
                aligned.loc[cpg] = mean_map.loc[cpg]
        else:
            fallback_mean = float(np.nanmean(betas_matrix.to_numpy(dtype=float)))
            LOGGER.warning(
                "Training means unavailable; imputing %d missing CpGs with global dataset mean %.6f",
                len(missing),
                fallback_mean,
            )
            for cpg in missing:
                aligned.loc[cpg] = fallback_mean

    if aligned.isna().any().any():
        raise ValueError("Aligned beta matrix contains NaN values after imputation")

    weights = coeffs.set_index("cpg").loc[cpg_ids, "coefficient"].astype(float)
    # Avoid pandas.DataFrame.dot here: explicit multiply+sum keeps the tiny
    # 353-CpG calculation simple, deterministic, and independent of platform
    # BLAS/DataFrame-dot behavior.
    beta_values = aligned.to_numpy(dtype=np.float64, copy=False)
    weight_values = weights.to_numpy(dtype=np.float64, copy=False)
    linear_predictor = pd.Series(
        intercept + (beta_values.T * weight_values).sum(axis=1),
        index=aligned.columns,
        name="linear_predictor",
    )
    predicted = inverse_horvath_age_transform(linear_predictor)
    return pd.Series(predicted, index=betas_matrix.columns, name="predicted_age")
