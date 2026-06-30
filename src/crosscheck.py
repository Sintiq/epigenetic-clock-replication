"""Lightweight independent cross-check against Biolearn's Horvathv1 data.

Biolearn's public Horvathv1 model definition uses:

    Horvath1.csv
    anti_trafo(sum + 0.696)

Importing `biolearn.model` currently requires torch in this environment. This
cross-check therefore reads the Biolearn coefficient CSV directly and applies
the documented transform in a separate implementation, avoiding our committed
coefficient file and our `predict_age` internals.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

try:
    from .download_data import load_cached_data
    from .horvath_clock import predict_age
except ImportError:  # Allow running this file directly as `python src/crosscheck.py`.
    from download_data import load_cached_data
    from horvath_clock import predict_age


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / "results"
CROSSCHECK_PATH = RESULTS_DIR / "crosscheck_biolearn_lightweight.csv"


def _biolearn_horvath1_path() -> Path:
    spec = importlib.util.find_spec("biolearn")
    if spec is None or spec.origin is None:
        raise ImportError("biolearn package is not installed")
    package_dir = Path(spec.origin).resolve().parent
    path = package_dir / "data" / "Horvath1.csv"
    if not path.exists():
        raise FileNotFoundError(f"Biolearn Horvath1 coefficient file not found: {path}")
    return path


def _biolearn_anti_trafo(values: np.ndarray, adult_age: float = 20.0) -> np.ndarray:
    """Biolearn/Horvath inverse transform convention."""
    return np.where(
        values < 0,
        (1.0 + adult_age) * np.exp(values) - 1.0,
        (1.0 + adult_age) * values + adult_age,
    )


def _predict_with_biolearn_csv(betas: pd.DataFrame, biolearn_csv: Path) -> pd.Series:
    coeffs = pd.read_csv(biolearn_csv)
    expected_columns = {"CpGmarker", "CoefficientTraining"}
    if set(coeffs.columns) != expected_columns:
        raise ValueError(f"Unexpected Biolearn Horvath1 columns: {list(coeffs.columns)}")
    if len(coeffs) != 353:
        raise ValueError(f"Expected 353 Biolearn Horvath1 CpGs, found {len(coeffs)}")
    if coeffs["CpGmarker"].duplicated().any():
        raise ValueError("Biolearn Horvath1 coefficient file contains duplicate CpGs")

    cpgs = coeffs["CpGmarker"].astype(str).tolist()
    missing = [cpg for cpg in cpgs if cpg not in betas.index]
    if missing:
        raise ValueError(f"Processed beta matrix missing Biolearn Horvath1 CpGs: {missing[:10]}")

    aligned = betas.loc[cpgs]
    beta_values = aligned.to_numpy(dtype=np.float64, copy=False)
    weight_values = coeffs["CoefficientTraining"].to_numpy(dtype=np.float64, copy=False)
    transformed = (beta_values.T * weight_values).sum(axis=1) + 0.696
    predicted = _biolearn_anti_trafo(transformed)
    return pd.Series(predicted, index=aligned.columns, name="biolearn_lightweight_predicted_age")


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    _, betas = load_cached_data()
    our_predictions = predict_age(betas)
    biolearn_path = _biolearn_horvath1_path()
    biolearn_predictions = _predict_with_biolearn_csv(betas, biolearn_path)

    aligned = pd.concat([our_predictions, biolearn_predictions], axis=1, join="inner")
    if len(aligned) != betas.shape[1]:
        raise ValueError(
            f"Prediction sample mismatch: compared {len(aligned)} of {betas.shape[1]} samples"
        )

    diff = aligned.iloc[:, 0] - aligned.iloc[:, 1]
    mean_abs_diff = float(np.mean(np.abs(diff)))
    max_abs_diff = float(np.max(np.abs(diff)))
    correlation = float(pearsonr(aligned.iloc[:, 0], aligned.iloc[:, 1]).statistic)

    result = pd.DataFrame(
        [
            {"metric": "n_samples_compared", "value": int(len(aligned))},
            {"metric": "mean_abs_difference_years", "value": mean_abs_diff},
            {"metric": "max_abs_difference_years", "value": max_abs_diff},
            {"metric": "pearson_r_between_predictions", "value": correlation},
            {"metric": "biolearn_coefficients_rows", "value": 353},
            {"metric": "biolearn_intercept_from_model_definition", "value": 0.696},
        ]
    )
    result.to_csv(CROSSCHECK_PATH, index=False)
    print(result.to_string(index=False))
    print(f"Biolearn coefficient source: {biolearn_path}")
    print(f"Wrote {CROSSCHECK_PATH}")


if __name__ == "__main__":
    main()
