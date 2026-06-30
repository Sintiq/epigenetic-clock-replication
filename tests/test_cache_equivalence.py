"""Prove derived cache matches a fresh raw-file build.

This is intentionally slower than a normal unit test. It reparses the local GEO
SOFT metadata and re-extracts the Horvath 353-CpG beta subset from the split raw
beta files, then compares that fresh result with the cached parquet files.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from download_data import build_fresh_data, cache_present, load_cached_data  # noqa: E402


BETA_ATOL = 0.0
BETA_RTOL = 0.0


def assert_cache_equivalence() -> None:
    if not cache_present():
        raise FileNotFoundError("Cache is missing. Run `make data` before cache equivalence test.")

    print("Building fresh data from raw GEO files for cache equivalence proof...")
    fresh_samples, fresh_betas = build_fresh_data()
    print("Loading derived cache...")
    cached_samples, cached_betas = load_cached_data()

    assert_frame_equal(
        fresh_samples.reset_index(drop=True),
        cached_samples.reset_index(drop=True),
        check_dtype=True,
        check_exact=True,
    )
    assert list(fresh_betas.index) == list(cached_betas.index)
    assert list(fresh_betas.columns) == list(cached_betas.columns)
    assert fresh_betas.shape == cached_betas.shape
    assert str(fresh_betas.dtypes.unique().tolist()) == str(cached_betas.dtypes.unique().tolist())
    if not np.allclose(
        fresh_betas.to_numpy(dtype=np.float64),
        cached_betas.to_numpy(dtype=np.float64),
        rtol=BETA_RTOL,
        atol=BETA_ATOL,
        equal_nan=False,
    ):
        diff = np.abs(
            fresh_betas.to_numpy(dtype=np.float64) - cached_betas.to_numpy(dtype=np.float64)
        )
        raise AssertionError(f"Beta cache differs from fresh build; max_abs_diff={diff.max()}")

    print("Cache equivalence passed.")
    print(f"samples_shape={fresh_samples.shape}")
    print(f"betas_shape={fresh_betas.shape}")
    print(f"beta_tolerance_atol={BETA_ATOL}")
    print(f"beta_tolerance_rtol={BETA_RTOL}")


def test_cache_equivalence() -> None:
    assert_cache_equivalence()


if __name__ == "__main__":
    assert_cache_equivalence()
