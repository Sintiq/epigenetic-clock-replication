"""Data preparation for GSE40279.

This module is intentionally conservative: it should fail loudly if sample age
metadata cannot be parsed from GEO metadata. Do not guess ages from row order or
file names.
"""

from __future__ import annotations

import argparse
import gzip
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

try:
    from .horvath_clock import DEFAULT_COEFFICIENTS_PATH, load_coefficients
except ImportError:  # Allow running this file directly as `python src/download_data.py`.
    from horvath_clock import DEFAULT_COEFFICIENTS_PATH, load_coefficients


ACCESSION = "GSE40279"

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data" / "raw"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
SOFT_FILE = RAW_DIR / f"{ACCESSION}_family.soft.gz"
HORVATH_BETAS_PATH = PROCESSED_DIR / f"{ACCESSION}_horvath353_betas.parquet"
SAMPLES_PATH = PROCESSED_DIR / f"{ACCESSION}_samples.csv"
BETAS_CACHE_PATH = RAW_DIR / "betas_subset.parquet"
SAMPLES_CACHE_PATH = RAW_DIR / "samples_metadata.parquet"


CHARACTERISTIC_KEY_VALUE_PATTERN = re.compile(
    r"^\s*(?P<key>[^:=]+?)\s*[:=]\s*(?P<value>.*?)\s*$"
)
AGE_VALUE_PATTERN = re.compile(r"^\s*(?P<age>[-+]?\d+(?:\.\d+)?)\s*(?:y|yr|yrs|years?)?\s*$", re.I)
SAMPLE_START_PATTERN = re.compile(r"^\^SAMPLE\s*=\s*(?P<sample_id>GSM\d+)")
CHARACTERISTIC_PATTERN = re.compile(r"^!Sample_characteristics_ch1\s*=\s*(?P<value>.*)$")
TITLE_PATTERN = re.compile(r"^!Sample_title\s*=\s*(?P<value>.*)$")
SOURCE_NAME_PATTERN = re.compile(r"^!Sample_source_name_ch1\s*=\s*(?P<value>.*)$")


@dataclass(frozen=True)
class SampleMetadata:
    sample_id: str
    age: float
    title: str | None = None
    source_name: str | None = None


def _parse_age_characteristic(value: str) -> float | None:
    """Return a numeric age if the characteristic is clearly an age field."""
    key_value = CHARACTERISTIC_KEY_VALUE_PATTERN.match(value)
    if not key_value:
        return None
    key = re.sub(r"[^a-z0-9]+", "", key_value.group("key").lower())
    if key not in {"age", "agey", "ageyears"}:
        return None
    age_match = AGE_VALUE_PATTERN.match(key_value.group("value"))
    if not age_match:
        raise ValueError(f"Found age characteristic but could not parse numeric age: {value!r}")
    return float(age_match.group("age"))


def parse_sample_metadata_from_soft(soft_path: Path = SOFT_FILE) -> pd.DataFrame:
    """Parse sample_id and age from a GEO family SOFT file.

    Raises:
        FileNotFoundError: if the SOFT file is missing.
        ValueError: if no samples or any missing/ambiguous ages are found.
    """
    if not soft_path.exists():
        raise FileNotFoundError(f"Missing GEO SOFT file: {soft_path}")

    samples: list[SampleMetadata] = []
    current_sample: str | None = None
    current_title: str | None = None
    current_source_name: str | None = None
    current_ages: list[float] = []

    def flush_current() -> None:
        nonlocal current_sample, current_title, current_source_name, current_ages
        if current_sample is None:
            return
        if len(current_ages) != 1:
            raise ValueError(
                f"Expected exactly one age for {current_sample}, found {len(current_ages)}. "
                "Refusing to guess from GEO metadata."
            )
        if not current_source_name:
            raise ValueError(f"Missing source_name for {current_sample}; cannot map beta columns safely.")
        samples.append(SampleMetadata(current_sample, current_ages[0], current_title, current_source_name))
        current_sample = None
        current_title = None
        current_source_name = None
        current_ages = []

    with gzip.open(soft_path, "rt", encoding="utf-8", errors="replace") as handle:
        in_sample_table = False
        for line in handle:
            line = line.rstrip("\n")
            if in_sample_table:
                if line.startswith("!sample_table_end"):
                    in_sample_table = False
                continue

            if line.startswith("!sample_table_begin"):
                in_sample_table = True
                continue

            sample_match = SAMPLE_START_PATTERN.match(line)
            if sample_match:
                flush_current()
                current_sample = sample_match.group("sample_id")
                continue

            if current_sample is None:
                continue

            title_match = TITLE_PATTERN.match(line)
            if title_match:
                current_title = title_match.group("value").strip()
                continue

            source_name_match = SOURCE_NAME_PATTERN.match(line)
            if source_name_match:
                current_source_name = source_name_match.group("value").strip()
                continue

            characteristic_match = CHARACTERISTIC_PATTERN.match(line)
            if characteristic_match:
                maybe_age = _parse_age_characteristic(characteristic_match.group("value"))
                if maybe_age is not None:
                    current_ages.append(maybe_age)

    flush_current()

    if not samples:
        raise ValueError(f"No samples with age metadata parsed from {soft_path}")

    frame = pd.DataFrame([sample.__dict__ for sample in samples])
    frame = frame.sort_values("sample_id").reset_index(drop=True)
    if frame["age"].isna().any():
        raise ValueError("Parsed sample table contains missing ages")
    if frame["source_name"].isna().any() or frame["source_name"].duplicated().any():
        raise ValueError("Parsed sample table contains missing or duplicate source_name values")
    return frame


def extract_horvath_beta_matrix(
    samples: pd.DataFrame,
    coefficients_path: Path = DEFAULT_COEFFICIENTS_PATH,
    raw_dir: Path = RAW_DIR,
    out_path: Path | None = HORVATH_BETAS_PATH,
    chunksize: int = 50_000,
) -> pd.DataFrame:
    """Extract only Horvath clock CpGs from the split GSE40279 beta files.

    The raw GEO beta files are large and split by sample columns. For the
    Horvath replication, keeping only the 353 clock CpGs is enough for the
    analysis and avoids materializing the full 473k x 656 matrix.
    """
    _, coeffs = load_coefficients(coefficients_path)
    wanted_cpgs = coeffs["cpg"].astype(str).tolist()
    wanted_set = set(wanted_cpgs)

    source_to_sample = dict(zip(samples["source_name"].astype(str), samples["sample_id"].astype(str)))
    raw_files = sorted(raw_dir.glob(f"{ACCESSION}_average_beta_GSM*.txt.gz"))
    if not raw_files:
        raise FileNotFoundError(f"No split beta files found under {raw_dir}")

    split_frames: list[pd.DataFrame] = []
    observed_sample_columns: set[str] = set()

    for raw_file in raw_files:
        print(f"Reading {raw_file.name}")
        matched_chunks: list[pd.DataFrame] = []
        for chunk in pd.read_csv(raw_file, sep="\t", compression="gzip", chunksize=chunksize):
            if "ID_REF" not in chunk.columns:
                raise ValueError(f"{raw_file} missing ID_REF column")
            matched = chunk[chunk["ID_REF"].isin(wanted_set)].copy()
            if matched.empty:
                continue
            matched_chunks.append(matched)

        if not matched_chunks:
            raise ValueError(f"No Horvath CpGs found in {raw_file}")

        split = pd.concat(matched_chunks, ignore_index=True)
        if split["ID_REF"].duplicated().any():
            dupes = split.loc[split["ID_REF"].duplicated(), "ID_REF"].tolist()
            raise ValueError(f"Duplicate CpGs in {raw_file}: {dupes[:10]}")

        rename_map = {"ID_REF": "cpg"}
        sample_columns = [column for column in split.columns if column != "ID_REF"]
        unknown_columns = [column for column in sample_columns if column not in source_to_sample]
        if unknown_columns:
            raise ValueError(f"{raw_file} has beta columns not found in SOFT source_name metadata: {unknown_columns[:10]}")
        for column in sample_columns:
            rename_map[column] = source_to_sample[column]
        observed_sample_columns.update(rename_map[column] for column in sample_columns)

        split = split.rename(columns=rename_map).set_index("cpg")
        split = split.apply(pd.to_numeric, errors="raise").astype("float32")
        split_frames.append(split)
        print(f"  matched_cpgs={len(split)} sample_columns={len(sample_columns)}")

    betas = pd.concat(split_frames, axis=1)
    if betas.columns.duplicated().any():
        dupes = betas.columns[betas.columns.duplicated()].tolist()
        raise ValueError(f"Duplicate sample columns after merging split beta files: {dupes[:10]}")

    missing_samples = sorted(set(samples["sample_id"].astype(str)) - set(betas.columns.astype(str)))
    extra_samples = sorted(set(betas.columns.astype(str)) - set(samples["sample_id"].astype(str)))
    if missing_samples or extra_samples:
        raise ValueError(
            "Beta/sample metadata mismatch: "
            f"missing_samples={missing_samples[:10]} extra_samples={extra_samples[:10]}"
        )

    betas = betas.reindex(index=wanted_cpgs, columns=samples["sample_id"].astype(str).tolist())
    missing_cpgs = [cpg for cpg in wanted_cpgs if cpg not in betas.index or betas.loc[cpg].isna().all()]
    if missing_cpgs:
        print(f"Missing Horvath CpGs in beta files: {len(missing_cpgs)}")

    if betas.isna().any().any():
        na_count = int(betas.isna().sum().sum())
        raise ValueError(f"Horvath beta subset contains {na_count} missing beta values")

    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        betas.to_parquet(out_path)
        print(f"Wrote Horvath beta subset: {out_path}")
    print(f"rows_cpg={betas.shape[0]} columns_samples={betas.shape[1]}")
    return betas


def _nonempty_file(path: Path) -> bool:
    return path.exists() and path.is_file() and path.stat().st_size > 0


def cache_present() -> bool:
    """Return true when both derived cache files exist and are non-empty."""
    return _nonempty_file(BETAS_CACHE_PATH) and _nonempty_file(SAMPLES_CACHE_PATH)


def validate_samples(samples: pd.DataFrame) -> None:
    required_columns = {"sample_id", "age", "source_name"}
    missing = required_columns - set(samples.columns)
    if missing:
        raise ValueError(f"Sample metadata missing required columns: {sorted(missing)}")
    if samples.empty:
        raise ValueError("Sample metadata is empty")
    if samples["sample_id"].isna().any() or samples["sample_id"].duplicated().any():
        raise ValueError("Sample metadata contains missing or duplicate sample_id values")
    if samples["source_name"].isna().any() or samples["source_name"].duplicated().any():
        raise ValueError("Sample metadata contains missing or duplicate source_name values")
    if samples["age"].isna().any():
        raise ValueError("Sample metadata contains missing ages")


def validate_betas(betas: pd.DataFrame, samples: pd.DataFrame) -> None:
    _, coeffs = load_coefficients(DEFAULT_COEFFICIENTS_PATH)
    wanted_cpgs = coeffs["cpg"].astype(str).tolist()
    wanted_samples = samples["sample_id"].astype(str).tolist()
    if betas.empty:
        raise ValueError("Beta matrix is empty")
    if list(betas.index.astype(str)) != wanted_cpgs:
        raise ValueError("Beta matrix CpG index does not match Horvath coefficient order")
    if list(betas.columns.astype(str)) != wanted_samples:
        raise ValueError("Beta matrix sample columns do not match sample metadata order")
    if betas.isna().any().any():
        na_count = int(betas.isna().sum().sum())
        raise ValueError(f"Beta matrix contains {na_count} missing values")


def load_cached_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load derived cache files and validate their basic contract."""
    if not cache_present():
        raise FileNotFoundError(
            "Derived cache is missing. Run `make data` first. "
            f"Expected {SAMPLES_CACHE_PATH} and {BETAS_CACHE_PATH}."
        )
    samples = pd.read_parquet(SAMPLES_CACHE_PATH)
    betas = pd.read_parquet(BETAS_CACHE_PATH)
    validate_samples(samples)
    validate_betas(betas, samples)
    return samples, betas


def write_processed_copies(samples: pd.DataFrame, betas: pd.DataFrame) -> None:
    """Write data/processed working copies for inspection and compatibility."""
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    samples.to_csv(SAMPLES_PATH, index=False)
    betas.to_parquet(HORVATH_BETAS_PATH)


def build_fresh_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build metadata and Horvath beta subset from raw GEO files."""
    samples = parse_sample_metadata_from_soft()
    validate_samples(samples)
    betas = extract_horvath_beta_matrix(samples, out_path=None)
    validate_betas(betas, samples)
    return samples, betas


def write_cache(samples: pd.DataFrame, betas: pd.DataFrame) -> None:
    """Persist derived cache files in data/raw."""
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    samples.to_parquet(SAMPLES_CACHE_PATH, index=False)
    betas.to_parquet(BETAS_CACHE_PATH)
    print(f"Wrote samples cache: {SAMPLES_CACHE_PATH}")
    print(f"Wrote beta cache: {BETAS_CACHE_PATH}")


def ensure_cached_data(force_rebuild: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load cache if present; otherwise build it from raw files.

    The cache is derived. The raw GEO files and Horvath coefficient list remain
    the source of truth.
    """
    if cache_present() and not force_rebuild:
        print("cache present, skipping parse")
        samples, betas = load_cached_data()
    else:
        if force_rebuild:
            print("forcing fresh parse and cache rebuild")
        else:
            print("cache missing, parsing raw GEO files")
        samples, betas = build_fresh_data()
        write_cache(samples, betas)
        samples, betas = load_cached_data()

    write_processed_copies(samples, betas)
    return samples, betas


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare GSE40279 metadata and Horvath beta subset.")
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="Parse sample metadata only. Beta matrix conversion is a later step.",
    )
    parser.add_argument(
        "--horvath-betas-only",
        action="store_true",
        help="Build only the 353-CpG Horvath beta subset from already parsed metadata.",
    )
    parser.add_argument(
        "--force-rebuild-cache",
        action="store_true",
        help="Ignore derived cache and rebuild metadata/beta subset from raw GEO files.",
    )
    args = parser.parse_args()

    if args.force_rebuild_cache:
        ensure_cached_data(force_rebuild=True)
        return

    if args.metadata_only:
        samples = parse_sample_metadata_from_soft()
        validate_samples(samples)
        PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
        samples.to_csv(SAMPLES_PATH, index=False)
        print(f"Wrote {len(samples)} sample metadata rows to {SAMPLES_PATH}")
        return

    if args.horvath_betas_only:
        if SAMPLES_PATH.exists():
            samples = pd.read_csv(SAMPLES_PATH)
            validate_samples(samples)
        else:
            samples = parse_sample_metadata_from_soft()
            validate_samples(samples)
            PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
            samples.to_csv(SAMPLES_PATH, index=False)
        extract_horvath_beta_matrix(samples)
        return

    ensure_cached_data(force_rebuild=False)


if __name__ == "__main__":
    main()
