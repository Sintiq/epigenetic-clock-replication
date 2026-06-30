"""Run the Horvath clock replication analysis.

This script intentionally refuses to fabricate a result. It requires processed
sample metadata, processed beta values, and validated coefficients.
"""

from __future__ import annotations

from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

try:
    from .download_data import load_cached_data
    from .horvath_clock import predict_age
except ImportError:  # Allow running this file directly as `python src/analyze.py`.
    from download_data import load_cached_data
    from horvath_clock import predict_age


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / "results"


def _fit_slope_intercept(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    """Fit y = slope*x + intercept without numpy.linalg."""
    x_values = x.to_numpy(dtype=float)
    y_values = y.to_numpy(dtype=float)
    x_mean = float(x_values.mean())
    y_mean = float(y_values.mean())
    denominator = float(((x_values - x_mean) ** 2).sum())
    if denominator == 0.0:
        raise ValueError("Cannot fit slope: chronological age has zero variance")
    slope = float(((x_values - x_mean) * (y_values - y_mean)).sum() / denominator)
    intercept = float(y_mean - slope * x_mean)
    return slope, intercept


def _write_svg_scatter(merged: pd.DataFrame, slope: float, intercept: float, out_path: Path) -> None:
    """Write a simple SVG scatter plot without Matplotlib native rendering."""
    width = 720
    height = 720
    margin = 72
    plot_width = width - 2 * margin
    plot_height = height - 2 * margin

    x_values = merged["age"].to_numpy(dtype=float)
    y_values = merged["predicted_age"].to_numpy(dtype=float)
    low = float(min(x_values.min(), y_values.min()))
    high = float(max(x_values.max(), y_values.max()))
    pad = max(2.0, (high - low) * 0.05)
    low -= pad
    high += pad

    def sx(value: float) -> float:
        return margin + (value - low) / (high - low) * plot_width

    def sy(value: float) -> float:
        return height - margin - (value - low) / (high - low) * plot_height

    yx_line = (sx(low), sy(low), sx(high), sy(high))
    fit_line = (sx(low), sy(slope * low + intercept), sx(high), sy(slope * high + intercept))
    points = "\n".join(
        f'<circle cx="{sx(float(x)):.2f}" cy="{sy(float(y)):.2f}" r="2.2" fill="#2f6f9f" fill-opacity="0.62" />'
        for x, y in zip(x_values, y_values, strict=True)
    )

    tick_start = np.ceil(low / 20.0) * 20.0
    tick_end = np.floor(high / 20.0) * 20.0
    tick_values = np.linspace(tick_start, tick_end, 5)
    ticks = []
    for value in tick_values:
        value = float(value)
        x = sx(value)
        y = sy(value)
        ticks.append(f'<line x1="{x:.2f}" y1="{height-margin:.2f}" x2="{x:.2f}" y2="{height-margin+6:.2f}" stroke="#333" />')
        ticks.append(f'<text x="{x:.2f}" y="{height-margin+24:.2f}" text-anchor="middle" font-size="12">{value:.0f}</text>')
        ticks.append(f'<line x1="{margin-6:.2f}" y1="{y:.2f}" x2="{margin:.2f}" y2="{y:.2f}" stroke="#333" />')
        ticks.append(f'<text x="{margin-12:.2f}" y="{y+4:.2f}" text-anchor="end" font-size="12">{value:.0f}</text>')

    title = escape("Horvath 2013 clock on GSE40279")
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect width="100%" height="100%" fill="white" />
  <text x="{width/2:.0f}" y="34" text-anchor="middle" font-family="Arial, sans-serif" font-size="20" font-weight="700">{title}</text>
  <line x1="{margin}" y1="{height-margin}" x2="{width-margin}" y2="{height-margin}" stroke="#222" />
  <line x1="{margin}" y1="{margin}" x2="{margin}" y2="{height-margin}" stroke="#222" />
  {''.join(ticks)}
  <line x1="{yx_line[0]:.2f}" y1="{yx_line[1]:.2f}" x2="{yx_line[2]:.2f}" y2="{yx_line[3]:.2f}" stroke="#111" stroke-width="1.2" stroke-dasharray="6 5" />
  <line x1="{fit_line[0]:.2f}" y1="{fit_line[1]:.2f}" x2="{fit_line[2]:.2f}" y2="{fit_line[3]:.2f}" stroke="#c43d3d" stroke-width="2.0" />
  {points}
  <text x="{width/2:.0f}" y="{height-18}" text-anchor="middle" font-family="Arial, sans-serif" font-size="14">Chronological age</text>
  <text x="20" y="{height/2:.0f}" text-anchor="middle" font-family="Arial, sans-serif" font-size="14" transform="rotate(-90 20 {height/2:.0f})">Predicted DNAmAge</text>
  <text x="{width-margin-120}" y="{margin+18}" font-family="Arial, sans-serif" font-size="13" fill="#111">dashed: y=x</text>
  <text x="{width-margin-120}" y="{margin+38}" font-family="Arial, sans-serif" font-size="13" fill="#c43d3d">red: fitted line</text>
</svg>
'''
    out_path.write_text(svg, encoding="utf-8")


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    samples, betas = load_cached_data()
    if not {"sample_id", "age"}.issubset(samples.columns):
        raise ValueError("Sample table must contain sample_id and age")

    predicted = predict_age(betas)

    merged = samples.set_index("sample_id").join(predicted, how="inner")
    if merged.empty:
        raise ValueError("No overlap between sample metadata and predicted ages")

    residual = merged["predicted_age"] - merged["age"]
    mae = float(np.median(np.abs(residual)))
    rmse = float(np.sqrt(np.mean(np.square(residual.to_numpy(dtype=float)))))
    r, p_value = pearsonr(merged["age"], merged["predicted_age"])
    slope, intercept = _fit_slope_intercept(merged["age"], merged["predicted_age"])

    metrics = pd.DataFrame(
        [
            {"metric": "median_absolute_error_years", "value": mae, "published_reference": "~3.6"},
            {"metric": "pearson_r", "value": float(r), "published_reference": ">0.9"},
            {"metric": "pearson_p_value", "value": float(p_value), "published_reference": ""},
            {"metric": "rmse_years", "value": rmse, "published_reference": ""},
            {"metric": "slope", "value": float(slope), "published_reference": ""},
            {"metric": "intercept", "value": float(intercept), "published_reference": ""},
            {"metric": "n_samples", "value": int(len(merged)), "published_reference": "~656"},
        ]
    )
    metrics.to_csv(RESULTS_DIR / "metrics.csv", index=False)

    _write_svg_scatter(merged, slope, intercept, RESULTS_DIR / "predicted_vs_chronological_age.svg")

    replicated = mae <= 4.6 and r > 0.9
    note = [
        "# Replication note",
        "",
        f"Median absolute error: {mae:.3f} years",
        f"Pearson r: {r:.4f}",
        f"RMSE: {rmse:.3f} years",
        f"Slope: {slope:.4f}",
        "",
        "Pre-declared success criterion: MAE within approximately 1 year of the",
        "published blood reference (~3.6 years) and Pearson r above 0.9.",
        "",
        f"Conclusion: {'replicated under this criterion' if replicated else 'did not replicate under this criterion'}.",
        "",
        "If not replicated, inspect normalization differences, missing CpGs, the",
        "age transform, tissue/platform assumptions, and batch effects.",
    ]
    (RESULTS_DIR / "REPLICATION_NOTE.md").write_text("\n".join(note) + "\n", encoding="utf-8")
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
