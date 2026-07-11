#!/usr/bin/env python3
"""
distance_profile.py

Compute distance-dependent homopolymer enrichment profiles:

    E_k(d) = P(run length >= k | distance from nearest boundary lies in bin d)

This script is intended as the next step after boundary_enrichment_v3.py.
It uses the `run_table.csv` file produced by that script.

Input
-----
run_table.csv must contain at least:

- accession
- base
- run_length
- nearest_boundary_distance

Optional but useful columns:

- organism
- sequence_length
- in_boundary_window

Outputs
-------
1. distance_profile.csv
   Binned estimates of E_k(d), all bases combined and base-specific.

2. distance_profile_by_genome.csv
   Per-genome binned estimates, useful for checking whether the aggregate
   profile is dominated by a few large genomes.

3. distance_profile_summary.csv
   Compact summary of near/far enrichment ratios.

4. PDF plots:
   - distance_profile_all_bases.pdf
   - distance_profile_by_base_k*.pdf

Example
-------
python distance_profile.py \
    --run-table analysis_initial/boundary_enrichment_v3/run_table.csv \
    --outdir analysis_initial/boundary_distance_profile \
    --thresholds 5,10,15 \
    --bins 0,25,50,100,200,400,800,1600,inf

Notes
-----
- Distances are in nucleotides/bases.
- Runs with no nearest boundary are dropped.
- Use --normalise-by-genome to compute an unweighted mean across genomes,
  rather than pooling all runs. This is often a good sensitivity check.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import pandas as pd


BASES = ["A", "C", "G", "T"]


def parse_thresholds(text: str) -> list[int]:
    return [int(x.strip()) for x in text.split(",") if x.strip()]


def parse_bins(text: str) -> list[float]:
    out = []
    for x in text.split(","):
        x = x.strip().lower()
        if not x:
            continue
        if x in {"inf", "infinity", "∞"}:
            out.append(float("inf"))
        else:
            out.append(float(x))
    if len(out) < 2:
        raise ValueError("At least two bin edges are required")
    if out != sorted(out):
        raise ValueError("Bin edges must be sorted increasingly")
    return out


def bin_label(left: float, right: float) -> str:
    if math.isinf(right):
        return f"{int(left)}+"
    return f"{int(left)}-{int(right)}"


def assign_distance_bins(df: pd.DataFrame, edges: Sequence[float]) -> pd.DataFrame:
    labels = [bin_label(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]
    out = df.copy()
    out["distance_bin"] = pd.cut(
        out["nearest_boundary_distance"],
        bins=edges,
        labels=labels,
        right=False,
        include_lowest=True,
    )
    out["distance_bin"] = out["distance_bin"].astype(str)
    return out[out["distance_bin"] != "nan"].copy()


def load_run_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)

    required = {"accession", "base", "run_length", "nearest_boundary_distance"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns from run table: {sorted(missing)}")

    df = df.copy()
    df["run_length"] = pd.to_numeric(df["run_length"], errors="coerce")
    df["nearest_boundary_distance"] = pd.to_numeric(
        df["nearest_boundary_distance"],
        errors="coerce",
    )

    df = df.dropna(subset=["run_length", "nearest_boundary_distance", "base"])
    df = df[df["base"].isin(BASES)]
    df["run_length"] = df["run_length"].astype(int)
    df["nearest_boundary_distance"] = df["nearest_boundary_distance"].astype(float)

    return df


def compute_profile_pooled(
    df: pd.DataFrame,
    thresholds: Sequence[int],
    edges: Sequence[float],
) -> pd.DataFrame:
    """
    Compute pooled distance profiles.

    Each run is one observation. Large genomes therefore contribute more runs.
    """

    binned = assign_distance_bins(df, edges)
    labels = [bin_label(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]

    rows = []

    for base_group in ["all"] + BASES:
        if base_group == "all":
            sub = binned
        else:
            sub = binned[binned["base"] == base_group]

        for threshold in thresholds:
            for i, label in enumerate(labels):
                left = edges[i]
                right = edges[i + 1]
                g = sub[sub["distance_bin"] == label]
                n = len(g)
                n_ge = int((g["run_length"] >= threshold).sum())
                probability = n_ge / n if n else float("nan")

                rows.append(
                    {
                        "profile_type": "pooled_runs",
                        "base": base_group,
                        "threshold": threshold,
                        "bin_index": i,
                        "distance_bin": label,
                        "distance_left": left,
                        "distance_right": right,
                        "n_runs": n,
                        "n_runs_ge_threshold": n_ge,
                        "E": probability,
                    }
                )

    return pd.DataFrame(rows)


def compute_profile_by_genome(
    df: pd.DataFrame,
    thresholds: Sequence[int],
    edges: Sequence[float],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compute per-genome distance profiles and unweighted genome means.

    The per-genome table has one row per genome/base/threshold/bin.
    The summary table averages E over genomes with at least one run in the bin.
    """

    binned = assign_distance_bins(df, edges)
    labels = [bin_label(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]

    genome_rows = []

    for accession, genome_df in binned.groupby("accession"):
        organism = genome_df["organism"].iloc[0] if "organism" in genome_df.columns else ""

        for base_group in ["all"] + BASES:
            if base_group == "all":
                sub = genome_df
            else:
                sub = genome_df[genome_df["base"] == base_group]

            for threshold in thresholds:
                for i, label in enumerate(labels):
                    left = edges[i]
                    right = edges[i + 1]
                    g = sub[sub["distance_bin"] == label]
                    n = len(g)
                    n_ge = int((g["run_length"] >= threshold).sum())
                    E = n_ge / n if n else float("nan")

                    genome_rows.append(
                        {
                            "accession": accession,
                            "organism": organism,
                            "base": base_group,
                            "threshold": threshold,
                            "bin_index": i,
                            "distance_bin": label,
                            "distance_left": left,
                            "distance_right": right,
                            "n_runs": n,
                            "n_runs_ge_threshold": n_ge,
                            "E": E,
                        }
                    )

    per_genome = pd.DataFrame(genome_rows)

    summary = (
        per_genome
        .dropna(subset=["E"])
        .groupby(
            [
                "base",
                "threshold",
                "bin_index",
                "distance_bin",
                "distance_left",
                "distance_right",
            ],
            as_index=False,
        )
        .agg(
            profile_type=("E", lambda _: "genome_unweighted_mean"),
            n_genomes=("accession", "nunique"),
            mean_E=("E", "mean"),
            median_E=("E", "median"),
            sd_E=("E", "std"),
            total_runs=("n_runs", "sum"),
            total_runs_ge_threshold=("n_runs_ge_threshold", "sum"),
        )
    )

    summary = summary.rename(columns={"mean_E": "E"})

    return per_genome, summary


def near_far_summary(profile: pd.DataFrame, near_max: float, far_min: float) -> pd.DataFrame:
    """
    Summarise E(d) by comparing near-boundary and far-from-boundary bins.

    For each base and threshold:
    - near_E: weighted pooled probability for bins with distance_right <= near_max
    - far_E: weighted pooled probability for bins with distance_left >= far_min
    """

    rows = []

    # Use pooled counts where available.
    pooled = profile[profile["profile_type"] == "pooled_runs"].copy()

    for (base, threshold), sub in pooled.groupby(["base", "threshold"]):
        near = sub[sub["distance_right"] <= near_max]
        far = sub[sub["distance_left"] >= far_min]

        near_n = int(near["n_runs"].sum())
        near_ge = int(near["n_runs_ge_threshold"].sum())
        far_n = int(far["n_runs"].sum())
        far_ge = int(far["n_runs_ge_threshold"].sum())

        near_E = near_ge / near_n if near_n else float("nan")
        far_E = far_ge / far_n if far_n else float("nan")
        ratio = near_E / far_E if far_E and not math.isnan(far_E) else float("nan")
        difference = near_E - far_E if not math.isnan(near_E) and not math.isnan(far_E) else float("nan")

        # Haldane-Anscombe corrected odds ratio.
        a = near_ge
        b = near_n - near_ge
        c = far_ge
        d = far_n - far_ge
        odds_ratio = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)) if near_n and far_n else float("nan")

        rows.append(
            {
                "base": base,
                "threshold": threshold,
                "near_max": near_max,
                "far_min": far_min,
                "near_n_runs": near_n,
                "near_n_runs_ge_threshold": near_ge,
                "near_E": near_E,
                "far_n_runs": far_n,
                "far_n_runs_ge_threshold": far_ge,
                "far_E": far_E,
                "near_far_ratio": ratio,
                "near_far_difference": difference,
                "near_far_odds_ratio": odds_ratio,
            }
        )

    return pd.DataFrame(rows)


def plot_all_bases(profile: pd.DataFrame, outdir: Path) -> None:
    df = profile[
        (profile["profile_type"] == "pooled_runs") &
        (profile["base"] == "all")
    ].copy()

    plt.figure()
    for threshold, sub in df.groupby("threshold"):
        sub = sub.sort_values("bin_index")
        plt.plot(sub["bin_index"], sub["E"], marker="o", label=f"k ≥ {threshold}")

    labels = df.sort_values("bin_index")["distance_bin"].drop_duplicates().tolist()
    plt.xticks(range(len(labels)), labels, rotation=45, ha="right")
    plt.xlabel("Distance to nearest boundary / bases")
    plt.ylabel(r"$E_k(d)=P(\lambda \geq k \mid d)$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "distance_profile_all_bases.pdf", bbox_inches="tight")
    plt.close()


def plot_by_base(profile: pd.DataFrame, outdir: Path) -> None:
    df = profile[profile["profile_type"] == "pooled_runs"].copy()
    df = df[df["base"].isin(BASES)]

    for threshold, threshold_df in df.groupby("threshold"):
        plt.figure()
        for base, sub in threshold_df.groupby("base"):
            sub = sub.sort_values("bin_index")
            plt.plot(sub["bin_index"], sub["E"], marker="o", label=base)

        labels = threshold_df.sort_values("bin_index")["distance_bin"].drop_duplicates().tolist()
        plt.xticks(range(len(labels)), labels, rotation=45, ha="right")
        plt.xlabel("Distance to nearest boundary / bases")
        plt.ylabel(r"$E_k(d)=P(\lambda \geq k \mid d)$")
        plt.title(f"Threshold k ≥ {threshold}")
        plt.legend(title="Base")
        plt.tight_layout()
        plt.savefig(outdir / f"distance_profile_by_base_k{threshold}.pdf", bbox_inches="tight")
        plt.close()


def plot_near_far(summary: pd.DataFrame, outdir: Path) -> None:
    df = summary.copy()
    df = df[df["base"].isin(["all"] + BASES)]

    for threshold, sub in df.groupby("threshold"):
        sub = sub.set_index("base").loc[[b for b in ["all"] + BASES if b in set(sub["base"])]].reset_index()

        plt.figure()
        plt.bar(sub["base"], sub["near_far_odds_ratio"])
        plt.axhline(1.0, linestyle="--")
        plt.xlabel("Base")
        plt.ylabel("Near/far odds ratio")
        plt.title(f"Near/far enrichment, k ≥ {threshold}")
        plt.tight_layout()
        plt.savefig(outdir / f"near_far_odds_ratio_k{threshold}.pdf", bbox_inches="tight")
        plt.close()


def write_markdown_report(
    outdir: Path,
    profile: pd.DataFrame,
    summary: pd.DataFrame,
    thresholds: Sequence[int],
    bins: Sequence[float],
) -> None:
    report = []
    report.append("# Distance-dependent boundary enrichment report\n")
    report.append("This report estimates:\n")
    report.append(r"\[")
    report.append(r"E_k(d)=P(\lambda\ge k\mid d),")
    report.append(r"\]")
    report.append("where \\(d\\) is the distance from a homopolymer run to the nearest annotated feature boundary.\n")

    report.append("## Parameters\n")
    report.append(f"- Thresholds: {', '.join(str(k) for k in thresholds)}")
    report.append(f"- Distance bins: {', '.join('inf' if math.isinf(x) else str(int(x)) for x in bins)}\n")

    report.append("## Near/far summary\n")
    cols = [
        "base",
        "threshold",
        "near_E",
        "far_E",
        "near_far_ratio",
        "near_far_difference",
        "near_far_odds_ratio",
        "near_n_runs",
        "far_n_runs",
    ]
    report.append(summary[cols].to_markdown(index=False))
    report.append("\n")

    report.append("## Output files\n")
    report.append("- `distance_profile.csv`")
    report.append("- `distance_profile_by_genome.csv`")
    report.append("- `distance_profile_summary.csv`")
    report.append("- `distance_profile_all_bases.pdf`")
    report.append("- `distance_profile_by_base_k*.pdf`")
    report.append("- `near_far_odds_ratio_k*.pdf`\n")

    (outdir / "distance_profile_report.md").write_text("\n".join(report), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute E_k(d)=P(run length >= k | distance to nearest boundary)."
    )

    parser.add_argument(
        "--run-table",
        type=Path,
        required=True,
        help="run_table.csv from boundary_enrichment_v3.py.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Output directory.",
    )
    parser.add_argument(
        "--thresholds",
        default="5,10,15",
        help="Comma-separated run-length thresholds. Default: 5,10,15.",
    )
    parser.add_argument(
        "--bins",
        default="0,25,50,100,200,400,800,1600,inf",
        help="Comma-separated distance-bin edges. Use inf for infinity.",
    )
    parser.add_argument(
        "--near-max",
        type=float,
        default=100,
        help="Maximum distance defining near-boundary region for summary. Default: 100.",
    )
    parser.add_argument(
        "--far-min",
        type=float,
        default=400,
        help="Minimum distance defining far-from-boundary region for summary. Default: 400.",
    )
    parser.add_argument(
        "--normalise-by-genome",
        action="store_true",
        help="Also write genome-unweighted mean profile as the main interpreted sensitivity check.",
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    thresholds = parse_thresholds(args.thresholds)
    bins = parse_bins(args.bins)

    run_table = load_run_table(args.run_table)

    pooled = compute_profile_pooled(run_table, thresholds, bins)
    per_genome, genome_mean = compute_profile_by_genome(run_table, thresholds, bins)

    # Combine pooled and genome-unweighted profiles into one table.
    genome_mean_for_export = genome_mean.copy()
    if not genome_mean_for_export.empty:
        genome_mean_for_export["n_runs"] = genome_mean_for_export["total_runs"]
        genome_mean_for_export["n_runs_ge_threshold"] = genome_mean_for_export["total_runs_ge_threshold"]

    common_cols = [
        "profile_type",
        "base",
        "threshold",
        "bin_index",
        "distance_bin",
        "distance_left",
        "distance_right",
        "n_runs",
        "n_runs_ge_threshold",
        "E",
    ]

    profile = pd.concat(
        [
            pooled[common_cols],
            genome_mean_for_export[[c for c in common_cols if c in genome_mean_for_export.columns]],
        ],
        ignore_index=True,
    )

    summary = near_far_summary(profile, args.near_max, args.far_min)

    profile.to_csv(args.outdir / "distance_profile.csv", index=False)
    per_genome.to_csv(args.outdir / "distance_profile_by_genome.csv", index=False)
    summary.to_csv(args.outdir / "distance_profile_summary.csv", index=False)

    plot_all_bases(profile, args.outdir)
    plot_by_base(profile, args.outdir)
    plot_near_far(summary, args.outdir)

    try:
        write_markdown_report(args.outdir, profile, summary, thresholds, bins)
    except Exception as exc:
        print(f"WARNING: could not write Markdown report: {exc}")

    print(f"Wrote distance-profile outputs to {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
