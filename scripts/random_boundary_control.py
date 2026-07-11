#!/usr/bin/env python3
"""
random_boundary_control.py

Random-boundary control for distance-dependent homopolymer enrichment.

This script tests whether the observed distance profile

    E_k(d) = P(run length >= k | distance from nearest boundary is in bin d)

is stronger than expected if feature-boundary positions were randomised within
each genome.

It uses the `run_table.csv` produced by boundary_enrichment_v3.py.  That file
already contains all homopolymer runs and the actual boundary hits/nearest
boundaries.  To build the null model, this script reconstructs the observed
boundary positions from the `boundary_hits` and `nearest_boundary` columns,
then repeatedly randomises the same number of boundary positions within each
genome.

No NCBI download is needed.

Inputs
------
Required:
- run_table.csv from boundary_enrichment_v3.py

Expected columns:
- accession
- sequence_length
- base
- run_length
- run_start_1based
- run_end_1based
- nearest_boundary_distance
- boundary_hits
- nearest_boundary

Outputs
-------
- observed_distance_profile.csv
- random_boundary_null_profiles.csv
- random_boundary_control_summary.csv
- random_boundary_control_report.md
- random_boundary_control_all_bases.pdf
- random_boundary_control_by_base_k*.pdf

Example
-------
python random_boundary_control.py \
    --run-table analysis_initial/boundary_enrichment_v3/run_table.csv \
    --outdir analysis_initial/random_boundary_control \
    --thresholds 5,10,15 \
    --bins 0,25,50,100,200,400,800,1600,inf \
    --n-random 500 \
    --seed 12345

Interpretation
--------------
For each threshold/base/bin, the script reports:

- observed_E
- null_mean_E
- null_sd_E
- z_score
- empirical_p_upper

where empirical_p_upper is the fraction of randomisations with E >= observed_E,
with a +1 correction.
"""

from __future__ import annotations

import argparse
import bisect
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASES = ["A", "C", "G", "T"]


BOUNDARY_HIT_RE = re.compile(r":(?:start|end):(\d+)(?:;|$)")
NEAREST_BOUNDARY_RE = re.compile(r":(?:start|end):(\d+):feature=")


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


def load_run_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)

    required = {
        "accession",
        "sequence_length",
        "base",
        "run_length",
        "run_start_1based",
        "run_end_1based",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    df = df.copy()
    for col in ["sequence_length", "run_length", "run_start_1based", "run_end_1based"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["sequence_length", "run_length", "run_start_1based", "run_end_1based", "base"])
    df = df[df["base"].isin(BASES)]

    df["sequence_length"] = df["sequence_length"].astype(int)
    df["run_length"] = df["run_length"].astype(int)
    df["run_start_1based"] = df["run_start_1based"].astype(int)
    df["run_end_1based"] = df["run_end_1based"].astype(int)

    for col in ["boundary_hits", "nearest_boundary", "nearest_boundary_distance"]:
        if col not in df.columns:
            df[col] = ""

    return df


def extract_boundary_positions_from_text(text: str) -> list[int]:
    if not isinstance(text, str) or not text:
        return []
    out = []
    out.extend(int(m.group(1)) for m in BOUNDARY_HIT_RE.finditer(text))
    out.extend(int(m.group(1)) for m in NEAREST_BOUNDARY_RE.finditer(text))
    return out


def infer_boundaries_from_run_table(df: pd.DataFrame) -> dict[str, list[int]]:
    """
    Infer observed boundary positions per accession from textual fields in run_table.

    The run_table does not store a separate boundary list, but each run contains
    boundary hits and nearest-boundary descriptions.  Across all runs this
    recovers the boundary positions well enough for the random-position control.
    """

    boundaries: dict[str, set[int]] = defaultdict(set)

    for row in df[["accession", "boundary_hits", "nearest_boundary"]].itertuples(index=False):
        accession = str(row.accession)

        for pos in extract_boundary_positions_from_text(str(row.boundary_hits)):
            boundaries[accession].add(pos)

        for pos in extract_boundary_positions_from_text(str(row.nearest_boundary)):
            boundaries[accession].add(pos)

    return {acc: sorted(pos) for acc, pos in boundaries.items()}


def distance_to_nearest_boundary_for_runs(
    starts: np.ndarray,
    ends: np.ndarray,
    boundaries: Sequence[int],
) -> np.ndarray:
    """
    Compute distance from each run interval to the nearest boundary point.

    Distance is zero when the boundary point lies inside the run.
    """

    if len(boundaries) == 0:
        return np.full(len(starts), np.inf)

    b = np.array(sorted(boundaries), dtype=np.int64)
    distances = np.empty(len(starts), dtype=float)

    for i, (s, e) in enumerate(zip(starts, ends)):
        idx = bisect.bisect_left(b, int((s + e) // 2))
        candidates = []
        if idx > 0:
            candidates.append(b[idx - 1])
        if idx < len(b):
            candidates.append(b[idx])
        # Also check insertion around start/end in case a boundary lies inside a long run.
        idx_s = bisect.bisect_left(b, int(s))
        if idx_s < len(b):
            candidates.append(b[idx_s])
        if idx_s > 0:
            candidates.append(b[idx_s - 1])

        best = math.inf
        for p in candidates:
            if s <= p <= e:
                best = 0
                break
            if p < s:
                d = s - p
            else:
                d = p - e
            if d < best:
                best = d
        distances[i] = best

    return distances


def assign_bins(distances: np.ndarray, edges: Sequence[float]) -> np.ndarray:
    """
    Assign distances to bin indices.

    Returns -1 for distances outside bins.
    """

    out = np.full(len(distances), -1, dtype=int)
    for i in range(len(edges) - 1):
        left = edges[i]
        right = edges[i + 1]
        mask = (distances >= left) & (distances < right)
        out[mask] = i
    return out


def compute_profile_from_distances(
    df: pd.DataFrame,
    distances: np.ndarray,
    thresholds: Sequence[int],
    edges: Sequence[float],
) -> pd.DataFrame:
    """
    Compute E_k(d) for all bases and base-specific groups from distances.
    """

    labels = [bin_label(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]
    bin_indices = assign_bins(distances, edges)

    work = df[["accession", "base", "run_length"]].copy()
    work["distance"] = distances
    work["bin_index"] = bin_indices
    work = work[work["bin_index"] >= 0].copy()

    rows = []

    for base_group in ["all"] + BASES:
        if base_group == "all":
            sub = work
        else:
            sub = work[work["base"] == base_group]

        for threshold in thresholds:
            for i, label in enumerate(labels):
                g = sub[sub["bin_index"] == i]
                n = len(g)
                n_ge = int((g["run_length"] >= threshold).sum())
                rows.append(
                    {
                        "base": base_group,
                        "threshold": threshold,
                        "bin_index": i,
                        "distance_bin": label,
                        "distance_left": edges[i],
                        "distance_right": edges[i + 1],
                        "n_runs": n,
                        "n_runs_ge_threshold": n_ge,
                        "E": n_ge / n if n else np.nan,
                    }
                )

    return pd.DataFrame(rows)


def compute_observed_distances(df: pd.DataFrame, boundaries_by_acc: dict[str, list[int]]) -> np.ndarray:
    """
    Compute observed nearest-boundary distances.

    Prefer the numeric nearest_boundary_distance column if valid.  If absent or
    partially missing, recompute from inferred boundary positions.
    """

    numeric = pd.to_numeric(df.get("nearest_boundary_distance", ""), errors="coerce")
    if numeric.notna().all():
        return numeric.to_numpy(dtype=float)

    distances = np.empty(len(df), dtype=float)

    offset = 0
    for accession, sub in df.groupby("accession", sort=False):
        n = len(sub)
        starts = sub["run_start_1based"].to_numpy()
        ends = sub["run_end_1based"].to_numpy()
        boundaries = boundaries_by_acc.get(accession, [])
        distances[offset:offset+n] = distance_to_nearest_boundary_for_runs(starts, ends, boundaries)
        offset += n

    return distances


def randomise_boundaries_for_accession(
    sequence_length: int,
    n_boundaries: int,
    rng: np.random.Generator,
) -> list[int]:
    """
    Generate random boundary positions within a genome.

    Sampling is without replacement when feasible.
    """

    if n_boundaries <= 0:
        return []
    n_boundaries = min(n_boundaries, sequence_length)

    positions = rng.choice(
        np.arange(1, sequence_length + 1, dtype=np.int64),
        size=n_boundaries,
        replace=False,
    )
    positions.sort()
    return positions.tolist()


def compute_random_distances(
    df: pd.DataFrame,
    boundary_counts: dict[str, int],
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Randomise boundary positions within each genome and compute run distances.
    """

    distances = np.empty(len(df), dtype=float)
    offset = 0

    for accession, sub in df.groupby("accession", sort=False):
        n = len(sub)
        seq_len = int(sub["sequence_length"].iloc[0])
        n_boundaries = boundary_counts.get(accession, 0)
        random_boundaries = randomise_boundaries_for_accession(seq_len, n_boundaries, rng)

        starts = sub["run_start_1based"].to_numpy()
        ends = sub["run_end_1based"].to_numpy()
        distances[offset:offset+n] = distance_to_nearest_boundary_for_runs(starts, ends, random_boundaries)

        offset += n

    return distances


def summarise_null(
    observed: pd.DataFrame,
    null_profiles: pd.DataFrame,
) -> pd.DataFrame:
    """
    Compare observed E against random-boundary null distribution.
    """

    group_cols = ["base", "threshold", "bin_index", "distance_bin", "distance_left", "distance_right"]

    null_summary = (
        null_profiles
        .groupby(group_cols, as_index=False)
        .agg(
            null_mean_E=("E", "mean"),
            null_sd_E=("E", "std"),
            null_median_E=("E", "median"),
            null_q025_E=("E", lambda x: np.nanquantile(x, 0.025)),
            null_q975_E=("E", lambda x: np.nanquantile(x, 0.975)),
            null_mean_n_runs=("n_runs", "mean"),
        )
    )

    merged = observed.merge(null_summary, on=group_cols, how="left")
    merged = merged.rename(
        columns={
            "E": "observed_E",
            "n_runs": "observed_n_runs",
            "n_runs_ge_threshold": "observed_n_runs_ge_threshold",
        }
    )

    # Empirical upper-tail p-value: P(null >= observed), +1 corrected.
    p_rows = []
    for key, obs_row in merged.groupby(group_cols, sort=False):
        if not isinstance(key, tuple):
            key = (key,)
        cond = np.ones(len(null_profiles), dtype=bool)
        for col, val in zip(group_cols, key):
            cond &= (null_profiles[col].to_numpy() == val)
        vals = null_profiles.loc[cond, "E"].to_numpy(dtype=float)
        obs = float(obs_row["observed_E"].iloc[0])
        vals = vals[~np.isnan(vals)]
        if len(vals) == 0 or np.isnan(obs):
            p_upper = np.nan
            p_lower = np.nan
        else:
            p_upper = (1 + np.sum(vals >= obs)) / (len(vals) + 1)
            p_lower = (1 + np.sum(vals <= obs)) / (len(vals) + 1)
        p_rows.append((*key, p_upper, p_lower))

    p_df = pd.DataFrame(p_rows, columns=group_cols + ["empirical_p_upper", "empirical_p_lower"])
    merged = merged.merge(p_df, on=group_cols, how="left")

    merged["z_score"] = (
        (merged["observed_E"] - merged["null_mean_E"]) / merged["null_sd_E"]
    )
    merged.loc[merged["null_sd_E"] == 0, "z_score"] = np.nan
    merged["observed_minus_null"] = merged["observed_E"] - merged["null_mean_E"]
    merged["observed_over_null"] = merged["observed_E"] / merged["null_mean_E"]

    return merged


def plot_observed_vs_null(summary: pd.DataFrame, outdir: Path) -> None:
    all_df = summary[summary["base"] == "all"].copy()

    plt.figure()
    for threshold, sub in all_df.groupby("threshold"):
        sub = sub.sort_values("bin_index")
        x = np.arange(len(sub))
        plt.plot(x, sub["observed_E"], marker="o", label=f"obs k≥{threshold}")
        plt.plot(x, sub["null_mean_E"], marker="x", linestyle="--", label=f"null k≥{threshold}")

    labels = all_df.sort_values("bin_index")["distance_bin"].drop_duplicates().tolist()
    plt.xticks(range(len(labels)), labels, rotation=45, ha="right")
    plt.xlabel("Distance to nearest boundary / bases")
    plt.ylabel(r"$E_k(d)$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "random_boundary_control_all_bases.pdf", bbox_inches="tight")
    plt.close()

    for threshold, threshold_df in summary[summary["base"].isin(BASES)].groupby("threshold"):
        plt.figure()
        for base, sub in threshold_df.groupby("base"):
            sub = sub.sort_values("bin_index")
            x = np.arange(len(sub))
            plt.plot(x, sub["observed_E"], marker="o", label=f"{base} observed")
            plt.plot(x, sub["null_mean_E"], linestyle="--", label=f"{base} null")

        labels = threshold_df.sort_values("bin_index")["distance_bin"].drop_duplicates().tolist()
        plt.xticks(range(len(labels)), labels, rotation=45, ha="right")
        plt.xlabel("Distance to nearest boundary / bases")
        plt.ylabel(r"$E_k(d)$")
        plt.title(f"Random-boundary control, k≥{threshold}")
        plt.legend(ncol=2)
        plt.tight_layout()
        plt.savefig(outdir / f"random_boundary_control_by_base_k{threshold}.pdf", bbox_inches="tight")
        plt.close()


def write_report(
    outdir: Path,
    summary: pd.DataFrame,
    n_random: int,
    thresholds: Sequence[int],
    bins: Sequence[float],
) -> None:
    report = []
    report.append("# Random-boundary control report\n")
    report.append("This analysis randomises annotated boundary positions within each genome while preserving the number of boundaries per genome.\n")
    report.append("It tests whether the observed distance-dependent homopolymer profile exceeds the random-boundary expectation.\n")

    report.append("## Parameters\n")
    report.append(f"- Randomisations: {n_random}")
    report.append(f"- Thresholds: {', '.join(str(k) for k in thresholds)}")
    report.append(f"- Distance bins: {', '.join('inf' if math.isinf(x) else str(int(x)) for x in bins)}\n")

    report.append("## Key 0-25 bp bin results\n")
    key = summary[summary["distance_bin"] == "0-25"].copy()
    cols = [
        "base",
        "threshold",
        "observed_E",
        "null_mean_E",
        "observed_minus_null",
        "observed_over_null",
        "z_score",
        "empirical_p_upper",
        "observed_n_runs",
    ]
    report.append(key[cols].to_markdown(index=False))
    report.append("\n")

    report.append("## Output files\n")
    report.append("- `observed_distance_profile.csv`")
    report.append("- `random_boundary_null_profiles.csv`")
    report.append("- `random_boundary_control_summary.csv`")
    report.append("- `random_boundary_control_all_bases.pdf`")
    report.append("- `random_boundary_control_by_base_k*.pdf`\n")

    (outdir / "random_boundary_control_report.md").write_text("\n".join(report), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Random-boundary control for distance-dependent homopolymer enrichment."
    )

    parser.add_argument(
        "--run-table",
        type=Path,
        required=True,
        help="run_table.csv produced by boundary_enrichment_v3.py.",
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
        "--n-random",
        type=int,
        default=500,
        help="Number of random boundary realisations. Default: 500.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=12345,
        help="Random seed. Default: 12345.",
    )
    parser.add_argument(
        "--save-null-profiles",
        action="store_true",
        help="Save every random profile. If omitted, the file is still saved, but this option is retained for clarity.",
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    thresholds = parse_thresholds(args.thresholds)
    bins = parse_bins(args.bins)

    print(f"Reading {args.run_table}")
    run_table = load_run_table(args.run_table)

    print("Inferring observed boundary positions from run_table")
    boundaries_by_acc = infer_boundaries_from_run_table(run_table)
    boundary_counts = {acc: len(pos) for acc, pos in boundaries_by_acc.items()}

    n_missing = sum(1 for acc in run_table["accession"].unique() if boundary_counts.get(acc, 0) == 0)
    if n_missing:
        print(f"WARNING: {n_missing} accession(s) have no inferred boundary positions")

    print("Computing observed distance profile")
    observed_distances = compute_observed_distances(run_table, boundaries_by_acc)
    observed = compute_profile_from_distances(run_table, observed_distances, thresholds, bins)
    observed.to_csv(args.outdir / "observed_distance_profile.csv", index=False)

    rng = np.random.default_rng(args.seed)

    null_profiles = []
    print(f"Running {args.n_random} random-boundary controls")
    for r in range(1, args.n_random + 1):
        if r == 1 or r % 25 == 0 or r == args.n_random:
            print(f"  randomisation {r}/{args.n_random}", flush=True)

        random_distances = compute_random_distances(run_table, boundary_counts, rng)
        profile = compute_profile_from_distances(run_table, random_distances, thresholds, bins)
        profile.insert(0, "randomisation", r)
        null_profiles.append(profile)

    null_profiles_df = pd.concat(null_profiles, ignore_index=True)
    null_profiles_df.to_csv(args.outdir / "random_boundary_null_profiles.csv", index=False)

    print("Summarising observed versus random-boundary null")
    summary = summarise_null(observed, null_profiles_df)
    summary.to_csv(args.outdir / "random_boundary_control_summary.csv", index=False)

    plot_observed_vs_null(summary, args.outdir)
    try:
        write_report(args.outdir, summary, args.n_random, thresholds, bins)
    except Exception as exc:
        print(f"WARNING: could not write Markdown report: {exc}")

    print(f"Wrote random-boundary control outputs to {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
