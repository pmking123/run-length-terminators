#!/usr/bin/env python3
"""
boundary_enrichment_v3.py

Fast systematic boundary-enrichment analysis for homopolymer runs.

This version fixes the processing bottleneck seen with records such as
NC_115876.1.  The earlier script checked feature overlaps using a nested loop
over all runs and all features, which can become extremely slow for large or
feature-rich GenBank records.

Version 3 changes
-----------------
- Uses boundary positions plus binary search for boundary-window lookup.
- Uses a sweep-line interval index for feature-overlap lookup.
- Adds --max-sequence-length to skip very large records if desired.
- Adds per-record timing diagnostics.
- Keeps timeout/caching/resume functionality from v2.

Typical command
---------------
python boundary_enrichment_v3.py \
    --features analysis_initial/results/features.csv \
    --email your.name@example.com \
    --outdir analysis_initial/boundary_enrichment_v3 \
    --boundary-window 250 \
    --cache-dir analysis_initial/boundary_enrichment_v3/gb_cache \
    --resume

Optional, while diagnosing:
---------------------------
python boundary_enrichment_v3.py \
    --features analysis_initial/results/features.csv \
    --email your.name@example.com \
    --outdir analysis_initial/boundary_enrichment_debug \
    --accessions NC_115876.1 \
    --boundary-window 250 \
    --verbose
"""

from __future__ import annotations

import argparse
import bisect
import csv
import math
import sys
import time
from collections import Counter
from dataclasses import dataclass
from multiprocessing import Process, Queue
from pathlib import Path
from queue import Empty
from typing import Iterable, Iterator, Sequence

from Bio import Entrez, SeqIO
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature


VALID_BASES = frozenset("ACGT")


@dataclass(frozen=True)
class Run:
    base: str
    start0: int
    end0: int

    @property
    def length(self) -> int:
        return self.end0 - self.start0

    @property
    def start1(self) -> int:
        return self.start0 + 1

    @property
    def end1(self) -> int:
        return self.end0

    @property
    def midpoint1(self) -> float:
        return (self.start1 + self.end1) / 2.0


@dataclass(frozen=True)
class Boundary:
    accession: str
    position1: int
    feature_type: str
    feature_label: str
    boundary_kind: str
    feature_start1: int
    feature_end1: int
    strand: str


@dataclass(frozen=True)
class FeatureSpan:
    start1: int
    end1: int
    strand: str


@dataclass(frozen=True)
class FeatureInterval:
    start1: int
    end1: int
    feature_type: str
    feature_label: str


@dataclass(frozen=True)
class Candidate:
    accession: str
    source_length: int | None
    source_longest_run: int | None
    source_gc_frequency: float | None


def print_flush(*args, **kwargs) -> None:
    print(*args, **kwargs, flush=True)


def unique_preserving_order(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for v in values:
        v = str(v).strip()
        if v and v not in seen:
            seen.add(v)
            out.append(v)
    return out


def detect_column(fieldnames: Sequence[str], candidates: Sequence[str]) -> str | None:
    lowered = {f.lower(): f for f in fieldnames}
    for c in candidates:
        if c.lower() in lowered:
            return lowered[c.lower()]
    return None


def parse_optional_int(value: str | None) -> int | None:
    if value is None or str(value).strip() == "":
        return None
    return int(float(str(value)))


def parse_optional_float(value: str | None) -> float | None:
    if value is None or str(value).strip() == "":
        return None
    return float(str(value))


def read_existing_accessions(path: Path) -> set[str]:
    if not path.exists() or path.stat().st_size == 0:
        return set()
    try:
        with path.open("r", newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None or "accession" not in reader.fieldnames:
                return set()
            return {row["accession"] for row in reader if row.get("accession")}
    except Exception:
        return set()


def read_candidates(
    features_csv: Path,
    *,
    id_column: str | None,
    length_column: str | None,
    longest_run_column: str | None,
    gc_column: str | None,
    accessions: Sequence[str] | None,
    limit: int | None,
) -> list[Candidate]:
    with features_csv.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"{features_csv} has no header row")

        fields = reader.fieldnames
        id_col = id_column or detect_column(fields, ["id", "accession", "seq_id", "record_id"])
        len_col = length_column or detect_column(fields, ["length", "sequence_length", "genome_length"])
        lr_col = longest_run_column or detect_column(fields, ["longest_run", "max_run", "overall_longest_run"])
        gc_col = gc_column or detect_column(fields, ["GC_frequency", "gc_frequency", "gc", "gc_fraction"])

        if id_col is None:
            raise ValueError("Could not detect accession column. Use --id-column.")

        rows = list(reader)

    by_acc: dict[str, Candidate] = {}
    for row in rows:
        acc = (row.get(id_col) or "").strip()
        if not acc:
            continue
        by_acc[acc] = Candidate(
            accession=acc,
            source_length=parse_optional_int(row.get(len_col)) if len_col else None,
            source_longest_run=parse_optional_int(row.get(lr_col)) if lr_col else None,
            source_gc_frequency=parse_optional_float(row.get(gc_col)) if gc_col else None,
        )

    if accessions:
        selected = [by_acc.get(a, Candidate(a, None, None, None)) for a in unique_preserving_order(accessions)]
    else:
        selected = list(by_acc.values())

    if limit is not None:
        selected = selected[:limit]

    return selected


def cache_path_for_accession(cache_dir: Path, accession: str) -> Path:
    safe = accession.replace("/", "_").replace("\\", "_")
    return cache_dir / f"{safe}.gb"


def _fetch_worker(accession: str, email: str, api_key: str | None, output_path: str, queue: Queue) -> None:
    try:
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as handle:
            text = handle.read()
        if not text.strip():
            raise RuntimeError("empty GenBank record")
        Path(output_path).write_text(text, encoding="utf-8")
        queue.put(("ok", output_path))
    except Exception as exc:
        queue.put(("error", repr(exc)))


def fetch_genbank_with_timeout(
    accession: str,
    *,
    email: str,
    api_key: str | None,
    cache_dir: Path,
    timeout: int,
    max_retries: int,
    delay: float,
    force_refetch: bool,
):
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_path_for_accession(cache_dir, accession)

    if path.exists() and path.stat().st_size > 0 and not force_refetch:
        try:
            return SeqIO.read(str(path), "genbank")
        except Exception:
            path.unlink(missing_ok=True)

    tmp = path.with_suffix(".tmp.gb")
    last_error = ""

    for attempt in range(1, max_retries + 1):
        q: Queue = Queue()
        p = Process(target=_fetch_worker, args=(accession, email, api_key, str(tmp), q))
        p.start()
        p.join(timeout)

        if p.is_alive():
            p.terminate()
            p.join(5)
            last_error = f"fetch timeout after {timeout}s on attempt {attempt}"
        else:
            try:
                status, payload = q.get_nowait()
            except Empty:
                status, payload = "error", "worker exited without status"

            if status == "ok":
                tmp.replace(path)
                if delay > 0:
                    time.sleep(delay)
                return SeqIO.read(str(path), "genbank")
            last_error = f"attempt {attempt}: {payload}"

        tmp.unlink(missing_ok=True)
        if delay > 0:
            time.sleep(delay)

    raise RuntimeError(last_error or "fetch failed")


def iter_runs(sequence: str) -> Iterator[Run]:
    sequence = sequence.upper()
    n = len(sequence)
    i = 0
    while i < n:
        b = sequence[i]
        if b not in VALID_BASES:
            i += 1
            continue
        j = i + 1
        while j < n and sequence[j] == b:
            j += 1
        yield Run(base=b, start0=i, end0=j)
        i = j


def strand_to_string(strand: int | None) -> str:
    if strand == 1:
        return "+"
    if strand == -1:
        return "-"
    return "."


def feature_spans(feature: SeqFeature) -> list[FeatureSpan]:
    loc = feature.location
    parts = loc.parts if isinstance(loc, CompoundLocation) else [loc]
    spans = []
    for p in parts:
        if isinstance(p, FeatureLocation):
            spans.append(FeatureSpan(start1=int(p.start) + 1, end1=int(p.end), strand=strand_to_string(p.strand)))
    return spans


def qualifier_first(feature: SeqFeature, key: str) -> str:
    vals = feature.qualifiers.get(key, [])
    return str(vals[0]) if vals else ""


def feature_label(feature: SeqFeature) -> str:
    for key in ("gene", "product", "locus_tag", "protein_id", "note"):
        val = qualifier_first(feature, key)
        if val:
            return val
    return feature.type


def source_metadata(record) -> dict[str, str]:
    meta = {
        "description": record.description,
        "organism": "",
        "taxonomy": "; ".join(record.annotations.get("taxonomy", [])),
        "molecule_type": record.annotations.get("molecule_type", ""),
        "topology": record.annotations.get("topology", ""),
    }
    for f in record.features:
        if f.type == "source":
            for k in ("organism", "host", "isolate", "country", "collection_date"):
                meta[k] = qualifier_first(f, k)
            break
    return meta


def interval_distance_to_point(start1: int, end1: int, point1: int) -> int:
    if start1 <= point1 <= end1:
        return 0
    if point1 < start1:
        return start1 - point1
    return point1 - end1


def intervals_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start <= b_end and b_start <= a_end


def extract_boundaries_and_intervals(
    record,
    accession: str,
    boundary_types: set[str],
) -> tuple[list[Boundary], list[FeatureInterval]]:
    boundaries: list[Boundary] = []
    intervals: list[FeatureInterval] = []
    seq_len = len(record.seq)

    for feature in record.features:
        if feature.type == "source":
            continue

        label = feature_label(feature)

        for span in feature_spans(feature):
            if span.start1 < 1 or span.end1 > seq_len:
                continue

            intervals.append(
                FeatureInterval(
                    start1=span.start1,
                    end1=span.end1,
                    feature_type=feature.type,
                    feature_label=label,
                )
            )

            if boundary_types and feature.type not in boundary_types:
                continue

            boundaries.append(
                Boundary(
                    accession=accession,
                    position1=span.start1,
                    feature_type=feature.type,
                    feature_label=label,
                    boundary_kind="start",
                    feature_start1=span.start1,
                    feature_end1=span.end1,
                    strand=span.strand,
                )
            )
            boundaries.append(
                Boundary(
                    accession=accession,
                    position1=span.end1,
                    feature_type=feature.type,
                    feature_label=label,
                    boundary_kind="end",
                    feature_start1=span.start1,
                    feature_end1=span.end1,
                    strand=span.strand,
                )
            )

    # Deduplicate boundaries.
    seen = set()
    unique_boundaries = []
    for b in boundaries:
        key = (b.position1, b.feature_type, b.feature_label, b.boundary_kind, b.feature_start1, b.feature_end1, b.strand)
        if key not in seen:
            seen.add(key)
            unique_boundaries.append(b)

    unique_boundaries.sort(key=lambda b: b.position1)
    intervals.sort(key=lambda x: x.start1)

    return unique_boundaries, intervals


def prepare_boundary_index(boundaries: Sequence[Boundary]) -> tuple[list[int], list[Boundary]]:
    sorted_boundaries = sorted(boundaries, key=lambda b: b.position1)
    return [b.position1 for b in sorted_boundaries], list(sorted_boundaries)


def boundary_window_lookup(
    run: Run,
    positions: Sequence[int],
    boundaries: Sequence[Boundary],
    window: int,
) -> tuple[bool, str, int | None, str]:
    """
    Fast lookup of boundary hits within window and nearest boundary.

    Boundary hit criterion: boundary point lies within [run.start-window, run.end+window].
    """

    if not positions:
        return False, "", None, ""

    left = run.start1 - window
    right = run.end1 + window

    lo = bisect.bisect_left(positions, left)
    hi = bisect.bisect_right(positions, right)
    hits = boundaries[lo:hi]

    hit_text = ";".join(
        f"{b.feature_type}:{b.feature_label}:{b.boundary_kind}:{b.position1}"
        for b in hits
    )

    # nearest point using insertion position around run midpoint/start.
    idx = bisect.bisect_left(positions, int(run.midpoint1))
    candidate_indices = [idx - 1, idx, lo, hi - 1]
    best = None
    for j in candidate_indices:
        if 0 <= j < len(boundaries):
            b = boundaries[j]
            d = interval_distance_to_point(run.start1, run.end1, b.position1)
            if best is None or d < best[0]:
                best = (d, b)

    if best is None:
        return bool(hits), hit_text, None, ""

    dist, b = best
    desc = (
        f"{b.feature_type}:{b.feature_label}:{b.boundary_kind}:{b.position1}:"
        f"feature={b.feature_start1}-{b.feature_end1}:strand={b.strand}"
    )
    return bool(hits), hit_text, dist, desc


def feature_overlap_sweep(
    runs: Sequence[Run],
    intervals: Sequence[FeatureInterval],
) -> list[tuple[str, str]]:
    """
    Return overlapping feature types/labels for each run using a sweep-line scan.

    This replaces the expensive O(n_runs * n_features) nested loop.
    """

    out: list[tuple[str, str]] = []
    active: list[FeatureInterval] = []
    j = 0
    intervals = list(intervals)

    for run in runs:
        # Add intervals that have started by run.end1.
        while j < len(intervals) and intervals[j].start1 <= run.end1:
            active.append(intervals[j])
            j += 1

        # Keep intervals that have not ended before run.start1.
        active = [iv for iv in active if iv.end1 >= run.start1]

        types = []
        labels = []
        for iv in active:
            if intervals_overlap(run.start1, run.end1, iv.start1, iv.end1):
                types.append(iv.feature_type)
                labels.append(f"{iv.feature_type}:{iv.feature_label}:{iv.start1}-{iv.end1}")

        out.append((";".join(sorted(set(types))), ";".join(sorted(set(labels)))))

    return out


def shannon_entropy(lengths: Sequence[int]) -> float:
    if not lengths:
        return 0.0
    c = Counter(lengths)
    n = sum(c.values())
    return -sum((v / n) * math.log2(v / n) for v in c.values())


def mean(values: Sequence[int | float]) -> float:
    return float(sum(values) / len(values)) if values else 0.0


def median(values: Sequence[int | float]) -> float:
    if not values:
        return 0.0
    s = sorted(values)
    m = len(s) // 2
    return float(s[m]) if len(s) % 2 else float((s[m - 1] + s[m]) / 2)


def summarise_run_subset(runs: Sequence[Run], prefix: str, thresholds: Sequence[int]) -> dict[str, int | float]:
    lengths = [r.length for r in runs]
    out: dict[str, int | float] = {
        f"{prefix}_n_runs": len(runs),
        f"{prefix}_total_bases_in_runs": sum(lengths),
        f"{prefix}_mean_run_length": mean(lengths),
        f"{prefix}_median_run_length": median(lengths),
        f"{prefix}_max_run_length": max(lengths) if lengths else 0,
        f"{prefix}_run_length_entropy": shannon_entropy(lengths),
    }

    for b in sorted(VALID_BASES):
        bruns = [r for r in runs if r.base == b]
        blens = [r.length for r in bruns]
        out[f"{prefix}_{b}_n_runs"] = len(bruns)
        out[f"{prefix}_{b}_mean_run_length"] = mean(blens)
        out[f"{prefix}_{b}_max_run_length"] = max(blens) if blens else 0

    for k in thresholds:
        n_ge = sum(1 for x in lengths if x >= k)
        out[f"{prefix}_n_runs_ge_{k}"] = n_ge
        out[f"{prefix}_fraction_runs_ge_{k}"] = n_ge / len(lengths) if lengths else 0.0

    return out


def odds_ratio(a: int, b: int, c: int, d: int, pseudocount: float = 0.5) -> float:
    return ((a + pseudocount) * (d + pseudocount)) / ((b + pseudocount) * (c + pseudocount))


def analyse_record(
    *,
    candidate: Candidate,
    record,
    boundary_types: set[str],
    boundary_window: int,
    thresholds: Sequence[int],
    verbose: bool,
) -> tuple[list[dict], dict, list[dict]]:
    t0 = time.perf_counter()
    acc = candidate.accession
    seq = str(record.seq).upper()
    seq_len = len(seq)
    meta = source_metadata(record)

    if verbose:
        print_flush(f"    extracting runs...", end="")
    runs = list(iter_runs(seq))
    if verbose:
        print_flush(f" {len(runs)}")

    if verbose:
        print_flush(f"    extracting boundaries/features...", end="")
    boundaries, intervals = extract_boundaries_and_intervals(record, acc, boundary_types)
    positions, sorted_boundaries = prepare_boundary_index(boundaries)
    if verbose:
        print_flush(f" boundaries={len(boundaries)}, intervals={len(intervals)}")

    if verbose:
        print_flush("    computing feature overlaps...")
    overlap_info = feature_overlap_sweep(runs, intervals)

    run_rows = []
    boundary_runs = []
    nonboundary_runs = []

    if verbose:
        print_flush("    classifying runs...")
    for idx, run in enumerate(runs, start=1):
        in_window, boundary_hits, nearest_dist, nearest_desc = boundary_window_lookup(
            run,
            positions,
            sorted_boundaries,
            boundary_window,
        )

        ftypes, flabels = overlap_info[idx - 1]

        if in_window:
            boundary_runs.append(run)
        else:
            nonboundary_runs.append(run)

        run_rows.append(
            {
                "accession": acc,
                "organism": meta.get("organism", ""),
                "molecule_type": meta.get("molecule_type", ""),
                "sequence_length": seq_len,
                "run_index": idx,
                "base": run.base,
                "run_length": run.length,
                "run_start_1based": run.start1,
                "run_end_1based": run.end1,
                "run_midpoint_1based": run.midpoint1,
                "in_boundary_window": int(in_window),
                "boundary_window": boundary_window,
                "boundary_hits": boundary_hits,
                "nearest_boundary_distance": "" if nearest_dist is None else nearest_dist,
                "nearest_boundary": nearest_desc,
                "overlap_feature_types": ftypes,
                "overlap_feature_labels": flabels,
            }
        )

    genome_row = {
        "accession": acc,
        "description": meta.get("description", ""),
        "organism": meta.get("organism", ""),
        "taxonomy": meta.get("taxonomy", ""),
        "molecule_type": meta.get("molecule_type", ""),
        "topology": meta.get("topology", ""),
        "host": meta.get("host", ""),
        "country": meta.get("country", ""),
        "collection_date": meta.get("collection_date", ""),
        "sequence_length": seq_len,
        "source_length": candidate.source_length,
        "source_gc_frequency": candidate.source_gc_frequency,
        "source_longest_run": candidate.source_longest_run,
        "n_features": len(record.features),
        "n_feature_intervals": len(intervals),
        "n_boundaries": len(boundaries),
        "boundary_window": boundary_window,
        "n_runs": len(runs),
        "n_boundary_runs": len(boundary_runs),
        "n_nonboundary_runs": len(nonboundary_runs),
        "fraction_runs_boundary": len(boundary_runs) / len(runs) if runs else 0.0,
        "processing_seconds": round(time.perf_counter() - t0, 3),
    }
    genome_row.update(summarise_run_subset(boundary_runs, "boundary", thresholds))
    genome_row.update(summarise_run_subset(nonboundary_runs, "nonboundary", thresholds))

    for k in thresholds:
        a = int(genome_row[f"boundary_n_runs_ge_{k}"])
        b = int(genome_row["boundary_n_runs"]) - a
        c = int(genome_row[f"nonboundary_n_runs_ge_{k}"])
        d = int(genome_row["nonboundary_n_runs"]) - c
        genome_row[f"odds_ratio_boundary_ge_{k}"] = odds_ratio(a, b, c, d)

    # Fast boundary summary. For each boundary, only inspect runs whose start/end could be nearby.
    # Because runs are ordered by position, use midpoints for a compact approximation.
    run_midpoints = [r.midpoint1 for r in runs]
    boundary_rows = []
    for b in sorted_boundaries:
        lo = bisect.bisect_left(run_midpoints, b.position1 - boundary_window)
        hi = bisect.bisect_right(run_midpoints, b.position1 + boundary_window)
        nearby = runs[lo:hi]
        # Include any long run spanning into the window but with midpoint outside, rare but possible.
        # This bounded local correction avoids global scans.
        lo2 = max(0, lo - 5)
        hi2 = min(len(runs), hi + 5)
        nearby = [
            r for r in runs[lo2:hi2]
            if interval_distance_to_point(r.start1, r.end1, b.position1) <= boundary_window
        ]
        lens = [r.length for r in nearby]
        boundary_rows.append(
            {
                "accession": acc,
                "organism": meta.get("organism", ""),
                "sequence_length": seq_len,
                "boundary_position_1based": b.position1,
                "boundary_kind": b.boundary_kind,
                "feature_type": b.feature_type,
                "feature_label": b.feature_label,
                "feature_start_1based": b.feature_start1,
                "feature_end_1based": b.feature_end1,
                "feature_strand": b.strand,
                "boundary_window": boundary_window,
                "n_runs_in_window": len(nearby),
                "max_run_in_window": max(lens) if lens else 0,
                "mean_run_in_window": mean(lens),
                "n_runs_ge_5_in_window": sum(1 for x in lens if x >= 5),
                "n_runs_ge_10_in_window": sum(1 for x in lens if x >= 10),
                "n_runs_ge_15_in_window": sum(1 for x in lens if x >= 15),
            }
        )

    return run_rows, genome_row, boundary_rows


def write_csv(rows: Sequence[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0].keys())
    for r in rows[1:]:
        for k in r.keys():
            if k not in fields:
                fields.append(k)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def append_csv_rows(rows: Sequence[dict], path: Path) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists() and path.stat().st_size > 0
    fields = list(rows[0].keys())
    if exists:
        with path.open("r", newline="", encoding="utf-8-sig") as handle:
            reader = csv.reader(handle)
            try:
                fields = next(reader)
            except StopIteration:
                exists = False
    with path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        if not exists:
            writer.writeheader()
        writer.writerows(rows)


def read_existing_csv_rows(path: Path) -> list[dict]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def aggregate_enrichment(run_rows: Sequence[dict], thresholds: Sequence[int]) -> list[dict]:
    out = []
    total_boundary = sum(1 for r in run_rows if int(r["in_boundary_window"]) == 1)
    total_nonboundary = sum(1 for r in run_rows if int(r["in_boundary_window"]) == 0)

    for k in thresholds:
        a = sum(1 for r in run_rows if int(r["in_boundary_window"]) == 1 and int(r["run_length"]) >= k)
        c = sum(1 for r in run_rows if int(r["in_boundary_window"]) == 0 and int(r["run_length"]) >= k)
        b = total_boundary - a
        d = total_nonboundary - c
        out.append(
            {
                "threshold": k,
                "base": "",
                "boundary_runs_ge_threshold": a,
                "boundary_runs_lt_threshold": b,
                "nonboundary_runs_ge_threshold": c,
                "nonboundary_runs_lt_threshold": d,
                "boundary_tail_probability": a / total_boundary if total_boundary else 0.0,
                "nonboundary_tail_probability": c / total_nonboundary if total_nonboundary else 0.0,
                "tail_probability_difference": (a / total_boundary if total_boundary else 0.0)
                - (c / total_nonboundary if total_nonboundary else 0.0),
                "odds_ratio": odds_ratio(a, b, c, d),
            }
        )

    for base in sorted(VALID_BASES):
        rows = [r for r in run_rows if r["base"] == base]
        tb = sum(1 for r in rows if int(r["in_boundary_window"]) == 1)
        tn = sum(1 for r in rows if int(r["in_boundary_window"]) == 0)
        for k in thresholds:
            a = sum(1 for r in rows if int(r["in_boundary_window"]) == 1 and int(r["run_length"]) >= k)
            c = sum(1 for r in rows if int(r["in_boundary_window"]) == 0 and int(r["run_length"]) >= k)
            b = tb - a
            d = tn - c
            out.append(
                {
                    "threshold": k,
                    "base": base,
                    "boundary_runs_ge_threshold": a,
                    "boundary_runs_lt_threshold": b,
                    "nonboundary_runs_ge_threshold": c,
                    "nonboundary_runs_lt_threshold": d,
                    "boundary_tail_probability": a / tb if tb else 0.0,
                    "nonboundary_tail_probability": c / tn if tn else 0.0,
                    "tail_probability_difference": (a / tb if tb else 0.0)
                    - (c / tn if tn else 0.0),
                    "odds_ratio": odds_ratio(a, b, c, d),
                }
            )
    return out


def parse_csv_list(text: str, cast=str) -> list:
    if text is None or text.strip() == "":
        return []
    return [cast(x.strip()) for x in text.split(",") if x.strip()]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Fast boundary-enrichment analysis for homopolymer runs.")
    p.add_argument("--features", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--email", required=True)
    p.add_argument("--api-key", default=None)
    p.add_argument("--delay", type=float, default=0.34)

    p.add_argument("--fetch-timeout", type=int, default=60)
    p.add_argument("--max-retries", type=int, default=2)
    p.add_argument("--cache-dir", type=Path, default=None)
    p.add_argument("--force-refetch", action="store_true")
    p.add_argument("--resume", action="store_true")

    p.add_argument("--id-column", default=None)
    p.add_argument("--length-column", default=None)
    p.add_argument("--longest-run-column", default=None)
    p.add_argument("--gc-column", default=None)
    p.add_argument("--accessions", nargs="*", default=None)
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--max-sequence-length", type=int, default=None, help="Skip records longer than this after fetching.")

    p.add_argument("--boundary-window", type=int, default=250)
    p.add_argument("--boundary-types", default="CDS,gene,mat_peptide")
    p.add_argument("--min-run-lengths", default="5,10,15,20,30")
    p.add_argument("--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    cache_dir = args.cache_dir or (args.outdir / "gb_cache")

    boundary_types = set(parse_csv_list(args.boundary_types, str))
    thresholds = parse_csv_list(args.min_run_lengths, int)

    candidates = read_candidates(
        args.features,
        id_column=args.id_column,
        length_column=args.length_column,
        longest_run_column=args.longest_run_column,
        gc_column=args.gc_column,
        accessions=args.accessions,
        limit=args.limit,
    )

    if args.resume:
        done = read_existing_accessions(args.outdir / "genome_summary.csv")
        before = len(candidates)
        candidates = [c for c in candidates if c.accession not in done]
        print_flush(f"Resume mode: skipping {before - len(candidates)} completed accession(s)")

    print_flush(f"Analysing {len(candidates)} accession(s)")
    print_flush(f"Boundary window: ±{args.boundary_window}")
    print_flush(f"Boundary types: {', '.join(sorted(boundary_types))}")
    print_flush(f"Cache: {cache_dir}")

    for i, cand in enumerate(candidates, start=1):
        total_t0 = time.perf_counter()
        print_flush(f"[{i}/{len(candidates)}] {cand.accession}")

        try:
            fetch_t0 = time.perf_counter()
            record = fetch_genbank_with_timeout(
                cand.accession,
                email=args.email,
                api_key=args.api_key,
                cache_dir=cache_dir,
                timeout=args.fetch_timeout,
                max_retries=args.max_retries,
                delay=args.delay,
                force_refetch=args.force_refetch,
            )
            fetch_seconds = time.perf_counter() - fetch_t0
            seq_len = len(record.seq)

            if args.max_sequence_length is not None and seq_len > args.max_sequence_length:
                raise RuntimeError(
                    f"sequence length {seq_len} exceeds --max-sequence-length {args.max_sequence_length}"
                )

            run_rows, genome_row, boundary_rows = analyse_record(
                candidate=cand,
                record=record,
                boundary_types=boundary_types,
                boundary_window=args.boundary_window,
                thresholds=thresholds,
                verbose=args.verbose,
            )
            genome_row["fetch_seconds"] = round(fetch_seconds, 3)
            genome_row["total_seconds"] = round(time.perf_counter() - total_t0, 3)

            append_csv_rows(run_rows, args.outdir / "run_table.csv")
            append_csv_rows([genome_row], args.outdir / "genome_summary.csv")
            append_csv_rows(boundary_rows, args.outdir / "boundary_summary.csv")

            print_flush(
                f"  ok: length={seq_len}; runs={genome_row['n_runs']}; "
                f"features={genome_row['n_features']}; boundaries={genome_row['n_boundaries']}; "
                f"time={genome_row['total_seconds']}s"
            )

        except Exception as exc:
            row = {"accession": cand.accession, "error": str(exc)}
            append_csv_rows([row], args.outdir / "failures.csv")
            print_flush(f"  SKIPPED: {exc}")

    rows = read_existing_csv_rows(args.outdir / "run_table.csv")
    if rows:
        write_csv(aggregate_enrichment(rows, thresholds), args.outdir / "enrichment_summary.csv")

    print_flush("\nWrote/updated output directory:")
    print_flush(f"  {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
