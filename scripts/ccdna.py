#!/usr/bin/env python3
"""
Core routines for coloured-composition analysis of DNA/RNA sequences.

A sequence is mapped to a valid coloured composition
    (lambda_1^b1, ..., lambda_k^bk)
where lambda_i is a maximal run length and b_i is A/C/G/T/U.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from math import log2, sqrt
from pathlib import Path
import csv
import random
import re
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple

BASES_DNA = set("ACGT")
BASES_RNA = set("ACGU")
IUPAC_DNA = set("ACGTRYSWKMBDHVN")

Run = Tuple[str, int]


def read_fasta(path: str | Path) -> Iterator[Tuple[str, str]]:
    """Yield (header, sequence) records from a FASTA file."""
    header = None
    chunks: List[str] = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).upper()
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks).upper()


def clean_sequence(seq: str, alphabet: str = "DNA", keep_ambiguous: bool = False) -> str:
    """Clean a sequence to A/C/G/T or A/C/G/U.

    If keep_ambiguous=False, ambiguous symbols are removed. This keeps run
    statistics focused on observed canonical bases. For analyses where N-runs
    are meaningful, set keep_ambiguous=True.
    """
    seq = seq.upper().replace(" ", "").replace("\n", "")
    if alphabet.upper() == "RNA":
        allowed = BASES_RNA | ({"N"} if keep_ambiguous else set())
        seq = seq.replace("T", "U")
    else:
        allowed = BASES_DNA | ({"N"} if keep_ambiguous else set())
        seq = seq.replace("U", "T")
    return "".join(ch for ch in seq if ch in allowed)


def coloured_composition(seq: str) -> List[Run]:
    """Return maximal runs as [(base, length), ...]."""
    if not seq:
        return []
    runs: List[Run] = []
    current = seq[0]
    length = 1
    for ch in seq[1:]:
        if ch == current:
            length += 1
        else:
            runs.append((current, length))
            current = ch
            length = 1
    runs.append((current, length))
    return runs


def base_counts(seq: str) -> Dict[str, int]:
    c = Counter(seq)
    return {b: c.get(b, 0) for b in sorted(set(seq) | BASES_DNA)}


def shannon_entropy_from_counts(counts: Sequence[int]) -> float:
    total = sum(counts)
    if total == 0:
        return 0.0
    ent = 0.0
    for n in counts:
        if n:
            p = n / total
            ent -= p * log2(p)
    return ent


def gini(values: Sequence[float]) -> float:
    """Gini coefficient for non-negative values."""
    vals = sorted(v for v in values if v >= 0)
    n = len(vals)
    if n == 0:
        return 0.0
    total = sum(vals)
    if total == 0:
        return 0.0
    weighted = sum((i + 1) * v for i, v in enumerate(vals))
    return (2 * weighted) / (n * total) - (n + 1) / n


def durfee_size(partition: Sequence[int]) -> int:
    """Durfee square size of a partition given as unsorted or sorted parts."""
    parts = sorted(partition, reverse=True)
    d = 0
    for i, part in enumerate(parts, start=1):
        if part >= i:
            d = i
        else:
            break
    return d


def run_length_distribution(runs: Sequence[Run]) -> Counter:
    return Counter((base, length) for base, length in runs)


def transition_counts(runs: Sequence[Run]) -> Counter:
    return Counter((runs[i][0], runs[i + 1][0]) for i in range(len(runs) - 1))


def sequence_id_from_header(header: str) -> str:
    """Extract a compact accession-like identifier from a FASTA header."""
    token = header.split()[0]
    token = token.replace("|", "_")
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", token)


def features_for_sequence(header: str, seq: str, alphabet: str = "DNA") -> Dict[str, object]:
    seq = clean_sequence(seq, alphabet=alphabet)
    runs = coloured_composition(seq)
    lengths = [length for _, length in runs]
    bases = sorted(BASES_RNA if alphabet.upper() == "RNA" else BASES_DNA)
    n = len(seq)
    k = len(runs)
    c = Counter(seq)

    row: Dict[str, object] = {
        "id": sequence_id_from_header(header),
        "description": header,
        "length": n,
        "n_runs": k,
        "mean_run_length": n / k if k else 0.0,
        "longest_run": max(lengths) if lengths else 0,
        "run_length_entropy": shannon_entropy_from_counts(list(Counter(lengths).values())),
        "run_length_gini": gini(lengths),
        "overall_durfee": durfee_size(lengths),
    }

    for b in bases:
        blens = [length for base, length in runs if base == b]
        row[f"count_{b}"] = c.get(b, 0)
        row[f"freq_{b}"] = c.get(b, 0) / n if n else 0.0
        row[f"runs_{b}"] = len(blens)
        row[f"mean_run_{b}"] = sum(blens) / len(blens) if blens else 0.0
        row[f"max_run_{b}"] = max(blens) if blens else 0
        row[f"durfee_{b}"] = durfee_size(blens)
        row[f"entropy_run_lengths_{b}"] = shannon_entropy_from_counts(list(Counter(blens).values())) if blens else 0.0

    # Adjacent run colour transitions, e.g. A->C.
    tc = transition_counts(runs)
    denom = max(k - 1, 1)
    for b1 in bases:
        for b2 in bases:
            if b1 == b2:
                continue
            row[f"transition_{b1}{b2}"] = tc.get((b1, b2), 0) / denom

    return row


def write_features_csv(records: Iterable[Tuple[str, str]], output_csv: str | Path, alphabet: str = "DNA") -> None:
    rows = [features_for_sequence(header, seq, alphabet=alphabet) for header, seq in records]
    if not rows:
        raise ValueError("No FASTA records found.")
    fieldnames = list(rows[0].keys())
    with open(output_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def iid_shuffle(seq: str, rng: random.Random) -> str:
    chars = list(seq)
    rng.shuffle(chars)
    return "".join(chars)


def random_iid_sequence(length: int, probs: Dict[str, float], rng: random.Random) -> str:
    bases = list(probs.keys())
    weights = [probs[b] for b in bases]
    return "".join(rng.choices(bases, weights=weights, k=length))


def js_distance(p: Dict[Tuple[str, int], float], q: Dict[Tuple[str, int], float]) -> float:
    """Jensen-Shannon distance between two discrete distributions."""
    keys = set(p) | set(q)
    def kl(a, b):
        s = 0.0
        for key in keys:
            av = a.get(key, 0.0)
            bv = b.get(key, 0.0)
            if av > 0 and bv > 0:
                s += av * log2(av / bv)
        return s
    m = {key: 0.5 * (p.get(key, 0.0) + q.get(key, 0.0)) for key in keys}
    return sqrt(0.5 * kl(p, m) + 0.5 * kl(q, m))


def empirical_run_distribution(seq: str, alphabet: str = "DNA") -> Dict[Tuple[str, int], float]:
    runs = coloured_composition(clean_sequence(seq, alphabet=alphabet))
    counts = run_length_distribution(runs)
    total = sum(counts.values())
    return {key: value / total for key, value in counts.items()} if total else {}
