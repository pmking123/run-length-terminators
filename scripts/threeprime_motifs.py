#!/usr/bin/env python3
"""
threeprime_motifs.py

Mechanism check for the long-A/T-run enrichment at viral gene 3' ends.

For every CDS in the GenBank cache it:
  1. builds a single-base-resolution composition profile around the 5' and 3'
     ends, oriented in the CODING direction so that positive offset = OUTSIDE
     the gene (downstream of the stop for 3' ends, upstream of the start for 5');
  2. tallies recognisable gene-end motifs in the 3' downstream flank
     (poly-A tract, poly-T/U tract, canonical poly-A signal AATAAA/ATTAAA);
  3. compares those motif rates against a matched set of CDS-interior windows.

If the long A/T runs at 3' ends are termination / polyadenylation signals, the
3' profile will show an A/T spike just downstream of the stop and the motif
rates will be enriched several-fold over the interior control.

Input  : --gb-cache  directory of *.gb GenBank files
           (e.g. scripts/analysis_initial/boundary_enrichment/gb_cache)
Outputs: <outdir>/threeprime_positional_profile.csv
         <outdir>/threeprime_motif_summary.csv
         <outdir>/threeprime_flanks.fasta
         <outdir>/threeprime_profile.png
         <outdir>/threeprime_motif_enrichment.png

Needs: biopython, numpy, pandas, matplotlib.  No network access.
"""
from __future__ import annotations
import argparse, random, re
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO

COMP = str.maketrans("ACGTacgtNn", "TGCATGCANN")
BASES = ["A", "C", "G", "T"]


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def cds_features(record):
    """Yield (terminal5_idx0, terminal3_idx0, strand) for each CDS, 0-based."""
    for feat in record.features:
        if feat.type != "CDS":
            continue
        strand = feat.location.strand
        if strand not in (1, -1):
            continue
        start = int(feat.location.start)        # 0-based, inclusive
        end = int(feat.location.end) - 1        # 0-based, inclusive
        if end <= start:
            continue
        if strand == 1:
            yield start, end, 1                 # 5' = start, 3' = end
        else:
            yield end, start, -1                # 5' = end, 3' = start (genomic)


def flank_bases(seq: str, terminal0: int, strand: int, W: int):
    """
    Return list of bases for offsets -W..+W in CODING direction.
    offset 0 = terminal coding base; positive offset = OUTSIDE the gene flank.
    Missing positions (off the genome) are returned as ''.
    """
    L = len(seq)
    out = []
    for o in range(-W, W + 1):
        if strand == 1:
            g = terminal0 + o      # outside (o>0) -> higher coord (downstream)
            b = seq[g] if 0 <= g < L else ""
        else:
            g = terminal0 - o      # outside (o>0) -> lower coord; complement
            b = seq[g].translate(COMP) if 0 <= g < L else ""
        out.append(b.upper())
    return out


def longest_run(s: str, base: str) -> int:
    best = cur = 0
    for ch in s:
        if ch == base:
            cur += 1; best = max(best, cur)
        else:
            cur = 0
    return best


def motif_flags(flank: str, min_run: int):
    """flank = coding-strand bases immediately downstream of the stop (3' UTR side)."""
    return {
        f"polyA_run_ge{min_run}": int(longest_run(flank, "A") >= min_run),
        f"polyT_run_ge{min_run}": int(longest_run(flank, "T") >= min_run),
        "polyA_signal_AATAAA": int("AATAAA" in flank or "ATTAAA" in flank),
        "AT_frac_ge_0.7": int(sum(flank.count(b) for b in "AT") >= 0.7 * len(flank)) if flank else 0,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gb-cache", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--window", type=int, default=50)
    ap.add_argument("--min-run", type=int, default=5)
    ap.add_argument("--control-per-genome", type=int, default=60)
    ap.add_argument("--min-control-dist", type=int, default=200)
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    W = args.window

    gb_files = sorted(p for p in args.gb_cache.glob("*.gb") if not p.name.endswith(".tmp.gb"))
    if not gb_files:
        raise SystemExit(f"No .gb files in {args.gb_cache}")

    # offset -> base counts, for 5' and 3'
    prof = {end: {b: np.zeros(2 * W + 1) for b in BASES} for end in ("5prime", "3prime")}
    prof_n = {end: np.zeros(2 * W + 1) for end in ("5prime", "3prime")}
    motif_3p = []      # list of dicts (one per CDS 3' flank)
    motif_ctrl = []    # control windows
    flank_records = []

    n_cds = 0
    for gb in gb_files:
        try:
            recs = list(SeqIO.parse(str(gb), "genbank"))
        except Exception:
            continue
        for record in recs:
            seq = str(record.seq).upper()
            L = len(seq)
            ends = list(cds_features(record))
            if not ends:
                continue
            term3_positions = []
            for t5, t3, strand in ends:
                n_cds += 1
                # ---- positional profiles ----
                for end_name, term in (("5prime", t5), ("3prime", t3)):
                    fb = flank_bases(seq, term, strand, W)
                    for i, b in enumerate(fb):
                        if b in prof[end_name]:
                            prof[end_name][b][i] += 1
                            prof_n[end_name][i] += 1
                # ---- 3' downstream flank motif tally ----
                fb3 = flank_bases(seq, t3, strand, W)
                flank = "".join(fb3[W + 1:])           # offsets +1..+W = outside (3' UTR side)
                if flank:
                    motif_3p.append(motif_flags(flank, args.min_run))
                    acc = record.id
                    flank_records.append(f">{acc}_CDS3p_strand{strand}\n{flank}")
                term3_positions.append(t3)
            # ---- interior controls for this genome ----
            term3_positions = np.array(term3_positions)
            tries = 0
            got = 0
            while got < args.control_per_genome and tries < args.control_per_genome * 20:
                tries += 1
                p = rng.randint(W, L - W - 1)
                if len(term3_positions) and np.min(np.abs(term3_positions - p)) < args.min_control_dist:
                    continue
                strand = rng.choice([1, -1])
                fb = flank_bases(seq, p, strand, W)
                flank = "".join(fb[W + 1:])
                if flank:
                    motif_ctrl.append(motif_flags(flank, args.min_run))
                    got += 1

    # ----- positional profile table + plot -----
    rows = []
    for end_name in ("5prime", "3prime"):
        for i in range(2 * W + 1):
            offset = i - W
            n = prof_n[end_name][i]
            row = {"end": end_name, "offset": offset, "n": int(n)}
            for b in BASES:
                row[f"frac_{b}"] = prof[end_name][b][i] / n if n else np.nan
            rows.append(row)
    prof_df = pd.DataFrame(rows)
    prof_df.to_csv(args.outdir / "threeprime_positional_profile.csv", index=False)

    # ----- motif summary -----
    m3 = pd.DataFrame(motif_3p); mc = pd.DataFrame(motif_ctrl)
    summ = []
    for col in m3.columns:
        r3 = m3[col].mean(); rc = mc[col].mean() if col in mc else np.nan
        # Haldane-style continuity correction keeps the ratio finite at rc=0
        eps3 = 0.5 / max(len(m3), 1); epsc = 0.5 / max(len(mc), 1)
        enr = (r3 + eps3) / (rc + epsc) if np.isfinite(rc) else np.nan
        summ.append(dict(motif=col, threeprime_rate=r3, control_rate=rc,
                         enrichment=enr, n_3prime=len(m3), n_control=len(mc)))
    summ_df = pd.DataFrame(summ)
    summ_df.to_csv(args.outdir / "threeprime_motif_summary.csv", index=False)
    (args.outdir / "threeprime_flanks.fasta").write_text("\n".join(flank_records) + "\n")

    # ----- plots -----
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.3), sharey=True)
    colors = {"A": "#2980b9", "T": "#c0392b", "G": "#27ae60", "C": "#e67e22"}
    for j, end_name in enumerate(("5prime", "3prime")):
        sub = prof_df[prof_df["end"] == end_name]
        for b in BASES:
            ax[j].plot(sub["offset"], sub[f"frac_{b}"], color=colors[b], label=b, lw=1.6)
        ax[j].axvline(0, ls="--", color="#555", lw=1)
        ax[j].set_title(f"{end_name} end  (offset>0 = outside gene)", fontsize=10)
        ax[j].set_xlabel("offset from terminal coding base (bp)")
        ax[j].grid(alpha=.3)
    ax[0].set_ylabel("base fraction"); ax[0].legend(ncol=4, fontsize=8)
    plt.tight_layout(); plt.savefig(args.outdir / "threeprime_profile.png", dpi=130); plt.close()

    sd = summ_df.dropna(subset=["enrichment"])
    fig, ax = plt.subplots(figsize=(7, 4))
    y = np.arange(len(sd))[::-1]
    ax.barh(y, sd["enrichment"], color="#8e44ad")
    ax.axvline(1.0, ls="--", color="#555")
    ax.set_yticks(y); ax.set_yticklabels(sd["motif"], fontsize=8)
    ax.set_xlabel("3' flank rate / interior-control rate")
    ax.set_title("Gene-end motif enrichment in 3' flanks")
    plt.tight_layout(); plt.savefig(args.outdir / "threeprime_motif_enrichment.png", dpi=130); plt.close()

    print(f"CDS analysed: {n_cds:,}   3' flanks: {len(m3):,}   control windows: {len(mc):,}")
    print("\n=== 3' downstream-flank motif rates vs interior control ===")
    print(summ_df.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(f"\nOutputs in {args.outdir}/")


if __name__ == "__main__":
    main()
