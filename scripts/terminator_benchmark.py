#!/usr/bin/env python3
"""
terminator_benchmark.py

Experiment 1 + 2: does an annotation-free run-length signal flag real intrinsic
terminators, including ones a dedicated predictor (TransTermHP) misses?

The caller is deliberately trivial: scan the top strand for maximal homopolymer
runs of T (forward-oriented U-tract) or A (reverse-oriented U-tract) of length
>= k. That single threshold k is the only knob. Each candidate is then
corroborated two independent ways, neither of which the caller used:

  * STRUCTURAL: fold the window immediately 5' of the U-tract (in RNA sense) and
    test for a hairpin (ViennaRNA MFE if available; inverted-repeat score
    otherwise). Intrinsic terminators have a GC-rich stem-loop there.
  * POSITIONAL: is the U-tract just downstream of a same-orientation CDS 3' end?

If a TransTermHP report (or a generic BED of terminator intervals) is supplied,
the caller's candidates are overlapped with it, giving the four quadrants
(both / run-only / predictor-only). The headline number is how many RUN-ONLY
candidates are corroborated -- i.e. plausible terminators the predictor missed.

Inputs
------
--gb-cache DIR     GenBank files (sequence + CDS features). Required.
--transterm FILE   TransTermHP .tt/.bag output (optional).
--reference-bed F  Generic terminator intervals: TSV 'accession start end strand'
                   (optional; alternative to --transterm).
--outdir DIR

Key options
-----------
--k-tract 6        U-tract length threshold (the one knob; also sweeps around it)
--hairpin-window 30   bp folded upstream of the tract
--mfe-threshold -3.0  kcal/mol; hairpin called if MFE <= this
--gene-end-window 50  bp; U-tract within this of a same-strand CDS 3' end = proximal
--overlap-tol 12      bp padding for candidate<->reference overlap

Outputs
-------
terminator_candidates.csv        every run candidate + corroboration flags
terminator_benchmark_summary.csv overlap quadrants + corroboration rates
terminator_sensitivity_sweep.csv k vs recall / run-only counts
terminator_benchmark_report.md

Needs: biopython, numpy, pandas. ViennaRNA optional (falls back to inverted-repeat).
"""
from __future__ import annotations
import argparse, re, shutil, subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO

COMP = str.maketrans("ACGTacgtNn", "TGCATGCANN")
try:
    import RNA  # ViennaRNA python bindings
    _HAVE_RNA = True
except Exception:
    _HAVE_RNA = False


class _RNAfoldProc:
    """Persistent RNAfold subprocess: fold windows via the module's binary when
    the Python bindings aren't importable (common on HPC module setups)."""
    def __init__(self, exe):
        self.p = subprocess.Popen([exe, "--noPS"], stdin=subprocess.PIPE,
                                  stdout=subprocess.PIPE, text=True, bufsize=1)
        self._e = re.compile(r"\(\s*(-?\d+\.\d+)\)\s*$")

    def fold(self, seq):
        s = seq.replace("T", "U").replace("t", "u").upper()
        if len(s) < 4:
            return 0.0, ""
        self.p.stdin.write(s + "\n"); self.p.stdin.flush()
        self.p.stdout.readline()                       # echoed sequence
        line = self.p.stdout.readline().strip()
        if not line:
            return 0.0, ""
        struct = line.split()[0]
        m = self._e.search(line)
        return (float(m.group(1)) if m else 0.0), struct


_FOLDER = None
if not _HAVE_RNA:
    _exe = shutil.which("RNAfold")
    if _exe:
        try:
            _FOLDER = _RNAfoldProc(_exe)
        except Exception:
            _FOLDER = None


def _fold_mode():
    if _HAVE_RNA:
        return "ViennaRNA python bindings"
    if _FOLDER is not None:
        return "RNAfold binary (module)"
    return "OFF (inverted-repeat fallback)"


def revcomp(s): return s.translate(COMP)[::-1]


# ---------------------------------------------------------------- structure
def best_inverted_repeat(seq):
    """Fold-free fallback: longest near-perfect stem of a hairpin in `seq`.
    Returns (stem_len, gc_frac_of_stem). Crude, used only without ViennaRNA."""
    n = len(seq); best = 0; best_gc = 0.0
    pair = {"A": "T", "T": "A", "G": "C", "C": "G"}
    for c in range(1, n - 1):
        for loop in range(3, 11):
            li, ri = c - 1, c + loop
            stem = 0; gc = 0
            while li >= 0 and ri < n and pair.get(seq[li]) == seq[ri]:
                stem += 1
                if seq[li] in "GC": gc += 1
                li -= 1; ri += 1
            if stem > best:
                best = stem; best_gc = gc / stem if stem else 0.0
    return best, best_gc


def hairpin_score(window):
    """Return (mfe, has_hairpin, gc_frac). Uses ViennaRNA (bindings or binary) if present."""
    gc = (window.count("G") + window.count("C")) / len(window) if window else 0.0
    if _HAVE_RNA and len(window) >= 8:
        struct, mfe = RNA.fold(window.replace("T", "U"))
        return mfe, bool((mfe <= HAIRPIN_MFE) and ("((" in struct)), gc
    if _FOLDER is not None and len(window) >= 8:
        mfe, struct = _FOLDER.fold(window)
        return mfe, bool((mfe <= HAIRPIN_MFE) and ("((" in struct)), gc
    stem, stem_gc = best_inverted_repeat(window)
    pseudo = -(stem * (1.0 + stem_gc))
    return pseudo, bool(stem >= 4), gc


# ---------------------------------------------------------------- genome IO
def load_genome(record):
    seq = str(record.seq).upper()
    cds3p = []   # (three_prime_1based, strand)
    for feat in record.features:
        if feat.type != "CDS":
            continue
        st = feat.location.strand
        if st not in (1, -1):
            continue
        s0 = int(feat.location.start); e0 = int(feat.location.end) - 1
        if e0 <= s0:
            continue
        cds3p.append((e0 + 1 if st == 1 else s0 + 1, st))
    return seq, cds3p


def iter_runs(seq, base, k):
    """Yield (start1, end1, length) for maximal runs of `base` with length>=k."""
    n = len(seq); i = 0
    while i < n:
        if seq[i] == base:
            j = i
            while j < n and seq[j] == base:
                j += 1
            if j - i >= k:
                yield (i + 1, j, j - i)
            i = j
        else:
            i += 1


def _gene_end_prox(s1, e1, base, plus3p, minus3p, gew, k, mode):
    """Is the tract just downstream of a CDS 3' end, under the chosen strand model?
    plus_down  = tract lies within the window 3' (higher coord) of a + CDS end.
    minus_down = tract lies within the window 3' (lower coord)  of a - CDS end.
      same  (intrinsic terminator): T<->plus_down,  A<->minus_down
      polyA (Mononegavirales gene-end / poly-A signal): A<->plus_down, T<->minus_down
      any   : either, regardless of base
    """
    win = gew + k
    plus_down = bool(len(plus3p) and ((s1 - plus3p >= 0) & (s1 - plus3p <= win)).any())
    minus_down = bool(len(minus3p) and ((minus3p - e1 >= 0) & (minus3p - e1 <= win)).any())
    if mode == "any":
        return plus_down or minus_down
    if mode == "polyA":
        return plus_down if base == "A" else minus_down
    return plus_down if base == "T" else minus_down          # "same" (default)


def call_candidates(seq, cds3p, k, H, gew, mode="same"):
    plus3p  = np.array(sorted(p for p, s in cds3p if s == 1), dtype=np.int64)
    minus3p = np.array(sorted(p for p, s in cds3p if s == -1), dtype=np.int64)
    L = len(seq); rows = []
    for base, orient in (("T", "+"), ("A", "-")):
        for s1, e1, ln in iter_runs(seq, base, k):
            if orient == "+":                         # hairpin upstream (lower coord)
                w = seq[max(0, s1 - 1 - H): s1 - 1]
            else:                                     # A-run: hairpin downstream, revcomp
                w = revcomp(seq[e1: min(L, e1 + H)])
            prox = _gene_end_prox(s1, e1, base, plus3p, minus3p, gew, ln, mode)
            mfe, hp, gc = hairpin_score(w)
            rows.append(dict(orient=orient, tract_start=s1, tract_end=e1,
                             tract_len=ln, hairpin_mfe=mfe, has_hairpin=hp,
                             hairpin_gc=gc, gene_end_proximal=prox))
    cols = ["orient", "tract_start", "tract_end", "tract_len", "hairpin_mfe",
            "has_hairpin", "hairpin_gc", "gene_end_proximal"]
    return pd.DataFrame(rows, columns=cols)


def null_corroboration(seq, cds3p, n, H, gew, k, rng, mode="same"):
    """Background corroboration rate under the chosen gene-end mode: apply the
    identical hairpin + gene-end tests at n random positions with a random tract
    base, so run-only rates can be judged against chance. Returns
    (hairpin_rate, proximal_rate, corroborated_rate, n)."""
    plus3p  = np.array(sorted(p for p, s in cds3p if s == 1), dtype=np.int64)
    minus3p = np.array(sorted(p for p, s in cds3p if s == -1), dtype=np.int64)
    L = len(seq)
    if L <= 2 * H + 2 or n <= 0:
        return (np.nan, np.nan, np.nan, 0)
    hp = pr = corr = 0
    for _ in range(n):
        p = int(rng.integers(H + 1, L - H))
        base = "T" if rng.random() < 0.5 else "A"
        if base == "T":
            w = seq[p - 1 - H: p - 1]
        else:
            w = revcomp(seq[p: p + H])
        prox = _gene_end_prox(p, p, base, plus3p, minus3p, gew, k, mode)
        _, has, _ = hairpin_score(w)
        hp += int(has); pr += int(prox); corr += int(has or prox)
    return (hp / n, pr / n, corr / n, n)


# ---------------------------------------------------------------- reference
def parse_transterm(path):
    """Parse TransTermHP terminator lines -> list of (start,end,strand)."""
    rx = re.compile(r"TERM\s+\d+\s+(\d+)\s+-\s+(\d+)\s+([+-])")
    out = []
    cur = None
    for line in Path(path).read_text(errors="ignore").splitlines():
        m = re.match(r"SEQUENCE\s+(\S+)", line.strip())
        if m:
            cur = m.group(1); continue
        mm = rx.search(line)
        if mm:
            a, b, st = int(mm.group(1)), int(mm.group(2)), mm.group(3)
            out.append((cur, min(a, b), max(a, b), st))
    return out


def parse_bed(path):
    df = pd.read_csv(path, sep=r"\s+", header=None,
                     names=["accession", "start", "end", "strand"])
    return [(r.accession, int(min(r.start, r.end)), int(max(r.start, r.end)),
             str(r.strand)) for r in df.itertuples()]


def overlap(cand_df, refs, tol):
    """Mark each candidate matched if within tol of a reference interval; return
    matched-candidate mask and the count of reference intervals recovered."""
    matched = np.zeros(len(cand_df), dtype=bool)
    ref_hit = np.zeros(len(refs), dtype=bool)
    cs = cand_df["tract_start"].to_numpy(); ce = cand_df["tract_end"].to_numpy()
    for ri, (_, rs, re_, _st) in enumerate(refs):
        ov = (cs - tol <= re_) & (ce + tol >= rs)
        if ov.any():
            matched |= ov; ref_hit[ri] = True
    return matched, ref_hit


HAIRPIN_MFE = -3.0  # set in main from args


def main():
    global HAIRPIN_MFE
    ap = argparse.ArgumentParser()
    ap.add_argument("--gb-cache", required=True, type=Path)
    ap.add_argument("--transterm", type=Path, default=None)
    ap.add_argument("--reference-bed", type=Path, default=None)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--k-tract", type=int, default=6)
    ap.add_argument("--hairpin-window", type=int, default=30)
    ap.add_argument("--mfe-threshold", type=float, default=-3.0)
    ap.add_argument("--gene-end-window", type=int, default=50)
    ap.add_argument("--gene-end-mode", choices=["same", "polyA", "any"], default="same",
                    help="strand model for gene-end proximity: same=intrinsic terminator (phages); polyA=Mononegavirales gene-end/poly-A signal; any=either strand")
    ap.add_argument("--overlap-tol", type=int, default=12)
    ap.add_argument("--null-per-genome", type=int, default=300,
                    help="random positions per genome for the background "
                         "corroboration rate (0 disables the null)")
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--sweep", default="4,5,6,7,8,10")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    HAIRPIN_MFE = args.mfe_threshold

    gbs = sorted(p for p in args.gb_cache.glob("*.gb") if not p.name.endswith(".tmp.gb"))
    if not gbs:
        raise SystemExit(f"No .gb files in {args.gb_cache}")

    # reference terminators grouped by accession
    refs_all = []
    if args.transterm:
        refs_all = parse_transterm(args.transterm)
    elif args.reference_bed:
        refs_all = parse_bed(args.reference_bed)
    refs_by_acc = {}
    for acc, s, e, st in refs_all:
        refs_by_acc.setdefault(acc, []).append((acc, s, e, st))

    print(f"RNA folding: {_fold_mode()}")
    all_cand = []; sweep_rows = []
    null_rng = np.random.default_rng(args.seed)
    null_hp = []; null_pr = []; null_corr = []   # per-genome rates, weighted later
    null_w = []
    ks = sorted({args.k_tract, *[int(x) for x in args.sweep.split(",")]})

    # main pass at the chosen k, full corroboration + overlap
    n_ref_total = 0; n_ref_recovered = 0
    for gb in gbs:
        for record in SeqIO.parse(str(gb), "genbank"):
            seq, cds3p = load_genome(record)
            if not seq:
                continue
            cand = call_candidates(seq, cds3p, args.k_tract,
                                   args.hairpin_window, args.gene_end_window, args.gene_end_mode)
            cand.insert(0, "accession", record.id)
            refs = refs_by_acc.get(record.id, [])
            if refs:
                matched, ref_hit = overlap(cand, refs, args.overlap_tol)
                cand["matches_reference"] = matched
                n_ref_total += len(refs); n_ref_recovered += int(ref_hit.sum())
            else:
                cand["matches_reference"] = pd.NA
            all_cand.append(cand)
            if args.null_per_genome > 0:
                hpr, prr, cor, nn = null_corroboration(
                    seq, cds3p, args.null_per_genome, args.hairpin_window,
                    args.gene_end_window, args.k_tract, null_rng, args.gene_end_mode)
                if nn:
                    null_hp.append(hpr); null_pr.append(prr); null_corr.append(cor); null_w.append(nn)
    C = pd.concat(all_cand, ignore_index=True) if all_cand else pd.DataFrame()
    C["corroborated"] = C["has_hairpin"] | C["gene_end_proximal"]
    C.to_csv(args.outdir / "terminator_candidates.csv", index=False)

    # sensitivity sweep over k (recall of reference + run-only volume)
    if refs_all:
        for k in ks:
            tot = rec = run_only = 0
            for gb in gbs:
                for record in SeqIO.parse(str(gb), "genbank"):
                    seq, cds3p = load_genome(record)
                    if not seq: continue
                    cand = call_candidates(seq, cds3p, k, args.hairpin_window, args.gene_end_window, args.gene_end_mode)
                    refs = refs_by_acc.get(record.id, [])
                    if not refs:
                        continue
                    matched, ref_hit = overlap(cand, refs, args.overlap_tol)
                    tot += len(refs); rec += int(ref_hit.sum())
                    run_only += int((~matched).sum())
            sweep_rows.append(dict(k_tract=k, ref_total=tot, ref_recovered=rec,
                                   recall=rec / tot if tot else np.nan,
                                   run_only_candidates=run_only))
        sweep = pd.DataFrame(sweep_rows)
        sweep.to_csv(args.outdir / "terminator_sensitivity_sweep.csv", index=False)

    # summary
    summ = {}
    summ["n_candidates"] = len(C)
    summ["frac_with_hairpin"] = float(C["has_hairpin"].mean()) if len(C) else np.nan
    summ["frac_gene_end_proximal"] = float(C["gene_end_proximal"].mean()) if len(C) else np.nan
    summ["frac_corroborated"] = float(C["corroborated"].mean()) if len(C) else np.nan
    if refs_all:
        mref = C["matches_reference"]
        ro = C[mref.eq(False).fillna(False)]
        both = C[mref.eq(True).fillna(False)]
        summ["reference_terminators"] = n_ref_total
        summ["reference_recovered_by_runs"] = n_ref_recovered
        summ["reference_recall"] = n_ref_recovered / n_ref_total if n_ref_total else np.nan
        summ["candidates_matching_reference"] = int(len(both))
        summ["run_only_candidates"] = int(len(ro))
        summ["run_only_with_hairpin"] = float(ro["has_hairpin"].mean()) if len(ro) else np.nan
        summ["run_only_gene_end_proximal"] = float(ro["gene_end_proximal"].mean()) if len(ro) else np.nan
        summ["run_only_corroborated"] = float(ro["corroborated"].mean()) if len(ro) else np.nan
        summ["run_only_corroborated_count"] = int(ro["corroborated"].sum()) if len(ro) else 0
        # ---- matched null: background corroboration rate + enrichment ----
        if null_w:
            w = np.array(null_w, dtype=float); wsum = w.sum()
            nh = float(np.nansum(np.array(null_hp) * w) / wsum)
            npx = float(np.nansum(np.array(null_pr) * w) / wsum)
            ncr = float(np.nansum(np.array(null_corr) * w) / wsum)
            summ["null_hairpin_rate"] = nh
            summ["null_gene_end_proximal_rate"] = npx
            summ["null_corroborated_rate"] = ncr
            ro_hp = float(ro["has_hairpin"].mean()) if len(ro) else np.nan
            ro_px = float(ro["gene_end_proximal"].mean()) if len(ro) else np.nan
            ro_cr = float(ro["corroborated"].mean()) if len(ro) else np.nan
            summ["enrich_run_only_hairpin"] = ro_hp / nh if nh else np.nan
            summ["enrich_run_only_gene_end"] = ro_px / npx if npx else np.nan
            summ["enrich_run_only_corroborated"] = ro_cr / ncr if ncr else np.nan
            # excess count = run-only corroborated beyond background expectation
            summ["run_only_corroborated_excess"] = int(round(len(ro) * (ro_cr - ncr))) if len(ro) else 0
    pd.DataFrame([summ]).to_csv(args.outdir / "terminator_benchmark_summary.csv", index=False)

    # report
    lines = ["# Terminator-benchmark report\n",
             f"- RNA folding: **{_fold_mode()}**",
             f"- gene-end mode: {args.gene_end_mode}\n"
             f"- U-tract threshold k = {args.k_tract}; hairpin window {args.hairpin_window} bp; "
             f"MFE<= {args.mfe_threshold}; gene-end window {args.gene_end_window} bp\n",
             "## Summary\n"]
    for k, v in summ.items():
        lines.append(f"- {k}: {v:.4f}" if isinstance(v, float) else f"- {k}: {v}")
    if refs_all:
        lines.append("\n## Interpretation\n")
        lines.append("`reference_recall` = fraction of TransTermHP terminators the run-length "
                     "caller also flags. `run_only_corroborated_count` = run-only candidates with "
                     "an independent hairpin and/or a gene-end position -- plausible terminators the "
                     "predictor missed. High run-only-corroborated => the run signal adds real calls; "
                     "near-zero => it mostly re-finds (or over-calls) what the predictor already has.")
    (args.outdir / "terminator_benchmark_report.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
