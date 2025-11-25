#!/usr/bin/env python3
"""
Reimplementation of legacy primer.pl in Python.

Key features (mirrors original intent):
- Read FASTA with one or more sequences.
- Non-tiled mode: find up to N primer pairs near each end of the sequence.
- Tiled mode: walk across the sequence producing overlapping amplicons.
- Primer filters: length, GC%, Tm (64.9 + 41*(G+C-16.4)/N), homopolymer runs,
  3'-end clamp, simple hairpin (palindromic ends), simple dimer check vs partner.
- Optional BLAST specificity check with e-value threshold.
- Outputs a report (.primers.txt) and the amplicon fragments (.fragments.fasta).

NOTE: The original Perl file is partially truncated, so a few behaviors (esp. tiling
      details and BLAST DB path) are inferred and modernized conservatively.
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
from collections import OrderedDict
from typing import Dict, Tuple, Optional

# ---------------------- FASTA utilities ----------------------

def read_fasta(path: str) -> Dict[str, str]:
    seqs: Dict[str, str] = OrderedDict()
    name: Optional[str] = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].strip()
                # de-dup headers if repeated
                if name in seqs:
                    idx = 2
                    base = name
                    while f"{base}_{idx}" in seqs:
                        idx += 1
                    name = f"{base}_{idx}"
                seqs[name] = []
            else:
                if name is None:
                    raise ValueError("FASTA missing header before sequence lines")
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def revcomp(s: str) -> str:
    table = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return s.translate(table)[::-1]

# ---------------------- Primer metrics ----------------------

def tm(seq: str) -> float:
    """Melting temp (Perl formula): 64.9 + 41*((G+C-16.4)/N)."""
    s = seq.upper()
    A, T, G, C = s.count("A"), s.count("T"), s.count("G"), s.count("C")
    n = A + T + G + C
    if n == 0:
        return 0.0
    return round(64.9 + 41.0 * ((G + C - 16.4) / n), 2)


def gc_percent(seq: str) -> float:
    s = seq.upper()
    n = len(s)
    if n == 0:
        return 0.0
    return 100.0 * (s.count("G") + s.count("C")) / n


def has_long_homopolymer(seq: str, max_run: int) -> bool:
    """Return True if there's a homopolymer run longer than max_run."""
    s = seq.upper()
    patt = re.compile(rf"(A{{{max_run+1},}}|T{{{max_run+1},}}|G{{{max_run+1},}}|C{{{max_run+1},}})")
    return bool(patt.search(s))


def palindrome_ends(seq: str, k: int = 4) -> bool:
    """Simple hairpin proxy: first k vs reverse-complement of last k."""
    s = seq.upper()
    if len(s) < 2 * k:
        return False
    return s[:k] == revcomp(s[-k:])


def bad_3prime(seq: str) -> bool:
    """Penalize weak clamps and long G/C clamps at 3' end.
    - Disallow ending with GGG or CCC.
    - Disallow last two bases both A/T (weak clamp).
    """
    s = seq.upper()
    if len(s) >= 3 and (s.endswith("GGG") or s.endswith("CCC")):
        return True
    if len(s) >= 2 and set(s[-2:]).issubset({"A", "T"}):
        return True
    return False


def dimer_with(primer: str, other: Optional[str], k: int = 3) -> bool:
    if not other:
        return False
    p3 = primer[-k:].upper()
    o3 = revcomp(other)[-k:].upper()
    return p3 == o3

# ---------------------- Selection logic ----------------------

def primer_is_valid(primer: str, args: argparse.Namespace, partner: Optional[str] = None) -> bool:
    if args.require_upper and primer != primer.upper():
        return False
    if "N" in primer.upper():
        return False
    L = len(primer)
    if L < args.min_len or L > args.max_len:
        return False
    if palindrome_ends(primer, 4):
        return False
    if has_long_homopolymer(primer, args.max_homopoly):
        return False
    if bad_3prime(primer):
        return False
    gc = gc_percent(primer)
    if not (args.min_gc <= gc <= args.max_gc):
        return False
    T = tm(primer)
    if not (args.min_tm <= T <= args.max_tm):
        return False
    if dimer_with(primer, partner, 3):
        return False
    return True


def choose_primer(window: str, full_seq: str, args: argparse.Namespace, partner: Optional[str] = None) -> Optional[str]:
    """Scan a window for the first primer passing all filters.
    The original Perl loop incremented a base offset and tried lengths from min..max.
    We emulate that and cap the offset search ~25 bases like the legacy script.
    """
    max_offset = min(25, max(0, len(window) - args.min_len))
    for offset in range(args.start_offset, max_offset + 1):
        max_len_here = min(args.max_len, len(window) - offset)
        for L in range(args.min_len, max_len_here + 1):
            primer = window[offset : offset + L]
            if not primer_is_valid(primer, args, partner):
                continue
            if args.skip_blast or not args.blast_db:
                return primer
            if primer_passes_blast(primer, args):
                return primer
    return None


def primer_passes_blast(primer: str, args: argparse.Namespace) -> bool:
    """Return True if BLAST finds ≤ 1 hit (self only). Uses blastn-short if available.
    If BLAST fails or isn't installed, we WARN and return True (do not block).
    """
    try:
        with tempfile.TemporaryDirectory() as td:
            q = os.path.join(td, "q.fa")
            with open(q, "w") as fh:
                fh.write(">primer\n" + primer + "\n")
            cmd = [
                "blastn",
                "-task",
                "blastn-short",
                "-db",
                args.blast_db,
                "-query",
                q,
                "-evalue",
                str(args.evalue),
                "-outfmt",
                "6 qseqid sseqid pident length evalue bitscore",
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode != 0:
                sys.stderr.write(f"[warn] blastn failed ({res.returncode}): {res.stderr.strip()}\n")
                return True
            hits = [ln for ln in res.stdout.splitlines() if ln.strip()]
            return len(hits) <= 1
    except FileNotFoundError:
        sys.stderr.write("[warn] blastn not found; skipping BLAST specificity check.\n")
        return True

# ---------------------- Workflows ----------------------

def non_tiled(sequence: str, name: str, args: argparse.Namespace, out_fh, path_fh) -> None:
    """Find up to N primer pairs near the ends, moving the windows in 30-bp steps."""
    success = 0
    step = 30
    fwd = sequence
    rev = revcomp(sequence)
    start = 0
    rstart = 0
    while success < args.num_pairs and start <= args.f_start and rstart <= args.r_start:
        fwin = fwd[start : start + args.f_start]
        rwin = rev[rstart : rstart + args.r_start]
        fprimer = choose_primer(fwin, fwd, args, partner=None)
        if not fprimer:
            start += step
            continue
        rprimer = choose_primer(rwin, rev, args, partner=fprimer)
        if not rprimer:
            rstart += step
            continue
        # map reverse primer to forward coordinates
        rev_rprimer = revcomp(rprimer)
        fpos = fwd.find(fprimer, start)
        rpos_end = fwd.rfind(rev_rprimer) + len(rev_rprimer)
        if fpos == -1 or rpos_end <= 0:
            # fallback conservative mapping
            fpos = fwd.find(fprimer)
            rpos_end = fwd.find(rev_rprimer, fpos + len(fprimer)) + len(rev_rprimer)
            if fpos == -1 or rpos_end <= 0:
                start += step
                rstart += step
                continue
        frag = fwd[fpos:rpos_end]
        global_f = args.global_offset + fpos + 1  # 1-based inclusive
        global_r = args.global_offset + rpos_end   # 1-based inclusive
        out_fh.write(f"{name}_{global_f}-{global_r} Exp.Frag.size {len(frag)} bp of sequence {len(sequence)} bp\n")
        out_fh.write(f"Forward Primer:\t{fprimer}\tTm={tm(fprimer)}\n")
        out_fh.write(f"Reverse Primer:\t{rprimer}\tTm={tm(rprimer)}\n\n")
        path_fh.write(f">{name}_{global_f}-{global_r}\n{frag}\n\n")
        success += 1
        start += step
        rstart += step


def tiled(sequence: str, name: str, args: argparse.Namespace, out_fh, path_fh) -> None:
    """Tile amplicons of ~tile_len with overlap across the sequence.
    For each tile, pick a forward primer near the left edge and a reverse primer near the right edge.
    """
    step = max(1, args.tile_len - args.overlap)
    N = len(sequence)
    rev = revcomp(sequence)
    for pos in range(0, N, step):
        right = min(pos + args.tile_len, N)
        # left window on forward strand
        fregion = sequence[pos : min(pos + args.f_start, N)]
        # right window on reverse strand: start around the right edge
        rev_start = max(0, N - right)
        rwin = rev[rev_start : rev_start + args.r_start]
        fprimer = choose_primer(fregion, sequence, args, partner=None)
        rprimer = choose_primer(rwin, rev, args, partner=fprimer) if fprimer else None
        if not (fprimer and rprimer):
            continue
        # map positions
        fpos = sequence.find(fprimer, pos)
        rev_rprimer = revcomp(rprimer)
        # search near the right edge for the reverse primer
        search_start = max(fpos + len(fprimer), right - args.r_start - args.max_len)
        rpos = sequence.find(rev_rprimer, max(pos, search_start), right)
        if rpos == -1:
            rpos = sequence.rfind(rev_rprimer, pos, right)
        if rpos == -1:
            continue
        rpos_end = rpos + len(rev_rprimer)
        frag = sequence[fpos:rpos_end]
        global_f = args.global_offset + fpos + 1
        global_r = args.global_offset + rpos_end
        out_fh.write(f"{name}_{global_f}-{global_r} Exp.Frag.size {len(frag)} bp of sequence {N} bp\n")
        out_fh.write(f"Forward Primer:\t{fprimer}\tTm={tm(fprimer)}\n")
        out_fh.write(f"Reverse Primer:\t{rprimer}\tTm={tm(rprimer)}\n\n")
        path_fh.write(f">{name}_{global_f}-{global_r}\n{frag}\n\n")

# ---------------------- CLI ----------------------

def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Legacy primer.pl reimplementation in Python")
    ap.add_argument("-i", "--input", required=True, help="Input FASTA file")
    ap.add_argument("-n", "--num-pairs", type=int, default=5, help="Number of primer pairs (non-tiled mode)")
    ap.add_argument("-f", "--f-start", type=int, default=100, help="Search window from 5' end for forward primer")
    ap.add_argument("-r", "--r-start", type=int, default=100, help="Search window from 3' end (on reverse) for reverse primer")
    ap.add_argument("-e", "--evalue", type=float, default=1.0, help="BLAST e-value threshold")
    ap.add_argument("-b", "--skip-blast", action="store_true", help="Skip BLAST specificity check")
    ap.add_argument("-s", "--min-len", type=int, default=20, help="Minimum primer length")
    ap.add_argument("-l", "--max-len", type=int, default=30, help="Maximum primer length")
    ap.add_argument("-t", "--tile-len", type=int, default=None, help="If set, tile across sequence with this approximate amplicon size")
    ap.add_argument("-o", "--overlap", type=int, default=100, help="Overlap between tiles (ignored if not tiling)")
    ap.add_argument("-m", "--min-tm", type=float, default=55.0, help="Minimum Tm (°C)")
    ap.add_argument("-x", "--max-tm", type=float, default=62.0, help="Maximum Tm (°C)")
    ap.add_argument("-M", "--min-gc", type=float, default=40.0, help="Minimum GC%")
    ap.add_argument("-X", "--max-gc", type=float, default=60.0, help="Maximum GC%")
    ap.add_argument("-p", "--max-homopoly", type=int, default=4, help="Max successive identical bases allowed")
    ap.add_argument("-O", "--out-prefix", default=None, help="Output prefix (default: basename of input)")
    ap.add_argument("-q", "--require-upper", action="store_true", help="Require primers to be all uppercase (no masked bases)")
    ap.add_argument("-a", "--start-offset", type=int, default=0, help="Offset within window to begin searching")
    ap.add_argument("-D", "--blast-db", default=None, help="BLAST database basename (makeblastdb -dbtype nucl)")
    ap.add_argument("--global-offset", type=int, default=0, help="Coordinate offset to add to reported positions")
    return ap


def main() -> None:
    ap = build_argparser()
    args = ap.parse_args()

    seqs = read_fasta(args.input)
    out_prefix = args.out_prefix or os.path.splitext(os.path.basename(args.input))[0]
    out_txt = out_prefix + ".primers.txt"
    out_fa = out_prefix + ".fragments.fasta"

    with open(out_txt, "w") as out_fh, open(out_fa, "w") as path_fh:
        for name, seq in seqs.items():
            sys.stderr.write(f"Processing {name} (len={len(seq)})\n")
            if args.tile_len and args.tile_len > 0:
                tiled(seq, name, args, out_fh, path_fh)
            else:
                non_tiled(seq, name, args, out_fh, path_fh)


if __name__ == "__main__":
    main()
