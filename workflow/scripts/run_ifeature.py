"""Sequence signatures via iFeature. Builds a FASTA from peptide sequences,
runs each descriptor once (batch over all sequences), joins into one wide table.

Paper's 1248-feature vector = 9 descriptors EXCLUDING CTriad (the 10 listed sum
to 1591; 1591-343=1248). CTriad is optional here (include_ctriad).

Peptides with fewer than min_signature_len STANDARD residues are excluded:
iFeature ignores non-standard residues (e.g. X), and DDE/CTDT divide by
(len-1)/(pair count), which raises ZeroDivisionError once a sequence collapses
to <2 standard residues (e.g. 'LXX' -> 'L'). Excluded IDs are logged.
"""
import concurrent.futures as cf
import csv
import hashlib
import os
import subprocess
import sys
import tempfile

import checkpoint

VERSION = "1"

# 9 descriptors that sum to the paper's 1248 features; CTriad optional (+343).
CORE_9 = ["AAC", "DPC", "DDE", "GAAC", "GDPC", "GTPC", "CTDC", "CTDT", "CTDD"]
STD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def build_fasta(pairs, seq_col, path, min_len=3):
    seen, excluded = set(), []
    with open(path, "w") as fh:
        for row in pairs:
            sid = row["id"]
            # keep ONLY standard residues (iFeature ignores X anyway; make it explicit)
            seq = "".join(c for c in row[seq_col].upper() if c in STD_AA)
            if sid in seen:
                continue
            if len(seq) < min_len:            # length of STANDARD residues
                excluded.append((sid, len(seq)))
                continue
            seen.add(sid)
            fh.write(f">{sid}\n{seq}\n")
    return len(seen), excluded


def run_descriptor(ifeature_dir, fasta, desc, out_tsv):
    proc = subprocess.run(
        ["python3", "iFeature.py", "--file", fasta, "--type", desc, "--out", out_tsv],
        cwd=ifeature_dir, capture_output=True, text=True, timeout=1800)
    if proc.returncode != 0 or not os.path.exists(out_tsv):
        print(f"  {desc} STDERR: {proc.stderr.strip()[:300]}", file=sys.stderr)
        return False
    return True


def ensure_descriptor(args):
    """Top-level (picklable) worker: compute descriptor into its cache path, or
    reuse it if a prior run already produced it. Returns (desc, out_tsv, ok)."""
    desc, ifeature_dir, fasta, workdir = args
    out_tsv = os.path.join(workdir, f"{desc}.tsv")
    if os.path.exists(out_tsv):
        return desc, out_tsv, True             # resumed
    ok = run_descriptor(ifeature_dir, os.path.abspath(fasta), desc,
                        os.path.abspath(out_tsv))
    return desc, out_tsv, ok


def load_matrix(path):
    rows = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")[1:]  # drop leading '#'
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows[f[0]] = f[1:]
    return header, rows


def main():
    p = snakemake.params                                   # noqa: F821
    ifeature_dir = os.path.expanduser(p.ifeature_dir)
    min_len = getattr(p, "min_signature_len", 3)

    descriptors = list(CORE_9)
    if getattr(p, "include_ctriad", False):
        descriptors.append("CTriad")

    pairs = list(csv.DictReader(open(snakemake.input.pairs), delimiter="\t"))  # noqa: F821
    tmpdir = tempfile.mkdtemp()
    fasta = os.path.join(tmpdir, "peptides.fasta")
    n, excluded = build_fasta(pairs, "pep_seq", fasta, min_len=min_len)
    print(f"FASTA: {n} unique peptides (excluded {len(excluded)} shorter than "
          f"{min_len} standard residues)", file=sys.stderr)

    with open(snakemake.output.excluded, "w") as ef:       # noqa: F821
        ef.write("id\tlength\n")
        for sid, ln in excluded:
            ef.write(f"{sid}\t{ln}\n")

    # Cache each descriptor's output, content-addressed by the FASTA, so an
    # interrupted run resumes (skips completed descriptors) and a changed peptide
    # set recomputes. Descriptors run in parallel (each is a batch over the FASTA).
    threads = getattr(snakemake, "threads", 1) or 1        # noqa: F821
    with open(fasta, "rb") as _fh:
        fasta_hash = hashlib.sha1(_fh.read()).hexdigest()[:16]
    workdir = checkpoint.namespace(p.ckpt, VERSION,
                                   {"desc": descriptors, "min_len": min_len,
                                    "fasta": fasta_hash})

    jobs = [(desc, ifeature_dir, fasta, workdir) for desc in descriptors]
    pool = min(threads, len(descriptors))
    if pool > 1:
        with cf.ProcessPoolExecutor(max_workers=pool) as ex:
            done = list(ex.map(ensure_descriptor, jobs))
    else:
        done = [ensure_descriptor(j) for j in jobs]

    all_headers, all_rows = [], {}
    order = None
    for desc in descriptors:                    # deterministic column order
        out_tsv, ok = next((o, k) for (d, o, k) in done if d == desc)
        if not ok or not os.path.exists(out_tsv):
            print(f"  {desc}: FAILED", file=sys.stderr)
            continue
        hdr, rows = load_matrix(out_tsv)
        all_headers += [f"{desc}_{h}" for h in hdr]
        if order is None:
            order = list(rows.keys())
        for sid, vals in rows.items():
            all_rows.setdefault(sid, []).extend(vals)
        print(f"  {desc}: {len(hdr)} features", file=sys.stderr)

    with open(snakemake.output.signatures, "w") as fh:     # noqa: F821
        fh.write("id\t" + "\t".join(all_headers) + "\n")
        for sid in (order or []):
            if len(all_rows.get(sid, [])) == len(all_headers):
                fh.write(sid + "\t" + "\t".join(all_rows[sid]) + "\n")

    print(f"DONE {len(all_headers)} features x {len(order or [])} peptides",
          file=sys.stderr)


if __name__ == "__main__":
    main()