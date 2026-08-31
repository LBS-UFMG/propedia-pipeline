#!/usr/bin/env python3
"""Verify v15-style X padding on the known-divergent peptides BEFORE the full
re-extraction. Runs the updated extract_pairs sequence logic (polymer_set +
modeled_aa + seq_of) on a few PDBs and compares the peptide sequence to propedia26's
stored value. A match means the padding is restored and it's safe to re-extract.

Usage (on the box with Biopython, e.g. mioglobina in the venv):
    python tests/verify_x_padding.py
CIFs are downloaded from RCSB into a temp dir (needs network); nothing else written.
"""
import gzip
import os
import sys
import tempfile
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "workflow", "scripts"))
import extract_pairs as ep  # noqa: E402
from Bio.PDB.MMCIFParser import MMCIFParser  # noqa: E402

# (pdb, peptide_chain, propedia26 PEPTIDE_SEQ).  1A1M-C-A is an all-standard control
# that must stay unchanged; the rest carry terminal/non-standard residues v15 kept as X.
CASES = [
    ("1A3R", "P", "VKAETRLNPDLQPTEX"),
    ("1A7C", "B", "XTVASSX"),
    ("1APT", "I", "XVVX"),
    ("1ABI", "I", "XPXGGGGGNGDXEEIPEEYL"),
    ("1A1M", "C", "TPYDINQML"),          # control: no X, must be unchanged
]


def main():
    cif_dir = tempfile.mkdtemp(prefix="xpad_")
    parser = MMCIFParser(QUIET=True)
    print(f"{'entry':10} {'ours':24} {'propedia26':24} match?")
    print("-" * 70)
    ok = 0
    for pid, pep, exp in CASES:
        dest = ep.shard_path(cif_dir, pid)
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        try:
            with urllib.request.urlopen(
                    f"https://files.rcsb.org/download/{pid}.cif.gz", timeout=60) as r:
                open(dest, "wb").write(r.read())
            with gzip.open(dest, "rt") as fh:
                fh.read(64)
            model, poly = ep.load_first_model(pid, cif_dir, parser)
            res = ep.modeled_aa(model[pep], poly)
            seq = ep.seq_of(res)
        except Exception as e:                              # noqa: BLE001
            print(f"{pid}-{pep:8} ERROR: {e}")
            continue
        hit = seq == exp
        ok += hit
        print(f"{pid}-{pep:8} {seq:24} {exp:24} {'OK' if hit else 'DIFF'}")
    print("-" * 70)
    print(f"{ok}/{len(CASES)} peptide sequences match propedia26")
    return 0 if ok == len(CASES) else 1


if __name__ == "__main__":
    sys.exit(main())
