"""Multipro extraction (Phase 1): group pep-pro pairs by (PDB, peptide chain);
an entry is Multipro when the peptide contacts >=2 protein chains. Emits the
grouping plus the per-protein-chain columns that are DERIVABLE from the pep-pro
pairs (colon-joined in item order). Recomputed features (surface, PRODIGY on the
multi-chain complex) are Phase 2.

Master set = pairs that pass BSA>0 (present in surface), matching the pep-pro CSV.
"""
import csv
import os
import sys
from collections import defaultdict


def colon(vals):
    return ":".join(str(v) for v in vals)


OUT_COLS = ["cluster_id", "count", "items", "PDB_ID", "PEPTIDE_CHAIN",
            "PROTEIN_CHAIN", "PEPTIDE_SIZE", "PEPTIDE_SEQ",
            "PROTEIN_SIZE", "PROTEIN_SEQ"]


def main():
    pairs = list(csv.DictReader(open(snakemake.input.pairs), delimiter="\t"))  # noqa: F821
    surface = {r["id"] for r in csv.DictReader(open(snakemake.input.surface),  # noqa: F821
                                               delimiter="\t")}
    kept = [r for r in pairs if r["id"] in surface]

    groups = defaultdict(list)
    for r in kept:
        groups[(r["pdb"], r["pep_chain"])].append(r)

    n = 0
    with open(snakemake.output.multipro, "w") as fh:       # noqa: F821
        fh.write("\t".join(OUT_COLS) + "\n")
        for (pdb, pep), members in sorted(groups.items()):
            members = sorted(members, key=lambda r: r["prot_chain"])
            prot_chains = [r["prot_chain"] for r in members]
            if len(set(prot_chains)) < 2:
                continue
            row = {
                "cluster_id": f"{pdb}-{pep}",
                "count": len(members),
                "items": colon(r["id"] for r in members),
                "PDB_ID": pdb,
                "PEPTIDE_CHAIN": pep,
                "PROTEIN_CHAIN": colon(prot_chains),
                "PEPTIDE_SIZE": members[0]["pep_size"],
                "PEPTIDE_SEQ": members[0]["pep_seq"],
                "PROTEIN_SIZE": colon(r["prot_size"] for r in members),
                "PROTEIN_SEQ": colon(r["prot_seq"] for r in members),
            }
            fh.write("\t".join(str(row[c]) for c in OUT_COLS) + "\n")
            n += 1
    print(f"DONE {n} multipro entries", file=sys.stderr)


if __name__ == "__main__":
    main()