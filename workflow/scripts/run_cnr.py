"""CNR clustering: group extracted pairs by identical peptide sequence (100%
identity), matching v15's seq100_clusters grouping. Representative = smallest
entry id in the group present in this run (v15's rep is DB-global; grouping is
reproduced, label may differ on a sample)."""
import csv
import sys
from collections import defaultdict


def main():
    pairs = list(csv.DictReader(open(snakemake.input.pairs), delimiter="\t"))  # noqa: F821
    groups = defaultdict(list)
    for r in pairs:
        groups[r["pep_seq"]].append(r["id"])

    with open(snakemake.output.clusters, "w") as fh:       # noqa: F821
        fh.write("id\tpep_seq\tcnr_cluster\tcnr_size\tis_representative\n")
        n_groups = 0
        for seq, ids in groups.items():
            rep = sorted(ids)[0]
            n_groups += 1
            for eid in sorted(ids):
                fh.write(f"{eid}\t{seq}\t{rep}\t{len(ids)}\t"
                         f"{'yes' if eid == rep else 'no'}\n")
    print(f"DONE {len(pairs)} entries -> {n_groups} CNR clusters", file=sys.stderr)


if __name__ == "__main__":
    main()