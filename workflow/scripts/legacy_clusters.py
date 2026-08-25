"""Legacy cluster columns, inherited from the precomputed v15 cluster files:
sequence-cluster, interface-cluster, binding-cluster (label on the cluster's FIRST
member only -- the representative), and is_leader / leader_id (redundancy grouping).

NOT recomputed: the three tool clusters need Hammock/MUSTANG/ProBiS, and the
redundancy leader is not a simple alphabetical pick (only ~87.5%); so these are
read from clusters_v15/. New entries in a future update are unclustered until the
tools are re-run (documented limitation; matches v15's inherited-from-v1 design).
"""
import csv
import os
import sys


def load_rep_labels(path):
    """{first_member: label} from `label \\t size \\t m1, m2, ...`."""
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            members = [m.strip() for m in parts[2].split(",") if m.strip()]
            if members:
                out[members[0]] = parts[0]
    return out


def load_redundant(path):
    """{member: (leader, is_leader)} from `leader \\t size \\t members`."""
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            leader = parts[0]
            for m in (x.strip() for x in parts[2].split(",") if x.strip()):
                out[m] = (leader, "yes" if m == leader else "no")
    return out


OUT_COLS = ["id", "sequence-cluster", "interface-cluster", "binding-cluster",
            "is_leader", "leader_id"]


def main():
    p = snakemake.params                                   # noqa: F821
    cdir = os.path.expanduser(p.cluster_dir)
    seq = load_rep_labels(os.path.join(cdir, "sequence.tsv"))
    itf = load_rep_labels(os.path.join(cdir, "interface.tsv"))
    bnd = load_rep_labels(os.path.join(cdir, "binding.tsv"))
    red = load_redundant(os.path.join(cdir, "redundant.tsv"))

    pairs = list(csv.DictReader(open(snakemake.input.pairs), delimiter="\t"))  # noqa: F821
    with open(snakemake.output.clusters, "w") as fh:       # noqa: F821
        fh.write("\t".join(OUT_COLS) + "\n")
        for r in pairs:
            eid = r["id"]
            leader, isl = red.get(eid, ("", ""))
            fh.write("\t".join([eid, seq.get(eid, ""), itf.get(eid, ""),
                                bnd.get(eid, ""), isl, leader]) + "\n")
    print(f"DONE {len(pairs)} rows", file=sys.stderr)


if __name__ == "__main__":
    main()