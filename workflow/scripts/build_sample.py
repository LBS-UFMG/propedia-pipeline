"""Pick a sample of oracle (v15) PDB IDs that are already downloaded,
so we can validate extraction against the existing database right now."""
import csv
import os
import sys


def shard_path(cif_dir, pid):
    sub = pid[1:3].lower() if len(pid) >= 3 else "_"
    return os.path.join(cif_dir, sub, f"{pid}.cif.gz")


def main():
    p = snakemake.params                                   # noqa: F821
    # oracle is the ;-delimited v15 table; PDB IDs live in the PDB_ID column
    seen, oracle = set(), []
    with open(snakemake.input.oracle) as fh:               # noqa: F821
        for row in csv.DictReader(fh, delimiter=";"):
            pid = (row.get("PDB_ID") or "").strip().upper()
            if pid and pid not in seen:
                seen.add(pid)
                oracle.append(pid)
    sample = []
    for pid in oracle:
        if os.path.exists(shard_path(p.cif_dir, pid)):
            sample.append(pid)
            if len(sample) >= p.sample_size:
                break
    with open(snakemake.output.sample, "w") as fh:         # noqa: F821
        fh.write("\n".join(sample) + "\n")
    print(f"sample: {len(sample)} of {p.sample_size} requested "
          f"(from {len(oracle)} oracle IDs)", file=sys.stderr)


if __name__ == "__main__":
    main()