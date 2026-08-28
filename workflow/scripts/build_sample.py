"""Pick a REPRESENTATIVE, seeded-random sample of oracle (v15) PDB IDs that are
already downloaded, so we validate against the existing database on an unbiased
slice spanning all eras/sizes — not the alphabetically-first (oldest, smallest)
entries the previous first-N logic returned. Reproducible via `sample.seed`."""
import csv
import os
import random
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
    # restrict to downloaded structures, then take a seeded-random sample across
    # the whole downloaded set (unbiased w.r.t. deposition era / structure size)
    present = [pid for pid in oracle if os.path.exists(shard_path(p.cif_dir, pid))]
    rng = random.Random(p.seed)
    k = min(p.sample_size, len(present))
    sample = sorted(rng.sample(present, k))
    with open(snakemake.output.sample, "w") as fh:         # noqa: F821
        fh.write("\n".join(sample) + "\n")
    print(f"sample: {len(sample)} random of {len(present)} downloaded oracle IDs "
          f"(seed={p.seed}, requested {p.sample_size})", file=sys.stderr)


if __name__ == "__main__":
    main()
