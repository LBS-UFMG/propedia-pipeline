"""Structural signatures via SIGNA (aCSM-ALL). Runs on peptide-only PDBs.
Propedia 26 params: cutoff_limit=10, cutoff_step=0.2, cumulative=True -> 1800 features.
SIGNA writes a CSV keyed by file PATH; we normalize to entry ID and rejoin."""
import csv
import os
import sys


def main():
    p = snakemake.params                                   # noqa: F821
    signa_dir = os.path.expanduser(p.signa_dir)
    sys.path.insert(0, signa_dir)
    import signa                                           # noqa: E402

    pep_dir = os.path.abspath(p.pep_pdb_dir)
    raw_out = os.path.abspath(p.raw_output)
    os.makedirs(os.path.dirname(raw_out), exist_ok=True)

    # SIGNA prints a lot; run it and let it write raw_out
    signa.read_folder(folder=pep_dir, signa_type="acsm-all",
                      cumulative=p.cumulative, output=raw_out,
                      cutoff_limit=p.cutoff_limit, cutoff_step=p.cutoff_step,
                      format="pdb")

    # Normalize: first column is a path like /.../1A1M-C-A.pdb -> entry id
    n = 0
    with open(raw_out) as fin, open(snakemake.output.signatures, "w") as fout:  # noqa: F821
        for line in fin:
            parts = line.rstrip("\n").split(",")
            if not parts or not parts[0]:
                continue
            eid = os.path.splitext(os.path.basename(parts[0]))[0]
            fout.write(eid + "," + ",".join(parts[1:]) + "\n")
            n += 1
    print(f"DONE {n} signatures normalized", file=sys.stderr)


if __name__ == "__main__":
    main()