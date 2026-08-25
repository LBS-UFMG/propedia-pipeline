"""Write peptide-only PDB files from the extracted two-chain pairs, for SIGNA.
Propedia 26 computes structural signatures on the peptide structure alone."""
import csv
import os
import sys


def main():
    p = snakemake.params                                   # noqa: F821
    pairs = list(csv.DictReader(open(snakemake.input.pairs), delimiter="\t"))  # noqa: F821
    os.makedirs(p.pep_pdb_dir, exist_ok=True)
    n = 0
    for row in pairs:
        eid = row["id"]
        pep_chain = row["pep_chain"]
        src = os.path.join(p.pair_pdb_dir, f"{eid}.pdb")
        if not os.path.exists(src):
            continue
        out = os.path.join(p.pep_pdb_dir, f"{eid}.pdb")
        with open(src) as fh, open(out, "w") as ofh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM", "TER")):
                    # chain ID is column 22 (index 21) in PDB format
                    if len(line) > 21 and line[21] == pep_chain:
                        ofh.write(line)
            ofh.write("END\n")
        n += 1

    # create the marker Snakemake declared as this rule's output
    with open(snakemake.output.marker, "w") as fh:   # noqa: F821
        fh.write(f"wrote {n} peptide-only PDBs\n")

    print(f"wrote {n} peptide-only PDBs to {p.pep_pdb_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()