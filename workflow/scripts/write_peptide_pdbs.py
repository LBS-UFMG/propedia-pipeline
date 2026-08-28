"""Write peptide-only PDB files from the extracted two-chain mmCIF pairs, for SIGNA.
Propedia 26 computes structural signatures on the peptide structure alone.

SIGNA has only a weak hand-rolled mmCIF reader, so we keep a PDB shim JUST here:
the peptide is <=50 residues and always fits PDB losslessly. The peptide chain is
copied into a fresh single-chain structure and its ID forced to 'A', so even a
multi-character author chain ID (which PDB cannot store) serializes cleanly. aCSM
signatures are computed on coordinates and ignore the chain label.
"""
import csv
import os
import shutil
import sys

from Bio.PDB import MMCIFParser, PDBIO
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure


def main():
    p = snakemake.params                                   # noqa: F821
    pairs = list(csv.DictReader(open(snakemake.input.pairs), delimiter="\t"))  # noqa: F821
    # clear stale peptide PDBs first — SIGNA scans this whole folder, so orphans
    # from a previous (larger) run would be counted as extra signatures.
    if os.path.isdir(p.pep_pdb_dir):
        shutil.rmtree(p.pep_pdb_dir)
    os.makedirs(p.pep_pdb_dir, exist_ok=True)
    parser = MMCIFParser(QUIET=True)
    io = PDBIO()
    n = 0
    for row in pairs:
        eid = row["id"]
        pep_chain = row["pep_chain"]
        src = os.path.join(p.pair_cif_dir, f"{eid}.cif")
        if not os.path.exists(src):
            continue
        try:
            model = next(iter(parser.get_structure(eid, src)))
        except Exception:                                  # noqa: BLE001
            continue
        if pep_chain not in model:
            continue
        # fresh single-chain structure; force a PDB-legal single-char chain id
        st = Structure(eid)
        m = Model(0)
        st.add(m)
        pep = model[pep_chain].copy()
        pep.id = "A"
        m.add(pep)
        io.set_structure(st)
        io.save(os.path.join(p.pep_pdb_dir, f"{eid}.pdb"))
        n += 1

    # create the marker Snakemake declared as this rule's output
    with open(snakemake.output.marker, "w") as fh:   # noqa: F821
        fh.write(f"wrote {n} peptide-only PDBs\n")

    print(f"wrote {n} peptide-only PDBs to {p.pep_pdb_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
