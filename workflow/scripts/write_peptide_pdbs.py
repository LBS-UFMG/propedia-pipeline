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
import sys

from Bio.PDB import MMCIFParser, PDBIO
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure


def main():
    p = snakemake.params                                   # noqa: F821
    pairs = list(csv.DictReader(open(snakemake.input.pairs), delimiter="\t"))  # noqa: F821
    os.makedirs(p.pep_pdb_dir, exist_ok=True)              # resume: do NOT clear
    parser = MMCIFParser(QUIET=True)
    io = PDBIO()
    n = skipped = 0
    want = set()
    for row in pairs:
        eid = row["id"]
        pep_chain = row["pep_chain"]
        out = os.path.join(p.pep_pdb_dir, f"{eid}.pdb")
        want.add(f"{eid}.pdb")
        if os.path.exists(out):        # resume: already written on a prior run
            skipped += 1
            continue
        src = os.path.join(p.pair_cif_dir, f"{eid}.cif")
        if not os.path.exists(src):
            want.discard(f"{eid}.pdb")
            continue
        try:
            model = next(iter(parser.get_structure(eid, src)))
        except Exception:                                  # noqa: BLE001
            want.discard(f"{eid}.pdb")
            continue
        if pep_chain not in model:
            want.discard(f"{eid}.pdb")
            continue
        # fresh single-chain structure; force a PDB-legal single-char chain id
        st = Structure(eid)
        m = Model(0)
        st.add(m)
        pep = model[pep_chain].copy()
        pep.id = "A"
        m.add(pep)
        io.set_structure(st)
        # atomic: write to tmp then rename, so a killed run never leaves a partial
        # PDB that resume would treat as complete
        tmp = out + ".tmp"
        io.save(tmp)
        os.replace(tmp, out)
        n += 1

    # prune stale peptide PDBs (smaller/different sample) — SIGNA scans this folder
    for f in os.listdir(p.pep_pdb_dir):
        if f.endswith(".pdb") and f not in want:
            os.remove(os.path.join(p.pep_pdb_dir, f))

    with open(snakemake.output.marker, "w") as fh:   # noqa: F821
        fh.write(f"wrote {n} (+{skipped} resumed) peptide-only PDBs\n")
    print(f"wrote {n} new, {skipped} resumed peptide-only PDBs to {p.pep_pdb_dir}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
