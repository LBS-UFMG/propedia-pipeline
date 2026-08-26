"""Multipro Phase 2a: write the multi-chain PDB for each Multipro entry
(peptide chain + all its protein chains), from the CIF. Reuses extract_pairs for
CIF parsing / chain selection so atom content matches the pep-pro PDBs. The model
is cached per PDB (several Multipro entries can share one structure).
"""
import csv
import os
import sys

from Bio.PDB import MMCIFParser, PDBIO

import extract_pairs as ep


def main():
    p = snakemake.params                                   # noqa: F821
    os.makedirs(p.pdb_out_dir, exist_ok=True)
    rows = list(csv.DictReader(open(snakemake.input.multipro), delimiter="\t"))  # noqa: F821
    parser = MMCIFParser(QUIET=True)
    io = PDBIO()
    cache = {}
    n = 0
    for r in rows:
        pid = r["PDB_ID"]
        chains = {r["PEPTIDE_CHAIN"]} | set(r["PROTEIN_CHAIN"].split(":"))
        if pid not in cache:
            try:
                cache[pid] = ep.load_first_model(pid, p.cif_dir, parser)
            except Exception:                              # noqa: BLE001
                cache[pid] = None
        model = cache[pid]
        if model is None:
            continue
        io.set_structure(model)
        io.save(os.path.join(p.pdb_out_dir, f'{r["cluster_id"]}.pdb'),
                ep.PairSelect(chains))
        n += 1
        if n % 50 == 0:
            print(f"{n} multipro pdbs", file=sys.stderr)

    open(os.path.join(p.pdb_out_dir, ".written"), "w").write(f"{n}\n")
    print(f"DONE {n} multipro pdbs", file=sys.stderr)


if __name__ == "__main__":
    main()