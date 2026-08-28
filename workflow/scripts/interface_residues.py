"""Interface residues (Propedia 6 A definition): protein residues having any atom
within `cutoff` A of any peptide atom. Reads the extracted two-chain mmCIF written
by extract_pairs. Output: id, interface_residues = comma-joined protein residue
sequence numbers (author numbering), ascending, each residue once.

NB: v15 leaves ~5.9% of entries blank despite a valid <=6 A contact (likely a
defect). This script emits the correct non-empty list in those cases; it
reproduces v15's list exactly wherever v15 has one (validate to confirm).
"""
import csv
import os
import sys

from Bio.PDB import MMCIFParser, NeighborSearch

import checkpoint

VERSION = "2"            # bump on logic change -> invalidates old checkpoints


def interface_resseqs(cif_path, pep_chain, prot_chain, cutoff):
    parser = MMCIFParser(QUIET=True)
    model = next(iter(parser.get_structure("x", cif_path)))
    if pep_chain not in model or prot_chain not in model:
        return None
    pep_atoms = list(model[pep_chain].get_atoms())
    prot_atoms = list(model[prot_chain].get_atoms())
    if not pep_atoms or not prot_atoms:
        return None
    ns = NeighborSearch(prot_atoms)
    hits = set()
    for atom in pep_atoms:
        for res in ns.search(atom.coord, cutoff, level="R"):
            hits.add(res.id[1])          # author residue sequence number (resseq)
    return sorted(hits)


def worker(item):
    """worker(item) -> (resseq_list|None, status)."""
    eid, cif_path, pep, prot, cutoff = item
    try:
        res = interface_resseqs(cif_path, pep, prot, cutoff)
    except Exception as exc:                               # noqa: BLE001
        return None, f"retry:exc:{exc}"
    return (res, "ok") if res is not None else (None, "no_interface")


def main():
    p = snakemake.params                                   # noqa: F821
    threads = getattr(snakemake, "threads", 1) or 1        # noqa: F821
    pairs = list(csv.DictReader(open(snakemake.input.pairs), delimiter="\t"))  # noqa: F821

    items = [(r["id"], os.path.join(p.cif_dir, f"{r['id']}.cif"),
              r["pep_chain"], r["prot_chain"], p.cutoff) for r in pairs
             if os.path.exists(os.path.join(p.cif_dir, f"{r['id']}.cif"))]
    workdir = checkpoint.namespace(p.ckpt, VERSION, {"cutoff": p.cutoff})
    results = checkpoint.run(items, worker, workdir, threads=threads,
                             id_of=lambda it: it[0], stage="interface")

    written = 0
    with open(snakemake.output.interface, "w") as fh:      # noqa: F821
        fh.write("id\tinterface_residues\n")
        for row in pairs:
            res = results.get(row["id"])
            if not res or res["status"] != "ok" or res["record"] is None:
                continue
            fh.write(f"{row['id']}\t{','.join(str(x) for x in res['record'])}\n")
            written += 1
    print(f"DONE {written} interface lists written", file=sys.stderr)


if __name__ == "__main__":
    main()