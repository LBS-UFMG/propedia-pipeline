"""Interface residues (Propedia 6 A definition): protein residues having any atom
within `cutoff` A of any peptide atom. Reads the extracted two-chain PDB written
by extract_pairs. Output: id, interface_residues = comma-joined protein residue
sequence numbers (author numbering), ascending, each residue once.

NB: v15 leaves ~5.9% of entries blank despite a valid <=6 A contact (likely a
defect). This script emits the correct non-empty list in those cases; it
reproduces v15's list exactly wherever v15 has one (validate to confirm).
"""
import csv
import os
import sys

from Bio.PDB import PDBParser, NeighborSearch


def interface_resseqs(pdb_path, pep_chain, prot_chain, cutoff):
    parser = PDBParser(QUIET=True)
    model = next(iter(parser.get_structure("x", pdb_path)))
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


def main():
    p = snakemake.params                                   # noqa: F821
    pairs = list(csv.DictReader(open(snakemake.input.pairs), delimiter="\t"))  # noqa: F821
    n = written = 0
    with open(snakemake.output.interface, "w") as fh:      # noqa: F821
        fh.write("id\tinterface_residues\n")
        for row in pairs:
            eid = row["id"]
            path = os.path.join(p.pdb_dir, f"{eid}.pdb")
            if not os.path.exists(path):
                continue
            n += 1
            res = interface_resseqs(path, row["pep_chain"], row["prot_chain"], p.cutoff)
            if res is None:
                continue
            fh.write(f"{eid}\t{','.join(str(x) for x in res)}\n")
            written += 1
            if n % 100 == 0:
                print(f"{n} processed", file=sys.stderr)
    print(f"DONE {written}/{n} interface lists written", file=sys.stderr)


if __name__ == "__main__":
    main()