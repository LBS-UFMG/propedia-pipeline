"""Rule 3 (sample): from each candidate CIF, classify chains by MODELED
amino-acid count, find peptide(2-50)-protein(>50) pairs in contact (<=cutoff A),
write the two-chain complex as a clean ATOM-only PDB, and record sequences.
"""
import gzip
import os
import sys
import tempfile

from Bio.PDB import MMCIFParser, NeighborSearch, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.IUPACData import protein_letters_3to1_extended


def one(resname):
    code = protein_letters_3to1_extended.get(resname.capitalize(), "X")
    return code if len(code) == 1 else "X"


def modeled_aa(chain):
    # Count any amino-acid residue that is part of the polymer (het flag 'H_' or ' '),
    # standard or modified, matching v15's convention of keeping modified residues as X.
    return [r for r in chain if is_aa(r, standard=False)]


def seq_of(residues):
    return "".join(one(r.resname) for r in residues)


class PairSelect(Select):
    def __init__(self, chain_ids):
        self.chain_ids = set(chain_ids)

    def accept_model(self, model):
        return True

    def accept_chain(self, chain):
        return chain.id in self.chain_ids

    def accept_residue(self, residue):
        return is_aa(residue, standard=False)


def shard_path(cif_dir, pid):
    sub = pid[1:3].lower() if len(pid) >= 3 else "_"
    return os.path.join(cif_dir, sub, f"{pid}.cif.gz")


def load_first_model(pid, cif_dir, parser):
    path = shard_path(cif_dir, pid)
    with gzip.open(path, "rt") as fh:
        data = fh.read()
    tmp = tempfile.NamedTemporaryFile("w", suffix=".cif", delete=False)
    try:
        tmp.write(data)
        tmp.close()
        structure = parser.get_structure(pid, tmp.name)
    finally:
        os.unlink(tmp.name)
    return next(iter(structure))            # first model (first NMR conformer)


def process(pid, p, io, parser):
    rows = []
    try:
        model = load_first_model(pid, p.cif_dir, parser)
    except Exception as exc:                # noqa: BLE001
        return rows, f"parse_error"
    chain_res = {c.id: modeled_aa(c) for c in model}
    chain_res = {c: r for c, r in chain_res.items() if r}
    peptides = {c: r for c, r in chain_res.items()
                if p.pep_min <= len(r) <= p.pep_max}
    proteins = {c: r for c, r in chain_res.items() if len(r) >= p.prot_min}
    if not peptides or not proteins:
        return rows, "no_pair"

    ns = NeighborSearch(list(model.get_atoms()))
    for pep_id, pep_res in peptides.items():
        for prot_id, prot_res in proteins.items():
            if pep_id == prot_id:
                continue
            contact = False
            for atom in (a for r in pep_res for a in r.get_atoms()):
                for near in ns.search(atom.coord, p.cutoff, level="A"):
                    if near.get_parent().get_parent().id == prot_id:
                        contact = True
                        break
                if contact:
                    break
            if not contact:
                continue
            entry_id = f"{pid}-{pep_id}-{prot_id}"
            io.set_structure(model)
            io.save(os.path.join(p.pdb_out_dir, f"{entry_id}.pdb"),
                    PairSelect([prot_id, pep_id]))
            rows.append({"id": entry_id, "pdb": pid,
                         "pep_chain": pep_id, "prot_chain": prot_id,
                         "pep_size": len(pep_res), "prot_size": len(prot_res),
                         "pep_seq": seq_of(pep_res),
                         "prot_seq": seq_of(prot_res)})
    return rows, "ok"


def main():
    p = snakemake.params                                   # noqa: F821
    os.makedirs(p.pdb_out_dir, exist_ok=True)
    sample = [l.strip() for l in open(snakemake.input.sample) if l.strip()]  # noqa: F821
    parser = MMCIFParser(QUIET=True)
    io = PDBIO()
    all_rows, stats = [], {}
    for i, pid in enumerate(sample, 1):
        rows, status = process(pid, p, io, parser)
        all_rows.extend(rows)
        stats[status] = stats.get(status, 0) + 1
        if i % 50 == 0:
            print(f"{i}/{len(sample)} pairs={len(all_rows)} {stats}",
                  file=sys.stderr)
    cols = ["id", "pdb", "pep_chain", "prot_chain",
            "pep_size", "prot_size", "pep_seq", "prot_seq"]
    with open(snakemake.output.pairs, "w") as fh:          # noqa: F821
        fh.write("\t".join(cols) + "\n")
        for r in all_rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"DONE pairs={len(all_rows)} {stats}", file=sys.stderr)


if __name__ == "__main__":
    main()