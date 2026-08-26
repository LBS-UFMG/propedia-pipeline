"""Multipro Phase 2b: surface area on the assembled multi-chain complex.
Same FreeSASA convention as run_freesasa (validated BSA r=1.000 vs v15), but
'protein' = ALL protein chains together and 'peptide' = the single peptide chain:
    ASA_Complex = whole multi-chain PDB
    ASA_Protein = all protein chains (no peptide)
    ASA_Peptide = peptide chain alone
    BSA = (ASA_Protein + ASA_Peptide - ASA_Complex) / 2
Single values per Multipro entry (not per chain).
"""
import csv
import os
import sys

import run_freesasa as fs   # reuse asa_of + freesasa silent setup


def process(pdb_path, pep_chain, prot_chains):
    with open(pdb_path) as fh:
        atom_lines = [l for l in fh if l.startswith(("ATOM", "HETATM"))]
    pep = [l for l in atom_lines if len(l) > 21 and l[21] == pep_chain]
    prot = [l for l in atom_lines if len(l) > 21 and l[21] in prot_chains]
    if not pep or not prot:
        return None
    asa_complex, per_chain = fs.areas_of(atom_lines)
    asa_pep = fs.asa_of(pep)
    asa_prot = fs.asa_of(prot)
    bsa = (asa_prot + asa_pep - asa_complex) / 2.0
    b_pep = asa_pep - per_chain.get(pep_chain, 0.0)
    b_pro = asa_prot - sum(a for c, a in per_chain.items() if c in prot_chains)
    bpp = (100.0 * b_pep / asa_pep) if asa_pep > 0 else 0.0
    return {
        "ASA_Complex": round(asa_complex, 2), "ASA_Peptide": round(asa_pep, 2),
        "ASA_Protein": round(asa_prot, 2), "BSA": round(bsa, 2),
        "BPepA": round(b_pep, 2), "BProA": round(b_pro, 2), "BPP%": round(bpp, 2),
    }


def main():
    p = snakemake.params                                   # noqa: F821
    rows = list(csv.DictReader(open(snakemake.input.multipro), delimiter="\t"))  # noqa: F821
    cols = ["cluster_id", "ASA_Complex", "ASA_Protein", "ASA_Peptide",
            "BProA", "BPepA", "BPP%", "BSA"]
    n = 0
    with open(snakemake.output.surface, "w") as fh:        # noqa: F821
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            cid = row["cluster_id"]
            path = os.path.join(p.pdb_dir, f"{cid}.pdb")
            if not os.path.exists(path):
                continue
            r = process(path, row["PEPTIDE_CHAIN"], set(row["PROTEIN_CHAIN"].split(":")))
            if r is None:
                continue
            r["cluster_id"] = cid
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
            n += 1
            if n % 50 == 0:
                print(f"{n} surfaces", file=sys.stderr)
    print(f"DONE {n} multipro surfaces", file=sys.stderr)


if __name__ == "__main__":
    main()