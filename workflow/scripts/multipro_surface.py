"""Multipro Phase 2b: surface area on the assembled multi-chain complex.
Same FreeSASA convention as run_freesasa (validated BSA r=1.000 vs v15), but
'protein' = ALL protein chains together and 'peptide' = the single peptide chain:
    ASA_Complex = whole multi-chain structure
    ASA_Protein = all protein chains (no peptide)
    ASA_Peptide = peptide chain alone
    BSA = (ASA_Protein + ASA_Peptide - ASA_Complex) / 2
Single values per Multipro entry (not per chain). Reads the multi-chain mmCIF.
"""
import csv
import os
import sys

import run_freesasa as fs   # reuse sub_structure/areas_of/asa_of + freesasa setup


def process(model, pep_chain, prot_chains):
    present = {c.id for c in model}
    if pep_chain not in present or not (prot_chains & present):
        return None
    asa_complex, per_chain = fs.areas_of(fs.sub_structure(model, {pep_chain} | prot_chains))
    asa_pep = fs.asa_of(fs.sub_structure(model, {pep_chain}))
    asa_prot = fs.asa_of(fs.sub_structure(model, prot_chains))
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
            path = os.path.join(p.cif_dir, f"{cid}.cif")
            if not os.path.exists(path):
                continue
            try:
                r = process(fs.load_model(path), row["PEPTIDE_CHAIN"],
                            set(row["PROTEIN_CHAIN"].split(":")))
            except Exception:                              # noqa: BLE001
                r = None
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
