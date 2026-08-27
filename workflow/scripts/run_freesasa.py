"""Surface area via FreeSASA (open-source NACCESS replacement — reproducible,
pip-installable). Computes the seven v15 ASA/BSA columns per pair and applies
the BSA>0 filter. Runs on the extracted two-chain mmCIF (complex), splitting into
isolated protein and peptide to get ΔASA/BSA per the paper's equation 1:
    BSA = (ASA_protein + ASA_peptide - ASA_complex) / 2

mmCIF path: the structure is parsed with Biopython (auth chain IDs) and handed to
FreeSASA via `structureFromBioPDB`, so multi-character chain IDs and large atom
counts work where the old PDB-text path could not. This path is numerically
identical to the previous PDB-file path on structures both can represent
(verified: Δ=0.000% on ASA_complex/peptide/protein), so the validated r≈1.0 vs
v15/NACCESS is preserved.

NOTE: FreeSASA != NACCESS numerically (different algorithm/radii); this is a
reproducible re-implementation, expected to correlate closely, not match exactly.
"""
import csv
import os
import sys

import freesasa
from Bio.PDB import MMCIFParser
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure

# quiet FreeSASA's stderr chatter
freesasa.setVerbosity(freesasa.silent)

# Match NACCESS as closely as FreeSASA allows: NACCESS radii + Lee-Richards
# (NACCESS's algorithm) at high slice count. Reproducible, still pip-only.
_NACCESS = freesasa.Classifier.getStandardClassifier("naccess")
_PARAMS = freesasa.Parameters({"algorithm": freesasa.LeeRichards,
                               "probe-radius": 1.4, "n-slices": 100})

_PARSER = MMCIFParser(QUIET=True)


def load_model(cif_path):
    """First model of a mmCIF file (first NMR conformer)."""
    return next(iter(_PARSER.get_structure("x", cif_path)))


def sub_structure(model, chain_ids):
    """A Bio.PDB Structure holding copies of only the named chains."""
    st = Structure("x")
    m = Model(0)
    st.add(m)
    for ch in model:
        if ch.id in chain_ids:
            m.add(ch.copy())
    return st


def areas_of(struct):
    """(total_area, {chain_id: area}) from ONE FreeSASA calc (NACCESS-matched)."""
    result = freesasa.calc(freesasa.structureFromBioPDB(struct, _NACCESS), _PARAMS)
    per = {ch: sum(r.total for r in resd.values())
           for ch, resd in result.residueAreas().items()}
    return result.totalArea(), per


def asa_of(struct):
    """Total ASA of a Bio.PDB Structure."""
    return areas_of(struct)[0]


def process(model, pep_chain, prot_chain):
    present = {c.id for c in model}
    if pep_chain not in present or prot_chain not in present:
        return None
    asa_complex, per_chain = areas_of(sub_structure(model, {pep_chain, prot_chain}))
    asa_pep = asa_of(sub_structure(model, {pep_chain}))
    asa_prot = asa_of(sub_structure(model, {prot_chain}))
    bsa = (asa_prot + asa_pep - asa_complex) / 2.0
    # true ΔASA: isolated area minus area still exposed within the complex
    b_pep = asa_pep - per_chain.get(pep_chain, 0.0)     # peptide area buried
    b_pro = asa_prot - sum(a for c, a in per_chain.items() if c != pep_chain)
    bpp = (100.0 * b_pep / asa_pep) if asa_pep > 0 else 0.0
    return {
        "ASA_Complex": round(asa_complex, 2),
        "ASA_Peptide": round(asa_pep, 2),
        "ASA_Protein": round(asa_prot, 2),
        "BSA": round(bsa, 2),
        "BPepA": round(b_pep, 2),
        "BProA": round(b_pro, 2),
        "BPP%": round(bpp, 2),
    }


def main():
    p = snakemake.params                                   # noqa: F821
    pairs = list(csv.DictReader(open(snakemake.input.pairs), delimiter="\t"))  # noqa: F821
    bsa_min = p.bsa_threshold
    cols = ["id", "ASA_Complex", "ASA_Peptide", "ASA_Protein",
            "BSA", "BPepA", "BProA", "BPP%"]
    n = kept = 0
    with open(snakemake.output.surface, "w") as fh:        # noqa: F821
        fh.write("\t".join(cols) + "\n")
        for row in pairs:
            eid = row["id"]
            path = os.path.join(p.cif_dir, f"{eid}.cif")
            if not os.path.exists(path):
                continue
            n += 1
            try:
                r = process(load_model(path), row["pep_chain"], row["prot_chain"])
            except Exception:                              # noqa: BLE001
                r = None
            if r is None:
                continue
            if r["BSA"] <= bsa_min:        # paper's BSA>0 selection
                continue
            kept += 1
            r["id"] = eid
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
            if n % 100 == 0:
                print(f"{n} processed, {kept} kept (BSA>{bsa_min})", file=sys.stderr)
    print(f"DONE {n} processed, {kept} passed BSA>{bsa_min}", file=sys.stderr)


if __name__ == "__main__":
    main()
