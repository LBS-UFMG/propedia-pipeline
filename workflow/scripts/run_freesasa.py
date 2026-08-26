"""Surface area via FreeSASA (open-source NACCESS replacement — reproducible,
pip-installable). Computes the seven v15 ASA/BSA columns per pair and applies
the BSA>0 filter. Runs on the extracted two-chain PDB (complex), splitting into
isolated protein and peptide to get ΔASA/BSA per the paper's equation 1:
    BSA = (ASA_protein + ASA_peptide - ASA_complex) / 2
NOTE: FreeSASA != NACCESS numerically (different algorithm/radii); this is a
reproducible re-implementation, expected to correlate closely, not match exactly.
"""
import csv
import os
import sys
import tempfile

import freesasa

# quiet FreeSASA's stderr chatter
freesasa.setVerbosity(freesasa.silent)

# Match NACCESS as closely as FreeSASA allows: NACCESS radii + Lee-Richards
# (NACCESS's algorithm) at high slice count. Reproducible, still pip-only.
_NACCESS = freesasa.Classifier.getStandardClassifier("naccess")
_PARAMS = freesasa.Parameters({"algorithm": freesasa.LeeRichards,
                               "probe-radius": 1.4, "n-slices": 100})

def asa_of(lines):
    """Total ASA of a set of PDB ATOM lines written to a temp file."""
    tmp = tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False)
    try:
        tmp.writelines(l if l.endswith("\n") else l + "\n" for l in lines)
        tmp.write("END\n")
        tmp.close()
        s = freesasa.Structure(tmp.name, _NACCESS)
        return freesasa.calc(s, _PARAMS).totalArea()
    finally:
        os.unlink(tmp.name)

def areas_of(lines):
    """(total_area, {chain_id: area}) from ONE FreeSASA calc (NACCESS-matched)."""
    tmp = tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False)
    try:
        tmp.writelines(l if l.endswith("\n") else l + "\n" for l in lines)
        tmp.write("END\n")
        tmp.close()
        result = freesasa.calc(freesasa.Structure(tmp.name, _NACCESS), _PARAMS)
        per = {ch: sum(r.total for r in resd.values())
               for ch, resd in result.residueAreas().items()}
        return result.totalArea(), per
    finally:
        os.unlink(tmp.name)

def process(pdb_path, pep_chain, prot_chain):
    with open(pdb_path) as fh:
        atom_lines = [l for l in fh if l.startswith(("ATOM", "HETATM"))]
    pep = [l for l in atom_lines if len(l) > 21 and l[21] == pep_chain]
    prot = [l for l in atom_lines if len(l) > 21 and l[21] == prot_chain]
    if not pep or not prot:
        return None
    asa_complex, per_chain = areas_of(atom_lines)
    asa_pep = asa_of(pep)
    asa_prot = asa_of(prot)
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
            path = os.path.join(p.pdb_dir, f"{eid}.pdb")
            if not os.path.exists(path):
                continue
            n += 1
            r = process(path, row["pep_chain"], row["prot_chain"])
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