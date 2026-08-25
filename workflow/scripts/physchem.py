"""ProtParam physicochemical properties for extracted pairs, matching v15 columns.
Non-standard residues (X etc.) are stripped before analysis (ProtParam requires
the 20 canonical amino acids)."""
import sys
from Bio.SeqUtils.ProtParam import ProteinAnalysis

STD = set("ACDEFGHIKLMNPQRSTVWY")

# Residue atomic composition = free amino acid MINUS one water (peptide-bond
# condensation). v15's Formula/TotalAtoms are the SUM of residue compositions
# with NO added terminal water (verified exactly against v15 for peptide+protein).
RES = {
    'G': (2, 3, 1, 1, 0), 'A': (3, 5, 1, 1, 0), 'S': (3, 5, 1, 2, 0),
    'P': (5, 7, 1, 1, 0), 'V': (5, 9, 1, 1, 0), 'T': (4, 7, 1, 2, 0),
    'C': (3, 5, 1, 1, 1), 'L': (6, 11, 1, 1, 0), 'I': (6, 11, 1, 1, 0),
    'N': (4, 6, 2, 2, 0), 'D': (4, 5, 1, 3, 0), 'Q': (5, 8, 2, 2, 0),
    'K': (6, 12, 2, 1, 0), 'E': (5, 7, 1, 3, 0), 'M': (5, 9, 1, 1, 1),
    'H': (6, 7, 3, 1, 0), 'F': (9, 9, 1, 1, 0), 'R': (6, 12, 4, 1, 0),
    'Y': (9, 9, 1, 2, 0), 'W': (11, 10, 2, 1, 0),
}  # (C, H, N, O, S)


def clean(seq):
    return "".join(c for c in seq.upper() if c in STD)


def formula_atoms(seq):
    """Return (formula_string, total_atom_count) for a cleaned sequence."""
    C = H = N = O = S = 0
    for aa in seq:
        c, h, n, o, s = RES[aa]
        C += c; H += h; N += n; O += o; S += s
    formula = f"C{C}H{H}N{N}O{O}" + (f"S{S}" if S else "")
    return formula, C + H + N + O + S


def aliphatic_index(seq):
    # Ikai (1980): X_Ala + 2.9*X_Val + 3.9*(X_Ile + X_Leu), mole% *100
    n = len(seq)
    if n == 0:
        return 0.0
    a = seq.count("A") / n * 100
    v = seq.count("V") / n * 100
    i = seq.count("I") / n * 100
    l = seq.count("L") / n * 100
    return round(a + 2.9 * v + 3.9 * (i + l), 3)


def hydrophobic_percent(seq):
    # percent of residues that are hydrophobic (A,V,L,I,P,F,M,W)
    n = len(seq)
    if n == 0:
        return 0.0
    h = sum(seq.count(x) for x in "AVLIPFMW")
    return round(h / n * 100, 3)


def props(seq):
    c = clean(seq)
    if len(c) == 0:
        return None
    pa = ProteinAnalysis(c)
    pos = sum(c.count(x) for x in "KR")   # simple charged-residue tallies
    neg = sum(c.count(x) for x in "DE")
    ext_reduced, ext_cystine = pa.molar_extinction_coefficient()  # (NoDisulfide, Disulfide)
    formula, atoms = formula_atoms(c)
    return {
        "clean_len": len(c),
        "stripped": len(seq) - len(c),
        "MW": round(pa.molecular_weight(), 3),
        "pI": round(pa.isoelectric_point(), 3),
        "GRAVY": round(pa.gravy(), 3),
        "Instability": round(pa.instability_index(), 3),
        "Aliphatic": aliphatic_index(c),
        "Hydrophobic": hydrophobic_percent(c),
        "Pos": pos,
        "Neg": neg,
        "ExtCoeff_Disulfide": ext_cystine,
        "ExtCoeff_NoDisulfide": ext_reduced,
        "Formula": formula,
        "TotalAtoms": atoms,
    }


def main():
    rows = list(open(snakemake.input.pairs))              # noqa: F821
    header = rows[0].rstrip("\n").split("\t")
    idx = {h: i for i, h in enumerate(header)}
    out_cols = ["id", "chain_type", "clean_len", "stripped", "MW", "pI",
                "GRAVY", "Instability", "Aliphatic", "Hydrophobic", "Pos", "Neg",
                "ExtCoeff_Disulfide", "ExtCoeff_NoDisulfide", "Formula", "TotalAtoms"]
    with open(snakemake.output.physchem, "w") as fh:      # noqa: F821
        fh.write("\t".join(out_cols) + "\n")
        for line in rows[1:]:
            f = line.rstrip("\n").split("\t")
            eid = f[idx["id"]]
            for ctype, seqcol in (("peptide", "pep_seq"), ("protein", "prot_seq")):
                p = props(f[idx[seqcol]])
                if not p:
                    continue
                p.update(id=eid, chain_type=ctype)
                fh.write("\t".join(str(p[c]) for c in out_cols) + "\n")
    print("physchem written", file=sys.stderr)


if __name__ == "__main__":
    main()