"""ProtParam physicochemical properties for extracted pairs, matching v15 columns.
Non-standard residues (X etc.) are stripped before analysis (ProtParam requires
the 20 canonical amino acids)."""
import sys
from Bio.SeqUtils.ProtParam import ProteinAnalysis

STD = set("ACDEFGHIKLMNPQRSTVWY")


def clean(seq):
    return "".join(c for c in seq.upper() if c in STD)


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
    }


def main():
    rows = list(open(snakemake.input.pairs))              # noqa: F821
    header = rows[0].rstrip("\n").split("\t")
    idx = {h: i for i, h in enumerate(header)}
    out_cols = ["id", "chain_type", "clean_len", "stripped", "MW", "pI",
                "GRAVY", "Instability", "Aliphatic", "Hydrophobic", "Pos", "Neg"]
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