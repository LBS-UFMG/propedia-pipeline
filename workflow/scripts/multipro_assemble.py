"""Multipro Phase 2d: assemble the multipro CSV (exact v4 order + additive columns) by
joining Phase 1 grouping + Phase 2b surface + Phase 2c PRODIGY with the pep-pro
outputs (physchem/metadata/interface/legacy/pisa), joined via `items`. Single-value
columns come from items[0]; per-protein-chain columns are colon-joined in items
order. One source of truth -- nothing is recomputed here.
"""
import csv
import sys

DEG = "\u02da"

V4_ORDER = [
    "cluster_id", "items", "PDB_ID", "TITLE", "RESOLUTION", "CLASSIFICATION",
    "DEPOSITION_DATE", "STRUCTURE_METHOD", "PROTEIN_CHAIN", "PEPTIDE_CHAIN",
    "PROTEIN_SIZE", "PEPTIDE_SIZE", "PROTEIN_DESC", "PEPTIDE_DESC",
    "PROTEIN_SEQ", "PEPTIDE_SEQ", "leader_id",
    "peptide_Length", "peptide_MW", "peptide_pI", "peptide_InstabilityIndex",
    "peptide_AliphaticIndex", "peptide_GRAVY", "peptide_HydrophobicPercent",
    "peptide_PositiveResidues", "peptide_NegativeResidues", "peptide_Formula",
    "peptide_TotalAtoms", "peptide_ExtCoeff_Disulfide", "peptide_ExtCoeff_NoDisulfide",
    "protein_Length", "protein_MW", "protein_pI", "protein_InstabilityIndex",
    "protein_AliphaticIndex", "protein_GRAVY", "protein_HydrophobicPercent",
    "protein_PositiveResidues", "protein_NegativeResidues", "protein_Formula",
    "protein_TotalAtoms", "protein_ExtCoeff_Disulfide", "protein_ExtCoeff_NoDisulfide",
    "No. of intermolecular contacts", "No. of charged-charged contacts",
    "No. of charged-polar contacts", "No. of charged-apolar contacts",
    "No. of polar-polar contacts", "No. of apolar-polar contacts",
    "No. of apolar-apolar contacts", "Percentage of apolar NIS residues",
    "Percentage of charged NIS residues", "Predicted binding affinity (kcal.mol-1)",
    f"Predicted dissociation constant (M) at 25.0{DEG}C", "Interface Residues",
    "count", "cluster_id", "ASA_Complex", "ASA_Protein", "ASA_Peptide",
    "BProA", "BPepA", "BPP%", "BSA",
]
# Appended after the sacred v4 order. The PISA columns are the multipro read-out of the
# biological-vs-crystal annotation (crystal-packing reviewer comment): a peptide in a
# multipro assembly touches several protein chains, so we report the per-chain PISA class
# and CSS (colon-joined in items order, parity with PROTEIN_CHAIN) plus ONE aggregate call
# for the entry. Aggregated from the pep-pro pisa.tsv via `items` -- no PISA re-run.
ADDITIVE = ["FIRST_RELEASE", "pisa_interface_class",
            "pisa_interface_class_per_chain", "pisa_css_per_chain"]


def aggregate_pisa_class(classes):
    """One biological/crystal call for a multipro entry from its per-protein-chain classes.
    A peptide that forms >=1 biological interface within the assembly is biological; if all
    scored interfaces are crystal-packing it is crystal-packing; non-X-ray => not_applicable;
    otherwise indeterminate; blank only when PISA could not be assessed for any constituent."""
    cs = list(classes)
    if any(c == "biological" for c in cs):
        return "biological"
    if cs and all(c == "not_applicable" for c in cs):
        return "not_applicable"
    nonblank = [c for c in cs if c and c != "not_applicable"]
    if not nonblank:
        return "not_applicable" if any(c == "not_applicable" for c in cs) else ""
    if all(c == "crystal-packing" for c in nonblank):
        return "crystal-packing"
    return "indeterminate"          # a mix incl. indeterminate, none biological

# multipro physchem suffix -> pep-pro physchem long-field
PHYS = {"MW": "MW", "pI": "pI", "InstabilityIndex": "Instability",
        "AliphaticIndex": "Aliphatic", "GRAVY": "GRAVY",
        "HydrophobicPercent": "Hydrophobic", "PositiveResidues": "Pos",
        "NegativeResidues": "Neg", "Formula": "Formula", "TotalAtoms": "TotalAtoms",
        "ExtCoeff_Disulfide": "ExtCoeff_Disulfide",
        "ExtCoeff_NoDisulfide": "ExtCoeff_NoDisulfide"}
PROD_COLS = {
    "No. of intermolecular contacts", "No. of charged-charged contacts",
    "No. of charged-polar contacts", "No. of charged-apolar contacts",
    "No. of polar-polar contacts", "No. of apolar-polar contacts",
    "No. of apolar-apolar contacts", "Percentage of apolar NIS residues",
    "Percentage of charged NIS residues", "Predicted binding affinity (kcal.mol-1)",
    f"Predicted dissociation constant (M) at 25.0{DEG}C",
}
SURF_COLS = {"ASA_Complex", "ASA_Protein", "ASA_Peptide", "BProA", "BPepA", "BPP%", "BSA"}
META_SINGLE = {"TITLE", "RESOLUTION", "CLASSIFICATION", "DEPOSITION_DATE",
               "STRUCTURE_METHOD", "PEPTIDE_DESC"}


def load_tsv(path, key="id"):
    with open(path) as fh:
        return {r[key]: r for r in csv.DictReader(fh, delimiter="\t")}  


def load_physchem(path):
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out.setdefault(r["id"], {})[r["chain_type"]] = r
    return out


def clean(v):
    return str(v).replace(";", ",").replace("\t", " ").replace("\n", " ").replace("\r", "")


def colon(vals):
    return ":".join(str(v) for v in vals)


def main():
    inp = snakemake.input                                  # noqa: F821
    mp = list(csv.DictReader(open(inp.multipro), delimiter="\t"))
    surface = load_tsv(inp.surface, key="cluster_id")
    prodigy = load_tsv(inp.prodigy, key="cluster_id")
    metadata = load_tsv(inp.metadata)
    interface = load_tsv(inp.interface)
    legacy = load_tsv(inp.legacy)
    physchem = load_physchem(inp.physchem)
    provenance = load_tsv(inp.provenance)
    pisa = load_tsv(inp.pisa)          # pep-pro PISA, keyed by pair id (== multipro item)

    def phys(eid, ctype, field):
        return (physchem.get(eid, {}).get(ctype) or {}).get(field, "")

    with open(snakemake.output.csv, "w") as fh:            # noqa: F821
        fh.write(";".join(V4_ORDER + ADDITIVE) + "\n")
        for r in mp:
            cid = r["cluster_id"]
            items = r["items"].split(":")
            i0 = items[0]
            surf = surface.get(cid, {})
            prod = prodigy.get(cid, {})
            row = []
            for col in V4_ORDER:
                if col in ("cluster_id",):
                    v = cid
                elif col in ("items", "PDB_ID", "PROTEIN_CHAIN", "PEPTIDE_CHAIN",
                             "PROTEIN_SIZE", "PEPTIDE_SIZE", "PROTEIN_SEQ",
                             "PEPTIDE_SEQ", "count"):
                    v = r.get(col, "")
                elif col == "peptide_Length":
                    v = r.get("PEPTIDE_SIZE", "")
                elif col == "protein_Length":
                    v = r.get("PROTEIN_SIZE", "")
                elif col == "leader_id":
                    v = (legacy.get(i0) or {}).get("leader_id", "")
                elif col in META_SINGLE:
                    v = (metadata.get(i0) or {}).get(col, "")
                elif col == "PROTEIN_DESC":
                    v = colon((metadata.get(i) or {}).get("PROTEIN_DESC", "") for i in items)
                elif col == "Interface Residues":
                    v = colon((interface.get(i) or {}).get("interface_residues", "") for i in items)
                elif col.startswith("peptide_"):
                    v = phys(i0, "peptide", PHYS[col[len("peptide_"):]])
                elif col.startswith("protein_"):
                    field = PHYS[col[len("protein_"):]]
                    v = colon(phys(i, "protein", field) for i in items)
                elif col in PROD_COLS:
                    v = prod.get(col, "")
                elif col in SURF_COLS:
                    v = surf.get(col, "")
                else:
                    v = ""
                row.append(clean(v))
            # additive: FIRST_RELEASE of the entry (from items[0]'s pep-pro provenance)
            row.append(clean((provenance.get(i0) or {}).get("FIRST_RELEASE", "")))
            # additive: PISA biological-vs-crystal read-out, aggregated from the
            # per-protein-chain pep-pro interfaces (colon-joined in items order).
            per_chain_class = [(pisa.get(i) or {}).get("pisa_interface_class", "") for i in items]
            per_chain_css = [(pisa.get(i) or {}).get("pisa_css", "") for i in items]
            row.append(clean(aggregate_pisa_class(per_chain_class)))
            row.append(clean(colon(per_chain_class)))
            row.append(clean(colon(per_chain_css)))
            fh.write(";".join(row) + "\n")
    print(f"DONE assembled {len(mp)} multipro rows x {len(V4_ORDER + ADDITIVE)} cols",
          file=sys.stderr)


if __name__ == "__main__":
    main()