"""Assemble the final Propedia CSV (;-delimited) from the per-feature TSVs.

Master row set = extracted pairs that pass the BSA>0 filter (i.e. ids present in
surface.tsv), matching v15's inclusion criteria. Left-join every feature onto it;
missing cells are blank. Columns are emitted in the exact v15 order, followed by
the two additive organism columns.
"""
import csv
import sys

DEG = "\u02da"   # U+02DA RING ABOVE, as used in the v15 header

# ---- exact v15 column order (71) ----
V15_ORDER = [
    "id", "AAP", "ABP", "ACP", "AIP",
    "ASA_Complex", "ASA_Peptide", "ASA_Protein", "BPP%", "BPepA", "BProA", "BSA",
    "CLASSIFICATION", "DEPOSITION_DATE", "Interface Residues",
    "No. of apolar-apolar contacts", "No. of apolar-polar contacts",
    "No. of charged-apolar contacts", "No. of charged-charged contacts",
    "No. of charged-polar contacts", "No. of intermolecular contacts",
    "No. of polar-polar contacts",
    "PDB_ID", "PEPTIDE_CHAIN", "PEPTIDE_DESC", "PEPTIDE_SEQ", "PEPTIDE_SIZE",
    "PROTEIN_CHAIN", "PROTEIN_DESC", "PROTEIN_SEQ", "PROTEIN_SIZE",
    "Percentage of apolar NIS residues", "Percentage of charged NIS residues",
    "Predicted binding affinity (kcal.mol-1)",
    f"Predicted dissociation constant (M) at 25.0{DEG}C",
    "QSP", "RESOLUTION", "SBP", "STRUCTURE_METHOD", "TITLE",
    "binding-cluster", "interface-cluster", "is_leader", "leader_id", "organism",
    "peptide_AliphaticIndex", "peptide_ExtCoeff_Disulfide",
    "peptide_ExtCoeff_NoDisulfide", "peptide_Formula", "peptide_GRAVY",
    "peptide_HydrophobicPercent", "peptide_InstabilityIndex", "peptide_MW",
    "peptide_NegativeResidues", "peptide_PositiveResidues", "peptide_TotalAtoms",
    "peptide_pI",
    "protein_AliphaticIndex", "protein_ExtCoeff_Disulfide",
    "protein_ExtCoeff_NoDisulfide", "protein_Formula", "protein_GRAVY",
    "protein_HydrophobicPercent", "protein_InstabilityIndex", "protein_MW",
    "protein_NegativeResidues", "protein_PositiveResidues", "protein_TotalAtoms",
    "protein_pI",
    "seq100_clusters", "sequence-cluster",
]
PISA_COLS = ["pisa_status", "pisa_interface_class",
             "pisa_assembly_done", "pisa_n_interfaces",
             "pisa_interface_id", "pisa_chain_1", "pisa_chain_2", "pisa_css",
             "pisa_area", "pisa_solv_en", "pisa_pvalue", "pisa_tipo",
             "pisa_n_hbonds", "pisa_n_saltbridges"]
ADDITIVE = ["PEPTIDE_ORGANISM", "PROTEIN_ORGANISM", "FIRST_RELEASE"] + PISA_COLS

# physchem long-field -> v15 suffix
PHYS_SUFFIX = {
    "Aliphatic": "AliphaticIndex", "ExtCoeff_Disulfide": "ExtCoeff_Disulfide",
    "ExtCoeff_NoDisulfide": "ExtCoeff_NoDisulfide", "Formula": "Formula",
    "GRAVY": "GRAVY", "Hydrophobic": "HydrophobicPercent",
    "Instability": "InstabilityIndex", "MW": "MW", "Neg": "NegativeResidues",
    "Pos": "PositiveResidues", "TotalAtoms": "TotalAtoms", "pI": "pI",
}
SUFFIX_PHYS = {f"{c}_{v}": (c, k) for c in ("peptide", "protein")
               for k, v in PHYS_SUFFIX.items()}   # e.g. "peptide_MW" -> ("peptide","MW")

PAIRS_MAP = {"PDB_ID": "pdb", "PEPTIDE_CHAIN": "pep_chain", "PEPTIDE_SEQ": "pep_seq",
             "PEPTIDE_SIZE": "pep_size", "PROTEIN_CHAIN": "prot_chain",
             "PROTEIN_SEQ": "prot_seq", "PROTEIN_SIZE": "prot_size"}
SURFACE_COLS = {"ASA_Complex", "ASA_Peptide", "ASA_Protein", "BPP%", "BPepA", "BProA", "BSA"}
PRODIGY_MAP = {
    "No. of apolar-apolar contacts": "n_apolar_apolar",
    "No. of apolar-polar contacts": "n_apolar_polar",
    "No. of charged-apolar contacts": "n_charged_apolar",
    "No. of charged-charged contacts": "n_charged_charged",
    "No. of charged-polar contacts": "n_charged_polar",
    "No. of intermolecular contacts": "n_intermolecular",
    "No. of polar-polar contacts": "n_polar_polar",
    "Percentage of apolar NIS residues": "nis_apolar",
    "Percentage of charged NIS residues": "nis_charged",
    "Predicted binding affinity (kcal.mol-1)": "dg",
    f"Predicted dissociation constant (M) at 25.0{DEG}C": "kd",
}
META_COLS = {"CLASSIFICATION", "DEPOSITION_DATE", "RESOLUTION", "STRUCTURE_METHOD",
             "TITLE", "organism", "PEPTIDE_DESC", "PROTEIN_DESC"}
THERAP_COLS = {"AAP", "ABP", "ACP", "AIP", "QSP", "SBP"}
LEGACY_COLS = {"sequence-cluster", "interface-cluster", "binding-cluster",
               "is_leader", "leader_id"}
BLANK_COLS = set()   # every v15 column is now sourced


def load_tsv(path):
    with open(path) as fh:
        return {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def load_physchem(path):
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out.setdefault(r["id"], {})[r["chain_type"]] = r
    return out


def clean(v):
    return str(v).replace(";", ",").replace("\t", " ").replace("\n", " ").replace("\r", "")


def value(col, eid, pair, surface, prodigy, interface, metadata, cnr,
          physchem, therapeutic, legacy, provenance, pisa):
    if col == "id":
        return eid
    if col == "FIRST_RELEASE":
        return (provenance.get(eid) or {}).get("FIRST_RELEASE", "")
    if col.startswith("pisa_"):
        return (pisa.get(eid) or {}).get(col, "")
    if col in BLANK_COLS:
        return ""
    if col in THERAP_COLS:
        return (therapeutic.get(eid) or {}).get(col, "")
    if col in LEGACY_COLS:
        return (legacy.get(eid) or {}).get(col, "")
    if col in PAIRS_MAP:
        return pair.get(PAIRS_MAP[col], "")
    if col in SURFACE_COLS:
        return (surface.get(eid) or {}).get(col, "")
    if col in PRODIGY_MAP:
        return (prodigy.get(eid) or {}).get(PRODIGY_MAP[col], "")
    if col == "Interface Residues":
        return (interface.get(eid) or {}).get("interface_residues", "")
    if col in META_COLS or col in ADDITIVE:
        return (metadata.get(eid) or {}).get(col, "")
    if col == "seq100_clusters":
        return (cnr.get(eid) or {}).get("cnr_cluster", "")
    if col in SUFFIX_PHYS:
        ctype, field = SUFFIX_PHYS[col]
        return (physchem.get(eid, {}).get(ctype) or {}).get(field, "")
    return ""


def main():
    inp = snakemake.input                                  # noqa: F821
    pairs = list(csv.DictReader(open(inp.pairs), delimiter="\t"))
    surface = load_tsv(inp.surface)
    prodigy = load_tsv(inp.prodigy)
    interface = load_tsv(inp.interface)
    metadata = load_tsv(inp.metadata)
    cnr = load_tsv(inp.clusters)
    physchem = load_physchem(inp.physchem)
    therapeutic = load_tsv(inp.therapeutic)
    legacy = load_tsv(inp.legacy)
    provenance = load_tsv(inp.provenance)
    pisa = load_tsv(inp.pisa)

    # master = extracted pairs that passed BSA>0 (present in surface), sorted by id
    kept = [r for r in pairs if r["id"] in surface]
    kept.sort(key=lambda r: r["id"])

    cols = V15_ORDER + ADDITIVE
    with open(snakemake.output.csv, "w") as fh:            # noqa: F821
        fh.write(";".join(cols) + "\n")
        for r in kept:
            eid = r["id"]
            row = [clean(value(c, eid, r, surface, prodigy, interface,
                               metadata, cnr, physchem, therapeutic, legacy,
                               provenance, pisa))
                   for c in cols]
            fh.write(";".join(row) + "\n")
    print(f"DONE assembled {len(kept)} rows x {len(cols)} cols", file=sys.stderr)


if __name__ == "__main__":
    main()