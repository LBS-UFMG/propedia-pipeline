"""PDB metadata columns from the CIF header, matching v15:
CLASSIFICATION, DEPOSITION_DATE, RESOLUTION, STRUCTURE_METHOD, TITLE, organism,
and per-chain PEPTIDE_DESC / PROTEIN_DESC (entity description).

organism reproduces v15's single per-pair value (the PROTEIN chain's entity
source; v15's rule is the receptor's organism, not the first entity). We ALSO
emit PEPTIDE_ORGANISM and PROTEIN_ORGANISM, resolving each chain's organism
separately — information v15 collapsed to one value (see reproduction_notes).
A PDB may appear in many pairs, so the CIF is parsed once per pdb id and cached.
"""
import csv
import gzip
import os
import sys

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

import checkpoint

VERSION = "2"            # bump on logic change -> invalidates old checkpoints


def _first(d, key):
    v = d.get(key)
    if v is None:
        return ""
    if isinstance(v, list):
        for x in v:
            if x not in ("", "?", "."):
                return x
        return ""
    return v if v not in ("?", ".") else ""


def _resolution(d):
    for key in ("_refine.ls_d_res_high",
                "_reflns.d_resolution_high",
                "_em_3d_reconstruction.resolution"):
        raw = _first(d, key)
        if raw:
            try:
                return str(float(raw))       # 2.30 -> "2.3"
            except ValueError:
                return raw
    return ""


def _aslist(x):
    if x is None:
        return []
    return x if isinstance(x, list) else [x]


def _chain_entity_map(d):
    """auth chain id -> entity_id, via _entity_poly.pdbx_strand_id."""
    ents = _aslist(d.get("_entity_poly.entity_id"))
    strands = _aslist(d.get("_entity_poly.pdbx_strand_id"))
    out = {}
    for eid, s in zip(ents, strands):
        for ch in str(s).split(","):
            ch = ch.strip()
            if ch:
                out[ch] = eid
    return out


def _entity_desc_map(d):
    ids = _aslist(d.get("_entity.id"))
    descs = _aslist(d.get("_entity.pdbx_description"))
    return {i: ("" if x in ("?", ".") else x) for i, x in zip(ids, descs)}


def _entity_org_map(d):
    """entity_id -> source organism (scientific), first available category."""
    out = {}
    for eid_key, name_key in (
            ("_entity_src_gen.entity_id", "_entity_src_gen.pdbx_gene_src_scientific_name"),
            ("_entity_src_nat.entity_id", "_entity_src_nat.pdbx_organism_scientific"),
            ("_pdbx_entity_src_syn.entity_id", "_pdbx_entity_src_syn.organism_scientific")):
        for eid, name in zip(_aslist(d.get(eid_key)), _aslist(d.get(name_key))):
            if name not in ("", "?", ".") and eid not in out:
                out[eid] = name
    return out


def _global_organism(d):
    for key in ("_entity_src_gen.pdbx_gene_src_scientific_name",
                "_entity_src_nat.pdbx_organism_scientific",
                "_pdbx_entity_src_syn.organism_scientific"):
        v = _first(d, key)
        if v:
            return v
    return ""


def parse_cif(path):
    d = MMCIF2Dict(gzip.open(path, "rt"))
    return {
        "CLASSIFICATION": _first(d, "_struct_keywords.pdbx_keywords"),
        "TITLE": _first(d, "_struct.title"),
        "STRUCTURE_METHOD": _first(d, "_exptl.method"),
        "DEPOSITION_DATE": _first(d, "_pdbx_database_status.recvd_initial_deposition_date"),
        "RESOLUTION": _resolution(d),
        "_chain_ent": _chain_entity_map(d),
        "_ent_desc": _entity_desc_map(d),
        "_ent_org": _entity_org_map(d),
        "_global_org": _global_organism(d),
    }


def organism_for(meta, chain):
    """Organism of a chain's entity; global source organism as fallback."""
    eid = meta["_chain_ent"].get(chain)
    org = meta["_ent_org"].get(eid, "") if eid else ""
    return (org or meta["_global_org"]).upper()


def desc_for(meta, chain):
    eid = meta["_chain_ent"].get(chain)
    return meta["_ent_desc"].get(eid, "") if eid else ""


def shard_path(cif_dir, pid):
    sub = pid[1:3].lower() if len(pid) >= 3 else "_"
    return os.path.join(cif_dir, sub, f"{pid}.cif.gz")


OUT_COLS = ["id", "CLASSIFICATION", "DEPOSITION_DATE", "RESOLUTION",
            "STRUCTURE_METHOD", "TITLE", "organism",
            "PEPTIDE_ORGANISM", "PROTEIN_ORGANISM",
            "PEPTIDE_DESC", "PROTEIN_DESC"]


def worker(item):
    """worker(item) -> (meta_dict|None, status). One CIF header parse per PDB id.
    A parse failure is deterministic, so it is checkpointed (not retried)."""
    pid, cif_path = item
    try:
        return parse_cif(cif_path), "ok"
    except Exception as exc:                               # noqa: BLE001
        return None, f"parse_error:{exc}"


def main():
    p = snakemake.params                                   # noqa: F821
    threads = getattr(snakemake, "threads", 1) or 1        # noqa: F821
    pairs = list(csv.DictReader(open(snakemake.input.pairs), delimiter="\t"))  # noqa: F821

    pids = sorted({row["pdb"] for row in pairs})
    items = [(pid, shard_path(p.cif_dir, pid)) for pid in pids]
    workdir = checkpoint.namespace(p.ckpt, VERSION, {})
    results = checkpoint.run(items, worker, workdir, threads=threads,
                             id_of=lambda it: it[0], stage="metadata")

    n = 0
    with open(snakemake.output.metadata, "w") as fh:       # noqa: F821
        fh.write("\t".join(OUT_COLS) + "\n")
        for row in pairs:
            res = results.get(row["pdb"])
            if not res or res["status"] != "ok" or not res["record"]:
                continue
            meta = res["record"]
            pep, prot = row["pep_chain"], row["prot_chain"]
            rec = {
                "id": row["id"],
                "CLASSIFICATION": meta["CLASSIFICATION"],
                "DEPOSITION_DATE": meta["DEPOSITION_DATE"],
                "RESOLUTION": meta["RESOLUTION"],
                "STRUCTURE_METHOD": meta["STRUCTURE_METHOD"],
                "TITLE": meta["TITLE"],
                "organism": organism_for(meta, prot),          # v15-compatible (protein's)
                "PEPTIDE_ORGANISM": organism_for(meta, pep),
                "PROTEIN_ORGANISM": organism_for(meta, prot),
                "PEPTIDE_DESC": desc_for(meta, pep),
                "PROTEIN_DESC": desc_for(meta, prot),
            }
            fh.write("\t".join(str(rec[c]).replace("\t", " ").replace("\n", " ")
                               for c in OUT_COLS) + "\n")
            n += 1
    print(f"DONE {n} metadata rows", file=sys.stderr)


if __name__ == "__main__":
    main()