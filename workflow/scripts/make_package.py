"""Package the pipeline outputs into the propedia26 website file tree.

Emits the per-entry layout the web app serves (see docs/propedia26_publishing.md):
  <web>/data/<mode>/csv/<id[0]>/<id>.csv            per-entry row (;-delimited, no header)
  <web>/data/<mode>/contacts/<id>/<PDB>_contacts.csv per-entry COCaDA contacts
  <web>/data/<mode>/multipro/csv/<id[0]>/<id>.csv    per-entry multipro row
  <web>/data/<mode>/pdb/<id[0]>/<id>.cif             per-entry complex structure (mmCIF)
  <web>/data/<mode>/multipro/pdb/<id[0]>/<id>.cif    per-entry multipro complex (mmCIF)
  <web>/data/clusters/*.tsv                          cluster tables (copied)
  <web>/columns_peppro.txt / columns_multipro.txt    the column ORDER (manifest)

COLUMN ORDER: we emit our own documented order (from assemble.py). The website's
final order is a pure PERMUTATION of these named columns, applied later with
reorder_columns.py once the target order is known (a live v17 file or a defined
site schema) — so this step does not depend on that decision.
"""
import csv
import os
import shutil
import sys


def _clear_make(d):
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.makedirs(d, exist_ok=True)


def _load_order(path):
    """Target column order = one column name per line (blank lines / # ignored)."""
    if not path or not os.path.exists(os.path.expanduser(path)):
        return None
    with open(os.path.expanduser(path)) as fh:
        return [l.strip() for l in fh if l.strip() and not l.startswith("#")]


def split_per_entry(csv_path, out_dir, id_col="id", target=None):
    """Split a ;-delimited master CSV into per-entry files sharded by id[0], header
    dropped (the site reads by position). If `target` (a list of column names) is
    given, each row is PERMUTED to that order (missing cols -> blank) — this is the
    one hook needed to match the site's final layout later. Returns (n, header_out).

    Target->source column matching is case-INSENSITIVE, so the site's v17 PISA columns
    (`PISA_status`, `PISA_area`, …) resolve to our lowercase `pisa_*` columns. Target
    columns with no source match (e.g. the richer v17 PISA fields we don't compute:
    PISA_nres_1, PISA_diss_energy, …) fall through to blank — intentionally, for now."""
    _clear_make(out_dir)
    with open(csv_path) as fh:
        src = fh.readline().rstrip("\n").split(";")
        pos = {c: i for i, c in enumerate(src)}
        pos_ci = {c.lower(): i for i, c in enumerate(src)}   # case-insensitive fallback
        idx = pos[id_col]
        header_out = target if target else src

        def col_idx(c):
            i = pos.get(c)
            return i if i is not None else pos_ci.get(c.lower())

        n = 0
        for line in fh:
            f = line.rstrip("\n").split(";")
            eid = f[idx]
            if target:
                out = []
                for c in header_out:
                    i = col_idx(c)
                    out.append(f[i] if i is not None and i < len(f) else "")
            else:
                out = f
            shard = os.path.join(out_dir, eid[0])
            os.makedirs(shard, exist_ok=True)
            with open(os.path.join(shard, f"{eid}.csv"), "w") as ofh:
                ofh.write(";".join(out) + "\n")
            n += 1
    return n, header_out


# Explore-page summary TSV (data/propedia26_v17.tsv): TAB-delimited, NO header, one row
# per pep-pro entry. The web app (Home.php EXPLORE_ARQUIVO / Entry.php getPisaCss) reads it
# BY FIXED COLUMN INDEX: [0-7] are the displayed columns, [8-24] drive the Explore filters,
# and [22] is the PISA CSS the Entry page shows. Each item is the source column in
# propedia.csv to project into that TAB position (case-insensitive, so pisa_* -> PISA_*).
# Indices/semantics are taken verbatim from Home.php's filter code and EXPLORE_CLASSES.
EXPLORE_SPEC = [
    "id",                    # 0  displayed / entry key (Home.php EXPLORE_COLUNAS_EXIBIDAS=8)
    "PROTEIN_SIZE",          # 1  displayed
    "PEPTIDE_SIZE",          # 2  displayed + size filter (c[2])
    "PEPTIDE_SEQ",           # 3  displayed + canonical-only filter (c[3], tests for 'x')
    "TITLE",                 # 4  displayed
    "CLASSIFICATION",        # 5  displayed + classification filter (c[5])
    "is_leader",             # 6  displayed + redundancy filter (c[6] == 'yes')
    "leader_id",             # 7  displayed
    "pisa_n_hbonds",         # 8  min H-bonds filter (c[8])          [CONFIRM: PISA vs contacts]
    "pisa_n_saltbridges",    # 9  salt-bridge presence filter (c[9]) [CONFIRM: PISA vs contacts]
    "BSA",                   # 10 min BSA filter (c[10])
    "BPP%",                  # 11 min buried-peptide% filter (c[11])
    "RESOLUTION",            # 12 max resolution filter (c[12])
    "STRUCTURE_METHOD",      # 13 method filter (c[13])
    "peptide_HydrophobicPercent",  # 14 min hydrophobic% filter (c[14])
    "peptide_PositiveResidues",    # 15 min positive-residues filter (c[15])
    "AAP",                   # 16 therapeutic score (EXPLORE_CLASSES)
    "ABP",                   # 17
    "ACP",                   # 18
    "AIP",                   # 19
    "QSP",                   # 20
    "SBP",                   # 21
    "pisa_css",              # 22 PISA CSS (Entry getPisaCss reads exactly col 22)
    "Predicted binding affinity (kcal.mol-1)",              # 23 max binding-affinity filter (c[23])
    "Predicted dissociation constant (M) at 25.0˚C",   # 24 min dissociation filter (c[24])
]


def build_explore_tsv(csv_path, out_path):
    """Project propedia.csv into the Explore summary TSV (TAB, no header) in EXPLORE_SPEC
    order. Returns the row count. Source columns are matched case-insensitively; any spec
    column absent from propedia.csv is written blank."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(csv_path) as fh:
        src = fh.readline().rstrip("\n").split(";")
        pos = {c: i for i, c in enumerate(src)}
        pos_ci = {c.lower(): i for i, c in enumerate(src)}
        idxs = [pos.get(c, pos_ci.get(c.lower())) for c in EXPLORE_SPEC]
        n = 0
        with open(out_path, "w") as ofh:
            for line in fh:
                f = line.rstrip("\n").split(";")
                row = [(f[i] if i is not None and i < len(f) else "") for i in idxs]
                ofh.write("\t".join(row) + "\n")
                n += 1
    return n


def build_blast_fastas(csv_path, pep_path, rec_path):
    """Write the two BLAST subject FASTAs the site's Blast.php searches
    (data/peptides.fasta, data/receptors.fasta) from propedia.csv. Header format
    matches the shipped files: `>{id} | {DESC}` — BLAST takes {id} as sseqid (up to the
    space) and the blast view links /entry/{id} via explode('|')[0]. One record per
    pep-pro entry. Returns (n_peptides, n_receptors)."""
    os.makedirs(os.path.dirname(pep_path), exist_ok=True)
    os.makedirs(os.path.dirname(rec_path), exist_ok=True)
    npep = nrec = 0

    def hdr(desc, default):
        d = (desc or "").replace("\r", " ").replace("\n", " ").strip() or default
        return d

    with open(csv_path) as fh, open(pep_path, "w") as pf, open(rec_path, "w") as rf:
        for r in csv.DictReader(fh, delimiter=";"):
            eid = r["id"]
            pseq = (r.get("PEPTIDE_SEQ") or "").strip()
            rseq = (r.get("PROTEIN_SEQ") or "").strip()
            if pseq:
                pf.write(f">{eid} | {hdr(r.get('PEPTIDE_DESC'), 'Peptide')}\n{pseq}\n")
                npep += 1
            if rseq:
                rf.write(f">{eid} | {hdr(r.get('PROTEIN_DESC'), 'Protein')}\n{rseq}\n")
                nrec += 1
    return npep, nrec


def package_contacts(cocada_dir, out_dir, ids):
    """Place each entry's COCaDA contacts at contacts/<id>/<PDB>_contacts.csv
    (renaming our <id>_contacts.csv -> <PDB>_contacts.csv, as the site expects)."""
    _clear_make(out_dir)
    n = 0
    for eid in ids:
        src = os.path.join(cocada_dir, eid, f"{eid}_contacts.csv")
        if not os.path.exists(src):
            continue
        pdb = eid.split("-")[0]
        dst_dir = os.path.join(out_dir, eid)
        os.makedirs(dst_dir, exist_ok=True)
        shutil.copy2(src, os.path.join(dst_dir, f"{pdb}_contacts.csv"))
        n += 1
    return n


def package_multipro_contacts(mp_cocada_dir, out_dir, cluster_ids):
    """Place each Multipro entry's COCaDA contacts (run on the multi-chain complex)
    at multipro/contacts/<cluster_id>/<PDB>_contacts.csv."""
    _clear_make(out_dir)
    n = 0
    for cid in cluster_ids:
        pdb = cid.split("-")[0]
        src = os.path.join(mp_cocada_dir, cid, f"{pdb}_contacts.csv")
        if not os.path.exists(src):
            continue
        dst_dir = os.path.join(out_dir, cid)
        os.makedirs(dst_dir, exist_ok=True)
        shutil.copy2(src, os.path.join(dst_dir, f"{pdb}_contacts.csv"))
        n += 1
    return n


def layout_structures(cif_dir, out_dir, ids, ext_out=".cif"):
    """Shard the per-entry complex structures into the site's structure tree:
        <out_dir>/<id[0]>/<id><ext_out>
    Source is the pipeline's lossless mmCIF (`<cif_dir>/<id>.cif`). The web app
    serves these under data/<mode>/pdb/ (pep-pro) and data/<mode>/multipro/pdb/
    (multipro). We emit mmCIF (`.cif`) — the site's structure viewer reads mmCIF.

    No silent loss (mirrors make_zips.py): every id whose source CIF is absent is
    collected and returned so the caller can WARN and record it. Returns
    (n_copied, missing_ids)."""
    _clear_make(out_dir)
    cif_dir = os.path.expanduser(cif_dir)
    n = 0
    missing = []
    for eid in ids:
        src = os.path.join(cif_dir, f"{eid}.cif")
        if not os.path.exists(src):
            missing.append(eid)
            continue
        shard = os.path.join(out_dir, eid[0])
        os.makedirs(shard, exist_ok=True)
        shutil.copy2(src, os.path.join(shard, f"{eid}{ext_out}"))
        n += 1
    return n, missing


def copy_clusters(legacy_dir, out_dir):
    """Copy the cluster tables the site's Clusters page reads. What we have comes
    from the inherited legacy clusters; therapeutic-class lists (AAP..SBP) and
    pdb_classes are TODO (see docs)."""
    _clear_make(out_dir)
    copied = []
    src = os.path.expanduser(legacy_dir)
    if os.path.isdir(src):
        for f in os.listdir(src):
            if f.endswith((".tsv", ".csv")):
                shutil.copy2(os.path.join(src, f), os.path.join(out_dir, f))
                copied.append(f)
    return copied


def main():
    p = snakemake.params                                   # noqa: F821
    web = os.path.expanduser(p.web_dir)
    mode = p.mode
    base = os.path.join(web, "data", mode)

    # optional target column orders (for the final permutation to the site layout)
    target = _load_order(getattr(p, "column_order", None))
    target_mp = _load_order(getattr(p, "column_order_multipro", None))

    # --- pep-pro: per-entry csv + contacts ---
    n_csv, header = split_per_entry(snakemake.input.propedia,            # noqa: F821
                                    os.path.join(base, "csv"), target=target)
    ids = [r["id"] for r in csv.DictReader(open(snakemake.input.propedia),  # noqa: F821
                                           delimiter=";")]
    n_con = package_contacts(os.path.expanduser(p.cocada_dir),
                             os.path.join(base, "contacts"), ids)

    # --- multipro: per-entry csv + contacts (COCaDA on the multi-chain complex) ---
    n_mp, mp_header = split_per_entry(snakemake.input.multipro,          # noqa: F821
                                      os.path.join(base, "multipro", "csv"),
                                      id_col="cluster_id", target=target_mp)
    mp_ids = [r["cluster_id"] for r in csv.DictReader(open(snakemake.input.multipro),  # noqa: F821
                                                      delimiter=";")]
    n_mpcon = package_multipro_contacts(os.path.expanduser(p.multipro_cocada_dir),
                                        os.path.join(base, "multipro", "contacts"),
                                        mp_ids)

    # --- per-entry structure trees (pep-pro + multipro), the site's 3D viewer ---
    # Emit only when a source CIF dir is configured (keeps package usable without
    # shipping structures). Structures are the lossless mmCIF from state/<mode>/.
    struct_subdir = getattr(p, "structure_subdir", "pdb")
    n_pdb = n_mppdb = 0
    miss_pdb, miss_mppdb = [], []
    cif_dir = getattr(p, "cif_dir", "") or ""
    mp_cif_dir = getattr(p, "multipro_cif_dir", "") or ""
    if cif_dir:
        n_pdb, miss_pdb = layout_structures(
            cif_dir, os.path.join(base, struct_subdir), ids)
    if mp_cif_dir:
        n_mppdb, miss_mppdb = layout_structures(
            mp_cif_dir, os.path.join(base, "multipro", struct_subdir), mp_ids)

    # --- Explore summary TSV (pep-pro): data/<explore_tsv_name>, read by fixed index ---
    explore_name = getattr(p, "explore_tsv_name", "propedia26_v17.tsv")
    n_expl = build_explore_tsv(snakemake.input.propedia,                 # noqa: F821
                               os.path.join(web, "data", explore_name))

    # --- master bulk CSVs (Download page): the whole tables, unsplit ---
    # propedia.csv -> data/propedia26_v17.csv ; multipro_final.csv -> data/multipro_v4.csv
    peppro_csv_name = getattr(p, "peppro_csv_name", "propedia26_v17.csv")
    multipro_csv_name = getattr(p, "multipro_csv_name", "multipro_v4.csv")
    os.makedirs(os.path.join(web, "data"), exist_ok=True)
    shutil.copy2(snakemake.input.propedia,                              # noqa: F821
                 os.path.join(web, "data", peppro_csv_name))
    shutil.copy2(snakemake.input.multipro,                             # noqa: F821
                 os.path.join(web, "data", multipro_csv_name))

    # --- BLAST subject FASTAs (pep-pro): data/peptides.fasta, data/receptors.fasta ---
    n_pepf, n_recf = build_blast_fastas(snakemake.input.propedia,        # noqa: F821
                                        os.path.join(web, "data", "peptides.fasta"),
                                        os.path.join(web, "data", "receptors.fasta"))

    # --- clusters (shared across modes) ---
    clusters = copy_clusters(p.legacy_dir, os.path.join(web, "data", "clusters"))

    # --- column manifests (drive the later permutation to the site's order) ---
    with open(os.path.join(web, "columns_peppro.txt"), "w") as fh:
        for i, c in enumerate(header):
            fh.write(f"{i}\t{c}\n")
    with open(os.path.join(web, "columns_multipro.txt"), "w") as fh:
        for i, c in enumerate(mp_header):
            fh.write(f"{i}\t{c}\n")

    with open(snakemake.output.marker, "w") as fh:         # noqa: F821
        fh.write(f"peppro_csv={n_csv} contacts={n_con} multipro_csv={n_mp} "
                 f"multipro_contacts={n_mpcon} explore_tsv={n_expl} "
                 f"peptides_fasta={n_pepf} receptors_fasta={n_recf} "
                 f"pdb={n_pdb} multipro_pdb={n_mppdb} "
                 f"clusters={len(clusters)}\n")
        if miss_pdb:
            fh.write(f"WARNING: {len(miss_pdb)} pep-pro entries have NO complex "
                     f"CIF in {cif_dir} (structure tree incomplete)\n")
        if miss_mppdb:
            fh.write(f"WARNING: {len(miss_mppdb)} multipro entries have NO complex "
                     f"CIF in {mp_cif_dir} (structure tree incomplete)\n")
    print(f"PACKAGED -> {base}\n  per-entry csv: {n_csv}\n  contacts: {n_con}\n"
          f"  multipro csv: {n_mp}\n  multipro contacts: {n_mpcon}\n"
          f"  structures ({struct_subdir}/.cif): pep-pro={n_pdb}"
          + (f" (MISSING {len(miss_pdb)})" if miss_pdb else "")
          + f", multipro={n_mppdb}"
          + (f" (MISSING {len(miss_mppdb)})" if miss_mppdb else "") + "\n"
          f"  master csv: data/{peppro_csv_name}, data/{multipro_csv_name}\n"
          f"  explore tsv: {n_expl} rows -> data/{explore_name}\n"
          f"  blast fasta: peptides={n_pepf} receptors={n_recf} -> data/*.fasta\n"
          f"  clusters copied: {len(clusters)} {clusters}\n"
          f"  column manifests: columns_peppro.txt ({len(header)} cols), "
          f"columns_multipro.txt ({len(mp_header)} cols)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
