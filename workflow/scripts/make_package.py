"""Package the pipeline outputs into the propedia26 website file tree.

Emits the per-entry layout the web app serves (see docs/propedia26_publishing.md):
  <web>/data/<mode>/csv/<id[0]>/<id>.csv            per-entry row (;-delimited, no header)
  <web>/data/<mode>/contacts/<id>/<PDB>_contacts.csv per-entry COCaDA contacts
  <web>/data/<mode>/multipro/csv/<id[0]>/<id>.csv    per-entry multipro row
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
    one hook needed to match the site's final layout later. Returns (n, header_out)."""
    _clear_make(out_dir)
    with open(csv_path) as fh:
        src = fh.readline().rstrip("\n").split(";")
        pos = {c: i for i, c in enumerate(src)}
        idx = pos[id_col]
        header_out = target if target else src
        n = 0
        for line in fh:
            f = line.rstrip("\n").split(";")
            eid = f[idx]
            out = [(f[pos[c]] if c in pos and pos[c] < len(f) else "") for c in header_out] \
                if target else f
            shard = os.path.join(out_dir, eid[0])
            os.makedirs(shard, exist_ok=True)
            with open(os.path.join(shard, f"{eid}.csv"), "w") as ofh:
                ofh.write(";".join(out) + "\n")
            n += 1
    return n, header_out


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
                 f"multipro_contacts={n_mpcon} clusters={len(clusters)}\n")
    print(f"PACKAGED -> {base}\n  per-entry csv: {n_csv}\n  contacts: {n_con}\n"
          f"  multipro csv: {n_mp}\n  multipro contacts: {n_mpcon}\n"
          f"  clusters copied: {len(clusters)} {clusters}\n"
          f"  column manifests: columns_peppro.txt ({len(header)} cols), "
          f"columns_multipro.txt ({len(mp_header)} cols)\n"
          f"  TODO (docs): Explore summary tsv, final column permutation to the site order",
          file=sys.stderr)


if __name__ == "__main__":
    main()
