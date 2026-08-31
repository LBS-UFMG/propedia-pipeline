"""Build the website's bulk-download ZIP bundles from EXISTING pipeline outputs.

READ-ONLY over `state/<mode>/` and `results/<mode>/`: it only reads the finished
structures/tables and writes zip files. It never recomputes, deletes, or modifies
any pipeline output, so it is safe to run against a completed build without any risk
of re-running the pipeline.

Bundles (names configurable via config `package.zip_names`):
  propedia.zip            pep-pro complex structures (mmCIF), one per propedia.csv entry
  peptides_pdb.zip        peptide-only PDBs, one per propedia.csv entry
  multipro.zip            multipro complex structures (mmCIF), the whole multipro set
  sequence_signature.zip  iFeature sequence signatures (seq_signatures.tsv [+ excluded])
  structural_signature.zip SIGNA structural signatures (struct_signatures.csv)
  clusters.zip            the cluster tables (from the legacy clusters dir)

Safety features (so "no information is lost" is VISIBLE, never silent):
  - each zip is written to <name>.tmp then atomically renamed (a kill never leaves a
    corrupt/partial zip that looks complete);
  - for the per-entry structure zips, the expected id set comes from the final CSV;
    any id whose structure file is missing on disk is COUNTED and written to
    <out>/<name>.missing.txt, and a WARNING is printed — nothing is dropped silently;
  - a summary + a `.zipped` marker record every count.

NOTE ON FORMAT: complexes are shipped as mmCIF (`.cif`) — the pipeline's lossless
intermediate. Peptide-only structures are PDB (`.pdb`, as produced for SIGNA).
"""
import csv
import os
import sys
import zipfile

DEFAULT_NAMES = {
    "peppro":     "propedia.zip",
    "peptides":   "peptides_pdb.zip",
    "multipro":   "multipro.zip",
    "seq_sig":    "sequence_signature.zip",
    "struct_sig": "structural_signature.zip",
    "clusters":   "clusters.zip",
}


def _ids(csv_path, id_col, delim=";"):
    with open(csv_path) as fh:
        return [r[id_col] for r in csv.DictReader(fh, delimiter=delim)]


def _new_zip(path, level):
    return zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED, compresslevel=level)


def _zip_by_id(out_dir, name, ids, src_dir, ext, arc_prefix, level):
    """Zip src_dir/<id><ext> for each id. Returns (added, missing_ids). Atomic."""
    src_dir = os.path.expanduser(src_dir)
    final = os.path.join(out_dir, name)
    tmp = final + ".tmp"
    added, missing = 0, []
    with _new_zip(tmp, level) as z:
        for eid in ids:
            src = os.path.join(src_dir, f"{eid}{ext}")
            if os.path.exists(src):
                z.write(src, f"{arc_prefix}/{eid}{ext}")
                added += 1
            else:
                missing.append(eid)
        n_members = len(z.namelist())
    assert n_members == added, f"{name}: zip member count {n_members} != added {added}"
    os.replace(tmp, final)
    if missing:
        with open(os.path.join(out_dir, f"{name}.missing.txt"), "w") as fh:
            fh.write("\n".join(missing) + "\n")
    return added, missing


def _zip_dir(out_dir, name, src_dir, exts, arc_prefix, level):
    """Zip every file in src_dir matching exts. Returns count. Atomic."""
    src_dir = os.path.expanduser(src_dir)
    final = os.path.join(out_dir, name)
    tmp = final + ".tmp"
    added = 0
    with _new_zip(tmp, level) as z:
        if os.path.isdir(src_dir):
            for f in sorted(os.listdir(src_dir)):
                if f.startswith(".") or not f.endswith(exts):
                    continue
                z.write(os.path.join(src_dir, f), f"{arc_prefix}/{f}")
                added += 1
        n_members = len(z.namelist())
    assert n_members == added, f"{name}: zip member count {n_members} != added {added}"
    os.replace(tmp, final)
    return added


def _zip_files(out_dir, name, files, level):
    """Zip a fixed list of individual files (flat). Returns count. Atomic."""
    final = os.path.join(out_dir, name)
    tmp = final + ".tmp"
    added = 0
    with _new_zip(tmp, level) as z:
        for src in files:
            if os.path.exists(src):
                z.write(src, os.path.basename(src))
                added += 1
    os.replace(tmp, final)
    return added


def main():
    p = snakemake.params                                       # noqa: F821
    inp = snakemake.input                                      # noqa: F821
    out_dir = os.path.expanduser(p.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    names = {**DEFAULT_NAMES, **(dict(getattr(p, "names", {}) or {}))}
    level = int(getattr(p, "compresslevel", 6))

    peppro_ids = _ids(inp.propedia, "id")
    n_mp_rows = max(0, sum(1 for _ in open(inp.multipro)) - 1)  # ; header
    warnings = []

    # 1. pep-pro complexes (dataset entries only)
    n_cx, miss_cx = _zip_by_id(out_dir, names["peppro"], peppro_ids,
                               p.cif_dir, ".cif", "propedia", level)
    # 2. peptide-only PDBs (dataset entries; X-only peptides may legitimately lack one)
    n_pep, miss_pep = _zip_by_id(out_dir, names["peptides"], peppro_ids,
                                 p.pep_pdb_dir, ".pdb", "peptides", level)
    # 3. multipro complexes (the whole produced set)
    n_mp = _zip_dir(out_dir, names["multipro"], p.multipro_cif_dir,
                    (".cif",), "multipro", level)
    # 4. sequence signatures
    seq_files = [inp.seq_sig]
    excl = os.path.join(os.path.dirname(inp.seq_sig), "seq_signatures_excluded.tsv")
    if os.path.exists(excl):
        seq_files.append(excl)
    n_seq = _zip_files(out_dir, names["seq_sig"], seq_files, level)
    # 5. structural signatures
    n_str = _zip_files(out_dir, names["struct_sig"], [inp.struct_sig], level)
    # 6. cluster tables
    n_clu = _zip_dir(out_dir, names["clusters"], p.legacy_dir,
                     (".tsv", ".csv"), "clusters", level)

    # ---- integrity reporting (no silent loss) ----
    if miss_cx:
        warnings.append(f"{names['peppro']}: {len(miss_cx)} dataset entries have NO "
                        f"complex CIF on disk (see {names['peppro']}.missing.txt)")
    if miss_pep:
        warnings.append(f"{names['peptides']}: {len(miss_pep)} dataset entries have NO "
                        f"peptide PDB (X-only peptides can legitimately lack one; "
                        f"see {names['peptides']}.missing.txt)")
    if n_mp != n_mp_rows:
        warnings.append(f"{names['multipro']}: zipped {n_mp} CIFs but multipro_final "
                        f"has {n_mp_rows} rows (mismatch — investigate before shipping)")

    lines = [
        f"propedia entries (propedia.csv) : {len(peppro_ids)}",
        f"  {names['peppro']:24s}: {n_cx} complexes"
        + (f"  (MISSING {len(miss_cx)})" if miss_cx else ""),
        f"  {names['peptides']:24s}: {n_pep} peptide PDBs"
        + (f"  (missing {len(miss_pep)})" if miss_pep else ""),
        f"multipro_final rows            : {n_mp_rows}",
        f"  {names['multipro']:24s}: {n_mp} complexes",
        f"  {names['seq_sig']:24s}: {n_seq} file(s)",
        f"  {names['struct_sig']:24s}: {n_str} file(s)",
        f"  {names['clusters']:24s}: {n_clu} cluster table(s)",
    ]
    report = "\n".join(lines)
    print("ZIP BUNDLES ->", out_dir, file=sys.stderr)
    print(report, file=sys.stderr)
    for w in warnings:
        print("WARNING:", w, file=sys.stderr)

    with open(snakemake.output.marker, "w") as fh:             # noqa: F821
        fh.write(report + "\n")
        for w in warnings:
            fh.write("WARNING: " + w + "\n")


if __name__ == "__main__":
    main()
