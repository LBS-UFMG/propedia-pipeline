"""Structural signatures via SIGNA (aCSM-ALL). Runs on peptide-only PDBs.
Propedia 26 params: cutoff_limit=10, cutoff_step=0.2, cumulative=True -> 1800 features.

To make SIGNA resumable and parallel, the peptide PDBs are split into CHUNKS; each
chunk computes signatures (calling SIGNA's per-file ``read()``) and writes a durable
chunk CSV. Chunk files are content-addressed by a hash of their member ids, so an
interrupted run resumes (existing chunks skipped) and a changed sample auto-invalidates
only the chunks that changed. The final CSV is the concatenation of the current chunks."""
import concurrent.futures as cf
import contextlib
import hashlib
import os
import sys

import checkpoint

VERSION = "2"           # bump: per-file read() (SystemExit-safe) instead of read_folder
CHUNK = 500             # peptides per chunk (checkpoint granularity)


def _chunk_key(names):
    return hashlib.sha1(",".join(sorted(names)).encode()).hexdigest()[:16]


def run_chunk(item):
    """Compute aCSM signatures for one chunk of peptide PDBs -> chunk CSV (id,feat...).

    Calls SIGNA's ``read()`` PER FILE (which is exactly what ``read_folder`` does
    internally) so we can catch the ``exit()`` SIGNA calls on a degenerate peptide
    ("No match found." -> bare ``exit()``, i.e. SystemExit). read_folder is all-or-
    nothing: one such peptide aborts the whole batch AND kills the host process. Here a
    bad peptide is skipped and counted; the rest of the chunk still gets written. Returns
    (out_csv, n_skipped)."""
    signa_dir, pep_dir, names, out_csv, cl, cs, cum = item
    if signa_dir not in sys.path:
        sys.path.insert(0, signa_dir)
    import signa                                           # noqa: E402
    tmp_out = out_csv + ".tmp"
    skipped = 0
    with open(tmp_out, "w") as fout, open(os.devnull, "w") as devnull:
        for nm in names:
            path = os.path.abspath(os.path.join(pep_dir, nm))
            eid = os.path.splitext(nm)[0]
            try:
                # mirror read_folder's positional call; verbose off; hush SIGNA's prints
                with contextlib.redirect_stdout(devnull):
                    sig = signa.read(path, "acsm-all", cl, cs, True, "ALL",
                                     False, cum, ",", "AMBER")
            except SystemExit:            # SIGNA exit() on empty/degenerate peptide
                skipped += 1
                continue
            except Exception:             # any other SIGNA failure on one file -> skip it
                skipped += 1
                continue
            if sig:
                fout.write(eid + "," + str(sig).rstrip("\n") + "\n")
    os.replace(tmp_out, out_csv)
    return out_csv, skipped


def main():
    p = snakemake.params                                   # noqa: F821
    signa_dir = os.path.expanduser(p.signa_dir)
    pep_dir = os.path.abspath(p.pep_pdb_dir)
    threads = getattr(snakemake, "threads", 1) or 1        # noqa: F821

    workdir = checkpoint.namespace(os.path.abspath(p.ckpt), VERSION,
                                   {"cl": p.cutoff_limit, "cs": p.cutoff_step,
                                    "cum": p.cumulative})

    pdbs = sorted(f for f in os.listdir(pep_dir) if f.endswith(".pdb"))
    chunks = [pdbs[i:i + CHUNK] for i in range(0, len(pdbs), CHUNK)]
    chunk_files, todo, want = [], [], set()
    for names in chunks:
        key = _chunk_key(names)
        out_csv = os.path.join(workdir, f"chunk_{key}.csv")
        chunk_files.append(out_csv)
        want.add(f"chunk_{key}.csv")
        if not os.path.exists(out_csv):    # resume: skip already-computed chunks
            todo.append((signa_dir, pep_dir, names, out_csv,
                         p.cutoff_limit, p.cutoff_step, p.cumulative))
    print(f"[signa] {len(chunks)} chunks, {len(todo)} to compute (threads={threads})",
          file=sys.stderr)

    skipped = 0
    if threads > 1 and todo:
        with cf.ProcessPoolExecutor(max_workers=threads) as ex:
            for i, (_, sk) in enumerate(ex.map(run_chunk, todo), 1):
                skipped += sk
                if i % 5 == 0:
                    print(f"[signa] {i}/{len(todo)} chunks done", file=sys.stderr)
    else:
        for i, it in enumerate(todo, 1):
            _, sk = run_chunk(it)
            skipped += sk
            if i % 5 == 0:
                print(f"[signa] {i}/{len(todo)} chunks done", file=sys.stderr)
    if skipped:
        print(f"[signa] skipped {skipped} peptides SIGNA could not process "
              f"(empty/degenerate; 'No match found.')", file=sys.stderr)

    for f in os.listdir(workdir):          # prune stale chunks (changed sample)
        if f.startswith("chunk_") and f not in want:
            os.remove(os.path.join(workdir, f))

    n = 0
    with open(snakemake.output.signatures, "w") as fout:   # noqa: F821
        for cf_path in chunk_files:
            with open(cf_path) as fin:
                for line in fin:
                    fout.write(line)
                    n += 1
    print(f"DONE {n} signatures normalized", file=sys.stderr)


if __name__ == "__main__":
    main()
