"""Structural signatures via SIGNA (aCSM-ALL). Runs on peptide-only PDBs.
Propedia 26 params: cutoff_limit=10, cutoff_step=0.2, cumulative=True -> 1800 features.

SIGNA is a folder tool (one call over a whole folder). To make it resumable and
parallel, the peptide PDBs are split into CHUNKS; each chunk runs SIGNA over a
temp dir of symlinks and writes a durable chunk CSV. Chunk files are
content-addressed by a hash of their member ids, so an interrupted run resumes
(existing chunks skipped) and a changed sample auto-invalidates only the chunks
that changed. The final CSV is the concatenation of the current chunks."""
import concurrent.futures as cf
import csv
import hashlib
import os
import shutil
import sys
import tempfile

import checkpoint

VERSION = "1"
CHUNK = 500              # peptides per SIGNA invocation


def _chunk_key(names):
    return hashlib.sha1(",".join(sorted(names)).encode()).hexdigest()[:16]


def run_chunk(item):
    """Run SIGNA over one chunk of peptide PDBs -> normalized chunk CSV (id,feat...)."""
    signa_dir, pep_dir, names, out_csv, cl, cs, cum = item
    sys.path.insert(0, signa_dir)
    import signa                                           # noqa: E402
    tmp = tempfile.mkdtemp()
    try:
        for nm in names:
            os.symlink(os.path.abspath(os.path.join(pep_dir, nm)),
                       os.path.join(tmp, nm))
        raw = out_csv + ".raw"
        signa.read_folder(folder=tmp, signa_type="acsm-all", cumulative=cum,
                          output=raw, cutoff_limit=cl, cutoff_step=cs, format="pdb")
        # normalize the leading file PATH to the entry id, write atomically
        tmp_out = out_csv + ".tmp"
        with open(raw) as fin, open(tmp_out, "w") as fout:
            for line in fin:
                parts = line.rstrip("\n").split(",")
                if not parts or not parts[0]:
                    continue
                eid = os.path.splitext(os.path.basename(parts[0]))[0]
                fout.write(eid + "," + ",".join(parts[1:]) + "\n")
        os.replace(tmp_out, out_csv)
        os.remove(raw)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return out_csv


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

    if threads > 1 and todo:
        with cf.ProcessPoolExecutor(max_workers=threads) as ex:
            for i, _ in enumerate(ex.map(run_chunk, todo), 1):
                if i % 5 == 0:
                    print(f"[signa] {i}/{len(todo)} chunks done", file=sys.stderr)
    else:
        for i, it in enumerate(todo, 1):
            run_chunk(it)
            if i % 5 == 0:
                print(f"[signa] {i}/{len(todo)} chunks done", file=sys.stderr)

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
