"""Run COCaDA on the downloaded CIF for each sampled entry, with per-entry
checkpointing (resumes after interruption) and a bounded process pool.

Invocation: -c pep,prot -inter -ph 7.4 -s -m 1. Runs on the ORIGINAL CIF (not
our extracted CIF), so contacts are independent of extraction convention.

CHAIN-ID TRANSLATION (important): the rest of Propedia uses auth_asym_id chain
IDs, but COCaDA selects on label_asym_id. These diverge in ~90% of PDB entries.
We translate each pair's auth chain ID to its POLYMER label_asym_id (built from
the CIF's _atom_site ATOM records; 1:1 for polymers) before calling COCaDA.
Chains with no polymer label are reported explicitly, not silently zeroed.

Each COCaDA call is -m 1 (single-core); Snakemake/the pool owns parallelism."""
import csv
import gzip
import os
import subprocess
import sys
import tempfile

import checkpoint

VERSION = "2"            # bump on logic change -> invalidates old checkpoints


def shard_path(cif_dir, pid):
    sub = pid[1:3].lower() if len(pid) >= 3 else "_"
    return os.path.join(cif_dir, sub, f"{pid}.cif.gz")


def auth_to_polymer_label(cif_path):
    """Map auth_asym_id -> polymer label_asym_id from a CIF's _atom_site block."""
    cols = []
    in_loop = False
    mapping = {}
    with open(cif_path, "rt") as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                cols.append(line.strip())
                in_loop = True
                continue
            if in_loop and (line.startswith("ATOM") or line.startswith("HETATM")):
                idx = {c: i for i, c in enumerate(cols)}
                parts = line.split()
                try:
                    if parts[idx["_atom_site.group_PDB"]] != "ATOM":
                        continue
                    lab = parts[idx["_atom_site.label_asym_id"]]
                    aut = parts[idx["_atom_site.auth_asym_id"]]
                except (KeyError, IndexError):
                    continue
                mapping.setdefault(aut, lab)          # first polymer label wins
            elif in_loop and line.startswith("#"):
                break
    return mapping


def run_one(args):
    """worker(item) -> (record, status). 'retry:*' statuses are not checkpointed."""
    entry_id, pid, pep, prot, cif_dir, cocada_dir, ph, out_root = args
    src = shard_path(cif_dir, pid)
    if not os.path.exists(src):
        return None, "missing_cif"
    entry_out = os.path.abspath(os.path.join(out_root, entry_id))
    # Clear any prior output so a stale *_contacts.csv can't be picked up.
    if os.path.isdir(entry_out):
        for _f in os.listdir(entry_out):
            os.remove(os.path.join(entry_out, _f))
    os.makedirs(entry_out, exist_ok=True)
    tmp = tempfile.NamedTemporaryFile("w", suffix=".cif", delete=False)
    try:
        with gzip.open(src, "rt") as fh:
            tmp.write(fh.read())
        tmp.close()
        a2l = auth_to_polymer_label(tmp.name)
        pep_lab, prot_lab = a2l.get(pep), a2l.get(prot)
        if pep_lab is None or prot_lab is None:
            missing = [c for c, lab in ((pep, pep_lab), (prot, prot_lab)) if lab is None]
            return None, f"no_polymer_label:{','.join(missing)}"
        # Pin numpy's OpenBLAS to ONE thread per subprocess. COCaDA does small
        # per-structure linear algebra; without this, OpenBLAS fans each cocada.py
        # out to ~3 cores, so a 32-worker pool demands ~96 cores and oversubscribes
        # the box (load >> ncpu). Parallelism comes from the POOL across pairs, not
        # from BLAS within one pair. Same numeric result, exact core accounting.
        env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
                   MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
        proc = subprocess.run(
            ["python3", "cocada.py", "-f", tmp.name,
             "-c", f"{pep_lab},{prot_lab}", "-inter", "-m", "1",
             "-o", entry_out, "-ph", str(ph), "-s"],
            cwd=cocada_dir, capture_output=True, text=True, timeout=300, env=env)
        produced = [f for f in os.listdir(entry_out) if f.endswith("_contacts.csv")]
        if not produced:
            return None, "no_output"
        canonical = os.path.join(entry_out, f"{entry_id}_contacts.csv")
        src_csv = os.path.join(entry_out, produced[0])
        if src_csv != canonical:
            os.replace(src_csv, canonical)
        with open(canonical) as fh:
            n = sum(1 for _ in fh) - 1     # minus header
        return {"pdb": pid, "pep": pep, "prot": prot, "contacts": n}, "ok"
    except subprocess.TimeoutExpired:
        return None, "retry:timeout"
    except Exception as exc:                       # noqa: BLE001
        return None, f"retry:exc:{exc}"
    finally:
        os.unlink(tmp.name)


def main():
    p = snakemake.params                                   # noqa: F821
    cocada_dir = os.path.expanduser(p.cocada_dir)
    out_root = os.path.abspath(p.out_root)
    threads = getattr(snakemake, "threads", 1) or 1        # noqa: F821
    os.makedirs(out_root, exist_ok=True)
    pairs = list(csv.DictReader(open(snakemake.input.pairs), delimiter="\t"))  # noqa: F821

    items = [(r["id"], r["pdb"], r["pep_chain"], r["prot_chain"],
              p.cif_dir, cocada_dir, p.ph, out_root) for r in pairs]
    workdir = checkpoint.namespace(p.ckpt, VERSION, {"ph": p.ph})
    results = checkpoint.run(items, run_one, workdir, threads=threads,
                             id_of=lambda it: it[0], stage="cocada")

    ok, reasons = 0, {}
    with open(snakemake.output.summary, "w") as summ:      # noqa: F821
        summ.write("id\tpdb\tpep\tprot\ttotal_contacts\n")
        for r in pairs:
            eid = r["id"]
            res = results.get(eid)
            if res and res["status"] == "ok" and res["record"]:
                v = res["record"]
                summ.write(f"{eid}\t{v['pdb']}\t{v['pep']}\t{v['prot']}\t{v['contacts']}\n")
                ok += 1
            else:
                status = res["status"] if res else "missing"
                reasons[status] = reasons.get(status, 0) + 1
    print(f"DONE ok={ok} reasons={reasons}", file=sys.stderr)


if __name__ == "__main__":
    main()
