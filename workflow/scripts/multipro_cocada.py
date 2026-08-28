"""COCaDA on the multi-chain complex (Multipro), per-entry checkpointed.

Runs on the extracted multipro CIF (state/<mode>/multipro_cif/<cluster_id>.cif),
which contains only the peptide + its protein chains, relabeled single-char by
MMCIFIO. So, unlike the pep-pro stage, there is no auth->label translation and no
`-c` selection: COCaDA analyses all interchain contacts in the complex (`-inter`).

Output: one <PDB>_contacts.csv per Multipro entry under out_root/<cluster_id>/
(same format the website serves), plus a small summary. The packager copies these
into multipro/contacts/<cluster_id>/<PDB>_contacts.csv.
"""
import csv
import os
import subprocess
import sys

import checkpoint

VERSION = "1"


def run_one(args):
    """worker(item) -> (record, status). 'retry:*' statuses are not checkpointed."""
    cluster_id, cif_path, cocada_dir, ph, out_root = args
    if not os.path.exists(cif_path):
        return None, "missing_cif"
    entry_out = os.path.abspath(os.path.join(out_root, cluster_id))
    if os.path.isdir(entry_out):
        for _f in os.listdir(entry_out):
            os.remove(os.path.join(entry_out, _f))
    os.makedirs(entry_out, exist_ok=True)
    try:
        proc = subprocess.run(
            ["python3", "cocada.py", "-f", os.path.abspath(cif_path),
             "-inter", "-m", "1", "-o", entry_out, "-ph", str(ph), "-s"],
            cwd=cocada_dir, capture_output=True, text=True, timeout=600)
        produced = [f for f in os.listdir(entry_out) if f.endswith("_contacts.csv")]
        if not produced:
            return None, "no_output"
        pdb = cluster_id.split("-")[0]
        canonical = os.path.join(entry_out, f"{pdb}_contacts.csv")
        src = os.path.join(entry_out, produced[0])
        if src != canonical:
            os.replace(src, canonical)
        with open(canonical) as fh:
            n = sum(1 for _ in fh) - 1
        return {"contacts": n}, "ok"
    except subprocess.TimeoutExpired:
        return None, "retry:timeout"
    except Exception as exc:                       # noqa: BLE001
        return None, f"retry:exc:{exc}"


def main():
    p = snakemake.params                                   # noqa: F821
    cocada_dir = os.path.expanduser(p.cocada_dir)
    out_root = os.path.abspath(p.out_root)
    threads = getattr(snakemake, "threads", 1) or 1        # noqa: F821
    os.makedirs(out_root, exist_ok=True)
    rows = list(csv.DictReader(open(snakemake.input.multipro), delimiter="\t"))  # noqa: F821

    items = [(r["cluster_id"], os.path.join(p.cif_dir, f'{r["cluster_id"]}.cif'),
              cocada_dir, p.ph, out_root) for r in rows]
    workdir = checkpoint.namespace(p.ckpt, VERSION, {"ph": p.ph})
    results = checkpoint.run(items, run_one, workdir, threads=threads,
                             id_of=lambda it: it[0], stage="multipro_cocada")

    ok, reasons = 0, {}
    with open(snakemake.output.summary, "w") as fh:        # noqa: F821
        fh.write("cluster_id\ttotal_contacts\n")
        for r in rows:
            cid = r["cluster_id"]
            res = results.get(cid)
            if res and res["status"] == "ok" and res["record"]:
                fh.write(f"{cid}\t{res['record']['contacts']}\n")
                ok += 1
            else:
                st = res["status"] if res else "missing"
                reasons[st] = reasons.get(st, 0) + 1
    print(f"DONE multipro cocada: ok={ok} reasons={reasons}", file=sys.stderr)


if __name__ == "__main__":
    main()
