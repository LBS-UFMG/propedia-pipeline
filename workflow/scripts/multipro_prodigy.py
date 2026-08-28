"""Multipro Phase 2c: PRODIGY per protein chain, run on the multi-chain complex
so values reproduce v15's multi-chain-context numbers (which differ from the
isolated pep-pro pair -- other chains affect NIS). Colon-joined per chain in
PROTEIN_CHAIN order. Reuses run_prodigy's output parser.
"""
import csv
import os
import subprocess
import sys

import run_prodigy as rp
import checkpoint

VERSION = "2"            # bump on logic change -> invalidates old checkpoints

DEG = "\u02da"
# v15 multipro column -> run_prodigy.parse key
COLMAP = [
    ("No. of intermolecular contacts", "n_intermolecular"),
    ("No. of charged-charged contacts", "n_charged_charged"),
    ("No. of charged-polar contacts", "n_charged_polar"),
    ("No. of charged-apolar contacts", "n_charged_apolar"),
    ("No. of polar-polar contacts", "n_polar_polar"),
    ("No. of apolar-polar contacts", "n_apolar_polar"),
    ("No. of apolar-apolar contacts", "n_apolar_apolar"),
    ("Percentage of apolar NIS residues", "nis_apolar"),
    ("Percentage of charged NIS residues", "nis_charged"),
    ("Predicted binding affinity (kcal.mol-1)", "dg"),
    (f"Predicted dissociation constant (M) at 25.0{DEG}C", "kd"),
]


def run_chain(pdb_path, pep, prot, cutoff, temp):
    try:
        proc = subprocess.run(["prodigy", pdb_path, "--selection", pep, prot,
                               "--distance-cutoff", str(cutoff),
                               "--temperature", str(temp)],
                              capture_output=True, text=True, timeout=rp.TIMEOUT)
        return rp.parse(proc.stdout + proc.stderr)
    except Exception:                                       # noqa: BLE001
        return {}


def worker(item):
    """worker(item) -> (record, status). Runs PRODIGY per protein chain on the
    multi-chain complex; colon-joins per chain in PROTEIN_CHAIN order."""
    cid, cif_path, pep, prots, cutoff, temp = item
    try:
        per_chain = [run_chain(cif_path, pep, pc, cutoff, temp) for pc in prots]
    except Exception as exc:                               # noqa: BLE001
        return None, f"retry:exc:{exc}"
    rec = {col: ":".join(str(v.get(key, "")) for v in per_chain)
           for col, key in COLMAP}
    return rec, "ok"


def main():
    p = snakemake.params                                   # noqa: F821
    threads = getattr(snakemake, "threads", 1) or 1        # noqa: F821
    rows = list(csv.DictReader(open(snakemake.input.multipro), delimiter="\t"))  # noqa: F821
    out_cols = ["cluster_id"] + [c for c, _ in COLMAP]

    items = [(r["cluster_id"], os.path.join(p.cif_dir, f'{r["cluster_id"]}.cif'),
              r["PEPTIDE_CHAIN"], r["PROTEIN_CHAIN"].split(":"),
              p.distance_cutoff, p.temperature) for r in rows
             if os.path.exists(os.path.join(p.cif_dir, f'{r["cluster_id"]}.cif'))]
    workdir = checkpoint.namespace(p.ckpt, VERSION,
                                   {"cutoff": p.distance_cutoff, "temp": p.temperature})
    results = checkpoint.run(items, worker, workdir, threads=threads,
                             id_of=lambda it: it[0], stage="multipro_prodigy")

    n = 0
    with open(snakemake.output.prodigy, "w") as fh:        # noqa: F821
        fh.write("\t".join(out_cols) + "\n")
        for r in rows:
            res = results.get(r["cluster_id"])
            if not res or res["status"] != "ok" or not res["record"]:
                continue
            rec = dict(res["record"], cluster_id=r["cluster_id"])
            fh.write("\t".join(rec[c] for c in out_cols) + "\n")
            n += 1
    print(f"DONE {n} multipro prodigy", file=sys.stderr)


if __name__ == "__main__":
    main()