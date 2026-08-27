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


def main():
    p = snakemake.params                                   # noqa: F821
    rows = list(csv.DictReader(open(snakemake.input.multipro), delimiter="\t"))  # noqa: F821
    out_cols = ["cluster_id"] + [c for c, _ in COLMAP]
    n = 0
    with open(snakemake.output.prodigy, "w") as fh:        # noqa: F821
        fh.write("\t".join(out_cols) + "\n")
        for r in rows:
            cid = r["cluster_id"]
            path = os.path.join(p.cif_dir, f"{cid}.cif")
            if not os.path.exists(path):
                continue
            pep = r["PEPTIDE_CHAIN"]
            prots = r["PROTEIN_CHAIN"].split(":")
            per_chain = [run_chain(path, pep, pc, p.distance_cutoff, p.temperature)
                         for pc in prots]
            rec = {"cluster_id": cid}
            for col, key in COLMAP:
                rec[col] = ":".join(str(v.get(key, "")) for v in per_chain)
            fh.write("\t".join(rec[c] for c in out_cols) + "\n")
            n += 1
            if n % 20 == 0:
                print(f"{n} multipro prodigy", file=sys.stderr)
    print(f"DONE {n} multipro prodigy", file=sys.stderr)


if __name__ == "__main__":
    main()