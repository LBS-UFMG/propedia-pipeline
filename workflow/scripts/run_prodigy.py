"""Throttled PRODIGY runner: bounded parallelism, per-call timeout, real error
logging. Safe to run without starving the machine."""
import csv, os, re, subprocess, sys
from concurrent.futures import ProcessPoolExecutor, as_completed

LABELS = {
    "No. of intermolecular contacts": "n_intermolecular",
    "No. of charged-charged contacts": "n_charged_charged",
    "No. of charged-polar contacts": "n_charged_polar",
    "No. of charged-apolar contacts": "n_charged_apolar",
    "No. of polar-polar contacts": "n_polar_polar",
    "No. of apolar-polar contacts": "n_apolar_polar",
    "No. of apolar-apolar contacts": "n_apolar_apolar",
    "Percentage of apolar NIS residues": "nis_apolar",
    "Percentage of charged NIS residues": "nis_charged",
    "Predicted binding affinity (kcal.mol-1)": "dg",
    "Predicted dissociation constant (M) at 25.0": "kd",
}
OUT_COLS = ["id","dg","kd","n_intermolecular","n_charged_charged","n_charged_polar",
            "n_charged_apolar","n_polar_polar","n_apolar_polar","n_apolar_apolar",
            "nis_apolar","nis_charged"]
NUM = re.compile(r"(-?\d+\.?\d*(?:e-?\d+)?)")
PDB_DIR = "results/sample/pdb"
MAX_WORKERS = 4          # <= the throttle that prevents the crash
TIMEOUT = 120


def parse(text):
    vals = {}
    for line in text.splitlines():
        for label, key in LABELS.items():
            if label in line and ":" in line:
                m = NUM.search(line.split(":", 1)[1])
                if m: vals[key] = m.group(1)
                break
    return vals


def run_one(args):
    eid, pep, prot = args
    path = os.path.join(PDB_DIR, f"{eid}.pdb")
    if not os.path.exists(path):
        return eid, None, "missing_pdb"
    try:
        p = subprocess.run(["prodigy", path, "--selection", pep, prot,
                            "--distance-cutoff", "6.0",
                            "--temperature", "25"],
                           capture_output=True, text=True, timeout=TIMEOUT)
        v = parse(p.stdout + p.stderr)
        if "dg" not in v:
            reason = "no_contacts" if "No contacts" in (p.stdout+p.stderr) else "no_dg"
            return eid, None, reason
        v["id"] = eid
        return eid, v, "ok"
    except subprocess.TimeoutExpired:
        return eid, None, "timeout"
    except Exception as exc:                     # noqa: BLE001
        return eid, None, f"exc:{exc}"


def main():
    pairs = list(csv.DictReader(open("results/sample/pairs.tsv"), delimiter="\t"))
    jobs = [(r["id"], r["pep_chain"], r["prot_chain"]) for r in pairs]
    ok = 0
    reasons = {}
    with open("results/sample/prodigy.tsv", "w") as fh, \
         open("results/sample/prodigy_errors.tsv", "w") as ef:
        fh.write("\t".join(OUT_COLS) + "\n")
        ef.write("id\treason\n")
        with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
            futs = {ex.submit(run_one, j): j[0] for j in jobs}
            for i, fut in enumerate(as_completed(futs), 1):
                eid, v, status = fut.result()
                if status == "ok":
                    fh.write("\t".join(v.get(c, "") for c in OUT_COLS) + "\n"); ok += 1
                else:
                    ef.write(f"{eid}\t{status}\n")
                    reasons[status] = reasons.get(status, 0) + 1
                if i % 50 == 0:
                    print(f"{i}/{len(jobs)} ok={ok} reasons={reasons}", file=sys.stderr)
    print(f"DONE ok={ok} reasons={reasons}", file=sys.stderr)


if __name__ == "__main__":
    main()