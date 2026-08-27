"""Throttled PRODIGY runner: bounded parallelism, per-call timeout, real error
logging. Paths, distance cutoff and temperature come from Snakemake, so this
works unchanged in both sample and full mode. The worker/timeout throttle stays
local (it's the safeguard that prevented the OOM/IO crash), tune if needed."""
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
OUT_COLS = ["id", "dg", "kd", "n_intermolecular", "n_charged_charged",
            "n_charged_polar", "n_charged_apolar", "n_polar_polar",
            "n_apolar_polar", "n_apolar_apolar", "nis_apolar", "nis_charged"]
NUM = re.compile(r"(-?\d+\.?\d*(?:e-?\d+)?)")

MAX_WORKERS = 4          # throttle that prevents the OOM/IO crash; tune if needed
TIMEOUT = 120            # per-call seconds


def parse(text):
    vals = {}
    for line in text.splitlines():
        for label, key in LABELS.items():
            if label in line and ":" in line:
                m = NUM.search(line.split(":", 1)[1])
                if m:
                    vals[key] = m.group(1)
                break
    return vals


def run_one(args):
    # everything run_one needs is passed in — no module globals, no snakemake,
    # so it is safe under both 'fork' and 'spawn' process-pool start methods.
    eid, pep, prot, cif_dir, cutoff, temperature = args
    path = os.path.join(cif_dir, f"{eid}.cif")
    if not os.path.exists(path):
        return eid, None, "missing_cif"
    try:
        proc = subprocess.run(["prodigy", path, "--selection", pep, prot,
                               "--distance-cutoff", str(cutoff),
                               "--temperature", str(temperature)],
                              capture_output=True, text=True, timeout=TIMEOUT)
        out = proc.stdout + proc.stderr
        v = parse(out)
        if "dg" not in v:
            return eid, None, "no_contacts" if "No contacts" in out else "no_dg"
        v["id"] = eid
        return eid, v, "ok"
    except subprocess.TimeoutExpired:
        return eid, None, "timeout"
    except Exception as exc:                     # noqa: BLE001
        return eid, None, f"exc:{exc}"


def main():
    pairs_path = snakemake.input.pairs           # noqa: F821
    out_path   = snakemake.output.prodigy        # noqa: F821
    err_path   = snakemake.output.errors         # noqa: F821
    cif_dir     = snakemake.params.cif_dir          # noqa: F821
    cutoff      = snakemake.params.distance_cutoff  # noqa: F821
    temperature = snakemake.params.temperature      # noqa: F821

    pairs = list(csv.DictReader(open(pairs_path), delimiter="\t"))
    jobs = [(r["id"], r["pep_chain"], r["prot_chain"], cif_dir, cutoff, temperature)
            for r in pairs]

    ok, reasons = 0, {}
    with open(out_path, "w") as fh, open(err_path, "w") as ef:
        fh.write("\t".join(OUT_COLS) + "\n")
        ef.write("id\treason\n")
        with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
            futs = {ex.submit(run_one, j): j[0] for j in jobs}
            for i, fut in enumerate(as_completed(futs), 1):
                eid, v, status = fut.result()
                if status == "ok":
                    fh.write("\t".join(v.get(c, "") for c in OUT_COLS) + "\n")
                    ok += 1
                else:
                    ef.write(f"{eid}\t{status}\n")
                    reasons[status] = reasons.get(status, 0) + 1
                if i % 50 == 0:
                    print(f"{i}/{len(jobs)} ok={ok} reasons={reasons}", file=sys.stderr)
    print(f"DONE ok={ok} reasons={reasons}", file=sys.stderr)


if __name__ == "__main__":
    main()