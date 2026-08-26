"""Run COCaDA on the downloaded CIF for each sampled entry, matching the lab's
invocation: -c pep,prot -inter -ph 7.4 -s. Runs on the ORIGINAL CIF (not our
extracted PDB), so contacts are independent of extraction convention.
Each entry writes to its own output subdir to avoid filename collisions.
Single-core per call; Snakemake owns parallelism."""
import csv
import gzip
import os
import subprocess
import sys
import tempfile


def shard_path(cif_dir, pid):
    sub = pid[1:3].lower() if len(pid) >= 3 else "_"
    return os.path.join(cif_dir, sub, f"{pid}.cif.gz")


def run_one(entry_id, pid, pep, prot, cif_dir, cocada_dir, ph, out_root):
    src = shard_path(cif_dir, pid)
    if not os.path.exists(src):
        return "missing_cif"
    # Absolute: COCaDA runs with cwd=cocada_dir, so a relative -o would land there.
    entry_out = os.path.abspath(os.path.join(out_root, entry_id))
    os.makedirs(entry_out, exist_ok=True)
    tmp = tempfile.NamedTemporaryFile("w", suffix=".cif", delete=False)
    try:
        with gzip.open(src, "rt") as fh:
            tmp.write(fh.read())
        tmp.close()
        proc = subprocess.run(
            ["python3", "cocada.py", "-f", tmp.name,
             "-c", f"{pep},{prot}", "-inter", "-m", "1",
             "-o", entry_out, "-ph", str(ph), "-s"],
            cwd=cocada_dir, capture_output=True, text=True, timeout=300)
        produced = [f for f in os.listdir(entry_out) if f.endswith("_contacts.csv")]
        if not produced:
            return "no_output"
        canonical = os.path.join(entry_out, f"{entry_id}_contacts.csv")
        src_csv = os.path.join(entry_out, produced[0])
        if src_csv != canonical:
            os.replace(src_csv, canonical)
        with open(canonical) as fh:
            n = sum(1 for _ in fh) - 1     # minus header
        return n
    except subprocess.TimeoutExpired:
        return "timeout"
    except Exception as exc:                       # noqa: BLE001
        return f"exc:{exc}"
    finally:
        os.unlink(tmp.name)


def main():
    p = snakemake.params                                   # noqa: F821
    cocada_dir = os.path.expanduser(p.cocada_dir)
    out_root = p.out_root
    os.makedirs(out_root, exist_ok=True)
    pairs = list(csv.DictReader(open(snakemake.input.pairs), delimiter="\t"))  # noqa: F821

    ok = 0
    reasons = {}
    with open(snakemake.output.summary, "w") as summ:      # noqa: F821
        summ.write("id\tpdb\tpep\tprot\ttotal_contacts\n")
        for i, row in enumerate(pairs, 1):
            eid = row["id"]
            res = run_one(eid, row["pdb"], row["pep_chain"], row["prot_chain"],
                          p.cif_dir, cocada_dir, p.ph, out_root)
            if isinstance(res, int):
                ok += 1
                summ.write(f"{eid}\t{row['pdb']}\t{row['pep_chain']}\t"
                           f"{row['prot_chain']}\t{res}\n")
            else:
                reasons[res] = reasons.get(res, 0) + 1
            if i % 50 == 0:
                print(f"{i}/{len(pairs)} ok={ok} reasons={reasons}", file=sys.stderr)
    print(f"DONE ok={ok} reasons={reasons}", file=sys.stderr)


if __name__ == "__main__":
    main()