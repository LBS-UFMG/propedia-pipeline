"""COCaDA on the multi-chain complex (Multipro), per-entry checkpointed.

Runs on the ORIGINAL downloaded CIF (not the MMCIFIO-extracted multipro CIF):
COCaDA cannot parse Biopython MMCIFIO output ('NoneType' object is not iterable),
so — exactly like the pep-pro `run_cocada.py` — we feed COCaDA the original RCSB CIF
and select the Multipro chain set (peptide + its protein chains). Chain IDs are
translated auth_asym_id -> polymer label_asym_id (COCaDA selects on label IDs),
then COCaDA analyses all interchain contacts among the selected chains (`-inter`).

Output: one <PDB>_contacts.csv per Multipro entry under out_root/<cluster_id>/
(same format the website serves). The packager copies these into
multipro/contacts/<cluster_id>/<PDB>_contacts.csv.
"""
import csv
import gzip
import os
import subprocess
import sys
import tempfile

import checkpoint

VERSION = "2"            # bump: run on the ORIGINAL CIF with -c (MMCIFIO CIF broke COCaDA)


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
    cluster_id, pdb, pep, prots, cif_dir, cocada_dir, ph, out_root = args
    src = shard_path(cif_dir, pdb)
    if not os.path.exists(src):
        return None, "missing_cif"
    entry_out = os.path.abspath(os.path.join(out_root, cluster_id))
    if os.path.isdir(entry_out):                      # clear any stale output
        for _f in os.listdir(entry_out):
            os.remove(os.path.join(entry_out, _f))
    os.makedirs(entry_out, exist_ok=True)
    tmp = tempfile.NamedTemporaryFile("w", suffix=".cif", delete=False)
    try:
        with gzip.open(src, "rt") as fh:
            tmp.write(fh.read())
        tmp.close()
        a2l = auth_to_polymer_label(tmp.name)
        labels, missing, seen = [], [], set()
        for c in [pep] + [x for x in prots if x]:     # peptide first, then protein chains
            lab = a2l.get(c)
            if lab is None:
                missing.append(c)
            elif lab not in seen:
                seen.add(lab)
                labels.append(lab)
        if missing:
            return None, f"no_polymer_label:{','.join(sorted(set(missing)))}"
        if len(labels) < 2:                           # need >=2 chains for interchain
            return None, "too_few_chains"
        proc = subprocess.run(
            ["python3", "cocada.py", "-f", tmp.name,
             "-c", ",".join(labels), "-inter", "-m", "1",
             "-o", entry_out, "-ph", str(ph), "-s"],
            cwd=cocada_dir, capture_output=True, text=True, timeout=600)
        produced = [f for f in os.listdir(entry_out) if f.endswith("_contacts.csv")]
        if not produced:
            err = (proc.stderr or proc.stdout or "").strip().splitlines()
            return None, "no_output:" + (err[-1][:100] if err else "")
        canonical = os.path.join(entry_out, f"{pdb}_contacts.csv")
        src_csv = os.path.join(entry_out, produced[0])
        if src_csv != canonical:
            os.replace(src_csv, canonical)
        with open(canonical) as fh:
            n = sum(1 for _ in fh) - 1                 # minus header
        return {"contacts": n}, "ok"
    except subprocess.TimeoutExpired:
        return None, "retry:timeout"
    except Exception as exc:                           # noqa: BLE001
        return None, f"retry:exc:{exc}"
    finally:
        os.unlink(tmp.name)


def main():
    p = snakemake.params                                   # noqa: F821
    cocada_dir = os.path.expanduser(p.cocada_dir)
    out_root = os.path.abspath(p.out_root)
    threads = getattr(snakemake, "threads", 1) or 1        # noqa: F821
    os.makedirs(out_root, exist_ok=True)
    rows = list(csv.DictReader(open(snakemake.input.multipro), delimiter="\t"))  # noqa: F821

    items = [(r["cluster_id"], r["PDB_ID"], r["PEPTIDE_CHAIN"],
              r["PROTEIN_CHAIN"].split(":"), p.cif_dir, cocada_dir, p.ph, out_root)
             for r in rows]
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
                # collapse parameterized statuses (no_output:<msg>, no_polymer_label:<c>)
                key = st.split(":", 1)[0] if ":" in st else st
                reasons[key] = reasons.get(key, 0) + 1
    print(f"DONE multipro cocada: ok={ok} reasons={reasons}", file=sys.stderr)


if __name__ == "__main__":
    main()
