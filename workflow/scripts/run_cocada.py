"""Run COCaDA on the downloaded CIF for each sampled entry.

Invocation: -c pep,prot -inter -ph 7.4 -s -m 1. Runs on the ORIGINAL CIF (not
our extracted PDB), so contacts are independent of extraction convention.

CHAIN-ID TRANSLATION (important): the rest of Propedia uses auth_asym_id chain
IDs, but COCaDA selects on label_asym_id. These diverge in ~90% of PDB entries.
Passing auth IDs straight to COCaDA silently selects the wrong chains (or none),
which is a latent defect in the published Propedia contact column. We therefore
translate each pair's auth chain ID to its POLYMER label_asym_id (built from the
CIF's _atom_site ATOM records; the mapping is 1:1 for polymers) before calling
COCaDA. Chains with no polymer label are reported explicitly, not silently zeroed.

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


def auth_to_polymer_label(cif_path):
    """Map auth_asym_id -> polymer label_asym_id from a CIF's _atom_site block.

    Only ATOM (polymer) records are used, which makes the mapping 1:1: an auth
    chain's het groups get their own label IDs, but its single polymer entity
    has exactly one label ID. Returns {auth_id: label_id}.
    """
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
                # first polymer label wins (1:1 for polymers)
                mapping.setdefault(aut, lab)
            elif in_loop and line.startswith("#"):
                break
    return mapping


def run_one(entry_id, pid, pep, prot, cif_dir, cocada_dir, ph, out_root):
    src = shard_path(cif_dir, pid)
    if not os.path.exists(src):
        return "missing_cif"
    # Absolute: COCaDA runs with cwd=cocada_dir, so a relative -o would land there.
    entry_out = os.path.abspath(os.path.join(out_root, entry_id))
    # Clear any prior output so a stale *_contacts.csv can't be picked up
    # (COCaDA names its file by PDB id, so re-runs would otherwise collide).
    if os.path.isdir(entry_out):
        for _f in os.listdir(entry_out):
            os.remove(os.path.join(entry_out, _f))
    os.makedirs(entry_out, exist_ok=True)
    tmp = tempfile.NamedTemporaryFile("w", suffix=".cif", delete=False)
    try:
        with gzip.open(src, "rt") as fh:
            tmp.write(fh.read())
        tmp.close()
        # Translate auth chain IDs (Propedia convention) -> polymer label IDs
        # (COCaDA convention). See module docstring.
        a2l = auth_to_polymer_label(tmp.name)
        pep_lab = a2l.get(pep)
        prot_lab = a2l.get(prot)
        if pep_lab is None or prot_lab is None:
            missing = [c for c, lab in ((pep, pep_lab), (prot, prot_lab)) if lab is None]
            return f"no_polymer_label:{','.join(missing)}"
        proc = subprocess.run(
            ["python3", "cocada.py", "-f", tmp.name,
             "-c", f"{pep_lab},{prot_lab}", "-inter", "-m", "1",
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