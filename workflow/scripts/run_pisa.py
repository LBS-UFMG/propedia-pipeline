"""PISA interface analysis (CCP4 `pisa`), per-entry checkpointed, then crossed to
the pep-pro pairs.

Runs on the COMPLETE mmCIF (our cif_dir has full structures with _cell/_symmetry,
which PISA needs — the extracted per-pair CIFs lack CRYST1 and can't be used).
X-RAY only: cryo-EM/NMR have no crystal lattice, so PISA can't assess assemblies
there; those are filtered out via our own metadata (STRUCTURE_METHOD) and left
blank. The PISA-running logic (analyse -> xml interfaces -> parse) is adapted from
the lab's pisa_css.py.

Crossing rule (validated on 1WRZ-B-A vs the live site): for a pep-pro pair, take
the PISA interface whose unordered chain set == {pep_chain, prot_chain} with the
LARGEST area (PISA sorts by area, so the smallest matching interface_id). Ion/
ligand "interfaces" (chain like '[CA]A:151') don't match a clean chain set and are
excluded for free.

Output: pisa.tsv keyed by pair id, with named pisa_* columns.
"""
import csv
import gzip
import itertools
import os
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

import checkpoint

VERSION = "1"
XRAY = "X-RAY DIFFRACTION"
_CTR = itertools.count()

# fields captured per interface (mirrors the lab pisa_css.py extraction)
IFACE_FIELDS = ["interface_id", "tipo", "chain_1", "chain_2", "css", "area",
                "solv_en", "pvalue", "n_hbonds", "n_saltbridges"]
# pair-level output columns
OUT_COLS = ["id", "pisa_status", "pisa_assembly_done", "pisa_n_interfaces",
            "pisa_interface_id", "pisa_chain_1", "pisa_chain_2", "pisa_css",
            "pisa_area", "pisa_solv_en", "pisa_pvalue", "pisa_tipo",
            "pisa_n_hbonds", "pisa_n_saltbridges"]


def shard_path(cif_dir, pid):
    sub = pid[1:3].lower() if len(pid) >= 3 else "_"
    return os.path.join(cif_dir, sub, f"{pid}.cif.gz")


def _texto(node, *tags):
    for tag in tags:
        x = node.find(tag) if node is not None else None
        if x is not None and x.text:
            return x.text.strip()
    return ""


def _write_cfg(template, data_root, dest):
    """Copy pisa.cfg replacing DATA_ROOT, to isolate parallel sessions."""
    lines = open(template, encoding="utf-8", errors="replace").read().splitlines()
    out_lines, swap = [], False
    for l in lines:
        if swap and l.strip() and not l.lstrip().startswith("#"):
            out_lines.append(data_root); swap = False; continue
        out_lines.append(l)
        if l.strip() == "DATA_ROOT":
            swap = True
    open(dest, "w", encoding="utf-8").write("\n".join(out_lines) + "\n")
    return dest


def _decompress(src, folder):
    base = os.path.basename(src)
    if base.endswith(".gz"):
        target = os.path.join(folder, base[:-3])
        with gzip.open(src, "rb") as fin, open(target, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        return target
    return src


def run_one(args):
    """worker(item) -> (interfaces|None, status). Runs PISA on one full structure and
    returns every interface. Timeouts/transient errors use 'retry:' so they re-run."""
    pid, cif_path, ccp4_dir, cfg_modelo, timeout = args
    if not os.path.exists(cif_path):
        return None, "missing_cif"
    pisa_bin = os.path.join(ccp4_dir, "bin", "pisa")
    data_root = tempfile.mkdtemp(prefix="pisa_dr_")
    tmpdir = tempfile.mkdtemp(prefix="pisa_tmp_")
    session = "s%d_%d" % (os.getpid(), next(_CTR))
    try:
        cfg = _write_cfg(cfg_modelo, data_root, os.path.join(data_root, "pisa.cfg"))
        struct = _decompress(cif_path, tmpdir)
        r = subprocess.run([pisa_bin, session, "-analyse", struct, cfg],
                           capture_output=True, text=True, timeout=timeout)
        out = (r.stdout or "") + (r.stderr or "")
        assembly = "yes" if "assembly analysis: done" in out else "no"
        if r.returncode != 0 or "quit" in out:
            return {"assembly_done": assembly}, "error"
        x = subprocess.run([pisa_bin, session, "-xml", "interfaces", cfg],
                           capture_output=True, text=True, timeout=timeout)
        if x.returncode != 0:
            return {"assembly_done": assembly}, "error"
        root = ET.fromstring(x.stdout)
        interfaces = list(root.iter("interface"))
        rows = []
        for i in interfaces:
            mols = i.findall("molecule")
            ch = [_texto(m, "chain_id") for m in mols[:2]]
            rows.append({
                "interface_id": _texto(i, "id"), "tipo": _texto(i, "type"),
                "chain_1": ch[0] if ch else "", "chain_2": ch[1] if len(ch) > 1 else "",
                "css": _texto(i, "css"), "area": _texto(i, "int_area"),
                "solv_en": _texto(i, "int_solv_en"), "pvalue": _texto(i, "pvalue"),
                "n_hbonds": _texto(i.find("h-bonds"), "n_bonds") if i.find("h-bonds") is not None else "",
                "n_saltbridges": _texto(i.find("salt-bridges"), "n_bonds") if i.find("salt-bridges") is not None else "",
            })
        return {"assembly_done": assembly, "n_interfaces": len(interfaces),
                "interfaces": rows}, ("ok" if rows else "no_interface")
    except subprocess.TimeoutExpired:
        return None, "retry:timeout"
    except ET.ParseError as e:
        return None, "xml_invalid:%s" % str(e)[:120]
    except Exception as e:                                 # noqa: BLE001
        return None, "retry:exc:%s" % str(e)[:120]
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
        for name in os.listdir(data_root):
            if name.endswith(session):
                shutil.rmtree(os.path.join(data_root, name), ignore_errors=True)
        shutil.rmtree(data_root, ignore_errors=True)


def _fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return -1.0


def cross(pep, prot, record):
    """Pick the pep-pro interface: chain set == {pep,prot}, largest area."""
    want = {pep, prot}
    matches = [f for f in (record or {}).get("interfaces", [])
               if {f.get("chain_1"), f.get("chain_2")} == want]
    if not matches:
        return None
    return max(matches, key=lambda f: _fnum(f.get("area")))


def main():
    p = snakemake.params                                   # noqa: F821
    ccp4_dir = os.path.expanduser(p.ccp4_dir)
    cfg_modelo = os.path.join(ccp4_dir, "share", "pisa", "pisa.cfg")
    timeout = getattr(p, "timeout", 900)
    threads = getattr(snakemake, "threads", 1) or 1        # noqa: F821

    # CCP4 environment (inherited by the fork pool workers)
    os.environ.setdefault("CCP4", ccp4_dir)
    os.environ.setdefault("CLIBD", os.path.join(ccp4_dir, "lib", "data"))
    os.environ.setdefault("CCP4_SCR", tempfile.mkdtemp(prefix="ccp4_scr_"))

    pairs = list(csv.DictReader(open(snakemake.input.pairs), delimiter="\t"))  # noqa: F821
    method = {r["id"]: r.get("STRUCTURE_METHOD", "")
              for r in csv.DictReader(open(snakemake.input.metadata), delimiter="\t")}  # noqa: F821
    # X-ray PDBs only (per our metadata); one entry per unique PDB id
    xray_pdbs = sorted({r["pdb"] for r in pairs if method.get(r["id"]) == XRAY})
    items = [(pid, shard_path(p.cif_dir, pid), ccp4_dir, cfg_modelo, timeout)
             for pid in xray_pdbs]
    workdir = checkpoint.namespace(p.ckpt, VERSION, {"ccp4": ccp4_dir})
    results = checkpoint.run(items, run_one, workdir, threads=threads,
                             id_of=lambda it: it[0], stage="pisa")

    n = ok = 0
    with open(snakemake.output.pisa, "w") as fh:           # noqa: F821
        fh.write("\t".join(OUT_COLS) + "\n")
        for r in pairs:
            eid, pdb = r["id"], r["pdb"]
            res = results.get(pdb)
            rec = res["record"] if res else None
            status = res["status"] if res else ("not_xray" if method.get(eid) != XRAY else "missing")
            iface = cross(r["pep_chain"], r["prot_chain"], rec) if status == "ok" else None
            row = {c: "" for c in OUT_COLS}
            row["id"] = eid
            row["pisa_status"] = status
            row["pisa_assembly_done"] = (rec or {}).get("assembly_done", "")
            row["pisa_n_interfaces"] = (rec or {}).get("n_interfaces", "")
            if iface:
                for k in IFACE_FIELDS:
                    row["pisa_" + k] = iface.get(k, "")
                ok += 1
            fh.write("\t".join(str(row[c]) for c in OUT_COLS) + "\n")
            n += 1
    print(f"DONE pisa: {n} pairs, {ok} with a matched interface "
          f"({len(xray_pdbs)} X-ray PDBs analysed)", file=sys.stderr)


if __name__ == "__main__":
    main()
