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
import signal
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

import checkpoint

VERSION = "2"            # bumped: extraction now also captures per-molecule interface
                         # stats + assembly energetics -> checkpoints must recompute
XRAY = "X-RAY DIFFRACTION"
_CTR = itertools.count()

# fields captured per interface. The _1/_2 pairs are the two interacting molecules'
# per-molecule interface stats, in the same order as chain_1/chain_2.
IFACE_FIELDS = ["interface_id", "tipo", "chain_1", "chain_2", "css", "area",
                "solv_en", "pvalue", "n_hbonds", "n_saltbridges",
                "nres_1", "natoms_1", "area_1", "solv_en_1",
                "nres_2", "natoms_2", "area_2", "solv_en_2"]
# assembly-level energetics, taken from the assembly whose composition == the pep-pro
# chain set (see find_assembly). PISA_bsa here is the assembly's buried area; note the
# per-molecule area_1+area_2 also sums to it for a 2-body assembly.
ASM_FIELDS = ["diss_energy", "entropy", "int_energy", "asa", "bsa", "diss_area"]
# pair-level output columns.  pisa_interface_class is the DERIVED biological-vs-crystal
# annotation (the reviewer-requested read-out): biological / crystal-packing when PISA
# gave a CSS for the pep-pro interface, indeterminate when it ran but could not score
# this interface, not_applicable for non-X-ray methods (no crystal lattice => the
# crystal-packing question does not arise), blank when PISA could not be assessed.
OUT_COLS = ["id", "pisa_status", "pisa_interface_class",
            "pisa_assembly_done", "pisa_n_interfaces",
            "pisa_interface_id", "pisa_chain_1", "pisa_chain_2", "pisa_css",
            "pisa_area", "pisa_solv_en", "pisa_pvalue", "pisa_tipo",
            "pisa_n_hbonds", "pisa_n_saltbridges",
            "pisa_nres_1", "pisa_natoms_1", "pisa_area_1", "pisa_solv_en_1",
            "pisa_nres_2", "pisa_natoms_2", "pisa_area_2", "pisa_solv_en_2",
            "pisa_diss_energy", "pisa_entropy", "pisa_int_energy",
            "pisa_asa", "pisa_bsa", "pisa_diss_area"]


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


class _Timeout(Exception):
    pass


def _pisa_run(cmd, timeout):
    """Run a pisa command in its own process group and, on timeout, SIGKILL the whole
    group (CCP4 pisa spawns children that outlive a plain kill -> orphan pile-up)."""
    # errors="replace": PISA emits occasional non-UTF-8 bytes (e.g. 0x81, CP-1252) in
    # text fields; strict decoding raised UnicodeDecodeError -> caught as a transient
    # 'retry:exc' that never self-heals (deterministic) and silently blanked PISA for
    # ~3.7k X-ray entries. A replaced byte only ever lands in a text field, never in the
    # numeric interface/assembly values we parse.
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, encoding="utf-8", errors="replace",
                            start_new_session=True)
    try:
        out, err = proc.communicate(timeout=timeout)
        return proc.returncode, (out or "") + (err or "")
    except subprocess.TimeoutExpired:
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            pass
        proc.communicate()
        raise _Timeout()


_PISA_KEEP = {"_entry", "_cell", "_symmetry", "_atom_sites", "_atom_site", "_atom_type"}


def _minimal_cif(path):
    """Coordinates+crystal-only mmCIF for PISA. CCP4 PISA's mmCIF reader chokes on some
    full-RCSB metadata loops (observed: `_struct_conf` HELX_P rows -> rc=19 'problem
    reading the coordinate file'), which silently blanked PISA for ~2.4k X-ray entries.
    PISA only needs cell/symmetry/atom_site, so keep just those categories — verified
    numerically identical (CSS unchanged on 1A1M). Used only as a fallback after a
    first-try parse failure, so the common path is never touched."""
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    out, i, n = [], 0, len(lines)
    while i < n:
        s = lines[i].strip()
        if s.startswith("data_") or s == "" or s == "#":
            out.append(lines[i]); i += 1; continue
        if s == "loop_":
            j = i + 1; hdr = []
            while j < n and lines[j].strip().startswith("_"):
                hdr.append(lines[j]); j += 1
            cat = hdr[0].strip().split(".")[0] if hdr else ""
            k = j
            while (k < n and lines[k].strip() != "loop_"
                   and not lines[k].strip().startswith("_")
                   and not lines[k].strip().startswith("data_")
                   and lines[k].strip() != "#"):
                k += 1
            if cat in _PISA_KEEP:
                out.append("loop_"); out += hdr; out += lines[j:k]
            i = k; continue
        if s.startswith("_"):
            cat = s.split(".")[0]; blk = [lines[i]]; j = i + 1
            if len(s.split()) == 1 and j < n and lines[j].startswith(";"):
                blk.append(lines[j]); j += 1
                while j < n and not lines[j].startswith(";"):
                    blk.append(lines[j]); j += 1
                if j < n:
                    blk.append(lines[j]); j += 1
            if cat in _PISA_KEEP:
                out += blk
            i = j; continue
        i += 1
    return "\n".join(out) + "\n"


def run_one(args):
    """worker(item) -> (interfaces|None, status). Runs PISA on one full structure and
    returns every interface. A timeout is a PERSISTED status (recorded, not retried) so
    a pathological structure never wedges future runs."""
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
        rc, out = _pisa_run([pisa_bin, session, "-analyse", struct, cfg], timeout)
        if rc != 0 or "quit" in out:
            # PISA's mmCIF reader trips on some full-RCSB metadata loops (e.g.
            # _struct_conf) -> retry on a coordinates+crystal-only CIF (identical
            # numbers). Only on failure, so the common path pays nothing.
            try:
                mini = os.path.join(tmpdir, "min.cif")
                with open(mini, "w", encoding="utf-8") as fh:
                    fh.write(_minimal_cif(struct))
                rc, out = _pisa_run([pisa_bin, session, "-analyse", mini, cfg], timeout)
            except Exception:                                  # noqa: BLE001
                pass
        assembly = "yes" if "assembly analysis: done" in out else "no"
        if rc != 0 or "quit" in out:
            return {"assembly_done": assembly}, "error"
        xrc, xout = _pisa_run([pisa_bin, session, "-xml", "interfaces", cfg], timeout)
        if xrc != 0:
            return {"assembly_done": assembly}, "error"
        root = ET.fromstring(xout)
        interfaces = list(root.iter("interface"))
        rows = []
        for i in interfaces:
            mols = i.findall("molecule")
            ch = [_texto(m, "chain_id") for m in mols[:2]]
            m1, m2 = (mols[0] if mols else None), (mols[1] if len(mols) > 1 else None)
            rows.append({
                "interface_id": _texto(i, "id"), "tipo": _texto(i, "type"),
                "chain_1": ch[0] if ch else "", "chain_2": ch[1] if len(ch) > 1 else "",
                "css": _texto(i, "css"), "area": _texto(i, "int_area"),
                "solv_en": _texto(i, "int_solv_en"), "pvalue": _texto(i, "pvalue"),
                "n_hbonds": _texto(i.find("h-bonds"), "n_bonds") if i.find("h-bonds") is not None else "",
                "n_saltbridges": _texto(i.find("salt-bridges"), "n_bonds") if i.find("salt-bridges") is not None else "",
                # per-molecule interface stats (mol order == chain_1/chain_2 order)
                "nres_1": _texto(m1, "int_nres"), "natoms_1": _texto(m1, "int_natoms"),
                "area_1": _texto(m1, "int_area"), "solv_en_1": _texto(m1, "int_solv_en"),
                "nres_2": _texto(m2, "int_nres"), "natoms_2": _texto(m2, "int_natoms"),
                "area_2": _texto(m2, "int_area"), "solv_en_2": _texto(m2, "int_solv_en"),
            })
        # assembly-level energetics: run the assembly xml too (already computed by
        # -analyse) and keep every assembly; the pep-pro row later picks the assembly
        # whose composition == {pep_chain, prot_chain} (validated on 1A1M-C-A).
        assemblies = []
        arc, aout = _pisa_run([pisa_bin, session, "-xml", "assemblies", cfg], timeout)
        if arc == 0:
            try:
                aroot = ET.fromstring(aout)
                for a in aroot.iter("assembly"):
                    assemblies.append({
                        "composition": _texto(a, "composition"),
                        "diss_energy": _texto(a, "diss_energy"),
                        "entropy": _texto(a, "entropy"),
                        "asa": _texto(a, "asa"), "bsa": _texto(a, "bsa"),
                        "diss_area": _texto(a, "diss_area"),
                        "int_energy": _texto(a, "int_energy"),
                    })
            except ET.ParseError:
                assemblies = []          # richer fields blank; interface fields still fine
        return {"assembly_done": assembly, "n_interfaces": len(interfaces),
                "interfaces": rows, "assemblies": assemblies}, ("ok" if rows else "no_interface")
    except _Timeout:
        return None, "timeout"          # persisted: recorded, never retried/re-hung
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


def find_assembly(assemblies, pep, prot):
    """Pick the assembly for the pep-pro pair: the one whose chain composition equals
    {pep, prot}. PISA's `composition` concatenates chain ids (e.g. 'AC'); matched by
    sorted characters, which is exact for single-char chains (>99% of the PDB) and
    degrades to no-match (blank richer fields) for the rare multi-char-chain case.
    Returns the assembly dict, or None."""
    want = "".join(sorted(pep + prot))
    for a in assemblies or []:
        if "".join(sorted(a.get("composition", ""))) == want:
            return a
    return None


def _css_val(iface):
    """CSS as float, or None if PISA emitted no parseable score for this interface."""
    try:
        return float((iface or {}).get("css", ""))
    except (TypeError, ValueError):
        return None


def classify(status, iface, threshold):
    """Derived biological-vs-crystal label for the pep-pro interface.

    biological / crystal-packing  -> PISA scored this interface (css >=/< threshold)
    indeterminate                 -> PISA ran but could not score this pep-pro interface
    not_applicable                -> non-X-ray method: no crystal lattice, so there is no
                                     crystal-packing artifact to flag
    "" (blank)                    -> PISA could not be assessed (missing CIF / error / timeout)
    """
    if status == "not_xray":
        return "not_applicable"
    if status != "ok":
        return ""                       # missing / error / timeout / xml_invalid
    css = _css_val(iface)
    if css is None:
        return "indeterminate"          # ran, but no CSS for the matched interface (or none matched)
    return "biological" if css >= threshold else "crystal-packing"


def _write_blank(path, pairs, status):
    """Emit a well-formed pisa.tsv where every pair carries `status` and blank PISA
    fields. Lets the DAG complete (propedia.csv just has empty PISA columns) when PISA
    cannot run, instead of erroring or perpetually retrying."""
    with open(path, "w") as fh:
        fh.write("\t".join(OUT_COLS) + "\n")
        for r in pairs:
            row = {c: "" for c in OUT_COLS}
            row["id"] = r["id"]
            row["pisa_status"] = status
            fh.write("\t".join(str(row[c]) for c in OUT_COLS) + "\n")


def main():
    p = snakemake.params                                   # noqa: F821
    ccp4_dir = os.path.expanduser(p.ccp4_dir)
    cfg_modelo = os.path.join(ccp4_dir, "share", "pisa", "pisa.cfg")
    pisa_bin = os.path.join(ccp4_dir, "bin", "pisa")
    on_missing = getattr(p, "on_missing", "skip")   # "skip" (turnkey) | "error" (enforce)
    timeout = getattr(p, "timeout", 600)   # per pisa call; a timed-out structure is
                                           # recorded (blank PISA), not retried
    css_threshold = float(getattr(p, "css_biological_threshold", 0.5))
    threads = getattr(snakemake, "threads", 1) or 1        # noqa: F821

    pairs = list(csv.DictReader(open(snakemake.input.pairs), delimiter="\t"))  # noqa: F821

    # Turnkey degradation: CCP4/PISA is a heavy external dependency. If it is absent,
    # don't wedge the whole build (or perpetually retry a missing binary) — either skip
    # PISA with a clear per-entry status, or hard-fail if the run declared it required.
    if not (os.path.exists(pisa_bin) and os.path.exists(cfg_modelo)):
        msg = (f"PISA unavailable: expected `pisa` at {pisa_bin} and config at "
               f"{cfg_modelo}. Install CCP4 and set machine.ccp4_dir to enable the "
               f"biological-vs-crystal interface annotation.")
        if on_missing == "error":
            sys.exit("ERROR: " + msg + " (pisa.on_missing='error')")
        print("=" * 78 + f"\nWARNING: {msg}\n"
              "Proceeding WITHOUT PISA — every entry gets pisa_status='pisa_unavailable' "
              "and blank PISA columns (pisa.on_missing='skip').\n" + "=" * 78,
              file=sys.stderr)
        _write_blank(snakemake.output.pisa, pairs, "pisa_unavailable")  # noqa: F821
        return

    # CCP4 environment (inherited by the fork pool workers)
    os.environ.setdefault("CCP4", ccp4_dir)
    os.environ.setdefault("CLIBD", os.path.join(ccp4_dir, "lib", "data"))
    os.environ.setdefault("CCP4_SCR", tempfile.mkdtemp(prefix="ccp4_scr_"))

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
    classes = {}
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
            if status == "ok":
                asm = find_assembly((rec or {}).get("assemblies"),
                                    r["pep_chain"], r["prot_chain"])
                if asm:
                    for k in ASM_FIELDS:
                        row["pisa_" + k] = asm.get(k, "")
            row["pisa_interface_class"] = classify(status, iface, css_threshold)
            classes[row["pisa_interface_class"] or "blank"] = \
                classes.get(row["pisa_interface_class"] or "blank", 0) + 1
            fh.write("\t".join(str(row[c]) for c in OUT_COLS) + "\n")
            n += 1
    breakdown = ", ".join(f"{k}={classes[k]}" for k in sorted(classes))
    print(f"DONE pisa: {n} pairs, {ok} with a matched interface "
          f"({len(xray_pdbs)} X-ray PDBs analysed); "
          f"interface_class[css>={css_threshold}]: {breakdown}", file=sys.stderr)


if __name__ == "__main__":
    main()
