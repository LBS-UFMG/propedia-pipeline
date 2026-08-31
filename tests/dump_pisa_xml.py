#!/usr/bin/env python3
"""Dump the RAW PISA XML for one structure so we can see the exact schema before
extending run_pisa.py to capture the richer interface/assembly fields.

Reuses the pipeline's own PISA plumbing (run_pisa._write_cfg/_decompress/_pisa_run)
so it exercises the same code path as the real stage.

Usage (on the CCP4 box, e.g. mioglobina):
    python tests/dump_pisa_xml.py --ccp4-dir /opt/xtal/ccp4-9 --pid 1A1M
writes  pisa_1A1M_interfaces.xml  and  pisa_1A1M_assemblies.xml  to the cwd and
prints the tag names found under the first <interface> and its <molecule> elements
(a quick map of the available fields). Nothing else is written.
"""
import argparse
import gzip
import os
import sys
import tempfile
import urllib.request
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "workflow", "scripts"))
import run_pisa  # noqa: E402

RCSB = "https://files.rcsb.org/download/{pid}.cif.gz"


def _fetch(pid, cif_dir):
    dest = run_pisa.shard_path(cif_dir, pid)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    with urllib.request.urlopen(RCSB.format(pid=pid), timeout=60) as r:
        data = r.read()
    with open(dest, "wb") as fh:
        fh.write(data)
    with gzip.open(dest, "rt") as fh:
        fh.read(64)                         # sanity: really gzip
    return dest


def _child_tags(node):
    """tag -> sample text for the direct children of `node` (one level)."""
    out = {}
    for c in node:
        out.setdefault(c.tag, (c.text or "").strip()[:24])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccp4-dir", required=True)
    ap.add_argument("--pid", default="1A1M")
    ap.add_argument("--timeout", type=int, default=600)
    args = ap.parse_args()

    ccp4 = os.path.expanduser(args.ccp4_dir)
    pisa_bin = os.path.join(ccp4, "bin", "pisa")
    cfg_modelo = os.path.join(ccp4, "share", "pisa", "pisa.cfg")
    if not (os.path.exists(pisa_bin) and os.path.exists(cfg_modelo)):
        print(f"CCP4/PISA not found: need {pisa_bin} and {cfg_modelo}", file=sys.stderr)
        return 2
    os.environ.setdefault("CCP4", ccp4)
    os.environ.setdefault("CLIBD", os.path.join(ccp4, "lib", "data"))
    os.environ.setdefault("CCP4_SCR", tempfile.mkdtemp(prefix="ccp4_scr_"))

    cif_dir = tempfile.mkdtemp(prefix="pisa_dump_cif_")
    data_root = tempfile.mkdtemp(prefix="pisa_dr_")
    tmpdir = tempfile.mkdtemp(prefix="pisa_tmp_")
    session = "dump%d" % os.getpid()
    cfg = run_pisa._write_cfg(cfg_modelo, data_root, os.path.join(data_root, "pisa.cfg"))
    cif = _fetch(args.pid, cif_dir)
    struct = run_pisa._decompress(cif, tmpdir)

    rc, out = run_pisa._pisa_run([pisa_bin, session, "-analyse", struct, cfg], args.timeout)
    print(f"-analyse rc={rc}; assembly done: {'assembly analysis: done' in out}")

    for kind in ("interfaces", "assemblies"):
        xrc, xout = run_pisa._pisa_run([pisa_bin, session, "-xml", kind, cfg], args.timeout)
        path = f"pisa_{args.pid}_{kind}.xml"
        with open(path, "w") as fh:
            fh.write(xout)
        print(f"-xml {kind}: rc={xrc} -> {path} ({len(xout)} bytes)")
        if xrc == 0:
            try:
                root = ET.fromstring(xout)
                if kind == "interfaces":
                    ifc = next(root.iter("interface"), None)
                    if ifc is not None:
                        print("  <interface> child tags:", _child_tags(ifc))
                        mol = ifc.find("molecule")
                        if mol is not None:
                            print("  <molecule> child tags:", _child_tags(mol))
                else:
                    asm = next(root.iter("assembly"), None)
                    if asm is not None:
                        print("  <assembly> child tags:", _child_tags(asm))
            except ET.ParseError as e:
                print("  (parse error:", e, ")")
    print("\nSend me the two pisa_*.xml files (or the printed tag maps).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
