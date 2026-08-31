#!/usr/bin/env python3
"""Verify the EXTENDED PISA extraction on one structure before the full re-run.

Runs the updated run_pisa.run_one on 1A1M, applies the same crossing +
assembly-matching the stage does, and prints the pep-pro row's richer PISA fields
next to propedia26's stored values for 1A1M-C-A. A match means the parser is right
and it's safe to bump the stage and re-run PISA on the full build.

Usage (on the CCP4 box):
    python tests/verify_pisa_richer.py --ccp4-dir /opt/xtal/ccp4-9
"""
import argparse
import gzip
import os
import sys
import tempfile
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "workflow", "scripts"))
import run_pisa  # noqa: E402

# propedia26's stored values for 1A1M-C-A (peptide C = _1, protein A = _2)
EXPECT = {
    "pisa_nres_1": "9", "pisa_natoms_1": "67", "pisa_area_1": "977", "pisa_solv_en_1": "-6.142",
    "pisa_nres_2": "35", "pisa_natoms_2": "122", "pisa_area_2": "665", "pisa_solv_en_2": "-3.243",
    "pisa_diss_energy": "8.067", "pisa_entropy": "7.39", "pisa_int_energy": "-9.384",
    "pisa_asa": "15346", "pisa_bsa": "1643", "pisa_diss_area": "821",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccp4-dir", required=True)
    ap.add_argument("--pid", default="1A1M")
    ap.add_argument("--pep", default="C")
    ap.add_argument("--prot", default="A")
    args = ap.parse_args()

    ccp4 = os.path.expanduser(args.ccp4_dir)
    pisa_bin = os.path.join(ccp4, "bin", "pisa")
    cfg = os.path.join(ccp4, "share", "pisa", "pisa.cfg")
    if not (os.path.exists(pisa_bin) and os.path.exists(cfg)):
        print(f"CCP4/PISA not found: need {pisa_bin} and {cfg}", file=sys.stderr)
        return 2
    os.environ.setdefault("CCP4", ccp4)
    os.environ.setdefault("CLIBD", os.path.join(ccp4, "lib", "data"))
    os.environ.setdefault("CCP4_SCR", tempfile.mkdtemp(prefix="ccp4_scr_"))

    cif_dir = tempfile.mkdtemp(prefix="pisa_verify_")
    dest = run_pisa.shard_path(cif_dir, args.pid)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    with urllib.request.urlopen(
            f"https://files.rcsb.org/download/{args.pid}.cif.gz", timeout=60) as r:
        open(dest, "wb").write(r.read())
    with gzip.open(dest, "rt") as fh:
        fh.read(64)

    rec, status = run_pisa.run_one((args.pid, dest, ccp4, cfg, 600))
    print(f"status={status}  n_interfaces={(rec or {}).get('n_interfaces')}  "
          f"n_assemblies={len((rec or {}).get('assemblies') or [])}")
    if status != "ok":
        return 1

    iface = run_pisa.cross(args.pep, args.prot, rec)
    asm = run_pisa.find_assembly(rec.get("assemblies"), args.pep, args.prot)
    row = {}
    if iface:
        for k in run_pisa.IFACE_FIELDS:
            row["pisa_" + k] = iface.get(k, "")
    if asm:
        for k in run_pisa.ASM_FIELDS:
            row["pisa_" + k] = asm.get(k, "")

    print(f"\n{'field':20} {'ours':>16}   {'propedia26':>12}   match?")
    print("-" * 60)
    ok = 0
    for k, exp in EXPECT.items():
        got = row.get(k, "")
        try:                                    # numeric compare, 1% tol (version drift)
            hit = abs(float(got) - float(exp)) <= max(0.01, abs(float(exp)) * 0.01)
        except ValueError:
            hit = str(got) == str(exp)
        ok += hit
        print(f"{k:20} {got:>16}   {exp:>12}   {'OK' if hit else 'DIFF'}")
    print("-" * 60)
    print(f"{ok}/{len(EXPECT)} within 1% of propedia26 "
          f"(diss_energy is expected to drift a little — PISA version).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
