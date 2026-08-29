#!/usr/bin/env python3
"""Smoke-test the PISA stage on a few known X-ray protein-peptide complexes BEFORE a
full run. It exercises the real pipeline code path (run_pisa.run_one -> cross ->
classify), so a pass means CCP4 is wired correctly and the biological-vs-crystal
annotation is being produced.

What it checks, per test entry (PDB id + expected pep/prot chain):
  1. CCP4 `pisa` + pisa.cfg are present and callable
  2. PISA runs on the full mmCIF and returns parseable interfaces with CSS values
  3. the pep-pro crossing finds an interface for the {pep,prot} chain set
  4. classify() yields a sensible label (biological/crystal-packing) with a CSS

Usage (on the Linux box with CCP4 activated):
    python tests/smoke_test_pisa.py --ccp4-dir /opt/xtal/ccp4-9
    python tests/smoke_test_pisa.py --ccp4-dir ~/ccp4-9 --threshold 0.5 \
        --entry 1WRZ:B:A --entry 1A1M:C:A
CIFs are downloaded from RCSB into a temp dir (needs network); nothing else is written.
Exit code 0 = all pass, 1 = a check failed, 2 = CCP4/PISA not found.
"""
import argparse
import gzip
import os
import sys
import tempfile
import urllib.request

# make workflow/scripts importable (run_pisa + its `import checkpoint`)
HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, "..", "workflow", "scripts")
sys.path.insert(0, SCRIPTS)

import run_pisa  # noqa: E402

RCSB = "https://files.rcsb.org/download/{pid}.cif.gz"
# (pid, pep_chain, prot_chain) — small X-ray complexes; 1WRZ-B-A is the entry the
# crossing logic was originally validated against on the live site.
DEFAULT_ENTRIES = [("1WRZ", "B", "A"), ("1A1M", "C", "A")]


def fetch_cif(pid, cif_dir):
    """Download <pid>.cif.gz into the sharded layout run_pisa.shard_path expects."""
    dest = run_pisa.shard_path(cif_dir, pid)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if not os.path.exists(dest):
        with urllib.request.urlopen(RCSB.format(pid=pid), timeout=60) as r:
            data = r.read()
        with open(dest, "wb") as fh:
            fh.write(data)
    # sanity: it must really be gzip
    with gzip.open(dest, "rt") as fh:
        fh.read(64)
    return dest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccp4-dir", required=True,
                    help="CCP4 root (has bin/pisa and share/pisa/pisa.cfg)")
    ap.add_argument("--threshold", type=float, default=0.5,
                    help="CSS biological threshold (matches config pisa.css_biological_threshold)")
    ap.add_argument("--timeout", type=int, default=600)
    ap.add_argument("--entry", action="append", default=[],
                    help="PID:PEP:PROT (repeatable); overrides the defaults")
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

    entries = ([tuple(e.split(":")) for e in args.entry] if args.entry
               else DEFAULT_ENTRIES)
    cif_dir = tempfile.mkdtemp(prefix="pisa_smoke_cif_")
    print(f"PISA smoke test — ccp4={ccp4}  threshold={args.threshold}  "
          f"entries={len(entries)}\n" + "-" * 70)

    failures = 0
    for pid, pep, prot in entries:
        try:
            cif = fetch_cif(pid, cif_dir)
        except Exception as e:                              # noqa: BLE001
            print(f"FAIL {pid}: could not fetch CIF ({e})")
            failures += 1
            continue
        record, status = run_pisa.run_one(
            (pid, cif, ccp4, cfg, args.timeout))
        if status != "ok":
            print(f"FAIL {pid}: PISA status={status} (expected ok)")
            failures += 1
            continue
        n_iface = record.get("n_interfaces", 0)
        iface = run_pisa.cross(pep, prot, record)
        cls = run_pisa.classify(status, iface, args.threshold)
        if iface is None:
            print(f"FAIL {pid}-{pep}-{prot}: no interface matched chain set "
                  f"{{{pep},{prot}}} among {n_iface} interfaces")
            failures += 1
            continue
        css = iface.get("css", "")
        area = iface.get("area", "")
        ok_css = run_pisa._css_val(iface) is not None
        flag = "OK " if ok_css and cls in ("biological", "crystal-packing") else "WARN"
        if flag != "OK ":
            failures += 1
        print(f"{flag} {pid}-{pep}-{prot}: interfaces={n_iface}  css={css}  "
              f"area={area}  class={cls}")

    print("-" * 70)
    print(f"{'PASS' if failures == 0 else 'FAIL'}: "
          f"{len(entries) - failures}/{len(entries)} entries OK")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
