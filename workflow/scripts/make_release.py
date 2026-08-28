"""Snapshot the finished deliverables into a dated, immutable release directory with
a provenance manifest, so each Propedia release is fully described (snapshot date,
git commit, tool config, row/entry counts) and reproducible.

The DURABLE state that produced it (per-entry checkpoints + extracted structures)
lives in `state/<mode>/` and is preserved between releases — an update bumps
`pdb_snapshot_date`, re-runs, and only the new PDB entries are computed."""
import csv
import json
import os
import shutil
import subprocess
import sys
import time


def _rows(path, delim):
    with open(path) as fh:
        return max(0, sum(1 for _ in fh) - 1)


def _git(*args):
    try:
        return subprocess.check_output(["git", *args], text=True,
                                       stderr=subprocess.DEVNULL).strip()
    except Exception:                                      # noqa: BLE001
        return None


def main():
    p = snakemake.params                                   # noqa: F821
    inp = snakemake.input                                  # noqa: F821
    os.makedirs(p.release_dir, exist_ok=True)

    # copy the deliverables verbatim into the release dir
    copied = {}
    for key, src in (("propedia.csv", inp.propedia),
                     ("multipro_final.csv", inp.multipro),
                     ("reproduction_report.txt", inp.report)):
        shutil.copy2(src, os.path.join(p.release_dir, key))
        copied[key] = os.path.getsize(src)

    manifest = {
        "release": os.path.basename(p.release_dir),
        "snapshot_date": p.snapshot_date,
        "mode": p.mode,
        "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "git_commit": _git("rev-parse", "HEAD"),
        "git_branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
        "git_dirty": bool(_git("status", "--porcelain")),
        "state_dir": p.state_dir,           # the durable memory that produced this
        "counts": {
            "propedia_entries": _rows(inp.propedia, ";"),
            "multipro_entries": _rows(inp.multipro, ";"),
        },
        "files": copied,
    }
    with open(snakemake.output.manifest, "w") as fh:       # noqa: F821
        json.dump(manifest, fh, indent=2)

    print(f"RELEASE {manifest['release']}: "
          f"{manifest['counts']['propedia_entries']} propedia + "
          f"{manifest['counts']['multipro_entries']} multipro entries "
          f"(commit {manifest['git_commit']}) -> {p.release_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
