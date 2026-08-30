"""Capture the exact software environment for a reproducible, self-describing release.

`requirements.txt` deliberately carries version RANGES (the dev spec). For deposition
(the NAR "deposit source code" requirement) a release must instead pin the EXACT
versions that produced it. This module records, at release time:
  - the Python interpreter version,
  - the installed versions of every Python dependency (a full pip freeze -> a
    `requirements-lock.txt` shipped inside the release), and
  - the versions / git commits of the external tools (COCaDA, SIGNA, iFeature,
    CCP4/PISA, prodigy-prot).

It is a plain importable module (no Snakemake dependency) so it can also be run
standalone for a quick environment audit:  `python workflow/scripts/env_capture.py`.
Every probe is best-effort: anything that cannot be determined is recorded as null
rather than failing the release.
"""
import os
import subprocess
import sys

# Python packages whose exact version is worth surfacing in the manifest summary
# (the full set is in the shipped requirements-lock.txt). Keep in sync with
# requirements.txt.
KEY_PACKAGES = [
    "biopython", "freesasa", "prodigy-prot", "scikit-learn",
    "snakemake", "numpy", "pandas", "tqdm",
]


def _run(cmd, cwd=None):
    """Run a command, return stripped stdout, or None on any failure."""
    try:
        out = subprocess.check_output(
            cmd, cwd=cwd, text=True, stderr=subprocess.DEVNULL, timeout=30)
        return out.strip() or None
    except Exception:                                          # noqa: BLE001
        return None


def python_version():
    return sys.version.split()[0]


def pip_freeze():
    """Full `pip freeze` as a single string (for requirements-lock.txt). None if
    pip is unavailable."""
    return _run([sys.executable, "-m", "pip", "freeze", "--all"])


def package_versions(names=KEY_PACKAGES):
    """Map package name -> installed version (None if not importable)."""
    try:
        from importlib.metadata import PackageNotFoundError, version
    except ImportError:                                        # pragma: no cover
        return {n: None for n in names}
    out = {}
    for n in names:
        try:
            out[n] = version(n)
        except PackageNotFoundError:
            out[n] = None
    return out


def git_ref(path):
    """Describe a git checkout: {commit, ref, dirty}. None if not a git repo."""
    path = os.path.expanduser(path or "")
    if not path or not os.path.isdir(path):
        return None
    commit = _run(["git", "rev-parse", "HEAD"], cwd=path)
    if commit is None:
        return None                       # present but not a git checkout
    # nearest tag, else branch/short-sha
    ref = (_run(["git", "describe", "--tags", "--always", "--dirty"], cwd=path)
           or _run(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=path))
    dirty = bool(_run(["git", "status", "--porcelain"], cwd=path))
    return {"commit": commit, "ref": ref, "dirty": dirty}


def ccp4_version(ccp4_dir):
    """Best-effort CCP4 version string. None if not resolvable."""
    ccp4_dir = os.path.expanduser(ccp4_dir or "")
    if not ccp4_dir or not os.path.isdir(ccp4_dir):
        return None
    # CCP4 ships a version stamp in a few known locations across releases.
    for rel in ("lib/ccp4/MASTER_version", "include/ccp4/ccp4_version.h"):
        p = os.path.join(ccp4_dir, rel)
        if os.path.isfile(p):
            try:
                with open(p) as fh:
                    txt = fh.read().strip()
                if txt:
                    return txt.splitlines()[0][:200]
            except OSError:
                pass
    # fall back to the directory basename (e.g. "ccp4-9")
    return os.path.basename(os.path.normpath(ccp4_dir)) or None


def tool_versions(machine):
    """External-tool provenance from the config `machine:` block (paths).

    machine: the config['machine'] dict (external tool locations).
    Returns {tool: {path, commit, ref, dirty} | {path, version} | None}.
    """
    machine = machine or {}
    tools = {}
    for key in ("cocada_dir", "signa_dir", "ifeature_dir"):
        path = machine.get(key)
        name = key[:-4]                                        # strip "_dir"
        ref = git_ref(path)
        tools[name] = ({"path": os.path.expanduser(path or ""), **ref}
                       if ref else {"path": path, "commit": None})
    # CCP4 (for PISA): not a git checkout — version stamp / dir name
    ccp4 = machine.get("ccp4_dir")
    tools["ccp4"] = {"path": ccp4, "version": ccp4_version(ccp4)}
    # prodigy-prot: installed as a Python package
    tools["prodigy-prot"] = {
        "version": package_versions(["prodigy-prot"])["prodigy-prot"]}
    return tools


def capture(machine=None):
    """Full environment snapshot as a JSON-serializable dict."""
    return {
        "python": python_version(),
        "packages": package_versions(),
        "tools": tool_versions(machine),
    }


def main():
    import json
    print(json.dumps(capture(), indent=2))


if __name__ == "__main__":
    main()
