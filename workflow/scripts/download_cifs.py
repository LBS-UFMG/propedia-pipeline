"""Rule 2: download gzipped mmCIF for every candidate PDB ID.

Idempotent and resumable: files already present are skipped, downloads are
validated (full gzip decompress) in memory before an atomic rename into a
sharded directory, so any file on disk is known-good. HTTP 404 is treated as
'not found' (obsolete/withdrawn) and recorded rather than retried.
"""
import gzip
import os
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

BASE_URL = "https://files.rcsb.org/download/{pid}.cif.gz"
UA = "propedia-pipeline/1.0 (academic; contact: raquelcm@dcc.ufmg.br)"


def shard_path(cif_dir, pid):
    sub = pid[1:3].lower() if len(pid) >= 3 else "_"
    return os.path.join(cif_dir, sub, f"{pid}.cif.gz")


def download_one(pid, cif_dir, retries, timeout):
    final = shard_path(cif_dir, pid)
    if os.path.exists(final) and os.path.getsize(final) > 0:
        return pid, "skipped"
    os.makedirs(os.path.dirname(final), exist_ok=True)
    url = BASE_URL.format(pid=pid)
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
            gzip.decompress(data)                      # integrity check
            tmp = f"{final}.tmp.{os.getpid()}"
            with open(tmp, "wb") as fh:
                fh.write(data)
            os.replace(tmp, final)                     # atomic within shard dir
            return pid, "downloaded"
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return pid, "not_found"
            last = f"http_{exc.code}"
        except Exception as exc:                       # noqa: BLE001
            last = str(exc)
        time.sleep(2 * (attempt + 1))
    return pid, f"failed:{last}"


def main():
    ids_file = snakemake.input.ids                     # noqa: F821
    cif_dir = snakemake.params.cif_dir                 # noqa: F821
    retries = snakemake.params.max_retries             # noqa: F821
    timeout = snakemake.params.timeout_seconds         # noqa: F821
    workers = snakemake.threads                        # noqa: F821
    done_out = snakemake.output.done                   # noqa: F821
    fail_out = snakemake.output.failures               # noqa: F821

    with open(ids_file) as fh:
        ids = [ln.strip() for ln in fh if ln.strip()]

    counts = {"downloaded": 0, "skipped": 0, "not_found": 0, "failed": 0}
    failures = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(download_one, pid, cif_dir, retries, timeout): pid
                   for pid in ids}
        for fut in tqdm(as_completed(futures), total=len(ids),
                        desc="Downloading CIFs", unit="file"):
            pid, status = fut.result()
            key = "failed" if status.startswith("failed") else status
            counts[key] += 1
            if key in ("failed", "not_found"):
                failures.append(f"{pid}\t{status}")

    os.makedirs(os.path.dirname(fail_out), exist_ok=True)
    with open(fail_out, "w") as fh:
        fh.write("\n".join(failures) + ("\n" if failures else ""))
    with open(done_out, "w") as fh:
        fh.write(f"snapshot download complete\n{counts}\n")
    print(f"DONE {counts}", file=sys.stderr)
    # Only truly failed (not 404) should be considered an error worth rerunning.
    if counts["failed"]:
        print(f"WARNING: {counts['failed']} transient failures — "
              f"rerun to retry (see {fail_out})", file=sys.stderr)


if __name__ == "__main__":
    main()