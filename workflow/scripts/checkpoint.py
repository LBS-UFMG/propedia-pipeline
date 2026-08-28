"""Per-entry checkpointing + bounded parallelism for the heavy per-entry stages,
so a long run resumes after interruption instead of recomputing from zero.

Contract
--------
A stage supplies a top-level ``worker(item) -> (record, status)`` and a list of
self-contained ``items``. Each item is processed independently; the PARENT process
persists one JSON file per id into a DURABLE work dir that snakemake does NOT track,
so partial progress survives a kill. On restart, items whose file exists are skipped.

``status`` convention: a status starting with ``"retry:"`` is a transient failure —
it is NOT persisted, so the item is retried on the next run. Any other status
(``"ok"``, ``"no_contacts"``, ``"missing"``, …) IS persisted (with its record, which
may be ``None`` for a legitimate no-result), so it is never recomputed.

The stage rebuilds its declared output from the returned records at the very end.
If killed before that merge, the declared output never appeared → snakemake reruns
the rule → this helper resumes from the checkpoint files.

Workers must be picklable (define them at module top level) and items must contain
only picklable values (str/num/tuple), so the pool works under both fork and spawn.
"""
import concurrent.futures as cf
import hashlib
import json
import os
import sys


def namespace(base_dir, version, params):
    """Return (and create) a checkpoint work dir namespaced by a signature of
    ``version`` + ``params``. A code change (bump ``version``) or a param change
    yields a fresh dir, so stale results are never silently reused."""
    h = hashlib.sha1()
    h.update(str(version).encode())
    h.update(repr(sorted((str(k), str(v)) for k, v in dict(params).items())).encode())
    workdir = os.path.join(base_dir, h.hexdigest()[:12])
    os.makedirs(workdir, exist_ok=True)
    return workdir


def _load(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except Exception:                                       # noqa: BLE001
        return None


def _atomic_write(path, obj):
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w") as fh:
        json.dump(obj, fh)
    os.replace(tmp, path)                                   # atomic within the dir


def run(items, worker, workdir, threads=1, id_of=None, stage="", log_every=200):
    """Process ``items`` with per-id checkpointing. Returns {id: {"record","status"}}
    for ALL items (those resumed from disk + those computed this run)."""
    os.makedirs(workdir, exist_ok=True)
    id_of = id_of or (lambda it: it[0])

    results, todo = {}, []
    for it in items:
        iid = str(id_of(it))
        cached = _load(os.path.join(workdir, f"{iid}.json"))
        if cached is not None:
            results[iid] = cached
        else:
            todo.append((iid, it))
    print(f"[{stage}] resume: {len(results)} cached, {len(todo)} to compute "
          f"(threads={threads})", file=sys.stderr)

    done = 0

    def handle(iid, record, status):
        nonlocal done
        if not str(status).startswith("retry:"):
            _atomic_write(os.path.join(workdir, f"{iid}.json"),
                          {"record": record, "status": status})
        results[iid] = {"record": record, "status": status}
        done += 1
        if done % log_every == 0:
            print(f"[{stage}] {done}/{len(todo)} computed", file=sys.stderr)

    if threads and threads > 1 and todo:
        with cf.ProcessPoolExecutor(max_workers=threads) as ex:
            futs = {ex.submit(worker, it): iid for iid, it in todo}
            for fut in cf.as_completed(futs):
                iid = futs[fut]
                try:
                    record, status = fut.result()
                except Exception as exc:                    # noqa: BLE001
                    record, status = None, f"retry:exc:{exc}"
                handle(iid, record, status)
    else:
        for iid, it in todo:
            try:
                record, status = worker(it)
            except Exception as exc:                        # noqa: BLE001
                record, status = None, f"retry:exc:{exc}"
            handle(iid, record, status)

    print(f"[{stage}] DONE: {len(results)} total ({done} this run)", file=sys.stderr)
    return results
