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
from concurrent.futures.process import BrokenProcessPool

# A native crash (segfault) in a worker kills its process and BREAKS the whole
# ProcessPoolExecutor — every still-pending future then raises BrokenProcessPool.
# We process in fresh-pool CHUNKS so a crash only affects its own chunk, and the
# unfinished items of a crashed chunk are re-run one-per-subprocess to isolate the
# culprit. `SEGFAULT_STATUS` is PERSISTED (never retried); a one-retry guard keeps a
# merely transient crash from permanently condemning an innocent item.
SEGFAULT_STATUS = "segfault"


def _one_isolated(worker, it):
    """Run ONE item in its own single-worker pool so a native crash is contained to
    it. Returns (record, status). A crash -> (None, SEGFAULT_STATUS); a normal Python
    exception -> transient 'retry:exc:...'."""
    ex = cf.ProcessPoolExecutor(max_workers=1)
    try:
        fut = ex.submit(worker, it)
        try:
            return fut.result()
        except BrokenProcessPool:
            return None, SEGFAULT_STATUS
        except Exception as exc:                               # noqa: BLE001
            return None, f"retry:exc:{exc}"
    finally:
        # wait=True DRAINS the pool's worker process + management thread before we
        # move on. With wait=False they linger undrained, and the interpreter's
        # concurrent.futures atexit handler then blocks forever joining them at exit
        # (a hang that only shows AFTER the stage's work is done). A broken pool's
        # workers are already dead, so this still returns promptly.
        ex.shutdown(wait=True, cancel_futures=True)


def _isolate(worker, pending):
    """Process `pending` [(iid, it), ...] one-per-subprocess, retrying a crash ONCE
    before condemning it as SEGFAULT_STATUS (guards against a transient native crash).
    Yields (iid, record, status)."""
    for iid, it in pending:
        rec, st = _one_isolated(worker, it)
        if st == SEGFAULT_STATUS:                              # confirm before persisting
            rec, st = _one_isolated(worker, it)
        yield iid, rec, st


def _run_chunk(worker, chunk, threads):
    """Run one chunk through a fresh pool. Returns (done: {iid: (rec, status)},
    crashed: bool). `done` holds only items that COMPLETED before any pool break."""
    done = {}
    ex = cf.ProcessPoolExecutor(max_workers=threads)
    try:
        futs = {ex.submit(worker, it): iid for iid, it in chunk}
        for fut in cf.as_completed(futs):
            iid = futs[fut]
            try:
                done[iid] = fut.result()
            except BrokenProcessPool:
                return done, True                              # pool is dead
            except Exception as exc:                           # noqa: BLE001
                done[iid] = (None, f"retry:exc:{exc}")
        return done, False
    except BrokenProcessPool:
        return done, True
    finally:
        # wait=True so no undrained worker/management thread survives to the
        # interpreter's atexit join (see _one_isolated). A crashed pool is already
        # dead here, so draining is immediate.
        ex.shutdown(wait=True, cancel_futures=True)


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

    crashed_total = 0
    if threads and threads > 1 and todo:
        # Process in fresh-pool chunks so a native crash is contained to its chunk;
        # isolate the crashed chunk's unfinished items to catch the culprit.
        chunk_size = max(threads * 32, 64)
        for i in range(0, len(todo), chunk_size):
            chunk = todo[i:i + chunk_size]
            got, crashed = _run_chunk(worker, chunk, threads)
            for iid, (record, status) in got.items():
                handle(iid, record, status)
            if crashed:
                leftover = [(iid, it) for iid, it in chunk if iid not in got]
                print(f"[{stage}] worker crash in chunk; isolating "
                      f"{len(leftover)} item(s) to find the culprit", file=sys.stderr)
                for iid, record, status in _isolate(worker, leftover):
                    if status == SEGFAULT_STATUS:
                        crashed_total += 1
                        print(f"[{stage}] SEGFAULT on {iid} -> persisted "
                              f"'{SEGFAULT_STATUS}' (skipped)", file=sys.stderr)
                    handle(iid, record, status)
    elif todo:
        # Single-threaded: still isolate each item, because an in-process native crash
        # would kill the whole run.
        for iid, record, status in _isolate(worker, todo):
            if status == SEGFAULT_STATUS:
                crashed_total += 1
            handle(iid, record, status)

    tail = f", {crashed_total} segfaulted" if crashed_total else ""
    print(f"[{stage}] DONE: {len(results)} total ({done} this run{tail})",
          file=sys.stderr)
    return results
