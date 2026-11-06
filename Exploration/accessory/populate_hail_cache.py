#!/usr/bin/env python
"""Pre-populate the persistent gnomAD/ClinVar Hail cache (hail_cache.sqlite).

Purpose
    Real annotation jobs only query Hail for variants that are missing from the
    persistent cache.  This script fills that cache ahead of time for a list of
    variants, so later jobs are served from the cache.  It writes nothing except
    the cache itself (no user data); all progress is logged to stderr.

Accepted inputs (one per line; the first tab-separated column is used; blank
lines and lines starting with '#' are ignored; same rules as the web form)
    rsID                    e.g. rs429358
    chr:pos:ref:alt         e.g. 19:44908684:T:C
    chr:pos                 e.g. 19:44908684
    chr:start-end interval  <= 1 Mb; only variants with MAF > 0.01 are used

    rsID and chr:pos resolve to ONE allele, exactly as the annotation pipeline
    does (same functions), so the cache keys match those of real jobs.  Use
    --build to say which genome build the coordinates are in (hg19/hg38, or the
    aliases grch37/grch38).  Only variants found in gnomAD are cached; variants
    not in gnomAD are counted and queried again next time.  LD partners do not
    need to be cached.

How to run (inside the container, the input must be on a mounted path)
    The container is named "snpxplorer" (see info_for_running.sh) and the app
    uses `python`:

    docker exec -d snpxplorer sh -c 'python /Exploration/accessory/populate_hail_cache.py --input /Annotation/cache/inputs.txt --build hg38 >> /Annotation/LOGS/populate_hail_cache.log 2>&1'

    Use --dry-run first to see how many variants are resolved, already cached
    and planned for querying (no Hail, no writes).

Behaviour
    - Resumable: variants already in the cache are skipped, so just rerun.
    - Ctrl-C safe: the Hail slot is released; finished chunks stay cached.
    - Takes ONE Hail slot per chunk (shared with web jobs, HAIL_MAX_CONCURRENT).
    - Speed: roughly 1 hour per ~1000 new variants.  Larger sorted chunks
      amortise the partition reads (default --chunk-size 2000).
    - Safeguard: if the existing cache was built with different gnomAD/ClinVar
      versions, writing would drop and recreate the whole cache.  The script
      then aborts (exit 4) unless --allow-rebuild is given.

Exit codes: 0 ok, 1 usage/unreadable input, 2 nothing resolvable,
3 at least one chunk failed, 4 cache version mismatch, 130 interrupted.

Automatic mode: --add N (instead of --input)
    Caches the next N not-yet-cached variants, without preparing an input file,
    and remembers where it stopped, so every run continues where the last ended.

    docker exec snpxplorer python /Exploration/accessory/populate_hail_cache.py --add 1000 --dry-run
    docker exec -d snpxplorer sh -c 'python /Exploration/accessory/populate_hail_cache.py --add 1000 >> /Annotation/LOGS/populate_hail_cache.log 2>&1'

    Sources and order (--source auto, the default): first the GWAS catalog
    variants (--gwas-file, hg38, in file order), and when that list is
    exhausted the remainder is taken from the whole variant_info table in
    rowid order (--source gwas / --source genome select one of them).
    MAF filter: only variants with min(AF, 1-AF) >= --min-maf (default 0.01)
    are picked.  --build is ignored with --add (everything is hg38).
    Bookmark file: --state-file, default /Annotation/cache/populate_state.json
    (one entry per source and MAF threshold; written atomically after every
    successful chunk).  Variants that are not in gnomAD are never cached, but
    the bookmark moves past them, so they are not queried again.  A changed
    GWAS file (fingerprint) restarts the GWAS list from its beginning; a new
    --min-maf starts a new bookmark from the beginning (the old one is kept).
    Single-run lock: a second --add run on the same state file exits with 1
    (lock file <state-file>.lock); --dry-run takes no lock and writes nothing.
    All paths must be container paths (/Annotation/..., /Data/...), because the
    script runs inside the container.  Exit 2 = every requested source is
    exhausted and nothing was picked; exit 3 = a chunk failed (bookmark stays
    after the last good chunk; rerun to retry).

Parallel mode: --parallel N (1..6, default 1; ignored with --dry-run)
    Runs N chunks at the same time in N separate OS processes (spawned, not
    forked), each taking ONE Hail slot, exactly like a web job.  Every Hail
    process needs about 1 GB RAM, so N is limited to 6: at least 2 of the 8
    Hail slots (HAIL_MAX_CONCURRENT) stay free for user jobs.  If fewer slots
    are free, the extra chunks simply wait for a slot.  Recommended:

    docker exec -d snpxplorer sh -c 'python /Exploration/accessory/populate_hail_cache.py --add 12000 --parallel 6 >> /Annotation/LOGS/populate_hail_cache.log 2>&1'

    Chunks may finish out of order.  With --add the bookmark only ever moves to
    the end of the contiguous prefix of finished chunks (chunk 1..k all ok), so
    a crash or Ctrl-C never skips an unfinished chunk; chunks finished beyond a
    gap are cached but not bookmarked and are skipped on the rerun because they
    are cached.  With --add the first failed chunk stops new submissions, the
    running chunks are finished, and the exit code is 3.  Ctrl-C / SIGTERM
    terminate all workers and release their slots (exit 130).
    Every chunk runs in a FRESH process, and after each chunk (parallel and
    --parallel 1) leftover child processes named "java" (Hail's JVM, which
    hl.stop() does not end) are terminated, so JVMs never pile up.
"""
import argparse
import gzip
import hashlib
import importlib
import json
import multiprocessing as mpc
import os
import re
import signal
import sqlite3
import sys
import tempfile
import threading
import time
import types
import urllib.parse
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from pathlib import Path

try:  # POSIX only; the container is Linux
    import fcntl
except ImportError:  # pragma: no cover
    fcntl = None

DEFAULT_CHUNK = 2000
DEFAULT_STATE_FILE = "/Annotation/cache/populate_state.json"
DEFAULT_GWAS_RELATIVE = "databases/GWAS_catalog/Gwas_catalog_hg38.txt.gz"
DEFAULT_BATCH = 50000
DEFAULT_MIN_MAF = 0.01
CACHE_CHECK_BATCH = 1000
PROGRESS_EVERY = 1000000
MAF_EPSILON = 1e-12
STATE_VERSION = 1
MAX_PARALLEL = 6
DEPS_ENV = "POPULATE_HAIL_DEPS_MODULE"  # test hook: module with load() -> deps for the workers
VALID_CHROMS = {str(c) for c in range(1, 23)} | {"X", "Y", "MT"}
EXIT_OK, EXIT_USAGE, EXIT_NOTHING, EXIT_CHUNK_FAILED, EXIT_VERSION, EXIT_INTERRUPT = 0, 1, 2, 3, 4, 130


def _stderr(message):
    print(message, file=sys.stderr, flush=True)


def read_inputs(path):
    """Return the unique, non-comment first-column entries of path, in order."""
    values = []
    with open(path, encoding="utf-8") as source:
        for line in source:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            first = line.split("\t")[0].strip()
            if first:
                values.append(first)
    return list(dict.fromkeys(values))


def normalize_build(build):
    """Same normalisation as the pipeline's main(): grch38 -> hg38, grch37 -> hg19."""
    value = str(build).strip().lower()
    value = {"grch38": "hg38", "grch37": "hg19"}.get(value, value)
    if value not in ("hg19", "hg38"):
        raise ValueError("invalid build %r (use hg19, hg38, grch37 or grch38)" % (build,))
    return value


def variant_id(key):
    return "%s-%s-%s-%s" % (key[0], key[1], key[2], key[3])


def resolve_all(lines, build, deps):
    """Resolve inputs to hg38 variant keys.

    Returns (sorted unique keys, [(input, reason)]).  Mirrors
    prepare_annotation_queries: explicit inputs resolve to one variant;
    intervals expand to all variants with MAF > INTERVAL_MAF_CUTOFF.
    """
    keys = set()
    unresolved = []
    for line in lines:
        entry = deps.classify_submission_line(line)
        kind = entry["type"]
        if kind == "invalid":
            unresolved.append((line, entry.get("reason") or "invalid_input_format"))
        elif kind == "interval":
            try:
                variants = deps.fetch_interval_variants(entry, build, deps.DB_FILE)
            except Exception as error:  # ValueError for failed liftover, DB errors
                _stderr("[POPULATE] interval %s failed: %s" % (line, error))
                unresolved.append((line, "interval_error"))
                continue
            records = variants.to_dict(orient="records") if hasattr(variants, "to_dict") else list(variants)
            for record in records:
                try:
                    maf = float(record.get("AF"))
                except (TypeError, ValueError):
                    continue
                if maf != maf or maf <= deps.INTERVAL_MAF_CUTOFF:
                    continue
                keys.add(tuple(deps.variant_key(record)))
        else:
            resolved = deps.resolve_explicit_variant(entry, build, deps.DB_FILE)
            if resolved is None:
                reason = "not_in_db"
                if build == "hg19" and kind in ("coord", "coord_alleles"):
                    lifted = deps.liftover_hg19_to_hg38("chr%s" % entry["chrom"], entry["pos"])
                    if lifted is None:
                        reason = "liftover_failed"
                unresolved.append((line, reason))
            else:
                keys.add(tuple(deps.variant_key(resolved)))
    return sorted(keys), unresolved


def plan_chunks(keys, cache_read, size):
    """Drop keys already cached; return sorted chunks of at most size keys."""
    ordered = sorted({tuple(key) for key in keys})
    cached = cache_read([variant_id(key) for key in ordered]) if ordered else {}
    misses = [key for key in ordered if variant_id(key) not in cached]
    return [misses[start:start + size] for start in range(0, len(misses), size)]


def _verify_persisted(deps, remaining, fresh):
    """Raise RuntimeError unless every queried variant was returned and, when found, re-reads from the cache."""
    ids = [variant_id(key) for key in remaining]
    if any(i not in fresh for i in ids):
        raise RuntimeError("cache write not persisted")
    found = [i for i in ids if fresh[i].get("total_af") is not None]
    if found:
        back = deps._hail_cache_read(found)
        if any(i not in back for i in found):
            raise RuntimeError("cache write not persisted")


def _proc_info(pid):
    """Return (comm, ppid, state) of a Linux process from /proc, or None."""
    try:
        with open("/proc/%d/stat" % pid) as handle:
            data = handle.read()
        comm = data[data.index("(") + 1:data.rindex(")")]
        fields = data.rsplit(")", 1)[1].split()
        return comm, int(fields[1]), fields[0]
    except (OSError, ValueError, IndexError):
        return None


def _child_pids(parent=None, comm=None):
    """PIDs of the (non-zombie) direct children of parent (default: this process), optionally by comm."""
    parent = os.getpid() if parent is None else parent
    found = []
    try:
        names = os.listdir("/proc")
    except OSError:
        return found
    for name in names:
        if not name.isdigit():
            continue
        info = _proc_info(int(name))
        if info is not None and info[1] == parent and info[2] != "Z" and (comm is None or info[0] == comm):
            found.append(int(name))
    return found


def _reap(pid):
    try:
        os.waitpid(pid, os.WNOHANG)
    except OSError:
        pass


def _kill_children(comm="java", parent=None, grace=10.0):
    """SIGTERM the direct children of this process whose comm is `comm` (None = all), wait up to
    `grace` s, then SIGKILL the survivors.  Hail's hl.stop() leaves its Java gateway running.
    Never raises; returns the number of children found.  Linux only (no /proc: returns 0)."""
    try:
        pids = _child_pids(parent, comm)
        for pid in pids:
            try:
                os.kill(pid, signal.SIGTERM)
            except OSError:
                pass
        deadline = time.time() + grace
        pending = list(pids)
        while pending:
            for pid in pending:
                _reap(pid)
            pending = [pid for pid in pending if (_proc_info(pid) or (None, None, "Z"))[2] != "Z"]
            if not pending or time.time() >= deadline:
                break
            time.sleep(0.05)
        for pid in pending:
            try:
                os.kill(pid, signal.SIGKILL)
            except OSError:
                pass
        end = time.time() + 5
        for pid in pending:
            while time.time() < end:  # SIGKILL cannot be ignored; give the kernel a moment
                _reap(pid)
                if (_proc_info(pid) or (None, None, "Z"))[2] == "Z":
                    break
                time.sleep(0.02)
        return len(pids)
    except Exception:
        return 0


def _cleanup_java(comm="java"):
    """Kill leftover Hail JVMs among our children and log how many there were."""
    count = _kill_children(comm)
    if count:
        _stderr("[HAIL] cleaned up %d leftover java process(es)" % count)
    return count


def populate(chunks, deps, log, clock, on_success=None, stop_on_failure=False, verify_write=False):
    """Query and cache each chunk under one Hail slot; return statistics.

    on_success(index) is called (1-based chunk index) after each successful
    chunk; stop_on_failure stops at the first failed chunk; verify_write turns
    a silently failed cache write into a chunk failure.
    """
    stats = {"chunks": len(chunks), "done": 0, "failed_chunks": 0,
             "queried": 0, "cached": 0, "not_in_gnomad": 0}
    for index, chunk in enumerate(chunks, 1):
        started = clock()
        queried = cached = missing = 0
        gate = deps.acquire_slot(deps.HAIL_SLOT_PREFIX, deps.HAIL_MAX_CONCURRENT,
                                 poll_seconds=10, label="HAIL", what="Hail query slot")
        try:
            try:
                # Another job may have filled the cache while we waited.
                present = deps._hail_cache_read([variant_id(key) for key in chunk])
                remaining = [key for key in chunk if variant_id(key) not in present]
                if remaining:
                    fresh = deps._query_hail(remaining)
                    queried = len(remaining)
                    deps._hail_cache_write(fresh.values())
                    if verify_write:
                        _verify_persisted(deps, remaining, fresh)
                    for annotation in fresh.values():
                        if annotation.get("total_af") is None:
                            missing += 1
                        else:
                            cached += 1
                stats["done"] += 1
            except Exception as error:  # includes GnomadLookupError; KeyboardInterrupt propagates
                stats["failed_chunks"] += 1
                log("chunk %d/%d FAILED (%s): %s" % (index, len(chunks), type(error).__name__, error))
                if stop_on_failure:
                    break
                continue
        finally:
            gate.close()
            _cleanup_java()  # hl.stop() leaves the Java gateway running; only "java" children are touched
        stats["queried"] += queried
        stats["cached"] += cached
        stats["not_in_gnomad"] += missing
        log("chunk %d/%d: queried %d, cached %d, not in gnomAD %d, elapsed %.0fs"
            % (index, len(chunks), queried, cached, missing, clock() - started))
        if on_success is not None:
            on_success(index)
    return stats


# ---------------------------------------------------------------- parallel mode (--parallel)
_WORKER_DEPS = None


def _worker_sigterm(signum, frame):
    raise SystemExit(143)


def _worker_init(deps_module):
    """Runs once in every spawned worker: signals, then load the dependencies."""
    global _WORKER_DEPS
    signal.signal(signal.SIGINT, signal.SIG_IGN)  # Ctrl-C is handled by the parent
    signal.signal(signal.SIGTERM, _worker_sigterm)  # finally blocks release the slot
    if deps_module:
        _WORKER_DEPS = importlib.import_module(deps_module).load()
    else:
        _WORKER_DEPS = load_pipeline()


def _chunk_worker(index, keys, verify_write):
    """One chunk under one Hail slot (mirrors populate()); never raises Exception."""
    deps = _WORKER_DEPS
    result = {"index": index, "queried": 0, "cached": 0, "missing": 0,
              "t_start": time.time(), "t_end": None, "pid": os.getpid(), "error": None}
    try:
        gate = None
        try:
            gate = deps.acquire_slot(deps.HAIL_SLOT_PREFIX, deps.HAIL_MAX_CONCURRENT,
                                     poll_seconds=10, label="HAIL", what="Hail query slot")
            try:
                present = deps._hail_cache_read([variant_id(key) for key in keys])
                remaining = [key for key in keys if variant_id(key) not in present]
                if remaining:
                    fresh = deps._query_hail(remaining)
                    result["queried"] = len(remaining)
                    deps._hail_cache_write(fresh.values())
                    if verify_write:
                        _verify_persisted(deps, remaining, fresh)
                    for annotation in fresh.values():
                        if annotation.get("total_af") is None:
                            result["missing"] += 1
                        else:
                            result["cached"] += 1
            except Exception as error:  # Hail exceptions may not pickle: report text only
                result["error"] = "%s: %s" % (type(error).__name__, error)
                result["queried"] = result["cached"] = result["missing"] = 0
        finally:
            if gate is not None:
                gate.close()
            _cleanup_java()  # runs after the slot is released; also on SIGTERM (SystemExit)
    except SystemExit:
        # SIGTERM: the slot is released; a pool worker would otherwise survive the exception
        os._exit(143)
    result["t_end"] = time.time()
    return result


def _terminate_pool(executor, futures):
    """Kill every worker of the executor (terminate, wait 10 s, then kill); never raises."""
    for future in futures:
        future.cancel()
    processes = list((getattr(executor, "_processes", None) or {}).values())
    manager = (getattr(executor, "_executor_manager_thread", None)
               or getattr(executor, "_queue_management_thread", None))
    # Flag the shutdown first: while the manager thread is alive its wakeup pipe is still valid;
    # once the workers are dead it closes that pipe itself.
    try:
        executor.shutdown(wait=False)
    except Exception:
        pass
    for process in processes:
        try:
            process.terminate()
        except Exception:
            pass
    deadline = time.time() + 10
    for process in processes:
        try:
            process.join(max(0.0, deadline - time.time()))
        except Exception:
            pass
    for process in processes:
        try:
            if process.is_alive():
                for child in _child_pids(process.pid, "java"):  # worker never got to clean up
                    try:
                        os.kill(child, signal.SIGKILL)
                    except OSError:
                        pass
                os.kill(process.pid, signal.SIGKILL)
                process.join(5)
        except Exception:
            pass
    # let the executor's manager thread end before the interpreter exits (bounded wait)
    if manager is not None:
        try:
            manager.join(10)
        except Exception:
            pass


def _raise_interrupt(signum, frame):
    raise KeyboardInterrupt


def populate_parallel(chunks, n, deps_module, log, clock, on_success=None, stop_on_failure=False,
                      verify_write=False):
    """Like populate(), with up to n chunks in flight in n separate worker processes.

    Chunks complete in any order; on_success(index) is still called in index
    order and only for the contiguous prefix of successful chunks.
    """
    stats = {"chunks": len(chunks), "done": 0, "failed_chunks": 0,
             "queried": 0, "cached": 0, "not_in_gnomad": 0}
    if not chunks:
        return stats
    workers = max(1, min(n, len(chunks)))
    started = clock()
    context = mpc.get_context("spawn")

    def new_executor():
        # One single-worker pool per chunk = one fresh process per chunk, so no Hail/JVM state
        # survives from an earlier chunk (ProcessPoolExecutor has no max_tasks_per_child before 3.11).
        return ProcessPoolExecutor(max_workers=1, mp_context=context,
                                   initializer=_worker_init, initargs=(deps_module,))

    running = {}  # future -> chunk index
    pools = {}  # future -> its single-worker executor
    done = {}  # index -> succeeded?
    next_submit, next_commit = 1, 1
    stopping = False
    try:
        old_term = signal.signal(signal.SIGTERM, _raise_interrupt)
    except ValueError:  # not the main thread
        old_term = None
    try:
        try:
            while True:
                while not stopping and next_submit <= len(chunks) and len(running) < workers:
                    keys = [tuple(key) for key in chunks[next_submit - 1]]
                    pool = new_executor()
                    future = pool.submit(_chunk_worker, next_submit, keys, verify_write)
                    running[future] = next_submit
                    pools[future] = pool
                    next_submit += 1
                if not running:
                    break
                finished, _ = wait(list(running), return_when=FIRST_COMPLETED)
                for future in sorted(finished, key=lambda f: running[f]):
                    index = running.pop(future)
                    pool = pools.pop(future)
                    try:
                        result = future.result()
                    except Exception as error:  # e.g. BrokenProcessPool
                        result = {"index": index, "error": "%s: %s" % (type(error).__name__, error)}
                    pool.shutdown(wait=True)  # the worker process exits here
                    if result["error"] is not None:
                        stats["failed_chunks"] += 1
                        done[index] = False
                        kind, _, message = result["error"].partition(": ")
                        log("chunk %d/%d FAILED (%s): %s" % (index, len(chunks), kind, message))
                        if stop_on_failure:
                            stopping = True
                        continue
                    done[index] = True
                    stats["done"] += 1
                    stats["queried"] += result["queried"]
                    stats["cached"] += result["cached"]
                    stats["not_in_gnomad"] += result["missing"]
                    log("chunk %d/%d: queried %d, cached %d, not in gnomAD %d, elapsed %.0fs"
                        % (index, len(chunks), result["queried"], result["cached"], result["missing"],
                           result["t_end"] - result["t_start"]))
                while next_commit in done and done[next_commit]:
                    if on_success is not None:
                        on_success(next_commit)
                    next_commit += 1
        except KeyboardInterrupt:
            killers = [threading.Thread(target=_terminate_pool, args=(pool, [future]))
                       for future, pool in list(pools.items())]  # in parallel: each may wait for its JVM
            for killer in killers:
                killer.start()
            for killer in killers:
                killer.join()
            raise
    finally:
        if old_term is not None:
            signal.signal(signal.SIGTERM, old_term)
    elapsed = clock() - started
    log("parallel %d: elapsed %.0fs, %.2f variants/s"
        % (workers, elapsed, stats["queried"] / elapsed if elapsed > 0 else 0.0))
    return stats


# ---------------------------------------------------------------- automatic mode (--add)
class StateError(Exception):
    """Unreadable/corrupt state file or lock problem."""


def state_key(source, maf):
    return "%s|maf=%.6g" % (source, maf)


def load_state(path):
    """Read the bookmark file; absent -> empty state; corrupt -> StateError."""
    try:
        with open(path, encoding="utf-8") as handle:
            data = json.load(handle)
    except FileNotFoundError:
        return {"version": STATE_VERSION, "entries": {}}
    except (OSError, ValueError) as error:  # ValueError covers JSON and Unicode decode errors
        raise StateError("state file %s unreadable (%s): fix or delete it" % (path, error))
    if (not isinstance(data, dict) or data.get("version") != STATE_VERSION
            or not isinstance(data.get("entries"), dict)
            or not all(isinstance(v, dict) for v in data["entries"].values())):
        raise StateError("state file %s has an unexpected structure: fix or delete it" % path)
    return data


def save_state(path, state):
    """Atomic save: temp file in the same directory, fsync, os.replace."""
    directory = os.path.dirname(os.path.abspath(path))
    descriptor, temp = tempfile.mkstemp(dir=directory, prefix=".populate_state.", suffix=".tmp")
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(state, handle, indent=1, sort_keys=True)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp, path)
    except BaseException:
        try:
            os.unlink(temp)
        except OSError:
            pass
        raise


def acquire_lock(state_path):
    """Take a non-blocking exclusive lock on <state_path>.lock; return the open handle."""
    lock_path = str(state_path) + ".lock"
    try:
        os.makedirs(os.path.dirname(os.path.abspath(lock_path)), exist_ok=True)
        handle = open(lock_path, "a")
    except OSError as error:
        raise StateError("cannot create lock %s: %s" % (lock_path, error))
    try:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except (OSError, AttributeError):
        handle.close()
        raise StateError("another --add run holds the lock %s" % lock_path)
    return handle


def load_gwas_entries(path):
    """Ordered unique entries of the hg38 GWAS catalogue: all rsIDs of the SNP column
    (lowercase), or 'Chr:Pos(hg38)' for rows without an rsID."""
    entries = []
    with gzip.open(str(path), "rt", encoding="utf-8", errors="replace") as source:
        for line in source:
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) < 4:
                continue
            try:
                pos = int(cols[2].strip())
            except ValueError:
                continue
            found = [value.lower() for value in re.findall(r"rs\d+", cols[3], flags=re.IGNORECASE)]
            if found:
                entries.extend(found)
            else:
                chrom = cols[0].strip()
                if chrom.lower().startswith("chr"):
                    chrom = chrom[3:]
                entries.append("%s:%d" % (chrom.upper(), pos))
    return list(dict.fromkeys(entries))


def _chrom_ok(chrom):
    return str(chrom).upper() in VALID_CHROMS


def _new_counters():
    return {"scanned": 0, "filtered": 0, "cached": 0, "unresolved": 0,
            "duplicate": 0, "picked": 0, "exhausted": False}


def pick_gwas(entries, start, n, maf, deps, seen, log=None):
    """Pick up to n uncached, MAF-passing variants from entries[start:].

    Returns (picked [(key, entry_index)], scan_end, counters); scan_end is the
    exclusive index of the last examined entry.
    """
    counters = _new_counters()
    picked = []
    buffer = []  # (index, outcome, key)

    def flush():
        """Process buffered entries in order; return the scan_end if n was reached."""
        candidates = [item for item in buffer if item[1] == "cand"]
        ids = [variant_id(item[2]) for item in candidates]
        cached = deps._hail_cache_read(ids) if ids else {}
        for index, outcome, key in buffer:
            counters["scanned"] += 1
            if outcome != "cand":
                counters[outcome] += 1
                continue
            if key in seen:
                counters["duplicate"] += 1
                continue
            seen.add(key)
            if variant_id(key) in cached:
                counters["cached"] += 1
                continue
            picked.append((key, index))
            counters["picked"] += 1
            if len(picked) >= n:
                return index + 1
        return None

    index = start
    total = len(entries)
    while index < total:
        entry = entries[index]
        outcome, key = "cand", None
        parsed = deps.classify_submission_line(entry)
        resolved = None
        if parsed.get("type") in ("rsid", "coord", "coord_alleles"):
            resolved = deps.resolve_explicit_variant(parsed, "hg38", deps.DB_FILE)
        if resolved is None:
            outcome = "unresolved"
        else:
            try:
                af = float(resolved.get("AF"))
            except (TypeError, ValueError):
                af = None
            if af is None or af != af or min(af, 1.0 - af) < maf - MAF_EPSILON:
                outcome = "filtered"
            else:
                try:
                    key = tuple(deps.variant_key(resolved))
                except (TypeError, ValueError, KeyError):
                    outcome = "filtered"
                else:
                    if not _chrom_ok(key[0]):
                        outcome = "filtered"
        buffer.append((index, outcome, key))
        index += 1
        if sum(1 for item in buffer if item[1] == "cand") >= CACHE_CHECK_BATCH or index == total:
            end = flush()
            buffer = []
            if end is not None:
                return picked, end, counters
            if log and (index - start) % 100000 < CACHE_CHECK_BATCH * 2 and index != total:
                log("[POPULATE] gwas scanned %d/%d entries, picked %d" % (index, total, len(picked)))
    counters["exhausted"] = True
    return picked, max(index, start), counters


def pick_genome(start_rowid, n, maf, deps, batch, seen, log=None):
    """Pick up to n uncached, MAF-passing variants from variant_info with rowid > start_rowid.

    Returns (picked [(key, rowid)], scan_end rowid, counters).
    """
    counters = _new_counters()
    picked = []
    last = start_rowid
    next_log = (start_rowid // PROGRESS_EVERY + 1) * PROGRESS_EVERY
    uri = "file:" + urllib.parse.quote(str(deps.DB_FILE)) + "?mode=ro"
    connection = sqlite3.connect(uri, uri=True, timeout=30)
    try:
        sql = ("SELECT rowid, Chromosome, Position_hg38, REF, ALT, AF FROM variant_info "
               "WHERE rowid > ? AND AF >= ? AND AF <= ? ORDER BY rowid LIMIT ?")
        low, high = maf - MAF_EPSILON, 1.0 - maf + MAF_EPSILON
        while len(picked) < n:
            rows = connection.execute(sql, (last, low, high, batch)).fetchall()
            if not rows:
                counters["exhausted"] = True
                break
            items = []
            for rowid, chrom, pos, ref, alt, af in rows:
                record = {"Chromosome": chrom, "Position_hg38": pos, "REF": ref, "ALT": alt, "AF": af}
                try:
                    key = tuple(deps.variant_key(record))
                except (TypeError, ValueError, KeyError):
                    key = None
                items.append((rowid, key if key is not None and _chrom_ok(key[0]) else None))
            ids = [variant_id(key) for _, key in items if key is not None]
            cached = deps._hail_cache_read(ids) if ids else {}
            for rowid, key in items:
                last = rowid
                counters["scanned"] += 1
                if key is None:
                    counters["filtered"] += 1
                elif key in seen:
                    counters["duplicate"] += 1
                else:
                    seen.add(key)
                    if variant_id(key) in cached:
                        counters["cached"] += 1
                    else:
                        picked.append((key, rowid))
                        counters["picked"] += 1
                        if len(picked) >= n:
                            break
            if log and last >= next_log:
                log("[POPULATE] genome scan at rowid %d, picked %d" % (last, len(picked)))
                next_log = (last // PROGRESS_EVERY + 1) * PROGRESS_EVERY
    finally:
        connection.close()
    return picked, last, counters


def read_cache_meta(path):
    """Read hail_cache_meta read-only; return a dict, or None if absent/unreadable. Never raises."""
    connection = None
    try:
        uri = "file:" + urllib.parse.quote(str(path)) + "?mode=ro"
        connection = sqlite3.connect(uri, uri=True, timeout=30)
        tables = {name for (name,) in connection.execute(
            "SELECT name FROM sqlite_master WHERE type='table'")}
        if "hail_cache_meta" not in tables:
            return None
        return dict(connection.execute("SELECT key, value FROM hail_cache_meta"))
    except Exception:
        return None
    finally:
        if connection is not None:
            try:
                connection.close()
            except Exception:
                pass


def load_pipeline():
    """Import the annotation pipeline lazily (needs the container environment)."""
    sys.path.insert(0, "/Annotation/BIN")
    import standalone_annotation_optimized_mp as mp
    names = ["classify_submission_line", "resolve_explicit_variant", "variant_key",
             "fetch_interval_variants", "liftover_hg19_to_hg38", "INTERVAL_MAF_CUTOFF",
             "DB_FILE", "_hail_cache_read", "_hail_cache_write", "_hail_cache_versions",
             "HAIL_CACHE_PATH", "_query_hail", "GnomadLookupError", "acquire_slot",
             "HAIL_SLOT_PREFIX", "HAIL_MAX_CONCURRENT", "DATA_PATH"]
    return types.SimpleNamespace(**{name: getattr(mp, name) for name in names})


class _Parser(argparse.ArgumentParser):
    def error(self, message):
        self.print_usage(sys.stderr)
        _stderr("%s: error: %s" % (self.prog, message))
        raise SystemExit(EXIT_USAGE)


def _build_type(value):
    try:
        return normalize_build(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error))


def _chunk_type(value):
    try:
        number = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be an integer >= 1")
    if number < 1:
        raise argparse.ArgumentTypeError("must be an integer >= 1")
    return number


def _parallel_type(value):
    message = "must be an integer between 1 and %d" % MAX_PARALLEL
    try:
        number = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(message)
    if not 1 <= number <= MAX_PARALLEL:
        raise argparse.ArgumentTypeError(message)
    return number


def _maf_type(value):
    try:
        number = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be a number between 0 and 0.5")
    if not 0.0 <= number <= 0.5:
        raise argparse.ArgumentTypeError("must be a number between 0 and 0.5")
    return number


def _parse(argv):
    parser = _Parser(description="Pre-populate the gnomAD/ClinVar Hail cache.")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--input", help="file with one variant/interval per line")
    mode.add_argument("--add", type=_chunk_type, metavar="N",
                      help="automatically cache the next N uncached variants (GWAS catalogue, then whole genome)")
    parser.add_argument("--build", default=None, type=_build_type, metavar="{hg19,hg38,grch37,grch38}",
                        help="genome build of --input (default hg38; ignored with --add)")
    parser.add_argument("--chunk-size", type=_chunk_type, default=DEFAULT_CHUNK)
    parser.add_argument("--parallel", type=_parallel_type, default=1, metavar="N",
                        help="run N chunks at once in N processes, one Hail slot each (1-%d, default 1; "
                             "~1 GB RAM per process; ignored with --dry-run)" % MAX_PARALLEL)
    parser.add_argument("--dry-run", action="store_true",
                        help="resolve and report only; no Hail, no slot, no writes")
    parser.add_argument("--allow-rebuild", action="store_true",
                        help="allow dropping/recreating a cache with different versions")
    parser.add_argument("--source", choices=("auto", "gwas", "genome"), default="auto",
                        help="--add only: GWAS catalogue, whole variant_info, or GWAS then genome")
    parser.add_argument("--min-maf", type=_maf_type, default=DEFAULT_MIN_MAF,
                        help="--add only: minimum minor allele frequency (0-0.5, default 0.01)")
    parser.add_argument("--state-file", default=DEFAULT_STATE_FILE, help="--add only: bookmark file")
    parser.add_argument("--gwas-file", default=None,
                        help="--add only: hg38 GWAS catalogue (default <DATA_PATH>/%s)" % DEFAULT_GWAS_RELATIVE)
    parser.add_argument("--batch", type=_chunk_type, default=DEFAULT_BATCH,
                        help="--add only: rows per variant_info read")
    return parser.parse_args(argv)


def check_cache_versions(deps, dry_run, allow_rebuild):
    """Version safeguard shared by both modes; return EXIT_VERSION to abort, else None."""
    cache_path = deps.HAIL_CACHE_PATH
    if not os.path.exists(str(cache_path)):
        return None
    stored = read_cache_meta(cache_path)
    try:
        current = deps._hail_cache_versions()
    except Exception as error:
        _stderr("[POPULATE] cannot determine cache versions: %s" % error)
        current = None
    if stored is None or stored != current:
        message = ("existing cache meta %s differs from current versions %s; "
                   "writing would drop and recreate the whole cache" % (stored, current))
        if dry_run:
            print("WARNING: " + message + " (a real run would abort without --allow-rebuild)")
        elif not allow_rebuild:
            _stderr("[POPULATE] ABORT: %s. Use --allow-rebuild to proceed." % message)
            return EXIT_VERSION
        else:
            _stderr("[POPULATE] WARNING: %s; --allow-rebuild given" % message)
    return None


def _fingerprint(entries):
    return hashlib.sha256("\n".join(entries).encode("utf-8")).hexdigest()


def _bookmark_line(state, maf):
    genome = state["entries"].get(state_key("genome", maf), {})
    gwas = state["entries"].get(state_key("gwas", maf), {})
    return "bookmark now at genome rowid %s / gwas %s/%s (maf=%.6g)" % (
        genome.get("rowid", 0), gwas.get("index", 0), gwas.get("total", "?"), maf)


def _print_counters(source, c):
    print("%s: scanned %d, filtered %d, already cached %d, unresolved %d, duplicates %d, picked %d%s" % (
        source, c["scanned"], c["filtered"], c["cached"], c["unresolved"], c["duplicate"], c["picked"],
        " (source exhausted)" if c["exhausted"] else ""))


def run_add(args, deps, started):
    maf = args.min_maf
    if args.build is not None:
        _stderr("[POPULATE] --build %s ignored with --add (GWAS file and variant_info are hg38)" % args.build)
    sources = ["gwas", "genome"] if args.source == "auto" else [args.source]
    lock = None
    try:
        if not args.dry_run:
            try:
                lock = acquire_lock(args.state_file)
            except StateError as error:
                _stderr("[POPULATE] %s" % error)
                return EXIT_USAGE
        try:
            state = load_state(args.state_file)
        except StateError as error:
            _stderr("[POPULATE] %s" % error)
            return EXIT_USAGE

        entries = None
        if "gwas" in sources:
            gwas_path = args.gwas_file or (Path(str(deps.DATA_PATH)) / DEFAULT_GWAS_RELATIVE)
            try:
                entries = load_gwas_entries(gwas_path)
            except (OSError, EOFError, ValueError) as error:
                _stderr("[POPULATE] cannot read GWAS file %s: %s" % (gwas_path, error))
                if args.source == "gwas":
                    return EXIT_USAGE
                _stderr("[POPULATE] WARNING: continuing with the genome source only")
                sources.remove("gwas")

        existing = state["entries"]

        def ensure_entry(source):
            key = state_key(source, maf)
            if key not in existing:
                if existing:
                    _stderr("[POPULATE] WARNING: no bookmark for %s; starting from the beginning "
                            "(existing bookmarks: %s)" % (key, ", ".join(sorted(existing))))
                existing[key] = {"rowid": 0} if source == "genome" else {"index": 0}
            return existing[key]

        gwas_entry = genome_entry = None
        if "gwas" in sources:
            gwas_entry = ensure_entry("gwas")
            fingerprint = _fingerprint(entries)
            if "fingerprint" in gwas_entry and gwas_entry["fingerprint"] != fingerprint:
                _stderr("[POPULATE] WARNING: GWAS file fingerprint changed; restarting the GWAS list at index 0")
                gwas_entry["index"] = 0
            gwas_entry["fingerprint"] = fingerprint
            gwas_entry["total"] = len(entries)
            gwas_entry["index"] = int(gwas_entry.get("index", 0))

        # ---- pick
        seen = set()
        remaining = args.add
        plan = []  # (source, picked, scan_end, start)
        try:
            if gwas_entry is not None:
                start = gwas_entry["index"]
                picked, end, counters = pick_gwas(entries, start, remaining, maf, deps, seen, _stderr)
                plan.append(("gwas", picked, end, start))
                _print_counters("gwas", counters)
                remaining -= len(picked)
            if "genome" in sources and remaining > 0:
                genome_entry = ensure_entry("genome")
                genome_entry["rowid"] = int(genome_entry.get("rowid", 0))
                start = genome_entry["rowid"]
                picked, end, counters = pick_genome(start, remaining, maf, deps, args.batch, seen, _stderr)
                plan.append(("genome", picked, end, start))
                _print_counters("genome", counters)
        except KeyboardInterrupt:
            _stderr("[POPULATE] interrupted while scanning; bookmark unchanged since the last good chunk")
            return EXIT_INTERRUPT
        sys.stdout.flush()

        chunks, tags = [], []
        for source, picked, end, start in plan:
            for first in range(0, len(picked), args.chunk_size):
                part = picked[first:first + args.chunk_size]
                last_chunk = first + args.chunk_size >= len(picked)
                cover = end if last_chunk else (part[-1][1] + 1 if source == "gwas" else part[-1][1])
                chunks.append(sorted(key for key, _ in part))
                tags.append((source, cover))
        total_picked = sum(len(item[1]) for item in plan)
        print("picked %d variants in %d chunks (chunk size %d)" % (total_picked, len(chunks), args.chunk_size))

        def set_bookmark(source, cover):
            if source == "gwas":
                gwas_entry["index"] = cover
            else:
                genome_entry["rowid"] = cover

        if args.dry_run:
            for key in [k for chunk in chunks for k in chunk][:10]:
                print("  %s" % variant_id(key))
            for source, picked, end, start in plan:
                covers = [cover for tag_source, cover in tags if tag_source == source]
                target = covers[-1] if covers else end
                print("would set bookmark %s: %s" % (
                    source, target if source == "genome" else "%d/%d" % (target, len(entries))))
            check_cache_versions(deps, True, args.allow_rebuild)
            return EXIT_OK

        if chunks:
            code = check_cache_versions(deps, False, args.allow_rebuild)
            if code is not None:
                return code
        # nothing picked but rows scanned: advance without Hail
        advanced = False
        for source, picked, end, start in plan:
            if not picked and end > start:
                set_bookmark(source, end)
                advanced = True
        if advanced:
            try:
                save_state(args.state_file, state)
            except OSError as error:
                _stderr("[POPULATE] cannot save state %s: %s" % (args.state_file, error))
                return EXIT_USAGE
        if not chunks:
            _stderr("[POPULATE] nothing to cache: every requested source is exhausted or fully cached")
            _stderr("[POPULATE] " + _bookmark_line(state, maf))
            return EXIT_NOTHING

        def on_success(index):
            source, cover = tags[index - 1]
            set_bookmark(source, cover)
            save_state(args.state_file, state)

        try:
            if args.parallel == 1:
                stats = populate(chunks, deps, _stderr, time.time, on_success=on_success,
                                 stop_on_failure=True, verify_write=True)
            else:
                stats = populate_parallel(chunks, args.parallel, os.environ.get(DEPS_ENV), _stderr,
                                          time.time, on_success=on_success,
                                          stop_on_failure=True, verify_write=True)
        except KeyboardInterrupt:
            _stderr("[POPULATE] interrupted; finished chunks are cached and bookmarked, rerun to resume")
            _stderr("[POPULATE] " + _bookmark_line(state, maf))
            return EXIT_INTERRUPT
        except OSError as error:
            _stderr("[POPULATE] cannot save state %s: %s" % (args.state_file, error))
            return EXIT_USAGE
        _stderr("[POPULATE] done: chunks %d/%d ok, failed %d, queried %d, cached %d, "
                "not in gnomAD %d, elapsed %.0fs"
                % (stats["done"], stats["chunks"], stats["failed_chunks"], stats["queried"],
                   stats["cached"], stats["not_in_gnomad"], time.time() - started))
        _stderr("[POPULATE] " + _bookmark_line(state, maf))
        return EXIT_CHUNK_FAILED if stats["failed_chunks"] else EXIT_OK
    finally:
        if lock is not None:
            lock.close()


def main(argv=None, deps=None):
    started = time.time()
    try:
        args = _parse(argv)
    except SystemExit as exit_:
        return exit_.code if isinstance(exit_.code, int) else EXIT_USAGE
    if args.add is not None:
        if deps is None:
            deps = load_pipeline()
        return run_add(args, deps, started)
    build = args.build or "hg38"
    try:
        lines = read_inputs(args.input)
    except (OSError, UnicodeDecodeError) as error:
        _stderr("[POPULATE] cannot read input %s: %s" % (args.input, error))
        return EXIT_USAGE
    if deps is None:
        deps = load_pipeline()

    keys, unresolved = resolve_all(lines, build, deps)
    if not keys:
        _stderr("[POPULATE] nothing resolvable (%d inputs, %d unresolved)" % (len(lines), len(unresolved)))
        return EXIT_NOTHING
    chunks = plan_chunks(keys, deps._hail_cache_read, args.chunk_size)
    misses = sum(len(chunk) for chunk in chunks)

    reasons = {}
    for _, reason in unresolved:
        reasons[reason] = reasons.get(reason, 0) + 1
    print("inputs: %d" % len(lines))
    print("resolved keys: %d" % len(keys))
    print("unresolved: %d %s" % (len(unresolved), dict(sorted(reasons.items()))))
    for item, reason in unresolved[:20]:
        print("  %s\t%s" % (item, reason))
    print("cache hits: %d" % (len(keys) - misses))
    print("cache misses: %d" % misses)
    print("planned chunks: %d (chunk size %d)" % (len(chunks), args.chunk_size))
    sys.stdout.flush()

    code = check_cache_versions(deps, args.dry_run, args.allow_rebuild)
    if code is not None:
        return code

    if args.dry_run or not chunks:
        if not args.dry_run:
            _stderr("[POPULATE] everything is already cached")
        return EXIT_OK

    try:
        if args.parallel == 1:
            stats = populate(chunks, deps, _stderr, time.time)
        else:
            stats = populate_parallel(chunks, args.parallel, os.environ.get(DEPS_ENV), _stderr,
                                      time.time, stop_on_failure=False, verify_write=False)
    except KeyboardInterrupt:
        _stderr("[POPULATE] interrupted; finished chunks are cached, rerun to resume")
        return EXIT_INTERRUPT
    _stderr("[POPULATE] done: chunks %d/%d ok, failed %d, queried %d, cached %d, "
            "not in gnomAD %d, elapsed %.0fs"
            % (stats["done"], stats["chunks"], stats["failed_chunks"], stats["queried"],
               stats["cached"], stats["not_in_gnomad"], time.time() - started))
    return EXIT_CHUNK_FAILED if stats["failed_chunks"] else EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
