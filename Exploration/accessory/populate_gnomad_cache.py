#!/usr/bin/env python3
"""Pre-populate a local gnomAD v4.1 cache (gnomad_cache.sqlite) from the public VCFs.

Purpose
    For a list of variants (GWAS catalogue, or the whole variant_info table) look
    each one up in the public gnomAD v4.1 exomes and genomes sites VCFs, using
    remote tabix queries (pysam) with LOCAL .tbi index files, and store allele
    frequencies, predictions and one VEP transcript per variant in a SQLite
    file.  Only the cache file is written (plus a bookmark file and a lock
    file); progress is logged to stderr.  Standalone: stdlib + (lazily) pysam.

Files
    {base-url}/{exomes|genomes}/gnomad.{dataset}.v4.1.sites.chr{C}.vcf.bgz
    {tbi-dir}/{exomes|genomes}/gnomad.{dataset}.v4.1.sites.chr{C}.vcf.bgz.tbi
    The default base URL is
    https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf
    A missing local .tbi for a needed chromosome is a usage error (exit 1).
    Chromosome MT is skipped (counted as filtered): v4.1 has no MT sites VCF.

Smoke test (40 variants, 2 chunks of 20)
    python3 populate_gnomad_cache.py --gwas GWAS.txt.gz --variant-db variant_info.db --tbi-dir TBI_DIR --out gnomad_cache.sqlite --add 40 --chunk-size 20 --threads 8 --dry-run
    python3 populate_gnomad_cache.py --gwas GWAS.txt.gz --variant-db variant_info.db --tbi-dir TBI_DIR --out gnomad_cache.sqlite --add 40 --chunk-size 20 --threads 8
    python3 populate_gnomad_cache.py --compare-with hail_cache.sqlite 200 --tbi-dir TBI_DIR --out gnomad_cache.sqlite

Whole genome source instead of the GWAS catalogue
    python3 populate_gnomad_cache.py --variant-info --variant-db variant_info.db --tbi-dir TBI_DIR --out gnomad_cache.sqlite --add 100000

Speed: about 3.5 variants/s with --threads 8 (measured), i.e. roughly 1 hour
per 12000 variants.

Cached values (table gnomad_cache, key variant_key = CHR-POS-REF-ALT, hg38)
    total_af and af_<african|european_non_finnish|east_asian|middle_east|finnish|
    south_asian|admixed_american>: AC/AN summed over the exomes and genomes
    records that exist (a missing key counts as 0/0; AN 0 -> NULL), the same
    rule as the annotation pipeline.  ac, an, nhomalt (summed), in_exomes,
    in_genomes, cadd_phred, revel_max, spliceai_ds_max, pangolin_largest_ds,
    phylop, sift_max, polyphen_max (exome value if present, else genome) and
    one VEP transcript (MANE select > canonical > first; exome VEP if the
    exome record has one, else genome): gene_symbol, gene_id, transcript_id,
    consequence, impact, biotype, hgvsc, hgvsp, mane_select, canonical, lof.
    Variants found in neither file are not stored (counted "not in gnomAD").

Behaviour
    - Resumable bookmark (--state-file, default <out>.state.json, atomic
      writes): after every successful chunk the bookmark moves to the end of
      the contiguous prefix of finished chunks.  A changed GWAS file
      (fingerprint) restarts the GWAS list; a new --min-maf starts a new
      bookmark.  Variants already in --out are skipped.
    - One run at a time: lock file <state-file>.lock (second run exits 1).
      --dry-run takes no lock, writes nothing and makes no remote call.
    - Lookups are retried up to 3 times per file (1 s, 4 s backoff, handle
      reopened).  A variant still failing makes its chunk fail: the successful
      variants of the chunk are stored, the bookmark does not move, exit 3.
    - Ctrl-C / SIGTERM: the running chunk is not committed, earlier chunks and
      bookmarks are kept, exit 130.
    - --compare-with HAIL_SQLITE N: look up N random cached variants (--seed)
      and compare the AF columns with the Hail cache (no writes).

Exit codes: 0 ok, 1 usage/unreadable input/lock/missing .tbi, 2 nothing to do
(source exhausted), 3 a chunk failed (or compare found a difference),
4 --out was built with another schema/gnomAD version, 130 interrupted.
"""
import argparse
import gzip
import hashlib
import json
import math
import os
import queue
import random
import re
import signal
import sqlite3
import sys
import tempfile
import threading
import time
import urllib.parse

try:  # POSIX only
    import fcntl
except ImportError:  # pragma: no cover
    fcntl = None

DEFAULT_BASE_URL = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf"
DEFAULT_CHUNK = 2000
DEFAULT_BATCH = 50000
DEFAULT_MIN_MAF = 0.01
DEFAULT_THREADS = 8
MAX_THREADS = 32
CACHE_CHECK_BATCH = 1000
PROGRESS_EVERY = 1000000
MAF_EPSILON = 1e-12
STATE_VERSION = 1
SCHEMA_VERSION = "1"
GNOMAD_VERSION = "4.1"
ATTEMPTS = 3
BACKOFF = (1, 4, 10)
DATASETS = ("exomes", "genomes")
VALID_CHROMS = {str(c) for c in range(1, 23)} | {"X", "Y", "MT"}
CACHEABLE_CHROMS = VALID_CHROMS - {"MT"}  # v4.1 has no MT sites VCF
EXIT_OK, EXIT_USAGE, EXIT_NOTHING, EXIT_CHUNK_FAILED, EXIT_VERSION, EXIT_INTERRUPT = 0, 1, 2, 3, 4, 130

HAIL_POPULATIONS = {
    "afr": "african", "nfe": "european_non_finnish",
    "eas": "east_asian", "mid": "middle_east", "fin": "finnish",
    "sas": "south_asian", "amr": "admixed_american",
}
AF_COLUMNS = ["total_af"] + ["af_%s" % label for label in HAIL_POPULATIONS.values()]
PREDICTIONS = ["cadd_phred", "revel_max", "spliceai_ds_max", "pangolin_largest_ds",
               "phylop", "sift_max", "polyphen_max"]
COLUMN_TYPES = (
    [("variant_key", "TEXT PRIMARY KEY")]
    + [(c, "REAL") for c in AF_COLUMNS]
    + [("ac", "INTEGER"), ("an", "INTEGER"), ("nhomalt", "INTEGER"),
       ("in_exomes", "INTEGER"), ("in_genomes", "INTEGER")]
    + [(c, "REAL") for c in PREDICTIONS]
    + [("gene_symbol", "TEXT"), ("gene_id", "TEXT"), ("transcript_id", "TEXT"), ("consequence", "TEXT"),
       ("impact", "TEXT"), ("biotype", "TEXT"), ("hgvsc", "TEXT"), ("hgvsp", "TEXT"),
       ("mane_select", "TEXT"), ("canonical", "INTEGER"), ("lof", "TEXT")]
)
COLUMNS = [name for name, _ in COLUMN_TYPES]
META = [("schema_version", SCHEMA_VERSION), ("gnomad_version", GNOMAD_VERSION),
        ("vep_transcript_rule", "mane_select>canonical>first;exome>genome"),
        ("created_by", "populate_gnomad_cache.py")]


def _stderr(message):
    print(message, file=sys.stderr, flush=True)


# ---------------------------------------------------------------- helpers (copied from the pipeline)
def normalize_chromosome(chrom):
    """Normalize a chromosome label to the representation used in variant_info."""
    normalized = str(chrom).strip().upper().replace("CHR", "")
    return "MT" if normalized == "M" else normalized


def variant_key(record):
    """Return a stable identity that preserves multiallelic variants."""
    return (
        normalize_chromosome(record["Chromosome"]),
        int(record["Position_hg38"]),
        str(record["REF"]).upper(),
        str(record["ALT"]).upper(),
    )


def variant_id(key):
    return "%s-%s-%s-%s" % (key[0], key[1], key[2], key[3])


def _frequency_parts(info, suffix):
    """(AC, AN) of one record for '' (total) or '_afr' etc.; missing/invalid -> 0, 0."""
    if info is None:
        return 0, 0
    try:
        ac = int(info.get("AC" + suffix) or 0)
    except (TypeError, ValueError):
        ac = 0
    try:
        an = int(info.get("AN" + suffix) or 0)
    except (TypeError, ValueError):
        an = 0
    return ac, an


def _combined_counts(exome, genome, suffix):
    """Sum AC and AN over the records that exist (missing key -> 0/0)."""
    counts = [_frequency_parts(item, suffix) for item in (exome, genome) if item is not None]
    return sum(c[0] for c in counts), sum(c[1] for c in counts)


def _combined_af(exome, genome, suffix):
    if exome is None and genome is None:
        return None
    ac, an = _combined_counts(exome, genome, suffix)
    return ac / an if an > 0 else None


# ---------------------------------------------------------------- VCF parsing
def parse_info(text):
    """INFO column -> dict; flags -> True, '.' -> None."""
    result = {}
    if not text or text == ".":
        return result
    for item in text.split(";"):
        if not item:
            continue
        name, sep, value = item.partition("=")
        if not sep:
            result[name] = True
        else:
            result[name] = None if value == "." else value
    return result


def parse_vep_format(header_lines):
    """Field names of the VEP INFO field from the ##INFO=<ID=vep ...> header line, or None."""
    for text in header_lines:
        if str(text).startswith("##INFO=<ID=vep,"):
            match = re.search(r'Format: ([^"]+)', text)
            if match:
                return [name.strip() for name in match.group(1).strip().split("|")]
    return None


def pick_transcript(entries):
    """First entry with MANE_SELECT, else first CANONICAL == YES, else the first; None if empty."""
    for entry in entries:
        if entry.get("MANE_SELECT"):
            return entry
    for entry in entries:
        if entry.get("CANONICAL") == "YES":
            return entry
    return entries[0] if entries else None


def _number(value):
    if value is None or value is True:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(number) else number


def _text(value):
    return value if value not in (None, "", ".") else None


def _first(*values):
    for value in values:
        if value is not None:
            return value
    return None


def _vep_transcript(record):
    """Chosen VEP entry dict of a record {'info', 'vep_format'}, or None."""
    if record is None or not record["vep_format"]:
        return None
    raw = record["info"].get("vep")
    if not raw or raw is True:
        return None
    fmt = record["vep_format"]
    entries = [dict(zip(fmt, part.split("|"))) for part in raw.split(",") if part]
    return pick_transcript(entries)


def build_row(key, exome, genome):
    """Cache row (dict of COLUMNS) from the exome/genome records ({'info', 'vep_format'} or None)."""
    ex = exome["info"] if exome else None
    ge = genome["info"] if genome else None
    row = {"variant_key": variant_id(key)}
    row["total_af"] = _combined_af(ex, ge, "")
    for ancestry, label in HAIL_POPULATIONS.items():
        row["af_%s" % label] = _combined_af(ex, ge, "_" + ancestry)
    row["ac"], row["an"] = _combined_counts(ex, ge, "")
    nhomalt = 0
    for info in (ex, ge):
        if info is not None:
            try:
                nhomalt += int(info.get("nhomalt") or 0)
            except (TypeError, ValueError):
                pass
    row["nhomalt"] = nhomalt
    row["in_exomes"] = 1 if exome else 0
    row["in_genomes"] = 1 if genome else 0
    for name in PREDICTIONS:
        row[name] = _number(_first((ex or {}).get(name), (ge or {}).get(name)))
    chosen = None
    for record in (exome, genome):
        if record is not None and record["info"].get("vep") not in (None, "", True):
            chosen = _vep_transcript(record)
            break
    names = {"gene_symbol": "SYMBOL", "gene_id": "Gene", "transcript_id": "Feature", "consequence": "Consequence",
             "impact": "IMPACT", "biotype": "BIOTYPE", "hgvsc": "HGVSc", "hgvsp": "HGVSp",
             "mane_select": "MANE_SELECT", "lof": "LoF"}
    for column, field in names.items():
        value = _text(chosen.get(field)) if chosen else None
        if value is not None and column in ("hgvsc", "hgvsp") and "%" in value:
            value = urllib.parse.unquote(value)
        row[column] = value
    row["canonical"] = (1 if chosen.get("CANONICAL") == "YES" else 0) if chosen else None
    return row


# ---------------------------------------------------------------- remote lookup
def default_opener(url, index):
    """Open a remote bgzipped VCF with a local tabix index (pysam is imported lazily)."""
    import pysam
    return pysam.TabixFile(url, index=index)


class LookupFailed(Exception):
    """A lookup still failed after all retries."""


class GnomadSource:
    """Looks up one variant in the exomes and genomes VCFs; one handle per (dataset, chrom) per thread."""

    def __init__(self, tbi_dir, opener, base_url):
        self.tbi_dir = str(tbi_dir)
        self.opener = opener
        self.base_url = str(base_url).rstrip("/")
        self._local = threading.local()
        self._lock = threading.Lock()
        self._formats = {}

    @staticmethod
    def file_name(dataset, chrom):
        return "gnomad.%s.v4.1.sites.chr%s.vcf.bgz" % (dataset, chrom)

    def url(self, dataset, chrom):
        return "%s/%s/%s" % (self.base_url, dataset, self.file_name(dataset, chrom))

    def index_path(self, dataset, chrom):
        return os.path.join(self.tbi_dir, dataset, self.file_name(dataset, chrom) + ".tbi")

    def _handle(self, dataset, chrom):
        handles = getattr(self._local, "handles", None)
        if handles is None:
            handles = self._local.handles = {}
        handle = handles.get((dataset, chrom))
        if handle is None:
            handle = handles[(dataset, chrom)] = self.opener(self.url(dataset, chrom), self.index_path(dataset, chrom))
        return handle

    def _drop(self, dataset, chrom):
        handle = getattr(self._local, "handles", {}).pop((dataset, chrom), None)
        if handle is not None:
            try:
                handle.close()
            except Exception:
                pass

    def _vep_format(self, dataset, chrom, handle):
        url = self.url(dataset, chrom)
        with self._lock:
            if url not in self._formats:
                self._formats[url] = parse_vep_format(list(handle.header))
            return self._formats[url]

    def _find(self, dataset, key):
        chrom, pos, ref, alt = key
        last = None
        for attempt in range(ATTEMPTS):
            try:
                handle = self._handle(dataset, chrom)
                if "chr" + chrom not in list(handle.contigs):
                    return None
                fmt = self._vep_format(dataset, chrom, handle)
                for text in handle.fetch("chr" + chrom, pos - 1, pos):
                    fields = text.split("\t")
                    if (int(fields[1]) == pos and fields[3].upper() == ref and fields[4].upper() == alt):
                        return {"info": parse_info(fields[7] if len(fields) > 7 else "."), "vep_format": fmt}
                return None
            except (OSError, IOError, ValueError) as error:
                last = error
                self._drop(dataset, chrom)
                if attempt + 1 < ATTEMPTS:
                    time.sleep(BACKOFF[attempt])
        raise LookupFailed("%s %s: %s: %s" % (dataset, variant_id(key), type(last).__name__, last))

    def lookup(self, key):
        """Row dict for key, or None when the variant is in neither file."""
        exome = self._find("exomes", key)
        genome = self._find("genomes", key)
        if exome is None and genome is None:
            return None
        return build_row(key, exome, genome)


# ---------------------------------------------------------------- output database
class VersionMismatch(Exception):
    pass


def _ro_connect(path):
    uri = "file:" + urllib.parse.quote(str(path)) + "?mode=ro"
    return sqlite3.connect(uri, uri=True, timeout=30)


def check_output_meta(path):
    """Raise VersionMismatch unless path is absent/empty or has the current schema and gnomAD version."""
    if not os.path.exists(str(path)):
        return
    try:
        connection = _ro_connect(path)
        try:
            tables = {n for (n,) in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
            if not tables:
                return
            if "gnomad_cache_meta" not in tables:
                raise VersionMismatch("%s has no gnomad_cache_meta table" % path)
            meta = dict(connection.execute("SELECT key, value FROM gnomad_cache_meta"))
        finally:
            connection.close()
    except sqlite3.Error as error:
        raise VersionMismatch("%s is not a readable cache (%s)" % (path, error))
    if meta.get("schema_version") != SCHEMA_VERSION or meta.get("gnomad_version") != GNOMAD_VERSION:
        raise VersionMismatch("%s has schema_version=%s gnomad_version=%s, expected %s / %s" % (
            path, meta.get("schema_version"), meta.get("gnomad_version"), SCHEMA_VERSION, GNOMAD_VERSION))


def count_rows(path):
    """Number of variants in the output database, or None if it can't be read. Never raises."""
    if not os.path.exists(str(path)):
        return 0
    try:
        connection = _ro_connect(path)
        try:
            return connection.execute("SELECT count(*) FROM gnomad_cache").fetchone()[0]
        finally:
            connection.close()
    except sqlite3.Error:
        return None


def cache_has(path, ids):
    """Subset of ids present in the output cache (empty when it does not exist yet). Never raises."""
    ids = list(ids)
    if not ids or not os.path.exists(str(path)):
        return set()
    found = set()
    try:
        connection = _ro_connect(path)
        try:
            tables = {n for (n,) in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
            if "gnomad_cache" not in tables:
                return set()
            for start in range(0, len(ids), 500):
                part = ids[start:start + 500]
                marks = ",".join("?" * len(part))
                found.update(r[0] for r in connection.execute(
                    "SELECT variant_key FROM gnomad_cache WHERE variant_key IN (%s)" % marks, part))
        finally:
            connection.close()
    except sqlite3.Error as error:
        _stderr("[POPULATE] cache read failed: %s" % error)
        return set()
    return found


def store_rows(path, rows):
    """One transaction: INSERT OR REPLACE rows, then read them back; raises RuntimeError if not persisted."""
    if not rows:
        return
    directory = os.path.dirname(os.path.abspath(str(path)))
    os.makedirs(directory, exist_ok=True)
    connection = sqlite3.connect(str(path), timeout=30, isolation_level=None)
    try:
        connection.execute("PRAGMA busy_timeout=30000")
        connection.execute("PRAGMA journal_mode=WAL")
        connection.execute("BEGIN IMMEDIATE")
        try:
            connection.execute("CREATE TABLE IF NOT EXISTS gnomad_cache (%s) WITHOUT ROWID" % ", ".join(
                "%s %s" % item for item in COLUMN_TYPES))
            connection.execute("CREATE TABLE IF NOT EXISTS gnomad_cache_meta (key TEXT PRIMARY KEY, value TEXT)")
            connection.executemany("INSERT OR IGNORE INTO gnomad_cache_meta(key, value) VALUES (?, ?)", META)
            connection.executemany(
                "INSERT OR REPLACE INTO gnomad_cache(%s) VALUES (%s)" % (", ".join(COLUMNS), ",".join("?" * len(COLUMNS))),
                [tuple(row[c] for c in COLUMNS) for row in rows])
            connection.execute("COMMIT")
        except BaseException:
            connection.execute("ROLLBACK")
            raise
        keys = [row["variant_key"] for row in rows]
        count = 0
        for start in range(0, len(keys), 500):
            part = keys[start:start + 500]
            count += connection.execute("SELECT COUNT(*) FROM gnomad_cache WHERE variant_key IN (%s)" % ",".join("?" * len(part)),
                                        part).fetchone()[0]
        if count != len(set(keys)):
            raise RuntimeError("cache write not persisted")
    finally:
        connection.close()


# ---------------------------------------------------------------- state, lock, GWAS (copied from populate_hail_cache)
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
    except (OSError, ValueError) as error:
        raise StateError("state file %s unreadable (%s): fix or delete it" % (path, error))
    if (not isinstance(data, dict) or data.get("version") != STATE_VERSION
            or not isinstance(data.get("entries"), dict)
            or not all(isinstance(v, dict) for v in data["entries"].values())):
        raise StateError("state file %s has an unexpected structure: fix or delete it" % path)
    return data


def save_state(path, state):
    """Atomic save: temp file in the same directory, fsync, os.replace."""
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
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
        raise StateError("another run holds the lock %s" % lock_path)
    return handle


def load_gwas_entries(path):
    """Ordered unique entries of the hg38 GWAS catalogue: all rsIDs of the SNP column
    (lowercase), or 'Chr:Pos(hg38)' for rows without an rsID."""
    entries = []
    with gzip.open(str(path), "rt", encoding="utf-8", errors="replace") as source:
        for text in source:
            cols = text.rstrip("\r\n").split("\t")
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


def _fingerprint(entries):
    return hashlib.sha256("\n".join(entries).encode("utf-8")).hexdigest()


# ---------------------------------------------------------------- picking
def resolve(entry, conn):
    """Resolve a GWAS entry to a variant_info record (dict) like resolve_explicit_variant, or None."""
    if re.fullmatch(r"rs\d+", entry):
        query = ("SELECT Chromosome, Position_hg38, REF, ALT, AF FROM variant_info WHERE rsID = ? "
                 "ORDER BY Chromosome, Position_hg38, REF, ALT LIMIT 1")
        params = (entry,)
    else:
        match = re.fullmatch(r"([0-9]+|X|Y|MT):(\d+)", entry)
        if not match or normalize_chromosome(match.group(1)) not in VALID_CHROMS:
            return None
        query = ("SELECT Chromosome, Position_hg38, REF, ALT, AF FROM variant_info "
                 "WHERE Chromosome = ? AND Position_hg38 = ? ORDER BY REF, ALT LIMIT 1")
        params = (normalize_chromosome(match.group(1)), int(match.group(2)))
    row = conn.execute(query, params).fetchone()
    if row is None:
        return None
    return {"Chromosome": row[0], "Position_hg38": row[1], "REF": row[2], "ALT": row[3], "AF": row[4]}


def _new_counters():
    return {"scanned": 0, "filtered": 0, "cached": 0, "unresolved": 0,
            "duplicate": 0, "picked": 0, "exhausted": False}


def _maf_ok(af, maf):
    try:
        af = float(af)
    except (TypeError, ValueError):
        return False
    return af == af and min(af, 1.0 - af) >= maf - MAF_EPSILON


def pick_gwas(entries, start, n, maf, conn, has, seen, log=None):
    """Pick up to n uncached, MAF-passing variants from entries[start:].

    Returns (picked [(key, entry_index)], scan_end, counters); scan_end is the
    exclusive index of the last examined entry.  has(ids) -> set of cached ids.
    """
    counters = _new_counters()
    picked = []
    buffer = []  # (index, outcome, key)

    def flush():
        candidates = [item for item in buffer if item[1] == "cand"]
        cached = has([variant_id(item[2]) for item in candidates]) if candidates else set()
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
        outcome, key = "cand", None
        resolved = resolve(entries[index], conn)
        if resolved is None:
            outcome = "unresolved"
        elif not _maf_ok(resolved.get("AF"), maf):
            outcome = "filtered"
        else:
            try:
                key = variant_key(resolved)
            except (TypeError, ValueError, KeyError):
                outcome = "filtered"
            else:
                if key[0] not in CACHEABLE_CHROMS:
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


def pick_genome(start_rowid, n, maf, conn, has, batch, seen, log=None):
    """Pick up to n uncached, MAF-passing variants from variant_info with rowid > start_rowid.

    Returns (picked [(key, rowid)], scan_end rowid, counters).
    """
    counters = _new_counters()
    picked = []
    last = start_rowid
    next_log = (start_rowid // PROGRESS_EVERY + 1) * PROGRESS_EVERY
    sql = ("SELECT rowid, Chromosome, Position_hg38, REF, ALT, AF FROM variant_info "
           "WHERE rowid > ? AND AF >= ? AND AF <= ? ORDER BY rowid LIMIT ?")
    low, high = maf - MAF_EPSILON, 1.0 - maf + MAF_EPSILON
    while len(picked) < n:
        rows = conn.execute(sql, (last, low, high, batch)).fetchall()
        if not rows:
            counters["exhausted"] = True
            break
        items = []
        for rowid, chrom, pos, ref, alt, af in rows:
            record = {"Chromosome": chrom, "Position_hg38": pos, "REF": ref, "ALT": alt, "AF": af}
            try:
                key = variant_key(record)
            except (TypeError, ValueError, KeyError):
                key = None
            items.append((rowid, key if key is not None and key[0] in CACHEABLE_CHROMS else None))
        ids = [variant_id(key) for _, key in items if key is not None]
        cached = has(ids) if ids else set()
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
    return picked, last, counters


# ---------------------------------------------------------------- chunk runner
class ChunkRunner:
    """Daemon worker threads that look variants up; run_chunk() handles one chunk at a time."""

    def __init__(self, source, threads):
        self.source = source
        self.stop = threading.Event()
        self.tasks = queue.Queue()
        self.results = queue.Queue()
        self.workers = [threading.Thread(target=self._work, daemon=True) for _ in range(threads)]
        for worker in self.workers:
            worker.start()

    def _work(self):
        while not self.stop.is_set():
            key = self.tasks.get()
            if key is None or self.stop.is_set():
                return
            try:
                row = self.source.lookup(key)
                self.results.put(("ok" if row is not None else "none", key, row))
            except Exception as error:  # LookupFailed and anything unexpected: fail this variant only
                self.results.put(("fail", key, "%s: %s" % (type(error).__name__, error)))

    def run_chunk(self, keys):
        """Return (rows, not_found_count, [(key, error)]); KeyboardInterrupt propagates."""
        for key in keys:
            self.tasks.put(key)
        rows, missing, failed = [], 0, []
        for _ in keys:
            while True:
                try:
                    status, key, value = self.results.get(timeout=0.2)
                    break
                except queue.Empty:
                    continue
            if status == "ok":
                rows.append(value)
            elif status == "none":
                missing += 1
            else:
                failed.append((key, value))
        return rows, missing, failed

    def close(self):
        self.stop.set()
        for _ in self.workers:
            self.tasks.put(None)


def _raise_interrupt(signum, frame):
    raise KeyboardInterrupt


# ---------------------------------------------------------------- CLI
class _Parser(argparse.ArgumentParser):
    def error(self, message):
        self.print_usage(sys.stderr)
        _stderr("%s: error: %s" % (self.prog, message))
        raise SystemExit(EXIT_USAGE)


def _int_min1(value):
    try:
        number = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be an integer >= 1")
    if number < 1:
        raise argparse.ArgumentTypeError("must be an integer >= 1")
    return number


def _threads_type(value):
    message = "must be an integer between 1 and %d" % MAX_THREADS
    try:
        number = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(message)
    if not 1 <= number <= MAX_THREADS:
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
    parser = _Parser(description="Pre-populate a local gnomAD v4.1 cache from the public VCFs.")
    parser.add_argument("--gwas", metavar="FILE", help="hg38 GWAS catalogue (.txt.gz) as variant source")
    parser.add_argument("--variant-info", action="store_true", help="use the whole variant_info table (rowid order)")
    parser.add_argument("--variant-db", metavar="PATH", help="variant_info SQLite database (read-only)")
    parser.add_argument("--add", type=_int_min1, metavar="N", help="cache the next N uncached variants")
    parser.add_argument("--min-maf", type=_maf_type, default=DEFAULT_MIN_MAF,
                        help="minimum minor allele frequency on variant_info.AF (default 0.01, inclusive)")
    parser.add_argument("--threads", type=_threads_type, default=DEFAULT_THREADS, help="lookup threads (1-32, default 8)")
    parser.add_argument("--chunk-size", type=_int_min1, default=DEFAULT_CHUNK)
    parser.add_argument("--batch", type=_int_min1, default=DEFAULT_BATCH, help="rows per variant_info read")
    parser.add_argument("--tbi-dir", metavar="DIR", help="directory with exomes/ and genomes/ .tbi files")
    parser.add_argument("--base-url", default=DEFAULT_BASE_URL)
    parser.add_argument("--out", metavar="FILE", help="output SQLite cache")
    parser.add_argument("--state-file", default=None, help="bookmark file (default <out>.state.json)")
    parser.add_argument("--dry-run", action="store_true", help="no lock, no writes, no remote calls")
    parser.add_argument("--compare-with", nargs=2, metavar=("HAIL_SQLITE", "N"),
                        help="compare N random variants with a Hail cache (no writes)")
    parser.add_argument("--seed", type=int, default=1, help="random seed for --compare-with")
    args = parser.parse_args(argv)
    if not args.tbi_dir:
        parser.error("--tbi-dir is required")
    if args.compare_with:
        try:
            args.compare_n = _int_min1(args.compare_with[1])
        except argparse.ArgumentTypeError:
            parser.error("--compare-with N must be an integer >= 1")
        return args
    if bool(args.gwas) == bool(args.variant_info):
        parser.error("exactly one of --gwas FILE or --variant-info is required")
    if not args.variant_db:
        parser.error("--variant-db is required")
    if not args.out:
        parser.error("--out is required")
    if args.add is None:
        parser.error("--add is required (unless --compare-with)")
    if args.state_file is None:
        args.state_file = args.out + ".state.json"
    return args


def missing_tbi(source, chroms):
    return [source.index_path(ds, c) for c in sorted(chroms) for ds in DATASETS
            if not os.path.isfile(source.index_path(ds, c))]


def _print_counters(source, c):
    print("%s: scanned %d, filtered %d, already cached %d, unresolved %d, duplicates %d, picked %d%s" % (
        source, c["scanned"], c["filtered"], c["cached"], c["unresolved"], c["duplicate"], c["picked"],
        " (source exhausted)" if c["exhausted"] else ""))


def _bookmark_line(entry, source, maf):
    if source == "gwas":
        return "bookmark now at gwas %s/%s (maf=%.6g)" % (entry.get("index", 0), entry.get("total", "?"), maf)
    return "bookmark now at genome rowid %s (maf=%.6g)" % (entry.get("rowid", 0), maf)


def run_add(args, source, started):
    maf = args.min_maf
    kind = "gwas" if args.gwas else "genome"
    lock = None
    conn = None
    try:
        try:
            conn = _ro_connect(args.variant_db)
            conn.execute("SELECT 1 FROM variant_info LIMIT 1").fetchall()
        except sqlite3.Error as error:
            _stderr("[POPULATE] cannot read variant_info in %s: %s" % (args.variant_db, error))
            return EXIT_USAGE
        entries = None
        if kind == "gwas":
            try:
                entries = load_gwas_entries(args.gwas)
            except (OSError, EOFError, ValueError) as error:
                _stderr("[POPULATE] cannot read GWAS file %s: %s" % (args.gwas, error))
                return EXIT_USAGE
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

        existing = state["entries"]
        key = state_key(kind, maf)
        if key not in existing:
            if existing:
                _stderr("[POPULATE] WARNING: no bookmark for %s; starting from the beginning "
                        "(existing bookmarks: %s)" % (key, ", ".join(sorted(existing))))
            existing[key] = {"index": 0} if kind == "gwas" else {"rowid": 0}
        entry = existing[key]
        if kind == "gwas":
            fingerprint = _fingerprint(entries)
            if "fingerprint" in entry and entry["fingerprint"] != fingerprint:
                _stderr("[POPULATE] WARNING: GWAS file fingerprint changed; restarting the GWAS list at index 0")
                entry["index"] = 0
            entry["fingerprint"] = fingerprint
            entry["total"] = len(entries)
            entry["index"] = int(entry.get("index", 0))
        else:
            entry["rowid"] = int(entry.get("rowid", 0))

        def has(ids):
            return cache_has(args.out, ids)

        seen = set()
        try:
            if kind == "gwas":
                start = entry["index"]
                picked, end, counters = pick_gwas(entries, start, args.add, maf, conn, has, seen, _stderr)
            else:
                start = entry["rowid"]
                picked, end, counters = pick_genome(start, args.add, maf, conn, has, args.batch, seen, _stderr)
        except KeyboardInterrupt:
            _stderr("[POPULATE] interrupted while scanning; bookmark unchanged since the last good chunk")
            return EXIT_INTERRUPT
        except sqlite3.Error as error:
            _stderr("[POPULATE] variant_info read failed: %s" % error)
            return EXIT_USAGE
        _print_counters(kind, counters)
        sys.stdout.flush()

        chunks, covers = [], []
        for first in range(0, len(picked), args.chunk_size):
            part = picked[first:first + args.chunk_size]
            last_chunk = first + args.chunk_size >= len(picked)
            covers.append(end if last_chunk else (part[-1][1] + 1 if kind == "gwas" else part[-1][1]))
            chunks.append(sorted(k for k, _ in part))
        print("picked %d variants in %d chunks (chunk size %d)" % (len(picked), len(chunks), args.chunk_size))

        missing = missing_tbi(source, {k[0] for chunk in chunks for k in chunk})
        if missing:
            _stderr("[POPULATE] missing local tbi index file(s): %s" % ", ".join(missing[:6]))
            return EXIT_USAGE

        def set_bookmark(cover):
            entry["index" if kind == "gwas" else "rowid"] = cover

        if args.dry_run:
            for k in [k for chunk in chunks for k in chunk][:10]:
                print("  %s" % variant_id(k))
            target = covers[-1] if covers else end
            print("would set bookmark %s: %s" % (kind, "%d/%d" % (target, len(entries)) if kind == "gwas" else target))
            try:
                check_output_meta(args.out)
            except VersionMismatch as error:
                print("WARNING: %s (a real run would abort with exit 4)" % error)
            return EXIT_OK

        if chunks:
            try:
                check_output_meta(args.out)
            except VersionMismatch as error:
                _stderr("[POPULATE] ABORT: %s" % error)
                return EXIT_VERSION
        if not picked and end > start:
            set_bookmark(end)
            try:
                save_state(args.state_file, state)
            except OSError as error:
                _stderr("[POPULATE] cannot save state %s: %s" % (args.state_file, error))
                return EXIT_USAGE
        if not chunks:
            _stderr("[POPULATE] nothing to cache: the source is exhausted or fully cached")
            _stderr("[POPULATE] " + _bookmark_line(entry, kind, maf))
            return EXIT_NOTHING

        runner = ChunkRunner(source, args.threads)
        totals = {"ok": 0, "looked": 0, "stored": 0, "missing": 0, "failed": 0}
        exit_code = EXIT_OK
        try:
            for index, chunk in enumerate(chunks, 1):
                t0 = time.time()
                rows, missing_count, failed = runner.run_chunk(chunk)
                try:
                    store_rows(args.out, rows)
                except (sqlite3.Error, OSError, RuntimeError) as error:
                    _stderr("[POPULATE] chunk %d/%d store failed: %s" % (index, len(chunks), error))
                    totals["failed"] += len(chunk)
                    exit_code = EXIT_CHUNK_FAILED
                    break
                elapsed = time.time() - t0
                totals["looked"] += len(chunk)
                totals["stored"] += len(rows)
                totals["missing"] += missing_count
                totals["failed"] += len(failed)
                _stderr("chunk %d/%d: looked up %d, stored %d, not in gnomAD %d, failed %d, elapsed %.0fs, %.1f variants/s" % (
                    index, len(chunks), len(chunk), len(rows), missing_count, len(failed), elapsed,
                    len(chunk) / elapsed if elapsed > 0 else 0.0))
                if failed:
                    for failed_key, message in failed[:5]:
                        _stderr("[POPULATE] FAILED %s: %s" % (variant_id(failed_key), message))
                    exit_code = EXIT_CHUNK_FAILED
                    break
                totals["ok"] += 1
                set_bookmark(covers[index - 1])
                save_state(args.state_file, state)
        except KeyboardInterrupt:
            runner.close()
            _stderr("[POPULATE] interrupted; finished chunks are cached and bookmarked, rerun to resume")
            _stderr("[POPULATE] " + _bookmark_line(entry, kind, maf))
            return EXIT_INTERRUPT
        except OSError as error:
            runner.close()
            _stderr("[POPULATE] cannot save state %s: %s" % (args.state_file, error))
            return EXIT_USAGE
        runner.close()
        _stderr("[POPULATE] done: chunks %d/%d ok, looked up %d, stored %d, not in gnomAD %d, failed %d, elapsed %.0fs" % (
            totals["ok"], len(chunks), totals["looked"], totals["stored"], totals["missing"], totals["failed"],
            time.time() - started))
        _stderr("[POPULATE] " + _bookmark_line(entry, kind, maf))
        total = count_rows(args.out)
        _stderr("[POPULATE] database now holds %s variants" % ("?" if total is None else total))
        return exit_code
    finally:
        if conn is not None:
            conn.close()
        if lock is not None:
            lock.close()


# ---------------------------------------------------------------- compare mode
def run_compare(args, source):
    hail_path, n = args.compare_with[0], args.compare_n
    transcript_cols = ["gene_symbol", "major_consequence", "transcript_id"]
    try:
        connection = _ro_connect(hail_path)
        try:
            keys = sorted(r[0] for r in connection.execute("SELECT variant_key FROM hail_cache WHERE total_af IS NOT NULL"))
            sample = random.Random(args.seed).sample(keys, min(n, len(keys)))
            columns = ["variant_key"] + AF_COLUMNS + transcript_cols
            hail = {}
            for start in range(0, len(sample), 500):
                part = sample[start:start + 500]
                for row in connection.execute("SELECT %s FROM hail_cache WHERE variant_key IN (%s)" % (
                        ", ".join(columns), ",".join("?" * len(part))), part):
                    hail[row[0]] = dict(zip(columns, row))
        finally:
            connection.close()
    except sqlite3.Error as error:
        _stderr("[COMPARE] cannot read hail cache %s: %s" % (hail_path, error))
        return EXIT_USAGE
    lookup_keys, skipped = [], 0
    for vid in sample:
        chrom, pos, ref, alt = vid.split("-", 3)
        if chrom in CACHEABLE_CHROMS:
            lookup_keys.append((chrom, int(pos), ref, alt))
        else:
            skipped += 1
    missing = missing_tbi(source, {k[0] for k in lookup_keys})
    if missing:
        _stderr("[COMPARE] missing local tbi index file(s): %s" % ", ".join(missing[:6]))
        return EXIT_USAGE
    runner = ChunkRunner(source, args.threads)
    old_term = None
    try:
        try:
            old_term = signal.signal(signal.SIGTERM, _raise_interrupt)
        except ValueError:
            pass
        try:
            rows, not_found, failed = runner.run_chunk(lookup_keys)
        except KeyboardInterrupt:
            _stderr("[COMPARE] interrupted")
            return EXIT_INTERRUPT
    finally:
        runner.close()
        if old_term is not None:
            signal.signal(signal.SIGTERM, old_term)
    print("sampled %d of %d hail variants with total_af (seed %d)" % (len(sample), len(keys), args.seed))
    print("looked up %d, found %d, failed %d, skipped (MT) %d" % (len(lookup_keys), len(rows), len(failed), skipped))
    print("not found in gnomAD: %d" % not_found)
    mismatch = bool(failed)
    for column in AF_COLUMNS:
        match, worst = 0, 0.0
        for row in rows:
            a, b = row[column], hail[row["variant_key"]][column]
            if a == b:
                match += 1
            elif a is not None and b is not None:
                worst = max(worst, abs(a - b))
        print("%s: %d/%d match, worst abs diff %.3g" % (column, match, len(rows), worst))
        mismatch = mismatch or match != len(rows)
    for column, mine in (("gene_symbol", "gene_symbol"), ("major_consequence", "consequence"), ("transcript_id", "transcript_id")):
        total = match = 0
        for row in rows:
            theirs = hail[row["variant_key"]][column]
            if theirs is None:
                continue
            total += 1
            value = row[mine]
            if column == "major_consequence" and value is not None:
                value = value.split("&")[0]
            match += value == theirs
        print("%s: %d/%d match (transcript field, informational)" % (column, match, total))
    for key, message in failed[:5]:
        _stderr("[COMPARE] FAILED %s: %s" % (variant_id(key), message))
    return EXIT_CHUNK_FAILED if mismatch else EXIT_OK


def main(argv=None, opener=None):
    started = time.time()
    try:
        args = _parse(argv)
    except SystemExit as exit_:
        return exit_.code if isinstance(exit_.code, int) else EXIT_USAGE
    source = GnomadSource(args.tbi_dir, opener or default_opener, args.base_url)
    if args.compare_with:
        return run_compare(args, source)
    old_term = None
    if not args.dry_run:
        try:
            old_term = signal.signal(signal.SIGTERM, _raise_interrupt)
        except ValueError:  # not the main thread
            pass
    try:
        return run_add(args, source, started)
    finally:
        if old_term is not None:
            signal.signal(signal.SIGTERM, old_term)


if __name__ == "__main__":
    sys.exit(main())
