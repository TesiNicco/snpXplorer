# Handoff: fix silent false "not in gnomAD" in populate_gnomad_cache.py

For: the Claude session working on the compute server where `populate_gnomad_cache.py` runs.
From: the session that wrote the script (on the web server, which has no pysam).
Status: **not fixed yet**. The user asked for this note so you can make targeted changes there.

## What the script is

`populate_gnomad_cache.py` (same folder) fills a standalone SQLite cache (`gnomad_cache`, key `"C-P-R-A"` on hg38) with gnomAD v4.1 data. It reads the public sites VCFs remotely via `pysam.TabixFile(url, index=<local .tbi>)`, using one handle per (dataset, chrom) per worker thread (`GnomadSource._handle`, `threading.local`). It stores:

- **Frequencies:** combined exome+genome AC/AN/AF, plus 7 ancestry AFs.
- **Other columns:** nhomalt, `in_exomes`/`in_genomes`, 7 prediction scores, and VEP for the MANE/canonical transcript.

Variants found in neither file are **not stored**. The bookmark (`<out>.state.json`) moves past them, so they are never retried. That is by design for genuinely absent variants.

Constraints to keep:

- Python ≥3.8, stdlib plus pysam only.
- Same schema, key format and meta (`schema_version=1`, `gnomad_version=4.1`).
- Same bookmark semantics (contiguous prefix, advanced only after a chunk is committed).
- Same exit codes (0/1/2/3/4/130).

## What happened (user's run: `--gwas ... --add 10000 --threads 32`)

- Chunks 1–3 behaved normally: 2000/2000 stored, **0 not in gnomAD**, about 24–34 variants/s.
- Chunk 4 slowed to 11 variants/s. Then:
  - libcurl reported `error 92 (Stream error in the HTTP/2 framing layer)` (several times),
  - followed by hundreds of `[E::hts_itr_next] Failed to seek to offset <huge number>: Illegal seek`,
  - and `[E::bgzf_read_block] Failed to read BGZF header at offset ...`.
- Chunk 5 reported `stored 1797, not in gnomAD 203, failed 0`. Earlier GWAS chunks had 0 not-found.
- At exit: `Libcurl reported error 43 (A libcurl function was given a bad argument)` ×5. This came while closing handles and is harmless.

## Diagnosis

After the HTTP/2 stream error, some per-thread handles became corrupt. Every later `handle.fetch(...)` on them logs an htslib error to stderr but **raises no Python exception and yields no records**.

In `GnomadSource._find` (around line 346), "no matching record" returns `None` (not found). The retry/reopen path only runs on `OSError/IOError/ValueError`, so a corrupt handle produces **silent false negatives**:

1. Both datasets empty → variant counted "not in gnomAD", not stored, and skipped forever by the bookmark. These are the 203.
2. **Worse, and invisible:** one dataset's handle is corrupt and the other is fine → the row is stored with that dataset missing. `in_exomes`/`in_genomes` is wrong, and the combined AC/AN/AF is computed from one dataset only, so the **AF is wrong**.

Variants looked up before the error are fine. A corrupt handle yields nothing; it does not yield wrong records.

## Fix requested

1. **Confirm every per-dataset miss with a second, independent handle.**
   - Keep a second persistent handle per (dataset, chrom) per thread (a "verify" handle) and re-fetch on a miss.
   - If the verify handle finds the record, the primary handle was corrupt: drop it (close and reopen on next use) and use the record.
   - Only treat it as a miss when both agree.
   - Reasons for a second handle rather than reopening per miss:
     - Exome misses are normal for non-coding variants, so they are frequent.
     - Opening a handle costs a remote header read.
   - Optional: also replace both handles every N fetches as hygiene.
2. **Safety net per chunk.** If the fraction of keys "not in gnomAD" in a chunk exceeds a threshold, treat the chunk as **failed** (exit 3; bookmark not moved; rows still committed is fine). Reasonable default: 5%, as a CLI option.
   - Observed rates: GWAS mode 0%; `--variant-info` mode 2.1% genuine (21/1000, identical to the old Hail pipeline); the broken chunk 10%.
   - Also log how many misses were "recovered" by the verify handle.
3. **Keep threads moderate.** Suggest the default stays 8 and the user runs 16 rather than 32.
   - Measured on this compute server: 32 → 17–22 variants/s, 64 → 18.8 (no real gain).
   - The HTTP/2 errors appeared at 32.
4. **Tests.** Extend the existing fake-opener tests with a `CorruptTabix` that, after k fetches, returns no records without raising, for one dataset or both. Assert:
   - (a) no false "not found" is stored or skipped;
   - (b) no row is stored with a dataset wrongly missing;
   - (c) the threshold fails the chunk with exit 3 and the bookmark unchanged.

   The existing test suite lives on the web server at `/tmp/.../tests_step_gnomadcache/test_gnomad_cache.py` (not available to you). Write new tests next to the script or in a scratch dir, using the same injectable `opener` argument of `main(argv, opener=...)`.

## Data recovery after the fix

The database is small (22,116 rows at the time of writing), and rows from the faulty run can't be told apart: there is no timestamp column, and exome misses are normal. **Recommendation: rebuild from scratch** with the fixed script. Use a new `--out` file and a new `--state-file` (or delete both old files). At ~15–20 variants/s that is roughly 20–30 minutes for the same GWAS range. Then continue as before:

```bash
python populate_gnomad_cache.py --gwas Gwas_catalog_hg38.txt.gz --variant-db /path/variant_info.db \
    --tbi-dir gnomad_v4.1_tbi/ --out db/gnomad_cache_v2.sqlite --add 30000 --threads 16
```

Alternative if the user prefers to keep the file: delete rows from the affected GWAS range, then rerun with a fresh `--state-file`. Already-stored rows are skipped by the cache check, so only missing ones are looked up. But rows with a wrongly missing dataset (point 2 of the diagnosis) can't be identified reliably, so a rebuild is cleaner.

## Other facts you may need

- Remote VCFs: `https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/{exomes,genomes}/gnomad.{exomes,genomes}.v4.1.sites.chr{1..22,X,Y}.vcf.bgz`. Local `.tbi` files live in `--tbi-dir/{exomes,genomes}/`.
- The combined AF rule matched the old Hail-based cache exactly (480/480 variants tested on the web server).
- `--variant-db` must be the same `variant_info.db` the web server uses (`Data/databases/Genes/variant_info.db`, 12.6 GB). Otherwise rsID → allele resolution can differ from real annotation jobs. The user is running against `/project/holstegelab/.../bravo_db/variant_info.db`; check it's identical (size or md5).
- Org policy for this project: no personally identifiable data in outputs or logs. The cache holds only public variant data.
