#!/usr/bin/env python3
import sqlite3
import re
from pathlib import Path
import pandas as pd
import io
import subprocess
import math
from typing import Optional, List
from liftover import get_lifter
import argparse
from gprofiler import GProfiler
import plotly.graph_objects as go
import sys
import pickle
import random
from functools import reduce
from pygosemsim import download
from pygosemsim import graph
from pygosemsim import similarity
import networkx as nx

# ---------------------------------------------------------
# Configuration
# ---------------------------------------------------------
DATA_PATH = Path("/Users/nicco/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/snpXplorer/Data")
DB_FILE = DATA_PATH / "databases/Genes/variants_info.sqlite"

# Liftover object (hg19 → hg38)
lifter = get_lifter("hg19", "hg38")

# ---------------------------------------------------------
# Helper to sanitize (kept in case you reuse)
# ---------------------------------------------------------
def sanitize_for_json(obj):
    """
    Recursively convert NaN/Inf to None.
    Not strictly needed in CLI mode, but kept for completeness.
    """
    if isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            return None
        return obj
    if isinstance(obj, dict):
        return {k: sanitize_for_json(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [sanitize_for_json(v) for v in obj]
    return obj

# ---------------------------------------------------------
# DB helper: lookup by rsID
# ---------------------------------------------------------
def lookup_by_rsid(rsid: str):
    conn = sqlite3.connect(str(DB_FILE))
    cur = conn.cursor()
    cur.execute(
        """
        SELECT chr_hg38, position, position_hg19, marker_id, maf
        FROM rsids
        WHERE rsid = ?
        LIMIT 1
        """,
        (rsid,),
    )
    row = cur.fetchone()
    conn.close()
    if row is None:
        return None
    chr_hg38, pos38, pos19, marker_id, maf = row
    ref, alt = marker_id.split(":")[-2::]
    return {
        "query_type": "rsid",
        "rsid": rsid,
        "chr_hg38": chr_hg38,
        "pos_hg38": pos38,
        "pos_hg19": pos19,
        "ref": ref,
        "alt": alt,
        "maf": maf,
    }

# ---------------------------------------------------------
# DB helper: lookup by CHR + POS (hg38)
# ---------------------------------------------------------
def lookup_by_coord_hg38(chrom: str, pos38: int):
    conn = sqlite3.connect(str(DB_FILE))
    cur = conn.cursor()
    cur.execute(
        """
        SELECT rsid, position_hg19, marker_id, maf
        FROM rsids
        WHERE chr_hg38 = ? AND position = ?
        LIMIT 1
        """,
        (chrom, pos38),
    )
    row = cur.fetchone()
    conn.close()
    if row is None:
        return None
    rsid, pos19, marker_id, maf = row
    ref, alt = marker_id.split(":")[-2::]
    return {
        "query_type": "coord",
        "input_chr": chrom,
        "input_pos_hg38": pos38,
        "rsid": rsid,
        "chr_hg38": chrom,
        "pos_hg38": pos38,
        "pos_hg19": pos19,
        "ref": ref,
        "alt": alt,
        "maf": maf,
    }

# ---------------------------------------------------------
# Liftover wrapper using get_lifter()
# ---------------------------------------------------------
def liftover_hg19_to_hg38(chrom: str, pos19: int):
    """
    Convert hg19 → hg38 using SNPxplorer's internal liftover tools.
    lifter.convert_coordinate returns a list of (chr, pos, strand).
    """
    try:
        lifted = lifter.convert_coordinate(chrom, int(pos19))
    except Exception:
        return None
    if not lifted:
        return None
    chr38, pos38, _strand = lifted[0]
    return chr38, int(pos38)

# ---------------------------------------------------------
# Variant type detection
# ---------------------------------------------------------
def parse_variant_query(q: str):
    q = q.strip()
    snpinfo = []
    # file input
    if Path(q).is_file():
        with open(q, "r") as f:
            lines = f.readlines()
        queries = [line.strip() for line in lines if line.strip()]
    else:
        queries = [q]
    # for each variant, validate format
    for snp in queries:
        if snp.lower().startswith("rs") and snp[2:].isdigit():
            snpinfo.append({"type": "rsid", "rsid": snp})
        else:
            m = re.match(r"^(chr)?(\w+)[\:\s]+(\d+)$", snp, flags=re.IGNORECASE)
            if m:
                chrom = m.group(2).upper()
                pos = int(m.group(3))
                if 'chr' not in chrom.lower():
                    chrom = 'chr' + chrom
                snpinfo.append({"type": "coord", "chrom": chrom, "pos": pos, 'query': snp})
            else:
                snpinfo.append({"type": "invalid", "query": snp})
    return snpinfo

# ---------------------------------------------------------
# Query CADD score
# ---------------------------------------------------------
def query_cadd_score(chrom, pos38):
    """
    Query CADD annotation for a given hg38 chromosome and position.
    Returns a list of dictionaries (JSON-serializable).
    """
    try:
        chrom = str(chrom).upper().replace("CHR", "")
        CADD_DB_FILE = DATA_PATH / f"databases/CADD/subset_chr{chrom}.db"
        conn = sqlite3.connect(CADD_DB_FILE)
        df = pd.read_sql_query(
            "SELECT * FROM variants_combined WHERE pos = ?;",
            conn,
            params=(pos38,),
        )
        conn.close()
    except Exception:
        return []
    if df.empty:
        return []
    df = df.reset_index(drop=True)
    df = df.sort_values(by="phred_max", ascending=False)
    return df.to_dict(orient="records")

# ---------------------------------------------------------
# Query eQTLs
# ---------------------------------------------------------
def query_eqtls(chrom, pos38):
    """
    Query eQTL annotation for a given hg38 chromosome and position.
    Returns a list of dictionaries (JSON-serializable).
    """
    try:
        chrom_clean = str(chrom).upper().replace("CHR", "")
        EQTL_DB_FILE = DATA_PATH / f"databases/eQTLs/chr{chrom_clean}_eQTL.txt.gz"
        region = f"{chrom_clean}:{pos38}-{pos38}"
        result = subprocess.run(
            ["tabix", str(EQTL_DB_FILE), region],
            capture_output=True,
            check=True,
            text=True
        )
        output = result.stdout.strip()
        if output == "":
            return []
        df = pd.read_csv(io.StringIO(output), sep="\t", header=None)
        df.columns = [
            "chrom",
            "pos",
            "ref",
            "alt",
            "ensemble",
            "gene",
            "tissue",
            "tss_distance",
            "maf",
            "pval_nominal",
            "slope",
        ]
        return df.to_dict(orient="records")
    except Exception:
        return []

# ---------------------------------------------------------
# Query sQTLs
# ---------------------------------------------------------
def query_sqtls(chrom, pos38):
    """
    Query sQTL annotation for a given hg38 chromosome and position.
    Returns a list of dictionaries (JSON-serializable).
    """
    try:
        chrom_clean = str(chrom).upper().replace("CHR", "")
        SQTL_DB_FILE = DATA_PATH / f"databases/sQTLs/chr{chrom_clean}_sQTL.txt.gz"
        region = f"{chrom_clean}:{pos38}-{pos38}"
        result = subprocess.run(
            ["tabix", str(SQTL_DB_FILE), region],
            capture_output=True,
            check=True,
            text=True
        )
        output = result.stdout.strip()
        if output == "":
            return []
        df = pd.read_csv(io.StringIO(output), sep="\t", header=None)
        df.columns = [
            "chrom",
            "pos",
            "ref",
            "alt",
            "ensemble",
            "gene",
            "tissue",
            "tss_distance",
            "maf",
            "pval_nominal",
            "slope",
        ]
        return df.to_dict(orient="records")
    except Exception:
        return []

# ---------------------------------------------------------
# LD helpers
# ---------------------------------------------------------
def _open_db_for_query(db_path: str) -> sqlite3.Connection:
    MMAP_BYTES = 1 << 30  # 1 GiB mmap
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=NORMAL;")
    conn.execute(f"PRAGMA mmap_size={MMAP_BYTES};")
    return conn

def partners_by_position(db_path: str, pos: int, r2_min: float = 0.0, limit: Optional[int] = None):
    """
    Resolve LD partners by exact match on 'pos' column of the variants table.
    Returns rows (partner_uniq, r2, dist_bp).
    """
    R2_SCALE = 1000
    r2m_min = int(round(max(0.0, min(1.0, r2_min)) * R2_SCALE))
    sql = f"""
    WITH v AS (
      SELECT id FROM variants WHERE pos = ?
    ),
    hits AS (
      SELECT
        CASE WHEN ld.v1 = v.id THEN ld.v2 ELSE ld.v1 END AS partner_id,
        ld.r2_milli, ld.dist_bp
      FROM ld, v
      WHERE (ld.v1 = v.id OR ld.v2 = v.id)
        AND ld.r2_milli >= ?
    )
    SELECT variants.uniq AS partner_uniq,
           CAST(hits.r2_milli AS REAL)/{R2_SCALE} AS r2,
           hits.dist_bp
    FROM hits
    JOIN variants ON variants.id = hits.partner_id
    GROUP BY partner_uniq
    ORDER BY r2 DESC, dist_bp ASC
    {f"LIMIT {int(limit)}" if limit else ""}
    """
    with _open_db_for_query(str(db_path)) as conn:
        return conn.execute(sql, (pos, r2m_min)).fetchall()

# ---------------------------------------------------------
# Query LD
# ---------------------------------------------------------
def query_ld(chrom, pos38):
    """
    Query LD partners for given chromosome and position.
    Returns:
      ld_records: list[dict]
      all_partner_positions: list[int]
    """
    try:
        chrom_clean = str(chrom).upper().replace("CHR", "")
        db_path = DATA_PATH / f"databases/LD_db/ld_chr{chrom_clean}.sqlite"
        ld_partners = partners_by_position(db_path, int(pos38), r2_min=0.2)
        if not ld_partners:
            return [], []
        # partner_uniq, r2, dist_bp
        df_ld = pd.DataFrame(ld_partners, columns=["partner_uniq", "r2", "dist_bp"])
        df_ld["chr"] = chrom_clean
        df_ld["target"] = int(pos38)
        # Extract numeric partner positions safely
        df_ld["partner_pos"] = (
            df_ld["partner_uniq"]
            .astype(str)
            .str.split(":")
            .str[0]
            .apply(lambda x: pd.to_numeric(x, errors="coerce"))
        )
        # Drop rows where partner_pos is NaN
        df_ld = df_ld.dropna(subset=["partner_pos"])
        df_ld["partner_pos"] = df_ld["partner_pos"].astype(int)
        # Collect partner positions
        all_partner_positions = df_ld["partner_pos"].tolist()
        # Lookup rsIDs for partners
        records = []
        for v in all_partner_positions:
            info = lookup_by_coord_hg38(chrom, v)
            if info is not None:
                records.append({"partner_pos": int(info["pos_hg38"]), "rsid": info["rsid"]})
        if records:
            df_rs = pd.DataFrame(records).drop_duplicates(subset=["partner_pos"])
            # Ensure same dtype
            df_rs["partner_pos"] = df_rs["partner_pos"].astype(int)
            df_ld["partner_pos"] = df_ld["partner_pos"].astype(int)
            # Now merging is safe
            df_ld = df_ld.merge(df_rs, on="partner_pos", how="left")
        else:
            df_ld["rsid"] = None
        # Keep only useful columns
        df_ld = df_ld[["chr", "partner_pos", "rsid", "partner_uniq", "r2", "dist_bp"]]
        df_ld = df_ld[df_ld["r2"] >= 0.2]
        # Add query position
        df_ld['query_pos'] = int(pos38)
        if df_ld.empty:
            return [], []
        return df_ld.to_dict(orient="records"), all_partner_positions
    except Exception as e:
        print("Error in query_ld:", e, file=sys.stderr, flush=True)
        return [], []

# ---------------------------------------------------------
# Query GWAS associations
# ---------------------------------------------------------
def query_gwas_associations(chrom, pos38):
    """
    Query GWAS associations for a given hg38 chromosome and position.
    Returns a list of dictionaries (JSON-serializable).
    """
    try:
        chrom_clean = str(chrom).upper().replace("CHR", "")
        GWAS_DB_FILE = DATA_PATH / "databases/haplotypes/All_indep_gwas_sumstats_AI_hg38.txt.gz"
        region = f"{chrom_clean}:{pos38}-{pos38}"
        result = subprocess.run(
            ["tabix", str(GWAS_DB_FILE), region],
            capture_output=True,
            check=True,
            text=True
        ).stdout.strip()
        if result == "":
            return []
        df = pd.read_csv(io.StringIO(result), sep="\t", header=None)
        df.columns = [
            'gwas_id', 'trait', 'chrom', 'pos_hg19', 'rsid',
            'ea', 'nea', 'eaf', 'beta', 'se', 'p', 'n',
            'study_id', 'position_hg38'
        ]
        # subset columns to return
        df = df[['trait', 'ea', 'nea', 'eaf', 'beta', 'se', 'p', 'n', 'gwas_id', 'rsid']]
        # 2 decimal places for beta and se and eaf
        df['beta'] = df['beta'].round(2)
        df['se'] = df['se'].round(2)
        df['eaf'] = df['eaf'].round(3)
        # 2 decimal places for p-value in scientific notation
        df['p'] = df['p'].apply(lambda x: f"{x:.2e}")
    except Exception:
        return []
    if df.empty:
        return []
    df = df.reset_index(drop=True)
    return df.to_dict(orient="records")

# ---------------------------------------------------------
# Helper for annotating LD partners with CADD, eQTL, sQTL
# ---------------------------------------------------------
def annotate_ld_partners(chrom, ld_rows, pos38):
    """
    For each LD partner, query CADD / eQTL / sQTL and flatten results.
    Returns three lists of dicts: (ld_cadd, ld_eqtl, ld_sqtl).
    Each row is augmented with LD metadata (partner_pos, rsid, r2, dist_bp).
    """
    ld_cadd = []
    ld_eqtl = []
    ld_sqtl = []
    chrom_clean = str(chrom).upper().replace("CHR", "")
    for rec in ld_rows:
        partner_pos = int(rec["partner_pos"])
        rsid = rec.get("rsid")
        r2 = rec.get("r2")
        dist = rec.get("dist_bp")
        # CADD for partner
        cadd_rows = query_cadd_score(chrom_clean, partner_pos)
        for row in cadd_rows:
            r = dict(row)
            r["ld_partner_pos"] = partner_pos
            r["ld_rsid"] = rsid
            r["ld_r2"] = r2
            r["ld_dist_bp"] = dist
            r['ld_query_pos'] = pos38
            ld_cadd.append(r)
        # eQTL for partner
        eqtl_rows = query_eqtls(chrom_clean, partner_pos)
        for row in eqtl_rows:
            r = dict(row)
            r["ld_partner_pos"] = partner_pos
            r["ld_rsid"] = rsid
            r["ld_r2"] = r2
            r["ld_dist_bp"] = dist
            r['ld_query_pos'] = pos38
            ld_eqtl.append(r)
        # sQTL for partner
        sqtl_rows = query_sqtls(chrom_clean, partner_pos)
        for row in sqtl_rows:
            r = dict(row)
            r["ld_partner_pos"] = partner_pos
            r["ld_rsid"] = rsid
            r["ld_r2"] = r2
            r["ld_dist_bp"] = dist
            r['ld_query_pos'] = pos38
            ld_sqtl.append(r)
    return ld_cadd, ld_eqtl, ld_sqtl

# ---------------------------------------------------------
# Core logic for querying a variant
# ---------------------------------------------------------
def run_variant_query(q: str, build: str = "hg38"):
    print(f"Running variant query for: {q} (build={build})", file=sys.stderr, flush=True)
    # Parse query
    parsed = parse_variant_query(q)
    # Define output
    output = {}
    # Iterate over parsed queries -- it's always a list
    for p in parsed:
        print(f"Processing parsed query: {p}", file=sys.stderr, flush=True)
        # Invalid query
        if p["type"] == "invalid":
            output[p['query']] = 'Error: invalid query format'
            continue
        # rsID path
        if p["type"] == "rsid":
            info = lookup_by_rsid(p["rsid"])
            if info is None:
                output[p['rsid']] = 'Error: rsID not found in DB'
            try:
                chr38 = info["chr_hg38"]
                pos38 = info["pos_hg38"]
                # Add CADD annotation
                cadd_dic = query_cadd_score(chr38, pos38)
                # Add eQTL annotation
                eqtl_dic = query_eqtls(chr38, pos38)
                # Add sQTL annotation
                sqtl_dic = query_sqtls(chr38, pos38)
                # Add LD annotation
                ld_dic, all_vars = query_ld(chr38, pos38)
                # Add GWAS associations
                gwas_dic = query_gwas_associations(chr38, pos38)
                # CADD / eQTL / sQTL for all LD partners
                if ld_dic:
                    ld_cadd, ld_eqtl, ld_sqtl = annotate_ld_partners(chr38, ld_dic, pos38)
                else:
                    ld_cadd, ld_eqtl, ld_sqtl = None, None, None
            except Exception:
                cadd_dic = None
                eqtl_dic = None
                sqtl_dic = None
                ld_dic = None
                gwas_dic = None
            # Convert to dfs
            output[p['rsid']] = []
            output[p['rsid']].append([info] if info is not None else [])
            output[p['rsid']].append(cadd_dic if cadd_dic is not None else [])
            output[p['rsid']].append(eqtl_dic if eqtl_dic is not None else [])
            output[p['rsid']].append(sqtl_dic if sqtl_dic is not None else [])
            output[p['rsid']].append(ld_dic if ld_dic is not None else [])
            output[p['rsid']].append(ld_cadd if ld_cadd is not None else [])
            output[p['rsid']].append(ld_eqtl if ld_eqtl is not None else [])
            output[p['rsid']].append(ld_sqtl if ld_sqtl is not None else [])
            output[p['rsid']].append(gwas_dic if gwas_dic is not None else [])
            continue
        # coordinate path
        chrom = p["chrom"]
        pos = p["pos"]
        if build == "hg38":
            info = lookup_by_coord_hg38(chrom, pos)
            if info is None:
                output[p['query']] = 'Error: coordinate not found in DB'
                continue
            try:
                chr38 = info["chr_hg38"]
                pos38 = info["pos_hg38"]
                # Add CADD annotation
                cadd_dic = query_cadd_score(chr38, pos38)
                # Add eQTL annotation
                eqtl_dic = query_eqtls(chr38, pos38)
                # Add sQTL annotation
                sqtl_dic = query_sqtls(chr38, pos38)
                # Add LD annotation
                ld_dic, all_vars = query_ld(chr38, pos38)
                # Add GWAS associations
                gwas_dic = query_gwas_associations(chr38, pos38)
                # CADD / eQTL / sQTL for all LD partners
                if ld_dic:
                    ld_cadd, ld_eqtl, ld_sqtl = annotate_ld_partners(chr38, ld_dic, pos38)
                else:
                    ld_cadd, ld_eqtl, ld_sqtl = None, None, None
            except Exception:
                cadd_dic = None
                eqtl_dic = None
                sqtl_dic = None
                ld_dic = None
                gwas_dic = None
                ld_cadd = None
                ld_eqtl = None
                ld_sqtl = None
            # Convert to dfs
            output[p['query']] = []
            output[p['query']].append([info] if info is not None else [])
            output[p['query']].append(cadd_dic if cadd_dic is not None else [])
            output[p['query']].append(eqtl_dic if eqtl_dic is not None else [])
            output[p['query']].append(sqtl_dic if sqtl_dic is not None else [])
            output[p['query']].append(ld_dic if ld_dic is not None else [])
            output[p['query']].append(ld_cadd if ld_cadd is not None else [])
            output[p['query']].append(ld_eqtl if ld_eqtl is not None else [])
            output[p['query']].append(ld_sqtl if ld_sqtl is not None else [])
            output[p['query']].append(gwas_dic if gwas_dic is not None else [])
            continue
        if build == "hg19":
            lifted = liftover_hg19_to_hg38(chrom, pos)
            if lifted is None:
                output[p['query']] = 'Error: liftover failed from hg19 to hg38'
                continue
            chr38, pos38 = lifted
            info = lookup_by_coord_hg38(chr38, pos38)
            if info is None:
                output[p['query']] = 'Error: coordinate not found in DB'
                continue
            info["liftover"] = {"from_hg19": f"{chrom}:{pos}", "to_hg38": f"{chr38}:{pos38}"}
            try:
                # Add CADD annotation based on hg38 location
                cadd_dic = query_cadd_score(chr38, pos38)
                # Add eQTL annotation based on hg38 location
                eqtl_dic = query_eqtls(chr38, pos38)
                # Add sQTL annotation based on hg38 location
                sqtl_dic = query_sqtls(chr38, pos38)
                # Add LD annotation based on hg38 location
                ld_dic, all_vars = query_ld(chr38, pos38)
                # Add GWAS associations
                gwas_dic = query_gwas_associations(chr38, pos38)
                # CADD / eQTL / sQTL for all LD partners
                if ld_dic:
                    ld_cadd, ld_eqtl, ld_sqtl = annotate_ld_partners(chr38, ld_dic, pos38)
                else:
                    ld_cadd, ld_eqtl, ld_sqtl = None, None, None
            except Exception:
                cadd_dic = None
                eqtl_dic = None
                sqtl_dic = None
                ld_dic = None
                gwas_dic = None
            # Convert to dfs
            output[p['query']] = []
            output[p['query']].append([info] if info is not None else [])
            output[p['query']].append(cadd_dic if cadd_dic is not None else [])
            output[p['query']].append(eqtl_dic if eqtl_dic is not None else [])
            output[p['query']].append(sqtl_dic if sqtl_dic is not None else [])
            output[p['query']].append(ld_dic if ld_dic is not None else [])
            output[p['query']].append(ld_cadd if ld_cadd is not None else [])
            output[p['query']].append(ld_eqtl if ld_eqtl is not None else [])
            output[p['query']].append(ld_sqtl if ld_sqtl is not None else [])
            output[p['query']].append(gwas_dic if gwas_dic is not None else [])
            continue
    return output

# ---------------------------------------------------------
# Build list of DataFrames from info dict
# ---------------------------------------------------------
def info_to_dataframes(info: dict) -> List[pd.DataFrame]:
    """
    Convert the 'info' dict returned by run_variant_query into
    a list of DataFrames in a fixed order.
    """
    # 1) variant summary as single-row DF
    variant_fields = [
        "query_type", "rsid", "chr_hg38", "pos_hg38", "pos_hg19",
        "ref", "alt", "maf"
    ]
    base_row = {k: info.get(k) for k in variant_fields}
    if "liftover" in info:
        # optionally attach liftover info as extra columns
        base_row["liftover_from_hg19"] = info["liftover"].get("from_hg19")
        base_row["liftover_to_hg38"] = info["liftover"].get("to_hg38")
    df_variant = pd.DataFrame([base_row])
    def to_df(key: str) -> pd.DataFrame:
        records = info.get(key, [])
        if not records:
            return pd.DataFrame()
        return pd.DataFrame(records)
    df_cadd = to_df("cadd")
    df_eqtl = to_df("eqtl")
    df_sqtl = to_df("sqtl")
    df_ld = to_df("ld")
    df_gwas = to_df("gwas")
    df_ld_cadd = to_df("ld_cadd")
    df_ld_eqtl = to_df("ld_eqtl")
    df_ld_sqtl = to_df("ld_sqtl")
    return [
        df_variant,
        df_cadd,
        df_eqtl,
        df_sqtl,
        df_ld,
        df_gwas,
        df_ld_cadd,
        df_ld_eqtl,
        df_ld_sqtl,
    ]

# ---------------------------------------------------------
# Identify closest gene
# ---------------------------------------------------------
def closest_gene(query_info):
    # get chromosome and position
    chrom = query_info["chrom"].values[0].lower()
    if 'chr' not in chrom:
        chrom = 'chr' + chrom
    start_pos = query_info["pos"].values[0] - 500000
    end_pos = query_info["pos"].values[0] + 500000
    # tabix command
    GENE_DB_FILE = DATA_PATH / "databases/Genes/genes_hg38.txt.gz"
    result = subprocess.run(["tabix", str(GENE_DB_FILE), f"{chrom}:{start_pos}-{end_pos}"], capture_output=True, check=True, text=True)
    output = result.stdout.strip()
    if output == "":
        return []
    df = pd.read_csv(io.StringIO(output), sep="\t", header=None)
    df.columns = ["transcript_id", "chrom", "strand", "tx_start", "tx_end", "cds_start", "cds_end", "exon_num", "exon_start", "exon_end", "gene_symbol"]
    # add gene size
    df["gene_size"] = df["tx_end"] - df["tx_start"]
    # remove duplicates based on gene_symbol, keeping the longest gene
    df = df.sort_values(by="gene_size", ascending=False).drop_duplicates(subset=["gene_symbol"])
    # add distance to variant from tx_start and tx_end, and keep minimum
    variant_pos = query_info["pos"].values[0]
    df["dist_to_variant"] = df.apply(lambda row: min(abs(row["tx_start"] - variant_pos), abs(row["tx_end"] - variant_pos)), axis=1)
    # sort by distance
    df = df.sort_values(by="dist_to_variant", ascending=True)
    # return closest gene symbol
    return [df.iloc[0]["gene_symbol"]]

# ---------------------------------------------------------
# Identify most likely gene
# ---------------------------------------------------------
def identify_most_likely_gene(info: dict):
    # Set container for genes
    likely_genes = {}
    # Iterate over snps in info
    for snp in info.keys():
        # Extract dataframes
        info_df, cadd_df, eqtl_df, sqtl_df, ld_df, gwas_df, ld_cadd_df, ld_eqtl_df, ld_sqtl_df = convert_info_dict_to_dfs({snp: info[snp]})
        # First check if variant is coding in CADD
        coding_rows = cadd_df[cadd_df["annotypes"].str.contains("coding", case=False, na=False)]
        coding_rows = coding_rows[~coding_rows["annotypes"].str.contains("noncoding", case=False, na=False)]
        cadd_genes = []
        if not coding_rows.empty:
            genes = [x.split(';') for x in coding_rows["genes"].unique().tolist()]
            genes = [gene for sublist in genes for gene in sublist]  # flatten
            # unique genes
            genes = list(set(genes))
            likely_genes[snp] = {'genes': genes, 'source': 'coding'}
            continue
        else:
            cadd_genes = cadd_df["genes"].dropna().unique().tolist()
        # Check variants in LD for coding impact
        coding_ld_rows = ld_cadd_df[ld_cadd_df["annotypes"].str.contains("coding", case=False, na=False)]
        coding_ld_rows = coding_ld_rows[~coding_ld_rows["annotypes"].str.contains("noncoding", case=False, na=False)]
        cadd_ld_genes = []
        if not coding_ld_rows.empty:
            genes = [x.split(';') for x in coding_ld_rows["genes"].unique().tolist()]
            genes = [gene for sublist in genes for gene in sublist]  # flatten
            # unique genes
            genes = list(set(genes))
            likely_genes[snp] = {'genes': genes, 'source': 'coding_ld'}
            continue
        else:
            cadd_ld_genes = ld_cadd_df["genes"].dropna().unique().tolist()
        # Next check QTLs
        if not eqtl_df.empty or not sqtl_df.empty:
            eqtl_genes = []
            sqtl_genes = []
            if not eqtl_df.empty:
                eqtl_genes = eqtl_df["gene"].dropna().unique().tolist()
            if not sqtl_df.empty:
                sqtl_genes = sqtl_df["gene"].dropna().unique().tolist()
            qtl_combined = list(set(eqtl_genes + sqtl_genes + cadd_genes))
            likely_genes[snp] = {'genes': qtl_combined, 'source': 'qtl_query'}
            continue
        # Next check LD QTLs
        if not ld_eqtl_df.empty or not ld_sqtl_df.empty:
            ld_eqtl_genes = []
            ld_sqtl_genes = []
            if not ld_eqtl_df.empty:
                ld_eqtl_genes = ld_eqtl_df["gene"].dropna().unique().tolist()
            if not ld_sqtl_df.empty:
                ld_sqtl_genes = ld_sqtl_df["gene"].dropna().unique().tolist()
            ld_qtl_combined = list(set(ld_eqtl_genes + ld_sqtl_genes + cadd_ld_genes))
            likely_genes[snp] = {'genes': ld_qtl_combined, 'source': 'qtl_ld'}
            continue
        # If none of the above, return closest gene (placeholder)
        likely_genes[snp] = {'genes': closest_gene(cadd_df), 'source': 'closest_gene'}
    # Convert to dataframe
    likely_genes_df = pd.DataFrame.from_dict(likely_genes, orient='index')
    likely_genes_df['query'] = likely_genes_df.index
    likely_genes_df = likely_genes_df.reset_index(drop=True)
    return likely_genes_df

# ---------------------------------------------------------
# Convert info dict to DataFrames
# ---------------------------------------------------------
def convert_info_dict_to_dfs(info: dict):
    combined = {}
    # number of evidence types (CADD, eQTL, sQTL, etc.)
    n_types = len(next(iter(info.values())))
    # Prepare empty lists for each evidence type
    for i in range(n_types):
        combined[i] = []
    # Now collect by index
    for snp, entries in info.items():
        for i, edict_list in enumerate(entries):
            combined[i].extend(edict_list)
    # Convert to DataFrames
    info_df = pd.DataFrame(combined[0])
    cadd_df = pd.DataFrame(combined[1])
    eqtl_df = pd.DataFrame(combined[2])
    sqtl_df = pd.DataFrame(combined[3])
    ld_df = pd.DataFrame(combined[4])
    gwas_df = pd.DataFrame(combined[8])
    ld_cadd_df = pd.DataFrame(combined[5])
    ld_eqtl_df = pd.DataFrame(combined[6])
    ld_sqtl_df = pd.DataFrame(combined[7])
    return info_df, cadd_df, eqtl_df, sqtl_df, ld_df, gwas_df, ld_cadd_df, ld_eqtl_df, ld_sqtl_df

# ---------------------------------------------------------
# Create bootstrap datasets of most likely genes
# ---------------------------------------------------------
def create_bootstrap_datasets(most_likely_gene: pd.DataFrame, n_iterations: int = 10):
    bootstrap_datasets = []
    for i in range(n_iterations):
        # Sample from genes list with replacement
        sampled_genes = most_likely_gene['genes'].apply(lambda genes: random.choices(genes, k=1)[0] if genes else None)
        bootstrap_datasets.append(list(set(sampled_genes.tolist())))
    return bootstrap_datasets

# ---------------------------------------------------------
# Gene-set enrichment analysis
# ---------------------------------------------------------
def gene_set_enrichment_analysis(most_likely_gene: pd.DataFrame, n_iterations: int = 10):
    # Create client
    gp = GProfiler(return_dataframe=True)
    # Set bootstrap iterations
    bootstrap_iterations = n_iterations
    # Create bootstrap datasets
    bootstrap_datasets = create_bootstrap_datasets(most_likely_gene, n_iterations=bootstrap_iterations)
    # Run gene-set enrichment for each bootstrap dataset
    enrichment_results = []
    for b in bootstrap_datasets:
        print(f"Running GSEA for bootstrap dataset with {len(b)} genes.", file=sys.stderr, flush=True)
        if not b:
            continue
        res = gp.profile(organism='hsapiens', query=b, sources=['GO:BP', 'KEGG', 'REAC', 'WP'], significance_threshold_method='fdr', no_evidences=False, user_threshold=1)
        enrichment_results.append(res)
    # average results across bootstraps
    avg_enrichment = average_enrichment_results(enrichment_results)
    return avg_enrichment

# ---------------------------------------------------------
# Average enrichment results across bootstraps
# ---------------------------------------------------------
def average_enrichment_results(enrichment_results: List[pd.DataFrame]):
    if not enrichment_results:
        return pd.DataFrame()
    else:
        # Concatenate all results
        concat_df = pd.concat(enrichment_results, ignore_index=True)
        # Convert intersections column to lists (sometimes they come as strings)
        if "intersections" in concat_df.columns:
            concat_df["intersections"] = concat_df["intersections"].apply(lambda x: x if isinstance(x, list) else [])
        else:
            concat_df["intersections"] = [[] for _ in range(len(concat_df))]
        # Group by unique term identifier (e.g., 'term_id')
        group_cols = ['native']
        # Numerical metrics to average
        mean_cols = ['p_value', 'precision', 'recall']
        # Columns where values are constant (take first)
        keep_first_cols = ['term_size', 'query_size', 'source', 'name', 'description']
        # Group and aggregate
        grouped = concat_df.groupby(group_cols)
        aggregated_rows = []
        for keys, sub in grouped:
            row = {}
            # Unpack keys
            for col, val in zip(group_cols, keys):
                row[col] = val
            # Average numerical metrics
            for col in mean_cols:
                if col in sub.columns:
                    row[col] = sub[col].mean()
                else:
                    row[col] = None
            # Take first values for consistent metadata
            for col in keep_first_cols:
                if col in sub.columns:
                    row[col] = sub[col].iloc[0]
                else:
                    row[col] = None
            # Union intersections
            if "intersections" in sub.columns:
                union_set = set()
                for lst in sub["intersections"]:
                    if isinstance(lst, list):
                        union_set.update(lst)
                row["intersections"] = sorted(list(union_set))
            else:
                row["intersections"] = []
            aggregated_rows.append(row)
        # Create final DataFrame
        out_df = pd.DataFrame(aggregated_rows)
    return out_df

# ---------------------------------------------------------
# Semantic similarity analysis -- pygosemsim
# ---------------------------------------------------------
def semantic_pygosemsim(enrichment_df: pd.DataFrame, p_threshold: float = 0.05):
    # Set Path
    go_path = DATA_PATH / "../Annotation/INPUTS_OTHER/20220510_go"
    # Load Go graph
    G = graph.from_resource(go_path)
    # Prepare list of GO terms
    term_list = list(G)
    # Take precalculated lower bounds
    similarity.precalc_lower_bounds(G)
    # Parse enrichment terms -- take only p<threshold
    sb = enrichment_df[(enrichment_df['p_value'] < p_threshold) & (enrichment_df['source'].str.contains('GO'))]
    # sort by p_value
    sb = sb.sort_values(by='p_value', ascending=True)
    # remove duplicates
    sb = sb.drop_duplicates(subset=['native'])
    # get list of GO terms
    go_list = sb['native'].unique().tolist()
    # get a dictionary of native (key) and pvalue (value)
    go_pval_dict = dict(zip(sb['native'], sb['p_value']))
    # Make sure terms are in the graph
    go_list = [term for term in go_list if term in term_list]
    go_pval_dict = {k: v for k, v in go_pval_dict.items() if k in go_list}
    # Calculate pairwise semantic similarity
    dist_mt = []
    for i in range(len(go_list)):
        row = []
        for j in range(len(go_list)):
            if i == j:
                row.append(1.0)
            else:
                sim = similarity.lin(G, go_list[i], go_list[j])
                row.append(sim)
        dist_mt.append(row)
    return dist_mt, go_list

# ---------------------------------------------------------
# CLI entrypoint
# ---------------------------------------------------------
def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Variant annotation query (standalone CLI). Outputs a pickled list of pandas DataFrames to stdout.")
    parser.add_argument("-q", "--query", required=True, help="Variant query (rsID like rs7412 or coordinate like 19:45411941 or chr19 45411941). Can be a file with one variant per line.")
    parser.add_argument("-b", "--build", default="hg38", choices=["hg19", "hg38"], help="Reference genome build of the input query (default: hg38).")
    parser.add_argument("-r", "--random", required=False, help="Random number to create folder for results (Default is AnnotateMe_results.")
    args = parser.parse_args()
    # Place arguments into variables
    query = args.query
    build = args.build
    random_number = args.random
    # query = 'rs7412'
    # build = 'hg38'
    # random = 123456
    # Run query
    info = run_variant_query(query, build)
    # Identify most likely gene
    most_likely_gene = identify_most_likely_gene(info)
    # Combine dictionaries by index across dictionary keys into dataframes (eg info[all_queries][0] + info[all_queries][1] + ...)
    info_df, cadd_df, eqtl_df, sqtl_df, ld_df, gwas_df, ld_cadd_df, ld_eqtl_df, ld_sqtl_df = convert_info_dict_to_dfs(info)
    # Gene-set enrichment analysis
    enrichment_df = gene_set_enrichment_analysis(most_likely_gene, n_iterations=10)
    # Semantic similarity analysis with pygosemsim
    dist_mt, go_list = semantic_pygosemsim(enrichment_df, p_threshold=0.05)
    
if __name__ == "__main__":
    main()