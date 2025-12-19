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
import numpy as np
import os
import time
from collections import defaultdict
import pickle
import random
from functools import reduce
from pygosemsim import download
from pygosemsim import graph
from pygosemsim import similarity
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from sklearn.manifold import MDS
import networkx as nx
from wordcloud import WordCloud
import matplotlib.pyplot as plt
from collections import Counter
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
import psutil

# ---------------------------------------------------------
# Configuration
# ---------------------------------------------------------
# Local data path
DATA_PATH = Path("/Users/nicco/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/snpXplorer/Data")
PATH_SCRIPT = Path("/Users/nicco/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/snpXplorer/Annotation/BIN/")
# Server data path
#DATA_PATH = Path("/Data")
#PATH_SCRIPT = '/Annotation/BIN/'

DB_FILE = DATA_PATH / "databases/Genes/variant_info.db"

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
def fetch_by_rsid(DB_FILE, rsid):
    """
    Query the SQLite database for a single rsID.
    Returns a pandas DataFrame (can have 0, 1, or multiple rows).
    """
    conn = sqlite3.connect(str(DB_FILE))
    query = """
        SELECT *
        FROM variant_info
        WHERE rsID = ?
    """
    info = pd.read_sql_query(query, conn, params=(rsid,)).to_dict(orient="records")[0]
    info['query_type'] = 'rsid'
    try:
        info['AF'] = round(float(info.get('AF')), 2)
    except Exception:
        info['AF'] = None
    return {"query_type": "rsid", "rsid": rsid, "chr_hg38": info['Chromosome'], "pos_hg38": info['Position_hg38'], "pos_hg19": info['Position_hg19'], "ref": info['REF'], "alt": info['ALT'], "af": info['AF']}

# ---------------------------------------------------------
# DB helper: lookup by CHR + POS (hg38)
# ---------------------------------------------------------
def fetch_by_position(DB_FILE, chrom, pos38):
    """
    Query the SQLite database for a variant using chromosome and hg38 position.
    Returns a pandas DataFrame.
    """
    conn = sqlite3.connect(str(DB_FILE))
    chrom = str(chrom).upper().replace("CHR", "")
    query = """
        SELECT *
        FROM variant_info
        WHERE Chromosome = ?
          AND Position_hg38 = ?
    """
    info = pd.read_sql_query(query, conn, params=(str(chrom), int(pos38))).to_dict(orient="records")[0]
    try:
        info['AF'] = round(float(info.get('AF')), 2)
    except Exception:
        info['AF'] = None
    return {"query_type": "coord", "input_chr": chrom, "input_pos_hg38": pos38, "rsid": info['rsID'], "chr_hg38": chrom, "pos_hg38": pos38, "pos_hg19": info['Position_hg19'], "ref": info['REF'], "alt": info['ALT'], "af": info['AF']}

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
def query_eqtls(chrom, pos38, qtl_tissues):
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
        if 'all' not in qtl_tissues:
            df = df[df['tissue'].isin(qtl_tissues)]
        return df.to_dict(orient="records")
    except Exception:
        return []

# ---------------------------------------------------------
# Query sQTLs
# ---------------------------------------------------------
def query_sqtls(chrom, pos38, qtl_tissues):
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
        if 'all' not in qtl_tissues:
            df = df[df['tissue'].isin(qtl_tissues)]
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
def query_ld(chrom, pos38, ld_threshold):
    """
    Query LD partners for given chromosome and position.
    Returns:
      ld_records: list[dict]
      all_partner_positions: list[int]
    """
    try:
        chrom_clean = str(chrom).upper().replace("CHR", "")
        db_path = DATA_PATH / f"databases/LD_db/ld_chr{chrom_clean}.sqlite"
        ld_partners = partners_by_position(db_path, int(pos38), r2_min=ld_threshold)
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
            info = fetch_by_position(DB_FILE, chrom, v)
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
        df_ld = df_ld[df_ld["r2"] >= ld_threshold]
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
def annotate_ld_partners(chrom, ld_rows, pos38, qtl_tissues):
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
        eqtl_rows = query_eqtls(chrom_clean, partner_pos, qtl_tissues)
        for row in eqtl_rows:
            r = dict(row)
            r["ld_partner_pos"] = partner_pos
            r["ld_rsid"] = rsid
            r["ld_r2"] = r2
            r["ld_dist_bp"] = dist
            r['ld_query_pos'] = pos38
            ld_eqtl.append(r)
        # sQTL for partner
        sqtl_rows = query_sqtls(chrom_clean, partner_pos, qtl_tissues)
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
def run_variant_query(q: str, build: str = "hg38", qtl_tissues: List[str] = ["all"], email: str = "", output_folder: str = "", ld_threshold: float = 0.2) -> dict:
    print(f"Running variant query for: {q} (build={build})", file=sys.stderr, flush=True)
    # Parse query
    parsed = parse_variant_query(q)
    if len(parsed) == 0:
        return {}
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
            try:
                info = fetch_by_rsid(DB_FILE, p["rsid"])
            except Exception:
                info = None
            if info is None:
                output[p['rsid']] = 'Error: rsID not found in DB'
                cadd_dic, eqtl_dic, sqtl_dic, ld_dic, gwas_dic = None, None, None, None, None
                continue
            try:
                chr38 = info["chr_hg38"]
                pos38 = info["pos_hg38"]
                # Add CADD annotation
                cadd_dic = query_cadd_score(chr38, pos38)
                # Add eQTL annotation
                eqtl_dic = query_eqtls(chr38, pos38, qtl_tissues)
                # Add sQTL annotation
                sqtl_dic = query_sqtls(chr38, pos38, qtl_tissues)
                # Add LD annotation
                ld_dic, all_vars = query_ld(chr38, pos38, ld_threshold)
                # Add GWAS associations
                gwas_dic = query_gwas_associations(chr38, pos38)
                # Add SV annotation
                sv_list = extract_sv(DATA_PATH, chr38, pos38, build, ['all'])
                # CADD / eQTL / sQTL for all LD partners
                if ld_dic:
                    ld_cadd, ld_eqtl, ld_sqtl = annotate_ld_partners(chr38, ld_dic, pos38, qtl_tissues)
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
                sv_list = None
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
            output[p['rsid']].append(sv_list if sv_list is not None else [])
            continue
        # coordinate path
        chrom = p["chrom"]
        pos = p["pos"]
        if build == "hg38":
            try:
                info = fetch_by_position(DB_FILE, chrom, pos)
            except Exception:
                info = None
            if info is None:
                output[p['query']] = 'Error: coordinate not found in DB'
                cadd_dic, eqtl_dic, sqtl_dic, ld_dic, gwas_dic = None, None, None, None, None
                continue
            try:
                chr38 = info["chr_hg38"]
                pos38 = info["pos_hg38"]
                # Add CADD annotation
                cadd_dic = query_cadd_score(chr38, pos38)
                # Add eQTL annotation
                eqtl_dic = query_eqtls(chr38, pos38, qtl_tissues)
                # Add sQTL annotation
                sqtl_dic = query_sqtls(chr38, pos38, qtl_tissues)
                # Add LD annotation
                ld_dic, all_vars = query_ld(chr38, pos38, ld_threshold)
                # Add GWAS associations
                gwas_dic = query_gwas_associations(chr38, pos38)
                # Add SV annotation
                sv_list = extract_sv(DATA_PATH, chr38, pos38, build, ['all'])
                # CADD / eQTL / sQTL for all LD partners
                if ld_dic:
                    ld_cadd, ld_eqtl, ld_sqtl = annotate_ld_partners(chr38, ld_dic, pos38, qtl_tissues)
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
                sv_list = None
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
            output[p['query']].append(sv_list if sv_list is not None else [])
            continue
        if build == "hg19":
            try:
                lifted = liftover_hg19_to_hg38(chrom, pos)
            except Exception:
                lifted = None
            if lifted is None:
                output[p['query']] = 'Error: liftover failed from hg19 to hg38'
                cadd_dic, eqtl_dic, sqtl_dic, ld_dic, gwas_dic = None, None, None, None, None
                continue
            chr38, pos38 = lifted
            info = fetch_by_position(DB_FILE, chr38, pos38)
            if info is None:
                output[p['query']] = 'Error: coordinate not found in DB'
                continue
            info["liftover"] = {"from_hg19": f"{chrom}:{pos}", "to_hg38": f"{chr38}:{pos38}"}
            try:
                # Add CADD annotation based on hg38 location
                cadd_dic = query_cadd_score(chr38, pos38)
                # Add eQTL annotation based on hg38 location
                eqtl_dic = query_eqtls(chr38, pos38, qtl_tissues)
                # Add sQTL annotation based on hg38 location
                sqtl_dic = query_sqtls(chr38, pos38, qtl_tissues)
                # Add LD annotation based on hg38 location
                ld_dic, all_vars = query_ld(chr38, pos38, ld_threshold)
                # Add GWAS associations
                gwas_dic = query_gwas_associations(chr38, pos38)
                # Add SV annotation
                sv_list = extract_sv(DATA_PATH, chr38, pos38, build, ['all'])
                # CADD / eQTL / sQTL for all LD partners
                if ld_dic:
                    ld_cadd, ld_eqtl, ld_sqtl = annotate_ld_partners(chr38, ld_dic, pos38, qtl_tissues)
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
                sv_list = None
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
            output[p['query']].append(sv_list if sv_list is not None else [])
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
    # remove genes starting with LIN or LOC
    df = df[~df["gene_symbol"].str.startswith(("LINC", "LOC"))]
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
    invalid_queries = []
    # Iterate over snps in info
    for snp in info.keys():
        if 'Error' in info[snp]:
            invalid_queries.append(snp)
            continue
        # Extract dataframes
        info_df, cadd_df, eqtl_df, sqtl_df, ld_df, gwas_df, ld_cadd_df, ld_eqtl_df, ld_sqtl_df, sv_df = convert_info_dict_to_dfs({snp: info[snp]})
        # Get rsid
        rsid = info_df["rsid"].values[0]
        # First check if variant is coding in CADD
        coding_rows = cadd_df[cadd_df["annotypes"].str.contains("coding", case=False, na=False)]
        coding_rows = coding_rows[~coding_rows["annotypes"].str.contains("noncoding", case=False, na=False)]
        cadd_genes = []
        if not coding_rows.empty:
            try:
                genes = [x.split(';') for x in coding_rows["genes"].unique().tolist()]
                genes = [gene for sublist in genes for gene in sublist]  # flatten
            except Exception:
                genes = []
            # unique genes
            genes = list(set(genes))
            if len(genes) >0:
                likely_genes[rsid] = {'genes': genes, 'source': 'coding'}
                continue
            else:
                cadd_genes = cadd_df["genes"].dropna().unique().tolist()
        else:
            cadd_genes = cadd_df["genes"].dropna().unique().tolist()
        # Check variants in LD for coding impact
        coding_ld_rows = pd.DataFrame()
        if not ld_cadd_df.empty:
            coding_ld_rows = ld_cadd_df[ld_cadd_df["annotypes"].str.contains("coding", case=False, na=False)]
            coding_ld_rows = coding_ld_rows[~coding_ld_rows["annotypes"].str.contains("noncoding", case=False, na=False)]
            if not coding_ld_rows.empty:
                try:
                    genes = [x.split(';') for x in coding_ld_rows["genes"].unique().tolist()]
                    genes = [gene for sublist in genes for gene in sublist]  # flatten
                except Exception:
                    genes = []
                # unique genes
                genes = list(set(genes))
                if len(genes) >0:
                    likely_genes[rsid] = {'genes': genes, 'source': 'coding_ld'}
                    continue
                else:
                    cadd_ld_genes = ld_cadd_df["genes"].dropna().unique().tolist()
            else:
                cadd_ld_genes = ld_cadd_df["genes"].dropna().unique().tolist()
        else:
            cadd_ld_genes = []
        # Next check QTLs
        if not eqtl_df.empty or not sqtl_df.empty:
            eqtl_genes = []
            sqtl_genes = []
            if not eqtl_df.empty:
                eqtl_genes = eqtl_df["gene"].dropna().unique().tolist()
            if not sqtl_df.empty:
                sqtl_genes = sqtl_df["gene"].dropna().unique().tolist()
            qtl_combined = list(set(eqtl_genes + sqtl_genes + cadd_genes))
            if len(qtl_combined) >0:
                likely_genes[rsid] = {'genes': qtl_combined, 'source': 'qtl_query'}
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
            if len(ld_qtl_combined) >0:
                likely_genes[rsid] = {'genes': ld_qtl_combined, 'source': 'qtl_ld'}
                continue
        # If none of the above, return closest gene (placeholder)
        likely_genes[rsid] = {'genes': closest_gene(cadd_df), 'source': 'closest_gene'}
    # Convert to dataframe
    likely_genes_df = pd.DataFrame.from_dict(likely_genes, orient='index')
    likely_genes_df['query'] = likely_genes_df.index
    likely_genes_df = likely_genes_df.reset_index(drop=True)
    # Remove invalid queries from info dictionary -- to implement
    filtered = {k: v for k, v in info.items() if k not in invalid_queries}
    return likely_genes_df, filtered, invalid_queries

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
    sv_df = pd.DataFrame(combined[9], columns=['Repeat Class', 'Chromosome', 'Start Position (hg38)', 'End Position (hg38)', 'Length', 'Repeat Name', 'Repeat Family', 'Color'])
    return info_df, cadd_df, eqtl_df, sqtl_df, ld_df, gwas_df, ld_cadd_df, ld_eqtl_df, ld_sqtl_df, sv_df

# ---------------------------------------------------------
# Merge tables
# ---------------------------------------------------------
def merge_info(info_df, cadd_df, eqtl_df, sqtl_df, ld_df, gwas_df, ld_cadd_df, ld_eqtl_df, ld_sqtl_df, most_likely_gene, enrichment_df, sv_df):
    # Return empty DF if info_df is empty
    if info_df.empty:
        return pd.DataFrame()
    # Rename info_df columns
    info_df_renamed = info_df.rename(columns={"query_type": "Query Type", "rsid": "RsID", "chr_hg38": "Chromosome", "pos_hg38": "Position (hg38)", "pos_hg19": "Position (hg19)", "ref": "Reference", "alt": "Alternative", "af": "Alternative Allele Frequency"}).reset_index(drop=True)
    # Add locus ID as chrom:pos:alleles (alleles are alphabetically ordered)
    info_df_renamed['Locus ID'] = info_df_renamed['Chromosome'].astype(str) + ":" + info_df_renamed['Position (hg38)'].astype(str) + ":" + info_df_renamed[['Reference', 'Alternative']].apply(lambda row: ':'.join(sorted([str(row['Reference']), str(row['Alternative'])])), axis=1)
    # Check if cadd_df is empty
    if cadd_df.empty:
        cadd_df_renamed = pd.DataFrame(columns=["Locus ID", "CADD Reference", "CADD Alternative", "CADD Score", "CADD Genes", "CADD Annotation Types", "CADD Consequences", "CADD CpG Max", "CADD GC Max", "CADD H3K27ac Max", "CADD H3K4me3 Max", "CADD DNase Max", "CADD SpliceAI Acc Loss Max", "CADD SpliceAI Don Loss Max", "CADD SIFT Max", "CADD PolyPhen Max", "CADD GERPRS Max", "CADD PhyloP Max"])
    else:
        # Rename cadd_df columns
        cadd_df_renamed = cadd_df.rename(columns={"ref": "CADD Reference", "alt": "CADD Alternative", "phred_max": "CADD Score", "genes": "CADD Genes", "annotypes": "CADD Annotation Types", "consequences": "CADD Consequences", "cpg_max": "CADD CpG Max", "gc_max": "CADD GC Max", "h3k27ac_max": "CADD H3K27ac Max", "h3k4me3_max": "CADD H3K4me3 Max", "dnase_max": "CADD DNase Max", "spliceai_acc_loss_max": "CADD SpliceAI Acc Loss Max", "spliceai_don_loss_max": "CADD SpliceAI Don Loss Max", "siftval_max": "CADD SIFT Max", "polyphenval_max": "CADD PolyPhen Max", "gerprs_max": "CADD GERPRS Max", "priphyloP_max": "CADD PhyloP Max"}).reset_index(drop=True)
        # Add locus ID to cadd_df
        cadd_df_renamed['Locus ID'] = cadd_df_renamed['chrom'].astype(str) + ":" + cadd_df_renamed['pos'].astype(str) + ":" + cadd_df_renamed[['CADD Reference', 'CADD Alternative']].apply(lambda row: ':'.join(sorted([str(row['CADD Reference']), str(row['CADD Alternative'])])), axis=1)
        cadd_df_renamed = cadd_df_renamed.drop(columns=['chrom', 'pos'])
        # Sort by phred score descending
        cadd_df_renamed = cadd_df_renamed.sort_values(by='CADD Score', ascending=False).reset_index(drop=True)
        # Remove duplicate locus IDs, keeping highest CADD score
        cadd_df_renamed = cadd_df_renamed.drop_duplicates(subset=['Locus ID'])
    # Merge info_df and cadd_df on Locus ID
    info_cadd_df = pd.merge(info_df_renamed, cadd_df_renamed, on='Locus ID', how='left')
    # Check if eqtl_df is empty
    if eqtl_df.empty:
        eqtl_grouped = pd.DataFrame(columns=["Locus ID", "eQTL Reference", "eQTL Alternative", "eQTL ensemble", "eQTL Gene", "eQTL Tissue", "eQTL TSS Distance", "eQTL P-value", "eQTL Slope", "eQTL MAF"])
        eqtl_df_renamed = eqtl_grouped
    else:
        # Rename eqtl_df columns
        eqtl_df_renamed = eqtl_df.rename(columns={"ref": "eQTL Reference", "alt": "eQTL Alternative", "ensemble": "eQTL ensemble", "gene": "eQTL Gene", "tissue": "eQTL Tissue", "tss_distance": "eQTL TSS Distance", "pval_nominal": "eQTL P-value", "slope": "eQTL Slope", "maf": "eQTL MAF"}).reset_index(drop=True)
        # Add locus ID to eqtl_df
        eqtl_df_renamed['Locus ID'] = eqtl_df_renamed['chrom'].astype(str) + ":" + eqtl_df_renamed['pos'].astype(str) + ":" + eqtl_df_renamed[['eQTL Reference', 'eQTL Alternative']].apply(lambda row: ':'.join(sorted([str(row['eQTL Reference']), str(row['eQTL Alternative'])])), axis=1)
        # Compress multiple eQTLs per locus into lists
        eqtl_grouped = eqtl_df_renamed.groupby('Locus ID').agg({'eQTL Reference': lambda x: ','.join(x.dropna().astype(str)), 'eQTL Alternative': lambda x: ','.join(x.dropna().astype(str)), 'eQTL ensemble': lambda x: ','.join(x.dropna().astype(str)), 'eQTL Gene': lambda x: ','.join(x.dropna().astype(str)), 'eQTL Tissue': lambda x: ','.join(x.dropna().astype(str)), 'eQTL TSS Distance': lambda x: ','.join(x.dropna().astype(str)), 'eQTL P-value': lambda x: ','.join(x.dropna().astype(str)), 'eQTL Slope': lambda x: ','.join(x.dropna().astype(str)), 'eQTL MAF': lambda x: ','.join(x.dropna().astype(str))}).reset_index()
        # Drop chrom and pos columns from eqtl_df
        eqtl_grouped = eqtl_grouped.drop(columns=['chrom', 'pos'], errors='ignore')
    # Merge info_cadd_df and eqtl_grouped on Locus ID
    info_cadd_eqtl_df = pd.merge(info_cadd_df, eqtl_grouped, on='Locus ID', how='left')
    # Check if sqtl_df is empty
    if sqtl_df.empty:
        sqtl_grouped = pd.DataFrame(columns=["Locus ID", "sQTL Reference", "sQTL Alternative", "sQTL ensemble", "sQTL Gene", "sQTL Tissue", "sQTL TSS Distance", "sQTL P-value", "sQTL Slope", "sQTL MAF"])
        sqtl_df_renamed = sqtl_grouped
    else:
        # Similarly process sqtl_df
        sqtl_df_renamed = sqtl_df.rename(columns={"ref": "sQTL Reference", "alt": "sQTL Alternative", "ensemble": "sQTL ensemble", "gene": "sQTL Gene", "tissue": "sQTL Tissue", "tss_distance": "sQTL TSS Distance", "pval_nominal": "sQTL P-value", "slope": "sQTL Slope", "maf": "sQTL MAF"}).reset_index(drop=True)
        # Add locus ID to sqtl_df
        sqtl_df_renamed['Locus ID'] = sqtl_df_renamed['chrom'].astype(str) + ":" + sqtl_df_renamed['pos'].astype(str) + ":" + sqtl_df_renamed[['sQTL Reference', 'sQTL Alternative']].apply(lambda row: ':'.join(sorted([str(row['sQTL Reference']), str(row['sQTL Alternative'])])), axis=1)
        # Compress multiple sQTL
        sqtl_grouped = sqtl_df_renamed.groupby('Locus ID').agg({'sQTL Reference': lambda x: ','.join(x.dropna().astype(str)), 'sQTL Alternative': lambda x: ','.join(x.dropna().astype(str)), 'sQTL ensemble': lambda x: ','.join(x.dropna().astype(str)), 'sQTL Gene': lambda x: ','.join(x.dropna().astype(str)), 'sQTL Tissue': lambda x: ','.join(x.dropna().astype(str)), 'sQTL TSS Distance': lambda x: ','.join(x.dropna().astype(str)), 'sQTL P-value': lambda x: ','.join(x.dropna().astype(str)), 'sQTL Slope': lambda x: ','.join(x.dropna().astype(str)), 'sQTL MAF': lambda x: ','.join(x.dropna().astype(str))}).reset_index()
        # Drop chrom and pos columns from sqtl_df
        sqtl_grouped = sqtl_grouped.drop(columns=['chrom', 'pos'], errors='ignore')
    # Merge with previous
    merged_df = pd.merge(info_cadd_eqtl_df, sqtl_grouped, on='Locus ID', how='left')
    # Check if most_likely_gene is empty
    if most_likely_gene.empty:
        most_likely_gene_renamed = pd.DataFrame(columns=["RsID", "Most Likely Genes", "Source"])
    else:
        most_likely_gene_renamed = most_likely_gene.rename(columns={"query": "RsID", "genes": "Most Likely Genes", "source": "Source"}).reset_index(drop=True)
        most_likely_gene_renamed['Most Likely Genes'] = most_likely_gene_renamed['Most Likely Genes'].apply(lambda x: ','.join(x) if isinstance(x, list) else x)
    # Merge with previous
    merged_df = pd.merge(merged_df, most_likely_gene_renamed, on='RsID', how='left')
    # Check if ld_cadd_df is empty
    if ld_cadd_df.empty:
        ld_cadd_df_renamed = pd.DataFrame(columns=["Chromosome", "Position (hg38)", "LD CADD Reference", "LD CADD Alternative", "LD CADD Score", "LD CADD Genes", "LD CADD Annotation Types", "LD CADD Consequences", "LD CADD CpG Max", "LD CADD GC Max", "LD CADD H3K27ac Max", "LD CADD H3K4me3 Max", "LD CADD DNase Max", "LD CADD SpliceAI Acc Loss Max", "LD CADD SpliceAI Don Loss Max", "LD CADD SIFT Max", "LD CADD PolyPhen Max", "LD CADD GERPRS Max", "LD CADD PhyloP Max"])
    else:
        # Rename ld_cadd_df columns
        ld_cadd_df_renamed = ld_cadd_df.rename(columns={"ref": "LD CADD Reference", "alt": "LD CADD Alternative", "phred_max": "LD CADD Score", "genes": "LD CADD Genes", "annotypes": "LD CADD Annotation Types", "consequences": "LD CADD Consequences", "cpg_max": "LD CADD CpG Max", "gc_max": "LD CADD GC Max", "h3k27ac_max": "LD CADD H3K27ac Max", "h3k4me3_max": "LD CADD H3K4me3 Max", "dnase_max": "LD CADD DNase Max", "spliceai_acc_loss_max": "LD CADD SpliceAI Acc Loss Max", "spliceai_don_loss_max": "LD CADD SpliceAI Don Loss Max", "siftval_max": "LD CADD SIFT Max", "polyphenval_max": "LD CADD PolyPhen Max", "gerprs_max": "LD CADD GERPRS Max", "priphyloP_max": "LD CADD PhyloP Max", "chrom": "Chromosome", "pos": "Position (hg38)"}).reset_index(drop=True)
    # Check if ld_eqtl_df is empty
    if ld_eqtl_df.empty:
        ld_eqtl_df_renamed = pd.DataFrame(columns=["Chromosome", "Position (hg38)", "LD eQTL Reference", "LD eQTL Alternative", "LD eQTL ensemble", "LD eQTL Gene", "LD eQTL Tissue", "LD eQTL TSS Distance", "LD eQTL P-value", "LD eQTL Slope", "LD eQTL MAF"])
    else:
        # Rename ld_eqtl_df columns
        ld_eqtl_df_renamed = ld_eqtl_df.rename(columns={"ref": "LD eQTL Reference", "alt": "LD eQTL Alternative", "ensemble": "LD eQTL ensemble", "gene": "LD eQTL Gene", "tissue": "LD eQTL Tissue", "tss_distance": "LD eQTL TSS Distance", "pval_nominal": "LD eQTL P-value", "slope": "LD eQTL Slope", "maf": "LD eQTL MAF", "chrom": "Chromosome", "pos": "Position (hg38)"}).reset_index(drop=True)
    # Check if ld_sqtl_df is empty
    if ld_sqtl_df.empty:
        ld_sqtl_df_renamed = pd.DataFrame(columns=["Chromosome", "Position (hg38)", "LD sQTL Reference", "LD sQTL Alternative", "LD sQTL ensemble", "LD sQTL Gene", "LD sQTL Tissue", "LD sQTL TSS Distance", "LD sQTL P-value", "LD sQTL Slope", "LD sQTL MAF"])
    else:    
        # Rename ld_sqtl_df columns
        ld_sqtl_df_renamed = ld_sqtl_df.rename(columns={"ref": "LD sQTL Reference", "alt": "LD sQTL Alternative", "ensemble": "LD sQTL ensemble", "gene": "LD sQTL Gene", "tissue": "LD sQTL Tissue", "tss_distance": "LD sQTL TSS Distance", "pval_nominal": "LD sQTL P-value", "slope": "LD sQTL Slope", "maf": "LD sQTL MAF", "chrom": "Chromosome", "pos": "Position (hg38)"}).reset_index(drop=True)
    # Check if ld_df is empty
    if ld_df.empty:
        ld_df_renamed = pd.DataFrame(columns=["Query position", "LD Partner Position (hg38)", "LD Partner RsID", "LD R2", "LD Distance (bp)", "Chromosome", "Position (hg38)"])
    else:
        # Rename ld_df columns
        ld_df_renamed = ld_df.rename(columns={"query_pos": "Query position", "partner_pos": "LD Partner Position (hg38)", "rsid": "LD Partner RsID", "r2": "LD R2", "dist_bp": "LD Distance (bp)", "chr": "Chromosome", "pos": "Position (hg38)"}).reset_index(drop=True)
    # Check if gwas_df is empty
    if gwas_df.empty:
        gwas_df_renamed = pd.DataFrame(columns=["GWAS Trait", "Effect Allele", "Non-effect Allele", "GWAS Beta", "GWAS SE", "GWAS P-value", "GWAS ID", "Sample size", "RsID"])
    else:
        # Rename gwas_df columns
        gwas_df_renamed = gwas_df.rename(columns={"trait": "GWAS Trait", "ea": "Effect Allele", "nea": "Non-effect Allele", "beta": "GWAS Beta", "se": "GWAS SE", "pval": "GWAS P-value", "gwas_id": "GWAS ID", "n": "Sample size", "rsid": "RsID"}).reset_index(drop=True)
    # Check if enrichment_df is empty
    if enrichment_df.empty:
        enrichment_df_renamed = pd.DataFrame(columns=["Enrichment Term ID", "Enrichment Term Name", "Enrichment Source", "Enrichment P-value", "Enrichment Term Size", "Enrichment Query Size", "Enrichment Precision", "Enrichment Recall", "Enrichment Term Description", "Enrichment intersections"])
    else:
        enrichment_df_renamed = enrichment_df.rename(columns={"native": "Enrichment Term ID", "name": "Enrichment Term Name", "source": "Enrichment Source", "p_value": "Enrichment P-value", "term_size": "Enrichment Term Size", "query_size": "Enrichment Query Size", "precision": "Enrichment Precision", "recall": "Enrichment Recall", "description": "Enrichment Term Description", "intersections": "Enrichment intersections"}).reset_index(drop=True)
        enrichment_df_renamed['Enrichment intersections'] = enrichment_df_renamed['Enrichment intersections'].apply(lambda x: ','.join(x) if isinstance(x, list) else x)
    # Check sv_df
    if sv_df.empty:
        sv_df_renamed = pd.DataFrame(columns=['Repeat Class', 'Chromosome', 'Start Position (hg38)', 'End Position (hg38)', 'Length', 'Repeat Name', 'Repeat Family', 'Color'])
    else:
        sv_df = sv_df.drop(columns=['Color'])
    return cadd_df_renamed, eqtl_df_renamed, sqtl_df_renamed, merged_df, ld_cadd_df_renamed, ld_eqtl_df_renamed, ld_sqtl_df_renamed, ld_df_renamed, gwas_df_renamed, enrichment_df_renamed, sv_df

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
def gene_set_enrichment_analysis(most_likely_gene: pd.DataFrame, n_iterations: int = 10, gsea_sets = ['GO:BP', 'KEGG', 'REAC', 'WP']):
    try:
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
            res = gp.profile(organism='hsapiens', query=b, sources=gsea_sets, significance_threshold_method='fdr', no_evidences=False, user_threshold=1)
            enrichment_results.append(res)
        # average results across bootstraps
        avg_enrichment = average_enrichment_results(enrichment_results)
        # sort by p_value
        avg_enrichment = avg_enrichment.sort_values(by='p_value', ascending=True).reset_index(drop=True)
    except Exception as e:
        print(f"Error during gene-set enrichment analysis: {e}", file=sys.stderr, flush=True)
        avg_enrichment = pd.DataFrame()
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
def semantic_pygosemsim(enrichment_df: pd.DataFrame, p_threshold: float = 0.05, output_folder: str = ""):
    try:
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
        # Convert to matrix with GO terms as row and column names
        dist_mt = pd.DataFrame(dist_mt, index=go_list, columns=go_list)
        # Save distance matrix
        output_folder_gset = f"{output_folder}/gene_set_enrichment"
        os.makedirs(output_folder_gset, exist_ok=True)
        dist_mt.to_csv(f"{output_folder_gset}/semantic_similarity_matrix.tsv", sep="\t", index=True, header=True)
    except Exception as e:
        print(f"Error during semantic similarity analysis: {e}", file=sys.stderr, flush=True)
        dist_mt = pd.DataFrame()
        go_list = []
    return dist_mt, go_list

# ---------------------------------------------------------
# Function count frequency of all GO terms words
# ---------------------------------------------------------
def count_go_term_words(percentace=0.02):
    # Load GO terms
    all_go_bp = pd.read_csv(DATA_PATH / "../Annotation/BIN/go_terms_BP_all.txt", sep="\t", header=None)
    # Put all in a list
    text_all = ' '.join(all_go_bp[0].dropna().tolist()).split()
    # Count frequency of each word
    word_counts = Counter(text_all)
    # Take top X% most common words
    most_common_count = int(len(word_counts) * percentace)
    word_counts = dict(word_counts.most_common(most_common_count))
    return word_counts

# ---------------------------------------------------------
# Function to guide wordclouds after clustering
# ---------------------------------------------------------
def guide_wordclouds(output_folder, enrichment_df):
    try:
        # Calculate word frequencies across all GO terms
        word_counts = count_go_term_words(percentace=0.02)
        # List files in output folder
        files = [x for x in os.listdir(f"{output_folder}/gene_set_enrichment/") if x.startswith('clustering') and x.endswith('.tsv')]
        # Container for clustered GO terms
        clustered_go_terms_list = []
        # Iterate over files
        for file in files:
            # Load clustered GO terms
            clustered_go_terms = pd.read_csv(f"{output_folder}/gene_set_enrichment/{file}", sep="\t")
            # Get max cluster ID
            max_cluster_id = clustered_go_terms['cluster'].max()
            # Create folder for wordclouds
            os.makedirs(f"{output_folder}/gene_set_enrichment/wordclouds_{max_cluster_id}_clusters", exist_ok=True)
            # Draw wordclouds
            clustered_go_terms = draw_wordcloud(clustered_go_terms, enrichment_df, word_counts, f"{output_folder}/gene_set_enrichment/wordclouds_{max_cluster_id}_clusters")
            # Save updated clustered_go_terms with GO term names
            clustered_go_terms.to_csv(f"{output_folder}/gene_set_enrichment/{file}", sep="\t", index=False)
            clustered_go_terms_list.append(clustered_go_terms)
    except Exception as e:
        print(f"Error during wordcloud generation: {e}", file=sys.stderr, flush=True)
        clustered_go_terms_list = pd.DataFrame()
    return clustered_go_terms_list

# ---------------------------------------------------------
# Function to draw wordclouds
# ---------------------------------------------------------
def draw_wordcloud(clustered_go_terms, enrichment_df, word_counts, outpath):
    # Merge to get GO term names
    clustered_go_terms = pd.merge(clustered_go_terms, enrichment_df[['native', 'name']], left_on='term', right_on='native', how='left')
    # Iterate over clusters
    for cluster_id, group in clustered_go_terms.groupby('cluster'):
        # Combine all GO term names in the cluster
        text = ' '.join(group['name'].dropna().tolist()).split()
        # Exclude common words
        filtered_text = [word for word in text if word not in word_counts]
        # Recreate text
        text = ' '.join(filtered_text)
        # Generate word cloud
        wc = WordCloud(width=800, height=400, background_color="white", max_words=20).generate(text)
        plt.figure(figsize=(10, 5))
        plt.imshow(wc, interpolation="bilinear")
        plt.axis("off")
        plt.tight_layout()
        # Save word cloud
        plt.savefig(f"{outpath}/wordcloud_cluster_{cluster_id}.png")
        plt.close()
    return clustered_go_terms

# ---------------------------------------------------------
# Function to guide weights for pathway-PRS
# ---------------------------------------------------------
def prepare_pathway_prs_weights(clustered_go_terms_list, merged_df, output_folder, enrichment_df):
    # Folder to save weights
    os.makedirs(f"{output_folder}/pathway_prs_weights", exist_ok=True)
    # Iterate over clustered_go_terms_list
    for clustered_go_terms in clustered_go_terms_list:
        # Derive weights
        weights_df = derive_pathway_prs_weights(clustered_go_terms, merged_df, enrichment_df)
        # Save weights
        max_cluster_id = clustered_go_terms['cluster'].max()
        weights_df.to_csv(f"{output_folder}/pathway_prs_weights/pathway_prs_weights_{max_cluster_id}_clusters.tsv", sep="\t", index=True, header=True)

# ---------------------------------------------------------
# Function to derive weights for pathway-PRS
# ---------------------------------------------------------
def derive_pathway_prs_weights(clustered_go_terms, merged_df, enrichment_df):
    # Merge with enrichment_df to get genes in each GO term
    merged_enrichment = pd.merge(clustered_go_terms, enrichment_df[['Enrichment Term ID', 'Enrichment intersections']], left_on='term', right_on='Enrichment Term ID', how='left')
    # Container for weights
    weights = {}
    # Iterate over snps
    for idx, row in merged_df.iterrows():
        rsid = row['RsID']
        gene_list = row['Most Likely Genes'].split(',') if pd.notna(row['Most Likely Genes']) else []
        # Iterate over genes and grep all rows from merged_enrichment where the gene is in Enrichment intersections
        total_counts_per_snp = 0
        dict_per_snp = {}
        for gene in gene_list:
            gene_rows = merged_enrichment[merged_enrichment['Enrichment intersections'].str.contains(gene, na=False)]
            # Count number of hits per cluster
            cluster_counts = gene_rows['cluster'].value_counts().to_dict()
            # Derive weights by dividing cluster counts by the total cluster counts
            total_counts_per_snp += sum(cluster_counts.values())
            for cluster_id, count in cluster_counts.items():
                dict_per_snp[cluster_id] = dict_per_snp.get(cluster_id, 0) + count
        # Normalize weights
        if total_counts_per_snp > 0:
            for cluster_id in dict_per_snp:
                dict_per_snp[cluster_id] = dict_per_snp[cluster_id] / total_counts_per_snp
        # Otherwise set weights to zero
        else:
            for cluster_id in merged_enrichment['cluster'].unique():
                dict_per_snp[cluster_id] = 0.0
        # Save weights per rsid
        weights[rsid] = dict_per_snp
    # Convert to DataFrame
    weights_df = pd.DataFrame.from_dict(weights, orient='index').fillna(0)
    return weights_df

# ---------------------------------------------------------
# Function to send email at start of analysis
# ---------------------------------------------------------
def email_at_start(email: str, output_folder: str, query: str, analysis_type: str, build: str, qtl_tissues: List[str], gsea_sets: List[str], random_number: str = None):
    # First read configuration file values for the emails
    cfg = pd.read_table(f"{DATA_PATH}/../Annotation/config_email.txt")
    sender = str(cfg['username'].values[0])
    port = int(cfg['port'].values[0])
    psw = str(cfg['psw'].values[0])
    host = str(cfg['host'].values[0])
    cc_add = 'snpxplorer@gmail.com'
    # Send email to myself to notify a new request has been received
    message_email = f"Dear user, \nsnpXplorer received an annotation request from you. \n You receive thie email to confirm that your request is being processed. A typical job takes about 30 minutes to complete, however, due to the high number of requests, jobs may be delayed. \n\n The following settings were requested: \n input --> {query} \n analysis_type --> {analysis_type} \n analysis mode --> {gsea_sets} \n interest_tissue --> {qtl_tissues} \n ref_version --> {build} \n run_ID --> {random_number} \n\n Thanks for using snpXplorer! \n\n snpXplorer Team"
    msg = MIMEMultipart()
    msg["From"] = sender
    msg["To"] = ", ".join(email) if isinstance(email, list) else email
    msg["Cc"] = ", ".join(cc_add) if isinstance(cc_add, list) else cc_add
    msg["Subject"] = 'snpXplorer request'
    # Add body
    msg.attach(MIMEText(message_email, "plain"))
    # Combine "to" and "cc" for sending
    all_emails = []
    if isinstance(email, list):
        all_emails.extend(email)
    else:
        all_emails.append(email)
    if cc_add:
        if isinstance(cc_add, list):
            all_emails.extend(cc_add)
        else:
            all_emails.append(cc_add)
    # Send email using SSL
    with smtplib.SMTP_SSL(host, port) as server:
        server.login(sender, psw)
        server.sendmail(sender, all_emails, msg.as_string())
    print("Email sent successfully!")

# ---------------------------------------------------------
# Function to send email at end of analysis
# ---------------------------------------------------------
def email_at_end(email: str, random_number: str = None):    
    # First read configuration file values for the emails
    cfg = pd.read_table(f"{DATA_PATH}/../Annotation/config_email.txt")
    sender = str(cfg['username'].values[0])
    port = int(cfg['port'].values[0])
    psw = str(cfg['psw'].values[0])
    host = str(cfg['host'].values[0])
    cc_add = 'snpxplorer@gmail.com'
    # Send email to myself to notify a new request has been received
    message_email = f"Dear user, \n\n Thanks so much for using snpXplorer and its annotation pipeline. \n We hope you find the tool useful. \n\n We now implemented a new way to download your results directly from the web-server. Please open https://snpxplorer.net/download/ and follow the instructions to get your results. \n\nYour Run ID --> {random_number}\n\nBest wishes, \n snpXplorer Team."
    msg = MIMEMultipart()
    msg["From"] = sender
    msg["To"] = ", ".join(email) if isinstance(email, list) else email
    msg["Cc"] = ", ".join(cc_add) if isinstance(cc_add, list) else cc_add
    msg["Subject"] = 'snpXplorer request'
    # Add body
    msg.attach(MIMEText(message_email, "plain"))
    # Combine "to" and "cc" for sending
    all_emails = []
    if isinstance(email, list):
        all_emails.extend(email)
    else:
        all_emails.append(email)
    if cc_add:
        if isinstance(cc_add, list):
            all_emails.extend(cc_add)
        else:
            all_emails.append(cc_add)
    # Send email using SSL
    with smtplib.SMTP_SSL(host, port) as server:
        server.login(sender, psw)
        server.sendmail(sender, all_emails, msg.as_string())
    print("Email sent successfully!")

# ---------------------------------------------------------
# Function to send email for errors
# ---------------------------------------------------------
def send_email_error(email, q, output_folder):
    # First read configuration file values for the emails
    cfg = pd.read_table(f"{DATA_PATH}/../Annotation/config_email.txt")
    sender = str(cfg['username'].values[0])
    port = int(cfg['port'].values[0])
    psw = str(cfg['psw'].values[0])
    host = str(cfg['host'].values[0])
    cc_add = 'snpxplorer@gmail.com'
    # Send email to myself to notify a new request has been received
    message_email = f"Dear user, \nsnpXplorer encountered an error while processing your annotation request. \n Please check your input query: {q} \n\n If the problem persists, please contact us at {sender}. \n\n Best wishes, \n snpXplorer Team"
    msg = MIMEMultipart()
    msg["From"] = sender
    msg["To"] = ", ".join(email) if isinstance(email, list) else email
    msg["Cc"] = ", ".join(cc_add) if isinstance(cc_add, list) else cc_add
    msg["Subject"] = 'snpXplorer request - ERROR'
    # Add body
    msg.attach(MIMEText(message_email, "plain"))
    # Combine "to" and "cc" for sending
    all_emails = []
    if isinstance(email, list):
        all_emails.extend(email)
    else:
        all_emails.append(email)
    if cc_add:
        if isinstance(cc_add, list):
            all_emails.extend(cc_add)
        else:
            all_emails.append(cc_add)
    # Send email using SSL
    with smtplib.SMTP_SSL(host, port) as server:
        server.login(sender, psw)
        server.sendmail(sender, all_emails, msg.as_string())
    # Move query file to error folder
    cmd = f"mv {q} {output_folder}/"
    os.system(cmd)
    # Compress output folder
    cmd = f"cd {os.path.dirname(output_folder)} && zip -r {os.path.basename(output_folder)}.zip {os.path.basename(output_folder)}"
    os.system(cmd)
    # Remove output folder
    cmd = f"rm -rf {output_folder}"
    os.system(cmd)
    print("Error email sent successfully!")

# ---------------------------------------------------------
# Check resources before starting the process
# ---------------------------------------------------------
def wait_for_memory(analysis_type: str, poll_interval: int = 60):
    """
    Block until at least `min_free_gb` GB of RAM are available.
    """
    MIN_FREE_GB = 2 if analysis_type == "annotation" else 4
    while True:
        mem = psutil.virtual_memory()
        free_gb = mem.available / (1024 ** 3)
        if free_gb >= MIN_FREE_GB:
            print(f"[MEM] {free_gb:.2f} GB free – starting job.", file=sys.stderr, flush=True)
            return
        print(f"[MEM] Only {free_gb:.2f} GB free (< {MIN_FREE_GB} GB). Waiting...",
              file=sys.stderr, flush=True)
        time.sleep(poll_interval)

# ---------------------------------------------------------
# Function to query SVs
# ---------------------------------------------------------
def extract_sv(DATA_PATH, chrom, start_pos, refGen, svtypes):
    try:
        # define end_pos
        end_pos = start_pos + 1000
        start_pos = start_pos - 1000
        # set reference prefix
        refPrefix = 'hg19' if refGen == 'GRCh37' else 'hg38'
        # use tabix to find genes -- enlarge window by 50kb up and down
        cmd = 'tabix %s/databases/Structural_variants/harmonized_svs_%s.bed.gz chr%s:%s-%s' %(str(DATA_PATH).replace(' ', '\ '), refPrefix, str(chrom), str(start_pos), str(end_pos))
        svs = [x.rstrip().split('\t') for x in os.popen(cmd)]
        # select based on the input selected
        if 'all' in svtypes:
            pass
        else:
            svs = [x for x in svs if x[0] in svtypes]
    except Exception as e:
        svs = []
    return svs

# ---------------------------------------------------------
# Helper to append to disk
# ---------------------------------------------------------
def append_df(path: str, df: pd.DataFrame):
    """
    Append df to a TSV file, creating it with header if it doesn't exist.
    """
    if df is None or df.empty:
        return
    write_header = not os.path.exists(path)
    df.to_csv(path, sep="\t", index=False, mode="a", header=write_header, na_rep="NA")

# ---------------------------------------------------------
# Function to clean redundant files
# ---------------------------------------------------------
def cleanRedundantFiles(output_folder):
    # Container
    pnginfo = {}
    # List files .png and their sizes
    for p in Path(f"{output_folder}/gene_set_enrichment/").glob("*.png"):
        size = p.stat().st_size
        pnginfo[p] = size
    # Group files by size
    by_size = defaultdict(list)
    for path, size in pnginfo.items():
        by_size[size].append(path)
    # Helper to extract the number
    def extract_max_n(path: Path) -> int:
        m = re.search(r'clustering_max_(\d+)', path.name)
        return int(m.group(1)) if m else float("inf")
    # Select files to keep
    files_to_keep = []
    for size, paths in by_size.items():
        if len(paths) == 1:
            files_to_keep.append(paths[0])
        else:
            # keep the one with the smallest clustering_max_N
            keep = min(paths, key=extract_max_n)
            files_to_keep.append(keep)
    # Remove the other files
    for path in pnginfo:
        if path not in files_to_keep:
            tspath = Path(str(path).replace('dendrogram_clusters.png', 'clusters.tsv'))
            path.unlink()
            tspath.unlink()    

# ---------------------------------------------------------
# CLI entrypoint
# ---------------------------------------------------------
def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Variant annotation query (standalone CLI). Outputs a pickled list of pandas DataFrames to stdout.")
    parser.add_argument("-q", "--query", required=True, help="Variant query (rsID like rs7412 or coordinate like 19:45411941 or chr19 45411941). Can be a file with one variant per line.")
    parser.add_argument("-b", "--build", default="hg38", choices=["hg19", "hg38", 'grch38', 'grch37'], help="Reference genome build of the input query (default: hg38).")
    parser.add_argument("-o", "--output", default=".", help="Output directory for results (Default is current working directory).")
    parser.add_argument("-r", "--random", required=False, help="Random number to create folder for results (Default is AnnotateMe_results).")
    parser.add_argument("-t", "--type", default="annotation", choices=["annotation", "enrichment"], help="Type of analysis to perform (default: annotation).")
    parser.add_argument("-ts", "--qtl-tissues", default="all", help="Comma-separated list of tissues to consider for eQTL/sQTL annotation (default: all).")
    parser.add_argument("-gs", "--gsea-sets", default="GO:BP", help="Comma-separated list of gene sets to use for gene-set enrichment analysis (default: GO:BP). Options: GO:BP, KEGG, REAC, WP.")
    parser.add_argument("-e", "--email", required=True, help="Email address to send results to.")
    args = parser.parse_args()
    # Place arguments into variables
    query = args.query
    build = args.build.lower()
    output_folder = args.output
    random_number = args.random
    analysis_type = args.type
    qtl_tissues = args.qtl_tissues.split(",") if args.qtl_tissues != "all" else ["all"]
    gsea_sets = args.gsea_sets.split(",")
    email = args.email
    # query = '/Users/nicco/Downloads/for_snpxplorer/snpXplorer_input_86821.txt'
    # build = 'hg38'
    # random_number = 86821
    # output_folder = '/Users/nicco/Downloads'
    # analysis_type = 'enrichment'
    # qtl_tissues = ['Whole_Blood']
    # gsea_sets = ['GO:BP', 'KEGG', 'REAC', 'Wiki']
    # email = 'tesinicco@gmail.com'
    
    # Fix reference built
    if build == 'grch38':
        build = 'hg38'
    elif build == 'grch37':
        build = 'hg19'
    
    # Compile output folder name
    if random_number:
        output_folder = os.path.join(output_folder, f"snpXplorer_annotation_{random_number}")
    else:
        output_folder = os.path.join(output_folder, "snpXplorer_annotation")
    
    # Check if output folder exists and in case stop execution
    if not os.path.exists(output_folder):
        # Create output folder
        os.makedirs(output_folder)
        # Add documentation to the output folder
        cmd = f"cp {PATH_SCRIPT}/snpXplorer_output_description.pdf {output_folder}/"
        os.system(cmd)
    else:
        print(f"Output folder {output_folder} already exists. Please remove it or provide a different random number.", file=sys.stderr, flush=True)
        sys.exit(1)
    
    # Email notification
    email_at_start(email, output_folder, query, analysis_type, build, qtl_tissues, gsea_sets, random_number)
    
    # Decide whether to run the process or not
    wait_for_memory(analysis_type)

    # Initialize output files
    most_likely_gene_list = []
    invalid_queries = []
    variant_combined_path = f"{output_folder}/variant_annotation_combined.tsv"
    cadd_path = f"{output_folder}/target_annotations/cadd_annotation_target.tsv"
    eqtl_path = f"{output_folder}/target_annotations/eqtl_annotation_target.tsv"
    sqtl_path = f"{output_folder}/target_annotations/sqtl_annotation_target.tsv"
    gwas_path = f"{output_folder}/target_annotations/gwas_annotation_target.tsv"
    ld_path = f"{output_folder}/ld_partner_annotations/ld_annotation_target.tsv"
    ld_cadd_path = f"{output_folder}/ld_partner_annotations/cadd_annotation_LD_partners.tsv"
    ld_eqtl_path = f"{output_folder}/ld_partner_annotations/eqtl_annotation_LD_partners.tsv"
    ld_sqtl_path = f"{output_folder}/ld_partner_annotations/sqtl_annotation_LD_partners.tsv"
    sv_path = f"{output_folder}/target_annotations/sv_annotation_target.tsv"
    invalid_queries_path = f"{output_folder}/invalid_queries.txt"
    os.makedirs(f"{output_folder}/target_annotations", exist_ok=True)
    os.makedirs(f"{output_folder}/ld_partner_annotations", exist_ok=True)

    # Buffering settings to write to disk every N queries
    FLUSH_EVERY = 30
    processed_valid = 0
    buf_variant = []
    buf_cadd = []
    buf_eqtl = []
    buf_sqtl = []
    buf_gwas = []
    buf_ld = []
    buf_ld_cadd = []
    buf_ld_eqtl = []
    buf_ld_sqtl = []
    buf_sv = []
    
    # If query is a file: process line-by-line
    if Path(query).is_file():
        with open(query) as f:
            queries = [l.strip() for l in f if l.strip()]
            queries = list(set(queries))  # Remove duplicates
    else:
        queries = [query]

    # Set ld threshold
    ld_threshold = 0.4
    
    # Iterate over elements in queries
    for q_line in queries:
        # Run query
        info_single = run_variant_query(q_line, build, qtl_tissues, email, output_folder, ld_threshold)
        # Identify most likely gene
        most_likely_gene_single, info_single_filtered, invalid_queries_single = identify_most_likely_gene(info_single)
        if not most_likely_gene_single.empty:
            # Append to list
            most_likely_gene_list.append(most_likely_gene_single)
            # Combine dictionaries by index across dictionary keys into dataframes (eg info[all_queries][0] + info[all_queries][1] + ...)
            info_df, cadd_df, eqtl_df, sqtl_df, ld_df, gwas_df, ld_cadd_df, ld_eqtl_df, ld_sqtl_df, sv_df = convert_info_dict_to_dfs(info_single_filtered)
            # Merge tables
            cadd_df_m, eqtl_df_m, sqtl_df_m, merged_df_single, ld_cadd_df_m, ld_eqtl_df_m, ld_sqtl_df_m, ld_df_m, gwas_df_m, _, sv_df_m = merge_info(info_df, cadd_df, eqtl_df, sqtl_df, ld_df, gwas_df, ld_cadd_df, ld_eqtl_df, ld_sqtl_df, most_likely_gene_single, pd.DataFrame(), sv_df)
            # Store in buffers
            buf_variant.append(merged_df_single)
            buf_cadd.append(cadd_df_m)
            buf_eqtl.append(eqtl_df_m)
            buf_sqtl.append(sqtl_df_m)
            buf_gwas.append(gwas_df_m)
            buf_ld.append(ld_df_m)
            buf_ld_cadd.append(ld_cadd_df_m)
            buf_ld_eqtl.append(ld_eqtl_df_m)
            buf_ld_sqtl.append(ld_sqtl_df_m)
            buf_sv.append(sv_df_m)
            processed_valid += 1
            # Append partial results to disk
            if processed_valid % FLUSH_EVERY == 0:
                if buf_variant:
                    append_df(variant_combined_path, pd.concat(buf_variant, ignore_index=True))
                    buf_variant = []
                if buf_cadd:
                    append_df(cadd_path, pd.concat(buf_cadd, ignore_index=True))
                    buf_cadd = []
                if buf_eqtl:
                    append_df(eqtl_path, pd.concat(buf_eqtl, ignore_index=True))
                    buf_eqtl = []
                if buf_sqtl:
                    append_df(sqtl_path, pd.concat(buf_sqtl, ignore_index=True))
                    buf_sqtl = []
                if buf_gwas:
                    append_df(gwas_path, pd.concat(buf_gwas, ignore_index=True))
                    buf_gwas = []
                if buf_ld:
                    append_df(ld_path, pd.concat(buf_ld, ignore_index=True))
                    buf_ld = []
                if buf_ld_cadd:
                    append_df(ld_cadd_path, pd.concat(buf_ld_cadd, ignore_index=True))
                    buf_ld_cadd = []
                if buf_ld_eqtl:
                    append_df(ld_eqtl_path, pd.concat(buf_ld_eqtl, ignore_index=True))
                    buf_ld_eqtl = []
                if buf_ld_sqtl:
                    append_df(ld_sqtl_path, pd.concat(buf_ld_sqtl, ignore_index=True))
                    buf_ld_sqtl = []
                if buf_sv:
                    append_df(sv_path, pd.concat(buf_sv, ignore_index=True))
                    buf_sv = []
        else:
            invalid_queries.extend(invalid_queries_single)

    # Final flush to disk
    if buf_variant:
        append_df(variant_combined_path, pd.concat(buf_variant, ignore_index=True))
    if buf_cadd:
        append_df(cadd_path, pd.concat(buf_cadd, ignore_index=True))
    if buf_eqtl:
        append_df(eqtl_path, pd.concat(buf_eqtl, ignore_index=True))
    if buf_sqtl:    
        append_df(sqtl_path, pd.concat(buf_sqtl, ignore_index=True))
    if buf_gwas:
        append_df(gwas_path, pd.concat(buf_gwas, ignore_index=True))
    if buf_ld:
        append_df(ld_path, pd.concat(buf_ld, ignore_index=True))
    if buf_ld_cadd:
        append_df(ld_cadd_path, pd.concat(buf_ld_cadd, ignore_index=True))
    if buf_ld_eqtl:
        append_df(ld_eqtl_path, pd.concat(buf_ld_eqtl, ignore_index=True))
    if buf_ld_sqtl:
        append_df(ld_sqtl_path, pd.concat(buf_ld_sqtl, ignore_index=True))
    if buf_sv:
        append_df(sv_path, pd.concat(buf_sv, ignore_index=True))

    # Save invalid queries
    if invalid_queries:
        with open(invalid_queries_path, "w") as f:
            for iq in invalid_queries:
                f.write(f"{iq}\n")
    
    # If all queries were invalid, send error email and exit
    if len(invalid_queries) == len(queries):
        send_email_error(email, query, output_folder)
        sys.exit(1)
    
    # At this point, combine all most likely genes into a single DataFrame
    if most_likely_gene_list:
        most_likely_gene = pd.concat(most_likely_gene_list, ignore_index=True)
    else:
        most_likely_gene = pd.DataFrame()

    # Gene-set enrichment analysis
    if analysis_type == "enrichment" and not most_likely_gene.empty:
        enrichment_df = gene_set_enrichment_analysis(most_likely_gene, n_iterations=200, gsea_sets=gsea_sets)
        dist_mt, go_list = semantic_pygosemsim(enrichment_df, p_threshold=0.05, output_folder=output_folder)

        clustered_go_terms_list = []
        if not dist_mt.empty:
            cmd = f"Rscript {PATH_SCRIPT}/standalone_annotation.R -d {output_folder}/gene_set_enrichment/semantic_similarity_matrix.tsv -o {output_folder}/gene_set_enrichment"
            os.system(cmd)
            clustered_go_terms_list = guide_wordclouds(output_folder, enrichment_df)
            cleanRedundantFiles(output_folder)
            
        # Rename columns of enrichment_df for consistency
        enrichment_df = enrichment_df.rename(columns={"native": "Enrichment Term ID", "name": "Enrichment Term Name", "source": "Enrichment Source", "p_value": "Enrichment P-value", "term_size": "Enrichment Term Size", "query_size": "Enrichment Query Size", "precision": "Enrichment Precision", "recall": "Enrichment Recall", "description": "Enrichment Term Description", "intersections": "Enrichment intersections"}).reset_index(drop=True)
        # Convert intersections to comma-separated strings
        enrichment_df['Enrichment intersections'] = enrichment_df['Enrichment intersections'].apply(lambda x: ','.join(x) if isinstance(x, list) else x)
        # Write to disk
        enrichment_df.to_csv(f"{output_folder}/gene_set_enrichment/gene_set_enrichment_results.tsv", sep="\t", index=False)
    else:
        enrichment_df = pd.DataFrame()
        dist_mt = pd.DataFrame()
        clustered_go_terms_list = []
    
    # Read again the merged_df from disk
    merged_df = pd.read_csv(variant_combined_path, sep="\t")
    
    # Prepare pathway-PRS weights if clustering was performed
    if analysis_type == "enrichment" and clustered_go_terms_list:
        prepare_pathway_prs_weights(clustered_go_terms_list, merged_df, output_folder, enrichment_df)
        
    # Plots with R
    cmd = f"Rscript {PATH_SCRIPT}/plot_annotation.R -i {output_folder}/variant_annotation_combined.tsv -g {output_folder}/target_annotations/gwas_annotation_target.tsv -o {output_folder}/plots_variant_annotation"
    os.system(cmd)
    
    # Move query file in output folder
    cmd = f"mv {query} {output_folder}/"
    os.system(cmd)
    # Compress output folder
    cmd = f"cd {os.path.dirname(output_folder)} && zip -r {os.path.basename(output_folder)}.zip {os.path.basename(output_folder)}"
    os.system(cmd)
    # Remove output folder after compression
    cmd = f"rm -rf {output_folder}"
    os.system(cmd)
    
    # Send email notification at the end
    email_at_end(email, random_number)    
        
if __name__ == "__main__":
    main()

