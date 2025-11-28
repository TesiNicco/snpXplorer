import sqlite3
import re
from pathlib import Path
import pandas as pd
from flask import Flask, request, jsonify, render_template, Blueprint, redirect, url_for, session
from liftover import get_lifter
import io
import subprocess
import math
from typing import Optional

# Blueprint setup (if needed)
snpbot_bp = Blueprint('snpbot', __name__, template_folder='templates')

# ---------------------------------------------------------
# Configuration
# ---------------------------------------------------------
DATA_PATH = Path("/Users/nicco/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/snpXplorer/Data")
#DATA_PATH = Path('/Data')
DB_FILE = DATA_PATH / "databases/Genes/variants_info.sqlite"

# Liftover object (hg19 → hg38)
lifter = get_lifter("hg19", "hg38")

#app = Flask(__name__)

# ---------------------------------------------------------
# Helpter to sanitize
# ---------------------------------------------------------
def sanitize_for_json(obj):
    """
    Recursively convert NaN/Inf to None so Flask's JSON encoder doesn't choke.
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
def lookup_by_rsid(rsid):
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
def lookup_by_coord_hg38(chrom, pos38):
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
def liftover_hg19_to_hg38(chrom, pos19):
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
def parse_variant_query(q):
    q = q.strip()

    # rsID
    if q.lower().startswith("rs") and q[2:].isdigit():
        return {"type": "rsid", "rsid": q}

    # chr:pos or chr pos (with optional "chr")
    m = re.match(r"^(chr)?(\w+)[\:\s]+(\d+)$", q, flags=re.IGNORECASE)
    if m:
        chrom = m.group(2).upper()
        pos = int(m.group(3))
        return {"type": "coord", "chrom": chrom, "pos": pos}

    return {"type": "invalid"}

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
        print(str(EQTL_DB_FILE))
        print(region)
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

def _open_db_for_query(db_path: str) -> sqlite3.Connection:
    MMAP_BYTES = 1 << 30  # 1 GiB mmap
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=NORMAL;")
    conn.execute(f"PRAGMA mmap_size={MMAP_BYTES};")
    return conn

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
        if df_ld.empty:
            return [], []
        return df_ld.to_dict(orient="records"), all_partner_positions
    except Exception as e:
        print("Error in query_ld:", e, flush=True)
        return [], []

#---------------------------------------------------------
# Query GWAS associations
#---------------------------------------------------------
def query_gwas_associations(chrom, pos38):
    """
    Query GWAS associations for a given hg38 chromosome and position.
    Returns a list of dictionaries (JSON-serializable).
    """
    try:
        chrom_clean = str(chrom).upper().replace("CHR", "")
        GWAS_DB_FILE = DATA_PATH / f"databases/haplotypes/All_indep_gwas_sumstats_AI_hg38.txt.gz"
        region = f"{chrom_clean}:{pos38}-{pos38}"
        result = subprocess.run(["tabix", str(GWAS_DB_FILE), region], capture_output=True, check=True, text=True).stdout.strip()
        if result == "":
            return []
        df = pd.read_csv(io.StringIO(result), sep="\t", header=None)
        df.columns = ['gwas_id', 'trait', 'chrom', 'pos_hg19', 'rsid', 'ea', 'nea', 'eaf', 'beta', 'se', 'p', 'n', 'study_id', 'position_hg38']
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
def annotate_ld_partners(chrom, ld_rows):
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
            ld_cadd.append(r)
        # eQTL for partner
        eqtl_rows = query_eqtls(chrom_clean, partner_pos)
        for row in eqtl_rows:
            r = dict(row)
            r["ld_partner_pos"] = partner_pos
            r["ld_rsid"] = rsid
            r["ld_r2"] = r2
            r["ld_dist_bp"] = dist
            ld_eqtl.append(r)
        # sQTL for partner
        sqtl_rows = query_sqtls(chrom_clean, partner_pos)
        for row in sqtl_rows:
            r = dict(row)
            r["ld_partner_pos"] = partner_pos
            r["ld_rsid"] = rsid
            r["ld_r2"] = r2
            r["ld_dist_bp"] = dist
            ld_sqtl.append(r)
    return ld_cadd, ld_eqtl, ld_sqtl

# ---------------------------------------------------------
# Core logic for querying a variant
# ---------------------------------------------------------
def run_variant_query(q, build="hg38"):
    print(f"Running variant query for: {q} (build={build})", flush=True)
    parsed = parse_variant_query(q)

    if parsed["type"] == "invalid":
        return {"error": "Query not recognized. Use rsID (rs123) or chr:pos."}, 400

    # rsID path
    if parsed["type"] == "rsid":
        info = lookup_by_rsid(parsed["rsid"])
        if info is None:
            return {"error": f"{parsed['rsid']} not found"}, 404

        try:
            chr38 = info["chr_hg38"]
            pos38 = info["pos_hg38"]
            # Add CADD annotation
            info["cadd"] = query_cadd_score(chr38, pos38)
            # Add eQTL annotation
            info["eqtl"] = query_eqtls(chr38, pos38)
            # Add sQTL annotation
            info["sqtl"] = query_sqtls(chr38, pos38)
            # Add LD annotation
            info["ld"], all_vars = query_ld(chr38, pos38)
            # Add GWAS associations
            info["gwas"] = query_gwas_associations(chr38, pos38)

            # CADD / eQTL / sQTL for all LD partners
            if info["ld"]:
                ld_cadd, ld_eqtl, ld_sqtl = annotate_ld_partners(chr38, info["ld"])
                info["ld_cadd"] = ld_cadd
                info["ld_eqtl"] = ld_eqtl
                info["ld_sqtl"] = ld_sqtl
            else:
                info["ld_cadd"] = []
                info["ld_eqtl"] = []
                info["ld_sqtl"] = []
        except Exception:
            info.setdefault("cadd", [])
            info.setdefault("eqtl", [])
            info.setdefault("sqtl", [])
            info.setdefault("ld", [])
            info.setdefault("ld_cadd", [])
            info.setdefault("ld_eqtl", [])
            info.setdefault("ld_sqtl", [])
            info.setdefault("gwas", [])

        # record in session
        info["query"] = q
        info["build"] = build
        session['single_variant_query'] = info
        return info, 200

    # coordinate path
    chrom = parsed["chrom"]
    pos = parsed["pos"]

    if build == "hg38":
        info = lookup_by_coord_hg38(chrom, pos)
        if info is None:
            return {"error": f"{chrom}:{pos} not found in hg38"}, 404

        try:
            chr38 = info["chr_hg38"]
            pos38 = info["pos_hg38"]
            # Add CADD annotation
            info["cadd"] = query_cadd_score(chr38, pos38)
            # Add eQTL annotation
            info["eqtl"] = query_eqtls(chr38, pos38)
            # Add sQTL annotation
            info["sqtl"] = query_sqtls(chr38, pos38)
            # Add LD annotation
            info["ld"], all_vars = query_ld(chr38, pos38)
            # Add GWAS associations
            info["gwas"] = query_gwas_associations(chr38, pos38)
            # CADD / eQTL / sQTL for all LD partners
            if info["ld"]:
                ld_cadd, ld_eqtl, ld_sqtl = annotate_ld_partners(chr38, info["ld"])
                info["ld_cadd"] = ld_cadd
                info["ld_eqtl"] = ld_eqtl
                info["ld_sqtl"] = ld_sqtl
            else:
                info["ld_cadd"] = []
                info["ld_eqtl"] = []
                info["ld_sqtl"] = []

        except Exception:
            info.setdefault("cadd", [])
            info.setdefault("eqtl", [])
            info.setdefault("sqtl", [])
            info.setdefault("ld", [])
            info.setdefault("ld_cadd", [])
            info.setdefault("ld_eqtl", [])
            info.setdefault("ld_sqtl", [])
            info.setdefault("gwas", [])

        # record in session
        info["query"] = q
        info["build"] = build
        session['single_variant_query'] = info
        return info, 200

    if build == "hg19":
        lifted = liftover_hg19_to_hg38(chrom, pos)
        if lifted is None:
            return {"error": f"Cannot liftover {chrom}:{pos} (hg19→hg38)"}, 400

        chr38, pos38 = lifted
        info = lookup_by_coord_hg38(chr38, pos38)
        if info is None:
            return {"error": f"{chrom}:{pos}→{chr38}:{pos38} not in DB"}, 404

        info["liftover"] = {
            "from_hg19": f"{chrom}:{pos}",
            "to_hg38": f"{chr38}:{pos38}",
        }
        try:
            # Add CADD annotation based on hg38 location
            info["cadd"] = query_cadd_score(chr38, pos38)
            # Add eQTL annotation based on hg38 location
            info["eqtl"] = query_eqtls(chr38, pos38)
            # Add sQTL annotation based on hg38 location
            info["sqtl"] = query_sqtls(chr38, pos38)
            # Add LD annotation based on hg38 location
            info["ld"], all_vars = query_ld(chr38, pos38)
            # Add GWAS associations
            info["gwas"] = query_gwas_associations(chr38, pos38)
            # CADD / eQTL / sQTL for all LD partners
            if info["ld"]:
                ld_cadd, ld_eqtl, ld_sqtl = annotate_ld_partners(chr38, info["ld"])
                info["ld_cadd"] = ld_cadd
                info["ld_eqtl"] = ld_eqtl
                info["ld_sqtl"] = ld_sqtl
            else:
                info["ld_cadd"] = []
                info["ld_eqtl"] = []
                info["ld_sqtl"] = []
        except Exception:
            info.setdefault("cadd", [])
            info.setdefault("eqtl", [])
            info.setdefault("sqtl", [])
            info.setdefault("ld", [])
            info.setdefault("ld_cadd", [])
            info.setdefault("ld_eqtl", [])
            info.setdefault("ld_sqtl", [])
            info.setdefault("gwas", [])
            
        # record in session
        info["query"] = q
        info["build"] = build
        session['single_variant_query'] = info
        return info, 200
    return {"error": "build must be hg19 or hg38"}, 400

# ---------------------------------------------------------
# API endpoint (JSON)
# ---------------------------------------------------------
@snpbot_bp.route("/api/variant/query")
def variant_query():
    q = request.args.get("query", "").strip()
    build = request.args.get("build", "hg38").lower()

    if not q:
        return jsonify({"error": "Missing ?query= parameter"}), 400

    info, status = run_variant_query(q, build)
    info = sanitize_for_json(info)
    return jsonify(info), status

# ---------------------------------------------------------
# Session endpoint to retrieve last query
# ---------------------------------------------------------
@snpbot_bp.route("/api/variant/last")
def last_variant():
    info = session.get("single_variant_query")
    if not info:
        return jsonify({"error": "no previous query"}), 404
    return jsonify(sanitize_for_json(info)), 200

# ---------------------------------------------------------
# Main page (search UI)
# ---------------------------------------------------------
@snpbot_bp.route("/")
def index():
    return render_template("variant_search.html")

# ---------------------------------------------------------
# Run
# ---------------------------------------------------------
#if __name__ == "__main__":
#    app.run(port=5001, debug=True)