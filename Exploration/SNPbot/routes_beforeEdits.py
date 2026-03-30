import sqlite3
import re
from pathlib import Path
import pandas as pd
from flask import Flask, request, jsonify, render_template, Blueprint, redirect, url_for, session
from liftover import get_lifter
import io
import os
import subprocess
from datetime import timedelta, datetime
import math
from typing import Optional

# Blueprint setup (if needed)
snpbot_bp = Blueprint('snpbot', __name__, template_folder='templates')

# ---------------------------------------------------------
# Configuration
# ---------------------------------------------------------
#DATA_PATH = Path("/Users/nicco/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/snpXplorer/Data")
DATA_PATH = Path('/Data')
DB_FILE = DATA_PATH / "databases/Genes/variant_info.db"

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
    return info

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
    return info

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
        # sort by pval_nominal
        df = df.sort_values(by="pval_nominal", ascending=True)
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
            info = fetch_by_position(DB_FILE, chrom, v)
            if info is not None:
                records.append({"partner_pos": int(info["Position_hg38"]), "rsid": info["rsID"]})
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
        # otherwise add also query info
        ld_query_alleles = query_ld_query(chrom_clean, pos38)
        df_ld["query"] = chrom_clean + ':' + str(pos38) + ':' + ld_query_alleles
        return df_ld.to_dict(orient="records"), all_partner_positions
    except Exception as e:
        print("Error in query_ld:", e, flush=True)
        return [], []

def query_ld_query(chrom, pos38):
    """
    Query query-SNP for given chromosome and position.
    """
    try:
        db_path = DATA_PATH / f"databases/LD_db/ld_chr{chrom}.sqlite"
        with _open_db_for_query(db_path) as conn:
            rows = conn.execute("SELECT id, uniq FROM variants WHERE pos = ? ORDER BY uniq",(pos38,)).fetchall()
        if not rows:
            return ''
        elif len(rows) == 1:
            return rows[0][1].split(':')[1] + ':' + rows[0][1].split(':')[2]
        else:
            # return the first one (should not happen)
            return rows[0][1].split(':')[1] + ':' + rows[0][1].split(':')[2]
    except:
        return ''

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
        a1, a2 = rec["partner_uniq"].split(":")[1::]
        # CADD for partner
        cadd_rows = query_cadd_score(chrom_clean, partner_pos)
        # keep only hits where the alleles match the query
        cadd_rows = [x for x in cadd_rows if ( (x['ref'] == a1 and x['alt'] == a2) or (x['ref'] == a2 and x['alt'] == a1) )]
        for row in cadd_rows:
            r = dict(row)
            r["ld_partner_pos"] = partner_pos
            r["ld_rsid"] = rsid
            r["ld_r2"] = r2
            r["ld_dist_bp"] = dist
            ld_cadd.append(r)
        # eQTL for partner
        eqtl_rows = query_eqtls(chrom_clean, partner_pos)
        # keep only hits where the alleles match the query
        eqtl_rows = [x for x in eqtl_rows if ( (x['ref'] == a1 and x['alt'] == a2) or (x['ref'] == a2 and x['alt'] == a1) )]
        for row in eqtl_rows:
            r = dict(row)
            r["ld_partner_pos"] = partner_pos
            r["ld_rsid"] = rsid
            r["ld_r2"] = r2
            r["ld_dist_bp"] = dist
            ld_eqtl.append(r)
        # sQTL for partner
        sqtl_rows = query_sqtls(chrom_clean, partner_pos)
        # keep only hits where the alleles match the query
        sqtl_rows = [x for x in sqtl_rows if ( (x['ref'] == a1 and x['alt'] == a2) or (x['ref'] == a2 and x['alt'] == a1) )]
        for row in sqtl_rows:
            r = dict(row)
            r["ld_partner_pos"] = partner_pos
            r["ld_rsid"] = rsid
            r["ld_r2"] = r2
            r["ld_dist_bp"] = dist
            ld_sqtl.append(r)
    return ld_cadd, ld_eqtl, ld_sqtl

# ---------------------------------------------------------
# Function to query SVs
# ---------------------------------------------------------
def extract_sv(chrom, start_pos, refGen, svtypes):
    try:
        # define end_pos
        end_pos = start_pos + 10000
        start_pos = start_pos - 10000
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
        # Convert to dataframe
        svs = pd.DataFrame(svs, columns=["Repeat Class", "Chromosome", "Start (hg38)", "End (hg38)", "Length", "Repeat Name", "Repeat Family", "Color"])
        # Drop Color column
        svs = svs.drop(columns=["Color"])
        # Dataframe to dictionary
        svs = svs.to_dict(orient="records")
    except Exception as e:
        svs = pd.DataFrame(columns=["Repeat Class", "Chromosome", "Start (hg38)", "End (hg38)", "Length", "Repeat Name", "Repeat Family"])
        svs = svs.to_dict(orient="records")
    return svs

# ---------------------------------------------------------
# Identify most likely gene
# ---------------------------------------------------------
def identify_most_likely_gene(info: dict):
    # Set container for genes
    likely_genes = {}
    invalid_queries = []
    # Iterate over snps in info
    try:
        # Extract dataframes
        cadd_df, eqtl_df, sqtl_df, ld_df, gwas_df, ld_cadd_df, ld_eqtl_df, ld_sqtl_df, sv_df = pd.DataFrame(info['cadd_top']), pd.DataFrame(info['eqtl']), pd.DataFrame(info['sqtl']), pd.DataFrame(info['ld']), pd.DataFrame(info['gwas']), pd.DataFrame(info['ld_cadd']), pd.DataFrame(info['ld_eqtl']), pd.DataFrame(info['ld_sqtl']), pd.DataFrame(info['svs'])
        # First check if variant is coding in CADD
        coding_rows = cadd_df[cadd_df["annotypes"].str.contains("coding", case=False, na=False)]
        coding_rows = coding_rows[~coding_rows["annotypes"].str.contains("noncoding", case=False, na=False)]
        cadd_genes = []
        if not coding_rows.empty:
            genes = [x.split(';') for x in coding_rows["genes"].unique().tolist()]
            genes = [gene for sublist in genes for gene in sublist]  # flatten
            # unique genes
            genes = list(set(genes))
            if len(genes) >0:
                likely_genes = {'genes': genes, 'source': 'Coding'}
                # return most likely gene
                info['likely_gene'] = likely_genes
                return info
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
                genes = [x.split(';') for x in coding_ld_rows["genes"].unique().tolist()]
                genes = [gene for sublist in genes for gene in sublist]  # flatten
                # unique genes
                genes = list(set(genes))
                if len(genes) >0:
                    likely_genes = {'genes': genes, 'source': 'LD with Coding'}
                    # return most likely gene
                    info['likely_gene'] = likely_genes
                    return info
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
                likely_genes = {'genes': qtl_combined, 'source': 'QTL'}
                # return most likely gene
                info['likely_gene'] = likely_genes
                return info
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
                likely_genes = {'genes': ld_qtl_combined, 'source': 'LD with QTL'}
                # return most likely gene
                info['likely_gene'] = likely_genes
                return info
        # If none of the above, return closest gene (placeholder)
        likely_genes = {'genes': closest_gene(cadd_df), 'source': 'Closest Gene'}
        # return most likely gene
        info['likely_gene'] = likely_genes
        return info
    except Exception as e:
        likely_genes = {'genes': [], 'source': 'Error'}
        info['likely_gene'] = likely_genes
        return info

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
# Monitor single variant query
# ---------------------------------------------------------
def add_search_to_file(q, build, DATA_PATH):
    log_file = f"{DATA_PATH}/monitor/single_query_logs.txt"
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_entry = f"{timestamp}\t{q}\t{build}\n"
    with open(log_file, "a") as f:
        f.write(log_entry)

# ---------------------------------------------------------
# Core logic for querying a variant
# ---------------------------------------------------------
def run_variant_query(q, build="hg38"):
    print(f"Running variant query for: {q} (build={build})", flush=True)
    parsed = parse_variant_query(q)

    if parsed["type"] == "invalid":
        return {"error": "Query not recognized. Use rsID (rs123) or chr:pos."}, 400

    # Monitoring: record search
    add_search_to_file(q, build, DATA_PATH)
    
    # rsID path
    if parsed["type"] == "rsid":
        info = fetch_by_rsid(DB_FILE, parsed["rsid"])
        if not info:
            return {"error": f"{parsed['rsid']} not found"}, 404

        try:
            chr38 = info["Chromosome"]
            pos38 = info["Position_hg38"]
            # Add CADD annotation
            info["cadd"] = query_cadd_score(chr38, pos38)
            # Only top 100 for the UI
            info["cadd_top"] = sorted(info["cadd"], key=lambda x: x.get("phred_max", 0), reverse=True)[:100]
            # Add eQTL annotation
            info["eqtl"] = query_eqtls(chr38, pos38)
            # Only top 100 for the UI
            info["eqtl_top"] = sorted(info["eqtl"], key=lambda x: x.get("pval_nominal", 1))[:100]
            # Add sQTL annotation
            info["sqtl"] = query_sqtls(chr38, pos38)
            # Only top 100 for the UI
            info["sqtl_top"] = sorted(info["sqtl"], key=lambda x: x.get("pval_nominal", 1))[:100]
            # Add LD annotation
            info["ld"], all_vars = query_ld(chr38, pos38)
            # Add GWAS associations
            info["gwas"] = query_gwas_associations(chr38, pos38)
            # Add structural variant (SV) annotations
            info["svs"] = extract_sv(chr38, pos38, refGen="GRCh38", svtypes=["all"])
            
            # CADD / eQTL / sQTL for all LD partners
            if info["ld"]:
                # If there are >100 LD partners, limit to top 100 by r2
                #if len(info["ld"]) > 100:
                #    info["ld"] = sorted(info["ld"], key=lambda x: x.get("r2", 0), reverse=True)[:100]
                info["ld"] = sorted(info["ld"], key=lambda x: x.get("r2", 0), reverse=True)
                ld_cadd, ld_eqtl, ld_sqtl = annotate_ld_partners(chr38, info["ld"])
                info["ld_cadd"] = ld_cadd
                info["ld_eqtl"] = ld_eqtl
                info["ld_sqtl"] = ld_sqtl
                # Only top 100 for the UI
                info["ld_cadd_top"] = sorted(info["ld_cadd"], key=lambda x: x.get("phred_max", 0), reverse=True)[:100]
                info["ld_eqtl_top"] = sorted(info["ld_eqtl"], key=lambda x: x.get("pval_nominal", 1))[:100]
                info["ld_sqtl_top"] = sorted(info["ld_sqtl"], key=lambda x: x.get("pval_nominal", 1))[:100]
            else:
                info["ld_cadd"] = []
                info["ld_eqtl"] = []
                info["ld_sqtl"] = []
                info["ld_cadd_top"] = []
                info["ld_eqtl_top"] = []
                info["ld_sqtl_top"] = []
            
            # Identify most likely gene
            info = identify_most_likely_gene(info)
    
        except Exception:
            info.setdefault("cadd", [])
            info.setdefault("eqtl", [])
            info.setdefault("sqtl", [])
            info.setdefault("ld", [])
            info.setdefault("ld_cadd", [])
            info.setdefault("ld_eqtl", [])
            info.setdefault("ld_sqtl", [])
            info.setdefault("gwas", [])
            info.setdefault("svs", [])
            info.setdefault("likely_gene", [])

        # record in session
        info["query"] = q
        info["build"] = build
        session['single_variant_query'] = info
        return info, 200

    # coordinate path
    chrom = parsed["chrom"].upper().replace("CHR", "")
    pos = parsed["pos"]

    if build == "hg38":
        info = fetch_by_position(DB_FILE, chrom, pos)
        if info is None:
            return {"error": f"{chrom}:{pos} not found in hg38"}, 404

        try:
            chr38 = info["Chromosome"]
            pos38 = info["Position_hg38"]
            # Add CADD annotation
            info["cadd"] = query_cadd_score(chr38, pos38)
            # Only top 100 for the UI
            info["cadd_top"] = sorted(info["cadd"], key=lambda x: x.get("phred_max", 0), reverse=True)[:100]
            # Add eQTL annotation
            info["eqtl"] = query_eqtls(chr38, pos38)
            # Only top 100 for the UI
            info["eqtl_top"] = sorted(info["eqtl"], key=lambda x: x.get("pval_nominal", 1))[:100]
            # Add sQTL annotation
            info["sqtl"] = query_sqtls(chr38, pos38)
            # Only top 100 for the UI
            info["sqtl_top"] = sorted(info["sqtl"], key=lambda x: x.get("pval_nominal", 1))[:100]
            # Add LD annotation
            info["ld"], all_vars = query_ld(chr38, pos38)
            # Add GWAS associations
            info["gwas"] = query_gwas_associations(chr38, pos38)
            # Add structural variant (SV) annotations
            info["svs"] = extract_sv(chr38, pos38, refGen="GRCh38", svtypes=["all"])
            
            # CADD / eQTL / sQTL for all LD partners
            if info["ld"]:
                # If there are >100 LD partners, limit to top 100 by r2
                if len(info["ld"]) > 100:
                    info["ld"] = sorted(info["ld"], key=lambda x: x.get("r2", 0), reverse=True)[:100]
                ld_cadd, ld_eqtl, ld_sqtl = annotate_ld_partners(chr38, info["ld"])
                info["ld_cadd"] = ld_cadd
                info["ld_eqtl"] = ld_eqtl
                info["ld_sqtl"] = ld_sqtl
                # Only top 100 for the UI
                info["ld_cadd_top"] = sorted(info["ld_cadd"], key=lambda x: x.get("phred_max", 0), reverse=True)[:100]
                info["ld_eqtl_top"] = sorted(info["ld_eqtl"], key=lambda x: x.get("pval_nominal", 1))[:100]
                info["ld_sqtl_top"] = sorted(info["ld_sqtl"], key=lambda x: x.get("pval_nominal", 1))[:100]
            else:
                info["ld_cadd"] = []
                info["ld_eqtl"] = []
                info["ld_sqtl"] = []
                info["ld_cadd_top"] = []
                info["ld_eqtl_top"] = []
                info["ld_sqtl_top"] = []
            
            # Identify most likely gene
            info = identify_most_likely_gene(info)

        except Exception:
            info.setdefault("cadd", [])
            info.setdefault("eqtl", [])
            info.setdefault("sqtl", [])
            info.setdefault("ld", [])
            info.setdefault("ld_cadd", [])
            info.setdefault("ld_eqtl", [])
            info.setdefault("ld_sqtl", [])
            info.setdefault("gwas", [])
            info.setdefault("svs", [])
            info.setdefault("likely_gene", [])

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
        info = fetch_by_position(DB_FILE, chr38, pos38)
        if info is None:
            return {"error": f"{chrom}:{pos}→{chr38}:{pos38} not in DB"}, 404

        info["liftover"] = {
            "from_hg19": f"{chrom}:{pos}",
            "to_hg38": f"{chr38}:{pos38}",
        }
        try:
            # Add CADD annotation based on hg38 location
            info["cadd"] = query_cadd_score(chr38, pos38)
            # Only top 100 for the UI
            info["cadd_top"] = sorted(info["cadd"], key=lambda x: x.get("phred_max", 0), reverse=True)[:100]
            # Add eQTL annotation based on hg38 location
            info["eqtl"] = query_eqtls(chr38, pos38)
            # Only top 100 for the UI
            info["eqtl_top"] = sorted(info["eqtl"], key=lambda x: x.get("pval_nominal", 1))[:100]
            # Add sQTL annotation based on hg38 location
            info["sqtl"] = query_sqtls(chr38, pos38)
            # Only top 100 for the UI
            info["sqtl_top"] = sorted(info["sqtl"], key=lambda x: x.get("pval_nominal", 1))[:100]
            # Add LD annotation based on hg38 location
            info["ld"], all_vars = query_ld(chr38, pos38)
            # Add GWAS associations
            info["gwas"] = query_gwas_associations(chr38, pos38)
            # Add structural variant (SV) annotations
            info["svs"] = extract_sv(chr38, pos38, refGen="GRCh38", svtypes=["all"])
            
            # CADD / eQTL / sQTL for all LD partners
            if info["ld"]:
                # If there are >100 LD partners, limit to top 100 by r2
                if len(info["ld"]) > 100:
                    info["ld"] = sorted(info["ld"], key=lambda x: x.get("r2", 0), reverse=True)[:100]
                ld_cadd, ld_eqtl, ld_sqtl = annotate_ld_partners(chr38, info["ld"])
                info["ld_cadd"] = ld_cadd
                info["ld_eqtl"] = ld_eqtl
                info["ld_sqtl"] = ld_sqtl
                # Only top 100 for the UI
                info["ld_cadd_top"] = sorted(info["ld_cadd"], key=lambda x: x.get("phred_max", 0), reverse=True)[:100]
                info["ld_eqtl_top"] = sorted(info["ld_eqtl"], key=lambda x: x.get("pval_nominal", 1))[:100]
                info["ld_sqtl_top"] = sorted(info["ld_sqtl"], key=lambda x: x.get("pval_nominal", 1))[:100]
            else:
                info["ld_cadd"] = []
                info["ld_eqtl"] = []
                info["ld_sqtl"] = []
                info["ld_cadd_top"] = []
                info["ld_eqtl_top"] = []
                info["ld_sqtl_top"] = []
            
            # Identify most likely gene
            info = identify_most_likely_gene(info)
            
        except Exception:
            info.setdefault("cadd", [])
            info.setdefault("eqtl", [])
            info.setdefault("sqtl", [])
            info.setdefault("ld", [])
            info.setdefault("ld_cadd", [])
            info.setdefault("ld_eqtl", [])
            info.setdefault("ld_sqtl", [])
            info.setdefault("gwas", [])
            info.setdefault("svs", [])
            info.setdefault("likely_gene", [])
            
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