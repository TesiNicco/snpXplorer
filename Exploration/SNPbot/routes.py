import sqlite3
import re
from pathlib import Path
import pandas as pd
from flask import Flask, request, jsonify, render_template, Blueprint, redirect, url_for, session
from liftover import get_lifter
import io
import os
import pickle
import subprocess
import sys
from datetime import timedelta, datetime
import math
from typing import Optional

try:
    from .gtex_alphagenome_mapping import gtex_to_alphagenome_terms
except ImportError:
    from gtex_alphagenome_mapping import gtex_to_alphagenome_terms

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

def normalize_tissue_selection(tissues):
    """
    Normalize tissue selection to a clean list of strings.
    Accepts None, a single string, or a list/tuple.
    """
    if tissues is None:
        return ["all_tissues"]
    if isinstance(tissues, str):
        tissues = [x.strip() for x in tissues.split(",") if x.strip()]
    elif isinstance(tissues, (list, tuple, set)):
        tissues = [str(x).strip() for x in tissues if str(x).strip()]
    else:
        return ["all_tissues"]
    if not tissues:
        return ["all_tissues"]
    tissue_keys = {t.lower() for t in tissues}
    if "all_tissues" in tissue_keys:
        return ["all_tissues"]
    return tissues

def map_gtex_to_alphagenome_terms(tissues):
    """
    Convert GTEx tissue identifiers to AlphaGenome biosample keyword terms.
    """
    tissues = normalize_tissue_selection(tissues)
    if "all_tissues" in {t.lower() for t in tissues}:
        return []

    keywords = []
    seen = set()
    for tissue in tissues:
        for term in gtex_to_alphagenome_terms.get(tissue, []):
            term_clean = str(term).strip()
            if not term_clean:
                continue
            key = term_clean.lower()
            if key in seen:
                continue
            seen.add(key)
            keywords.append(term_clean)
    return keywords

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
def refseq_accession_to_chrom(accession):
    accession = accession.upper()
    if accession == "NC_012920.1":
        return "MT"
    m = re.match(r"^NC_0*(\d+)\.\d+$", accession)
    if not m:
        return None
    chrom_num = int(m.group(1))
    if 1 <= chrom_num <= 22:
        return str(chrom_num)
    if chrom_num == 23:
        return "X"
    if chrom_num == 24:
        return "Y"
    return None

def parse_variant_query(q):
    q = q.strip()
    # rsID
    if q.lower().startswith("rs") and q[2:].isdigit():
        return {"type": "rsid", "rsid": q}
    # SPDI accession:position:deleted:inserted
    m = re.match(r"^(NC_[0-9]+\.[0-9]+):(\d+):([^:]*):([^:]*)$", q, flags=re.IGNORECASE)
    if m:
        chrom = refseq_accession_to_chrom(m.group(1))
        if chrom is not None:
            pos = int(m.group(2)) + 1
            return {"type": "coord", "chrom": chrom, "pos": pos}
    # gnomAD-like chr-pos-ref-alt or chr-pos
    m = re.match(r"^(chr)?([A-Za-z0-9]+)-(\d+)(?:-[^-]+(?:-[^-]+)?)?$", q, flags=re.IGNORECASE)
    if m:
        chrom = m.group(2).upper()
        pos = int(m.group(3))
        return {"type": "coord", "chrom": chrom, "pos": pos}
    # chr:pos or chr pos (with optional "chr")
    m = re.match(r"^(chr)?(\w+)[\:\s]+(\d+)$", q, flags=re.IGNORECASE)
    if m:
        chrom = m.group(2).upper()
        pos = int(m.group(3))
        return {"type": "coord", "chrom": chrom, "pos": pos}
    return {"type": "invalid"}

# ---------------------------------------------------------
# Query gnomAD information -- frequency and ClinVar annotations
# ---------------------------------------------------------
def query_gnomad_info(info):
    """
    Query gnomAD for cohort frequencies and ClinVar annotations for a variant.
    Expects `info` to contain Chromosome, Position_hg38, REF, and ALT.
    """
    import json
    import urllib.request

    chrom = str(info.get("Chromosome", "")).upper().replace("CHR", "")
    pos38 = info.get("Position_hg38")
    ref = info.get("REF")
    alt = info.get("ALT")

    result = {
        "variant_id": None,
        "dataset": "gnomad_r4",
        "total_af": None,
        "ancestry_af": {
            "african": None,
            "european_non_finnish": None,
            "east_asian": None,
            "middle_east": None,
            "finnish": None,
            "south_asian": None,
            "admixed_american": None,
        },
        "clinvar": {
            "variation_id": None,
            "ref": None,
            "alt": None,
            "gene_symbol": None,
            "conditions": [],
            "germline_classification": None,
            "clinical_significance": None,
            "major_consequence": None,
            "transcript_id": None,
            "hgvsc": None,
            "hgvsp": None,
            "last_evaluated": None,
            "review_status": None,
        },
        "links": {
            "gnomad": None,
            "clinvar": None,
        },
    }

    if not chrom or pos38 is None or not ref or not alt:
        return result

    variant_id = f"{chrom}-{int(pos38)}-{ref}-{alt}"
    result["variant_id"] = variant_id
    result["links"]["gnomad"] = f"https://gnomad.broadinstitute.org/variant/{variant_id}?dataset=gnomad_r4"

    population_map = {
        "afr": "african",
        "nfe": "european_non_finnish",
        "eas": "east_asian",
        "mid": "middle_east",
        "fin": "finnish",
        "sas": "south_asian",
        "amr": "admixed_american",
    }

    def compute_af(ac, an):
        try:
            ac = int(ac)
            an = int(an)
        except (TypeError, ValueError):
            return None
        if an <= 0:
            return None
        return ac / an

    def population_totals(payload):
        totals = {key: {"ac": 0, "an": 0} for key in population_map}
        if not isinstance(payload, dict):
            return totals
        for pop in payload.get("populations") or []:
            pop_id = pop.get("id")
            if pop_id not in population_map:
                continue
            totals[pop_id]["ac"] += int(pop.get("ac") or 0)
            totals[pop_id]["an"] += int(pop.get("an") or 0)
        return totals

    try:
        graphql_query = """
        query VariantSummary($variantId: String!, $dataset: DatasetId!, $referenceGenome: ReferenceGenomeId!) {
          variant(variantId: $variantId, dataset: $dataset) {
            variant_id
            transcript_consequences {
              gene_symbol
              transcript_id
              hgvsc
              hgvsp
              major_consequence
            }
            exome {
              ac
              an
              af
              populations {
                id
                ac
                an
              }
            }
            genome {
              ac
              an
              af
              populations {
                id
                ac
                an
              }
            }
          }
          clinvar_variant(variant_id: $variantId, reference_genome: $referenceGenome) {
            clinvar_variation_id
            ref
            alt
            clinical_significance
            last_evaluated
            review_status
            submissions {
              conditions {
                name
                medgen_id
              }
            }
          }
        }
        """
        graphql_payload = json.dumps(
            {
                "query": graphql_query,
                "variables": {
                    "variantId": variant_id,
                    "dataset": "gnomad_r4",
                    "referenceGenome": "GRCh38",
                },
            }
        ).encode("utf-8")
        graphql_request = urllib.request.Request(
            "https://gnomad.broadinstitute.org/api",
            data=graphql_payload,
            headers={"Content-Type": "application/json"},
            method="POST",
        )
        with urllib.request.urlopen(graphql_request, timeout=15) as response:
            graphql_response = json.loads(response.read().decode("utf-8"))
    except Exception:
        return result

    data = graphql_response.get("data") or {}
    variant_data = data.get("variant") or {}
    clinvar_data = data.get("clinvar_variant") or {}
    transcript_consequences = variant_data.get("transcript_consequences") or []
    transcript_consequence = transcript_consequences[0] if transcript_consequences else {}

    exome = variant_data.get("exome") or {}
    genome = variant_data.get("genome") or {}
    total_ac = int(exome.get("ac") or 0) + int(genome.get("ac") or 0)
    total_an = int(exome.get("an") or 0) + int(genome.get("an") or 0)
    result["total_af"] = compute_af(total_ac, total_an)

    exome_pops = population_totals(exome)
    genome_pops = population_totals(genome)
    for pop_id, label in population_map.items():
        combined_ac = exome_pops[pop_id]["ac"] + genome_pops[pop_id]["ac"]
        combined_an = exome_pops[pop_id]["an"] + genome_pops[pop_id]["an"]
        result["ancestry_af"][label] = compute_af(combined_ac, combined_an)

    clinvar_variation_id = clinvar_data.get("clinvar_variation_id")
    if clinvar_variation_id:
        result["clinvar"]["variation_id"] = clinvar_variation_id
        result["links"]["clinvar"] = f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar_variation_id}/"
        result["clinvar"]["ref"] = clinvar_data.get("ref")
        result["clinvar"]["alt"] = clinvar_data.get("alt")
        condition_names = []
        seen_conditions = set()
        for submission in clinvar_data.get("submissions") or []:
            for condition in submission.get("conditions") or []:
                condition_name = condition.get("name")
                if not condition_name or condition_name in seen_conditions:
                    continue
                seen_conditions.add(condition_name)
                condition_names.append(condition_name)

        result["clinvar"]["conditions"] = condition_names
        result["clinvar"]["gene_symbol"] = transcript_consequence.get("gene_symbol")
        result["clinvar"]["germline_classification"] = clinvar_data.get("clinical_significance")
        result["clinvar"]["clinical_significance"] = clinvar_data.get("clinical_significance")
        result["clinvar"]["major_consequence"] = transcript_consequence.get("major_consequence")
        result["clinvar"]["transcript_id"] = transcript_consequence.get("transcript_id")
        result["clinvar"]["hgvsc"] = transcript_consequence.get("hgvsc")
        result["clinvar"]["hgvsp"] = transcript_consequence.get("hgvsp")
        result["clinvar"]["last_evaluated"] = clinvar_data.get("last_evaluated")
        result["clinvar"]["review_status"] = clinvar_data.get("review_status")

    return result

# ---------------------------------------------------------
# Query alphagenome annotation
# ---------------------------------------------------------
def alphagenome_annotation(info, tissues=None):
    """
    Run AlphaGenome scoring for a single variant and return a pandas DataFrame.
    """
    try:
        chrom = str(info["Chromosome"]).upper().replace("CHR", "")
        pos38 = int(info["Position_hg38"])
        ref = str(info["REF"]).upper()
        alt = str(info["ALT"]).upper()
        variant_id = f"chr{chrom}:{pos38}:{ref}:{alt}"
        alphagenome_dir = DATA_PATH.parent / "Exploration" / "SNPbot"
        script_path = alphagenome_dir / "alphagenome_score.py"
        cmd = [
            sys.executable,
            str(script_path),
            variant_id,
            "--sequence-length", "16KB",
            "--organism", "human",
            "--output-root", "./",
            "--genes",
            "--stdout",
        ]
        alpha_terms = map_gtex_to_alphagenome_terms(tissues)
        if alpha_terms:
            cmd.extend(["--tissue", ",".join(alpha_terms)])
        result = subprocess.run(
            cmd,
            capture_output=True,
            check=True,
            cwd=str(alphagenome_dir),
        )
        if not result.stdout:
            return pd.DataFrame()
        df = pickle.loads(result.stdout)
        if isinstance(df, pd.DataFrame):
            return df
        return pd.DataFrame(df)
    except Exception as e:
        print(f"Error in alphagenome_annotation: {e}", flush=True)
        return pd.DataFrame()

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
def query_eqtls(chrom, pos38, tissues=None):
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
        tissues = normalize_tissue_selection(tissues)
        if "all_tissues" not in {t.lower() for t in tissues}:
            df = df[df["tissue"].isin(tissues)]
        # sort by pval_nominal
        df = df.sort_values(by="pval_nominal", ascending=True)
        return df.to_dict(orient="records")
    except Exception:
        return []

# ---------------------------------------------------------
# Query sQTLs
# ---------------------------------------------------------
def query_sqtls(chrom, pos38, tissues=None):
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
        tissues = normalize_tissue_selection(tissues)
        if "all_tissues" not in {t.lower() for t in tissues}:
            df = df[df["tissue"].isin(tissues)]
        df = df.sort_values(by="pval_nominal", ascending=True)
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
def annotate_ld_partners(chrom, ld_rows, tissues=None):
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
        eqtl_rows = query_eqtls(chrom_clean, partner_pos, tissues=tissues)
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
        sqtl_rows = query_sqtls(chrom_clean, partner_pos, tissues=tissues)
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
        clinvar_gene = (
            info.get("gnomad_clinvar", {})
            .get("clinvar", {})
            .get("gene_symbol")
        )
        if clinvar_gene:
            likely_genes = {'genes': [clinvar_gene], 'source': 'ClinVar'}
            info['likely_gene'] = likely_genes
            return info
        # Extract dataframes
        cadd_df, eqtl_df, sqtl_df, ld_df, gwas_df, ld_cadd_df, ld_eqtl_df, ld_sqtl_df, sv_df = (
            pd.DataFrame(info.get('cadd_top', [])),
            pd.DataFrame(info.get('eqtl', [])),
            pd.DataFrame(info.get('sqtl', [])),
            pd.DataFrame(info.get('ld', [])),
            pd.DataFrame(info.get('gwas', [])),
            pd.DataFrame(info.get('ld_cadd', [])),
            pd.DataFrame(info.get('ld_eqtl', [])),
            pd.DataFrame(info.get('ld_sqtl', [])),
            pd.DataFrame(info.get('svs', [])),
        )
        # First check if variant is coding in CADD
        if "annotypes" in cadd_df.columns:
            coding_rows = cadd_df[cadd_df["annotypes"].str.contains("coding", case=False, na=False)]
            coding_rows = coding_rows[~coding_rows["annotypes"].str.contains("noncoding", case=False, na=False)]
        else:
            coding_rows = pd.DataFrame()
        cadd_genes = []
        if not coding_rows.empty:
            genes_col = coding_rows["genes"] if "genes" in coding_rows.columns else pd.Series([], dtype="object")
            genes = [x.split(';') for x in genes_col.dropna().astype(str).unique().tolist()]
            genes = [gene for sublist in genes for gene in sublist]  # flatten
            # unique genes
            genes = list(set(genes))
            if len(genes) >0:
                likely_genes = {'genes': genes, 'source': 'Coding'}
                # return most likely gene
                info['likely_gene'] = likely_genes
                return info
            else:
                if "genes" in cadd_df.columns:
                    cadd_genes = cadd_df["genes"].dropna().astype(str).unique().tolist()
        else:
            if "genes" in cadd_df.columns:
                cadd_genes = cadd_df["genes"].dropna().astype(str).unique().tolist()
        # Check variants in LD for coding impact
        coding_ld_rows = pd.DataFrame()
        if not ld_cadd_df.empty:
            if "ld_r2" in ld_cadd_df.columns:
                ld_cadd_df = ld_cadd_df[pd.to_numeric(ld_cadd_df["ld_r2"], errors="coerce") >= 0.6]
            if "annotypes" in ld_cadd_df.columns:
                coding_ld_rows = ld_cadd_df[ld_cadd_df["annotypes"].str.contains("coding", case=False, na=False)]
                coding_ld_rows = coding_ld_rows[~coding_ld_rows["annotypes"].str.contains("noncoding", case=False, na=False)]
            else:
                coding_ld_rows = pd.DataFrame()
            if not coding_ld_rows.empty:
                genes_col = coding_ld_rows["genes"] if "genes" in coding_ld_rows.columns else pd.Series([], dtype="object")
                genes = [x.split(';') for x in genes_col.dropna().astype(str).unique().tolist()]
                genes = [gene for sublist in genes for gene in sublist]  # flatten
                # unique genes
                genes = list(set(genes))
                if len(genes) >0:
                    likely_genes = {'genes': genes, 'source': 'LD with Coding'}
                    # return most likely gene
                    info['likely_gene'] = likely_genes
                    return info
                else:
                    if "genes" in ld_cadd_df.columns:
                        cadd_ld_genes = ld_cadd_df["genes"].dropna().astype(str).unique().tolist()
                    else:
                        cadd_ld_genes = []
            else:
                if "genes" in ld_cadd_df.columns:
                    cadd_ld_genes = ld_cadd_df["genes"].dropna().astype(str).unique().tolist()
                else:
                    cadd_ld_genes = []
        else:
            cadd_ld_genes = []
        # Next check QTLs
        if not eqtl_df.empty or not sqtl_df.empty:
            eqtl_genes = qtl_target_ids(eqtl_df)
            sqtl_genes = qtl_target_ids(sqtl_df)
            qtl_combined = list(set(eqtl_genes + sqtl_genes + cadd_genes))
            if len(qtl_combined) >0:
                likely_genes = {'genes': qtl_combined, 'source': 'QTL'}
                # return most likely gene
                info['likely_gene'] = likely_genes
                return info
        # Next check LD QTLs
        if not ld_eqtl_df.empty or not ld_sqtl_df.empty:
            ld_eqtl_genes = qtl_target_ids(ld_eqtl_df)
            ld_sqtl_genes = qtl_target_ids(ld_sqtl_df)
            ld_qtl_combined = list(set(ld_eqtl_genes + ld_sqtl_genes + cadd_ld_genes))
            if len(ld_qtl_combined) >0:
                likely_genes = {'genes': ld_qtl_combined, 'source': 'LD with QTL'}
                # return most likely gene
                info['likely_gene'] = likely_genes
                return info
        # Next check AlphaGenome recurrence among top absolute scores
        alpha_genes = alphagenome_likely_genes(info.get("alphagenome", []))
        if len(alpha_genes) > 0:
            likely_genes = {'genes': alpha_genes, 'source': 'AlphaGenome'}
            info['likely_gene'] = likely_genes
            return info
        # Next check LD-CADD genes (including noncoding) before position fallback
        if len(cadd_ld_genes) > 0:
            genes = [x.split(';') for x in cadd_ld_genes if isinstance(x, str)]
            genes = [gene.strip() for sublist in genes for gene in sublist if gene and gene.strip()]
            genes = list(dict.fromkeys(genes))  # unique, preserve order
            if len(genes) > 0:
                likely_genes = {'genes': genes, 'source': 'LD (CADD)'}
                info['likely_gene'] = likely_genes
                return info
        # If none of the above, return closest gene (placeholder)
        likely_genes = {'genes': closest_gene(info), 'source': 'Closest Gene'}
        # return most likely gene
        info['likely_gene'] = likely_genes
        return info
    except Exception as e:
        likely_genes = {'genes': [], 'source': 'Error'}
        info['likely_gene'] = likely_genes
        return info

def alphagenome_likely_genes(alpha_rows, top_n=100, min_fraction=0.20):
    """
    Identify recurrent AlphaGenome-supported genes from the top absolute scores.
    """
    if not alpha_rows:
        return []

    try:
        alpha_df = pd.DataFrame(alpha_rows)
        if alpha_df.empty or "gene_name" not in alpha_df.columns:
            return []

        alpha_df = alpha_df.copy()
        alpha_df["gene_name"] = alpha_df["gene_name"].astype(str).str.strip()
        alpha_df = alpha_df[alpha_df["gene_name"] != ""]
        alpha_df = alpha_df[alpha_df["gene_name"].str.lower() != "nan"]
        if alpha_df.empty:
            return []

        if "score" in alpha_df.columns:
            alpha_df["score_numeric"] = pd.to_numeric(alpha_df["score"], errors="coerce")
            alpha_df["abs_score"] = alpha_df["score_numeric"].abs()
            alpha_df = alpha_df.sort_values(by="abs_score", ascending=False, na_position="last")

        alpha_top = alpha_df.head(top_n)
        if alpha_top.empty:
            return []

        gene_counts = alpha_top["gene_name"].value_counts(normalize=True)
        supported = gene_counts[gene_counts >= min_fraction]
        return supported.index.tolist()
    except Exception:
        return []

def qtl_target_ids(qtl_df):
    """
    Prefer gene symbols for QTL targets, but fall back to Ensembl IDs when the
    symbol is missing.
    """
    if qtl_df is None or qtl_df.empty:
        return []

    targets = []
    seen = set()

    for _, row in qtl_df.iterrows():
        gene = row.get("gene")
        ensemble = row.get("ensemble")

        gene = "" if pd.isna(gene) else str(gene).strip()
        ensemble = "" if pd.isna(ensemble) else str(ensemble).strip()

        value = gene or ensemble
        if not value:
            continue

        key = value.lower()
        if key in seen:
            continue
        seen.add(key)
        targets.append(value)

    return targets

# ---------------------------------------------------------
# Identify closest gene
# ---------------------------------------------------------
def closest_gene(query_info):
    # get chromosome and position
    chrom = str(query_info["Chromosome"])
    if 'chr' not in chrom:
        chrom = 'chr' + chrom
    start_pos = query_info["Position_hg38"] - 500000
    end_pos = query_info["Position_hg38"] + 500000
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
    variant_pos = query_info["Position_hg38"]
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

def resolve_variant_info(q, build="hg38"):
    """
    Resolve a query to the base variant info record with hg38 alleles/position.
    """
    parsed = parse_variant_query(q)
    if parsed["type"] == "invalid":
        return {"error": "Query not recognized. Use rsID (rs123), chr:pos, chr pos, gnomAD-style chr-pos-ref-alt, or SPDI accession:position:deleted:inserted."}, 400

    if parsed["type"] == "rsid":
        info = fetch_by_rsid(DB_FILE, parsed["rsid"])
        if not info:
            return {"error": f"{parsed['rsid']} not found"}, 404
        return info, 200

    chrom = parsed["chrom"].upper().replace("CHR", "")
    pos = parsed["pos"]

    if build == "hg38":
        info = fetch_by_position(DB_FILE, chrom, pos)
        if info is None:
            return {"error": f"{chrom}:{pos} not found in hg38"}, 404
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
        return info, 200

    return {"error": "build must be hg19 or hg38"}, 400

# ---------------------------------------------------------
# Core logic for querying a variant
# ---------------------------------------------------------
def run_variant_query(q, build="hg38", tissues=None):
    print(f"Running variant query for: {q} (build={build})", flush=True)
    tissues = normalize_tissue_selection(tissues)

    # Monitoring: record search
    add_search_to_file(q, build, DATA_PATH)

    info, status = resolve_variant_info(q, build)
    if status != 200:
        return info, status

    # rsID path
    parsed = parse_variant_query(q)
    if parsed["type"] == "rsid":

        try:
            chr38 = info["Chromosome"]
            pos38 = info["Position_hg38"]
            # Add gnomAD/ClinVar annotation
            info["gnomad_clinvar"] = query_gnomad_info(info)
            # Add CADD annotation
            info["cadd"] = query_cadd_score(chr38, pos38)
            # Only top 100 for the UI
            info["cadd_top"] = sorted(info["cadd"], key=lambda x: x.get("phred_max", 0), reverse=True)[:100]
            # Add eQTL annotation
            info["eqtl"] = query_eqtls(chr38, pos38, tissues=tissues)
            # Only top 100 for the UI
            info["eqtl_top"] = sorted(info["eqtl"], key=lambda x: x.get("pval_nominal", 1))[:100]
            # Add sQTL annotation
            info["sqtl"] = query_sqtls(chr38, pos38, tissues=tissues)
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
                ld_cadd, ld_eqtl, ld_sqtl = annotate_ld_partners(chr38, info["ld"], tissues=tissues)
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
            info.setdefault("gnomad_clinvar", {})
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
        info["tissues"] = tissues
        session['single_variant_query'] = info
        return info, 200

    # coordinate path
    chrom = parsed["chrom"].upper().replace("CHR", "")
    pos = parsed["pos"]

    if build == "hg38":
        try:
            chr38 = info["Chromosome"]
            pos38 = info["Position_hg38"]
            # Add gnomAD/ClinVar annotation
            info["gnomad_clinvar"] = query_gnomad_info(info)
            # Add CADD annotation
            info["cadd"] = query_cadd_score(chr38, pos38)
            # Only top 100 for the UI
            info["cadd_top"] = sorted(info["cadd"], key=lambda x: x.get("phred_max", 0), reverse=True)[:100]
            # Add eQTL annotation
            info["eqtl"] = query_eqtls(chr38, pos38, tissues=tissues)
            # Only top 100 for the UI
            info["eqtl_top"] = sorted(info["eqtl"], key=lambda x: x.get("pval_nominal", 1))[:100]
            # Add sQTL annotation
            info["sqtl"] = query_sqtls(chr38, pos38, tissues=tissues)
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
                ld_cadd, ld_eqtl, ld_sqtl = annotate_ld_partners(chr38, info["ld"], tissues=tissues)
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
            info.setdefault("gnomad_clinvar", {})
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
        info["tissues"] = tissues
        session['single_variant_query'] = info
        return info, 200

    if build == "hg19":
        chr38 = info["Chromosome"]
        pos38 = info["Position_hg38"]
        try:
            # Add gnomAD/ClinVar annotation based on hg38 location
            info["gnomad_clinvar"] = query_gnomad_info(info)
            # Add CADD annotation based on hg38 location
            info["cadd"] = query_cadd_score(chr38, pos38)
            # Only top 100 for the UI
            info["cadd_top"] = sorted(info["cadd"], key=lambda x: x.get("phred_max", 0), reverse=True)[:100]
            # Add eQTL annotation based on hg38 location
            info["eqtl"] = query_eqtls(chr38, pos38, tissues=tissues)
            # Only top 100 for the UI
            info["eqtl_top"] = sorted(info["eqtl"], key=lambda x: x.get("pval_nominal", 1))[:100]
            # Add sQTL annotation based on hg38 location
            info["sqtl"] = query_sqtls(chr38, pos38, tissues=tissues)
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
                ld_cadd, ld_eqtl, ld_sqtl = annotate_ld_partners(chr38, info["ld"], tissues=tissues)
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
            info.setdefault("gnomad_clinvar", {})
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
        info["tissues"] = tissues
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
    tissues = request.args.getlist("tissues")
    if not tissues:
        tissues_param = request.args.get("tissues", "").strip()
        tissues = tissues_param if tissues_param else None

    if not q:
        return jsonify({"error": "Missing ?query= parameter"}), 400

    info, status = run_variant_query(q, build, tissues=tissues)
    info = sanitize_for_json(info)
    return jsonify(info), status

# ---------------------------------------------------------
# AlphaGenome endpoint (JSON)
# ---------------------------------------------------------
@snpbot_bp.route("/api/variant/alphagenome")
def variant_alphagenome():
    q = request.args.get("query", "").strip()
    build = request.args.get("build", "hg38").lower()
    tissues = request.args.getlist("tissues")
    if not tissues:
        tissues_param = request.args.get("tissues", "").strip()
        tissues = tissues_param if tissues_param else None
    tissues = normalize_tissue_selection(tissues)

    if not q:
        return jsonify({"error": "Missing ?query= parameter"}), 400

    session_info = session.get("single_variant_query")
    if session_info and session_info.get("query") == q and session_info.get("build") == build:
        info = dict(session_info)
    else:
        info, status = resolve_variant_info(q, build)
        if status != 200:
            return jsonify(sanitize_for_json(info)), status

    alpha_df = alphagenome_annotation(info, tissues=tissues)
    alpha_records = alpha_df.to_dict(orient="records") if not alpha_df.empty else []
    info["alphagenome"] = alpha_records

    if {"cadd_top", "eqtl", "sqtl", "ld", "gwas", "ld_cadd", "ld_eqtl", "ld_sqtl", "svs"}.issubset(info.keys()):
        info = identify_most_likely_gene(info)

    if session_info and session_info.get("query") == q and session_info.get("build") == build:
        session["single_variant_query"] = info

    payload = {
        "query": q,
        "build": build,
        "tissues": tissues,
        "alphagenome": alpha_records,
        "likely_gene": info.get("likely_gene"),
    }
    return jsonify(sanitize_for_json(payload)), 200

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
