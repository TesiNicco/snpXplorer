from flask import Flask, request, render_template, send_file, session, make_response, redirect, url_for, jsonify, Response, stream_with_context, Blueprint
from queue import Queue, Empty
from uuid import uuid4
import json
import requests
from urllib.parse import urljoin
from flask_sqlalchemy import SQLAlchemy
from flask_session import Session
from werkzeug.security import generate_password_hash, check_password_hash
import redis
import matplotlib
matplotlib.use('Agg')
import logging
import traceback
import polars as pl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_hex
import io
from contextlib import redirect_stdout
from pathlib import Path
from datetime import datetime, timezone
from io import StringIO
import subprocess
import math
import textwrap
import pickle
import threading
import pandas as pd
from typing import Optional, Sequence, Tuple, Any
import re
import numpy as np
import webcolors
import random
import jinja2
import matplotlib.colors as mcolors
import base64
import statsmodels.api as sm
from liftover import get_lifter
import os
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
with redirect_stdout(io.StringIO()):
    from ieugwaspy import api, query, config as ieugwas_config
from datetime import timedelta
import sqlite3
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import plotly.io as pio
import plotly.colors as pc
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.pairwise import cosine_similarity
from scipy.cluster.hierarchy import dendrogram, linkage

# Import Blueprint
from SNPbot.routes import snpbot_bp
from PRSpage.prs import prs_bp

# Import all functions from snpxplorer_functions.py
from snpxplorer_functions import *

# differences between local and server
# local
#data_path = '/Users/nicco/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/snpXplorer/Data'
#path_annot = '/Users/nicco/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/snpXplorer'
# server
data_path = '/Data'
path_annot = ''
from pandasgwas import get_variants_by_variant_id

# Initialize the App
app = Flask(__name__)
# Add configuration for the SQLAlchemy database (for usernames and passwords)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///db.sqlite'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)
# Add configuration for the Session to save to redis
app.config['SESSION_TYPE'] = 'redis'
app.config['SESSION_REDIS'] = redis.from_url('redis://localhost:6379')
app.config['SESSION_USE_SIGNER'] = True
app.config['SECRET_KEY'] = 'secret'
# ✔ Set session timeout: 1 hour
app.config['PERMANENT_SESSION_LIFETIME'] = timedelta(minutes=10)
app.config['SESSION_PERMANENT'] = True
Session(app)

# Ensure session is permanent
@app.before_request
def make_session_permanent():
    session.permanent = True

# Initialize Queues
CONSOLE_QUEUES = {}
# Register Blueprint
app.register_blueprint(snpbot_bp, url_prefix='/snpbot')
app.register_blueprint(prs_bp, url_prefix='/prs')

# Function to manage the publishing of messages
def publish(msg, console_id="default"):
    q = CONSOLE_QUEUES.get(console_id)
    if q:
        q.put_nowait(str(msg))

# function to make a search row for the IEU Open GWAS database
def build_search_row(meta: pd.DataFrame, fields=('id','trait','author','population', 'pmid')) -> None:
    """Create a new row '_search' that concatenates colname + selected rows per column."""
    pieces = [pd.Series(meta.columns, index=meta.columns, name='colname')]
    for rn in fields:
        if rn in meta.index:
            pieces.append(meta.loc[rn].rename(rn))
    # concat per column and join into a single lowercase string
    tmp = pd.concat(pieces, axis=1).astype(str)
    search_series = tmp.apply(lambda r: " ".join(x for x in r.values if x and x.lower() != "nan"), axis=1).str.lower()
    meta.loc['_search'] = search_series

# Cache the OpenGWAS catalogue so an unavailable API cannot block every start.
GWAS_META_CACHE = Path(data_path) / "monitor" / "opengwas_metadata.pkl.gz"
GWAS_META_CACHE_MAX_AGE = timedelta(days=7)
GWAS_META_LOCK = threading.RLock()
GWAS_META_REFRESH_STARTED = False


def _validate_gwas_metadata(candidate: pd.DataFrame) -> pd.DataFrame:
    if not isinstance(candidate, pd.DataFrame) or candidate.empty:
        raise ValueError("OpenGWAS returned an empty metadata catalogue.")
    if "id" not in candidate.index or "trait" not in candidate.index:
        raise ValueError("OpenGWAS metadata is missing the required id/trait rows.")
    return candidate


def _fetch_gwas_metadata() -> pd.DataFrame:
    headers = {
        "X-API-SOURCE": f"ieugwaspy/{ieugwas_config.__version__}",
    }
    if ieugwas_config.env.get("jwt"):
        headers["Authorization"] = f"Bearer {ieugwas_config.env['jwt']}"
    if ieugwas_config.env.get("test_mode_key"):
        headers["X-TEST-MODE-KEY"] = ieugwas_config.env["test_mode_key"]

    response = requests.get(
        urljoin(ieugwas_config.env["base_url"], "gwasinfo"),
        headers=headers,
        timeout=(5, 60),
    )
    response.raise_for_status()
    return _validate_gwas_metadata(pd.DataFrame(response.json()))


def _write_gwas_metadata_cache(candidate: pd.DataFrame) -> None:
    GWAS_META_CACHE.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = GWAS_META_CACHE.with_name(
        f".{GWAS_META_CACHE.name}.{os.getpid()}.tmp"
    )
    try:
        candidate.drop(index="_search", errors="ignore").to_pickle(
            temporary_path,
            compression="gzip",
        )
        os.replace(temporary_path, GWAS_META_CACHE)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def _load_initial_gwas_metadata():
    cache_is_stale = True
    if GWAS_META_CACHE.exists():
        try:
            cached = _validate_gwas_metadata(
                pd.read_pickle(GWAS_META_CACHE, compression="gzip")
            )
            cache_age = datetime.now(timezone.utc) - datetime.fromtimestamp(
                GWAS_META_CACHE.stat().st_mtime,
                tz=timezone.utc,
            )
            cache_is_stale = cache_age > GWAS_META_CACHE_MAX_AGE
            build_search_row(cached)
            logging.info(
                "Loaded OpenGWAS metadata from %s (age: %s).",
                GWAS_META_CACHE,
                cache_age,
            )
            return cached, cache_is_stale
        except Exception:
            logging.exception("Could not read the OpenGWAS metadata cache; refreshing it.")

    fresh = _fetch_gwas_metadata()
    _write_gwas_metadata_cache(fresh)
    build_search_row(fresh)
    return fresh, False


def get_gwas_metadata() -> pd.DataFrame:
    with GWAS_META_LOCK:
        return meta


def _refresh_gwas_metadata_cache() -> None:
    global meta
    try:
        refreshed = _fetch_gwas_metadata()
        _write_gwas_metadata_cache(refreshed)
        build_search_row(refreshed)
        with GWAS_META_LOCK:
            meta = refreshed
            app.config["GWAS_META"] = refreshed
        logging.info("OpenGWAS metadata cache refreshed successfully.")
    except Exception:
        logging.exception(
            "OpenGWAS metadata refresh failed; continuing with the cached catalogue."
        )


meta, GWAS_META_CACHE_IS_STALE = _load_initial_gwas_metadata()
app.config["GWAS_META"] = meta


@app.before_request
def refresh_stale_gwas_metadata():
    """Refresh stale metadata after Flask is already ready to serve requests."""
    global GWAS_META_REFRESH_STARTED
    if not GWAS_META_CACHE_IS_STALE or GWAS_META_REFRESH_STARTED:
        return
    with GWAS_META_LOCK:
        if GWAS_META_REFRESH_STARTED:
            return
        GWAS_META_REFRESH_STARTED = True
    threading.Thread(
        target=_refresh_gwas_metadata_cache,
        name="opengwas-metadata-refresh",
        daemon=True,
    ).start()

# Read browsing options for genes and rsids
genes_map = pd.read_csv('%s/databases/Genes/gene_names.txt' %(Path(data_path)), sep='\t')
# remove duplicates
genes_map = genes_map.drop_duplicates()
# remove genes starting with loc, linc, mir
gene_col = genes_map.columns[0]
genes_map = genes_map[~genes_map[gene_col].str.lower().str.startswith(('loc', 'linc', 'mir'))]
genes_lower = genes_map[gene_col].astype(str).str.lower()

# Helper for the SQlite database
rsid_db_path = '%s/databases/Genes/variant_info.db' %(Path(data_path))
def get_db():
    return sqlite3.connect(str(rsid_db_path))

# Create a user class that contains id, username and password
class User(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(20), unique=True, nullable=False)
    password = db.Column(db.String(64), nullable=False)
    email = db.Column(db.String(120), unique=True, nullable=True)

# Create database
with app.app_context():
    db.create_all()

# Function to generate IDs for each page
@app.get("/page")
def page():
    console_id = str(uuid4())
    return render_template("page.html", console_id=console_id)

# SSE backend
@app.get("/events")
def sse_events():
    cid = request.args.get("cid") or "default"
    q = CONSOLE_QUEUES.setdefault(cid, Queue())

    @stream_with_context
    def gen():
        # greet
        yield 'data: {"message":"__hello__"}\n\n'
        while True:
            try:
                msg = q.get(timeout=25)
                yield f'data: {json.dumps({"message": msg})}\n\n'
            except Empty:
                # heartbeat to keep proxies alive
                yield 'data: {"message": ""}\n\n'

    headers = {
        "Cache-Control": "no-cache",
        "X-Accel-Buffering": "no",
        "Connection": "keep-alive",
    }
    return Response(gen(), mimetype="text/event-stream", headers=headers)

# Function to clear the database -- to do so, connext to http://localhost:5000/clear_database
@app.route('/clear_database')
def clear_database():
    db.reflect() # Reflect all tables to Metadata
    db.drop_all() # Drop all those tables
    return "Database Cleared!"

# Registration route
@app.route('/register', methods=['GET', 'POST'])
def register():
    success = None  # set default value for success message
    if request.method == 'POST':
        # getting form data
        username = request.form.get('username')
        password = request.form.get('password')
        email = request.form.get('email')
        confirm_password = request.form.get('confirm_password')

        # validating the form data
        if len(username) < 5:
            error = 'Username must have at least 5 characters.'
            return render_template('user.html', error=error)
        elif len(password) < 8:
            error = 'Password must have at least 8 characters.'
            return render_template('user.html', error=error)
        elif password != confirm_password:
            error = 'Passwords do not match.'
            return render_template('user.html', error=error)

        # hashing and salting the password
        hashed_password = generate_password_hash(password, method='sha256')

        # creating a new user object and adding to the db
        new_user = User(username=username, password=hashed_password, email=email)
        try:
            db.session.add(new_user)
            db.session.commit()
            success = "You have successfully registered. You can Log In now."
        except:
            db.session.rollback()
            error = "Something went wrong during registration. Please try again."
            return render_template('user.html', error_register=error)

        # return success message and render user.html again with success message
        return render_template('user.html', success_register=success)

    return render_template('user.html')

# Login route
@app.route('/login', methods=['GET', 'POST'])
def login():
    # Check if the use already logged in
    if 'user_id' in session:
        return redirect(url_for('dashboard'))
    
    success = None
    error = None
    
    if request.method == 'POST':
        # getting form data
        username = request.form.get('username')
        password = request.form.get('password')

        # querying user from the db
        user = User.query.filter_by(username=username).first()

        # checking if user exists and if password is correct
        if user and check_password_hash(user.password, password):
            success = "Welcome back %s!" %(username)
            # Save the user's ID in the session
            session['user_id'] = user.id
            session['success_login'] = success
            return redirect(url_for('dashboard'))
        else:
            error = "Invalid username or password"

    return render_template('user.html', error_login=error, success_login=success)

# Redirect to user's dashboard
@app.route('/dashboard')
def dashboard():
   return render_template('dashboard.html')

# Index tab
@app.route("/")
def index():
    return render_template("index.html")

# API endpoint to search GWAS traits
@app.get("/api/gwas")
def api_gwas():
    meta_snapshot = get_gwas_metadata()
    q = (request.args.get("q") or "").strip().lower()
    page = int(request.args.get("page", 1))
    page_size = 50
    if not q:
        return jsonify({"results": [], "pagination": {"more": False}})
    # build _search row if missing
    mask = meta_snapshot.loc["_search"].str.contains(q, na=False, regex=False)
    cols = meta_snapshot.columns[mask]
    total = len(cols)
    start = (page - 1) * page_size
    end = min(start + page_size, total)
    results = []
    for cid in cols[start:end]:
        trait = str(meta_snapshot.at["trait", cid]) if "trait" in meta_snapshot.index else ""
        author = str(meta_snapshot.at["author", cid]) if "author" in meta_snapshot.index else ""
        year = str(meta_snapshot.at["year", cid]) if "year" in meta_snapshot.index else ""
        ancestry = gwas_ancestry_code(meta_snapshot, cid)
        label = f"[{ancestry}] {trait} — {author} ({year}) [{cid}]".strip()
        results.append({"id": label, "text": label})
    return jsonify({"results": results, "pagination": {"more": end < total}})

# API endpoint to search genes and rsids
@app.get("/api/locus")
def api_locus():
    q = (request.args.get("q") or "").strip()
    page = int(request.args.get("page", 1))
    page_size = 50
    if not q:
        return jsonify({"results": [], "pagination": {"more": False}})
    # rsID search if it starts with 'rs' (case-insensitive)
    if q.lower().startswith("rs"):
        like = q + '%'
        offset = (page - 1) * page_size
        conn = get_db()
        cur = conn.cursor()
        cur.execute("""
            SELECT rsID, Chromosome, Position_hg38
            FROM variant_info
            WHERE rsID LIKE ?
            ORDER BY (rsID = ?) DESC,     -- exact match first
                    LENGTH(rsID),        -- shorter IDs next
                    rsID
            LIMIT ? OFFSET ?
        """, (like, q, page_size + 1, offset))
        rows = cur.fetchall()
        conn.close()
        results = []
        for rsid, chrom, pos in rows[:page_size]:
            chrom_label = chrom if str(chrom).lower().startswith("chr") else f"chr{chrom}"
            results.append({"id": rsid, "text": f"{rsid} — {chrom_label}:{pos}"})
        return jsonify({"results": results, "pagination": {"more": len(rows) > page_size}})
    # Otherwise: gene name search (case-insensitive substring)
    # Filter in-memory; fast for ~27k names
    mask = genes_lower.str.contains(q.lower(), na=False)
    # Simple pagination
    matches = genes_map.loc[mask, gene_col].astype(str).tolist()
    start = (page - 1) * page_size
    end = start + page_size
    slice_ = matches[start:end]
    results = [{"id": g, "text": g} for g in slice_]
    more = end < len(matches)
    return jsonify({"results": results, "pagination": {"more": more}})

# API endpoint to search traits for haplotypes
@app.get("/api/traits")
def api_traits():
    """
    Return trait suggestions for Select2.
    Uses the same trait file as in haplotypes().
    """
    q = (request.args.get("q") or "").strip()
    page = int(request.args.get("page", 1))
    page_size = 50
    if not q:
        return jsonify({"results": [], "pagination": {"more": False}})
    trait_names = combined_trait_names[combined_trait_names['name'].notna()].copy()
    q_lower = q.lower()
    matches = trait_names[
        trait_names['name'].astype(str).str.lower().str.contains(q_lower, regex=False)
        | trait_names['gwas_id'].astype(str).str.lower().str.contains(q_lower, regex=False)
    ]
    start = (page - 1) * page_size
    end = start + page_size
    slice_ = matches.iloc[start:end]
    more = end < len(matches)
    results = [
        {
            "id": row.selection_id,
            "text": f"[{SEMANTIC_ANCESTRY_CODES.get(row.ancestry, row.ancestry.upper())}] {row.name}",
            "ancestry": row.ancestry,
        }
        for row in slice_.itertuples(index=False)
    ]
    return jsonify({"results": results, "pagination": {"more": more}})

# Function to run the pipeline in a separate thread
# Global, in-process result cache for SSE runs
RESULTS = {}  # { console_id: {template_context...} }
def run_pipeline(console_id, formdata):
    """Background runner that mirrors your current POST body and
    stores the final render context in RESULTS[console_id]."""
    try:
        meta_snapshot = get_gwas_metadata()
        # --- read inputs (same as your POST branch) ---
        refGen = formdata.get("refGenome")
        browse = formdata.get("browse", "")
        gwas = formdata.getlist('gwas_data')
        window = int(formdata.get("window", 5000) or 5000)
        recomb = formdata.get("recomb", "No")
        exons = formdata.get("exons", "No")
        axestype = formdata.get("axestype", "manh")
        plotype = formdata.get("plotype", "Scatter")
        qtl_tissues = formdata.getlist('QTLtissuesExplo')
        ld = formdata.get("ld", "No")
        svtypes = formdata.getlist('sv_source')
        yrange_type = formdata.get("y_range", "Auto")
        yrange = [int(formdata.get("y_min")), int(formdata.get("y_max"))] if yrange_type == "Manual" else "Auto"
        plotHaplo = formdata.get("plotHaplo", "No")
        reference_ancestry = reference_ancestry_for_studies(gwas, meta_snapshot)

        # --- your original steps, with logs ---
        publish("Finding genomic location...", console_id)
        chrom, start_pos, end_pos, browse_type = readBrowseOption(data_path, browse, window, refGen)
        logging.info("All good with readBrowseOption")
        
        # Save search to file
        add_search_to_file(data_path, chrom, start_pos, end_pos, refGen, gwas, browse_type)

        publish("Gathering data of interest...", console_id)
        df, df_info = get_data_plot(data_path=data_path, gwas=gwas, chrom=chrom, start_pos=start_pos, end_pos=end_pos, refGen=refGen, meta=meta_snapshot, window=window)
        genes = extract_genes(data_path, chrom, start_pos, end_pos, refGen)
        svs, svs_df = extract_sv(data_path, chrom, start_pos, end_pos, refGen, svtypes)
        recomb_data = extract_recomb(data_path, chrom, start_pos, end_pos, refGen) if recomb == 'Yes' else "None"
        gtex_df = get_gtex(data_path, genes, refGen)
        ld_df = extract_ld(data_path, chrom, df, refGen, browse_type, browse, reference_ancestry) if ld == "Yes" else "None"
        haplo_df = extract_haplo(data_path, df, chrom, refGen, reference_ancestry) if plotHaplo == "Yes" else "None"
        logging.info("All good with data extraction")

        publish("Rendering plots...", console_id)
        fig = scatterplot_plotly(df=df, chrom=chrom, start_pos=start_pos, end_pos=end_pos, gwas=gwas, genes=genes, svs=svs, browse_type=browse_type, refGen=refGen, recomb_data=recomb_data, exons=exons, browse=browse, plotype=plotype, ld=ld, ld_df=ld_df, yrange=yrange, haplo_df=haplo_df, plotHaplo=plotHaplo, axestype=axestype)
        fig.update_layout(autosize=True, margin=dict(l=40, r=20, t=90, b=40))
        plot_url = pio.to_html(fig, include_plotlyjs='cdn', full_html=False, default_width='100%', default_height='100%', config={'responsive': True})
        img_beta = beta_maf_plot(df, 'Beta', 'EAF', title="Beta vs MAF")
        plot_beta = pio.to_html(img_beta, include_plotlyjs='cdn', full_html=False, default_width='100%', default_height='100%', config={'responsive': True})
        logging.info("All good with scatterplot")
        img_gtex = gtex_heatmap(gtex_df)
        plot_gtex = pio.to_html(img_gtex, include_plotlyjs='cdn', full_html=False, default_width='100%', default_height='100%', config={'responsive': True})
        logging.info("All good with gtex heatmap")
        publish("Preparing tables...", console_id)

        # SNPs table
        df_sorted = df.sort_values(by='Pvalue', ascending=False).copy()
        df_sorted['Locus'] = df_sorted.apply(lambda x: f"{x['Chrom']}:{x['Position']}", axis=1)
        df_sorted['Beta'] = df_sorted['Beta'].round(2)
        df_sorted['SE'] = df_sorted['SE'].round(2)
        df_sorted['Beta_SE'] = df_sorted.apply(lambda x: f"{x['Beta']} ({x['SE']})", axis=1)
        df_sorted['Alleles'] = df_sorted.apply(lambda x: f"{x['EA']}/{x['NEA']}", axis=1)
        df_sorted = df_sorted[['Locus', 'Rsid', 'Gwas', 'Alleles', 'EAF', 'Beta_SE', 'Pvalue']]
        df_sorted['Pvalue'] = df_sorted['Pvalue'].round(2)
        snps_table = df_sorted.to_dict(orient='records')
        logging.info("All good with SNPs table")
        # GWAS table
        df_info_sub = df_info[['trait', 'population', 'id', 'pmid', 'author', 'year', 'sample_size', 'build']].drop_duplicates()
        df_info_sub['Ancestry'] = df_info_sub['population'].apply(normalize_gwas_ancestry)
        df_info_sub = df_info_sub[['trait', 'Ancestry', 'id', 'pmid', 'author', 'year', 'sample_size', 'build']]
        df_info_sub = df_info_sub.rename(columns={'trait':'Trait', 'id':'ID', 'pmid':'PMID', 'author':'Author', 'year':'Year', 'sample_size':'Sample Size', 'build':'Build'})
        df_info_table = df_info_sub.to_dict(orient='records')
        logging.info("All good with GWAS table")
        # SVs table
        svs_df = svs_df.copy()
        svs_df['Region'] = svs_df.apply(lambda x: f"{x['chrom']}:{x['start']}-{x['end']}", axis=1)
        svs_df = svs_df[['Region', 'len', 'repName', 'repClass', 'repFamily']]
        svs_table = svs_df.to_dict(orient='records')
        logging.info("All good with SVs table")
        # GWAS Catalog
        gwas_table, gwascat_df = extract_gwascatalog(data_path, chrom, start_pos, end_pos, refGen)
        logging.info("All good with GWAS Catalog table")
        # Persist for download routes if you still want them after reload
        # (we can’t touch flask.session in a thread, so stash CSVs here too)
        ctx = {
            "plot_url": plot_url,
            "plot_gtex": plot_gtex,
            "plot_beta": plot_beta,
            "table_snps": snps_table,
            "table_svs": svs_table,
            "table_gwas": gwas_table,
            "table_info": df_info_table,
            # anything needed for later or downloads:
            "_csv_payloads": {
                "gtex_df": gtex_df.to_csv(index=False),
                "df": df.to_csv(index=False),
                "svs_df": svs_df.to_csv(index=False),
                "gwascat_df": gwascat_df.to_csv(index=False),
                "df_info_sub": df_info_sub.to_csv(index=False),
            },
            "_state": {
                "chrom": chrom, 
                "start_pos": start_pos, 
                "end_pos": end_pos,
                "gwas": gwas, 
                "genes": genes, 
                "svs": svs, 
                "browse_type": browse_type,
                "refGen": refGen, 
                "recomb": recomb,
                "recomb_data": (recomb_data.to_csv(index=False) if recomb == "Yes" else "None"),
                "ld_df": (ld_df.to_csv(index=False) if ld == "Yes" else "None"),
                "ld": ld,
                "exons": exons, 
                "browse": browse, 
                "axestype": axestype,
                "plotype": plotype,
                "window": window,
                "svtypes": svtypes,
                "yrange_type": yrange_type,
                "yrange": (yrange if yrange == "Auto" else json.dumps(yrange)),
                "plotHaplo": plotHaplo,
                "reference_ancestry": reference_ancestry,
                "haplo_df": (haplo_df.to_csv(index=False) if plotHaplo == "Yes" else "None"),
            }
        }

        RESULTS[console_id] = ctx
        publish("[done]", console_id)

    except Exception as e:
        # Log the error
        logging.exception("Error in run_pipeline")
        try:
            publish(f"[error] {type(e).__name__}: {e}", console_id)
        except:
            logging.exception("publish() failed while reporting an error")

@app.post('/exploration/run')
def exploration_run():
    console_id = request.form.get('console_id', 'default')
    CONSOLE_QUEUES.setdefault(console_id, Queue())  # ensure queue exists
    formdata = request.form.copy()  # pass a safe copy into the thread
    threading.Thread(
        target=run_pipeline,
        args=(console_id, formdata),
        daemon=True
    ).start()

    return jsonify(ok=True), 202

# Added to try to solve the stucking issue
@app.get('/exploration/status')
def exploration_status():
    cid = request.args.get('cid')
    if not cid:
        return jsonify(done=False), 200

    done = cid in RESULTS
    return jsonify(done=done), 200

# Exploration tab
@app.route('/exploration/', methods=["GET", "POST"])
def exploration():
    # 1) If redirected after background run
    cid = request.args.get("cid")
    if cid and cid in RESULTS:
        ctx = RESULTS.pop(cid)
        # Persist into the session (so coming back later still works)
        for k, v in ctx.get("_csv_payloads", {}).items():
            session[k] = v
        for k, v in ctx.get("_state", {}).items():
            session[k] = v

        return render_template(
            "exploration.html",
            plot_url=ctx["plot_url"],
            plot_beta=ctx["plot_beta"],
            plot_gtex=ctx["plot_gtex"],
            table_snps=ctx["table_snps"],
            table_svs=ctx["table_svs"],
            table_gwas=ctx["table_gwas"],
            table_info=ctx["table_info"],
            browse_value=session.get("browse"),
            gwas=session.get("gwas", []),
            console_id=cid,
            # added
            refGenome=session.get("refGen", "GRCh37"),
            window=session.get("window", 5000),
            axestype=session.get("axestype", "manh"),
            plotype=session.get("plotype", "Scatter"),
            exons=session.get("exons", "No"),
            recomb=session.get("recomb", "No"),
            ld=session.get("ld", "No"),
            plotHaplo=session.get("plotHaplo", "No"),
            reference_ancestry=session.get("reference_ancestry", "European"),
            svtypes=session.get("svtypes", ['tr']),
            y_range = session.get("yrange_type", "Auto"),
            y_min = session.get("yrange", [0, 10])[0] if session.get("yrange_type", "Auto") == "Manual" else "",
            y_max = session.get("yrange", [0, 10])[1] if session.get("yrange_type", "Auto") == "Manual" else "",
        )

    # 3) Rebuild from session when returning later without cid
    if 'df' in session:
        df = pd.read_csv(StringIO(session['df']))
        gtex_df = pd.read_csv(StringIO(session['gtex_df']))
        svs_df = pd.read_csv(StringIO(session['svs_df']))
        gwascat = pd.read_csv(StringIO(session['gwascat_df']))
        table_info = pd.read_csv(StringIO(session['df_info_sub']))

        chrom = session['chrom']
        start_pos = session['start_pos']
        end_pos = session['end_pos']
        gwas = session['gwas']
        genes = session['genes']
        svs = session['svs']
        browse_type = session['browse_type']
        refGen = session['refGen']
        recomb_data = session['recomb_data']
        ld = session['ld']
        exons = session['exons']
        browse = session['browse']
        axestype = session['axestype']
        plotype = session['plotype']
        recomb = session['recomb']
        window = session['window']
        svtypes = session['svtypes']
        ld_df = pd.read_csv(StringIO(session['ld_df'])) if 'ld_df' in session and session['ld_df'] != "None" else "None"
        plotHaplo = session['plotHaplo']
        reference_ancestry = session.get('reference_ancestry', 'European')
        haplo_df = pd.read_csv(StringIO(session['haplo_df'])) if 'haplo_df' in session and session['haplo_df'] != "None" else "None"
        yrange_type = session['yrange_type']
        yrange = [int(x) for x in json.loads(session['yrange'])] if yrange_type == "Manual" else "Auto"

        # tables
        df_sorted = df.sort_values(by='Pvalue', ascending=False).copy()
        df_sorted['Locus'] = df_sorted.apply(lambda x: f"{x['Chrom']}:{x['Position']}", axis=1)
        df_sorted['Beta'] = df_sorted['Beta'].round(2)
        df_sorted['SE'] = df_sorted['SE'].round(2)
        df_sorted['Beta_SE'] = df_sorted.apply(lambda x: f"{x['Beta']} ({x['SE']})", axis=1)
        df_sorted['Alleles'] = df_sorted.apply(lambda x: f"{x['EA']}/{x['NEA']}", axis=1)
        df_sorted = df_sorted[['Locus', 'Rsid', 'Gwas', 'Alleles', 'EAF', 'Beta_SE', 'Pvalue']]
        df_sorted['Pvalue'] = df_sorted['Pvalue'].round(2)
        snps_table = df_sorted.to_dict(orient='records')
        svs_table = svs_df.to_dict(orient='records')
        gwas_table = gwascat.to_dict(orient='records')
        info_table = table_info.to_dict(orient='records')

        # plots
        fig = scatterplot_plotly(df=df, chrom=chrom, start_pos=start_pos, end_pos=end_pos, gwas=gwas, genes=genes, svs=svs, browse_type=browse_type, refGen=refGen, recomb_data=recomb_data, exons=exons, browse=browse, plotype=plotype, ld=ld, ld_df=ld_df, yrange=yrange, haplo_df=haplo_df, plotHaplo=plotHaplo, axestype=axestype)
        fig.update_layout(autosize=True, margin=dict(l=40, r=20, t=90, b=40))
        plot_url = pio.to_html(fig, include_plotlyjs='cdn', full_html=False, default_width='100%', default_height='100%', config={'responsive': True})

        img_gtex = gtex_heatmap(gtex_df)
        plot_gtex = pio.to_html(img_gtex, include_plotlyjs='cdn', full_html=False, default_width='100%', default_height='100%', config={'responsive': True})
        
        img_beta = beta_maf_plot(df, 'Beta', 'EAF', title="Beta vs MAF")
        plot_beta = pio.to_html(img_beta, include_plotlyjs='cdn', full_html=False, default_width='100%', default_height='100%', config={'responsive': True})

        return render_template(
            "exploration.html",
            plot_url=plot_url, plot_gtex=plot_gtex, plot_beta=plot_beta,
            table_snps=snps_table, table_svs=svs_table, table_gwas=gwas_table, table_info=info_table,
            browse_value=browse, gwas=gwas,
            console_id=str(uuid4()),
            refGenome=session.get("refGen", "GRCh37"),
            window=session.get("window", 5000),
            axestype=session.get("axestype", "manh"),
            plotype=session.get("plotype", "Scatter"),
            exons=session.get("exons", "No"),
            recomb=session.get("recomb", "No"),
            ld=session.get("ld", "No"),
            plotHaplo=session.get("plotHaplo", "No"),
            reference_ancestry=reference_ancestry,
            svtypes=session.get("svtypes", ['tr']),
            y_range = session.get("yrange_type", "Auto"),
            y_min = session.get("yrange", [0, 10])[0] if session.get("yrange_type", "Auto") == "Manual" else "",
            y_max = session.get("yrange", [0, 10])[1] if session.get("yrange_type", "Auto") == "Manual" else "",
        )
    
    # 4) first load
    return render_template("exploration.html", console_id=str(uuid4()))

# Download SNP table
@app.route('/downloadSNPs')
def download_snp_table():
    if 'df' in session:
        # retrieve table data from session object
        csv_data = session['df']
        # convert CSV string to file-like object
        file_like = StringIO(csv_data)
        # create response with file-like object as download
        response = make_response(file_like.getvalue())
        response.headers["Content-Disposition"] = "attachment; filename=SNPs_table.csv"
        response.headers["Content-type"] = "text/csv"    
        return response
    return 'Data not available. Please go back and generate a plot before downloading.'

# Download SV table
@app.route('/downloadSVs')
def download_sv_table():
    if 'svs_df' in session:
        # retrieve table data from session object
        csv_data = session['svs_df']
        # convert CSV string to file-like object
        file_like = StringIO(csv_data)
        # create response with file-like object as download
        response = make_response(file_like.getvalue())
        response.headers["Content-Disposition"] = "attachment; filename=SVs_table.csv"
        response.headers["Content-type"] = "text/csv"    
        return response
    return 'Data not available. Please go back and generate a plot before downloading.'

# Download GWAS-Catalog table
@app.route('/downloadGWAScat')
def download_gwascat_table():
    if 'gwascat_df' in session:
        # retrieve table data from session object
        csv_data = session['gwascat_df']
        # convert CSV string to file-like object
        file_like = StringIO(csv_data)
        # create response with file-like object as download
        response = make_response(file_like.getvalue())
        response.headers["Content-Disposition"] = "attachment; filename=GWASCatalog_table.csv"
        response.headers["Content-type"] = "text/csv"    
        return response
    return 'Data not available. Please go back and generate a plot before downloading.'

# Download GTEx Expression
@app.route('/downloadGTExExpression')
def download_gtex_expression():
    if 'gtex_df' in session:
        # retrieve table data from session object
        csv_data = session['gtex_df']
        # convert CSV string to file-like object
        file_like = StringIO(csv_data)
        # create response with file-like object as download
        response = make_response(file_like.getvalue())
        response.headers["Content-Disposition"] = "attachment; filename=GTEx_expression.csv"
        response.headers["Content-type"] = "text/csv"    
        return response
    return 'Data not available. Please go back and generate a plot before downloading.'

# Download SNP plot
@app.route('/download_plot/SNP_plot.png')
def download_snp_plot():
    if 'df' in session:
        chrom = session['chrom']
        start_pos = session['start_pos']
        end_pos = session['end_pos']
        gwas = session['gwas']
        genes = session['genes']
        svs = session['svs']
        browse_type = session['browse_type']
        refGen = session['refGen']
        recomb_data = session['recomb_data']
        exons = session['exons']
        browse = session['browse']
        plotype = session['plotype']
        df = pd.read_csv(StringIO(session.get('df')))
        gtex_df = pd.read_csv(StringIO(session.get('gtex_df')))
        svs_df = pd.read_csv(StringIO(session.get('svs_df')))
        gwascat = pd.read_csv(StringIO(session.get('gwascat_df')))
        filename = 'SNP_plot.png'
        # Plot
        img = scatterplot(df = df, chrom = chrom, start_pos = start_pos, end_pos = end_pos, gwas = gwas, genes = genes, svs = svs, browse_type = browse_type, refGen = refGen, recomb_data = recomb_data, exons = exons, browse = browse, plotype = plotype)
        return send_file(img, mimetype='image/png', as_attachment=True, download_name=filename)
    return 'Data not available. Please go back and generate a plot before downloading.'

# Download SNP plot
@app.route('/download_gtex/GTEx_plot.png')
def download_gtex_plot():
    if 'df' in session:
        gtex_df = pd.read_csv(StringIO(session.get('gtex_df')))
        filename = 'GTEx_plot.png'
        # Plot
        img_gtex = gtex_heatmap(gtex_df)
        return send_file(img_gtex, mimetype='image/png', as_attachment=True, download_name=filename)
    return 'Data not available. Please go back and generate a plot before downloading.'

# About tab
@app.route('/about/')
def about():
    return render_template('about.html')

# Submit feedback
@app.route("/submit_feedback", methods=["POST"])
def submit_feedback():
    FEEDBACK_FILE = Path(f"{data_path}/monitor/feedback.jsonl")  # adjust path if needed
    obj = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "name": (request.form.get("name") or "Anonymous").strip() or "Anonymous",
        "rating": int(request.form["rating"]) if request.form.get("rating") else None,
        "message": (request.form.get("message") or "").strip(),
    }
    FEEDBACK_FILE.parent.mkdir(parents=True, exist_ok=True)
    with FEEDBACK_FILE.open("a", encoding="utf-8") as f:
        f.write(json.dumps(obj, ensure_ascii=False) + "\n")
    return redirect(url_for("about", sent=1))

# User tab
@app.route('/user/')
def user():
    if 'user_id' in session:
        return redirect(url_for('dashboard'))
    if request.method == 'POST':
        if request.form['action'] == 'login':
            # Handle login form submission
            return 'Login button pressed'
        elif request.form['action'] == 'register':
            #Handle registration form submission
            return 'Register button pressed'
    else:
        # Render user.html template with GET request
        return render_template('user.html')

# Logout
@app.route('/logout')
def logout():
    # remove the username from the session if it is there
    session.pop('user_id', None)
    return redirect(url_for('login'))

# Annotation section
@app.route('/annotation/', methods=["GET", "POST"])
def annotation():
    #if 'user_id' in session:
    #    message = '%s You will find the results in your dashboard.' %(session['success_login'])
    #else:
    #    message = 'Please enter your email address below. The results will be sent to you:'
    #    email = '<textarea name="emailAnnot" style="width: 100%; height: 30px;" placeholder="Enter your email here"></textarea>'
    # check if the user has logged in -- in this case, show a messages, otherwise, show a textarea for the user to input the email address
    # email
    email = request.form.get('email', '')
    email = email.replace(' ', '').rstrip()
    # check if email is valid
    emailResponse = is_valid_email(email)
    # modify the message linked to submission
    messageSubmission = ''
    run_id = None
    limit_message = ('Input limits: Annotation allows up to 10000 selected variants and '
                     'Gene-set Enrichment Analysis allows up to 1000 explicitly entered variants. '
                     'Please provide one rsID, coordinate, or interval per line.')
    if request.method == "POST":
        # Read inputs
        # list of variants
        textarea_input = request.form.get('SNPlist', '')
        my_list = textarea_input.replace(' ', '').replace(';', '\n').replace(',', '\n').split('\n')
        my_list = [x for x in my_list if x]
        # input type
        #inpType = request.form["inputType"]
        inpType = 'rsid'  # for now, only rsid
        # if the input type is not rsid, then read also the reference genome
        refGeno = request.form["annotRefGen"]
        ancestry = (
            'African'
            if request.form.get('annotAncestry', 'European').strip().lower() == 'african'
            else 'European'
        )
        # analysis type
        analType = request.form["analType"]
        # if the analysis type is GSEA, then look at the enrichment sources
        gsea_source = request.form.getlist("EnrichSource") if analType == 'GSEA' else []
        # QTL tissues
        qtl_tissues = request.form.getlist('QTLtissuesAnnot')
        # Intervals are expanded by the standalone annotation pipeline, where
        # the definitive post-expansion limit is enforced. This route only
        # provides immediate validation for explicitly entered variants.
        interval_pattern = re.compile(
            r'^(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y|MT|M):\d+-\d+$',
            re.IGNORECASE,
        )
        interval_inputs = [value for value in my_list if interval_pattern.fullmatch(value)]
        explicit_inputs = [value for value in my_list if value not in interval_inputs]
        snp_limit = 10000 if analType == 'Annot' else 1000
        if not my_list:
            messageSubmission = f'Please enter at least one variant or interval. {limit_message}'
        elif analType == 'GSEA' and interval_inputs:
            messageSubmission = 'Genomic intervals are available for Annotation only. Please remove the interval or select Annotation.'
        elif len(set(explicit_inputs)) > snp_limit:
            messageSubmission = f'Too many explicit variants provided for the selected analysis. {limit_message}'
        elif emailResponse == True:
            # run the script here
            command, run_id = run_annotation(
                my_list,
                inpType,
                refGeno,
                analType,
                gsea_source,
                qtl_tissues,
                email,
                ancestry,
            )
            print(command)
            subprocess.Popen(command, shell=True)
            messageSubmission = f'Your job has been submitted. Save Run ID {run_id} to download the results when ready.'
        else:
            messageSubmission = 'Email is not correct. Please check!'
    return render_template('annotation.html', validEmail=emailResponse, messageSubmission=messageSubmission, run_id=run_id)

# Help tab
@app.route('/help/')
def help():
    return render_template('help.html')

# Download tab
@app.route('/download/', methods=["GET", "POST"])
def download():
    run_id = None
    messageError = None
    messageToUser = None
    if request.method == 'POST':
        run_id = request.form["run_id"]
        messageError, messageToUser = check_runID(run_id)
    return render_template('download.html', run_id=run_id, messageError=messageError, messageToUser=messageToUser)

# Actual download button for the annotation results
@app.route('/downloadResults/<run_id>', methods=["GET"])
def download_results(run_id):
    # check here the run ID
    filename = f"snpXplorer_annotation_{run_id}.zip"
    filepath = os.path.join(path_annot + "/Annotation/RUNS", filename)
    return send_file(filepath, as_attachment=True, download_name=filename)

# Haplotype tab
@app.route('/haplotype/', methods=["GET", "POST"])
def haplotypes():
    meta_snapshot = get_gwas_metadata()
    trait_names = combined_trait_names.copy()
    trait_list = trait_names['name'].values.tolist()
    if request.method == "POST":
        # Clear previous detail state on new query
        session.pop('haplo_last_detail', None)
        # Read inputs
        browse = (request.form.get("browse") or "").strip()
        browse_type = request.form.get("browse_type")
        region_ancestry = (
            'African'
            if (request.form.get('region_ancestry') or '').strip().lower() == 'african'
            else 'European'
        )
        window = 50000 if browse_type == 'gene' else 100000
        refGen = 'GRCh38'
        # Monitor searches
        add_haplo_to_file(data_path, browse, refGen, browse_type)
        if browse_type == 'trait':
            # get info to plot umap and heatmap
            clusters_of_interest, all_indices, trait_list, cluster_index_map, umap, cluster_representatives, browse_resolved, browse_label, sim_mat, names_sub, cluster_labels_for_index, cluster_traits, traits_interest, idx_trait, all_cluster_representatives, semantic_artifact, selected_trait, combined_context = guide_haplotypes_traits(data_path, browse, window, refGen, trait_names)
            ancestry = selected_trait['ancestry']
            # get table of trait info and gwas associations
            traits_info_table, gwas_assoc_df, all_haplo_df = get_traits_info(data_path, traits_interest, cluster_representatives, meta_snapshot, ancestry=ancestry)
            traits_info_csv = traits_info_table.drop(columns=['__is_representative'], errors='ignore')
            gwas_assoc_df_all = gwas_assoc_df.copy()
            # get plot
            plot_url, gwas_assoc_df, plot_haplo_df = plot_haplotype_traits(clusters_of_interest, all_indices, trait_list, cluster_index_map, umap, cluster_representatives, browse_label, sim_mat, names_sub, cluster_labels_for_index, cluster_traits, idx_trait, all_cluster_representatives, gwas_assoc_df, all_haplo_df, semantic_artifact, selected_trait, combined_context)
            # show only top rows in tables
            MAX_ROWS = 100
            gwas_assoc_df_show = top_gwas_rows(gwas_assoc_df, MAX_ROWS)
            all_haplo_df_show = all_haplo_df.head(MAX_ROWS).copy()
            # store full/subset tsvs in session for download
            session["haplo_trait_full_csv"] = {
                "gwas_assoc_suggestive": gwas_assoc_df_all.to_csv(index=False, sep="\t"),
                "gwas_assoc_established": gwas_assoc_df.to_csv(index=False, sep="\t"),
                "all_haplo_suggestive": all_haplo_df.to_csv(index=False, sep="\t"),
                "all_haplo_established": plot_haplo_df.to_csv(index=False, sep="\t"),
                "traits_info": traits_info_csv.to_csv(index=False, sep="\t"),
            }
            session.modified = True
            # store last state in session
            session['haplo_last'] = {"mode": "trait", "browse": browse_resolved, "browse_label": browse_label, "ancestry": ancestry, "browse_type": browse_type, "refGen": refGen, "window": window, "plot_url": plot_url, "ld_plot_url": None, "fig_url": None, "haplo_summary": [], "chrom": None, "start_pos": None, "end_pos": None, "hap_id_interest": None, "table_traits": traits_info_table.to_dict(orient="records"), "gwas_assoc_table": gwas_assoc_df_show.to_dict(orient="records"), "all_haplo_table": all_haplo_df_show.to_dict(orient="records"), "gwas_assoc_total": int(len(gwas_assoc_df)), "all_haplo_total": int(len(all_haplo_df))}
            return render_template("haplotypes.html", browse_value=browse_resolved, browse_label=browse_label, trait_ancestry=ancestry, traits=trait_list, plot_url=plot_url, ld_plot_url=None, fig_url=None, haplo_summary=[], chrom=None, start_pos=None, end_pos=None, hap_id_interest=None, table_traits=traits_info_table.to_dict(orient="records"), gwas_assoc_table=gwas_assoc_df_show.to_dict(orient="records"), all_haplo_table=all_haplo_df_show.to_dict(orient="records"), gwas_assoc_total=len(gwas_assoc_df), all_haplo_total=len(all_haplo_df))
        elif browse_type == 'chromosome':
            session['haplo_last'] = {"mode": "chromosome", "browse": browse, "browse_type": browse_type, "refGen": refGen, "window": window, "plot_url": None, "ld_plot_url": None, "fig_url": None, "haplo_summary": [], "chrom": None, "start_pos": None, "end_pos": None, "hap_id_interest": None, "table_traits": None, "gwas_assoc_table": None, "all_haplo_table": None, "gwas_assoc_total": None, "all_haplo_total": None}
            return render_template("haplotypes.html", browse_value=browse, traits=trait_list, plot_url=None, fig_url=None, haplo_summary=[])
        else:
            browse_resolved, plot_url, ld_plot_url, haplo_summary, chrom, start_pos, end_pos, hap_id_interest = guide_haplotypes_snps_genes(
                data_path,
                browse,
                window,
                refGen,
                ancestry=region_ancestry,
            )
            haplo_summary_records = haplo_summary.to_dict(orient="records")
            fig_url = extract_haplotype_effect_on_traits(
                data_path,
                haplo_summary,
                chrom,
                ancestry=region_ancestry,
            )
            session['haplo_last'] = {"mode": "region", "browse": browse_resolved, "ancestry": region_ancestry, "browse_type": browse_type, "refGen": refGen, "window": window, "plot_url": plot_url, "ld_plot_url": ld_plot_url, "fig_url": fig_url, "haplo_summary": haplo_summary_records, "chrom": chrom, "start_pos": start_pos, "end_pos": end_pos, "hap_id_interest": hap_id_interest, "table_traits": None, "gwas_assoc_table": None, "all_haplo_table": None, "gwas_assoc_total": None,"all_haplo_total": None}
            return render_template("haplotypes.html", browse_value=browse_resolved, region_browse_type=browse_type, region_ancestry=region_ancestry, traits=trait_list, plot_url=plot_url, ld_plot_url=ld_plot_url, fig_url=fig_url, haplo_summary=haplo_summary_records, chrom=chrom, start_pos=start_pos, end_pos=end_pos, hap_id_interest=hap_id_interest,table_traits=None, gwas_assoc_table=None, all_haplo_table=None)
    # ---- Restore from session ----
    last = session.get('haplo_last')
    if last:
        return render_template("haplotypes.html", traits=trait_list, browse_value=last.get("browse"), browse_label=last.get("browse_label"), trait_ancestry=last.get("ancestry") if last.get('mode') == 'trait' else None, region_ancestry=last.get("ancestry") if last.get('mode') == 'region' else None, region_browse_type=last.get("browse_type") if last.get('mode') == 'region' else None, plot_url=last.get("plot_url"), ld_plot_url=last.get("ld_plot_url"), fig_url=last.get("fig_url"), haplo_summary=last.get("haplo_summary") or [], chrom=last.get("chrom"), start_pos=last.get("start_pos"), end_pos=last.get("end_pos"), hap_id_interest=last.get("hap_id_interest"), table_traits=last.get("table_traits"), gwas_assoc_table=last.get("gwas_assoc_table"), all_haplo_table=last.get("all_haplo_table"), gwas_assoc_total=last.get("gwas_assoc_total"), all_haplo_total=last.get("all_haplo_total"))
    # First load
    return render_template("haplotypes.html", traits=trait_list, plot_url=None, ld_plot_url=None, fig_url=None, haplo_summary=[], chrom=None, start_pos=None, end_pos=None, hap_id_interest=None, table_traits=None, gwas_assoc_table=None, all_haplo_table=None, gwas_assoc_total=None, all_haplo_total=None)

@app.get("/haplotype/download/<which>")
def haplo_download(which):
    blob = session.get("haplo_trait_full_csv", {})
    table_text = blob.get(which)
    if not table_text:
        return ("No table available for download.", 404)
    filename = f"{which}.tsv"
    return Response(table_text, mimetype="text/tab-separated-values", headers={"Content-Disposition": f"attachment; filename={filename}"})



@app.post('/haplotypes/ancestry-overview')
def haplotype_ancestry_overview_plot():
    """Build overview panels for the active region and selected ancestries."""
    last = session.get('haplo_last') or {}
    if last.get('mode') != 'region':
        return jsonify(error='No active regional haplotype query.'), 400

    payload = request.get_json(silent=True)
    if not isinstance(payload, dict):
        return jsonify(error='Invalid overview request.'), 400
    requested = payload.get('ancestries')
    if (not isinstance(requested, list) or not requested or len(requested) > 2
            or any(not isinstance(value, str) or value not in ('European', 'African')
                   for value in requested)
            or len(set(requested)) != len(requested)):
        return jsonify(error='Choose European and/or African ancestry.'), 400
    limit = payload.get('limit', 10)
    if type(limit) is not int or limit not in (10, 20):
        return jsonify(error='Choose a top 10 or top 20 view.'), 400
    chrom = last.get('chrom')
    start_pos = last.get('start_pos')
    end_pos = last.get('end_pos')
    if chrom is None or start_pos is None or end_pos is None:
        return jsonify(error='The active regional query is incomplete.'), 400

    primary = last.get('ancestry', 'European')
    ancestries = [value for value in (primary, 'African' if primary == 'European' else 'European')
                  if value in requested]
    figure = haplotype_ancestry_overview(
        data_path, chrom, int(start_pos), int(end_pos), ancestries,
        refGen=last.get('refGen', 'GRCh38'), limit=limit,
    )
    return Response(pio.to_json(figure), mimetype='application/json')


@app.post("/haplotypes/regional-ld")
def haplotype_regional_ld():
    """Rebuild the regional LD roof for a user-selected haplotype subset."""
    last = session.get('haplo_last') or {}
    if last.get('mode') != 'region':
        return jsonify(error='No active regional haplotype query.'), 400

    payload = request.get_json(silent=True) or {}
    requested_ids = payload.get('haplotypes')
    if not isinstance(requested_ids, list):
        return jsonify(error='Haplotypes must be provided as a list.'), 400

    summary_records = last.get('haplo_summary') or []
    label_by_id = {
        str(record.get('ID')): str(record.get('haplo_index'))
        for record in summary_records
        if record.get('ID') is not None and record.get('haplo_index') is not None
    }
    requested_set = {str(hap_id) for hap_id in requested_ids}
    selected_ids = [hap_id for hap_id in label_by_id if hap_id in requested_set]
    if not selected_ids:
        return jsonify(error='Select at least one available haplotype.'), 400
    if len(selected_ids) != len(requested_set):
        return jsonify(error='One or more selected haplotypes are not available in this region.'), 400

    chrom = last.get('chrom')
    start_pos = last.get('start_pos')
    end_pos = last.get('end_pos')
    if chrom is None or start_pos is None or end_pos is None:
        return jsonify(error='The active regional query is incomplete.'), 400

    ancestry = last.get('ancestry', 'European')
    haplo_df, _ = extract_haplo_data(
        data_path, chrom, start_pos, end_pos, ancestry=ancestry
    )
    selected_map = {hap_id: label_by_id[hap_id] for hap_id in selected_ids}
    selected_df = haplo_df[haplo_df['ID'].astype(str).isin(selected_ids)].copy()
    figure = build_haplotype_ld_roof(
        data_path,
        chrom,
        selected_df,
        selected_map,
        ancestry=ancestry,
    )
    return Response(pio.to_json(figure), mimetype='application/json')

@app.route("/haplotypes/detail")
def haplotype_detail():
    hap_id = request.args.get("hap_id")
    # Derive chromosome and position
    chrom = request.args.get("chrom")
    start_pos = int(request.args.get("start_pos"))
    end_pos = int(request.args.get("end_pos"))
    ancestry = (
        'African'
        if (request.args.get('ancestry') or '').strip().lower() == 'african'
        else 'European'
    )
    refGen = 'GRCh38'
    # get data of the region of interest
    haplo_df, haplo_dict = extract_haplo_data(
        data_path,
        chrom,
        start_pos-1_000_000,
        end_pos+1_000_000,
        ancestry=ancestry,
    )
    haplo_df = haplo_df[haplo_df['ID'] == hap_id]   
    # get snps
    snps = ' '.join([str(chrom).upper().replace('CHR', '') + ':' + x.split(':')[1] + '-' + x.split(':')[1] for x in haplo_dict[hap_id]])
    # get ld between snps
    ld_df = get_ld_between_snps(snps, data_path, chrom, ancestry=ancestry)
    
    # get genes
    genes = extract_genes(data_path, chrom, start_pos-250_000, end_pos+250_000, refGen)
    # gather snp association
    association_file = (
        '20260917_All_indep_gwas_sumstats_AI_hg38_African.txt.gz'
        if ancestry == 'African'
        else 'All_indep_gwas_sumstats_AI_hg38.txt.gz'
    )
    cmd = "tabix %s/databases/haplotypes/%s %s" %(
        data_path.replace(' ', '\ '),
        association_file,
        snps,
    )
    snps_data = [x.rstrip().split('\t') for x in os.popen(cmd)]
    # convert to dataframe
    snps_df = pd.DataFrame(snps_data, columns=["id", "trait", "chr", "position", "rsid", "ea", "nea", "eaf", "beta", "se", "p", "n", "study_id", "position_hg38"])
    # also gather cadd consequences of the snps
    cadd_df = query_cadd_score(snps, data_path)
    cadd_df['pos'] = cadd_df['pos'].astype(int)
    snps_df['position_hg38'] = snps_df['position_hg38'].astype(int)
    cadd_df = cadd_df.merge(snps_df[['rsid', 'position_hg38', 'ea', 'nea']], left_on='pos', right_on='position_hg38', how='outer')
    # restrict to rows where the ref and alt matches with ea and nea
    cadd_df = cadd_df[((cadd_df['ref'] == cadd_df['ea']) & (cadd_df['alt'] == cadd_df['nea'])) | ((cadd_df['ref'] == cadd_df['nea']) & (cadd_df['alt'] == cadd_df['ea']))]

    # create plotly figure
    fig = make_subplots(
        rows=3,
        cols=2,
        shared_xaxes=True,              # still share X for col=1 (rows 2 & 3)
        vertical_spacing=0.08,
        horizontal_spacing=0.08,
        row_heights=[0.40, 0.40, 0.20],
        column_widths=[0.4, 0.6],
        specs=[
            [{"type": "heatmap", "colspan": 2}, None],                    # Row 1: top heatmap spans both columns
            [{"secondary_y": True}, {"type": "heatmap", "rowspan": 2}],   # Row 2: SNPs + trait-trait heatmap (rows 2–3)
            [{},                      None]                               # Row 3: gene track
        ]
    )

    # Overall figure height
    fig.update_layout(height=1700)

    # Grid only for the SNP/gene panels (rows 2 and 3, col 1)
    for r in [2, 3]:
        fig.update_xaxes(showgrid=True, gridcolor="#e0e0e0", zeroline=False, row=r, col=1)
        fig.update_yaxes(showgrid=True, gridcolor="#e0e0e0", zeroline=False, row=r, col=1)

    fig.update_layout(
        plot_bgcolor="white",
        paper_bgcolor="white",
        title_font=dict(size=22),
        font=dict(family="Arial", size=14)
    )

    # ------------------------------------------------------------------
    # (0) NEW TOP-ROW HEATMAP (placeholder)
    # ------------------------------------------------------------------
    # LD df with clustering-based ordering (no dendrogram shown)
    ld_plot = ld_df.copy()
    ld_plot = ld_plot.fillna(0.0)
    # enforce symmetry for stable clustering without shrinking LD values:
    # some pairs are present only in one direction, so use max(A, A^T), not average.
    m = ld_plot.values.astype(float)
    m = np.maximum(m, m.T)
    np.fill_diagonal(m, 1.0)
    ld_plot = pd.DataFrame(m, index=ld_plot.index, columns=ld_plot.columns)
    if ld_plot.shape[0] > 1:
        dist = 1 - ld_plot.values
        np.fill_diagonal(dist, 0.0)
        dist = np.clip(dist, 0.0, 1.0)
        dist_condensed = squareform(dist, checks=False)
        Z = linkage(dist_condensed, method="average")
        order = leaves_list(Z)
        ld_plot = ld_plot.iloc[order, order]
    corr_vals = ld_plot.values
    traits = list(ld_plot.columns)
    custom_ld_colorscale = [
        [0.0,   "navy"],   # 0.0 → base color (≤0.2)
        [0.2,   "navy"],   # stop until 0.2
        [0.20001, "lightblue"], # >0.2
        [0.4,   "lightblue"],   # up to 0.4
        [0.40001, "green"],     # >0.4
        [0.6,   "green"],       # up to 0.6
        [0.60001, "orange"],    # >0.6
        [0.8,   "orange"],      # up to 0.8
        [0.80001, "red"],       # >0.8
        [1.0,   "red"],         # max
    ]
    overall_heatmap = go.Heatmap(
        z=corr_vals,
        x=traits,
        y=traits,
        zmin=0,
        zmax=1,
        colorscale=custom_ld_colorscale,
        colorbar=dict(title="LD R2", orientation="v", x=1.02, y=0.83, len=0.3, thickness=15),
        hovertemplate=(
            "Trait 1: %{y}<br>"  
            "Trait 2: %{x}<br>"  
            "r = %{z:.2f}<extra></extra>"
        ),
        showscale=True,
    )
    fig.add_trace(overall_heatmap, row=1, col=1)
        
    # ------------------------------------------------------------------
    # (1) SNP ASSOCIATIONS (row 2, col 1)
    # ------------------------------------------------------------------
    # count associations per trait
    trait_counts = snps_df['trait'].value_counts().to_dict()
    # palette of 20 colors - tab20
    cmap = plt.get_cmap("tab20", len(trait_counts))
    palette = [mcolors.to_hex(cmap(i)) for i in range(len(trait_counts))]

    # iterate over traits
    for i, (trait, count) in enumerate(trait_counts.items()):
        # subset dataframe
        sub_df = snps_df[snps_df['trait'] == trait]
        # positions in Mb
        positions = sub_df['position'].astype(int) / 1_000_000
        # pvalues in -log10
        pvalues = -np.log10(sub_df['p'].astype(float))
        # create scatter trace
        fig.add_trace(
            go.Scatter(
                x=positions,
                y=pvalues,
                mode='markers',
                marker=dict(
                    color=palette[i % len(palette)],
                    size=6,
                    line=dict(width=1, color="black")
                ),
                name=f'{trait} (n={count})',
                text=[
                    f'Trait: {trait}<br>'
                    f'SNP: {row["rsid"]}<br>'
                    f'Position: {row["chr"]}:{row["position"]}<br>'
                    f'P-value: {row["p"]}<br>'
                    f'Beta: {row["beta"]} (SE: {row["se"]})'
                    for _, row in sub_df.iterrows()
                ],
                hoverinfo='text'
            ),
            row=2,
            col=1
        )

    # update layout for SNP associations panel
    fig.update_layout(
        title=f'Haplotype {hap_id}',
        xaxis_title='',
        yaxis_title='',
        showlegend=True,
        # place legend below the SNP panel
        legend=dict(orientation='h', y=-0.10, x=-0.05, xanchor='left', yanchor='top', itemwidth=30, font=dict(size=9), tracegroupgap=2),
        margin=dict(t=120, b=150),
        # increase overall text sizes
        title_font=dict(size=24, family="Arial", color="black", weight="bold"),
        font=dict(size=14),
        legend_font=dict(size=14)
    )

    # increase marker size for marker traces and enlarge annotation text
    fig.update_traces(marker=dict(size=14), selector=dict(mode='markers'))
    fig.update_annotations(font=dict(size=14))

    # make axis title and tick fonts larger for SNP panel (row 2, col 1)
    fig.update_yaxes(
        title_text="-log10(P-value)",
        row=2,
        col=1,
        title_font=dict(size=20, family="Arial", color="black", weight="bold"),
        tickfont=dict(size=16, family="Arial", color="black"),
    )
    fig.update_xaxes(
        row=2,
        col=1,
        title_font=dict(size=18, family="Arial", color="black", weight="bold"),
        tickfont=dict(size=16, family="Arial", color="black"),
    )

    # ------------------------------------------------------------------
    # (2) GENE TRACK (row 3, col 1)
    # ------------------------------------------------------------------
    if genes:
        index_pos = 4 if refGen == 'GRCh37' else 5
        index_text = 0 if refGen == 'GRCh37' else -2
        index_exon = 9 if refGen == 'GRCh37' else 8
        index_strand = 3 if refGen == 'GRCh37' else 2
        index_transc = 1 if refGen == 'GRCh37' else 0
        index_exonN = 8 if refGen == 'GRCh37' else 7
        max_y = max(g[-1] for g in genes) if len(genes) else 1

        # set parameters for gene names depending on the number of genes
        if len(genes) < 10:
            fsize, yoff, gw = 16, max_y * 0.10, 6
        elif len(genes) < 20:
            fsize, yoff, gw = 10, max_y * 0.06, 4
        else:
            fsize, yoff, gw = 5, max_y * 0.02, 2
        if yoff == 0:
            yoff = 0.10

        for g in genes:
            # create 2d array for hovertemplate
            row = np.array([
                str(g[index_text]),
                str(g[index_transc]),
                str(g[index_strand]),
                int(g[index_pos]),
                int(g[index_pos + 1]),
                int(g[index_exonN])
            ], dtype=object)
            customdata = np.vstack([row, row])   # shape (2, 6)
            x0 = int(g[index_pos])   / 1_000_000
            x1 = int(g[index_pos+1]) / 1_000_000
            ylane = int(g[-1])
            color = 'grey' if g[index_strand] == '+' else 'black'
            fig.add_trace(
                go.Scatter(
                    x=[x0, x1], y=[ylane, ylane],
                    mode='lines',
                    hoverinfo='skip',
                    line=dict(color=normalize_color(color), width=gw),
                    showlegend=False
                ),
                row=3,
                col=1
            )
            fig.update_layout(hovermode='closest')

            # add arrow for strand
            right = (color == 'grey')
            symbol = 'triangle-right' if right else 'triangle-left'
            x_tip = x1 if right else x0
            fig.add_trace(
                go.Scatter(
                    x=[x_tip], y=[ylane],
                    mode='markers',
                    marker=dict(symbol=symbol, size=10, color=color, line=dict(width=0)),
                    hoverinfo='skip', showlegend=False
                ),
                row=3,
                col=1
            )

            # add gene name and hovering information
            label = f"<b><i>{g[index_text]}</i></b>"
            hovertext = (
                f"<b>{g[index_text]}</b><br>"
                f"Transcript: {g[index_transc]}<br>"
                f"Strand: {g[index_strand]}<br>"
                f"TxStart: {int(g[index_pos])}<br>"
                f"TxEnd: {int(g[index_pos+1])}<br>"
                f"Exons: {int(g[index_exonN])}"
            )
            fig.add_annotation(
                x=x0 + (x1 - x0) / 2,
                y=ylane + yoff,
                text=label,
                showarrow=False,
                font=dict(size=fsize),
                hovertext=hovertext,
                row=3,
                col=1,
            )

        fig.update_yaxes(range=[-max_y*0.20, max_y+max_y*0.20], row=3, col=1)
    else:
        fig.add_annotation(
            x=(start_pos + (end_pos-start_pos)/2)/1_000_000,
            y=0.5,
            text="No Genes to plot",
            showarrow=False,
            row=3,
            col=1,
        )

    # increase gene-track axis fonts (row 3, col 1)
    fig.update_yaxes(
        title_text="Gene track",
        row=3,
        col=1,
        title_font=dict(size=20, family="Arial", color="black", weight="bold"),
        tickfont=dict(size=16, family="Arial", color="black"),
    )
    fig.update_xaxes(
        row=3,
        col=1,
        title_font=dict(size=18, family="Arial", color="black", weight="bold"),
        tickfont=dict(size=16, family="Arial", color="black"),
    )

    # ------------------------------------------------------------------
    # (3) X-RANGE FOR THE REGION (rows 2 & 3, col 1)
    # ------------------------------------------------------------------
    min_x_pos = snps_df['position'].astype(int).min() / 1_000_000
    max_x_pos = snps_df['position'].astype(int).max() / 1_000_000
    region_size = max_x_pos - min_x_pos
    # extend by 40% on each side
    min_x_pos -= region_size * 0.40
    max_x_pos += region_size * 0.40

    fig.update_xaxes(range=[min_x_pos, max_x_pos], row=2, col=1)
    fig.update_xaxes(range=[min_x_pos, max_x_pos], row=3, col=1)

    # ------------------------------------------------------------------
    # (4) TRAIT-TRAIT CORRELATION HEATMAP (RIGHT PANEL, rows 2-3, col 2)
    # ------------------------------------------------------------------
    # Correlate trait-effect patterns across SNPs. Missing SNP-trait
    # associations remain missing and are excluded pairwise rather than being
    # interpreted as zero effects.
    snps_df['beta'] = snps_df['beta'].astype(float)
    snps_df['p'] = snps_df['p'].astype(float)
    snps_df = (snps_df.loc[snps_df.groupby(["rsid", "trait"])["p"].idxmin()].reset_index(drop=True)).copy()
    beta_matrix = snps_df.pivot(index='rsid', columns='trait', values='beta')
    min_overlapping_snps = 3
    pearson_corr = beta_matrix.corr(method='pearson', min_periods=min_overlapping_snps)
    spearman_corr = beta_matrix.corr(method='spearman', min_periods=min_overlapping_snps)
    overlap_counts = beta_matrix.notna().astype(int).T.dot(beta_matrix.notna().astype(int))
    traits = list(beta_matrix.columns)

    # Establish one fixed ordering from the Pearson matrix so switching the
    # correlation method does not move traits around. Missing correlations are
    # neutral only for this ordering step and stay missing in both heatmaps.
    if len(traits) > 1:
        corr_for_ordering = pearson_corr.to_numpy(dtype=float)
        corr_for_ordering = np.nan_to_num(corr_for_ordering, nan=0.0)
        corr_for_ordering = (corr_for_ordering + corr_for_ordering.T) / 2
        corr_for_ordering = np.clip(corr_for_ordering, -1.0, 1.0)
        np.fill_diagonal(corr_for_ordering, 1.0)
        dist = 1 - corr_for_ordering
        np.fill_diagonal(dist, 0.0)
        dist_condensed = squareform(dist, checks=False)
        Z = linkage(dist_condensed, method="average")
        order = leaves_list(Z)
    else:
        order = np.arange(len(traits))

    pearson_clustered = pearson_corr.to_numpy(dtype=float)[np.ix_(order, order)]
    spearman_clustered = spearman_corr.to_numpy(dtype=float)[np.ix_(order, order)]
    overlap_clustered = overlap_counts.to_numpy(dtype=int)[np.ix_(order, order)]
    traits_clustered = [traits[i] for i in order]
    traits_wrapped = [t[:30] for t in traits_clustered]

    def json_safe_correlation_matrix(matrix):
        return [
            [float(value) if np.isfinite(value) else None for value in row]
            for row in matrix
        ]

    pearson_plot_values = json_safe_correlation_matrix(pearson_clustered)
    spearman_plot_values = json_safe_correlation_matrix(spearman_clustered)
    pearson_na_text = np.where(np.isfinite(pearson_clustered), "", "/").tolist()
    spearman_na_text = np.where(np.isfinite(spearman_clustered), "", "/").tolist()

    # define palette
    custom_corr_colorscale = [[0.0,  "#4B0082"], [0.5,  "#FFFFFF"], [1.0,  "#B22222"]]

    # Use separate traces for the two methods. The buttons only change the
    # visibility of these traces, keeping their payloads JSON-safe and leaving
    # every other subplot trace untouched.
    pearson_heatmap = go.Heatmap(
        z=pearson_plot_values,
        x=traits_wrapped,
        y=traits_wrapped,
        customdata=overlap_clustered.tolist(),
        text=pearson_na_text,
        texttemplate="%{text}",
        textfont=dict(color="darkgrey", size=18, family="Arial Black"),
        colorscale=custom_corr_colorscale,
        zmin=-1,
        zmax=1,
        zmid=0,
        colorbar=dict(title="Pearson r", orientation="v", x=1.02, y=0.30, len=0.3, thickness=15),
        hovertemplate=(
            "Trait 1: %{y}<br>"
            "Trait 2: %{x}<br>"
            "Pearson r: %{z:.2f}<br>"
            "Overlapping SNPs: %{customdata}<extra></extra>"
        ),
        visible=True,
    )
    pearson_trace_index = len(fig.data)
    fig.add_trace(pearson_heatmap, row=2, col=2)

    spearman_heatmap = go.Heatmap(
        z=spearman_plot_values,
        x=traits_wrapped,
        y=traits_wrapped,
        customdata=overlap_clustered.tolist(),
        text=spearman_na_text,
        texttemplate="%{text}",
        textfont=dict(color="darkgrey", size=18, family="Arial Black"),
        colorscale=custom_corr_colorscale,
        zmin=-1,
        zmax=1,
        zmid=0,
        colorbar=dict(title="Spearman ρ", orientation="v", x=1.02, y=0.30, len=0.3, thickness=15),
        hovertemplate=(
            "Trait 1: %{y}<br>"
            "Trait 2: %{x}<br>"
            "Spearman ρ: %{z:.2f}<br>"
            "Overlapping SNPs: %{customdata}<extra></extra>"
        ),
        visible=False,
    )
    spearman_trace_index = len(fig.data)
    fig.add_trace(spearman_heatmap, row=2, col=2)

    fig.update_layout(
        updatemenus=[dict(
            type="buttons",
            direction="right",
            active=0,
            x=1.0,
            xanchor="right",
            y=0.54,
            yanchor="bottom",
            buttons=[
                dict(
                    label="Pearson",
                    method="restyle",
                    args=[{"visible": [True, False]}, [pearson_trace_index, spearman_trace_index]],
                ),
                dict(
                    label="Spearman",
                    method="restyle",
                    args=[{"visible": [False, True]}, [pearson_trace_index, spearman_trace_index]],
                ),
            ],
        )]
    )

    # shrink the heatmap panel inside column 2 (rows 2-3)
    fig.update_xaxes(row=2, col=2, domain=[0.48, 1.00], tickangle=45, title_text="", tickfont=dict(size=10))
    fig.update_yaxes(row=2, col=2, autorange="reversed", title_text="", tickfont=dict(size=10))

    # Customize layout for the heatmap
    # --- Row 1 (both columns) ---
    fig.update_yaxes(domain=[0.60, 1.00], row=1, col=1)
    fig.update_yaxes(domain=[0.60, 1.00], row=1, col=2)
    # --- Row 2 (both columns; col 2 spans rows 2–3 but we control yaxis2 domain) ---
    fig.update_yaxes(domain=[0.18, 0.50], row=2, col=1)
    fig.update_yaxes(domain=[0.00, 0.50], row=2, col=2)
    # --- Row 3 (only col 1) ---
    fig.update_yaxes(domain=[0.00, 0.18], row=3, col=1)
    
    # Subtitles
    fig.add_annotation(
        text="<b></b>",
        x=0.5, y=1.05,
        xref="paper", yref="paper",
        showarrow=False,
        font=dict(size=18, color="darkgrey")
    )
    fig.add_annotation(
        text="<b></b>",
        x=0.15, y=0.55,
        xref="paper", yref="paper",
        showarrow=False,
        font=dict(size=18, color="darkgrey")
    )
    fig.add_annotation(
        text="<b></b>",
        x=0.75, y=0.55,
        xref="paper", yref="paper",
        showarrow=False,
        font=dict(size=18, color="darkgrey")
    )    
    
    # Haplotype table
    haplo_df = haplo_df[['ID', 'BP1', 'BP2', 'KB', 'NSNPS', 'assoc_count', 'n_traits', 'rsids', 'traits']].copy()
    haplo_df.rename(columns={'ID':'Haplotype ID', 'BP1':'Start (hg38)', 'BP2':'End (hg38)', 'KB':'Size (kb)', 'NSNPS':'# SNPs', 'assoc_count':'# Associations', 'n_traits':'# Traits', 'rsids':'SNPs (RsIDs)', 'traits':'Traits'}, inplace=True)
    # return tables as html
    haplo_table_html = haplo_df.to_html(
        classes="detail-table",
        index=False,
        border=0,
        escape=False
    )
    
    # SNPs table
    snps_df = snps_df[['rsid', 'position', 'position_hg38', 'trait', 'ea', 'nea', 'eaf', 'beta', 'se', 'p', 'n', 'id']].copy()
    snps_df.rename(columns={'position':'Pos (hg19)', 'position_hg38':'Pos (hg38)', 'rsid':'RsID', 'trait':'Trait', 'ea':'Effect', 'nea':'Other', 'eaf':'EAF', 'beta':'Beta', 'se':'SE', 'p':'P-value', 'n':'N', 'id':'Study ID'}, inplace=True)
    # sort by pvalue
    snps_df['P-value'] = snps_df['P-value'].astype(float)
    snps_df = snps_df.sort_values(by='P-value', ascending=True)
    snps_table_html = snps_df.to_html(
        classes="detail-table",
        index=False,
        border=0,
        escape=False
    )
    
    # CADD table
    # reorder columns
    cadd_df = cadd_df[['rsid', 'pos', 'position_hg38', 'ref', 'alt', 'phred_max', 'genes', 'annotypes', 'consequences', 'cpg_max', 'gc_max', 'h3k27ac_max', 'h3k4me3_max', 'dnase_max', 'spliceai_acc_loss_max', 'spliceai_don_loss_max', 'siftval_max', 'polyphenval_max', 'gerprs_max', 'priphyloP_max']].copy()
    cadd_df.rename(columns={'pos':'Pos (hg19)', 'position_hg38':'Pos (hg38)', 'rsid':'RsID', 'ref':'Ref', 'alt':'Alt', 'phred_max':'Phred score', 'genes':'Gene', 'annotypes':'Annotypes', 'consequences':'Consequences', 'cpg_max':'CpG', 'gc_max':'GC', 'h3k27ac_max':'H3K27ac', 'h3k4me3_max':'H3K4me3', 'dnase_max':'DNase', 'spliceai_acc_loss_max':'SpliceAI (acceptor)', 'spliceai_don_loss_max':'SpliceAI (donor)', 'siftval_max':'SIFT', 'polyphenval_max':'PolyPhen', 'gerprs_max':'GERP', 'priphyloP_max':'PhyloP'}, inplace=True)
    cadd_df = cadd_df.drop_duplicates()
    # sort by phred score
    cadd_df['Phred score'] = cadd_df['Phred score'].astype(float)
    cadd_df = cadd_df.sort_values(by='Phred score', ascending=False)
    snps_cadd_table_html = cadd_df.to_html(
        classes="detail-table",
        index=False,
        border=0,
        escape=False
    )
    
    # Remember which detail was last opened
    session['haplo_last_detail'] = {"hap_id": hap_id, "chrom": chrom, "start_pos": start_pos, "end_pos": end_pos, "ancestry": ancestry}

    # Return as Plotly JSON dict
    return jsonify(fig=fig.to_dict(), haplo_table=haplo_table_html, snps_table=snps_table_html, snps_cadd_table=snps_cadd_table_html)

if __name__ == "__main__":
    app.run()
