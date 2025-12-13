from flask import Flask, request, render_template, send_file, session, make_response, redirect, url_for, jsonify, Response, stream_with_context, Blueprint
from queue import Queue, Empty
from uuid import uuid4
import json
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
from pathlib import Path
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
from ieugwaspy import api, query
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

# Import all functions from snpxplorer_functions.py
from snpxplorer_functions import *

# differences between local and server
# local
#data_path = '/Users/nicco/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/snpXplorer/Data'
# server
data_path = '/Data'
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

# Function to manage the publishing of messages
def publish(msg, console_id="default"):
    q = CONSOLE_QUEUES.get(console_id)
    if q:
        q.put_nowait(str(msg))

# Read IEU Open GWAS database
meta = pd.DataFrame(query.gwasinfo())
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
# Build a search row once to allow fast searching
build_search_row(meta)

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
    q = (request.args.get("q") or "").strip().lower()
    page = int(request.args.get("page", 1))
    page_size = 50
    if not q:
        return jsonify({"results": [], "pagination": {"more": False}})
    # build _search row if missing
    if "_search" not in meta.index:
        build_search_row(meta)
    mask = meta.loc["_search"].str.contains(q, na=False, regex=False)
    cols = meta.columns[mask]
    total = len(cols)
    start = (page - 1) * page_size
    end = min(start + page_size, total)
    results = []
    for cid in cols[start:end]:
        trait = str(meta.at["trait", cid]) if "trait" in meta.index else ""
        author = str(meta.at["author", cid]) if "author" in meta.index else ""
        year = str(meta.at["year", cid]) if "year" in meta.index else ""
        label = f"{trait} — {author} ({year}) [{cid}]".strip()
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
    # Read traits (you can optimize by loading once globally if you want)
    trait_names = pd.read_csv(f"{data_path}/databases/haplotypes/trait_names_and_ids.txt", sep="\t")
    trait_list = trait_names['name'].dropna().astype(str).tolist()
    q_lower = q.lower()
    matches = [t for t in trait_list if q_lower in t.lower()]
    start = (page - 1) * page_size
    end = start + page_size
    slice_ = matches[start:end]
    more = end < len(matches)
    results = [{"id": t, "text": t} for t in slice_]
    return jsonify({"results": results, "pagination": {"more": more}})

# Function to run the pipeline in a separate thread
# Global, in-process result cache for SSE runs
RESULTS = {}  # { console_id: {template_context...} }
def run_pipeline(console_id, formdata):
    """Background runner that mirrors your current POST body and
    stores the final render context in RESULTS[console_id]."""
    try:
        # --- read inputs (same as your POST branch) ---
        refGen = formdata.get("refGenome")
        browse = formdata.get("browse", "")
        gwas = formdata.getlist('gwas_data')
        window = int(formdata.get("window", 25000) or 25000)
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

        # --- your original steps, with logs ---
        publish("Finding genomic location...", console_id)
        chrom, start_pos, end_pos, browse_type = readBrowseOption(data_path, browse, window, refGen)
        logging.info("All good with readBrowseOption")

        publish("Gathering data of interest...", console_id)
        df, df_info = get_data_plot(data_path=data_path, gwas=gwas, chrom=chrom, start_pos=start_pos, end_pos=end_pos, refGen=refGen, meta=meta, window=window)
        genes = extract_genes(data_path, chrom, start_pos, end_pos, refGen)
        svs, svs_df = extract_sv(data_path, chrom, start_pos, end_pos, refGen, svtypes)
        recomb_data = extract_recomb(data_path, chrom, start_pos, end_pos, refGen) if recomb == 'Yes' else "None"
        gtex_df = get_gtex(data_path, genes, refGen)
        ld_df = extract_ld(data_path, chrom, df, refGen, browse_type, browse) if ld == "Yes" else "None"
        haplo_df = extract_haplo(data_path, df, chrom, refGen) if plotHaplo == "Yes" else "None"
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
        df_info_sub = df_info[['trait', 'id', 'pmid', 'author', 'year', 'sample_size', 'build']].drop_duplicates()
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
            window=session.get("window", 25000),
            axestype=session.get("axestype", "manh"),
            plotype=session.get("plotype", "Scatter"),
            exons=session.get("exons", "No"),
            recomb=session.get("recomb", "No"),
            ld=session.get("ld", "No"),
            plotHaplo=session.get("plotHaplo", "No"),
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
            window=session.get("window", 25000),
            axestype=session.get("axestype", "manh"),
            plotype=session.get("plotype", "Scatter"),
            exons=session.get("exons", "No"),
            recomb=session.get("recomb", "No"),
            ld=session.get("ld", "No"),
            plotHaplo=session.get("plotHaplo", "No"),
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
    if request.method == "POST":
        # Read inputs
        # list of variants
        textarea_input = request.form.get('SNPlist', '')
        my_list = textarea_input.replace(' ', '').replace(';', '\n').replace(',', '\n').split('\n')
        # input type
        #inpType = request.form["inputType"]
        inpType = 'rsid'  # for now, only rsid
        # if the input type is not rsid, then read also the reference genome
        refGeno = request.form["annotRefGen"]
        # analysis type
        analType = request.form["analType"]
        # if the analysis type is GSEA, then look at the enrichment sources
        gsea_source = request.form.getlist("EnrichSource") if analType == 'GSEA' else []
        # QTL tissues
        qtl_tissues = request.form.getlist('QTLtissuesAnnot')
        # update message
        if emailResponse == True:
            messageSubmission = 'Your job has been submitted! Check your email shortly.'
            # run the script here
            command = run_annotation(my_list, inpType, refGeno, analType, gsea_source, qtl_tissues, email)
            print(command)
            command_run = subprocess.Popen(command, shell=True)
        else:
            messageSubmission = 'Email is not correct. Please check!'
    return render_template('annotation.html', validEmail=emailResponse, messageSubmission=messageSubmission) 

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
    filepath = os.path.join("/Annotation/RUNS", filename)
    return send_file(filepath, as_attachment=True, download_name=filename)

# Haplotype tab
@app.route('/haplotype/', methods=["GET", "POST"])
def haplotypes():
    trait_names = pd.read_csv(f"{data_path}/databases/haplotypes/trait_names_and_ids.txt", sep="\t")
    trait_list = trait_names['name'].values.tolist()
    if request.method == "POST":
        # whenever the user submits the form, clear previous session data
        session.pop('haplo_last_detail', None)
        # get query
        browse = (request.form.get("browse") or "").strip()
        browse_type = request.form.get("browse_type")
        window = 50000 if browse_type == 'gene' else 100000
        refGen = 'GRCh38'
        if browse_type == 'trait':
            # get info to plot umap and heatmap
            clusters_of_interest, all_indices, trait_list, cluster_index_map, umap, cluster_representatives, browse_resolved, sim_mat, names_sub, cluster_labels_for_index, cluster_traits, traits_interest, idx_trait, all_cluster_representatives = guide_haplotypes_traits(data_path, browse, window, refGen, trait_names)
            # get table of traits info
            traits_info_table, gwas_assoc_df, all_haplo_df = get_traits_info(data_path, traits_interest, cluster_representatives, meta)
            # get plot
            plot_url, gwas_assoc_df = plot_haplotype_traits(clusters_of_interest, all_indices, trait_list, cluster_index_map, umap, cluster_representatives, browse, sim_mat, names_sub, cluster_labels_for_index, cluster_traits, idx_trait, all_cluster_representatives, gwas_assoc_df, all_haplo_df)
            
            # persist in session for download details later
            session['haplo_last'] = {"mode": "trait", "browse": browse_resolved, "browse_type": browse_type, "refGen": refGen, "window": window, "plot_url": plot_url, "haplo_summary": [], "chrom": None, "start_pos": None, "end_pos": None, "hap_id_interest": None, "table_traits": traits_info_table.to_dict(orient="records"), "gwas_assoc_table": gwas_assoc_df.to_dict(orient="records"), "all_haplo_table": all_haplo_df.to_dict(orient="records")}
                        
            return render_template("haplotypes.html", browse_value=browse_resolved, traits=trait_list, plot_url=plot_url, haplo_summary=[], chrom=None, start_pos=None, end_pos=None, hap_id_interest=None, table_traits=traits_info_table.to_dict(orient="records"), gwas_assoc_table=gwas_assoc_df.to_dict(orient="records"), all_haplo_table=all_haplo_df.to_dict(orient="records"))
        elif browse_type == 'chromosome':
            # nothing heavy but do remember the browse for later
            session['haplo_last'] = {"mode": "chromosome", "browse": browse, "browse_type": browse_type, "refGen": refGen, "window": window, "plot_url": None, "haplo_summary": [], "chrom": None, "start_pos": None, "end_pos": None, "hap_id_interest": None, "table_traits": None, "gwas_assoc_table": None, "all_haplo_table": None}
            
            return render_template("haplotypes.html", browse_value=browse, traits=trait_list, plot_url=None, haplo_summary=[])
        else:
            # get data and plots
            browse_resolved, plot_url, haplo_summary, chrom, start_pos, end_pos, hap_id_interest = guide_haplotypes_snps_genes(data_path, browse, window, refGen)
            haplo_summary_records=haplo_summary.to_dict(orient="records")
            
            # persist in session for download details later
            session['haplo_last'] = {"mode": "region", "browse": browse_resolved, "browse_type": browse_type, "refGen": refGen, "window": window, "plot_url": plot_url, "haplo_summary": haplo_summary_records,"chrom": chrom, "start_pos": start_pos, "end_pos": end_pos, "hap_id_interest": hap_id_interest, "table_traits": None, "gwas_assoc_table": None, "all_haplo_table": None}
            
            # return render
            return render_template("haplotypes.html", browse_value=browse_resolved, traits=trait_list, plot_url=plot_url, haplo_summary=haplo_summary_records, chrom=chrom, start_pos=start_pos, end_pos=end_pos, hap_id_interest=hap_id_interest, table_traits=None, gwas_assoc_table=None, all_haplo_table=None)
    
    # Get restored from session if available
    last = session.get('haplo_last')
    if last:
        return render_template("haplotypes.html", traits=trait_list, browse_value=last.get("browse"), plot_url=last.get("plot_url"), haplo_summary=last.get("haplo_summary") or [], chrom=last.get("chrom"), start_pos=last.get("start_pos"), end_pos=last.get("end_pos"), hap_id_interest=last.get("hap_id_interest"), table_traits=last.get("table_traits"), gwas_assoc_table=last.get("gwas_assoc_table"), all_haplo_table=last.get("all_haplo_table"))

    # First load
    return render_template("haplotypes.html", traits=trait_list, plot_url=None, haplo_summary=[], chrom=None, start_pos=None, end_pos=None, hap_id_interest=None, table_traits=None, gwas_assoc_table=None, all_haplo_table=None)

@app.route("/haplotypes/detail")
def haplotype_detail():
    hap_id = request.args.get("hap_id")
    # Derive chromosome and position
    chrom = request.args.get("chrom")
    start_pos = int(request.args.get("start_pos"))
    end_pos = int(request.args.get("end_pos"))
    refGen = 'GRCh38'
    # get data of the region of interest
    haplo_df, haplo_dict = extract_haplo_data(data_path, chrom, start_pos-1_000_000, end_pos+1_000_000)
    haplo_df = haplo_df[haplo_df['ID'] == hap_id]   
    # get snps
    snps = ' '.join([str(chrom).upper().replace('CHR', '') + ':' + x.split(':')[1] + '-' + x.split(':')[1] for x in haplo_dict[hap_id]])
    # get ld between snps
    ld_df = get_ld_between_snps(snps, data_path, chrom)
    
    # get genes
    genes = extract_genes(data_path, chrom, start_pos-250_000, end_pos+250_000, refGen)
    # gather snp association
    cmd = "tabix %s/databases/haplotypes/All_indep_gwas_sumstats_AI_hg38.txt.gz %s" %(data_path.replace(' ', '\ '), snps)
    snps_data = [x.rstrip().split('\t') for x in os.popen(cmd)]
    # convert to dataframe
    snps_df = pd.DataFrame(snps_data, columns=["id", "trait", "chr", "position", "rsid", "ea", "nea", "eaf", "beta", "se", "p", "n", "study_id", "position_hg38"])
    # also gather cadd consequences of the snps
    cadd_df = query_cadd_score(snps, data_path)
    cadd_df['pos'] = cadd_df['pos'].astype(int)
    snps_df['position_hg38'] = snps_df['position_hg38'].astype(int)
    cadd_df = cadd_df.merge(snps_df[['rsid', 'position_hg38']], left_on='pos', right_on='position_hg38', how='outer')

    # create plotly figure
    fig = fig = make_subplots(
        rows=3,
        cols=2,
        shared_xaxes=True,              # still share X for col=1 (rows 2 & 3)
        vertical_spacing=0.04,
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
    # LD df
    corr_vals = ld_df.values
    traits = list(ld_df.columns)
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
        colorbar=dict(title="r", orientation="v", x=1.02, y=0.83, len=0.3, thickness=15),
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
                    size=10,
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
    # correlation of the beta values between traits
    snps_df['beta'] = snps_df['beta'].astype(float)
    snps_df['p'] = snps_df['p'].astype(float)
    snps_df = (snps_df.loc[snps_df.groupby(["rsid", "trait"])["p"].idxmin()].reset_index(drop=True)).copy()
    beta_matrix = snps_df.pivot(index='rsid', columns='trait', values='beta').fillna(0)
    beta_corr = beta_matrix.corr().round(2)
    if isinstance(beta_corr, pd.DataFrame):
        corr_vals = beta_corr.values
        traits = list(beta_corr.columns)
    else:
        corr_vals = beta_corr
        traits = [f"T{i+1}" for i in range(corr_vals.shape[0])]
    # clustering
    if len(corr_vals) >1:
        dist = 1 - corr_vals
        # condensed form for linkage()
        dist_condensed = squareform(dist, checks=False)
        Z = linkage(dist_condensed, method="average")
        # Order of items according to the dendrogram leaves
        order = leaves_list(Z)
        corr_clustered = corr_vals[np.ix_(order, order)]
        # original trait labels (before wrapping)
        traits_clustered = [traits[i] for i in order]
        # shorten names and wrap long trait names
        traits_wrapped = [t[:30] for t in traits_clustered]
    else:
        corr_clustered = corr_vals
        traits_wrapped = [t[:30] for t in traits]

    # define palette
    custom_corr_colorscale = [[0.0,  "#4B0082"], [0.5,  "#FFFFFF"], [1.0,  "#B22222"]]

    # create heatmap
    heatmap = go.Heatmap(
        z=corr_clustered,
        x=traits_wrapped,
        y=traits_wrapped,
        colorscale=custom_corr_colorscale,
        zmin=-1,
        zmax=1,
        zmid=0,
        colorbar=dict(title="r", orientation="v", x=1.02, y=0.30, len=0.3, thickness=15),
        hovertemplate=(
            "Trait 1: %{y}<br>"
            "Trait 2: %{x}<br>"
            "r = %{z:.2f}<extra></extra>"
        ),
    )
    fig.add_trace(heatmap, row=2, col=2)

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
        text="<b>LD Heatmap</b>",
        x=0.5, y=1.05,
        xref="paper", yref="paper",
        showarrow=False,
        font=dict(size=18, color="darkgrey")
    )
    fig.add_annotation(
        text="<b>SNP Associations</b>",
        x=0.15, y=0.55,
        xref="paper", yref="paper",
        showarrow=False,
        font=dict(size=18, color="darkgrey")
    )
    fig.add_annotation(
        text="<b>Trait–Trait Correlation</b>",
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
    snps_cadd_table_html = cadd_df.to_html(
        classes="detail-table",
        index=False,
        border=0,
        escape=False
    )
    
    # Remember which detail was last opened
    session['haplo_last_detail'] = {"hap_id": hap_id, "chrom": chrom, "start_pos": start_pos, "end_pos": end_pos}

    # Return as Plotly JSON dict
    return jsonify(fig=fig.to_dict(), haplo_table=haplo_table_html, snps_table=snps_table_html, snps_cadd_table=snps_cadd_table_html)

if __name__ == "__main__":
    app.run()