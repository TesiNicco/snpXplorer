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

# differences between local and server
# local
#from pandasgwas.get_SNPs import get_variants_by_variant_id
#data_path = '/Users/nicco/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/snpXplorer/Data'
# server
data_path = '/Data'
from pandasgwas import get_variants_by_variant_id

# Initialize the App
app = Flask(__name__)
# Register Blueprint
app.register_blueprint(snpbot_bp, url_prefix='/snpbot')
# Add configuration for the SQLAlchemy database (for usernames and passwords)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///db.sqlite'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)
# Add configuration for the Session to save to redis
app.config['SESSION_TYPE'] = 'redis'
app.config['SESSION_REDIS'] = redis.from_url('redis://localhost:6379')
app.config['SESSION_USE_SIGNER'] = True
app.config['SECRET_KEY'] = 'secret'
Session(app)
# Initialize Queues
CONSOLE_QUEUES = {}
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
gene_col = genes_map.columns[0]
genes_lower = genes_map[gene_col].astype(str).str.lower()

# Helper for the SQlite database
rsid_db_path = '%s/databases/Genes/variants_info.sqlite' %(Path(data_path))
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
            SELECT rsid, chr_hg38, position
            FROM rsids
            WHERE rsid LIKE ?
            ORDER BY (rsid = ?) DESC,     -- exact match first
                     LENGTH(rsid),        -- shorter (closer) IDs next (e.g. rs7412 before rs7412123)
                     rsid
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

        publish("Gathering data of interest...", console_id)
        print(gwas)
        df, df_info = get_data_plot(data_path=data_path, gwas=gwas, chrom=chrom, start_pos=start_pos, end_pos=end_pos, refGen=refGen)
        genes = extract_genes(data_path, chrom, start_pos, end_pos, refGen)
        svs, svs_df = extract_sv(data_path, chrom, start_pos, end_pos, refGen, svtypes)
        recomb_data = extract_recomb(data_path, chrom, start_pos, end_pos, refGen) if recomb == 'Yes' else "None"
        gtex_df = get_gtex(data_path, genes, refGen)
        ld_df = extract_ld(data_path, chrom, df, refGen, browse_type, browse) if ld == "Yes" else "None"
        haplo_df = extract_haplo(data_path, df, chrom, refGen) if plotHaplo == "Yes" else "None"

        publish("Rendering plots...", console_id)
        fig = scatterplot_plotly(df=df, chrom=chrom, start_pos=start_pos, end_pos=end_pos, gwas=gwas, genes=genes, svs=svs, browse_type=browse_type, refGen=refGen, recomb_data=recomb_data, exons=exons, browse=browse, plotype=plotype, ld=ld, ld_df=ld_df, yrange=yrange, haplo_df=haplo_df, plotHaplo=plotHaplo)
        fig.update_layout(autosize=True, margin=dict(l=40, r=20, t=90, b=40))
        plot_url = pio.to_html(fig, include_plotlyjs='cdn', full_html=False, default_width='100%', default_height='100%', config={'responsive': True})

        img_gtex = gtex_heatmap(gtex_df)
        plot_gtex = pio.to_html(img_gtex, include_plotlyjs='cdn', full_html=False, default_width='100%', default_height='100%', config={'responsive': True})

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
        # GWAS table
        df_info_sub = df_info[['trait', 'id', 'pmid', 'author', 'year', 'sample_size', 'build']].drop_duplicates()
        df_info_sub = df_info_sub.rename(columns={'trait':'Trait', 'id':'ID', 'pmid':'PMID', 'author':'Author', 'year':'Year', 'sample_size':'Sample Size', 'build':'Build'})
        df_info_table = df_info_sub.to_dict(orient='records')
        # SVs table
        svs_df = svs_df.copy()
        svs_df['Region'] = svs_df.apply(lambda x: f"{x['chrom']}:{x['start']}-{x['end']}", axis=1)
        svs_df = svs_df[['Region', 'len', 'repName', 'repClass', 'repFamily']]
        svs_table = svs_df.to_dict(orient='records')
        # GWAS Catalog
        gwas_table, gwascat_df = extract_gwascatalog(data_path, chrom, start_pos, end_pos, refGen)
        # Persist for download routes if you still want them after reload
        # (we can’t touch flask.session in a thread, so stash CSVs here too)
        ctx = {
            "plot_url": plot_url,
            "plot_gtex": plot_gtex,
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
        publish(f"[error] {type(e).__name__}: {e}", console_id)

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
        fig = scatterplot_plotly(df=df, chrom=chrom, start_pos=start_pos, end_pos=end_pos, gwas=gwas, genes=genes, svs=svs, browse_type=browse_type, refGen=refGen, recomb_data=recomb_data, exons=exons, browse=browse, plotype=plotype, ld=ld, ld_df=ld_df, yrange=yrange, haplo_df=haplo_df, plotHaplo=plotHaplo)
        fig.update_layout(autosize=True, margin=dict(l=40, r=20, t=90, b=40))
        plot_url = pio.to_html(fig, include_plotlyjs='cdn', full_html=False, default_width='100%', default_height='100%', config={'responsive': True})

        img_gtex = gtex_heatmap(gtex_df)
        plot_gtex = pio.to_html(img_gtex, include_plotlyjs='cdn', full_html=False, default_width='100%', default_height='100%', config={'responsive': True})

        return render_template(
            "exploration.html",
            plot_url=plot_url, plot_gtex=plot_gtex,
            table_snps=snps_table, table_svs=svs_table, table_gwas=gwas_table, table_info=info_table,
            browse_value=browse, gwas=gwas,
            console_id=str(uuid4()),
            refGenome=session.get("refGen", "GRCh37"),
            window=session.get("window", 25000),
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
        inpType = request.form["inputType"]
        # if the input type is not rsid, then read also the reference genome
        refGeno = request.form["annotRefGen"] if inpType != 'RsID' else 'GRCh37'
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
    filename = f"AnnotateMe_results_{run_id}.tar.gz"
    filepath = os.path.join("/Annotation/RUNS", filename)
    return send_file(filepath, as_attachment=True, download_name=filename)

# Haplotype tab
@app.route('/haplotype/', methods=["GET", "POST"])
def haplotypes():
    trait_names = pd.read_csv(f"{data_path}/databases/haplotypes/trait_names_and_ids.txt", sep="\t")
    trait_list = trait_names['name'].values.tolist()
    if request.method == "POST":
        browse = (request.form.get("browse") or "").strip()
        browse_type = request.form.get("browse_type")
        window = 50000 if browse_type == 'gene' else 100000
        refGen = 'GRCh38'
        if browse_type == 'trait':
            # get info to plot umap and heatmap
            clusters_of_interest, all_indices, trait_list, cluster_index_map, umap, cluster_representatives, browse, sim_mat, names_sub, cluster_labels_for_index, cluster_traits, traits_interest, idx_trait, all_cluster_representatives = guide_haplotypes_traits(data_path, browse, window, refGen, trait_names)
            # get table of traits info
            traits_info_table, gwas_assoc_df = get_traits_info(data_path, traits_interest, cluster_representatives)
            # get plot
            plot_url, gwas_assoc_df = plot_haplotype_traits(clusters_of_interest, all_indices, trait_list, cluster_index_map, umap, cluster_representatives, browse, sim_mat, names_sub, cluster_labels_for_index, cluster_traits, idx_trait, all_cluster_representatives, gwas_assoc_df)
            return render_template("haplotypes.html", browse_value=browse, traits=trait_list, plot_url=plot_url, haplo_summary=[], chrom=None, start_pos=None, end_pos=None, hap_id_interest=None, table_traits=traits_info_table.to_dict(orient="records"), gwas_assoc_table=gwas_assoc_df.to_dict(orient="records"))
        elif browse_type == 'chromosome':
            return render_template("haplotypes.html", browse_value=browse, traits=trait_list, plot_url=None, haplo_summary=[])
        else:
            # get data and plots
            browse, plot_url, haplo_summary, chrom, start_pos, end_pos, hap_id_interest = guide_haplotypes_snps_genes(data_path, browse, window, refGen)
            # return render
            return render_template("haplotypes.html", browse_value=browse, traits=trait_list, plot_url=plot_url, haplo_summary=haplo_summary.to_dict(orient="records"), chrom=chrom, start_pos=start_pos, end_pos=end_pos, hap_id_interest=hap_id_interest, table_traits=None, gwas_assoc_table=None)
    # GET
    return render_template("haplotypes.html",
                           traits=trait_list,
                           plot_url=None,
                           haplo_summary=[])

@app.route("/haplotypes/detail")
def haplotype_detail():
    hap_id = request.args.get("hap_id")
    # Derive chromosome and position
    chrom = request.args.get("chrom")
    start_pos = int(request.args.get("start_pos"))
    end_pos = int(request.args.get("end_pos"))
    refGen = 'GRCh38'
    # get data of the region of interest
    haplo_df, haplo_dict = extract_haplo_data(data_path, chrom, start_pos, end_pos)
    haplo_df = haplo_df[haplo_df['ID'] == hap_id]   
    # get snps
    snps = ' '.join([str(chrom) + ':' + x.split(':')[1] + '-' + x.split(':')[1] for x in haplo_dict[hap_id]])
    # get ld between snps
    ld_df = get_ld_between_snps(snps, data_path)
    # get genes
    genes = extract_genes(data_path, chrom, start_pos, end_pos, refGen)
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
    fig.update_layout(height=1500)

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
    beta_matrix = snps_df.pivot(index='rsid', columns='trait', values='beta').fillna(0)
    beta_corr = beta_matrix.corr().round(2)
    if isinstance(beta_corr, pd.DataFrame):
        corr_vals = beta_corr.values
        traits = list(beta_corr.columns)
    else:
        corr_vals = beta_corr
        traits = [f"T{i+1}" for i in range(corr_vals.shape[0])]
    # clustering
    dist = 1 - corr_vals
    # condensed form for linkage()
    dist_condensed = squareform(dist, checks=False)
    Z = linkage(dist_condensed, method="average")
    # Order of items according to the dendrogram leaves
    order = leaves_list(Z)
    corr_clustered = corr_vals[np.ix_(order, order)]
    # original trait labels (before wrapping)
    traits_clustered = [traits[i] for i in order]
    # wrap long trait names
    traits_wrapped = [wrap_label(t, max_len=20) for t in traits_clustered]

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

    # Return as Plotly JSON dict
    return jsonify(fig=fig.to_dict(), haplo_table=haplo_table_html, snps_table=snps_table_html, snps_cadd_table=snps_cadd_table_html)

# Global functions
# function to get trait info
def get_traits_info(data_path, traits_interest, representative):
    all_gw_info = []
    gwas_assoc_df = pd.DataFrame()
    # iterate over df rows
    for index, row in traits_interest.iterrows():
        gwas_id = row['id'].replace('_hits_p5e-5.txt.gz', '')
        # get gwas info
        gwas_info = meta[gwas_id]
        # convert to dictionary
        gwas_info = dict(gwas_info)
        # make a dictionary
        all_gw_info.append(gwas_info)
        # check if this is the representative gwas
        if row['id'] == representative['Final Representative file'].values[0]:
            # in that case, get associations for manhattan plot
            gwas_associations = [x.split('\t') for x in subprocess.run(["tabix", data_path + "/databases/haplotypes/All_indep_gwas_sumstats_AI_hg38_byTrait.txt.gz", row['id'].replace('_hits_p5e-5.txt.gz', '')], capture_output=True,check=True, text=True).stdout.strip().split('\n')]
            # convert to dataframe
            gwas_assoc_df = pd.DataFrame(gwas_associations, columns=["id", "trait", "chr", "position", "rsid", "ea", "nea", "eaf", "beta", "se", "p", "n", "study_id", "position_hg38"])
    # convert all_gw_info to dataframe
    all_gw_info = pd.DataFrame(all_gw_info)
    # subset of columns to keep
    all_gw_info = all_gw_info[["id", "trait", "pmid", "author", "year", "sample_size", "population"]]
    # pmid as integer NA proof
    all_gw_info['pmid'] = all_gw_info['pmid'].fillna(0).astype(int)
    # rename columns
    all_gw_info.rename(columns={"id":"GWAS ID", "trait":"Trait", "pmid":"PubMed ID", "author":"First Author", "year":"Year", "sample_size":"Sample Size", "population":"Population"}, inplace=True)
    return all_gw_info, gwas_assoc_df

# function to extract number from gwas id
def extract_numeric_id(fname):
    try:
        # keep only the part after the second dash
        part = fname.split('-')[2].split('_')[0]
        # extract numeric prefix
        return int(num)
    except:
        return None

# function to plot haplotype traits
def plot_haplotype_traits(clusters_of_interest, all_indices, trait_list, cluster_index_map, umap, cluster_representatives, browse, sim_mat, names_sub, cluster_labels_for_index, cluster_trait, idx_trait, all_cluster_representatives, gwas_assoc_df):
    # ------------------------------------------------
    # Shared color map for clusters (UMAP + annotation)
    # ------------------------------------------------
    palette = [
        "#000000",  # black
        "#E69F00",  # orange
        "#56B4E9",  # sky blue
        "#009E73",  # bluish green
        "#0072B2",  # blue
        "#D55E00",  # vermilion
    ]
    color_map = {cid: palette[j % len(palette)] for j, cid in enumerate(clusters_of_interest)}

    # -------------------------------
    # Create FIGURE with 2 SUBPLOTS
    # -------------------------------
    fig = make_subplots(
        rows=2,
        cols=2,
        column_widths=[0.6, 0.40],
        row_heights=[0.6, 0.4],         # adjust as you want
        horizontal_spacing=0.04,
        vertical_spacing=0.25,         # spacing between row1 and row2
        specs=[
            [{"type": "xy"}, {"type": "xy"}],      # Row 1: UMAP | Heatmap
            [{"colspan": 2, "type": "xy"}, None]   # Row 2: NEW PLOT ACROSS BOTH COLUMNS
        ],
        subplot_titles=(
            "UMAP Trait Projection",
            "Trait Similarity Heatmap",
            "Manhattan plot of representative GWAS"                     # Title for row 2
        )
    )
    # increase subplot subtitle font size
    fig.update_annotations(font=dict(size=20, family="Arial", color="darkgrey", weight="bold"))

    # -------------------------------
    # LEFT subplot: UMAP
    # -------------------------------
    # Use the SAME indices used for the heatmap for colored clusters
    all_cluster_idxs = np.unique(all_indices)
    all_points = np.arange(len(trait_list))
    background_idxs = np.setdiff1d(all_points, all_cluster_idxs)
    # Grey background points (everything NOT in the clusters of interest subset)
    if background_idxs.size > 0:
        fig.add_trace(
            go.Scatter(
                x=umap[background_idxs, 0],
                y=umap[background_idxs, 1],
                mode="markers",
                marker=dict(size=6, color="lightgrey", opacity=0.3),
                showlegend=False,
                hoverinfo="skip",
            ),
            row=1,
            col=1
        )
    # Colored clusters (but capped to the same subset as the heatmap)
    for j, cid in enumerate(clusters_of_interest):
        # Intersect to keep at most the terms we actually used (max_terms_per_cluster)
        idxs_full = np.array(cluster_index_map[cid])
        idxs = np.intersect1d(idxs_full, all_cluster_idxs)
        if idxs.size == 0:
            continue
        fig.add_trace(
            go.Scatter(
                x=umap[idxs, 0],
                y=umap[idxs, 1],
                mode="markers",
                marker=dict(size=9, color=color_map[cid]),
                name=all_cluster_representatives.loc[all_cluster_representatives['Cluster Label'] == cid, 'Final Representative'].values[0],  # show representative name in legend
                text=[trait_list[i] for i in idxs],
                hoverinfo="text",
            ),
            row=1,
            col=1
        )
    # Highlight queried trait with a star
    fig.add_trace(
        go.Scatter(
            x=[umap[idx_trait, 0]],
            y=[umap[idx_trait, 1]],
            mode="markers+text",
            text=["★"],
            textposition="top center",
            marker=dict(
                size=15,
                color="black",
                symbol="star",
                line=dict(width=2, color="yellow"),
            ),
            name=browse,
            hovertext=browse,
            hoverinfo="text",
        ),
        row=1,
        col=1
    )
    fig.update_xaxes(title_text="UMAP-1", row=1, col=1, showticklabels=False, title_font=dict(size=20, family="Arial", color="black", weight="bold"))
    fig.update_yaxes(title_text="UMAP-2", row=1, col=1, showticklabels=False, title_font=dict(size=20, family="Arial", color="black", weight="bold"))
    # Legend below UMAP
    fig.update_layout(
        legend=dict(
            orientation="h",
            x=0,
            y=0.50,
            xanchor="left",
            font=dict(size=16, family="Arial"),
            itemsizing="constant"
        ),
    )

    # -------------------------------
    # RIGHT subplot: Heatmap
    # -------------------------------
    n_traits = sim_mat.shape[0]
    # Heatmap with flipped colors (blue = -1, red = +1) and no axis labels
    hover_text = [
        [
            f"{names_sub[i]} ↔ {names_sub[j]}"
            for j in range(n_traits)
        ]
        for i in range(n_traits)
    ]
    fig.add_trace(
        go.Heatmap(
            z=sim_mat,
            x=list(range(n_traits)),
            y=list(range(n_traits)),
            text=hover_text,
            hovertemplate=(
                "Trait pair: %{text}<br>"
                "Sim=%{z:.3f}<extra></extra>"
            ),
            zmin=0,
            zmax=1,
            colorscale=[
                [0.0, "blue"],   # -1
                [0.5, "white"],  #  0
                [1.0, "red"],    # +1
            ],
            colorbar=dict(title="Cosine similarity", orientation="v", x=1.02, y=0.83, len=0.3, thickness=15),
            showscale=True,
        ),
        row=1,
        col=2
    )
    # Cluster annotation strip on the right
    cluster_colors_for_traits = [
        color_map.get(cid, "#d3d3d3") for cid in cluster_labels_for_index
    ]
    fig.add_trace(
        go.Scatter(
            x=[n_traits + 0.5] * n_traits,   # vertical strip to the right
            y=list(range(n_traits)),
            mode="markers",
            marker=dict(
                color=cluster_colors_for_traits,
                symbol="square",
                size=8,
                line=dict(width=0),
            ),
            showlegend=False,
            hoverinfo="skip",
        ),
        row=1,
        col=2
    )
    # Hide axis labels / ticks for heatmap panel
    fig.update_xaxes(
        showticklabels=False,
        title_text="",
        range=[-0.5, n_traits + 1.0],  # make room for annotation strip
        row=1,
        col=2
    )
    fig.update_yaxes(
        showticklabels=False,
        title_text="",
        row=1,
        col=2
    )

    # -------------------------------
    # BOTTOM subplot: Manhattan of representative GWAS
    # -------------------------------
    chrom_lengths = {1 : 248956422, 2: 242193529, 3: 198295559, 4: 190214555, 5: 181538259, 6 : 170805979, 7: 159345973, 8: 145138636, 9: 138394717, 10: 133797422, 11 : 135086622, 12: 133275309, 13: 114364328, 14: 107043718, 15: 101991189, 16 : 90338345, 17: 83257441, 18: 80373285, 19: 58617616, 20: 64444167, 21 : 46709983, 22: 50818468}
    if not gwas_assoc_df.empty:
        # Prepare data
        gwas_assoc_df['position_hg38'] = gwas_assoc_df['position_hg38'].astype(int)
        gwas_assoc_df['neg_log10_p'] = -np.log10(gwas_assoc_df['p'].astype(float))
        # Sort by pvalue
        gwas_assoc_df['p'] = gwas_assoc_df['p'].astype(float)
        gwas_assoc_df = gwas_assoc_df.sort_values(by='p').reset_index(drop=True)
        # Take top 10000 associations if more than that
        if gwas_assoc_df.shape[0] > 10000:
            gwas_assoc_df = gwas_assoc_df.iloc[:10000, :].copy().reset_index(drop=True)
        # Sort by chromosome and position
        gwas_assoc_df['chr'] = gwas_assoc_df['chr'].astype(int)
        gwas_assoc_df['position_hg38'] = gwas_assoc_df['position_hg38'].astype(int)
        gwas_assoc_df = gwas_assoc_df.sort_values(by=['chr', 'position_hg38']).reset_index(drop=True)
        # Compute cumulative positions -- for chr 1 it's the same, for chr 2 onwards, add the lengths of previous chromosomes
        cum_length = 0
        cum_positions = []
        for chrom in range(1, 23):
            chrom_assoc = gwas_assoc_df[gwas_assoc_df['chr'] == chrom]
            for pos in chrom_assoc['position_hg38']:
                cum_positions.append(pos + cum_length)
            cum_length += chrom_lengths[chrom]
        gwas_assoc_df['cum_position'] = cum_positions
        # Alternating chromosome colors
        colors = ["#1f77b4", "#ff7f0e"]
        gwas_assoc_df["color"] = gwas_assoc_df["chr"] % 2
        gwas_assoc_df["color"] = gwas_assoc_df["color"].map({0: colors[0], 1: colors[1]})
        # Scatter plot
        fig.add_trace(
            go.Scatter(
                x=gwas_assoc_df['cum_position'],
                y=gwas_assoc_df['neg_log10_p'],
                mode='markers',
                marker=dict(
                    color=gwas_assoc_df['color'],
                    size=7,
                    line=dict(width=1, color=gwas_assoc_df['color'])
                ),
                name='Associations',
                text=[
                    f'Position: {row["chr"]}:{row["position_hg38"]}<br>'
                    f'RsID: {row["rsid"]}<br>'
                    f'Effect Allele: {row["ea"]}<br>'
                    f'Other Allele: {row["nea"]}<br>'
                    f'P-value: {row["p"]}<br>'
                    f'Beta: {row["beta"]} (SE: {row["se"]})'
                    for _, row in gwas_assoc_df.iterrows()
                ],
                hoverinfo='text'
            ),
            row=2,
            col=1
        )
        # Genome-wide significance line
        fig.add_hline(y=-np.log10(5e-08), line=dict(color="red", width=2, dash="dash"), annotation_text="Genome-wide significant", annotation_position="top right", row=2, col=1)
        fig.add_hline(y=-np.log10(5e-05), line=dict(color="navy", width=2, dash="dash"), annotation_text="Suggestive significant", annotation_position="top right", row=2, col=1)
        # X-axis ticks at chromosome centers
        ticks = gwas_assoc_df.groupby("chr")["cum_position"].mean()
        tick_labels = gwas_assoc_df.groupby("chr")["chr"].first()
        fig.update_xaxes(title_text="Chromosomes (hg38)", row=2, col=1, title_font=dict(size=20, family="Arial", color="black", weight="bold"), tickmode="array", tickvals=ticks, ticktext=tick_labels)
        fig.update_yaxes(title_text="-log10(P-value)", row=2, col=1, title_font=dict(size=20, family="Arial", color="black", weight="bold"))
    else:
        fig.add_annotation(
            x=0.5,
            y=0.5,
            xref="x2 domain",
            yref="y2 domain",
            text="No GWAS associations available for the representative trait.",
            showarrow=False,
            font=dict(size=16, family="Arial", color="darkgrey"),
        )
    
    # Overall layout
    fig.update_layout(
        template="plotly_white",
        height=1100,
        title={
            "text": f"<b>Semantic Neighborhood of '{browse}'</b>",
            "font": {"size": 32},
        }
    )
    # return as HTML
    plot_url = pio.to_html(fig, include_plotlyjs='cdn', full_html=False)
    # subset of the gwas_assoc_df to return
    gwas_assoc_df = gwas_assoc_df[['id', 'trait', 'chr', 'position_hg38', 'rsid', 'ea', 'nea', 'eaf', 'beta', 'se', 'p', 'n']].copy()
    # rename columns
    gwas_assoc_df.rename(columns={'id':'GWAS ID', 'trait':'Trait', 'chr':'Chromosome', 'position_hg38':'Position (hg38)', 'rsid':'RsID', 'ea':'Effect Allele', 'nea':'Other Allele', 'eaf':'EAF', 'beta':'Beta', 'se':'SE', 'p':'P-value', 'n':'N'}, inplace=True)
    # and sort by p-value
    gwas_assoc_df = gwas_assoc_df.sort_values(by='P-value').reset_index(drop=True)
    return plot_url, gwas_assoc_df

# function to guide haplotype plots with traits as inputs
def guide_haplotypes_traits(data_path, browse, window, refGen, trait_names):
    # Put traits into a list
    trait_list = trait_names['name'].values.tolist()
    # Class to handle numpy unpickling issues
    class NumpyCompatUnpickler(pickle.Unpickler):
        def find_class(self, module, name):
            # Remap old internal numpy paths to the current ones
            if module == "numpy._core":
                module = "numpy.core"
            elif module == "numpy._core.multiarray":
                module = "numpy.core.multiarray"
            return super().find_class(module, name)
    # import pickle object with embeddings
    with open(f"{data_path}/databases/haplotypes/trait_embeddings_umap_thr05.pkl", "rb") as f:
        embeddings_obj = NumpyCompatUnpickler(f).load()
    # Also read cluster representatives
    cluster_representatives = pd.read_csv(f"{data_path}/databases/haplotypes/trait_clusters_summary_AI_thresholds.csv", sep=",")
    cluster_representatives = cluster_representatives[cluster_representatives['Threshold'] == 0.5].copy().reset_index(drop=True)

    # Derive umap, embeddings, labels, cluster_map
    umap = embeddings_obj["umap"]
    embeddings = embeddings_obj["embeddings"]
    labels = embeddings_obj["labels"]
    cluster_map = embeddings_obj["cluster_map"]
    
    # Find trait index
    idx_trait = find_trait_index(browse, trait_list)
    # Find cluster id
    cluster_trait = int(labels[idx_trait])
    # Build cluster to index map
    cluster_index_map = build_cluster_index_map(labels)
    # Compute centroids
    centroids = {}
    for cid, idxs in cluster_index_map.items():
        c = embeddings[idxs].mean(axis=0)
        c /= np.linalg.norm(c) + 1e-9
        centroids[cid] = c

    # Compute similarity of all clusters to trait's cluster
    c0 = centroids[cluster_trait].reshape(1, -1)
    other_clusters = [cid for cid in centroids if cid != cluster_trait]
    sims = cosine_similarity(c0, np.vstack([centroids[c] for c in other_clusters]))[0]
    order = np.argsort(-sims)
    # Pick 5 closest clusters
    k = min(5, len(other_clusters))
    neighbor_clusters = [other_clusters[i] for i in order[:k]]
    clusters_of_interest = [cluster_trait] + neighbor_clusters

    # ------------------------------------------------
    # Build similarity heatmap between traits in clusters
    # ------------------------------------------------
    all_indices = []
    cluster_labels_for_index = []
    max_terms_per_cluster = 100
    for cid in clusters_of_interest:
        idxs = cluster_index_map[cid]
        if len(idxs) > max_terms_per_cluster:
            idxs = idxs[:max_terms_per_cluster]
        all_indices.extend(idxs)
        cluster_labels_for_index.extend([cid] * len(idxs))
    all_indices = np.array(all_indices)
    cluster_labels_for_index = np.array(cluster_labels_for_index)
    emb_sub = embeddings[all_indices]
    names_sub = [trait_list[i] for i in all_indices]
    # Order by cluster ID, then by similarity to the queried trait
    order = np.argsort(cluster_labels_for_index)
    sim_to_trait = cosine_similarity(emb_sub, embeddings[idx_trait].reshape(1, -1))[:, 0]
    new_order = []
    for cid in clusters_of_interest:
        mask = cluster_labels_for_index == cid
        idxs = np.where(mask)[0]
        if len(idxs) > 0:
            local = idxs[np.argsort(-sim_to_trait[idxs])]
            new_order.extend(local.tolist())
    order = np.array(new_order)
    emb_sub = emb_sub[order]
    names_sub = [names_sub[i] for i in order]
    cluster_labels_for_index = cluster_labels_for_index[order]
    sim_mat = cosine_similarity(emb_sub)

    # Also get the list of trait in the cluster of interest
    traits_in_same_cluster = [trait_list[i] for i in cluster_index_map[cluster_trait]]
    traits_sub = trait_names[trait_names['name'].isin(traits_in_same_cluster)].copy().reset_index(drop=True)
    representative = cluster_representatives[cluster_representatives['Cluster Label'] == cluster_trait]
    return clusters_of_interest, all_indices, trait_list, cluster_index_map, umap, representative, browse, sim_mat, names_sub, cluster_labels_for_index, cluster_trait, traits_sub, idx_trait, cluster_representatives
    
# function to build cluster to index map
def build_cluster_index_map(labels):
    """
    labels: np.array shape (n_traits,)
    return: dict cluster_id -> list of indices
    """
    cluster_index_map = {}
    for i, lab in enumerate(labels):
        cluster_index_map.setdefault(int(lab), []).append(i)
    return cluster_index_map

# function to find trait index
def find_trait_index(trait_query, trait_names):
    """
    trait_query: string
    trait_names: list-like of trait names (len = n_traits, same order as embeddings/labels)
    Returns index of best match (exact, then case-insensitive, then substring).
    """
    names = pd.Series(trait_names)
    # 1) exact
    exact = names[names == trait_query]
    if len(exact) == 1:
        return exact.index[0]
    elif len(exact) > 1:
        return exact.index[0]
    # 2) case-insensitive exact
    ci = names.str.lower() == trait_query.lower()
    if ci.any():
        return np.where(ci)[0][0]
    # 3) substring (case-insensitive)
    sub = names.str.lower().str.contains(trait_query.lower())
    if sub.any():
        return np.where(sub)[0][0]
    raise ValueError(f"Trait '{trait_query}' not found in trait_names.")

# function to guide haplotype plots with snps/genes as inputs
def guide_haplotypes_snps_genes(data_path, browse, window, refGen):
    # Derive chromosome and position
    chrom, start_pos, end_pos, browse_type = readBrowseOption(data_path, browse, window, refGen)
    # get data of the region of interest
    haplo_df, haplo_dict = extract_haplo_data(data_path, chrom, start_pos, end_pos)
    # identify the haplotype when an individual SNP is searched
    hap_id_interest = None
    if browse_type == 'RsID' or browse_type == 'Single position':
        if browse_type == 'RsID':
            # grep rsid in haplo_df rsids column and take the ID; treat NA as False and ignore case
            hap_id_interest = haplo_df[haplo_df['rsids'].str.contains(browse, na=False, case=False)]['ID'].unique().tolist()
        else:
            # grep locus in haplo_df SNPS column and take the ID; treat NA as False
            locus_pos = ':' + browse.split(':')[1] + ':'
            hap_id_interest = haplo_df[haplo_df['SNPS'].str.contains(locus_pos, na=False)]['ID'].unique().tolist()
    # get genes
    genes = extract_genes(data_path, chrom, start_pos, end_pos, refGen)
    haplo_summary = (haplo_df.groupby("ID").agg(n_snps=('NSNPS', 'max'), start_bp=('BP1', 'min'), end_bp=('BP2', 'max'), n_traits=('n_traits', 'max'), traits=('traits', 'max')).reset_index())
    # sort by position and assign an index
    haplo_summary = haplo_summary.sort_values(by='start_bp').reset_index(drop=True)
    haplo_summary['haplo_index'] = haplo_summary.index
    # sort haplo summary by number of traits descending
    haplo_summary = haplo_summary.sort_values(by='n_traits', ascending=False).reset_index(drop=True)
    # reorder columns
    haplo_summary = haplo_summary[['haplo_index', 'ID', 'n_snps', 'start_bp', 'end_bp', 'n_traits', 'traits']]
    # if hap_id_interest is not empty, move that haplotype on top
    if hap_id_interest == []:
        hap_id_interest = 'Not found'
    else:
        haplo_summary['is_interest'] = haplo_summary['ID'].apply(lambda x: x in hap_id_interest)
        haplo_summary = haplo_summary.sort_values(by='is_interest', ascending=False).drop(columns=['is_interest']).reset_index(drop=True)
        hap_id_interest = hap_id_interest[0]  # take the first one for highlighting in the plot
    # create plot
    fig = haplotype_plotly_regional(haplo_df, haplo_dict, genes)
    plot_url = pio.to_html(fig, include_plotlyjs='cdn', full_html=False)
    return browse, plot_url, haplo_summary, chrom, start_pos, end_pos, hap_id_interest

# function to query cadd scores
def query_cadd_score(snps, data_path):
    """
    Query CADD annotation for a given hg38 chromosome and position.
    Returns a list of dictionaries (JSON-serializable).
    """
    try:
        cadd_info = []
        for snp in snps.split(" "):
            chrom = snp.split(":")[0].replace("chr", "")
            pos38 = int(snp.split(":")[1].split("-")[0])
            chrom = str(chrom).upper().replace("CHR", "")
            CADD_DB_FILE = f"{data_path}/databases/CADD/subset_chr{chrom}.db"
            conn = sqlite3.connect(CADD_DB_FILE)
            df = pd.read_sql_query(
                "SELECT * FROM variants_combined WHERE pos = ?;",
                conn,
                params=(pos38,),
            )
            if not df.empty:
                cadd_info.append(df)
        conn.close()
        # concatenate all dataframes in the list
        if cadd_info:
            df = pd.concat(cadd_info, ignore_index=True)
            df = df.sort_values(by="phred_max", ascending=False)
            return df
        else:
            return []
    except Exception:
        return []

# function to get rsid from position
def lookup_by_coord_hg38(chrom, pos38):
    DB_FILE = "%s/databases/Genes/variants_info.sqlite" % (data_path)
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

# function to get LD between SNPs
def get_ld_between_snps(snps, data_path):
    snp_list = snps.split(' ')
    snp_positions = [int(snp.split(':')[1].split('-')[0]) for snp in snp_list]
    chr = snp_list[0].split(':')[0]
    # get rsids for the SNPs
    rsids_list = {}
    for snp in snp_positions:
        rsid_info = lookup_by_coord_hg38('chr' + chr, snp)
        rsids_list[snp] = rsid_info['rsid'] if rsid_info else 'NA'
    # convert to dataframe
    rsids_df = pd.DataFrame.from_dict(rsids_list, orient='index', columns=['rsid']).reset_index().rename(columns={'index': 'position'})
    # define interval
    min_pos = min(snp_positions)
    max_pos = max(snp_positions)
    # run function to get LD by interval
    db_path = f"{data_path}/databases/LD_db/ld_chr{chr}.sqlite"
    ld_results = partners_by_interval(db_path, min_pos, max_pos, r2_min=0.2)
    # create a dataframe
    if ld_results:
        ld_df = pd.DataFrame(ld_results, columns=['anchor_uniq', 'partner_uniq', 'r2', 'dist_bp'])
        # extract positions for anchor and partner
        ld_df['anchor_pos'] = ld_df['anchor_uniq'].apply(lambda x: int(x.split(':')[0]))
        ld_df['partner_pos'] = ld_df['partner_uniq'].apply(lambda x: int(x.split(':')[0]))
        # filter to keep only the requested SNPs
        ld_df = ld_df[ld_df['anchor_pos'].isin(snp_positions) & ld_df['partner_pos'].isin(snp_positions)]
        # deduplicate
        # create unordered pair columns
        ld_df["p1"] = ld_df[["anchor_pos", "partner_pos"]].min(axis=1)
        ld_df["p2"] = ld_df[["anchor_pos", "partner_pos"]].max(axis=1)
        # if you want to keep the "best" record per pair (e.g. highest r2, then smallest distance)
        ld_df = (ld_df.sort_values(["r2", "dist_bp"], ascending=[False, True]).drop_duplicates(subset=["p1", "p2"], keep="first"))
        # optional: drop helper cols
        ld_df = ld_df.drop(columns=["p1", "p2"])
        # merge with rsids
        ld_df = ld_df.merge(rsids_df.rename(columns={'position': 'anchor_pos', 'rsid': 'anchor_rsid'}), on='anchor_pos', how='left')
        ld_df = ld_df.merge(rsids_df.rename(columns={'position': 'partner_pos', 'rsid': 'partner_rsid'}), on='partner_pos', how='left')
        # when anchor_rsid is NA, fill with anchor_uniq
        ld_df['anchor_rsid'] = ld_df.apply(lambda row: row['anchor_uniq'] if row['anchor_rsid'] == 'NA' else row['anchor_rsid'], axis=1)
        ld_df['partner_rsid'] = ld_df.apply(lambda row: row['partner_uniq'] if row['partner_rsid'] == 'NA' else row['partner_rsid'], axis=1)
        # finally, convert to a df with same rows and columns -- colnames and rownames are anchor_rsid and partner_rsid -- value is r2 -- NaN if no LD
        ld_df = ld_df.pivot(index='anchor_rsid', columns='partner_rsid', values='r2')
        # Build a unified ordered list of SNP IDs
        snps = sorted(set(ld_df.index) | set(ld_df.columns))
        # Reindex both the rows and columns using the same order
        ld_df = ld_df.reindex(index=snps, columns=snps)
        for i in range(len(snps)):
            ld_df.iat[i, i] = 1.0
        ld_df = ld_df.fillna(0.0)   # missing pairs → 0
        return ld_df
    return None

# function to get LD by interval
def partners_by_interval(db_path: str, start_pos: int, end_pos: int, r2_min: float = 0.2, limit: Optional[int] = None) -> Sequence[Any]:
    """
    Return partner variants in LD with any variant whose `pos` is in [start_pos, end_pos].
    Results are grouped by partner_uniq (uniq of the *partner* variant).
    """
    R2_SCALE = 1000
    r2m_min = int(round(max(0.0, min(1.0, r2_min)) * R2_SCALE))
    sql = f"""
    WITH v AS (
      SELECT id AS anchor_id
      FROM variants
      WHERE pos BETWEEN ? AND ?
    ),
    hits AS (
      SELECT
        v.anchor_id,
        CASE WHEN ld.v1 = v.anchor_id THEN ld.v2 ELSE ld.v1 END AS partner_id,
        ld.r2_milli,
        ld.dist_bp
      FROM ld
      JOIN v ON ld.v1 = v.anchor_id OR ld.v2 = v.anchor_id
      WHERE ld.r2_milli >= ?
    )
    SELECT
      anchor.uniq  AS anchor_uniq,
      partner.uniq AS partner_uniq,
      CAST(hits.r2_milli AS REAL) / {R2_SCALE} AS r2,
      hits.dist_bp
    FROM hits
    JOIN variants AS anchor  ON anchor.id  = hits.anchor_id
    JOIN variants AS partner ON partner.id = hits.partner_id
    ORDER BY r2 DESC, dist_bp ASC
    {f"LIMIT {int(limit)}" if limit else ""}
    """
    with _open_db_for_query(db_path) as conn:
        return conn.execute(sql, (start_pos, end_pos, r2m_min)).fetchall()

# function to wrap long trait names
def wrap_label(s, max_len=20):
    s = str(s)
    if len(s) <= max_len:
        return s
    # split into chunks of ~max_len and join with <br>
    return "<br>".join(textwrap.wrap(s, max_len))

# function to set height of the figure based on the number of haplotypes
def set_haplotype_plot_height(haplo_df):
    if not haplo_df.empty:
        n_haplotypes = haplo_df.shape[0]
        if n_haplotypes <= 5:
            px = 600
        elif n_haplotypes <= 10:
            px = 700
        elif n_haplotypes <= 20:
            px = 800
        elif n_haplotypes <= 30:
            px = 1000
        else:
            px = 1200
    else:
        px = 400
    return px

# function to plot haplotype regional plot
def haplotype_plotly_regional(haplo_df, haplo_dict, genes, refGen='GRCh38'):
    # sort haplo_df by position
    haplo_df = haplo_df.sort_values(by='BP1', ascending=True).reset_index(drop=True)
    # set figure height based on number of haplotypes
    height_px = set_haplotype_plot_height(haplo_df)
    # set up figure -- 2 panels associations and genes
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.02, row_heights=[0.45, 0.275], specs=[[{"secondary_y": True}], [{}]])
    fig.update_layout(height=height_px)
    for r in [1, 2]:
        fig.update_xaxes(showgrid=True, gridcolor="#e0e0e0", zeroline=False, row=r, col=1)
        fig.update_yaxes(showgrid=True, gridcolor="#e0e0e0", zeroline=False, row=r, col=1)
    fig.update_layout(plot_bgcolor="white", paper_bgcolor="white")
    # define color palette, violet to white to red with 100 values
    cmap = LinearSegmentedColormap.from_list("violet_white_red", ['violet', 'white', 'red'])
    colors = [to_hex(cmap(i/99)) for i in range(100)]
    # find maximum and minimum number of traits in haplo_df
    max_n_traits = haplo_df['n_traits'].max() if not haplo_df.empty else 1
    min_n_traits = haplo_df['n_traits'].min() if not haplo_df.empty else 0
    # iterate over haplo_df rows and plot each haplotype block
    for index, row in haplo_df.iterrows():
        haplo_id = row['ID']
        snps = row['SNPS'].split('|')
        snps_positions = [int(x.split(':')[1]) / 1_000_000 for x in snps]
        n_traits = int(row['n_traits'])
        n_traits_list = [n_traits for _ in snps]
        trait_list = row['traits'].split('___') if pd.notna(row['traits']) else []
        index_list = [index for _ in snps]
        # the color depends on the number of traits, normalized to 0-100 based on the maximum number of traits in the haplo_df
        counter = int(((n_traits - min_n_traits) / (max_n_traits - min_n_traits)) * 99) if max_n_traits != min_n_traits else 50
        # create a scatter trace for the haplotype block
        fig.add_trace(go.Scatter(
            x=snps_positions,
            y=index_list,
            mode='markers+lines',
            line=dict(color="grey", width=2),  # <- explicitly set line color
            marker=dict(
                color=n_traits_list,
                size=10,
                line=dict(width=1, color="black"),
                colorscale="Spectral",
                colorbar=dict(title="r", orientation="v", x=1.02, y=0.83, len=0.3, thickness=15),
                cmin=min_n_traits,
                cmax=max_n_traits,
                showscale=(index == 0)  # show only once
            ),
            name=f'Haplotype {haplo_id}',
            text=[
                f'Haplotype: {haplo_id}<br>SNP: {snp}<br>Traits: {nt}' 
                for snp, nt in zip(snps, n_traits_list)
            ],
            hoverinfo='text'
        ), row=1, col=1)

    # update layout
    fig.update_layout(title='', xaxis_title='', yaxis_title='', showlegend=False)

    # ensure y-axis ticks are integers only (0,1,2,...) and set a small padding around the range
    max_idx = (haplo_df.shape[0] - 1) if not haplo_df.empty else 1
    fig.update_yaxes(
        range=[-0.2, max_idx + 0.2],
        tick0=0,
        dtick=1,
        tickmode='linear',
        title_text="Haplotype number",
        title_font=dict(size=18, family="Arial", color="black", weight="bold"),
        tickfont=dict(size=14, family="Arial", color="black"),
        row=1,
        col=1
    )
    # then the second panel with the genes
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
            row = np.array([str(g[index_text]), str(g[index_transc]), str(g[index_strand]), int(g[index_pos]), int(g[index_pos + 1]), int(g[index_exonN])], dtype=object)
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
                    showlegend=False), row=2, col=1)
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
                ), row=2, col=1)

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
                x=x0 + (x1 - x0) / 2, y=ylane + yoff,
                text=label, showarrow=False, font=dict(size=fsize), hovertext=hovertext, row=2, col=1)

        fig.update_yaxes(range=[-max_y*0.20, max_y+max_y*0.20], row=2, col=1)
    else:
        fig.add_annotation(
            x=(start_pos + (end_pos-start_pos)/2)/1_000_000,
            y=0.5, text="No Genes to plot", showarrow=False,
            row=2, col=1
        )
    fig.update_yaxes(title_text="Gene track", row=2, col=1, title_font=dict(size=18, family="Arial", color="black", weight="bold"), tickfont=dict(size=14, family="Arial", color="white"))
    return fig

# function to get haplotype data of interest
def extract_haplo_data(data_path, chrom, start_pos, end_pos):
    # define data to look at
    haplo_path = f"{data_path}/databases/haplotypes/ld_clusters_with_gwas_AI_latest.tsv.gz".replace(' ', '\ ')
    # tabix query
    chrom = str(chrom).upper().replace('CHR', '')
    tabix_cmd = f"tabix {haplo_path} {chrom}:{start_pos}-{end_pos}"
    data_list = [x for x in os.popen(tabix_cmd)]
    # convert list to df
    haplo_df = pd.read_csv(StringIO(''.join(data_list)), sep='\t', header=None)
    haplo_df.columns = ['CHR', 'BP1', 'BP2', 'KB', 'NSNPS', 'SNPS', 'assoc_count', 'rsids', 'traits']
    # add column for the number of traits
    haplo_df['n_traits'] = haplo_df['traits'].apply(
        lambda x: len(set(t.strip() for t in x.split('___') if t and t.strip().lower() != 'nan')) 
        if pd.notnull(x) else 0
    )
    # sort by n_traits descending
    haplo_df = haplo_df.sort_values(by='n_traits', ascending=False)
    # add haplotype ID as chrom:BP1-BP2
    haplo_df['ID'] = haplo_df.apply(lambda x: f"{x['CHR']}:{x['BP1']}-{x['BP2']}", axis=1)
    # get snps as well
    haplo_dict = {}
    for index, row in haplo_df.iterrows():
        haplo_id = row['ID']
        snps = row['SNPS'].split('|')
        haplo_dict[haplo_id] = snps
    return haplo_df, haplo_dict

# function to get variants in LD
def extract_ld(data_path, chrom, df, refGen, browse_type, browse):
    if browse_type != 'RsID':
        # take most significant variant in the region, create a string like chr:pos:ea:nea
        top_snp = df.loc[df['Pvalue'].idxmax()]     
    else:
        # take the rsID from the browse field, create a string like chr:pos:ea:nea
        rsid = browse.replace(' ', '').split()[0].lower()
        # get the relative information from the df
        top_snp = df[df['Rsid'].str.lower() == rsid]
    if refGen == 'GRCh37':
        uniq1 = f"{top_snp['Position_hg38']}:{top_snp['EA']}:{top_snp['NEA']}"
        uniq2 = f"{top_snp['Position_hg38']}:{top_snp['NEA']}:{top_snp['EA']}"
    else:
        uniq1 = f"{top_snp['Position']}:{top_snp['EA']}:{top_snp['NEA']}"
        uniq2 = f"{top_snp['Position']}:{top_snp['NEA']}:{top_snp['EA']}"       
    # define database path
    db_path = f"{data_path}/databases/LD_db/ld_chr{chrom}.sqlite"
    # query the database
    res1 = partners_by_uniq(db_path, uniq1, r2_min=0.2, limit=100)
    res2 = partners_by_uniq(db_path, uniq2, r2_min=0.2, limit=100)
    # combine results
    res = res1 + res2
    # convert to df
    ld_df = pd.DataFrame(res, columns=['partner_uniq', 'r2', 'dist_bp'])
    # add the original variant
    ld_df = pd.concat([ld_df, pd.DataFrame({'partner_uniq': [uniq1], 'r2': [1.0], 'dist_bp': [0]})], ignore_index=True)
    # extract position from partner_uniq
    ld_df['pos'] = ld_df['partner_uniq'].apply(lambda x: int(x.split(':')[0]))
    # assign colors for ld: >=0.8 red, >=0.6 orange, >=0.4 green, >=0.2 lightblue, <0.2 grey
    ld_df['color'] = ld_df['r2'].apply(assign_color)
    return ld_df

# function to assign colors for ld
def assign_color(r2):
    if r2 >= 0.8:
        return 'red'
    elif r2 >= 0.6:
        return 'orange'
    elif r2 >= 0.4:
        return 'green'
    elif r2 >= 0.2:
        return 'lightblue'
    else:
        return 'grey'

# function to get variants in LD by id
def partners_by_uniq(db_path: str, uniq: str, r2_min: float = 0.0, limit: Optional[int] = None):
    """
    Return partner uniq IDs, r2 (float), and distance for a focal uniq ID.
    """
    DIST_MAX     = 500_000    # keep LD pairs with distance <= 500kb
    BATCH_SIZE   = 100_000    # rows per batch insert (tweak 50k-200k)
    VAR_FLUSH_EVERY = 200_000 # how often to flush new variants to DB
    PAGE_SIZE    = 32768      # 32KiB pages help large DBs (VACUUM required)
    CACHE_PAGES  = 200_000    # ~200k pages -> ~6.4GB cache at 32K per page (adjust)
    MMAP_BYTES   = 1<<30      # 1GiB mmap (adjust or set 0 if OS limits)
    R2_SCALE     = 1000       # store r2 as int(r2 * R2_SCALE)
    r2m_min = int(round(max(0.0, min(1.0, r2_min)) * R2_SCALE))
    sql = f"""
    WITH v AS (
      SELECT id FROM variants WHERE uniq = ?
    ),
    hits AS (
      SELECT
        CASE WHEN ld.v1 = v.id THEN ld.v2 ELSE ld.v1 END AS partner_id,
        ld.r2_milli, ld.dist_bp
      FROM ld, v
      WHERE ld.v1 = v.id OR ld.v2 = v.id
        AND ld.r2_milli >= ?
    )
    SELECT variants.uniq AS partner_uniq,
           CAST(hits.r2_milli AS REAL)/{R2_SCALE} AS r2,
           hits.dist_bp
    FROM hits
    JOIN variants ON variants.id = hits.partner_id
    GROUP BY partner_uniq   -- collapse any accidental dup entries
    ORDER BY r2 DESC, dist_bp ASC
    {f"LIMIT {int(limit)}" if limit else ""}
    """
    with _open_db_for_query(db_path) as conn:
        cur = conn.cursor()
        cur.execute(sql, (uniq, r2m_min))
        return cur.fetchall()

# function to get variants in LD by position
def partners_by_position(db_path: str, pos: int, r2_min: float = 0.0, limit: Optional[int] = None):
    """
    Same as above but resolve uniq by position (exact match of 'pos' column).
    """
    DIST_MAX     = 500_000    # keep LD pairs with distance <= 500kb
    BATCH_SIZE   = 100_000    # rows per batch insert (tweak 50k-200k)
    VAR_FLUSH_EVERY = 200_000 # how often to flush new variants to DB
    PAGE_SIZE    = 32768      # 32KiB pages help large DBs (VACUUM required)
    CACHE_PAGES  = 200_000    # ~200k pages -> ~6.4GB cache at 32K per page (adjust)
    MMAP_BYTES   = 1<<30      # 1GiB mmap (adjust or set 0 if OS limits)
    R2_SCALE     = 1000       # store r2 as int(r2 * R2_SCALE)
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
      WHERE ld.v1 = v.id OR ld.v2 = v.id
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
    with _open_db_for_query(db_path) as conn:
        return conn.execute(sql, (pos, r2m_min)).fetchall()

# helper to open a database connection with safer defaults
def _open_db_for_query(db_path: str) -> sqlite3.Connection:
    DIST_MAX     = 500_000    # keep LD pairs with distance <= 500kb
    BATCH_SIZE   = 100_000    # rows per batch insert (tweak 50k-200k)
    VAR_FLUSH_EVERY = 200_000 # how often to flush new variants to DB
    PAGE_SIZE    = 32768      # 32KiB pages help large DBs (VACUUM required)
    CACHE_PAGES  = 200_000    # ~200k pages -> ~6.4GB cache at 32K per page (adjust)
    MMAP_BYTES   = 1<<30      # 1GiB mmap (adjust or set 0 if OS limits)
    R2_SCALE     = 1000       # store r2 as int(r2 * R2_SCALE)
    conn = sqlite3.connect(db_path)
    # Safer defaults for querying; WAL helps concurrency
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=NORMAL;")
    conn.execute(f"PRAGMA mmap_size={MMAP_BYTES};")
    return conn

# function to get haplotypes in the region
def extract_haplo(data_path, df, chrom, refGen):
    # take min and max position wrt hg38
    if refGen == 'GRCh37':
        min_pos = df['Position_hg38'].min()
        max_pos = df['Position_hg38'].max()
    else:
        min_pos = df['Position'].min()
        max_pos = df['Position'].max()
    # take positions of interest
    cmd = 'tabix %s/databases/HaploBlocks/chrAll_haploblocks.txt.gz chr%s:%s-%s' %(data_path.replace(' ', '\ '), str(chrom).replace('chr', ''), str(min_pos - 1000), str(max_pos + 1000))
    haplo = [x.rstrip().split('\t') for x in os.popen(cmd)]
    # iterate over haplotypes and flag the variants in df
    df_haplo = assign_haplo_to_variants(df, haplo, refGen)
    return df_haplo

# function to assign haplotype to variants
def assign_haplo_to_variants(df, haplo, refGen):
    # normalize df chrom as 'chrN'
    pos_col = 'Position_hg38' if refGen == 'GRCh37' else 'Position'
    df_ = df.copy()
    df_["chrom"] = "chr" + df_["Chrom"].astype(str).str.replace("^chr", "", regex=True)
    df_["pos"]   = df_[pos_col].astype(int)
    df_["EA"]    = df_["EA"].str.upper()
    df_["NEA"]   = df_["NEA"].str.upper()
    # build (chrom, pos) -> block_id map from haplo list
    site2block = {}
    block_id2name = {}
    for i, b in enumerate(haplo):
        chrom, start, end, length_mb, nsites, sites = b
        # a readable name (optional)
        bid_name = f"{chrom}:{start}-{end}"
        block_id2name[i] = bid_name
        for s in sites.split("|"):
            # robust split in case of indels like G:GC
            c, p, a1, a2 = s.split(":", 3)  # "16:81905859:C:G"
            ck = "chr" + c.lstrip("chr")
            p  = int(p)
            a1 = a1.upper()
            a2 = a2.upper()
            # store both directions so lookup is O(1) without a second pass
            site2block[(ck, p, a1, a2)] = i
            site2block[(ck, p, a2, a1)] = i
    # assign block IDs to df
    keys = list(zip(df_["chrom"], df_["pos"], df_["EA"], df_["NEA"]))
    block_idx = pd.Series(keys).map(site2block)  # returns block index or NaN
    # if you prefer a readable label instead of numeric index:
    block_lbl = block_idx.map(block_id2name)
    out = pd.DataFrame({
        "variant": df_.get("Rsid", df_["chrom"] + ":" + df_["pos"].astype(str)),
        "block": block_lbl  # or use 'block_idx' if you prefer integers
    })
    return out

# function to check whether the run_id is correct
def check_runID(run_id):
    run_id = run_id.replace(' ', '')
    if run_id == '':
        messageError = True
        messageToUser = 'Please insert a valid Run ID'
    else:
        # check whether the file exists
        if os.path.exists('/Annotation/RUNS/AnnotateMe_results_%s.tar.gz' %(run_id)):
            messageError = False
            messageToUser = 'Valid Run ID. Download will start soon'
        else:
            messageError = True
            messageToUser = 'Not valid Run ID. Try again.'
    return messageError, messageToUser

# function to extract coordinates based on input type
def readBrowseOption(data_path, browse, window, refGen):
    # check if input is in the form 1:100000
    if ':' in browse and '-' not in browse:
        chrom, start_pos, end_pos, browse_type = int(browse.split(':')[0]), int(browse.split(':')[1]) - window, int(browse.split(':')[1]) + window, 'Single position'
    elif ':' in browse and '-' in browse:
        chrom, start_pos, end_pos, browse_type = int(browse.split(':')[0]), int(browse.split(':')[1].split('-')[0]) - window, int(browse.split(':')[1].split('-')[1]) + window, 'Interval'
    elif browse.lower().startswith("rs"):   # user picked an rsID
        conn = sqlite3.connect(f"{data_path}/databases/Genes/variants_info.sqlite")
        cur = conn.cursor()
        cur.execute("""
            SELECT chr_hg38, position, position_hg19
            FROM rsids
            WHERE rsid = ?
            LIMIT 1
        """, (browse,))
        row = cur.fetchone()
        conn.close()
        chrom, pos_hg38, pos_hg19 = row[0], row[1], row[2]
        if refGen == 'GRCh37':
            start_pos, end_pos = int(pos_hg19) - window, int(pos_hg19) + window
        else:
            start_pos, end_pos = int(pos_hg38) - window, int(pos_hg38) + window
        browse_type = 'RsID'
    else:
        genes_path = '%s/databases/Genes/genes_hg19.txt.gz' %(Path(data_path)) if refGen == 'GRCh37' else '%s/databases/Genes/genes_hg38.txt.gz' %(Path(data_path))
        line = list(os.popen('zgrep -i -w %s %s' % (browse.replace(' ', ''), genes_path.replace(' ', '\ '))))[0].strip()
        genes = re.split(r'[\t ]+', line)  # split on 1+ tabs or spaces
        if refGen == "GRCh37":
            chrom, start_pos, end_pos, browse_type = genes[2].replace('chr', ''), int(genes[4]) - window, int(genes[5]) + window, 'Gene'
        else:
            chrom, start_pos, end_pos, browse_type = genes[1].replace('chr', ''), int(genes[3]) - window, int(genes[4]) + window, 'Gene'
    return chrom, start_pos, end_pos, browse_type

# function to run annotation analysis
def run_annotation(my_list, inpType, refGeno, analType, gsea_source, qtl_tissues, email):
    # Sample a number for randomization
    random_number = random.randint(1, 10000)
    # Take snps and save them
    fpath_server = '/Annotation/RUNS/annotateMe_input_%s.txt' %(random_number)
    filename = 'annotateMe_input_%s.txt' %(random_number)
    with open(fpath_server, 'w') as finp:
        for x in my_list:
            finp.write('%s\n' %(x))
    # Define log file
    log_filename_server = '/Annotation/RUNS/annotateMe_run_%s.log' %(random_number)
    # Take input type
    ftype = str(inpType)
    if ftype == 'Colon':
        ftype = 1
    elif ftype == 'Tab':
        ftype = 2
    else:
        ftype = 3
    # Take analysis type
    analysis_type = str(analType)
    analysis_type = 'mapping' if analysis_type == 'Annot' else 'enrichment'
    print(gsea_source)
    analysis_mode = ','.join(list(gsea_source)) if analysis_type == 'enrichment' else 'None'
    # Take gtex tissues
    gtex_tissues = ','.join(list(qtl_tissues))
    # Then run annotate me externally in background -- this depends on the analysis_type requested
    command = "Rscript /Annotation/BIN/MAIN.R %s %s %s %s %s %s %s %s > %s" %(filename, ftype, email, analysis_type, analysis_mode, gtex_tissues, refGeno, random_number, log_filename_server)
    return command

# function to extract GTEx information
def get_gtex(data_path, genes, refGen):
    genename_index = 0 if refGen == 'GRCh37' else -2
    gtex_genes = []
    if len(genes) >0:
        for x in genes:
            try:
                tmp = [x.rstrip().split('\t') for x in list(os.popen('zgrep %s %s/databases/GTEx_Analysis_2017_06_05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz' %(x[genename_index], data_path.replace(' ', '\ '))))][0]
                gtex_genes.append(tmp)
            except:
                pass
    # convert to dataframe
    colnames = ["Name", "Description", "Adipose - Subcutaneous", "Adipose - Visceral (Omentum)", "Adrenal Gland", "Artery - Aorta", "Artery - Coronary", "Artery - Tibial", "Bladder", "Brain - Amygdala", "Brain - Anterior cingulate cortex (BA24)", "Brain - Caudate (basal ganglia)", "Brain - Cerebellar Hemisphere", "Brain - Cerebellum", "Brain - Cortex", "Brain - Frontal Cortex (BA9)", "Brain - Hippocampus", "Brain - Hypothalamus", "Brain - Nucleus accumbens (basal ganglia)", "Brain - Putamen (basal ganglia)", "Brain - Spinal cord (cervical c-1)", "Brain - Substantia nigra", "Breast - Mammary Tissue", "Cells - Cultured fibroblasts", "Cells - EBV-transformed lymphocytes", "Cervix - Ectocervix", "Cervix - Endocervix", "Colon - Sigmoid", "Colon - Transverse", "Esophagus - Gastroesophageal Junction", "Esophagus - Mucosa", "Esophagus - Muscularis", "Fallopian Tube", "Heart - Atrial Appendage", "Heart - Left Ventricle", "Kidney - Cortex", "Kidney - Medulla", "Liver", "Lung", "Minor Salivary Gland", "Muscle - Skeletal", "Nerve - Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)", "Small Intestine - Terminal Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole Blood"]
    gtex = pd.DataFrame(gtex_genes, columns=colnames)
    return gtex

# function to liftover sumstats
def liftover_sumstats(sumstats, from_build, to_build):
    if from_build == 'GRCh37':
        lifter = get_lifter('hg19', 'hg38')
    else:
        lifter = get_lifter('hg38', 'hg19')
    # helped function
    def lift_row(row):
        chrom = f"chr{row['chr']}"   # liftover expects 'chrN'
        pos   = row['position']
        res = lifter[chrom][pos]
        if res:   # non-empty result
            return res[0][1]   # lifted position
        else:
            return None 
    # apply to all rows
    sumstats['positions_lifted'] = sumstats.apply(lift_row, axis=1)
    return sumstats

# function to get data to plot given the gwas name, chromosome, start and end position
def get_data_plot(data_path, gwas, chrom, start_pos, end_pos, refGen):
    # define region to look at
    if refGen == 'GRCh37':
        # set liftover
        liftover_info = get_lifter('hg19', 'hg38')
        region = '%s:%s-%s' %(str(chrom).replace('chr', ''), str(start_pos), str(end_pos))
        region_hg38 = '%s:%s-%s' %(str(chrom).replace('chr', ''), str(liftover_info[chrom][start_pos][0][1]), str(liftover_info[chrom][end_pos][0][1]))
    else:
        # set liftover
        liftover_info = get_lifter('hg38', 'hg19')
        region_hg38 = '%s:%s-%s' %(str(chrom), str(start_pos), str(end_pos))
        region = '%s:%s-%s' %(str(chrom), str(liftover_info[chrom][start_pos][0][1]), str(liftover_info[chrom][end_pos][0][1]))
    # get gwas info
    data_list = pd.DataFrame()
    data_list_info = []
    sumstats_build = ''
    # iterate over gwas
    for gw in gwas:
        # extract gwas info from label -- label is "Alzheimer's disease — Bellenguez C (2022) [ebi-a-GCST90027158]" --> extract ebi-a-GCST90027158
        gwas_id = gw.split('[')[-1].replace(']', '').strip()
        # get gwas info
        gwas_info = meta[gwas_id]
        # check whether the gwas is mapped to hg19 or hg38
        if 'hg19' in gwas_info.loc['build'].lower():
            sumstats_build = 'GRCh37'
            # get region of interest with hg19 coordinates
            res = query.associations(id=[gwas_id], variant=[region], proxies=0)
        else:
            sumstats_build = 'GRCh38'
            # get region of interest with hg38 coordinates
            res = query.associations(id=[gwas_id], variant=[region_hg38], proxies=0)
        # put results into a dataframe
        sumstats = pd.DataFrame(res)
        # do liftover -- it's fast and useful to have both positions
        sumstats = liftover_sumstats(sumstats, sumstats_build, refGen)
        # update column names
        if sumstats_build == 'GRCh37':
            if refGen == 'GRCh37':
                sumstats = sumstats.rename(columns={"chr": "Chrom", "position": "Position", "p": "Pvalue", "rsid": "Rsid", "eaf": "EAF", "ea": "EA", "nea": "NEA", "positions_lifted": "Position_hg38", "id": "ID", "trait": "Gwas", "beta": "Beta", "se": "SE", 'n': 'N'}, errors='ignore')
            else:
                sumstats = sumstats.rename(columns={"chr": "Chrom", "position": "Position_hg19", "p": "Pvalue", "rsid": "Rsid", "eaf": "EAF", "ea": "EA", "nea": "NEA", "positions_lifted": "Position", "id": "ID", "trait": "Gwas", "beta": "Beta", "se": "SE", 'n': 'N'}, errors='ignore')
        else:
            if refGen == 'GRCh37':
                sumstats = sumstats.rename(columns={"chr": "Chrom", "position": "Position_hg38", "p": "Pvalue", "rsid": "Rsid", "eaf": "EAF", "ea": "EA", "nea": "NEA", "positions_lifted": "Position", "id": "ID", "trait": "Gwas", "beta": "Beta", "se": "SE", 'n': 'N'}, errors='ignore')
            else:
                sumstats = sumstats.rename(columns={"chr": "Chrom", "position": "Position", "p": "Pvalue", "rsid": "Rsid", "eaf": "EAF", "ea": "EA", "nea": "NEA", "positions_lifted": "Position_hg19", "id": "ID", "trait": "Gwas", "beta": "Beta", "se": "SE", 'n': 'N'}, errors='ignore')
        # append to data_list
        data_list = pd.concat([data_list, sumstats], ignore_index=True) if data_list.shape[0] >0 else sumstats
        data_list_info.append(gwas_info)
    # make position and pvalue numeric
    data_list["Position"] = pd.to_numeric(data_list['Position'], errors='coerce')
    data_list["Pvalue"] = pd.to_numeric(data_list['Pvalue'], errors='coerce')
    # find smallest nonzero p-value
    min_nonzero = data_list.loc[data_list["Pvalue"] > 0, "Pvalue"].min()
    # replace zeros with slightly smaller random values
    n_zeros = (data_list["Pvalue"] == 0).sum()
    if n_zeros > 0:
        data_list.loc[data_list["Pvalue"] == 0, "Pvalue"] = min_nonzero * np.random.uniform(0.1, 1.0, n_zeros)
    # then take the -log10 of the pvalue
    data_list['Pvalue'] = -np.log10(data_list['Pvalue'])
    # then take the position value as megabases
    data_list["Position_plot"] = data_list["Position"] / 1000000
    # drop duplicates based on position and gwas
    data_list = data_list.drop_duplicates(subset=['ID', 'Position'])
    # combine data_list_info into a transposed dataframe
    data_list_info_df = pd.DataFrame(data_list_info)
    return data_list, data_list_info_df

# function to extract genes to plot
def extract_genes(data_path, chrom, start_pos, end_pos, refGen):
    # set prefix of the reference used
    refPrefix = 'hg19' if refGen == 'GRCh37' else 'hg38'
    # use tabix to find genes -- enlarge window by 50kb up and down
    cmd = 'tabix %s/databases/Genes/genes_%s.txt.gz chr%s:%s-%s' %(data_path.replace(' ', '\ '), refPrefix, str(chrom).replace('chr', ''), str(start_pos - 1000), str(end_pos + 1000))
    genes = [x.rstrip().split('\t') for x in os.popen(cmd)]
    # in case of GRCh38, there may be multiple entries for the same gene -- keep the one where the gene is largest
    if refGen == 'GRCh38':
        unique_genes = {}
        for g in genes:
            key = g[-1]
            if key not in unique_genes:
                unique_genes[key] = g
            else:
                if g[0] > unique_genes[key][0]:
                    unique_genes[key] = g
        genes = list(unique_genes.values())
    # exclude long non coding stuff
    if refGen == 'GRCh37':
        genes = [x for x in genes if 'LOC' not in x[0]]
    else:
        genes = [x for x in genes if 'LOC' not in x[-1]]
    # add y axis
    genes = [[*sublist, i] for i, sublist in enumerate(genes)]
    return genes

# function to fix svs names for subsetting
def fix_svs_names(svtypes):
    svtypes_fixed = []
    for sv in svtypes:
        if sv == 'tr':
            svtypes_fixed.append('Tandem Repeat')
            svtypes_fixed.append('Satellite')
        elif sv == 'sine':
            svtypes_fixed.append('SINE')
        elif sv == 'line':
            svtypes_fixed.append('LINE')
        elif sv == 'ltr':
            svtypes_fixed.append('LTR')
    return svtypes_fixed            

# function to extract genes to plot
def extract_sv(data_path, chrom, start_pos, end_pos, refGen, svtypes):
    # fix sv names for subsetting
    svtypes = fix_svs_names(svtypes)
    # set reference prefix
    refPrefix = 'hg19' if refGen == 'GRCh37' else 'hg38'
    # use tabix to find genes -- enlarge window by 50kb up and down
    cmd = 'tabix %s/databases/Structural_variants/harmonized_svs_%s.bed.gz chr%s:%s-%s' %(data_path.replace(' ', '\ '), refPrefix, str(chrom), str(start_pos), str(end_pos))
    svs = [x.rstrip().split('\t') for x in os.popen(cmd)]
    # select based on the input selected
    if 'all' in svtypes:
        pass
    else:
        svs = [x for x in svs if x[0] in svtypes]
    # add y axis
    svs = [[*sublist, i] for i, sublist in enumerate(svs)]
    # convert to dataframe as well
    colnames = ['repClass', 'chrom', 'start', 'end', 'len', 'repName', 'repFamily', 'color', 'y']
    svs_df = pd.DataFrame(svs, columns=colnames)
    return svs, svs_df

# function to extract recombination rates
def extract_recomb(data_path, chrom, start_pos, end_pos, refGen):
    # set reference prefix
    refPrefix = 'hg19' if refGen == 'GRCh37' else 'hg38'
    # use tabix to extract recombination rates
    cmd = 'tabix %s/databases/Recombination_rates/recombination_rates_%s.txt.gz chr%s:%s-%s' %(data_path.replace(' ', '\ '), refPrefix, str(chrom), str(start_pos), str(end_pos))
    rec_rates = [x.rstrip().split('\t') for x in os.popen(cmd)]
    # convert to dataframe
    rec_rates_df = pd.DataFrame(rec_rates, columns=['chrom', 'position', 'rate', 'map'])
    rec_rates_df['position'] = pd.to_numeric(rec_rates_df['position'])/1000000
    rec_rates_df['rate'] = pd.to_numeric(rec_rates_df['rate'])
    return rec_rates_df

# function to extract GWAS catalog information
def extract_gwascatalog(data_path, chrom, start_pos, end_pos, refGen):
    # set reference prefix
    refPrefix = 'hg19' if refGen == 'GRCh37' else 'hg38'
    # use tabix to extract recombination rates
    cmd = 'tabix %s/databases/GWAS_catalog/Gwas_catalog_%s.txt.gz %s:%s-%s' %(data_path.replace(' ', '\ '), refPrefix, str(chrom), str(start_pos), str(end_pos))
    gwascat = [x.rstrip().split('\t') for x in os.popen(cmd)]
    # convert to dataframe
    colnames = ['chr', 'pos', 'pos_hg38', 'snp', 'gene', 'trait', 'pubmed_id', 'pvalue'] if refGen == 'GRCh37' else ['chr', 'pos_hg37', 'pos', 'snp', 'gene', 'trait', 'pubmed_id', 'pvalue']
    gwascat = pd.DataFrame(gwascat, columns=colnames)
    # prepare table for output: add locus, remove redundant columns
    gwascat['Locus'] = gwascat.apply(lambda x: str(x['chr']) + ':' + str(x['pos']), axis=1)
    gwascat = gwascat[['Locus', 'snp', 'gene', 'trait', 'pubmed_id', 'pvalue']]
    gwas_table = gwascat.to_dict(orient='records')
    return gwas_table, gwascat

# function to draw the scatterplot given a dataframe and positions to plot
def scatterplot(df, chrom, start_pos, end_pos, gwas, genes, svs, browse_type, refGen, recomb_data, exons, browse, plotype):
    # define structure for the plot
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14, 12), gridspec_kw={'height_ratios': [2, 1.5, 1.5]})
    ##################
    # Scatterplot here
    # define colors and map them to datasets
    all_colors = ['red', 'blue', 'orange', 'yellow', 'grey', 'green']
    #all_colors = [x for x in list(webcolors.CSS3_NAMES_TO_HEX.keys()) if x != 'white']
    random.shuffle(all_colors)
    gwas_unique = df['Gwas'].unique().tolist()
    # sometimes using the GWAS-Cat data, there are too many colors needed and it crashes. In these cases, take the top 100 traits based on significance
    if len(gwas_unique) > 70:
        df_gwascat = df[df['Gwas'].isin(gwas_unique)]
        df_gwascat = df_gwascat[~df_gwascat['Gwas'].isin(gwas)]
        df_NOgwascat = df[df['Gwas'].isin(gwas)]
        df_gwascat = df[df['Gwas'].isin(gwas_unique)]
        # sort the dataframe by the numeric column in descending order
        df_sorted = df_gwascat.sort_values(by='Pvalue', ascending=False)
        # select the top 100 unique values from the character column
        unique_values = df_sorted['Gwas'].unique()[:70]
        # subset the dataframe based on the selected values
        subset = df_sorted[df_sorted['Gwas'].isin(unique_values)]
        # re-join data
        df = pd.concat([df_NOgwascat, subset])
        gwas_unique = df['Gwas'].unique().tolist()
    color_map = {val: all_colors[index] for index, val in enumerate(gwas_unique)}
    # Set subplot background color to grey
    ax1.set_facecolor('#f0f0f0')
    ax2.set_facecolor('#f0f0f0')
    ax3.set_facecolor('#f0f0f0')
    # Then put a white vertical grid
    minimum_value_pos = df['Position_plot'].min()
    maximum_value_pos = df['Position_plot'].max()
    num_values = 10     # define grid size
    equally_spaced_values = np.linspace(minimum_value_pos, maximum_value_pos, num_values)
    for val in equally_spaced_values:
        ax1.axvline(x=val, color='white', linestyle='-', linewidth=1.5, zorder = 1)
        ax2.axvline(x=val, color='white', linestyle='-', linewidth=1.5, zorder = 1)
        ax3.axvline(x=val, color='white', linestyle='-', linewidth=1.5, zorder = 1)
    # plot data depending on plot type
    if plotype == 'Scatter':
        colors = df['Gwas'].map(color_map)
        # scatter plot
        ax1.scatter(df['Position_plot'], df['Pvalue'], c=colors, zorder = 2)
    else:
        # smooth plot - need to iterate over gwas types
        for gw in gwas:
            if gw == 'GWAS-Cat':
                # if the whole GWAS-Catalog was selected, still use scatter as it is too sparse to do densities
                tmp = df.loc[~df['Gwas'].isin(gwas)].copy()
                colors = tmp['Gwas'].map(color_map)
                ax1.scatter(tmp['Position_plot'], tmp['Pvalue'], c=colors, zorder = 2)
            else:
                tmp = df.loc[df['Gwas'] == gw].copy()
                # use rolling function to select local maxima
                tmp['smoothed'] = tmp['Pvalue'].rolling(window=20).max()
                # fit lowess
                lowess = sm.nonparametric.lowess(tmp['smoothed'], tmp['Position_plot'], frac=0.02)
                # fill the area under the lowess line
                ax1.fill_between(lowess[:, 0], lowess[:, 1], color=color_map[gw], alpha=0.4)
                ax1.plot(lowess[:, 0], lowess[:, 1], color=color_map[gw], linewidth=3)
    # set y-lim -- max pvalue + 10% of the max value
    max_p_assoc = int(max(df['Pvalue'])) + int(max(df['Pvalue']))*0.08
    ax1.set_ylim(0, max(max_p_assoc, 6))
    # if rsid was requested, the corresponding marker will be in a different color (black)
    if browse_type == 'RsID':
        subset_rsid = df.loc[df['Rsid'] == browse.lower()].copy()
        # add colors
        colors_rsid = subset_rsid['Gwas'].map(color_map)
        # also label snp of interest
        subset_rsid['Position_label_x'] = subset_rsid['Position_plot'] - (end_pos - start_pos)/1000000*0.25
        subset_rsid['Position_label_y'] = subset_rsid['Pvalue']
        subset_rsid['Labels'] = 'Input SNP ~ ' + subset_rsid['Gwas']
        ax1.scatter(subset_rsid['Position_plot'], subset_rsid['Pvalue'], c=colors_rsid, marker = '^')
        subset_rsid_lst = subset_rsid.values.tolist()
        if len(gwas_unique) <5:
            for x in subset_rsid_lst:
                ax1.annotate(x[-1], xy=(x[-4], x[-2]), xytext=(x[-3], x[-2]), arrowprops=dict(facecolor='black', shrink=0.002, width=1, mutation_scale=10))
    # add text for rsid of the most significant variant: the selection is different in case of the GWAS-Cat: in this case we take only the top 3 traits
    top_snp = df.loc[df.groupby('Gwas')['Pvalue'].idxmax(), ['Position_plot', 'Pvalue', 'Gwas', 'Rsid']]
    if 'GWAS-Cat' in gwas:
        top_snp_noCat = top_snp[top_snp['Gwas'].isin(gwas)]
        top_snp_Cat = top_snp.sort_values('Pvalue', ascending=False)
        top_snp_Cat = top_snp_Cat.head(3)
        top_snp = pd.concat([top_snp_Cat, top_snp_noCat], axis=0)
    top_snp['Position_label_x'] = top_snp['Position_plot'] + (end_pos - start_pos)/1000000*0.04; top_snp['Position_label_y'] = top_snp['Pvalue']
    top_snp['Labels'] = top_snp['Rsid'] + ' ~ ' + top_snp['Gwas']; top_snp_lst = top_snp.values.tolist()
    for x in top_snp_lst:
        if browse_type != 'Rsid':
            ax1.annotate(x[-1], xy=(x[0], x[1]), xytext=(x[4], x[1]), arrowprops=dict(facecolor='black', shrink=0.002, width=1, mutation_scale=10))
        else:
            if x[3] not in str(subset_rsid['Rsid']):
                ax1.annotate(x[-1], xy=(x[0], x[1]), xytext=(x[4], x[1]), arrowprops=dict(facecolor='black', shrink=0.002, width=1, mutation_scale=10))
    # if requested, plot recombination rates
    if isinstance(recomb_data, pd.DataFrame):
        ax1_2 = ax1.twinx()
        ax1_2.plot(recomb_data['position'], recomb_data['rate'], c='green')
        ax1_2.set_ylabel('Recombination Rates (cM/Mb)', color='green')
        ax1_2.set_ylim(-0.5, 100)
    # add title
    title = '%s:%s-%s Input:%s (%s)' %(str(chrom), str(start_pos), str(end_pos), browse_type, refGen)
    ax1.set_title(title)
    # axes names
    ax1.set_xlabel("")
    ax1.set_ylabel("-log10(P-value)")
    # disable x-axis values
    ax1.set_xticks([])
    # Add the legend -- if the number of characters in the legend name is >20, then put a . and truncate it
    for key in list(color_map.keys()):
        if len(key) > 20:
            new_key = key[:20] + "."
            color_map[new_key] = color_map.pop(key)
    handles = [plt.plot([],[],color=color_map[color], marker="o", ls="", label=color)[0] for color in color_map.keys()]
    # plot the legend outside the plot (on top). Use columns containing 7 values. If more that 28 elements need to be plotted, make it extra small
    if len(color_map.keys()) >30:
        legend = ax1.legend(handles=handles, loc = 'upper center', bbox_to_anchor=(0.5, 1.45), fontsize=5.4, ncol = math.ceil(len(color_map.keys())/10))
    elif len(color_map.keys()) >15:
        legend = ax1.legend(handles=handles, loc = 'upper center', bbox_to_anchor=(0.5, 1.45), fontsize=7, ncol = math.ceil(len(color_map.keys())/8))
    elif len(color_map.keys()) >=10:
        legend = ax1.legend(handles=handles, loc = 'upper center', bbox_to_anchor=(0.5, 1.45), fontsize=8.5, ncol = math.ceil(len(color_map.keys())/8))
    else:
        legend = ax1.legend(handles=handles, loc = 'upper center', bbox_to_anchor=(0.5, 1.35), fontsize=10, ncol = math.ceil(len(color_map.keys())/8))
    ##################
    # Gene track here
    max_y = max(g[-1] for g in genes) if len(genes) >0 else 1
    index_pos = 4 if refGen == 'GRCh37' else 5
    index_text = 0 if refGen == 'GRCh37' else -2
    index_exon = 9 if refGen == 'GRCh37' else 8
    if len(genes) >0:
        for g in genes:
            ax2.plot([int(g[index_pos])/1000000, int(g[index_pos + 1])/1000000], [int(g[-1]), int(g[-1])], 'r-', linewidth=3)
            if len(genes) < 10:
                ax2.text(int(g[index_pos])/1000000 + (int(g[index_pos + 1])/1000000 - int(g[index_pos])/1000000)/2, int(g[-1])+max_y*0.06, g[index_text], ha='center', fontsize = 10)
            elif len(genes) <20:
                ax2.text(int(g[index_pos])/1000000 + (int(g[index_pos + 1])/1000000 - int(g[index_pos])/1000000)/2, int(g[-1])+max_y*0.03, g[index_text], ha='center', fontsize = 7)
            # show exons if requested
            if exons == 'Yes':
                # extract exon starts and end
                exons_start = [int(x)/1000000 for x in g[index_exon].split(',') if x != '']
                exons_end = [int(x)/1000000 for x in g[index_exon+1].split(',') if x != '']
                ax2.plot([exons_start, exons_end], [int(g[-1]), int(g[-1])], 'r-', linewidth=5)
    else:
        ax2.text((start_pos + (end_pos-start_pos)/2)/1000000, 0.5, "No Genes to plot", ha='center')
    # axes names
    ax2.set_xlabel("")
    ax2.set_ylabel("Gene track")
    # disable axes values
    ax2.set_xticks([])
    ax2.set_yticks([])
    # set x-limit
    ax2.set_xlim(int(start_pos)/1000000, int(end_pos)/1000000)
    ax2.set_ylim(0-max_y*0.20, max_y+max_y*0.20)
    ##################
    # Structural variants
    max_y = max(sv[-1] for sv in svs) if len(svs) >0 else 1
    index_col = 4 if refGen == 'GRCh37' else 5
    index_sv = 3 if refGen == 'GRCh37' else 4
    index_diff = -3 if refGen == 'GRCh37' else -5
    if len(svs) >0:
        for sv in svs:
            # if the size is too small (<5% of the window), then plot a square
            if int(sv[index_diff]) < (end_pos - start_pos)*0.001:
                ax3.plot(int(sv[1])/1000000, sv[-1], color=sv[index_col], marker = 's', markersize = 5)
            else:
                ax3.plot([int(sv[1])/1000000, int(sv[2])/1000000], [sv[-1], sv[-1]], color=sv[index_col], linewidth=5)
    else:
        ax3.text((start_pos + (end_pos-start_pos)/2)/1000000, 0.5, "No Structural variants to plot", ha='center')
    # axes names
    ax3.set_xlabel("Genomic Position (Mb)")
    ax3.set_ylabel("Structural variants")
    # disable axes values
    ax3.set_yticks([])
    # set x-limit
    ax3.set_xlim(int(start_pos)/1000000, int(end_pos)/1000000)
    ax3.set_ylim(0-max_y*0.10, max_y+max_y*0.20)
    # Add the legend
    color_map_sv = {x[index_sv]: x[index_col] for x in svs}
    handles_sv = [plt.plot([],[],color=color_map_sv[color], marker="s", ls="", label=color)[0] for color in color_map_sv.keys()]
    if len(svs) >0:
        ax3.legend(handles=handles_sv, loc = 'upper center', ncol=6, fontsize=8, markerscale=0.75)
    #plt.show()         # enable this for debuggin
    # save image
    img = io.BytesIO()
    fig.savefig(img, format="png")
    img.seek(0)
    return img

# function to draw the scatterplot plotly given a dataframe and positions to plot
def scatterplot_plotly(df, chrom, start_pos, end_pos, gwas, genes, svs, browse_type, refGen, recomb_data, exons, browse, plotype, ld, ld_df, yrange, haplo_df, plotHaplo):
    # ------- colors like your Matplotlib version -------
    all_colors = ['red', 'blue', 'orange', 'grey', 'green', 'purple', 'brown', 'pink', 'olive', 'cyan', 'magenta', 'teal', 'lavender', 'turquoise', 'tan', 'salmon', 'gold', 'lightgreen', 'lightblue']
    rng = np.random.default_rng(42)
    all_colors = list(all_colors)
    rng.shuffle(all_colors)

    gwas_unique = df['Gwas'].unique().tolist()
    if len(gwas_unique) > 70:
        df_NOgwascat = df[df['Gwas'].isin(gwas)]
        df_gwascat = df[df['Gwas'].isin(gwas_unique)]
        df_sorted_tmp = df_gwascat.sort_values(by='Pvalue', ascending=False)
        unique_values = df_sorted_tmp['Gwas'].unique()[:70]
        subset = df_sorted_tmp[df_sorted_tmp['Gwas'].isin(unique_values)]
        df = pd.concat([df_NOgwascat, subset], ignore_index=True)
        gwas_unique = df['Gwas'].unique().tolist()
    
    # color map per dataset (cycle if more datasets than colors)
    if len(gwas_unique) > len(all_colors):
        reps = math.ceil(len(gwas_unique) / len(all_colors))
        palette = (all_colors * reps)[:len(gwas_unique)]
    else:
        palette = all_colors[:len(gwas_unique)]
    color_map = dict(zip(gwas_unique, palette))

    # ---------- figure with 4 rows ----------
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.02, row_heights=[0.45, 0.275, 0.275], specs=[[{"secondary_y": True}], [{}], [{}]])

    # background & white-ish grid feel
    for r in (1, 2, 3):
        fig.update_xaxes(showgrid=True, gridcolor="#e0e0e0", zeroline=False, row=r, col=1)
        fig.update_yaxes(showgrid=True, gridcolor="#e0e0e0", zeroline=False, row=r, col=1)
    fig.update_layout(plot_bgcolor="white", paper_bgcolor="white")

    # check if LD should be plotted -- ld should be Yes, and len(gwas_unique) should be 1
    if ld == 'Yes' and len(gwas_unique) == 1 and not ld_df.empty:
        plot_ld = True
    else:
        plot_ld = False
    
    # check if haplotype blocks should be plotted
    if plotHaplo == "Yes":
        print('not none')
        plot_haplo = True
        # merge with df
        df = df.merge(haplo_df, how='left', left_on='Rsid', right_on='variant')
    else:
        plot_haplo = False

    # ---------- panel 1: associations ----------
    if plotype == 'Scatter':
        hover_cut = -math.log10(5e-3)
        for gw in gwas_unique:
            # set symbol
            symbol = 'circle'
            tmp = df[df['Gwas'] == gw]
            if tmp.empty:
                continue
            # if LD is requested, then color by r2
            if plot_ld:
                # combine ld info with tmp
                if refGen == 'GRCh37':
                    tmp = tmp.merge(ld_df[['pos', 'r2', 'color']], how='left', left_on='Position_hg38', right_on='pos')
                else:
                    tmp = tmp.merge(ld_df[['pos', 'r2', 'color']], how='left', left_on='Position', right_on='pos')
                tmp['base_color'] = list(map(normalize_color, tmp['color'].fillna('grey')))
            else:
                tmp['base_color'] = normalize_color(color_map[gw])

            # if haplotype blocks are requested, then color by haplotype block
            if plot_haplo:
                # set all symbols
                all_symbols = ['circle', 'square', 'diamond', 'cross', 'x', 'triangle-up', 'triangle-down', 'pentagon', 'hexagon', 'star']
                tmp["base_color"] = "#CCCCCC"
                # assign a color per haplotype block, otherwise grey
                haplo_colors = {}
                unique_haplos = tmp['block'].dropna().unique().tolist()
                for i, hb in enumerate(unique_haplos):
                    haplo_colors[hb] = all_colors[i % len(all_colors)]
                # apply colors
                tmp['base_color'] = tmp.apply(lambda x: normalize_color(haplo_colors[x['block']]) if pd.notna(x['block']) else x['base_color'], axis=1)
                # define symbol
                symbol = all_symbols[gwas_unique.index(gw)]

            # hover points (significant)
            sig = tmp[tmp['Pvalue'] > hover_cut]
            if not sig.empty:
                fig.add_trace(
                    go.Scattergl(
                        x=sig['Position_plot'], 
                        y=sig['Pvalue'],
                        mode='markers',
                        marker=dict(color=sig['base_color'], size=12, symbol=symbol),
                        name=gw,
                        legendgroup=gw,
                        hovertemplate=(
                            "<b>%{customdata[0]}</b><br>"
                            "Gwas: %{customdata[1]}<br>"
                            "Pos: %{customdata[2]}<br>"
                            "-log10(p): %{y:.2f}<br>"
                            "Beta: %{customdata[3]}<br>"
                            "EA: %{customdata[4]}<br>"
                            "NEA: %{customdata[5]}<extra></extra>"
                        ),
                        customdata=np.stack([
                            sig['Rsid'], 
                            sig['Gwas'], 
                            sig['Position'], 
                            sig['Beta'], 
                            sig['EA'], 
                            sig['NEA']
                        ], axis=1),
                    ),
                    row=1, col=1
                )

            # non-hover points
            nonsig = tmp[tmp['Pvalue'] <= hover_cut]
            if not nonsig.empty:
                fig.add_trace(
                    go.Scattergl(
                        x=nonsig['Position_plot'], y=nonsig['Pvalue'],
                        mode='markers',
                        marker=dict(color=nonsig['base_color'], size=8, opacity=0.8),
                        name=gw, legendgroup=gw, showlegend=False, hoverinfo='skip'
                    ), row=1, col=1)
    else:
        # Smooth mode: GWAS-Cat (others) as scatter, selected GWAS as LOWESS
        if 'GWAS-Cat' in gwas:
            tmp = df.loc[~df['Gwas'].isin(gwas)].copy()
            if not tmp.empty:
                hover_cut = -math.log10(5e-3)
                sig = tmp[tmp['Pvalue'] > hover_cut]
                if not sig.empty:
                    colors = normalize_color(list(sig['Gwas'].map(color_map)))
                    fig.add_trace(
                        go.Scattergl(
                            x=sig['Position_plot'], y=sig['Pvalue'],
                            mode='markers',
                            marker=dict(color=colors, size=6),
                            name='GWAS-Cat (sig)', legendgroup='GWAS-Cat',
                            hovertemplate="<b>%{customdata[0]}</b><br>"
                                          "Gwas: %{customdata[1]}<br>"
                                          "Pos: %{x}<br>"
                                          "-log10(p): %{y:.2f}<extra></extra>",
                            customdata=np.stack([sig['Rsid'], sig['Gwas']], axis=1),
                        ),
                        row=1, col=1
                    )
                nonsig = tmp[tmp['Pvalue'] <= hover_cut]
                if not nonsig.empty:
                    colors = normalize_color(list(nonsig['Gwas'].map(color_map)))
                    fig.add_trace(
                        go.Scattergl(
                            x=nonsig['Position_plot'], y=nonsig['Pvalue'],
                            mode='markers',
                            marker=dict(color=colors, size=5, opacity=0.8),
                            name='GWAS-Cat', legendgroup='GWAS-Cat',
                            showlegend=False, hoverinfo='skip'
                        ),
                        row=1, col=1
                    )

        for gw in gwas_unique:
            tmp = df[df['Gwas'] == gw].copy()
            if tmp.empty:
                continue
            tmp['smoothed'] = tmp['Pvalue'].rolling(window=20).max()
            low = sm.nonparametric.lowess(tmp['smoothed'], tmp['Position_plot'], frac=0.02, return_sorted=True)
            color = normalize_color(color_map.get(gw, 'blue'))
            rgb = mcolors.to_rgb(color)
            fillcolor = f"rgba({rgb[0]},{rgb[1]},{rgb[2]},0.5)"
            if len(low) > 1:
                xlow, ylow = low[:, 0], low[:, 1]
                fig.add_trace(
                    go.Scatter(
                        x=xlow, y=ylow,
                        mode='lines',
                        line=dict(color=color, width=3),
                        fill='tozeroy',
                        fillcolor=fillcolor,
                        opacity=0.25,
                        name=gw, legendgroup=gw,
                        hovertemplate=(
                            "Position: %{x}<br>"
                            "Significance: %{y:.2f}<extra></extra>"
                        )
                    ),
                    row=1, col=1
                )

    # y-range
    if isinstance(yrange, list):
        min_y = yrange[0]
        max_y = yrange[1] + math.ceil(yrange[1]) * 0.10
        fig.update_yaxes(range=[min_y, max_y], row=1, col=1, title_text="-log10(P-value)", title_font=dict(size=18, family="Arial", color="black", weight="bold"), tickfont=dict(size=14, family="Arial", color="black"))
    else:
        max_p = float(np.nanmax(df['Pvalue'])) if len(df) else 6.0
        max_p_assoc = math.ceil(max_p) + math.ceil(max_p) * 0.10
        fig.update_yaxes(range=[0, max(max_p_assoc, 6)], row=1, col=1, title_text="-log10(P-value)", title_font=dict(size=18, family="Arial", color="black", weight="bold"), tickfont=dict(size=14, family="Arial", color="black"))

    # rsID highlight
    if browse_type == 'RsID':
        subset_rsid = df[df['Rsid'].str.lower() == str(browse).lower()].copy()
        if not subset_rsid.empty:
            colors = normalize_color(list(subset_rsid['Gwas'].map(color_map)))
            fig.add_trace(
                go.Scatter(
                    x=subset_rsid['Position_plot'],
                    y=subset_rsid['Pvalue'],
                    mode='markers',
                    marker=dict(size=12, color=normalize_color('black')),
                    name="Input SNP",
                    legendgroup="Input SNP",
                    hovertemplate=(
                        "<b>%{customdata[0]}</b><br>"
                        "Gwas: %{customdata[1]}<br>"
                        "Pos: %{customdata[2]}<br>"
                        "-log10(p): %{y:.2f}<br>"
                        "Beta: %{customdata[3]}<br>"
                        "EA: %{customdata[4]}<br>"
                        "NEA: %{customdata[5]}<extra></extra>"
                    ),
                    customdata=np.stack([
                        subset_rsid['Rsid'],
                        subset_rsid['Gwas'],
                        subset_rsid['Position'],
                        subset_rsid['Beta'],
                        subset_rsid['EA'],
                        subset_rsid['NEA']
                    ], axis=1),
                ),
                row=1, col=1
            )

    # top SNP annotations
    if not df.empty:
        top_idx = df.groupby('Gwas')['Pvalue'].idxmax()
        top_snp = df.loc[top_idx, ['Position_plot', 'Pvalue', 'Gwas', 'Rsid']].copy()
        if 'GWAS-Cat' in gwas:
            top_snp_noCat = top_snp[top_snp['Gwas'].isin(gwas)]
            top_snp_Cat = top_snp.sort_values('Pvalue', ascending=False).head(3)
            top_snp = pd.concat([top_snp_Cat, top_snp_noCat], axis=0)

        for _, row in top_snp.iterrows():
            if browse_type == 'RsID' and str(row['Rsid']).lower() == str(browse).lower():
                continue
            fig.add_annotation(
                x=row['Position_plot'], y=row['Pvalue'],
                xanchor='left', yanchor='middle',
                ax=40, ay=0, showarrow=True, arrowwidth=1, arrowhead=2,
                text=f"{row['Rsid']} ~ {row['Gwas']}",
                font=dict(size=14, family="Arial", color="black", weight="bold"),
                row=1, col=1
            )

    # recombination secondary y-axis
    if isinstance(recomb_data, pd.DataFrame) and not recomb_data.empty:
        fig.add_trace(
            go.Scatter(
                x=recomb_data['position'], y=recomb_data['rate'],
                mode='lines', line=dict(width=2),
                name="Recombination", legendgroup="Recombination",
                hovertemplate="Pos: %{x}<br>Rate: %{y:.2f} cM/Mb<extra></extra>"
            ),
            row=1, col=1, secondary_y=True
        )
        fig.update_yaxes(title_text="Recomb (cM/Mb)", range=[-0.5, 100], secondary_y=True, row=1, col=1, title_font=dict(size=18, family="Arial", color="green", weight="bold"), tickfont=dict(size=14, family="Arial", color="green"))

    # title & layout
    #title = f"{chrom}:{start_pos}-{end_pos} Input:{browse_type} ({refGen})"
    fig.update_layout(
        legend_orientation="h",
        legend_yanchor="bottom",
        legend_y=1.18,
        legend_x=0.5,
        legend_xanchor="center",
        margin=dict(l=40, r=20, t=90, b=40),
        hovermode="x unified",
        autosize=True
    )

    # ---------- panel 2: gene track ----------
    if genes:
        index_pos  = 4 if refGen == 'GRCh37' else 5
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
            row = np.array([str(g[index_text]), str(g[index_transc]), str(g[index_strand]), int(g[index_pos]), int(g[index_pos + 1]), int(g[index_exonN])], dtype=object)
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
                    showlegend=False), row=2, col=1)
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
                ), row=2, col=1)

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
                x=x0 + (x1 - x0) / 2, y=ylane + yoff,
                text=label, showarrow=False, font=dict(size=fsize), hovertext=hovertext, row=2, col=1)

            if exons == 'Yes':
                starts = [int(x)/1_000_000 for x in str(g[index_exon]).split(',') if x]
                ends   = [int(x)/1_000_000 for x in str(g[index_exon+1]).split(',') if x]
                for xs, xe in zip(starts, ends):
                    fig.add_trace(
                        go.Scatter(
                            x=[xs, xe], y=[ylane, ylane],
                            mode='lines',
                            line=dict(color=normalize_color(color), width=12),
                            hoverinfo='skip', showlegend=False
                        ),
                        row=2, col=1
                    )
        fig.update_yaxes(range=[-max_y*0.20, max_y+max_y*0.20], row=2, col=1)

    else:
        fig.add_annotation(
            x=(start_pos + (end_pos-start_pos)/2)/1_000_000,
            y=0.5, text="No Genes to plot", showarrow=False,
            row=2, col=1
        )
    fig.update_yaxes(title_text="Gene track", row=2, col=1, title_font=dict(size=18, family="Arial", color="black", weight="bold"), tickfont=dict(size=14, family="Arial", color="white"))

    # ---------- panel 3: structural variants ----------
    color_map_sv = {}
    if svs:
        index_start = 2
        index_end = 3
        index_col = 7
        index_sv = 0
        index_name = 5
        index_fam = 6
        index_diff = 4
        max_y_sv = max(sv[-1] for sv in svs) if len(svs) else 1

        # based on the size, split into points and segments
        by_type_points, by_type_segments = {}, {}
        for sv in svs:
            svtype = sv[index_sv]
            color  = normalize_color(sv[index_col])
            color_map_sv[svtype] = color
            ylane = sv[-1]
            size_bp = int(sv[index_diff])
            x0 = int(sv[index_start]) / 1_000_000
            x1 = int(sv[index_end]) / 1_000_000
            if size_bp < (end_pos - start_pos) * 0.01:
                by_type_points.setdefault(svtype, {"x": [], "y": [], "len": [], "class": [], "name": [], "fam": []})
                by_type_points[svtype]["x"].append(x0)
                by_type_points[svtype]["y"].append(ylane)
                by_type_points[svtype]["len"].append(size_bp)
                by_type_points[svtype]["class"].append(sv[index_sv])
                by_type_points[svtype]["name"].append(sv[index_name])
                by_type_points[svtype]["fam"].append(sv[index_fam])
            else:
                by_type_segments.setdefault(svtype, {"x": [], "y": [], "len": [], "class": [], "name": [], "fam": []})
                by_type_segments[svtype]["x"].extend([x0, x1, None])
                by_type_segments[svtype]["y"].extend([ylane, ylane, None])
                by_type_segments[svtype]["len"].append(size_bp)
                by_type_segments[svtype]["class"].append(sv[index_sv])
                by_type_segments[svtype]["name"].append(sv[index_name])
                by_type_segments[svtype]["fam"].append(sv[index_fam])

        # add traces
        for svtype, data in by_type_points.items():
            fig.add_trace(
                go.Scatter(
                    x=data["x"], y=data["y"], mode='markers',
                    marker=dict(symbol='square', size=8, color=normalize_color(color_map_sv[svtype])),
                    name=svtype, legendgroup=svtype,
                    hovertemplate=(
                        "<b>Class:</b> %{customdata[0]}<br>"
                        "<b>Name:</b> %{customdata[1]}<br>"
                        "<b>Family:</b> %{customdata[2]}<br>"
                        "<b>Length (bp):</b> %{customdata[3]}<extra></extra>"
                    ),
                    customdata=np.stack([
                        data["class"],  # SV type/class
                        data["name"],   # SV name
                        data["fam"],    # SV family
                        data["len"],    # SV length
                    ], axis=-1) if "class" in data else None,
                ),
                row=3, col=1
            )
        for svtype, data in by_type_segments.items():
            x = data["x"]
            y = data["y"]

            # Build per-point customdata aligned to x/y
            cd = []
            for cls, name, fam, ln in zip(data["class"], data["name"], data["fam"], data["len"]):
                row = [cls, name, fam, ln]
                cd.extend([row, row, [None, None, None, None]])  # start, end, separator

            fig.add_trace(
                go.Scatter(
                    x=x, y=y, mode='lines',
                    line=dict(color=normalize_color(color_map_sv[svtype]), width=6),
                    name=svtype, legendgroup=svtype,
                    showlegend=(svtype not in by_type_points),
                    customdata=np.array(cd, dtype=object),
                    hovertemplate=(
                        "<b>Class:</b> %{customdata[0]}<br>"
                        "<b>Name:</b> %{customdata[1]}<br>"
                        "<b>Family:</b> %{customdata[2]}<br>"
                        "<b>Length (bp):</b> %{customdata[3]}<extra></extra>"
                    ),
                ),
                row=3, col=1
            )

        fig.update_yaxes(range=[-max_y_sv*0.10, max_y_sv+max_y_sv*0.20], row=3, col=1)
    else:
        fig.add_annotation(
            x=(start_pos + (end_pos-start_pos)/2)/1_000_000,
            y=0.5,
            text="<b><span style='font-size:20px'>No Structural variants to plot</span></b>",
            showarrow=False,
            font=dict(size=20, color="black", family="Arial",),
            row=3, col=1
        )

    fig.update_xaxes(title_text="Genomic Position (Mb)", row=3, col=1, title_font=dict(size=18, family="Arial", color="black", weight="bold"), tickfont=dict(size=16, family="Arial", color="black"), showticklabels=True, ticks="outside", showline=True, automargin=True)
    fig.update_yaxes(title_text="Structural variants", row=3, col=1, title_font=dict(size=18, family="Arial", color="black", weight="bold"), tickfont=dict(size=14, family="Arial", color="white"))

    # x-range across panels
    xr = [start_pos/1_000_000, end_pos/1_000_000]
    fig.update_xaxes(range=xr, row=2, col=1)
    fig.update_xaxes(range=xr, row=3, col=1)

    # truncate legend labels (purely cosmetic)
    truncated = {}
    for k in list(color_map.keys()):
        if len(k) > 20:
            truncated[k[:20] + "."] = color_map.pop(k)
    color_map.update(truncated)
    
    # increase legend size
    fig.update_layout(legend=dict(font=dict(size=18), itemsizing='constant', yanchor="top", y=1.12, xanchor="left", x=0), margin=dict(t=10))
   
    return fig

# function to normalize colors
def normalize_color(c):
    """Return a Plotly-safe color.
       Accepts '#RRGGBB', '#RRGGBBAA', '#RGB', bare hex, named CSS, or lists/Series."""
    if isinstance(c, (list, tuple, np.ndarray, pd.Series)):
        return [normalize_color(x) for x in list(c)]
    if not isinstance(c, str):
        return c
    s = c.strip()
    # allow bare hex
    if re.fullmatch(r'[0-9a-fA-F]{8}', s): s = '#' + s
    if re.fullmatch(r'[0-9a-fA-F]{6}', s): s = '#' + s
    if re.fullmatch(r'[0-9a-fA-F]{3}', s): s = '#' + s
    # #RRGGBBAA -> rgba(...)
    m = re.fullmatch(r'#([0-9a-fA-F]{8})', s)
    if m:
        h = m.group(1)
        r, g, b, a = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16), int(h[6:8], 16) / 255.0
        return f'rgba({r},{g},{b},{a:.3f})'
    # #RGB -> #RRGGBB
    m = re.fullmatch(r'#([0-9a-fA-F]{3})', s)
    if m:
        h = m.group(1)
        s = '#' + ''.join(ch*2 for ch in h)
    return s  # '#RRGGBB' or named CSS

# function to generate the GTEx heatmap
def gtex_heatmap(gtex_info):
    cols_heatmap = gtex_info.columns[2:]
    zdf = gtex_info[cols_heatmap].apply(pd.to_numeric, errors="coerce")
    x_labels = list(cols_heatmap)
    y_labels = gtex_info['Description'].astype(str).tolist()
    n_rows, n_cols = zdf.shape

    # y as numeric indices so rows are unique
    y_idx = list(range(n_rows))

    # CUSTOMDATA must match z's shape (n_rows x n_cols)
    cd = np.tile(np.array(y_labels, dtype=object)[:, None], (1, n_cols))

    fig = go.Figure(
        data=go.Heatmap(
            z=zdf.values,
            x=x_labels,
            y=y_idx,
            colorscale='Viridis',
            customdata=cd,  # 2D
            hovertemplate="<b>%{customdata}</b><br>Tissue: %{x}<br>Value: %{z:.3g}<extra></extra>",
        )
    )

    fig.update_yaxes(
        tickmode='array',
        tickvals=y_idx,
        ticktext=y_labels,
        autorange='reversed',
        automargin=True,
        tickfont=dict(size=12)
    )
    fig.update_xaxes(tickangle=45, tickfont=dict(size=12), automargin=True)

    per_row = 22
    fig.update_layout(height=max(400, min(140 + per_row * n_rows, 1200)),
                      margin=dict(l=110, r=40, t=40, b=110))
    return fig

# function to check if an email is valid
def is_valid_email(email):
    pattern = r'^[\w\.-]+@[\w\.-]+\.\w+$'
    return re.match(pattern, email) is not None

if __name__ == "__main__":
    app.run()