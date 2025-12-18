import sqlite3
import re
from pathlib import Path
import pandas as pd
from flask import Flask, request, jsonify, render_template, Blueprint, redirect, url_for, session, current_app, send_from_directory, abort
from liftover import get_lifter
import io
import os
import subprocess
from datetime import timedelta
import math
from liftover import get_lifter
import os
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
from ieugwaspy import api, query
from datetime import timedelta, datetime
import sqlite3
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import plotly.io as pio
import plotly.colors as pc
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.pairwise import cosine_similarity
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np
import pickle
from typing import Optional

# Blueprint setup (if needed)
prs_bp = Blueprint('prs', __name__, template_folder='templates')

# ---------------------------------------------------------
# Configuration
# ---------------------------------------------------------
#DATA_PATH = Path("/Users/nicco/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/snpXplorer/Data")
DATA_PATH = Path('/Data')
DB_FILE = DATA_PATH / "databases/Genes/variant_info.db"
trait_names = pd.read_csv(DATA_PATH / "databases/haplotypes/trait_names_and_ids.txt", sep="\t")
TRAIT_LIST = trait_names["name"].astype(str).tolist()
PRS_DIR = DATA_PATH / "databases/PRS"

# ---------------------------------------------------------
# Functions
# ---------------------------------------------------------
# Function to get info about trait and umap
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

# Function to plot UMAP + heatmap
def plot_haplotype_traits(clusters_of_interest, all_indices, trait_list, cluster_index_map, umap, cluster_representatives, browse, sim_mat, names_sub, cluster_labels_for_index, cluster_trait, idx_trait, all_cluster_representatives):
    # ------------------------------------------------
    # Shared color map for clusters (UMAP + annotation)
    # ------------------------------------------------
    palette = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00"]
    color_map = {cid: palette[j % len(palette)] for j, cid in enumerate(clusters_of_interest)}

    # -------------------------------
    # Create FIGURE with 2 SUBPLOTS
    # -------------------------------
    fig = make_subplots(rows=1, cols=2, column_widths=[0.6, 0.40], horizontal_spacing=0.04, subplot_titles=("UMAP Trait Projection", "Trait Similarity Heatmap"))
    # increase subplot subtitle font size
    fig.update_annotations(font=dict(size=20, family="Arial", color="darkgrey", weight="bold"))

    # -------------------------------
    # LEFT subplot: UMAP
    # -------------------------------
    # Use the SAME indices used for the heatmap for colored clusters
    all_cluster_idxs = np.unique(all_indices)
    all_points = np.arange(len(trait_list))
    background_idxs = np.setdiff1d(all_points, all_cluster_idxs)
    # Colored clusters (but capped to the same subset as the heatmap)
    for j, cid in enumerate(clusters_of_interest):
        # Intersect to keep at most the terms we actually used (max_terms_per_cluster)
        idxs_full = np.array(cluster_index_map[cid])
        idxs = np.intersect1d(idxs_full, all_cluster_idxs)
        if idxs.size == 0:
            continue
        fig.add_trace(go.Scatter(x=umap[idxs, 0], y=umap[idxs, 1], mode="markers", marker=dict(size=9, color=color_map[cid]), name=all_cluster_representatives.loc[all_cluster_representatives['Cluster Label'] == cid, 'Final Representative'].values[0], text=[trait_list[i] for i in idxs], hoverinfo="text",), row=1, col=1)
    # Highlight queried trait with a star
    fig.add_trace(go.Scatter(x=[umap[idx_trait, 0]], y=[umap[idx_trait, 1]], mode="markers+text", text=["★"], textposition="top center", marker=dict(size=15, color="black", symbol="star", line=dict(width=2, color="yellow")), name=browse, hovertext=browse, hoverinfo="text",), row=1, col=1)
    # Update axes
    fig.update_xaxes(title_text="UMAP-1", row=1, col=1, showticklabels=False, title_font=dict(size=20, family="Arial", color="black", weight="bold"))
    fig.update_yaxes(title_text="UMAP-2", row=1, col=1, showticklabels=False, title_font=dict(size=20, family="Arial", color="black", weight="bold"))
    # Legend below UMAP
    fig.update_layout(legend=dict(orientation="h", x=0, y=-0.05, xanchor="left", font=dict(size=16, family="Arial"), itemsizing="constant"))

    # -------------------------------
    # RIGHT subplot: Heatmap
    # -------------------------------
    n_traits = sim_mat.shape[0]
    # Heatmap with flipped colors (blue = -1, red = +1) and no axis labels
    hover_text = [[f"{names_sub[i]} ↔ {names_sub[j]}" for j in range(n_traits)] for i in range(n_traits)]
    fig.add_trace(go.Heatmap(z=sim_mat, x=list(range(n_traits)), y=list(range(n_traits)), text=hover_text, hovertemplate=("Trait pair: %{text}<br>" "Sim=%{z:.3f}<extra></extra>"), zmin=0, zmax=1, colorscale=[[0.0, "blue"], [0.5, "white"], [1.0, "red"]], colorbar=dict(title="Cosine similarity", orientation="v", x=1.02, y=0.83, len=0.3, thickness=15), showscale=True), row=1, col=2)
    # Cluster annotation strip on the right
    cluster_colors_for_traits = [color_map.get(cid, "#d3d3d3") for cid in cluster_labels_for_index]
    fig.add_trace(go.Scatter(x=[n_traits + 0.5] * n_traits, y=list(range(n_traits)), mode="markers", marker=dict(color=cluster_colors_for_traits, symbol="square", size=8, line=dict(width=0)), showlegend=False, hoverinfo="skip"), row=1, col=2)
    # Hide axis labels / ticks for heatmap panel
    fig.update_xaxes(showticklabels=False, title_text="", range=[-0.5, n_traits + 1.0], row=1, col=2)
    fig.update_yaxes(showticklabels=False, title_text="", row=1, col=2)
    
    # Overall layout
    fig.update_layout(template="plotly_white", height=650, title={"text": f"<b>Semantic Neighborhood of '{browse}'</b>", "font": {"size": 32}})
    # return as HTML
    plot_url = pio.to_html(fig, include_plotlyjs='cdn', full_html=False)
    return plot_url

# Function to retrieve PRS info
def getPRSinfo(meta_sb, cluster_rep, browse_resolved):
    # Define path to PRS data
    prs_path = DATA_PATH / "databases/PRS/"
    # Get files of interest
    prs_interest = meta_sb['ID'].values.tolist()
    # Container for results
    prs_info = []
    # Define pvalue thresholds 
    pvals = [5e-5, 5e-6, 5e-7, 5e-8]
    # Iterate over PRS files and extract relevant info
    for prs in prs_interest:
        try:
            tmp_path = prs_path / f"{prs}_hits_p5e-5_clumped_reformatted.txt.gz"
            tmp_dic = {'id' : prs, 'file_name' : tmp_path.name}
            tmp_prs = pd.read_csv(tmp_path, sep="\t")
            tmp_prs['P'] = tmp_prs['P'].astype(float)
            # Save all snps
            tmp_dic['all'] = tmp_prs.shape[0]
            # Find snps at pvalue thresholds
            for pv in pvals:
                tmp_prs_pv = tmp_prs[tmp_prs['P'] < pv].copy()
                tmp_dic[str(pv)] = tmp_prs_pv.shape[0]
            prs_info.append(tmp_dic)
        except:
            continue
    # Convert to dataframe
    prs_df = pd.DataFrame(prs_info)
    # Flag the representative
    prs_df['representative'] = prs_df['id'].apply(lambda x: 'Yes' if x == cluster_rep else 'No')
    # Flag non-curatedness
    prs_df['curated'] = 'No'
    # Merge with meta_sb
    prs_df = prs_df.merge(meta_sb[['Trait', 'ID', 'Author', 'Year', 'PMID', 'Sample Size']], left_on='id', right_on='ID')
    # Remove ID column
    prs_df = prs_df.drop(columns=['ID'])
    # Take care of curated PRS
    curated_prs = curatedPRS()
    # Check if curated PRS should be considered
    if 'alzheimer' in browse_resolved.lower():
        # Merge with curated
        prs_df = pd.concat([prs_df, curated_prs[curated_prs['Trait'] == "Alzheimer's Disease - Bellenguez et al. - 2022"]], ignore_index=True)
    elif 'parkinson' in browse_resolved.lower():
        # Merge with curated
        prs_df = pd.concat([prs_df, curated_prs[curated_prs['Trait'] == "Parkinson's Disease - Nalls et al. - 2025"]], ignore_index=True)
    # Sort by curatedness first, then representativeness
    prs_df = prs_df.sort_values(by=['curated', 'representative'], ascending=[False, False])
    # Rename columns
    prs_df = prs_df.rename(columns = {'id': 'ID', 'all': 'All SNPs', '5e-05': 'P<5e-05', '5e-06': 'P<5e-06', '5e-07': 'P<5e-07', '5e-08': 'P<5e-08', 'representative': 'Suggested', 'curated': 'Curated'})
    return prs_df

# Function to manage manually curated PRSs
def curatedPRS():
    curated_paths = {"Alzheimer's Disease - Bellenguez et al. - 2022": "Bellenguez_2022_forJordan.txt;Bellenguez et al.;2022;35379992",
                     "Parkinson's Disease - Nalls et al. - 2025": "PDsnps_info_largestNalls_jordan.txt;Nalls et al.;2025;NA"}
    # Container for results
    curated = []
    # Iterate through files and read them
    counter = 1
    for tr in curated_paths.keys():
        # Get info
        tr_path, tr_author, tr_year, tr_pmid = curated_paths[tr].split(';')
        # Get file
        tmp_path = DATA_PATH / f"databases/PRS/{tr_path}"
        # Read PRS
        tmp_prs = pd.read_csv(tmp_path, sep="\t")
        curated.append({'id': f"curated-{counter}", "file_name": tmp_path.name, "Author": tr_author, "Year": tr_year, "PMID": tr_pmid, "Trait": tr, 'all': tmp_prs.shape[0], '5e-05': 0, '5e-06': 0, '5e-07': 0, '5e-08': tmp_prs.shape[0], 'representative': 'Yes', 'curated': 'Yes'})
        counter += 1
    return pd.DataFrame(curated)

# Function to find trait index
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

# Function to build cluster to index map
def build_cluster_index_map(labels):
    """
    labels: np.array shape (n_traits,)
    return: dict cluster_id -> list of indices
    """
    cluster_index_map = {}
    for i, lab in enumerate(labels):
        cluster_index_map.setdefault(int(lab), []).append(i)
    return cluster_index_map

# Function to monitor trait searches
def add_search_to_file(DATA_PATH, q):
    log_file = f"{DATA_PATH}/monitor/prs_logs.txt"
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_entry = f"{timestamp}\t{q}\n"
    with open(log_file, "a") as f:
        f.write(log_entry)

# ---------------------------------------------------------
# Routes
# ---------------------------------------------------------
# Route for PRS page
@prs_bp.route("/")
def index():
    return render_template("prs.html", trait_list=TRAIT_LIST)

# Route to handle trait search
@prs_bp.route("/traits", methods=["GET", "POST"])
def trait_search():
    meta = current_app.config["GWAS_META"]
    # Remember last selection ---
    if request.method == "POST":
        query = (request.form.get("trait") or "").strip()
        if query:
            session["prs_last_trait"] = query
        else:
            session.pop("prs_last_trait", None)
    else:
        query = (session.get("prs_last_trait") or "").strip()

    # Initialize variables
    plot_html = None
    meta_rows = []
    prs_rows = []

    # If we have a query (POST or remembered GET), compute page content
    if query:
        q = query.lower()
        # Monitoring: record search
        add_search_to_file(DATA_PATH, q)
        # Get all relevant info
        clusters_of_interest, all_indices, trait_list, cluster_index_map, umap, cluster_representatives, browse_resolved, sim_mat, names_sub, cluster_labels_for_index, cluster_traits, traits_interest, idx_trait, all_cluster_representatives = guide_haplotypes_traits(DATA_PATH, q, 100000, "hg38", trait_names)
        # Generate plots
        plot_html = plot_haplotype_traits(clusters_of_interest, all_indices, trait_list, cluster_index_map, umap, cluster_representatives, browse_resolved, sim_mat, names_sub, cluster_labels_for_index, cluster_traits, idx_trait, all_cluster_representatives)
        # Subset metadata
        meta_sb = meta.loc[:, meta.loc["trait"].isin(traits_interest["name"])].copy()
        meta_subset_long = meta_sb.T.reset_index(drop=False)
        meta_subset_long = meta_subset_long.rename(columns={"trait": "Trait", "author": "Author", "year": "Year", "pmid": "PMID", "sample_size": "Sample Size", "population": "Population", "doi": "DOI", "id": "ID"})
        meta_sb = meta_subset_long[["Trait", "ID", "Author", "Year", "PMID", "Sample Size", "Population"]].copy()
        # Get cluster representative
        cluster_rep = cluster_representatives["Final Representative file"].replace("_hits_p5e-5.txt.gz", "", regex=True).values[0]
        # Flag cluster representative
        meta_sb["Cluster Representative"] = meta_sb["ID"].apply(lambda x: "Yes" if x == cluster_rep else "No")
        # Get PRS info
        prs_info = getPRSinfo(meta_sb, cluster_rep, browse_resolved)
        prs_info["download_url"] = prs_info["file_name"].apply(lambda fn: url_for("prs.download_prs_file", filename=fn) if fn else "")
        prs_info = prs_info.drop(columns=["file_name"], errors="ignore")
        # Prepare tables
        meta_rows = meta_sb.to_dict(orient="records")
        prs_rows = prs_info.to_dict(orient="records")

    return render_template("prs.html", trait_list=TRAIT_LIST, query=query, plot_html=plot_html, meta_rows=meta_rows, prs_rows=prs_rows)

# Serve PRS files
@prs_bp.route("/download/<path:filename>")
def download_prs_file(filename):
    if not filename:
        abort(404)
    filename = Path(filename).name
    return send_from_directory(PRS_DIR, filename, as_attachment=True)