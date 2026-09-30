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
from functools import lru_cache

# Blueprint setup (if needed)
prs_bp = Blueprint('prs', __name__, template_folder='templates')

# ---------------------------------------------------------
# Configuration
# ---------------------------------------------------------
#DATA_PATH = Path("/Users/nicco/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/snpXplorer/Data")
DATA_PATH = Path('/Data')
DB_FILE = DATA_PATH / "databases/Genes/variant_info.db"
TRAIT_FILE = DATA_PATH / "databases/haplotypes/combined_traits_and_ids.txt"
COMBINED_EMBEDDINGS_FILE = DATA_PATH / "databases/haplotypes/20260918_combined_umap_embeddings.pkl"
CLUSTER_SUMMARY_FILES = {
    "European": DATA_PATH / "databases/haplotypes/trait_clusters_summary_AI_thresholds.csv",
    "African": DATA_PATH / "databases/haplotypes/20260916_trait_clusters_summary_AI_thresholds_African.csv",
}

trait_names = pd.read_csv(TRAIT_FILE, sep="\t")
trait_names = trait_names[
    ~(
        trait_names["id"].astype(str).eq("trait_association_counts.txt.gz")
        & trait_names["name"].astype(str).eq("File Name")
    )
].copy()
trait_names["ancestry"] = trait_names["ancestry"].astype(str).str.strip().str.title()
trait_names["dataset_index"] = trait_names.groupby("ancestry", sort=False).cumcount()
trait_names["gwas_id"] = trait_names["id"].astype(str).str.replace(
    "_hits_p5e-5.txt.gz", "", regex=False
)
trait_names["selection_id"] = trait_names["ancestry"] + "|||" + trait_names["gwas_id"]

ANCESTRY_CODES = {"European": "EUR", "African": "AFR"}
SEMANTIC_CLUSTER_PALETTE = [
    "#000000",
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#0072B2",
    "#D55E00",
]
searchable_traits = trait_names[trait_names["name"].notna()].copy()
TRAIT_OPTIONS = [
    {
        "id": row.selection_id,
        "text": f"[{ANCESTRY_CODES.get(row.ancestry, row.ancestry.upper())}] {row.name}",
    }
    for row in searchable_traits.itertuples(index=False)
]
PRS_DIR = DATA_PATH / "databases/PRS"

# ---------------------------------------------------------
# Functions
# ---------------------------------------------------------
@lru_cache(maxsize=1)
def load_semantic_artifacts():
    """Load and validate the joint European/African semantic artifacts."""
    class NumpyCompatUnpickler(pickle.Unpickler):
        def find_class(self, module, name):
            if module == "numpy._core":
                module = "numpy.core"
            elif module == "numpy._core.multiarray":
                module = "numpy.core.multiarray"
            return super().find_class(module, name)
    with open(COMBINED_EMBEDDINGS_FILE, "rb") as handle:
        artifact = NumpyCompatUnpickler(handle).load()

    dataset_order = list(artifact["datasets"])
    summaries = {}
    for ancestry in dataset_order:
        dataset = artifact["datasets"][ancestry]
        ancestry_rows = trait_names[trait_names["ancestry"] == ancestry]
        expected = len(dataset["labels"])
        if len(ancestry_rows) != expected:
            raise ValueError(
                f"{ancestry} trait rows ({len(ancestry_rows)}) do not match "
                f"the semantic dataset ({expected})."
            )
        summary_all = pd.read_csv(CLUSTER_SUMMARY_FILES[ancestry])
        label_set = {int(value) for value in dataset["cluster_map"]}
        matching_summary = None
        for threshold in sorted(summary_all["Threshold"].dropna().unique()):
            candidate = summary_all[summary_all["Threshold"] == threshold].copy()
            if set(candidate["Cluster Label"].astype(int)) == label_set:
                matching_summary = candidate.reset_index(drop=True)
                break
        if matching_summary is None:
            raise ValueError(f"No cluster summary matches the {ancestry} semantic labels.")
        summaries[ancestry] = matching_summary

    # The current joint UMAP was generated European-first even though the
    # serialized datasets dictionary is African-first. New artifacts can make
    # this explicit by including a top-level ``coordinate_order`` list.
    coordinate_order = artifact.get("coordinate_order")
    if coordinate_order is None:
        known_order = ["European", "African"]
        coordinate_order = [name for name in known_order if name in artifact["datasets"]]
        coordinate_order.extend(name for name in dataset_order if name not in coordinate_order)
    if set(coordinate_order) != set(dataset_order):
        raise ValueError("Joint UMAP coordinate order does not match the ancestry datasets.")

    offset = 0
    joint_offsets = {}
    for ancestry in coordinate_order:
        joint_offsets[ancestry] = offset
        offset += len(artifact["datasets"][ancestry]["labels"])
    if offset != len(artifact["coordinates"]):
        raise ValueError("Joint UMAP coordinates do not match the ancestry datasets.")
    artifact["coordinate_order"] = coordinate_order
    artifact["joint_offsets"] = joint_offsets
    artifact["summaries"] = summaries
    return artifact


def resolve_trait_selection(selection):
    """Resolve a Select2 value to exactly one ancestry-specific study row."""
    selection = str(selection or "").strip()
    match = trait_names[trait_names["selection_id"] == selection]
    if not match.empty:
        return match.iloc[0]

    # Backward-compatible fallback for a remembered pre-ancestry trait value.
    named = trait_names[trait_names["name"].notna()].copy()
    exact = named[named["name"].astype(str).str.lower() == selection.lower()]
    if exact.empty:
        exact = named[named["name"].astype(str).str.lower().str.contains(selection.lower(), regex=False)]
    if exact.empty:
        raise ValueError(f"Trait '{selection}' was not found.")
    european = exact[exact["ancestry"] == "European"]
    return european.iloc[0] if not european.empty else exact.iloc[0]


def _cluster_indices(labels):
    result = {}
    for index, label in enumerate(labels):
        result.setdefault(int(label), []).append(index)
    return result


def build_semantic_context(artifact, selected_row, scope):
    """Build a query-centred ancestry-specific or combined semantic neighborhood."""
    query_ancestry = selected_row["ancestry"]
    query_index = int(selected_row["dataset_index"])
    query_dataset = artifact["datasets"][query_ancestry]
    query_embedding = query_dataset["embeddings"][query_index]
    query_cluster = int(query_dataset["labels"][query_index])

    ancestries = [scope] if scope != "Combined" else list(artifact["datasets"])
    cluster_indices = {}
    centroids = {}
    for ancestry in ancestries:
        dataset = artifact["datasets"][ancestry]
        by_cluster = _cluster_indices(dataset["labels"])
        for cluster_id, indices in by_cluster.items():
            key = (ancestry, cluster_id)
            cluster_indices[key] = indices
            centroid = dataset["embeddings"][indices].mean(axis=0)
            centroid /= np.linalg.norm(centroid) + 1e-9
            centroids[key] = centroid

    query_key = (query_ancestry, query_cluster)
    selected_clusters = []
    if query_key in centroids:
        selected_clusters.append(query_key)
    remaining = [key for key in centroids if key != query_key]
    if remaining:
        similarities = cosine_similarity(
            query_embedding.reshape(1, -1),
            np.vstack([centroids[key] for key in remaining]),
        )[0]
        ranked_clusters = [remaining[index] for index in np.argsort(-similarities)]
        selected_clusters.extend(ranked_clusters[:6 - len(selected_clusters)])

        # A combined neighborhood must visibly contain each available ancestry.
        # If the globally closest clusters all come from one ancestry, replace
        # its least-similar extra cluster with the closest missing ancestry.
        if scope == "Combined":
            for ancestry in ancestries:
                if any(key[0] == ancestry for key in selected_clusters):
                    continue
                replacement = next(
                    (key for key in ranked_clusters if key[0] == ancestry),
                    None,
                )
                if replacement is None:
                    continue
                ancestry_counts = {
                    name: sum(key[0] == name for key in selected_clusters)
                    for name in ancestries
                }
                replace_index = next(
                    (
                        index
                        for index in range(len(selected_clusters) - 1, -1, -1)
                        if selected_clusters[index] != query_key
                        and ancestry_counts[selected_clusters[index][0]] > 1
                    ),
                    None,
                )
                if replace_index is None:
                    selected_clusters.append(replacement)
                else:
                    selected_clusters[replace_index] = replacement

    points = []
    heatmap_rows = []
    for cluster_key in selected_clusters:
        ancestry, cluster_id = cluster_key
        dataset = artifact["datasets"][ancestry]
        indices = np.asarray(cluster_indices[cluster_key], dtype=int)
        similarities = cosine_similarity(
            dataset["embeddings"][indices], query_embedding.reshape(1, -1)
        )[:, 0]
        indices = indices[np.argsort(-similarities)[:100]]
        ancestry_rows = trait_names[trait_names["ancestry"] == ancestry].set_index("dataset_index")
        names = [
            ancestry_rows.at[int(index), "name"]
            if int(index) in ancestry_rows.index and pd.notna(ancestry_rows.at[int(index), "name"])
            else f"Unnamed {ANCESTRY_CODES.get(ancestry, ancestry)} study"
            for index in indices
        ]
        display_names = [
            f"[{ANCESTRY_CODES.get(ancestry, ancestry)}] {name}"
            for name in names
        ]
        if scope == "Combined":
            coordinates = artifact["coordinates"][artifact["joint_offsets"][ancestry] + indices]
        else:
            coordinates = dataset["umap"][indices]
        points.append({
            "key": cluster_key,
            "indices": indices,
            "coordinates": coordinates,
            "names": display_names,
        })
        for index, name in zip(indices, display_names):
            heatmap_rows.append({
                "key": cluster_key,
                "embedding": dataset["embeddings"][int(index)],
                "name": name,
            })

    # Different GWAS can carry the same displayed trait name. Keep only the
    # first (highest-ranked) occurrence so heatmap axes do not repeat labels.
    unique_heatmap_rows = []
    seen_heatmap_names = set()
    for row in heatmap_rows:
        normalized_name = row["name"].strip().casefold()
        if normalized_name in seen_heatmap_names:
            continue
        seen_heatmap_names.add(normalized_name)
        unique_heatmap_rows.append(row)
    heatmap_rows = unique_heatmap_rows

    embeddings = np.vstack([row["embedding"] for row in heatmap_rows])
    names = [row["name"] for row in heatmap_rows]
    similarity_matrix = cosine_similarity(embeddings)
    if scope == "Combined":
        query_coordinate = artifact["coordinates"][
            artifact["joint_offsets"][query_ancestry] + query_index
        ]
    else:
        query_coordinate = query_dataset["umap"][query_index]
    return {
        "scope": scope,
        "points": points,
        "heatmap_names": names,
        "heatmap_keys": [row["key"] for row in heatmap_rows],
        "similarity_matrix": similarity_matrix,
        "query_coordinate": query_coordinate,
        "query_key": query_key,
    }


def _cluster_name(artifact, cluster_key):
    ancestry, cluster_id = cluster_key
    summary = artifact["summaries"][ancestry]
    row = summary[summary["Cluster Label"].astype(int) == int(cluster_id)]
    if not row.empty:
        representative = row.iloc[0].get("Final Representative")
        if pd.notna(representative) and str(representative).strip():
            return f"[{ANCESTRY_CODES.get(ancestry, ancestry)}] {representative}"
    fallback = artifact["datasets"][ancestry]["cluster_map"].get(cluster_id, {}).get(
        "representative", f"Cluster {cluster_id}"
    )
    return f"[{ANCESTRY_CODES.get(ancestry, ancestry)}] {fallback}"


def plot_haplotype_traits(artifact, selected_row):
    """Plot ancestry-specific and combined query-centred semantic neighborhoods."""
    query_ancestry = selected_row["ancestry"]
    contexts = [
        build_semantic_context(artifact, selected_row, query_ancestry),
        build_semantic_context(artifact, selected_row, "Combined"),
    ]
    figure = make_subplots(
        rows=1,
        cols=2,
        column_widths=[0.6, 0.40],
        horizontal_spacing=0.04,
        subplot_titles=("UMAP semantic trait projection", "Semantic trait similarity"),
    )
    figure.update_annotations(font=dict(size=20, family="Arial", color="darkgrey", weight="bold"))

    umap_traces_by_context = []
    context_color_maps = []
    for context_index, context in enumerate(contexts):
        visible = context_index == 0
        trace_indices = []
        color_map = {
            point_group["key"]: SEMANTIC_CLUSTER_PALETTE[index % len(SEMANTIC_CLUSTER_PALETTE)]
            for index, point_group in enumerate(context["points"])
        }
        context_color_maps.append(color_map)
        for point_index, point_group in enumerate(context["points"]):
            cluster_key = point_group["key"]
            color = color_map[cluster_key]
            figure.add_trace(
                go.Scatter(
                    x=point_group["coordinates"][:, 0],
                    y=point_group["coordinates"][:, 1],
                    mode="markers",
                    marker=dict(size=9, color=color),
                    name=_cluster_name(artifact, cluster_key),
                    text=point_group["names"],
                    hoverinfo="text",
                    legendrank=10 + point_index,
                    visible=visible,
                ),
                row=1,
                col=1,
            )
            trace_indices.append(len(figure.data) - 1)

        query_label = f"[{ANCESTRY_CODES.get(selected_row['ancestry'], selected_row['ancestry'])}] {selected_row['name']}"
        figure.add_trace(
            go.Scatter(
                x=[context["query_coordinate"][0]],
                y=[context["query_coordinate"][1]],
                mode="markers",
                marker=dict(size=15, color="black", symbol="star", line=dict(width=2, color="yellow")),
                name="Queried trait",
                hovertext=query_label,
                hoverinfo="text",
                legendrank=0,
                visible=visible,
            ),
            row=1,
            col=1,
        )
        trace_indices.append(len(figure.data) - 1)
        umap_traces_by_context.append(trace_indices)

    # The heatmap always reflects the ancestry of the selected study, even when
    # the UMAP is switched to the combined ancestry view.
    heatmap_context = contexts[0]
    names = heatmap_context["heatmap_names"]
    n_traits = len(names)
    hover_text = [
        [f"{names[row]} ↔ {names[column]}" for column in range(n_traits)]
        for row in range(n_traits)
    ]
    figure.add_trace(
        go.Heatmap(
            z=heatmap_context["similarity_matrix"],
            x=list(range(n_traits)),
            y=list(range(n_traits)),
            text=hover_text,
            hovertemplate="Trait pair: %{text}<br>Sim=%{z:.3f}<extra></extra>",
            zmin=0,
            zmax=1,
            colorscale=[[0.0, "blue"], [0.5, "white"], [1.0, "red"]],
            colorbar=dict(
                title="Cosine similarity",
                orientation="v",
                x=1.02,
                y=0.83,
                len=0.3,
                thickness=15,
            ),
            showscale=True,
        ),
        row=1,
        col=2,
    )
    heatmap_trace_index = len(figure.data) - 1

    annotation_colors = [
        context_color_maps[0].get(key, "#d3d3d3")
        for key in heatmap_context["heatmap_keys"]
    ]
    unique_annotation_colors = list(dict.fromkeys(annotation_colors))
    color_to_index = {
        color: index for index, color in enumerate(unique_annotation_colors)
    }
    strip_width = 4
    strip_column = np.array(
        [color_to_index[color] for color in annotation_colors],
        dtype=float,
    ).reshape(-1, 1)
    strip_values = np.repeat(strip_column, strip_width, axis=1)
    if len(unique_annotation_colors) == 1:
        strip_colorscale = [
            [0.0, unique_annotation_colors[0]],
            [1.0, unique_annotation_colors[0]],
        ]
    else:
        denominator = float(len(unique_annotation_colors) - 1)
        strip_colorscale = [
            [index / denominator, color]
            for index, color in enumerate(unique_annotation_colors)
        ]
    figure.add_trace(
        go.Heatmap(
            z=strip_values,
            x=list(range(n_traits, n_traits + strip_width)),
            y=list(range(n_traits)),
            zmin=0,
            zmax=max(1, len(unique_annotation_colors) - 1),
            colorscale=strip_colorscale,
            showscale=False,
            hoverinfo="skip",
        ),
        row=1,
        col=2,
    )
    annotation_trace_index = len(figure.data) - 1

    buttons = []
    other_ancestries = [name for name in artifact["datasets"] if name != query_ancestry]
    combined_label = " + ".join([query_ancestry] + other_ancestries)
    labels = [f"{query_ancestry} only", combined_label]
    for context_index, label in enumerate(labels):
        visibility = [False] * len(figure.data)
        for trace_index in umap_traces_by_context[context_index]:
            visibility[trace_index] = True
        visibility[heatmap_trace_index] = True
        visibility[annotation_trace_index] = True
        buttons.append(dict(
            label=label,
            method="update",
            args=[
                {"visible": visibility},
                {
                    "title.text": f"<b>Semantic Neighborhood of '{selected_row['name']}'</b>",
                    "xaxis.autorange": True,
                    "yaxis.autorange": True,
                },
            ],
        ))

    figure.update_xaxes(
        title_text="UMAP-1",
        showticklabels=False,
        title_font=dict(size=20, family="Arial", color="black", weight="bold"),
        row=1,
        col=1,
    )
    figure.update_yaxes(
        title_text="UMAP-2",
        showticklabels=False,
        title_font=dict(size=20, family="Arial", color="black", weight="bold"),
        row=1,
        col=1,
    )
    figure.update_xaxes(
        showticklabels=False,
        title_text="",
        range=[-0.5, n_traits + strip_width - 0.5],
        row=1,
        col=2,
    )
    figure.update_yaxes(showticklabels=False, title_text="", row=1, col=2)
    figure.update_layout(
        template="plotly_white",
        height=650,
        title={
            "text": f"<b>Semantic Neighborhood of '{selected_row['name']}'</b>",
            "font": {"size": 32},
            "x": 0.5,
            "xanchor": "center",
            "y": 0.99,
            "yanchor": "top",
        },
        margin=dict(t=145, b=110),
        updatemenus=[dict(
            buttons=buttons,
            direction="down",
            showactive=True,
            active=0,
            x=0.0,
            xanchor="left",
            y=1.10,
            yanchor="top",
        )],
        legend=dict(
            orientation="h",
            x=0,
            y=-0.12,
            xanchor="left",
            font=dict(size=16, family="Arial"),
            itemsizing="constant",
        ),
    )
    return pio.to_html(figure, include_plotlyjs="cdn", full_html=False)

# Function to retrieve PRS info
def getPRSinfo(meta_sb, study_ids, cluster_rep, browse_resolved, ancestry):
    # Define path to PRS data
    prs_path = PRS_DIR / "African" if ancestry == "African" else PRS_DIR
    # Get files of interest
    prs_interest = list(dict.fromkeys(study_ids))
    # Container for results
    prs_info = []
    # Define pvalue thresholds 
    pvals = [5e-5, 5e-6, 5e-7, 5e-8]
    # Iterate over PRS files and extract relevant info
    for prs in prs_interest:
        try:
            tmp_path = prs_path / f"{prs}_hits_p5e-5_clumped_reformatted.txt.gz"
            relative_name = (
                str(Path("African") / tmp_path.name)
                if ancestry == "African"
                else tmp_path.name
            )
            tmp_dic = {'id' : prs, 'file_name' : relative_name}
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
    if prs_df.empty:
        return pd.DataFrame()
    # Flag the representative
    prs_df['representative'] = prs_df['id'].apply(lambda x: 'Yes' if x == cluster_rep else 'No')
    # Flag non-curatedness
    prs_df['curated'] = 'No'
    # Merge with meta_sb
    metadata_columns = ['Trait', 'ID', 'Author', 'Year', 'PMID', 'Sample Size']
    available_metadata = meta_sb.reindex(columns=metadata_columns)
    prs_df = prs_df.merge(
        available_metadata,
        left_on='id',
        right_on='ID',
        how='left',
    )
    # Remove ID column
    prs_df = prs_df.drop(columns=['ID'])
    # Take care of curated PRS
    curated_prs = curatedPRS()
    # Check if curated PRS should be considered
    if ancestry == "European" and 'alzheimer' in browse_resolved.lower():
        # Merge with curated
        prs_df = pd.concat([prs_df, curated_prs[curated_prs['Trait'] == "Alzheimer's Disease - Bellenguez et al. - 2022"]], ignore_index=True)
    elif ancestry == "European" and 'parkinson' in browse_resolved.lower():
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
    return render_template("prs.html", trait_options=TRAIT_OPTIONS)

# Route to handle trait search
@prs_bp.route("/traits", methods=["GET", "POST"])
def trait_search():
    meta = current_app.config["GWAS_META"]
    # Remember last selection ---
    if request.method == "POST":
        selection = (request.form.get("trait") or "").strip()
        if selection:
            session["prs_last_trait"] = selection
        else:
            session.pop("prs_last_trait", None)
    else:
        selection = (session.get("prs_last_trait") or "").strip()

    # Initialize variables
    plot_html = None
    meta_rows = []
    prs_rows = []

    # If we have a query (POST or remembered GET), compute page content
    query_label = ""
    if selection:
        selected_row = resolve_trait_selection(selection)
        query_label = f"[{ANCESTRY_CODES.get(selected_row['ancestry'], selected_row['ancestry'])}] {selected_row['name']}"
        # Monitoring: record search
        add_search_to_file(DATA_PATH, query_label.lower())
        artifact = load_semantic_artifacts()
        # Generate plots
        plot_html = plot_haplotype_traits(artifact, selected_row)

        ancestry = selected_row["ancestry"]
        dataset = artifact["datasets"][ancestry]
        selected_index = int(selected_row["dataset_index"])
        selected_cluster = int(dataset["labels"][selected_index])
        cluster_indices = np.where(dataset["labels"] == selected_cluster)[0]
        ancestry_rows = trait_names[trait_names["ancestry"] == ancestry].set_index("dataset_index")
        cluster_rows = ancestry_rows.loc[cluster_indices].copy()
        cluster_study_ids = cluster_rows["gwas_id"].astype(str).tolist()

        # Subset live OpenGWAS metadata using the ancestry-specific study IDs.
        available_ids = [study_id for study_id in cluster_study_ids if study_id in meta.columns]
        if available_ids:
            meta_subset_long = meta.loc[:, available_ids].T.reset_index(drop=False)
            meta_subset_long = meta_subset_long.rename(columns={
                "trait": "Trait",
                "author": "Author",
                "year": "Year",
                "pmid": "PMID",
                "sample_size": "Sample Size",
                "population": "Population",
                "doi": "DOI",
                "id": "ID",
            })
            meta_sb = meta_subset_long.reindex(columns=[
                "Trait", "ID", "Author", "Year", "PMID", "Sample Size", "Population"
            ])
        else:
            meta_sb = pd.DataFrame(columns=[
                "Trait", "ID", "Author", "Year", "PMID", "Sample Size", "Population"
            ])

        # Get the representative from the summary matching the stored labels.
        cluster_summary = artifact["summaries"][ancestry]
        representative = cluster_summary[
            cluster_summary["Cluster Label"].astype(int) == selected_cluster
        ]
        cluster_rep = ""
        if not representative.empty:
            cluster_rep = re.sub(
                r"_hits_p5e-5\.txt(?:\.gz)?$", "",
                str(representative.iloc[0]["Final Representative file"]),
            )
        # Flag cluster representative
        meta_sb["Cluster Representative"] = meta_sb["ID"].apply(lambda x: "Yes" if x == cluster_rep else "No")
        # Get PRS info
        prs_info = getPRSinfo(
            meta_sb,
            cluster_study_ids,
            cluster_rep,
            str(selected_row["name"]),
            ancestry,
        )
        if not prs_info.empty:
            prs_info["download_url"] = prs_info["file_name"].apply(
                lambda filename: url_for("prs.download_prs_file", filename=filename)
                if filename else ""
            )
            prs_info = prs_info.drop(columns=["file_name"], errors="ignore")
        # Prepare tables
        meta_rows = meta_sb.to_dict(orient="records")
        prs_rows = prs_info.to_dict(orient="records")

    return render_template(
        "prs.html",
        trait_options=TRAIT_OPTIONS,
        query=selection,
        query_label=query_label,
        plot_html=plot_html,
        meta_rows=meta_rows,
        prs_rows=prs_rows,
    )

# Serve PRS files
@prs_bp.route("/download/<path:filename>")
def download_prs_file(filename):
    if not filename:
        abort(404)
    return send_from_directory(PRS_DIR, filename, as_attachment=True)
