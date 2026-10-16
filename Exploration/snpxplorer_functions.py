from flask import Flask, request, render_template, send_file, session, make_response, redirect, url_for, jsonify, Response, stream_with_context, Blueprint
from queue import Queue, Empty
from uuid import uuid4
import json
from flask_sqlalchemy import SQLAlchemy
from flask_session import Session
from werkzeug.security import generate_password_hash, check_password_hash
import redis
import matplotlib
from datetime import timedelta
matplotlib.use('Agg')
import logging
import traceback
import polars as pl
from datetime import datetime
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
from scipy.ndimage import rotate as scipy_rotate
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
from pandasgwas import get_variants_by_variant_id
import plotly.express as px
from PRSpage.prs import (
    ANCESTRY_CODES as SEMANTIC_ANCESTRY_CODES,
    trait_names as combined_trait_names,
    load_semantic_artifacts,
    resolve_trait_selection,
    build_semantic_context,
    _cluster_name,
)

# differences between local and server
# local
#from pandasgwas.get_SNPs import get_variants_by_variant_id
#data_path = '/Users/nicco/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/snpXplorer/Data'
#path_server = '/Users/nicco/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/snpXplorer'
# server
data_path = '/Data'
path_server = ''


def normalize_gwas_ancestry(population):
    """Return a compact ancestry code from OpenGWAS population metadata."""
    if population is None or pd.isna(population):
        return 'UNK'

    value = str(population).strip().lower()
    if not value or value in {'nan', 'none', 'na', 'n/a', 'nr', 'unknown'}:
        return 'UNK'

    exact_codes = {
        'afr': 'AFR',
        'eur': 'EUR',
        'eas': 'EAS',
        'sas': 'SAS',
        'amr': 'AMR',
        'mixed': 'MIX',
        'multi-ancestry': 'MIX',
        'multi ancestry': 'MIX',
    }
    if value in exact_codes:
        return exact_codes[value]

    ancestry_markers = {
        'AFR': ('african', 'afro-caribbean', 'afro caribbean', 'black'),
        'EUR': ('european', 'white british', 'white'),
        'EAS': ('east asian', 'chinese', 'japanese', 'korean'),
        'SAS': ('south asian', 'indian', 'pakistani', 'bangladeshi'),
        'AMR': ('hispanic', 'latino', 'latin american', 'ad-mixed american'),
    }
    detected = {
        code
        for code, markers in ancestry_markers.items()
        if any(marker in value for marker in markers)
    }
    if len(detected) > 1 or any(
        marker in value
        for marker in ('mixed', 'multi-ancestry', 'multi ancestry', 'multiple ancestries')
    ):
        return 'MIX'
    if len(detected) == 1:
        return detected.pop()
    return 'OTHER'


def gwas_id_from_selection(selection):
    """Extract an OpenGWAS ID from the label stored by the study selector."""
    return str(selection).rsplit('[', 1)[-1].replace(']', '').strip()


def gwas_ancestry_code(meta, gwas_id):
    """Look up and normalize the population metadata for one OpenGWAS study."""
    if gwas_id not in meta.columns or 'population' not in meta.index:
        return 'UNK'
    return normalize_gwas_ancestry(meta.at['population', gwas_id])


def reference_ancestry_for_studies(studies, meta):
    """Use AFR resources only for exclusively African-ancestry selections."""
    ancestry_codes = {
        gwas_ancestry_code(meta, gwas_id_from_selection(study))
        for study in studies
    }
    return 'African' if ancestry_codes == {'AFR'} else 'European'


def ancestry_resource_path(data_path, resource, ancestry, filename):
    """Resolve ancestry-specific LD/haplotype files with European as fallback."""
    base_path = Path(data_path) / 'databases' / resource
    if str(ancestry).strip().upper() in {'AFR', 'AFRICAN'}:
        base_path = base_path / 'African'
    return base_path / filename

# Global functions
# function to get all haplotypes associated with a trait
def get_haplotypes_by_trait(row, ancestry='European'):
    # read file of interest with polars
    filename = (
        "20260917_ld_clusters_with_gwas_AI_latest_African.tsv.gz"
        if ancestry == 'African'
        else "ld_clusters_with_gwas_AI_latest.tsv.gz"
    )
    df = (
        pl.scan_csv(f"{data_path}/databases/haplotypes/{filename}", separator="\t")
        .filter(pl.col("traits").str.contains(row['name'], literal=True))
        .collect()
    )
    # convert to dataframe
    df_dic = df.to_dicts()
    df_df = pd.DataFrame(df_dic, columns=df.columns)
    return df_df

# function to get trait info
def get_traits_info(data_path, traits_interest, representative, meta, ancestry='European'):
    all_gw_info = []
    all_haplo_df = pd.DataFrame()
    gwas_assoc_df = pd.DataFrame()
    representative_file = str(representative['Final Representative file'].values[0])
    representative_id = re.sub(r'_hits_p5e-5\.txt(?:\.gz)?$', '', representative_file)
    # iterate over df rows
    for index, row in traits_interest.iterrows():
        gwas_id = str(row.get('gwas_id') or re.sub(r'_hits_p5e-5\.txt(?:\.gz)?$', '', str(row['id'])))
        # get gwas info
        gwas_info = meta[gwas_id]
        # convert to dictionary
        gwas_info = dict(gwas_info)
        gwas_info['__is_representative'] = (gwas_id == representative_id)
        # make a dictionary
        all_gw_info.append(gwas_info)
        # check if this is the representative gwas
        if gwas_id == representative_id:
            # in that case, get associations for manhattan plot
            association_file = (
                "20260917_All_indep_gwas_sumstats_AI_hg38_byTrait_African.txt.gz"
                if ancestry == 'African'
                else "All_indep_gwas_sumstats_AI_hg38_byTrait.txt.gz"
            )
            output = subprocess.run(
                ["tabix", f"{data_path}/databases/haplotypes/{association_file}", gwas_id],
                capture_output=True,
                check=True,
                text=True,
            ).stdout.strip()
            gwas_associations = [x.split('\t') for x in output.split('\n') if x]
            # convert to dataframe
            gwas_assoc_df = pd.DataFrame(gwas_associations, columns=["id", "trait", "chr", "position", "rsid", "ea", "nea", "eaf", "beta", "se", "p", "n", "study_id", "position_hg38"])
            # get haplotype info for this gwas
            all_haplo_df = get_haplotypes_by_trait(row, ancestry=ancestry)
    # convert all_gw_info to dataframe
    all_gw_info = pd.DataFrame(all_gw_info)
    # subset of columns to keep
    all_gw_info = all_gw_info[["id", "trait", "pmid", "author", "year", "sample_size", "population", "__is_representative"]]
    all_gw_info.insert(1, 'ancestry', ancestry)
    # pmid as integer NA proof
    all_gw_info['pmid'] = all_gw_info['pmid'].fillna(0).astype(int)
    # sort by sample size (descending) but keep representative trait at the very top
    all_gw_info['sample_size'] = pd.to_numeric(all_gw_info['sample_size'], errors='coerce')
    all_gw_info = all_gw_info.sort_values(
        by=['__is_representative', 'sample_size'],
        ascending=[False, False],
        na_position='last'
    ).reset_index(drop=True)
    # rename columns
    all_gw_info.rename(columns={"id":"GWAS ID", "ancestry":"Ancestry", "trait":"Trait", "pmid":"PubMed ID", "author":"First Author", "year":"Year", "sample_size":"Sample Size", "population":"Population"}, inplace=True)
    ancestry_code = SEMANTIC_ANCESTRY_CODES.get(ancestry, ancestry.upper())
    all_gw_info['Trait'] = all_gw_info['Trait'].apply(
        lambda value: f"[{ancestry_code}] {value}"
    )
    all_haplo_df.rename(columns={'CHR': 'Chromosome', 'BP1':'Start (hg38)', 'BP2':'End (hg38)', 'KB':'Size (kb)', 'NSNPS':'# SNPs', 'SNPS':'SNPs', 'assoc_count':'# Associations', 'rsids':'SNPs (RsIDs)', 'traits':'Traits'}, inplace=True)
    return all_gw_info, gwas_assoc_df, all_haplo_df


def get_semantic_neighborhood_table(artifact, selected_row, scope, meta):
    """Export the studies represented by each cluster in the selected UMAP view."""
    context = build_semantic_context(artifact, selected_row, scope)
    ancestry_rows = {
        ancestry: rows.set_index('dataset_index')
        for ancestry, rows in combined_trait_names.groupby('ancestry')
    }
    records = []
    for cluster_order, group in enumerate(context['points']):
        ancestry, cluster_id = group['key']
        summary = artifact['summaries'][ancestry]
        representative_rows = summary[summary['Cluster Label'].astype(int) == int(cluster_id)]
        representative_id = ''
        if not representative_rows.empty:
            representative_file = representative_rows.iloc[0].get('Final Representative file')
            if pd.notna(representative_file):
                representative_id = re.sub(r'_hits_p5e-5\.txt(?:\.gz)?$', '', str(representative_file))
        representative_trait = _cluster_name(artifact, group['key'])
        code = SEMANTIC_ANCESTRY_CODES.get(ancestry, ancestry.upper())
        # Export every study in each displayed cluster; the UMAP itself caps points at 100.
        cluster_indices = np.flatnonzero(artifact['datasets'][ancestry]['labels'] == cluster_id)
        for index in cluster_indices:
            study = ancestry_rows[ancestry].loc[int(index)]
            gwas_id = str(study['gwas_id'])
            details = meta[gwas_id].to_dict() if gwas_id in meta.columns else {}
            trait = details.get('trait')
            if pd.isna(trait) or not str(trait).strip():
                trait = study['name']
            records.append({
                'GWAS ID': gwas_id,
                'Ancestry': ancestry,
                'Cluster': int(cluster_id),
                'Cluster Representative Trait': representative_trait,
                'Trait': f'[{code}] {trait}',
                'PubMed ID': details.get('pmid'),
                'First Author': details.get('author'),
                'Year': details.get('year'),
                'Sample Size': details.get('sample_size'),
                'Population': details.get('population'),
                '__cluster_order': cluster_order,
                '__is_representative': gwas_id == representative_id,
            })

    columns = [
        'GWAS ID', 'Ancestry', 'Cluster', 'Cluster Representative Trait',
        'Trait', 'PubMed ID', 'First Author', 'Year', 'Sample Size', 'Population',
    ]
    table = pd.DataFrame.from_records(records)
    if table.empty:
        return pd.DataFrame(columns=columns)
    table['Sample Size'] = pd.to_numeric(table['Sample Size'], errors='coerce')
    table = table.sort_values(
        ['__cluster_order', '__is_representative', 'Sample Size'],
        ascending=[True, False, False],
        na_position='last',
    )
    return table[columns].reset_index(drop=True)

# function to extract number from gwas id
def extract_numeric_id(fname):
    try:
        # keep only the part after the second dash
        part = fname.split('-')[2].split('_')[0]
        # extract numeric prefix
        return int(num)
    except:
        return None

def subset_haplotype_rows_for_gwas_overlap(all_haplo_df, gwas_assoc_df):
    if all_haplo_df is None or all_haplo_df.empty or gwas_assoc_df is None or gwas_assoc_df.empty:
        return all_haplo_df.iloc[0:0].copy() if isinstance(all_haplo_df, pd.DataFrame) else pd.DataFrame()

    gwas_tmp = gwas_assoc_df.copy()
    gwas_tmp['chr'] = pd.to_numeric(gwas_tmp['chr'], errors='coerce')
    gwas_tmp['position_hg38'] = pd.to_numeric(gwas_tmp['position_hg38'], errors='coerce')
    gwas_tmp = gwas_tmp.dropna(subset=['chr', 'position_hg38']).copy()
    gwas_tmp['chr'] = gwas_tmp['chr'].astype(int)
    gwas_tmp['position_hg38'] = gwas_tmp['position_hg38'].astype(int)
    gwas_pos_by_chr = (
        gwas_tmp.groupby('chr')['position_hg38']
        .apply(lambda s: set(s.tolist()))
        .to_dict()
    )

    def haplotype_overlaps_gwas(row):
        snps_value = row.get('SNPs')
        if pd.isna(snps_value):
            return False
        snps = [s.strip() for s in str(snps_value).split('|') if s and s.strip()]
        if not snps:
            return False
        first_parts = snps[0].split(':')
        if len(first_parts) < 2:
            return False
        chrom_token = first_parts[0].strip().upper().replace('CHR', '')
        try:
            chrom = int(chrom_token)
        except ValueError:
            return False
        if chrom not in gwas_pos_by_chr:
            return False
        gwas_positions = gwas_pos_by_chr[chrom]
        for snp in snps:
            parts = snp.split(':')
            if len(parts) < 2:
                continue
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            if pos in gwas_positions:
                return True
        return False

    return all_haplo_df[all_haplo_df.apply(haplotype_overlaps_gwas, axis=1)].copy()

# function to plot haplotype traits
def plot_haplotype_traits(clusters_of_interest, all_indices, trait_list, cluster_index_map, umap, cluster_representatives, browse_label, sim_mat, names_sub, cluster_labels_for_index, cluster_trait, idx_trait, all_cluster_representatives, gwas_assoc_df, all_haplo_df, artifact, selected_row, combined_context):
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
            "UMAP semantic trait projection",
            "Semantic trait similarity",
            "Manhattan plot of representative GWAS"                     # Title for row 2
        )
    )
    # increase subplot subtitle font size
    fig.update_annotations(font=dict(size=20, family="Arial", color="darkgrey", weight="bold"))
    ancestry = selected_row['ancestry']

    # -------------------------------
    # LEFT subplot: UMAP
    # -------------------------------
    primary_context = build_semantic_context(artifact, selected_row, ancestry)
    semantic_contexts = [primary_context, combined_context]
    semantic_traces_by_context = []
    for context_index, context in enumerate(semantic_contexts):
        visible = context_index == 0
        trace_indices = []
        context_color_map = {
            point_group['key']: palette[index % len(palette)]
            for index, point_group in enumerate(context['points'])
        }
        for point_index, point_group in enumerate(context['points']):
            cluster_key = point_group['key']
            fig.add_trace(
                go.Scatter(
                    x=point_group['coordinates'][:, 0],
                    y=point_group['coordinates'][:, 1],
                    mode='markers',
                    marker=dict(size=9, color=context_color_map[cluster_key]),
                    name=_cluster_name(artifact, cluster_key),
                    text=point_group['names'],
                    hoverinfo='text',
                    legendrank=10 + point_index,
                    visible=visible,
                ),
                row=1,
                col=1,
            )
            trace_indices.append(len(fig.data) - 1)

        fig.add_trace(
            go.Scatter(
                x=[context['query_coordinate'][0]],
                y=[context['query_coordinate'][1]],
                mode='markers',
                marker=dict(size=15, color='black', symbol='star', line=dict(width=2, color='yellow')),
                name='Queried trait',
                hovertext=browse_label,
                hoverinfo='text',
                legendrank=0,
                visible=visible,
            ),
            row=1,
            col=1,
        )
        trace_indices.append(len(fig.data) - 1)
        semantic_traces_by_context.append(trace_indices)
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
    cluster_colors_for_traits = [color_map.get(cid, "#d3d3d3") for cid in cluster_labels_for_index]
    # Build a 1-column heatmap strip so width can be increased in X only
    unique_strip_colors = list(dict.fromkeys(cluster_colors_for_traits))
    color_to_idx = {c: i for i, c in enumerate(unique_strip_colors)}
    strip_width_cols = 4
    strip_col = np.array([color_to_idx[c] for c in cluster_colors_for_traits], dtype=float).reshape(-1, 1)
    strip_vals = np.repeat(strip_col, strip_width_cols, axis=1)
    if len(unique_strip_colors) == 1:
        strip_colorscale = [[0.0, unique_strip_colors[0]], [1.0, unique_strip_colors[0]]]
    else:
        strip_colorscale = []
        denom = float(len(unique_strip_colors) - 1)
        for i, c in enumerate(unique_strip_colors):
            p = i / denom
            strip_colorscale.append([p, c])
    fig.add_trace(
        go.Heatmap(
            z=strip_vals,
            # first strip column centered at n_traits: left edge is n_traits-0.5 (touches heatmap edge)
            x=list(range(n_traits, n_traits + strip_width_cols)),
            y=list(range(n_traits)),
            zmin=0,
            zmax=max(1, len(unique_strip_colors) - 1),
            colorscale=strip_colorscale,
            showscale=False,
            hoverinfo="skip",
        ),
        row=1,
        col=2
    )
    # Hide axis labels / ticks for heatmap panel
    fig.update_xaxes(
        showticklabels=False,
        title_text="",
        # include the strip columns exactly; no right-side gap
        range=[-0.5, n_traits + strip_width_cols - 0.5],
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
    plot_haplo_df = all_haplo_df.iloc[0:0].copy() if isinstance(all_haplo_df, pd.DataFrame) else pd.DataFrame()
    if not gwas_assoc_df.empty:
        # Prepare data
        gwas_assoc_df['position_hg38'] = gwas_assoc_df['position_hg38'].astype(int)
        # Replace zero p-values with the smallest non-zero p in the dataframe (fallback to 1e-300)
        p_floats = gwas_assoc_df['p'].astype(float)
        zeros_mask = p_floats == 0
        if zeros_mask.any():
            nonzero_vals = p_floats[~zeros_mask]
            min_nonzero = nonzero_vals.min() if not nonzero_vals.empty else 1e-300
            gwas_assoc_df.loc[zeros_mask, 'p'] = min_nonzero
        gwas_assoc_df['neg_log10_p'] = -np.log10(gwas_assoc_df['p'].astype(float))
        # Sort by pvalue
        gwas_assoc_df['p'] = gwas_assoc_df['p'].astype(float)
        gwas_assoc_df = gwas_assoc_df.sort_values(by='p').reset_index(drop=True)
        # Take only p<5e-08
        gwas_assoc_df = gwas_assoc_df[gwas_assoc_df['p'] < 5e-08].copy().reset_index(drop=True)
        if gwas_assoc_df.shape[0] > 10000:
            gwas_assoc_df = gwas_assoc_df.iloc[:10000, :].copy().reset_index(drop=True)
        # Sort by chromosome and position
        gwas_assoc_df['chr'] = gwas_assoc_df['chr'].astype(int)
        gwas_assoc_df['position_hg38'] = gwas_assoc_df['position_hg38'].astype(int)
        gwas_assoc_df = gwas_assoc_df.sort_values(by=['chr', 'position_hg38']).reset_index(drop=True)
        # Keep only haplotypes containing at least one SNP overlapping a plotted GWAS association
        plot_haplo_df = subset_haplotype_rows_for_gwas_overlap(all_haplo_df, gwas_assoc_df)
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
        # Now prepare the haplotype associations
        plot_haplo_df['Chromosome'] = plot_haplo_df['Chromosome'].astype(int)
        plot_haplo_df['Start (hg38)'] = plot_haplo_df['Start (hg38)'].astype(int)
        plot_haplo_df['End (hg38)'] = plot_haplo_df['End (hg38)'].astype(int)
        plot_haplo_df['midpoint'] = (plot_haplo_df['Start (hg38)'] + plot_haplo_df['End (hg38)']) // 2
        # Calculate cumulative positions for haplotypes
        cum_length = 0
        haplo_cum_positions = []
        for chrom in range(1, 23):
            chrom_haplo = plot_haplo_df[plot_haplo_df['Chromosome'] == chrom]
            for pos in chrom_haplo['midpoint']:
                haplo_cum_positions.append(pos + cum_length)
            cum_length += chrom_lengths[chrom]
        # sort haplotypes by chromosome and Start (hg38)
        plot_haplo_df = plot_haplo_df.sort_values(by=['Chromosome', 'Start (hg38)']).reset_index(drop=True)
        plot_haplo_df['cum_position'] = haplo_cum_positions
        # calculate distribution -- choose bin
        bin_size = 10_000_000  # 10 Mb bins
        plot_haplo_df['bin'] = (plot_haplo_df['cum_position'] // bin_size).astype(int)
        # Add haplo ID
        plot_haplo_df['haplo_id'] = plot_haplo_df['Chromosome'].astype(str) + ':' + plot_haplo_df['Start (hg38)'].astype(str) + '-' + plot_haplo_df['End (hg38)'].astype(str)
        # Count haplotypes per bin
        haplo_bin_counts = (
            plot_haplo_df
            .groupby(["bin", "Chromosome"])
            .agg(
                haplo_count = ("haplo_id", "count"),
                haplo_ids   = ("haplo_id", list)
            )
            .reset_index()
        )
        # Center bars on bin midpoints so they align with SNP towers in the same bins
        haplo_bin_counts['genomic_pos'] = (haplo_bin_counts['bin'] * bin_size) + (bin_size / 2.0)
        # Add number of haplo_ids per bin
        haplo_bin_counts['n_haplo_ids'] = haplo_bin_counts['haplo_ids'].apply(len)
        # Add colors
        haplo_bin_counts['color'] = haplo_bin_counts['Chromosome'] % 2
        haplo_bin_counts['color'] = haplo_bin_counts['color'].map({0: colors[0], 1: colors[1]})
        # Scale the y otherwise nothing is visible
        max_haplo_count = haplo_bin_counts['haplo_count'].max() if not haplo_bin_counts.empty else 1
        # scale -negative log10 pvalues to max haplo count
        gwas_assoc_df['neg_log10_p'] = gwas_assoc_df['neg_log10_p'] * (max_haplo_count / gwas_assoc_df['neg_log10_p'].abs().max())
        # Negate the -log10 pvalues for better visualization
        gwas_assoc_df['neg_log10_p'] = -gwas_assoc_df['neg_log10_p']
        gwas_assoc_df['orig_neg_log10_p'] = -np.log10(gwas_assoc_df['p'].astype(float))  # keep original for hovertext
        # Bar plot
        fig.add_trace(
            go.Bar(
                x=haplo_bin_counts['genomic_pos'],
                y=haplo_bin_counts['haplo_count'],
                width=bin_size * 0.9,
                marker_color=haplo_bin_counts['color'],
                marker_line=dict(width=1, color=haplo_bin_counts['color']),
                hovertext=[
                    f'Number of Haplotypes: {row["haplo_count"]}'
                    for _, row in haplo_bin_counts.iterrows()
                ],
                hoverinfo='text',
                showlegend=False,
            ), row=2, col=1
        )
        # thin line at y=0
        fig.add_hline(y=0, line=dict(color="black", width=1), row=2, col=1)
        # Add GWAS associations
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
                hoverinfo='text',
                showlegend=False,
            ),
            row=2,
            col=1
        )
        # X-axis ticks at chromosome centers
        ticks = gwas_assoc_df.groupby("chr")["cum_position"].mean()
        tick_labels = gwas_assoc_df.groupby("chr")["chr"].first()
        fig.update_xaxes(title_text="Chromosomes (hg38)", row=2, col=1, title_font=dict(size=20, family="Arial", color="black", weight="bold"), tickmode="array", tickvals=ticks, ticktext=tick_labels)
        # vertical lines between chromosomes
        cum_length = 0
        for chrom in range(1, 23):
            cum_length += chrom_lengths[chrom]
            fig.add_vline(x=cum_length, line=dict(color="lightgrey", width=1, dash="dash"), row=2, col=1)
        # custom title for top y-axis at y=max_haplo_count/2
        fig.add_annotation(
            xref="x3", yref="y3",
            x=0,                      # left of plot
            y=max_haplo_count/2,            # halfway up positive space
            text="<b>Haplotype count</b>",
            textangle=90,
            showarrow=False,
            font=dict(size=14),
            align="left"
        )
        fig.add_annotation(
            xref="x3", yref="y3",
            x=0,
            y=-max_haplo_count/2,            # halfway down negative space
            text="<b>-log10(pvalue)</b>",
            showarrow=False,
            textangle=90,
            font=dict(size=14),
            align="left"
        )
        # y axis
        p_ticks = np.linspace(-max_haplo_count, max_haplo_count, 6).tolist()
        scale = (-gwas_assoc_df["neg_log10_p"] / gwas_assoc_df["orig_neg_log10_p"]).mean()
        p_names = [(-val / scale) if val < 0 else val for val in p_ticks]
        fig.update_yaxes(
            showline=True,
            ticks="outside",
            side="right",
            tickvals = p_ticks,
            ticktext = [f"{val:.0f}" for val in p_names],
            row=2, col=1,    # your bottom plot
        )
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
        height=1300,
        title={
            "text": f"<b>Semantic Neighborhood of '{selected_row['name']}'</b>",
            "font": {"size": 32},
            "x": 0.5,
            "xanchor": "center",
            "y": 0.99,
            "yanchor": "top",
        },
        margin=dict(t=150, b=80),
    )
    semantic_trace_indices = set(
        index
        for trace_indices in semantic_traces_by_context
        for index in trace_indices
    )
    visibility_by_scope = {}
    for scope, active_indices in zip(
        ('primary', 'combined'), semantic_traces_by_context
    ):
        visibility = [index not in semantic_trace_indices for index in range(len(fig.data))]
        for index in active_indices:
            visibility[index] = True
        visibility_by_scope[scope] = visibility
    fig.update_layout(meta={'semantic_scope_visibility': visibility_by_scope})
    # return as HTML
    plot_url = pio.to_html(fig, include_plotlyjs='cdn', full_html=False)
    # subset of the gwas_assoc_df to return
    gwas_assoc_df = gwas_assoc_df[['id', 'trait', 'chr', 'position_hg38', 'rsid', 'ea', 'nea', 'eaf', 'beta', 'se', 'p', 'n']].copy()
    # rename columns
    gwas_assoc_df.rename(columns={'id':'GWAS ID', 'trait':'Trait', 'chr':'Chromosome', 'position_hg38':'Position (hg38)', 'rsid':'RsID', 'ea':'Effect Allele', 'nea':'Other Allele', 'eaf':'EAF', 'beta':'Beta', 'se':'SE', 'p':'P-value', 'n':'N'}, inplace=True)
    # and sort by p-value
    gwas_assoc_df = gwas_assoc_df.sort_values(by='P-value').reset_index(drop=True)
    return plot_url, gwas_assoc_df, plot_haplo_df

# function to extract haplotype effect on traits
def extract_haplotype_effect_on_traits(data_path, haplo_summary, chrom, ancestry='European'):
    # subset to haplotypes with >1 trait associations
    haplo_tmp = haplo_summary[haplo_summary['n_traits'] >= 1].copy().reset_index(drop=True)
    if haplo_tmp.empty:
        return None
    # if more than 10, restrict to top 10 by number of traits
    if haplo_tmp.shape[0] > 10:
        haplo_tmp = haplo_tmp.sort_values(by='n_traits', ascending=False).iloc[:10, :].copy().reset_index(drop=True)
    # extract all haplo data
    start_pos = min(haplo_tmp['start_bp'].astype(int))
    end_pos = max(haplo_tmp['end_bp'].astype(int))
    haplo_df, haplo_dict = extract_haplo_data(
        data_path,
        chrom,
        start_pos-1_000_000,
        end_pos+1_000_000,
        ancestry=ancestry,
    )
    association_file = (
        '20260917_All_indep_gwas_sumstats_AI_hg38_African.txt.gz'
        if str(ancestry).strip().upper() in {'AFR', 'AFRICAN'}
        else 'All_indep_gwas_sumstats_AI_hg38.txt.gz'
    )
    # iterate over haplotypes
    haplo_effects = []
    for idx, row in haplo_tmp.iterrows():
        haplo_id = row['ID']
        # subset to haplotype of interest
        haplo_row = haplo_df[haplo_df['ID'] == haplo_id].copy().reset_index(drop=True)
        # extract snps
        snps = ' '.join([str(chrom).upper().replace('CHR', '') + ':' + x.split(':')[1] + '-' + x.split(':')[1] for x in haplo_dict[haplo_id]])
        # extract effects
        cmd = "tabix %s/databases/haplotypes/%s %s" %(
            data_path.replace(' ', '\ '),
            association_file,
            snps,
        )
        snps_data = [x.rstrip().split('\t') for x in os.popen(cmd)]
        snps_df = pd.DataFrame(snps_data, columns=["id", "trait", "chr", "position", "rsid", "ea", "nea", "eaf", "beta", "se", "p", "n", "study_id", "position_hg38"])
        # Make effect estimates numeric and derive a signed association Z-score.
        # Invalid or zero standard errors cannot yield a finite Z-score and are
        # therefore treated as missing in the standardized view.
        snps_df['beta'] = pd.to_numeric(snps_df['beta'], errors='coerce')
        snps_df['se'] = pd.to_numeric(snps_df['se'], errors='coerce')
        valid_se = snps_df['se'].notna() & np.isfinite(snps_df['se']) & (snps_df['se'] != 0)
        snps_df['z'] = np.nan
        snps_df.loc[valid_se, 'z'] = snps_df.loc[valid_se, 'beta'] / snps_df.loc[valid_se, 'se']
        snps_df.loc[~np.isfinite(snps_df['z']), 'z'] = np.nan
        # Average raw effects and signed Z-scores across the SNP associations
        # available for each haplotype-trait pair. Both are unweighted means.
        trait_effects = (
            snps_df.groupby('trait')
            .agg(
                beta=('beta', 'mean'),
                z=('z', 'mean'),
                n_beta=('beta', 'count'),
                n_z=('z', 'count'),
            )
            .reset_index()
        )
        # add haplo id
        trait_effects['haplo_id'] = haplo_id
        # append to list
        haplo_effects.append(trait_effects)
    # combine into a single dataframe
    haplo_effects = pd.concat(haplo_effects, ignore_index=True)
    # add haplotype index by merging with haplo_tmp
    haplo_effects = haplo_effects.merge(haplo_tmp[['ID', 'haplo_index']], left_on='haplo_id', right_on='ID', how='left')
    # Build matching matrices so the figure can switch between raw effects and
    # signed association Z-scores without changing the displayed traits.
    heatmap_beta = haplo_effects.pivot(index="haplo_index", columns="trait", values="beta")
    heatmap_z = haplo_effects.pivot(index="haplo_index", columns="trait", values="z")
    heatmap_n_beta = haplo_effects.pivot(index="haplo_index", columns="trait", values="n_beta")
    heatmap_n_z = haplo_effects.pivot(index="haplo_index", columns="trait", values="n_z")
    # preserve haplotype order from the original dataframe
    haplo_order = haplo_effects["haplo_index"].drop_duplicates()
    heatmap_beta = heatmap_beta.reindex(haplo_order)
    heatmap_z = heatmap_z.reindex(index=haplo_order, columns=heatmap_beta.columns)
    heatmap_n_beta = heatmap_n_beta.reindex(index=haplo_order, columns=heatmap_beta.columns)
    heatmap_n_z = heatmap_n_z.reindex(index=haplo_order, columns=heatmap_beta.columns)
    # If >40 traits, limit both views to the same top 40 by number of
    # associated haplotypes.
    if heatmap_beta.shape[1] > 40:
        trait_counts = haplo_effects.groupby("trait")["haplo_id"].nunique().reset_index(name="haplo_count")
        top_traits = trait_counts.sort_values(by="haplo_count", ascending=False).head(40)["trait"].tolist()
        heatmap_beta = heatmap_beta.loc[:, top_traits]
        heatmap_z = heatmap_z.loc[:, top_traits]
        heatmap_n_beta = heatmap_n_beta.loc[:, top_traits]
        heatmap_n_z = heatmap_n_z.loc[:, top_traits]

    trait_labels = [col if len(col) <= 40 else col[:40] + "..." for col in heatmap_beta.columns]

    def symmetric_limit(matrix):
        values = matrix.to_numpy(dtype=float)
        finite_values = values[np.isfinite(values)]
        if finite_values.size == 0:
            return 1.0
        limit = float(np.max(np.abs(finite_values)))
        return limit if limit > 0 else 1.0

    max_abs_beta = symmetric_limit(heatmap_beta)
    max_abs_z = symmetric_limit(heatmap_z)
    beta_values = heatmap_beta.to_numpy(dtype=float)
    z_values = heatmap_z.to_numpy(dtype=float)
    customdata = np.dstack((
        beta_values,
        z_values,
        heatmap_n_beta.to_numpy(dtype=float),
        heatmap_n_z.to_numpy(dtype=float),
    ))
    beta_na_text = np.where(heatmap_beta.isna(), "/", "")
    z_na_text = np.where(heatmap_z.isna(), "/", "")

    fig = go.Figure(
        data=go.Heatmap(
            z=z_values,
            x=trait_labels,
            y=heatmap_z.index.tolist(),
            customdata=customdata,
            colorscale="RdBu_r",
            zmin=-max_abs_z,
            zmax=max_abs_z,
            zmid=0,
            xgap=1,
            ygap=1,
            text=z_na_text,
            texttemplate="%{text}",
            textfont=dict(color="lightgrey", size=12),
            colorbar=dict(
                title="Mean signed Z-score",
                orientation="h",
                x=0.5,
                xanchor="center",
                y=1.08,
                yanchor="bottom",
                len=0.3,
                thickness=15,
            ),
            hovertemplate=(
                "Haplotype: %{y}<br>"
                "Trait: %{x}<br>"
                "Mean raw effect: %{customdata[0]:.4g}<br>"
                "Mean signed Z-score: %{customdata[1]:.4g}<br>"
                "Associations contributing to raw mean: %{customdata[2]:.0f}<br>"
                "Associations contributing to Z-score mean: %{customdata[3]:.0f}"
                "<extra></extra>"
            ),
        )
    )
    fig.update_layout(
        height=max(500, 35 * heatmap_beta.shape[0]),
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(t=150),
        xaxis=dict(tickangle=45, automargin=True),
        yaxis=dict(
            automargin=True,
            tickmode="array",
            tickvals=heatmap_beta.index.tolist(),
            ticktext=heatmap_beta.index.tolist(),
        ),
        updatemenus=[dict(
            type="buttons",
            direction="right",
            active=0,
            x=1.0,
            xanchor="right",
            y=1.22,
            yanchor="bottom",
            buttons=[
                dict(
                    label="Association Z-score",
                    method="restyle",
                    args=[{
                        "z": [z_values],
                        "zmin": [-max_abs_z],
                        "zmax": [max_abs_z],
                        "text": [z_na_text],
                        "colorbar.title.text": ["Mean signed Z-score"],
                    }],
                ),
                dict(
                    label="Raw effect",
                    method="restyle",
                    args=[{
                        "z": [beta_values],
                        "zmin": [-max_abs_beta],
                        "zmax": [max_abs_beta],
                        "text": [beta_na_text],
                        "colorbar.title.text": ["Mean raw effect"],
                    }],
                ),
            ],
        )],
    )
    fig.update_yaxes(
        autorange="reversed",
        title_text="Haplotype",
        title_font=dict(size=20, family="Arial", color="black", weight="bold"),
        tickfont=dict(size=16, family="Arial", color="black"),
    )
    # convert to HTML
    plot_url = pio.to_html(fig, include_plotlyjs='cdn', full_html=False)
    return plot_url

# function to guide haplotype plots with traits as inputs
def guide_haplotypes_traits(data_path, browse, window, refGen, trait_names):
    artifact = load_semantic_artifacts()
    selected_row = resolve_trait_selection(browse)
    ancestry = selected_row['ancestry']
    ancestry_code = SEMANTIC_ANCESTRY_CODES.get(ancestry, ancestry.upper())
    browse_label = f"[{ancestry_code}] {selected_row['name']}"
    browse_resolved = selected_row['selection_id']

    ancestry_traits = (
        combined_trait_names[combined_trait_names['ancestry'] == ancestry]
        .sort_values('dataset_index')
        .reset_index(drop=True)
    )
    trait_list = ancestry_traits['name'].values.tolist()
    dataset = artifact['datasets'][ancestry]
    umap = dataset['umap']
    embeddings = dataset['embeddings']
    labels = dataset['labels']
    idx_trait = int(selected_row['dataset_index'])
    # Find cluster id
    cluster_trait = int(labels[idx_trait])
    # Build cluster to index map
    cluster_index_map = build_cluster_index_map(labels)
    primary_context = build_semantic_context(artifact, selected_row, ancestry)
    clusters_of_interest = [group['key'][1] for group in primary_context['points']]
    # ------------------------------------------------
    # Build similarity heatmap between traits in clusters
    # ------------------------------------------------
    all_indices = []
    for point_group in primary_context['points']:
        all_indices.extend(point_group['indices'].tolist())
    all_indices = np.array(all_indices)
    # Also get the list of trait in the cluster of interest
    cluster_member_indices = set(cluster_index_map[cluster_trait])
    traits_sub = trait_names[
        (trait_names['ancestry'] == ancestry)
        & (trait_names['dataset_index'].isin(cluster_member_indices))
    ].copy().reset_index(drop=True)
    cluster_representatives = artifact['summaries'][ancestry]
    representative = cluster_representatives[
        cluster_representatives['Cluster Label'].astype(int) == cluster_trait
    ]
    # Use the same de-duplicated, ancestry-specific heatmap construction as
    # the PRS page. Only the UMAP is allowed to switch to the joint projection.
    sim_mat = primary_context['similarity_matrix']
    names_sub = primary_context['heatmap_names']
    cluster_labels_for_index = np.array([
        cluster_key[1] for cluster_key in primary_context['heatmap_keys']
    ])
    combined_context = build_semantic_context(artifact, selected_row, 'Combined')
    return (
        clusters_of_interest,
        all_indices,
        trait_list,
        cluster_index_map,
        umap,
        representative,
        browse_resolved,
        browse_label,
        sim_mat,
        names_sub,
        cluster_labels_for_index,
        cluster_trait,
        traits_sub,
        idx_trait,
        cluster_representatives,
        artifact,
        selected_row,
        combined_context,
    )
    
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

# function to convert a zero-based index to spreadsheet-like letters
def index_to_letters(idx):
    # 0 -> A, 25 -> Z, 26 -> AA, etc.
    idx = int(idx)
    letters = ""
    while True:
        idx, rem = divmod(idx, 26)
        letters = chr(65 + rem) + letters
        if idx == 0:
            break
        idx -= 1
    return letters

# function to summarize trait repetitions within a haplotype
def summarize_traits_with_counts(traits_value):
    """
    Convert strings like:
    A___B___A___C
    into:
    A (2 SNPs); B (1 SNPs); C (1 SNPs)
    Preserves first-seen order.
    """
    if pd.isna(traits_value):
        return ""
    raw = [t.strip() for t in str(traits_value).split("___")]
    raw = [t for t in raw if t and t.lower() != "nan"]
    if not raw:
        return ""
    counts = {}
    order = []
    for t in raw:
        if t not in counts:
            counts[t] = 0
            order.append(t)
        counts[t] += 1
    return "; ".join([f"{t} ({counts[t]} SNPs)" for t in order])


def build_haplotype_ld_matrix(
    data_path,
    chrom,
    haplo_df,
    haplo_label_map,
    max_variants=1000,
    ancestry='European',
):
    """Build an allele-aware LD matrix ordered by haplotype and position."""
    ordered_haplotypes = list(haplo_label_map.items())
    haplotype_rank = {hap_id: rank for rank, (hap_id, _) in enumerate(ordered_haplotypes)}
    label_by_haplotype = dict(ordered_haplotypes)

    variants = {}
    selected = haplo_df[haplo_df['ID'].isin(haplotype_rank)].copy()
    selected['__rank'] = selected['ID'].map(haplotype_rank)
    selected = selected.sort_values(['__rank', 'BP1'])
    for _, row in selected.iterrows():
        hap_id = row['ID']
        for variant_id in str(row['SNPS']).split('|'):
            parts = variant_id.split(':')
            if len(parts) != 4:
                continue
            variant_chrom, position, ref, alt = parts
            try:
                position = int(position)
            except ValueError:
                continue
            key = (str(variant_chrom).upper().replace('CHR', ''), position, ref.upper(), alt.upper())
            if key not in variants:
                variants[key] = {
                    'variant_id': variant_id,
                    'position': position,
                    'ref': ref.upper(),
                    'alt': alt.upper(),
                    'haplotype_id': hap_id,
                    'haplotype_label': label_by_haplotype[hap_id],
                    'haplotype_rank': haplotype_rank[hap_id],
                    'db_id': None,
                }

    ordered_variants = sorted(
        variants.values(),
        key=lambda item: (
            item['haplotype_rank'],
            item['position'],
            item['ref'],
            item['alt'],
        ),
    )
    if len(ordered_variants) > max_variants:
        raise ValueError(
            f"Regional LD contains {len(ordered_variants):,} variants; "
            f"the interactive limit is {max_variants:,}."
        )
    if len(ordered_variants) < 2:
        raise ValueError("At least two variants are required for the regional LD plot.")

    db_path = ancestry_resource_path(
        data_path,
        'LD_db',
        ancestry,
        f"ld_chr{str(chrom).upper().replace('CHR', '')}.sqlite",
    )
    if not db_path.is_file():
        raise ValueError(f"LD database is unavailable for chromosome {chrom}.")

    positions = sorted({variant['position'] for variant in ordered_variants})
    rows_by_position = {}
    connection = sqlite3.connect(f"file:{db_path.resolve()}?mode=ro", uri=True)
    try:
        for offset in range(0, len(positions), 900):
            position_chunk = positions[offset:offset + 900]
            placeholders = ','.join('?' for _ in position_chunk)
            query_sql = f"SELECT id, uniq, pos FROM variants WHERE pos IN ({placeholders})"
            for db_id, uniq, position in connection.execute(query_sql, position_chunk):
                rows_by_position.setdefault(int(position), []).append((int(db_id), str(uniq)))

        for variant in ordered_variants:
            candidates = []
            for db_id, uniq in rows_by_position.get(variant['position'], []):
                uniq_parts = uniq.split(':')
                if len(uniq_parts) != 3:
                    continue
                observed_alleles = (uniq_parts[1].upper(), uniq_parts[2].upper())
                expected = (variant['ref'], variant['alt'])
                if observed_alleles in {expected, expected[::-1]}:
                    candidates.append((db_id, uniq))
            if candidates:
                variant['db_id'] = min(candidates, key=lambda item: item[0])[0]

        matrix_size = len(ordered_variants)
        matrix = np.full((matrix_size, matrix_size), np.nan, dtype=float)
        np.fill_diagonal(matrix, 1.0)
        matrix_index_by_db_id = {
            variant['db_id']: index
            for index, variant in enumerate(ordered_variants)
            if variant['db_id'] is not None
        }
        selected_db_ids = set(matrix_index_by_db_id)
        seen_pairs = set()
        for db_id in sorted(selected_db_ids):
            ld_rows = connection.execute(
                "SELECT v1, v2, r2_milli FROM ld WHERE v1 = ? "
                "UNION ALL "
                "SELECT v1, v2, r2_milli FROM ld WHERE v2 = ?",
                (db_id, db_id),
            )
            for v1, v2, r2_milli in ld_rows:
                v1 = int(v1)
                v2 = int(v2)
                if v1 not in selected_db_ids or v2 not in selected_db_ids or v1 == v2:
                    continue
                pair = (min(v1, v2), max(v1, v2))
                if pair in seen_pairs:
                    continue
                seen_pairs.add(pair)
                row_index = matrix_index_by_db_id[v1]
                column_index = matrix_index_by_db_id[v2]
                r2 = int(r2_milli) / 1000.0
                matrix[row_index, column_index] = r2
                matrix[column_index, row_index] = r2
    finally:
        connection.close()

    return matrix, ordered_variants


def haplotype_ld_roof_plotly(matrix, ordered_variants):
    """Rotate a lower LD matrix into a Plotly roof with a haplotype strip."""
    observed_matrix = np.isfinite(matrix)
    display_matrix = np.nan_to_num(matrix, nan=0.0)
    lower_matrix = np.where(
        np.tri(matrix.shape[0], matrix.shape[1], dtype=bool),
        display_matrix,
        np.nan,
    )
    observed_lower = np.where(
        np.tri(matrix.shape[0], matrix.shape[1], dtype=bool),
        observed_matrix.astype(float),
        np.nan,
    )
    rotated = scipy_rotate(
        lower_matrix,
        45,
        reshape=True,
        order=0,
        mode='constant',
        cval=np.nan,
        prefilter=False,
    )
    rotated_observed = scipy_rotate(
        observed_lower,
        45,
        reshape=True,
        order=0,
        mode='constant',
        cval=np.nan,
        prefilter=False,
    )

    diagonal_marker = np.full(matrix.shape, np.nan, dtype=float)
    np.fill_diagonal(diagonal_marker, 1.0)
    rotated_diagonal = scipy_rotate(
        diagonal_marker,
        45,
        reshape=True,
        order=0,
        mode='constant',
        cval=np.nan,
        prefilter=False,
    )
    diagonal_counts = np.sum(np.isfinite(rotated_diagonal), axis=1)
    diagonal_row = int(np.argmax(diagonal_counts))
    roof = rotated[diagonal_row:, :]
    roof_observed = rotated_observed[diagonal_row:, :]
    diagonal_columns = np.where(np.isfinite(rotated_diagonal[diagonal_row]))[0]
    if len(diagonal_columns) < 2:
        raise ValueError("The regional LD matrix could not be transformed for plotting.")

    if not np.any(np.isfinite(roof)):
        raise ValueError("No LD values are available for the selected haplotypes.")
    roof_height, roof_width = roof.shape
    diagonal_start = int(diagonal_columns.min())
    diagonal_end = int(diagonal_columns.max())
    diagonal_x = np.linspace(diagonal_start, diagonal_end, len(ordered_variants))

    haplotype_labels = [variant['haplotype_label'] for variant in ordered_variants]
    label_order = list(dict.fromkeys(haplotype_labels))
    label_rank = {label: index for index, label in enumerate(label_order)}
    annotation = np.full((1, roof_width), np.nan, dtype=float)
    group_centers = []
    group_boundaries = []
    group_start = 0
    for index in range(1, len(haplotype_labels) + 1):
        if index == len(haplotype_labels) or haplotype_labels[index] != haplotype_labels[group_start]:
            start_x = diagonal_x[group_start]
            end_x = diagonal_x[index - 1]
            left = int(round(start_x if group_start == 0 else (diagonal_x[group_start - 1] + start_x) / 2))
            right = int(round(end_x if index == len(haplotype_labels) else (end_x + diagonal_x[index]) / 2))
            annotation[0, max(0, left):min(roof_width, right + 1)] = label_rank[haplotype_labels[group_start]]
            group_centers.append((start_x + end_x) / 2)
            if index < len(haplotype_labels):
                group_boundaries.append((end_x + diagonal_x[index]) / 2)
            group_start = index

    ld_colorscale = [
        [0.0, 'navy'], [0.2, 'navy'],
        [0.20001, 'lightblue'], [0.4, 'lightblue'],
        [0.40001, 'green'], [0.6, 'green'],
        [0.60001, 'orange'], [0.8, 'orange'],
        [0.80001, 'red'], [1.0, 'red'],
    ]
    haplotype_colors = pc.qualitative.Plotly + pc.qualitative.Safe + pc.qualitative.Set3
    n_labels = max(1, len(label_order))
    if n_labels == 1:
        annotation_colorscale = [
            [0.0, haplotype_colors[0]],
            [1.0, haplotype_colors[0]],
        ]
    else:
        annotation_colorscale = []
        for index in range(n_labels):
            start = max(0.0, (index - 0.5) / (n_labels - 1))
            end = min(1.0, (index + 0.5) / (n_labels - 1))
            color = haplotype_colors[index % len(haplotype_colors)]
            annotation_colorscale.extend([[start, color], [end, color]])

    hover_text = np.full(roof.shape, '', dtype=object)
    plotted_cells = np.isfinite(roof)
    observed_cells = plotted_cells & (roof_observed >= 0.5)
    missing_cells = plotted_cells & ~observed_cells
    hover_text[observed_cells] = [
        f'LD r²: {value:.3f}' for value in roof[observed_cells]
    ]
    hover_text[missing_cells] = 'LD r²: NA (not stored)'

    figure = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.01,
        row_heights=[0.94, 0.06],
    )
    figure.add_trace(
        go.Heatmap(
            z=roof,
            zmin=0,
            zmax=1,
            colorscale=ld_colorscale,
            showscale=True,
            colorbar=dict(
                title=dict(text='LD r²', font=dict(size=17)),
                orientation='h',
                x=0.5,
                xanchor='center',
                y=1.08,
                yanchor='bottom',
                len=0.52,
                thickness=20,
                tickmode='array',
                tickvals=[0.1, 0.3, 0.5, 0.7, 0.9],
                ticktext=['NA and 0–0.2', '>0.2–0.4', '>0.4–0.6', '>0.6–0.8', '>0.8–1.0'],
                tickfont=dict(size=14),
            ),
            text=hover_text,
            hovertemplate='%{text}<extra></extra>',
        ),
        row=1,
        col=1,
    )
    figure.add_trace(
        go.Heatmap(
            z=annotation,
            zmin=0,
            zmax=max(1, len(label_order) - 1),
            colorscale=annotation_colorscale,
            showscale=False,
            hoverinfo='skip',
        ),
        row=2,
        col=1,
    )
    for boundary in group_boundaries:
        figure.add_vline(x=boundary, line_width=1, line_color='white', row=1, col=1)
    figure.update_xaxes(
        tickmode='array',
        tickvals=group_centers,
        ticktext=[f'<b>{label}</b>' for label in label_order],
        tickfont=dict(size=20),
        ticks='',
        showgrid=False,
        zeroline=False,
        row=2,
        col=1,
    )
    figure.update_xaxes(showticklabels=False, showgrid=False, zeroline=False, row=1, col=1)
    figure.update_yaxes(
        showticklabels=False,
        showgrid=False,
        zeroline=False,
        scaleanchor='x',
        scaleratio=1,
        constrain='domain',
        row=1,
        col=1,
    )
    figure.update_yaxes(showticklabels=False, showgrid=False, zeroline=False, row=2, col=1)
    figure.update_layout(
        height=850,
        margin=dict(l=30, r=30, t=115, b=65),
        plot_bgcolor='white',
        paper_bgcolor='white',
        showlegend=False,
        hovermode='closest',
    )
    return figure


def build_haplotype_ld_roof(data_path, chrom, haplo_df, haplo_label_map, ancestry='European'):
    """Return the regional haplotype LD roof, or an explanatory Plotly panel."""
    try:
        matrix, ordered_variants = build_haplotype_ld_matrix(
            data_path,
            chrom,
            haplo_df,
            haplo_label_map,
            ancestry=ancestry,
        )
        return haplotype_ld_roof_plotly(matrix, ordered_variants)
    except Exception as exc:
        logging.exception("Unable to create regional haplotype LD plot")
        figure = go.Figure()
        figure.add_annotation(
            text=f"Regional LD plot unavailable: {exc}",
            x=0.5,
            y=0.5,
            xref='paper',
            yref='paper',
            showarrow=False,
            font=dict(size=16, color='#666666'),
        )
        figure.update_layout(
            height=260,
            margin=dict(l=20, r=20, t=20, b=20),
            plot_bgcolor='white',
            paper_bgcolor='white',
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
        )
        return figure

def summarize_region_haplotypes(haplo_df):
    """Build the complete ranked table used for a regional haplotype query."""
    columns = ['haplo_index', 'ID', 'n_snps', 'start_bp', 'end_bp', 'n_traits', 'traits']
    if haplo_df.empty:
        return pd.DataFrame(columns=columns)
    summary = (
        haplo_df.groupby('ID')
        .agg(
            n_snps=('NSNPS', 'max'),
            start_bp=('BP1', 'min'),
            end_bp=('BP2', 'max'),
            n_traits=('n_traits', 'max'),
            traits=('traits', 'max'),
        )
        .reset_index()
    )
    summary['traits'] = summary['traits'].apply(summarize_traits_with_counts)
    summary = summary.sort_values(
        by=['n_traits', 'start_bp'], ascending=[False, True]
    ).reset_index(drop=True)
    summary['haplo_index'] = summary.index.map(index_to_letters)
    return summary[columns]


# function to guide haplotype plots with snps/genes as inputs
def guide_haplotypes_snps_genes(data_path, browse, window, refGen, ancestry='European'):
    # Derive chromosome and position
    chrom, start_pos, end_pos, browse_type = readBrowseOption(data_path, browse, window, refGen)
    # get data of the region of interest
    haplo_df, _ = extract_haplo_data(
        data_path, chrom, start_pos, end_pos, ancestry=ancestry
    )
    # identify the haplotype when an individual SNP is searched
    hap_id_interest = []
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
    haplo_summary = summarize_region_haplotypes(haplo_df)
    # Start the regional LD panel with the six highest-ranked haplotypes. The
    # interactive selector can subsequently request any haplotypes in the
    # complete regional summary.
    ld_plot_df = haplo_summary.head(6).copy()
    ld_haplo_label_map = dict(zip(ld_plot_df['ID'], ld_plot_df['haplo_index']))
    # keep first matching haplotype for highlighting if available
    if hap_id_interest == []:
        hap_id_interest = 'Not found'
    else:
        hap_id_interest = hap_id_interest[0]  # take the first one for highlighting in the plot
    # Use the same layout for the initial single-ancestry view and comparisons.
    fig = haplotype_ancestry_overview(
        data_path, chrom, start_pos, end_pos, [ancestry], refGen=refGen,
        limit=10, preloaded={ancestry: haplo_df}, genes=genes,
    )
    plot_url = pio.to_html(
        fig, include_plotlyjs='cdn', full_html=False,
        config={'responsive': True, 'displaylogo': False},
    )
    ld_fig = build_haplotype_ld_roof(
        data_path,
        chrom,
        haplo_df[haplo_df['ID'].isin(ld_haplo_label_map)].copy(),
        ld_haplo_label_map,
        ancestry=ancestry,
    )
    ld_plot_url = pio.to_html(
        ld_fig,
        include_plotlyjs=False,
        full_html=False,
        config={'responsive': True, 'displaylogo': False},
    )
    return (
        browse,
        plot_url,
        ld_plot_url,
        haplo_summary,
        chrom,
        start_pos,
        end_pos,
        hap_id_interest,
    )


def haplotype_ancestry_overview(data_path, chrom, start_pos, end_pos, ancestries,
                               refGen='GRCh38', limit=10, preloaded=None, genes=None):
    """Build aligned ancestry panels for the same genomic interval."""
    if not ancestries or any(a not in ('European', 'African') for a in ancestries):
        raise ValueError('Choose European and/or African ancestry.')
    if limit not in (10, 20):
        raise ValueError('Overview limit must be 10 or 20.')

    panels = []
    for ancestry in ancestries:
        if preloaded is not None and ancestry in preloaded:
            haplo_df = preloaded[ancestry]
        else:
            haplo_df, _ = extract_haplo_data(
                data_path, chrom, start_pos, end_pos, ancestry=ancestry
            )
        if haplo_df.empty:
            panels.append((ancestry, haplo_df))
            continue
        ranked = (haplo_df.groupby('ID', as_index=False)
                  .agg(n_traits=('n_traits', 'max'), start_bp=('BP1', 'min'))
                  .sort_values(['n_traits', 'start_bp'], ascending=[False, True])
                  .head(limit))
        rank_by_id = dict(zip(ranked['ID'], range(len(ranked))))
        visible = haplo_df[haplo_df['ID'].isin(rank_by_id)].copy()
        visible = visible.drop_duplicates('ID')
        visible['plot_rank'] = visible['ID'].map(rank_by_id)
        visible = visible.sort_values(['plot_rank', 'BP1'])
        panels.append((ancestry, visible))

    if genes is None:
        genes = extract_genes(data_path, chrom, start_pos, end_pos, refGen)
    n_panels = len(panels)
    heights = [max(3, len(frame)) for _, frame in panels] + [9]
    figure = make_subplots(
        rows=n_panels + 1, cols=1, shared_xaxes=True,
        row_heights=heights, vertical_spacing=0.025,
    )
    trait_values = [int(value) for _, frame in panels for value in frame['n_traits']]
    min_traits = min(trait_values, default=0)
    max_traits = max(trait_values, default=1)
    if min_traits == max_traits:
        max_traits = min_traits + 1
    scale_shown = False

    for row_number, (ancestry, frame) in enumerate(panels, start=1):
        code = 'EUR' if ancestry == 'European' else 'AFR'
        ticks, labels = [], []
        for _, haplo in frame.iterrows():
            snps = str(haplo['SNPS']).split('|')
            positions = [int(snp.split(':')[1]) / 1_000_000 for snp in snps]
            rank = int(haplo['plot_rank'])
            label = f'{code} {index_to_letters(rank)}'
            ticks.append(rank)
            labels.append(label)
            count = int(haplo['n_traits'])
            figure.add_trace(go.Scatter(
                x=positions, y=[rank] * len(positions), mode='markers+lines',
                line=dict(color='#8c8c8c', width=1.5),
                marker=dict(size=7, color=[count] * len(positions),
                            colorscale='Spectral', cmin=min_traits, cmax=max_traits,
                            showscale=not scale_shown,
                            colorbar=dict(title='Traits', x=1.01, len=0.45),
                            line=dict(color='black', width=1)),
                name=label, showlegend=False,
                meta=dict(kind='haplotype', ancestry=ancestry, snps=snps,
                          haplotype_id=str(haplo['ID']), label=label),
                customdata=[[snp] for snp in snps],
                hovertemplate=(f'{label}<br>Haplotype: {haplo["ID"]}'
                               '<br>SNP: %{customdata[0]}'
                               f'<br>Associated traits: {count}<extra></extra>'),
            ), row=row_number, col=1)
            scale_shown = True
        figure.update_yaxes(
            tickmode='array', tickvals=ticks, ticktext=labels,
            range=[max(len(ticks) - 0.5, 0.5), -0.5],
            title_text=ancestry, tickfont=dict(size=11), row=row_number, col=1,
        )
        if not ticks:
            figure.add_annotation(
                x=(start_pos + end_pos) / 2_000_000, y=0,
                text=f'No {ancestry.lower()} haplotypes in this region',
                showarrow=False, row=row_number, col=1,
            )

    gene_row = n_panels + 1
    if genes:
        index_pos = 4 if refGen == 'GRCh37' else 5
        index_name = 0 if refGen == 'GRCh37' else -2
        index_strand = 3 if refGen == 'GRCh37' else 2
        for gene in genes:
            x0 = int(gene[index_pos]) / 1_000_000
            x1 = int(gene[index_pos + 1]) / 1_000_000
            lane = int(gene[-1])
            figure.add_trace(go.Scatter(
                x=[x0, x1], y=[lane, lane], mode='lines',
                line=dict(color='#666' if gene[index_strand] == '+' else '#222', width=5),
                text=[str(gene[index_name])] * 2, hovertemplate='%{text}<extra></extra>',
                showlegend=False,
            ), row=gene_row, col=1)
            figure.add_annotation(
                x=(x0 + x1) / 2, y=lane + 0.15,
                text=str(gene[index_name]), showarrow=False,
                font=dict(size=14), row=gene_row, col=1,
            )
        figure.update_yaxes(range=[-0.5, max(int(g[-1]) for g in genes) + 0.7],
                            showticklabels=False, title_text='Genes', row=gene_row, col=1)
    else:
        figure.add_annotation(
            x=(start_pos + end_pos) / 2_000_000, y=0,
            text='No genes to plot', showarrow=False, row=gene_row, col=1,
        )
        figure.update_yaxes(showticklabels=False, title_text='Genes', row=gene_row, col=1)

    figure.update_xaxes(range=[start_pos / 1_000_000, end_pos / 1_000_000],
                        showgrid=True, gridcolor='#e0e0e0')
    figure.update_xaxes(title_text=f'Chromosome {chrom} position (Mb)', row=gene_row, col=1)
    figure.update_layout(
        height=round((280 + 14 * limit) * (1.35 if n_panels == 2 else 1)),
        margin=dict(l=85, r=80, t=30, b=55),
        plot_bgcolor='white', paper_bgcolor='white',
        hovermode='closest', showlegend=False,
    )
    return figure


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
def fetch_by_position(chrom, pos38):
    """
    Query the SQLite database for a variant using chromosome and hg38 position.
    Returns a pandas DataFrame.
    """
    DB_FILE = "%s/databases/Genes/variant_info.db" % (data_path)
    conn = sqlite3.connect(str(DB_FILE))
    chrom = str(chrom).upper().replace("CHR", "")
    query = """
        SELECT *
        FROM variant_info
        WHERE Chromosome = ?
          AND Position_hg38 = ?
    """
    rows = pd.read_sql_query(query, conn, params=(str(chrom), int(pos38))).to_dict(orient="records")
    conn.close()
    if not rows:
        return None
    info = rows[0]
    try:
        info['AF'] = round(float(info.get('AF')), 2)
    except Exception:
        info['AF'] = None
    return {
        "query_type": "coord",
        "input_chr": chrom,
        "input_pos_hg38": pos38,
        "rsid": info['rsID'],
        "chr_hg38": chrom,
        "pos_hg38": pos38,
        "pos_hg19": info['Position_hg19'],
        "ref": info['REF'],
        "alt": info['ALT'],
        "maf": info['AF'],
    }

# function to get LD between SNPs
def get_ld_between_snps(snps, data_path, chrom, ancestry='European'):
    snp_list = snps.split(' ')
    snp_positions = [int(snp.split(':')[1].split('-')[0]) for snp in snp_list]
    chr = snp_list[0].split(':')[0].upper().replace("CHR", "")
    chrom = chrom.upper().replace("CHR", "")
    # get rsids for the SNPs
    rsids_list = {}
    for snp in snp_positions:
        rsid_info = fetch_by_position('chr' + chr, snp)
        rsids_list[snp] = rsid_info['rsid'] if rsid_info else 'NA'
    # convert to dataframe
    rsids_df = pd.DataFrame.from_dict(rsids_list, orient='index', columns=['rsid']).reset_index().rename(columns={'index': 'position'})
    # define interval
    min_pos = min(snp_positions)
    max_pos = max(snp_positions)
    # run function to get LD by interval
    db_path = ancestry_resource_path(
        data_path,
        'LD_db',
        ancestry,
        f"ld_chr{chrom}.sqlite",
    )
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
        elif n_haplotypes <= 40:
            px = 1200
        else:
            px = 1400
    else:
        px = 400
    return px

# function to plot haplotype regional plot
def haplotype_plotly_regional(haplo_df, haplo_dict, genes, refGen='GRCh38', haplo_label_map=None):
    # Keep only haplotypes in the provided ranking map (already top-50 by #traits)
    # and force the plot order to follow the table order (A, B, C, ...).
    if haplo_label_map:
        rank_map = {hid: i for i, hid in enumerate(haplo_label_map.keys())}
        haplo_df = haplo_df[haplo_df['ID'].isin(rank_map)].copy()
        haplo_df['plot_rank'] = haplo_df['ID'].map(rank_map)
        haplo_df = haplo_df.sort_values(by=['plot_rank', 'BP1'], ascending=[True, True]).reset_index(drop=True)
    else:
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
    y_labels = []
    for index, row in haplo_df.iterrows():
        haplo_id = row['ID']
        snps = row['SNPS'].split('|')
        snps_positions = [int(x.split(':')[1]) / 1_000_000 for x in snps]
        n_traits = int(row['n_traits'])
        n_traits_list = [n_traits for _ in snps]
        trait_list = row['traits'].split('___') if pd.notna(row['traits']) else []
        y_val = row['plot_rank'] if 'plot_rank' in row.index else index
        index_list = [y_val for _ in snps]
        y_labels.append((y_val, haplo_label_map.get(haplo_id, str(y_val)) if haplo_label_map else str(y_val)))
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
    if y_labels:
        y_ticks = [x[0] for x in y_labels]
        y_ticktext = [x[1] for x in y_labels]
    else:
        y_ticks = [0]
        y_ticktext = [""]
    fig.update_yaxes(
        range=[-0.2, max_idx + 0.2],
        tickmode='array',
        tickvals=y_ticks,
        ticktext=y_ticktext,
        autorange='reversed',
        title_text="Haplotype",
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
def extract_haplo_data(data_path, chrom, start_pos, end_pos, ancestry='European'):
    # define data to look at
    haplotype_filename = (
        '20260917_ld_clusters_with_gwas_AI_latest_African.tsv.gz'
        if str(ancestry).strip().upper() in {'AFR', 'AFRICAN'}
        else 'ld_clusters_with_gwas_AI_latest.tsv.gz'
    )
    haplo_path = f"{data_path}/databases/haplotypes/{haplotype_filename}"
    # tabix query
    chrom = str(chrom).upper().replace('CHR', '')
    tabix_result = subprocess.run(
        ['tabix', haplo_path, f'{chrom}:{start_pos}-{end_pos}'],
        capture_output=True,
        check=True,
        text=True,
    )
    data_list = tabix_result.stdout.splitlines()
    # convert list to df
    columns = ['CHR', 'BP1', 'BP2', 'KB', 'NSNPS', 'SNPS', 'assoc_count', 'rsids', 'traits']
    if not data_list:
        return pd.DataFrame(columns=columns + ['n_traits', 'ID']), {}
    haplo_df = pd.read_csv(StringIO('\n'.join(data_list)), sep='\t', header=None, names=columns)
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
def extract_ld(data_path, chrom, df, refGen, browse_type, browse, ancestry='European'):
    try:
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
        db_path = ancestry_resource_path(
            data_path,
            'LD_db',
            ancestry,
            f"ld_chr{str(chrom).replace('chr', '')}.sqlite",
        )
        # query the database
        res1 = partners_by_uniq(str(db_path), uniq1, r2_min=0.2, limit=100)
        res2 = partners_by_uniq(str(db_path), uniq2, r2_min=0.2, limit=100)
        # combine results
        res = res1 + res2
        # convert to df
        ld_df = pd.DataFrame(res, columns=['partner_uniq', 'r2', 'dist_bp'])
        # add the original variant
        ld_df = pd.concat([ld_df, pd.DataFrame({'partner_uniq': [uniq1], 'r2': [1.0], 'dist_bp': [0]})], ignore_index=True)
        # extract position from partner_uniq
        ld_df['pos'] = ld_df['partner_uniq'].apply(lambda x: int(float(str(x).split(':')[0])))
        # assign colors for ld: >=0.8 red, >=0.6 orange, >=0.4 green, >=0.2 lightblue, <0.2 grey
        ld_df['color'] = ld_df['r2'].apply(assign_color)
    except Exception:
        ld_df = pd.DataFrame(columns=['partner_uniq', 'r2', 'dist_bp', 'pos', 'color'])
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
      WHERE (ld.v1 = v.id OR ld.v2 = v.id)
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
def extract_haplo(data_path, df, chrom, refGen, ancestry='European'):
    # take min and max position wrt hg38
    if refGen == 'GRCh37':
        min_pos = df['Position_hg38'].min()
        max_pos = df['Position_hg38'].max()
    else:
        min_pos = df['Position'].min()
        max_pos = df['Position'].max()
    # take positions of interest
    is_african = str(ancestry).strip().upper() in {'AFR', 'AFRICAN'}
    haplotype_filename = (
        '20260914_thr0.2_haplotypes.tsv.gz'
        if is_african
        else '20251122_ld_clusters_low_recombination_intervals.tsv.gz'
    )
    haplotype_path = ancestry_resource_path(
        data_path,
        'haplotypes',
        ancestry,
        haplotype_filename,
    )
    region = f"{str(chrom).replace('chr', '')}:{int(min_pos - 500000)}-{int(max_pos + 500000)}"
    tabix_result = subprocess.run(
        ['tabix', str(haplotype_path), region],
        capture_output=True,
        check=True,
        text=True,
    )
    haplo = [line.split('\t') for line in tabix_result.stdout.splitlines() if line]
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
        if os.path.exists(path_server + '/Annotation/RUNS/snpXplorer_annotation_%s.zip' %(run_id)):
            messageError = False
            messageToUser = 'Valid Run ID. Download will start soon'
        else:
            messageError = True
            messageToUser = 'Not valid Run ID. Try again.'
    return messageError, messageToUser

# function to fetch by rsid
def fetch_by_rsid(rsid):
    """
    Query the SQLite database for a single rsID.
    Returns a pandas DataFrame (can have 0, 1, or multiple rows).
    """
    DB_FILE = "%s/databases/Genes/variant_info.db" % (data_path)
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

# function to extract coordinates based on input type
def readBrowseOption(data_path, browse, window, refGen):
    # check if input is in the form 1:100000
    if ':' in browse and '-' not in browse:
        chrom, start_pos, end_pos, browse_type = int(browse.split(':')[0]), int(browse.split(':')[1]) - window, int(browse.split(':')[1]) + window, 'Single position'
    elif ':' in browse and '-' in browse:
        chrom, start_pos, end_pos, browse_type = int(browse.split(':')[0]), int(browse.split(':')[1].split('-')[0]) - window, int(browse.split(':')[1].split('-')[1]) + window, 'Interval'
    elif browse.lower().startswith("rs"):   # user picked an rsID
        tmp_info = fetch_by_rsid(browse.replace(' ', '').split()[0])
        chrom, pos_hg38, pos_hg19 = tmp_info['Chromosome'], tmp_info['Position_hg38'], tmp_info['Position_hg19']
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
def run_annotation(
    my_list,
    inpType,
    refGeno,
    analType,
    gsea_source,
    qtl_tissues,
    email,
    ancestry='European',
):
    # Sample a number for randomization
    random_number = random.randint(1, 100000)
    # Take snps and save them
    fpath_server = path_server + '/Annotation/RUNS/snpXplorer_input_%s.txt' %(random_number)
    filename = 'snpXplorer_input_%s.txt' %(random_number)
    with open(fpath_server, 'w') as finp:
        for x in my_list:
            finp.write('%s\n' %(x))
    # Define log file
    #log_filename_server = '/Annotation/RUNS/snpXplorer_run_%s.log' %(random_number)
    # Take analysis type
    analysis_type = str(analType)
    analysis_type = 'annotation' if analysis_type == 'Annot' else 'enrichment'
    analysis_mode = ','.join(list(gsea_source)) if analysis_type == 'enrichment' else 'None'
    analysis_mode = analysis_mode.replace('Default', 'GO:BP')
    # Take gtex tissues
    gtex_tissues = ','.join(list(qtl_tissues))
    refGeno = str(refGeno).lower()
    ancestry = 'African' if str(ancestry).strip().lower() == 'african' else 'European'
    # Then run annotate me externally in background -- this depends on the analysis_type requested
    #command = "Rscript /Annotation/BIN/MAIN.R %s %s %s %s %s %s %s %s > %s" %(filename, ftype, email, analysis_type, analysis_mode, gtex_tissues, refGeno, random_number, log_filename_server)
    command = f"python {path_server}/Annotation/BIN/standalone_annotation_optimized_mp.py -q {fpath_server} -t {analysis_type} -o {path_server}/Annotation/RUNS/ -r {random_number} -b {refGeno} -a {ancestry} -ts {gtex_tissues} -gs {analysis_mode} -e {email} --workers 2 --gsea-workers 2 --gnomad-source hail "
    print(command)
    return command, str(random_number)

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
def get_data_plot(data_path, gwas, chrom, start_pos, end_pos, refGen, meta, window):
    # define region to look at
    if refGen == 'GRCh37':
        # set liftover
        liftover_info = get_lifter('hg19', 'hg38')
        region = '%s:%s-%s' %(str(chrom).replace('chr', ''), str(start_pos), str(end_pos))
        try:
            region_hg38 = '%s:%s-%s' %(str(chrom).replace('chr', ''), str(liftover_info[chrom][start_pos][0][1]), str(liftover_info[chrom][end_pos][0][1]))
        except:
            region_hg38 = '%s:%s-%s' %(str(chrom).replace('chr', ''), str(liftover_info[chrom][start_pos][0][1]), str(liftover_info[chrom][start_pos][0][1] + window))
    else:
        # set liftover
        liftover_info = get_lifter('hg38', 'hg19')
        region_hg38 = '%s:%s-%s' %(str(chrom), str(start_pos), str(end_pos))
        try:
            region = '%s:%s-%s' %(str(chrom), str(liftover_info[chrom][start_pos][0][1]), str(liftover_info[chrom][end_pos][0][1]))
        except:
            region = '%s:%s-%s' %(str(chrom), str(liftover_info[chrom][start_pos][0][1]), str(liftover_info[chrom][start_pos][0][1] + window))
        region = region.upper().replace('CHR', '')
    # get gwas info
    data_list = pd.DataFrame()
    data_list_info = []
    sumstats_build = ''
    # iterate over gwas
    for gw in gwas:
        # extract gwas info from label -- label is "Alzheimer's disease — Bellenguez C (2022) [ebi-a-GCST90027158]" --> extract ebi-a-GCST90027158
        gwas_id = gwas_id_from_selection(gw)
        # get gwas info
        gwas_info = meta[gwas_id]
        ancestry_code = gwas_ancestry_code(meta, gwas_id)
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
        if 'Gwas' in sumstats.columns:
            sumstats['Gwas'] = sumstats['Gwas'].apply(
                lambda trait: f"[{ancestry_code}] {trait}"
            )
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
        else:
            svtypes_fixed.append(sv)
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
def scatterplot_plotly(df, chrom, start_pos, end_pos, gwas, genes, svs, browse_type, refGen, recomb_data, exons, browse, plotype, ld, ld_df, yrange, haplo_df, plotHaplo, axestype):
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
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.02, row_heights=[0.45, 0.30, 0.25], specs=[[{"secondary_y": True}], [{}], [{}]])

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

    # Define y-axis variable
    yaxis_interest = 'Pvalue' if axestype == 'manh' else 'Beta'
    yaxis_label = '-log10(P-value)' if axestype == 'manh' else 'Beta'

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
                    go.Scatter(
                        x=sig['Position_plot'], 
                        y=sig[yaxis_interest],
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
                    go.Scatter(
                        x=nonsig['Position_plot'], y=nonsig[yaxis_interest],
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
                        go.Scatter(
                            x=sig['Position_plot'], y=sig[yaxis_interest],
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
                        go.Scatter(
                            x=nonsig['Position_plot'], y=nonsig[yaxis_interest],
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
            tmp['smoothed'] = tmp[yaxis_interest].rolling(window=20).max()
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
        fig.update_yaxes(range=[min_y, max_y], row=1, col=1, title_text=yaxis_label, title_font=dict(size=18, family="Arial", color="black", weight="bold"), tickfont=dict(size=14, family="Arial", color="black"))
    else:
        if (axestype == 'manh'):
            max_p = float(np.nanmax(df['Pvalue'])) if len(df) else 6.0
            max_p_assoc = math.ceil(max_p) + math.ceil(max_p) * 0.10
            fig.update_yaxes(range=[0, max(max_p_assoc, 6)], row=1, col=1, title_text=yaxis_label, title_font=dict(size=18, family="Arial", color="black", weight="bold"), tickfont=dict(size=14, family="Arial", color="black"))
        else:
            max_p = float(np.nanmax(df['Beta'])) if len(df) else 1
            max_p_assoc = math.ceil(max_p) + math.ceil(max_p) * 0.10
            fig.update_yaxes(range=[0, max(max_p_assoc, 1)], row=1, col=1, title_text=yaxis_label, title_font=dict(size=18, family="Arial", color="black", weight="bold"), tickfont=dict(size=14, family="Arial", color="black"))

    # rsID highlight
    if browse_type == 'RsID':
        subset_rsid = df[df['Rsid'].str.lower() == str(browse).lower()].copy()
        if not subset_rsid.empty:
            colors = normalize_color(list(subset_rsid['Gwas'].map(color_map)))
            fig.add_trace(
                go.Scatter(
                    x=subset_rsid['Position_plot'],
                    y=subset_rsid[yaxis_interest],
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
        top_snp = df.loc[top_idx, ['Position_plot', 'Pvalue', 'Gwas', 'Rsid', 'Beta']].copy()
        if 'GWAS-Cat' in gwas:
            top_snp_noCat = top_snp[top_snp['Gwas'].isin(gwas)]
            top_snp_Cat = top_snp.sort_values('Pvalue', ascending=False).head(3)
            top_snp = pd.concat([top_snp_Cat, top_snp_noCat], axis=0)

        for _, row in top_snp.iterrows():
            if browse_type == 'RsID' and str(row['Rsid']).lower() == str(browse).lower():
                continue
            fig.add_annotation(
                x=row['Position_plot'], y=row[yaxis_interest],
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
        if len(genes) < 8:
            fsize, yoff, gw = 16, max_y * 0.10, 6
        elif len(genes) < 10:
            fsize, yoff, gw = 14, max_y * 0.08, 5
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

# function to draw beta-maf
def beta_maf_plot(df, xcol, ycol, title=""):
    # Define marker types based on Gwas
    markers = ['circle', 'square', 'diamond', 'cross', 'x', 'triangle-up', 'triangle-down', 'pentagon', 'hexagon', 'star']
    gwas_unique = df['Gwas'].unique().tolist()
    # limit number of GWAS to number of markers
    max_gwas = len(markers)
    if len(gwas_unique) > max_gwas:
        gwas_unique = gwas_unique[:max_gwas]
    # map each GWAS to a marker
    marker_map = {gw: markers[i] for i, gw in enumerate(gwas_unique)}
    df['marker'] = df['Gwas'].map(marker_map).fillna(markers[0])
    # sort by Pvalue descending so that most significant are plotted last
    df = df.sort_values(by='Pvalue', ascending=True)
    # create plotly figure
    fig = go.Figure()
    fig.update_xaxes(showgrid=True, gridcolor="#e0e0e0", zeroline=True)
    fig.update_yaxes(showgrid=True, gridcolor="#e0e0e0", zeroline=False)
    fig.add_trace(go.Scatter(
        x=df[xcol],
        y=df[ycol],
        mode="markers",
        marker=dict(size=11, opacity=0.8, color=df['Pvalue'], symbol=df['marker'], colorscale='Viridis', colorbar=dict(title='-log10(p)'), showscale=True),
        hovertemplate=("<b>%{customdata[0]}</b><br>"
                    "Gwas: %{customdata[1]}<br>"
                    "Pos: %{customdata[2]}<br>"
                    "-log10(p): %{customdata[3]}<br>"
                    "Beta: %{x:.2f}<br>"
                    "EA: %{customdata[4]}<br>"
                    "NEA: %{customdata[5]}<br>"
                    "EAF: %{y:.2f}<extra></extra>"
                    ),
                        customdata=np.stack([
                            df['Rsid'], 
                            df['Gwas'], 
                            df['Position'],
                            df['Pvalue'], 
                            df['EA'], 
                            df['NEA']
                        ], axis=1)
    ))
    fig.update_layout(legend=dict(x=1.02, y=1, xanchor="left", yanchor="top"))
    fig.update_yaxes(title_text="Frequency (Effect allele)", title_font=dict(size=18, family="Arial", color="black", weight="bold"), tickfont=dict(size=14, family="Arial", color="white"))
    fig.update_xaxes(title_text="Beta (Effect allele)", title_font=dict(size=18, family="Arial", color="black", weight="bold"), tickfont=dict(size=14, family="Arial", color="black"))
    fig.update_layout(plot_bgcolor="white", paper_bgcolor="white")  
    return fig

# function to check if an email is valid
def is_valid_email(email):
    pattern = r'^[\w\.-]+@[\w\.-]+\.\w+$'
    return re.match(pattern, email) is not None

# function to monitor exploration jobs
def add_search_to_file(data_path, chrom, start_pos, end_pos, refGen, gwas, browse_type):
    log_file = f"{data_path}/monitor/exploration_logs.txt"
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    gwas_str = ";".join(gwas) if isinstance(gwas, list) else str(gwas)
    log_entry = f"{timestamp}\t{chrom}\t{start_pos}\t{end_pos}\t{refGen}\t{gwas_str}\t{browse_type}\n"
    with open(log_file, "a") as f:
        f.write(log_entry)

# function to monitor exploration jobs
def add_haplo_to_file(data_path, browse, refGen, browse_type):
    log_file = f"{data_path}/monitor/haplotypes_logs.txt"
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_entry = f"{timestamp}\t{browse}\t{browse_type}\t{refGen}\n"
    with open(log_file, "a") as f:
        f.write(log_entry)

# function to get top n GWAS rows
def top_gwas_rows(df: pd.DataFrame, n: int = 100) -> pd.DataFrame:
    """Return top n most significant rows from GWAS assoc table."""
    if df is None or df.empty:
        return df
    out = df.copy()
    # Prefer raw p if present
    if "p" in out.columns:
        out["p"] = pd.to_numeric(out["p"], errors="coerce")
        out = out.sort_values("p", ascending=True, na_position="last")
        return out.head(n)
    # Else prefer original -log10(p) if present (bigger = more significant)
    if "orig_neg_log10_p" in out.columns:
        out["orig_neg_log10_p"] = pd.to_numeric(out["orig_neg_log10_p"], errors="coerce")
        out = out.sort_values("orig_neg_log10_p", ascending=False, na_position="last")
        return out.head(n)
    # Else fall back to neg_log10_p (depending on your scaling it might be negative)
    if "neg_log10_p" in out.columns:
        out["neg_log10_p"] = pd.to_numeric(out["neg_log10_p"], errors="coerce")
        # If your scaled values are negative, "more significant" may be more negative.
        # Safer: sort by absolute magnitude descending.
        out = out.sort_values("neg_log10_p", ascending=True, na_position="last")
        return out.head(n)
    # Otherwise just take first n
    return out.head(n)
