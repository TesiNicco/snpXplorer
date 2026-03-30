# AI-based way to do semantic similarity

# # If needed, install deps:
# conda activate nlpy

from sentence_transformers import SentenceTransformer
from sklearn.preprocessing import normalize
from sklearn.metrics.pairwise import cosine_distances
from sklearn.cluster import AgglomerativeClustering
from collections import defaultdict
import numpy as np
import re
import pandas as pd
import os
import pickle
import umap.umap_ as umap
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# 1) Preprocessing (lightweight; keep meaning, remove boilerplate)
# ------------------------------------------------------------
STOPPHRASES = [r'\bgwas (of|on|for)\b', r'\bgenome[-\s]?wide association (study|analysis)\b', r'\bmeta[-\s]?analysis\b', r'\bstudy of\b', r'\bassociation with\b']
STOPWORDS = set(['the','and','of','for','to','with','in','on','a','an','by','from'])

def normalize_trait(t: str) -> str:
    s = t.lower()
    s = re.sub(r'[\(\)\[\]\{\},:;]', ' ', s)
    for pat in STOPPHRASES:
        s = re.sub(pat, ' ', s)
    s = re.sub(r'\s+', ' ', s).strip()
    # drop leading/trailing stopwords; keep acronyms (e.g., CAD, LDL)
    tokens = [w for w in s.split() if w not in STOPWORDS]
    return ' '.join(tokens) if tokens else t.lower().strip()

# ------------------------------------------------------------
# 2) Load models (general + biomedical)
#    You can swap models or add more (and adjust weights)
# ------------------------------------------------------------
MODEL_NAMES = [
    'sentence-transformers/all-MiniLM-L6-v2',
    'amrothemich/sapbert-sentence-transformers',  # biomedical concept similarity
    # or: 'pritamdeka/SapBERT-mnli-snli-scinli-scitail-mednli-stsb'
]
WEIGHTS = np.array([0.5, 0.5], dtype=np.float32)                # tweakable: e.g., [0.3, 0.7]

models = [SentenceTransformer(m) for m in MODEL_NAMES]

# ------------------------------------------------------------
# 3) Encode + Ensemble (weighted, L2-normalized)
# ------------------------------------------------------------
def encode_ensemble(texts, models, weights, batch_size=256, show_progress_bar=True):
    all_embs = []
    for model in models:
        embs = model.encode(
            texts,
            batch_size=batch_size,
            show_progress_bar=show_progress_bar,
            convert_to_numpy=True,
            normalize_embeddings=True,   # ensures cosine sim = dot product
        )
        all_embs.append(embs.astype(np.float32))
    # stack: [n_models, n_items, dim] -> weighted sum across models
    max_dim = max(e.shape[1] for e in all_embs)
    # pad if dims differ (rare, but MiniLM=384, SapBERT=768) -> project to common space via simple zero-pad
    padded = []
    for e in all_embs:
        if e.shape[1] < max_dim:
            pad = np.zeros((e.shape[0], max_dim - e.shape[1]), dtype=np.float32)
            e = np.hstack([e, pad])
        padded.append(e)
    weights = weights / weights.sum()
    ens = np.tensordot(weights, np.stack(padded, axis=0), axes=(0,0))  # [n_items, max_dim]
    ens = normalize(ens)  # final L2 norm
    return ens

# ------------------------------------------------------------
# 4) Clustering (cosine distance, hierarchical; threshold is key)
# ------------------------------------------------------------
def cluster_embeddings(embeddings, distance_threshold=0.20):
    # cosine distance in [0,2]; similar items -> small distance
    dist = cosine_distances(embeddings)
    # Agglomerative on a precomputed distance matrix. 'average' handles cosine well.
    clustering = AgglomerativeClustering(
        n_clusters=None,
        metric='precomputed',
        linkage='average',
        distance_threshold=distance_threshold
    )
    labels = clustering.fit_predict(dist)
    return labels, dist

# ------------------------------------------------------------
# 5) Representative selection per cluster
#    Option A: pick closest to cluster centroid in ensemble space
#    Option B: shortest normalized string (often nice for labels)
# ------------------------------------------------------------
def pick_representatives(names, embeddings, labels, method='centroid'):
    reps = {}
    for lab in np.unique(labels):
        idx = np.where(labels == lab)[0]
        if method == 'shortest':
            # Use normalized strings length as proxy
            sub = [(i, len(normalize_trait(names[i]))) for i in idx]
            rep = min(sub, key=lambda x: x[1])[0]
        else:
            # centroid: pick item with max cosine similarity to mean
            centroid = normalize(np.mean(embeddings[idx], axis=0).reshape(1, -1))[0]
            sims = embeddings[idx] @ centroid
            rep = idx[int(np.argmax(sims))]
        reps[lab] = rep
    return reps

# ------------------------------------------------------------
# 6) Driver: de-duplicate trait names
# ------------------------------------------------------------
def deduplicate_traits(trait_names,
                       weights=WEIGHTS,
                       distance_threshold=0.20,
                       pick='centroid'):
    # keep original + normalized for reporting
    normed = [normalize_trait(t) for t in trait_names]
    ens = encode_ensemble(normed, models, weights, show_progress_bar=True)
    labels, dist = cluster_embeddings(ens, distance_threshold=distance_threshold)
    reps = pick_representatives(trait_names, ens, labels, method=pick)
    # build outputs: mapping and reduced list
    clusters = defaultdict(list)
    for i, lab in enumerate(labels):
        clusters[lab].append(i)
    representative_names = [trait_names[reps[lab]] for lab in sorted(clusters)]
    mapping = {
        int(lab): {
            'representative': trait_names[reps[lab]],
            'members': [trait_names[i] for i in clusters[lab]]
        } for lab in clusters
    }
    return representative_names, mapping, labels, reps, ens, dist

# ------------------------------------------------------------
# Helper to clean GO definitions
# ------------------------------------------------------------
def clean_go_description(text: str) -> str:
    GO_PREFIX_PATTERNS = [r'^"?', r'^any process that modulates the (frequency|rate|extent) of ', r'^any process that modulates the (rate|frequency),? or extent of ', r'^any process that modulates ', r'^the (series of )?molecular signals initiated by ', r'^the chemical reactions and pathways resulting in the formation of ', r'^the chemical reactions and pathways involving ', r'^a process that modulates ']
    # Patterns for trailing noise
    GO_SUFFIX_PATTERNS = [r'\[[^\]]*\]', r'"$']
    s = text.strip()
    # 1. Remove suffixes first (metadata blocks)
    for pat in GO_SUFFIX_PATTERNS:
        s = re.sub(pat, '', s).strip()
    # 2. Lowercase for normalization operations
    s = s.lower()
    # 3. Remove prefixes and boilerplate phrases
    for pat in GO_PREFIX_PATTERNS:
        s = re.sub(pat, '', s)
    # 4. Remove repeated boilerplate if present twice
    for pat in GO_PREFIX_PATTERNS:
        s = re.sub(pat, '', s)
    # 5. Clean punctuation + normalize spaces
    s = re.sub(r'[",:;{}\[\]\(\)]', ' ', s)
    s = re.sub(r'\s+', ' ', s).strip()
    return s

def clean_do_names(term: str) -> str:
    REGULATION_PATTERNS = [r'positive regulation of', r'negative regulation of', r'regulation of', r'modulation of', r'process involved in', r'involved in', r'response to', r'mediated by', r'dependent on']
    REPLACE_PATTERNS = {r'biosynthetic process': 'biosynthesis', r'catabolic process': 'catabolism', r'metabolic process': 'metabolism', r'transport process': 'transport'}
    s = term.lower()
    # Remove regulation prefixes
    for pat in REGULATION_PATTERNS:
        s = re.sub(pat, '', s)
    # Replace long biology phrases
    for pat, repl in REPLACE_PATTERNS.items():
        s = re.sub(pat, repl, s)
    # Cleanup leftover "of", "by", "to"
    s = re.sub(r'\b(of|by|to|the|a)\b', ' ', s)
    # Normalize punctuation
    s = re.sub(r'[,:\(\)]', ' ', s)
    # Collapse spaces
    s = re.sub(r'\s+', ' ', s).strip()
    return s
# ------------------------------------------------------------

# Main execution
# ------------------------------------------------------------
# Get the trait names - this comes from the other script
inp_file = pd.read_csv('enrichment_terms_temp.txt', sep="\t")
# Description-based
#trait_names = inp_file['description'].tolist()
#trait_names = [clean_go_description(t) for t in trait_names]
# Name-based
trait_names = inp_file['name'].tolist()
trait_names = [clean_do_names(t) for t in trait_names]

# Run deduplication
# distance_treshold can be tuned; lower = more clusters, higher = fewer clusters [0.3-0.4 is aggressive]
reps_list = []
cluster_map_list = []
labels_list = []
rawreps_list = []
ens_list = []
dist_list = []
for thr in [0.2, 0.4, 0.6, 0.7, 0.75, 0.8]:
    print(f"Processing threshold: {thr}")
    reps, cluster_map, labels, reps_raw, ens, dist = deduplicate_traits(trait_names, weights=np.array([0.4, 0.6], dtype=np.float32), distance_threshold=thr, pick='centroid')
    reps_list.append(reps)
    cluster_map_list.append(cluster_map)
    labels_list.append(labels)
    rawreps_list.append(reps_raw)
    ens_list.append(ens)
    dist_list.append(dist)

for i, thr in enumerate([0.2, 0.4, 0.6, 0.7, 0.75, 0.8]):
    print(f"Threshold: {thr}")
    reps = reps_list[i]
    cluster_map = cluster_map_list[i]
    labels = labels_list[i]
    print(f"Original: {len(trait_names)}  ->  Unique reps: {len(reps)}")

# Convert cluster_map to a DataFrame for easier viewing
all_dfs = pd.DataFrame()
for i, thr in enumerate([0.6, 0.7, 0.75, 0.8]):
    cluster_map = cluster_map_list[i]
    rows = []
    for lab, info in cluster_map.items():
        rep = info['representative']
        members = "___".join(info['members'])
        n_members = len(info['members'])
        rows.append({'Cluster Label': lab, 'Representative': rep, 'Members': members, 'Number of Members': n_members})
    df = pd.DataFrame(rows)
    df['Threshold'] = thr
    all_dfs = pd.concat([all_dfs, df], ignore_index=True)

# We need to pick the representative -- basically the study with most associations
# Reset index
#all_dfs = all_dfs.reset_index(drop=True)
# Iterate over rows of the df
#all_dfs['Final Representative'] = ''
#all_dfs['Final Representative associations'] = 0
#all_dfs['Final Representative file'] = ''
#for idx, row in all_dfs.iterrows():
    print(f"Processing cluster {idx+1} out of {len(all_dfs)}        ", end='\r')
    # get all members
    all_members = row['Members'].split('___')
    # get traits info
    traits_info = assoc_counts_df[assoc_counts_df['Trait Name'].isin(all_members)]
    if traits_info.empty:
        all_dfs.loc[idx, 'Final Representative'] = 0
        all_dfs.loc[idx, 'Final Representative associations'] = 0
        all_dfs.loc[idx, 'Final Representative file'] = ''
        continue
    # sort by association count
    traits_info = traits_info.sort_values(by='Association Count', ascending=False)
    # pick the top one
    all_dfs.loc[idx, 'Final Representative'] = traits_info.iloc[0]['Trait Name']
    all_dfs.loc[idx, 'Final Representative associations'] = traits_info.iloc[0]['Association Count']
    all_dfs.loc[idx, 'Final Representative file'] = traits_info.iloc[0]['File Name']
# Save final summary
#all_dfs.to_csv("trait_clusters_summary_AI_thresholds.csv", index=False)

# ------------------------------------------------------------
# ------------------------------------------------------------
# Visualize embeddings with UMAP (optional)
# threshold 0.75 chosen arbitrarily
ens = ens_list[2]  # pick one threshold's embeddings
labels = labels_list[2]
clusters = cluster_map_list[2]
u = umap.UMAP(n_neighbors=50, min_dist=0.1, metric='cosine', random_state=42).fit_transform(ens)
# Plot
plt.figure(figsize=(12, 10))
scatter = plt.scatter(u[:, 0], u[:, 1], c=labels, cmap='tab20', s=80, alpha=0.95, edgecolors='k', linewidths=0.25, marker='o')
plt.title("UMAP of Trait Embeddings (semantic clusters)")
# Build legend from cluster dictionary 'clusters' (cluster_map_list[2])
uniq_labels = np.unique(labels)
cmap = plt.get_cmap('tab20')
norm = plt.Normalize(vmin=labels.min(), vmax=labels.max())
handles = []
legend_texts = []
for lab in sorted(uniq_labels):
    color = cmap(norm(lab))
    rep = clusters.get(int(lab), {}).get('representative', str(int(lab)))
    h = plt.Line2D([0], [0],
                   marker='o',
                   color='w',
                   markerfacecolor=color,
                   markeredgecolor='k',
                   markersize=8,
                   linewidth=0)
    handles.append(h)
    legend_texts.append(f"{int(lab)}: {rep}")

plt.legend(handles, legend_texts, title='Cluster: representative', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xlabel("UMAP 1")
plt.ylabel("UMAP 2")
plt.tight_layout()
plt.show()
plt.savefig("umap_trait_embeddings.png", dpi=300, bbox_inches='tight')

# save umap and embeddings
with open("trait_embeddings_umap_thr05.pkl", "wb") as f:
    pickle.dump({"umap": u, "embeddings": ens, "labels": labels, 'cluster_map': cluster_map}, f)
