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

# Function to plot cluster counts and cluster-size distributions by threshold. It takes a DataFrame containing cluster information and generates two plots: one for the number of clusters across different thresholds and another for the distribution of cluster sizes.
def plot_cluster_threshold_summary(all_dfs):
    """Plot cluster counts and cluster-size distributions by threshold."""
    required_columns = {'Threshold', 'Number of Members'}
    missing_columns = required_columns.difference(all_dfs.columns)
    if missing_columns:
        raise ValueError(
            f"all_dfs is missing required columns: {sorted(missing_columns)}"
        )
    plot_df = all_dfs.dropna(subset=list(required_columns)).copy()
    if plot_df.empty:
        raise ValueError("all_dfs has no non-missing values to plot")
    thresholds = sorted(plot_df['Threshold'].unique())
    positions = np.arange(len(thresholds))
    threshold_labels = [f'{threshold:g}' for threshold in thresholds]
    cluster_counts = plot_df.groupby('Threshold').size().reindex(thresholds)
    member_counts = [
        plot_df.loc[plot_df['Threshold'] == threshold, 'Number of Members'].values
        for threshold in thresholds
    ]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].bar(positions, cluster_counts.values, color='#4C78A8')
    axes[0].set(
        xlabel='Threshold',
        ylabel='Number of clusters',
        title='Number of clusters across thresholds',
        xticks=positions,
        xticklabels=threshold_labels,
    )
    violins = axes[1].violinplot(
        member_counts,
        positions=positions,
        widths=0.8,
        showmeans=False,
        showmedians=False,
        showextrema=False,
    )
    for body in violins['bodies']:
        body.set_facecolor('#72B7B2')
        body.set_edgecolor('#2F4F4F')
        body.set_alpha(0.8)
    axes[1].boxplot(
        member_counts,
        positions=positions,
        widths=0.18,
        patch_artist=True,
        showfliers=False,
        boxprops={'facecolor': 'white', 'edgecolor': '#333333'},
        medianprops={'color': '#E45756', 'linewidth': 1.5},
        whiskerprops={'color': '#333333'},
        capprops={'color': '#333333'},
    )
    axes[1].set(
        xlabel='Threshold',
        ylabel='Number of members (log scale)',
        title='Cluster size across thresholds',
        xticks=positions,
        xticklabels=threshold_labels,
    )
    axes[1].set_yscale('log')
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    return fig, axes

# ------------------------------------------------------------
ancestry = "African"
# 0) Get the trait names - takes a little bit of time
if ancestry == "European":
    data_path = '/project/holstegelab/Share/nicco/haplotypeBlocks/R2cluster/OpenGWAS_sumstats/European'
elif ancestry == "African":
    data_path = '/project/holstegelab/Share/nicco/haplotypeBlocks/R2cluster/OpenGWAS_sumstats/African'
# List files in the directory
files = os.listdir(data_path)
# Keep only .gz that do not start with eqtl
if ancestry == "European":
    files = [f for f in files if f.endswith('.gz') and not f.startswith('eqtl')]
elif ancestry == "African":
    files = [f for f in files if f.endswith('.txt') and not f.startswith('eqtl')]
# Then iterate over files
trait_names = []
trait_names_ids = {}
for f in files:
    if f in ['trait_names.txt', 'trait_names_and_ids.txt']:
        continue
    print(f"Done {files.index(f)+1} out of {len(files)}        ", end='\r')
    tool = 'zcat' if ancestry == 'European' else 'cat'
    cmd = f"{tool} {os.path.join(data_path, f)} | head -n 2"
    fname = [x for x in os.popen(cmd) if not x.startswith('id')][0].rstrip().split('\t')[1]
    trait_names.append(fname)
    trait_names_ids[f.replace('.txt', '')] = fname
# Write to file for future use
with open(f"{data_path}/trait_names.txt", 'w') as out:
    for t in trait_names:
        out.write(t + '\n')
# Write the trait names and their corresponding file IDs to a separate file
outf = open(f"{data_path}/trait_names_and_ids.txt", 'w')
outf.write('id\tname\n')
for k, v in trait_names_ids.items():
    key = k + '.txt.gz'
    outf.write(f"{key}\t{v}\n")
outf.close()
# Read from file
trait_names = []
with open(f"{data_path}/trait_names.txt", 'r') as infile:
   for line in infile:
       trait_names.append(line.rstrip())

# ------------------------------------------------------------
# 1) Preprocessing (lightweight; keep meaning, remove boilerplate)
# ------------------------------------------------------------
STOPPHRASES = [r'\bgwas (of|on|for)\b', r'\bgenome[-\s]?wide association (study|analysis)\b', r'\bmeta[-\s]?analysis\b', r'\bstudy of\b', r'\bassociation with\b']
STOPWORDS = set(['the','and','of','for','to','with','in','on','a','an','by','from'])

# Function to normalize trait names: lowercase, remove punctuation, remove stopwords, and collapse whitespace
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
# 2) Load biomedical semantic models
#    BioLORD captures broader biomedical/ontological relationships;
#    SapBERT contributes concept and synonym similarity.
# ------------------------------------------------------------
MODEL_NAMES = [
    'FremyCompany/BioLORD-2023',
    'amrothemich/sapbert-sentence-transformers',  # biomedical concept similarity
]
WEIGHTS = np.array([0.7, 0.3], dtype=np.float32)  # BioLORD, SapBERT

models = [SentenceTransformer(m) for m in MODEL_NAMES]

# ------------------------------------------------------------
# 3) Encode + Ensemble (weighted, L2-normalized)
# ------------------------------------------------------------
# Function to encode a list of texts using multiple models and combine their embeddings into an ensemble representation.
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
    weights = weights / weights.sum()
    # Keep independently trained latent spaces separate. Scaling each normalized
    # embedding by sqrt(weight) makes cosine similarity in the concatenated space
    # equal to the weighted mean of the individual models' cosine similarities.
    weighted_embs = [np.sqrt(weight) * emb for weight, emb in zip(weights, all_embs)]
    ens = normalize(np.concatenate(weighted_embs, axis=1))
    return ens

# ------------------------------------------------------------
# 4) Clustering (cosine distance, hierarchical; threshold is key)
# ------------------------------------------------------------
# Function to cluster embeddings using Agglomerative Clustering based on cosine distance. The distance_threshold parameter controls the granularity of clustering.
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
# 5) Conservatively assign singleton clusters to a larger cluster
# ------------------------------------------------------------
# Function to reassign singleton clusters (clusters with only one member) to the nearest larger cluster based on cosine similarity. This helps reduce isolated traits forming arbitrary new groups.
def reassign_singletons(embeddings, labels, min_similarity=0.40):
    labels = labels.copy()
    unique_labels, counts = np.unique(labels, return_counts=True)
    singleton_labels = unique_labels[counts == 1]
    target_labels = unique_labels[counts > 1]
    if len(singleton_labels) == 0 or len(target_labels) == 0:
        return labels
    # Calculate centroids only from the original non-singleton clusters. This
    # prevents isolated traits from forming arbitrary new groups with each other.
    target_centroids = np.vstack([
        np.mean(embeddings[labels == lab], axis=0) for lab in target_labels
    ])
    target_centroids = normalize(target_centroids)
    reassigned = 0
    for singleton_label in singleton_labels:
        singleton_idx = np.where(labels == singleton_label)[0][0]
        similarities = target_centroids @ embeddings[singleton_idx]
        best_target = int(np.argmax(similarities))
        if similarities[best_target] >= min_similarity:
            labels[singleton_idx] = target_labels[best_target]
            reassigned += 1
    print(
        f"Reassigned {reassigned} of {len(singleton_labels)} singleton clusters "
        f"at minimum cosine similarity {min_similarity:.2f}"
    )
    return labels

# ------------------------------------------------------------
# 6) Representative selection per cluster
#    Option A: pick closest to cluster centroid in ensemble space
#    Option B: shortest normalized string (often nice for labels)
# ------------------------------------------------------------
# Function to pick a representative trait name for each cluster. The method parameter allows choosing between selecting the item closest to the cluster centroid or the shortest normalized string.
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
# 7) Driver: de-duplicate trait names
# ------------------------------------------------------------
# Function to deduplicate trait names by clustering their embeddings and selecting representative names. Returns the list of representative names, a mapping of clusters, labels, representatives, embeddings, and distance matrix.
def deduplicate_traits(trait_names,
                       weights=WEIGHTS,
                       distance_threshold=0.20,
                       singleton_min_similarity=0.40,
                       pick='centroid'):
    # keep original + normalized for reporting
    normed = [normalize_trait(t) for t in trait_names]
    ens = encode_ensemble(normed, models, weights, show_progress_bar=True)
    labels, dist = cluster_embeddings(ens, distance_threshold=distance_threshold)
    labels = reassign_singletons(
        ens, labels, min_similarity=singleton_min_similarity
    )
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

# Main execution
# Run deduplication
# distance_treshold can be tuned; lower = more clusters, higher = fewer clusters [0.3-0.4 is aggressive]
reps_list = []
cluster_map_list = []
labels_list = []
rawreps_list = []
ens_list = []
dist_list = []
# before it was [0.1, 0.2, 0.3, 0.4, 0.5], but now we only do 0.5
for thr in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]:
    print(f"Processing threshold: {thr}")
    reps, cluster_map, labels, reps_raw, ens, dist = deduplicate_traits(trait_names, weights=WEIGHTS, distance_threshold=thr, pick='centroid')
    reps_list.append(reps)
    cluster_map_list.append(cluster_map)
    labels_list.append(labels)
    rawreps_list.append(reps_raw)
    ens_list.append(ens)
    dist_list.append(dist)

# before it was [0.1, 0.2, 0.3, 0.4, 0.5], but now we only do 0.5
for i, thr in enumerate([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]):
    print(f"Threshold: {thr}")
    reps = reps_list[i]
    cluster_map = cluster_map_list[i]
    labels = labels_list[i]
    print(f"Original: {len(trait_names)}  ->  Unique reps: {len(reps)}")

# Thresholds - European
# Threshold: 0.1
# Original: 10604  ->  Unique reps: 1842
# Threshold: 0.2
# Original: 10604  ->  Unique reps: 1662
# Threshold: 0.3
# Original: 10604  ->  Unique reps: 1343
# Threshold: 0.4
# Original: 10604  ->  Unique reps: 978
# Threshold: 0.5
# Original: 10604  ->  Unique reps: 604
# Threshold: 0.6
# Original: 10604  ->  Unique reps: 319
# Threshold: 0.7
# Original: 10604  ->  Unique reps: 136

# Thresholds - African
# Threshold: 0.1
# Original: 751  ->  Unique reps: 236
# Threshold: 0.2
# Original: 751  ->  Unique reps: 197
# Threshold: 0.3
# Original: 751  ->  Unique reps: 172
# Threshold: 0.4
# Original: 751  ->  Unique reps: 159
# Threshold: 0.5
# Original: 751  ->  Unique reps: 136
# Threshold: 0.6
# Original: 751  ->  Unique reps: 109
# Threshold: 0.7
# Original: 751  ->  Unique reps: 66

# Bundle them into one dictionary
data = {
    "reps": reps_list,
    "cluster_map": cluster_map_list,
    "labels": labels_list,
    "rawreps": rawreps_list,
    "ens": ens_list,
    "dist": dist_list
}

# Save to file
with open(f"20260916_trait_clusters_{ancestry}.pkl", "wb") as f:
    pickle.dump(data, f)
#class NumpyCompatUnpickler(pickle.Unpickler):
#    def find_class(self, module, name):
#        # Remap old internal numpy paths to the current ones
#        if module == "numpy._core":
#            module = "numpy.core"
#        elif module == "numpy._core.multiarray":
#            module = "numpy.core.multiarray"
#        return super().find_class(module, name)
#with open("for_snpxplorer/trait_clusters.pkl", "rb") as f:
#    data = NumpyCompatUnpickler(f).load()
#reps_list = data["reps"]
#cluster_map_list = data["cluster_map"]
#labels_list = data["labels"]

# Convert cluster_map to a DataFrame for easier viewing
all_dfs = pd.DataFrame()
# before it was [0.1, 0.2, 0.3, 0.4, 0.5], but now we only do 0.5
for i, thr in enumerate([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]):
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
# Save to CSV for inspection
all_dfs.to_csv(f"20260916_trait_clusters_summary_AI_thresholds_{ancestry}.csv", index=False, sep='\t')
# Read from CSV for inspection
all_dfs = pd.read_csv(f"20260916_trait_clusters_summary_AI_thresholds_{ancestry}.csv", sep='\t')
fig, axes = plot_cluster_threshold_summary(all_dfs)
fig.savefig(f"20260916_trait_clusters_summary_AI_thresholds_{ancestry}.png", dpi=300, bbox_inches='tight')
plt.close(fig)
#all_dfs = pd.read_csv("for_snpxplorer/trait_clusters_summary_AI_thresholds.csv", sep='\t')

# Grep alzh in Representatives
#alzh_df = all_dfs[all_dfs['Representative'].str.contains('alzh', case=False, na=False)]
# in this case, AD is only 1 cluster when using 0.4 and 0.5 threhsolds, while it's multiple clusters for lower thresholds

# Separately, only once, calculate how many associations per study for each trait
flist = os.listdir(data_path)
trait_assoc_counts = {}
for f in flist:
    print(f"Processing file {flist.index(f)+1} out of {len(flist)}        ", end='\r')
    if ancestry == 'European' and f.endswith('.gz') and not f.startswith('eqtl'):
        cmd = f"zcat {os.path.join(data_path, f)} | wc -l"
        count = int(os.popen(cmd).read().strip()) - 1  # subtract header
        fname = [x for x in os.popen(f"zcat {os.path.join(data_path, f)} | head -n 2") if not x.startswith('id')][0].rstrip().split('\t')[1]
        trait_assoc_counts[fname] = [count, f]
    elif ancestry == 'African' and f.endswith('.txt') and not f.startswith('eqtl') and not f.startswith('trait_names'):
        cmd = f"cat {os.path.join(data_path, f)} | wc -l"
        count = int(os.popen(cmd).read().strip()) - 1  # subtract header
        fname = [x for x in os.popen(f"cat {os.path.join(data_path, f)} | head -n 2") if not x.startswith('id')][0].rstrip().split('\t')[1]
        trait_assoc_counts[fname] = [count, f]
# Convert to DataFrame with trait name, association count, file name
assoc_counts_df = pd.DataFrame.from_dict(trait_assoc_counts, orient='index', columns=['Association Count', 'File Name'])
assoc_counts_df['Trait Name'] = assoc_counts_df.index
# Save to file
assoc_counts_df.to_csv(f"{data_path}/trait_association_counts.txt", sep='\t', index=False)
# Do rstrip to all Trait Name entries
assoc_counts_df['Trait Name'] = assoc_counts_df['Trait Name'].str.rstrip()
#assoc_counts_df = pd.read_csv(f"{data_path}/trait_association_counts.txt", sep='\t')

# we need to pick the representative -- basically the study with most associations
# reset index
# iterate over rows of the df
all_dfs['Final Representative'] = ''
all_dfs['Final Representative associations'] = 0
all_dfs['Final Representative file'] = ''
for idx, row in all_dfs.iterrows():
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
all_dfs.to_csv(f"20260916_trait_clusters_summary_AI_thresholds_{ancestry}.csv", index=False)

# ------------------------------------------------------------
# ------------------------------------------------------------
# Visualize embeddings with UMAP (optional)
# derive ens for the 0.5 threshold
ens = ens_list[4]  # index 4 corresponds to threshold 0.5
u = umap.UMAP(n_neighbors=50, min_dist=0.1, metric='cosine', random_state=42).fit_transform(ens)

# Plot
plt.figure(figsize=(24, 10))
plt.scatter(u[:,0], u[:,1], c=labels, cmap='tab20', s=20, alpha=0.7)
plt.title("", fontsize=20)
plt.xlabel("UMAP 1", fontsize=18)
plt.ylabel("UMAP 2", fontsize=18)
plt.tick_params(axis='both', labelsize=16)
plt.savefig(f"20260916_umap_trait_embeddings_{ancestry}.png", dpi=300, bbox_inches='tight')

# save umap and embeddings
with open(f"20260916_trait_embeddings_umap_thr05_{ancestry}.pkl", "wb") as f:
    pickle.dump({"umap": u, "embeddings": ens, "labels": labels, 'cluster_map': cluster_map}, f)

# DIAGNOSTICS & EVALUATION
# ------------------------------------------------------------
# Inspect the five closest distinct neighbors of Alzheimer's disease
# ------------------------------------------------------------
query_trait = "Multiple sclerosis"
query_norm = normalize_trait(query_trait)
query_matches = [i for i, trait in enumerate(trait_names) if normalize_trait(trait) == query_norm]
if not query_matches:
    print(f"Trait not found: {query_trait}")
else:
    query_idx = query_matches[0]
    similarities = ens @ ens[query_idx]  # embeddings are L2-normalized
    ranked_indices = np.argsort(similarities)[::-1]
    neighbors = []
    seen = {query_norm}  # exclude the query and duplicate copies of it
    for idx in ranked_indices:
        neighbor_norm = normalize_trait(trait_names[idx])
        if neighbor_norm in seen:
            continue
        seen.add(neighbor_norm)
        neighbors.append((trait_names[idx], similarities[idx]))
        if len(neighbors) == 5:
            break
    print(f"\nTop 5 closest neighbors of {query_trait!r}:")
    for rank, (trait, similarity) in enumerate(neighbors, start=1):
        print(f"{rank}. {trait} (cosine similarity: {similarity:.4f})")

# ------------------------------------------------------------
# Inspect the five closest clusters to Alzheimer's disease
# ------------------------------------------------------------
if query_matches:
    query_idx = query_matches[0]
    query_cluster = labels[query_idx]
    closest_clusters = []

    for cluster_label in np.unique(labels):
        if cluster_label == query_cluster:
            continue
        member_indices = np.where(labels == cluster_label)[0]
        centroid = normalize(
            np.mean(ens[member_indices], axis=0).reshape(1, -1)
        )[0]
        similarity = ens[query_idx] @ centroid
        closest_clusters.append((cluster_label, similarity, len(member_indices)))

    closest_clusters.sort(key=lambda item: item[1], reverse=True)

    print(f"\nTop 5 closest clusters to {query_trait!r}:")
    for rank, (cluster_label, similarity, size) in enumerate(closest_clusters[:5], start=1):
        representative = cluster_map[int(cluster_label)]['representative']
        print(
            f"{rank}. {representative} "
            f"(cluster {cluster_label}, n={size}, cosine similarity: {similarity:.4f})"
        )
