# Script to align semantic clustering results

# conda activate nlpy

# Libraries
import pickle
import numpy as np
import umap
import matplotlib.pyplot as plt

# Data paths
african_path = "20260916_trait_embeddings_umap_thr05_African.pkl"
european_path = "20260714_trait_embeddings_umap_thr05.pkl"

# Load data
# African data
datasets = {}
with open(african_path, "rb") as f:
    datasets["African"] = pickle.load(f)
# European data
with open(european_path, "rb") as f:
    datasets["European"] = pickle.load(f)

# Combine embeddings
combined_embeddings = np.vstack([datasets["European"]["embeddings"], datasets["African"]["embeddings"]])

# Combined umap
coordinates = umap.UMAP(n_neighbors=50, min_dist=0.1, metric="cosine", random_state=42).fit_transform(combined_embeddings)

# Separate cluster identifies for each ancestry
cluster_keys = [(ancestry, int(label)) for ancestry in ["European", "African"] for label in np.unique(datasets[ancestry]["labels"])]
palette = plt.get_cmap("turbo")(np.linspace(0, 1, len(cluster_keys)))
cluster_colors = dict(zip(cluster_keys, palette))

# Combined plot
fig, ax = plt.subplots(figsize=(14, 10))
offset = 0
for ancestry, marker in [("European", "o"), ("African", "^")]:
    labels = datasets[ancestry]["labels"]
    n = len(labels)
    xy = coordinates[offset:offset + n]
    colors = [cluster_colors[(ancestry, int(label))] for label in labels]
    ax.scatter(
        xy[:, 0], xy[:, 1],
        c=colors,
        marker=marker,
        s=25,
        alpha=0.7,
        linewidths=0,
    )
    ax.scatter([], [], color="gray", marker=marker, label=ancestry)
    offset += n

ax.set(xlabel="UMAP 1", ylabel="UMAP 2")
ax.legend(title="GWAS ancestry")
fig.tight_layout()
fig.savefig("20260918_combined_umap_clusters.png", dpi=300)

# Try to make a focus plot on hypertension
query = "triglycerides"
fig, ax = plt.subplots(figsize=(12, 9))
offset = 0
for ancestry, marker in [("European", "o"), ("African", "^")]:
    data = datasets[ancestry]
    labels = np.asarray(data["labels"])
    xy = coordinates[offset:offset + len(labels)]
    offset += len(labels)
    selected = {
        int(cluster_id): info
        for cluster_id, info in data["cluster_map"].items()
        if any(query.casefold() in name.casefold()
               for name in info["members"])
    }
    for cluster_id, info in selected.items():
        mask = labels == cluster_id
        representative = info["representative"]
        ax.scatter(
            xy[mask, 0], xy[mask, 1],
            color=cluster_colors[(ancestry, cluster_id)],
            marker=marker,
            s=45,
            alpha=0.8,
            label=f"{ancestry}: {representative}",
        )
        print(
            f"{ancestry} | cluster {cluster_id} | "
            f"{mask.sum()} records | {representative}"
        )
if ax.collections:
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 0.98), fontsize=14)
ax.set(xlabel="UMAP 1", ylabel="UMAP 2", title=f"Clusters containing '{query}'")
ax.set_xlabel("UMAP 1", fontsize=18)
ax.set_ylabel("UMAP 2", fontsize=18)
fig.tight_layout()
fig.savefig("20260918_combined_umap_triglycerides_clusters.png", dpi=300)

# Save pkl object with embeddings and cluster information
with open("20260918_combined_umap_embeddings.pkl", "wb") as f:
    pickle.dump({
        "coordinates": coordinates,
        "datasets": datasets,
        "cluster_colors": cluster_colors,
    }, f)
