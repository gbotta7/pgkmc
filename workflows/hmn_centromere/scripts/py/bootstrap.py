import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import cosine_distances
import sys

from tree_utils import *

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

joined_tsv = snakemake.input["joined_tsv"]
max_missing_frac = snakemake.params["max_missing_frac"]
maf = snakemake.params["maf"]
n_boot = snakemake.params["n_boot"]
seed = 0

def load_kc_matrix(tsv_path):
    df = pd.read_csv(tsv_path, sep="\t", index_col=0, dtype=str).fillna("0,0")
    kc = np.zeros((len(df), df.shape[1], 2), dtype=np.uint16)
    for j, sample in enumerate(df.columns):
        kc[:, j, :] = df[sample].str.split(",", expand=True).astype(int).to_numpy()
    return kc


def filter_samples(kc, max_missing_frac):
    missing = ((kc[:, :, 0] == 0) & (kc[:, :, 1] == 0)).mean(axis=0)
    return kc[:, missing <= max_missing_frac, :]


def filter_maf(kc, maf):
    return kc[(kc > 0).mean(axis=1).min(axis=1) > maf]


def build_tree(X3):
    """X3: (n SNP-mers, 2 alleles, m haplotypes) -> linkage, leaf order, cluster labels."""
    n, _, m = X3.shape
    X = X3.reshape(2 * n, m)
    Z, leaves = get_tree(X.T, cosine_distances)
    clusters, _, _, _ = assign_clusters_from_tree(Z, cosine_distances(X.T))
    return Z, leaves, np.asarray(clusters)


def retained(ref_idx, boot_labels):
    """Per haplotype: is it in the replicate cluster that best matches (Jaccard) its reference cluster?"""
    _, b_idx = np.unique(boot_labels, return_inverse=True)
    C = np.zeros((ref_idx.max() + 1, b_idx.max() + 1))
    np.add.at(C, (ref_idx, b_idx), 1)
    J = C / (C.sum(1)[:, None] + C.sum(0)[None, :] - C)
    return b_idx == J.argmax(1)[ref_idx]


kc = filter_maf(filter_samples(load_kc_matrix(joined_tsv), max_missing_frac), maf)
X3 = np.transpose(np.log1p(kc.astype(np.float64)), (0, 2, 1))   # (n, 2, m)
n, _, m = X3.shape

Z, leaves, labels = build_tree(X3)
_, ref_idx = np.unique(labels, return_inverse=True)

rng = np.random.default_rng(seed)
kept = np.zeros(m)
for rep in range(n_boot):
    _, _, cl_b = build_tree(X3[rng.integers(0, n, n)])
    kept += retained(ref_idx, cl_b)
kept /= n_boot

members = {i: np.array([i]) for i in range(m)}
splits = []
for i, (a, b, _, _) in enumerate(Z):
    a, b = int(a), int(b)
    members[m + i] = np.concatenate([members[a], members[b]])
    if len(np.unique(labels[members[m + i]])) >= 2:
        splits.append((m + i, kept[members[a]].mean(), kept[members[b]].mean()))
split_id = [s[0] for s in splits]
split_val = [(s[1] + s[2]) / 2 for s in splits]

fig, ax = plt.subplots(figsize=(14, 7))
x = {int(l): i for i, l in enumerate(leaves)}
y = {i: 0.0 for i in range(m)}
for i, (a, b, h, _) in enumerate(Z):
    a, b = int(a), int(b)
    ax.plot([x[a], x[a], x[b], x[b]], [y[a], h, h, y[b]], c="0.5", lw=0.6, zorder=1)
    x[m + i], y[m + i] = (x[a] + x[b]) / 2, h

sc = ax.scatter([x[i] for i in split_id], [y[i] for i in split_id], c=split_val, s=60,
                cmap="RdYlGn", vmin=0, vmax=1, edgecolor="k", linewidth=0.6, zorder=3)
for i, v in zip(split_id, split_val):
    ax.annotate(f"{v:.0%}", (x[i], y[i]), fontsize=8, xytext=(0, 6), textcoords="offset points", ha="center")

ymin = -0.05 * Z[:, 2].max()
ax.imshow(kept[leaves][None, :], cmap="RdYlGn", vmin=0, vmax=1, aspect="auto",
          extent=(-0.5, m - 0.5, ymin, 0), zorder=0)

ax.set_xlim(-1, m); ax.set_ylim(ymin, Z[:, 2].max() * 1.05)
ax.set_xticks([]); ax.set_ylabel("Cosine distance")
# ax.set_title(f"{CHR}: dots = splits between clusters (mean retention of the two sides)\n"
#              "strip = fraction of bootstrap replicates each haplotype stays in its cluster")
fig.colorbar(sc, ax=ax, fraction=0.025, pad=0.01, shrink=0.5, label="Retention")
fig.tight_layout()
fig.savefig(snakemake.output["boot_tree_plot"], dpi=300)