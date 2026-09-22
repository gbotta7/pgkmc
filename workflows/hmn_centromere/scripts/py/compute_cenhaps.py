import itertools
import heapq
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import re
import sys
from collections import deque
from scipy.cluster.hierarchy import fcluster, cophenet, linkage, dendrogram, to_tree
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import silhouette_score, normalized_mutual_info_score
from sklearn.metrics.pairwise import cosine_distances

from plot_utils import *
from tree_utils import *

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

META_COLS = ["Hap", "Sex_chrom", "Sample", "Sex", "Population", "Country", "Extra"]

def to_hap_key(s):
    sample, suf = re.match(r"\d+_(.+)\.(pat|mat|hap1|hap2)$", s).groups()
    return f"{sample}.{suffix_map[suf]}"

def load_kc_matrix(tsv_path):
    """
    Parse a joined snpmer TSV (one column of "a,b" counts per sample) into a
    3D numpy array. Returns kc[i, j, 0] = first allele count,
    kc[i, j, 1] = second allele count, for snpmer i, sample j.
    """
    df = pd.read_csv(tsv_path, sep="\t", index_col=0, dtype=str).fillna("0,0")

    sample_ids = list(df.columns)
    snpmers = list(df.index)

    kc = np.zeros((len(snpmers), len(sample_ids), 2), dtype=np.uint16)
    for j, sample in enumerate(sample_ids):
        for i, kc_str in enumerate(df[sample]):
            parts = kc_str.split(",")
            if len(parts) != 2:
                print(f"[warn] unexpected KC field {kc_str!r} at sample {sample}, snpmer {snpmers[i]}, skipping")
                continue
            kc[i, j, 0] = int(parts[0])
            kc[i, j, 1] = int(parts[1])

    variant_info = []
    for s in snpmers:
        left, rest = s.split("[")
        alleles, right = rest.split("]")
        ref, alt = alleles.split("/")
        variant_info.append({"snpmer": s, "ref": ref, "alt": alt,
                             "left": left, "right": right})

    return kc, variant_info, sample_ids

def filter_by_maf_binarized(kc, variant_info, maf=0.05):
    kc_bin = (kc > 0).astype(np.uint8)
    allele_freqs = kc_bin.mean(axis=1)
    af = allele_freqs.min(axis=1)
    keep_mask = af > maf

    kc_filt = kc[keep_mask]
    variant_info_filt = [v for v, keep in zip(variant_info, keep_mask) if keep]

    return kc_filt, variant_info_filt, keep_mask

def kc_get_features(X, dtype=np.float64):
    X = np.asarray(X, dtype=dtype)

    tot = X[:,:,0] + X[:,:, 1]
    comp = np.where(tot > 0, X[..., 0] / tot, np.nan)
    # comp = X[:,:,0] / (tot + eps)
    log_depth = np.log1p(tot)

    return comp, log_depth


#### PIPELINE
kc, variant_info, sample_ids = load_kc_matrix(snakemake.input["joined_tsv"])
kc_filt, variant_info_filt, keep_mask = filter_by_maf_binarized(kc, variant_info, maf=snakemake.params["maf"])
print(f"Kept {keep_mask.sum()} / {len(keep_mask)} variants (MAF > {snakemake.params["maf"]}, presence/absence)")

# Reshape the matrix and stack both alleles
n, m, c = kc_filt.shape
X = np.log1p(np.asarray(kc_filt, dtype=np.float64))
X = np.transpose(X, (0, 2, 1)).reshape(2 * n, m) 
X = X.reshape(2 * n, m)

Z_cols, col_leaves = get_tree(X.T, cosine_distances)
D = cosine_distances(X.T)

col_clusters, best_k, scores, cut_height = assign_clusters_from_tree(Z_cols, D)

# Extract comp and log_depth for plotting
comp, log_depth = kc_get_features(kc_filt)

meta = pd.read_csv(snakemake.input["meta"], sep="\t", header=None, names=META_COLS, dtype=str, na_values=[], keep_default_na=False)
anno = pd.read_csv(snakemake.input["anno"], sep="\t", names=["Hap", "Cenhaps_Karen"])
# Build the join key on df1: pat -> .1, mat -> .2
suffix_map = {
    "pat": "1", "mat": "2",
    "hap1": "1", "hap2": "2",
}
meta["Hap"] = meta["Sample"] + "." + meta["Hap"].str.extract(r"\.(pat|mat|hap1|hap2)$")[0].map(suffix_map)
meta_merged = meta.merge(anno, on="Hap", how="left").set_index("Hap")

sample_ids = [to_hap_key(s) for s in sample_ids]

kc_heatmap_clustered(comp,
                     log_depth,
                     Z_cols,
                     col_leaves,
                     sample_ids=sample_ids,
                     meta=meta_merged,
                     annotate=("Population", "Cenhaps_Karen"),
                     col_clusters=col_clusters,
                     save_path=snakemake.output["hm"]
                     )

out_df = pd.DataFrame({"sample": sample_ids, "cenhap": col_clusters})
out_df.to_csv(snakemake.output["cenhaps"], sep="\t", index=False)

# Compute the SNP-mer decision tree on all SNP-mers
comp_full, _ = kc_get_features(kc)
tree = build_snpmer_tree(Z_cols, comp_full, cut_height)

labels_walk, scores_walk, n_obs_walk = assign_clusters_by_tree_walk(tree, comp_full)

kc_heatmap_clustered(comp,
                     log_depth,
                     Z_cols,
                     col_leaves,
                     sample_ids=sample_ids,
                     meta=meta_merged,
                     annotate=("Population", "Cenhaps_Karen"),
                     col_clusters=labels_walk,
                     save_path=snakemake.output["hm_walk"]
                     )

# Write NMI
nmi = normalized_mutual_info_score(col_clusters, labels_walk)
with open(snakemake.output["nmi"], "w") as f:
    f.write("chr\tNMI\n")
    f.write(f"{snakemake.params['chr']}\t{nmi}\n")


from collections import deque

def write_tree_files(chr, tree, variant_info, snpmer_path, structure_path):
    """
    BFS over the tree, numbering internal (splitting) nodes 1..n in
    top-to-bottom order. Writes:
      - snpmer_path: SNPmer\tnode_id\tmove   (one row per SNPmer)
      - structure_path: node_id\tparent_id\tside\tleft_child\tright_child
        (child ids are 0 if that side is a leaf / has no further split)
    """
    node_counter = 0

    def is_leaf(node):
        return "leaf" in node or "leaves" in node

    node_counter += 1
    root_id = node_counter
    queue = deque([(tree, root_id, 0, "ROOT")])

    with open(snpmer_path, "w") as f_snp, open(structure_path, "w") as f_struct:
        f_snp.write("SNPmer\tnode_id\tmove\tchr\n")
        f_struct.write("node_id\tparent_id\tside\tleft_child\tright_child\tchr\n")

        while queue:
            node, node_id, parent_id, side = queue.popleft()

            left, right = node["left"], node["right"]
            left_id = 0
            right_id = 0

            if not is_leaf(left):
                node_counter += 1
                left_id = node_counter
                queue.append((left, left_id, node_id, "L"))

            if not is_leaf(right):
                node_counter += 1
                right_id = node_counter
                queue.append((right, right_id, node_id, "R"))

            for m in node["markers"]:
                snpmer = variant_info[m["marker"]]["snpmer"]
                move = m["direction"]
                f_snp.write(f"{snpmer}\t{node_id}\t{move}\t{chr}\n")

            f_struct.write(f"{node_id}\t{parent_id}\t{side}\t{left_id}\t{right_id}\t{chr}\n")

    return node_counter

n_nodes = write_tree_files(
    snakemake.params['chr'],
    tree,
    variant_info,
    snakemake.output["snpmers"],
    snakemake.output["tree_structure"],
)
print(f"Wrote {n_nodes} nodes")