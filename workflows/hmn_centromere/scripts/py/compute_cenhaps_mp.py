import itertools
import heapq
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import re
import sys
from itertools import combinations_with_replacement
from collections import deque
from scipy.cluster.hierarchy import fcluster, cophenet, linkage, dendrogram, to_tree
from scipy.spatial.distance import pdist, squareform
from scipy.stats import beta
from sklearn.metrics import normalized_mutual_info_score
from sklearn.metrics.pairwise import cosine_distances

from plot_utils import *
from tree_utils import *

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

META_COLS = ["Hap", "Sex_chrom", "Sample", "Sex", "Population", "Country", "Extra"]

def to_hap_key(s):
    sample, suf = re.match(r"\d+_(.+)\.(pat|mat|hap1|hap2)$", s).groups()
    return f"{sample}.{suffix_map[suf]}"

def load_kc_matrix(tsv_path, compression=None):
    """
    Parse a joined snpmer TSV (one column of "a,b" counts per sample) into a
    3D numpy array. Returns kc[i, j, 0] = first allele count,
    kc[i, j, 1] = second allele count, for snpmer i, sample j.
    """
    df = pd.read_csv(tsv_path, compression=compression, sep="\t", index_col=0, dtype=str).fillna("0,0")

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

def filter_samples_by_missingness(kc, sample_ids, max_missing_frac):
    is_missing = (kc[:, :, 0] == 0) & (kc[:, :, 1] == 0)
    missing_frac = is_missing.mean(axis=0)

    keep_mask = missing_frac <= max_missing_frac

    kc_filt = kc[:, keep_mask, :]
    sample_ids_filt = [s for s, keep in zip(sample_ids, keep_mask) if keep]

    n_removed = (~keep_mask).sum()
    print(f"[info] removed {n_removed}/{len(keep_mask)} samples with >{max_missing_frac:.0%} (0,0) across SNPs")

    return kc_filt, sample_ids_filt, keep_mask

# def select_snpmers(comp, n_select=200, k=10):
#     X = np.where(np.isnan(comp), np.nanmean(comp, axis=1, keepdims=True), comp)  # fill gaps with the SNP-mer mean
#     X = X - X.mean(axis=1, keepdims=True) # centre each SNP-mer
#     U, s, Vt = np.linalg.svd(X, full_matrices=False) # PCA
#     score = ((U[:, :k] * s[:k]) ** 2).sum(axis=1) # variance of each SNP-mer in the top-k PCs
#     sel = np.zeros(len(comp), bool)
#     sel[np.argsort(-score)[:n_select]] = True # keep the N highest
#     return sel
def knee(w, m=30):
    """Index (1-based) of the knee of a decreasing curve: the point furthest from the straight line
    joining the first and the m-th value."""
    m = min(m, len(w))
    y = w[:m] / w[0]
    x = np.arange(m) / (m - 1)
    d = np.abs((y[-1] - y[0]) * x - (x[-1] - x[0]) * (y - y[0])) / np.hypot(y[-1] - y[0], x[-1] - x[0])
    return int(np.argmax(d)) + 1

def select_snpmers(comp, max_select=500, k=None, alpha=1e-7, min_sd=0.10, min_call_rate=0.5, min_k=3, max_k=20):
    """
    Keep the SNP-mers whose variation across haplotypes is shared with the top-k PCs, up to max_select.
 
    A SNP-mer counts as informative if the share of its variance in the top-k PCs (its "communality")
    is higher than a structureless SNP-mer would reach by chance (the 1 - alpha quantile of a
    Beta(k/2, (H-1-k)/2) null), its sd across haplotypes is >= min_sd, and it is called in at least
    min_call_rate of the haplotypes. The number kept is the number of informative SNP-mers, capped at
    max_select; if there are more than max_select, the ones with the most variance in the top-k PCs are kept.
    alpha is deliberately tiny: the PCs are fitted on the same data, so structureless SNP-mers reach a higher
    communality than the Beta null predicts, and alpha = 1e-3 lets a lot of noise through.
    """
    H = comp.shape[1]
    valid = (~np.isnan(comp)).mean(axis=1) >= min_call_rate                       # enough calls
    Xv = comp[valid]
    X = np.where(np.isnan(Xv), np.nanmean(Xv, axis=1, keepdims=True), Xv)         # fill gaps with the SNP-mer mean
    X = X - X.mean(axis=1, keepdims=True)                                         # centre each SNP-mer
    U, s, Vt = np.linalg.svd(X, full_matrices=False)   
    # Knee plot
    if k is None:
        k = int(np.clip(knee(s ** 2), min_k, max(min_k, min(max_k, H - 2))))
    # PCA
    captured = ((U[:, :k] * s[:k]) ** 2).sum(axis=1)                              # variance of each SNP-mer in the top-k PCs
    total = (X ** 2).sum(axis=1)                                                  # its total variance
    communality = captured / np.maximum(total, 1e-12)                             # share of it that is shared structure
    thr = beta.ppf(1 - alpha, k / 2, max(H - 1 - k, 1) / 2)                       # communality reached by chance
    ok = (communality > thr) & (np.sqrt(total / H) >= min_sd)
    idx = np.where(ok)[0]
    n_pass = len(idx)
    idx = idx[np.argsort(-captured[idx])][:max_select]                            # keep the top max_select
    print(f"[info] select_snpmers: {int(valid.sum())}/{len(comp)} SNP-mers with enough calls, "
          f"{n_pass} above the chance threshold ({thr:.3f}), keeping {len(idx)} (max {max_select})")
    sel = np.zeros(len(comp), bool)
    sel[np.where(valid)[0][idx]] = True
    return sel

def load_t2t_status(t2t_path, chrom):
    """Per-haplotype T2T status for one chromosome, keyed by '<sample>.<1|2>'.
    Returns a dict: key -> 'contig' | 'scaffold'. Haplotypes missing from the
    table for this chromosome are incomplete."""
    t2t = pd.read_csv(t2t_path, sep="\t", dtype=str)
    chrom = "chr" + str(chrom).removeprefix("chr")   # accepts "1" or "chr1"
    t2t = t2t[t2t["chr_name"] == chrom].copy()
    t2t["key"] = t2t["hap"].map(to_hap_key)
    return dict(zip(t2t["key"], t2t["level"]))

#### PIPELINE
kc, variant_info, sample_ids = load_kc_matrix(snakemake.input["joined_tsv"], None)
all_sample_ids = list(sample_ids)

kc, sample_ids, keep_mask = filter_samples_by_missingness(kc, sample_ids, snakemake.params["max_missing_frac"])
maf = snakemake.params["maf"]
kc_filt1, variant_info_filt1, keep_mask = filter_by_maf_binarized(kc, variant_info, maf=maf)
print(f"Kept {keep_mask.sum()} / {len(keep_mask)} variants (MAF > {maf}, presence/absence)")

# Keep only SNP-mers that occur exactly once in the centromere and never outside the centromere
cts_cm = pd.read_csv(snakemake.input["snpmers_cts_cm"], sep="\t").set_index("snpmer")["occurrences"]
cts_no_cm = pd.read_csv(snakemake.input["snpmers_cts_no_cm"], sep="\t").set_index("snpmer")["occurrences"]

keep_spec = []
for v in variant_info_filt1:
    s = v["snpmer"]
    keep_spec.append(cts_cm.get(s, 0) == 1 and cts_no_cm.get(s, 0) == 0)
keep_spec = np.array(keep_spec)
kc_filt2 = kc_filt1[keep_spec]
variant_info_filt2 = [v for v, k in zip(variant_info_filt1, keep_spec) if k]

comp_filt2, log_depth_filt2 = kc_get_features(kc_filt2)

# Keep most variable SNP-mers
sel = select_snpmers(comp_filt2, max_select=snakemake.params["max_snpmers"])
kc_filt3 = kc_filt2[sel]
variant_info_filt3 = [v for v, k in zip(variant_info_filt2, sel) if k]
comp_filt3, log_depth_filt3 = kc_get_features(kc_filt3)

# Reshape the matrix and stack both alleles
n, m, c = kc_filt3.shape
X = np.log1p(np.asarray(kc_filt3, dtype=np.float64))
X = np.transpose(X, (0, 2, 1)).reshape(2 * n, m)
X = X.reshape(2 * n, m)

Z_cols, col_leaves = get_tree(X.T, cosine_distances)
D = cosine_distances(X.T)

col_clusters, best_k, scores, cut_height = assign_clusters_from_tree(Z_cols, D)

meta = pd.read_csv(snakemake.input["meta"], sep="\t", header=None, names=META_COLS, dtype=str, na_values=[], keep_default_na=False)
anno = pd.read_csv(snakemake.input["anno"], sep="\t", names=["Hap", "Cenhaps_Karen"])
# Build the join key on df1: pat -> .1, mat -> .2
suffix_map = {
    "pat": "1", "mat": "2",
    "hap1": "1", "hap2": "2",
}
meta["Hap"] = meta["Sample"] + "." + meta["Hap"].str.extract(r"\.(pat|mat|hap1|hap2)$")[0].map(suffix_map)
meta_merged = meta.merge(anno, on="Hap", how="left").set_index("Hap")

# Add T2T information
t2t_level = load_t2t_status(snakemake.input["t2t"], snakemake.params["chr"])
lvl = pd.Series(meta_merged.index.map(t2t_level), index=meta_merged.index)   # NaN if not in table
meta_merged["Complete"] = np.where(lvl.notna(), "complete", "incomplete")
meta_merged["Level"] = lvl.where(lvl.notna()) 

sample_ids = [to_hap_key(s) for s in sample_ids]

kc_heatmap_clustered(comp_filt3,
                     log_depth_filt3,
                     Z_cols,
                     col_leaves,
                     sample_ids=sample_ids,
                     meta=meta_merged,
                    #  annotate=("Population", "Cenhaps_Karen", "Complete", "Level"),
                    annotate=("Population", "Cenhaps_Karen"),
                     col_clusters=col_clusters,
                     save_path=snakemake.output["hm"]
                     )

# NMI between the clusters and the Cenhaps_Karen annotation
ann = meta_merged[~meta_merged.index.duplicated()]["Cenhaps_Karen"].reindex(sample_ids)  # same order as col_clusters
has_ann = ann.notna().to_numpy()
nmi = normalized_mutual_info_score(ann[has_ann].astype(str), np.asarray(col_clusters)[has_ann])
print(f"NMI(col_clusters, Cenhaps_Karen) = {nmi:.3f} on {has_ann.sum()}/{len(has_ann)} haplotypes with an annotation")

### Compute profiles of genotypes with all SNP-mers
col_clusters = np.asarray(col_clusters)
clusters = np.unique(col_clusters)
snpmer_names = [v["snpmer"] for v in variant_info_filt3]
cluster_mean = np.column_stack(
    [np.nanmean(comp_filt3[:, col_clusters == k], axis=1) for k in clusters]
)

pairs = list(combinations_with_replacement(range(len(clusters)), 2))
geno_labels = [f"{clusters[a]}/{clusters[b]}" for a, b in pairs]
geno_mean = np.column_stack(
    [np.nanmean(cluster_mean[:, [a, b]], axis=1) for a, b in pairs]
)
geno_df = pd.DataFrame(geno_mean, index=snpmer_names, columns=geno_labels)
geno_df["chr"] = str(snakemake.params["chr"])
geno_df.to_csv(snakemake.output["gt_profiles"], sep="\t", float_format="%.3f", index_label="snpmer")


### Genotype per sample, encoded as the column index in gt_profiles
# The profile columns are in the same order as `pairs`, so a sample's genotype
# index is the position of its sorted (a, b) cluster-position pair in `pairs`.
# 0-based; the trailing "chr" column is not counted. "?" = a haplotype is
# missing/filtered, so no pair can be assigned.
cluster_pos = {c: i for i, c in enumerate(clusters)}     # cluster label -> position in `clusters`
pair_to_idx = {pair: i for i, pair in enumerate(pairs)}  # (a, b) with a <= b -> column index

hap_to_cluster = dict(zip(sample_ids, col_clusters))     # sample_ids are already "<sample>.<1|2>" keys

all_hap_keys = [to_hap_key(s) for s in all_sample_ids]
all_samples = sorted({h.rsplit(".", 1)[0] for h in all_hap_keys})

records = []
for sample in all_samples:
    found = [hap_to_cluster.get(f"{sample}.{i}") for i in ("1", "2")]
    if any(c is None for c in found):
        records.append({"sample": sample, "genotype": "?"})
        continue
    a, b = sorted(cluster_pos[c] for c in found)
    records.append({"sample": sample, "genotype": pair_to_idx[(a, b)]})

pd.DataFrame(records).to_csv(snakemake.output["gt_clusters"], sep="\t", index=False)