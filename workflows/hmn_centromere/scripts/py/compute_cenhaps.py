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
from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score
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

def beam_select_snpmers(comp_s, labels, max_select=None, beam_width=3, tol=0.005):
    """
    Find a small set of SNP-mers that reproduces the genotype labels with the nearest-neighbour
    classifier used by `pgkmc clust`.
 
    comp_s : (n_snpmers x n_samples) mean composition of each sample, NaN = absent.
    labels : genotype label ("a/b") of each sample, from the clustering of all SNP-mers.
 
    Objective (same geometry as the haplotype tree): every SNP-mer gives two coordinates (c, 1-c), and
    (0, 0) when absent; the distance between two samples is the cosine distance over the selected SNP-mers.
    NOTE: this is not the distance used by the classifier in clust_mp.c if that one is RMSE.
    Score of a set = accuracy * (1 + mean margin), over the query samples (those whose label is shared with
    at least one other sample; a label that occurs once cannot be recovered):
      accuracy : fraction of queries whose nearest other sample (the sample itself is excluded:
                 leave-one-out) has the right label;
      margin_i : distance to the nearest sample with another label - distance to the nearest sample with
                 the same label, clipped to [-1, 1]. It is positive when the call is right and negative
                 when it is wrong, so a near miss costs less than a clear miss.
    (1 + mean margin) is in [0, 2], so the score is never negative and is 0 only at accuracy 0; accuracy
    dominates and the margin ranks sets with similar accuracy by how clearly they separate the labels.
    Search: every SNP-mer is a candidate at every step; there is no pre-ranking. Beam search grows the set
    one SNP-mer at a time, keeping the beam_width best sets at each size. There is no early stopping: the
    set grows until all SNP-mers are in it (or max_select, if given). The result is the smallest set
    whose score is within a fraction `tol` of the best score over all sizes.
    Distances are updated incrementally (dot products and norms are sums over SNP-mers).
    """
    n, N = comp_s.shape
    _, inv = np.unique(np.asarray(labels).astype(str), return_inverse=True)
    counts = np.bincount(inv)
    q = np.where(counts[inv] >= 2)[0]                       # query samples
    if len(q) < 2 or len(np.unique(inv[q])) < 2:
        raise RuntimeError("beam_select_snpmers needs at least two labels shared by two or more samples")
 
    P = ~np.isnan(comp_s)
    U = np.where(P, comp_s, 0.0).astype(np.float32)         # allele-1 coordinate, 0 if absent
    V = np.where(P, 1.0 - comp_s, 0.0).astype(np.float32)   # allele-2 coordinate, 0 if absent
    Pf = P.astype(np.float32)
    R = U * U + V * V                                       # squared-norm contribution of each SNP-mer
 
    inv_q = inv[q]
    nq = len(q)
    SELF = np.zeros((nq, N), bool)
    SELF[np.arange(nq), q] = True
    SAME = (inv_q[:, None] == inv[None, :]) & ~SELF
    OTHER = ~SAME & ~SELF
 
    def score_from(DOT, NRM, SH):
        denom = np.sqrt(np.outer(NRM[q], NRM))
        with np.errstate(divide="ignore", invalid="ignore"):
            D = 1.0 - DOT / denom
        D = np.where((SH > 0) & ~SELF, D, np.inf)           # no SNP-mer in common -> not a candidate
        nn = D.argmin(1)
        dmin = D[np.arange(nq), nn]
        correct = (inv[nn] == inv_q) & np.isfinite(dmin)
        Dm = np.where(np.isfinite(D), D, 1.0)
        d_same = np.where(SAME, Dm, 2.0).min(1)
        d_other = np.where(OTHER, Dm, 2.0).min(1)
        M = float(np.clip(d_other - d_same, -1.0, 1.0).mean())     # signed margin, all queries
        acc = float(correct.mean())
        return acc * (1.0 + M), acc, M
 
    def mats(idx):                                          # dot / norm / shared-count for a whole set
        idx = list(idx)
        DOT = U[idx][:, q].T @ U[idx] + V[idx][:, q].T @ V[idx]
        SH = Pf[idx][:, q].T @ Pf[idx]
        return DOT, R[idx].sum(0), SH
 
    def extend(m, i):                                       # add SNP-mer i to a set
        DOT, NRM, SH = m
        return (DOT + np.outer(U[i, q], U[i]) + np.outer(V[i, q], V[i]),
                NRM + R[i], SH + np.outer(Pf[i, q], Pf[i]))
 
    empty = (np.zeros((nq, N), np.float32), np.zeros(N, np.float32), np.zeros((nq, N), np.float32))
    beam = [((), empty)]
    history = []                                            # best set at every size: (size, score, acc, chosen)
    n_steps = n if max_select is None else min(max_select, n)
    for step in range(n_steps):
        pool = {}
        for pi, (chosen, m) in enumerate(beam):
            for i in range(n):
                if i in chosen:
                    continue
                key = tuple(sorted(chosen + (int(i),)))
                if key in pool:                             # same set reached through another path
                    continue
                sc, acc, mg = score_from(*extend(m, i))
                pool[key] = (sc, acc, pi, int(i))
        top = sorted(pool.values(), key=lambda t: -t[0])[:beam_width]
        beam = [(beam[pi][0] + (i,), extend(beam[pi][1], i)) for _, _, pi, i in top]
        history.append((step + 1, top[0][0], top[0][1], beam[0][0]))
        if (step + 1) % 10 == 0:
            print(f"[info] beam search: size {step + 1}, score {top[0][0]:.4f}, LOO accuracy {top[0][1]:.3f}")
 
    top_score = max(h[1] for h in history)
    size, score, acc, chosen = next(h for h in history if h[1] >= top_score - tol * abs(top_score))
    base = score_from(*mats(range(n)))
    print(f"[info] beam_select_snpmers: {len(q)}/{N} query samples, {n} SNP-mers; "
          f"kept {size} SNP-mers (LOO accuracy {acc:.3f}, score {score:.4f}); "
          f"all {n} SNP-mers: LOO accuracy {base[1]:.3f}, score {base[0]:.4f}")
    return np.array(sorted(chosen), dtype=int)

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
# 1. Filters (no feature selection): samples, MAF, centromere specificity
suffix_map = {
    "pat": "1", "mat": "2",
    "hap1": "1", "hap2": "2",
}

kc, variant_info, sample_ids = load_kc_matrix(snakemake.input["joined_tsv"], None)
all_sample_ids = list(sample_ids)

kc, sample_ids, keep_mask = filter_samples_by_missingness(kc, sample_ids, snakemake.params["max_missing_frac"])
maf = snakemake.params["maf"]
kc_filt1, variant_info_filt1, keep_mask = filter_by_maf_binarized(kc, variant_info, maf=maf)
print(f"Kept {keep_mask.sum()} / {len(keep_mask)} variants (MAF > {maf}, presence/absence)")

comp_filt1, log_depth_filt1 = kc_get_features(kc_filt1)

# Keep only SNP-mers that occur exactly in only one centromere and never outside the centromeres
cts_cm = pd.read_csv(snakemake.input["snpmers_cts_cm"], sep="\t").set_index("snpmer")["occurrences"]
cts_no_cm = pd.read_csv(snakemake.input["snpmers_cts_no_cm"], sep="\t").set_index("snpmer")["occurrences"]

keep_spec = []
for v in variant_info_filt1:
    s = v["snpmer"]
    keep_spec.append(cts_cm.get(s, 0) == 1 and cts_no_cm.get(s, 0) == 0)
keep_spec = np.array(keep_spec)
kc_filt2 = kc_filt1[keep_spec]
variant_info_filt2 = [v for v, k in zip(variant_info_filt1, keep_spec) if k]
print(f"Kept {len(variant_info_filt2)} / {len(keep_spec)} SNP-mers specific to the centromere")

comp_filt2, log_depth_filt2 = kc_get_features(kc_filt2)

# 2. Cluster the haplotypes on all the filtered SNP-mers, but also cross-chromosomal (log1p counts of both alleles, absent = 0)
n, m, c = kc_filt1.shape
X = np.log1p(np.asarray(kc_filt1, dtype=np.float64))
X = np.transpose(X, (0, 2, 1)).reshape(2 * n, m)

Z_cols, col_leaves = get_tree(X.T, cosine_distances)
D = cosine_distances(X.T)

col_clusters, best_k, scores, cut_height = assign_clusters_from_tree(Z_cols, D)

sample_ids = [to_hap_key(s) for s in sample_ids] # "<sample>.<1|2>"

# 3. Genotype samples with the clustering from all the SNP-mers
col_clusters = np.asarray(col_clusters)
clusters = np.unique(col_clusters)
cluster_pos = {c: i for i, c in enumerate(clusters)}

hap_col = {h: j for j, h in enumerate(sample_ids)}
hap_to_cluster = dict(zip(sample_ids, col_clusters))

all_hap_keys = [to_hap_key(s) for s in all_sample_ids]
all_samples = sorted({h.rsplit(".", 1)[0] for h in all_hap_keys})

mat_names, mat_cols, mat_labels, records = [], [], [], []
for sample in all_samples:
    h1, h2 = f"{sample}.1", f"{sample}.2"
    if h1 not in hap_col or h2 not in hap_col:
        records.append({"sample": sample, "genotype": "?"})
        continue
    a, b = sorted((cluster_pos[hap_to_cluster[h1]], cluster_pos[hap_to_cluster[h2]]))
    label = f"{clusters[a]}/{clusters[b]}"
    col = np.nanmean(comp_filt1[:, [hap_col[h1], hap_col[h2]]], axis=1)
    mat_names.append(sample)
    mat_cols.append(col)
    mat_labels.append(label)
    records.append({"sample": sample, "genotype": label})

comp_samples = np.column_stack(mat_cols) # (SNP-mers x samples)
comp_samples_spec = comp_samples[keep_spec]

# 4. Beam search for a small set of SNP-mers that reproduces the labels
sel_idx = beam_select_snpmers(
    comp_samples_spec, mat_labels,
    max_select=snakemake.params["max_snpmers"],
    beam_width=snakemake.params["beam_width"],
    tol=snakemake.params["tol"],
)
kc_filt3 = kc_filt2[sel_idx]
variant_info_filt3 = [variant_info_filt2[i] for i in sel_idx]
comp_filt3, log_depth_filt3 = kc_get_features(kc_filt3)
snpmer_names = [v["snpmer"] for v in variant_info_filt3]

meta = pd.read_csv(snakemake.input["meta"], sep="\t", header=None, names=META_COLS, dtype=str, na_values=[], keep_default_na=False)
anno = pd.read_csv(snakemake.input["anno"], sep="\t", names=["Hap", "Cenhaps_KHMiga"])
# Build the join key on df1: pat -> .1, mat -> .2
meta["Hap"] = meta["Sample"] + "." + meta["Hap"].str.extract(r"\.(pat|mat|hap1|hap2)$")[0].map(suffix_map)
meta_merged = meta.merge(anno, on="Hap", how="left").set_index("Hap")

# Add T2T information
t2t_level = load_t2t_status(snakemake.input["t2t"], snakemake.params["chr"])
lvl = pd.Series(meta_merged.index.map(t2t_level), index=meta_merged.index)   # NaN if not in table
meta_merged["Complete"] = np.where(lvl.notna(), "complete", "incomplete")
meta_merged["Level"] = lvl.where(lvl.notna())

# Heatmap: tree from all SNP-mers, all SNP-mers vs selected SNP-mers
kc_heatmap_clustered(comp_filt1,
                     log_depth_filt1,
                     Z_cols,
                     col_leaves,
                     sample_ids=sample_ids,
                     meta=meta_merged,
                     annotate=("Population", "Cenhaps_KHMiga"),
                     col_clusters=col_clusters,
                     save_path=snakemake.output["hm"]
                     )

kc_heatmap_clustered(comp_filt3,
                     log_depth_filt3,
                     Z_cols,
                     col_leaves,
                     sample_ids=sample_ids,
                     meta=meta_merged,
                     annotate=("Population", "Cenhaps_KHMiga"),
                     col_clusters=col_clusters,
                     save_path=snakemake.output["hm_sel"]
                     )

# NMI between the clusters and the Cenhaps_KHMiga annotation
ann = meta_merged[~meta_merged.index.duplicated()]["Cenhaps_KHMiga"].reindex(sample_ids)  # same order as col_clusters
has_ann = ann.notna().to_numpy()
nmi = normalized_mutual_info_score(ann[has_ann].astype(str), np.asarray(col_clusters)[has_ann])
print(f"NMI(col_clusters, Cenhaps_KHMiga) = {nmi:.3f} on {has_ann.sum()}/{len(has_ann)} haplotypes with an annotation")

### Outputs: matrix of the selected SNP-mers (SNP-mers x samples) + genotype label of each sample
# Only samples with a genotype are in the matrix, so every matrix column has a label.
gt_matrix = pd.DataFrame(comp_samples_spec[sel_idx], index=snpmer_names, columns=mat_names)
gt_matrix["chr"] = str(snakemake.params["chr"])
gt_matrix.to_csv(snakemake.output["gt_matrix"], sep="\t", float_format="%.3f", na_rep="NA", index_label="snpmer")

pd.DataFrame(records).to_csv(snakemake.output["gt_labels"], sep="\t", index=False)