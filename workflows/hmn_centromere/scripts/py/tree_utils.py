import itertools
import heapq
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import re
import sys
from scipy.cluster.hierarchy import fcluster, cophenet, linkage, dendrogram, to_tree
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import silhouette_score
from sklearn.metrics.pairwise import cosine_distances

def get_tree(X, distance_func, method="ward"):
    D = np.asarray(distance_func(X), dtype=np.float64)
    col_condensed = squareform(D, checks=False)

    Z_cols = linkage(col_condensed, method=method)
    col_dendro = dendrogram(Z_cols, no_plot=True)
    col_leaves = col_dendro["leaves"]

    return Z_cols, col_leaves

def height_for_k_clusters(Z, k):
    n = Z.shape[0] + 1
    if k <= 1:
        return Z[-1, 2]
    if k >= n:
        return 0.0
    return Z[n - 1 - k, 2]
    
def assign_clusters_from_tree(Z, dist_matrix=None, k_range=(3, 20), min_cluster_size=1, plot_diagnostics=False):
    n = Z.shape[0] + 1

    if dist_matrix is None:
        coph_condensed = cophenet(Z)
        dist_matrix = squareform(coph_condensed)

    lo, hi = k_range
    scores = {}
    labelings = {}
    for k in range(lo, hi + 1):
        labels_k = fcluster(Z, t=k, criterion="maxclust")
        _, counts = np.unique(labels_k, return_counts=True)
        if counts.min() < min_cluster_size:
            continue
        if len(np.unique(labels_k)) < 2:
            continue
        score = silhouette_score(dist_matrix, labels_k, metric="precomputed")
        scores[k] = score
        labelings[k] = labels_k

    if not scores:
        raise ValueError(
            "No candidate k in k_range produced valid clusters "
            f"(min_cluster_size={min_cluster_size}). Try widening k_range "
            "or lowering min_cluster_size."
        )

    best_k = max(scores, key=scores.get)
    labels = labelings[best_k]
    height = height_for_k_clusters(Z, best_k)

    if plot_diagnostics:
        fig, ax = plt.subplots(figsize=(5, 3))
        ks = sorted(scores)
        ax.plot(ks, [scores[k] for k in ks], marker="o")
        ax.axvline(best_k, color="red", linestyle="--", alpha=0.6)
        ax.set_xlabel("k (number of clusters)")
        ax.set_ylabel("Silhouette score")
        ax.set_title(f"Best k = {best_k} (silhouette = {scores[best_k]:.3f})")
        plt.tight_layout()
        plt.show()

    return labels, best_k, scores, height

def _leaf_sets(Z):
    """node_id -> set of original leaf indices under that node."""
    n = Z.shape[0] + 1
    leaves = {i: {i} for i in range(n)}
    for k in range(Z.shape[0]):
        a, b = int(Z[k, 0]), int(Z[k, 1])
        leaves[n + k] = leaves[a] | leaves[b]
    return leaves

def find_splits_above_height(Z, cut_height):
    n = Z.shape[0] + 1
    leaves = _leaf_sets(Z)
    out = []
    for k in range(Z.shape[0]):
        height = Z[k, 2]
        if height <= cut_height:
            continue
        a, b = int(Z[k, 0]), int(Z[k, 1])
        out.append({
            "node": n + k,
            "height": height,
            "left": np.array(sorted(leaves[a])),
            "right": np.array(sorted(leaves[b])),
        })
    return out

def find_perfect_markers(comp, left_idx, right_idx, exclude=None):
    L = comp[:, left_idx]
    R = comp[:, right_idx]

    left_present = np.isfinite(L)
    right_present = np.isfinite(R)

    has_left = left_present.any(axis=1)
    has_right = right_present.any(axis=1)

    left_all_present = left_present.all(axis=1)
    right_all_present = right_present.all(axis=1)
    left_all_na = ~left_present.any(axis=1)
    right_all_na = ~right_present.any(axis=1)

    left_all_zero = np.where(left_present, L == 0, True).all(axis=1)
    left_all_one  = np.where(left_present, L == 1, True).all(axis=1)
    right_all_zero = np.where(right_present, R == 0, True).all(axis=1)
    right_all_one  = np.where(right_present, R == 1, True).all(axis=1)

    # Case A: normal perfect split -- both sides have data and disagree
    dir_left_zero = left_all_zero & right_all_one & has_left & has_right
    dir_left_one  = left_all_one & right_all_zero & has_left & has_right

    # Case B: one side entirely NA, other side fully called (values can be mixed)
    dir_left_none  = left_all_na & right_all_present
    dir_right_none = right_all_na & left_all_present

    perfect = dir_left_zero | dir_left_one | dir_left_none | dir_right_none
    marker_idx = np.nonzero(perfect)[0]

    covered = {}
    directions = {}
    for i in marker_idx:
        if dir_left_none[i]:
            covered[i] = np.concatenate([left_idx, right_idx])
            directions[i] = "left_none"
        elif dir_right_none[i]:
            covered[i] = np.concatenate([left_idx, right_idx])
            directions[i] = "right_none"
        elif dir_left_zero[i]:
            covered[i] = np.concatenate([left_idx[left_present[i]], right_idx[right_present[i]]])
            directions[i] = "left_zero"
        else:  # dir_left_one[i]
            covered[i] = np.concatenate([left_idx[left_present[i]], right_idx[right_present[i]]])
            directions[i] = "left_one"

    return marker_idx, covered, directions

def _node_children(Z, node, n):
    """Return (a, b) child node ids for a merge node, or None if node is a leaf."""
    if node < n:
        return None
    k = node - n
    return int(Z[k, 0]), int(Z[k, 1])

def build_snpmer_tree(Z, comp, cut_height=0.0, marker_names=None):
    n = Z.shape[0] + 1
    leaves = _leaf_sets(Z)
    used = set()
    root = 2 * n - 2

    if marker_names is None:
        marker_names = np.arange(comp.shape[0])
    else:
        marker_names = np.asarray(marker_names)

    def pick_markers(left_idx, right_idx):
        """
        Greedy set cover: repeatedly pick the not-yet-used marker that covers
        the most currently-uncovered samples in this split, until every
        sample in left_idx/right_idx is covered by at least one chosen
        marker, or no more candidates help.
    
        Case A markers ("left_zero"/"left_one" -- an actual allelic difference)
        are preferred over "left_none"/"right_none". The none-markers are only
        used to mop up samples that no Case A marker could cover.
        """
        marker_idx, covered, directions = find_perfect_markers(comp, left_idx, right_idx)
        all_samples = set(left_idx.tolist()) | set(right_idx.tolist())
    
        def greedy_cover(candidate_pool, remaining):
            chosen = []
            candidates = {i for i in candidate_pool if i not in used}
            while remaining and candidates:
                best = max(candidates, key=lambda i: len(remaining & set(covered[i].tolist())))
                gain = len(remaining & set(covered[best].tolist()))
                if gain == 0:
                    break
                chosen.append({
                    "marker": int(best),
                    "direction": directions[best],
                    "covered": covered[best],
                    "n_covered": int(len(covered[best])),
                })
                used.add(best)
                candidates.discard(best)
                remaining -= set(covered[best].tolist())
            return chosen, remaining
    
        case_a_idx = [i for i in marker_idx if directions[i] in ("left_zero", "left_one")]
        none_idx = [i for i in marker_idx if directions[i] in ("left_none", "right_none")]
    
        remaining = set(all_samples)
        chosen, remaining = greedy_cover(case_a_idx, remaining)
    
        if remaining:
            more, remaining = greedy_cover(none_idx, remaining)
            chosen.extend(more)
    
        return chosen, sorted(remaining)

    def recurse(node):
        children = _node_children(Z, node, n)
        if children is None:
            return {"leaf": int(node)}

        a, b = children
        height = Z[node - n, 2]
        left_idx = np.array(sorted(leaves[a]))
        right_idx = np.array(sorted(leaves[b]))

        if height <= cut_height:
            return {"node": int(node), "height": float(height),
                    "leaves": left_idx.tolist() + right_idx.tolist()}

        markers, uncovered = pick_markers(left_idx, right_idx)
        for m in markers:
            m["marker_name"] = marker_names[m["marker"]].item()

        return {
            "node": int(node),
            "height": float(height),
            "markers": markers,              # list of {marker, marker_name, direction, n_covered}
            "uncovered_samples": uncovered,  # samples no available marker could resolve
            "left": recurse(a),
            "right": recurse(b),
        }

    return recurse(root)


def classify_sample_tree_walk(tree, genotype):
    """
    Walk the tree top-down like decision tree, following the markers'
    majority vote at each node. If a node's markers give no information for
    this sample, branch into both children rather than guessing.
    """
    results = []

    def recurse(node, matches, observed):
        if "leaf" in node:
            results.append((("leaf", node["leaf"]), matches, observed))
            return
        if not node.get("markers") and "leaves" in node:
            results.append((("flat", node["node"]), matches, observed))
            return

        left_votes = 0
        right_votes = 0
        local_observed = 0
        for m in node["markers"]:
            g = genotype[m["marker"]]
            direction = m["direction"]

            if direction == "left_none":
                # left side is all-NA, right side is fully called:
                # presence/absence of a call IS the signal -- always informative
                local_observed += 1
                if not np.isfinite(g):
                    left_votes += 1
                else:
                    right_votes += 1
                continue

            if direction == "right_none":
                local_observed += 1
                if not np.isfinite(g):
                    right_votes += 1
                else:
                    left_votes += 1
                continue

            # Case A markers ("left_zero" / "left_one"): only informative
            # when this sample actually has a call
            if not np.isfinite(g):
                continue
            local_observed += 1
            expected_left = 0 if direction == "left_zero" else 1
            if g == expected_left:
                left_votes += 1
            else:
                right_votes += 1

        if left_votes == right_votes:
            # no informative markers, or a genuine tie: explore both branches
            recurse(node["left"], matches, observed)
            recurse(node["right"], matches, observed)
        elif left_votes > right_votes:
            recurse(node["left"], matches + left_votes, observed + local_observed)
        else:
            recurse(node["right"], matches + right_votes, observed + local_observed)

    recurse(tree, 0, 0)
    return results


def classify_sample_by_tree_walk(tree, genotype):
    candidates = classify_sample_tree_walk(tree, genotype)
    if not candidates:
        return None, None, 0
    scored = []
    for key, matches, observed in candidates:
        score = matches / observed if observed > 0 else 0.0
        scored.append((score, observed, key))
    scored.sort(key=lambda r: (-r[0], -r[1]))
    best_score, best_observed, best_key = scored[0]
    return best_key, best_score, best_observed


def assign_clusters_by_tree_walk(tree, comp):
    n_samples = comp.shape[1]
    labels = np.zeros(n_samples, dtype=int)
    scores = np.full(n_samples, np.nan)
    n_observed = np.zeros(n_samples, dtype=int)
    key_to_label = {}
    next_id = itertools.count(1)

    for s in range(n_samples):
        genotype = comp[:, s]
        key, score, observed = classify_sample_by_tree_walk(tree, genotype)
        if key is None:
            key = ("unresolved", s)
            score = np.nan
        if key not in key_to_label:
            key_to_label[key] = next(next_id)
        labels[s] = key_to_label[key]
        scores[s] = score
        n_observed[s] = observed

    return labels, scores, n_observed