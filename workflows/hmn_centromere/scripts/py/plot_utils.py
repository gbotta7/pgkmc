import itertools
import heapq
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import sys
from scipy.cluster.hierarchy import fcluster, cophenet, linkage, dendrogram, to_tree
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import silhouette_score
from sklearn.metrics.pairwise import cosine_distances

def _categorical_track(labels, cmap_name, color_overrides=None, nan_color="0.75"):
    """Map a 1-D array of labels to an RGBA column plus legend material."""
    labels = np.asarray(labels, dtype=object)
    is_na = pd.isna(labels)
    levels = list(pd.unique(labels[~is_na]))
    levels.sort(key=str)
    n = len(levels)

    if color_overrides:
        base = plt.get_cmap(cmap_name, max(n, 1))
        colors = np.array([
            mpl.colors.to_rgba(color_overrides.get(lev, base(i)))
            for i, lev in enumerate(levels)
        ])
    else:
        colors = plt.get_cmap(cmap_name, max(n, 1))(np.arange(n))

    level_to_idx = {lev: i for i, lev in enumerate(levels)}
    nan_rgba = mpl.colors.to_rgba(nan_color)

    rgba = np.empty((len(labels), 4))
    for i, l in enumerate(labels):
        rgba[i] = nan_rgba if is_na[i] else colors[level_to_idx[l]]

    legend_levels = list(levels)
    legend_colors = list(colors)
    if is_na.any():
        legend_levels.append("NA")
        legend_colors.append(np.array(nan_rgba))

    return rgba[:, None, :], legend_levels, np.array(legend_colors)

def _cluster_link_colors(Z, leaf_clusters, cluster_colors, default="0.35"):
    """link_color_func: colour a link by its cluster, grey where clusters merge."""
    n = Z.shape[0] + 1
    node_cluster = {i: leaf_clusters[i] for i in range(n)}
    for k in range(Z.shape[0]):
        a, b = int(Z[k, 0]), int(Z[k, 1])
        ca, cb = node_cluster[a], node_cluster[b]
        node_cluster[n + k] = ca if ca is not None and ca == cb else None

    def _color(k):
        c = node_cluster.get(k)
        if c is None:
            return default
        return mpl.colors.to_hex(cluster_colors[c])

    return _color


def _spread_labels(desired, min_sep, lo, hi):
    """Push apart labels closer than min_sep, keeping each group centred."""
    desired = np.asarray(desired, float)
    n = len(desired)
    if n == 0:
        return desired
    order = np.argsort(desired)
    d = desired[order]

    def gpos(g):
        c = d[g].mean()
        start = c - (len(g) - 1) * min_sep / 2
        return start + np.arange(len(g)) * min_sep

    groups = [[k] for k in range(n)]
    merged = True
    while merged:
        merged = False
        i = 0
        while i < len(groups) - 1:
            if gpos(groups[i])[-1] + min_sep > gpos(groups[i + 1])[0] + 1e-12:
                groups[i:i + 2] = [groups[i] + groups[i + 1]]
                merged = True
            else:
                i += 1

    out = np.empty(n)
    for g in groups:
        p = gpos(g)
        # keep the whole group inside the axes
        if p[0] < lo:
            p += lo - p[0]
        if p[-1] > hi:
            p -= p[-1] - hi
        out[order[g]] = p
    return out


def _shade_clusters(ax, clusters_ord, cluster_colors, alpha=0.18,
                    fontsize=13):
    """Faded band behind each contiguous run of a cluster, plus a dark label."""
    runs, start = [], 0
    for i in range(1, len(clusters_ord) + 1):
        if i == len(clusters_ord) or clusters_ord[i] != clusters_ord[start]:
            runs.append((clusters_ord[start], start, i))
            start = i

    # a cluster can be split across the leaf order; label only its largest run
    longest = {}
    for lev, s, e in runs:
        if lev not in longest or (e - s) > (longest[lev][1] - longest[lev][0]):
            longest[lev] = (s, e)

    for lev, s, e in runs:
        ax.axhspan(10 * s, 10 * e, color=cluster_colors[lev], alpha=alpha,
                   lw=0, zorder=0)

    for artist in ax.collections + ax.lines:
        artist.set_zorder(2)

    x_root, x_leaf = ax.get_xlim()  # leaves sit at x_leaf (0 for "left")
    span_x = x_leaf - x_root
    f = lambda t: x_root + t * span_x   # 0 = root side, 1 = leaf side

    levs = list(longest)
    y_true = np.array([10 * (longest[l][0] + longest[l][1]) / 2 for l in levs])

    # minimum vertical separation, in data units, from the text height
    y0, y1 = sorted(ax.get_ylim())
    ax_h_pt = ax.get_position().height * ax.figure.get_figheight() * 72
    min_sep = fontsize / max(ax_h_pt, 1e-9) * (y1 - y0)

    y_lab = _spread_labels(y_true, min_sep, y0 + min_sep / 2, y1 - min_sep / 2)

    tol = 0.15 * min_sep
    for lev, yt, yl in zip(levs, y_true, y_lab):
        color = cluster_colors[lev]
        if abs(yl - yt) > tol:
            # elbow: out from the text, then down/up to the band it belongs to
            ax.plot([f(0.030), f(0.016), f(0.016), f(0.004)],
                    [yl, yl, yt, yt],
                    color=color, lw=1.0, solid_capstyle="round", zorder=9,
                    path_effects=[pe.withStroke(linewidth=2.2,
                                                foreground="white")])
            x_text = f(0.034)
        else:
            x_text = f(0.02)

        ax.text(x_text, yl, str(lev),
                color=color, ha="left", va="center",
                fontsize=fontsize, fontweight="bold", zorder=10,
                path_effects=[pe.withStroke(linewidth=2.5, foreground="white")])


def kc_heatmap_clustered(comp,
                         log_depth,
                         Z_cols,
                         col_leaves,
                         sample_ids=None,
                         meta=None,
                         annotate=("Population",),
                         col_clusters=None,
                         row_metric="euclidean",
                         row_method="average",
                         depth_norm_percentile=(2, 98),
                         cmap_comp="coolwarm",
                         figsize=(12, 10),
                         cluster_cmap="tab20",
                         cluster_alpha=0.18,
                         cluster_label_fontsize=12,
                         annot_cmaps=None,
                         annot_nan_color="0.75",
                         cbar_frac=0.25,
                         cbar_fontsize=10,
                         cbar_width=0.018,
                         cbar_gap=0.022,
                         legend_pad=None,
                         save_path=None):
    n_variants, n_samples = comp.shape

    # ---- cluster variants ---------------------------------------------------
    comp_col_ord = comp[:, col_leaves]
    var_condensed = nan_pdist(comp_col_ord, metric=row_metric)
    Z_vars = linkage(var_condensed, method=row_method)
    var_leaves = dendrogram(Z_vars, no_plot=True)["leaves"]

    comp_ord = comp_col_ord[var_leaves, :]
    depth_ord = log_depth[:, col_leaves][var_leaves, :]

    # ---- normalize depth ----------------------------------------------------
    lo, hi = np.percentile(depth_ord, depth_norm_percentile)
    alpha_frac = np.clip((depth_ord - lo) / max(hi - lo, 1e-9), 0, 1)
    alpha = 0.15 + 0.85 * alpha_frac
    nan_mask = np.isnan(comp_ord)
    alpha[nan_mask] = 1.0  # fully opaque where composition is NaN, so grey reads clearly

    cmap = plt.get_cmap(cmap_comp)
    cmap.set_bad(color="0.75")  # grey for NaN composition
    norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
    rgba = cmap(norm(comp_ord))
    rgba[..., 3] = alpha
    rgba = np.transpose(rgba, (1, 0, 2))  # -> (haplotypes, variants, 4)

    # ---- assemble annotation tracks ----------------------------------------
    annot_cmaps = {"Population": "plasma", **(annot_cmaps or {})}
    tracks = []

    if isinstance(annotate, str):  # a bare string would iterate per character
        annotate = (annotate,)

    if meta is not None and annotate:
        if isinstance(meta, str):
            meta = load_meta(meta, sample_ids=sample_ids)
        elif sample_ids is not None:
            meta = meta.loc[np.asarray(sample_ids, dtype=str)]
        elif len(meta) != n_samples:
            raise ValueError(
                "Pass sample_ids so metadata rows can be matched to rows."
            )

        meta_ord = meta.iloc[col_leaves]
        for col_name in annotate:
            col, levels, colors = _categorical_track(
                meta_ord[col_name].values,
                annot_cmaps.get(col_name, "tab20"),
                nan_color=annot_nan_color,
            )
            tracks.append((col_name, col, levels, colors))

    n_tracks = len(tracks)

    # ---- cluster palette (dendrogram bands + link colours) -----------------
    cluster_colors = None
    if col_clusters is not None:
        col_clusters = np.asarray(col_clusters)
        levels = sorted(pd.unique(col_clusters), key=str)
        base = plt.get_cmap(cluster_cmap, max(len(levels), 1))
        cluster_colors = {lev: base(i) for i, lev in enumerate(levels)}

    # ---- figure -------------------------------------------------------------xw
    if legend_pad is None:
        legend_pad = max(1.6, 0.7 * n_tracks)
    width_ratios = [1] + [0.06] * n_tracks + [3, legend_pad]
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(1, len(width_ratios), width_ratios=width_ratios,
                          wspace=0.08)

    ax_tree = fig.add_subplot(gs[0, 0])
    ax_hm = fig.add_subplot(gs[0, 1 + n_tracks])

    link_cf = (None if cluster_colors is None
               else _cluster_link_colors(Z_cols, col_clusters, cluster_colors))
    dendrogram(Z_cols, ax=ax_tree, no_labels=True, color_threshold=0,
               orientation="left", above_threshold_color="black",
               link_color_func=link_cf)
    # ax_tree.set_xlabel("Distance", fontsize=14)
    ax_tree.set_xticks([])
    ax_tree.set_yticks([])
    ax_tree.set_frame_on(False) # remove the borders of tree plot
    # scipy stacks leaves bottom-up; imshow draws row 0 at the top, so flip
    ax_tree.set_ylim(10 * n_samples, 0)

    if cluster_colors is not None:
        _shade_clusters(ax_tree, col_clusters[col_leaves], cluster_colors,
                        alpha=cluster_alpha, fontsize=cluster_label_fontsize)

    # annotation strips carry no titles: the legends already name them
    for i, (_, track_rgba, _, _) in enumerate(tracks, start=1):
        ax_t = fig.add_subplot(gs[0, i])
        ax_t.imshow(track_rgba, aspect="auto", interpolation="none")
        ax_t.set_xticks([])
        ax_t.set_yticks([])

    ax_hm.imshow(rgba, aspect="auto", interpolation="none")
    ax_hm.set_xlabel(f"SNP-mers (n={n_variants})", fontsize=12)
    ax_hm.set_ylabel(f"Haplotypes (m={n_samples})", fontsize=12, labelpad=8)
    ax_hm.yaxis.set_label_position("right")  # left side belongs to the tree
    ax_hm.set_xticks([])
    ax_hm.set_yticks([])

    # ---- legend block -------------------------------------------------------xw
    margin = gs[0, -1].get_position(fig)
    x_center = 0.5 * (margin.x0 + margin.x1)

    legs = []
    for title, _, levels, colors in tracks:
        if len(levels) > 20:
            continue
        handles = [mpl.patches.Patch(color=colors[i], label=str(lev))
                   for i, lev in enumerate(levels)]
        leg = fig.legend(
            handles=handles, title=title, loc="lower left",
            bbox_to_anchor=(0.0, 0.0), ncol=1, frameon=False,
            fontsize=8, title_fontsize=9, labelspacing=0.3,
            handlelength=1.0, handleheight=1.0, handletextpad=0.4,
            borderaxespad=0.0, borderpad=0.0,
        )
        leg._legend_box.align = "left"
        fig.add_artist(leg)
        legs.append(leg)

    fig.canvas.draw()  # needed before legend extents can be measured
    inv = fig.transFigure.inverted()
    boxes = [leg.get_window_extent().transformed(inv) for leg in legs]
    leg_w = [b.width for b in boxes]
    leg_h = max([b.height for b in boxes], default=0.0)

    gap_x = 0.012          # between legend columns
    gap_y = 0.035          # between the legends and the colour scales
    total_w = sum(leg_w) + gap_x * max(len(leg_w) - 1, 0)

    block_h = cbar_frac + (leg_h + gap_y if legs else 0.0)
    y_bottom = 0.5 - block_h / 2

    x = x_center - total_w / 2
    for leg, w in zip(legs, leg_w):
        leg.set_bbox_to_anchor((x, y_bottom + cbar_frac + gap_y),
                               transform=fig.transFigure)
        x += w + gap_x

    # the composition colour scale + dual blue/red depth swatch, side by side
    depth_width = 2 * cbar_width       # blue + red columns share one axes
    pair_w = cbar_width + cbar_gap + depth_width
    x0 = x_center - pair_w / 2
    ax_cbar = fig.add_axes([x0, y_bottom, cbar_width, cbar_frac])
    ax_alpha = fig.add_axes([x0 + cbar_width + cbar_gap, y_bottom,
                             depth_width, cbar_frac])

    # ---- comp colorbar ------------------------------------------------------
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, cax=ax_cbar)
    cbar.set_label("Allele composition", fontsize=cbar_fontsize)
    cbar.ax.tick_params(labelsize=cbar_fontsize - 1)
    # label/ticks on the left so they don't run into the depth strip
    cbar.ax.yaxis.set_ticks_position("left")
    cbar.ax.yaxis.set_label_position("left")

    # ---- depth legend: two swatches (blue & red), same palette as comp -----
    # Sampled directly from `cmap` (white centre -> blue, white centre -> red)
    # instead of alpha-blending pure endpoint colors onto white, so the hues
    # exactly match what the composition colorbar shows at those cmap values.
    depth_max = depth_ord.max()  # true max log1p(N_copies) among plotted cells

    n_steps = 256
    blue_vals = np.linspace(0.5, 0.0, n_steps)  # bottom (0 depth) -> top (max depth)
    red_vals = np.linspace(0.5, 1.0, n_steps)

    swatch_img = np.ones((n_steps, 2, 4))
    swatch_img[:, 0, :3] = cmap(blue_vals)[:, :3]
    swatch_img[:, 1, :3] = cmap(red_vals)[:, :3]
    # alpha stays 1.0 throughout: these are real cmap colors, not blended

    ax_alpha.imshow(swatch_img, aspect="auto", extent=[0, 2, 0, 1], origin="lower")
    ax_alpha.set_facecolor("white")
    ax_alpha.set_xticks([])
    ax_alpha.set_yticks([0.0, 1.0])
    ax_alpha.set_yticklabels(["0", f"{depth_max:.1f}"], fontsize=cbar_fontsize - 1)
    ax_alpha.yaxis.tick_right()
    ax_alpha.set_ylabel("Log1p(N_copies)", fontsize=cbar_fontsize)
    ax_alpha.yaxis.set_label_position("right")

    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
    # plt.show()

def nan_pdist(X, metric="euclidean"):
    X = np.asarray(X, dtype=float)
    n, m = X.shape
    mask = np.isfinite(X)
    Xf = np.where(mask, X, 0.0)

    out = np.empty(n * (n - 1) // 2)
    idx = 0
    for i in range(n):
        xi, mi = Xf[i], mask[i]
        for j in range(i + 1, n):
            xj, mj = Xf[j], mask[j]
            both = mi & mj
            k = both.sum()
            if k == 0:
                out[idx] = np.nan
                idx += 1
                continue

            if metric == "euclidean":
                diff = (xi - xj)[both]
                d = np.sqrt(np.sum(diff ** 2) * (m / k))
            elif metric == "cityblock":
                diff = (xi - xj)[both]
                d = np.sum(np.abs(diff)) * (m / k)
            elif metric == "cosine":
                a, b = xi[both], xj[both]
                na, nb = np.linalg.norm(a), np.linalg.norm(b)
                d = 1.0 - (a @ b) / (na * nb) if na > 0 and nb > 0 else np.nan
            else:
                raise ValueError(f"Unsupported metric for nan_pdist: {metric}")

            out[idx] = d
            idx += 1

    bad = ~np.isfinite(out)
    if bad.any():
        finite = out[~bad]
        fill = finite.max() * 1.5 if finite.size and finite.max() > 0 else 1.0
        print(f"nan_pdist: {bad.sum()}/{out.size} pair(s) share no finite "
              f"overlap; filling with {fill:.4g} (1.5x max observed distance) "
              f"so linkage() gets an all-finite matrix.", file=sys.stderr)
        out[bad] = fill

    return out