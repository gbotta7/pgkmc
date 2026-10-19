import os
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

def parse_positions(pos_field):
    """'SAMPLE#HAP#CONTIG:ALLELE:POS,...' -> list of (full_contig_id, pos)."""
    out = []
    if not isinstance(pos_field, str) or not pos_field.strip():
        return out
    for entry in pos_field.split(","):
        entry = entry.strip()
        contig_id, _allele, pos = entry.rsplit(":", 2)
        out.append((contig_id, int(pos)))
    return out

def strip_suffix(path, suffixes):
    base = os.path.basename(path)
    for suf in suffixes:
        if base.endswith(suf):
            return base[: -len(suf)]
    return base

cts_paths = snakemake.input["cts"]
bed_paths = snakemake.input["bed"]
index_path = snakemake.input["t2t"]
chrom = snakemake.params["chr"]


index_df = pd.read_csv(index_path, sep="\t")
valid_haps = set(index_df.loc[index_df["chr_name"] == chrom, "hap"])
print(f"{chrom}: {len(valid_haps)} haplotypes have a T2T assembled contig for this chromosome (of {index_df['hap'].nunique()} total)")

rng = np.random.default_rng(42)

records = []  # (haplotype_label, relative_position)
for cts_path, bed_path in zip(cts_paths, bed_paths):
    hap_label = strip_suffix(cts_path, ["_unfilt.centromeric.snpmers.tsv"])
    if hap_label not in valid_haps:
        print(f"{hap_label}: skipped (no complete assembly for {chrom} in index)")
        continue

    if os.path.getsize(bed_path) == 0:
        print(f"{hap_label}: skipped (empty bed file, no centromeric intervals for {chrom})")
        continue
    bed = pd.read_csv(bed_path, sep="\t", header=None, skiprows=1)
    hor_mask = bed[3].str.match(r"^active_hor", case=False, na=False)
    hor_bed = bed[hor_mask]
    bed_spans = {}
    for contig, group in hor_bed.groupby(0):
        bed_spans[contig] = (int(group[1].min()), int(group[2].max()))

    df = pd.read_csv(cts_path, sep="\t")
    pos_col = df.columns[-1]  # "positions"

    hap_positions = []
    for pos_field in df[pos_col].fillna(""):
        for contig_id, pos in parse_positions(pos_field):
            span = bed_spans.get(contig_id)
            if span is None:
                continue
            start, end = span
            if start <= pos <= end:
                length = end - start
                if length <= 0:
                    continue
                rel = (pos - start) / length
                if 0.0 <= rel <= 1.0:
                    hap_positions.append(rel)

    n_points = len(hap_positions)
    idx = rng.choice(n_points, size=int(n_points * snakemake.params["downsample_freq"]), replace=False)
    hap_positions = [hap_positions[i] for i in idx]
    records.extend((hap_label, rel) for rel in hap_positions)

    print(f"{hap_label}: {n_points} SNPmers found\n")

os.makedirs(os.path.dirname(snakemake.output["density_plot"]), exist_ok=True)
rel_df = pd.DataFrame(records, columns=["haplotype", "relative_position"])

# Plot density of SNPmers across the centromere
fig, ax = plt.subplots(figsize=(10, 6))
xs = np.linspace(0, 1, 500)

haplotypes = sorted(rel_df["haplotype"].unique())
for hap in haplotypes:
    vals = rel_df.loc[rel_df["haplotype"] == hap, "relative_position"].values
    if len(vals) < 2:
        continue
    kde = gaussian_kde(vals, bw_method=0.01)
    ax.plot(xs, kde(xs), color="steelblue", alpha=0.25, linewidth=1)

all_vals = rel_df["relative_position"].values
kde_all = gaussian_kde(all_vals, bw_method=0.01)
ax.plot(xs, kde_all(xs), color="navy", alpha=1.0, linewidth=2, label=f"Avg over haplotypes (n={len(haplotypes)})")

ax.set_xlabel("Relative position across centromere (0 = start, 1 = end)")
ax.set_ylabel("Density")
ax.set_title(f"SNPmer density across the {chrom} centromere")
ax.legend()
fig.tight_layout()
fig.savefig(snakemake.output["density_plot"], dpi=300)