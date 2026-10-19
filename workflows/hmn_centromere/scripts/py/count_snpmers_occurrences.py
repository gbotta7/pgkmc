import pandas as pd
import numpy as np
import sys

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

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

chrom_counts = {}
for tsv_path in snakemake.input["joined_tsvs"]:
    print(f"Loading {tsv_path}\n")
    kc, variant_info, sample_ids = load_kc_matrix(tsv_path, None)

    allele0_sum = kc[:, :, 0].sum(axis=1)
    allele1_sum = kc[:, :, 1].sum(axis=1)

    for i, info in enumerate(variant_info):
        if allele0_sum[i] > 0 or allele1_sum[i] > 0:
            snpmer = info["snpmer"]
            chrom_counts[snpmer] = chrom_counts.get(snpmer, 0) + 1

with open(snakemake.output["snpmers_cts"], "w") as f:
    f.write("snpmer\toccurrences\n")
    for snpmer, count in sorted(chrom_counts.items(), key=lambda x: -x[1]):
        f.write(f"{snpmer}\t{count}\n")