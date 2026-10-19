import os
import sys
import pandas as pd

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

def add_kc(a, b):
    sa = a.stack().str.split(",", expand=True).astype("int64")
    sb = b.stack().str.split(",", expand=True).astype("int64")
    if sa.shape[1] != sb.shape[1]:
        raise ValueError(f"KC field widths differ: {sa.shape[1]} vs {sb.shape[1]}")
    s = (sa + sb).astype(str)
    joined = s[0]
    for c in s.columns[1:]:
        joined = joined + "," + s[c]
    return joined.unstack().reindex(index=a.index, columns=a.columns)

def merge_pair(a, b, name_a, name_b):
    if list(a.columns) != list(b.columns):
        if set(a.columns) != set(b.columns):
            raise ValueError(f"{name_a} and {name_b} have different columns")
        b = b[a.columns]
    shared = a.index.intersection(b.index)
    only_a = a.index.difference(shared, sort=False)
    only_b = b.index.difference(shared, sort=False)
    summed = add_kc(a.loc[shared], b.loc[shared])
    merged = pd.concat([a.loc[only_a], b.loc[only_b], summed])
    print(f"{name_a}: {len(a)} snpmers | {name_b}: {len(b)} snpmers | "
          f"shared: {len(shared)} | merged: {len(merged)}")
    assert len(merged) == len(a) + len(b) - len(shared)
    return merged, shared


def decrement(occ, snpmers, label):
    """Merging two chromosomes into one: a snpmer present in both was
    counted twice, so remove one occurrence. Modifies `occ` in place."""
    hit = occ.index.intersection(snpmers)
    print(f"  {label}: {len(hit)}/{len(snpmers)} shared snpmers found in occurrence table")
    occ.loc[hit] -= 1

# Map chromosome / pair name -> path using the parent directory name
in_by_chr = {os.path.basename(os.path.dirname(p)): p
             for p in snakemake.input["ac_tsvs"]}
out_by_pair = {os.path.basename(os.path.dirname(p)): p
               for p in snakemake.output["merged"]}

cts_cm = pd.read_csv(snakemake.input["snpmers_cts_cm"], sep="\t").set_index("snpmer")["occurrences"]
cts_no_cm = pd.read_csv(snakemake.input["snpmers_cts_no_cm"], sep="\t").set_index("snpmer")["occurrences"]

for pair, (chr_a, chr_b) in snakemake.params["pairs"].items():
    a = pd.read_csv(in_by_chr[chr_a], sep="\t", index_col=0)
    b = pd.read_csv(in_by_chr[chr_b], sep="\t", index_col=0)
    merged, shared = merge_pair(a, b, chr_a, chr_b)
    merged.to_csv(out_by_pair[pair], sep="\t")
    decrement(cts_cm, shared, f"{pair} cm occ")
    decrement(cts_no_cm, shared, f"{pair} no_cm occ")
    del a, b, merged

cts_cm.rename_axis("snpmer").to_csv(snakemake.output["snpmers_cts_cm"], sep="\t", header=["occurrences"])
cts_no_cm.rename_axis("snpmer").to_csv(snakemake.output["snpmers_cts_no_cm"], sep="\t", header=["occurrences"])