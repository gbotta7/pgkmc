import os
import sys
import pandas as pd

sys.stdout = sys.stderr = open(snakemake.log[0], "w")


def zero_fill(df):
    """Zero string matching the KC field width of df (e.g. '0,0')."""
    if df.shape[0] == 0 or df.shape[1] == 0:
        return "0,0"
    width = str(df.iloc[0, 0]).count(",") + 1
    return ",".join(["0"] * width)


def merge_pair(a, b, name_a, name_b):
    # Suffix columns with the chromosome name so the two tables never collide
    a = a.add_suffix(f"_{name_a}")
    b = b.add_suffix(f"_{name_b}")

    overlap = a.columns.intersection(b.columns)
    if len(overlap) > 0:
        raise ValueError(
            f"{name_a} and {name_b} still share column names after suffixing: "
            f"{list(overlap)[:5]}"
        )

    shared = a.index.intersection(b.index)
    only_b = b.index.difference(a.index, sort=False)
    union_idx = a.index.append(only_b)

    fill_a = zero_fill(a)
    fill_b = zero_fill(b)

    # Keep all columns from both; snpmers missing from one side get 0,0
    a_full = a.reindex(union_idx).fillna(fill_a)
    b_full = b.reindex(union_idx).fillna(fill_b)
    merged = pd.concat([a_full, b_full], axis=1)

    print(f"{name_a}: {len(a)} snpmers, {a.shape[1]} cols | "
          f"{name_b}: {len(b)} snpmers, {b.shape[1]} cols | "
          f"shared: {len(shared)} | merged: {merged.shape[0]} snpmers, "
          f"{merged.shape[1]} cols")
    assert len(merged) == len(a) + len(b) - len(shared)
    assert merged.shape[1] == a.shape[1] + b.shape[1]
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
    a = pd.read_csv(in_by_chr[chr_a], sep="\t", index_col=0, dtype=str)
    b = pd.read_csv(in_by_chr[chr_b], sep="\t", index_col=0, dtype=str)
    merged, shared = merge_pair(a, b, chr_a, chr_b)
    merged.to_csv(out_by_pair[pair], sep="\t")
    decrement(cts_cm, shared, f"{pair} cm occ")
    decrement(cts_no_cm, shared, f"{pair} no_cm occ")
    del a, b, merged

cts_cm.rename_axis("snpmer").to_csv(snakemake.output["snpmers_cts_cm"], sep="\t", header=["occurrences"])
cts_no_cm.rename_axis("snpmer").to_csv(snakemake.output["snpmers_cts_no_cm"], sep="\t", header=["occurrences"])