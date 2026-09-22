import pandas as pd
import sys

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

def join_snpmer_counts(tsv_list):
    cols = []
    for path in sorted(tsv_list):
        df = pd.read_csv(path, sep="\t", usecols=[0, 1], dtype=str)
        cols.append(df.set_index(df.columns[0])[df.columns[1]])
    return pd.concat(cols, axis=1)

tsv_list = snakemake.input["tsv_list"]
df = join_snpmer_counts(tsv_list)
df.to_csv(snakemake.output["joined_tsv"], sep="\t", compression="gzip")
