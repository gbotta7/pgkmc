# import pandas as pd
# import sys

# sys.stdout = sys.stderr = open(snakemake.log[0], "w")

# def join_snpmer_counts(tsv_list):
#     result = None
#     for path in sorted(tsv_list):
#         df = pd.read_csv(path, sep="\t", usecols=[0, 1], dtype=str, compression="gzip")
#         df = df.set_index(df.columns[0])
#         if result is None:
#             result = df
#         else:
#             result = result.join(df, how="outer")
#     return result

# tsv_list = snakemake.input["tsv_list"]
# df = join_snpmer_counts(tsv_list)
# df.to_csv(snakemake.output["joined_tsv"], sep="\t", compression="gzip")

# import gzip
import pandas as pd
import sys

sys.stdout = sys.stderr = open(snakemake.log[0], "w", buffering=1)

CHUNKSIZE = 100_000  # tune for memory/speed tradeoff

def join_snpmer_counts_streaming(tsv_list, out_path, chunksize=CHUNKSIZE):
    paths = sorted(tsv_list)
    readers = [
        pd.read_csv(p, sep="\t", usecols=[0, 1], dtype=str, compression="gzip", chunksize=chunksize)
        for p in paths
    ]

    i = 0
    with open(out_path, "wt") as out:
        header_written = False
        chunk_num = 0
        while True:
            i += 1
            if i % 1_000_000 == 0:
                print(f"Processed {i} SNP-mers\n", flush=True)
            
            chunks = []
            for r in readers:
                try:
                    chunks.append(next(r))
                except StopIteration:
                    chunks.append(None)
            
            if all(c is None for c in chunks):
                break

            chunk_num += 1
            id_col = chunks[0].columns[0]
            ids = chunks[0][id_col].reset_index(drop=True)

            merged = pd.concat(
                [chunks[0].set_index(id_col)] +
                [c.set_index(c.columns[0]) for c in chunks[1:]],
                axis=1,
            )

            # drop rows where every sample's value is the literal "0,0"
            all_zero = (merged == "0,0").all(axis=1)
            merged = merged[~all_zero]

            merged.to_csv(out, sep="\t", header=not header_written)
            header_written = True

tsv_list = snakemake.input["tsv_list"]
join_snpmer_counts_streaming(tsv_list, snakemake.output["joined_tsv"])