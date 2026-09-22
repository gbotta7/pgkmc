import csv
import sys

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

trees = snakemake.input["trees"]
snpmers = snakemake.input["snpmers"]
tree_outpath = snakemake.output["tree"]
snpmers_outpath = snakemake.output["snpmers"]

header = None
rows = []

# Merge trees
for t in trees:
    with open(t, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        file_header = next(reader)

        if header is None:
            header = file_header

        for row in reader:
            rows.append(row)

with open(tree_outpath, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerows(rows)

rows = []
# Merge snpmers
for s in snpmers:
    with open(s, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        file_header = next(reader)

        if header is None:
            header = file_header

        for row in reader:
            rows.append(row)

with open(snpmers_outpath, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerows(rows)