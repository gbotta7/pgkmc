import csv
import sys

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

gt_profiles_in = snakemake.input["gt_profiles"]
gt_profiles_out = snakemake.output["gt_profiles"]

header = None
rows = []

# Merge genotype profiles
for gt in gt_profiles_in:
    with open(gt, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        file_header = next(reader)

        if header is None:
            header = file_header

        for row in reader:
            rows.append(row)

with open(gt_profiles_out, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerows(rows)