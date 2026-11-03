#!/bin/bash

# Subsample centromere sequences from fasta file
seqtk subseq input.fa.gz cm_regions.bed \
    | gzip -nc > test.fa.gz

# Subsample 1000 sequences from short-read file
samtools view -@4 -u -T ref.fa -s 9.00001 input.cram \
    | samtools collate -@4 -O -u - \
    | samtools fastq -@4 -n -0 /dev/null -s /dev/null - \
    | head -n 4000 \
    | gzip -nc > test.fq.gz