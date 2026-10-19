import glob

# Count SNP-mers and rare k-mers in centromeres
# rule detect_snpmers:
#     input:
#         fa_files = expand(
#             os.path.join(WORKDIR, "data/asm/{hap}.fa.gz"),
#             hap=hap_list
#         )
#     output:
#         snpmers = os.path.join(WORKDIR, "results/snps/genomes/chrALL/hmn462_unfilt.genomic.snpmers.txt")
#     params:
#         k = config["kmer_length"],
#         msf = 0,
#         f = 0
#     threads:
#         config["pgkmc_detect_threads"]
#     log:
#         "logs/snps/genomes/chrALL/detect_centromeric_snpmers_hmn462_unfilt.log"
#     shell:
#         """
#         pgkmc detect -v \
#             --snp \
#             -k {params.k} \
#             --msf {params.msf} \
#             -t {threads} \
#             -f {params.f} \
#             -o {output.snpmers} \
#             {input.fa_files} \
#             2> {log}
#         """

# rule count_snpmers_no_cm:
#     input:
#         fa = os.path.join(WORKDIR, "data/asm/{hap}.fa.gz"),
#         bed = os.path.join(WORKDIR, "results/chrBed/{chr}/{hap}.no_cm.bed"),
#         snpmers = os.path.join(WORKDIR, "results/snps/genomes/chrALL/hmn462_unfilt.genomic.snpmers.txt")
#     output:
#         cts = os.path.join(WORKDIR, "results/snps/genomes/{chr}/{hap}_unfilt.non_centromeric.snpmers.tsv.gz")
#     params:
#         k = config["kmer_length"]
#     threads:
#         config["pgkmc_count_threads"]
#     log:
#         "logs/snps/genomes/{chr}/count_non_centromeric_snpmers_{hap}_unfilt.log"
#     shell:
#         """
#         pgkmc count -v -w \
#             --snp \
#             -k {params.k} \
#             -t {threads} \
#             -b {input.bed} \
#             --kmers {input.snpmers} \
#             {input.fa} \
#             2> {log} | \
#             gzip -c > {output.cts}
#         """

# rule filterge_counts_no_cm:
#     input:
#         tsv_list = expand(
#             os.path.join(WORKDIR, "results/snps/genomes/{chr}/{hap}_unfilt.non_centromeric.snpmers.tsv.gz"),
#             chr="{chr}",
#             hap=hap_list
#         )
#     output:
#         joined_tsv = os.path.join(WORKDIR, "results/snps/genomes/{chr}/unfilt.non_centromeric.snpmers.tsv")
#     log:
#         "logs/snps/genomes/{chr}/filterge_non_centromeric_counts_unfilt.log"
#     script:
#         "../scripts/py/filterge_cts.py"

# rule compute_snpmers_occurrences_no_cm:
#     input:
#         joined_tsvs = expand(
#             os.path.join(WORKDIR, "results/snps/genomes/{chr}/unfilt.non_centromeric.snpmers.tsv"),
#             chr=config["chroms"]
#         )
#     output:
#         snpmers_cts = os.path.join(WORKDIR, "results/snps/genomes/chrALL/unfilt.non_centromeric.snpmers_occ.tsv")
#     log:
#         "logs/snps/genomes/chrALL/compute_non_centromeric_snpmers_occurrences.log"
#     script:
#         "../scripts/py/count_snpmers_occurrences.py"

# rule count_snpmers_cm:
#     input:
#         fa = os.path.join(WORKDIR, "data/asm/{hap}.fa.gz"),
#         bed = os.path.join(WORKDIR, "results/cenSat/{chr}/{hap}.centromeric.sf.bed"),
#         snpmers = os.path.join(WORKDIR, "results/snps/genomes/chrALL/hmn462_unfilt.genomic.snpmers.txt")
#     output:
#         cts = os.path.join(WORKDIR, "results/snps/centromeres/{chr}/{hap}_unfilt.centromeric.snpmers.tsv.gz")
#     params:
#         k = config["kmer_length"]
#     threads:
#         config["pgkmc_count_threads"]
#     log:
#         "logs/snps/centromeres/{chr}/count_centromeric_snpmers_{hap}_unfilt.log"
#     shell:
#         """
#         pgkmc count -v -w \
#             --snp \
#             -k {params.k} \
#             -t {threads} \
#             -b {input.bed} \
#             --kmers {input.snpmers} \
#             {input.fa} \
#             2> {log} | \
#             gzip -c > {output.cts}
#         """

# rule filterge_counts_cm:
#     input:
#         tsv_list = expand(
#             os.path.join(WORKDIR, "results/snps/centromeres/{chr}/{hap}_unfilt.centromeric.snpmers.tsv.gz"),
#             chr="{chr}",
#             hap=hap_list
#         )
#     output:
#         joined_tsv = os.path.join(WORKDIR, "results/snps/centromeres/{chr}/unfilt.centromeric.snpmers.tsv")
#     log:
#         "logs/snps/centromeres/{chr}/filterge_centromeric_counts_unfilt.log"
#     script:
#         "../scripts/py/filterge_cts.py"

# rule compute_snpmers_occurrences_cm:
#     input:
#         joined_tsvs = expand(
#             os.path.join(WORKDIR, "results/snps/centromeres/{chr}/unfilt.centromeric.snpmers.tsv"),
#             chr=config["chroms"]
#         )
#     output:
#         snpmers_cts = os.path.join(WORKDIR, "results/snps/centromeres/chrALL/unfilt.centromeric.snpmers_occ.tsv")
#     log:
#         "logs/snps/centromeres/chrALL/compute_centromeric_snpmers_occurrences.log"
#     script:
#         "../scripts/py/count_snpmers_occurrences.py"

# rule merge_acrocentric_cm:
#     input:
#         ac_tsvs = lambda wc: expand(
#             os.path.join(WORKDIR, "results/snps/centromeres/{chr}/unfilt.centromeric.snpmers.tsv"),
#             chr=sorted(ACRO_SINGLE)
#         ),
#         snpmers_cts_cm = os.path.join(WORKDIR, "results/snps/centromeres/chrALL/unfilt.centromeric.snpmers_occ.tsv"),
#         snpmers_cts_no_cm = os.path.join(WORKDIR, "results/snps/genomes/chrALL/unfilt.non_centromeric.snpmers_occ.tsv")
#     wildcard_constraints:
#         pair = ACRO_PAIR_RE
#     output:
#         merged = expand(
#             os.path.join(WORKDIR, "results/snps/centromeres/{pair}/unfilt.centromeric.snpmers.tsv"),
#             pair=ACRO_MERGED,
#         ),
#         snpmers_cts_cm = os.path.join(WORKDIR, "results/snps/centromeres/chrALL/unfilt.centromeric.snpmers_occ.refined.tsv"),
#         snpmers_cts_no_cm = os.path.join(WORKDIR, "results/snps/genomes/chrALL/unfilt.non_centromeric.snpmers_occ.refined.tsv")
#     params:
#         pairs = ACRO_PAIRS
#     log:
#         "logs/snps/centromeres/chrALL/merge_acrocentric.log"
#     script:
#         "../scripts/py/merge_acrocentric.py"

rule evaluate_tree_stability:
    input:
        joined_tsv = os.path.join(WORKDIR, "results/snps/centromeres/{chr}/unfilt.centromeric.snpmers.tsv")
    output:
        boot_tree_plot = os.path.join(WORKDIR, "plots/tree_stability/{chr}/bootstrap_tree.png")
    params:
        n_boot = 1000,
        maf = 0.05,
        max_missing_frac = 1.0
    log:
        "logs/snps/centromeres/{chr}/evaluate_tree_stability.log"
    conda:
        "../envs/py/ml.yaml"
    script:
        "../scripts/py/bootstrap.py"

rule compute_cenhaps:
    input:
        joined_tsv = os.path.join(WORKDIR, "results/snps/centromeres/{chr}/unfilt.centromeric.snpmers.tsv"),
        meta = os.path.join(WORKDIR, "data/metadata/human462.meta.tsv"),
        anno = os.path.join(WORKDIR, "data/cenhaps/assignments/{chr}.cenhap_predictions.tsv"),
        snpmers_cts_cm = os.path.join(WORKDIR, "results/snps/centromeres/chrALL/unfilt.centromeric.snpmers_occ.refined.tsv"),
        snpmers_cts_no_cm = os.path.join(WORKDIR, "results/snps/genomes/chrALL/unfilt.non_centromeric.snpmers_occ.refined.tsv"),
        t2t = os.path.join(WORKDIR, "results/annotations/t2t.tsv")
    # wildcard_constraints:
    #     chr = CENHAP_CHR_RE
    output:
        hm = os.path.join(WORKDIR, "plots/cenhaps/hm/{chr}_hm.png"),
        hm_sel = os.path.join(WORKDIR, "plots/cenhaps/hm_selected/{chr}_hm_selected.png"),
        gt_labels = os.path.join(WORKDIR, "results/snps/centromeres/{chr}/cm_gt_labels.tsv"),
        gt_matrix = os.path.join(WORKDIR, "results/snps/centromeres/{chr}/cm_gt_matrix.tsv")
    params:
        maf = 0.05,
        chr = "{chr}",
        max_missing_frac = 1.0,
        max_snpmers = 100,
        beam_width = 3,
        tol = 0.005
    log:
        "logs/snps/centromeres/{chr}/compute_cenhaps.log"
    conda:
        "../envs/py/ml.yaml"
    script:
        "../scripts/py/compute_cenhaps.py"

# Union of the SNP-mers used in the matrices of all chromosomes
rule merge_snpmers:
    input:
        gt_matrix = expand(
            os.path.join(WORKDIR, "results/snps/centromeres/{chr}/cm_gt_matrix.tsv"),
            # chr=CENHAP_CHRS
            chr=config["chroms"]
        )
    output:
        snpmers = os.path.join(WORKDIR, "results/snps/centromeres/chrALL/decision_snpmers.txt")
    log:
        "logs/snps/centromeres/chrALL/merge_decision_snpmers.log"
    shell:
        """
        cut -f1 {input.gt_matrix} | sort | uniq | grep -vx 'snpmer' > {output.snpmers}
        """


def find_cram(wildcards):
    matches = glob.glob(f"/1000gen/data/aligned_BAMs/ERR324/*/{wildcards.sample}.final.cram")
    if len(matches) == 0:
        raise ValueError(f"No CRAM found for sample {wildcards.sample}")
    if len(matches) > 1:
        raise ValueError(f"Multiple CRAMs found for sample {wildcards.sample}: {matches}")
    return matches[0]


# Genotype a sample: nearest neighbour of its SNP-mer composition among the reference samples
rule genotype_samples:
    input:
        gt_matrix = expand(
            os.path.join(WORKDIR, "results/snps/centromeres/{chr}/cm_gt_matrix.tsv"),
            # chr=CENHAP_CHRS
            chr=config["chroms"]
        ),
        gt_labels = expand(
            os.path.join(WORKDIR, "results/snps/centromeres/{chr}/cm_gt_labels.tsv"),
            # chr=CENHAP_CHRS
            chr=config["chroms"]
        ),
        snpmers = os.path.join(WORKDIR, "results/snps/centromeres/chrALL/decision_snpmers.txt"),
        ref = "/1000gen/refs/1000genomes_GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        cram = find_cram
    output:
        gtype = os.path.join(WORKDIR, "results/gtype/centromeres/{sample}/genotype.tsv")
    params:
        k = config["kmer_length"],
        matrix_list = os.path.join(WORKDIR, "results/snps/centromeres/chrALL/cm_gt_matrix_paths.txt"),
        labels_list = os.path.join(WORKDIR, "results/snps/centromeres/chrALL/cm_gt_labels_paths.txt")
    threads:
        3
    log:
        "logs/gtype/centromeres/{sample}/genotype.log"
    conda:
        "../envs/standalone/sambedtools.yaml"
    shell:
        """
        printf '%s\\n' {input.gt_matrix} > {params.matrix_list};
        printf '%s\\n' {input.gt_labels} > {params.labels_list};
        samtools fastq --reference {input.ref} {input.cram} | \
        pgkmc clust -k{params.k} -t{threads} -K 0.1g -o {output.gtype} -v \
            --matrices {params.matrix_list} --clusters {params.labels_list} \
            --kmers {input.snpmers} - 2> {log}
        """














# def find_fa(wildcards):
#     d = os.path.join(WORKDIR, "data/asm")
#     s = wildcards.sample
#     for a, b in (("pat", "mat"), ("hap1", "hap2")):
#         fa = glob.glob(f"{d}/*{s}.{a}.fa.gz")
#         fb = glob.glob(f"{d}/*{s}.{b}.fa.gz")
#         if fa and fb:
#             return [fa[0], fb[0]]

# rule genotype_samples:
#     input:
#         gt_profiles = expand(
#             os.path.join(WORKDIR, "results/snps/centromeres/{chr}/cenhaps_gt_profiles.tsv"),
#             chr=config["chroms"]
#         ),
#         snpmers = os.path.join(WORKDIR, "results/snps/centromeres/chrALL/decision_snpmers.txt"),
#         fa_files = find_fa
#     output:
#         gtype = os.path.join(WORKDIR, "results/gtype/centromeres/{sample}/genotype.tsv")
#     params:
#         k = config["kmer_length"],
#         matrix_list = os.path.join(WORKDIR, "results/snps/centromeres/chrALL/cenhaps_gt_profiles_paths.txt"),
#     threads:
#         3
#     log:
#         "logs/gtype/centromeres/{sample}/genotype.log"
#     shell:
#         """
#         printf '%s\\n' {input.gt_profiles} > {params.matrix_list};
#         zcat -f {input.fa_files} | \
#         pgkmc clust -k{params.k} -t{threads} -K 0.1g -o {output.gtype} -v \
#             --matrices {params.matrix_list} --kmers {input.snpmers} - 2> {log}
#         """

# rule plot_snp_density_cm:
#     input:
#         cts = expand(
#             os.path.join(WORKDIR, "results_chrom_spec/snps/centromeres/{chr}/{hap}_unfilt.centromeric.snpmers.tsv"),
#             chr = "{chr}",
#             hap = hap_list
#         ),
#         bed = expand(
#             os.path.join(WORKDIR, "data/annotations/cenSat/{hap}.cenSat.bed"),
#             chr = "{chr}",
#             hap = hap_list
#         ),
#         t2t = os.path.join(WORKDIR, "results/annotations/t2t.tsv")
#     output:
#         density_plot = os.path.join(WORKDIR, "plots_chrom_spec/snpmers_density/{chr}/density.png")
#     params:
#         chr = "{chr}",
#         downsample_freq = 0.001
#     log:
#         "logs/snps/centromeres/{chr}/plot_snpmers_density.log"
#     conda:
#         "../envs/py/ml.yaml"
#     script:
#         "../scripts/py/plot_snpmers_density.py"



# rule count_kmers_cm:
#     input:
#         fa = os.path.join(WORKDIR, "data/asm/{hap}.fa.gz"),
#         bed = os.path.join(WORKDIR, "results/cenSat/{chr}/{hap}.centromeric.sf.bed")
#     output:
#         kmers = os.path.join(WORKDIR, "results/rare_kmers/centromeres/{chr}/{hap}.centromeric.kmers.tsv.gz")
#     params:
#         kmers = os.path.join(WORKDIR, "results/rare_kmers/centromeres/{chr}/{hap}.centromeric.kmers.tsv"),
#         k = config["kmer_length"]
#     threads:
#         config["pgkmc_threads"]
#     log:
#         "logs/rare_kmers/centromeres/{chr}/count_centromeric_kmers_{hap}.log"
#     shell:
#         """
#         pgkmc count -v -w \
#             -k {params.k} \
#             -t {threads} \
#             -b {input.bed} \
#             -o {params.kmers} \
#             {input.fa} \
#             2> {log};
#         gzip {params.kmers}
#         """