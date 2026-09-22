import glob

# Count SNP-mers and rare k-mers in centromeres
rule detect_snpmers_cm:
    input:
        fa_files = expand(
            os.path.join(DATADIR, "{hap}.fa.gz"),
            hap=hap_list
        ),
        bed_list = os.path.join(WORKDIR, "results/annotations/bed_lists/{chr}/bed_list.centromeric.txt")
    output:
        snpmers = os.path.join(WORKDIR, "results/snps/centromeres/{chr}/hmn462_unfilt.centromeric.snpmers.txt")
    params:
        k = config["kmer_length"],
        msf = 0,
        f = 0
    threads:
        config["pgkmc_threads"]
    log:
        "logs/snps/centromeres/{chr}/detect_centromeric_snpmers_hmn462_unfilt.log"
    shell:
        """
        pgkmc detect -v \
            --snp \
            -k {params.k} \
            --msf {params.msf} \
            -t {threads} \
            -f {params.f} \
            -b {input.bed_list} \
            -o {output.snpmers} \
            {input.fa_files} \
            2> {log}
        """

rule count_snpmers_cm:
    input:
        fa = os.path.join(DATADIR, "{hap}.fa.gz"),
        bed = os.path.join(WORKDIR, "results/cenSat/{chr}/{hap}.centromeric.sf.bed"),
        snpmers = os.path.join(WORKDIR, "results/snps/centromeres/{chr}/hmn462_unfilt.centromeric.snpmers.txt")
    output:
        cts = os.path.join(WORKDIR, "results/snps/centromeres/{chr}/{hap}_unfilt.centromeric.snpmers.tsv")
    params:
        k = config["kmer_length"]
    threads:
        config["pgkmc_threads"]
    log:
        "logs/snps/centromeres/{chr}/count_centromeric_snpmers_{hap}_unfilt.log"
    shell:
        """
        pgkmc count -v -w \
            --snp \
            -k {params.k} \
            -t {threads} \
            -b {input.bed} \
            --kmers {input.snpmers} \
            -o {output.cts} \
            {input.fa} \
            2> {log}
        """

rule merge_counts_cm:
    input:
        tsv_list = expand(
            os.path.join(WORKDIR, "results/snps/centromeres/{chr}/{hap}_unfilt.centromeric.snpmers.tsv"),
            chr="{chr}",
            hap=hap_list
        )
    output:
        joined_tsv = os.path.join(WORKDIR, "results/snps/centromeres/{chr}/unfilt.centromeric.snpmers.tsv")
    log:
        "logs/snps/centromeres/{chr}/merge_centromeric_counts_unfilt.log"
    script:
        "../scripts/py/merge_cts.py"

# rule compute_cenhaps:
#     input:
#         joined_tsv = os.path.join(WORKDIR, "results/snps/centromeres/{chr}/unfilt.centromeric.snpmers.tsv.gz"),
#         meta = os.path.join(WORKDIR, "data/metadata/human462.meta.tsv"),
#         anno = os.path.join(WORKDIR, "data/cenhaps/assignments/{chr}.cenhap_predictions.tsv")
#     output:
#         hm = os.path.join(WORKDIR, "plots/cenhaps/hm/{chr}_hm.png"),
#         hm_walk = os.path.join(WORKDIR, "plots/cenhaps/hm_walk/{chr}_hm_walk.png"),
#         cenhaps = os.path.join(WORKDIR, "results/snps/centromeres/{chr}/cenhaps_clusters.tsv"),
#         nmi = os.path.join(WORKDIR, "results/snps/centromeres/{chr}/cenhaps_nmi.tsv"),
#         snpmers = os.path.join(WORKDIR, "results/snps/centromeres/{chr}/cenhaps_snpmers.tsv"),
#         tree_structure = os.path.join(WORKDIR, "results/snps/centromeres/{chr}/cenhaps_tree_structure.tsv")
#     params:
#         maf = 0.05,
#         chr = "{chr}"
#     log:
#         "logs/snps/centromeres/{chr}/compute_cenhaps.log"
#     conda:
#         "../envs/py/ml.yaml"
#     script:
#         "../scripts/py/compute_cenhaps.py"


# rule count_kmers_cm:
#     input:
#         fa = os.path.join(DATADIR, "{hap}.fa.gz"),
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

rule merge_snpmers:
    input:
        snpmers = expand(
            os.path.join(WORKDIR, "results/snps/centromeres/{chr}/cenhaps_snpmers.tsv"),
            chr=config["chroms"]
        ),
        trees = expand(
            os.path.join(WORKDIR, "results/snps/centromeres/{chr}/cenhaps_tree_structure.tsv"),
            chr=config["chroms"]
        )
    output:
        snpmers = os.path.join(WORKDIR, "results/snps/centromeres/chr_all/cenhaps_snpmers.tsv"),
        tree = os.path.join(WORKDIR, "results/snps/centromeres/chr_all/cenhaps_tree_structure.tsv")
    log:
        "logs/snps/centromeres/chr_all/merge_snpmers.log"
    script:
        "../scripts/py/merge_gt_snpmers.py"

def find_cram(wildcards):
    matches = glob.glob(f"/1000gen/data/aligned_BAMs/ERR324/*/{wildcards.sample}.final.cram")
    if len(matches) == 0:
        raise ValueError(f"No CRAM found for sample {wildcards.sample}")
    if len(matches) > 1:
        raise ValueError(f"Multiple CRAMs found for sample {wildcards.sample}: {matches}")
    return matches[0]

rule genotype_samples:
    input:
        snpmer_list = os.path.join(WORKDIR, "results/snps/centromeres/chr_all/cenhaps_snpmers.tsv"),
        snpmer_tree = os.path.join(WORKDIR, "results/snps/centromeres/chr_all/cenhaps_tree_structure.tsv"),
        ref = "/1000gen/refs/1000genomes_GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        cram = find_cram
    output:
        gtype = os.path.join(WORKDIR, "results/gtype/centromeres/{sample}/genotype.tsv")
    params:
        k = config["kmer_length"]
    threads:
        3
    log:
        "logs/gtype/centromeres/{sample}/genotype.log"
    conda:
        "../envs/standalone/sambedtools.yaml"
    shell:
        """
        samtools fastq --reference {input.ref} {input.cram} | \
        pgkmc clust -k{params.k} -t{threads} -K 0.1g -o {output.gtype} -v --kmers {input.snpmer_list} --tree {input.snpmer_tree} - 2> {log}
        """