# ### Download and merge chromosome aliases
# rule download_chrom_alias:
#     output:
#         chromalias = os.path.join(WORKDIR, "data/annotations/chromalias/{hap}.chromAlias.txt")
#     params:
#         hprc_url = config["hprc_url"]
#     log:
#         "logs/download/chromalias/download_{hap}_chromalias.log"
#     shell:
#         r"""
#         exec > {log} 2>&1
#         S3="{params.hprc_url}"

#         BASE={wildcards.hap}
#         IFS='_' read -r F1 F2 _ <<< "$BASE"
#         IFS='.' read -r S1 S2 _ <<< "$F2"
#         SAMPLE="$S1"
#         SAMPLE_HAP="${{S1}}_${{S2}}"
#         SAMPLE_FULL="${{F1}}_${{S1}}.${{S2}}"

#         PFX="submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/${{SAMPLE}}/assemblies/freeze_2/annotation/chrom_assignment/"
#         KEY=$(curl -fsS "$S3/?list-type=2&prefix=$PFX" \
#             | grep -o '<Key>[^<]*</Key>' | sed 's|</*Key>||g' \
#             | grep "/${{SAMPLE_HAP}}_.*\.chromAlias\.txt$" | head -1)

#         curl -fsSL -o {output.chromalias} "$S3/$KEY"
#         """

# rule create_chromalias_table:
#     input:
#         chromalias = expand(
#             os.path.join(WORKDIR, "data/annotations/chromalias/{hap}.chromAlias.txt"),
#             hap = hap_list
#         )
#     output:
#         chrom_table = os.path.join(WORKDIR, "results/annotations/chromalias_table.txt")
#     log:
#         "logs/annotations/chromalias/create_chromalias_table.txt"
#     shell:
#         """
#         printf 'assembly\\tucsc\\tgenbank\\n' > {output.chrom_table} 2> {log}
#         awk 'FNR>1' {input.chromalias} >> {output.chrom_table} 2>> {log}
#         """

# ### Download centromeric coordinates
# rule download_cenSat_bed:
#     output:
#         bed = os.path.join(WORKDIR, "data/annotations/cenSat/{hap}.cenSat.bed")
#     params:
#         hprc_url = config["hprc_url"]
#     log:
#         "logs/download/cenSat/download_{hap}_cenSat_anno.log"
#     shell:
#         r"""
#         exec > {log} 2>&1
#         S3="{params.hprc_url}"

#         BASE={wildcards.hap}
#         IFS='_' read -r F1 F2 _ <<< "$BASE"
#         IFS='.' read -r S1 S2 _ <<< "$F2"
#         SAMPLE="$S1"
#         SAMPLE_HAP="${{S1}}_${{S2}}"
#         SAMPLE_FULL="${{F1}}_${{S1}}.${{S2}}"

#         PFX="submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/${{SAMPLE}}/assemblies/freeze_2/annotation/censat/"
#         KEY=$(curl -fsS "$S3/?list-type=2&prefix=$PFX" \
#             | grep -o '<Key>[^<]*</Key>' | sed 's|</*Key>||g' \
#             | grep "/${{SAMPLE_HAP}}_.*\.cenSat\.bed$" | head -1)

#         curl -fsSL -o {output.bed} "$S3/$KEY"
#         """

# rule process_cenSat_bed:
#     input:
#         bed = os.path.join(WORKDIR, "data/annotations/cenSat/{hap}.cenSat.bed")
#     output:
#         cm_bed = os.path.join(WORKDIR, "results/cenSat/chrALL/{hap}.centromeric.bed"),
#         pcm_bed = os.path.join(WORKDIR, "results/cenSat/chrALL/{hap}.pericentromeric.bed")
#     log:
#         "logs/annotations/cenSat/process_{hap}_cenSat_bed.log"
#     shell:
#         r"""
#         awk -F'\t' '!/^track|^#|^browser/ {{
#             f4 = tolower($4)
#             if (f4 ~ /^active_hor\(/) print
#         }}' {input.bed} > {output.cm_bed} 2> {log}

#         awk -F'\t' '!/^track|^#|^browser/ {{
#             f4 = tolower($4)
#             if (f4 !~ /^active_hor\(/) print
#         }}' {input.bed} > {output.pcm_bed} 2> {log}
#         """

# rule fix_bed_sf_names:
#     input:
#         bed = os.path.join(WORKDIR, "results/cenSat/chrALL/{hap}.centromeric.bed"),
#         sf = os.path.join(WORKDIR, "data/sfalias/sf_alias_table.txt")
#     output:
#         bed = os.path.join(WORKDIR, "results/cenSat/chrALL/{hap}.centromeric.sf.bed")
#     log:
#         "logs/annotations/centromeres/chrALL/add_{hap}_centromeric_sf_names_bed.txt"
#     script:
#         "../scripts/py/add_bed_sf.py"

# rule split_cenSat_bed_chr_names:
#     input:
#         bed = os.path.join(WORKDIR, "results/cenSat/chrALL/{hap}.centromeric.sf.bed"),
#         chrom_table = os.path.join(WORKDIR, "results/annotations/chromalias_table.txt")
#     output:
#         bed = os.path.join(WORKDIR, "results/cenSat/{chr}/{hap}.centromeric.sf.bed")
#     params:
#         chrom = "{chr}"
#     log:
#         "logs/annotations/centromeres/{chr}/split_{hap}_centromeric_chr_names_bed.txt"
#     script:
#         "../scripts/py/split_bed_chr.py"

# rule create_bed_list_cm_single_chr:
#     input:
#         bed_paths = expand(
#             os.path.join(WORKDIR, "results/cenSat/{chr}/{hap}.centromeric.sf.bed"),
#             chr="{chr}",
#             hap=hap_list
#         )
#     output:
#         bed_list = os.path.join(WORKDIR, "results/annotations/bed_lists/{chr}/bed_list.centromeric.txt")
#     log:
#         "logs/annotations/centromeres/{chr}/create_centromeric_bed_list.txt"
#     shell:
#         """
#         printf '%s\\n' {input.bed_paths} > {output.bed_list} 2> {log}
#         """

# rule create_bed_list_cm_all_chr:
#     input:
#         bed_paths = expand(
#             os.path.join(WORKDIR, "results/cenSat/chrALL/{hap}.centromeric.sf.bed"),
#             chr="{chr}",
#             hap=hap_list
#         )
#     output:
#         bed_list = os.path.join(WORKDIR, "results/annotations/bed_lists/chrALL/bed_list.centromeric.txt")
#     log:
#         "logs/annotations/centromeres/chrALL/create_centromeric_bed_list.txt"
#     shell:
#         """
#         printf '%s\\n' {input.bed_paths} > {output.bed_list} 2> {log}
#         """


# ### Download and merge T2T chromosomes
# rule download_t2t_chr_anno:
#     output:
#         t2t = os.path.join(WORKDIR, "data/annotations/t2t/{hap}.t2t.tsv")
#     params:
#         hprc_url = config["hprc_url"]
#     log:
#         "logs/download/t2t/download_{hap}_t2t.log"
#     shell:
#         r"""
#         exec > {log} 2>&1
#         S3="{params.hprc_url}"

#         BASE={wildcards.hap}
#         IFS='_' read -r F1 F2 _ <<< "$BASE"
#         IFS='.' read -r S1 S2 _ <<< "$F2"
#         SAMPLE="$S1"
#         SAMPLE_HAP="${{S1}}_${{S2}}"
#         SAMPLE_FULL="${{F1}}_${{S1}}.${{S2}}"

#         PFX="submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/${{SAMPLE}}/assemblies/freeze_2/annotation/chrom_assignment/"
#         KEY=$(curl -fsS "$S3/?list-type=2&prefix=$PFX" \
#             | grep -o '<Key>[^<]*</Key>' | sed 's|</*Key>||g' \
#             | grep "/${{SAMPLE_HAP}}_.*\.t2t_chromosomes\.tsv$" | head -1)

#         curl -fsSL -o {output.t2t} "$S3/$KEY"
#         """

# rule create_t2t_table:
#     input:
#         t2t = expand(
#             os.path.join(WORKDIR, "data/annotations/t2t/{hap}.t2t.tsv"),
#             hap = hap_list
#         )
#     output:
#         t2t_table = os.path.join(WORKDIR, "results/annotations/t2t.tsv")
#     params:
#         hprc_url = config["hprc_url"]
#     log:
#         "logs/annotations/t2t/create_t2t_table.txt"
#     shell:
#         r"""
#         printf 'hap\tsequence_id\tchr_name\tlevel\tsequence_length\tnum_gaps\ttotal_gap_length\n' > {output.t2t_table} 2> {log}
#         for f in {input.t2t}; do
#             hap=$(basename "$f" .t2t.tsv)
#             awk -F'\t' -v OFS='\t' -v h="$hap" 'FNR>1 {{print h, $0}}' "$f" >> {output.t2t_table} 2>> {log}
#         done
#         """

rule create_full_chr_bed:
    input:
        fa = os.path.join(WORKDIR, "data/asm/{hap}.fa.gz")
    output:
        bed = os.path.join(WORKDIR, "results/chrBed/chrALL/{hap}.bed")
    log:
        "logs/chrBed/chrALL/get_{hap}_chr_full_bed.log"
    conda:
        "../envs/standalone/sambedtools.yaml"
    shell:
        r"""
        samtools faidx {input.fa} 2> {log};
        awk -v OFS='\t' '{{print $1, 0, $2}}' {input.fa}.fai > {output.bed} 2> {log}
        """

rule split_full_chr_bed_chr_names:
    input:
        bed = os.path.join(WORKDIR, "results/chrBed/chrALL/{hap}.bed"),
        chrom_table = os.path.join(WORKDIR, "results/annotations/chromalias_table.txt")
    output:
        bed = os.path.join(WORKDIR, "results/chrBed/{chr}/{hap}.bed")
    params:
        chrom = "{chr}"
    log:
        "logs/chrBed/{chr}/split_{hap}_chr_full_bed.log"
    script:
        "../scripts/py/split_bed_chr.py"

rule subtract_cm2full_chr_bed_chr_names:
    input:
        full_bed = os.path.join(WORKDIR, "results/chrBed/{chr}/{hap}.bed"),
        cm_bed = os.path.join(WORKDIR, "results/cenSat/{chr}/{hap}.centromeric.sf.bed")
    output:
        bed = os.path.join(WORKDIR, "results/chrBed/{chr}/{hap}.no_cm.bed")
    log:
        "logs/chrBed/{chr}/subtract_{hap}_chr_full_bed.log"
    conda:
        "../envs/standalone/sambedtools.yaml"
    shell:
        """
        bedtools subtract -a {input.full_bed} -b {input.cm_bed} > {output.bed}
        """