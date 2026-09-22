### Compute fai idx
rule index_fasta:
    input:
        fa = os.path.join(DATADIR, "{sample}.fa.gz")
    output:
        fai_idx = os.path.join(WORKDIR, "data/asm/fa_idx/{sample}.fa.gz.fai"),
        gzi_idx = os.path.join(WORKDIR, "data/asm/fa_idx/{sample}.fa.gz.gzi")
    log:
        "logs/annotations/asm/index_{sample}_fa.log"
    conda:
        "../envs/standalone/sambedtools.yaml"
    shell:
        r"""
        samtools faidx {input.fa} --fai-idx {output.fai_idx} --gzi-idx {output.gzi_idx} 2> {log}
        """

### Download and merge chromosome aliases
rule download_chrom_alias:
    output:
        chromalias = os.path.join(WORKDIR, "data/annotations/chromalias/{sample}.chromAlias.txt")
    params:
        hprc_url = config["hprc_url"]
    log:
        "logs/download/chromalias/download_{sample}_chromalias.log"
    shell:
        r"""
        exec > {log} 2>&1
        S3="{params.hprc_url}"

        BASE={wildcards.sample}
        IFS='_' read -r F1 F2 _ <<< "$BASE"
        IFS='.' read -r S1 S2 _ <<< "$F2"
        SAMPLE="$S1"
        SAMPLE_HAP="${{S1}}_${{S2}}"
        SAMPLE_FULL="${{F1}}_${{S1}}.${{S2}}"

        PFX="submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/${{SAMPLE}}/assemblies/freeze_2/annotation/chrom_assignment/"
        KEY=$(curl -fsS "$S3/?list-type=2&prefix=$PFX" \
            | grep -o '<Key>[^<]*</Key>' | sed 's|</*Key>||g' \
            | grep "/${{SAMPLE_HAP}}_.*\.chromAlias\.txt$" | head -1)

        curl -fsSL -o {output.chromalias} "$S3/$KEY"
        """

rule create_chromalias_table:
    input:
        chromalias = expand(
            os.path.join(WORKDIR, "data/annotations/chromalias/{sample}.chromAlias.txt"),
            sample = sample_list
        )
    output:
        chrom_table = os.path.join(WORKDIR, "results/annotations/chromalias_table.txt")
    log:
        "logs/annotations/chromalias/create_chromalias_table.txt"
    shell:
        """
        printf 'assembly\\tucsc\\tgenbank\\n' > {output.chrom_table} 2> {log}
        awk 'FNR>1' {input.chromalias} >> {output.chrom_table} 2>> {log}
        """

### Download centromeric coordinates
rule download_cenSat_bed:
    output:
        bed = os.path.join(WORKDIR, "data/annotations/cenSat/{sample}.cenSat.bed")
    params:
        hprc_url = config["hprc_url"]
    log:
        "logs/download/cenSat/download_{sample}_cenSat_anno.log"
    shell:
        r"""
        exec > {log} 2>&1
        S3="{params.hprc_url}"

        BASE={wildcards.sample}
        IFS='_' read -r F1 F2 _ <<< "$BASE"
        IFS='.' read -r S1 S2 _ <<< "$F2"
        SAMPLE="$S1"
        SAMPLE_HAP="${{S1}}_${{S2}}"
        SAMPLE_FULL="${{F1}}_${{S1}}.${{S2}}"

        PFX="submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/${{SAMPLE}}/assemblies/freeze_2/annotation/censat/"
        KEY=$(curl -fsS "$S3/?list-type=2&prefix=$PFX" \
            | grep -o '<Key>[^<]*</Key>' | sed 's|</*Key>||g' \
            | grep "/${{SAMPLE_HAP}}_.*\.cenSat\.bed$" | head -1)

        curl -fsSL -o {output.bed} "$S3/$KEY"
        """

rule process_cenSat_bed:
    input:
        bed = os.path.join(WORKDIR, "data/annotations/cenSat/{sample}.cenSat.bed")
    output:
        cm_bed = os.path.join(WORKDIR, "results/cenSat/chrALL/{sample}.centromeric.bed"),
        pcm_bed = os.path.join(WORKDIR, "results/cenSat/chrALL/{sample}.pericentromeric.bed")
    log:
        "logs/annotations/cenSat/process_{sample}_cenSat_bed.log"
    shell:
        r"""
        awk -F'\t' '!/^track|^#|^browser/ {{
            f4 = tolower($4)
            if (f4 ~ /^active_hor\(/) print
        }}' {input.bed} > {output.cm_bed} 2> {log}

        awk -F'\t' '!/^track|^#|^browser/ {{
            f4 = tolower($4)
            if (f4 !~ /^active_hor\(/) print
        }}' {input.bed} > {output.pcm_bed} 2> {log}
        """

rule fix_bed_sf_names:
    input:
        bed = os.path.join(WORKDIR, "results/cenSat/chrALL/{sample}.{region}.bed"),
        sf = os.path.join(WORKDIR, "data/sfalias/sf_alias_table.txt")
    output:
        bed = os.path.join(WORKDIR, "results/cenSat/chrALL/{sample}.{region}.sf.bed")
    log:
        "logs/annotations/centromeres/chrALL/add_{sample}_{region}_sf_names_bed.txt"
    script:
        "../scripts/py/add_bed_sf.py"

rule split_cenSat_bed_chr_names:
    input:
        bed = os.path.join(WORKDIR, "results/cenSat/chrALL/{sample}.{region}.sf.bed"),
        chrom_table = os.path.join(WORKDIR, "results/annotations/chromalias_table.txt"),
        # mashmap = os.path.join(WORKDIR, "data/annotations/mashmap/{sample}.pi95.paf"),
        # fai_idx = os.path.join(WORKDIR, "data/asm/fa_idx/{sample}.fa.gz.fai")
    output:
        bed = os.path.join(WORKDIR, "results/cenSat/{chr}/{sample}.{region}.sf.bed")
    params:
        chrom = "{chr}"
    log:
        "logs/annotations/centromeres/{chr}/split_{sample}_{region}_chr_names_bed.txt"
    script:
        "../scripts/py/split_bed_chr.py"

rule create_bed_list_cm:
    input:
        bed_paths = expand(
            os.path.join(WORKDIR, "results/cenSat/{chr}/{sample}.{region}.sf.bed"),
            chr="{chr}",
            sample=sample_list,
            region="{region}"
        )
    output:
        bed_list = os.path.join(WORKDIR, "results/annotations/bed_lists/{chr}/bed_list.{region}.txt")
    log:
        "logs/annotations/centromeres/{chr}/create_{region}_bed_list.txt"
    shell:
        """
        printf '%s\\n' {input.bed_paths} > {output.bed_list} 2> {log}
        """