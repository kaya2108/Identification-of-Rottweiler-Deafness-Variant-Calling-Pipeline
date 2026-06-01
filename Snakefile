"""
Rottweiler WES/WGS Pipeline
Reproduces: Hytonen et al. Human Genetics 2021
DOI: 10.1007/s00439-021-02286-z

WES: bwa mem + GATK MarkDuplicates + HaplotypeCaller
WGS: bwa mem + SAMBLASTER + Sambamba + HaplotypeCaller
"""

configfile: "config.yaml"

WES = config["samples"]["WES"]  # SRR13743383
WGS = config["samples"]["WGS"]  # SRR13743384
REF     = config["ref"]
DBSNP   = config["dbsnp"]
FASTQ   = config["fastq_dir"]
OUTDIR  = config["outdir"]
VEP_EXE = config["vep_exe"]
VEP_PERL = config["vep_perl"]

rule all:
    input:
        # QC
        f"{OUTDIR}/qc/{WES}_1_fastqc.html",
        f"{OUTDIR}/qc/{WES}_2_fastqc.html",
        f"{OUTDIR}/qc/{WGS}_1_fastqc.html",
        f"{OUTDIR}/qc/{WGS}_2_fastqc.html",
        f"{OUTDIR}/qc/multiqc_report.html",
        f"{OUTDIR}/qc/trimmed/{WES}_1P_fastqc.html",
        f"{OUTDIR}/qc/trimmed/{WGS}_1P_fastqc.html",
        # Final annotated VCF
        f"{OUTDIR}/final/joint.filtered.vep.vcf",
        # Candidate variants (homozygous missense)
        f"{OUTDIR}/final/candidate_variants.txt",


# SHARED STEPS (same for WES and WGS)


# FastWC
rule fastqc:
    input:
        r1 = f"{FASTQ}/{{sample}}_1.fastq.gz",
        r2 = f"{FASTQ}/{{sample}}_2.fastq.gz",
    output:
        html1 = f"{OUTDIR}/qc/{{sample}}_1_fastqc.html",
        html2 = f"{OUTDIR}/qc/{{sample}}_2_fastqc.html",
        zip1  = f"{OUTDIR}/qc/{{sample}}_1_fastqc.zip",
        zip2  = f"{OUTDIR}/qc/{{sample}}_2_fastqc.zip",
    threads: 2
    conda: "envs/qc.yaml"
    shell:
        """
        fastqc {input.r1} {input.r2} \
            --outdir {OUTDIR}/qc \
            --threads {threads}
        """

# MultiQC
rule multiqc:
    input:
        f"{OUTDIR}/qc/{WES}_1_fastqc.zip",
        f"{OUTDIR}/qc/{WES}_2_fastqc.zip",
        f"{OUTDIR}/qc/{WGS}_1_fastqc.zip",
        f"{OUTDIR}/qc/{WGS}_2_fastqc.zip",
    output:
        f"{OUTDIR}/qc/multiqc_report.html"
    conda: "envs/qc.yaml"
    shell:
        "multiqc {OUTDIR}/qc -o {OUTDIR}/qc"

# Trimmomatic 
rule trim:
    input:
        r1 = f"{FASTQ}/{{sample}}_1.fastq.gz",
        r2 = f"{FASTQ}/{{sample}}_2.fastq.gz",
    output:
        r1p = f"{OUTDIR}/trimmed/{{sample}}_1P.fastq.gz",
        r2p = f"{OUTDIR}/trimmed/{{sample}}_2P.fastq.gz",
        r1u = f"{OUTDIR}/trimmed/{{sample}}_1U.fastq.gz",
        r2u = f"{OUTDIR}/trimmed/{{sample}}_2U.fastq.gz",
        log = f"{OUTDIR}/trimmed/{{sample}}_trim.log",
    threads: 4
    conda: "envs/trim.yaml"
    shell:
        """
        trimmomatic PE -threads {threads} \
            {input.r1} {input.r2} \
            {output.r1p} {output.r1u} \
            {output.r2p} {output.r2u} \
            ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
            LEADING:3 TRAILING:3 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36 \
            2> {output.log}
        """
# FastQC on trimmed reads 
rule fastqc_trimmed:
    input:
        r1 = f"{OUTDIR}/trimmed/{{sample}}_1P.fastq.gz",
        r2 = f"{OUTDIR}/trimmed/{{sample}}_2P.fastq.gz",
    output:
        html1 = f"{OUTDIR}/qc/trimmed/{{sample}}_1P_fastqc.html",
        html2 = f"{OUTDIR}/qc/trimmed/{{sample}}_2P_fastqc.html",
        zip1  = f"{OUTDIR}/qc/trimmed/{{sample}}_1P_fastqc.zip",
        zip2  = f"{OUTDIR}/qc/trimmed/{{sample}}_2P_fastqc.zip",
    threads: 2
    conda: "envs/qc.yaml"
    shell:
        """
        fastqc {input.r1} {input.r2} \
            --outdir {OUTDIR}/qc/trimmed \
            --threads {threads}
        """


# WES BRANCH
# bwa mem + GATK MarkDuplicates (standard GATK best practices)

# 4a. WES Alignment 
rule align_WES:
    input:
        r1  = f"{OUTDIR}/trimmed/{WES}_1P.fastq.gz",
        r2  = f"{OUTDIR}/trimmed/{WES}_2P.fastq.gz",
        ref = REF,
    output:
        bam = f"{OUTDIR}/aligned/{WES}.sorted.bam",
        bai = f"{OUTDIR}/aligned/{WES}.sorted.bam.bai",
    params:
        rg = f"@RG\\tID:{WES}\\tSM:{WES}\\tPL:ILLUMINA\\tLB:{WES}\\tPU:{WES}"
    threads: 8
    conda: "envs/align.yaml"
    shell:
        """
        bwa mem -t {threads} -R '{params.rg}' {input.ref} {input.r1} {input.r2} \
            | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

# WES Mark Duplicates (GATK)
rule markdup_WES:
    input:
        bam = f"{OUTDIR}/aligned/{WES}.sorted.bam",
    output:
        bam     = f"{OUTDIR}/markdup/{WES}.markdup.bam",
        bai     = f"{OUTDIR}/markdup/{WES}.markdup.bai",
        metrics = f"{OUTDIR}/markdup/{WES}.markdup.metrics",
    conda: "envs/gatk.yaml"
    shell:
        """
        gatk MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --CREATE_INDEX true
        """

# WES BQSR 
rule bqsr_WES:
    input:
        bam   = f"{OUTDIR}/markdup/{WES}.markdup.bam",
        ref   = REF,
        dbsnp = DBSNP,
    output:
        table = f"{OUTDIR}/bqsr/{WES}.recal.table",
    conda: "envs/gatk.yaml"
    shell:
        """
        gatk BaseRecalibrator \
            -I {input.bam} \
            -R {input.ref} \
            --known-sites {input.dbsnp} \
            -O {output.table}
        """

rule apply_bqsr_WES:
    input:
        bam   = f"{OUTDIR}/markdup/{WES}.markdup.bam",
        table = f"{OUTDIR}/bqsr/{WES}.recal.table",
        ref   = REF,
    output:
        bam = f"{OUTDIR}/bqsr/{WES}.recal.bam",
        bai = f"{OUTDIR}/bqsr/{WES}.recal.bai",
    conda: "envs/gatk.yaml"
    shell:
        """
        gatk ApplyBQSR \
            -I {input.bam} \
            -R {input.ref} \
            --bqsr-recal-file {input.table} \
            -O {output.bam}
        """

# WES HaplotypeCaller
rule haplotypecaller_WES:
    input:
        bam = f"{OUTDIR}/bqsr/{WES}.recal.bam",
        ref = REF,
    output:
        gvcf = f"{OUTDIR}/gvcf/{WES}.g.vcf.gz",
    threads: 4
    conda: "envs/gatk.yaml"
    shell:
        """
        gatk HaplotypeCaller \
            -I {input.bam} \
            -R {input.ref} \
            -ERC GVCF \
            -O {output.gvcf}
        """


# WGS BRANCH
# bwa mem + SAMBLASTER (dedup) 
# Paper: SpeedSeq = bwa + SAMBLASTER 

# WGS Alignment + SAMBLASTER dedup 
rule align_WGS:
    input:
        r1  = f"{OUTDIR}/trimmed/{WGS}_1P.fastq.gz",
        r2  = f"{OUTDIR}/trimmed/{WGS}_2P.fastq.gz",
        ref = REF,
    output:
        bam = f"{OUTDIR}/aligned/{WGS}.sorted.bam",
        bai = f"{OUTDIR}/aligned/{WGS}.sorted.bam.bai",
    params:
        rg = f"@RG\\tID:{WGS}\\tSM:{WGS}\\tPL:ILLUMINA\\tLB:{WGS}\\tPU:{WGS}"
    threads: 2
    conda: "envs/align_wgs.yaml"
    shell:
        """
        bwa mem -t {threads} -R '{params.rg}' {input.ref} {input.r1} {input.r2} \
            | samblaster \
            | samtools view -bS - \
            | samtools sort -@ {threads} -o {output.bam}

        samtools index {output.bam}
        """

# WGS BQSR 
rule bqsr_WGS:
    input:
        bam   = f"{OUTDIR}/aligned/{WGS}.sorted.bam",
        ref   = REF,
        dbsnp = DBSNP,
    output:
        table = f"{OUTDIR}/bqsr/{WGS}.recal.table",
    conda: "envs/gatk.yaml"
    shell:
        """
        gatk BaseRecalibrator \
            -I {input.bam} \
            -R {input.ref} \
            --known-sites {input.dbsnp} \
            -O {output.table}
        """

rule apply_bqsr_WGS:
    input:
        bam   = f"{OUTDIR}/aligned/{WGS}.sorted.bam",
        table = f"{OUTDIR}/bqsr/{WGS}.recal.table",
        ref   = REF,
    output:
        bam = f"{OUTDIR}/bqsr/{WGS}.recal.bam",
        bai = f"{OUTDIR}/bqsr/{WGS}.recal.bai",
    conda: "envs/gatk.yaml"
    shell:
        """
        gatk ApplyBQSR \
            -I {input.bam} \
            -R {input.ref} \
            --bqsr-recal-file {input.table} \
            -O {output.bam}
        """

# WGS HaplotypeCaller
rule haplotypecaller_WGS:
    input:
        bam = f"{OUTDIR}/bqsr/{WGS}.recal.bam",
        ref = REF,
    output:
        gvcf = f"{OUTDIR}/gvcf/{WGS}.g.vcf.gz",
    threads: 4
    conda: "envs/gatk.yaml"
    shell:
        """
        gatk HaplotypeCaller \
            -I {input.bam} \
            -R {input.ref} \
            -ERC GVCF \
            -O {output.gvcf}
        """


# JOINT GENOTYPING (WES + WGS combined)

rule combine_gvcfs:
    input:
        wes = f"{OUTDIR}/gvcf/{WES}.g.vcf.gz",
        wgs = f"{OUTDIR}/gvcf/{WGS}.g.vcf.gz",
        ref = REF,
    output:
        f"{OUTDIR}/gvcf/combined.g.vcf.gz"
    conda: "envs/gatk.yaml"
    shell:
        """
        gatk CombineGVCFs \
            -R {input.ref} \
            -V {input.wes} \
            -V {input.wgs} \
            -O {output}
        """

rule genotype_gvcfs:
    input:
        gvcf  = f"{OUTDIR}/gvcf/combined.g.vcf.gz",
        ref   = REF,
        dbsnp = DBSNP,
    output:
        f"{OUTDIR}/genotyped/joint.vcf.gz"
    conda: "envs/gatk.yaml"
    shell:
        """
        gatk GenotypeGVCFs \
            -R {input.ref} \
            -V {input.gvcf} \
            --dbsnp {input.dbsnp} \
            -O {output}
        """


# HARD FILTERING

rule hard_filter:
    input:
        vcf = f"{OUTDIR}/genotyped/joint.vcf.gz",
        ref = REF,
    output:
        f"{OUTDIR}/filtered/joint.filtered.vcf.gz"
    conda: "envs/gatk.yaml"
    shell:
        """
        # SNPs
        gatk SelectVariants -R {input.ref} -V {input.vcf} \
            --select-type-to-include SNP \
            -O {OUTDIR}/filtered/snps.vcf.gz

        gatk VariantFiltration -R {input.ref} \
            -V {OUTDIR}/filtered/snps.vcf.gz \
            --filter-expression "QD < 2.0"              --filter-name "QD2" \
            --filter-expression "FS > 60.0"             --filter-name "FS60" \
            --filter-expression "MQ < 40.0"             --filter-name "MQ40" \
            --filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \
            --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            -O {OUTDIR}/filtered/snps.filtered.vcf.gz

        # Indels
        gatk SelectVariants -R {input.ref} -V {input.vcf} \
            --select-type-to-include INDEL \
            -O {OUTDIR}/filtered/indels.vcf.gz

        gatk VariantFiltration -R {input.ref} \
            -V {OUTDIR}/filtered/indels.vcf.gz \
            --filter-expression "QD < 2.0"               --filter-name "QD2" \
            --filter-expression "FS > 200.0"             --filter-name "FS200" \
            --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
            -O {OUTDIR}/filtered/indels.filtered.vcf.gz

        # Merge
        gatk MergeVcfs \
            -I {OUTDIR}/filtered/snps.filtered.vcf.gz \
            -I {OUTDIR}/filtered/indels.filtered.vcf.gz \
            -O {output}
        """

rule remove_chr_prefix:
    input:
        f"{OUTDIR}/filtered/joint.filtered.vcf.gz"
    output:
        f"{OUTDIR}/filtered/joint.filtered.nochr.vcf"
    shell:
        """
        gunzip -c {input} | sed 's/^chr//' > {output}
        """


# VEP ANNOTATION
# Note: VEP uses external x86 executable via config["vep_exe"] and config["vep_perl"]

rule vep:
    input:
        vcf = f"{OUTDIR}/filtered/joint.filtered.nochr.vcf"
    output:
        vcf   = f"{OUTDIR}/final/joint.filtered.vep.vcf",
        stats = f"{OUTDIR}/final/joint.filtered.vep.stats.html",
    params:
        cache_dir = config["vep_cache"],
    shell:
        """
        {VEP_PERL} {VEP_EXE} \
            --input_file {input.vcf} \
            --output_file {output.vcf} \
            --stats_file {output.stats} \
            --format vcf \
            --vcf \
            --species canis_lupus_familiaris \
            --assembly CanFam3.1 \
            --cache \
            --dir_cache {params.cache_dir} \
            --offline \
            --symbol \
            --numbers \
            --biotype \
            --canonical \
            --sift b \
            --force_overwrite \
            --cache_version 104
        """


# CANDIDATE VARIANT FILTERING
# Filter for homozygous missense variants using python script (candidate_variants.py)

rule candidate_variants:
    input:
        f"{OUTDIR}/final/joint.filtered.vep.vcf"
    output:
        csv = f"{OUTDIR}/final/candidate_variants.csv",
        txt = f"{OUTDIR}/final/candidate_variants.txt"

    script:
        "scripts/candidate_variants.py"
