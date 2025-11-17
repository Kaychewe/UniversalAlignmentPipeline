###############################################
# Snakefile — Full Workflow: FASTQ → QC BAM
###############################################
# TODO publish version 1 of pipeline 
import os

configfile: "config/config.yaml"

# Load config
SAMPLES = config["samples"]
OUT_ROOT = config["output_root"]
LOG_ROOT = config.get("log_root", "logs")

SKIP_SUBSAMPLE = config.get("skip_subsample", False)
SUBSAMPLE_FRACTION = float(config.get("subsample_fraction", 0.01))

SKIP_TRIMMING = config.get("skip_trimming", False)
THREADS = int(config.get("threads", 8))

os.makedirs(LOG_ROOT, exist_ok=True)

###############################################
# FINAL TARGET
###############################################

rule all:
    input:
        #
        # Trimmed FASTQ
        #
        expand(f"{OUT_ROOT}/{{sample}}/trimmed/{{sample}}_R1_clean.fastq.gz",
               sample=SAMPLES.keys()),
        expand(f"{OUT_ROOT}/{{sample}}/trimmed/{{sample}}_R2_clean.fastq.gz",
               sample=SAMPLES.keys()),

        #
        # Alignment (unsorted BAM)
        #
        expand(f"{OUT_ROOT}/{{sample}}/aligned/{{sample}}.unsorted.bam",
               sample=SAMPLES.keys()),

        #
        # Sort + index
        #
        expand(f"{OUT_ROOT}/{{sample}}/sorted/{{sample}}.sorted.bam",
               sample=SAMPLES.keys()),
        expand(f"{OUT_ROOT}/{{sample}}/sorted/{{sample}}.sorted.bam.bai",
               sample=SAMPLES.keys()),

        #
        # Dedup + metrics
        #
        expand(f"{OUT_ROOT}/{{sample}}/dedup/{{sample}}.dedup.bam",
               sample=SAMPLES.keys()),
        # expand(f"{OUT_ROOT}/{{sample}}/logs/markdup.metrics.txt",
        #        sample=SAMPLES.keys()),

        #
        # Readgroup-fixed BAM
        #
        expand(f"{OUT_ROOT}/{{sample}}/rgfixed/{{sample}}.rg.bam",
               sample=SAMPLES.keys()),

        #
        # QC outputs
        #
        expand(f"{OUT_ROOT}/{{sample}}/qc/{{sample}}.stats.txt",
               sample=SAMPLES.keys()),
        expand(f"{OUT_ROOT}/{{sample}}/qc/{{sample}}.flagstat.txt",
               sample=SAMPLES.keys()),
        expand(f"{OUT_ROOT}/{{sample}}/qc/{{sample}}.idxstats.txt",
               sample=SAMPLES.keys()),
        expand(f"{OUT_ROOT}/{{sample}}/qc/{{sample}}.depth.txt",
               sample=SAMPLES.keys()),
        expand(f"{OUT_ROOT}/{{sample}}/qc/{{sample}}.depth_summary.txt",
               sample=SAMPLES.keys())

############################################################
# RULE 1: Subsample FASTQ
############################################################
# TODO desired behavior. if skip subsampling use the full files via symbolic link 
# rule subsample_fastq:
#     input:
#         r1 = lambda wc: SAMPLES[wc.sample]["fastq_R1"],
#         r2 = lambda wc: SAMPLES[wc.sample]["fastq_R2"]
#     output:
#         r1_out = f"{OUT_ROOT}/{{sample}}/subsampled/{{sample}}_R1_sub.fastq.gz",
#         r2_out = f"{OUT_ROOT}/{{sample}}/subsampled/{{sample}}_R2_sub.fastq.gz"
#     params:
#         script = "scripts/subsample_fastq.sh",
#         frac = lambda wc: SAMPLES[wc.sample].get("subsample_fraction", SUBSAMPLE_FRACTION)
#     log:
#         f"{OUT_ROOT}/{{sample}}/logs/{{sample}}.subsample.log"
#     message:
#         "Subsampling {wildcards.sample}"

#     run:
#         outdir = f"{OUT_ROOT}/{wildcards.sample}/subsampled"
#         os.makedirs(outdir, exist_ok=True)

#         has_r2 = input.r2 is not None

#         if SKIP_SUBSAMPLE:
#             # Handle SE or PE
#             shell("ln -sf {input.r1} {output.r1_out}")
#             if has_r2:
#                 shell("ln -sf {input.r2} {output.r2_out}")
#         else:
#             if has_r2:
#                 shell("bash {params.script} {params.frac} {outdir} {input.r1} {input.r2} {wildcards.sample} &>> {log}")
#             else:
#                 shell("bash {params.script} {params.frac} {outdir} {input.r1} '' {wildcards.sample} &>> {log}")


############################################################
# RULE 2: fastp trimming
############################################################

rule fastp_trim:
    input:
        r1 = f"{OUT_ROOT}/{{sample}}/subsampled/{{sample}}_R1_sub.fastq.gz",
        r2 = f"{OUT_ROOT}/{{sample}}/subsampled/{{sample}}_R2_sub.fastq.gz"
    output:
        r1_clean = f"{OUT_ROOT}/{{sample}}/trimmed/{{sample}}_R1_clean.fastq.gz",
        r2_clean = f"{OUT_ROOT}/{{sample}}/trimmed/{{sample}}_R2_clean.fastq.gz",
        html =     f"{OUT_ROOT}/{{sample}}/trimmed/fastp.html",
        json =     f"{OUT_ROOT}/{{sample}}/trimmed/fastp.json"
    params:
        script = "scripts/02_fastp.sh",
        threads = THREADS
    log:
        f"{OUT_ROOT}/{{sample}}/logs/fastp.log"
    message:
        "fastp trimming {wildcards.sample}"

    run:
        outdir = f"{OUT_ROOT}/{wildcards.sample}/trimmed"
        logdir = f"{OUT_ROOT}/{wildcards.sample}/logs"
        os.makedirs(outdir, exist_ok=True)
        os.makedirs(logdir, exist_ok=True)

        if SKIP_TRIMMING:
            shell("ln -sf {input.r1} {output.r1_clean}")
            shell("ln -sf {input.r2} {output.r2_clean}")
            shell("echo 'skipped' > {output.html}")
            shell("echo 'skipped' > {output.json}")
        else:
            shell("""
                bash {params.script} \
                    {OUT_ROOT} \
                    {logdir} \
                    {wildcards.sample} \
                    {params.threads} \
                    {input.r1} \
                    {input.r2} \
                    &>> {log}
            """)


############################################################
# RULE 3: align_reads — bwa-mem2 or bismark
############################################################
# memory hungry step.. scale down threads?
rule align_reads:
    input:
        r1 = f"{OUT_ROOT}/{{sample}}/trimmed/{{sample}}_R1_clean.fastq.gz",
        r2 = f"{OUT_ROOT}/{{sample}}/trimmed/{{sample}}_R2_clean.fastq.gz"
    output:
        bam = f"{OUT_ROOT}/{{sample}}/aligned/{{sample}}.unsorted.bam"
    params:
        script = "scripts/04_align.sh",
        aligner = lambda wc: SAMPLES[wc.sample]["aligner"],
        # ref = lambda wc: f'{config["reference_root"]}/{SAMPLES[wc.sample]["assembly"]}/genome.fa',
        ref = lambda wc: f'{config["reference_root"]}/{SAMPLES[wc.sample]["assembly"]}/{SAMPLES[wc.sample]["assembly"]}.fa',
        ref_dir = lambda wc: f'{config["reference_root"]}/{SAMPLES[wc.sample]["assembly"]}',
        threads = THREADS
    log:
        f"{OUT_ROOT}/{{sample}}/logs/align.log"
    message:
        "Aligning {wildcards.sample}"

    run:
        outdir = f"{OUT_ROOT}/{wildcards.sample}/aligned"
        logdir = f"{OUT_ROOT}/{wildcards.sample}/logs"
        os.makedirs(outdir, exist_ok=True)
        os.makedirs(logdir, exist_ok=True)

        shell("""
            bash {params.script} \
                {OUT_ROOT} \
                {logdir} \
                {wildcards.sample} \
                {params.aligner} \
                {params.threads} \
                {params.ref} \
                {params.ref_dir} \
                {input.r1} \
                {input.r2} \
                &>> {log}
        """)


############################################################
# RULE 4: sort BAM
############################################################

rule sort_bam:
    input:
        bam = f"{OUT_ROOT}/{{sample}}/aligned/{{sample}}.unsorted.bam"
    output:
        sorted = f"{OUT_ROOT}/{{sample}}/sorted/{{sample}}.sorted.bam"
    params:
        script = "scripts/05_sort.sh",
        threads = THREADS
    log:
        f"{OUT_ROOT}/{{sample}}/logs/sort.log"
    message:
        "Sorting BAM for {wildcards.sample}"

    run:
        outdir = f"{OUT_ROOT}/{wildcards.sample}/sorted"
        logdir = f"{OUT_ROOT}/{wildcards.sample}/logs"
        os.makedirs(outdir, exist_ok=True)
        os.makedirs(logdir, exist_ok=True)

        shell("""
            bash {params.script} \
                {OUT_ROOT} \
                {logdir} \
                {wildcards.sample} \
                {params.threads} \
                {input.bam} \
                &>> {log}
        """)


############################################################
# RULE 5: index BAM
############################################################

rule index_sorted_bam:
    input:
        bam = f"{OUT_ROOT}/{{sample}}/sorted/{{sample}}.sorted.bam"
    output:
        bai = f"{OUT_ROOT}/{{sample}}/sorted/{{sample}}.sorted.bam.bai"
    params:
        script = "scripts/05_index.sh",
        threads = THREADS
    log:
        f"{OUT_ROOT}/{{sample}}/logs/index.log"
    message:
        "Indexing {wildcards.sample}"

    run:
        logdir = f"{OUT_ROOT}/{wildcards.sample}/logs"
        os.makedirs(logdir, exist_ok=True)

        shell("""
            bash {params.script} \
                {OUT_ROOT} \
                {logdir} \
                {wildcards.sample} \
                {params.threads} \
                {input.bam} \
                &>> {log}
        """)


############################################################
# RULE 6: mark duplicates — Picard
############################################################

rule mark_duplicates:
    input:
        bam = f"{OUT_ROOT}/{{sample}}/sorted/{{sample}}.sorted.bam"
    output:
        dedup  = f"{OUT_ROOT}/{{sample}}/dedup/{{sample}}.dedup.bam",
    params:
        script = "scripts/07_markdup.sh",
        threads = THREADS
    log:
        f"{OUT_ROOT}/{{sample}}/logs/markdup.log"
    message:
        "Marking duplicates {wildcards.sample}"

    run:
        dedup_dir = f"{OUT_ROOT}/{wildcards.sample}/dedup"
        logdir = f"{OUT_ROOT}/{wildcards.sample}/logs"
        os.makedirs(dedup_dir, exist_ok=True)
        os.makedirs(logdir, exist_ok=True)

        shell("""
            bash {params.script} \
                {OUT_ROOT} \
                {logdir} \
                {wildcards.sample} \
                {params.threads} \
                {input.bam} \
                &>> {log}
        """)


############################################################
# RULE 7: add readgroups
############################################################

rule add_readgroup:
    input:
        bam = f"{OUT_ROOT}/{{sample}}/dedup/{{sample}}.dedup.bam"
    output:
        rg = f"{OUT_ROOT}/{{sample}}/rgfixed/{{sample}}.rg.bam"
    params:
        script = "scripts/08_add_readgroup.sh",
        threads = THREADS
    log:
        f"{OUT_ROOT}/{{sample}}/logs/readgroup.log"
    message:
        "Adding RG to {wildcards.sample}"

    run:
        rgdir = f"{OUT_ROOT}/{wildcards.sample}/rgfixed"
        logdir = f"{OUT_ROOT}/{wildcards.sample}/logs"
        os.makedirs(rgdir, exist_ok=True)
        os.makedirs(logdir, exist_ok=True)

        shell("""
            bash {params.script} \
                {OUT_ROOT} \
                {logdir} \
                {wildcards.sample} \
                {params.threads} \
                {input.bam} \
                &>> {log}
        """)


############################################################
# RULE 8: QC BAM — stats, flagstat, depth
############################################################
# TODO: make this step optional
rule qc_bam:
    input:
        bam = f"{OUT_ROOT}/{{sample}}/rgfixed/{{sample}}.rg.bam"
    output:
        stats   = f"{OUT_ROOT}/{{sample}}/qc/{{sample}}.stats.txt",
        flagstat= f"{OUT_ROOT}/{{sample}}/qc/{{sample}}.flagstat.txt",
        idxstats= f"{OUT_ROOT}/{{sample}}/qc/{{sample}}.idxstats.txt",
        depth   = f"{OUT_ROOT}/{{sample}}/qc/{{sample}}.depth.txt",
        summary = f"{OUT_ROOT}/{{sample}}/qc/{{sample}}.depth_summary.txt"
    params:
        script = "scripts/08_qc_bam.sh",
        threads = THREADS
    log:
        f"{OUT_ROOT}/{{sample}}/logs/qc.log"
    message:
        "QC: {wildcards.sample}"

    run:
        qcdir = f"{OUT_ROOT}/{wildcards.sample}/qc"
        logdir = f"{OUT_ROOT}/{wildcards.sample}/logs"
        os.makedirs(qcdir, exist_ok=True)
        os.makedirs(logdir, exist_ok=True)

        shell("""
            bash {params.script} \
                {OUT_ROOT} \
                {logdir} \
                {wildcards.sample} \
                {params.threads} \
                {input.bam} \
                &>> {log}
        """)
