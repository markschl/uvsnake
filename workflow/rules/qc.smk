from os.path import dirname

localrules:
    stats_paired,
    clean_qc,


rule fastqc:
    input:
        lambda wildcards: config["_input"][wildcards.sample],
    output:
        outdir=directory("qc/fastqc/{sample}"),
        outfiles=expand(
            "qc/fastqc/{{sample}}/{{sample}}_R{read}_fastqc.{ext}",
            read=[1, 2] if config["_layout"] == "paired" else [1],
            ext=("html", "zip"),
        ),
    log:
        "logs/qc/fastqc/{sample}.log",
    group:
        "sample"
    conda:
        "envs/qc.yaml"
    shell:
        """
        fastqc -q -f fastq -t 1 -o {output.outdir} {input} &> {log}
        """


rule multiqc:
    params:
        fastqc_dir=lambda w, input: dirname(dirname(input.fastqc[0])),
    input:
        fastqc=expand(
            "qc/fastqc/{sample}/{sample}_R{read}_fastqc.html",
            sample=config["_sample_names"],
            read=[1, 2],
        ),
    output:
        "qc/multiqc/multiqc_report.html",
    log:
        "logs/qc/multiqc.log",
    conda:
        "envs/qc.yaml"
    resources:
        mem_mb=mem_func(1000),
        runtime=time_func(120),
    shell:
        """
        multiqc -f -m fastqc -o $(dirname {output}) {params.fastqc_dir} 2> {log}
        """


rule multiqc_all:
    input:
        fastqc=rules.multiqc.input,
        cutadapt=expand(
            "workdir/prepare_paired/2_trim/{sample}/{sample}_{direction}.log",
            sample=config["_sample_names"],
            direction=["fwd", "rev"],
        ),
    output:
        "qc/multiqc_all/multiqc_report.html",
    log:
        "logs/qc/multiq_all.log",
    conda:
        "envs/qc.yaml"
    resources:
        mem_mb=mem_func(1000),
        runtime=time_func(120),
    shell:
        """
        indir=$(dirname $(dirname {input.fastqc[0]}))
        outdir=$(dirname {output})
        cutadapt_dir=$(dirname $(dirname {input.cutadapt[0]}))
        multiqc -f -m fastqc -m cutadapt -o $outdir $indir $cutadapt_dir 2> {log}
        """


rule stats_paired:
    params:
        primer_combinations=config["_primer_combinations"],
    input:
        merge=expand(
            "workdir/prepare_paired/1_merge/{sample}/{sample}_stats.txt",
            sample=config["_sample_names"],
        ),
        trim=expand(
            "workdir/prepare_paired/2_trim/{sample}/{sample}_stats.txt",
            sample=config["_sample_names"],
        ),
        filter=expand(
            "workdir/cluster/1_filter_derep/{primers}/{sample}/{sample}_stats.txt",
            primers=config["_primer_combinations"],
            sample=config["_sample_names"],
        ),
    output:
        "results/sample_report.tsv",
    log:
        "logs/sample_report.log",
    conda:
        "envs/qc.yaml"
    script:
        "scripts/stats_paired.py"


rule clean_qc:
    log:
        "logs/clean_qc.log",
    shell:
        "rm -Rf qc 2> {log}"
