from os.path import dirname


localrules:
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
        "../envs/qc.yaml"
    shell:
        r"""
        exec &> "{log}"
        fastqc -q -f fastq -t 1 -o {output.outdir} {input:q}
        # rename output
        i=1
        for f in {input:q}; do
          sample=$(basename $f | sed -E 's/(\.(fq|fastq))?\.gz$//g')
          orig={output.outdir}/$sample"_fastqc"
          out="{output.outdir}/{wildcards.sample}"_R"$((i++))"_fastqc
          if [ "$orig.html" != "$out.html" ]; then
            mv -f $orig.html $out.html
            mv -f $orig.zip $out.zip
          fi
        done
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
        "../envs/qc.yaml"
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
        "../envs/qc.yaml"
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


rule clean_qc:
    log:
        "logs/clean_qc.log",
    conda:
        "../envs/uvsnake.yaml"
    shell:
        "rm -Rf qc 2> {log}"
