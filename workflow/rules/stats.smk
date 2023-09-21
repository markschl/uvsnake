localrules:
    stats,


rule stats:
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
        tsv="results/sample_report.tsv",
        report=report("report/stats.txt", category="Sample report"),
    log:
        "logs/sample_report.log",
    conda:
        "envs/qc.yaml"
    script:
        "../scripts/stats_paired.py"
