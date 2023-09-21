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
        "../envs/qc.yaml"
    script:
        "../scripts/stats_paired.py"


rule combine_logs:
    input:
        merge=expand(rules.merge_paired.log[0], sample=config["_sample_names"]),
        trim=expand(rules.trim_primers_paired.log[0], sample=config["_sample_names"]),
        filter_derep=expand(
            rules.filter_derep.log[0],
            sample=config["_sample_names"],
            primers=config["_primer_combinations"],
        ),
        collect_derep=expand(
            rules.collect_derep.log[0], primers=config["_primer_combinations"]
        ),
        cluster=expand(
            "logs/cluster/3_cluster/{primers}_{{method}}.log",
            primers=config["_primer_combinations"],
        ),
        otutab=expand(
            "logs/cluster/4_otutab/{primers}_{{method}}.log",
            primers=config["_primer_combinations"],
        ),
    output:
        "logs/cluster_{method}_all.log",
    log:
        "logs/combine_logs_{method}.log",
    conda:
        "../envs/uvsnake.yaml"
    shell:
        """
        exec 2> "{log}"
        exec 1> "{output}"
        echo "Paired-end read merging"
        echo "======================="
        tail -n+1 {input.merge:q}
        printf "\n\n\n"
        echo "Primer trimming"
        echo "================"
        tail -n+1 {input.trim:q}
        printf "\n\n\n"
        echo "Filter & de-replicate (per-sample)"
        echo "=================================="
        tail -n+1 {input.filter_derep:q}
        printf "\n\n\n"
        echo "Combine and de-replicate globally"
        echo "================================="
        tail -n+1 {input.collect_derep:q}
        printf "\n\n\n"
        echo "Cluster"
        echo "======="
        tail -n+1 {input.cluster:q}
        printf "\n\n\n"
        echo "OTU tab"
        echo "======="
        tail -n+1 {input.otutab:q}
        """
