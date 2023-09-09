
localrules:
    dump_config,
    make_cutadapt_fasta,


rule dump_config:
    params:
        config=config,
    output:
        "results/config.yaml",
    log:
        "logs/dump_config.log",
    conda:
        "envs/basic.yaml"
    script:
        "scripts/dump_config.py"


rule make_cutadapt_fasta:
    params:
        primer_config=config["_primers"],
    output:
        fwd="workdir/primers/forward.fasta",
        rev="workdir/primers/reverse.fasta",
        reverse_rev="workdir/primers/reverse_rev.fasta",
    log:
        "logs/make_primer_fasta.log",
    script:
        "scripts/make_primer_fasta.py"


rule combine_logs:
    input:
        merge=expand(rules.merge_paired.log[0], sample=config["_sample_names"]),
        trim=expand(rules.trim_primers_paired.log[0], sample=config["_sample_names"]),
        filter_derep=expand(rules.filter_derep.log[0],
                            sample=config["_sample_names"],
                            primers=config["_primer_combinations"]),
        collect_derep=expand(rules.collect_derep.log[0],
                             primers=config["_primer_combinations"]),
        cluster=expand("logs/cluster/3_cluster/{primers}_{{method}}.log",
                       primers=config["_primer_combinations"]),
    output:
        "logs/cluster_{method}_all.log",
    shell:
        """
        exec 1> "{output}"
        echo "Paired-end read merging"
        echo "======================="
        cat {input.merge:q}
        printf "\n\n\n"
        echo "Primer trimming"
        echo "================"
        cat {input.trim:q}
        printf "\n\n\n"
        echo "Filter & de-replicate (per-sample)"
        echo "=================================="
        cat {input.filter_derep:q}
        printf "\n\n\n"
        echo "Combine and de-replicate globally"
        echo "================================="
        cat {input.collect_derep:q}
        printf "\n\n\n"
        echo "Cluster"
        echo "======="
        cat {input.cluster:q}
        """
