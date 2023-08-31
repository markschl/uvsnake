

# filter and de-replicate (uses VSEARCH)
rule filter_derep:
    params:
        maxee_rate=config["filter"]["max_error_rate"],
        minlen=config["filter"]["min_length"],  # not absolutely required (done after trimming)
    input:
        expand("workdir/prepare_{layout}/2_trim/{{sample}}/{{primers}}.fastq.zst", layout=config["_layout"]),
    output:
        filtered_tmp=temp("workdir/cluster/1_filter_derep/{primers}/{sample}/{sample}_filtered.fasta"),
        discarded_tmp=temp("workdir/cluster/1_filter_derep/{primers}/{sample}/{sample}_discarded.fasta"),
        all="workdir/cluster/1_filter_derep/{primers}/{sample}/{sample}_all_uniques.fasta.zst",
        good="workdir/cluster/1_filter_derep/{primers}/{sample}/{sample}_good_uniques.fasta.zst",
        stats="workdir/cluster/1_filter_derep/{primers}/{sample}/{sample}_stats.txt",
    log:
        "logs/cluster/1_filter_derep/{primers}/{sample}.log",
    conda:
        "envs/vsearch.yaml"
    group:
        "sample"
    resources:
        mem_mb=2000,
    script:
        "scripts/filter_derep.sh"


# Combine de-replicated sample files into one
rule collect_derep:
    input:
        good=expand(
            "workdir/cluster/1_filter_derep/{{primers}}/{sample}/{sample}_good_uniques.fasta.zst",
            sample=config["_sample_names"]
        ),
        all=expand(
            "workdir/cluster/1_filter_derep/{{primers}}/{sample}/{sample}_all_uniques.fasta.zst",
            sample=config["_sample_names"]
        ),
    output:
        good="workdir/cluster/2_unique_all/{primers}/good_uniques.fasta.zst",
        all="workdir/cluster/2_unique_all/{primers}/all_uniques.fasta.zst",
    log:
        "logs/cluster/2_unique_all/{primers}.log",
    group:
        "denoise"
    conda:
        "envs/vsearch.yaml"
    resources:
        mem_mb=5000,
    script:
        "scripts/collect_derep.sh"


rule cluster_unoise3:
    params:
        min_size=config["unoise3"]["min_size"],
        program=with_default("program", "unoise3"),
        usearch_bin=config["usearch_binary"],
        maxaccepts=with_default("maxaccepts", "unoise3"),
        maxrejects=with_default("maxrejects", "unoise3"),
    input:
        "workdir/cluster/2_unique_all/{primers}/good_uniques.fasta.zst",
    output:
        "results/{primers}/unoise3.fasta",
    log:
        "logs/cluster/3_cluster/{primers}_unoise3.log",
    conda:
        "envs/vsearch.yaml"
    group:
        "cluster"
    # threads:
    # VSEARCH works in parallel (although cores seem to be used only ~50%) while
    # USEARCH v11 does not appear to use more than 1 thread
    # TODO: further validate VSEARCH threads setting
    threads:
        int(workflow.cores * 1.5) if with_default("program", "unoise3") == "vsearch" else 1
    resources:
        mem_mb=10000,
        runtime=24 * 60,
    script:
        "scripts/unoise3.sh"


rule cluster_uparse:
    params:
        usearch_bin=config["usearch_binary"],
        min_size=config["uparse"]["min_size"],
        maxaccepts=with_default("maxaccepts", "uparse"),
        maxrejects=with_default("maxrejects", "uparse"),
    input:
        "workdir/cluster/2_unique_all/{primers}/good_uniques.fasta.zst",
    output:
        tmp_in=temp("results/{primers}/uparse_input_tmp.fasta"),
        fa="results/{primers}/uparse.fasta",
    log:
        "logs/cluster/3_cluster/{primers}_uparse.log",
    conda:
        "envs/basic.yaml"
    group:
        "cluster"
    resources:
        mem_mb=10000,
        runtime=24 * 60,
    shell:
        """
        zstd -dcq {input} > {output.tmp_in}
        "{params.usearch_bin}" \
            -cluster_otus {output.tmp_in} \
            -otus {output.fa} \
            -relabel Otu \
            -maxaccepts {params.maxaccepts} \
            -maxrejects {params.maxrejects} \
            -minsize {params.min_size} &> {log}
        """


rule usearch_make_otutab:
    params:
        ident_threshold=config["otutab"]["ident_threshold"]/100,
        program=with_default("program", "otutab"),
        usearch_bin=config["usearch_binary"],
        maxaccepts=with_default("maxaccepts", "otutab"),
        maxrejects=with_default("maxrejects", "otutab"),
        # optional diagnostic output
        extra="true" if config["otutab"].get("extra", False) else "false",
        bam_out="workdir/otutab_mapping/{primers}/{what}.bam",
        map_out="workdir/otutab_mapping/{primers}/{what}_search.txt.gz",
    input:
        otus="results/{primers}/{what}.fasta",
        uniques="workdir/cluster/2_unique_all/{primers}/all_uniques.fasta.zst",
    output:
        tab="results/{primers}/{what}_otutab.txt.gz",
        biom="results/{primers}/{what}.biom",
        notmatched="workdir/cluster/4_otutab/{primers}/{what}_otutab_notmatched.fasta.zst",
    threads: workflow.cores
    log:
        "logs/cluster/4_otutab/{primers}_{what}.log",
    conda:
        lambda wildcards: "envs/vsearch-samtools.yaml" \
            if config["otutab"].get("extra", False) \
            else "envs/vsearch.yaml"
    group:
        "otutab"
    resources:
        mem_mb=10000,
        runtime=24 * 60,
    script:
        "scripts/make_otutab.sh"


# rule post_cluster:
#     params:
#         threshold=config["post_cluster"]["ident_threshold"]/100,
#         program=with_default("program", "post_cluster"),
#         usearch_bin=config["usearch_binary"],
#         maxaccepts=with_default("maxaccepts", "post_cluster"),
#         maxrejects=with_default("maxrejects", "post_cluster"),
#     input:
#         seqs="results/{primers}/unoise3.fasta",
#         biom="results/{primers}/unoise3.biom",
#     output:
#         seqs="results/{primers}/unoise3_cluster.fasta",
#         map="results/{primers}/unoise3_cluster_map.txt.gz",
#     log:
#         "logs/cluster/3_cluster/{primers}_post_cluster.log",
#     conda:
#         "envs/vsearch.yaml"
#     group:
#         "cluster"
#     threads:
#         int(workflow.cores * 1.5) if with_default("program", "post_cluster") == "vsearch" else 1
#     script:
#         "scripts/post_cluster.py"
