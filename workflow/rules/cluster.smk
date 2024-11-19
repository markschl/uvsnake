# filter and de-replicate (uses VSEARCH)
rule filter_derep:
    params:
        maxee_rate=config["filter"]["max_error_rate"],
    input:
        expand(
            "workdir/prepare_{layout}/2_trim/{{sample}}/{{primers}}.fastq.zst",
            layout=config["_layout"],
        ),
    output:
        filtered_tmp=temp(
            "workdir/cluster/1_filter_derep/{primers}/{sample}/{sample}_filtered.fasta"
        ),
        discarded_tmp=temp(
            "workdir/cluster/1_filter_derep/{primers}/{sample}/{sample}_discarded.fasta"
        ),
        all="workdir/cluster/1_filter_derep/{primers}/{sample}/{sample}_all_uniques.fasta.zst",
        good="workdir/cluster/1_filter_derep/{primers}/{sample}/{sample}_good_uniques.fasta.zst",
        stats="workdir/cluster/1_filter_derep/{primers}/{sample}/{sample}_stats.txt",
    log:
        "logs/cluster/1_filter_derep/{primers}/{sample}.log",
    conda:
        "../envs/uvsnake.yaml"
    group:
        "sample"
    resources:
        # mem_mb=mem_func(300), # very variable 10 - 230MB
        runtime=30,  #time_func(5, f=0.2),  # 0.004min/MiB  
    script:
        "../scripts/filter_derep.sh"


# Combine de-replicated sample files into one
rule collect_derep:
    input:
        good=expand(
            "workdir/cluster/1_filter_derep/{{primers}}/{sample}/{sample}_good_uniques.fasta.zst",
            sample=config["_sample_names"],
        ),
        all=expand(
            "workdir/cluster/1_filter_derep/{{primers}}/{sample}/{sample}_all_uniques.fasta.zst",
            sample=config["_sample_names"],
        ),
    output:
        good="workdir/cluster/2_unique_all/{primers}/good_uniques.fasta.zst",
        all="workdir/cluster/2_unique_all/{primers}/all_uniques.fasta.zst",
    log:
        "logs/cluster/2_unique_all/{primers}.log",
    conda:
        "../envs/uvsnake.yaml"
    resources:
        mem_mb=mem_func(10, f=20),  #  ~11MB per compressed input MB
        runtime=time_func(5, f=0.08),  # 0.0014min/MiB  
    shell:
        """
        exec &> "{log}"
        set -xeuo pipefail

        # First, combine the filtered and de-replicated sample files,
        # and de-replicate these sequences again.
        # This is important, because UNOISE assumes unique sequences
        # with size annotations.
        # TODO: the output is sorted by size according to VSEARCH docs.
        #   -> Make sure that this stays the same when updating to
        #      future versions.
        zstd -dcq {input.good:q} | 
            vsearch -derep_fulllength - -sizein -sizeout -output - |
            zstd -cq > "{output.good}"

        # Second, combine unfiltered (actually, length filtered after trimming)
        # and de-replicated sequences into one file. These will be used
        # for mapping against the clustered sequences to create the OTU table.
        # It is important *not* to de-replicate them again here,
        # otherwise part of the sample labels will be lost, leading to 
        # incorrect results.
        zstd -dcq {input.all:q} |
            zstd -cq > "{output.all}"
        """


rule cluster_unoise3:
    params:
        program=lambda _: config["unoise3"]["program"],
        min_size=lambda _: config["unoise3"]["min_size"],
        maxaccepts=lambda _: config["unoise3"]["maxaccepts"],
        maxrejects=lambda _: config["unoise3"]["maxrejects"],
    input:
        "workdir/cluster/2_unique_all/{primers}/good_uniques.fasta.zst",
    output:
        "results/{primers}/unoise3.fasta",
    log:
        "logs/cluster/3_cluster/{primers}_unoise3.log",
    conda:
        "../envs/uvsnake.yaml"
    threads:
        cluster_cores(config["unoise3"]["program"], len(config["_primer_combinations"]), workflow.cores)
    resources:
        mem_mb=unoise3_memfunc(
            config["unoise3"]["program"],
            config["unoise3"]["min_size"],
            workflow.cores,
        ),
        runtime=unoise3_timefunc(
            config["unoise3"]["program"],
            config["unoise3"]["min_size"],
            config["unoise3"]["maxrejects"],
            workflow.cores,
        ),
    script:
        "../scripts/unoise3.sh"


rule cluster_uparse:
    params:
        program=lambda _: config["uparse"]["program"],
        min_size=lambda _: config["uparse"]["min_size"],
        maxaccepts=lambda _: config["uparse"]["maxaccepts"],
        maxrejects=lambda _: config["uparse"]["maxrejects"],
    input:
        "workdir/cluster/2_unique_all/{primers}/good_uniques.fasta.zst",
    output:
        tmp_in=temp("results/{primers}/uparse_input_tmp.fasta"),
        fa="results/{primers}/uparse.fasta",
    log:
        "logs/cluster/3_cluster/{primers}_uparse.log",
    conda:
        "../envs/uvsnake.yaml"
    resources:
        mem_mb=uparse_memfunc(config["uparse"]["min_size"]),
        runtime=uparse_timefunc(
            config["uparse"]["min_size"],
            config["uparse"]["maxaccepts"],
        ),
    shell:
        """
        exec &> {log[0]}
        set -xeuo pipefail
        zstd -dcq {input} > {output.tmp_in}
        "{params.program}" \
            -cluster_otus {output.tmp_in} \
            -otus {output.fa} \
            -relabel Otu \
            -maxaccepts {params.maxaccepts} \
            -maxrejects {params.maxrejects} \
            -minsize {params.min_size}
        """


rule make_otutab:
    params:
        program=lambda _: config["otutab"]["program"],
        ident_threshold=lambda _: config["otutab"]["ident_threshold"] / 100,
        maxaccepts=lambda _: config["otutab"]["maxaccepts"],
        maxrejects=lambda _: config["otutab"]["maxrejects"],
    input:
        otus="results/{primers}/{what}.fasta",
        uniques="workdir/cluster/2_unique_all/{primers}/all_uniques.fasta.zst",
    output:
        tab="results/{primers}/{what}_otutab.txt.gz",
        biom="results/{primers}/{what}.biom",
        notmatched="workdir/cluster/4_otutab/{primers}/{what}_otutab_notmatched.fasta.zst",
        extra=otutab_extra_files(bam=False),
    threads: otutab_cores(config["otutab"]["program"], len(config["_primer_combinations"]), workflow.cores)
    log:
        "logs/cluster/4_otutab/{primers}_{what}.log",
    conda:
        "../envs/uvsnake.yaml"
    resources:
        mem_mb=otutab_memfunc(config["otutab"]["program"], workflow.cores),
        runtime=otutab_timefunc(
            config["otutab"]["program"],
            config["otutab"]["maxaccepts"],
            config["otutab"]["maxrejects"],
            workflow.cores,
        ),
    script:
        "../scripts/make_otutab.sh"


rule otutab_sam2bam:
    input:
        otus="results/{primers}/{what}.fasta",
        sam="workdir/cluster/4_otutab/{primers}/{what}.sam",
    output:
        bam="workdir/cluster/4_otutab/{primers}/{what}.bam",
        bai="workdir/cluster/4_otutab/{primers}/{what}.bam.bai",
    log:
        "logs/cluster/4_otutab/{primers}_{what}_sam2bam.log",
    conda:
        "../envs/samtools.yaml"
    threads:
        workflow.cores,
    resources:
        mem_mb=mem_func(100, 2),
        runtime=time_func(10, 0.05 / workflow.cores),
    shell:
        """
        rm -f "{input.otus}".fai
        samtools view -T "{input.otus}" -b "{input.sam}" |
            samtools sort -m {resources.mem_mb}M -@ {threads} >{output.bam}
        rm -f "{input.otus}".fai
        samtools index "{output.bam}"
        """


# rule post_cluster:
#     params:
#       threshold=lambda _: config["post_cluster"]["ident_threshold"]/100,
#       program=lambda _: config["post_cluster"]["program"],
#       maxaccepts=lambda _: config["post_cluster"]["maxaccepts"],
#       maxrejects=lambda _: config["post_cluster"]["maxrejects"],
#     input:
#         seqs="results/{primers}/unoise3.fasta",
#         biom="results/{primers}/unoise3.biom",
#     output:
#         seqs="results/{primers}/unoise3_cluster.fasta",
#         map="results/{primers}/unoise3_cluster_map.txt.gz",
#     log:
#         "logs/cluster/3_cluster/{primers}_post_cluster.log",
#     conda:
#         "../envs/uvsnake.yaml"
#     group:
#         "cluster"
#     threads:
#         int(workflow.cores * 1.5) if config["post_cluster"]["program"] == "vsearch" else 1
#     script:
#         "../scripts/post_cluster.py"
