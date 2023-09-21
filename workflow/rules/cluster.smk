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
        "envs/uvsnake.yaml"
    group:
        "sample"
    resources:
        # mem_mb=mem_func(300), # very variable 10 - 230MB
        runtime=30,  #time_func(5, f=0.2),  # 0.004min/MiB  
    script:
        "scripts/filter_derep.sh"


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
        "envs/uvsnake.yaml"
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


# Memory and time constraints for HPC:
# Approximation based on benchmarks and non-linear regression
# (could always be inproved)
# Some safety margin is added by multiplying the memory and time constraints with a factor
# In addition, constraints are doubled at every attempt of re-running after a fail
_mem_f = 2
_time_f = 4


def unoise3_memfunc(program, min_size, threads):
    if program == "usearch":
        return lambda _, input, attempt: (
            _mem_f * (180 * input.size_mb * min_size**-0.9 + 50)
        ) * 2 ** (attempt - 1)
    return lambda _, input, attempt: (
        _mem_f * (90 * input.size_mb * min_size**-1.4 * (threads - 0.46) + 50)
    ) * 2 ** (attempt - 1)


def unoise3_timefunc(program, min_size, maxrejects, threads):
    if program == "usearch":
        return lambda _, input, attempt: (
            _time_f * (30 * input.size_mb**1.6 * min_size**-2 + 5)
        ) * 2 ** (attempt - 1)
    return lambda _, input, attempt: (
        _time_f
        * (
            0.4
            * input.size_mb**1.52
            * min_size**-2
            * (threads**-0.63 - 0.58)
            * (maxrejects + 4)
            + 5
        )
    ) * 2 ** (attempt - 1)


def uparse_memfunc(min_size):
    return lambda _, input, attempt: (
        _mem_f * (30 * input.size_mb**0.5 * min_size**0.2 + 50)
    ) * 2 ** (attempt - 1)


def uparse_timefunc(min_size, maxaccepts):
    return lambda _, input, attempt: (
        _time_f * (2 * input.size_mb * min_size**-2 * maxaccepts**0.75 + 5)
    ) * 2 ** (attempt - 1)


from os.path import getsize
from math import log


def otutab_memfunc(program, threads):
    if program == "usearch":
        return lambda _, input, attempt: (
            _mem_f
            * (
                4e-4
                * (log(getsize(input.uniques)) ** 8.3e-3 - 1)
                * (getsize(input.otus) + 2e6)
                * (threads + 5.1)
                + 50
            )
        ) * 2 ** (attempt - 1)
    return lambda _, input, attempt: (
        _mem_f * (2e-5 * getsize(input.otus) * (threads - 0.8) + 20)
    ) * 2 ** (attempt - 1)


def otutab_timefunc(program, maxaccepts, maxrejects, threads):
    if program == "usearch":
        return lambda _, input, attempt: (
            _time_f
            * (
                1.5e-13
                * threads**-1
                * getsize(input.uniques)
                * (getsize(input.otus) + 4.2e5)
                * maxrejects**0.65
                + 5
            )
        ) * 2 ** (attempt - 1)
    return lambda _, input, attempt: (
        _time_f
        * (
            4e-15
            * threads**-1
            * getsize(input.uniques)
            * (getsize(input.otus) + 8.4e6)
            * (maxrejects + 20)
            + 5
        )
    ) * 2 ** (attempt - 1)


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
        "envs/uvsnake.yaml"
    # threads:
    # VSEARCH works in parallel (although cores seem to be used only ~50%) while
    # USEARCH v11 does not appear to use more than 1 thread
    threads: workflow.cores if with_default("program", "unoise3") == "vsearch" else 1
    resources:
        mem_mb=unoise3_memfunc(
            with_default("program", "unoise3"),
            config["unoise3"]["min_size"],
            workflow.cores,
        ),
        runtime=unoise3_timefunc(
            with_default("program", "unoise3"),
            config["unoise3"]["min_size"],
            with_default("maxrejects", "unoise3"),
            workflow.cores,
        ),
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
        "envs/uvsnake.yaml"
    resources:
        mem_mb=uparse_memfunc(config["uparse"]["min_size"]),
        runtime=uparse_timefunc(
            config["uparse"]["min_size"], with_default("maxaccepts", "uparse")
        ),
    shell:
        """
        exec &> {log[0]}
        set -xeuo pipefail
        zstd -dcq {input} > {output.tmp_in}
        "{params.usearch_bin}" \
            -cluster_otus {output.tmp_in} \
            -otus {output.fa} \
            -relabel Otu \
            -maxaccepts {params.maxaccepts} \
            -maxrejects {params.maxrejects} \
            -minsize {params.min_size}
        """


rule usearch_make_otutab:
    params:
        ident_threshold=config["otutab"]["ident_threshold"] / 100,
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
        "envs/uvsnake.yaml"
    resources:
        mem_mb=otutab_memfunc(with_default("program", "otutab"), workflow.cores),
        runtime=otutab_timefunc(
            with_default("program", "otutab"),
            with_default("maxaccepts", "otutab"),
            with_default("maxrejects", "otutab"),
            workflow.cores,
        ),
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
#         "envs/uvsnake.yaml"
#     group:
#         "cluster"
#     threads:
#         int(workflow.cores * 1.5) if with_default("program", "post_cluster") == "vsearch" else 1
#     script:
#         "scripts/post_cluster.py"
