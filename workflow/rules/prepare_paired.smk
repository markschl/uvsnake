
localrules:
    make_primer_fasta,


rule merge_paired:
    params:
        program=lambda _: config["merge"]["program"],
        overlap_ident=config["merge"]["overlap_ident"],
        max_diffs=config["merge"]["max_diffs"],
    input:
        r1=lambda wildcards: config["_input"][wildcards.sample][0],
        r2=lambda wildcards: config["_input"][wildcards.sample][1],
    output:
        tempdir=temp(directory("workdir/prepare_paired/1_merge/{sample}/tmp")),
        merged=temp("workdir/prepare_paired/1_merge/{sample}/{sample}.fastq.zst"),
        nm_r1=temp("workdir/prepare_paired/1_merge/{sample}/{sample}_notmerged_R1.fastq.zst"),
        nm_r2=temp("workdir/prepare_paired/1_merge/{sample}/{sample}_notmerged_R2.fastq.zst"),
        stats=temp("workdir/prepare_paired/1_merge/{sample}/{sample}_stats.txt"),
    log:
        "logs/prepare_paired/1_merge/{sample}.log",
    conda:
        "../envs/uvsnake.yaml"
    group:
        "sample"
    resources:
        # TODO: cannot use lambda functions (https://github.com/snakemake/snakemake/issues/2154#issuecomment-1553538507)
        # mem_mb=mem_func(40 if config["merge"]["program"] == "vsearch" else 120),  # VSEARCH: ~16MB, USEARCH: <=78MB
        runtime=60,  # time_func(5, f=0.2),  # 0.003-0.006 min/MB of gzip input
    script:
        "../scripts/merge_paired.sh"


rule trim_primers_paired:
    params:
        max_error_rate=config["primers"]["trim_settings"]["max_error_rate"],
        min_overlap=config["primers"]["trim_settings"]["min_overlap"],
        min_length=config["primers"]["trim_settings"]["min_length"],
        primer_comb=config["_primer_combinations"],
    input:
        fprimers="workdir/primers/forward.fasta",
        rprimers_rev="workdir/primers/reverse_rev.fasta",
        seq="workdir/prepare_paired/1_merge/{sample}/{sample}.fastq.zst",
    output:
        fwd_log="workdir/prepare_paired/2_trim/{sample}/{sample}_fwd.log",
        rev_log="workdir/prepare_paired/2_trim/{sample}/{sample}_rev.log",
        stats=temp("workdir/prepare_paired/2_trim/{sample}/{sample}_stats.txt"),
        by_primers=temp(expand(
            "workdir/prepare_paired/2_trim/{{sample}}/{primers}.fastq.zst",
            primers=config["_primer_combinations"],
        )),
        short="workdir/prepare_paired/2_trim/{sample}/too_short.fastq.zst",
    log:
        "logs/prepare_paired/2_trim/{sample}.log",
    conda:
        "../envs/uvsnake.yaml"
    group:
        "sample"
    resources:
        # mem_mb=mem_func(100),  # 57MB
        runtime=60,  # time_func(5, 0.15),  # 0.06min/MB (zstd input)
    script:
        "../scripts/trim_primers_paired.sh"
