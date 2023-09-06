
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
