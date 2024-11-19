
rule import_taxdb:
    params:
        format=lambda wildcards: config["_taxdb_from_id"][wildcards.id][1],
    input:
        db=lambda wildcards: config["_taxdb_from_id"][wildcards.id][0],
    output:
        db="workdir/sintax/taxdb_{id}.fasta",
    log:
        "logs/sintax/import_taxdb_{id}.log",
    conda:
        "../envs/uvsnake.yaml"
    group:
        "taxonomy"
    script:
        "../scripts/import_taxdb.py"


rule assign_taxonomy_sintax:
    params:
        program=lambda _: config["sintax"]["program"],
        confidence=lambda wildcards: config["sintax"][wildcards.primers]["confidence"],
        strand=lambda wildcards: config["sintax"][wildcards.primers]["strand"],
        rand_seed=lambda wildcards: config["sintax"][wildcards.primers]["rand_seed"],
    input:
        fa="results/{primers}/{what}.fasta",
        db=lambda wildcards: expand(
            "workdir/sintax/taxdb_{id}.fasta",
            id=config["_taxdb_id_from_path"][config["sintax"][wildcards.primers]["db"]]
        ),
    output:
        sintax="results/{primers}/{what}_sintax.txt.gz",
        taxtab="results/{primers}/{what}_sintax_taxtab.txt.gz",
    log:
        "logs/sintax/{primers}/{what}_sintax.log",
    group:
        "taxonomy"
    conda:
        "../envs/uvsnake.yaml"
    threads: workflow.cores,
    resources:
        mem_mb=5000,
        runtime=240,
    script:
        "../scripts/taxonomy_sintax.sh"


rule make_tax_fasta:
    input:
        fa="results/{primers}/{what}.fasta",
        tax="results/{primers}/{what}_sintax_taxtab.txt.gz",
    output:
        "results/{primers}/{what}_sintax.fasta.gz",
    log:
        "logs/sintax/{primers}/{what}_make_tax_fasta.log",
    conda:
        "../envs/uvsnake.yaml"
    group:
        "taxonomy"
    shell:
        """
        tax="{input.tax}"
        st set -m <(gzip -dc "$tax") -d "{{meta(2)}}" "{input.fa}" |
          st replace -d '__' ':' |
          st replace -dr ' *; *' ' ' |
          gzip -nc > {output}
        """


rule make_tax_biom:
    input:
        biom="results/{primers}/{what}.biom",
        tax="results/{primers}/{what}_{method}_taxtab.txt.gz",
    output:
        tax_tmp=temp("workdir/{primers}/{what}_{method}_tax_tmp.txt"),
        biom="results/{primers}/{what}_{method}.biom.gz",
    log:
        "logs/sintax/{primers}/{what}_{method}_make_tax_biom.log",
    conda:
        "../envs/uvsnake.yaml"
    group:
        "taxonomy"
    script:
        "../scripts/make_tax_biom.sh"
