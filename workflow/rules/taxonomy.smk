
rule import_taxdb:
    params:
        format=lambda _: config["sintax"]["db_format"]
    input:
        db=lambda _: config["sintax"]["db"],
    output:
        db="workdir/sintax/taxdb.fasta",
    log:
        "logs/sintax/import_taxdb.log",
    conda:
        "envs/biopython.yaml"
    group:
        "taxonomy"
    script:
        "scripts/import_taxdb.py"


rule assign_taxonomy_sintax:
    params:
        confidence=lambda _: config["sintax"]["confidence"],
        program=lambda _: with_default("program", "sintax"),
        usearch_bin=lambda _: config["usearch_binary"],
        # TODO: not sure if USEARCH even uses these settings (VSEARCH doesnt')
        maxaccepts=lambda _: with_default("maxaccepts", "sintax"),
        maxrejects=lambda _: with_default("maxrejects", "sintax"),
    input:
        fa="results/{primers}/{what}.fasta",
        db="workdir/sintax/taxdb.fasta",
    output:
        sintax="results/{primers}/{what}_sintax.txt.gz",
        taxtab="results/{primers}/{what}_sintax_taxtab.txt.gz",
    log:
        "logs/sintax/{primers}/{what}_sintax.log",
    group:
        "taxonomy"
    conda:
        "envs/vsearch.yaml"
    # threads:
    # VSEARCH works in parallel (although cores seem to be used only ~50%) while
    # USEARCH v11 does not appear to use more than 1 thread
    threads:
        lambda _: workflow.cores \
        if with_default("program", "sintax") == "vsearch" \
        else 1
    resources:
        mem_mb=5000,
    script:
        "scripts/taxonomy_sintax.sh"


rule make_tax_fasta:
    input:
        fa="results/{primers}/{what}.fasta",
        tax="results/{primers}/{what}_sintax_taxtab.txt.gz",
    output:
        "results/{primers}/{what}_sintax.fasta.gz",
    log:
        "logs/sintax/{primers}/{what}_make_tax_fasta.log",
    conda:
        "envs/basic.yaml"
    group:
        "taxonomy"
    shell:
        """
        tax="{input.tax}"
        st set -ul <(gzip -dc "$tax") -d {{l:2}} "{input.fa}" |
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
        "logs/taxonomy/{primers}/{what}_{method}_make_tax_biom.log",
    conda:
        "envs/biom.yaml"
    group:
        "taxonomy"
    script:
        "scripts/make_tax_biom.sh"
