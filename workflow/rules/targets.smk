"""
Target rules
"""


rule config:
    input:
        rules.dump_config.output,


rule quality:
    input:
        rules.multiqc.output,


rule quality_all:
    input:
        rules.multiqc_all.output,


rule unoise3:
    input:
        rules.config.input,
        expand(
            "results/{primers}/{what}",
            primers=config["_primer_combinations"],
            what=["unoise3.fasta", "unoise3_otutab.txt.gz", "unoise3.biom"],
        ),
        rules.stats.output,
        expand(rules.combine_logs.output, method="unoise3"),


rule uparse:
    input:
        rules.config.input,
        expand(
            "results/{primers}/{what}",
            primers=config["_primer_combinations"],
            what=["uparse.fasta", "uparse_otutab.txt.gz", "uparse.biom"],
        ),
        rules.stats.output,
        expand(rules.combine_logs.output, method="uparse"),


# rule post_cluster:
#     input:
#         rules.config.output,
#         expand(
#             "results/{primers}/{what}",
#             primers=config["_primer_combinations"],
#             what=["unoise3_cluster.fasta", "unoise3_cluster_otutab.txt.gz", "unoise3_cluster.biom"]
#         ),


# This rule is only run on existing files and will thus do nothing if there
# are no clustered/denoised FASTA files.
rule sintax:
    input:
        lambda wildcards: expand_clustered(
            "results/{primers}/{seqs}_{suffix}",
            suffix=[
                "sintax.txt.gz",
                "sintax_taxtab.txt.gz",
                "sintax.biom.gz",
                "sintax.fasta.gz",
            ],
            **wildcards,
        ),
