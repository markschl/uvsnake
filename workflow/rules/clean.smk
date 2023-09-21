"""
Rules for cleaning up
"""


localrules:
    clean,
    clean_all,


rule clean:
    log:
        "logs/clean.log",
    conda:
        "envs/uvsnake.yaml"
    shell:
        "rm -Rf workdir logs 2> {log}"


rule clean_all:
    log:
        "logs/clean_all.log",
    conda:
        "envs/uvsnake.yaml"
    shell:
        "rm -Rf results qc workdir logs 2> {log}"
