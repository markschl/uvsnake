"""
Rules for cleaning up
"""


localrules:
    clean,
    clean_all,


rule clean:
    shell:
        "rm -Rf workdir logs"


rule clean_all:
    shell:
        "rm -Rf results qc workdir logs"
