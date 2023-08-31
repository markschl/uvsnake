import os
from os.path import dirname
from glob import glob

from scripts.utils.primers import get_primer_combinations
from scripts.utils.sample_list import SampleList


#### Prepare ####

from snakemake.workflow import srcdir

# Set environment variable of workflow root dir to allow post-deploy
# scripts to run other scripts stored in that directory
os.environ['PIPELINE_DIR'] = dirname(dirname(dirname(srcdir('.'))))

# further add these settings to config (_underscore indicates that they are special)
config["_primers"], config["_primer_combinations"] = get_primer_combinations(config)
l = SampleList(config["input"]["sample_file"])
sample_names = list(s for s, _ in l.samples())
assert len(sample_names) == len(set(sample_names)), "Duplicate sample names found"
config["_input"] = dict(l.samples())
config["_sample_names"] = sample_names
config["_layout"] = l.layout

def with_default(name, group):
    value = config[group].get(name)
    if value is None:
        return config["defaults"][name]
    return value


#### Other ####

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


def expand_clustered(path, **wildcards):
    for f in glob("results/*/*.fasta"):
        parts = f.split(os.sep)
        yield from expand(path, primers=parts[1], seqs=parts[2].split(".")[0], **wildcards)
