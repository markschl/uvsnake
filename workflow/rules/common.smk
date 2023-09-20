import os
from os.path import dirname
from glob import glob
from collections import OrderedDict

from utils.primers import get_primer_combinations
from utils.sample_list import SampleList


#### Initialize ####

from snakemake.workflow import srcdir

# initialize samples if corresponding sections present in configuration
# (otherwise we assume that uvsnake is used as Snakemake module and not all 
# functionality may be needed)
if "input" in config and "primers" in config:
    # further add these settings to config (_underscore indicates that they are special)
    config["_primers"], config["_primer_combinations"] = get_primer_combinations(config)
    l = SampleList(config["input"]["sample_file"])
    config["_input"] = OrderedDict(l.samples())
    config["_sample_names"] = list(config["_input"])
    assert len(config["_sample_names"]) == len(set(config["_sample_names"])), \
         "Duplicate sample names found"
    config["_layout"] = l.layout
# from pprint import pprint; pprint(config)


#### Helpers ####


def with_default(name, group):
    """
    Helper function that obtains a value for program/maxaccepts/maxrejects,
    or the default value from the 'defaults' section.
    This is necessary because validate() does not fill all defaults from the 
    JSON schema
    """
    value = config[group].get(name)
    if value is None:
        return config["defaults"][name]
    return value


def expand_clustered(path, **wildcards):
    for f in glob("results/*/*.fasta"):
        parts = f.split(os.sep)
        yield from expand(path, primers=parts[1], seqs=parts[2].split(".")[0], **wildcards)


def mem_func(mem=5, f=0, max_mem=50000):
    def _mem_func(wildcards, input, attempt):
        _mem = mem + f * input.size_mb
        return round(min(max_mem, _mem * 2**(attempt - 1)))
    return _mem_func


def time_func(time=1, f=0, max_time=24*60*20):
    def _time_func(wildcards, input, attempt):
        _time = time + f * input.size_mb
        return round(min(max_time, _time * 2**(attempt - 1)))
    return _time_func
