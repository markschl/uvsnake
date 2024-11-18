import os
from os.path import dirname, exists
from glob import glob
from collections import OrderedDict


include: "utils.smk"


#### Initialize ####


# initialize samples if corresponding sections present in configuration
# (otherwise we assume that uvsnake is used as Snakemake module and not all
# functionality may be needed)
if "primers" in config:
    config["_primers"], config["_primer_combinations"] = get_primer_combinations(config)

if "sample_file" in config and exists(config["sample_file"]):
    # further add these settings to config (_underscore indicates that they are special)
    l = SampleList(config["sample_file"])
    config["_input"] = OrderedDict(l.samples())
    config["_sample_names"] = list(config["_input"])
    assert len(config["_sample_names"]) == len(
        set(config["_sample_names"])
    ), "Duplicate sample names found"
    config["_layout"] = l.layout

# from pprint import pprint; pprint(config)


#### Helpers ####


def cfg_path(*keys, default=None):
    """
    Get a nested dict entry or default if non-existent
    """
    return _cfg_path(config, *keys)


def _cfg_path(cfg, *keys, default=None):
    if len(keys) == 0:
        return cfg
    try:
        sub = cfg[keys[0]]
        return _cfg_path(sub, *keys[1:])
    except KeyError:
        return default


def cfg_or_global_default(*keys, fallback=None):
    """
    Helper function that obtains a value for program/maxaccepts/maxrejects,
    or the default value from the 'defaults' section.
    This is necessary because validate() does not fill all defaults from the
    JSON schema
    """
    value = cfg_path(*keys, default=fallback)
    if value is fallback:
        value = config["defaults"].get(keys[-1], fallback)
        assert not value is fallback, "No fallback setting for '{}' in 'defaults' section".format(keys[-1])
    return value


def usearch_bin():
    return config.get("usearch_binary", "usearch")



def otutab_extra_files(bam, **wildcards):
    out = []
    if config["otutab"].get("extra", False):
        prefix = expand("workdir/cluster/4_otutab/{primers}/{what}", allow_missing=True, **wildcards)[0]
        out.append(prefix + "_search.txt.gz")
        if cfg_or_global_default("otutab", "program") == "vsearch":
            if bam:
                out.append(prefix + ".bam")
                out.append(prefix + ".bam.bai")
            else:
                out.append(temp(prefix + ".sam"))
    return out


def expand_clustered(path, **wildcards):
    for f in glob("results/*/*.fasta"):
        parts = f.split(os.sep)
        yield from expand(path, primers=parts[1], seqs=parts[2].split(".")[0], **wildcards)


def mem_func(mem=5, f=0, max_mem=50000):
    def _mem_func(wildcards, input, attempt):
        _mem = mem + f * input.size_mb
        return round(min(max_mem, _mem * 2 ** (attempt - 1)))

    return _mem_func


def time_func(time=1, f=0, max_time=24 * 60 * 20):
    def _time_func(wildcards, input, attempt):
        _time = time + f * input.size_mb
        return round(min(max_time, _time * 2 ** (attempt - 1)))

    return _time_func
