# frequently used imports made available in all rules
import sys
import os
from os.path import dirname, exists
from glob import glob
from collections import OrderedDict


__complements = bytes.maketrans(
    b"ATCGRYKMBVDHNSWatcgrykmbvdhnsw -",
    b"TAGCYRMKVBHDNSWtagcyrmkvbhdnsw -"
)


def reverse_complement(seq):
    return seq.translate(__complements)[::-1]


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


def expand_clustered(path, **wildcards):
    for f in glob("results/*/*.fasta"):
        parts = f.split(os.sep)
        yield from expand(path, primers=parts[1], seqs=parts[2].split(".")[0], **wildcards)


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


class SampleList(object):
    default_header = {"paired": ["id", "R1", "R2"], "single": ["id", "R1"]}
    qiime_header = {
        "paired": [
            "sample-id",
            "forward-absolute-filepath",
            "reverse-absolute-filepath",
        ],
        "single": ["sample-id", "absolute-filepath"],
    }

    def __init__(self, sample_file=None, layout=None, reserved_chars=None):
        self._samples = []
        if reserved_chars is None:
            self.reserved_re = None
        else:
            import re
            self.reserved_re = re.compile(f"[{reserved_chars}]")
        if sample_file is not None:
            assert layout is None
            self._read_samples(sample_file)
        else:
            assert layout in ("single", "paired")
            self.layout = layout
            self.ncol = 3 if layout == "paired" else 2
            self.header = None
        self.n_reads = self.ncol - 1

    def _read_samples(self, sample_file):
        import csv
        with open(sample_file) as f:
            rdr = csv.reader(f, delimiter="\t")
            self.header = next(rdr)
            self.ncol = len(self.header)
            if self.ncol == 3:
                self.layout = "paired"
            elif self.ncol == 2:
                self.layout = "single"
            else:
                raise AssertionError(
                    "Invalid number of columns in sample file. "
                    "Either two (single-end) or three (paired-end) are expected"
                )
            assert (
                self.header == self.default_header[self.layout]
                or self.header == self.qiime_header[self.layout]
            ), (
                "Unknown {}-end sample file header: '{}'. " "Valid are '{}' or '{}'"
            ).format(
                self.layout,
                ",".join(self.header),
                ",".join(self.default_header[self.layout]),
                ",".join(self.qiime_header[self.layout]),
            )
            for row in rdr:
                self.add(row[0], row[1:])

    def add(self, sample, reads):
        row = [sample] + list(reads)
        assert len(row) == self.ncol
        if self.reserved_re is not None:
            (row[0], n) = self.reserved_re.subn("_", row[0])
            if n > 0:
                print(
                    f"Reserved characters replaced in sample name: {row[0]}",
                    file=sys.stderr,
                )
        self._samples.append(row)

    def samples(self):
        for row in self._samples:
            yield row[0], row[1:]

