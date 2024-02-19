import csv
import re
import sys
from itertools import product
import copy


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


__complements = bytes.maketrans(
    b"ATCGRYKMBVDHNSWatcgrykmbvdhnsw -", b"TAGCYRMKVBHDNSWtagcyrmkvbhdnsw -"
)


def reverse_complement(seq):
    return seq.translate(__complements)[::-1]


def parse_primer(primer_dict):
    assert (
        isinstance(primer_dict, dict) and len(primer_dict) > 0
    ), "Primers must be defined in the form: 'name: sequence1,sequence2,...'"
    name, seqs = next(iter(primer_dict.items()))
    seqs = seqs.strip()
    assert (
        name != "no_adapter"
    ), "Primers should not be named 'no_adapter', since that indicates reads having no primer"
    anchor = False
    if seqs.startswith("^"):
        anchor = True
        seqs = seqs[1:]
    seqs = [s.strip() for s in seqs.split(",")]
    return {"name": name, "sequences": seqs, "anchor": anchor}


def get_primer_combinations(config):
    primers = {"forward": [], "reverse": [], "reverse_rev": []}
    for direction in ["forward", "reverse"]:
        for p in config["primers"][direction]:
            primers[direction].append(parse_primer(p))

    # reverse complement reverse sequences
    r_rev = primers["reverse_rev"] = []
    for p in primers["reverse"]:
        p = copy.copy(p)
        p["sequences"] = [reverse_complement(seq) for seq in p["sequences"]]
        r_rev.append(p)

    # obtain primer combinations
    combinations = primers.get("combinations", "default")
    if combinations == "default":
        combinations = []
        for fwd, rev in product(primers["forward"], primers["reverse"]):
            combinations.append("{}...{}".format(fwd["name"], rev["name"]))
    else:
        assert isinstance(combinations, list)
        for c in combinations:
            s = c.split("...")
            assert len(s) == 2, (
                "Primer combinations must be in the form 'forward...reverse'. "
                "Encountered '{}' instead.".format(c)
            )
            assert s[0] in primers["forward"], "Unknown forward primer: {}".format(s[0])
            assert s[1] in primers["reverse"], "Unknown reverse primer: {}".format(s[1])
    return primers, combinations
