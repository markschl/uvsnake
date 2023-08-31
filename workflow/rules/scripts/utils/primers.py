

from copy import copy
from itertools import product


__complements = bytes.maketrans(
    b"ATCGRYKMBVDHNSWatcgrykmbvdhnsw -",
    b"TAGCYRMKVBHDNSWtagcyrmkvbhdnsw -"
)


def reverse_complement(seq):
    return seq.translate(__complements)[::-1]


def parse_primer(primer_dict):
    assert isinstance(primer_dict, dict) and len(primer_dict) > 0, \
        "Primers must be defined in the form: 'name: sequence1,sequence2,...'"
    name, seqs = next(iter(primer_dict.items()))
    seqs = seqs.strip()
    assert name != "no_adapter", \
        "Primers should not be named 'no_adapter', since that indicates reads having no primer"
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
        p = copy(p)
        p["sequences"] = [reverse_complement(seq) for seq in p["sequences"]]
        r_rev.append(p)

    # obtain primer combinations
    combinations = primers.get("combinations", "default")
    if combinations == "default":
        combinations = []
        for fwd, rev in product(primers["forward"], primers["reverse"]):
            combinations.append(
                "{}...{}".format(fwd["name"], rev["name"]))
    else:
        assert isinstance(combinations, list)
        for c in combinations:
            s = c.split("...")
            assert (len(s) == 2), \
                "Primer combinations must be in the form 'forward...reverse'. " \
                "Encountered '{}' instead.".format(c)
            assert (s[0] in primers["forward"]), \
                "Unknown forward primer: {}".format(s[0])
            assert (s[1] in primers["reverse"]), \
                "Unknown reverse primer: {}".format(s[1])
    return primers, combinations
