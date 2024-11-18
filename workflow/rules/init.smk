
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
        seqs = seqs[1:].strip()
        assert not seqs.startswith("^"), "The anchoring symbol '^' cannot be repeated"
    seqs = [s.strip() for s in seqs.split(",")]
    assert not any(s.startswith("^") for s in seqs[1:]), (
        "Primer {}: The anchoring symbol '^' can only be used at the start of a "
        "comma-delimited oligo list. Configuring anchoring for mixed oligos "
        "individually is not possible."
    ).format(name)

    return name, seqs, anchor


def get_primer_combinations(primer_cfg):
    from itertools import product
    from copy import copy

    primers = {"forward": {}, "reverse": {}, "reverse_rev": {}}
    for direction in ["forward", "reverse"]:
        for p in primer_cfg[direction]:
            name, seqs, anchor = parse_primer(p)
            primers[direction][name] = (seqs, anchor)

    # reverse complement reverse sequences
    r_rev = primers["reverse_rev"]
    for name, p in primers["reverse"].items():
        seqs, anchor = p
        seqs = [reverse_complement(seq) for seq in seqs]
        r_rev[name] = (seqs, anchor)

    # parse primer combinations
    combinations = primer_cfg.get("combinations", None)
    if combinations is None:
        combinations = []
        for fwd, rev in product(primers["forward"], primers["reverse"]):
            combinations.append("{}...{}".format(fwd, rev))
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


# Initialize primers, taxonomy and input samples if corresponding sections 
# present in configuration
# (otherwise we assume that uvsnake is used as Snakemake module and not all
# functionality may be needed)
if "primers" in config:
    config["_primers"], config["_primer_combinations"] = get_primer_combinations(config["primers"])

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

