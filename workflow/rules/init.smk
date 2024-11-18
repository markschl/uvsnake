
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

