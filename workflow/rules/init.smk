
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


def init_sintax_config(config):
    import hashlib
    required_keys = ["db", "db_format", "confidence"]
    optional_keys = ["program", "strand", "rand_seed"]
    cfg = config["sintax"]
    db2hash = config['_taxdb_id_from_path'] = {}
    hash2db = config['_taxdb_from_id'] = {}
    # fill in missing configurations and make path hashes
    for comb in config["_primer_combinations"]:
        try:
            _cfg = cfg[comb]
        except KeyError:
            _cfg = cfg[comb] = {}
        # fill in missing keys from global sintax config
        for k in required_keys + optional_keys:
            if not k in _cfg:
                try:
                    _cfg[k] = cfg[k]
                except KeyError:
                    pass
        # validate required keys
        for key in required_keys:
            assert key in _cfg, (
                "Missing key '{key}' in the 'sintax' configuration section for "
                "the primer combination '{comb}'. The key needs to be present either "
                "globally in the 'sintax' section or in the sintax/{comb} section:"
                "\n\nsintax:\n  {key}: <value>\n\nOR:\n\nsintax:\n  {comb}:\n    {key}: <value>\n"
            ).format(key=key, comb=comb)
        # make path hash
        path = _cfg["db"]
        fmt = _cfg["db_format"]
        abs_path = os.path.abspath(path)
        hash_value = hashlib.sha256(abs_path.encode()).hexdigest()
        while hash_value in hash2db and hash2db[hash_value][0] != abs_path:
            # hash collision found: make the hash unique again.
            # Collision are detected across different files
            # within the same run, but not across repeated runs with
            # changing configurations
            hash_value += "$"
        db2hash[path] = hash_value
        hash2db[hash_value] = (abs_path, fmt)


# Initialize primers, taxonomy and input samples if corresponding sections 
# present in configuration
# (otherwise we assume that uvsnake is used as Snakemake module and not all
# functionality may be needed)
if "primers" in config:
    config["_primers"], config["_primer_combinations"] = get_primer_combinations(config["primers"])
    if "sintax" in config:
        init_sintax_config(config)

if "sample_file" in config and exists(config["sample_file"]):
    # further add these settings to config (_underscore indicates that they are special)
    l = SampleList(config["sample_file"])
    config["_input"] = OrderedDict(l.samples())
    config["_sample_names"] = list(config["_input"])
    assert len(config["_sample_names"]) == len(
        set(config["_sample_names"])
    ), "Duplicate sample names found"
    config["_layout"] = l.layout

# fill in the correct USEARCH/VSEARCH binary for each section
for section in ["merge", "unoise3", "uparse", "otutab", "sintax"]:
    if section in config:
        program = config[section].get("program", None)
        if program is None:
            program = config["program"]
        program = program.strip().lower()
        if program == "usearch":
            program = config.get("usearch_binary", "usearch")
        else:
            assert program == "vsearch", "Unknown program '{}', valid are 'usearch' or 'vsearch'".format(program)
        config[section]["program"] = program


from pprint import pprint; pprint(config)
