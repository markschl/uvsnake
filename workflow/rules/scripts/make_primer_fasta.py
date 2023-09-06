from utils import file_logging


def _make_fasta(primers, outfile, add_anchor):
    with open(outfile, 'w') as f:
        for d in primers:
            name = d["name"]
            for s in d['sequences']:
                if d['anchor']:
                    s = add_anchor(s)
                f.write(f">{name}\n{s}\n")


def make_primer_fasta(primers, fwd_out, rev_out, rev_rev_out):
    _make_fasta(primers["forward"], fwd_out, lambda s: '^' + s)
    _make_fasta(primers["reverse"], rev_out, lambda s: '^' + s)
    _make_fasta(primers["reverse_rev"], rev_rev_out, lambda s: s + '$')


with file_logging(snakemake.log[0]):
    make_primer_fasta(snakemake.params.primer_config,
                        snakemake.output.fwd,
                        snakemake.output.rev,
                        snakemake.output.reverse_rev)
