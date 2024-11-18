
from contextlib import contextmanager
import sys
import re


def _make_fasta(primers, outfile, add_anchor):
    with open(outfile, 'w') as f:
        for name, p in primers.items():
            seqs, anchor = p
            assert re.search(r"\s", name) is None, (
                "Primer {}: primer names must not contain spaces"
            ).format(name)
            for s in seqs:
                s = s.strip()
                assert re.search(r"\s", s) is None, (
                    "Primer {}: primer sequences must not contain spaces"
                ).format(name)
                if not s:
                    continue
                if anchor:
                    s = add_anchor(s)
                f.write(f">{name}\n{s}\n")


def make_primer_fasta(primers, fwd_out, rev_out, rev_rev_out):
    _make_fasta(primers["forward"], fwd_out, lambda s: '^' + s)
    _make_fasta(primers["reverse"], rev_out, lambda s: '^' + s)
    _make_fasta(primers["reverse_rev"], rev_rev_out, lambda s: s + '$')


@contextmanager
def file_logging(f):
    with open(f, "w") as handle:
        sys.stderr = sys.stdout = handle
        yield


with file_logging(snakemake.log[0]):
    make_primer_fasta(
        snakemake.params.primer_config,
        snakemake.output.fwd,
        snakemake.output.rev,
        snakemake.output.reverse_rev
    )
