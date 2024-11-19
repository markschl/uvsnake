"""
QIIME-style to UTAX-style database conversion.
Reference: https://www.drive5.com/usearch/manual/tax_annot.html
"""

from contextlib import contextmanager
import gzip
import os
import re
import shutil
import sys

from Bio.SeqIO.FastaIO import FastaIterator, FastaWriter


def is_gzip(filename):
    with open(filename, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


@contextmanager
def fasta_reader(filename):
    _open = gzip.open if is_gzip(filename) else open
    with _open(filename, mode="rt", encoding="ascii", errors="replace") as f:
        yield FastaIterator(f)


@contextmanager
def fasta_writer(filename):
    with open(filename, "w") as out:
        yield FastaWriter(out, wrap=None)


def convert_taxdb_utax(input, output):
    reserved_chars = re.compile("[ ,:]")
    rank_pat = re.compile(r"\s*?([a-z]+)__(.*?)\s*")
    with fasta_reader(input) as records, fasta_writer(output) as out:
        for rec in records:
            s = rec.description.split(" ", 1)
            assert len(s) == 2, "Invalid taxonomy file header: {}".format(
                rec.description)
            lineage = s[1]
            # some characters have a special meaning in the UTAX format, convert to '_'
            # (including spaces in names)
            lineage = reserved_chars.sub('_', lineage)
            # split into components
            try:
                lineage = [rank_pat.fullmatch(r).groups()
                           for r in lineage.split(';')]
            except AttributeError:
                raise Exception(
                    "Not a valid QIIME-formatted lineage: {}".format(lineage))
            # remove empty ranks, since the UTAX format does not require every rank
            # in every lineage
            lineage_out = [(rank, name)
                           for rank, name in lineage if name.strip() != ""]
            # final format
            rec.id = "{id};tax={tax};".format(
                id=rec.id,
                tax=",".join(rank + ':' + name for rank, name in lineage_out)
            )
            rec.description = ""
            out.write_record(rec)


def import_taxdb(input, output, format):
    if format == "qiime":
        convert_taxdb_utax(input, output)
    else:
        assert format == "utax", f"Unknown format: {format}"
        input = os.path.abspath(os.path.expanduser(input))
        if is_gzip(input):
            with gzip.open(input, "rb") as f, open(output, "wb") as o:
                shutil.copyfileobj(f, o)
        else:
            os.symlink(input, output)


@contextmanager
def file_logging(f):
    with open(f, "w") as handle:
        sys.stderr = sys.stdout = handle
        yield


with file_logging(snakemake.log[0]):
    import_taxdb(
        snakemake.input.db,
        snakemake.output.db,
        snakemake.params.format
    )
