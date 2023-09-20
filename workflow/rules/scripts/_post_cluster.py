from utils import file_logging
import json

from Bio.SeqIO.FastaIO import FastaIterator, FastaWriter


def post_cluster(seqs, biom, program="vsearch", usearch_bin=None, maxaccepts=1, maxrejects=1):
    with open(biom) as f:
        d = json.loads(f)
    print(d)

# with file_logging(snakemake.log[0]):
post_cluster(snakemake.input.seqs, snakemake.input.biom,
             snakemake.output.seqs, snakemake.output.map,
             **snakemake.params)
