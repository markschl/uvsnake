from contextlib import contextmanager
from os.path import dirname, basename
import csv
from collections import defaultdict
import sys


def write_stats(merge_stats, trim_stats, filter_stats, primer_combinations, outfile, report_file):
    # read data into dict with sample names as keys
    # TODO: depends on file paths
    merge = read_stats(merge_stats, key=lambda f: basename(
        f).replace("_stats.txt", ""))
    trim = read_stats(trim_stats, key=lambda f: basename(
        f).replace("_stats.txt", ""))
    flt = read_stats(
        filter_stats,
        key=lambda f: (basename(dirname(dirname(f))),
                       basename(f).replace("_stats.txt", "")),
    )
    flt2 = defaultdict(dict)
    for k, n in flt.items():
        flt2[k[1]][k[0]] = n

    # combine and write to output
    with open(outfile, 'w') as out, open(report_file, 'w') as rep:
        writer = csv.writer(out, delimiter='\t')
        rep_writer = csv.writer(rep, delimiter='\t')
        header = [
            'sample',
            'raw',
            'merged', 'percent-merged',
            'fwd-primer', 'fwd-percent-of-merged',
            'rev-primer', 'rev-percent-of-merged',
            'long-enough', 'long-percent-of-merged']
        rep_header = [
            'sample',
            'raw',
            'merged', '%',
            'fwd-primer', '%',
            'rev-primer', '%',
            'long-enough', '%']
        for p in primer_combinations:
            header += [p, "filtered", "filtered-percent-of-trimmed"]
            rep_header += [p, "filtered", "% of trimmed"]
        writer.writerow(header)
        rep_writer.writerow(rep_header)

        for sample in merge:
            n_raw, n_merged = merge[sample]
            m2, n_fwd, n_rev, n_long_enough = trim[sample]
            assert (n_merged == m2), (
                f"Number of sequences in logfiles from read merging and primer "
                "trimming does not match: {n_merged} vs. {m2}")
            row = [
                sample,
                n_raw,
                n_merged, percent(n_merged, n_raw),
                n_fwd, percent(n_fwd, n_merged),
                n_rev, percent(n_rev, n_merged),
                n_long_enough, percent(n_long_enough, n_merged)
            ]
            f = flt2[sample]
            for p in primer_combinations:
                kept, removed = f[p]
                total = kept + removed
                row += [total, kept, percent(kept, total)]
            writer.writerow(row)
            rep_writer.writerow(row)


def read_stats(files, key):
    out = {}
    for f in files:
        with open(f, "r") as h:
            out[key(f)] = [int(n) for n in next(csv.reader(h, delimiter="\t"))]
    return out


def percent(x, y):
    return round(100 * x / y, 2) if y > 0 else 0.0


@contextmanager
def file_logging(f):
    with open(f, "w") as handle:
        sys.stderr = sys.stdout = handle


with file_logging(snakemake.log[0]):
    write_stats(
        snakemake.input.merge,
        snakemake.input.trim,
        snakemake.input.filter,
        snakemake.params.primer_combinations,
        snakemake.output.tsv,
        snakemake.output.report
    )
