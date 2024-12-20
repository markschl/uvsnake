#!/usr/bin/env bash

exec &>"${snakemake_log[0]}"
set -xeuo pipefail

# we know that the output files are all in the same directory, therefore
# we can use the dirname of the first file
split_files=(${snakemake_output[by_primers]})
outdir="$(dirname "${split_files[0]}")"
mkdir -p "$outdir"

# First, recognize forward and reverse primers, which are supplied as
# FASTA files containing one or more sequences. If > 1 sequences,
# we expect a primer mixture, where all sequences start or end at the same
# position within the locus.
#
# The logs are stored in extra files, which can be parsed by MultiQC.
# TODO: MultiQC may accept --json output someday https://github.com/ewels/MultiQC/issues/1933
#
# Recognize the forward primers and write the primer name to the header
fwd_out="$outdir/__trim_fwd.fastq"
short_file="${snakemake_output[short]%.zst}"
zstd -dcq "${snakemake_input[seq]}" |
    cutadapt - \
        -g "file:${snakemake_input[fprimers]}" \
        --suffix ' fwd={name}' \
        --error-rate "${snakemake_params[max_error_rate]}" \
        --overlap "${snakemake_params[min_overlap]}" \
        --no-indels \
        -o "$fwd_out" \
        &>"${snakemake_output[fwd_log]}"

# Recognize the reverse primers and
# distribute the reads by primer combination into different files
cutadapt "$fwd_out" \
    -a "file:${snakemake_input[rprimers_rev]}" \
    --suffix ' rev={name}' \
    --minimum-length "${snakemake_params[min_length]}" \
    --too-short-output "$short_file" \
    --error-rate "${snakemake_params[max_error_rate]}" \
    --overlap "${snakemake_params[min_overlap]}" \
    --no-indels \
    2>"${snakemake_output[rev_log]}" | # use 'seqtool' to split by primer combination (names in header)
    st split --fq -o "$outdir/{attr(fwd)}...{attr(rev)}.fastq"

rm "$fwd_out"

# Compress the files
# We want one ZSTD archive for every requested primer combination,
# even the ones that don't have any reads.
shopt -s nullglob
for comb in ${snakemake_params[primer_comb]}; do
    if [ -f "$outdir/$comb.fastq" ]; then
        zstd --rm -qf "$outdir/$comb.fastq"
    else
        echo "No sequences with both forward and reverse primer ($comb) were found in sample '${snakemake_wildcards[sample]}'"
        echo -n | zstd -cq >"$outdir/$comb.fastq.zst"
    fi
done

# compress remaining files from non-specified primer combinations
zstd --rm -qf "$outdir"/*.fastq

# compress reads that are too short
if [ -e "$short_file" ]; then
    zstd -qf --rm "$short_file"
else
    echo -n | zstd -qc >"${snakemake_output[short]}"
fi
# also compress files with missing primers
# (which Cutadapt names 'no_adapter')
for f in "$outdir/"no_adapter...*.fastq; do zstd --rm -qf "$f"; done
for f in "$outdir/"*...no_adapter.fastq; do zstd --rm -qf "$f"; done

# Finally, handle trimming statistics:
# add sample name to cutadapt logs:
# this hack is needed to tell MultiQC the sample name
sed -i -E "s/(Command line parameters[^$]+$)/\1 ${snakemake_wildcards[sample]}.fastq.gz/g" "${snakemake_output[fwd_log]}"
sed -i -E "s/(Command line parameters[^$]+$)/\1 ${snakemake_wildcards[sample]}.fastq.gz/g" "${snakemake_output[rev_log]}"
# general read statistics
if grep -q "No reads processed!" "${snakemake_output[fwd_log]}"; then
    n=0
    n_trimmed_f=0
    n_trimmed_r=0
    n_long=0
else
    # parse logfiles to obtain total numbers
    extract_num() {
        sed -E 's/[^0-9]+([0-9,]+).*/\1/g' | tr -d ','
    }
    n=$(grep 'Total reads processed' "${snakemake_output[fwd_log]}" | extract_num)
    n_trimmed_f=$(grep 'Reads with adapters' "${snakemake_output[fwd_log]}" | extract_num)
    n_trimmed_r=$(grep 'Reads with adapters' "${snakemake_output[rev_log]}" | extract_num)
    n_long=$(grep 'Reads written' "${snakemake_output[rev_log]}" | extract_num)
    # TODO this reports only the reverse sequences that are long enough, not very intuitive
fi
printf "$n\t$n_trimmed_f\t$n_trimmed_r\t$n_long" >"${snakemake_output[stats]}"
