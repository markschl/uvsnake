#!/usr/bin/env bash

exec &>"${snakemake_log[0]}"
set -xeuo pipefail

# First, filter using max. error rate
# This is always done with VSEARCH, since there is no difference
# in how USEARCH and VSEARCH do it.
sample="${snakemake_wildcards[sample]}"
out=$({
    zstd -dcq "${snakemake_input[0]}" |  # decompress input file
        st set -i "$sample.{seq_num}" --fq | # Create ID from sample name + sequence number
        st del -d --fq |                 # delete header descriptions
        vsearch -fastq_filter - \
            --fastaout "${snakemake_output[filtered_tmp]}" \
            --fastaout_discarded "${snakemake_output[discarded_tmp]}" \
            --fastq_maxee_rate "${snakemake_params[maxee_rate]}" \
            --fastq_qmax 48
} 2>&1)

# output stats file based on parsing of intercepted output of VSEARCH command
grep 'sequences kept' <<<"$out" |
    sed -E 's/ *([0-9]+)[^,]+, *([0-9]+).*/\1 \2/g' |
    tr ' ' '\t' >"${snakemake_output[stats]}"

# forward output
cat <<<"$out" 1>&2

#### de-replicate ####

# "good" uniques only used for clustering
good="${snakemake_output[good]%.zst}"
vsearch -derep_fulllength "${snakemake_output[filtered_tmp]}" \
    -sizeout \
    -output "$good"

# all uniques: used for constructing the OTU table
# (note -sizein argument, since the filtered sequences are already
# de-replicated and thus have size annotations)

cat "$good" "${snakemake_output[discarded_tmp]}" |
    vsearch -derep_fulllength - \
        -sizein -sizeout \
        -output - |
    zstd -cq >"${snakemake_output[all]}"

# compress
zstd -fq --rm "$good"
