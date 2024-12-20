#!/usr/bin/env bash

exec &>"${snakemake_log[0]}"
set -xeuo pipefail

# this is where temporary files go
tempdir="${snakemake_output[tempdir]}"
mkdir -p "$tempdir"

# first, extract the input files
gzip -dc "${snakemake_input[r1]}" >"$tempdir"/R1.fastq
gzip -dc "${snakemake_input[r2]}" >"$tempdir"/R2.fastq

# handle merging options that differ between USEARCH and VSEARCH
if [[ "${snakemake_params[program]}" == "vsearch" ]]; then
    # USEARCH -fastq_pctid and VSEARCH -fastq_maxdiffpct are reciprocal
    max_diffpct=$(bc <<<"100 - ${snakemake_params[overlap_ident]}")
    # also, we add -fastq_allowmergestagger (which is default in USEARCH)
    extra="--fastq_maxdiffpct $max_diffpct --fastq_allowmergestagger"
    bin=vsearch
else
    extra="-fastq_pctid ${snakemake_params[overlap_ident]}"
    # Unfortunately, -fastq_mergepairs fails with spaces in paths with usearch v11
    # Therefore we ensure that the tempdir path is relative and never absolute
    # (we know that there are no spaces in relative paths)
    tempdir=$(python3 <<<"import os; print(os.path.relpath('$tempdir'))")
fi

# perform the merging, capturing the output
output=$(
    "${snakemake_params[program]}" \
        -fastq_mergepairs "$tempdir"/R1.fastq \
        --reverse "$tempdir"/R2.fastq \
        -fastqout "$tempdir"/merged.fastq \
        -fastqout_notmerged_fwd "$tempdir"/notmerged_R1.fastq \
        -fastqout_notmerged_rev "$tempdir"/notmerged_R2.fastq \
        -fastq_maxdiffs "${snakemake_params[max_diffs]}" \
        -fastq_qmax 48 \
        -threads 1 \
        $extra \
        2>&1
)

# output stats file based on parsing of intercepted usearch output
# (fortunately, the output of USEARCH and VSEARCH is similar enough for a common parsing routine)
n=$(grep -Ei "[0-9]+ +(read )?pairs" <<<"$output" | sed -E 's/ *([0-9]+) *(read )?pairs.*/\1/ig')
n_merged=$(grep " Merged" <<<"$output" | sed -E 's/ *([0-9]+) *Merged.*/\1/g')
printf "$n\t$n_merged" >"${snakemake_output[stats]}"

# finally, compress the files
zstd -qc "$tempdir"/merged.fastq >"${snakemake_output[merged]}"
zstd -qc "$tempdir"/notmerged_R1.fastq >"${snakemake_output[nm_r1]}"
zstd -qc "$tempdir"/notmerged_R2.fastq >"${snakemake_output[nm_r2]}"

# clean up
rm "$tempdir"/R1.fastq "$tempdir"/R2.fastq \
    "$tempdir"/merged.fastq \
    "$tempdir"/notmerged_R1.fastq "$tempdir"/notmerged_R2.fastq
