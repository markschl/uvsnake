#!/usr/bin/env bash

exec &> ${snakemake_log[0]}
set -xeuo pipefail


# First, combine the filtered and de-replicated sample files,
# and de-replicate these sequences again.
# This is important, because UNOISE assumes unique sequences
# with size annotations.
# TODO: the output is sorted by size according to VSEARCH docs.
#   -> Make sure that this stays the same when updating to
#      future versions.
zstd -dcq ${snakemake_input[good]} | 
  vsearch -derep_fulllength - -sizein -sizeout -output - |
  zstd -cq > "${snakemake_output[good]}"

# Second, combine unfiltered (actually, length filtered after trimming)
# and de-replicated sequences into one file. These will be used
# for mapping against the denoised sequences to create the OTU table.
# It is important *not* to de-replicate them again here,
# otherwise part of the sample labels will be lost, leading to 
# incorrect results.
zstd -dcq ${snakemake_input[all]} |
  zstd -cq > "${snakemake_output[all]}"
