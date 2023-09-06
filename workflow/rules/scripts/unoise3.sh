#!/usr/bin/env bash

exec &> "${snakemake_log[0]}"
set -xeuo pipefail


if [[ ${snakemake_params[program]} == "usearch" ]]; then
    zstd -dqf "${snakemake_input[0]}"
    input_uncompressed="${snakemake_input[0]%.zst}"

    "${snakemake_params[usearch_bin]}" \
      -unoise3 "$input_uncompressed" \
      -zotus "${snakemake_output[0]}" \
      -minsize "${snakemake_params[min_size]}" \
      -maxaccepts "${snakemake_params[maxaccepts]}" \
      -maxrejects "${snakemake_params[maxrejects]}"
    
    rm "$input_uncompressed"

elif [[ ${snakemake_params[program]} == "vsearch" ]]; then
    # following code from https://github.com/torognes/vsearch/pull/283
    # using 'stdbuf' to prevent mixing of the command messages
    zstd -dcqf "${snakemake_input[0]}" \
      |
      stdbuf -eL vsearch --cluster_unoise - \
          --minsize "${snakemake_params[min_size]}" \
          --sizein \
          --sizeout \
          --maxaccepts "${snakemake_params[maxaccepts]}" \
          --maxrejects "${snakemake_params[maxrejects]}" \
          --threads "${snakemake[threads]}" \
          --centroids - \
        |
        stdbuf -eL vsearch --sortbysize - --output - |
        stdbuf -eL vsearch --uchime3_denovo - \
          --nonchimeras - \
          --sizein \
          --relabel Zotu \
        |
        st upper --wrap 80 `# convert masked letters to uppercase (alternative: use --qmask none)` \
        > "${snakemake_output[0]}"
else
    echo "unknown program: ${snakemake_params[program]}"
    exit 1
fi
