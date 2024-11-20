#!/usr/bin/env bash

exec &>${snakemake_log[0]}
set -xeuo pipefail

taxtab_header="Feature ID\tTaxon\tConfidence"

# determine the binary and parameters that differ between USEARCH/VSEARCH
extra=""
if [[ "${snakemake_params[program]}" == "vsearch" ]]; then
    # sintax_random was introduced in version v2.28.0 and may be default in future version,
    # but it still makes sense to use it
    extra="-sintax_random"
else
    # USEARCH doesn't like empty input files
    if [ ! -s "${snakemake_input[fa]}" ]; then
        echo -e "$taxtab_header" | gzip -c > "${snakemake_output[taxtab]}"
        echo -n | gzip -c > "${snakemake_output[sintax]}"
        exit 0
    fi
fi

# random seed: not compatible with multithreading
if [[ ${snakemake_params[rand_seed]} > 0 ]]; then
   extra="$extra -randseed ${snakemake_params[rand_seed]}"
   snakemake[threads]=1
fi

sintax_out="${snakemake_output[sintax]%.*}"

"${snakemake_params[program]}" \
    -sintax "${snakemake_input[fa]}" \
    -db "${snakemake_input[db]}" \
    -tabbedout "$sintax_out" \
    -sintax_cutoff ${snakemake_params[confidence]} \
    -strand ${snakemake_params[strand]} \
    -threads ${snakemake[threads]} \
    $extra \
    1>&2

# convert to Qiime-like tabular format
echo -e "$taxtab_header" | gzip -c > "${snakemake_output[taxtab]}"
# TODO: confidence not extracted
awk -v IFS="\t" -v OFS="\t" '{print $1, $4, ""}' "$sintax_out" |
    sed -E 's/([a-z]):([^,\t]+),/\1__\2;/g' | # convert lineages to QIIME format
    sed -E 's/([a-z]):([^,\t]+)/\1__\2/g' |   # convert last entry
    gzip -c >>"${snakemake_output[taxtab]}"

# compress
gzip -fn "$sintax_out"
