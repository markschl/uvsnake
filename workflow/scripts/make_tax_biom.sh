#!/usr/bin/env bash

exec &>"${snakemake_log[0]}"
set -xeuo pipefail

mkdir -p "$(dirname "${snakemake_output[tax_tmp]}")"

# write header
gzip -dc "${snakemake_input[tax]}" |
    sed 's/Taxon/taxonomy/g' |
    sed 's/Feature ID/# Feature ID/g' >"${snakemake_output[tax_tmp]}"

# write rest
if [[ $(wc -l <"${snakemake_output[tax_tmp]}") -ge 2 ]]; then
    biom add-metadata -i "${snakemake_input[biom]}" \
        -o /dev/stdout \
        --observation-metadata-fp "${snakemake_output[tax_tmp]}" \
        --sc-separated taxonomy --float-fields Confidence --output-as-json |
        gzip -nc >"${snakemake_output[biom]}"
else
    # no taxa
    echo -n | gzip -nc >"${snakemake_output[biom]}"
fi
