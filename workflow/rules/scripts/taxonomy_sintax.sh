#!/usr/bin/env bash

exec &> ${snakemake_log[0]}
set -xeuo pipefail


# determine the binary and parameters that differ between USEARCH/VSEARCH
if [ ${snakemake_params[program]} = "usearch" ]; then
  bin="${snakemake_params[usearch_bin]}"
  extra="-maxaccepts ${snakemake_params[maxaccepts]} -maxrejects ${snakemake_params[maxrejects]}"
elif [ ${snakemake_params[program]} = "vsearch" ]; then
  bin=vsearch
  extra=""
else
  echo "unknown program: ${snakemake_params[program]}"
  exit 1
fi

sintax_out="${snakemake_output[sintax]%.*}"

"$bin" \
  -sintax "${snakemake_input[fa]}" \
  -db "${snakemake_input[db]}" \
  -tabbedout "$sintax_out" \
  -strand both \
  -sintax_cutoff ${snakemake_params[confidence]} \
  -threads ${snakemake[threads]} \
  $extra \
  1>&2

# convert to Qiime-like tabular format
# TODO: confidence not extracted
printf "Feature ID\tTaxon\tConfidence\n" | gzip -c > "${snakemake_output[taxtab]}"
awk -v IFS="\t" -v OFS="\t" '{print $1, $4, ""}' "$sintax_out" |
  sed -E 's/([a-z]):([^,\t]+),/\1__\2;/g' |  # convert lineages to QIIME format
  sed -E 's/([a-z]):([^,\t]+)/\1__\2/g' |    # convert last entry
  gzip -c >> "${snakemake_output[taxtab]}"

# compress
gzip -fn "$sintax_out"
