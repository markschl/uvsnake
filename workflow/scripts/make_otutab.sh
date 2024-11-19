exec &>"${snakemake_log[0]}"
set -xeuo pipefail

# prepare optional output
extra=
map_out=
sam_out=
for f in ${snakemake_output[extra]}; do
    if [[ "$f" == *_search.txt.gz ]]; then
        map_out="${f%.gz}"
        extra="$extra --userout $map_out -userfields query+target+id -maxhits 1"
    fi
    if [[ "$f" == *.sam ]]; then
        sam_out="$f"
        extra="$extra --samout $f"
    fi
done

# run the search
if [[ "${snakemake_params[program]}" == "vsearch" ]]; then
    # Parameter differences to USEARCH: --strand plus --sizein
    # -otutab already assumes -strand plus and will fail without size annotations
    # VSEARCH -usearch_global needs to be told these things
    zstd -dcq "${snakemake_input[uniques]}" |
        vsearch \
            --usearch_global - \
            --db "${snakemake_input[otus]}" \
            --id "${snakemake_params[ident_threshold]}" \
            --maxaccepts "${snakemake_params[maxaccepts]}" \
            --maxrejects "${snakemake_params[maxrejects]}" \
            --otutabout "${snakemake_output[tab]%.gz}" \
            --biomout "${snakemake_output[biom]}" \
            --strand plus \
            --sizein \
            --notmatched "${snakemake_output[notmatched]%.zst}" \
            --threads ${snakemake[threads]} \
            $extra

else
    # usearch -usearch_global doesn't work with size annotations in v11 (although it does in v12),
    # but there is a dedicated -otutab command
    # https://www.drive5.com/usearch/manual/cmd_otutab.html
    uniques=$(mktemp "${snakemake_input[uniques]%.*.*}".XXXXXX.fasta)
    zstd -dqf "${snakemake_input[uniques]}" -o "$uniques"
    "${snakemake_params[program]}" \
        -otutab "$uniques" \
        -otus "${snakemake_input[otus]}" \
        -id "${snakemake_params[ident_threshold]}" \
        -maxaccepts "${snakemake_params[maxaccepts]}" \
        -maxrejects "${snakemake_params[maxrejects]}" \
        -otutabout "${snakemake_output[tab]%.gz}" \
        -notmatched "${snakemake_output[notmatched]%.zst}" \
        -threads ${snakemake[threads]} \
        $extra \
        1>&2
    # -biomout "${snakemake_output[biom]}" \
    rm "$uniques"

    # Convert OTU tab to BIOM
    # This is done because USEARCH -biomout can return an
    # invalid trailing comma in the sparse count table,
    # which can break import routines
    biom convert -i "${snakemake_output[tab]%.gz}" \
        -o "${snakemake_output[biom]}" \
        --table-type 'OTU table' --to-json
    # Alternative: edit the JSON directly (not well tested)
    # end="$(tail -n3 "${snakemake_output[biom]}")"
    # truncate -s -$(wc -c <<< "$end") "${snakemake_output[biom]}"
    # echo $end | sed -E 's/\], *\] *\} *$/]]}/g' >> "${snakemake_output[biom]}"
fi

# compress the output
gzip -n "${snakemake_output[tab]%.gz}"
zstd --rm -qf "${snakemake_output[notmatched]%.zst}"

if [[ -n "$map_out" ]]; then
    gzip -f "$map_out"
fi
