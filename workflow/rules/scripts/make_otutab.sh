

exec &> "${snakemake_log[0]}"
set -xeuo pipefail


if [[ "${snakemake_params[program]}" == "vsearch" ]]; then
    # prepare optional parameters
    extra=""
    sam="${snakemake_params[bam_out]%.*}.sam"
    mkdir -p $(dirname $sam)
    if [ "${snakemake_params[extra]}" = "true" ]; then
        extra="--samout $sam --userout ${snakemake_output[map_out]%.gz} -userfields query+target+id -maxhits 1"
    fi

    # run the search
    # Parameter differences to USEARCH: --strand plus --sizein
    # -otutab already assumes -strand plus and will fail without size annotations
    # https://www.drive5.com/usearch/manual/cmd_otutab.html
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

    # SAM > BAM
    if [ -s "$sam" ]; then
        rm -f "${snakemake_input[otus]}".fai "${snakemake_params[bam_out]}".bai
        samtools view -T "${snakemake_input[otus]}" -b "$sam" |
            samtools sort -@ ${snakemake[threads]} > "${snakemake_params[bam_out]}"
        rm "$sam" "${snakemake_input[otus]}".fai
        samtools index "${snakemake_params[bam_out]}"
    fi

elif [[ "${snakemake_params[program]}" == "usearch" ]]; then
    uniques=$(mktemp "${snakemake_input[uniques]%.*.*}".XXXXXX.fasta)
    zstd -dqf "${snakemake_input[uniques]}" -o "$uniques"
    "${snakemake_params[usearch_bin]}" \
        -otutab "$uniques" \
        -otus "${snakemake_input[otus]}" \
        -id "${snakemake_params[ident_threshold]}" \
        -maxaccepts "${snakemake_params[maxaccepts]}" \
        -maxrejects "${snakemake_params[maxrejects]}" \
        -otutabout "${snakemake_output[tab]%.gz}" \
        -biomout "${snakemake_output[biom]}" \
        -notmatched "${snakemake_output[notmatched]%.zst}" \
        -threads ${snakemake[threads]} \
        1>&2
    rm "$uniques"

else
    echo "unknown program: ${snakemake_params[program]}"
    exit 1
fi

# compress the output
gzip -n "${snakemake_output[tab]%.gz}"
zstd --rm -qf "${snakemake_output[notmatched]%.zst}"

if [ "${snakemake_params[extra]}" = "true" ]; then
    "${snakemake_output[map_out]%.gz}"
fi
