# Test data

The directory `test/fastq` contains example files from an Illumina NextSeq run of two uneven (staggered) mock communities (Schlegel et al. 2018).

A random subset was taken from two published read files:

```sh
# number of sequences to keep from each sample
sample_size=10000
# random seed to ensure reproducibility of sample
random_seed=20386
# sample accessions
declare -A accessions
accessions[mock1]=SRR25280624
accessions[mock2]=SRR25280620

cd test
# download
for sample in ${!accessions[@]}; do
  acc=${accessions[$sample]}
  echo "$acc > $sample"
  fasterq-dump --split-files -O fastq $acc
  # take a random subset using a fixed seed (for reproducibility)
  # Used seqtool v0.3.0 (https://github.com/markschl/seqtool)
  st sample -n $sample_size -s $random_seed fastq/"$acc"_1.fastq | 
    gzip -c > fastq/"$sample"_R1.fastq.gz
  st sample -n $sample_size -s $random_seed fastq/"$acc"_2.fastq |
    gzip -c > fastq/"$sample"_R2.fastq.gz
  rm fastq/*.fastq
done
```

## Download taxonomy database

For taxonomic assignments, we can download the UNITE database. Here, we use the Eukaryota dataset, which is more strictly filtered (doesn't contain all singletons) (https://doi.org/10.15156/BIO/2483917):

```sh
cd test

url=https://files.plutof.ut.ee/public/orig/8F/FC/8FFCC8A730E50FEEF8CFFEEFEF02A22FBCF7E02B7FD31C6649754834D2CB0E6F.tgz
wget --no-check-certificate -O unite.tar.gz $url
mkdir -p unite
tar -C unite -xzf unite.tar.gz

# We clean the taxonomy a bit, converting undefined names to 'rank__' with an
# empty name. When importing, the pipeline will simply remove those empty ranks, as
# the SINTAX does not need these names.
cat unite/*_dynamic_*.txt |
  sed -E 's/([a-z])__[^;]+?_Incertae_sedis;/\1__;/g' |  # clean intermediate ranks
  sed -E 's/s__.+?_sp$/s__/g' `# clean undefined species names` \
  > unite/tax.txt

# Subsequently, we use https://github.com/markschl/seqtool to add the taxonomy to the
# FASTA headers.
# Also, we remove poorly annotated sequences, which have no known
# order name.
# In this example, 177138 of 234479 sequences are retained
st set -d '{l:2}' -ul unite/tax.txt unite/*_dynamic_*.fasta |  # combine
  st find --exclude --desc ';o__;f__;g__;s__' |
  gzip -c > unite_refs.fasta.gz

# clean up
rm -R unite unite.tar.gz
```

## Running the denoising/clustering pipeline

These commands run all the clustering workflows (on a local computer) and compare the results. This needs to be done in the 'uvsnake' directory, so if you ran the setup code above, make sure to `cd ..`.

```sh
conda activate snakemake

# Run the UNOISE3 and UPARSE pipelines.
# To make sure that the order of ASVs does not change between runs,
# we use only one core (-c1).
./uvsnake test unoise3 uparse

# (optional) remove working directories
./uvsnake test clean

# Now, we can assign the taxonomy using the 'sintax' rule
./uvsnake test sintax

# We may also compare the results with the expected sequences using
# VSEARCH -usearch_global
# This allows us to match the known isolates with the observed OTUs
out=test/results/ITS3-KYO2...ITS4
vsearch -usearch_global $out/unoise3.fasta -db test/mock/mock_ITS2.fasta \
    -userout $out/mock_cmp.txt \
    -userfields 'query+target+id' \
    -id 0.97

## render the example Rmd (requires pandoc in PATH or RSTUDIO_PANDOC set, here for Ubuntu)
# If this doesn't work, you can still directly run the document in RStudio
export RSTUDIO_PANDOC=/usr/lib/rstudio/resources/app/bin/quarto/bin/tools
Rscript -e "rmarkdown::render('test/R_example/example.Rmd', 'github_document')"

# the following command removes everything (INCLUDING the results/ directory)
./uvsnake test clean_all
```

This figure from the `R_example` analysis shows that the mixed relative genomic DNA concentration of the isolates corresponds well to the relative read abundance in the samples.

![mock comparison](R_example/example_files/figure-gfm/unnamed-chunk-5-1.png)


## Validation using example workflows

In order to further carefully validate this software, the test data was further analyzed using example scripts from the online documentation of the different tools, currently:

* [USEARCH pipeline](https://www.drive5.com/usearch/manual/ex_miseq_its.html) for MiSeq 2x300 fungal ITS. Results were found to be slightly different due to differences in sequence order, resulting in one OTU being different in the UPARSE pipeline (explanation see [the script itself](../scripts/simple/compare.sh)).
* [VSEARCH "alternative" pipeline](https://github.com/torognes/vsearch/wiki/Alternative-VSEARCH-pipeline/c4859786f05bba35d8c306de4a3d64fea40d9dbf) slightly modified to use UNOISE3 following the [this description](https://github.com/torognes/vsearch/pull/283). Like the USEARCH approach, the VSEARCH "alternative" pipeline contains an extra step of read mapping against the OTUs to obtain the count table, using quality filtered reads in this case. The workflow from this repository maps the raw/unfiltered reads instead (with a 97% identity threshold), [as recommended by the USEARCH author](https://www.drive5.com/usearch/manual/cmd_otutab.html). In the future, this should be configurable.

The following script runs this pipeline vs. the "simple" workflows and compares the results:

```sh
# if running for the first time, do this first:
# conda env create -f scripts/simple/uvsearch_env.yaml
scripts/simple/compare.sh
```

## Reference

Schlegel, M., Queloz, V., and Sieber, T. N. (2018). The endophytic mycobiome of European ash and sycamore maple leaves – geographic patterns, host specificity and influence of ash dieback. *Frontiers in Microbiology* 9. doi: 10.3389/fmicb.2018.02345.
