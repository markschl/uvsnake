#!/usr/bin/bash

# This script compares the test pipeline results with
# simple example pipelines published by the pipeline authors
# in order to validate our code

set -euo pipefail

# In order for results to be reproducible, we use only one core.
# This is especially important for USEARCH, where the order of 
# merged reads changes with multithreading (apparently not with VSEARCH).
threads=1

# this relies on the 'uvsearch' env
(conda env create -qf scripts/simple/uvsearch_env.yaml || true)

source $(conda info --base)/etc/profile.d/conda.sh

# now prepare input files
fq=mock_example/fastq
out=mock_example/simple
rm -rf $out
mkdir -p $out/fq
for f in $fq/*.fastq.gz; do
  zcat $f > $out/fq/$(basename ${f%.gz})
done
# primer mixes that were used for the specific mock communities in the test dataset
fprimers=$out/fprimers.fa
echo -n > $fprimers
for seq in GATGAAGAACGYAGYRAA GATGAAGRACGCMGCGAA GATGACGAACGCAGCGAA GATGAAGAACACAGYGAA GATGAAGAGCGYAGCRAA GATGAAGAACGCGGCGAA TCGATGAAGAMCGTWGC; do
  printf ">ITS3-KYO2\n$seq\n" >> $fprimers
done
rprimers=$out/rprimers.fa  # these are reverse complemented
echo -n > $rprimers
for seq in GCATATCAATAAGCGGAGGA GCATATTAATAAGCGGAGGA GCATATCABTAAGCGGAGGA GCATATCAATAAGTGGAGGA GCATATCAATAAGCAGAGGA GCATATCAATAAGCGCAGGA GCATATTAWTCAGCGGAGGA GCATATCAATAAGGCGAGGA; do
  printf ">ITS3-KYO2\n$seq\n" >> $rprimers
done

# run the pipeline once using VSEARCH
conda activate snakemake
./uvsnake -d mock_example clean_all
./uvsnake -d mock_example unoise3 --config 'defaults={program: vsearch}'

# Run "simple" VSEARCH analysis and compare results
pout=mock_example/results/ITS3-KYO2...ITS4
vout=mock_example/simple/vsearch
conda activate uvsearch
scripts/simple/vsearch.sh $vout mock_example/config/config.yaml file:$fprimers file:$rprimers $threads $out/fq/*_R1.fastq &> mock_example/simple/vsearch.log
if ! cmp -s $vout/unoise3.fasta $pout/unoise3.fasta; then
  echo "VSEARCH unoise3 otutabs differ ($vout/unoise3_otutab.txt.gz $pout/unoise3_otutab.txt.gz)" >&2
  exit 1
fi
if ! cmp -s $vout/unoise3_otutab.txt.gz $pout/unoise3_otutab.txt.gz; then
  echo "VSEARCH zOTUs differ ($vout/unoise3_otutab.txt.gz $pout/unoise3_otutab.txt.gz)" >&2
  exit 1
fi
echo "VSEARCH 'simple' pipeline does not differ from UVSnake"

# for comparing biom with Meld/json.tool:
# meld <(python -m json.tool mock_example/simple/vsearch/unoise3.biom) <(python -m json.tool mock_example/results/ITS3-KYO2...ITS4/unoise3.biom)
# ...or the whole directories
# meld mock_example/results/ITS3-KYO2...ITS4 mock_example/simple/vsearch

# Run USEARCH pipeline
conda activate snakemake
./uvsnake -d mock_example clean_all
./uvsnake -d mock_example unoise3 uparse

# Run simple USEARCH analysis
conda activate uvsearch
uout=mock_example/simple/usearch
scripts/simple/usearch.sh $uout mock_example/config/config.yaml file:$fprimers file:$rprimers $threads $out/fq/*_R1.fastq &> mock_example/simple/usearch.log

# for USEARCH, direct comparisons are more difficult. This has to do with the de-replication,
# which our pipeline does separately per sample before collecting the unique sequences
# from all samples and de-replicating again. The standard USEARCH pipeline does it for all 
# samples together at once. The output is sorted by size, but in case of ties (same size),
# the order will vary between the two strategies. Since the output of greedy U/VSEARCH
# clustering algorithms depends on the order of the input sequences, the results can be different.
# (also not sure if -threads setting is always respected by USEARCH, which may lead to reordering)
#
# Here we compare just the sorted sequences and OTU counts. Small differences are still
# not a 100% proof for a bug in our code.
#
# Needs 'seqtool' is installed
sort_seqs() {
  st . --to-tsv seq $1 | sort
}
sort_freqs() {
  zcat $1 | cut -d$'\t' -f2,3 | sort
}
for what in unoise3 uparse; do
  if ! cmp -s <(sort_seqs $uout/$what.fasta) <(sort_seqs $pout/$what.fasta); then
    echo "USEARCH $what OTUs differ: ($uout/$what.fasta $pout/$what.fasta)" >&2
    # exit 1
  fi
  if ! cmp -s <(sort_freqs $uout/"$what"_otutab.txt.gz) <(sort_freqs $pout/"$what"_otutab.txt.gz); then
    echo "USEARCH $what otutabs differ: ($uout/"$what"_otutab.txt.gz $pout/"$what"_otutab.txt.gz)" >&2
    # exit 1
  fi
done

# for comparing biom with Meld:
# meld <(python -m json.tool mock_example/simple/usearch/unoise3.biom) <(python -m json.tool mock_example/results/ITS3-KYO2...ITS4/unoise3.biom)
# ...or the whole directories
# meld mock_example/results/ITS3-KYO2...ITS4 mock_example/simple/usearch
