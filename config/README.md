# Configuration

The configuration lives in the `config` directory and is named `config.yaml`. It roughly looks like this:

```yaml
defaults:
  program: usearch
  maxaccepts: 1
  maxrejects: 1

primers:
  forward:
    - fwd_name: SEQUENCE
  reverse:
    - rev_name: SEQUENCE
  trim_settings:
    min_overlap: 15
    max_error_rate: 0.1
    min_length: 100

merge:
  overlap_ident: 75
  max_diffs: 1000

filter:
  max_error_rate: 0.002

uparse:
  min_size: 2

unoise3:
  min_size: 8

otutab:
  ident_threshold: 97

sintax:
  db: taxonomic_database.fasta.gz
  db_format: utax
  confidence: 0.8
```

The command `snakedeploy deploy-workflow https://github.com/markschl/uvsnake . --tag v0.1.1` will copy a simple configuration template to `config/config.yaml`. All its sections can be modified according to your needs. In the following, we describe the available options.

## Primers

The `primers` section needs to contain the correct primer sequences:

```yaml
primers:
  forward:
    - fwd_name: SEQUENCE
  reverse:
    - rev_name: SEQUENCE
```

Mixes of primer oligos and multiple primer combinations are possible, in which case all possible combinations are searched (see also [below](#all-options)):

```yaml
primers:
  forward:
    - fwd_1: SEQUENCE
    - fwd_2: SEQUENCE
  reverse:
    - rev_1: SEQUENCE
    - rev_2: SEQUENCE
```

To restrict the combinations, use `combinations`, with the names separated by `...`:

```yaml
primers:
  (...)
  combinations:
    - fwd_1...rev_1
    - fwd_2...rev_2
```

### Anchoring

To make sure that primers are only recognized at the start of the sequence, insert the anchoring symbol *`^`* at the start: `^SEQUENCE` (for both forward and reverse primers). This can be especially useful if different primers were used to amplify overlapping portions of the same locus, an therefore some primer sequences may be wrongly recognized at some internal position.

With anchoring, base insertions or deletions at the start are still possible up to the given error rate (`max_error_rate` setting) (see also [Cutadapt documentation](https://cutadapt.readthedocs.io/en/stable/guide.html#anchored-5adapters)).


### Primer mixes

Forward and/or reverse primers may actually be mixes of multiple oligos, which can be represented as a comma-delimited list: `SEQUENCE1, SEQUENCE2, ...`.

Possible use cases:

- making sure that specific taxonomic groups are well amplified, but adding degeneracies would introduce too many unnecessary primer sequence variants
- phased nucleotide inserts to increase nucleotide diversity for improved Illumina sequence quality

Specifying all phased primer variants is not actually necessary, unless [anchoring](#anchoring) is used. If anchoring, specify the anchoring symbol only once at the start. Example:

```yaml
primers:
  forward:
    - phased_primer: 
       ^   PRIMER,
          NPRIMER,
         NNPRIMER,
        NNNPRIMER
  (...)
```

> Note: spaces and line breaks are optional but can be useful for clarity


## Sample files

The sample file `config/samples.tsv` should contain a list of paths of all demultiplexed FASTQ files files to be processed. In the following, we use the tool [make_sample_tab](https://github.com/markschl/ngs-sample-tab), which was created exactly for this purpose:

```sh
make_sample_tab -d path/to/fastq_files/run1 -f simple
```

The message indicates that paired-end Illumina reads were found in the `run1` directory:

```
Automatically inferred sample pattern: 'illumina'
120 samples from 'run1' (paired-end) written to samples_run1_paired.tsv
```

The content of `samples_run1_paired.tsv` may look like this:

```
id	R1	R2
sample1	path/to/fastq_files/run1/sample1_R1.fastq.gz	path/to/fastq_files/run1/sample1_R2.fastq.gz
sample2	path/to/fastq_files/run1/sample2_R1.fastq.gz	path/to/fastq_files/run1/sample2_R2.fastq.gz
(...)
```

We rename the file to `config/samples.tsv`:

```sh
mv samples_run1_paired.tsv config/samples.tsv
```

Alternatively, specify a custom location for the sample file in `config/config.yaml`:

```yaml
sample_file: samples_run1_paired.tsv
```

## USEARCH

The pipeline is configured to use VSEARCH for all steps, but in case of using USEARCH (`program: usearch`), first obtain the software [here](https://www.drive5.com/usearch/download.html), and make sure that it is installed as `usearch` in `$PATH`. Alternatively, specify `usearch_binary: path/to/usearch`.

## All options

In the following, the structure of `config/config.yaml` is shown along with detailed comments. The following sections are mandatory (described below): `primers`, `merge`, `filter`, `otutab`. The following sections are optional: `usearch_binary`, `sample_file`, `uparse`, `unoise3` and `sintax`.

*uparse*, *unoise3* and *sintax* are needed if the corresponding target rule is actually run.

```yaml
defaults:
  # Program to use by default for read merging, denoising and OTU table construction 
  # The output should still be similar between USEARCH and VSEARCH, except for
  # the read merging, where VSEARCH is more conservative (see comment below)
  # This default can still be overridden by the individual program settings
  # (set them to 'usearch' or 'vsearch' instead of 'default')
  program: vsearch  # {vsearch, usearch}
  # termination options
  # https://drive5.com/usearch/manual/termination_options.html
  # Their values have a strong impact on speed and sensitivity
  # * setting both to 0 searches the whole database
  # * setting to values > 1 results in a more precise, but slower search/clustering
  # These are the default settings by USEARCH and VSEARCH
  # the unoise3, uparse, post_cluster, otutab and sintax commands all
  # accept these settings as well, which overrides these defaults
  maxaccepts: 1
  maxrejects: 1

# Default search settings

# Path to USEARCH binary
# (download here: https://drive5.com/usearch/)
# By default, we assume that binary is in PATH as 'usearch'
usearch_binary: usearch

# Sample file, which is a three-column file with the following header:
# id  R1  R2
# The R1/R2 sample paths must be absolute, or relative to the 
# target directory where the results go (not the Snakemake directory).
# (QIIME2 manifest files are accepted, too)
input:
  sample_file: manifest_paired_run1.tsv

# primers/primer combinations
primers:
  forward:
    # this is a list with one or several primers (depending on what you have in the run)
    - fwd_primer1: SEQUENCE1
    # Instead or in addition to using degenerate primers, there can also
    # be a mix of different oligos, which all start at the same position
    # in the locus
    - fwd_primer2: SEQUENCE2, SEQUENCE3
  reverse:
    - rev_primer1: SEQUENCE4
     # (optional) list of primer combinations
     # (if not supplied, all combinations are used)
     # Every combination will obtain its own results directory
  combinations:
      - ITS3-KYO2...ITS4
      - fITS7...ITS4
  # Primer trimming settings
  trim_settings:
    min_overlap: 15
    max_error_rate: 0.15
    # keep only amplicons of at least this length (shorter may be adapters)
    min_length: 100

# paired end merging
# Note on 'program' setting:
# VSEARCH is more conservative regarding more reads of potentially
# bad quality will remain unmerged with the message
# "alignment score too low, or score drop too high" regardless of the options
# specified (overlap_ident, max_diffs)
# https://groups.google.com/g/vsearch-forum/c/ojqZYbEyUw4/m/9RMhPQbXAwAJ
merge:
  overlap_ident: 75  # percent
  # max. nucleotide differences
  # here, we don't care about this (setting the number very high),
  # we only rely on overlap identity
  max_diffs: 1000
  # (optional) other settings, overriding above defaults
  # program: usearch

# read filtering before denoising/clustering
# (done by VSEARCH)
filter:
  # Maximum per-base error rate
  # (e.g. 0.002 per bp means 0.8 errors per 400 bp)
  max_error_rate: 0.002

# UPARSE options (https://doi.org/10.1038/nmeth.2604)
# USEARCH command: https://www.drive5.com/usearch/manual/cmd_cluster_otus.html
# The clustering threshold is fixed at 97%
# (uparse target rule, creates 'otus.fasta')
uparse:
  # min. OTU size (USEARCH has 2 as default, discarding singletons)
  min_size: 2
  # (optional) other settings, overriding above defaults
  # program: usearch
  # maxaccepts: 8
  # maxrejects: 8

# UNOISE3 settings
# USEARCH command: https://www.drive5.com/usearch/manual/cmd_unoise3.html
# (unoise3 target rule, creates 'denoised.fasta')
unoise3:
  # min. zOTU size (USEARCH has 8 as default)
  min_size: 8
  # (optional) other settings, overriding above defaults
  # program: usearch
  # maxaccepts: 8
  # maxrejects: 8

# Mapping options for creation of OTU table
# for USEARCH: https://www.drive5.com/usearch/manual/cmd_otutab.html
# with VSEARCH, we use the standard -usearch_global
otutab:
  # similarity threshold (USEARCH -otutab uses 97% as default)
  ident_threshold: 97
  # extra: true will create two additional output files (may be large!):
  # * BAM file for inspection of mapped reads (VSEARCH only)
  # * A TSV file (named ..._search.txt.gz) mapping each read to each cluster.
  # They are placed in workdir/otutab_mapping
  extra: true
  # (optional) other settings, overriding above defaults
  # program: usearch
  maxaccepts: 64
  maxrejects: 64

# Taxonomy assignment options
sintax:
  # path to database (can be GZIP-compressed)
  db: db.fasta.gz
  # database format: (example has only kingdom and order defined)
  # * QIIME-style: >id k__Kingdom;p__;o__Order;f__;g__;s__;
  # * UTAX-style headers: >id;tax=k:Kingdom,o:Order;
  db_format: qiime  # {qiime, utax}
  # bootstrap threshold
  confidence: 0.9
  # (optional) other settings, overriding above defaults
  # program: vsearch

  # confidence setting override for first primer combination
  # (database 'db.fasta.gz' is still used)
  - fwd_primer1...rev_primer1:
    confidence: 0.8
    # random seed to ensure reproducibility (but only single thread is used, which is slower)
    rand_seed: 42

  # use another database file for second primer combination
  # (overriding the global setting, which is only used if there is
  # no specific setting for a primer combination)
  - fwd_primer2...rev_primer2:
    db: db2_utax.fasta.gz
    db_format: utax
```
