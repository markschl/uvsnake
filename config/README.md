# Configuration

The configuration lives in the `config` directory and is named `config.yaml`. For a simple project, it may look like this:

```yaml
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

The command `snakedeploy deploy-workflow https://github.com/markschl/uvsnake . --tag v0.1.1` will copy a configuration template to `config/config.yaml`. All its sections can be modified according to your needs. In the following, we describe the available options.

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

## Program (USEARCH and VSEARCH)

The pipeline is configured to use the [open-sourced version USEARCH (v12 beta)](https://github.com/rcedgar/usearch12) for all steps, with the option to use your own preferred USEARCH binary (`usearch_binary` setting), which you may [obtain here](https://github.com/rcedgar/usearch_old_binaries).

Alternatively, VSEARCH can be used for all steps using this setting:


```yaml
program: vsearch
(...)
```

Or, you can configure the program (`usearch` or `vsearch`) at each step individually (see [all options](#all-options)) below.

### Which program to use?

[USEARCH](https://www.drive5.com/usearch) has existed for years as a closed-source binary distributed a free 32-bit version (handles files < 4 GiB) and a paid 64-bit version. [VSEARCH](https://github.com/torognes/vsearch) emerged as an open-source counterpart that is very similar to use and produces similar results. In 2024, an USEARCH version has been [released as open source](https://github.com/rcedgar/usearch12) on GitHub (with most of the core functionality) and the older binaries were [released to the publich domain](https://github.com/rcedgar/usearch_old_binaries). Before, the 4 GiB restriction affected *very* large projects, primarily at OTU table construction step. In USEARCH workflows [handling all samples at once](https://www.drive5.com/usearch/manual/ex_miseq_its.bash), medium-sized projects would hit the limit, too. But this pipeline starts with processing samples individually, so this has not been a problem. With USEARCH v12, these considerations are gone anyway, so the choice of the program is completely up to the user.

Here we list some technical details on every step that might make it easier to decide:

|  | USEARCH | VSEARCH |
|---|---|---|
| Read merging | [Edgar & Flyvbjerg (2015)](https://doi.org/10.1093/bioinformatics/btv401); [fastq_mergepairs](https://rcedgar.github.io/usearch12_documentation/cmd_fastq_mergepairs.html) command | More conservative: [more reads of potentially bad quality remain unmerged](https://groups.google.com/g/vsearch-forum/c/ojqZYbEyUw4/m/9RMhPQbXAwAJ) with the message "alignment score too low, or score drop too high" regardless of the options specified (overlap_ident, max_diffs). The speed is comparable with single-threaded (per-sample) processing. |
| Primer trimming | ([Cudatapt](https://cutadapt.readthedocs.io/en/stable/guide.html) is used) |  |
| Quality filtering | Currently always done with VSEARCH. The results are exactly the same between the two. |  |
| UNOISE3 clustering | See [Edgar & Flyvbjerg (2015)](https://doi.org/10.1093/bioinformatics/btv401) and [Edgar (2016)](https://doi.org/10.1101/081257). The [unoise3 command](https://rcedgar.github.io/usearch12_documentation/cmd_unoise3.html) does clustering and chimera checking in a single step. | The clustering result is very similar, but [chimera checking is done in a separate step](https://github.com/torognes/vsearch/pull/283) instead of being integrated into the command. The runtime is faster and the command runs multi-threaded. |
| UPARSE clustering | [Edgar (2013)](https://doi.org/10.1038/nmeth.2604); the [cluster_otus](https://rcedgar.github.io/usearch12_documentation/cmd_cluster_otus.html) function does fixed-threshold clustering at 97% and chimera checking in one command. | (not implemented) |
| OTU table construction by mapping the raw (unfiltered) reads against the OTU table at 97% similarity (default) | Using the [otutab](https://drive5.com/usearch/manual/cmd_otutab.html) command. A [heuristic method](https://rcedgar.github.io/usearch12_documentation/aln_heuristics.html) is used to approximately find the best alignment fast unless `-fulldp` is specified (not used in uvsnake, since it makes USEARCH dramatically slower). The `maxaccepts` and `maxrejects` options are important parameters for [balancing speed vs. sensitivity](https://rcedgar.github.io/usearch12_documentation/termination_options.html). We use the rather high and conservative defaults of USEARCH v11 (8 and 256) at the expense of speed. Lowering `maxrejects` speeds up the search. | Using the `usearch_global` command with parameters equivalent to `otutab`. VSEARCH always calculates the optimal alignment (equivalent to USEARCH's `-fulldp`). As a consequence, constructing the OTU table by sequence search is roughly ~4 x slower in our experience. To speed up searching, it is necessary to lower the value for `maxrejects` to < 256. |
| SINTAX taxonomy assignment | [Edgar (2016)](https://doi.org/10.1101/074161), using the [sintax command](https://rcedgar.github.io/usearch12_documentation/cmd_sintax.html). | Implemented according to Edgar (2016), but the exact implementation [may differ slightly](https://github.com/torognes/vsearch/issues/535#issuecomment-2079405760). VSEARCH usually runs faster. |


## All options

In the following, the structure of `config/config.yaml` is shown along with detailed comments. The following sections are mandatory (described below): `primers`, `merge`, `filter`, `otutab`. The following sections are optional: `program`, `usearch_binary`, `input`, `uparse`, `unoise3` and `sintax`.

*uparse*, *unoise3* and *sintax* are needed if the corresponding target rule is actually run.

```yaml
# (optional) 
# Program to use by default for read merging, clustering and OTU table construction 
# This default can still be overridden by the individual program settings
program: usearch  # {usearch, vsearch}

# (optional)
# Path to custom USEARCH binary in case you don't want to use v12 beta
# (download here: https://github.com/rcedgar/usearch_old_binaries)
usearch_binary: usearch

# Sample file, which is a three-column file with the following header:
# id  R1  R2
# The R1/R2 sample paths must be absolute, or relative to the 
# target directory where the results go (not the Snakemake directory).
# (QIIME2 manifest files are accepted, too)
input:
  sample_file: manifest_paired_run1.tsv

# Primers: specify one or more forward and reverse primers.
# All combinations are further considered, unless restricted in `combinations`.
# Clustering is done separately for every combination.
primers:
  forward:
    # single sequence or list of mixed oligos
    - fwd_primer1: SEQUENCE
    # Here, the second primer is a mix of two phased variants
    # (N inserted before second sequence), and the variants are anchored (^)
    - fwd_primer2: ^ SEQUENCE, NSEQUENCE
  reverse:
    - rev_primer1: SEQUENCE
    - rev_primer2: SEQUENCE
  combinations:
      # (optional) list of primer combinations
      - fwd_primer1...rev_primer1
      - fwd_primer2...rev_primer2
  # Primer trimming settings
  trim_settings:
    min_overlap: 15
    # rate of substitutions/InDels allowed in the alignment
    # (also allows ^anchored sequences to be shifted to a certain extent)
    max_error_rate: 0.15
    # length filter: discard shorter amplicons (may be adapters)
    min_length: 100

# Paired-end read merging
merge:
  overlap_ident: 75  # percent
  # max. nucleotide differences
  # here, we don't care about this (setting the number very high),
  # we only rely on overlap identity
  max_diffs: 1000
  # (optional) other settings, overriding above defaults
  # Note: VSEARCH is more conservative (less reads merged)
  program: usearch  # {usearch, vsearch}

# Read filtering before denoising/clustering
# (currently *always* done by VSEARCH)
filter:
  # Maximum per-base error rate
  # (e.g. 0.002 per bp means 0.8 errors per 400 bp)
  max_error_rate: 0.002

# UPARSE options (https://doi.org/10.1038/nmeth.2604)
# https://www.drive5.com/usearch/manual/cmd_cluster_otus.html
# The clustering threshold is fixed at 97%
uparse:
  # min. OTU size (USEARCH has 2 as default, discarding singletons)
  min_size: 2
  # (optional) custom termination options overriding the defaults
  maxaccepts: 1
  maxrejects: 32

# UNOISE3 settings
# https://www.drive5.com/usearch/manual/cmd_unoise3.html
unoise3:
  # min. zOTU size (USEARCH has 8 as default)
  min_size: 8
  # (optional settings)
  program: usearch  # {usearch, vsearch}
  maxaccepts: 1
  maxrejects: 32

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
  # (optional settings)
  program: usearch  # {usearch, vsearch}
  maxaccepts: 8
  maxrejects: 256

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

  # (optional settings)
  program: vsearch
  # Which strand to search: use 'plus' if all database sequences have 
  # the same orientation
  strand: both  # {both, plus}
  # Random number generator seed: initialize with a seed for reproducible outcome
  # (but only a single thread is used, which is slower)
  rand_seed: 0  # 0 means no random seed

  # Override global SINTAX settings for individual primer combinations
  # (here: confidence and rand_seed, the database 'db.fasta.gz' is still used)
  - fwd_primer1...rev_primer1:
    confidence: 0.8
    rand_seed: 42

  # Use another database file for second primer combination
  # (overriding the global setting, which is only used if there is
  # no specific setting for a primer combination)
  - fwd_primer2...rev_primer2:
    db: db2_utax.fasta.gz
    db_format: utax
```
