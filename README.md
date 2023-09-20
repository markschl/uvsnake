# USEARCH/VSEARCH-based pipeline for amplicon data processing

This [Snakemake](https://snakemake.github.io) workflow aims at providing a simple way for clustering/denoising paired-end Illumina sequences with [USEARCH](https://drive5.com/usearch) and/or [VSEARCH](https://github.com/torognes/vsearch).


## Features

* *Standard pipeline*: Paired-end read merging, primer trimming, UNOISE3/UPARSE clustering and SINTAX taxonomy assignment
* Both USEARCH and VSEARCH can be used (configurable individually at every step)
* *Multiple primer combinations* can be recognized using [Cutadapt](https://cutadapt.readthedocs.io), and subsequent clustering/taxonomy assignment is done separately for every amplicon
* Quality control using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc) and [MultiQC](https://multiqc.info)
* Optimized for large projects with many samples (samples are *processed in parallel*, beneficial if running on a cluster)

## Workflow

The workflow follows the [recommendations](https://drive5.com/usearch/manual/uparse_pipeline.html) by the author of USEARCH, but with the addition of trimming primers with Cutadatapt and adding some QC

![workflow](docs/workflow.png)


## Installation (UNIX)

Download the [latest release](https://github.com/markschl/uvsnake/releases/latest) and unpack it to some directory. In a terminal, change into the directory (`cd directory`). Also make sure to install [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (along with the Conda package manager).

The pipeline is configured to use VSEARCH for all steps, but in case of using USEARCH (`program: usearch`), first obtain the software [here](https://www.drive5.com/usearch/download.html), and make sure that it is installed as `usearch` in `$PATH`. Alternatively, specify the path to USEARCH in the *usearch_binary* section of [`config.yaml`](config/config.template.yaml).


## Input files

The following illustration shows an overview of the input files, commands and output files:

![input](docs/input.png)

More details on all available options in `config.yaml` are also available as comments in the template file [config/config.yaml](config/config.template.yaml).


## Running

To run FastQC/MultiQC, UNOISE3 and UPARSE (assuming `config/config.yaml` to be in `analysis_dir`).

```sh
conda activate snakemake
./uvsnake analysis_dir quality unoise3 uparse
```

The equivalent (a bit more complicated) Snakemake command would be:

```sh
snakemake -d analysis_dir -c1 --conda --conda-prefix ~/uvsnake/conda quality unoise3 uparse
```

Apart from setting a few defaults (not all shown here), `uvsnake` can be used like `snakemake`, any number of Snakemake arguments can be added. The full command is always printed when executing. For instance, the script also assists with setting resource constraints on computer clusters and grouping small sample processing jobs into batches.

![results](docs/results.png)

Assign taxonomy with SINTAX:

```sh
./uvsnake analysis_dir sintax
```

![results](docs/taxonomy.png)

### Notes on taxonomy

The *taxonomy* rule should not be run together with the *unoise3* and/or *uparse* rules, but in a separate command *afterwards*. Taxonomic assignments are only done for clusters that are already present. If run before clustering, nothing happens.

### Cluster execution

Executing on computer clusters requires [additional settings](https://snakemake.readthedocs.io/en/stable/executing/cluster.html). Uvsnake already takes care of setting resource constraints as well as possible, but it may make sense to allow restarting failed jobs with `-T/--retries`, in case they failed because of wrong time/memory constraints. Here an example running on a SLURM cluster with a total of three tries:

```sh
./uvsnake -T 3 --slurm analysis_dir quality unoise3 uparse
```

## Example

The [test](test) directory offers more details using an example dataset and shows how this pipeline is validated.
