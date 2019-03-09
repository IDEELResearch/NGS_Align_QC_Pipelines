# QC of paired-end whole-genome sequence reads to *Plasmodium* spp

## Dependencies

These are available on Longleaf by default or via `module load ...`, unless otherwise noted:

* `java >= 1.8`
* `samtools >= 1.5`
* `python >= 3.5`
* `bedtools >= 2.26`
* `pybedtools` (Python package, requires command-line `bedtools` above)
* `snakemake >= 4.0` (available on Bioconda or Longleaf)
* `Picard >= 2.0`
* `GATK == 3.8` (`CallableLoci` not available in GATK4 yet)
* `mosdepth >= 0.2` (available on [Github](https://github.com/brentp/mosdepth))
* `R >= 3.5` (with `knitr`, `ggplot2` and `tidyverse`; all available on Longleaf)

The `Picard` and `GATK` tools don't actually have to be on your `$PATH`; the path to a recent version of the jar file (on Longleaf) is hard-coded in the Snakemake script.

## Overview

We assume that WGS reads have been aligned with the `wgs_pe_improved` pipeline from this repository. As with that pipeline, this one is controlled by a **runs** file and assumes a 1-to-1 relationship between samples and BAM files. Operations are parallelized over samples and (optionally) genomic regions.

Steps:

0. Make **interval files** defining genomic regions to be processed together
1. Do bare-minimum validation of BAM files
2. Generate alignment summary metrics with `Picard`
3. Summarize flag counts in BAM files with `samtools flagstats`
4. Use `CallableLoci` to segment BAM files by coverage and mapping quality
5. Summarize coverage in windows with `mosdepth`
6. Create report with `R/knitr`.

Options modifiable by the user at runtime can be supplied either at the command line (as `--config key=value`) or (better) in a YAML-formatted configuration file (as `--configfile config.yaml`) like the one below.

```
# config file for QC pipeline
reference: /proj/ideel/julianog/users/apm/genomes/pf_3d7/pf_3d7.fa
runs: run_map.txt
regions: $PF/all_chroms.named.bed
aligned: $SCRATCH/cambodia/aln/merged
min_depth: 4
max_depth: 250
min_base_qual: 20
min_mapping_qual: 10
threads: 4
memory: 8
outdir: $SCRATCH/cambodia/qc
windows: $PF/all_chroms.w5kb.bed
```

Options are as follows. Note that all file paths will be interpreted by Snakemake relative to the working directory at runtime; it may be useful to specify absolute paths to files that live outside the project directory (eg. reference genome.)

* **reference**: path to reference genome fasta file, which must be indexed for use with `bwa`.
* **runs**: text file with tab-separated lines containing sample ID in first column and run ID (ERR* or SRR* for ENA or SRA data, respectively, or some unique lane identifier for locally-generated data) in second column. Should be same file given to `wgs_pe_improved` pipeline.
* **regions**: BED file of genomic regions over which to do the QC, most importantly the `CallableLoci` step. Intervals with the same `name` (column 4) are grouped together and parallelization happens over those groups. If no `name` specified, all regions assumed to be part of one big group.
* **windows**: BED file of genomic regions defining windows over which coverage profiles are generated with `mosdepth`; non-overlapping 5 kb windows is usually a good compromise between resolution and information-overload in figures
* **aligned**: directory where BAM files live
* **min_depth**: minimum read depth to consider a site callable (inclusive)
* **max_depth**: maximum read depth to consider a site callable (inclusive)
* **min_base_qual**: minimum base quality to consider a site callable (inclusive)
* **min_mapping_qual**: minimum root-mean-squared mapping quality (`MAPQ`) to consider a site callable (inclusive)
* **outdir**: directory where all the output goes
* **threads**: number of threads to use per job; also will be used as the limit for the number of Java garbage-collection threads (which is unlimited by default, sometimes leading to problems on the cluster.)
* **memory**: memory limit per job in gigabytes.

I show briefly how to invoke the pipeline and enable dispatch of jobs on the cluster by Slurm, then move on to more detailed discussion of each step.

## Running the pipeline
The pipeline is invoked with Snakemake via a command like the following. Independent steps can be run in parallel via the Slurm job scheduler. The helper script `launch.py` handles this, and passes stuff like number of threads and memory requirements on to Slurm.

```
snakemake --snakefile qc_wgs.snake \
	--configfile config.yaml \
	--directory {project_dir} \
	--printshellcmds \
	--cluster $PWD/launch.py \
	-j {max_concurrent_jobs} \
	--rerun-incomplete \
	--latency-wait {filesystem_delay}
```

The option `--directory` sets the working directory for Snakemake; all relative file paths will be interpreted relative to this root, so it's best to set it explicitly even if it's just `$PWD`.

A few extra options are supplied to improve stability on the cluster. First, the `--latency-wait` option adds a delay (measured in seconds) between the end of job execution and the check for output files, to allow the cluster's file system to catch up with itself. I have found 30 seconds to be reasonable. The value of `-j` controls how many jobs are dispatched at once. The number actually running at one time will be controlled by Slurm, but it is poor form to flood the scheduler with thousands of jobs that will pend for hours/days so I usually only let Snakemake send ~1000 a time.

The call to Snakemake is best wrapped in a shell script for submission with Slurm.

## Step 1: validation of BAM files

Just call `Picard ValidateSamFile` on each BAM. The result is one file per sample called `{outdir}/{sample}.isvalid`. Assuming no problems, the contents of that file should be a single line `No errors found`.

## Step 2: validation of BAM files

Generate some summaries over paired-end alignments with `Picard CollectAlignmentSummaryMetrics`. For more details, see the []`Picard` documentation](https://broadinstitute.github.io/picard/command-line-overview.html). Briefly these include proportion of reads mapped at all; proportion of reads mapping in proper pairs; and mismatch rates. Results go in `{outdir}/{sample}.AlignmentSummaryMetrics.txt`.

## Step 3: BAM flag summary

Just call `samtools flagstat` on each BAM. Results overlap with what's provided by `Picard CollectAlignmentSummaryMetrics`, but we include them for completeness. Results go in `{outdir}/{sample}.flagstats`.

## Step 4: define "callable" sites

Use the old `GATK CallableLoci` tool to segment the genome base-by-base into these bins (straight from GATK documentation):

* `REF_N`: reference genome had an `N`
* `PASS`: base quality, mapping quality and read depth all within defined bounds
* `NO_COVERAGE`: no reads at all, even bad ones
* `LOW_COVERAGE`: read depth below defined bounds after applying filters on base quality and mapping quality
* `EXCESSIVE_COVERAGE`: read depth above defined bounds after applying filters on base quality and mapping quality
* `POOR_MAPPING_QUALITY`: too many reads failing base or mapping quality filters

If QC is parallelized over genomic regions, this step is performed separately for each chunk and then aggregated at the sample level.

Two files are produced per sample: a summary of the total number of bases in each bin (`{outdir}/{sample}.callable_summary.txt`) and a BED file in with each interval labelled as above (`{outdir}/{sample}.callable.bed`). **This latter file is the key to defining the "genomic denominator", ie. the proportion of the genome in each sample that is accessible for variant-calling.**

## Step 5: coverage summary in windows

Use `mosdepth` for (really) quick but dirty estimation of read depth in windows defined in file passed as **windows** in the configuration file. This summary is used to make coverage plots in the final report. For any downstream work it is probably better to produce more thoughtful summaries. Result goes in `{outdir}/{sample}.regions.bed.gz`.

## Step 6: create report

The RMarkddown template `wgs_qc_report.Rmd` takes the output from previous steps and makes a nice report that includes:

* list of BAM files failing validation, if any
* table of alignment summary metrics
* plot of coverage quantiles across samples
* per-sample plots of coverage profiles, plus summaries per chromosome (for nuclear chromosomes only)

This step may be a little rickety as it needs the various players in the `tidyverse` to get along nicely.
