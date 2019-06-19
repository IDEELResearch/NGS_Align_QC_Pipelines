# Estimate complexity of infection and perform haplotype deconvolution in mixed _Plasmodium_ infections with `dEploid`

## Dependencies

These are available on Longleaf by default or via `module load ...`, unless otherwise noted:

* `dEploid >= 0.6-beta` (available from [McVean lab github](https://github.com/DEploid-dev/DEploid))
* `bcftools >= 1.19`
* `python >= 3.5`
* `bedtools >= 2.26`
* `pybedtools` (Python package, requires command-line `bedtools` above)
* `vcfdo` (available [from another IDEEL repository](https://github.com/IDEELResearch/vcfdo/))
* `snakemake >= 4.0` (available from Bioconda or on Longleaf)

## Overview

This pipeline performs joint estimation of complexity of infection (COI) and haplotype deconvolution in mixed _Plasmodium_ isolates with `dEploid`. We assume no panel of reference haplotypes and use population allele frequencies as a prior for the phasing algorithm. This means that samples should be split _a priori_ into clusters with homogeneous ancestry for analysis. Input is a VCF file of nominally-diploid calls at biallelic SNVs with the `FORMAT/AD` field (read depths on each allele) filled in.

Steps:

0. Filter input VCF on quality and global population-level allele frequency (PLAF)
1. Estimate cluster-specific PLAFs in filtered VCF
2. Make per-sample VCFs as required by `dEploid`
3. Estimate COI and deconvolve (ie. phase) haplotypes with `dEploid` using cluster-specific PLAFs
4. Aggregate COI and strain proportions across samples
5. Merge all deconvolved haplotypes into a single VCF
6. Identify clonal samples (COI=1) and create additional VCF for potential use as reference panel in future analyses

Options modifiable by the user at runtime can be supplied either at the command line (as `--config key=value`) or (better) in a YAML-formatted configuration file (as `--configfile config.yaml`) like the one below.

```
vcf: good_snps.vcf.gz
exclude: exclude.bed
outdir: onecluster_PLMAF01_VQSLOD8
samples: good_samples.txt
clusters: good_clusters.txt
min_prop: 0.01
min_vqslod: 0.0
coverage: coverage_summary.txt
method: noibd
min_reads: 10
min_plmaf: 0.01
```

Options are as follows. Note that all file paths will be interpreted by Snakemake relative to the working directory at runtime; it may be useful to specify absolute paths to files that live outside the project directory (eg. reference genome.)

* **vcf**: input cohort VCF, should carry `FORMAT/AD` and `INFO/VQSLOD` annotations
* **exclude**: BED file defining regions to mask from analysis (eg. regions with dubious SNV calls)
* **samples**: text file with list of sample IDs to include, one per line
* **clusters**: text file assigning samples to clusters, whose lines are (sample ID, cluster ID) pairs
* **coverage**: text file giving approximate sequencing depth for each sample, whose lines are (sample ID, coverage depth) pairs where depth value should be rounded to nearest integer
* **min_vqslod**: fitler sites with VQSLOD score below this threshhold
* **min_reads**: filter sites with fewer than this many total reads (across whole cohort) supporting non-reference allele
* **min_plmaf**: filter sites with _global_ PLMAF below this threshhold (within-cluster PLMAFs may still be lower)
* **min_prop**: ignore haplotypes that constitute less than this proportion of a mixed isolate
* **outdir**: directory where all output goes
* **mode**: either `noibd` for `dEploid` "classic" mode; or `ibd` for IBD-aware mode [**NOT FULLY IMPLEMENTED**]

I show briefly how to invoke the pipeline and enable dispatch of jobs on the cluster by Slurm, then move on to more detailed discussion of each step.

## Running the pipeline
The pipeline is invoked with Snakemake via a command like the following. Independent steps can be run in parallel via the Slurm job scheduler. The helper script `launch.py` handles this, and passes stuff like number of threads and memory requirements on to Slurm.

```
snakemake --snakefile call_gatk4.snake \
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

## Step 0: filter input VCF

Remove sites that overlap masked regions (specified in **exclude** file); trim unused alleles; keep only biallelic sites with filter set to `PASS`; remove sites with too few reads supporting non-reference allele or _global_ PLMAF below threshhold.

## Step 1: estimate cluster-specific PLAFs

Use `vcfdo` to estimate PLAFs within each cluster of samples. These are critical input for `dEploid` phasing algorithm. Results in `{outdir}/{cluster}.plaf`.

## Step 2: make per-sample VCFs

Results in `{outdir}/sample_vcfs/{cluster}/{sample}.vcf.gz`. Note that all samples, across all clusters, have the same number of sites in the VCF, so that phased haplotypes "line up" at the end of the analysis.

## Step 3: estimate COI and phase with `dEploid`

TODO

Dispersion parameter (`-c`) is set independently per sample to reflect heterogeneity of sequencing depth -- though in practice this doesn't seem to have much effect.

## Steps 4-6: TODO
