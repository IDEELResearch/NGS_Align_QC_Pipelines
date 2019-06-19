# Small variant (SNV, short indel) discovery in WGS data with GATK4

## Dependencies

These are available on Longleaf by default or via `module load ...`, unless otherwise noted:

* `GATK >= 4.0.3.0`
* `java >= 1.8`
* `bcftools >= 1.19`
* `python >= 3.5`
* `bedtools >= 2.26`
* `pybedtools` (Python package, requires command-line `bedtools` above)
* `snakemake >= 4.0` (available from Bioconda or on Longleaf)

## Overview

This pipeline implements the [GATK Best Practices Workflow for germline short variant discovery](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145). We assume that WGS reads have been aligned with the `wgs_pe_improved` pipeline from this repository, which enforces assumes a 1-to-1 relationship between samples and BAM files. Initial gVCF generation is parallelized over samples and genomic regions. Joint genotyping step is parallelized over genomic regions. A set of intervals to use as "calling targets" -- which can just be tiles across the genome -- is required.

Steps:

0. Make **interval files** defining genomic regions to be processed together
1. Call variants separately in each sample and region with `HaplotypeCaller`
2. Merge per-sample gVCFs by genomic regions with `CombineGVCFs`
3. Perform joint genotyping across gVCFs by genomic regions with `GenotypeGVCFs`
4. Apply hard filters to chunk-wise VCFs
5. Concatenate raw variant calls into genome-wide VCF
6. Concatenate hard-filtered calls into genome-wide VCF

Options modifiable by the user at runtime can be supplied either at the command line (as `--config key=value`) or (better) in a YAML-formatted configuration file (as `--configfile config.yaml`) like the one below.

```
## basic config options
ref: $PF/pf_3d7.fa
regions: all_chroms.big_chunks.core_only.bed
samples: all_bams.txt
max_alleles: 3
outdir: vcfs_gatk_joint
theta: 0.008
tmpdir: $SCRATCH/tmp
mem: 24
mode: novel

## sites to genotype if mode='known'
known_sites: $IDEEL/resources/pf3k_v5/v5.1.filtered.PASS_only.biallelic_snps.sites_only.vcf.gz

## hard-filter cutoffs
min_DP: 5
max_DP: 1000000
min_RPRS_SNV: -8.0
min_RPRS_indel: -20.0
min_QD: 2.0
max_FS_SNV: 60.0
max_SOR_SNV: 3.0
max_FS_indel: 200.0
max_SOR_indel: 10.0
min_MQRS: -12.5
min_MQ: 25.0
```

Options are as follows. Note that all file paths will be interpreted by Snakemake relative to the working directory at runtime; it may be useful to specify absolute paths to files that live outside the project directory (eg. reference genome.)

* **ref**: path to reference genome fasta file, which must be indexed for use with `bwa`.
* **regions**: BED file of genomic regions in which to attempt variant calling. Intervals with the same `name` (column 4) are grouped together and parallelization happens over those groups. If no `name` specified, all regions assumed to be part of one big group.
* **samples**: text file with list of BAM files to process, one per line, expected to have 1-to-1 relationship with samples and to be named with a globally-unique smaple ID (as in `wgs_pe_improved` pipeline.)
* **max_alleles**: maximum number of alternate alleles to genotype at each site
* **theta**: population-scaled mutation rate per base, used as prior for variant calling (not super important if lots of data)
* **outdir**: directory where all the output goes
* **tmpdir**: directory to use for temporary files, to avoid flooding system-wide `/tmp`
* **mem**: memory limit per job in gigabytes
* **mode**: either `novel`, for usual variant discovery; or `known`, to genotype at known variant sites (provided in VCF specified in parameter **known_sites**)
* **min_DP, max_DP, ...**: treshholds for hard filters, explained in [GATK guideline "Hard-filtering germline short variants"](https://software.broadinstitute.org/gatk/documentation/article?id=11069)

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

## Step 1: per-sample variant calling with `HaplotypeCaller`

Perform variant discovery via the local _de novo_ assembly method implemented in `HaplotypeCaller`. When **mode**=`novel`, the tool runs in default mode (ie. ascertains putative variants with respect to the reference genome). When **mode**=`known`, runs in `GENOTYPE_GIVEN_ALLELES` mode, which is really genotype assignment rather than variant discovery. In the latter case a VCF (which can be sites-only) must be provided to specify which alleles to type at which sites. This step is parallelized over samples and genomic regions and produces gVCF files in `{outdir}/chunks/{region}/{sample}.g.vcf.gz`.

## Step 2: merging of gVCF files

Perform a "horizontal" merge of sample-wise gVCFs into a cohort gVCF that contains all the information needed for joint genotyping of putative variants in a given genomic region. Parallelized over genomic regions. Result goes in `{outdir}/chunks/{region}/combined.g.vcf.gz` (in `novel` mode).

## Step 3: joint genotyping

Assign (diploid) genotypes for each sample at putative variant sites using `GenotypeGVCFs`. Maximum number of alternate alleles that appear in genotypes is constrained by **max_alleles**. Parallelized over genomic regions. Result goes in `{outdir}/chunks/{region}/firstpass.vcf.gz`.

## Step 4: apply hard filters

Although the hard-filter approach -- that is, specify some cutoffs for various site-level metrics of variant quality, and exclude sites that fall below any of them -- pretty clearly underperforms machine-learning methods for cleaning up raw variant calls, it is simple and fast and good for a rough guess at how many calls _might_ be spurious. This step uses `VariantFiltration` to apply both site- and genotype-level filters, and sends the result to `{outdir}/chunks/{region}/filtered.vcf.gz`.

## Steps 5,6: make genome-wide VCFs

Just concatenate results from steps 3 and 4 to get one big raw (`{outdir}/all_raw.vcf.gz`) and hard-filtered (`{outdir}/filtered.vcf.gz`) VCF, and index them.
