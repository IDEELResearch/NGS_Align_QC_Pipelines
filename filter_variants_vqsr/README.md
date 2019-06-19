# Small variant (SNV, short indel) discovery in WGS data with GATK4

## Dependencies

These are available on Longleaf by default or via `module load ...`, unless otherwise noted:

* `GATK >= 4.0.3.0`
* `java >= 1.8`
* `snakemake >= 4.0` (available from Bioconda or on Longleaf)

## Overview

This pipeline implements the [GATK Variant Quality Score Recalibration (VQSR) workflow](https://software.broadinstitute.org/gatk/documentation/article?id=11084). We assume that a VCF with raw variant calls from the `call_variants_gatk4` is available as input. We focus here on SNVs only -- the approach would be more or less the same for indels but good training sets for our species of interest (_P. falciparum_ and _P. vivax_) are hard to come by.

Steps:

1. Train VQSR model on own data augmented with external reference datasets.
2. Apply classifier to own data.
3. Filter variants by VQSLOD score.

Options modifiable by the user at runtime can be supplied either at the command line (as `--config key=value`) or (better) in a YAML-formatted configuration file (as `--configfile config.yaml`) like the one below.

```
## input files
raw_vcf: all_raw.vcf.gz
reference: $PF/pf_3d7/pf_3d7.fa
targets: all_chroms.big_chunks.core_only.bed

## known and truth sets for building VQSR model
snps_vcf_TP: $IDEEL/resources/pf_crosses_v1.0/crosses_merged.PASS.vcf.gz
snps_prior_TP: 30
snps_vcf_TPFP: $IDEEL/resources/pf3k_v5/v5.1.filtered.PASS_only.sites_only.vcf.gz
snps_prior_TPFP: 15
snps_tranches: [ 100.0, 99.9, 99.0, 90.0 ]
snps_ncomponents: 3
snps_input_annotations: [ QD, MQ, MQRankSum, ReadPosRankSum, FS, SOR ]
snps_filter_level: 90.0

## resource requirements
recal_memory: 16
```

Options are as follows. Note that all file paths will be interpreted by Snakemake relative to the working directory at runtime; it may be useful to specify absolute paths to files that live outside the project directory (eg. reference genome.)

* **reference**: path to reference genome fasta file, which must have an index (`*.fai` format)
* **raw_vcf**: a totally unfiltered VCF from your workflow of choice, but needs to carry (at least) the annotations specified in next parameter
* **snps_input_annotations**: site-level annotations (in `INFO` field of VCF) to use in training classifier
* **targets**: BED file of genomic regions in which variant calling was attempted.
* **snps_vcf_TP**: a VCF file (can be sites-only, needs requested site-level annotations) of _true positive_ sites. Ideally need >10K sites but these should be filtered as stringently as possible; err on side of conservatism. Sites ascertained in pedigrees or crosses with known founders are ideal, because Mendelian rules can be used to push error rates arbitrarily low.
* **snps_prior_TP**: Phred-scaled likelihood of sites in the _true positive_ set not being true variants
* **snps_vcf_TPFP**: a VCF file _true positive + false positive_ sites. Sites discovered in population-scale resequencing projects are appropriate here (eg. Pf3K for _P. falciparum_).
* **snps_prior_TP**: Phred-scaled likelihood of sites in the _true positive + false positive_ set not being true variants. Set this to lower value than the prior on _true positives_.
* **snps_tranches**: cutpoints _x_ for tranches defined by minimum VQSLOD score that achieves sensitivity _x_% in training sets
* **snps_ncomponents**: number of components in Gaussian mixture model that is the core of the VQSR algorithm
* **snps_filter_level**: mark as `PASS` any variants with VQSLOD that achieves _x_% sensitivity in training set (lower = higher specificity, fewer positives both false and true)
* **recal_memory**: memory limit per job in gigabytes

I show briefly how to invoke the pipeline and enable dispatch of jobs on the cluster by Slurm. This is a very simple pipeline so step-by-step explanations are not provided here. Interested readers should refer to the documentation on the GATK VQSR approach and related forum posts to understand what the parameters mean.

## Outputs

A file `filtered.vcf.gz` with the input variants annotated by `ApplyVQSR`. Each is assigned a VQSLOD score and sites meeting the sensitivity-based threshhold have filter set to `PASS`. Some diagnostic plots are also generated in `snps.plots.R.pdf` that give some idea of what annotation(s) are more important for the classifier and how obviously it discriminates between putative true and false positives.

## Running the pipeline
The pipeline is invoked with Snakemake via a command like the following. Independent steps can be run in parallel via the Slurm job scheduler. The helper script `launch.py` handles this, and passes stuff like number of threads and memory requirements on to Slurm.

```
snakemake --snakefile run_vqsr.snake \
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
