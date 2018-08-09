# Alignment of paired-end whole-genome sequence reads to *P. falciparum*

## Dependencies

These are available on Longleaf by default or via `module load ...`, unless otherwise noted:

* `java >= 1.8`
* `samtools >= 1.5`
* `bwa >= 0.7.17`
* `cutadapt >= 1.15`
* `python >= 3.5`
* `snakemake >= 4.0` (available on Bioconda or Longleaf)
* `samblaster` (available on Bioconda)
* `Picard >= 2.0`

The `Picard` tools don't actually have to be on your `$PATH`; the path to a recent version of the jar file (on Longleaf) is hard-coded in the Snakemake script.

## Overview

Steps:

0. Organize raw fastq files under a single project directory
1. Trim adapter sequences
2. Align to genome
3. Merge lane-level bam files
4. Index merged bam file

Step 0 is assumed to be done already (but see suggestions later for how to do it.) Remaining steps are organized into a Snakemake script `align_wgs.snake`. Steps 1 and 2 run independently (ie. is parallelizable) over each lane-sample combination; step 3 runs independently on each sample. Options modifiable by the user at runtime can be supplied either at the command line (as `--config key=value`) or (better) in a YAML-formatted configuration file (as `--configfile config.yaml`) like the one below.

```
reference: /proj/ideel/julianog/users/apm/genomes/pf_3d7/pf_3d7.fa

adapters: /nas02/apps/trimmomatic-0.36/Trimmomatic-0.36/adapters/TruSeq3-PE.fa
min_trimmed_length: 30

logs: jobs
tmpdir: .tmp
readbuffer: 1000000
threads: 6
memory: 12

runs: pf3k_asia_run_map.txt
aligned: aln
fastq: fastq
```

Options are as follows. Note that all file paths will be interpreted by Snakemake relative to the working directory at runtime; it may be useful to specify absolute paths to files that live outside the project directory (eg. reference genome.)

* **reference**: path to reference genome fasta file, which must be indexed for use with `bwa`.
* **adapters**: fasta file of adapter sequences used for read trimming.
* **min_trimmed_length**: omit trimmed reads (and their mates) shorter than this from alignment.
* **logs**: directory where Slurm job scheduler will write log files
* **tmpdir**: directory used by Picard for temporary files when merging and sorting bams. Should be in scratch space. Should not point to the system-wide `/tmp` directory, which tends to fill up easily.
* **readbuffer**: maximum number of reads to hold in RAM at a time when merging and sorting bams; about a million is usually reasonable. If too large relative to the memory allotted to these steps, Java will run out of heap space and/or spend most of its time on garbage collection.
* **threads**: number of threads to use per `bwa` job; also will be used as the limit for the number of Java garbage-collection threads (which is unlimited by default, sometimes leading to problems on the cluster.)
* **memory**: memory limit for `bwa` jobs in gigabytes. Scales linearly with number of threads.
* **runs**: text file with tab-separated lines containing sample ID in first column and run ID (ERR* or SRR* for ENA or SRA data, respectively, or some unique lane identifier for locally-generated data) in second column.
* **aligned**: directory where bam files will go.
* **fastq**: directory where fastq files (or symlinks to same) are expected to live.

I show briefly how to invoke the pipeline and enable dispatch of jobs on the cluster by Slurm, then move on to more detailed discussion of each step.

## Running the pipeline
The pipeline is invoked with Snakemake via a command like the following. Independent steps can be run in parallel via the Slurm job scheduler. The helper script `launch.py` handles this, and passes stuff like number of threads and memory requirements on to Slurm.

```
snakemake --snakefile align_wgs.snake \
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

The call to Snakemake is best wrapped in a shell script (example provided in `launch_wgs.sh`) for submission with Slurm.

## Step 0: Organizing read files
Briefly, raw reads are expected to be paried-end with forward and reverse reads in separate gzipped files (ie. NOT interleaved). Each sample gets its own directory, and in that directory live all the "runs" (independent units of data generation, probably lanes) for that sample. Forward read files should be named `_1.fastq.gz` and reverse reads `_2.fastq.gz`. So fastq file paths (relative to project directory) will be `{fastq_dir}/{sample_id}/{run_id}_1.fastq.gz` and `{fastq_dir}/{sample_id}/{run_id}_2.fastq.gz`, respectively. This implies that run IDs must be unique within samples; it's probably a good idea that they be globally unique to avoid potential mischief. This is certainly true for ENA and SRA data, but depending on barcoding and sample labelling it may not always be true for files generated by UNC HTSF.

The **runs** file is a tab-separated text file tuples of (sample ID, run ID) with one line per run. In the case that a sample has multiple runs, it should appear on multiple lines. For example:

```
H34	ERR019142
H34	ERR019173
H36	ERR019143
H36	ERR019166
```

Read files (or symlinks) for the above example would look like:

```
fastq/H34/ERR019142_1.fastq.gz
fastq/H34/ERR019142_2.fastq.gz
fastq/H34/ERR019173_1.fastq.gz
fastq/H34/ERR019173_2.fastq.gz
fastq/H36/ERR019143_1.fastq.gz
fastq/H36/ERR019143_2.fastq.gz
fastq/H36/ERR019166_1.fastq.gz
fastq/H36/ERR019166_2.fastq.gz
```

See companion script `rename_htsf_fastq.py` for help generating a *runs* file and compatible symlinks from pre-existing HTSF data.

## Step 1: Trim adapters
There are many tools for adapter trimming; `cutadapt` ([documentation here](http://cutadapt.readthedocs.io/en/stable/guide.html)) has good performance and will write paired-end reads to stdout so it can be used on a stream. This is a big advantage. The call to `cutadapt` looks for a standard set of adapters in the fasta-formatted **adapters** file (can use the Illumina TruSeq ones provided with `Trimmomatic`) at the 3' end of both forward and reverse reads. If I understand the standard Illumina chemistry correctly, this should be sufficient for vanilla WGS (see [Illumina's explanation here](https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html)). Default parameters for error tolerance etc. are used otherwise.

Reads shorter than **min_trimmed_length** will be thrown out, along with their mates.

Diagnostic output is written to `{fastq_dir}/{sample}/{run_id}_trimming.log`.

Note that I do NOT perform additional trimming based on quality scores. The benefits and harms of doing this are subject to long-running debates among bioinformatics people. My sense is that, for WGS on not-ancient Illumina instruments at least, it's better to be conservative and keep as much data as possible (to increase specificity of alignments) and let variant-callers sort things out downstream.

## Step 2: Alignment to reference genome
Output from `cutadapt` is piped directly to `bwa mem`. The reference genome fasta at **reference** is assumed to be already indexed. Speed is linearly proportional to number of threads; most nodes on Longleaf have 16 or less. In practice 6-12 threads gives a reasonable balance of performance and time spent waiting for resources on the cluster. Memory usage (specified with **memory**) is linearly proportional to the number of threads used, the vast majority of which is spent on the in-memory index. That index consumes about twice as many gigabytes as there are gigabases in the genome -- so for *Plasmodium* and other small genomes, 8 Gb is probably sufficient even for a large number of threads.

A read group tag is generated from the sample ID and run ID as shown below:
```
@RG	LB:{sample_id}	ID:{run_id}	SM:{sample_id}	PL:illumina"
```
Since runs are processed independently, there will be only one read group per bam file at time of alignment, named like `{aln_dir}/lane_bams/{sample_id}/{run_id}.aligned.bam`.

I use additional options `-Y` (use soft-clipping for supplementary alignments; preserves whole read sequence in bam file) and `-K100000000` (read fixed-size chunk of data regardless of number of threads, so results are reproducible.)

Alignments are piped through `samblaster` for duplicate-marking and fixing of mate tags, then to `samtools` for conversion to name-sorted bam.

> **NOTE**: Run-level bams are temporary. Once the merged bam file for a sample has been created successfully, Snakemake will automatically delete corresponding run-level bams. Peak disk usage can still be >5x the size of the raw fastq files in the worst case, if all jobs are moving at about the same speed, but final disk usage is kept to a minimum.

## Step 3: Merge bam files
The `Picard MergeSamFiles` utility is used to merge run-level bam files into sample-level files. I use this tool because I have found it slightly less annoying than equivalent utilities in `samtools`, and it will sort the final bam file at the same time as merging. Depending on the value of **readbuffer** and the total number of reads for a sample, this step can generate many temporary files (written to **tmpdir**.) Some extra options are provided to the Java call to ensure that the jobs are not too thread-hungry. Memory limit is hard-coded at a value that has worked well in the past on Longleaf.

Merged bam files (one per sample) go in `{aln_dir}/merged/{sample}.bam`.

## Step 4: Index merged bam files
Uses `Picard BuildBamIndex` rather than `samtools index` because it offers finer control of memory usage to avoid causing problems on the cluster.

Index files (one per sample) go in `{aln_dir}/merged/{sample}.bam.bai`.
