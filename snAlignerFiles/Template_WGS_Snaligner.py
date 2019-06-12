###############################################################################
# Purpose:  SnakeMake File to take FASTQs and make BAMs
# Authors: Nick Brazeau (inspired by Christian Parobek)
# NOTES: Updated on June 10, 2019 to incorporate Andrew Morgan's
# BWA set seed and read chunking (for reproducibility) and use of samblaster
#Given: FASTQ
#Return: BAM
#Remarks: Note, you will have to change the paths and project specifics as well as the tools paths (depending on user) for each project/each new directory
#         You will need to add the BamMergeCommander to your path
#         Input for BamMergeCommander is a tab delimited file that first column contains a header and sample names or sample IDs (some unique n.alphanumeric for your sample)
#Remarks: PLEASE DO ALL ALIGNMENTS IN YOUR SCRATCH SPACE AND THEN TRANSFER OVER THE FILE PRODUCT TO YOUR ANALYSIS DIRECTORY
###############################################################################

####### Working Directory and Project Specifics ############
workdir: '/pine/o/y/onyen/somefolder'
WRKDIR = '/pine/o/y/onyen/somefolder' # same as above, this is for string concatenation
readWD = '/proj/ideel/YOURDIRECTORY/'
SAMPLES, = glob_wildcards(readWD + 'symlinks/{samp}_R1.fastq.gz')
MERGEDSAMPS, = glob_wildcards(WRKDIR + 'aln/{ms}.merged.bam')
MTDT = '/your/directory/metadata.txt'

################   REFERENCE    ##############
ref = '/your/reffasta'

######## Tools to Call #########
PICARD = '/proj/ideel/apps/linuxbrew/Cellar/picard-tools/2.18.4/bin/picard'
TMPDIR = '/pine/scr/o/y/onyen/somefolder'

###########################
####### Target #############
############################
#Rule all checks to see if file is the same and follows directions up to specified point
rule all:
#	input: expand('aln/{samp}.matefixed.bam', samp = SAMPLES)
#   input: 'aln/CmdlineMerge_cmppnasaln.sh'
#	input: 'merge.log.file'                                                 # Run to here and then check file
#	input: expand('aln/{merge}.bam.bai', merge = MERGEDSAMPS)

###############################################################################

############################
######## Alignment #########
############################


rule index_merged:
	input: 'aln/{merge}.merged.bam'
	output: 'aln/{merge}.merged.bam.bai'
	shell: '{JAVA} -jar {PICARD} BuildBamIndex INPUT={input} OUTPUT={output} TMP_DIR={TMPDIR}'

rule merge_matefixed:
	input: 'aln/CmdlineMerge_cmppnasaln.sh'
	output: 'merge.log.file'
	shell: 'bash {input}; echo "The merge was completed using the custom Rscript BamMergeCommander -- good idea to check the CmdlineMerge table to confirm merge was appropriate" > {output}'

rule make_merge_commandLine:
	input:  metadata='/yourproject/metadata.txt',
	params: dedupbams='/yourproject/aln',
	    outdir = '/yourproject/aln/',
	output: 'aln/CmdlineMerge_cmppnasaln.sh'
	shell: 'BamMergeCommander --input {params.dedupbams} --metadata {input.metadata} --outdir {params.outdir} --output {output}'
# This is a tool by Nick Brazeau. It is in the ideel binary. Must have this in your path to run this
# Merge steps are only needed if you have multiple runs

rule fastq2bam:
	input: read1 = 'symlinks/{samp}_R1.fastq.gz',
		   read2 = 'symlinks/{samp}_R2.fastq.gz',
		   adapters = '/nas/longleaf/apps/trimmomatic/0.36/Trimmomatic-0.36/adapters/TruSeq3-PE.fa',
	params: minlength = '50',
	output: aligned = 'aln/{samp}.matefixed.bam', trimlog = 'logs/{samp}.adaptertrimlog.txt',
	shell:
		r"""
		cutadapt --interleaved \
			-a file:{input.adapters} \
			-A file:{input.adapters} \
			-m {params.minlength} \
			{input.read1} {input.read2} \
			2>{output.trimlog} | \
		bwa mem -t {params.threads} -YK100000000 \
			-H '@PG\tID:cutadapt\tCL:cutadapt -a file:{input.adapters} -A file:{input.adapters} -m {params.minlength} --interleaved {input.read1} {input.read2} 2>{output.trimlog}' \
			-R '{params.rg}' -p \
			{input.ref} - | \
		samblaster --addMateTags | \
		samtools view -bhS - >{output.aligned}
		"""
