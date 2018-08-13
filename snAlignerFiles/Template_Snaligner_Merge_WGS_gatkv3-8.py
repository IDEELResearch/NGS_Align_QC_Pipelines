###############################################################################
# Purpose:  SnakeMake File to take FASTQs and make BAMs
# Authors: Nick Brazeau (inspired by Christian Parobek)
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
REF = '/your/reffasta'

######## Tools to Call #########
PICARD = '/proj/ideel/apps/linuxbrew/Cellar/picard-tools/2.18.4/bin/picard'
GATK = '/nas/longleaf/apps/gatk/3.8-0/GenomeAnalysisTK.jar'
TRIMMOMATIC = '/proj/ideel/apps/brew/share/java/trimmomatic-0.36.jar'
JAVA = '/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.171-8.b10.el7_5.x86_64/jre/bin/java'
TMPDIR = '/pine/scr/o/y/onyen/somefolder'

###########################
####### Target #############
############################
#Rule all checks to see if file is the same and follows directions up to specified point
rule all:
#	input: expand('symlinks/pairedfastqs/{samp}_R1.PAIREDtrimmomatictrimmed.fastq.gz',  samp = SAMPLES)
#	input: expand('aln/{samp}.matefixed.bam', samp = SAMPLES)
#       input: 'aln/CmdlineMerge_cmppnasaln.sh'
#	input: 'merge.log.file'                                                 # Run to here and then check file
#	input: expand('aln/{merge}.realn.bam', merge = MERGEDSAMPS)

###############################################################################

############################
######## Alignment #########
############################

rule realn_indels:
	input: bam = 'aln/{merge}.merged.bam',  targets = 'aln/{merge}.realigner.intervals',
	output: 'aln/{merge}.realn.bam'
	shell: '{JAVA} -jar {GATK} -T IndelRealigner \
		-R {REF} -I {input.bam} \
		-targetIntervals {input.targets} \
		-o {output}'

rule find_indels:
	input: bam = 'aln/{merge}.merged.bam', index = 'aln/{merge}.merged.bam.bai',
	output: 'aln/{merge}.realigner.intervals'
	shell: '{JAVA} -jar {GATK} -T RealignerTargetCreator \
		-R {REF} -I {input.bam} \
		-o {output}'

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

rule matefix:
	input: 'aln/{samp}.dedup.bam'
	output: 'aln/{samp}.matefixed.bam'
	shell: '{JAVA} -jar {PICARD} FixMateInformation INPUT={input} OUTPUT={output} TMP_DIR={TMPDIR}'


rule mark_dups:
	input: 'aln/{samp}.sorted.bam'
	output:'aln/{samp}.dedup.bam','aln/{samp}.dedup.metrics'
	shell: '{JAVA} -jar {PICARD} MarkDuplicates \
		I={input} O={output[0]} \
		METRICS_FILE={output[1]} \
		TMP_DIR={TMPDIR} REMOVE_DUPLICATES=FALSE \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000'
    # Note have changed this to false as part of keeping more information and GATK best practices (https://gatkforums.broadinstitute.org/gatk/discussion/6747) and picard defaults

rule sort_bam:
	input: 'aln/{samp}.raw.bam'
	output: 'aln/{samp}.sorted.bam'
	shell: '{JAVA} -jar {PICARD} SortSam \
		I={input} O={output} \
		SO=coordinate TMP_DIR={TMPDIR}'



rule fastq_to_bam:
	input: R1='/yourpath/symlinks/{samp}_R1.fastq.gz', R2='/yourpath/symlinks/{samp}_R2.fastq.gz'
#	input: R1='/yourpath/symlinks/{samp}_R1.PAIREDtrimmomatictrimmed.fastq.gz', R2='/yourpath/symlinks/{samp}_R2.PAIREDtrimmomatictrimmed.fastq.gz'
	output: 'aln/{samp}.raw.bam'
	shell: 'bwa mem {REF} {input.R1} {input.R2} \
		-R "@RG\\tID:bwa\\tPL:illumina\\tLB:{wildcards.samp}_lib\\tSM:{wildcards.samp[0]}{wildcards.samp[1]}{wildcards.samp[2]}{wildcards.samp[3]}{wildcards.samp[4]}" \
		 | samtools view -Sb - > {output}'
		# calling the @RG ID: 'bwa' because this resolves a clash with @PG ID
		# You need to control how long the read sample name (i..e SM) is by the wild cards for merging bams. Merged bams need different library names but same sample names


rule trim_illumina_Adaptors_fastqs:
	 input: readWD + 'symlinks/{samp}_1.fastq.gz', readWD + 'symlinks/{samp}_2.fastq.gz',
	 output: 'symlinks/{samp}_R1.PAIREDtrimmomatictrimmed.fastq.gz', 'symlinks/{samp}_R1.UNPAIREDtrimmomatictrimmed.fastq.gz', 'symlinks/{samp}_R2.PAIREDtrimmomatictrimmed.fastq.gz', 'symlinks/{samp}_R2.UNPAIREDtrimmomatictrimmed.fastq.gz',
	 shell: 'trimmomatic PE -threads 12 -trimlog symlinks/trim_log.txt {input[0]} {input[1]} {output[0]} {output[1]} {output[2]} {output[3]} ILLUMINACLIP:/nas/longleaf/apps/trimmomatic/0.36/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'
     # The TRUE at the end keeps the paired end reads in R2
     # Want to align the PAIRED trimmed
	 # You should only trim if you have read through...Ideally you don't, seems to happen often in reality
