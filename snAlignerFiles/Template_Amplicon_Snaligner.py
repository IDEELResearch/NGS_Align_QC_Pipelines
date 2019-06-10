###############################################################################
# Purpose:  SnakeMake File to align amplicons
# Authors: Nick Brazeau
#Given: FASTQ
#Return: BAM
#Remarks: Whether or not you stitch, use the realigner step, and trim your adaptors all depends on your project/project specifics.
#Remarks: PLEASE DO ALL ALIGNMENTS IN YOUR SCRATCH SPACE AND THEN TRANSFER OVER THE FILE PRODUCT TO YOUR ANALYSIS DIRECTORY
###############################################################################


####### Working Directory and Project Specifics ############
workdir: '/pine/o/y/onyen/somefolder'
WRKDIR = '/pine/o/y/onyen/somefolder' # same as above, this is for string concatenation
readWD = '/proj/ideel/YOURDIRECTORY/'
SAMPS, = glob_wildcards(WRKDIR + 'symlinks/{samp}_R1.fastq.gz')

######## REFERENCE SEQ ########
ref = '/ref/directory_for_bowtie'


######## Tools to Call #########
GATK = '/nas/longleaf/apps/gatk/3.8-0/GenomeAnalysisTK.jar'
PICARD = '/nas/longleaf/apps/picard/2.10.3/picard-2.10.3/picard.jar'
TRIMMOMATIC = '/proj/ideel/apps/brew/share/java/trimmomatic-0.36.jar'
TMPDIR = '/pine/scr/n/f/nfb/PicardandGATKscratch/'
JAVA = '/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.171-8.b10.el7_5.x86_64/jre/bin/java'

##########################################################################################

############################
####### Target #############
############################
#Rule all checks to see if file is the same and follows directions up to specified point
rule all:
#    input: expand('symlinks/{samp}_R1.PAIREDtrimmomatictrimmed.fastq.gz', samp=SAMPS)
#    input: expand('symlinks/{samp}_Stitch.tab.txt', samp=SAMPS)
#	input: expand('aln/{samp}.bam.bai', samp=SAMPS)

rule index_bam:
	input: 'aln/{samp}.bam.bai'
	output: 'aln/{samp}.bam.bai'
	shell: 'samtools index -b {input} > {output}'
rule stitchedfastq_to_bam:
	input: 'symlinks/{samp}.extendedFrags.fastq.gz'
	output: 'aln/{samp}.bam'
	shell: 'bowtie2 -x {ref} -U {input}  \
	    --rg-id "bowtie" --rg "PL:illumina\\tLB:{wildcards.samp}\\tSM:{wildcards.samp[0]}{wildcards.samp[1]}{wildcards.samp[3]}{wildcards.samp[4]}{wildcards.samp[5]}{wildcards.samp[6]}{wildcards.samp[7]}{wildcards.samp[8]}" \
		| samtools sort - | samtools view -Sb - > {output}'

rule stitch_R1_R2_fastqs:
	input: R1 = 'symlinks/{samp}_R1.PAIREDtrimmomatictrimmed.fastq.gz',
	       R2 ='symlinks/{samp}_R2.PAIREDtrimmomatictrimmed.fastq.gz',
	params: name='{samp}'
	output: 'symlinks/{samp}_Stitch.tab.txt'
	shell: '{FLASH} --max-overlap 300 -z {readWD}{input.R1} {readWD}{input.R2} -o {params.name} -d symlinks | grep -A 5 " Read combination statistics" > {output}'

rule trim_illumina_Adaptors_fastqs:
	 input: 'symlinks/{samp}_R1.fastq.gz', 'symlinks/{samp}_R2.fastq.gz',
	 output: 'symlinks/{samp}_R1.PAIREDtrimmomatictrimmed.fastq.gz', 'symlinks/{samp}_R1.UNPAIREDtrimmomatictrimmed.fastq.gz', 'symlinks/{samp}_R2.PAIREDtrimmomatictrimmed.fastq.gz', 'symlinks/{samp}_R2.UNPAIREDtrimmomatictrimmed.fastq.gz',
	 shell: 'trimmomatic PE -threads 12 -trimlog symlinks/{wildcards.samp}_trim_log.txt {input[0]} {input[1]} {output[0]} {output[1]} {output[2]} {output[3]} ILLUMINACLIP:/nas/longleaf/apps/trimmomatic/0.36/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'
