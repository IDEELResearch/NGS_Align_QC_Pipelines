###############################################################################
# Purpose:  SnakeMake File to align amplicons
# Authors: NFB                                                 
#Given: FASTQ
#Return: BAM 
#Remarks:                                               
###############################################################################


####### Working Directory and Project Specifics ############
WRKDIR = '/pine/'
readWD = '/proj/ideel/YOURDIRECTORY'
SAMPS, = glob_wildcards(WRKDIR + 'symlinks/{samp}_R1.fastq.gz')

######## REFERENCE SEQ ########
ref = '/ref/directory_for_bowtie'


######## Always on #########
PICARD = '/proj/ideel/apps/brew/share/java/picard.jar'
GATK = '/proj/ideel/apps/brew/share/java/GenomeAnalysisTK.jar'
FLASH = '/proj/ideel/apps/brew/Cellar/flash/1.2.11/bin/flash'
TRIMMOMATIC = '/proj/ideel/apps/brew/share/java/trimmomatic-0.36.jar'
TMPDIR = '/pine/scr/n/f/nfb/PicardandGATKscratch'

##########################################################################################

############################
####### Target #############
############################
#Rule all checks to see if file is the same and follows directions up to specified point
rule all:
#    input: expand('symlinks/{samp}_R1.PAIREDtrimmomatictrimmed.fastq.gz', samp=SAMPS)
#    input: expand('symlinks/{samp}_Stitch.tab.txt', samp=SAMPS)
#	input: expand('aln/{samp}.sorted.bam.bai', samp=SAMPS)
	input: expand('aln/{samp}.realn.bam', samp=SAMPS)
###############################################################################




############################
######## Alignment #########
############################

rule realn_indels:
	input: bam = 'aln/{samp}.sorted.bam', chrs = 'intervals/cytb.intervals', targets = 'aln/{samp}.realigner.intervals', 
	output: 'aln/{samp}.realn.bam'
	shell: 'java -jar {GATK} -T IndelRealigner \
		-R {ref2} -I {input.bam} \
		-L {input.chrs} -targetIntervals {input.targets} \
		-o {output}' 
		# all_chrs.intervals includes just chrs and mito -- it is similar to a bed file for GATK Caller
		# -fixMisencodedQuals must be added for SRA data
		# Interval file from CMP direcotry

rule find_indels:
	input: bam = 'aln/{samp}.sorted.bam', index = 'aln/{samp}.sorted.bam.bai', chrs = 'intervals/cytb.intervals'
	output: 'aln/{samp}.realigner.intervals'
	shell: 'java -jar {GATK} -T RealignerTargetCreator \
		-R {ref2} -I {input.bam} \
		-L {input.chrs} -o {output}'
		# all_chrs.intervals includes just  chrs and mito

rule index_bam: 
	input: 'aln/{samp}.sorted.bam'
	output: 'aln/{samp}.sorted.bam.bai'
	shell: 'java -jar {PICARD} BuildBamIndex INPUT={input} OUTPUT={output} TMP_DIR={TMPDIR}'


rule sort_bam:
	input: 'aln/{samp}.bam'
	output: 'aln/{samp}.sorted.bam'
	shell: 'java -jar {PICARD} SortSam \
		I={input} O={output} \
		SO=coordinate TMP_DIR={TMPDIR}'

# rule matefix:
# 	input: 'aln/{samp}.bam'
# 	output: 'aln/{samp}.matefixed.bam'
# 	shell: 'java -jar {PICARD} FixMateInformation INPUT={input} OUTPUT={output} TMP_DIR={TMPDIR}'


rule stitchedfastq_to_bam:
	input: 'symlinks/{samp}.extendedFrags.fastq.gz'
	output: 'aln/{samp}.bam'
	shell: 'bowtie2 -x {ref} -U {input} -S {output} \
	    --rg-id "bowtie" --rg "PL:illumina\tLB:{wildcards.samp}\tSM:{wildcards.samp[0]}{wildcards.samp[1]}{wildcards.samp[3]}{wildcards.samp[4]}{wildcards.samp[5]}{wildcards.samp[6]}{wildcards.samp[7]}{wildcards.samp[8]}" '

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

	 