###############################################################################
# Purpose:  SnakeMake File to take FASTQs and make BAMs
# Authors: Christian Parobek & Nick Brazeau                                        
#Given: FASTQ
#Return: BAM 
#Remarks: Note, you will have to change the paths and project specifics as well as the tools paths (depending on user) for each project/each new directory                               
#         You will need to add the BamMergeCommander to your path
#         Input for BamMergeCommander is a tab delimited file that first column contains a header and sample names or sample IDs (some unique n.alphanumeric for your sample)   
###############################################################################


####### Working Directory and Project Specifics ############
####### Working Directory and Project Specifics ############
workdir: '/pine/'
WRKDIR = '/pine/'
readWD = '/proj/ideel/YOURDIRECTORY'
SAMPLES, = glob_wildcards(readWD + 'symlinks/{samp}_R1.fastq.gz')
MERGEDSAMPS, = glob_wildcards(WRKDIR + 'aln/{ms}.merged.bam')
MTDT = '/your/directory/metadata.txt'


################   REFERENCE    ##############
REF = '/your/reffasta'
GFF = '/your/refgff'

######## Tools to Call #########
######## Always on #########
PICARD = '/proj/ideel/apps/brew/share/java/picard.jar'
GATK = '/proj/ideel/apps/brew/share/java/GenomeAnalysisTK.jar'
FLASH = '/proj/ideel/apps/brew/Cellar/flash/1.2.11/bin/flash'
TRIMMOMATIC = '/proj/ideel/apps/brew/share/java/trimmomatic-0.36.jar'
TMPDIR = '/pine/scr/n/f/nfb/PicardandGATKscratch'

##########################################################################################

###########################
####### Target #############
############################
#Rule all checks to see if file is the same and follows directions up to specified point
rule all:
#   input: expand('symlinks/pairedfastqs/{samp}_R1.PAIREDtrimmomatictrimmed.fastq.gz',  samp = SAMPLES) 
#	input: expand('aln/{samp}.matefixed.bam', samp = SAMPLES) 
   input: 'aln/CmdlineMerge_cmppnasaln.sh'
#	input: 'merge.log.file'                                                 # Run to here and then check file
#	input: expand('aln/{merge}.realn.bam', merge = MERGEDSAMPS) 

###############################################################################




############################
######## Alignment #########
############################

rule realn_indels:
	input: bam = 'aln/{merge}.merged.bam',  targets = 'aln/{merge}.realigner.intervals', 
	output: 'aln/{merge}.realn.bam'
	shell: 'java -jar {GATK} -T IndelRealigner \
		-R {REF} -I {input.bam} \
		-targetIntervals {input.targets} \
		-o {output}' 

rule find_indels:
	input: bam = 'aln/{merge}.merged.bam', index = 'aln/{merge}.merged.bam.bai',
	output: 'aln/{merge}.realigner.intervals'
	shell: 'java -jar {GATK} -T RealignerTargetCreator \
		-R {REF} -I {input.bam} \
		-o {output}'

rule index_merged: 
	input: 'aln/{merge}.merged.bam'
	output: 'aln/{merge}.merged.bam.bai'
	shell: 'java -jar {PICARD} BuildBamIndex INPUT={input} OUTPUT={output} TMP_DIR={TMPDIR}'


rule merge_matefixed:
	input: 'aln/CmdlineMerge_cmppnasaln.sh'
	output: 'merge.log.file'
	shell: 'bash {input}; echo "The merge was completed using the custom Rscript BamMergeCommander -- good idea to check the CmdlineMerge table to confirm merge was appropriate" > {output}'

rule make_merge_commandLine:
	input:  metadata='/proj/ideel/meshnick/users/NickB/Projects/DRC_Vivax/DRC_Pv_Sequences/PvSal1_Aln/DRC_Vivax_mtdt.txt',
	params: dedupbams='/proj/ideel/meshnick/users/NickB/Projects/DRC_Vivax/DRC_Pv_Sequences/PvSal1_Aln/aln', 
	    outdir = '/proj/ideel/meshnick/users/NickB/Projects/DRC_Vivax/DRC_Pv_Sequences/PvSal1_Aln/aln/',
	output: 'aln/CmdlineMerge_cmppnasaln.sh'
	shell: 'BamMergeCommander --input {params.dedupbams} --metadata {input.metadata} --outdir {params.outdir} --output {output}'
# note outdir needs forward slash, dedupbams does not

rule matefix:
	input: 'aln/{samp}.dedup.bam'
	output: 'aln/{samp}.matefixed.bam'
	shell: 'java -jar {PICARD} FixMateInformation INPUT={input} OUTPUT={output} TMP_DIR={TMPDIR}'
	
	
rule mark_dups:
	input: 'aln/{samp}.sorted.bam'
	output:'aln/{samp}.dedup.bam','aln/{samp}.dedup.metrics'
	shell: 'java -jar {PICARD} MarkDuplicates \
		I={input} O={output[0]} \
		METRICS_FILE={output[1]} \
		TMP_DIR={TMPDIR} REMOVE_DUPLICATES=FALSE \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000'
    # Note have changed this to false as part of keeping more information and GATK best practices (https://gatkforums.broadinstitute.org/gatk/discussion/6747) and picard defaults
    
rule sort_bam:
	input: 'aln/{samp}.raw.bam'
	output: 'aln/{samp}.sorted.bam'
	shell: 'java -jar {PICARD} SortSam \
		I={input} O={output} \
		SO=coordinate TMP_DIR={TMPDIR}'



rule fastq_to_bam:
	input: 'symlinks/{samp}_R1.PAIREDtrimmomatictrimmed.fastq.gz'
	output: 'aln/{samp}.raw.bam'
	shell: 'bwa mem {REF} {input[0]} \
		-R "@RG\tID:bwa\tPL:illumina\tLB:{wildcards.samp}_lib\tSM:{wildcards.samp[0]}{wildcards.samp[1]}{wildcards.samp[2]}{wildcards.samp[3]}{wildcards.samp[4]}" \
		 | samtools view -Sb - > {output}'
		# calling the @RG ID: 'bwa' because this resolves a clash with @PG ID --> I updated this recently to make it more unique for MERGING
		# Can controls how long the read sample name is by wild cards for tSM, important if want to merge file later and need different library names but same sample names'


rule trim_illumina_Adaptors_fastqs:
	 input: readWD + 'symlinks/pairedfastqs/{samp}_1.fastq.gz', readWD + 'symlinks/pairedfastqs/{samp}_2.fastq.gz', 
	 output: 'symlinks/pairedfastqs/{samp}_R1.PAIREDtrimmomatictrimmed.fastq.gz', 'symlinks/pairedfastqs/{samp}_R1.UNPAIREDtrimmomatictrimmed.fastq.gz', 'symlinks/pairedfastqs/{samp}_R2.PAIREDtrimmomatictrimmed.fastq.gz', 'symlinks/pairedfastqs/{samp}_R2.UNPAIREDtrimmomatictrimmed.fastq.gz',  
	 shell: 'trimmomatic PE -threads 12 -trimlog symlinks/trim_log.txt {input[0]} {input[1]} {output[0]} {output[1]} {output[2]} {output[3]} ILLUMINACLIP:/nas/longleaf/apps/trimmomatic/0.36/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'
     # The TRUE at the end keeps the paired end reads in R2
     # Want to align the PAIRED trimmed
     #adapterremoval, https://github.com/MikkelSchubert/adapterremoval
