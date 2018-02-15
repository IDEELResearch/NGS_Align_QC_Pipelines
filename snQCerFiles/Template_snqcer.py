
###############################################################################
# Purpose:  SnakeMake File TEMPLATE to check for Coverage and QC on the Bam files
# Authors: Nick Brazeau  & Dayne Filer   & Christian Parobek                                               
#Given: Realn.Bam
#Return: coverage, qc       
#Remarks: As before, you will need to change your paths 
#      	  There is a clash issue with snakemake using python 3 and multiqc python 2. Need to run from command line      
#		  Note, the heatmap coverage plot and the bed map plot are not very generalizable and need to be updated. Don't use...depreciated and not up to date at this time   
###############################################################################


####### Working Directory and Project Specifics ############
WRKDIR = '/proj/ideel/YOURDIRECTORY'
readWD = '/proj/ideel/YOURDIRECTORY'
FQS, = glob_wildcards(WRKDIR + 'symlinks/{fq}.fastq.gz')
SAMPLES, = glob_wildcards(WRKDIR + 'aln/{sample}.realn.bam')

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

############################
####### Target #############
############################
#Rule all checks to see if file is the same and follows directions up to specified point
rule all:
# SUMMARYSTATS
#	input: 'symlinks/FASTQC/Log.txt'
	input: expand('SumSTATsandQC/ValidateBAM/{sample}.validatebams.txt', sample = SAMPLES)
#	input: expand('SumSTATsandQC/AlignSummary/{sample}.AlignSummary.Metrics', sample = SAMPLES)
#	input: expand('SumSTATsandQC/FlagStats/{sample}.samtools.flagstats', sample = SAMPLES)
#	input: expand('SumSTATsandQC/coverage/data/{sample}.cov', sample = SAMPLES)
#	input: 'SumSTATsandQC/coverage/{params.prefix}cov_heatmap.pdf'
#   input: expand('SumSTATsandQC/CallableLoci/{sample}_callable_status.bed', sample = SAMPLES)
###############################################################################



######################################
########   Quality Control   #########
#####################################

rule CallableLoci_By_SAMPLE:
	input: bams = 'aln/{sample}.realn.bam',
	output: bedfile = 'SumSTATsandQC/CallableLoci/{sample}_callable_status.bed', summarytable = 'SumSTATsandQC/CallableLoci/{sample}_summarytable.txt'
	shell: 'java -jar {GATK} -T CallableLoci \
		-R {REF} \
		-I {input.bams} \
		--minBaseQuality 20 \
		--minMappingQuality 10 \
		--minDepth 4 \
		-summary {output.summarytable} \
		-o {output.bedfile}'

rule plot_coverage:
	input:	expand('SumSTATsandQC/coverage/data/{sample}.cov', sample = SAMPLES)
	params:	idir = 'SumSTATsandQC/coverage/data/', 
		prefix = '1877_Plasmepsin_', 
		odir = 'SumSTATsandQC/coverage/'
	output:	'SumSTATsandQC/coverage/{params.prefix}cov_plot.pdf', 
		'SumSTATsandQC/coverage/{params.prefix}cov_heatmap.pdf'
	shell:	'covPlotter -I {params.idir} -O {params.odir} -F {params.prefix}'


rule calculate_cov_forheatmap:
	input:	'aln/{sample}.recal.realn.bam'
	output:	'SumSTATsandQC/coverage/data/{sample}.cov'
	shell:	'bedtools genomecov \
		-ibam {input} -g {GFF} -max 500 | grep "genome\|PFC10_API\|M76\|Pf3" \
		> {output}'
# Note for Pv switch out to: -ibam {input} -g {GFF} -max 500 | grep "genome\|chr\|NC_" \
	
rule summary_stats:
	input: 'aln/{sample}.recal.realn.bam'
	output: 'SumSTATsandQC/FlagStats/{sample}.samtools.flagstats'
	shell: 'samtools flagstat {input} > {output}'
# Provides alignment information
		
rule AlignSummaryMetrics:
	input: 'aln/{sample}.recal.realn.bam'
	output: 'SumSTATsandQC/AlignSummary/{sample}.AlignSummary.Metrics'
	shell: 'java -jar {PICARD} CollectAlignmentSummaryMetrics \
          R={REF} \
          I={input} \
          O={output}'
# Need this for read length averages, etc. 

###############################################
########   VALIDATE SAM FILE MODULE   #########
###############################################
rule ValidateSamFile:
	input:'aln/{sample}.realn.bam'
	output: 'SumSTATsandQC/ValidateBAM/{sample}.validatebams.txt'
	shell: 'java -jar {PICARD} ValidateSamFile \
	    MODE=SUMMARY \
		I={input} > {output}'
	# After this step you should run from the command line `cat 'SumSTATsandQC/ValidateBams.tab.txt | grep "ERROR"` -- if there are errors, STOP and figure out why	
		
######################################
########   FASTQC MODULE   #########
#####################################
rule fastqc:
	input: expand('symlinks/{fq}.fastq.qz', fq=FQS)
	output: outdir='symlinks/FASTQC/', log='symlinks/FASTQC/Log.txt',
	shell: 'fastqc {input} --outdir {output.outdir}; echo "Fastqc succesfully completed" > {output.log} '

	
