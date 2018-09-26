
###############################################################################
# Purpose:  SnakeMake File TEMPLATE to check for Coverage and QC on the Bam files
# Authors: Nick Brazeau  & Dayne Filer   & Christian Parobek
#Given: Realn.Bam
#Return: coverage, qc
#Remarks: As before, you will need to change your paths
#Remarks: PLEASE RUN ALL OF THESE STEPS IN YOUR SCRATCH SPACE (except the final report)
################################################################################


####### Working Directory and Project Specifics ############
workdir: '/proj/ideel/YOURDIRECTORY'
WRKDIR = '/proj/ideel/YOURDIRECTORY'
readWD = '/proj/ideel/YOURDIRECTORY'
FQS, = glob_wildcards(WRKDIR + 'symlinks/{fq}.fastq.gz')
SAMPLES, = glob_wildcards(WRKDIR + 'aln/{sample}.realn.bam')

################   REFERENCE    ##############
REF = '/your/reffasta'
GFF = '/your/refgff'

######## Tools to Call #########
######## Always on #########
PICARD = '/proj/ideel/apps/linuxbrew/Cellar/picard-tools/2.18.4/bin/picard'
GATK = '/nas/longleaf/apps/gatk/3.8-0/GenomeAnalysisTK.jar'
TRIMMOMATIC = '/proj/ideel/apps/brew/share/java/trimmomatic-0.36.jar'
JAVA = '/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.171-8.b10.el7_5.x86_64/jre/bin/java'
TMPDIR = '/pine/scr/o/n/onyen/'


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
# 	input: expand('SumSTATsandQC/coverage/data/{sample}.long.cov', sample=SAMPLES)
#	input: 'report.txt'
###############################################################################


######################################
######   Pull it all together   ######
#####################################
rule render_rmarkdown:
	output: 'report.txt'
	shell:
		"""r
		Rscript -e 'rmarkdown::render("WGS_coverage.Rmd", clean=TRUE)' ; \
			echo "Report finished" > {output}
		"""
# IMPORTANT -- this step requires you to downloada  WGS_Coverage_template.Rmd file and put your appropriate paths in it
# You should also add in your project details!! (highly recommended)
# It can live wherever you want/wherever your project is... i.e. maybe the report lives in your analysis directory and not your scratch space


######################################
########   Coverage Calc   #########
#####################################
rule calculate_cov:
	input:	'aln/{sample}.realn.bam'
	output:	'SumSTATsandQC/coverage/data/{sample}.long.cov'
	shell:	'bedtools genomecov \
		-ibam {input} -d  \
		> {output}'

# note you can use | grep "chromosome names" if you have a lot of contigs that you don't care about

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

rule summary_stats:
	input: 'aln/{sample}.realn.bam'
	output: 'SumSTATsandQC/FlagStats/{sample}.samtools.flagstats'
	shell: 'samtools flagstat {input} > {output}'
# Provides alignment information

rule AlignSummaryMetrics:
	input: 'aln/{sample}.realn.bam'
	output: 'SumSTATsandQC/AlignSummary/{sample}.AlignSummary.Metrics'
	shell: '{JAVA} -jar {PICARD} CollectAlignmentSummaryMetrics \
          R={REF} \
          I={input} \
          O={output}'
# Need this for read length averages, etc.

rule ValidateSamFile:
	input:'aln/{sample}.realn.bam'
	output: 'SumSTATsandQC/ValidateBAM/{sample}.validatebams.txt'
	shell: '{JAVA} -jar {PICARD} ValidateSamFile \
	    MODE=SUMMARY \
		I={input} > {output}'
	# After this step you should run from the command line `cat 'SumSTATsandQC/ValidateBams.tab.txt | grep "ERROR"` -- if there are errors, STOP and figure out why

rule fastqc:
	input: expand('symlinks/{fq}.fastq.qz', fq=FQS)
	output: outdir='symlinks/FASTQC/', log='symlinks/FASTQC/Log.txt',
	shell: 'fastqc {input} --outdir {output.outdir}; echo "Fastqc succesfully completed" > {output.log} '
