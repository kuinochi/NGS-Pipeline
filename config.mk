# ------------------------------------------------------ #
# --- Configuration makefile of user-editable variables 
# --- AWARE: All paths must be absolute or relative to 
# --- the directory where Makefile exists
# ------------------------------------------------------ #

# ------------------------------------------------------ #
# --- Paths to input files
# ------------------------------------------------------ #

# Paths to reference genomes files 
# e.g. hg19.fasta, hg38.fasta
# Must be in FASTA format

REF_DIR		= ./Ref
REF_FA 		= ${REF_DIR}/ucsc.hg19.mt.fasta

# Common name of genome (used to name files)
# e.g. hg19, hg38, mm9 or mm10 
REF_NAME 	= hg19

# Input .fastq reads files
# -- READ_DIR:	Define the path to .fastq files, 
# -- 			The file must be in FASTQ format and end in '.fastq'
# -- 			or in gzip'd FASTQ format and end in '.fastq.gz'
# --
# -- FastQC will not name the output file properly if ending is '.fq'

READ_DIR	= ./Read/

# ------------------------------------------------------ #
# --- Paths to knownvcf files
# ------------------------------------------------------ #

MILLS_KG_INDEL=/home/hpc/cychen/user/dungchi/reference/knownvcf/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
DBSNP_138=/home/hpc/cychen/user/dungchi/reference/knownvcf/dbsnp_138.hg19.vcf
KGPHASE1_INDEL=/home/hpc/cychen/user/dungchi/reference/knownvcf/1000G_phase1.indels.hg19.sites.vcf

# ------------------------------------------------------ #
# --- Paths to external programs     
# ------------------------------------------------------ #

BWA 			= /Users/dungchi/Work/bin/bwa-0.7.13/bwa
BOWTIE2 		= bowtie2 
BOWTIE2_BUILD	= bowtie2-build
PICARD 			= picard.jar
GATK 			= GenomeAnalysisTK.jar
SAMTOOL 		= samtools
BCFTOOL 		= bcftool
#BEDTOOLS		=
#VCFTOOLS		=
#TABIX			=

# ------------------------------------------------------ #
# --- Parameters for external programs
# ------------------------------------------------------ #

# Number of threads
CPU = 1

# BWA parameters
BWA_PARAM 		= mem -t ${CPU} -a -M -K 10000000

# Bowtie2 parameters
BOWTIE2_PARAM 	= --end-to-end --sensitive -p ${CPU} -x ${REF_FA} -q

# Number of node or dataset used by GATK
GATK_NMT = 1

# SAMtools mpileup parameters
SNP_MIN_COV	=	5
SNP_MAX_COV	=	100

