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
# -- READ_1: 	Read 1 file name
# -- READ_2: 	Read 2 file name
# -- sNAME:		Sample Name
# -- sRGID: 	Read group identifier,  each read group's ID must be unique
# -- sRGSM:		The Name of the Sample in RG
# -- sRGPL:		Sample platform, for GATK, Valid values: ILLUMINA, SOLID, LS454, HELICOS and PACBIO.


READ_DIR	= ./Read
READ_1 		= ${READ_DIR}/mito_yoruba.r1.fastq
READ_2		= ${READ_DIR}/mito_yoruba.r2.fastq
sNAME		= mito_yoruba
sRGID		= GP_${sNAME}
sRGSM		= SM_${sNAME}
sRGPL		= ILLUMINA
_Output_Folder = ./results

# ------------------------------------------------------ #
# --- Paths to knownvcf files
# ------------------------------------------------------ #

MILLS_KG_INDEL	= Ref/chrM.mills.vcf
DBSNP_138		= Ref/chrM.dbsnp138.vcf
KGPHASE1_INDEL	= Ref/chrM.1000g.indel.vcf


# ------------------------------------------------------ #
# --- Paths to external programs     
# ------------------------------------------------------ #

BWA 			= /Users/dungchi/Work/bin/bwa-0.7.13/bwa
BOWTIE2 		= bowtie2 
BOWTIE2_BUILD	= bowtie2-build
PICARD 			= picard.jar
GATK 			= GenomeAnalysisTK.jar
SAMTOOL 		= samtools
BCFTOOL 		= bcftools
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

