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

READ_DIR	= /home/hpc/cychen/user/dungchi/tools/NGS-Pipeline/NGS-Pipeline/Read/

#READ_SE 	=
#READ_1 		= $(wildcard ${READ_DIR}/*.r1.fq)
#READ_2 		= $(wildcard ${READ_DIR}/*.r2.fq)

# Paired-end or single-end analysis? Must be either PE or SE
#READ_TYPE 	= 

# Path to input .sam file
IN_SAM 		=

# Path to input .bam file
IN_BAM 		=

# Path to input .vcf file
IN_VCF 		=

# ------------------------------------------------------ #
# --- Paths to knownvcf files
# ------------------------------------------------------ #
MILLS_KG_INDEL=/home/hpc/cychen/user/dungchi/reference/knownvcf/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
DBSNP_138=/home/hpc/cychen/user/dungchi/reference/knownvcf/hg19/dbsnp_138.hg19.vcf.gz
KGPHASE1_INDEL=/home/hpc/cychen/user/dungchi/reference/knownvcf/hg19/1000G_phase1.indels.hg19.sites.vcf.gz

# ------------------------------------------------------ #
# --- Paths to external programs     
# ------------------------------------------------------ #

SAMTOOL 		= /home/hpc/cychen/user/dungchi/bin/samtools-1.2/samtools
BWA 			= /home/hpc/cychen/user/dungchi/bin/bwa-0.7.12/bwa
BOWTIE2 		= /home/hpc/cychen/user/dungchi/bin/bowtie2-2.2.6/bowtie2 
BOWTIE2_BUILD	= /home/hpc/cychen/user/dungchi/bin/bowtie2-2.2.6/bowtie2-build
PICARD 			=  /home/hpc/cychen/user/dungchi/bin/picard-tools-1.139/picard.jar
GATK 			=  /home/hpc/cychen/user/dungchi/bin/GATK/GenomeAnalysisTK.jar

#ANNOVAR 		= $BIN/annovar
BCFTOOL 		= /home/hpc/cychen/user/dungchi/bin/bcftool
#BEDTOOLS		=
#VCFTOOLS		=
#TABIX			=

# ------------------------------------------------------ #
# --- Parameters for external programs
# ------------------------------------------------------ #

# Number of threads
CPU = 1

# BWA parameters
BWA_PARAM 		= mem -t ${CPU} -a -M

# Bowtie2 parameters
BOWTIE2_PARAM 	= --end-to-end --sensitive -p ${CPU} -x ${REF_FA} -q

# Number of node or dataset used by GATK
GATK_NMT = 1

# SAMtools mpileup parameters
SNP_MIN_COV	=	5
SNP_MAX_COV	=	100

# BAMtools filter parameters
#MAPQUAL=20

# Should we mark duplicates? TRUE or FALSE
#MARK_DUPS=TRUE

# Max number of file handles to keep open when Picard's MarkDuplicates writes to disk.
# This should be a bit lower than the per-process max number of files that can be open.
# You can find that max using command 'ulimit -n'
# This avoids the "java.io.FileNotFoundException: (Too many open files)" exception
#PICARD_MARK_DUP_MAX_FILES=4000

# -------------------------------------------------------------------------------------- #
# --- Parameters for multi-sample SNP calling
# -------------------------------------------------------------------------------------- #
