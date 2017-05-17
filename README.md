# About
This repository contains a simplified pipeline for WGS resequencing analysis, written in Make.
The primary goal of this pipeline is to generate two VCF from raw FASTQ files.

# Usage 
Type `make` to see the usage.

# Steps: 
First, you need to copy the config.mk.default to config.mk, then edit the necessary variable in it.

## 1. Index the reference for following steps
To indexing the reference FASTA file for downstream analysis.
Type `make IndexRef`

## 2. Alingment

To align the reads to the reference FASTA file, generating a SAM files.
Type `make Alignment`

## 3. Preprocess

Processing the SAM files to an pre-analysis-ready BAM file, including sorting, markduplicates and indexing.
Type `make Preprocess`

## 4. GATK_preprocess

Further processing the pre-analysis-ready BAM file to a analysis-ready BAM file, including RTC, IR, BQSR, PR by GATK.
Type `make GATK_preprocess`

## 5. Calling variants

Variants discovery either by samtools mpileup or GATK HaplotypeCaller.
Type `make callvariants_sam` and `make callvariants_gatk` 


## Overall
You can also try to run all the steps by typing `make all`


---------

### Appendix:
#### Testing data:
If you want, you can generate small example data for test by the following steps:

##### 1. The reference fasta:
    
    # You can generated one small reference sequence from  ucsc.hg19.fasta
    # 
    # Ex, Mitochondrial chromosome only:
    #
    `head -n 333 /Path/to/ucsc.hg19.fasta > ucsc.hg19.mt.fa`

##### 2. The example reads:
    
    # You can download the small reads set from:
    # 
    # https://bioinf.comav.upv.es/courses/sequence_analysis/_downloads/mito_yoruba_reads_pe.20k.fastq.gz
    #
    # The reads are Illumina paired end and we want to assemble them:

    `wget https://bioinf.comav.upv.es/courses/sequence_analysis/_downloads/mito_yoruba_reads_pe.20k.fastq.gz --no-check-certificate`

    ### Extract the file

    `gunzip -c mito_yoruba_reads_pe.20k.fastq.gz  > test.fastq` 

    ### Split one pair-end into two files

    `sga-deinterleave.pl test.fastq mito_yoruba.r1.fastq mito_yoruba.r2.fastq`

##### Known VCF needed

    # You will also need some known VCF in this pipeline
    # You can download from GATK bundle.
    `wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.gz`
    `wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.idx.gz`
    `wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz`
    `wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx.gz`
    `wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.gz`
    `wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.idx.gz`
