# This is a sample-pipeline for WGS resequencing analysis.

# The reference fasta:

    # The reference sequence was generate by ucsc.hg19.fasta
    # 
    # Mitochondrial chromosome only:
    #
    
    head -n 333 /home/hpc/cychen/user/dungchi/TaiwanBiobank/t1000g/ucsc.hg19.fasta > test.fa

# The example reads:
    # https://bioinf.comav.upv.es/courses/sequence_analysis/final_practice.html
    #
    # We have recived 20000 human mithocondrial reads from a yoruba individual. 
    #
    # Link: https://bioinf.comav.upv.es/courses/sequence_analysis/_downloads/mito_yoruba_reads_pe.20k.fastq.gz
    #
    # The reads are Illumina paired end and we want to assemble them:

    wget https://bioinf.comav.upv.es/courses/sequence_analysis/_downloads/mito_yoruba_reads_pe.20k.fastq.gz --no-check-certificate

# Extract the file

    gunzip -c mito_yoruba_reads_pe.20k.fastq.gz  > test.fastq 

# Split one pair-end into two files

    sga-deinterleave.pl test.fastq mito_yoruba.r1.fastq mito_yoruba.r2.fastq

# Known VCF need

    wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.gz
    wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.idx.gz
    wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
    wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx.gz
    wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.gz
    wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.idx.gz
