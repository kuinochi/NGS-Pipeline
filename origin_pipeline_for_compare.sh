#!/bin/bash
#PBS -q large
#PBS -N T1000g.issac
#PBS -l walltime=900:00:00,nodes=1:ppn=30
#PBS -e /home/hpc/cychen/user/dungchi/TaiwanBiobank/ISSAC_BAM/MountDevice/log.err
#PBS -o /home/hpc/cychen/user/dungchi/TaiwanBiobank/ISSAC_BAM/MountDevice/log.out

## 程式執行檔或mpirun指令
## Set Path Variable
set -e
#BIN=/home/u00cyc01/bin
#BWA=$BIN/bwa
#BOWTIE2=$BIN/bowtie2
#BOWTIE2BUILD=$BIN/bowtie2-build
#SAMTOOL=$BIN/samtools
#PICARD=$BIN/picard.jar
#GATK=$BIN/GATK/gatk_3.5/GenomeAnalysisTK.jar
#BCFTOOL=$BIN/bcftools
### ANNOVAR=$BIN/annovar
### SNPEFF=$SRC/snpEff
#JAVA=/pkg/biology/java/jdk_1.8.0_92/bin/java
### Set Reference File
#REF=/work5/NRPB_user/u00cyc01/Reference/ucsc.hg19.fasta
#KNOWN=/work5/NRPB_user/u00cyc01/Knownvcf

## Path Variable
BIN=/home/hpc/cychen/user/dungchi/bin
BWA=$BIN/bwa-0.7.12/bwa
BOWTIE2=$BIN/bowtie2-2.2.6
SAMTOOL=$BIN/samtools-1.2/samtools
PICARD=$BIN/picard-tools-1.139/picard.jar
#ANNOVAR=$BIN/annovar
GATK=$BIN/GATK/GenomeAnalysisTK.jar
BCFTOOL=$BIN/bcftools/bcftools
#SNPEFF=$SRC/snpEff
KNOWN=/home/hpc/cychen/user/dungchi/reference/knownvcf
JAVA=java
REF=/home/hpc/cychen/user/dungchi/TaiwanBiobank/t1000g/ucsc.hg19.fasta

## Working directory
#WD=/work5/NRPB_user/u00cyc01/Workspace  # Change Working directory to YourFolder.
WD=/home/hpc/cychen/user/dungchi/TaiwanBiobank/ISSAC_BAM/MountDevice
cd $WD

echo -e "##############################################################################"
echo -e "###              Work started: $(date)                ###"
echo -e "##############################################################################\n"

## Indexing the reference 
   # echo -e "#-- Indexing the reference fasta file: ${REF}"

   # echo -e "#-- Start indexing by bwa: $(date)"
   # $BWA index -a bwtsw $REF
   # echo -e "#-- BWA indexing done: $(date)"
   # 
   # echo -e "#-- Start indexing by bowtie2: $(date)"
   # $BOWTIE2BUILD $REF $REF
   # echo -e "#-- Bowtie2 indexing done: $(date)"
   # 
   # echo -e "#-- Start indexing by samtools: $(date)"
   # $SAMTOOL faidx $REF
   # echo -e "#-- Samtools indexing done: $(date)"

   # echo -e "#-- Start indexing by picard: $(date)"
   # java -jar $PICARD CreateSequenceDictionary R=${REF} O=ucsc.hg19.dict2
   # echo -e "#-- Picard Indexing done: $(date) \n"
  
   # echo -e "#-- Indexing the reference fasta file ${REF} done"

## Specify Sample file
## ${BSUB_CALL}"NGS-pipeline.sh ${FOLDER} ${SAMPLE} ${FILE}"
#FOLDER=$1
#SAMPLE=$2
#FILENAME=$3
#echo ${FOLDER}
#echo ${SAMPLE}
#echo ${FILENAME}


#SAMPLE=NGS2015036A
FILENAME=NGS2015036A_S1.bam

SAMPLE=NGS2015036A.NAS
#FILENAME="/data5/BAM_iSAAC\ \(152人\)/NGS2015036A/NGS2015036A_S1.bam"

echo -e "Process sample name: ${SAMPLE}"
echo -e "BAMFile: ${FILENAME}"

## Retrive Sample from cloud
# rm -rf ${SAMPLE}
# mkdir -p ${SAMPLE}
# echo -e "#-- Retrive ${SAMPLE} BAM file: ${FILENAME}" 
# echo -e "#-- lftp-cloud folder: ${FOLDER}" 
# echo -e "#-- $(date) "
#lftp -c "set xfer:clobber on; set ftp:passive-mode true; open -u u00cyc01,++c4Lab65334ntuH st.nchc.tw:2122; cd ${FOLDER}; mget -O ${WD}/${SAMPLE} -c ${FILENAME}"
#lftp -c "set xfer:clobber on; set ftp:passive-mode true; open -u u00cyc01,++c4Lab65334ntuH st.nchc.tw:2122; cd ${FOLDER}; pget -n 30 ${FILENAME} -o ${WD}/${SAMPLE}/${FILENAME}"
# echo -e "#-- Retrive done: $(date)"
#cd ${SAMPLE}
#ln -s ${FILENAME} .

## Extract Fastq 
  echo -e "#-- Extract Fastq from BAM file by PICARD SamToFastq: $(date)"
  $JAVA -XX:ParallelGCThreads=30 -Xmx48g -jar $PICARD SamToFastq I=${FILENAME} FASTQ=${SAMPLE}.r1.fastq SECOND_END_FASTQ=${SAMPLE}.r2.fastq \
  VALIDATION_STRINGENCY=LENIENT UNPAIRED_FASTQ=${SAMPLE}.unpaired.fastq INCLUDE_NON_PF_READS=True MAX_RECORDS_IN_RAM=5000000
  echo -e "#-- Extract Fastq done: $(date)"
         #samtools view -f 4 -@ 20 $i > unmapped.bam
        #java -Xmx40g -jar $PICARD CollectAlignmentSummaryMetrics R=$REF I=$i O=$f.collectalignmnetsummary.output.txt
        #java -Xmx40g -jar $PICARD CollectInsertSizeMetrics I=$i O=$f.insert_size_metrics.txt H=$f.insert_size_histogram.pdf M=0.05


## Mapping the reads to the reference genome
 ## Mapping by BWA
 
  echo -e "#-- Mapping reads to REF by BWA-mem" 
  echo -e "#-- Start align sample ${SAMPLE} reads: $(date)"
  $BWA mem -t 30 -a -M $REF ${SAMPLE}.r1.fastq ${SAMPLE}.r2.fastq > ${SAMPLE}.bwa.sam
  echo -e "#-- Sample ${SAMPLE} alignment by BWA is done: $(date)"
  
  ## Process SAM file to BAM file by Picard (Sorting, Mark Duplicates, Edit @RG, Indexing)
  echo -e "#-- Process SAM/BAM file by PICARD start: $(date)"
  
  echo -e "#-- PICARD SortSam for ${SAMPLE}.bwa.sam start: $(date)"
  $JAVA -XX:ParallelGCThreads=30 -Xmx48g -jar $PICARD SortSam I=${SAMPLE}.bwa.sam O=${SAMPLE}.bwa.srt.bam SO=coordinate
  echo -e "#-- PICARD SortSam for ${SAMPLE}.bwa.sam done: $(date)"
  
  echo -e "#-- PICARD MarkDuplicates for ${SAMPLE}.bwa.srt.bam $(date)"
  $JAVA -XX:ParallelGCThreads=30 -Xmx48g -jar $PICARD MarkDuplicates I=${SAMPLE}.bwa.srt.bam O=${SAMPLE}.bwa.srt.dedup.bam M=${SAMPLE}.metrics.txt
  echo -e "#-- PICARD MarkDuplicates for ${SAMPLE}.bwa.srt.bam done: $(date)"
  
  echo -e "#-- PICARD AddOrReplaceReadGroups: for ${SAMPLE}.bwa.srt.dedup.bam start: $(date)" ## Edit @RG groups
  $JAVA -XX:ParallelGCThreads=30 -Xmx48g -jar $PICARD AddOrReplaceReadGroups I=${SAMPLE}.bwa.srt.dedup.bam O=${SAMPLE}.bwa.srt.dedup.rgmd.bam RGID=${SAMPLE} RGSM=${SAMPLE} RGLB=lib_${SAMPLE} RGPL=ILLUMINA RGPU=${SAMPLE}
  echo -e "#-- PICARD AddOrReplaceReadGroups for ${SAMPLE}.bwa.srt.dedup.bam done: $(date)"
  
  echo -e "#-- PICARD BuildBamIndex for ${SAMPLE}.bwa.srt.dedup.rgmd.bam start: $(date)"
  $JAVA -jar $PICARD BuildBamIndex I=${SAMPLE}.bwa.srt.dedup.rgmd.bam
  echo -e "#-- PICARD BuildBamIndex for ${SAMPLE}.bwa.srt.dedup.rgmd.bam done: $(date)"
  
  echo -e "#-- Process SAM/BAM file by PICARD done: $(date)"
 
## DEFINE KNOWNVCF
MILLS_KG_INDEL=$KNOWN/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
DBSNP_138=$KNOWN/dbsnp_138.hg19.vcf
KGPHASE1_INDEL=$KNOWN/1000G_phase1.indels.hg19.sites.vcf

## Calling variant by GATK
echo -e "#-- Calling variants by GATK: $(date)"
echo -e "#-- 1. GATK RealignerTargetCreator start:  $(date)"
$JAVA -jar $GATK -T RealignerTargetCreator -nt 30 -nct 1 -R $REF -I ${SAMPLE}.bwa.srt.dedup.rgmd.bam -known $KGPHASE1_INDEL -known $MILLS_KG_INDEL -o ${SAMPLE}.target_intervals.list
echo -e "#-- 1. GATK RealignerTargetCreator done: $(date)"

echo -e "#-- 2. GATK IndelRealigner start: $(date)"
$JAVA -jar $GATK -T IndelRealigner -R $REF -I ${SAMPLE}.bwa.srt.dedup.rgmd.bam -known $KGPHASE1_INDEL -known $MILLS_KG_INDEL -targetIntervals ${SAMPLE}.target_intervals.list -o ${SAMPLE}.realign.bam --filter_bases_not_stored
echo -e "#-- 2. GATK IndelRealigner is done: $(date)"

echo -e "#-- 3.1 GATK BaseRecalibrator start: $(date)"
$JAVA -jar $GATK -T BaseRecalibrator -R $REF -nct 30 -nt 1 -I ${SAMPLE}.realign.bam -knownSites $DBSNP_138 -knownSites $MILLS_KG_INDEL -knownSites $KGPHASE1_INDEL -o ${SAMPLE}.base.recali.table
echo -e "#-- 3.1 GATK BaseReclibrator done: $(date)"

echo -e "#-- 3.2 GATK PrintReads start: $(date)"
$JAVA -jar $GATK -T PrintReads -R $REF -nct 30 -I ${SAMPLE}.realign.bam -BQSR ${SAMPLE}.base.recali.table -o ${SAMPLE}.realign.recal.bam
echo -e "#-- 3.2 GATK PrintReads done: $(date)"

echo -e "#-- 4. GATK HaplotypeCaller start: $(date)"
$JAVA -jar $GATK -T HaplotypeCaller -R $REF -I ${SAMPLE}.realign.recal.bam -o ${SAMPLE}.realign.recal.HC.vcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -A AlleleBalance -rf BadCigar --dbsnp $DBSNP_138
echo -e "#-- 4. GATK HaplotypeCaller done: $(date)"



echo -e "\n##############################################################################"
echo -e "###                Work done: $(date)                 ###"
echo -e "##############################################################################"

