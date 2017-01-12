#---------------------------------------------------------------------------------#
# Makefile for running next generation sequencing workflow, 
# including mapping, variant calling
# Please check the needed config.mk before running the pipeline
#---------------------------------------------------------------------------------#

#SHELL = !/bin/bash

# Get user editable variables
include config.mk

# 	.marked.realigned.fixed.recal.indexed

# Default goal, print usage
help: 
	@printf "\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n\n" \
	"Please enter the steps you want to proceed" \
	"-------------------------------------------" \
	"IndexRef       --  Indexing a FASTA file by four tools:" \
	"                   bwa index, samtools faidx, bowtie2-build" \
	"                   and CreateSequenceDictionary of picard." \
	"bwa            --  Align the reads to the reference genome by bwa-mem" \
	"bowtie2        --  Align the reads to the reference genome by bowtie2" \
	"Process        --  Process SAM file to BAM file by Picard:" \
	"                   Sorting, Mark Duplicates, Edit @RG, Indexing" \
	"variants_sam   --  Variants Calling (SNP/INDELs) by samtools " \
	"vairants_gatk  --  Variants Calling (SNP/INDELs) by samtools"


bwa: 			results/alignment/bwa/*.sam
bowtie2: 		results/alignment/bowtie2/*.sam
Process: 		index_ref ${BAM_FILE}
variants_sam:		$
vairants_gatk:

#=================================================================================#
#--------------------       Indexing reference FASTA       -----------------------#
#=================================================================================#

# Output files of indexing.
_SAMTOOLS_INDEX		= $(addsuffix .fai, ${REF_FA})
_BWA_INDEX_END	 	= .amb .ann .bwt .pac .sa
_BWA_INDEX			= $(addsuffix ${_BWA_INDEX_END}, ${REF_FA})
_BOWTIE2_INDEX_END	= .1.bt2 .2.bt2 .3.bt2 .4.bt2 .rev.1.bt2 .rev.2.bt2
_BOWTIE2_INDEX		= $(addsuffix ${_BOWTIE2_INDEX_END}, ${REF_FA})
_PICARD_DICT		= $(addsuffix .dict, ${REF_FA})

IndexRef: ${_SAMTOOLS_INDEX} ${_BWA_INDEX} ${_BOWTIE2_INDEX} ${_PICARD_DICT}
	@printf "\n%s\n%s\n" \
	"The file:  ${REF_FA}  has been indexed." \
	"Check the folder: ${REF_DIR}"

${_SAMTOOLS_INDEX}: ${REF_FA}
	@printf "\n%s\n" \
	"#--- Start indexing by samtools..."
	${SAMTOOL} faidx $<
	@printf "%s\n" \
	"#--- Samtools Indexing is done!"
	@touch ${_SAMTOOLS_INDEX}

${_BWA_INDEX}: ${REF_FA}
	@printf "\n%s\n" \
	"#--- Start indexing by bwa..."
	${BWA} index -a bwtsw $<
	@printf "%s\n" \
	"#--- BWA Indexing is done! \n"
	@touch ${_BWA_INDEX}

${_BOWTIE2_INDEX}: ${REF_FA}
	@printf "\n%s\n" \
	"#--- Start indexing by bowtie2..."
	${BOWTIE2_BUILD} $< $< 
	@printf "%s\n" \
	"#--- Bowtie2 Indexing is done! \n"
	@touch ${_BOWTIE2_INDEX}	

${_PICARD_DICT}: ${REF_FA}
	@printf "\n%s\n" \
	"#--- Start indexing by picard tools..."
	${PICARD} CreateSequenceDictionary R=$< O=$@
	@printf "%s\n" \
	"#--- Picard Indexing is done! \n"
	@touch ${_PICARD_DICT}


#=================================================================================#
#----------------------       Mapping to reference       -------------------------#
#=================================================================================#

#--- Align reads to the ref. genome by BWA

results/alignment/bwa/*.sam: ${READ_DIR}
	@mkdir -p results/alignment/bwa/ 

	for i in ${READ_DIR}/*.r1.fq; do F=`basename $$i .r1.fq` ;\
	echo "Start align sample \"$$F\" by bwa-mem, `date`" ;\
	${BWA} ${BWA_PARAM} \
		${REF_FA} ${READ_DIR}/$$F.r1.fq ${READ_DIR}/$$F.r2.fq > results/alignment/bwa/$$F.${REF_NAME}.bwa.sam ;\
	echo "Sample \"$$F\" alignment by bwa-mem is done! `date`\n"; done\

	@echo "Mapping by bwa-mem of all samples in the folder \"${READ_DIR}\" is done!"


#--- Align reads to the ref. genome by bowite2
results/alignment/bowtie2/*.sam: ${READ_DIR} 
	@mkdir -p results/alignment/bowtie2/;
	for i in ${READ_DIR}/*.r1.fq; do F=`basename $$i .r1.fq` ;\
	echo "Start align sample \"$$F\" by bowtie2, `date`" ;\
	${BOWTIE2} ${BOWTIE2_PARAM} \
	-1 ${READ_DIR}/$$F.r1.fq -2 ${READ_DIR}/$$F.r2.fq -S results/alignment/bowtie2/$$F.${REF_NAME}.bowtie2.sam;\
	echo "Sample \"$$F\"" alignment by bowtie2 is done! `date`\n"; done\

	@echo "Mapping by bowtie2 of all samples in the folder \"${READ_DIR}\" is done!"


#=================================================================================#
#----------------------      Alignment file process      -------------------------#
#=================================================================================#


# for i in `ls *.sam`; do
#    F=`basename $i .sam`
#    echo -e "Sorting sample $F... start at `date`\n"
#    java -jar $PICARD SortSam I=$F.sam O=$F.srt.bam SO=coordinate
#    echo -e "Sample $F sroting by is done... at `date` \n"
 
#    echo -e "Marking duplicate of sample $F.... start at `date`\n"
#    java -jar $PICARD MarkDuplicates I=$F.srt.bam O=$F.srt.dedup.bam M=$F.metrics.txt
#    echo -e "Sample $F mark_duplicate is done... at `date` \n"
 
#    echo -e "Edit @RG label of sample $F... start at `date`\n"
#    java -jar $PICARD AddOrReplaceReadGroups I=$F.srt.dedup.bam O=$F.srt.dedup.rgmd.bam RGID=$F RGSM=Sample_$F RGLB=lib_$F RGPL=ILLUMINA RGPU=unit_$F
#    echo -e "Sample $F editing RG is done... at `date`\n"
 
#    echo -e "Building bam index of sample $F... start at `date`\n"
#    java -jar $PICARD BuildBamIndex I=$F.srt.dedup.rgmd.bam
#    echo -e "Sample $F buildbamindex is done... at `date`\n"
#  done



#=================================================================================#
#------------------------       Variants Calling       ---------------------------#
#=================================================================================#

#---- Local realignment, step 1: RealignTargetCreator
results/${ID}.bwa.target_interval.list: 
	${GATK} -T RealignerTargetCreator -nt 20 -nct 1 -R $REF -I $F.srt.dedup.rgmd.bam /
	-known $KNOWN/1000G_phase1.indels.b37.vcf.gz -known $KNOWN/Mills_and_1000G_gold_standard.indels.b37.vcf.gz -o $F.target_intervals.list
          

#---- Local realignment, step 2: IndelRealigner

results/${IND_ID_W_PE_SE}.bwa.${GENOME_NAME}.passed.realn.bam: 
	
#---- Base Recalibraiont, step 1: 



#---- Variant calling, HaplotypeCaller

#/pkg/java/jdk1.8.0/bin/java -Xmx16g -jar /home/u00yhl02/software/20150613/GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar -T HaplotypeCaller -R /work3/u00yhl02/references/ucsc.hg19.fasta -I /work3/u00yhl02/TSC_fastq/TSC002A/TSC002A_11132015_bwamem.marked.realigned.fixed.recal.indexed.bam --dbsnp /work3/u00yhl02/references/dbsnp_137.hg19.vcf -o /work3/u00yhl02/TSC_fastq/TSC002A/TSC002A_11132015_bwamem.haplotype.SnpIndel.vcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -A AlleleBalance -rf BadCigar
