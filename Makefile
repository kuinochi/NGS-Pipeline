#---------------------------------------------------------------------------------#
# Makefile for running next generation sequencing workflow, 
# including mapping, variant calling
# Please check the needed config.mk before running the pipeline
#---------------------------------------------------------------------------------#

# Get user editable variables
include config.mk

# 	.marked.realigned.fixed.recal.indexed

# Default goal, print usage

help: 
	@printf "\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n\n" \
	"Please enter the steps you want to proceed" \
	"-------------------------------------------" \
	"IndexRef           --  Indexing a FASTA file by four tools:" \
	"                       bwa index, samtools faidx, bowtie2-build and CreateSequenceDictionary of picard." \
	"alignment          --  Align the reads to the reference genome by both bwa-mem and bowtie2" \
	"bwa                --  Align the reads to the reference genome by bwa-mem" \
	"bowtie2            --  Align the reads to the reference genome by bowtie2" \
	"Preprocess         --  Process SAM file to BAM file by Picard:" \
	"                       Sorting, Mark Duplicates, Edit @RG, Indexing" \
	"Process_gatk       --  Process BAM file by GATK, including:" \
	"                       RealignerTargetCreator, IndelRealigner, BaseRecalibrator, PrintReads" \
	"callvariants_sam   --  Variants Calling (SNP/INDELs) by samtools " \
	"callvairants_gatk  --  Variants Calling (SNP/INDELs) by samtools"


#=================================================================================#
#--------------------       Indexing reference FASTA       -----------------------#
#=================================================================================#

# Output files of indexing.
_SAMTOOLS_INDEX		= $(addsuffix .fai, ${REF_FA})
_BWA_INDEX_END	 	= .amb .ann .bwt .pac .sa
_BWA_INDEX			= $(foreach bw, ${_BWA_INDEX_END}, ${REF_FA}${bw})
_BOWTIE2_INDEX_END	= .1.bt2 .2.bt2 .3.bt2 .4.bt2 .rev.1.bt2 .rev.2.bt2
_BOWTIE2_INDEX		= $(foreach bt, ${_BOWTIE2_INDEX_END}, ${REF_FA}${bt})
_PICARD_DICT		= $(addsuffix .dict, ${REF_FA})
_Output_Folder = ./results


IndexRef: ${_SAMTOOLS_INDEX} ${_BWA_INDEX} ${_BOWTIE2_INDEX} ${_PICARD_DICT}
	@printf "\n%s\n%s\n\n" "The file:  ${REF_FA}  has alreadey been indexed." "Check the folder: ${REF_DIR} "

${_SAMTOOLS_INDEX}: ${REF_FA}
	@printf "\n%s\n" "#--- Start indexing by samtools..."
	\
	${SAMTOOL} faidx $<
	\
	@printf "%s\n" "#--- Samtools Indexing is done!"
	@touch ${_SAMTOOLS_INDEX}

${_BWA_INDEX}: ${REF_FA}
	@printf "\n%s\n" "#--- Start indexing by bwa..."
	\
	${BWA} index -a bwtsw $<
	\
	@printf "%s\n" "#--- BWA Indexing is done!"
	@touch ${_BWA_INDEX}

${_BOWTIE2_INDEX}: ${REF_FA} ${REF_FA}
	@printf "\n%s\n" "#--- Start indexing by bowtie2..."
	\
	${BOWTIE2_BUILD} $< $<
	\
	@printf "%s\n" "#--- Bowtie2 Indexing is done!"
	@touch ${_BOWTIE2_INDEX}	

${_PICARD_DICT}: ${REF_FA}
	@printf "\n%s\n" "#--- Start indexing by picard tools..."
	\
	java -jar ${PICARD} CreateSequenceDictionary R=$< O=$@
	\
	@printf "%s\n" "#--- Picard Indexing is done!"
	@touch ${_PICARD_DICT}


#=================================================================================#
#----------------------       Mapping to reference       -------------------------#
#=================================================================================#

alignment:	bwa bowtie2
	@printf "\n%s\n%s\n\n" \
	"The pair-end reads file in ${READ_DIR}  have all been mapped by both bwa-mem and bowtie2." \
	"Check the folder: ${READ_DIR} and ${_Output_Folder}/alignment "


#--- Align reads to the ref. genome by BWA

bwa: ${_Output_Folder}/alignment/*.bwa.sam
	@printf "\n%s\n%s\n\n" \
	"The pair-end reads file in ${READ_DIR}  have all been mapped by bwa-mem." \
	"Check the folder: ${READ_DIR} and ${_Output_Folder}/alignment "

${_Output_Folder}/alignment/*.bwa.sam: ${READ_DIR}
	@mkdir -p ${_Output_Folder}/alignment/
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start bwa-mem alignment..."

	@for i in ${READ_DIR}/*.r1.fastq;\
		do F=`basename $$i .r1.fastq` ;\
		printf "\n[%s] %s\n\n" "$$(date +%Y\/%m\/%d\ %T)" "Start align sample \"$${F}\" ...";\
		\
		${BWA} ${BWA_PARAM} \
			${REF_FA} \
			${READ_DIR}/$${F}.r1.fastq \
			${READ_DIR}/$${F}.r2.fastq \
			> ${_Output_Folder}/alignment/$${F}.${REF_NAME}.bwa.sam ;\
		\
	printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Sample \"$${F}\" alignment by bwa-mem is done!"; done\

	@printf "\n[%s] %s\n\n" "$$(date +%Y\/%m\/%d\ %T)" "Mapping by bwa-mem of all samples in the folder \"${READ_DIR}\" is done!"


#--- Align reads to the ref. genome by bowite2

bowtie2: 		${_Output_Folder}/alignment/*.bowtie2.sam
	@printf "\n%s\n%s\n\n" \
	"The pair-end reads file in ${READ_DIR}  have all been mapped by bowtie2." \
	"Check the folder: ${READ_DIR} and ${_Output_Folder}/alignment "

${_Output_Folder}/alignment/*.bowtie2.sam: ${READ_DIR} 
	@mkdir -p ${_Output_Folder}/alignment/
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start bowtie2 alignment..."

	@for i in ${READ_DIR}/*.r1.fastq; \
		do F=`basename $$i .r1.fastq` ;\
		printf "\n[%s] %s\n\n" "$$(date +%Y\/%m\/%d\ %T)" "Start align sample \"$${F}\" ..." ;\
		\
		${BOWTIE2} ${BOWTIE2_PARAM} \
		-1 ${READ_DIR}/$${F}.r1.fastq \
		-2 ${READ_DIR}/$${F}.r2.fastq \
		-S ${_Output_Folder}/alignment/$${F}.${REF_NAME}.bowtie2.sam;\
		\
		printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Sample \"$${F}\" alignment by bowtie2 is done!"; done\

	@printf "\n[%s] %s\n\n" "$$(date +%Y\/%m\/%d\ %T)" "Mapping by bowtie2 of all samples in the folder \"${READ_DIR}\" is done!"


#=================================================================================#
#----------------------      Alignment file process      -------------------------#
#=================================================================================#

Preprocess:	${_Output_Folder}/alignment/*.srt.dedup.rgmd.bai

	@printf "\n%s\n\n" "All file in ${_Output_Folder}/alignment have already been processed."


${_Output_Folder}/alignment/*.srt.dedup.rgmd.bai:  ${_Output_Folder}/alignment/*.srt.dedup.rgmd.bam
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start indexing bam files by Picard BuildBamIndex"
	@for i in ${_Output_Folder}/alignment/*.srt.dedup.rgmd.bam; do \
		F=`basename $${i} .srt.dedup.rgmd.bam`; \
		printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start indexing Sample $${F} by Picard BuildBamIndex..."; \
		\
		java -jar ${PICARD} BuildBamIndex I=$${i} ; \
		\
		printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Indexing Sample $${F} by Picard BuildBamIndex is done!"; \
	done
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Indexing bam file by Picard BuildBamIndex for all samples in ${_Output_Folder}/alignment is done!"


${_Output_Folder}/alignment/*.srt.dedup.rgmd.bam: ${_Output_Folder}/alignment/*.srt.dedup.bam
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start edit @RG tag in bam files by Picard AddOrReplaceReadGroups"
	@for i in ${_Output_Folder}/alignment/*.srt.dedup.bam; do \
		F=`basename $${i} .srt.dedup.bam`; \
		printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start Edit @RG for Sample $${F} by Picard AddOrReplaceReadGroups..."; \
		\
		java -jar ${PICARD} AddOrReplaceReadGroups \
			I=$${i} \
			O=${_Output_Folder}/alignment/$${F}.srt.dedup.rgmd.bam \
			RGID=$${F} RGSM=S_$${F} RGLB=Lib_$${F} RGPL=Illumina RGPU=Unit_$${F};\
		\
		printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "@RG is edited Sample $${F} by Picard AddOrReplaceReadGroups is done!";\
	done
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "@RG is edited by Picard AddOrReplaceReadGroups for all sample in ${_Output_Folder}/alignment!"


${_Output_Folder}/alignment/*.srt.dedup.bam: ${_Output_Folder}/alignment/*.srt.bam
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start mark duplicate sam files by Picard MarkDuplicates"
	@for i in ${_Output_Folder}/alignment/*.srt.bam; do \
		F=`basename $${i} .srt.bam`; \
		printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start mark duplicate Sample $${F} by Picard MarkDuplicates..."; \
		\
		java -jar ${PICARD} MarkDuplicates \
			I=$${i} \
			O=${_Output_Folder}/alignment/$${F}.srt.dedup.bam \
			M=${_Output_Folder}/alignment/$${F}.dedup.metrics.txt ; \
		\
		printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Mark duplicate Sample $${F} by Picard MarkDuplicates is done!"; \
	done
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Mark duplicate by Picard MarkDuplicates for all sample in ${_Output_Folder}/alignment is done!"


${_Output_Folder}/alignment/*.srt.bam:	${_Output_Folder}/alignment/*.sam
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start sorting sam files by Picard SortSam"
	@for i in ${_Output_Folder}/alignment/*.sam; do \
		F=`basename $${i} .sam`; \
		printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start sorting Sample $${F} by Picard SortSam..."; \
		\
		java -jar ${PICARD} SortSam \
			I=$${i} \
			O=${_Output_Folder}/alignment/$${F}.srt.bam \
			SO=coordinate; \
		\
		printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Sorting Sample $${F} by Picard SortSam is done!"; \
	done
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Sorting by Picard SortSam for all samples in ${_Output_Folder}/alignment is done!"

#=================================================================================#
#-----------------------        GATK Preprocess        ---------------------------#
#=================================================================================#



#=================================================================================#
#------------------------       Variants Calling       ---------------------------#
#=================================================================================#

callvariant_sam: ${_Output_Folder}/vcf/*.samtools.vcf
	@printf "\n%s\n%s\n\n" \
	"The bam file in ${_Output_Folder}/alignment have all conducted the variants calling by samtools." \
	"Check the folder: ${_Output_Folder}/vcf/ and ${_Output_Folder}/alignment/ "

${_Output_Folder}/vcf/*.samtools.vcf: ${_Output_Folder}/alignment/*.srt.dedup.rgmd.bai
	@mkdir -p ${_Output_Folder}/vcf
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start calling variants by samtools"
	@for i in ${_Output_Folder}/alignment/*.srt.dedup.rgmd.bam; do \
		F=`basename $${i} .srt.dedup.rgmd.bam`; \
		printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start calling variants for sample $${F} by samtools..."; \
		\
		${SAMTOOL} mpileup \
			-t DP,DV \
			-Buf ${REF_FA} \
			${_Output_Folder}/alignment/$${F}.srt.dedup.rgmd.bam | \
		${BCFTOOL} call \
			-vmOv \
			-o - > ${_Output_Folder}/vcf/$${F}.samtools.raw.vcf; \
		\
		printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Calling variants for sample $${F} by samtools is done!"; \
	done
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Calling variants for all sample in ${_Output_Folder}/alignment by samtools is done!"





#---- Local realignment, step 1: RealignTargetCreator
results/${ID}.bwa.target_interval.list: 
	${GATK} -T RealignerTargetCreator -nt 20 -nct 1 -R $REF -I $F.srt.dedup.rgmd.bam /
	-known $KNOWN/1000G_phase1.indels.b37.vcf.gz -known $KNOWN/Mills_and_1000G_gold_standard.indels.b37.vcf.gz -o $F.target_intervals.list
          

#---- Local realignment, step 2: IndelRealigner

results/${IND_ID_W_PE_SE}.bwa.${GENOME_NAME}.passed.realn.bam: 
	
#---- Base Recalibraiont, step 1: 



#---- Variant calling, HaplotypeCaller

#/pkg/java/jdk1.8.0/bin/java -Xmx16g -jar /home/u00yhl02/software/20150613/GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar -T HaplotypeCaller -R /work3/u00yhl02/references/ucsc.hg19.fasta -I /work3/u00yhl02/TSC_fastq/TSC002A/TSC002A_11132015_bwamem.marked.realigned.fixed.recal.indexed.bam --dbsnp /work3/u00yhl02/references/dbsnp_137.hg19.vcf -o /work3/u00yhl02/TSC_fastq/TSC002A/TSC002A_11132015_bwamem.haplotype.SnpIndel.vcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -A AlleleBalance -rf BadCigar
