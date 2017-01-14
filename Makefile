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
	"GATK_preprocess    --  Process BAM file by GATK, including:" \
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
#_PICARD_DICT		= $(addsuffix .dict, ${REF_FA})
_PICARD_DICT   = $(patsubst %.fasta, %.dict, ${REF_FA})
_Output_Folder = ./results



IndexRef: ${_SAMTOOLS_INDEX} ${_BWA_INDEX} ${_BOWTIE2_INDEX} ${_PICARD_DICT}
	@printf "\n%s\n%s\n\n" "The file:  ${REF_FA}  has alreadey been indexed." "Check the folder: ${REF_DIR}" 

${_SAMTOOLS_INDEX}: ${REF_FA}
	@printf "\n%s\n" "#--- Start indexing by samtools..."
	
	${SAMTOOL} faidx $<
	
	@printf "%s\n" "#--- Samtools Indexing is done!"
	@touch ${_SAMTOOLS_INDEX}

${_BWA_INDEX}: ${REF_FA}
	@printf "\n%s\n" "#--- Start indexing by bwa..."
	
	${BWA} index -a bwtsw $<
	
	@printf "%s\n" "#--- BWA Indexing is done!"
	@touch ${_BWA_INDEX}

${_BOWTIE2_INDEX}: ${REF_FA} ${REF_FA}
	@printf "\n%s\n" "#--- Start indexing by bowtie2..."

	${BOWTIE2_BUILD} $< $<
	
	@printf "%s\n" "#--- Bowtie2 Indexing is done!"
	@touch ${_BOWTIE2_INDEX}	

${_PICARD_DICT}: ${REF_FA}
	@printf "\n%s\n" "#--- Start indexing by picard tools..."
	
	java -jar ${PICARD} CreateSequenceDictionary R=$< O=$@
	
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
	@printf "\n%s\n\n" "All file in ${_Output_Folder}/alignment have already been preprocessed."
	@touch ${_Output_Folder}/alignment/*.srt.dedup.rgmd.bai

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
	@touch ${_Output_Folder}/alignment/*.srt.dedup.rgmd.bai


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

GATK_preprocess: ${_Output_Folder}/alignment/*.realign.recal.bai
	@touch ${_Output_Folder}/alignment/*.realign.recal.bai
	@printf "\n%s\n%s\n\n" "The bam file in ${_Output_Folder}/alignment have all conducted the GATK processing:" \
	"RealignerTargetCreator --> IndelRealigner --> BaseRecalibrator --> PrintReads" \
	"You can go next step (e.g. variants calling) now."

${_Output_Folder}/alignment/*.realign.recal.bai: ${_Output_Folder}/alignment/*.base.recali.table
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start GATK PrintReads"
	@for i in ${_Output_Folder}/alignment/*.base.recali.table; do \
		F=`basename $${i} .base.recali.table`; \
		printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start GATK PrintReads for sample $${F}..."; \
		\
		java -jar ${GATK} -T PrintReads \
			-R ${REF_FA} \
			-nct ${GATK_NMT} \
			-I ${_Output_Folder}/alignment/$${F}.realign.bam \
			-BQSR ${_Output_Folder}/alignment/$${F}.base.recali.table \
			-o ${_Output_Folder}/alignment/$${F}.realign.recal.bam; \
		\
		printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "GATK PrintReads for sample $${F} is done!"; \
		done
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "GATK PrintReads for all sample in ${_Output_Folder} is done!"


${_Output_Folder}/alignment/*.base.recali.table: ${MILLS_KG_INDEL} ${DBSNP_138} ${KGPHASE1_INDEL} ${_Output_Folder}/alignment/*.realign.bai
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start GATK BaseRecalibrator"
	@for i in ${_Output_Folder}/alignment/*.realign.bam; do \
		F=`basename $${i} .realign.bam`; \
		printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start GATK BaseRecalibrator for sample $${F}..."; \
		\
		java -jar ${GATK} -T BaseRecalibrator \
			-R ${REF_FA} \
			-nct ${GATK_NMT} \
			-nt 1 \
			-I ${_Output_Folder}/alignment/$${F}.realign.bam \
			-knownSites ${DBSNP_138} \
			-knownSites ${MILLS_KG_INDEL} \
			-knownSites ${KGPHASE1_INDEL} \
			-o ${_Output_Folder}/alignment/$${F}.base.recali.table; \
		\
		printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "GATK BaseRecalibrator for sample $${F} is done!"; \
		done
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "GATK BaseRecalibrator for all sample in ${_Output_Folder} is done!"


${_Output_Folder}/alignment/*.realign.bai: ${MILLS_KG_INDEL} ${KGPHASE1_INDEL} ${_Output_Folder}/alignment/*.target_intervals.list
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start GATK IndelRealigner"
	@for i in ${_Output_Folder}/alignment/*.target_intervals.list; do \
		F=`basename $${i} .target_intervals.list`; \
		printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start GATK IndelRealigner for sample $${F}..."; \
		\
		java -jar ${GATK} -T IndelRealigner \
			-R ${REF_FA} \
			-I ${_Output_Folder}/alignment/$${F}.srt.dedup.rgmd.bam \
			-known ${MILLS_KG_INDEL} \
			-known ${KGPHASE1_INDEL} \
			-targetIntervals ${_Output_Folder}/alignment/$${F}.target_intervals.list \
			-o ${_Output_Folder}/alignment/$${F}.realign.bam \
			--filter_bases_not_stored; \
		\
		printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "GATK IndelRealigner for sample $${F} is done!"; \
		done
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "GATK IndelRealigner for all sample in ${_Output_Folder} is done!"


${_Output_Folder}/alignment/*.target_intervals.list: ${MILLS_KG_INDEL} ${KGPHASE1_INDEL} Preprocess
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start GATK RealignerTargetCreator"
	@for i in ${_Output_Folder}/alignment/*.srt.dedup.rgmd.bam; do \
		F=`basename $${i} .srt.dedup.rgmd.bam`; \
		printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start GATK RealignerTargetCreator for sample $${F}..."; \
		\
		java -jar ${GATK} -T RealignerTargetCreator \
			-nt ${GATK_NMT} \
			-nct 1 \
			-R ${REF_FA} \
			-I ${_Output_Folder}/alignment/$${F}.srt.dedup.rgmd.bam \
			-known ${KGPHASE1_INDEL} \
			-known ${MILLS_KG_INDEL} \
			-o ${_Output_Folder}/alignment/$${F}.target_intervals.list; \
		\
		printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "GATK RealignerTargetCreator for sample $${F} is done!"; \
		done
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "GATK RealignerTargetCreator for all sample in ${_Output_Folder} is done!"




#=================================================================================#
#------------------------       Variants Calling       ---------------------------#
#=================================================================================#

callvariant_sam: ${_Output_Folder}/vcf/*.samtools.raw.vcf
	@printf "\n%s\n%s\n\n" \
	"The bam file in ${_Output_Folder}/alignment have all conducted the variants calling by samtools." \
	"Check the folder: ${_Output_Folder}/vcf/ and ${_Output_Folder}/alignment/ "

${_Output_Folder}/vcf/*.samtools.raw.vcf: ${_Output_Folder}/alignment/*.${REF_NAME}.*.realign.recal.bai
	@mkdir -p ${_Output_Folder}/vcf
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start calling variants by samtools"
	@for i in ${_Output_Folder}/alignment/*.${REF_NAME}.*.realign.recal.bam; do \
		F=`basename $${i} .realign.recal.bam`; \
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


callvariant_gatk: ${_Output_Folder}/vcf/*.gatk.raw.vcf

${_Output_Folder}/vcf/*.gatk.raw.vcf: ${_Output_Folder}/alignment/*.${REF_NAME}.*.realign.recal.bai
	@mkdir -p ${_Output_Folder}/vcf
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start GATK HaplotypeCaller"
	@for i in ${_Output_Folder}/alignment/*.${REF_NAME}.*.realign.recal.bam; do \
		F=`basename $${i} .realign.recal.bam`; \
		printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start GATK HaplotypeCaller for sample $${F}..."; \
		\
		java -jar ${GATK} -T HaplotypeCaller \
		-R ${REF_FA} \
		-I ${_Output_Folder}/alignment/$${F}.srt.dedup.rgmd.bam \
		-o ${_Output_Folder}/vcf/$${F}.gatk.raw.vcf \
		-stand_call_conf 50.0 \
		-stand_emit_conf 10.0 \
		-A AlleleBalance \
		-rf BadCigar \
		--dbsnp ${DBSNP_138}; \
		\
		printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "GATK HaplotypeCaller for sample $${F} is done!"; \
	done
	@printf "\n[%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Calling variants for all sample in ${_Output_Folder}/alignment by GATK HaplotypeCaller is done!"
