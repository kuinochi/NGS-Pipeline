#---------------------------------------------------------------------------------#
# Makefile for running next generation sequencing workflow, 
# including mapping, variant calling
# Please check the needed config.mk before running the pipeline
#---------------------------------------------------------------------------------#

# Get user editable variables
include config.mk

#
# Default goal, print usage
# 

help: 
	@printf "\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n\n" \
	"Please enter the steps you want to proceed" \
	"-------------------------------------------" \
	"IndexRef           --  Indexing a FASTA file by four tools:" \
	"                       bwa index, samtools faidx, bowtie2-build and CreateSequenceDictionary of picard." \
	"Alignment          --  Align the reads to the reference genome by bwa-mem" \
	"Preprocess         --  Process SAM file to BAM file by Picard:" \
	"                       Sorting, Mark Duplicates, Edit @RG, Indexing" \
	"GATK_preprocess    --  Process BAM file by GATK, including:" \
	"                       RealignerTargetCreator, IndelRealigner, BaseRecalibrator, PrintReads" \
	"callvariants_sam   --  Variants Calling (SNP/INDELs) by samtools " \
	"callvariants_gatk  --  Variants Calling (SNP/INDELs) by samtools"
 	
SPFH = ${_Output_Folder}/${sNAME}.${REF_NAME}

 
all: ${SPFH}.samtools.raw.vcf ${SPFH}.gatk.raw.vcf
	@printf "[PIPELINE CMD][%s] checking VCF file for ${sNAME} \n\n%s\n%s\n%s\n" \
	"$$(date +%Y\/%m\/%d\ %T)" \
	"There are two VCF files named \
	${sNAME}.${REF_NAME}.gatk.vcf and \
	${sNAME}.${REF_NAME}.samtools.vcf in ${_Output_Folder}/" \
	"Please check again." 

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




IndexRef: ${_SAMTOOLS_INDEX} ${_BWA_INDEX} ${_BOWTIE2_INDEX} ${_PICARD_DICT}
	@printf "\n%s\n%s\n\n" "The file:  ${REF_FA}  has alreadey been indexed." "Check the folder: ${REF_DIR}" 

${_SAMTOOLS_INDEX}: ${REF_FA}
	@printf "\n%s\n" "#--- Start indexing by samtools..."
	
	${SAMTOOL} faidx $<
	
	@printf "%s\n" "#--- Samtools Indexing is done!"
	@touch ${_SAMTOOLS_INDEX}

${_BWA_INDEX}: ${REF_FA}
	@printf "\n%s\n" "#--- Start indexing by bwa... # -a bwtsw"
	
	${BWA} index $<
	 
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
#-------------		Align reads to the ref. genome by BWA		------------------#
#=================================================================================#



Alignment: ${SPFH}.bwa.sam
	@printf "[PIPELINE CMD][%s] checking BAM file for ${sNAME} \n\n%s\n%s\n%s\n" \
	"$$(date +%Y\/%m\/%d\ %T)" \
	"There is a BAM file named ${sNAME}.${REF_NAME}.bwa.sam in ${_Output_Folder}/" \
	"Please check again." 

${SPFH}.bwa.sam: ${READ_1} ${READ_2} ${REF_FA} ${_PICARD_DICT}
	@mkdir -p ${_Output_Folder}/
	@printf "[PIPELINE CMD][%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start bwa-mem alignment for Sample: ${sNAME}"
	${BWA} ${BWA_PARAM} \
	-R "@RG\tID:${sRGID}\tSM:${sRGSM}\tPL:${sRGPL}" \
	${REF_FA} \
	${READ_1} \
	${READ_2} \
	> ${_Output_Folder}/${sNAME}.${REF_NAME}.bwa.sam
	@printf "[PIPELINE CMD][%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "End bwa-mem alignment for ${sNAME}"


#=================================================================================#
#----------------------      Alignment file process      -------------------------#
#=================================================================================#

Preprocess: ${SPFH}.bwa.sorted.deduped.bai
	@printf "[PIPELINE CMD][%s] checking deduped BAM file for ${sNAME} \n\n%s\n%s\n%s\n" \
	"$$(date +%Y\/%m\/%d\ %T)" \
	"There is a BAM file named ${sNAME}.${REF_NAME}.bwa.sorted.deduped.bam in ${_Output_Folder}/" \
	"Please check again." ;\
	
${SPFH}.bwa.sorted.deduped.bai: ${SPFH}.bwa.sorted.deduped.bam
	@printf "[PIPELINE CMD][%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start Build BAM index for deduped ${sNAME}"
	java -jar ${PICARD} BuildBamIndex \
		I=$< 
	@printf "[PIPELINE CMD][%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "End Build BAM index for deduped ${sNAME}"


${SPFH}.bwa.sorted.deduped.bam:	${SPFH}.bwa.sorted.bam	
	@printf "[PIPELINE CMD][%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Strat MarkDup for ${sNAME}"
	java -jar ${PICARD} MarkDuplicates \
		I=$< \
		O=${SPFH}.bwa.sorted.deduped.bam \
		M=${SPFH}.dedup.metrics.txt;
	@printf "[PIPELINE CMD][%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "End MarkDup for ${sNAME}"

${SPFH}.bwa.sorted.bam:	${SPFH}.bwa.sam
	@printf "[PIPELINE CMD][%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start sorting for ${sNAME}"
	java -jar ${PICARD} SortSam \
		I=$< \
		O=${_Output_Folder}/${sNAME}.${REF_NAME}.bwa.sorted.bam \
		SO=coordinate; 
	@printf "[PIPELINE CMD][%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "End Sorting for ${sNAME}"

#=================================================================================#
#-----------------------        GATK Preprocess        ---------------------------#
#=================================================================================#

GATK_preprocess: ${SPFH}.realigned.recaled.bai
	@printf "[PIPELINE CMD][%s] checking realigned and recaled BAM file for ${sNAME} \n\n%s\n%s\n%s\n" \
	"$$(date +%Y\/%m\/%d\ %T)" \
	"There is a BAM file named ${sNAME}.${REF_NAME}.realigned.recaled.bam in ${_Output_Folder}/" \
	"Please check again." ;\

${SPFH}.realigned.recaled.bai: ${SPFH}.recali.table
	@printf "[PIPELINE CMD][%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start PrintReads for ${sNAME}"
	java -jar ${GATK} -T PrintReads \
		-R ${REF_FA} \
		-nct ${GATK_NMT} -nt 1 \
		-I ${SPFH}.realigned.bam \
		-BQSR ${SPFH}.recali.table \
		-o ${SPFH}.realigned.recaled.bam
	@printf "[PIPELINE CMD][%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "End PrintReads for ${sNAME}"


${SPFH}.recali.table:	${MILLS_KG_INDEL} ${DBSNP_138} ${KGPHASE1_INDEL} ${SPFH}.realigned.bai
	@printf "[PIPELINE CMD][%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start BaseRecalibrator for ${sNAME}"
	java -jar ${GATK} -T BaseRecalibrator \
		-nct ${GATK_NMT} -nt 1 \
		-R ${REF_FA} \
		-I ${SPFH}.realigned.bam \
		-knownSites ${MILLS_KG_INDEL} \
		-knownSites ${KGPHASE1_INDEL} \
		-knownSites ${DBSNP_138} \
		-o ${SPFH}.recali.table
	@printf "[PIPELINE CMD][%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "End BaseRecalibrator for ${sNAME}"

${SPFH}.realigned.bai: ${MILLS_KG_INDEL} ${KGPHASE1_INDEL} ${SPFH}.target_intervals.list
	@printf "[PIPELINE CMD][%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start IndelRealigner for ${sNAME}"
	java -jar ${GATK} -T IndelRealigner \
		-R ${REF_FA} \
		-I ${SPFH}.bwa.sorted.deduped.bam \
		-known ${MILLS_KG_INDEL} \
		-known ${KGPHASE1_INDEL} \
		-targetIntervals ${SPFH}.target_intervals.list \
		-o ${SPFH}.realigned.bam \
		--filter_bases_not_stored;
	@printf "[PIPELINE CMD][%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "End IndelRealigner for ${sNAME}"

${SPFH}.target_intervals.list: ${MILLS_KG_INDEL} ${KGPHASE1_INDEL} ${SPFH}.bwa.sorted.deduped.bai
	@printf "[PIPELINE CMD][%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start RealignerTargetCreator for ${sNAME}"
	java -jar ${GATK} -T RealignerTargetCreator \
		-nt ${GATK_NMT} \
		-nct 1 \
		-R ${REF_FA} \
		-I ${SPFH}.bwa.sorted.deduped.bam \
		-known ${MILLS_KG_INDEL} \
		-known ${KGPHASE1_INDEL} \
		-o ${SPFH}.target_intervals.list
	@printf "[PIPELINE CMD][%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "End RealignerTargetCreator for ${sNAME}"
	

#=================================================================================#
#------------------------       Variants Calling       ---------------------------#
#=================================================================================#


callvariants_sam: ${SPFH}.samtools.raw.vcf
	@printf "[PIPELINE CMD][%s] checking VCF by samtools for ${sNAME} \n\n%s\n%s\n%s\n" \
	"$$(date +%Y\/%m\/%d\ %T)" \
	"There is a VCF file named ${sNAME}.${REF_NAME}.samtools.raw.vcf in ${_Output_Folder}/" \
	"Please check again." ;\

${SPFH}.samtools.raw.vcf: ${SPFH}.realigned.recaled.bai
	@printf "[PIPELINE CMD][%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start calling by samtools for ${sNAME}"
	${SAMTOOL} mpileup -t DP,DV -Buf ${REF_FA} ${SPFH}.realigned.recaled.bam | \
	${BCFTOOL} call -vmOv -o - > ${SPFH}.samtools.raw.vcf
	@printf "[PIPELINE CMD][%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "End calling by samtools for ${sNAME}"

callvariants_gatk: ${SPFH}.gatk.raw.vcf
	@printf "[PIPELINE CMD][%s] checking VCF by GATK for ${sNAME} \n\n%s\n%s\n%s\n" \
	"$$(date +%Y\/%m\/%d\ %T)" \
	"There is a VCF file named ${sNAME}.${REF_NAME}.gatk.raw.vcf in ${_Output_Folder}/" \
	"Please check again." ;\

${SPFH}.gatk.raw.vcf: ${SPFH}.realigned.recaled.bai
	@printf "[PIPELINE CMD][%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "Start calling by HaplotypeCaller for ${sNAME}"
	java -jar ${GATK} -T HaplotypeCaller \
		-R ${REF_FA} \
		-I ${SPFH}.realigned.recaled.bam \
		-o ${SPFH}.gatk.raw.vcf \
		-stand_call_conf 50.0 -stand_emit_conf 10.0 \
		-A AlleleBalance -rf BadCigar \
		--dbsnp ${DBSNP_138}
	@printf "[PIPELINE CMD][%s] %s\n" "$$(date +%Y\/%m\/%d\ %T)" "End calling by HaplotypeCaller for ${sNAME}"