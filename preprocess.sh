#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --time=23:00:00
#SBATCH --mem=16g
#SBATCH -J aibdbb


module load biobuilds/2017.05
module load jdk/8.u181

REF=/mainfs/hgig/public/HUMAN_REFS/HG38/REF_HLA/GRCh38_full_analysis_set_plus_decoy_hla.fa
DBSITES=/mainfs/hgig/public/HUMAN_REFS/HG38/common_dbsnp_151.hg38.vcf

cd ~/pirate/Sanger_IBD_adult/Oct2020/
cd $1

#. variables

#ls -1 *L*.gz | grep '_1\.fq' | xargs gunzip -c | gzip > ${1}_1.fq.gz
#ls -1 *L*.gz | grep '_2\.fq' | xargs gunzip -c | gzip > ${1}_2.fq.gz

###########################
##Concatentate Multilane FQ Files if required
###########################
##pair1_fqs=($(grep 'pair1' variables | cut -d '=' -f 2))
##pair2_fqs=($(grep 'pair2' variables | cut -d '=' -f 2))
##
##for f in "${pair1_fqs[@]}"; do
##    zcat "${f//\"/}" >> ${1}_cat1.fq
##done
##
##for f in "${pair2_fqs[@]}"; do
##    zcat "${f//\"/}" >> ${1}_cat2.fq
##done
#
#
############################################
## aligning fastq files  with the reference genome using bwa-mem
## -t = Threads -M = flag shorter split hits as secondary
## -R = Readgroups -O = Gap open penalty -E =Gap extension penalty
############################################
#
#bwa mem \
#	-K 100000000 \
#    -M \
#	-R '@RG\tID:'${1}'_lane1\tSM:'${1}'\tPL:ILLUMINA\tLB:Library' \
#	$REF \
#	${1}_1.fq.gz ${1}_2.fq.gz \
#	> ${1}_aligned.sam
#
### Remove contactenated FastQ files if they were generated.
##rm ${1}_cat1.fq ${1}_cat2.fq
#
#
#######################
## convert SAM to BAM
#######################
#
samtools view -T $REF -bS -o ${1}_aligned.bam ${1}.cram

#rm ${1}.cram

module unload biobuilds/2017.05

######################
# Picard tools
######################

module load picard


# sort bam file
picard SortSam \
	INPUT=${1}_aligned.bam \
	OUTPUT=${1}_aligned_sorted.bam \
	SORT_ORDER=coordinate \
	TMP_DIR=tmp1 \
	VALIDATION_STRINGENCY=SILENT \
	MAX_RECORDS_IN_RAM=2500000

rm ${1}_aligned.bam

# mark duplicates
picard MarkDuplicates \
	INPUT=${1}_aligned_sorted.bam \
	METRICS_FILE=${1}_dup_metrics \
	OUTPUT=${1}_marked_dups_sorted.bam \
	TMP_DIR=tmp1 \
	VALIDATION_STRINGENCY=SILENT \
	MAX_RECORDS_IN_RAM=2500000

rm ${1}_aligned_sorted.bam
rm -fr tmp1

# Sort BAM file
picard SortSam \
	INPUT=${1}_marked_dups_sorted.bam \
	OUTPUT=${1}.DelDup.bam \
	SORT_ORDER=coordinate \
	TMP_DIR=tmp2 \
	VALIDATION_STRINGENCY=SILENT \
	MAX_RECORDS_IN_RAM=2500000

rm ${1}_marked_dups_sorted.bam
rm -fr tmp2

# Index BAM file
picard BuildBamIndex \
	INPUT=${1}.DelDup.bam \
	TMP_DIR=tmp3 \
	VALIDATION_STRINGENCY=SILENT

# Fix mate pair information by picard
picard FixMateInformation \
	INPUT=${1}.DelDup.bam \
	OUTPUT=${1}.GATK.fixedmateinfo.bam \
	SORT_ORDER=coordinate \
	TMP_DIR=tmp3 \
	VALIDATION_STRINGENCY=SILENT \
	MAX_RECORDS_IN_RAM=500000 \
	CREATE_INDEX=true



rm ${1}.DelDup.bam ${1}.DelDup.bai
rm -fr tmp3

module load jdk/8.u181
##########################################################################################
# GATK BQRS; NB: Indel realignment not required in using HaplotypeCaller downstream
##########################################################################################


# Recalibratiing base quality (longer: more than 60 minutes)
# if it fails, might be an old Illumina sample-> add -fixMisencodedQuals
java -Xmx16G -jar /local/software/GATK/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar BaseRecalibrator \
	-I ${1}.GATK.fixedmateinfo.bam \
	-R $REF \
	--known-sites $DBSITES \
	-O ${1}.recal_data.table

 Shorter: less than 15 minutes
java -Xmx16G -jar /local/software/GATK/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar ApplyBQSR \
	-R $REF \
	-I ${1}.GATK.fixedmateinfo.bam \
	--bqsr-recal-file ${1}.recal_data.table \
	-O ${1}.GATK.recal.bam \

#############################
 Remove reference files, Penultimate BAM. Picard tmp directories
#############################

rm ${1}.GATK.fixedmateinfo.bam ${1}.GATK.fixedmateinfo.bai

picard AddOrReplaceReadGroups \
    CREATE_INDEX=true \
    I=${1}.GATK.recal.bam \
    O=${1}.GATK.recal.reh.bam \
    RGID=99 \
    RGLB=sanger202010 \
    RGPL=ILLUMINA \
    RGPU=sureselect5 \
    RGSM=$1i \
    RGDS=AIBD
