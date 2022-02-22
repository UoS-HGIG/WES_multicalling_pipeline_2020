module load biobuilds
module load GATK/4.1.2
module load jdk/8


##collect variant quality metrics before filtration, after vqsr

INPUT=$1

##filtration
vcftools --gzvcf ${INPUT}.recalibrated.vcf.gz --minDP 8 --minGQ 20 --recode-INFO-all --recode --out ${INPUT}_GQ_DP_filt

# To define the correct treshold for HWE check total number of variants in the previous step and then divide 0.05 by it: E.G with 457003 variants > p.bonferroni 1.094*10-7
#Apply HWE filter for case/control cohort but may not be suitable for case-only study

#vcftools --vcf IBD_GQ_DP_filt.recode.vcf --recode --recode-INFO-all --hwe 1.094e-7 --out IBD_HWE


~/bin/meanGQ_filter.sh ${INPUT}_GQ_DP_filt.recode.vcf ${INPUT}_GQ_DP_meanGQ_filt.vcf


## Apply callrate filter, snp has to be genotyped in at least 88% of the cohort
vcftools --vcf ${INPUT}_HWE_meanGQ_filt.vcf --max-missing 0.88 --recode-INFO-all --recode --out ${INPUT}_postfilt

