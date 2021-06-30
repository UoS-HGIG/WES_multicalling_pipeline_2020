#!/bin/bash
#SBATCH --time=04:00:00



module load picard
module load GATK/4.1.2
module load jdk/8.u181

DBSITES=/mainfs/hgig/public/HUMAN_REFS/HG38/common_dbsnp_151.hg38.vcf
REF=/mainfs/hgig/public/HUMAN_REFS/HG38/REF_HLA/GRCh38_full_analysis_set_plus_decoy_hla.fa

cd /home/gc1a20/pirate/Sanger_IBD_adult/IBD_AD_joint_pad
##gatk beta version
#gatk VariantEval \
#        -R ${REF} \
#           -O $2 \
#           --eval $1

##picard
picard CollectVariantCallingMetrics DBSNP=${DBSITES} \
    I=$1 \
    O=$2

##peddy of pedigress, gender, and ethnicity estimation
##hereby installed locally conda

#peddy -p 4 --sites hg38 --plot --prefix IBD_202103 ~/pirate/Sanger_IBD_adult/IBD_AD_joint_pad/IBD.vqsr.GQ_DP_meanGQ_p150_intersect.recode.vcf.gz IBD_none.ped
