#!/bin/bash
#SBATCH --ntasks-per-node=10
#SBATCH --time=04:00:00

OMNI=/mainfs/hgig/public/HUMAN_REFS/HG38/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf
ONE_THOUSAND=/mainfs/hgig/public/HUMAN_REFS/HG38/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf
DBSNP=/mainfs/hgig/public/HUMAN_REFS/HG38/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf
MILLS=/mainfs/hgig/public/HUMAN_REFS/HG38/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf
HAPMAP=/mainfs/hgig/public/HUMAN_REFS/HG38/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf
AXIOM=/mainfs/hgig/public/HUMAN_REFS/HG38/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz




module load GATK/4.1.2
module load biobuilds

cd /home/gc1a20/pirate/Sanger_IBD_adult/joint

#picard MergeVcfs \
#    I=split_vcf.list \
#    O=pibd_subFrom.sorted.vcf.gz \
#    D=/mainfs/hgig/public/HUMAN_REFS/HG38/REF_HLA/GRCh38_full_analysis_set_plus_decoy_hla.dict


gatk \
    MakeSitesOnlyVcf \
    --INPUT ps_indi.sorted.vcf.gz \
    --OUTPUT ps_indi.sitesOnly.sorted.vcf.gz


#########################
gatk --java-options "-Xms120g -Xmx160g" \
        VariantRecalibrator \
            -V ps_indi.sitesOnly.sorted.vcf.gz \
            -O ps_indi.sitesOnly.indel.recal \
            --tranches-file ps_indi.sitesOnly.indel.tranches \
                --rscript-file ps_indi.sitesOnly.indel.R \
                        --trust-all-polymorphic \
                            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff \
                                -tranche 100 -tranche 99.9 -tranche 99.0 -tranche 95.0 \
                                    -mode INDEL \
                                            -resource:mills,known=false,training=true,truth=true,prior=12 ${MILLS} \
                                                -resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${AXIOM} \
                                                    -resource:dbsnp,known=true,training=false,truth=false,prior=2 ${DBSNP}
gatk --java-options "-Xmx160g -Xms120g" \
          VariantRecalibrator \
            -V ps_indi.sitesOnly.sorted.vcf.gz \
            -O ps_indi.sitesOnly.snp.recal \
            --rscript-file ps_indi.sitesOnly.snp.R \
                --tranches-file ps_indi.sitesOnly.snp.tranches \
                                  --trust-all-polymorphic \
                                    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff \
                                        -tranche 100 -tranche 99.9 -tranche 99.0 -tranche 95.0 \
                                                    -mode SNP \
                                                            -resource:hapmap,known=false,training=true,truth=true,prior=15 ${HAPMAP} \
                                                                    -resource:omni,known=false,training=true,truth=true,prior=12 ${OMNI} \
                                                                            -resource:1000G,known=false,training=true,truth=false,prior=10 ${ONE_THOUSAND} \
                                                                                              -resource:dbsnp,known=true,training=false,truth=false,prior=7 ${DBSNP}



gatk --java-options "-Xmx16g -Xms12g" \
        ApplyVQSR \
            -O tmp.indel.recalibrated.vcf \
                -V ps_indi.sorted.vcf.gz \
                    --recal-file ps_indi.sitesOnly.indel.recal \
                        --tranches-file ps_indi.sitesOnly.indel.tranches \
                            --truth-sensitivity-filter-level 99.0 \
                                --create-output-variant-index true \
                                    -mode INDEL

gatk --java-options "-Xmx16g -Xms12g" \
        ApplyVQSR \
            -V tmp.indel.recalibrated.vcf \
                -O ps_indi.recalibrated.vcf.gz \
                    --recal-file ps_indi.sitesOnly.snp.recal \
                        --tranches-file ps_indi.sitesOnly.snp.tranches \
                            --truth-sensitivity-filter-level 99.0 \
                                --create-output-variant-index true \
                                    -mode SNP
