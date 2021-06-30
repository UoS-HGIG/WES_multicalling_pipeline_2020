#!/bin/bash
#SBATCH --mem=16g
#SBATCH --time=16:00:00
#SBATCH -J hoardnshop




module load jdk/8.u181


REF=/mainfs/hgig/public/HUMAN_REFS/HG38/REF_HLA/GRCh38_full_analysis_set_plus_decoy_hla.fa
DBSITES=/mainfs/hgig/public/HUMAN_REFS/HG38/common_dbsnp_151.hg38.vcf


mkdir -p /home/gc1a20/pirate/Sanger_IBD_adult/pibd_joint
cd /home/gc1a20/pirate/Sanger_IBD_adult/pibd_joint
#########################
# GATK variants calling
#########################

java -Xms12G -Xmx16G -jar /local/software/GATK/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar \
               GenomicsDBImport \
                      --genomicsdb-workspace-path $2 \
                             --batch-size 30 \
                                    -L $1 \
                                    --merge-input-intervals true \
                                        --sample-name-map /home/gc1a20/bin/gatk-workflows/PIBD.sample_map \
                                                  --tmp-dir /ssdfs/users/gc1a20/tmp \
                                                         --reader-threads 10

java -Xms12g -Xmx16G -jar /local/software/GATK/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar \
        GenotypeGVCFs \
            -R ${REF} \
            -L $1 --only-output-calls-starting-in-intervals \
            --merge-input-intervals true \
                -O ${2}.vcf.gz \
                    -D ${DBSITES} \
                        -new-qual \
                                -V gendb:///home/gc1a20/pirate/Sanger_IBD_adult/pibd_joint/$2
