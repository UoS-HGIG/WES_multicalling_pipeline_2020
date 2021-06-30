#!/bin/bash
#SBATCH --mem=16g
#SBATCH --time=08:00:00




cd /home/gc1a20/pirate/iss1g18/FASTQS/ADULT_RELATIVES
cd $1

module load picard
module load jdk/8.u181
REF=/mainfs/hgig/public/HUMAN_REFS/HG38/REF_HLA/GRCh38_full_analysis_set_plus_decoy_hla.fa

java -Xmx16G -jar /local/software/GATK/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar HaplotypeCaller \
	-ERC GVCF \
	-R $REF \
	-I ${1}.GATK.recal.bam \
	-L ~/ref/target_ref/SSV5_6_union_gc_pad150bp.interval_list \
  -O ${1}.onTargetPad150bp.g.vcf.gz \
	--dont-use-soft-clipped-bases

##optional for individual calling without limit of target region
java -Xmx16G -jar /local/software/GATK/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar HaplotypeCaller \
	-ERC GVCF \
	-R $REF \
	-I ${1}.GATK.recal.bam \
  -O ${1}.g.vcf.gz \
	--dont-use-soft-clipped-bases
java -Xmx16G -jar /local/software/GATK/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar GenotypeGVCFs \
        -R $REF \
        -V ${1}.g.vcf.gz \
        -O ${1}.vcf.gz
