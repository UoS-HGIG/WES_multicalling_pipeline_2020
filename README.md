# WES_multicalling_pipeline_2020
##preprocess.sh -->caller.sh -->joint_genotyping.sh -->vqsr.sh-->qc_based_filtration.sh -->variantEvaluation.sh

##removal of duplicate calls 

awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' jc.vcf |uniq > jc_uniq.vcf
