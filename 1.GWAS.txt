#!/bin/bash

#############
# Run GWAS
#############
for phe in UKB_COVID_Inf UKB_COVID_Sev UKB_T2D_adjBMI UKB_BMI UKB_WHR UKB_WC;do
bolt \
    --bed=~/genotype/ukb_cal_chr{1:22}_v2.bed \
    --bim=~/genotype/ukb_snp_chr{1:22}_v2.bim \
    --fam=~/genotype/ukb_snp_chr1_v2.fam \
    --exclude=~/phenotype/autosome_missing_gt_0.1.txt \
    --exclude=~/phenotype/autosome_maf_lt_0.001.txt \
    --phenoFile=~/phenotype.txt \
    --remove=~/removeID_missingGENO_Asthma.txt \
    --phenoCol="$phe" \
    --covarFile=~/covariate/ukb4777.processed_and_post2.plinkPCs.tab \
    --covarCol=cov_ASSESS_CENTER \
    --covarCol=cov_GENO_ARRAY \
    --covarMaxLevels=30 \
    --covarCol=cov_SEX \
    --qCovarCol=cov_AGE \
    --qCovarCol=cov_AGE_SQ \
    --qCovarCol=PC{1:20} \
    --LDscoresFile=~/phenotype/LDSCORE.1000G_EUR.tab.gz \
    --geneticMapFile=~/phenotype/genetic_map_hg19.txt.gz \
    --lmmForceNonInf \
    --numThreads=10 \
    --statsFile=~/BOLT_"$phe".stats.gz \
    --bgenFile=~/imputation/ukb_imp_chr{1:22}_v3.bgen \
    --bgenMinMAF=0.01 \
    --bgenMinINFO=0.8 \
    --sampleFile=~/phenotype/ukb_imp_s487406.sample \
    --statsFileBgenSnps=~/1.GWAS/BOLT_"$phe".bgen.stats.gz \
    --verboseStats
done