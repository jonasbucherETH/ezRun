#!/bin/bash

source /usr/local/ngseq/etc/lmod_profile
module add Tools/PLINK/1.9beta6.21

plink --vcf ~/sushi_project_JB/data/test_vcf_dataset/ragi_highcov_sa0001_1k.vcf.gz --double-id --allow-extra-chr --cluster --mds-plot 5 --out plink_3101

plink --vcf ~/sushi_project_JB/data/test_vcf_dataset/ragi_highcov_sa0001_1k.vcf.gz --double-id --allow-extra-chr --distance square --out plink_3101

# plink --vcf ~/sushi_project_JB/data/test_vcf_dataset/ragi_highcov_sa0001_1k.vcf.gz --double-id --allow-extra-chr --cluster --mds-plot 2 --out mds_dim2
