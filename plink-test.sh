#!/bin/bash
set -eu  ## also fail if there is an unset variable
set -o pipefail
set -x 
set -o history -o histexpand



## run with
##  /srv/GT/analysis/course_sushi_testCases/nightlyTest.sh > /srv/GT/analysis/course_sushi_testCases/nightlyTest.out 2> /srv/GT/analysis/course_sushi_testCases/nightlyTest.err


# dieWithMail() {
#         echo -e "The failed commandline is:\n $1" | mail -s "ERROR: sushi test run script failed" $mailaddress
#         exit 1
# }

# sushiCmd="/usr/local/ngseq/opt/Ruby_Gems/ruby/1.9.1//bin/sushi_fabric -I /srv/GT/analysis/course_sushi"
# readDataset="/srv/GT/analysis/course_sushi/public/gstore/projects/p1000/yeast_10k/dataset.tsv"
# alignedDataset="/srv/GT/analysis/course_sushi/public/gstore/projects/p1000/yeast10k_STAR_622_2015-05-04--19-28-50/dataset.tsv"
# countDataset1="/srv/GT/analysis/course_sushi/public/gstore/projects/p1000/yeast_10k_RSEM_2015-05-04--19-15-08/dataset.tsv"
# countDataset2="/srv/GT/analysis/course_sushi/public/gstore/projects/p1000/yeast_10k_CountOverlaps_638_2015-05-04--19-34-02/dataset.tsv"
# cd /srv/GT/analysis/course_sushi

dataset <- read.delim("~/sushi_project_JB/data/test_vcf_dataset/dataset.tsv")
ragi_highcov_sa0001_1k.vcf <- read.table("~/sushi_project_JB/data/test_vcf_dataset/ragi_highcov_sa0001_1k.vcf.gz", quote="\"")
populations <- read.delim("~/sushi_project_JB/data/test_vcf_dataset/populations.tsv")
populations_txt <- read.delim("~/sushi_project_JB/data/test_vcf_dataset/populations.txt") # same as populations (txt instead of tsv)

plink --vcf ragi_highcov_sa0001_1k.vcf --double-id --allow-extra-chr --cluster --mds-plot 4


# $sushiCmd --class FastqcApp -p 1000 --run -d $readDataset --next_dataset_name testcase_fastqc || dieWithMail " !! !:* "