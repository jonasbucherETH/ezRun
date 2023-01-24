# for vcf app testing
require(ezRun)
# require(stringr)


# import test dataset
dataset <- read.delim("~/sushi_project_JB/data/test_vcf_dataset/dataset.tsv")
ragi_highcov_sa0001_1k.vcf <- read.table("~/sushi_project_JB/data/test_vcf_dataset/ragi_highcov_sa0001_1k.vcf.gz", quote="\"")
populations <- read.delim("~/sushi_project_JB/data/test_vcf_dataset/populations.tsv")
populations_txt <- read.delim("~/sushi_project_JB/data/test_vcf_dataset/populations.txt") # same as populations (txt instead of tsv)


output_dir <- "/home/jobucher/git/ezRun/output_data"

####### PCA
# library(gdsfmt)
# library(SNPRelate)
# 
# vcf_f <- file.path("/srv/gstore/projects/p1535/test_vcf_dataset/ragi_highcov_sa0001_1k.vcf.gz")
# 
# # convert vcf to gds   
# snpgdsVCF2GDS(vcf_f, file.path(output_dir, "snp.gds"),  method="biallelic.only")
# 
# # open gds
# genofile <- snpgdsOpen(file.path(output_dir, "snp.gds"))

library(vcfR) # to convert vcf to genind (ade4/adegenet)

# genind <- vcfR2genind(ragi_highcov_sa0001_1k.vcf)
#######


####### MDS
# file for mds
mds <- file.path(output_dir, "mds")

# run plink for distance matrix (mds)
prefix_mds <- file.path(output_dir, "plink")

# for testing
# cmd = paste("/usr/local/ngseq/packages/Tools/PLINK/2.00alpha/bin/plink2 -h")
# cmd = paste("module load Tools/PLINK/1.9beta6.21")
# cmd = paste("source /usr/local/ngseq/etc/lmod_profile")
# ezSystem(cmd)



cwd <- getwd()
setwdNew("/home/jobucher/sushi_project_JB/data")
on.exit(setwd(cwd), add = TRUE)
cmd <- paste("plink --vcf", "/srv/gstore/projects/p1535/test_vcf_dataset/ragi_highcov_sa0001_1k.vcf.gz", "--double-id", "--allow-extra-chr", "--cluster", "--mds-plot", 4 , "--out", prefix_mds)
# /srv/gstore/projects/p1535/test_vcf_dataset

mds_read <- read.csv(file.path(output_dir, "mds.mds"), sep="")


result <- ezSystem(cmd)
gc()

### cmdscale trial
# dist <- gen2dist(genind, biallelic = FALSE)
# dist <- dist(ragi_highcov_sa0001_1k.vcf)
# mds_cmdscale <- cmdscale(dist)


####### PCA

