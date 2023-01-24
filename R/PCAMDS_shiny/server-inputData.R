inputDataReactive <- reactive({

  
  cat("input")
  
  ### load data for testing
  vcf <- read.vcfR("~/sushi_project_JB/data/test_vcf_dataset/ragi_highcov_sa0001_1k.vcf.gz")
  genind <- vcfR2genind(vcf)
  grouping_vars <- read.delim("~/sushi_project_JB/data/test_vcf_dataset/populations.txt")
  # pop(genind) <- populations_txt$Population
  
  X <- scaleGen(genind, NA.method="mean")
  
  
  # pca <- dudi.pca(X, center = TRUE, scale = TRUE, scan = FALSE)
  # #
  # # mds <- read.csv("~/git/ezRun/output_data/plink_mds.mds", sep="")
  # #
  # # ### load data
  # # # pca <- readRDS("PCA.rds")
  # # # grouping_vars <- readRDS("grouping_vars.rds")
  # # # mds <- read.csv("plink.mds", sep="")
  # #
  # #
  # # ### both testing & normal setup
  # n_pcs <- pca$nf # number of (> 0) principal components
  # n_grouping <- ncol(grouping_vars)
  # 
  # eig_sum <- sum(pca$eig)
  # pca_varprop <- pca$eig/eig_sum
  # pca_varprop <- pca_varprop[1:n_pcs]
  # 
  # pca_tab <- data.frame(grouping_vars, pca$li, stringsAsFactors = FALSE, row.names = NULL)
  # 
  # tab_varprop <- as.data.frame(t(pca_varprop), stringsAsFactors = FALSE)
  # PC_indeces <- seq(1+n_grouping, ncol(pca_tab))
  # 
  # for (i in 1:n_pcs){
  #   colnames(pca_tab)[i+n_grouping] <- paste0("PC", i)
  #   colnames(tab_varprop)[i] <- paste0("PC", i)
  # }
  # 
  # # all PCs in array for selecting input
  # pc_list <- colnames(pca_tab)[-(1:n_grouping)]
  
  
  # 
  # n_dim <- ncol(mds)-ncol(grouping_vars)-1   # number of dimensions kept in mds
  # mds_tab <- data.frame(grouping_vars, mds[, (ncol(grouping_vars)+2):ncol(mds)], stringsAsFactors = FALSE, row.names = NULL)
  
  # n_pcs = "please print this"
  # print(n_pcs)
  
  return(list(
    # "pca" = pca,
    # "n_pcs" = n_pcs,
    # "n_grouping" = n_grouping,
    # "pc_list" = pc_list,
    # "pca_varprop" = pca_varprop,
    # "pca_tab" = pca_tab,
    # "tab_varprop" = tab_varprop,
    "vcf" = vcf,
    "genind" = genind,
    "grouping_vars" = grouping_vars,
    "X" = X
    )
  )
  
  
})