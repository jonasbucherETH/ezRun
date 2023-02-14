inputDataReactive <- reactive({
# inputDataReactive <- {
  

  # Load data ---------------------------
  # pca <- readRDS("PCA.rds")
  # groupingVariables <- readRDS("grouping_vars.rds")
  # vcfRaw <- read.vcfR("vcf.rds")
  # vcfGenind <- vcfR2genind(vcfRaw)
  # datasetScaled <- scaleGen(vcfGenind, NA.method="mean")
  # 
  # mds <- read.csv("plink.mds", sep="")
  # 
  # distanceMatrixTSNE <- read_tsv("plink.dist", col_names = F)
  
  ### for testing
  # vcfA <- read.vcfR("~/git/ezRun/R/PCAMDS_shiny/data/sample_dataset4jonas_20230206/A_Eall.vcf.gz")
  # vcfB <- read.vcfR("~/git/ezRun/R/PCAMDS_shiny/data/sample_dataset4jonas_20230206/B_Eall.vcf.gz")
  # vcfA_corrected <- read.vcfR("~/git/ezRun/R/PCAMDS_shiny/data/sample_dataset4jonas_20230206/A_Eall.pruned.phased.sample_name_corrected.landraces.vcf.gz")
  # vcfB_corrected <- read.vcfR("~/git/ezRun/R/PCAMDS_shiny/data/sample_dataset4jonas_20230206/B_Eall.pruned.phased.sample_name_corrected.landraces.vcf.gz")
  # vcfGenind <- vcfR2genind(vcfA_corrected)
  # pop(vcfGenind) <- c(groupingVariables[,2])
  
  # colnames(vcfRaw@gt)
  
  ## PCA
  vcfRaw <- read.vcfR("data/ragi_highcov_sa0001_1k.vcf.gz")
  # vcfR@gt has sample IDs as colnames (col 1 = FORMAT)
  vcfGenind <- vcfR2genind(vcfRaw)
  groupingVariables <- read.delim("data/populations.txt") # read.table or read.delim (used before) -> delim should be fine
  # groupingVariables <- groupingVariables[1:41,]
  rownames(groupingVariables) <- groupingVariables[,1]
  factors <- colnames(groupingVariables)
  groupingVariables[42, 2] <- "DipPop"

  # # get sample IDs from vcf
  # sampleIDs <- colnames(vcfRaw@gt)[-1]
  # 
  datasetScaled <- scaleGen(vcfGenind, NA.method="mean")
  # datasetScaled <- datasetScaled[1:41,]
  # 
  # 
  # ## MDS
  mds <- read.csv("data/plink_3101.mds", sep="")
  # 
  # ## t-SNE
  distanceMatrixTSNE <- read_tsv("data/plink_3101.dist", col_names = F)
  
  ## UMAP
  
  
  # cat(class(distanceMatrixTSNE))


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
  
  
  # If there's only one factor, duplicate it so everything that expects a 
  # second factor doesn't break: 
  
  # factorLevels <- NULL
  # for (i in seq_along(datasetScaled[factors])){
  #   factorLevels[[i]] <- paste0(
  #     colnames(datasetScaled[factors[[i]]]),
  #     ": ",
  #     levels(as.factor(datasetScaled[, factors[i]])))
  # }
  # factorLevels <- unlist(factorLevels)
  # 
  # 
  if (length(factors) == 1) {
    factors <- as.numeric(c(factors, factors))
  }

  colourList <- list()
  for (i in factors) {
    for (l in levels(as.factor(groupingVariables[, i]))) {
      colourList[l] <- NA
    }
  }

  return(list(
    # "pca" = pca,
    # "n_pcs" = n_pcs,
    # "n_grouping" = n_grouping,
    # "pc_list" = pc_list,
    # "pca_varprop" = pca_varprop,
    # "pca_tab" = pca_tab,
    # "tab_varprop" = tab_varprop,
    "vcfRaw" = vcfRaw,
    "vcfGenind" = vcfGenind,
    "datasetScaled" = datasetScaled,
    "groupingVariables" = groupingVariables,
    "distanceMatrixTSNE" = distanceMatrixTSNE,
    "colourList" = colourList
    # "factorLevels" = factorLevels
    # "mds" = mds,
    # "distanceMatrix" = distanceMatrix
    
    )
  )
  
# }
})
  