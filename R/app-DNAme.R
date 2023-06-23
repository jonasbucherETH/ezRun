###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodDNAme <- function(input = NA, output = NA, param = NA,
                           htmlFile = "00index.html") {
  
  ##### ----- libraries
  library("edgeR")
  library("dmrseq")
  library("tidyverse")
  library("rGREAT")
  library("KEGGREST")
  library("biomaRt")
  library("BioMartGOGeneSets")
  library("parallel")
  library("AnnotationHub")
  library("GenomicFeatures")
  library("biomartr") # for getGO
  library("marge")
  library(monaLisa)
  library(GenomicRanges)
  library(SummarizedExperiment)
  library(JASPAR2020)
  library(TFBSTools)
  # library(BSgenome.Mmusculus.UCSC.mm10)
  library(ComplexHeatmap)
  library(circlize)
  # from app-HOMER
  require(rtracklayer)
  require(ChIPpeakAnno)
  library("methylKit")
  library("genomation") # annotation of DML/R
  library(BiocParallel)
  
  ##### ----- input & output paths
  # #setwdNew(basename(output$getColumn("Report")))
  dataset <- input$meta
  ans4Report <- list() # a list of results for rmarkdown report
  ans4Report[["dataset"]] <- dataset
  
  report_dir <- basename(output$getColumn("Report"))
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add = TRUE)
  
  ###################### actual app ######################
  ### General preparation
  # sampleIDs <- input$getColumn("Name")
  sampleIDs <- c("a","b","c","d","e","f","g")
  
  if (param$cores > 1){
    BPPARAM <- MulticoreParam(workers = param$cores)
  } else {
    BPPARAM <- SerialParam()
  }
  
  
  ### DMR
  # TODO: check if COV contain all contexts
  mkdirDMR = paste("mkdir dmr")
  ezSystem(mkdirDMR)
  
  # testCovariateData <- input$getColumn(param$grouping)
  # bsseqColData <- as.data.frame(input$getColumn(param$grouping), col.names = param$grouping, row.names = input$getColumn("Name"))
  bsseqColData <- as.data.frame(input$getColumn(param$grouping), row.names = sampleIDs)
  colnames(bsseqColData) <- param$grouping
  covFilesBismark <- input$getFullPaths("COV")
  
  bsseq <- bsseq::read.bismark(files = covFilesBismark,
                               rmZeroCov = FALSE,
                               strandCollapse = FALSE,
                               verbose = FALSE,
                               colData = bsseqColData)
  
  ### test
  # sampleIDs <- c("a","b","c","d","e","f","g")
  # bsseqColData <- as.data.frame(c("40","40","40","40","60","60","60"), row.names = sampleIDs)
  # colnames(bsseqColData) <- "Treatment"
  # 
  # covFilesBismark <- list.files("/srv/gstore/projects/p1535/Bismark_JBmm_test3_2023-03-27--15-58-43/", pattern = "cov", full.names = T)
  # 
  # bsseq <- bsseq::read.bismark(files = covFilesBismark,
  #                              rmZeroCov = FALSE,
  #                              strandCollapse = FALSE,
  #                              verbose = FALSE,
  #                              colData = bsseqColData)
  # 
  # lociCoverage <- which(DelayedMatrixStats::rowSums2(getCoverage(bsseq, type="Cov")==0) == 0)
  # bsseqFiltered <- bsseq[lociCoverage, ]
  # 
  # dmRegions <- dmrseq(
  #   bs = bsseqFiltered,
  #   testCovariate = "Treatment", 
  #   # A continuous covariate is assumed if the data type in the 'testCovariate' slot is continuous,
  #   # with the exception of if there are only two unique values (then a two group comparison is carried out)
  #   # adjustCovariate = param$adjustCovariate,
  #   verbose = T # keep this
  # 
  # )
  ### test end
  
  lociCoverage <- which(DelayedMatrixStats::rowSums2(getCoverage(bsseq, type="Cov")==0) == 0)
  
  bsseqFiltered <- bsseq[lociCoverage, ]
  
  dmRegions <- dmrseq(
    bs = bsseqFiltered,
    testCovariate = param$grouping, 
    # A continuous covariate is assumed if the data type in the 'testCovariate' slot is continuous,
    # with the exception of if there are only two unique values (then a two group comparison is carried out)
    # adjustCovariate = param$adjustCovariate,
    verbose = FALSE, # keep this
    BPPARAM = BPPARAM
  )
  
  significantRegions <- dmRegions[dmRegions$qval < 0.10, ]
  
  ### setwd before saving results
  setwd("dmr")

  saveRDS(bsseq, file="bsseq.rds")
  saveRDS(bsseqFiltered, file="bsseqFiltered.rds")
  # saveRDS(bsseqColData, file="bsseqColData.rds")
  saveRDS(dmRegions, file="dmRegions.rds")
  saveRDS(significantRegions, file="significantRegions.rds")
  saveRDS(bsseqColData, file="bsseqColData.rds")
  
  ########################### DML ###########################
  # TODO: what pipeline to use for methRead?
  # mkdirDML = paste("mkdir dml")
  # ezSystem(mkdirDML)
  # 
  # treatmentMethylKit <- rep(0, length(sampleIDs))
  # treatmentMethylKit[input$getColumn(param$grouping) == param$sampleGroup, ] <- 1
  # 
  # # TODO: ask Deepak for mincov value (probably = 0, because filtering happening afterwards)
  # methylRawCpG <- methRead(input$getColumn("CpG_Context"),
  #                          sample.id=sampleIDs,
  #                          treatment=treatmentMethylKit, # 0 = control, 1 = test
  #                          context="CpG",
  #                          mincov = 10
  # )
  # 
  # methylRawCHG <- methRead(input$getColumn("CHG_Context"),
  #                          sample.id=sampleIDs,
  #                          treatment=treatmentMethylKit, # 0 = control, 1 = test
  #                          context="CHG",
  #                          mincov = 10
  # )
  # 
  # methylRawCHH <- methRead(input$getColumn("CHH_Context"),
  #                          sample.id=sampleIDs,
  #                          treatment=treatmentMethylKit, # 0 = control, 1 = test
  #                          context="CHH",
  #                          mincov = 10
  # )
  # 
  # for(i in seq_along(sampleIDs)){
  #   methylationStats <- getMethylationStats(methylRawCpG[[i]],plot=FALSE,both.strands=FALSE)
  #   assign(x = paste0("methylationStats_CpG_", sampleIDs[i]), value = methylationStats)
  #   methylationStats <- getMethylationStats(methylRawCHG[[i]],plot=FALSE,both.strands=FALSE)
  #   assign(x = paste0("methylationStats_CHG_", sampleIDs[i]), value = methylationStats)
  #   methylationStats <- getMethylationStats(methylRawCHH[[i]],plot=FALSE,both.strands=FALSE)
  #   assign(x = paste0("methylationStats_CHH_", sampleIDs[i]), value = methylationStats)
  # 
  #   coverageStats <- getCoverageStats(methylRawCpG[[i]],plot=FALSE,both.strands=FALSE)
  #   assign(x = paste0("coverageStats_CpG_", sampleIDs[i]), value = coverageStats)
  #   coverageStats <- getCoverageStats(methylRawCHG[[i]],plot=FALSE,both.strands=FALSE)
  #   assign(x = paste0("coverageStats_CHG_", sampleIDs[i]), value = coverageStats)
  #   coverageStats <- getCoverageStats(methylRawCHH[[i]],plot=FALSE,both.strands=FALSE)
  #   assign(x = paste0("coverageStats_CHH_", sampleIDs[i]), value = coverageStats)
  # }
  # 
  # 
  # 
  # 
  # # TODO: adapt values (ask Deepak)
  # filteredMethylRawCpG  <- filterByCoverage(
  #   methylRawCpG,
  #   lo.count=10, # Bases/regions having lower coverage than this count is discarded
  #   lo.perc=NULL, # Bases/regions having lower coverage than this percentile is discarded
  #   hi.count=NULL, # might want to filter out very high coverages as well (PCR bias)
  #   hi.perc=99.9 # Bases/regions having higher coverage than this percentile is discarded
  # )
  # 
  # filteredMethylRawCHG  <- filterByCoverage(
  #   methylRawCHG,
  #   lo.count=10, # Bases/regions having lower coverage than this count is discarded
  #   lo.perc=NULL, # Bases/regions having lower coverage than this percentile is discarded
  #   hi.count=NULL, # might want to filter out very high coverages as well (PCR bias)
  #   hi.perc=99.9 # Bases/regions having higher coverage than this percentile is discarded
  # )
  # 
  # filteredMethylRawCHH  <- filterByCoverage(
  #   methylRawCHH,
  #   lo.count=10, # Bases/regions having lower coverage than this count is discarded
  #   lo.perc=NULL, # Bases/regions having lower coverage than this percentile is discarded
  #   hi.count=NULL, # might want to filter out very high coverages as well (PCR bias)
  #   hi.perc=99.9 # Bases/regions having higher coverage than this percentile is discarded
  # )
  # 
  # # TODO: normalizeCoverage (?)
  # # TODO: 3.1 and 3.2 (?)
  # 
  # methylBase_CpG <- unite(filteredMethylRawCpG, destrand=TRUE) # destrand = T only for CpG
  # methylBase_CHG <- unite(filteredMethylRawCHG, destrand=FALSE) # destrand = T only for CpG
  # methylBase_CHH <- unite(filteredMethylRawCHH, destrand=FALSE) # destrand = T only for CpG
  # 
  # diffMeth_CpG <- calculateDiffMeth(
  #   methylBase_CpG,
  #   # covariates = bsseqColData, # data.frame, to separate from treatment effect
  #   mc.cores = param$cores
  #   )
  # 
  # diffMeth_CHG <- calculateDiffMeth(
  #   methylBase_CHG,
  #   # covariates = bsseqColData, # data.frame
  #   mc.cores = param$cores
  # )
  # 
  # diffMeth_CHH <- calculateDiffMeth(
  #   methylBase_CHH,
  #   # covariates = bsseqColData, # data.frame
  #   mc.cores = param$cores
  # )
  # 
  # # TODO: ask deepak whether to split into hypo-/hyper-methylated bases
  # diff25p_hyper_CpG <- getMethylDiff(diffMeth_CpG, difference=25, qvalue=0.01, type="hyper")
  # diff25p_hypo_CpG <- getMethylDiff(diffMeth_CpG, difference=25, qvalue=0.01, type="hypo")
  # diff25p_both_CpG <- getMethylDiff(diffMeth_CpG, difference=25, qvalue=0.01, type="all")
  # 
  # diff25p_hyper_CHG <- getMethylDiff(diffMeth_CHG, difference=25, qvalue=0.01, type="hyper")
  # diff25p_hypo_CHG <- getMethylDiff(diffMeth_CHG, difference=25, qvalue=0.01, type="hypo")
  # diff25p_both_CHG <- getMethylDiff(diffMeth_CHG, difference=25, qvalue=0.01, type="all")
  # 
  # diff25p_hyper_CHH <- getMethylDiff(diffMeth_CHH, difference=25, qvalue=0.01, type="hyper")
  # diff25p_hypo_CHH <- getMethylDiff(diffMeth_CHH, difference=25, qvalue=0.01, type="hypo")
  # diff25p_both_CHH <- getMethylDiff(diffMeth_CHH, difference=25, qvalue=0.01, type="all")
  # 
  # ### setwd before saving results
  # setwd("dml")
  # 
  # saveRDS(diffMeth_CpG, file=paste0(diffMeth_CpG, ".rds"))
  # saveRDS(diffMeth_CHG, file=paste0(diffMeth_CHG, ".rds"))
  # saveRDS(diffMeth_CHH, file=paste0(diffMeth_CHH, ".rds"))
  # saveRDS(diff25p_both_CpG, file=paste0(diff25p_both_CpG, ".rds"))
  # saveRDS(diff25p_both_CHG, file=paste0(diff25p_both_CHG, ".rds"))
  # saveRDS(diff25p_both_CHH, file=paste0(diff25p_both_CHH, ".rds"))
  # 
  # saveRDS(methylRawCpG, file=paste0(methylRawCpG, ".rds"))
  # saveRDS(methylRawCHG, file=paste0(methylRawCHG, ".rds"))
  # saveRDS(methylRawCHH, file=paste0(methylRawCHH, ".rds"))
  # 
  # transcriptFeatures <- readTranscriptFeatures(paste0(param$refBuild), "Genes/genes.bed") # bed file
  # # transcriptFeatures <- readTranscriptFeatures("/srv/GT/reference/Mus_musculus/UCSC/mm10/Annotation/Version-2012-05-23/Genes/genes.bed")
  # 
  # # TODO: think about if and how to do this (still split by hypo and hyper?)
  # diff25p_both_CpG_annotated <- annotateWithGeneParts(as(diff25p_both_CpG,"GRanges"), transcriptFeatures)
  # 
  # # read the shores and flanking regions and name the flanks as shores 
  # # and CpG islands as CpGi
  # cpg.obj=readFeatureFlank(system.file("extdata", "cpgi.hg18.bed.txt", 
  #                                      package = "methylKit"),
  #                          feature.flank.name=c("CpGi","shores"))
  # #
  # # convert methylDiff object to GRanges and annotate
  # diffCpGann=annotateWithFeatureFlank(as(myDiff25p,"GRanges"),
  #                                     cpg.obj$CpGi,cpg.obj$shores,
  #                                     feature.name="CpGi",flank.name="shores")
  # 
  # ######## HOMER Motif Analysis
  # # TODO: conversion to bed files (?)
  # setwd(report_dir)
  # 
  # lmr <- as(diffMethLoci,"GRanges")
  # lmr_25p <- as(diffMeth_25p,"GRanges")
  # 
  # indeces_25p <- which(lmr@ranges@start %in% lmr_25p@ranges@start)
  # 
  # # paste0("peakID_", seq(1:length(lmr)))
  # 
  # df <- data.frame(seqnames=seqnames(lmr),
  #                  starts=start(lmr)-1,
  #                  ends=end(lmr),
  #                  # names=c(rep(".", length(lmr))), # col 4 = unique Peak ID (?)
  #                  names=paste0("peakID_", seq(1:length(lmr))), # col 4 = unique Peak ID (?)
  #                  scores=c(rep(".", length(lmr))), # col 5 = not used
  #                  strands=strand(lmr))
  # 
  # df_25p <- df[indeces_25p, ]
  # 
  # mkdirHomer = paste("mkdir homer")
  # ezSystem(mkdirHomer)
  # setwd("./homer")
  # 
  # mkdirHomerDMR = paste("mkdir region")
  # ezSystem(mkdirHomerDMR)
  # 
  # mkdirHomerDML = paste("mkdir locus")
  # ezSystem(mkdirHomerDML)
  # 
  # region_size <- 200
  # motif_length <- "8,10,12"
  # # genome = fasta file
  # # genome <- "hg38" # fasta file
  # genomeHomer <- file.path(param$ref_selector, '../../Sequence/WholeGenomeFasta/genome.fa')
  # output_dir <- "homer"
  # ### needed: bed file DMR/L, bed file background x2
  # cmd <- paste("findMotifsGenome.pl", bedDMR, genomeHomer, "region", "-size", region_size, "-len", motif_length, "-bg", bedBGDMR)
  # ezSystem(cmd)
  # 
  # cmd <- paste("findMotifsGenome.pl", bedDML, genomeHomer, "locus", "-size", region_size, "-len", motif_length, "-bg", bedBGDML)
  # ezSystem(cmd)
  # 
  # ### Functional analysis (GREAT)
  # setwd(report_dir)
  # mkdirGreat = paste("mkdir great")
  # ezSystem(mkdirGreat)
  # 
  # greatResult_BP <- great(gr = significantRegions, gene_sets = "BP", biomart_dataset = param$biomart_selection,
  #                         background = dmRegions, cores = param$cores)
  # greatResult_CC <- great(gr = significantRegions, gene_sets = "CC", biomart_dataset = param$biomart_selection,
  #                         background = dmRegions, cores = param$cores)
  # greatResult_MF <- great(gr = significantRegions, gene_sets = "MF", biomart_dataset = param$biomart_selection,
  #                         background = dmRegions, cores = param$cores)
  # 
  # enrichmentTable_BP <- getEnrichmentTable(greatResult_BP)
  # enrichmentTable_CC <- getEnrichmentTable(greatResult_CC)
  # enrichmentTable_MF <- getEnrichmentTable(greatResult_MF)
  # 
  # if(param$biomart_selection=="athaliana_eg_gene") {
  #   reactome <- "https://plantreactome.gramene.org/download/current/gene_ids_by_pathway_and_species.tab"
  #   react <- data.frame(data.table::fread(input = reactome, header = F, nThread = 16))
  #   rdb <- react[grep(pattern = "^R-ATH", x = react$V1), ]
  #   reactome_pathways <- split(rdb$V4, paste(rdb$V1, rdb$V2, sep = ": "))
  #   
  #   ## KEGG pathways
  #   kg <- keggList("organism")
  #   
  #   pathway2gene <- keggLink("pathway", "ath")
  #   pathwayName <- keggList("pathway", "ath")
  #   df1 <- data.frame(
  #     gene = gsub("ath:", "", names(pathway2gene)),
  #     pathID = gsub("path:", "", pathway2gene)
  #   )
  #   
  #   df2 <- data.frame(
  #     pathID = gsub("path:", "", names(pathwayName)),
  #     name = pathwayName
  #   )
  #   
  #   df_kegg <- merge(df2, df1)
  #   
  #   kegg_pathways <- split(df_kegg$gene, paste(df_kegg$pathID, df_kegg$name,
  #                                              sep = ": "
  #   ))
  #   
  #   greatResult_RE <- great(gr = significantRegions, gene_sets = reactome_pathways, extended_tss = extendedTSS,
  #                           background = dmRegions, cores = cores)
  #   
  # 
  #   greatResult_RE <- great(gr = significantRegions, gene_sets = reactome_pathways, tss_source = "TxDb.Athaliana.BioMart.plantsmart51",
  #                           background = dmRegions, cores = param$cores)
  #   greatResult_KE <- great(gr = significantRegions, gene_sets = kegg_pathways, tss_source = "TxDb.Athaliana.BioMart.plantsmart51",
  #                           background = dmRegions, cores = param$cores)
  #   
  #   enrichmentTable_RE <- getEnrichmentTable(greatResult_RE)
  #   enrichmentTable_KE <- getEnrichmentTable(greatResult_KE)
  #   
  #   setwd("great")
  #   saveRDS(greatResult_RE, file = paste0("greatResult_RE", ".rds"))
  #   saveRDS(greatResult_KE, file = paste0("greatResult_KE", ".rds"))
  #   saveRDS(enrichmentTable_RE, file = paste0("enrichmentTable_RE", ".rds"))
  #   saveRDS(enrichmentTable_KE, file = paste0("enrichmentTable_KE", ".rds"))
  # }
  # 
  # setwd("great")
  # saveRDS(greatResult_BP, file = paste0("greatResult_BP", ".rds"))
  # saveRDS(greatResult_CC, file = paste0("greatResult_CC", ".rds"))
  # saveRDS(greatResult_MF, file = paste0("greatResult_MF", ".rds"))
  # saveRDS(enrichmentTable_BP, file = paste0("enrichmentTable_BP", ".rds"))
  # saveRDS(enrichmentTable_CC, file = paste0("enrichmentTable_CC", ".rds"))
  # saveRDS(enrichmentTable_MF, file = paste0("enrichmentTable_MF", ".rds"))
  

  ###################### testing ######################
  
  # filePathProject <- file.path("/srv/gstore/projects")
  # setwd(filePathProject)
  # 
  # input <- "p1535/Bismark_JB_test1_2023-03-07--15-47-45" # for testing
  # dataDir <- "/srv/gstore/projects/p1535/Bismark_JB_test1_2023-03-07--15-47-45"
  # 
  # ### input dataset
  # inputDataset <- read_tsv(file.path(dataDir, "input_dataset.tsv"))
  # dataset <- read_tsv(file.path(dataDir, "dataset.tsv"))
  # input <- dataset
  # sampleNames <- dataset$Name
  
  # ### M-bias
  # filesBias <- file.path(filePathProject, input$getColumn("M-bias"))
  # filesBias <- list.files(file.path(filePathProject, input), pattern = "M-bias")
  # ### CpG context
  # # filesCpgContext <- file.path(filePathProject, input$getColumn("cCpG_Context"))
  # filesCpgContext <- list.files(file.path(filePathProject, input), pattern = "CpG_context")
  # ### report
  # # filesReport <- file.path(filePathProject, input$getColumn("TxtReport"))
  # filesReport <- list.files(file.path(filePathProject, input), pattern = "report.txt")
  # ### reference
  # filesReference <- file.path(filePathProject, input$getColumn("refBuild"))
  # filesReference <- list.files(file.path(filePathProject, input), pattern = "report.txt")
  
  
  
  
  # cov <- read.delim(file.path("/srv/gstore/projects", input$getColumn("Grouping File")))
  # rownames(groupingVariables) <- groupingVariables[,1]
  
  
  ##### read input data
  # coverage
  # coverageFiles <- file.path("/srv/gstore/projects", input$getColumn("Filtered VCF"))
  # coverage <- readBismark2DGE(input$getColumn("COV"), sample.names=input$getColumn("Name"))
  # coverage <- readBismark2DGE(input$`COV [File]`, sample.names=sampleNames)
  # coverageEdgeR <- ReadBismark2DGE(input$`COV [File]`, sample.names=sampleNames, data.table = TRUE)
  
  ### dmrseq read
  # We strongly recommend you use the 'genome wide cytosine report' (*.CX_report.txt.gz)
  # same as report? same as CpG_report?
  # coverageFileList <- list.files(dataDir, pattern = "*.report.txt", full.names = T)
  # coverageDmrseq <- read.bismark(files = list.files(dataDir, pattern = "*.report.txt", full.names = T),
  #                              rmZeroCov = F,
  #                              strandCollapse = FALSE,
  #                              verbose = TRUE,
  #                              BPPARAM = 4,
  #                              # colData = input)
  #                              colData = NULL)
  
  # coverageFileList <- list.files(dataDir, pattern = "*bismark.cov.gz", full.names = T)
  # coverageDmrseq <- bsseq::read.bismark(files = coverageFileList,
  #                                rmZeroCov = F,
  #                                strandCollapse = FALSE,
  #                                verbose = TRUE,
  #                                colData = input$Name)
  
  ### mouse data
  # dataDir <- "/srv/gstore/projects/p1535/Bismark_JBmm_test3_2023-03-27--15-58-43"
  # inputDataset <- read_tsv(file.path(dataDir, "input_dataset.tsv"))
  # dataset <- read_tsv(file.path(dataDir, "dataset.tsv"))
  # sampleNames <- dataset$Name
  # 
  # coverageFileList_mm <- list.files(dataDir, pattern = "*bismark.cov.gz", full.names = T)
  # coverageDmrseq_mm <- bsseq::read.bismark(files = coverageFileList_mm,
  #                                       rmZeroCov = F,
  #                                       strandCollapse = FALSE,
  #                                       verbose = TRUE,
  #                                       colData = inputDataset)
                                        # colData = as.data.frame(inputDataset))

  
                                 # colData = input)
                                 # colData = inputDataset$Name)
  # 
  # reportFileList <- list.files(dataDir, pattern = "*CpG_context.txt", full.names = T)
  # coverageDmrseq <- bsseq::read.bismark(files = reportFileList,
  #                                       rmZeroCov = F,
  #                                       strandCollapse = FALSE,
  #                                       verbose = TRUE,
  #                                       BPPARAM = 4)
  #                                       # colData = input)
  #                                       # colData = inputDataset$Name)
  
  # infile <- system.file("extdata/test_data.fastq_bismark.bismark.cov.gz",
  #                       package = 'bsseq')
  # bismarkBSseq <- read.bismark(files = infile,
  #                              rmZeroCov = TRUE,
  #                              strandCollapse = FALSE,
  #                              verbose = TRUE)
                                 
  ######################################################################################################################## 

  
  ######################################################################################################################## 
  ##### ----- CNVkit - https://cnvkit.readthedocs.io/en/stable/pipeline.html#target
  # require: species’ reference genome sequence, in FASTA format
  # optional: Gene annotation database, via RefSeq or Ensembl, in BED or “RefFlat” format
  
  # ### target: Prepare a BED file of baited regions for use with CNVkit.
  # # cnvkit.py target my_baits.bed --annotate refFlat.txt --split -o my_targets.bed
  # cmd <- paste("cnvkit.py target", input$getColumn("BedGraph")), "--annotate", input$getColumn("refBuild")), "--split -o", output_dir)
  # result <- ezSystem(cmd)
  # gc()
  
  ### batch: Prepare a BED file of baited regions for use with CNVkit.
  # --method wgs - for whole-genome sequencing

  # cmd <- paste("cnvkit.py batch", input$getColumn("BAM"),
  #              "--targets", input$getColumn("BedGraph"), 
  #              # "--annotate", input$getColumn("refBuild"), # optional; Format: UCSC, refFlat.txt or ensFlat.txt file (preferred), or BED, interval list, GFF, or similar.
  #              "--fasta", file.path(input$getColumn("refBuild"), "Genes/"), # required
  #              "--access", input$getColumn("BedGraph"), # difference to other bed graph?
  #              "--output-reference reference.cnn",
  #              "--output-dir", output_dir,
  #              # --diagram
  #              # --scatter
  #              )
  # result <- ezSystem(cmd)
  # gc()
  
  ######################################################################################################################## 
  ##### ----- dmrseq

  
  ## Copy the style files and templates
  styleFiles <- file.path(
    system.file("templates", package = "ezRun"),
    c(
      "fgcz.css", "DNAme.Rmd",
      "fgcz_header.html", "banner.png"
    )
  )
  file.copy(from = styleFiles, to = ".", overwrite = TRUE)
  
  ### generate the main reports
  # rmarkdown::render(
  #   input = "DNAme.Rmd", envir = new.env(),
  #   output_dir = ".", output_file = htmlFile, quiet = TRUE
  # )
  
  html_files <- c("00index.html",  "banner.png",  "fgcz.css",  "fgcz_header.html")
  
  return("Success")
}


##' @template app-template
##' @templateVar method ezMethodMinimal(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run
##' @section Functions:
##' \itemize{
##'   \item{\code{plotReadCountToLibConc(dataset, colname): }}
##'   {Plots \code{colname} from \code{dataset} against read counts in millions.}
##'   \item{\code{getQualityMatrix(inputFile): }}
##'   {Gets a quality count matrix from a fastq or gziped fastq.gz file with dimensions read quality and read length.}
##'   \item{\code{plotQualityMatrixAsHeatmap(qualMatrixList, isR2=FALSE, xScale=1, yScale=1): }}
##'   {Returns a png table of quality matrices interpreted as heatmaps.}
##'   \item{\code{plotQualityHeatmap(result, name=NULL, colorRange=c(0,sqrt(40)), colors=gray((1:256)/256), main=NULL, pngFileName=NULL, xScale=1, yScale=1): }}
##'   {Creates and returns the images used by \code{plotQualityMatrixAsHeatmap()}.}
##' }

EzAppDNAme <-
  setRefClass("EzAppDNAme",
              contains = "EzApp",
              methods = list(
                initialize = function() {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodDNAme
                  name <<- "EzAppDNAme"
                  appDefaults <<- rbind(perLibrary = ezFrame(Type = "logical", DefaultValue = TRUE, Description = "DNAme brabra"))
                }
              )
  )

