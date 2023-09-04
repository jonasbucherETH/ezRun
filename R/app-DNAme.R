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
  library(valr) # GRanges to bed
  library(GenomeInfoDb)
  
  ##### ----- input & output paths
  # #setwdNew(basename(output$getColumn("Report")))
  dataset <- input$meta
  ans4Report <- list() # a list of results for rmarkdown report
  ans4Report[["dataset"]] <- dataset
  
  # report_dir <- basename(output$getColumn("Report"))
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add = TRUE)
  
  ###################### actual app ######################
  ### General preparation
  # sampleNames <- input$getColumn("Name")
  # sampleNames <- c("a","b","c","d","e","f","g")
  
  if (param$cores > 1){
    BPPARAM <- MulticoreParam(workers = param$cores)
  } else {
    BPPARAM <- SerialParam()
  }
  
  # significantRegions <- readRDS("/srv/gstore/projects/p1535/DNAme_fun_mm_test6_60--over--40_2023-07-20--13-14-42/DNAme/CpG/significantRegions.rds")
  # writeBedFileRegions(regions = significantRegions, nameBed = "significantRegions")
  # findMotifsGenome.pl ~/git/sushi/significantRegions.bed mm10 homer -nlen 0 -noweight -preparsedDir .
  
  # coverageFiles <- input$getFullPaths("COV")
  # sampleNames <- names(coverageFiles)
  sampleNames <- param$samples
  cat(sampleNames)
  bsseqColData <- as.data.frame(input$getColumn(param$grouping), row.names = sampleNames)
  colnames(bsseqColData) <- param$grouping
  
  # ezSystem("mkdir region loci")
  
  if (param$allCytosineContexts) {
    contexts <- c("CpG", "CHG", "CHH")
  } else {
    contexts <- c("CpG")
  }
  
  for (i in seq_along(contexts)) {
    # using coverage files
    # if (contexts[i] == "CpG") {
    #   covColumnName <- "COV"
    # } else {
    #   covColumnName <- paste0("COV_", contexts[i])
    # } 
    # using genome-wide cytosine report
    if (contexts[i] == "CpG") {
      covColumnName <- "COV"
    } else {
      covColumnName <- paste0("COV_", contexts[i])
    }
    cytosineReportColumnName <- paste0(contexts[i], "_report")

    setwd(cwd)
    setwd(basename(output$getColumn("Report")))
    setwdNew(contexts[i])
    # ezSystem(paste("mkdir", contexts[i]))
    # ezSystem("mkdir regions loci")
    # ezSystem("mkdir regions/hyper regions/hypo loci/hyper loci/hypo")
    ezSystem("mkdir hyper hypo")
    ezSystem("mkdir hyper/regions hyper/loci hypo/regions hypo/loci")
    
    # testDir <- "/srv/gstore/projects/p1535/Bismark_awk_test3_2023-07-24--15-14-57"
    coverageFiles <- input$getFullPaths(covColumnName)
    cytosineReportFiles <- input$getFullPaths(cytosineReportColumnName)
    
    cat(coverageFiles)
    cat(cytosineReportFiles)
    
    
    # cat(coverageFiles)
    
    # dd <- "/srv/gstore/projects/p1535/Bismark_mm_CpG_2023-07-24--16-13-37"
    # coverageFiles <- list.files(testDir, pattern = "*CHH.gz.bismark.cov.gz", full.names = T)
    # cytosineReportFiles <- list.files(testDir, pattern = "*CG_report.txt.gz", full.names = T)
    
    # input_datset <- read_tsv("/srv/gstore/projects/p1535/DNAme_CXreport_test1_BB--over--AN_2023-07-24--15-31-12/input_dataset.tsv")
    # bsseqColData <- data.frame(input_datset$`Treatment [Factor]`)
    # rownames(bsseqColData) <- input_datset$Name
    # colnames(bsseqColData) <- "Treatment"
    # bsseq_6 <- bsseq::read.bismark(files = coverageFiles[1], loci = loci, verbose = T, strandCollapse = FALSE)
    # strandCollapse?
    # lociCoverage <- which(DelayedMatrixStats::rowSums2(getCoverage(bsseq_6, type="M")==0) == 0)
    # bsseqFiltered <- bsseq_6[lociCoverage, ]
    
    # round(colMeans(getCoverage(bsseq_6)), 1)
    # cat(sampleNames)

    # loci <- bsseq:::.readBismarkAsFWGRanges(unname(cytosineReportFiles[1]), rmZeroCov = F, strandCollapse = F)
    bsseq <- bsseq::read.bismark(files = unname(coverageFiles),
                                    rmZeroCov = FALSE,
                                    strandCollapse = FALSE,
                                    verbose = TRUE,
                                    colData = bsseqColData,
                                    BPPARAM = BPPARAM
                                    # loci = loci,
                                    # nThread = 3
                                 )
    
    # loci <- readRDS("/scratch/jonas/62355/DNAme/CpG/loci.rds")
    
    seqlevelsStyle(loci) <- "UCSC"
    seqlevelsStyle(bsseq) <- "UCSC"
    
    lociCoverage <- which(DelayedMatrixStats::rowSums2(getCoverage(bsseq, type="Cov")==0) == 0)
    bsseqFiltered <- bsseq[lociCoverage, ]
    # cat(table(strand(bsseqFiltered)))
    
    # cat(length(bsseqFiltered))
    
    dmRegions <- dmrseq(
      bs = bsseqFiltered,
      testCovariate = param$grouping,
      # testCovariate = "Treatment",
      verbose = TRUE,
      BPPARAM = BPPARAM,
      cutoff = param$cutoffRegions
    )

    # if (length(dmRegions) > 0) {
    #   seqlevelsStyle(dmRegions) <- "UCSC"
    # }

    significantRegions <- dmRegions[dmRegions$qval < param$qvalRegions, ]
    # TODO: split into hypo- and hyper-methylated regions
    # note that for a two-group comparison dmrseq uses alphabetical order of the covariate of interest
    if (sort(c(param$sampleGroup, param$refGroup))[1] == param$refGroup) {
      significantRegions_hyper <- significantRegions[significantRegions$stat > 0, ]
      significantRegions_hypo <- significantRegions[significantRegions$stat < 0, ]
    } else {
      significantRegions_hyper <- significantRegions[significantRegions$stat < 0, ]
      significantRegions_hypo <- significantRegions[significantRegions$stat > 0, ]
    }
    
    # saveRDS(dmRegions, file=paste0("dmRegions", ".rds"))
    # saveRDS(significantRegions, file=paste0("significantRegions", ".rds"))
    # saveRDS(significantRegions_hyper, file=file.path("hyper", "regions", paste0("significantRegions", ".rds")))
    # saveRDS(significantRegions_hypo, file=file.path("hypo", "regions", paste0("significantRegions", ".rds")))
    

    # keepStandardChromosomes
    # significantRegions <- dmRegions[1:(length(dmRegions)/3),]
    # qvalCutoff <- 0.5
    # significantRegions <- dmRegions[dmRegions$qval < qvalCutoff, ]
    # sampleNames <- input_datset$Name
    treatmentMethylKit <- rep(0, length(sampleNames))
    treatmentMethylKit[input$getColumn(param$grouping) == param$sampleGroup] <- 1
    # cat(treatmentMethylKit)
    # treatmentMethylKit[bsseqColData == "BB"] <- 1
    

    methylRaw <- methRead(location = as.list(coverageFiles),
                          sample.id = as.list(sampleNames),
                          treatment = treatmentMethylKit,
                          pipeline = "bismarkCoverage",
                          # pipeline = "bismarkCytosineReport",
                          assembly = param$biomart_selection,
                          context = contexts[i],
                          mincov = 0,
                          skip = 0
    )
    
    # overlap <- findOverlaps(bsseq@rowRanges, loci)
    # strand(methylGR[ovl@to]) <- strand(bsseq@rowRanges[ovl@from])
    
    # methylRaw <- methRead(location = as.list(cytosineReportFiles),
    #                       sample.id = as.list(sampleNames),
    #                       treatment = treatmentMethylKit,
    #                       # pipeline = "bismarkCoverage",
    #                       pipeline = "bismarkCytosineReport",
    #                       assembly = "ath",
    #                       context = "CHH",
    #                       mincov = 1,
    #                       skip = 0
    # )
    # 
    # methylBase <- methylKit::unite(methylRaw)
    # methylGR <- as(methylBase, "GRanges")
    # seqlevelsStyle(methylGR) <- "UCSC"
    # 
    # ovl <- findOverlaps(bsseq@rowRanges, methylGR)
    # strand(methylGR[ovl@to]) <- strand(bsseq@rowRanges[ovl@from])
    # 
    # which(bsseq@rowRanges@ranges@start == methylBase$start)
    # class(bsseq@rowRanges)
    
    filteredMethylRaw  <- filterByCoverage(
      methylRaw,
      # lo.count=param$minCoverageBases, # Bases/regions having lower coverage than this count is discarded
      lo.count=NULL, # Bases/regions having lower coverage than this count is discarded
      lo.perc=0.1, # Bases/regions having lower coverage than this percentile is discarded
      # lo.perc=NULL, # Bases/regions having lower coverage than this percentile is discarded
      hi.count=NULL, # might want to filter out very high coverages as well (PCR bias)
      hi.perc=99.9 # Bases/regions having higher coverage than this percentile is discarded
      # hi.perc=NULL # Bases/regions having higher coverage than this percentile is discarded
    )
    methylBase <- methylKit::unite(filteredMethylRaw, destrand=FALSE, min.per.group = NULL, mc.cores = param$cores) # destrand = T only for CpG
    dmLoci <- calculateDiffMeth(methylBase, mc.cores = param$cores)
    significantLoci <- getMethylDiff(dmLoci, difference=param$cutoffLoci, qvalue=param$qvalLoci, type="all")
    significantLoci <- as(significantLoci,"GRanges")
    seqlevelsStyle(significantLoci) <- "UCSC"
    dmLoci <- as(dmLoci,"GRanges")
    seqlevelsStyle(dmLoci) <- "UCSC"

    significantLoci_hyper <- significantLoci[significantLoci$meth.diff > 0, ]
    significantLoci_hypo <- significantLoci[significantLoci$meth.diff < 0, ]

    saveRDS(bsseq, file=paste0("bsseq", ".rds"))
    # writeBedFileRegions(regions = dmRegions, nameBed = "regions/dmRegions")
    # saveRDS(dmRegions, file=paste0("regions/dmRegions", ".rds"))

    # methType <- "hyper"
    # dmType <- "regions"

    region_size <- 200
    motif_length <- "8,10,12"
    genomeHomer <- file.path("/srv/GT/reference", dirname(dirname(param$refBuild)), 'Sequence/WholeGenomeFasta/genome.fa')
    # genomeHomer <- file.path("/srv/GT/reference", "abc", 'Sequence/WholeGenomeFasta/genome.fa')

    # greatHomer(significantLoci_hyper, significantLoci_hypo, dmLoci, "loci")
    # greatHomer(significantRegions_hyper, significantRegions_hypo, dmRegions, "regions")
  }
  
####################################### --- FUNCTIONS START --- #######################################
  
  writeBedFileRegions <- function(regions, nameBed) {
    dfGR <- data.frame(chr=seqnames(regions),
                       starts=start(regions),
                       ends=end(regions),
                       names=paste0("peakID_", seq(1:length(regions))), # col 4 = unique Peak ID (?)
                       scores=c(rep(".", length(regions))), # col 5 = not used
                       strands=strand(regions))
    # strands=c(rep("+", length(regions))))
    write.table(dfGR, paste0(nameBed, ".bed"), sep = "\t", col.names = F, row.names = F, quote = F)
  }
  
  greatFun <- function(dmRegions, significantRegions, type) { # dmType = region / locus
    greatResult_BP <- great(gr = significantRegions, gene_sets = "BP", biomart_dataset = param$biomart_selection,
                            background = dmRegions, cores = param$cores)
    greatResult_CC <- great(gr = significantRegions, gene_sets = "CC", biomart_dataset = param$biomart_selection,
                            background = dmRegions, cores = param$cores)
    greatResult_MF <- great(gr = significantRegions, gene_sets = "MF", biomart_dataset = param$biomart_selection,
                            background = dmRegions, cores = param$cores)
    
    # enrichmentTable_BP <- getEnrichmentTable(greatResult_BP)
    # enrichmentTable_CC <- getEnrichmentTable(greatResult_CC)
    # enrichmentTable_MF <- getEnrichmentTable(greatResult_MF)
    
    if(param$biomart_selection=="athaliana_eg_gene") {
      reactome <- "https://plantreactome.gramene.org/download/current/gene_ids_by_pathway_and_species.tab"
      react <- data.frame(data.table::fread(input = reactome, header = F, nThread = 16))
      rdb <- react[grep(pattern = "^R-ATH", x = react$V1), ]
      reactome_pathways <- split(rdb$V4, paste(rdb$V1, rdb$V2, sep = ": "))
      
      ## KEGG pathways
      kg <- keggList("organism")
      
      pathway2gene <- keggLink("pathway", "ath")
      pathwayName <- keggList("pathway", "ath")
      df1 <- data.frame(
        gene = gsub("ath:", "", names(pathway2gene)),
        pathID = gsub("path:", "", pathway2gene)
      )
      
      df2 <- data.frame(
        pathID = gsub("path:", "", names(pathwayName)),
        name = pathwayName
      )
      
      df_kegg <- merge(df2, df1)
      kegg_pathways <- split(df_kegg$gene, paste(df_kegg$pathID, df_kegg$name,
                                                 sep = ": "
      ))
      
      greatResult_RE <- great(gr = significantRegions, gene_sets = reactome_pathways, tss_source = "TxDb.Athaliana.BioMart.plantsmart51",
                              background = dmRegions, cores = param$cores)
      greatResult_KE <- great(gr = significantRegions, gene_sets = kegg_pathways, tss_source = "TxDb.Athaliana.BioMart.plantsmart51",
                              background = dmRegions, cores = param$cores)
      
      # enrichmentTable_RE <- getEnrichmentTable(greatResult_RE)
      # enrichmentTable_KE <- getEnrichmentTable(greatResult_KE)
      saveRDS(greatResult_RE, file = file.path(type, paste0("greatResultRE", ".rds")))
      saveRDS(greatResult_KE, file = file.path(type, paste0("greatResultKE", ".rds")))
    }
    saveRDS(greatResult_BP, file = file.path(type, paste0("greatResultBP", ".rds")))
    saveRDS(greatResult_CC, file = file.path(type, paste0("greatResultCC", ".rds")))
    saveRDS(greatResult_MF, file = file.path(type, paste0("greatResultMF", ".rds")))
  }
  
  greatHomer <- function(significantHyper, significantHypo, differentialSet, dmType) { # methType = hypo/hyper | dmType = regions/loci
    if (length(differentialSet) > 0) {
      writeBedFileRegions(regions = differentialSet, nameBed = dmType)
      saveRDS(differentialSet, file=paste0(dmType, ".rds"))
      if (length(significantHyper) > 0 && length(differentialSet) > length(significantHyper)) {
        methType <- "hyper"
        writeBedFileRegions(regions = significantHyper, nameBed = file.path(methType, dmType, "significant"))
        saveRDS(significantHyper, file=file.path(methType, dmType, paste0("significant", ".rds")))
        greatFun(dmRegions = differentialSet, significantRegions = significantHyper, type = file.path(methType, dmType))
        # cmd <- paste("findMotifsGenome.pl", file.path(methType, dmType, "significant.bed"), genomeHomer, file.path(methType, dmType, "homer"), "-size 200", "-bg", file.path(methType, dmType, "full.bed"), "-len", motif_length, "-keepOverlappingBg", "-preparsedDir", file.path(methType, dmType))
        # ezSystem(cmd)
      } else {
        cat(paste0("no significant hypermethylated ", dmType, " found"))
      }
      if (length(significantHypo) > 0 && length(differentialSet) > length(significantHypo)) {
        methType <- "hypo"
        writeBedFileRegions(regions = significantHypo, nameBed = file.path(methType, dmType, "significant"))
        saveRDS(significantHypo, file=file.path(methType, dmType, paste0("significant", ".rds")))
        greatFun(dmRegions = differentialSet, significantRegions = significantHypo, type = file.path(methType, dmType))
        # cmd <- paste("findMotifsGenome.pl", file.path(methType, dmType, "significant.bed"), genomeHomer, file.path(methType, dmType, "homer"), "-size 200", "-bg", paste0(dmType, ".bed"), "-len", motif_length, "-keepOverlappingBg", "-preparsedDir .")
        # ezSystem(cmd)
      } else {
        cat(paste0("no significant hypomethylated ", dmType, " found"))
      } 
    } else { # remove subdirectory if no regions/loci for given methType and dmType combination
      cat(paste0("no methylated ", dmType, " found"))
      cmd <- paste("rm -r ", file.path(methType, dmType))
      ezSystem(cmd)
    }
  }
####################################### --- FUNCTIONS END --- #######################################
    
    # writeBedFileRegions(regions = dmRegions, nameBed = "regions/dmRegions")
    # saveRDS(dmRegions, file=paste0("regions/dmRegions", ".rds"))
    # writeBedFileRegions(regions = significantRegions, nameBed = "regions/significantRegions")
    # saveRDS(significantRegions, file=paste0("regions/significantRegions", ".rds"))
    
    # writeBedFileRegions(regions = significantRegions_hyper, nameBed = "regions/hyper/significantRegions")
    # writeBedFileRegions(regions = significantRegions_hypo, nameBed = "regions/hypo/significantRegions")
    # saveRDS(significantRegions_hyper, file=paste0("regions/hyper/significantRegions", ".rds"))
    # saveRDS(significantRegions_hypo, file=paste0("regions/hypo/significantRegions", ".rds"))
    # saveRDS(bsseq, file=paste0("regions/bsseq", ".rds"))

    
    # saveRDS(dmLoci, file=file.path("loci", paste0("dmLoci", ".rds")))
    # saveRDS(significantLoci, file=file.path("loci", paste0("significantLoci", ".rds")))
    # saveRDS(significantLoci_hyper, file=file.path("loci", "hyper", paste0("significantLoci", ".rds")))
    # saveRDS(significantLoci_hypo, file=file.path("loci", "hypo", paste0("significantLoci", ".rds")))
    # writeBedFileRegions(regions = dmLoci, nameBed = "loci/dmLoci")
    # writeBedFileRegions(regions = significantLoci, nameBed = "loci/significantLoci")
    # writeBedFileRegions(regions = significantLoci_hyper, nameBed = "loci/hyper/significantLoci")
    # writeBedFileRegions(regions = significantLoci_hypo, nameBed = "loci/hypo/significantLoci")

    # greatFun(dmRegions = dmRegions, significantRegions = significantRegions_hyper, type = "regions/hyper")
    # greatFun(dmRegions = dmRegions, significantRegions = significantRegions_hypo, type = "regions/hypo")
    # greatFun(dmRegions = dmLoci, significantRegions = significantLoci_hyper, type = "loci/hyper")
    # greatFun(dmRegions = dmLoci, significantRegions = significantLoci_hypo, type = "loci/hypo")
    
    # region_size <- 200
    # motif_length <- "8,10,12"
    # genomeHomer <- file.path("/srv/GT/reference", dirname(dirname(param$refBuild)), 'Sequence/WholeGenomeFasta/genome.fa')
    # # genomeHomer <- "mm10"
    # bedFile <- "significantRegions.bed"
    # bedFileBG <- "dmRegions.bed"
    # # cmd <- paste("findMotifsGenome.pl", bedFile, genomeHomer, "homer", "-size 200", "-bg", bedFileBG, "-len", motif_length, "-keepOverlappingBg", "-preparsedDir .")
    # # ezSystem(cmd)
    # cmd <- paste("findMotifsGenome.pl", "regions/hyper/significantRegions.bed", genomeHomer, "regions/hyper/homer", "-size 200", "-bg", "regions/dmRegions.bed", "-len", motif_length, "-keepOverlappingBg", "-preparsedDir regions/hyper")
    # ezSystem(cmd)
    # cmd <- paste("findMotifsGenome.pl", "regions/hypo/significantRegions.bed", genomeHomer, "regions/hypo/homer", "-size 200", "-bg", "regions/dmRegions.bed", "-len", motif_length, "-keepOverlappingBg", "-preparsedDir regions/hypo")
    # ezSystem(cmd)
    # 
    # cmd <- paste("findMotifsGenome.pl", "loci/hyper/significantLoci.bed", genomeHomer, "loci/hyper/homer", "-size 200", "-bg", "loci/dmLoci.bed", "-len", motif_length, "-keepOverlappingBg", "-preparsedDir loci/hyper")
    # ezSystem(cmd)
    # cmd <- paste("findMotifsGenome.pl", "loci/hypo/significantLoci.bed", genomeHomer, "loci/hypo/homer", "-size 200", "-bg", "loci/dmLoci.bed", "-len", motif_length, "-keepOverlappingBg", "-preparsedDir loci/hypo")
    # ezSystem(cmd)

  # print(sampleNames)
  
  # coverageFilesCHG <- input$getFullPaths("COV_CHG")
  # coverageFilesCHG <- unname(coverageFilesCHG)
  # coverageFilesCHH <- input$getFullPaths("COV_CHH")
  # coverageFilesCHH <- unname(coverageFilesCHH)
  
  # testDir <- "/srv/gstore/projects/p1535/Bismark_ath_full8_2023-07-12--17-27-56"
  # testDir <- "/srv/gstore/projects/p1535/Bismark_mm_full3_2023-07-13--11-55-16"
  # CpG_report_files <- list.files(testDir, pattern = "CG_report.txt", full.names = T)
  # CHG_report_files <- list.files(testDir, pattern = "CHG_report.txt", full.names = T)
  # CHH_report_files <- list.files(testDir, pattern = "CHH_report.txt", full.names = T)
  # bsseqColData <- as.data.frame(c("BB","BB","AN","AN","BB","AN","AN"), row.names = c('A04_BB','A06_BB','A07_AN','A08_AN','A09_BB','A10_AN','A15_AN'))
  # bsseqColData <- as.data.frame(c(40,40,50,50,60,60,60), row.names = c('SRR4105496','SRR4105497','SRR4105498','SRR4105499','SRR4105500','SRR4105501','SRR4105502'))
  # colnames(bsseqColData) <- "Treatment"
  
  # CpG_report_files <- input$getFullPaths("CpG_report")
  # CHG_report_files <- input$getFullPaths("CHG_report")
  # CHH_report_files <- input$getFullPaths("CHH_report")
  
  # coverageFilesAll <- list.files(testDir, pattern = ".gz.bismark.cov.gz", full.names = T)
  # coverageFilesCHG <- list.files(testDir, pattern = "CHG.gz.bismark.cov.gz", full.names = T)
  # coverageFilesCHH <- list.files(testDir, pattern = "CHH.gz.bismark.cov.gz", full.names = T)
  # coverageFilesCpG <- coverageFilesAll[!(coverageFilesAll %in% coverageFilesCHG)]
  # coverageFilesCpG <- coverageFilesCpG[!(coverageFilesCpG %in% coverageFilesCHH)]
  
  # saveRDS(coverageCpG, file="dmr/coverageCpG.rds")
  
  # coverageCpG <- readRDS("/scratch/Bismark_JBmm_test3_2023-03-27--15-58-43_temp18497/DNAme/dmr/coverageCpG.rds")
  
  
  # coverageCpG <- as.list(coverageCpG)
  # print(class(coverageCpG))
  # print(coverageCpG)
  # print(coverageCpG[[1]])
  
  # coverageCpG <- list.files("/srv/gstore/projects/p1535/Bismark_JBmm_test3_2023-03-27--15-58-43", pattern = "cov", full.names = T)
  ## cat(class(coverageCpG))
  
 # cat("1")
  
  # bsseq <- bsseq::read.bismark(files = CX_reportFiles,
  #                                 rmZeroCov = FALSE,
  #                                 strandCollapse = FALSE,
  #                                 verbose = TRUE,
  #                                 colData = bsseqColData)
  
  # bsseqCpG <- bsseq::read.bismark(files = unname(coverageFilesCpG),
  #                              rmZeroCov = FALSE,
  #                              strandCollapse = FALSE,
  #                              verbose = FALSE,
  #                              colData = bsseqColData)
  
  # bsseqCHG <- bsseq::read.bismark(files = unname(coverageFilesCpG),
  #                                 rmZeroCov = FALSE,
  #                                 strandCollapse = FALSE,
  #                                 verbose = FALSE,
  #                                 colData = bsseqColData)
  # 
  # bsseqCHH <- bsseq::read.bismark(files = unname(coverageFilesCpG),
  #                                 rmZeroCov = FALSE,
  #                                 strandCollapse = FALSE,
  #                                 verbose = FALSE,
  #                                 colData = bsseqColData)
  
  
  # 
  # 
  # bsseqColData <- as.data.frame(c("BB","BB","BB","AN","AN","AN","AN"), row.names = c('A04_BB','A06_BB','A09_BB','A07_AN','A08_AN','A10_AN','A15_AN'))
  # colnames(bsseqColData) <- "Treatment"
  
  # bsseqCX <- bsseq::read.bismark(files = CX_reportFiles,
  #                                rmZeroCov = FALSE,
  #                                strandCollapse = FALSE,
  #                                verbose = TRUE,
  #                                colData = bsseqColData,
  #                                loci = NULL,
  #                                nThread = 2)
                                 # BACKEND = "HDF5Array" ,
                                 # dir = "/scratch/jonas/bsseq_temp")
    

 # cat("2")
  
  # pData(bsseq)$Treatment <- input$getColumn(param$grouping)
  
  # bsseq <- readRDS("/scratch/Bismark_JBmm_test3_2023-03-27--15-58-43_temp18497/DNAme/dmr/bsseq.rds")
  # dmRegions <- readRDS("/scratch/Bismark_JBmm_test3_2023-03-27--15-58-43_temp18497/DNAme/dmr/dmRegions.rds")
  # significantRegions <- readRDS("/scratch/Bismark_JBmm_test3_2023-03-27--15-58-43_temp18497/DNAme/dmr/significantRegions.rds")
  
  
  
  ### test
  # sampleNames <- c("a","b","c","d","e","f","g")
  # bsseqColData <- as.data.frame(c("40","40","40","40","60","60","60"), row.names = sampleNames)
  # colnames(bsseqColData) <- "Treatment"
  # 
  # coverageCpG <- list.files("/srv/gstore/projects/p1535/Bismark_JBmm_test3_2023-03-27--15-58-43/", pattern = "cov", full.names = T)
  # 
  # bsseq <- bsseq::read.bismark(files = coverageCpG,
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
  # lociCoverage <- which(DelayedMatrixStats::rowSums2(getCoverage(bsseq, type="Cov")==0) == 0)
  
  # lociCoverageCpG <- which(DelayedMatrixStats::rowSums2(getCoverage(bsseqCpG, type="Cov")==0) == 0)
  # lociCoverageCHG <- which(DelayedMatrixStats::rowSums2(getCoverage(bsseqCHG, type="Cov")==0) == 0)
  # lociCoverageCHH <- which(DelayedMatrixStats::rowSums2(getCoverage(bsseqCHH, type="Cov")==0) == 0)
  
 # cat("3")
  # bsseqFiltered <- bsseq[lociCoverage, ]
  
  # bsseqFilteredCpG <- bsseqCpG[lociCoverageCpG, ]
  # bsseqFilteredCHG <- bsseqCHG[lociCoverageCHG, ]
  # bsseqFilteredCHH <- bsseqCHH[lociCoverageCHH, ]
  
 # cat("3")
  
  # dmRegions <- dmrseq(
  #   bs = bsseqFilteredCpG,
  #   # bs = bsseqCpG,
  #   testCovariate = param$grouping,
  #   verbose = TRUE, # keep this
  #   BPPARAM = BPPARAM
  # )
  
  # dmRegionsCHG <- dmrseq(
  #   bs = bsseqFilteredCHG,
  #   testCovariate = param$grouping, 
  #   verbose = TRUE, # keep this
  #   BPPARAM = BPPARAM
  # )
  # 
  # dmRegionsCHH <- dmrseq(
  #   bs = bsseqFilteredCHH,
  #   testCovariate = param$grouping, 
  #   verbose = TRUE, # keep this
  #   BPPARAM = BPPARAM
  # )
  
  # dmRegionsCpG <- dmrseq(
  #   bs = bsseqFilteredCpG,
  #   testCovariate = "Treatment",
  #   verbose = TRUE,
  #   cutoff = 0.001
  # )
  # dmRegionsCHG <- dmrseq(
  #   bs = bsseqFilteredCHG,
  #   testCovariate = "Treatment",
  #   verbose = TRUE,
  #   cutoff = 0.001
  # )
  # dmRegionsCHH <- dmrseq(
  #   bs = bsseqFilteredCHH,
  #   testCovariate = "Treatment",
  #   verbose = TRUE,
  #   cutoff = 0.001
  # )
  # seqlevelsStyle(dmRegionsCpG) <- "UCSC"
  # dmRegionsCpG <- readRDS("/scratch/jonas/17679/DNAme/dmRegionsCpG.rds")
  
  # significantRegionsCpG <- dmRegionsCpG[1:(length(dmRegionsCpG)/3),]
  # cat("4")
  
  # qvalCutoff <- 0.5
  # significantRegionsCpG <- dmRegionsCpG[dmRegionsCpG$qval < qvalCutoff, ]
  # significantRegionsCHG <- dmRegionsCHG[dmRegionsCHG$qval < qvalCutoff, ]
  # significantRegionsCHH <- dmRegionsCHH[dmRegionsCHH$qval < qvalCutoff, ]
  
  # TODO: transform to bedGraph for HOMER
  ### setwd before saving results
  # setwd("dmr")
 # cat("5")
  # seqlevelsStyle(dmRegionsCpG) <- "UCSC"
  # seqlevelsStyle(significantRegionsCpG) <- "UCSC"
  # keepStandardChromosomes
  
  # writeBedFileRegions <- function(regions, nameBed) {
  #   dfGR <- data.frame(chr=seqnames(regions),
  #                      starts=start(regions),
  #                      ends=end(regions),
  #                      names=paste0("peakID_", seq(1:length(regions))), # col 4 = unique Peak ID (?)
  #                      scores=c(rep(".", length(regions))), # col 5 = not used
  #                      # strands=strand(regions))
  #                      strands=c(rep("+", length(regions))))
  #   write.table(dfGR, paste0(nameBed, ".bed"), sep = "\t", col.names = F, row.names = F)
  # }
  # writeBedFileRegions(regions = dmRegionsCpG, nameBed = "dmRegionsCpG")
  # writeBedFileRegions(regions = dmRegionsCHG, nameBed = "dmRegionsCHG")
  # writeBedFileRegions(regions = dmRegionsCHH, nameBed = "dmRegionsCHH")
  # writeBedFileRegions(regions = significantRegionsCpG, nameBed = "significantRegionsCpG")
  # writeBedFileRegions(regions = significantRegionsCHG, nameBed = "significantRegionsCHG")
  # writeBedFileRegions(regions = significantRegionsCHH, nameBed = "significantRegionsCHH")

  # saveRDS(bsseqCpG, file=paste0("bsseqCpG", ".rds"))
  # saveRDS(bsseqCHG, file=paste0(bsseqCHG, ".rds"))
  # saveRDS(bsseqCHH, file=paste0(bsseqCHH, ".rds"))
  # saveRDS(dmRegionsCpG, file=paste0("dmRegionsCpG", ".rds"))
  # saveRDS(dmRegionsCHG, file=paste0(dmRegionsCHG, ".rds"))
  # saveRDS(dmRegionsCHH, file=paste0(dmRegionsCHH, ".rds"))
  ### maybe don't save significantRegions, as it is merely a subset of dmRegions
  # saveRDS(significantRegionsCpG, file=paste0("significantRegionsCpG", ".rds"))
  # saveRDS(significantRegionsCHG, file=paste0(significantRegionsCHG, ".rds"))
  # saveRDS(significantRegionsCHH, file=paste0(significantRegionsCHH, ".rds"))
  # saveRDS(bsseq, file=paste0(bsseq, ".rds"))
  
  # awk -F "\t" ' if($6==CHH) { print $0}' /srv/gstore/projects/p1535/Bismark_ath_full7_2023-07-10--13-01-23/A04_BB.CX_report.txt.gz > ./CHH_report.txt.gz
  # awk -F "\t" '$6==CG' /srv/gstore/projects/p1535/Bismark_ath_full7_2023-07-10--13-01-23/A04_BB.CX_report.txt.gz > ./CG_report.txt.gz
  # awk -F "\t" '$6==CHG' /srv/gstore/projects/p1535/Bismark_ath_full7_2023-07-10--13-01-23/A04_BB.CX_report.txt.gz > ./CHG_report.txt.gz
  # awk -F "\t" '{ if($6 == CHH) { print $0} }' /srv/gstore/projects/p1535/Bismark_ath_full7_2023-07-10--13-01-23/A04_BB.CX_report.txt.gz > ./CHH_report.txt.gz
  # awk -F "\t" '{ if($6 == "CG") { print $0} }' /srv/gstore/projects/p1535/Bismark_ath_full7_2023-07-10--13-01-23/A04_BB.CX_report.txt.gz > ./CG_report.txt.gz
  # awk -F '\t' ' $6 == "CG" ' /srv/gstore/projects/p1535/Bismark_ath_full7_2023-07-10--13-01-23/A04_BB.CX_report.txt.gz > ./CG_report.txt.gz
  # scp jobucher@fgcz-c-047.uzh.ch:/srv/gstore/projects/p1535/Bismark_ath_full7_2023-07-10--13-01-23/A04_BB.CX_report.txt.gz .
  # scp jobucher@fgcz-c-047.uzh.ch:/srv/gstore/projects/p1535/Bismark_ath_full7_2023-07-10--13-01-23/A04_BB.CHG.gz.bismark.cov.gz .
  # awk -F '\t' '{ if($6 == CHH) { print } }' A04_BB.CX_report.txt > CG_report.txt
  # awk -F'\t' '{print > ($6 ".txt")}' A04_BB.CX_report.txt
  
  ######################################## DML ########################################
  # TODO: what pipeline to use for methRead?
  # cat("5")
  
  # setwdNew(basename(output$getColumn("Report")))
  # mkdirDML = paste("mkdir dml")
  # ezSystem(mkdirDML)
  

  # sampleNames <- c("SRR4105496","SRR4105497","SRR4105498","SRR4105499","SRR4105500","SRR4105501","SRR4105502")
  # treatmentMethylKit <- c(rep(0, 4), rep(1, 3))
                          
  # treatmentMethylKit <- rep(0, length(sampleNames))
  # treatmentMethylKit[input$getColumn(param$grouping) == param$sampleGroup] <- 1
  # treatmentMethylKit[bsseqColData$"Treatment" == 40] <- 1
  
  
  # cat("6")

  # TODO: ask Deepak for mincov value (probably = 0, because filtering happening afterwards)
  # TODO: check out readBismarkCoverage: https://gist.github.com/al2na/4839e615e2401d73fe51
  # methylKitFiles <- unname(input$getFullPaths("COV"))
  
  # methylKitFiles <- list.files("/srv/gstore/projects/p1535/Bismark_JBmm_test3_2023-03-27--15-58-43", pattern = "cov", full.names = T)
  
  
  # methylKitFiles <- as.list(c("/srv/gstore/projects/p1535/Bismark_JBmm_test3_2023-03-27--15-58-43/SRR4105496.CpG_context.txt",
  #                     "/srv/gstore/projects/p1535/Bismark_JBmm_test3_2023-03-27--15-58-43/SRR4105497.CpG_context.txt",
  #                     "/srv/gstore/projects/p1535/Bismark_JBmm_test3_2023-03-27--15-58-43/SRR4105498.CpG_context.txt",
  #                     "/srv/gstore/projects/p1535/Bismark_JBmm_test3_2023-03-27--15-58-43/SRR4105499.CpG_context.txt",
  #                     "/srv/gstore/projects/p1535/Bismark_JBmm_test3_2023-03-27--15-58-43/SRR4105500.CpG_context.txt",
  #                     "/srv/gstore/projects/p1535/Bismark_JBmm_test3_2023-03-27--15-58-43/SRR4105501.CpG_context.txt",
  #                     "/srv/gstore/projects/p1535/Bismark_JBmm_test3_2023-03-27--15-58-43/SRR4105502.CpG_context.txt"))
  
  # CX_reportFiles <- input$getFullPaths("CX_report")
  # CX_reportFiles <- unname(CX_reportFiles)
  # CX_reportFiles <- list.files(testDir, pattern = "*CX_report.txt.gz", full.names = T)
  # COV_Files <- list.files(testDir, pattern = "*cov*", full.names = T)
  
  # methylRawL <- methRead(location = as.list(CpG_report_files),
  #                          sample.id = as.list(rownames(bsseqColData)),
  #                          treatment = c(1,1,1,0,0,0,0), # 0 = control, 1 = test
  #                          pipeline = list(fraction=FALSE,chr.col=1,start.col=2,end.col=2,
  #                                          coverage.col=4,strand.col=3,freqC.col=5),
  #                          assembly = "Ath",
  #                          context= "CpG",
  #                          mincov = 0,
  #                          skip = 0
  # )
  
  # sampleNames <- rownames(bsseqColData)
  # treatmentMethylKit <- rep(0, length(rownames(bsseqColData)))
  # treatmentMethylKit[bsseqColData$Treatment == "BB"] <- 1

                          
  # methylRawCpG <- methRead(location = as.list(coverageFilesCpG),
  #                          sample.id = as.list(sampleNames),
  #                          treatment = treatmentMethylKit,
  #                          pipeline = "bismarkCoverage",
  #                          assembly = param$biomart_selection,
  #                          # assembly = "Ath",
  #                          context= "CpG",
  #                          mincov = 0,
  #                          skip = 0
  # )
  
  # methylRawCHG <- methRead(location = as.list(coverageFilesCHG),
  #                          sample.id = as.list(sampleNames),
  #                          treatment = c(1,1,1,0,0,0,0),
  #                          pipeline = "bismarkCoverage",
  #                          assembly = param$biomart_selection,
  #                          # assembly = "Ath",
  #                          context= "CHG",
  #                          mincov = 0,
  #                          skip = 0
  # )
  # 
  # methylRawCHH <- methRead(location = as.list(coverageFilesCHH),
  #                          sample.id = as.list(sampleNames),
  #                          treatment = c(1,1,1,0,0,0,0), # 0 = control, 1 = test
  #                          pipeline = "bismarkCoverage",
  #                          assembly = param$biomart_selection,
  #                          # assembly = "Ath",
  #                          context= "CHH",
  #                          mincov = 0,
  #                          skip = 0
  # )
  
  
  # 
  # for(i in seq_along(sampleNames)){
  #   methylationStats <- getMethylationStats(methylRawCpG[[i]],plot=FALSE,both.strands=FALSE)
  #   assign(x = paste0("methylationStats_CpG_", sampleNames[i]), value = methylationStats)
  #   methylationStats <- getMethylationStats(methylRawCHG[[i]],plot=FALSE,both.strands=FALSE)
  #   assign(x = paste0("methylationStats_CHG_", sampleNames[i]), value = methylationStats)
  #   methylationStats <- getMethylationStats(methylRawCHH[[i]],plot=FALSE,both.strands=FALSE)
  #   assign(x = paste0("methylationStats_CHH_", sampleNames[i]), value = methylationStats)
  # 
  #   coverageStats <- getCoverageStats(methylRawCpG[[i]],plot=FALSE,both.strands=FALSE)
  #   assign(x = paste0("coverageStats_CpG_", sampleNames[i]), value = coverageStats)
  #   coverageStats <- getCoverageStats(methylRawCHG[[i]],plot=FALSE,both.strands=FALSE)
  #   assign(x = paste0("coverageStats_CHG_", sampleNames[i]), value = coverageStats)
  #   coverageStats <- getCoverageStats(methylRawCHH[[i]],plot=FALSE,both.strands=FALSE)
  #   assign(x = paste0("coverageStats_CHH_", sampleNames[i]), value = coverageStats)
  # }
  # 
  # 
  # 
  # 
  # # TODO: adapt values (ask Deepak)
  # filteredMethylRawCpG  <- filterByCoverage(
  #   methylRawCpG,
  #   lo.count=NULL, # Bases/regions having lower coverage than this count is discarded
  #   lo.perc=NULL, # Bases/regions having lower coverage than this percentile is discarded
  #   hi.count=NULL, # might want to filter out very high coverages as well (PCR bias)
  #   hi.perc=99.9 # Bases/regions having higher coverage than this percentile is discarded
  # )
  
  # filteredMethylRawCHG  <- filterByCoverage(
  #   methylRawCHG,
  #   lo.count=1, # Bases/regions having lower coverage than this count is discarded
  #   lo.perc=NULL, # Bases/regions having lower coverage than this percentile is discarded
  #   hi.count=NULL, # might want to filter out very high coverages as well (PCR bias)
  #   hi.perc=99.9 # Bases/regions having higher coverage than this percentile is discarded
  # )
  # 
  # filteredMethylRawCHH  <- filterByCoverage(
  #   methylRawCHH,
  #   lo.count=1, # Bases/regions having lower coverage than this count is discarded
  #   lo.perc=NULL, # Bases/regions having lower coverage than this percentile is discarded
  #   hi.count=NULL, # might want to filter out very high coverages as well (PCR bias)
  #   hi.perc=99.9 # Bases/regions having higher coverage than this percentile is discarded
  # )
  
  ## cat("8")
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
  # methylBaseCpG <- unite(filteredMethylRawCpG, destrand=FALSE) # destrand = T only for CpG
  # methylBaseCHG <- unite(filteredMethylRawCHG, destrand=FALSE) # destrand = T only for CpG
  # methylBaseCHH <- unite(filteredMethylRawCHH, destrand=FALSE) # destrand = T only for CpG

  # 
  # dmLociCpG <- calculateDiffMeth(methylBaseCpG, mc.cores = param$cores)
  # dmLociCHG <- calculateDiffMeth(methylBaseCHG, mc.cores = param$cores)
  # dmLociCHH <- calculateDiffMeth(methylBaseCHH, mc.cores = param$cores)
  
  # dmLociCpG <- calculateDiffMeth(methylBaseCpG)
  # dmLociCHG <- calculateDiffMeth(methylBaseCHG)
  # dmLociCHH <- calculateDiffMeth(methylBaseCHH)
  
  
  ## cat("9")
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
  # significantLociCpG <- getMethylDiff(dmLociCpG, difference=25, qvalue=0.01, type="all")
  # significantLociCpG_hyper <- getMethylDiff(dmLociCpG, difference=25, qvalue=0.01, type="hyper")
  # significantLociCpG_hypo <- getMethylDiff(dmLociCpG, difference=25, qvalue=0.01, type="hypo")

  # significantLociCHG <- getMethylDiff(dmLociCHG, difference=25, qvalue=0.01, type="all")
  # significantLociCHG_hyper <- getMethylDiff(dmLociCHG, difference=25, qvalue=0.01, type="hyper")
  # significantLociCHG_hypo <- getMethylDiff(dmLociCHG, difference=25, qvalue=0.01, type="hypo")

  # significantLociCHH <- getMethylDiff(dmLociCHH, difference=25, qvalue=0.01, type="all")
  # significantLociCHH_hyper <- getMethylDiff(dmLociCHH, difference=25, qvalue=0.01, type="hyper")
  # significantLociCHH_hypo <- getMethylDiff(dmLociCHH, difference=25, qvalue=0.01, type="hypo")

  # setwd("homer/dml_test/")
  # getwd()
  
  ### methylKit objects to bedGraph file / data.frame
  # bedgraph(
  #   file.name = paste0("significantLociCpG", ".bed"),
  #   methylObj = significantLociCpG,
  #   col.name = 'qvalue' # 'pvalue','qvalue', 'meth.diff'
  # )
  # bedgraph(
  #   file.name = paste0("significantLociCHG", ".bed"),
  #   methylObj = significantLociCHG,
  #   col.name = 'qvalue' # 'pvalue','qvalue', 'meth.diff'
  # )
  # bedgraph(
  #   file.name = paste0("significantLociCHH", ".bed"),
  #   methylObj = significantLociCHH,
  #   col.name = 'qvalue' # 'pvalue','qvalue', 'meth.diff'
  # )
  # bedgraph(
  #   file.name = paste0("dmLociCpG", ".bed"),
  #   methylObj = dmLociCpG,
  #   col.name = 'qvalue' # 'pvalue','qvalue', 'meth.diff'
  # )
  # bedgraph(
  #   file.name = paste0("dmLociCHG", ".bed"),
  #   methylObj = dmLociCHG,
  #   col.name = 'qvalue' # 'pvalue','qvalue', 'meth.diff'
  # )
  # bedgraph(
  #   file.name = paste0("dmLociCHH", ".bed"),
  #   methylObj = dmLociCHH,
  #   col.name = 'qvalue' # 'pvalue','qvalue', 'meth.diff'
  # )
  
  # gr_cpg <- as(dmLociCpG,"GRanges")
  # export(gr_cpg, "dm.bed", "bed")
  # gr_cpg_sig <- as(significantLociCpG,"GRanges")
  # export(gr_cpg_sig, "sl.bed", "bed")
  # gr_cpg$peakID <- paste0("peakID_", seq(1:length(gr_cpg)))
  
  # indeces_25p <- which(lmr@ranges@start %in% lmr_25p@ranges@start)
  
  # paste0("peakID_", seq(1:length(lmr)))
  
  # df <- data.frame(chr=seqnames(gr_cpg),
  #                  starts=start(gr_cpg),
  #                  ends=end(gr_cpg),
  #                  # names=c(rep(".", length(lmr))), # col 4 = unique Peak ID (?)
  #                  names=paste0("peakID_", seq(1:length(gr_cpg))), # col 4 = unique Peak ID (?)
  #                  scores=c(rep(".", length(gr_cpg))), # col 5 = not used
  #                  strands=strand(gr_cpg))
  



  # ### setwd before saving results
  # setwd("dml")
  
  # significantLociCpG_bedGraph <- bedgraph(
  #   methylObj = significantLociCHG,
  #   col.name = 'qvalue' # 'pvalue','qvalue', 'meth.diff'
  # )# saveRDS(methylBaseCpG, file=paste0(methylBaseCpG, ".rds"))
  # saveRDS(methylBaseCHG, file=paste0(methylBaseCHG, ".rds"))
  # saveRDS(methylBaseCHH, file=paste0(methylBaseCHH, ".rds"))
  # saveRDS(dmLociCpG, file=paste0("dmLociCpG", ".rds"))
  # saveRDS(dmLociCHG, file=paste0(dmLociCHG, ".rds"))
  # saveRDS(dmLociCHH, file=paste0(dmLociCHH, ".rds"))
  # saveRDS(significantLociCpG, file=paste0("significantLociCpG", ".rds"))
  # saveRDS(significantLociCHG, file=paste0(diffMeth_CpG, ".rds"))
  # saveRDS(significantLociCHH, file=paste0(diffMeth_CpG, ".rds"))

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
  ##################### HOMER Motif Analysis ########################
  # TODO: conversion to bed files (?)
  # setwdNew(basename(output$getColumn("Report")))
  # 
  # mkdirHomer = paste("mkdir homer")
  # ezSystem(mkdirHomer)
  # setwd("homer")
  # 
  # mkdirHomerDMR = paste("mkdir region")
  # ezSystem(mkdirHomerDMR)
  
  # mkdirHomerDML = paste("mkdir locus")
  # ezSystem(mkdirHomerDML)

  # region_size <- 200
  # motif_length <- "8,10,12"
  # if (param$biomart_selection == "mmusculus_gene_ensembl") {
  #   genomeHomer <- "mm10"
  # }
  # else if (param$biomart_selection == "hsapiens_gene_ensembl") {
  #   genomeHomer <- "hg38"
  # }
  # else {
  #   genomeHomer <- file.path("/srv/GT/reference", dirname(dirname(param$refBuild)), 'Sequence/WholeGenomeFasta/genome.fa')
  # }
  # genomeHomer <- file.path("/srv/GT/reference", dirname(dirname(param$refBuild)), 'Sequence/WholeGenomeFasta/genome.fa')
  
  # genomeHomer <- "/srv/GT/reference/Mus_musculus/GENCODE/GRCm39/Sequence/WholeGenomeFasta/genome.fa"
  # dirname(dirname(genomeHomer))
  # output_dir <- "homer"
  ### needed: bed file DMR/L, bed file background x2

  # bedDMR <- "significantRegionsCpG.bed"
  # bedBGDMR <- "dmRegionsCpG.bed"
  # cmd <- paste("findMotifsGenome.pl", bedDMR, genomeHomer, "region", "-size", region_size, "-len", motif_length, "-bg", bedBGDMR)
  # cmd <- paste("findMotifsGenome.pl", bedDMR, genomeHomer, "homer", "-size 200", "-bg", bedBGDMR, "-len", motif_length, "-keepOverlappingBg", "-preparsedDir .")
  # ezSystem(cmd)
  # module load Tools/HOMER/4.11
  # findMotifsGenome.pl significantRegionsCpG.bed /srv/GT/reference/Mus_musculus/GENCODE/GRCm39/Sequence/WholeGenomeFasta/genome.fa homer -size 200 -bg dmRegionsCpG.bed -len 8,10,12 -keepOverlappingBg
  # findMotifsGenome.pl significantRegionsCpG.bed /srv/GT/reference/Mus_musculus/GENCODE/GRCm39/Sequence/WholeGenomeFasta homer -size 200 -bg dmRegionsCpG.bed -len 8,10,12 -keepOverlappingBg
  # findMotifsGenome.pl significantRegionsCpG.bed /srv/GT/reference/Mus_musculus/GENCODE/GRCm39/Sequence/WholeGenomeFasta homer -size 200 -len 8,10,12 -preparsedDir .
  # findMotifsGenome.pl significantRegionsCpG.bed /srv/GT/reference/Mus_musculus/GENCODE/GRCm39/Sequence/WholeGenomeFasta/genome.fa homer -size 200 -len 8,10,12 -preparsedDir .
  # findMotifsGenome.pl significantRegionsCpG.bed /srv/GT/reference/Mus_musculus/GENCODE/GRCm39/Sequence/WholeGenomeFasta/genome.fa homer -len 8,10,12 -preparsedDir . -nlen 0
  # findMotifsGenome.pl dmRegionsCpG.bed mm10 homer
  
  # bsseqCpG <- readRDS(paste0("/scratch/jonas/17679/DNAme/", "bsseqCpG", ".rds"))
  # greatResultBP_CpG <- readRDS(paste0("/scratch/jonas/17679/DNAme/", "greatResultBP_CpG", ".rds"))
  # greatResultMF_CpG <- readRDS(paste0("/scratch/jonas/17679/DNAme/", "greatResultMF_CpG", ".rds"))
  # greatResultCC_CpG <- readRDS(paste0("/scratch/jonas/17679/DNAme/", "greatResultCC_CpG", ".rds"))
  # enrichmentTable_BP <- getEnrichmentTable(greatResultBP_CpG)
  # enrichmentTable_CC <- getEnrichmentTable(greatResultCC_CpG)
  # enrichmentTable_MF <- getEnrichmentTable(greatResultMF_CpG)
  
  

  # cmd <- paste("findMotifsGenome.pl", bedDML, genomeHomer, "locus", "-size", region_size, "-len", motif_length, "-bg", bedBGDML)
  # ezSystem(cmd)
  
  # findMotifsGenome.pl significantLociCpG.bed /srv/GT/reference/Arabidopsis_thaliana/TAIR/TAIR10/Sequence/WholeGenomeFasta/genome.fa CpG -size 200 -len 8,10,12 -bg dmLociCpG.bed
  # findMotifsGenome.pl significantLociCHG.bed /srv/GT/reference/Arabidopsis_thaliana/TAIR/TAIR10/Sequence/WholeGenomeFasta/genome.fa CHG -size 200 -len 8,10,12 -bg dmLociCHG.bed
  # 
  # findMotifsGenome.pl sl.bed /srv/GT/reference/Arabidopsis_thaliana/TAIR/TAIR10/Sequence/WholeGenomeFasta/genome.fa CHG -size 200 -len 8,10,12 -bg dm.bed
  # 
  # findMotifsGenome.pl dmRegionsCHG.bed /srv/GT/reference/Arabidopsis_thaliana/TAIR/TAIR10/Sequence/WholeGenomeFasta/genome.fa CHG -size 200 -len 8,10,12
  # findMotifsGenome.pl dmRegionsCHG.bed /srv/GT/reference/Arabidopsis_thaliana/TAIR/TAIR10/Sequence/WholeGenomeFasta CHG -size 200 -len 8,10,12
  
  
  ######################## Functional analysis (GREAT) ##############################
  # setwdNew(basename(output$getColumn("Report")))
  # mkdirGreat = paste("mkdir great")
  # ezSystem(mkdirGreat)
  # 
  # setwd("great")
  # mkdirGreatDMR = paste("mkdir region")
  # ezSystem(mkdirGreatDMR)
  
  # mkdirGreatDML = paste("mkdir locus")
  # ezSystem(mkdirGreatDML)
  
  # greatFun <- function(dmRegions, significantRegions, context) { # dmType = region / locus
  #   greatResult_BP <- great(gr = significantRegions, gene_sets = "BP", biomart_dataset = param$biomart_selection,
  #                           background = dmRegions, cores = param$cores)
  #   greatResult_CC <- great(gr = significantRegions, gene_sets = "CC", biomart_dataset = param$biomart_selection,
  #                           background = dmRegions, cores = param$cores)
  #   greatResult_MF <- great(gr = significantRegions, gene_sets = "MF", biomart_dataset = param$biomart_selection,
  #                           background = dmRegions, cores = param$cores)
  #   
  #   enrichmentTable_BP <- getEnrichmentTable(greatResult_BP)
  #   enrichmentTable_CC <- getEnrichmentTable(greatResult_CC)
  #   enrichmentTable_MF <- getEnrichmentTable(greatResult_MF)
  #   
  #   if(param$biomart_selection=="athaliana_eg_gene") {
  #     reactome <- "https://plantreactome.gramene.org/download/current/gene_ids_by_pathway_and_species.tab"
  #     react <- data.frame(data.table::fread(input = reactome, header = F, nThread = 16))
  #     rdb <- react[grep(pattern = "^R-ATH", x = react$V1), ]
  #     reactome_pathways <- split(rdb$V4, paste(rdb$V1, rdb$V2, sep = ": "))
  #     
  #     ## KEGG pathways
  #     kg <- keggList("organism")
  #     
  #     pathway2gene <- keggLink("pathway", "ath")
  #     pathwayName <- keggList("pathway", "ath")
  #     df1 <- data.frame(
  #       gene = gsub("ath:", "", names(pathway2gene)),
  #       pathID = gsub("path:", "", pathway2gene)
  #     )
  #     
  #     df2 <- data.frame(
  #       pathID = gsub("path:", "", names(pathwayName)),
  #       name = pathwayName
  #     )
  #     
  #     df_kegg <- merge(df2, df1)
  #     
  #     kegg_pathways <- split(df_kegg$gene, paste(df_kegg$pathID, df_kegg$name,
  #                                                sep = ": "
  #     ))
  #     
  #     greatResult_RE <- great(gr = significantRegions, gene_sets = reactome_pathways, extended_tss = extendedTSS,
  #                             background = dmRegions, cores = cores)
  #     
  #     
  #     greatResult_RE <- great(gr = significantRegions, gene_sets = reactome_pathways, tss_source = "TxDb.Athaliana.BioMart.plantsmart51",
  #                             background = dmRegions, cores = param$cores)
  #     greatResult_KE <- great(gr = significantRegions, gene_sets = kegg_pathways, tss_source = "TxDb.Athaliana.BioMart.plantsmart51",
  #                             background = dmRegions, cores = param$cores)
  #     
  #     enrichmentTable_RE <- getEnrichmentTable(greatResult_RE)
  #     enrichmentTable_KE <- getEnrichmentTable(greatResult_KE)
  #     
  #     # setwd("great")
  #     saveRDS(greatResult_RE, file = paste0("greatResultRE_", context, ".rds"))
  #     saveRDS(greatResult_KE, file = paste0("greatResultKE_", context, ".rds"))
  #     # saveRDS(enrichmentTable_RE, file = paste0(dmType, "/", "enrichmentTable_RE_", ".rds"))
  #     # saveRDS(enrichmentTable_KE, file = paste0(dmType, "/", "enrichmentTable_KE_", ".rds"))
  #   }
  #   
  #   # setwd("great")
  #   saveRDS(greatResult_BP, file = paste0("greatResultBP_", context, ".rds"))
  #   saveRDS(greatResult_CC, file = paste0("greatResultCC_", context, ".rds"))
  #   saveRDS(greatResult_MF, file = paste0("greatResultMF_", context, ".rds"))
  #   # saveRDS(enrichmentTable_BP, file = paste0(dmType, "/", "enrichmentTable_BP", ".rds"))
  #   # saveRDS(enrichmentTable_CC, file = paste0(dmType, "/", "enrichmentTable_CC", ".rds"))
  #   # saveRDS(enrichmentTable_MF, file = paste0(dmType, "/", "enrichmentTable_MF", ".rds"))
  #   
  # }
  # 
  # greatFun(dmRegions = dmRegionsCpG, significantRegions = significantRegionsCpG, context = "CpG")
  # greatFun(dmRegions = dmRegionsCHG, significantRegions = significantRegionsCHG, context = "CHG")
  # greatFun(dmRegions = dmRegionsCHH, significantRegions = significantRegionsCHH, context = "CHH")

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

  # setwd again before saving these
  # setwd(basename(output$getColumn("Report")))
  setwd(cwd)
  
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

