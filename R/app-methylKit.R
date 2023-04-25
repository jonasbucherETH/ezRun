###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMethylKit <- function(input = NA, output = NA, param = NA,
                          htmlFile = "00index.html") {
  
  ### libraries
  library("methylKit")
  library("parallel")
  library("genomation") # annotation of DML/R
  
  
  # #setwdNew(basename(output$getColumn("Report")))
  dataset <- input$meta
  ans4Report <- list() # a list of results for rmarkdown report
  ans4Report[["dataset"]] <- dataset
  
  cores <- detectCores(all.tests = FALSE, logical = TRUE) - 1
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add = TRUE)
  
  dataDirSave <- "/home/jobucher/data/MethylKit"
  
  ### read data (testing)
  file.list <- list(system.file("extdata", 
                              "test1.myCpG.txt", package = "methylKit"),
                  system.file("extdata",
                              "test2.myCpG.txt", package = "methylKit"),
                  system.file("extdata", 
                              "control1.myCpG.txt", package = "methylKit"),
                  system.file("extdata", 
                              "control2.myCpG.txt", package = "methylKit"))
  
  sampleIDs <- list("test1","test2","ctrl1","ctrl2")
  
  methylRaw <- methRead(file.list,
                 sample.id=list("test1","test2","ctrl1","ctrl2"),
                 assembly="hg18",
                 treatment=c(1,1,0,0),
                 context="CpG",
                 mincov = 10
  )
  
  ### read data (final)
  # processBismarkAln / methRead
  # eg SAM files from Bismark, but also CpG.txt (see above)
  
  
  ##### 2.4 Descriptive statistics on samples
  methylationStats <- getMethylationStats(myobj[[2]],plot=FALSE,both.strands=FALSE)
  # cannot be saved
  
  # getMethylationStats(myobj[[2]],plot=TRUE,both.strands=FALSE)
  
  # coverageStats <- getCoverageStats(myobj[[2]],plot=F,both.strands=FALSE)
  #### NOTE: this doesn't work, immediately prints and doesn't save
  # coverageStatsList <- c()
  # for(i in 1:length(myobj)) { # myobj = list
  #   name <- paste0("coverageStats_", sampleIDs[i])
  #   coverageStats <- getCoverageStats(myobj[[i]], plot=F, both.strands=FALSE)
  #   assign(name, coverageStats)
  #   coverageStatsList <- cbind(coverageStatsList, name)
  # }

  ##### 2.5 Filtering samples based on read coverage
  filteredMethylRaw  <- filterByCoverage(
      methylRaw,
      lo.count=10, # Bases/regions having lower coverage than this count is discarded
      lo.perc=NULL, # Bases/regions having lower coverage than this percentile is discarded
      hi.count=NULL, # might want to filter out very high coverages as well (PCR bias)
      hi.perc=99.9 # Bases/regions having higher coverage than this percentile is discarded
  )
  
  
  # normalizeCoverage() function to normalize coverage between samples
  normalizedMethylRaw <- normalizeCoverage(
    filteredMethylRaw,
    method = "median" # median / mean
    # chunk.size = 1e+06,
    # save.db = TRUE,
  )
  
  ##### 3.1 Merging samples
  # merge all samples to one object for base-pair locations that are covered in all samples
  methylAll <- unite(normalizedMethylRaw, destrand=FALSE) # destrand=TRUE (default=FALSE) merges reads on both strands of a CpG dinucleotide
  # head(methylAll) # methylBase object
  
  ##### 3.2 Sample Correlation
  sampleCorrelation <- getCorrelation(methylAll, plot=FALSE,
                                      method = "pearson") # default:"pearson", "kendall" and "spearman"
  ### not saved
  # getCorrelation(methylAll, plot=TRUE, method = "pearson")
  
  ##### 3.3 Clustering samples
  # clusterSamples(methylAll, dist="correlation", method="ward", plot=TRUE)
  clusterDendrogram <- clusterSamples(methylAll, dist="correlation", method="ward", plot=FALSE)
  # PCASamples(methylAll, screeplot=TRUE)
  # PCASamples(methylAll)
  
  # If you are convinced that some principal components are accounting for batch effects,
  # you can remove those principal components from your data using removeComp
  
  
  ##### 3.4 Batch effects
  sampleAnnotation <- data.frame(batch_id=c("a","a","b","b"),
                              age=c(19,34,23,40))
  
  batchAssociation <- assocComp(mBase = methylAll, sampleAnnotation)

  
  newMethylAll <- removeComp(methylAll, comp=1) # remove PC1
  
  ##### 3.5 Tiling windows analysis -> to get DMRs
  
  ##### 3.6 Finding differentially methylated bases or regions
  diffMethLoci <- calculateDiffMeth(
    methylAll,
    covariates = NULL,
    overdispersion = "none" # c("none", "MN", "shrinkMN")
    # adjust = "SLIM", # multiple testing adjustment; default = "SLIM"
    # c("SLIM", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none", "qvalue")
    # effect = c("wmean", "mean", "predicted"),
    # test = c("F", "Chisq", "fast.fisher", "midPval"), 
    # mc.cores = cores,
    # slim = TRUE, # earlier version
    # weighted.mean = TRUE, # earlier version
    # chunk.size = 1e+06, # Number of rows to be taken as a chunk for processing the methylBaseDB objects
    # save.db = FALSE,
  )
  
  # After q-value calculation, we can select the differentially methylated
  # regions/bases based on q-value and percent methylation difference cutoffs
  
  diffMeth_25p <- getMethylDiff(
      diffMethLoci,
      difference = 25, # cutoff for absolute value of methylation percentage change between test and control
      qvalue = 0.01, # cutoff for qvalue of differential methylation statistic
      type = "all", # all / hyper / hypo
      chunk.size = 1e+06, # Number of rows to be taken as a chunk for processing the methylDiffDB objects
      save.db = FALSE # default=T for methyl(Diff)DB objects; false for methylDiff
      # (...) optional args when save.db = T
  )
  
  # also possible per chromosome
  diffMethPerChr <- diffMethPerChr(diffMethLoci,plot=F,qvalue.cutoff=0.01, meth.cutoff=25)
  
  ##### 3.8 Accounting for covariates
  covariates <- data.frame(age=c(30,80,34,30,80,40))
  sim.methylBase <- dataSim(replicates=6,sites=1000,
                          treatment=c(rep(1,3),rep(0,3)),
                          covariates=covariates,
                          sample.ids=c(paste0("test",1:3),paste0("ctrl",1:3))
  )
  diffMeth_covariates <- calculateDiffMeth(sim.methylBase,
                                  covariates=covariates,
                                  overdispersion="MN",test="Chisq",mc.cores=1)
  
  ##### 4 Annotating differentially methylated bases or regions
  gene.obj=readTranscriptFeatures(system.file("extdata", "refseq.hg18.bed.txt", 
                                              package = "methylKit"))
  
  diffMeth_25p_annotated <- annotateWithGeneParts(as(diffMeth_25p,"GRanges"),gene.obj)
  # Formal class AnnotationByGeneParts
  
  # Similarly, we can read the CpG island annotation and
  # annotate our differentially methylated bases/regions with them.
  cpg.obj <- readFeatureFlank(system.file("extdata", "cpgi.hg18.bed.txt", 
                                       package = "methylKit"),
                           feature.flank.name=c("CpGi","shores"))
  
  # convert methylDiff object to GRanges and annotate
  diffCpGann <- annotateWithFeatureFlank(as(diffMeth_25p,"GRanges"),
                                      cpg.obj$CpGi,cpg.obj$shores,
                                      feature.name="CpGi",flank.name="shores")
  
  ##### 4.1 Regional analysis
  promoters <- regionCounts(myobj,gene.obj$promoters)
  
  ##### 4.2 Convenience functions for annotation objects
  diffAnn <- annotateWithGeneParts(as(diffMeth_25p,"GRanges"),gene.obj)
  # get the distance to TSS and nearest gene name
  associationWithTSS <- getAssociationWithTSS(diffAnn)
  
  # get percentage/number of differentially methylated regions that overlap with intron/exon/promoters
  targetAnnotationStats <- getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)
  # promoter       exon     intron intergenic 
  # 28.15       0.00      11.11      60.74
  
  # plot the percentage of differentially methylated bases overlapping with exon/intron/promoters
  plotTargetAnnotation(diffAnn,precedence=TRUE,
                       main="differential methylation annotation")
  
  # plot the CpG island annotation
  # shows what percentage of differentially methylated bases
  # are on CpG islands, CpG island shores and other regions
  plotTargetAnnotation(diffCpGann,col=c("green","gray","white"),
                       main="differential methylation annotation")
  
  # get percentage of intron/exon/promoters that overlap with differentially methylated bases
  featsWithTargetsStats <- getFeatsWithTargetsStats(diffAnn,percentage=TRUE)
  # promoter     exon   intron 
  # 0.29     0.03     0.17 
  
  ########## 5 methylKit convenience functions
  ##### 5.1 Coercing methylKit objects to GRanges
  
  
  ############## save objects
  # setwdNew("/home/jobucher/data/MethylKit/doc")
  # on.exit(setwd(cwd), add = TRUE)
  
  saveRDS(filteredMethylRaw, file = paste0("filteredMethylRaw", ".rds"))
  saveRDS(normalizedMethylRaw, file = paste0("normalizedMethylRaw", ".rds"))
  saveRDS(methylAll, file = paste0("methylAll", ".rds"))
  saveRDS(clusterDendrogram, file = paste0("clusterDendrogram", ".rds"))
  saveRDS(batchAssociation, file = paste0("batchAssociation", ".rds"))
  saveRDS(diffMethLoci, file = paste0("diffMethLoci", ".rds"))
  saveRDS(diffMeth_25p, file = paste0("diffMeth_25p", ".rds"))
  saveRDS(diffMethPerChr, file = paste0("diffMethPerChr", ".rds"))
  saveRDS(diffMeth_covariates, file = paste0("diffMeth_covariates", ".rds"))
  saveRDS(gene.obj, file = paste0("gene.obj", ".rds"))
  saveRDS(diffMeth_25p_annotated, file = paste0("diffMeth_25p_annotated", ".rds"))
  
  
  
  ## Copy the style files and templates
  styleFiles <- file.path(
    system.file("templates", package = "ezRun"),
    c(
      "fgcz.css", "MethylKit.Rmd",
      "fgcz_header.html", "banner.png"
    )
  )
  file.copy(from = styleFiles, to = ".", overwrite = TRUE)
  
  ### generate the main reports
  # rmarkdown::render(
  #   input = "MethylKit.Rmd", envir = new.env(),
  #   output_dir = ".", output_file = htmlFile, quiet = TRUE
  # )
  
  html_files <- c("00index.html",  "banner.png",  "fgcz.css",  "fgcz_header.html")
  
  ## try this, taken (and adapted) from app-Vpipe
  # dir.create(param[['name']])
  # file.copy(from = html_files, to = param[['name']])
  # cmd <- paste('mv rmarkdownLib', param[['name']])
  
  ## this was before
  # file.copy(from = html_files, to = "plink")
  # cmd <- "mv rmarkdownLib plink"
  # ezSystem(cmd)
  
  # file.copy(from = html_files, to = output_dir)
  # cmd <- paste("mv rmarkdownLib ", output_dir) 
  # ezSystem(cmd)
  
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

EzAppMethylKit <-
  setRefClass("EzAppMethylKit",
              contains = "EzApp",
              methods = list(
                initialize = function() {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMethylKit
                  name <<- "EzAppMethylKit"
                  appDefaults <<- rbind(perLibrary = ezFrame(Type = "logical", DefaultValue = TRUE, Description = "MethylKit brabra"))
                }
              )
  )

