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
  
  methylRaw <- methRead(file.list,
                 sample.id=list("test1","test2","ctrl1","ctrl2"),
                 assembly="hg18",
                 treatment=c(1,1,0,0),
                 context="CpG",
                 mincov = 10
  )
  
  ### read data (final)
  # processBismarkAln / methRead
  
  
  ##### 2.4 Descriptive statistics on samples
  methylationStats <- getMethylationStats(myobj[[2]],plot=FALSE,both.strands=FALSE)
  getMethylationStats(myobj[[2]],plot=TRUE,both.strands=FALSE)
  getCoverageStats(myobj[[2]],plot=TRUE,both.strands=FALSE)
  
  ##### 2.5 Filtering samples based on read coverage
  filteredMethylRaw  <- filterByCoverage(
      methylRaw,
      lo.count=10, # Bases/regions having lower coverage than this count is discarded
      lo.perc=NULL, # Bases/regions having lower coverage than this percentile is discarded
      hi.count=NULL, # might want to filter out very high coverages as well (PCR bias)
      hi.perc=99.9 # Bases/regions having higher coverage than this percentile is discarded
  )
  
  ##### 3.1 Merging samples
  # merge all samples to one object for base-pair locations that are covered in all samples
  methylAll <- unite(methylRaw, destrand=FALSE) # destrand=TRUE (default=FALSE) merges reads on both strands of a CpG dinucleotide
  # head(methylAll) # methylBase object
  
  ##### 3.2 Sample Correlation
  sampleCorrelation <- getCorrelation(methylAll, plot=FALSE,
                                      method = "pearson") # default:"pearson", "kendall" and "spearman"
  # getCorrelation(methylAll, plot=TRUE, method = "pearson")
  
  ##### 3.3 Clustering samples
  # clusterSamples(methylAll, dist="correlation", method="ward", plot=TRUE)
  clusterDendrogram <- clusterSamples(methylAll, dist="correlation", method="ward", plot=FALSE)
  PCASamples(methylAll, screeplot=TRUE)
  PCASamples(methylAll)
  # If you are convinced that some principal components are accounting for batch effects,
  # you can remove those principal components from your data using removeComp
  
  
  ##### 3.4 Batch effects
  sampleAnnotation <- data.frame(batch_id=c("a","a","b","b"),
                              age=c(19,34,23,40))
  
  as <- assocComp(mBase = methylAll, sampleAnnotation)
  as
  
  newMethylAll <- removeComp(methylAll, comp=1) # remove PC1
  
  ##### 3.5 Tiling windows analysis -> to get DMRs
  
  ##### 3.6 Finding differentially methylated bases or regions
  diffMethLoci <- calculateDiffMeth(methylAll)
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

