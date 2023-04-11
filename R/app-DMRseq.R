###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodDMRseq <- function(input = NA, output = NA, param = NA,
                           htmlFile = "00index.html") {
  # #setwdNew(basename(output$getColumn("Report")))
  dataset <- input$meta
  ans4Report <- list() # a list of results for rmarkdown report
  ans4Report[["dataset"]] <- dataset
  
  # output_dir <- basename(output$getColumn("Report"))
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add = TRUE)
  
  ### libraries
  library("dmrseq")
  library("BiocParallel")
  
  # if (param$cores > 1){
  #   BPPARAM <- MulticoreParam(workers = param$cores)
  # } else {
  #   BPPARAM <- SerialParam() 
  # }
  # register(BPPARAM)
  # require(future)
  # plan("multicore", workers = param$cores)
  
  # stopifnot(param$sampleGroup != param$refGroup)
  
  ###### grouping input
  # input <- cleanupTwoGroupsInput(input, param)
  # param$testCovariateName <- param$testCovariate
  # param$testCovariate <- input$getColumn(param$testCovariate)
  # bsseqColData <- param@testCovariate
  # 
  # if (ezIsSpecified(param$adjustCovariate) && length(param$adjustCovariate) == 1) {
  #   param$adjustCovariateName <- param$adjustCovariate
  #   param$adjustCovariate <- input$getColumn(param$adjustCovariate)
  #   bsseqColData <- cbind(bsseqColData, param$adjustCovariate)
  # }
  # 
  # if (ezIsSpecified(param$testCovariate) && length(param$testCovariate) == 1) {
  #   param$testCovariateName <- param$testCovariate
  #   param$testCovariate <- input$getColumn(param$testCovariate)
  #   bsseqColData <- cbind(bsseqColData, param$testCovariate)
  # }
  
  # bsseqColData <- as.data.frame(input$getColumn("Treatment"), col.names = "Treatment", row.names = param$samples)
  # print(input$getColumn("Treatment"))
  # print(type(input$getColumn("Treatment")))
  # print(dim)
  # print(input$getFullPaths("COV"))
  # bsseqColData <- as.data.frame(t(input$getColumn("Treatment")), col.names = "Treatment")
  # print(bsseqColData)
  

  # bsseqColData <- as.data.frame("Treatment")
  
  # bsseqColData <- c(param@testCovariate, param@adjustCovariate, param@matchCovariate)
  dataDirBismark <- "/srv/gstore/projects/p1535/Bismark_JBmm_test3_2023-03-27--15-58-43"
  input_dataset <- read_tsv(file.path(dataDirBismark, "input_dataset.tsv"))
  # bsseqColData <- base::as.data.frame(input_dataset$Treatment, row.names = input_dataset$Name, col.names = c("Treatment"))
  bsseqColData <- data.frame(input_dataset$Treatment, row.names = input_dataset$Name)
  colnames(bsseqColData) <- "Treatment"
  
  # covFilesBismark <- input$getFullPaths("COV")
  covFilesBismark <- list.files(dataDirBismark, pattern = "cov", full.names = T)
  
  bsseq <- bsseq::read.bismark(files = covFilesBismark,
                                           rmZeroCov = FALSE,
                                           strandCollapse = FALSE,
                                           verbose = FALSE,
                                           colData = bsseqColData)
  
  lociCoverage <- which(DelayedMatrixStats::rowSums2(getCoverage(bsseq, type="Cov")==0) == 0)
  
  bsseqFiltered <- bsseq[lociCoverage, ]
  
  # testCovariate <- param$testCovariate
  testCovariate <- "Treatment"
  
  dmRegions <- dmrseq(
                 bs = bsseqFiltered,
                 testCovariate = testCovariate, 
                 # A continuous covariate is assumed if the data type in the 'testCovariate' slot is continuous,
                 # with the exception of if there are only two unique values (then a two group comparison is carried out)
                 # adjustCovariate = param$adjustCovariate,
                 cutoff = param$cutoff,
                 minNumRegion = param$minNumRegion,
                 smooth = param$smooth,
                 bpSpan = param$bpSpan,
                 minInSpan = param$minInSpan,
                 maxGapSmooth = param$maxGapSmooth,
                 maxGap = param$maxGap,
                 verbose = FALSE, # keep this
                 maxPerms = param$maxPerms,
                 # matchCovariate = param$matchCovariate, # opt
                 # eg if samples from different gender, but not covariate of interest
                 # -> avoids to do permutation with all-male vs all-female
                 # BPPARAM = BPPARAM,
                 stat = param$stat,
                 block = param$block,
                 blockSize = param$blockSize,
                 chrsPerChunk = param$chrsPerChunk
  )
  
  
  # saveRDS(bsseq, file="bsseq.rds")
  significantRegions <- dmRegions[dmRegions$qval < 0.05, ]
  # sum(significantRegions$stat > 0) / length(significantRegions)
  
  saveRDS(bsseq, file="bsseq.rds")
  saveRDS(bsseqFiltered, file="bsseqFiltered.rds")
  saveRDS(bsseqColData, file="bsseqColData.rds")
  saveRDS(dmRegions, file="dmRegions.rds")
  saveRDS(significantRegions, file="significantRegions.rds")
  
  ## Copy the style files and templates
  styleFiles <- file.path(
    system.file("templates", package = "ezRun"),
    c(
      "fgcz.css", "dmrseq.Rmd",
      "fgcz_header.html", "banner.png"
    )
  )
  file.copy(from = styleFiles, to = ".", overwrite = TRUE)
  
  ### generate the main reports
  # rmarkdown::render(
  #   input = "dmrseq.Rmd", envir = new.env(),
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

EzAppDMRseq <-
  setRefClass("EzAppDMRseq",
              contains = "EzApp",
              methods = list(
                initialize = function() {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodDMRseq
                  name <<- "EzAppDMRseq"
                  appDefaults <<- rbind(perLibrary = ezFrame(Type = "logical", DefaultValue = TRUE, Description = "DMRseq brabra"))
                }
              )
  )

