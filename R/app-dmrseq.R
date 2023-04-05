###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethoddmrseq <- function(input = NA, output = NA, param = NA,
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
  
  if (param$cores > 1){
    BPPARAM <- MulticoreParam(workers = param$cores)
  } else {
    BPPARAM <- SerialParam() 
  }
  register(BPPARAM)
  require(future)
  plan("multicore", workers = param$cores)

  bsseq <- bsseq::read.bismark(files = input$getFullPaths("COV [File]"),
                                           rmZeroCov = FALSE,
                                           strandCollapse = FALSE,
                                           verbose = FALSE,
                                           colData = param@Factor)
  
  lociCoverage <- which(DelayedMatrixStats::rowSums2(getCoverage(coverageDmrseq_mm, type="Cov")==0) == 0)
  
  bsseqFiltered <- bsseq[lociCoverage, ]
  
  # testCovariate <- param$testCovariate
  
  dmRegions <- dmrseq(bs = bsseqFiltered,
                               testCovariate = param$testCovariate, 
                               cutoff = param$cutoff,
                               minNumRegion = 5,
                               smooth = TRUE,
                               bpSpan = 1000,
                               minInSpan = 30,
                               maxGapSmooth = 5000,
                               maxGap = 1000,
                               verbose = TRUE,
                               maxPerms = 10,
                               # matchCovariate = deparse(substitute(matchCovariate)), # opt
                               # eg if samples from different gender, but not covariate of interest
                               # -> avoids to do permutation with all-male vs all-female
                               BPPARAM = BPPARAM,
                               stat = "stat",
                               block = FALSE,
                               blockSize = 5000,
                               chrsPerChunk = 1
  )
  
  # saveRDS(bsseq, file="bsseq.rds")
  significantRegions <- dmRegions[dmRegions$qval < 0.05, ]
  # sum(significantRegions$stat > 0) / length(significantRegions)
  
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
  rmarkdown::render(
    input = "dmrseq.Rmd", envir = new.env(),
    output_dir = ".", output_file = htmlFile, quiet = TRUE
  )
  
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

EzAppdmrseq <-
  setRefClass("EzAppdmrseq",
              contains = "EzApp",
              methods = list(
                initialize = function() {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethoddmrseq
                  name <<- "EzAppdmrseq"
                  appDefaults <<- rbind(perLibrary = ezFrame(Type = "logical", DefaultValue = TRUE, Description = "dmrseq brabra"))
                }
              )
  )

