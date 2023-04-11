###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodGreat <- function(input = NA, output = NA, param = NA,
                           htmlFile = "00index.html") {
  # #setwdNew(basename(output$getColumn("Report")))
  dataset <- input$meta
  ans4Report <- list() # a list of results for rmarkdown report
  ans4Report[["dataset"]] <- dataset
  
  # output_dir <- basename(output$getColumn("Report"))
  dataDir <- "~/data/dmrseq"
  dmRegionsFilePath <- file.path(dataDir, "dmRegions.rds")
  # significantRegionsFilePath <- file.path(dataDir, "significantRegions.rds")
  dmRegions <- readRDS(dmRegionsFilePath)
  geneSets	<- "GO:BP"
  tssSource <- ""
  
  # mmusculus_gene_ensembl
  
  greatResult <- great(dmRegions, "GO:BP", biomart_dataset = "mmusculus_gene_ensembl")
  
  greatResult <- great(gr = dmRegions, gene_sets = geneSets, tss_source = tssSource, biomart_dataset = NULL,
        min_gene_set_size = 5, mode = "basalPlusExt", basal_upstream = 5000,
        basal_downstream = 1000, extension = 1000000,
        extended_tss = NULL, background = NULL, exclude = "gap",
        cores = 1)
  
  enrichmentTable <- getEnrichmentTable(greatResult)
  head(enrichmentTable)
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add = TRUE)
  
  ### libraries
  library("rGREAT")
  
  ## Copy the style files and templates
  styleFiles <- file.path(
    system.file("templates", package = "ezRun"),
    c(
      "fgcz.css", "Great.Rmd",
      "fgcz_header.html", "banner.png"
    )
  )
  file.copy(from = styleFiles, to = ".", overwrite = TRUE)
  
  ### generate the main reports
  # rmarkdown::render(
  #   input = "Great.Rmd", envir = new.env(),
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

EzAppGreat <-
  setRefClass("EzAppGreat",
              contains = "EzApp",
              methods = list(
                initialize = function() {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodGreat
                  name <<- "EzAppGreat"
                  appDefaults <<- rbind(perLibrary = ezFrame(Type = "logical", DefaultValue = TRUE, Description = "Great brabra"))
                }
              )
  )

