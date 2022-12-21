###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodPCAMDS <- function(input = NA, output = NA, param = NA,
                           htmlFile = "00index.html") {
  # #setwdNew(basename(output$getColumn("Report")))
  dataset <- input$meta
  ans4Report <- list() # a list of results for rmarkdown report
  ans4Report[["dataset"]] <- dataset

  output_dir <- basename(output$getColumn("Report"))
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add = TRUE)

  
  ### PCA
  library(gdsfmt)
  library(SNPRelate)
  library(adegenet)
  library(ade4)
  library(vcfR)
  
  # vcf_f <- file.path("/srv/gstore/projects", input$getColumn("Filtered VCF"))
  # 
  # grouping_vars <- file.path("/srv/gstore/projects", input$getColumn("Grouping File"))
  # 
  # # convert vcf to gds   
  # snpgdsVCF2GDS(vcf_f, file.path(output_dir, "snp.gds"),  method="biallelic.only")
  # 
  # # open gds
  # genofile <- snpgdsOpen(file.path(output_dir, "snp.gds"))
  # 
  # pca <- snpgdsPCA(genofile, autosome.only = F, verbose = F)
  
  grouping_vars <- read.delim(file.path("/srv/gstore/projects", input$getColumn("Grouping File")))
  
  vcf <- read.vcfR(file.path("/srv/gstore/projects", input$getColumn("Filtered VCF")))
  genind <- vcfR2genind(vcf)
  # pop(genind) <- populations_txt$Population
  
  X <- scaleGen(genind, NA.method="mean")
  pca <- dudi.pca(X, center = TRUE, scale = TRUE, scan = FALSE, nf = 5)
  
  saveRDS(pca, file="PCA.rds")

  saveRDS(grouping_vars, file="grouping_vars.rds")
  
  ### MDS
  # file for mds
  mds <- file.path(output_dir, "plink.mds")
  
  # plink MDS
  # prefix <- file.path(output_dir, "plink")
  # cmd <- paste("plink --vcf", file.path("/srv/gstore/projects", input$getColumn("Filtered VCF")), "--double-id", "--allow-extra-chr", "--cluster", "--mds-plot", 4 , "--out", prefix)
  cmd <- paste("plink --vcf", file.path("/srv/gstore/projects", input$getColumn("Filtered VCF")), "--double-id", "--allow-extra-chr", "--cluster", "--mds-plot", 5)
  # this saves it to plink.mds
  result <- ezSystem(cmd)
  gc()
  
  # plink distance matrix
  cmd <- paste("plink --vcf", file.path("/srv/gstore/projects", input$getColumn("Filtered VCF")), "--double-id", "--allow-extra-chr", "--distance square")
  # this saves it to plink.dist
  result <- ezSystem(cmd)
  gc()
  
  
  # mds <- read.csv(file.path(output_dir, "plink.mds"), sep="")
  
  # saveRDS(mds, file="mds.rds")

  ## Copy the style files and templates
  styleFiles <- file.path(
    system.file("templates", package = "ezRun"),
    c(
      "fgcz.css", "PCAMDS.Rmd",
      "fgcz_header.html", "banner.png"
    )
  )
  file.copy(from = styleFiles, to = ".", overwrite = TRUE)

  ### generate the main reports
  rmarkdown::render(
    input = "PCAMDS.Rmd", envir = new.env(),
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

EzAppPCAMDS <-
  setRefClass("EzAppPCAMDS",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodPCAMDS
        name <<- "EzAppPCAMDS"
        appDefaults <<- rbind(perLibrary = ezFrame(Type = "logical", DefaultValue = TRUE, Description = "PCAMDS brabra"))
      }
    )
  )

