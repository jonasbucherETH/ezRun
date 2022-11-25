###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodVcfStats <- function(input = NA, output = NA, param = NA,
                           htmlFile = "00index.html") {
  # #setwdNew(basename(output$getColumn("Report")))
  dataset <- input$meta
  ans4Report <- list() # a list of results for rmarkdown report
  ans4Report[["dataset"]] <- dataset

  output_dir <- basename(output$getColumn("Report"))
  prefix <- file.path(output_dir, "vcf_stats")

  # For Rmd
  # SNP counts
  snp_counts <- file.path(output_dir, "vcf_stats.snps")

  # InDel counts
  # ToDo

  # Private SNP counts
  private_snp_counts <- file.path(output_dir, "vcf_stats.private")

  # Shared SNP counts
  shared_snp_counts <- file.path(output_dir, "vcf_stats.shared")

  # Transions/Transversions
  tstv <- file.path(output_dir, "vcf_stats.samples-tstv")

  # run vcf-stats
  cmd <- paste("vcf-stats", file.path("/srv/gstore/projects", input$getColumn("Filtered VCF")), "-p", prefix)
  result <- ezSystem(cmd)
  gc()

  ## pca prep
  #library(adegenet)
  #library(ade4)
  library(gdsfmt)
  library(SNPRelate)

  vcf_f <- file.path("/srv/gstore/projects", input$getColumn("Filtered VCF"))
 
  # convert vcf to gds   
  snpgdsVCF2GDS(vcf_f, file.path(output_dir, "snp.gds"),  method="biallelic.only")

  # open gds
  genofile <- snpgdsOpen(file.path(output_dir, "snp.gds"))

  # file for mds
  mds <- file.path(output_dir, "mds")

  # run plink for distance matrix (mds) 
  cmd <- paste("plink --vcf", file.path("/srv/gstore/projects", input$getColumn("Filtered VCF")), "--cluster", "--mds-plot", 10 , "--out", prefix)
  result <- ezSystem(cmd)
  gc()

  #get.attr.gdsn(index.gdsn(genofile, "snp.chromosome"))

  # Get the attribute of genotype
  #get.attr.gdsn(index.gdsn(genofile, "genotype"))

  
  #ccm_pca <- snpgdsPCA(genofile)

  #pca_dat <- file.path(output_dir, "vcf_stats.pca")

  # turn SNP data into genind format
  #snp_genind <- df2genind(snp_df, ploidy=2)
   
  #snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
  #snpset.id <- unlist(unname(snpset))

  #pca <- snpgdsPCA(genofile)
  #pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)

  # replace NA
  #snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
  #snpset.id <- unlist(unname(snpset))
  
  #pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)

  # case: no prior population information
  
#  pca_tab <- data.frame(sample.id = pca$sample.id,
#		        EV1 = pca$eigenvect[,1],    # the first eigenvector
#			    EV2 = pca$eigenvect[,2],    # the second eigenvector
#			    stringsAsFactors = FALSE)


  ## Copy the style files and templates
  styleFiles <- file.path(
    system.file("templates", package = "ezRun"),
    c(
      "fgcz.css", "VcfStats.Rmd",
      "fgcz_header.html", "banner.png"
    )
  )
  file.copy(from = styleFiles, to = ".", overwrite = TRUE)

  ### generate the main reports
  rmarkdown::render(
    input = "VcfStats.Rmd", envir = new.env(),
    output_dir = ".", output_file = htmlFile, quiet = TRUE
  )

  html_files <- c("00index.html",  "banner.png",  "fgcz.css",  "fgcz_header.html")
  file.copy(from = html_files, to = "vcf_stats")
  cmd <- "mv rmarkdownLib vcf_stats"
  ezSystem(cmd)

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

EzAppVcfStats <-
  setRefClass("EzAppVcfStats",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodVcfStats
        name <<- "EzAppVcfStats"
        appDefaults <<- rbind(perLibrary = ezFrame(Type = "logical", DefaultValue = TRUE, Description = "VcfStats brabra"))
      }
    )
  )

