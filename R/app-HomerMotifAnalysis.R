###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodHomerMotifAnalysis <- function(input = NA, output = NA, param = NA,
                          htmlFile = "00index.html") {
  
  ### libraries
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
  require(parallel)
  require(rtracklayer)
  require(ChIPpeakAnno)
  
  ##### test data dmrseq
  # generated from bash script
  dataDirDMRseq <- "~/data/dmrseq"
  dmRegionsFilePath <- file.path(dataDirDMRseq, "dmRegions.rds")
  significantRegionsFilePath <- file.path(dataDirDMRseq, "significantRegions.rds")
  dmRegions <- readRDS(dmRegionsFilePath)
  significantRegions <- readRDS(significantRegionsFilePath)
  
  ##### test data methylKit
  # generated from RStudio
  dataDirMethylKit <- "~/data/MethylKit/doc"
  diffMethLoci <- readRDS(file.path(dataDirMethylKit, paste0("diffMethLoci", ".rds")))
  diffMeth_25p <- readRDS(file.path(dataDirMethylKit, paste0("diffMeth_25p", ".rds")))
  diffMethPerChr <- readRDS(file.path(dataDirMethylKit, paste0("diffMethPerChr", ".rds")))
  gene.obj <- readRDS(file.path(dataDirMethylKit, paste0("gene.obj", ".rds")))
  diffMeth_25p_annotated <- readRDS(file.path(dataDirMethylKit, paste0("diffMeth_25p_annotated", ".rds")))
  
  ########################################
  # configureHomer.pl, FASA, GTF
  
  # cmd <- paste("plink --vcf", file.path("/srv/gstore/projects", input$getColumn("Filtered VCF")), "--double-id", "--allow-extra-chr", "--cluster", "--mds-plot", 2)
  # # this saves it to plink.mds
  # result <- ezSystem(cmd)
  # gc()
  
  # find_motifs_genome(x, path, genome, motif_length = c(8, 10, 12),
  #                    scan_size = 100, optimize_count = 8, background = "automatic",
  #                    local_background = FALSE, only_known = FALSE, only_denovo = FALSE,
  #                    fdr_num = 0, cores = parallel::detectCores(),
  #                    cache = .calc_free_mem()/4, overwrite = FALSE,
  #                    keep_minimal = FALSE)
  
  # findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]
  
  ############ monaLisa using methylKit output data
  # peak_bins <- bin(x = peak_change, binmode = "equalN", nElement = 400)
  # se <- calcBinnedMotifEnrR(seqs = peak_seqs,
  #                           bins = peak_bins,
  #                           pwmL = pwms)
  
  # lmr <- as(diffMethLoci,"GRanges")
  # summary(lmr$meth.diff)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # -78.262  -3.934   0.000   1.210   4.054  64.266 
  #
  # lmr_25p <- as(diffMeth_25p,"GRanges")
  # # combine the two sets or genomic regions
  # lmrsel2 <- c(lmr, lmr_25p)
  # 
  # bins2 <- rep(c("unchanged", "up"), c(length(lmr), length(lmr_25p)))
  # bins2 <- factor(bins2)
  # table(bins2)
  # 
  # se3 <- calcBinnedMotifEnrR(seqs = lmrseqs3,
  #                            pwmL = pwms[rownames(seSel)],
  #                            background = "genome",
  #                            genome = BSgenome.Mmusculus.UCSC.mm10,
  #                            genome.regions = NULL, # sample from full genome
  #                            genome.oversample = 2, 
  #                            BPPARAM = BiocParallel::SerialParam(RNGseed = 42),
  #                            verbose = TRUE)
  #
  # hist(lmr$meth.diff, 100, col = "gray", main = "",
  #      xlab = "Methylation difference", ylab = "Number of DMLs")
  
  ######################################## using HOMER ########################################
  # cmd <- paste("plink --vcf", file.path("/srv/gstore/projects", input$getColumn("Filtered VCF")), "--double-id", "--allow-extra-chr", "--cluster", "--mds-plot", 2)
  # # this saves it to plink.mds
  # result <- ezSystem(cmd)
  # gc()
  
  ######### findMotifsGenome.pl
  # findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]
  # x = x
  # path = 
  #   genome = 
  #   motif_length = c(8, 10, 12)
  # scan_size = 100
  # optimize_count = 8
  # background = bg_n
  # local_background = FALSE
  
  ####
  # lmr <- as(diffMethLoci,"GRanges")
  # lmr_25p <- as(diffMeth_25p,"GRanges")
  # df <- data.frame(seqnames=seqnames(lmr),
  #                  starts=start(lmr)-1,
  #                  ends=end(lmr),
  #                  names=c(rep(".", length(lmr))), # col 4 = unique Peak ID (?)
  #                  scores=c(rep(".", length(lmr))), # col 5 = not used
  #                  strands=strand(lmr))
  # df_25p <- data.frame(seqnames=seqnames(lmr_25p),
  #                  starts=start(lmr_25p)-1,
  #                  ends=end(lmr_25p),
  #                  names=c(rep(".", length(lmr_25p))), # col 4 = unique Peak ID (?)
  #                  scores=c(rep(".", length(lmr_25p))), # col 5 = not used
  #                  strands=strand(lmr_25p))
  # 
  # dataDirMethylKit <- "~/data/MethylKit/doc"
  # setwd(dataDirMethylKit)
  # write.table(df, file="dml.bed", quote=F, sep="\t", row.names=F, col.names=F)
  # write.table(df_25p, file="dml_25p.bed", quote=F, sep="\t", row.names=F, col.names=F)
  ####
  # bed_file <- read.table("~/data/MethylKit/doc/dml.bed")
  # bed_file_25p <- read.table("~/data/MethylKit/doc/dml_25p.bed")
  bed_file <- file.path("~/data/MethylKit/doc/dml.bed")
  bed_file_25p <- file.path("~/data/MethylKit/doc/dml_25p.bed")
  
  region_size <- 200
  motif_length <- "8,10,12"
  genome <- "hg18"
  output_dir <- "~/data/homer/test1/"
  # findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]
  # findMotifsGenome.pl ~/data/MethylKit/doc/dml.bed hg18 . -size 100 
  
  cmd <- paste("findMotifsGenome.pl", bed_file_25p, genome, output_dir, "-size", region_size, "-len", motif_length, "-bg", bed_file)
  ezSystem(cmd)
  # gc()
  

  
  # findMotifs.pl: 3 mandatory arguments: A gene ID input file, the name of the promoter set (which is tied to an organism), and an output directory for all of the output files
  

  # setwd("/home/jobucher/data/HomerMotifAnalysis/mm")
  # saveRDS(HomerMotifAnalysisResult_BP, file = paste0("HomerMotifAnalysisResult_BP", ".rds"))

  ## Copy the style files and templates
  styleFiles <- file.path(
    system.file("templates", package = "ezRun"),
    c(
      "fgcz.css", "HomerMotifAnalysis.Rmd",
      "fgcz_header.html", "banner.png"
    )
  )
  file.copy(from = styleFiles, to = ".", overwrite = TRUE)
  
  ### generate the main reports
  # rmarkdown::render(
  #   input = "HomerMotifAnalysis.Rmd", envir = new.env(),
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

EzAppHomerMotifAnalysis <-
  setRefClass("EzAppHomerMotifAnalysis",
              contains = "EzApp",
              methods = list(
                initialize = function() {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodHomerMotifAnalysis
                  name <<- "EzAppHomerMotifAnalysis"
                  appDefaults <<- rbind(perLibrary = ezFrame(Type = "logical", DefaultValue = TRUE, Description = "HomerMotifAnalysis brabra"))
                }
              )
  )

