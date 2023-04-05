###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


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


ezMethodDNAme <- function(input = NA, output = NA, param = NA,
                           htmlFile = "00index.html") {
  
  ##### ----- libraries
  library("edgeR")
  library("dmrseq")
  library("tidyverse")
  

  ##### ----- input & output paths
  # #setwdNew(basename(output$getColumn("Report")))
  dataset <- input$meta
  ans4Report <- list() # a list of results for rmarkdown report
  ans4Report[["dataset"]] <- dataset
  
  output_dir <- basename(output$getColumn("Report"))
  
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add = TRUE)
  
  filePathProject <- file.path("/srv/gstore/projects")
  setwd(filePathProject)
  
  input <- "p1535/Bismark_JB_test1_2023-03-07--15-47-45" # for testing
  dataDir <- "/srv/gstore/projects/p1535/Bismark_JB_test1_2023-03-07--15-47-45"
  
  ### input dataset
  inputDataset <- read_tsv(file.path(dataDir, "input_dataset.tsv"))
  dataset <- read_tsv(file.path(dataDir, "dataset.tsv"))
  input <- dataset
  sampleNames <- dataset$Name
  
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
  coverageEdgeR <- ReadBismark2DGE(input$`COV [File]`, sample.names=sampleNames, data.table = TRUE)
  
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
  dataDir <- "/srv/gstore/projects/p1535/Bismark_JBmm_test3_2023-03-27--15-58-43"
  inputDataset <- read_tsv(file.path(dataDir, "input_dataset.tsv"))
  dataset <- read_tsv(file.path(dataDir, "dataset.tsv"))
  sampleNames <- dataset$Name

  coverageFileList_mm <- list.files(dataDir, pattern = "*bismark.cov.gz", full.names = T)
  coverageDmrseq_mm <- bsseq::read.bismark(files = coverageFileList_mm,
                                        rmZeroCov = F,
                                        strandCollapse = FALSE,
                                        verbose = TRUE,
                                        colData = inputDataset)
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
  
  infile <- system.file("extdata/test_data.fastq_bismark.bismark.cov.gz",
                        package = 'bsseq')
  bismarkBSseq <- read.bismark(files = infile,
                               rmZeroCov = TRUE,
                               strandCollapse = FALSE,
                               verbose = TRUE)
                                 
  ######################################################################################################################## 
  ##### ----- edgeR ----- Filter by counts & DML analysis
  inputDataset <- read_tsv(file.path(dataDir, "input_dataset.tsv"))
  dataset <- read_tsv(file.path(dataDir, "dataset.tsv"))
  input <- dataset
  sampleNames <- dataset$Name
  # ordering
  orderChromosomes <- order(coverage$genes$Chr, coverage$genes$Locus)
  coverage <- coverage[orderChromosomes,]
  
  # coverage$samples$group <- factor(input_dataset$Population)
  coverage$samples$group <- factor(input$Treatment)
  
  coverage <- coverage[keep,, keep.lib.sizes=FALSE] # large DEGList
  
  # We now annotate the CpG loci with the identity of the nearest gene.
  # We search for the gene transcriptional start site (TSS) closest to each our CpGs:
  TSS <- nearestTSS(coverage$genes$Chr, coverage$genes$Locus, species="Mm")
  
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

  cmd <- paste("cnvkit.py batch", input$getColumn("BAM"),
               "--targets", input$getColumn("BedGraph"), 
               # "--annotate", input$getColumn("refBuild"), # optional; Format: UCSC, refFlat.txt or ensFlat.txt file (preferred), or BED, interval list, GFF, or similar.
               "--fasta", file.path(input$getColumn("refBuild"), "Genes/"), # required
               "--access", data/access-5kb-mappable.hg19.bed, # difference to other bed graph?
               "--output-reference reference.cnn",
               "--output-dir" output_dir,
               # --diagram
               # --scatter
               )
  result <- ezSystem(cmd)
  gc()
  
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
  rmarkdown::render(
    input = "DNAme.Rmd", envir = new.env(),
    output_dir = ".", output_file = htmlFile, quiet = TRUE
  )
  
  html_files <- c("00index.html",  "banner.png",  "fgcz.css",  "fgcz_header.html")
  
  return("Success")
}

### from deepak, to read bismark files more efficiently
# note: capital R instead of lower case
ReadBismark2DGE <- function (files, sample.names = NULL, readr = FALSE, data.table = FALSE, nParallel = 8, verbose = TRUE) {
  nsamples <- length(files)
  if (is.null(sample.names)) 
    sample.names <- removeExt(removeExt(removeExt(files)))
  if (readr) {
    OK <- requireNamespace("readr", quietly = TRUE)
    if (!OK) 
      stop("readr package required but not installed (or can't be loaded)")
  }
  if (data.table) {
    OK <- requireNamespace("data.table", quietly = TRUE)
    if (!OK)
      stop("data.table package required but not installed (or can't be loaded)")
  }
  CountList <- list()
  ChrRleList <- list()
  LocusList <- list()
  ChrNames <- c()
  MaxLocus <- 1L
  for (i in 1:nsamples) {
    if (verbose) 
      cat("Reading", files[i], "\n")
    if (readr) {
      x <- as.data.frame(suppressWarnings(readr::read_tsv(files[i], 
                                                          col_names = FALSE, col_types = "ci__ii", progress = FALSE)))
    }else if(data.table){
      x <- data.table::fread(input = files[i], sep = "\t", header = T, stringsAsFactors = F, check.names = F,
                             data.table = F, nThread = nParallel, drop = 3:4)
    }else{
      x <- read.delim(files[i], header = FALSE, colClasses = c("character", "integer", "NULL", "NULL", "integer", "integer"))
    }
    ChrRleList[[i]] <- rle(x[, 1])
    LocusList[[i]] <- x[, 2]
    CountList[[i]] <- as.matrix(x[, 3:4])
    ChrNames <- unique(c(ChrNames, ChrRleList[[i]]$values))
  }
  if (verbose) 
    cat("Hashing ...\n")
  for (i in 1:nsamples) ChrRleList[[i]]$values <- match(ChrRleList[[i]]$values, 
                                                        ChrNames)
  HashBase <- length(ChrNames) + 1L
  HashList <- list()
  HashUnique <- c()
  for (i in 1:nsamples) {
    HashList[[i]] <- inverse.rle(ChrRleList[[i]])/HashBase + 
      LocusList[[i]]
    HashUnique <- unique(c(HashUnique, HashList[[i]]))
  }
  if (verbose) 
    cat("Collating counts ...\n")
  counts <- matrix(0L, length(HashUnique), nsamples * 2L)
  j <- 1:2
  for (i in 1:nsamples) {
    m <- match(HashList[[i]], HashUnique)
    counts[m, j] <- CountList[[i]]
    j <- j + 2L
  }
  Locus <- as.integer(HashUnique)
  Chr <- as.integer((HashUnique - Locus) * HashBase + 0.5)
  attr(Chr, "levels") <- ChrNames
  class(Chr) <- "factor"
  Sample2 <- rep(sample.names, each = 2L)
  Methylation <- rep.int(c("Me", "Un"), nsamples)
  colnames(counts) <- paste(Sample2, Methylation, sep = "-")
  y <- DGEList(counts, genes = data.frame(Chr, Locus))
  row.names(y) <- paste(Chr, Locus, sep = "-")
  invisible(gc())
  return(y)
}
