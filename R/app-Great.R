###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodGreat <- function(input = NA, output = NA, param = NA,
                           htmlFile = "00index.html") {
  
  ### libraries
  library("rGREAT")
  library("KEGGREST")
  library("biomaRt")
  
  # #setwdNew(basename(output$getColumn("Report")))
  dataset <- input$meta
  ans4Report <- list() # a list of results for rmarkdown report
  ans4Report[["dataset"]] <- dataset
  
  # output_dir <- basename(output$getColumn("Report"))
  dataDir <- "~/data/dmrseq"
  dmRegionsFilePath <- file.path(dataDir, "dmRegions.rds")
  significantRegionsFilePath <- file.path(dataDir, "significantRegions.rds")
  
  dmRegions <- readRDS(dmRegionsFilePath)
  significantRegions <- readRDS(significantRegionsFilePath)
  
  ## Reactome pathways
  if(param$reactome_kegg & param$biomart_dataset=="athaliana_eg_gene") {
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
    
    geneSetsAthaliana <- c("reactome_pathways" = reactome_pathways, "kegg_pathways" = kegg_pathways)
    
    # gs = read_gmt(url("http://dsigdb.tanlab.org/Downloads/D2_LINCS.gmt"), 
    #               from = "SYMBOL", to = "ENTREZ", orgdb = "org.Hs.eg.db")
    
    ### TODO: add built-in gene set(s) to these
    greatResult <- great(gr = dmRegions, gene_sets = geneSetsAthaliana, tss_source = "TxDb.Athaliana.BioMart.plantsmart51",
                         min_gene_set_size = param$min_gene_set_size, mode = param$mode, basal_upstream = param$basal_upstream,
                         basal_downstream = param$basal_downstream, extension = param$extension,
                         background = dmRegions, exclude = param$exclude,
                         cores = param$cores)
  }
  
  else {
    # tssSource <- ""
    geneSets <- param$gene_sets
    
    greatResult <- great(gr = significantRegions, gene_sets = geneSets, biomart_dataset = biomart_dataset,
                         min_gene_set_size = param$min_gene_set_size, mode = param$mode, basal_upstream = param$basal_upstream,
                         basal_downstream = param$basal_downstream, extension = param$extension,
                         background = dmRegions, exclude = param$exclude, # gap regions for corresponding organism will be removed from the analysis
                         cores = param$cores)
  }
  
  # set.seed(123)
  # gr = randomRegions(nr = 1000, genome = "hg19")
  # gres <- great(gr, "GO:MF", "hg19")
    
  
  # greatResult <- gres
  enrichmentTable <- getEnrichmentTable(greatResult)
  # head(enrichmentTable)
  
  saveRDS(greatResult, file="greatResult.rds")
  saveRDS(enrichmentTable, file="enrichmentTable.rds")
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add = TRUE)
  
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

