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
  library("BioMartGOGeneSets")
  library("parallel")
  library("AnnotationHub")
  library("genekitr") # transId: gene id conversion
  
  # #setwdNew(basename(output$getColumn("Report")))
  dataset <- input$meta
  ans4Report <- list() # a list of results for rmarkdown report
  ans4Report[["dataset"]] <- dataset
  
  cores <- detectCores(all.tests = FALSE, logical = TRUE) - 1
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add = TRUE)
  
  # output_dir <- basename(output$getColumn("Report"))
  # dataDir <- "/home/jobucher/data/dmrseq"
  # dmRegionsFilePath <- file.path(dataDir, "dmRegions.rds")
  # significantRegionsFilePath <- file.path(dataDir, "significantRegions.rds")
  # 
  # dmRegions <- readRDS(dmRegionsFilePath)
  # significantRegions <- readRDS(significantRegionsFilePath)
  nRegions <- 10000
  dmRegions <- randomRegionsFromBioMartGenome(param$biomart_dataset, nr = nRegions)
  # dmRegions <- randomRegionsFromBioMartGenome("mmusculus_gene_ensembl", nr = nRegions)
  randomSubset <- sample(nRegions, nRegions/10)
  significantRegions <- dmRegions[randomSubset]
  
  tableBiomart <- readRDS(system.file("extdata", "all_supported_organisms.rds", package = "BioMartGOGeneSets"))
  tableBiomart$genesets = paste0("BP (", tableBiomart$n_bp_genesets, "), CC (", tableBiomart$n_cc_genesets, "), MF (", tableBiomart$n_mf_genesets, ")")
  colnames(tableBiomart)[colnames(tableBiomart) == "n_gene"] = "genes"
  
  getGeneSetsFromBioMart = function(dataset, ontology) {
    BioMartGOGeneSets::getBioMartGOGeneSets(dataset, ontology)
  }
  # geneSetsBP <- getGeneSetsFromBioMart("athaliana_eg_gene", "BP")
  
  geneSetsBP <- getGeneSetsFromBioMart(param$biomart_dataset, "BP")
  geneSetsCC <- getGeneSetsFromBioMart(param$biomart_dataset, "CC")
  geneSetsMF <- getGeneSetsFromBioMart(param$biomart_dataset, "MF")
  geneSetsAll <- c("BP" = geneSetsBP, "CC" = geneSetsCC, "MF" = geneSetsMF)
  
  # geneSetsBP <- getGeneSetsFromBioMart("mmusculus_gene_ensembl", "BP")
  # geneSetsCC <- getGeneSetsFromBioMart("mmusculus_gene_ensembl", "CC")
  # geneSetsMF <- getGeneSetsFromBioMart("mmusculus_gene_ensembl", "MF")
  # geneSetsAll <- c("BP" = geneSetsBP, "CC" = geneSetsCC, "MF" = geneSetsMF)
  
  ## Reactome pathways
  if(param$biomart_dataset=="athaliana_eg_gene") {
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
    
    geneSetsAthaliana <- c("RP" = reactome_pathways, "KP" = kegg_pathways)
    geneSetsAll <- c(geneSetsAll, geneSetsAthaliana)

  }
  
  ah <- AnnotationHub()
  # hub <- subset(hub, hub$species=='Drosophila melanogaster')
  ensdb <- query(ah, c("EnsDb", param$species)) 
  # ensdb <- query(ah, c("EnsDb", "Mus musculus")) 
  # ensdb <- query(ah, c("GRCm38", "EnsDb")) # 102 = latest version
  ensdb <- rev(ensdb) # reverse (latest versions come first)
  # id <- ensdb$ah_id[grep(pattern = param$txdb_dataset, x = ensdb$title)]
  id <- ensdb$ah_id[1]
  # id <- ensdb$ah_id[grep(pattern = "Mus musculus", x = ensdb$title)]

  gs <- genes(ensdb[[id]], columns = c("tx_id", "gene_id", "gene_biotype", "symbol"))
  # gene = gene[seqnames(gene) %in% paste0("chr", c(1:50, "X", "Y"))]
  # gl = seqlengths(gene)[paste0("chr", c(1:22, "X", "Y"))]  # restrict to normal chromosomes
  # gl <- seqlengths(gene) # restrict to normal chromosomes

  # speciesNameBiomart <- tableBiomart[tableBiomart$dataset == param$biomart_dataset, ]
  # speciesNameBiomart <- tableBiomart[tableBiomart$dataset == "mmusculus_gene_ensembl", ]
  # speciesNameBiomart <- speciesNameBiomart$name
  # speciesToTransId <- c("hsapiens_gene_ensembl" ,"mmusculus_gene_ensembl")
  # if(param$biomart_dataset %in% speciesToTransId) {
  #   org <- substr(param$biomart_dataset, 1, 2) # turns it into hs / mm
  #   newGeneId <- transId(gene$gene_id, transTo = "ensembl", org = org, keepNA = TRUE, unique = TRUE) # df with 2 columns (input_id, ensembl)
  #   gene$gene_id <- newGeneId$ensembl
  #   gene[is.na(gene$gene_id), ]
  # }
  
  extendedTSS <- extendTSS(gs)
  
  greatResult <- great(gr = significantRegions, gene_sets = geneSetsAll, extended_tss = extendedTSS,
                       background = dmRegions)

  # set.seed(123)
  # gr = randomRegions(nr = 1000, genome = "hg19")
  # gres <- great(gr, "GO:MF", "hg19")
  # 
  # 
  # greatResult <- gres
  
  enrichmentTable <- getEnrichmentTable(greatResult, min_region_hits = 5)
  enrichmentTable$gene_set_collection <- substr(enrichmentTable$id, 1, 2)
  enrichmentTable$gene_set_id <- substr(enrichmentTable$id, 4, 1000000L)
  
  # regionGeneAssociations <- getRegionGeneAssociations(greatResult, term_id = NULL, by_middle_points = FALSE,
                            # use_symbols = TRUE)
  
  saveRDS(greatResult, file="greatResult.rds")
  saveRDS(enrichmentTable, file="enrichmentTable.rds")
  saveRDS(geneSetsAll, file="geneSetsAll.rds")

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

