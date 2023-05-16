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
  # library("genekitr") # transId: gene id conversion
  library("GenomicFeatures")
  library("biomartr") # for getGO
  
  # #setwdNew(basename(output$getColumn("Report")))
  dataset <- input$meta
  ans4Report <- list() # a list of results for rmarkdown report
  ans4Report[["dataset"]] <- dataset
  
  cores <- detectCores(all.tests = FALSE, logical = TRUE) - 1
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add = TRUE)
  
  # output_dir <- basename(output$getColumn("Report"))
  # dataDir<- "/home/jobucher/data/dmrseq"
  dataDirSave <- "/home/jobucher/data/great"
  # dmRegionsFilePath <- file.path(dataDir, "dmRegions.rds")
  # significantRegionsFilePath <- file.path(dataDir, "significantRegions.rds")
  # 
  # dmRegions <- readRDS(dmRegionsFilePath)
  # significantRegions <- readRDS(significantRegionsFilePath)
  nRegions <- 10000
  # dmRegions <- randomRegionsFromBioMartGenome(param$biomart_dataset, nr = nRegions)
  
  print(param$biomart_selection)
  
  dmRegions <- randomRegionsFromBioMartGenome(param$biomart_selection, nr = nRegions)
  # dmRegions <- randomRegionsFromBioMartGenome("mmusculus_gene_ensembl", nr = nRegions)
  randomSubset <- sample(nRegions, nRegions/10)
  significantRegions <- dmRegions[randomSubset]
  
  # tableBiomart <- readRDS(system.file("extdata", "all_supported_organisms.rds", package = "BioMartGOGeneSets"))
  # tableBiomart$genesets = paste0("BP (", tableBiomart$n_bp_genesets, "), CC (", tableBiomart$n_cc_genesets, "), MF (", tableBiomart$n_mf_genesets, ")")
  # colnames(tableBiomart)[colnames(tableBiomart) == "n_gene"] = "genes"
  # df_genes = tableBiomart[tableBiomart$mart == "genes_mart", c("dataset", "name", "version", "taxon_id", "genbank_accession", "genesets"), drop = FALSE]
  # df_plants = tableBiomart[tableBiomart$mart == "plants_mart", c("dataset", "name", "version", "taxon_id", "genbank_accession", "genesets"), drop = FALSE]
  # df_metazoa = tableBiomart[tableBiomart$mart == "metazoa_mart", c("dataset", "name", "version", "taxon_id", "genbank_accession", "genesets"), drop = FALSE]
  # df_fungi = tableBiomart[tableBiomart$mart == "fungi_mart", c("dataset", "name", "version", "taxon_id", "genbank_accession", "genesets"), drop = FALSE]
  # df_protists = tableBiomart[tableBiomart$mart == "protists_mart", c("dataset", "name", "version", "taxon_id", "genbank_accession", "genesets"), drop = FALSE]
  # 
  # saveRDS(tableBiomart, file = "tableBiomart.rds")
  
  # write_csv(tableBiomart, "/srv/GT/analysis/jonas/jonas_test_sushi_20221115/master/selector_tests/biomart.csv")

  getGeneSets = function(dataset, ontology) {
    BioMartGOGeneSets::getBioMartGOGeneSets(dataset, ontology)
  }
  # geneSetsBP <- getGeneSetsFromBioMart("athaliana_eg_gene", "BP")
  

  geneSetsBP <- getGeneSets(param$biomart_selection, "BP")
  geneSetsCC <- getGeneSets(param$biomart_selection, "CC")
  geneSetsMF <- getGeneSets(param$biomart_selection, "MF")
  geneSetsAll <- c("BP" = geneSetsBP, "CC" = geneSetsCC, "MF" = geneSetsMF)
  
  ## Reactome pathways
  # if(param$biomart_dataset=="athaliana_eg_gene") {
  if(param$biomart_selection=="athaliana_eg_gene") {
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
    
    # ah <- AnnotationHub()
    # ensdb <- query(ah, c("EnsDb"))
    # ensdb <- ensdb[ensdb$species==""]
    # 
    # ensdb <- query(ah, c("EnsDb", param$species))
    # ensdb <- query(ah, c("EnsDb", "athaliana_eg_gene"))
    # # taxonomyid, genome, description, coordinate_1_based, maintainer, rdatadateadded, preparerclass, tags, rdatapath, sourceurl, sourcetype 
    # head(ensdb$taxonomyid)
    # ensdb <- rev(ensdb) # reverse (latest versions come first)
    # id <- ensdb$ah_id[1]
    # gs <- genes(ensdb[[id]], columns = c("tx_id", "gene_id", "gene_biotype", "symbol"))
    # extendedTSS <- extendTSS(gs)

    greatResult_RE <- great(gr = significantRegions, gene_sets = reactome_pathways, extended_tss = extendedTSS,
                          background = dmRegions, cores = cores)
    
    # nRegions <- 10000
    # # dmRegions <- randomRegionsFromBioMartGenome(param$biomart_dataset, nr = nRegions)
    # dmRegions <- randomRegionsFromBioMartGenome("athaliana_eg_gene", nr = nRegions)
    # # dmRegions <- randomRegionsFromBioMartGenome("mmusculus_gene_ensembl", nr = nRegions)
    # randomSubset <- sample(nRegions, nRegions/10)
    # significantRegions <- dmRegions[randomSubset]
    greatResult_RE <- great(gr = significantRegions, gene_sets = reactome_pathways, tss_source = "TxDb.Athaliana.BioMart.plantsmart51",
                            background = dmRegions, cores = 5)
    greatResult_KE <- great(gr = significantRegions, gene_sets = kegg_pathways, tss_source = "TxDb.Athaliana.BioMart.plantsmart51",
                            background = dmRegions, cores = 5)
    
    enrichmentTable_RE <- getEnrichmentTable(greatResult_RE)
    enrichmentTable_KE <- getEnrichmentTable(greatResult_KE)
    
    setwd("/home/jobucher/data/great/ath")
    saveRDS(greatResult_RE, file = paste0("greatResult_RE", ".rds"))
    saveRDS(greatResult_KE, file = paste0("greatResult_KE", ".rds"))
    saveRDS(enrichmentTable_RE, file = paste0("enrichmentTable_RE", ".rds"))
    saveRDS(enrichmentTable_KE, file = paste0("enrichmentTable_KE", ".rds"))

  }
  
  # ah <- AnnotationHub()
  # ensdb <- query(ah, c("EnsDb", param$species))
  # ensdb <- rev(ensdb) # reverse (latest versions come first)
  # id <- ensdb$ah_id[1]
  # gs <- genes(ensdb[[id]], columns = c("tx_id", "gene_id", "gene_biotype", "symbol"))
  # extendedTSS <- extendTSS(gs)
  # geneSetCollectionsAll <- substr(names(geneSetsAll), 1, 2) # BP, CC, ...
  # goTermsAll <- substr(names(geneSetsAll), 4, 1000000L)
  # names(geneSetsAll) <- goTermsAll
  # greatResult <- great(gr = significantRegions, gene_sets = geneSetsAll, extended_tss = extendedTSS,
  #                      background = dmRegions, cores = 5)
  
  greatResult_BP <- great(gr = significantRegions, gene_sets = "BP", biomart_dataset = param$biomart_selection,
                       background = dmRegions, cores = cores)
  greatResult_CC <- great(gr = significantRegions, gene_sets = "CC", biomart_dataset = param$biomart_selection,
                          background = dmRegions, cores = cores)
  greatResult_MF <- great(gr = significantRegions, gene_sets = "MF", biomart_dataset = param$biomart_selection,
                          background = dmRegions, cores = cores)
  # greatResult_BP <- great(gr = significantRegions, gene_sets = "BP", biomart_dataset = "mmusculus_gene_ensembl",
  #                         background = dmRegions, cores = 5)
  # greatResult_CC <- great(gr = significantRegions, gene_sets = "CC", biomart_dataset = "mmusculus_gene_ensembl",
  #                         background = dmRegions, cores = 5)
  # greatResult_MF <- great(gr = significantRegions, gene_sets = "MF", biomart_dataset = "mmusculus_gene_ensembl",
  #                         background = dmRegions, cores = 5)
  
  enrichmentTable_BP <- getEnrichmentTable(greatResult_BP)
  enrichmentTable_CC <- getEnrichmentTable(greatResult_CC)
  enrichmentTable_MF <- getEnrichmentTable(greatResult_MF)

#   regionGeneAssociations <- getRegionGeneAssociations(greatResult, term_id = NULL, by_middle_points = FALSE,
#                                                       use_symbols = TRUE)
  
  # set.seed(123)
  # gr = randomRegions(nr = 1000, genome = "hg19")
  # gres <- great(gr, "GO:MF", "hg19")
  # 
  # 
  # greatResult <- gres
 
  # setwdNew("/home/jobucher/data/great/mm")
  # on.exit(setwd(cwd), add = TRUE)
  setwd("/home/jobucher/data/great/mm")
  saveRDS(greatResult_BP, file = paste0("greatResult_BP", ".rds"))
  saveRDS(greatResult_CC, file = paste0("greatResult_CC", ".rds"))
  saveRDS(greatResult_MF, file = paste0("greatResult_MF", ".rds"))
  saveRDS(enrichmentTable_BP, file = paste0("enrichmentTable_BP", ".rds"))
  saveRDS(enrichmentTable_CC, file = paste0("enrichmentTable_CC", ".rds"))
  saveRDS(enrichmentTable_MF, file = paste0("enrichmentTable_MF", ".rds"))
  
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

