### https://www.bioconductor.org/packages/release/bioc/vignettes/rGREAT/inst/doc/local-GREAT.html

# sampleNames <- c("imr90.r1", "imr90.r2", "h1.r1", "h1.r2")
# cellType <- c("imr90", "imr90", "h1", "h1")

library("dmrseq")
library("tidyverse")
library("rGREAT")

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

great_analysis <- function() {
  set.seed(123)
  gr = randomRegions(nr = 1000, genome = "mm10")
  # great = main function; can use built-in annot. or self-provided ones
  # res = great(gr, gene_sets, tss_source, biomart_dataset = NULL,
                    # min_gene_set_size = 5, mode = "basalPlusExt", basal_upstream = 5000,
                    # basal_downstream = 1000, extension = 1000000,
                    # extended_tss = NULL, background = NULL, exclude = "gap",
                    # cores = 1, verbose = great_opt$verbose
          # )
  
  # If you provide a self-defined gene_sets or extended_tss,
  # you need to make sure they two have the same gene ID types.
  
  # geneSets <- read_gmt(url("http://dsigdb.tanlab.org/Downloads/D2_LINCS.gmt"), 
  #               from = "SYMBOL", to = "ENTREZ", orgdb = "org.Hs.eg.db")
  
  geneSets <- read_gmt(refBuild, 
                       from = "SYMBOL", to = "ENTREZ", orgdb = NULL)
  df <- read.table(url("https://jokergoo.github.io/rGREAT_suppl/data/GREATv4.genes.hg19.tsv"))
  # note there must be a 'gene_id' column
  tss <- GRanges(seqnames = df[, 2], ranges = IRanges(df[, 3], df[, 3]), 
                strand = df[, 4], gene_id = df[, 5])
  # extendTSS accepts GRanges object of gene or TSS, returns new GRanges object
  tssSource <- extendTSS(tss, genome = "hg19", gene_id_type = "SYMBOL")
  
  
  res <- great(gr, gene_sets = geneSets, tss_source = tssSource)
  res
  tb <- getEnrichmentTable(res)
  head(tb)
  return(res)

  
}

gr = randomRegions(nr = 1000, genome = "hg19")
res <- great(gr, "GO:BP", "txdb:mm10")

plotVolcano(res)
plotRegionGeneAssociations(res)
getRegionGeneAssociations(res)
plotRegionGeneAssociations(res, term_id = "HALLMARK_APOPTOSIS")
shinyReport(res)

# values in res@table
colnames(res@table)
# [1] "id"                    "genome_fraction"       "observed_region_hits"  "fold_enrichment"      
# [5] "p_value"               "p_adjust"              "mean_tss_dist"         "observed_gene_hits"   
# [9] "gene_set_size"         "fold_enrichment_hyper" "p_value_hyper"         "p_adjust_hyper" 

### Global region-gene associations
res@table$observed_gene_hits
res@table$observed_region_hits
res@table$fold_enrichment_hyper
tbl = getEnrichmentTables(res)

# plotRegionGeneAssociations(object, ontology = NULL, term_id = NULL, which_plot = 1:3,
#                            request_interval = 10, max_tries = 100, verbose = great_opt$verbose)
plotRegionGeneAssociations(res, term_id = NULL, which_plot = 1:3)

gg_associations <- function(resultGreat) {
  p1 <- resultGreat %>%
    ggplot(aes(x=observed_gene_hits, y=-log10(p_value))) +
    geom_bar() +
    labs(title = "Number of associated genes per region",
         xlab = "Number of associated genes per region")
    theme_bw()
  p
}

### volcano plot
# https://github.com/jokergoo/rGREAT/blob/master/R/volcano.R

colnames(res@table)
# [1] "id"                    "genome_fraction"       "observed_region_hits"  "fold_enrichment"      
# [5] "p_value"               "p_adjust"              "mean_tss_dist"         "observed_gene_hits"   
# [9] "gene_set_size"         "fold_enrichment_hyper" "p_value_hyper"         "p_adjust_hyper" 

gg_volcanoPlot <- function(resultGreat, xAxisVolcano = "Fold enrichment", yAxisVolcano = "P-value") {
  if(xAxisVolcano=="Fold enrichment") {
    x <- log2(resultGreat$fold_enrichment)
    xlab <- "log2 fold enrichment: log2(obs/exp)"
  }
  else {
    x <- resultGreat$z_score
    xlab <- "z-score: (obs-exp)/sd"
  }
  if(yAxisVolcano=="P-value") {
    y <- -log10(resultGreat$p_value)
    ylab <- "-log10(p_value)"
  }
  else {
    y <- -log10(resultGreat$p_adjust)
    ylab <- "-log10(p_adjust)"
  }
  p <- resultGreat %>%
    ggplot(aes(x=x, y=y)) +
    geom_point() +
    geom_hline(aes(yintercept=-log10(0.05)), col = "grey", lty = 2) +
    labs(xlab = xlab,
         ylab = ylab)
    theme_bw() 
  p
}

gg_volcanoPlot(res@table)

plot_volcano = function(
    observed_region_hits,
    genome_fraction,
    fold_enrichment,
    z_score,
    p_value,
    p_adjust,
    min_region_hits = 5, 
    x_values = "fold_enrichment",
    y_values = "p_value",
    main = NULL) {
  
  col_fun = colorRamp2(seq(min_region_hits, min_region_hits + min(50, max(observed_region_hits)), length = 11), rev(brewer.pal(11, "Spectral")))
  size_fun = function(x) x^0.3*2 + 0.1
  
  if(x_values == "fold_enrichment") {
    x = log2(fold_enrichment)
    xlab = "log2 fold enrichment: log2(obs/exp)"
  } else {
    x = z_score
    xlab = "z-score: (obs-exp)/sd"
  }
  
  
  if(y_values == "p_value") {
    y = -log10(p_value)
    ylab = "-log10(p_value)"
  } else {
    y = -log10(p_adjust)
    ylab = "-log10(p_adjust)"
  }
  y[is.infinite(y)] = max(y[is.finite(y)])
  
  l = observed_region_hits >= min_region_hits
  if(!any(l)) {
    plot(NULL, xlim = c(0, 1), ylim = c(0, 1), ann = FALSE, axes = FALSE)
    text(0.5, 0.5, qq("No term left with min_region_hits >= @{min_region_hits}"))
    return(invisible(NULL))
  }
  
  plot(x[l], y[l], pch = 16,
       col = col_fun(observed_region_hits[l]),
       cex = size_fun(genome_fraction[l]),
       xlab = xlab, ylab = ylab)
  if(is.null(main)) {
    title("Volcano plot")
  } else {
    title(main)
  }
  abline(h = -log10(0.05), col = "grey", lty = 2)
  text(par("usr")[2], -log10(0.05), qq("@{y_values} = 0.05"), adj = c(1.05, 1.5), cex = 0.8)
  if(x_values == "fold_enrichment") abline(v = 1, col = "grey", lty = 2)
  
  max_region_hits = min_region_hits + min(50, max(observed_region_hits))
  breaks = seq(min_region_hits, max_region_hits, length = 5)
  breaks = round(breaks)
  size = legend("topleft", pch = 16, cex = 0.8, legend = breaks, col = col_fun(breaks), 
                title = "Region hits", bty = "n")
  
  rg = range(genome_fraction[l])
  size_rg = size_fun(rg)
  breaks = ((seq(size_rg[1], size_rg[2], length = 5) - 0.1)/2)^(1/0.3)
  
  size = legend(x = size$rect$left, y = size$rect$top - size$rect$h*1.05,
                pch = 16, pt.cex = size_fun(breaks), cex = 0.8,
                legend = paste0(sprintf("%.2f", breaks*100), "%"), 
                col = "black", 
                title = "Genome fraction", bty = "n")
  
  text(size$rect$left, size$rect$top - size$rect$h, qq("Region hits >= @{min_region_hits}"), adj = c(-0.05, 2), cex = 0.8)
  
}
