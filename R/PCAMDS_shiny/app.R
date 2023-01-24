library("shiny")
library("shinymanager")
library("shinyjs")
library("tidyverse")
library("ggpubr")
library("plotly")
library("DESeq2")
library("RColorBrewer")
library("ComplexHeatmap")
library("tximport")
library("vsn")
library("clusterProfiler")
library("DT")
library("colourpicker")
library("writexl")
library("circlize")
library("GO.db")
library("shinydashboard")
library("shinyBS")
library("pixiedust")
library("ezRun")
library("kableExtra")
library("ggrepel")
library("gplots")
library("sortable")
library("waiter")
library("shinycssloaders")
library("ggprism")
library("ggbeeswarm")
library("rstatix")
library("gridExtra")
library("shinytitle")
library("shinylogs")

# JB libraries
library("magrittr")
library("adegenet")
library("ade4")
library("gdsfmt")
library("SNPRelate")
library("ggplot2")
library("shinythemes")
library("shinytest")
library("vcfR")
library("testthat")

# console.error = function () {
#   require("system").stderr.write(Array.prototype.join.call(arguments, ' ') + '\n');
# };

reactiveConsole(TRUE)

# For secure login: 
# library(digest)
# digest("password1", algo = "md5")
# credentials <- data.frame(
#   user = c("user1", "user2"),
#   password = c(
#     "password1", "password2"),
#   admin = c(FALSE, TRUE),
#   comment = "Login Page.",
#   stringsAsFactors = FALSE
# )

############### load data for testing
# vcf <- read.vcfR("~/sushi_project_JB/data/test_vcf_dataset/ragi_highcov_sa0001_1k.vcf.gz")
# genind <- vcfR2genind(vcf)
# grouping_vars <- read.delim("~/sushi_project_JB/data/test_vcf_dataset/populations.txt")
# # pop(genind) <- populations_txt$Population
# 
# X <- scaleGen(genind, NA.method="mean")

# pca <- dudi.pca(X, center = TRUE, scale = TRUE, scan = FALSE)

# mds <- read.csv("~/git/ezRun/output_data/plink_mds.mds", sep="")
# 
# ### load data
# # pca <- readRDS("PCA.rds")
# # grouping_vars <- readRDS("grouping_vars.rds")
# # mds <- read.csv("plink.mds", sep="")
# 
# 
# ### both testing & normal setup
# n_pcs <- pca$nf # number of (> 0) principal components
# 
# eig_sum <- sum(pca$eig)
# pca_varprop <- pca$eig/eig_sum
# pca_varprop <- pca_varprop[1:n_pcs]
# 
# pca_tab <- data.frame(grouping_vars, pca$li, stringsAsFactors = FALSE, row.names = NULL)
# 
# tab_varprop <- as.data.frame(t(pca_varprop), stringsAsFactors = FALSE)
# PC_indeces <- seq(1+ncol(grouping_vars), ncol(tab))
# 
# for (i in 1:n_pcs){
#   colnames(pca_tab)[i+ncol(grouping_vars)] <- paste0("PC", i)
#   colnames(tab_varprop)[i] <- paste0("PC", i)
# }


# 
# n_dim <- ncol(mds)-ncol(grouping_vars)-1   # number of dimensions kept in mds
# mds_tab <- data.frame(grouping_vars, mds[, (ncol(grouping_vars)+2):ncol(mds)], stringsAsFactors = FALSE, row.names = NULL)

##############

# spinner <- tagList(
#   spin_chasing_dots(),
#   span("Loading stuff...", style="color:white;")
# )

ui = dashboardPage(
  skin = "blue",
  dashboardHeader(
    title = "Dimensionality Reduction" #,
    # tags$li(
    #   a(
    #     href = 'mailto:sequencing@fgcz.ethz.ch?subject=exploreDEG-shiny-app-feedback', 
    #     "Request Features/Report Bugs"), 
    #   class = "dropdown"
    # ),
    # tags$li(
    #   a(href = 'http://www.fgcz.ch', 
    #     target = "_blank",
    #     img(src = 'fgcz_logo.png', title = "FGCZ", height = "30px"),
    #     style = "padding-top:10px; padding-bottom:5px;"),
    #   class = "dropdown"),
    # tags$li(
    #   a(href = 'http://www.ethz.ch/en.html',
    #     target = "_blank",
    #     img(src = 'eth_logo.png', title = "FGCZ", height = "22px"),
    #     style = "padding-top:13px; padding-bottom:10px;"),
    #   class = "dropdown"),
    # tags$li(
    #   a(href = 'http://www.uzh.ch/en.html',
    #     target = "_blank",
    #     img(src = 'University_of_Zurich_Logo.png', title = "FGCZ", height = "30px"),
    #     style = "padding-top:10px; padding-bottom:5px;"),
    #   class = "dropdown")
  ),
  dashboardSidebar(
    shinyjs::useShinyjs(),
    sidebarMenu(
      id = "tabs",
      menuItem(
        text = "PCA Plots",
        tabName = "tab-PCA",
        icon = icon("meteor")
      )
    )
  ), 
  dashboardBody(
    # use_tracking(),
    # tags$head(tags$link(rel = "shortcut icon", href = "sushi.png")),
    # tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "main.css")),
    # use_waiter(),
    tabItems(
      source("ui-PCA.R", local = TRUE)$value
      # source("~/git/ezRun/R/PCAMDS_shiny/ui-PCA.R", local = F)$value
      
    )
  )
)

# For secure login:
# ui <- secure_app(ui)

cat("middle")


server = function(input, output, session) {
  cat("reached server")
  # Server login if required:
  # res_auth <- secure_server(
  #   check_credentials = check_credentials(credentials)
  # )
  # 
  # output$auth_output <- renderPrint({
  #   reactiveValuesToList(res_auth)
  # })
  # track_usage(storage_mode = store_rds(path = "/scratch/shiny_logs/"), app_name="PCAMDS_shiny")
  
  source("server-inputData.R", local = TRUE)
  source("server-PCA.R", local = TRUE)
  # source("~/git/ezRun/R/PCAMDS_shiny/server-inputData.R", local = F)
  # source("~/git/ezRun/R/PCAMDS_shiny/server-PCA.R", local = F)

  
}

shinyApp(ui = ui, server = server)


