#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

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
# library("shinylogs") # not yet installed

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


### load data for testing
vcf <- read.vcfR("~/sushi_project_JB/data/test_vcf_dataset/ragi_highcov_sa0001_1k.vcf.gz")
genind <- vcfR2genind(vcf)
grouping_vars <- read.delim("~/sushi_project_JB/data/test_vcf_dataset/populations.txt")
# pop(genind) <- populations_txt$Population

X <- scaleGen(genind, NA.method="mean")

# library(FactoMineR)
# library(Factoshiny)
# library(testthat)

# pca_f <- PCA(X)
# pca_shiny <- PCAshiny(pca_f)

pca <- dudi.pca(X, center = TRUE, scale = TRUE, scan = FALSE, nf = 5)
n_pcs <- pca$nf # number of (> 0) principal components
n_grouping <- ncol(grouping_vars)
eig_sum <- sum(pca$eig)
pca_varprop <- pca$eig/eig_sum
pca_varprop <- pca_varprop[1:n_pcs]
pca_tab <- data.frame(grouping_vars, pca$li, stringsAsFactors = FALSE, row.names = NULL)
tab_varprop <- as.data.frame(t(pca_varprop), stringsAsFactors = FALSE)
PC_indeces <- seq(1+n_grouping, ncol(pca_tab))
for (i in 1:n_pcs){
  colnames(pca_tab)[i+n_grouping] <- paste0("PC", i)
  colnames(tab_varprop)[i] <- paste0("PC", i)
}
# all PCs in array for selecting input
pc_list <- colnames(pca_tab)[-(1:n_grouping)]


# Define UI for application that draws a histogram
ui <- 

  
  tabItem(
    tabName = "tab-PCA",
    fluidRow(
      column(
        width = 3,
        box(
          title = "Display options",
          width = NULL,
          solidHeader = TRUE,
          status = "primary",
          collapsible = TRUE,
          collapsed = TRUE,
          
          checkboxInput("sample_labels", "Display sample labels", value = TRUE),
          
          ### grouping input (colors & shapes)
          # selectInput("color_by", "Select groups to color by", choices = c("", colnames(grouping_vars))),
          # selectInput("group_by", "Select groups for shapes",  choices = c("", colnames(grouping_vars))),
          # selectInput("color_by", "Select groups to color by", choices = colnames(grouping_vars)),
          # selectInput("group_by", "Select groups for shapes",  choices = colnames(grouping_vars)),
          
          # varSelectInput("color_group", "Select groups to color by", data = pca_tab[,-PC_indeces]),
          # varSelectInput("shape_group", "Select groups for shapes", data = pca_tab[,-PC_indeces]),
          # selectInput("color_group", "Select groups to color by", choices = list("Sample", "Population")),
          # selectInput("shape_group", "Select groups for shapes", choices = list("Sample", "Population"), selected = "Population"),
          
          selectInput(
            inputId = "color_group",
            label = "Select groups to color by",
            choices = "", 
            selected = ""
          ),
          selectInput(
            inputId = "shape_group",
            label = "Select groups for shapes",
            choices = "", 
            selected = ""
          ),
          
          
          ### PC input (axis)
          # selectInput(inputId = "pick_pc_x", label = "Select PC for x-axis", choices = list("PC1", "PC2", "PC3", "PC4")),
          # selectInput(inputId = "pick_pc_y", label = "Select PC for y-axis", choices = list("PC1", "PC2", "PC3", "PC4"), selected = "PC2")
          
          # selectInput(inputId = "pick_pc_x", label = "Select PC for x-axis", choices = colnames(pca_tab)[PC_indeces]),
          # selectInput(inputId = "pick_pc_y", label = "Select PC for y-axis", choices = colnames(pca_tab)[PC_indeces], selected = colnames(pca_tab)[PC_indeces[2]])
          
          # selectInput(
          #   inputId = "pick_pc_x",
          #   label = "Select PC for x-axis",
          #   choices = "PC1", 
          #   selected = "PC1"
          # ),
          # selectInput(
          #   inputId = "pick_pc_y",
          #   label = "Select PC for y-axis",
          #   choices = "PC2", 
          #   selected = "PC2"
          # ),
          
          selectInput(
            inputId = "pick_pc_x",
            label = "Select PC for x-axis",
            choices = list("PC1","PC2","PC3","PC4","PC5"), 
            selected = "PC1"
          ),
          selectInput(
            inputId = "pick_pc_y",
            label = "Select PC for y-axis",
            choices = list("PC1","PC2","PC3","PC4","PC5"), 
            selected = "PC2"
          ),
          
          # selectInput("pick_pc_x", "Select PC for x-axis", choices = colnames(pca_tab)[PC_indeces]),
          # selectInput("pick_pc_y", "Select PC for y-axis", choices = colnames(pca_tab)[PC_indeces], selected = "PC2")
        ) # close box
      ), # close column
      
      column(
        width = 8,
        box(
          title = "PCA Plots",
          width = NULL,
          solidHeader = TRUE,
          status = "primary",
          plotOutput(
            outputId = "pca_plot",
            inline = TRUE
          )
        ), # close box
        
        box(
          title = "PCA Table",
          width = NULL,
          solidHeader = TRUE,
          status = "primary",
          tableOutput(outputId = "pca_table")
        ) # close box
        ,
        box(
          title = "Test text",
          width = NULL,
          solidHeader = TRUE,
          status = "primary",
          textOutput(outputId = "test_text")
        ) # close box
        
      ) # close column
    ) # close fluidRow
  ) # close tabItem
  


# Define server logic required to draw a histogram
server <- function(input, output, session) {

  observe({
    # withProgress(message = "Generating PCA Plots. Please wait...", {
    # pca <- inputDataReactive()$pca
    n_pcs <- inputDataReactive()$n_pcs
    # pca_varprop <- inputDataReactive()$pca_varprop
    # pca_tab <- inputDataReactive()$pca_tab
    # tab_varprop <- inputDataReactive()$tab_varprop
    
    vcf <- inputDataReactive()$vcf
    genind <- inputDataReactive()$genind
    grouping_vars <- inputDataReactive()$grouping_vars
    X <- inputDataReactive()$X
    
    # cat(n_pcs)
    
    # observeEvent(inputDataReactive()$n_pcs, {
    #   print(paste0("Test: ", inputDataReactive()$n_pcs))
    # })
    
    # observeEvent(input()$pick_pc_x, {
    #   print(paste0("Test: ", input()$pick_pc_x))
    # })
    
    
    # n_grouping <- inputDataReactive()$n_grouping
    # pc_list <- inputDataReactive()$pc_list
    
    output$pca_table <- renderTable(pca_tab)
    
    output$test_text <- renderText({n_pcs})
    
    
    
    # observeEvent({
    #   updateSelectInput(session = session, inputId = "pick_pc_x", selected = "PC1")
    #   updateSelectInput(session = session, inputId = "pick_pc_y", selected = "PC2")
    # })
    
    # updateVarSelectInput(session = session, "pick_pc_x", data = pca_tab[, PC_indeces], selected = "PC1")
    # updateVarSelectInput(session = session, "pick_pc_x", data = pca_tab[, PC_indeces], selected = "PC2")
    
    # observeEvent({
    #   updateSelectInput(session = session, "pick_pc_x", choices = pc_list, selected = "PC1")
    #   updateSelectInput(session = session, "pick_pc_x", choices = pc_list, selected = "PC2") 
    # })
    
    # updateSelectInput(session = session, "pick_pc_x", choices = pc_list, selected = "PC1")
    # updateSelectInput(session = session, "pick_pc_x", choices = pc_list, selected = "PC2")
    
    
    # for testing
    # vcf <- read.vcfR("~/sushi_project_JB/data/test_vcf_dataset/ragi_highcov_sa0001_1k.vcf.gz")
    # genind <- vcfR2genind(vcf)
    # grouping_vars <- read.delim("~/sushi_project_JB/data/test_vcf_dataset/populations.txt")
    # # pop(genind) <- populations_txt$Population
    # 
    # X <- scaleGen(genind, NA.method="mean")
    
    observeEvent( # number 1
      {
        # input$pick_pc_x
        # input$pick_pc_y
        # input$color_group
        # input$shape_group
        input$sample_labels
      }, ignoreInit = TRUE,
      {
        vcf <- vcf
        genind <- genind
        grouping_vars <- grouping_vars
        X <- X
        
        # compute PCA
        pca <- dudi.pca(X, center = TRUE, scale = TRUE, scan = FALSE, nf = 5)
        
        n_pcs <- pca$nf # number of (> 0) principal components
        print(n_pcs)
        n_grouping <- ncol(grouping_vars)
        
        eig_sum <- sum(pca$eig)
        pca_varprop <- pca$eig/eig_sum
        pca_varprop <- pca_varprop[1:n_pcs]
        
        pca_tab <- data.frame(grouping_vars, pca$li, stringsAsFactors = FALSE, row.names = NULL)
        
        tab_varprop <- as.data.frame(t(pca_varprop), stringsAsFactors = FALSE)
        PC_indeces <- seq(1+n_grouping, ncol(tab))
        
        for (i in 1:n_pcs){
          colnames(pca_tab)[i+n_grouping] <- paste0("PC", i)
          colnames(tab_varprop)[i] <- paste0("PC", i)
        }
        
        # all PCs in array for selecting input
        pc_list <- colnames(pca_tab)[-(1:n_grouping)]
        
        updateSelectInput(session = session, "pick_pc_x", choices = pc_list, selected = "PC1")
        updateSelectInput(session = session, "pick_pc_x", choices = pc_list, selected = "PC2")
        
        observeEvent({ # number 2
          input$pick_pc_x
          input$pick_pc_y
          input$color_group
          input$shape_group
        }, {
          # { # aaa
          
          # pca_ggplot <- ggplot(data = pca_tab, aes(x = .data[[input$pick_pc_x]], y = .data[[input$pick_pc_y]])) +
          # pca_ggplot <- ggplot(data = pca_tab, aes_string(x = input$pick_pc_x, y = input$pick_pc_y)) +
          pca_ggplot <- pca_tab %>%
            ggplot(aes_string(x = input$pick_pc_x, y = input$pick_pc_y)) +
            geom_point() +
            theme_classic()
          
          output$pca_plot <- renderPlot({
            pca_ggplot
          })
          
          observeEvent(input()$pick_pc_x, {
            print(paste0("Test: ", input()$pick_pc_x))
          })
          
          # test_table <- xtable::xtable(pca_tab)
          
          # test_table <- matrix(data = NA, nrow = 2, ncol = 2, byrow = FALSE,
          #                      dimnames = NULL)
          # 
          # test_table <- xtable::xtable(test_table)
          # 
          # 
          # output$pca_table <- renderTable(test_table)
          # output$pca_table <- renderTable(pca_tab)
          
          
          # }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")} # close withProgress
        }) # close observeEvent number 2
        # } # aaa
        
        
      }) # close observeEvent number 1
    # }) # close withProgress
  }) # close observe
  
  
  
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
