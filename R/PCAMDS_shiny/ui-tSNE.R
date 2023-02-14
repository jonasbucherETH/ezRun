tabItem(
  tabName = "tab-tSNE",
  fluidRow( ### NOTE: 1 row has width = 12
    
    column(
      width = 2, 
      box(
        title = "t-SNE Parameters",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        # inputs
        # actionButton(
        #   inputId = "goButton",
        #   label = "Apply parameters",
        #   icon = NULL
        # ),
        numericInput(
          inputId = "perplexityTSNE",
          label = "Perplexity parameter",
          min = 5, max = 5,
          value = 5, step = 5
        ),
        numericInput(
          inputId = "max_iterTSNE",
          label = "Maximum number of iterations",
          min = 500, max = 3000,
          value = 1000, step = 500
        ),        
        numericInput(
          inputId = "thetaTSNE",
          label = "theta: Speed/accuracy trade-off. Set to 0.0 for exact TSNE",
          min = 0, max = 1,
          value = 0.5, step = 0.1
        ),        
        numericInput(
          inputId = "etaTSNE",
          label = paste0(expression(eta), ": learning rate"),
          min = 50, max = 800,
          value = 200, step = 50
        ),
        br(), br(),
        actionBttn(
          inputId = "paramButtonTSNE",
          label = "Apply parameters",
          icon = icon("check"),
          style = "simple",
          color = "primary",
          size = "md",
          block = FALSE,
          no_outline = TRUE
        ),
        # selectInput(
        #   inputId = "removeSamplesTSNE",
        #   label = "Select samples to exclude",
        #   choices = "",
        #   selected = ""
        # ),
        
      ) # close box 
    ), # close t-SNE parameters column
    
    column(
      width = 7,
      box(
        title = "TSNE Plots",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        downloadButton(
          outputId = "downloadTSNE",
          label = "Download t-SNE Plot (PDF)"
        ),
        br(), br(),
        plotOutput(
          outputId = "tsneStatic",
          height = "80vh",
          width = "100%",
          inline = F
        )
      ) # close box 
    ), # close plot(s) column
    
    column(
      width = 3,
      box(
        title = "Display options",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        collapsible = TRUE,
        collapsed = FALSE,
        
        
        actionBttn(
          inputId = "displayButtonTSNE",
          label = "Apply display options",
          icon = icon("text-size"),
          style = "simple",
          color = "primary",
          size = "md",
          block = FALSE,
          no_outline = TRUE
        ),
        br(), br(),
        # actionButton("reset", "Default colours", icon = icon("undo"))
        
        #### ----------------------------
        selectizeInput(
          inputId = "select",
          "Select:",
          choices = "",
          selected = "",
          multiple = TRUE
        ),
        
        # selectizeInput("select", "Select:", 
        #                choices = as.list(unique(groupingVariables$Population)), 
        #                # selected = "X1", 
        #                multiple = TRUE),
        
        uiOutput('myPanel'),
        
        #### ---------------------------
        
        checkboxInput(
          inputId = "sampleLabelsTSNE",
          label = "Display sample labels",
          value = FALSE
        ),
        
        ### no choices, selected = "" as default
        selectInput(
          inputId = "colorGroupTSNE",
          label = "Select groups to color by",
          choices = "",
          selected = ""
        ),
        selectInput(
          inputId = "shapeGroupTSNE",
          label = "Select groups for shapes",
          choices = "",
          selected = ""
        ),
        # selectInput(
        #   inputId = "pickFactor1TSNE",
        #   label = "Select Dimension for x-axis",
        #   choices = "X1",
        #   selected = "X1"
        # ),
        # selectInput(
        #   inputId = "pickFactor2TSNE",
        #   label = "Select Dimension for y-axis",
        #   choices = "X2",
        #   selected = "X2"
        # ),

        # sliderInput("TSNEPlotWidth", "Width of plot", min = 100, max = 2000, value = 800, step = 10),
        # sliderInput("TSNEPlotHeight", "Height of plot", min = 100, max = 2000, value = 600, step = 10),
        
        numericInput(
          inputId = "pointSizeTSNE",
          label = "Sizes of points in TSNE plot", min = 1, max = 6,
          value = 3, step = 0.5
        ),
        numericInput(
          inputId = "textSizeTSNE",
          label = "Figure Font Size", min = 4, max = 30,
          value = 12, step = 0.5
        ),

      ) # close box
    ), # close column
    
  
  ) # close fluidRow
) # close tabItem
