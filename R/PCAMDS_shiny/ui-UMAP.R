tabItem(
  tabName = "tab-UMAP",
  fluidRow( ### NOTE: 1 row has width = 12
    
    column(
      width = 2, 
      box(
        title = "Parameters",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        # inputs
        actionButton(
          inputId = "goButtonUMAP",
          label = "Apply parameters",
          icon = NULL
        ),
        br(), br(),

        selectInput(
          inputId = "excludeSamplesUMAP",
          label = "Select samples to exclude",
          choices = "",
          selected = ""
        
        ),
      ) # close box 
    ), # close  column
    
    column(
      width = 7,
      box(
        title = "UMAP Plots",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        downloadButton(
          outputId = "downloadUMAP",
          label = "Download UMAP Plot (PDF)"
        ),
        br(), br(),
        plotOutput(
          outputId = "UMAPStatic",
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
        checkboxInput(
          inputId = "displayTitleUMAP",
          label = "Display Title",
          value = FALSE
        ),
        textInput(
          inputId = "UMAPTitle",
          label = "Title of UMAP plot",
          value = ""
        ),
        
        checkboxInput(
          inputId = "sampleLabelsUMAP",
          label = "Display sample labels",
          value = FALSE
        ),
        
        ### no choices, selected = "" as default
        selectInput(
          inputId = "colorGroupUMAP",
          label = "Select groups to color by",
          choices = "",
          selected = ""
        ),
        selectInput(
          inputId = "shapeGroupUMAP",
          label = "Select groups for shapes",
          choices = "",
          selected = ""
        ),
        # selectInput(
        #   inputId = "pickFactor1UMAP",
        #   label = "Select Dimension for x-axis",
        #   choices = "X1",
        #   selected = "X1"
        # ),
        # selectInput(
        #   inputId = "pickFactor2UMAP",
        #   label = "Select Dimension for y-axis",
        #   choices = "X2",
        #   selected = "X2"
        # ),
        
        # sliderInput("UMAPPlotWidth", "Width of plot", min = 100, max = 2000, value = 800, step = 10),
        # sliderInput("UMAPPlotHeight", "Height of plot", min = 100, max = 2000, value = 600, step = 10),
        
        numericInput(
          inputId = "pointSizeUMAP",
          label = "Sizes of points in UMAP plot", min = 1, max = 6,
          value = 3, step = 0.5
        ),
        numericInput(
          inputId = "textSizeUMAP",
          label = "Figure Font Size", min = 4, max = 30,
          value = 12, step = 0.5
        ),
        
      ) # close box
    ), # close column
    
    
  ) # close fluidRow
) # close tabItem
