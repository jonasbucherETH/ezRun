tabItem(
  tabName = "tab-PCA",
  fluidRow( ### NOTE: 1 row has width = 12
    column(
      width = 3,
      box(
        title = "Plot Settings",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        collapsible = TRUE,
        collapsed = FALSE,
        checkboxInput(
          inputId = "displayTitlePCA",
          label = "Display Title",
          value = FALSE
        ),
        textInput(
          inputId = "pcaTitle",
          label = "Title of PCA plot",
          value = ""
        ),
        checkboxInput(
          inputId = "sampleLabelsPCA",
          label = "Display sample labels",
          value = FALSE
        ),
        tags$b("Keep PCA axes proportional to variance?"),
        checkboxInput(
          inputId = "pcaAxesProp",
          label = "Keep axes proportional",
          value = TRUE
        ),
        
        selectInput(
          inputId = "colorGroupPCA",
          label = "Select groups to color by",
          choices = "",
          selected = ""
        ),
        selectInput(
          inputId = "shapeGroupPCA",
          label = "Select groups for shapes",
          choices = "",
          selected = ""
        ),
        
        selectInput(
          inputId = "excludeSamplesPCA",
          label = "Select samples to exclude",
          choices = "",
          selected = ""
        ),
        selectInput(
          inputId = "pickFactor1PCA",
          label = "Select PC for x-axis",
          choices = "PC1",
          selected = "PC1"
        ),
        selectInput(
          inputId = "pickFactor2PCA",
          label = "Select PC for y-axis",
          choices = "PC2",
          selected = "PC2"
        ),
        # sliderInput("pcaPlotWidth", "Width of plot", min = 100, max = 2000, value = 800, step = 10),
        # sliderInput("pcaPlotHeight", "Height of plot", min = 100, max = 2000, value = 600, step = 10),
        
        numericInput("pointSizePCA", "Sizes of points in PCA plot", min = 1, max = 6, value = 3, step = 0.5),
        numericInput(
          inputId = "textSizePCA",
          label = "Figure Font Size", min = 4, max = 30,
          value = 12, step = 0.5
        ),
        
        
        # selectInput(
        #   inputId = "pickFactor1PCA",
        #   label = "Select PC for x-axis",
        #   choices = c("PC1","PC2","PC3","PC4","PC5"),
        #   selected = "PC1"
        # ),
        # selectInput(
        #   inputId = "pickFactor2PCA",
        #   label = "Select PC for y-axis",
        #   choices = c("PC1","PC2","PC3","PC4","PC5"),
        #   selected = "PC2"
        # ),
        
        # selectInput("pickFactor1PCA", "Select PC for x-axis", choices = colnames(pca_tab)[PC_indeces]),
        # selectInput("pickFactor2PCA", "Select PC for y-axis", choices = colnames(pca_tab)[PC_indeces], selected = "PC2")
      ) # close box
    ), # close column
    
    column(
      width = 9, # 3 + 9 = 12 to fill row
      box(
        title = "PCA Plots",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        downloadButton(
          outputId = "downloadPCA",
          label = "Download PCA Plot (PDF)"
        ),
        br(), br(),
        plotOutput(
          outputId = "pcaStatic",
          height = "80vh",
          width = "100%",
          inline = F,
          brush = brushOpts(id = "pcaBrush")
        )
      ), # close box
      
      column(
        width = 3,
        box(
          title = "Brush points",
          width = NULL,
          solidHeader = TRUE,
          status = "primary",
          h4("Brushed points"), 
          # verbatimTextOutput("brush_info")
          tableOutput("brush_info"),
          br(), br(),
          textInput(
            inputId = "removeThesePointsPCA",
            label = "Type sample(s) to remove",
            value = ""  
          ),
          tags$b("Remove selected data and redo PCA?"),
          actionButton(
            inputId = "excludePointsPCA",
            label = "Remove selected data",
            icon = NULL
          )
        )
      ),
      
      column(
        width = 9,
        box(
          title = "Scree Plot",
          width = 6, 
          solidHeader = TRUE,
          status = "primary",
          plotOutput("pcaScree", inline = F),
          # DT::dataTableOutput("pcaVars")
        ),
        box(
          title = "PCA Loadings",
          width = 6, 
          solidHeader = TRUE,
          status = "primary",
          DT::dataTableOutput("pcaLoadings")
        )
      ),
      
      
      
      
      # box(
      #   title = "PCA Table",
      #   width = NULL,
      #   solidHeader = TRUE,
      #   status = "primary", br(),
      #   DT::dataTableOutput(outputId = "pca_table")
      # ) # close box
      
      #### testing stuff
      # ,
      # box(
      #   title = "Test text",
      #   width = NULL,
      #   solidHeader = TRUE,
      #   status = "primary",
      #   textOutput(outputId = "test_text")
      # ) # close box
      
      # ,
      # box(
      #   title = "Test pickFactor1PCA",
      #   width = NULL,
      #   solidHeader = TRUE,
      #   status = "primary",
      #   textOutput(outputId = "result")
      # ) # close box
    ) # close column
  ) # close fluidRow
) # close tabItem
