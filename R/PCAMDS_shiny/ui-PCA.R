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
        
        # tags$head(
        #   tags$style(HTML('.shiny-split-layout>div {overflow: hidden;}')),
        # ),
        
        # checkboxInput(
        #   inputId = "displayTitlePCA",
        #   label = "Display Title",
        #   value = FALSE
        # ),
        # textInput(
        #   inputId = "pcaTitle",
        #   label = "Title of PCA plot",
        #   value = ""
        # ),
        fluidRow(
          column(
            width = 9,
            textInput(
              inputId = "pcaTitle",
              label = "Title of PCA plot",
              value = ""
            )
          ),
          column(
            width = 3,
            br(),
            checkboxInput(
              inputId = "displayTitlePCA",
              label = "Display Title",
              # label = "",
              value = FALSE
            )
          )
        ),
        
        checkboxInput(
          inputId = "sampleLabelsPCA",
          label = "Display sample labels",
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
        
        # selectInput(
        #   inputId = "excludeSamplesPCA",
        #   label = "Select samples to exclude",
        #   choices = "",
        #   selected = ""
        # ),
        
        # tags$div(sliderInput("slide1", "Slider1", min = 0, max=10, value=4),  style="display:inline-block"),
        tags$b("Select PCs"),
        # hr(),
        # br(),
        
        # tags$div(
        #   selectInput(
        #     inputId = "pickFactor1PCA",
        #     label = "x-axis",
        #     choices = "PC1",
        #     selected = "PC1",
        #     width = "100px"
        #     ),
        #   # width = "200%",
        #   style="display:inline-block"
        #   ),
        # tags$div(
        #   selectInput(
        #     inputId = "pickFactor2PCA",
        #     label = "y-axis",
        #     choices = "PC2",
        #     selected = "PC2",
        #     width = "100px"
        #   ),
        #   # width = "200%",
        #   style="display:inline-block"
        # ),
        
        fluidRow(
          column(
            selectInput(
              inputId = "pickFactor1PCA",
              label = "x-axis",
              choices = "PC1",
              selected = "PC1",
              # width = "100px"
            ),
            width = 6
          ),
          column(
            selectInput(
              inputId = "pickFactor2PCA",
              label = "y-axis",
              choices = "PC2",
              selected = "PC2",
              # width = "100px"
            ),
            width = 6
          )
        ),

        
        br(),
        
        # tags$b("Keep PCA axes proportional to variance?"),
        checkboxInput(
          inputId = "pcaAxesProp",
          label = "Keep axes proportional",
          # label = "",
          value = TRUE
        ),

        sliderInput(
          inputId = "pcaPlotWidth",
          label = "Width of plot",
          min = 100, max = 2000,
          value = 800, step = 10,
          ticks = FALSE
        ),
        sliderInput(
          inputId = "pcaPlotHeight",
          label = "Height of plot",
          min = 100, max = 2000,
          value = 600, step = 10,
          ticks = FALSE
        ),
        
        fluidRow(
          column(
            numericInput(
              inputId = "pointSizePCA",
              label = "Point size", min = 1, max = 6,
              value = 3, step = 0.5,
              # width = "100px"
            ),
            width = 6
          ),
          column(
            numericInput(
              inputId = "textSizePCA",
              label = "Font Size", min = 4, max = 30,
              value = 12, step = 0.5,
              # width = "100px"
            ),
            width = 6
          )
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
        # downloadButton(
        #   outputId = "downloadPCA",
        #   label = "Download PCA Plot (PDF)"
        # ),
        # div(
        #   style="display:inline-block",
        #   downloadBttn(
        #     outputId = "downloadPCA",
        #     label = "Download PCA Plot (PDF)",
        #     style = "unite",
        #     color = "primary",
        #     size = "sm",
        #     block = FALSE,
        #     no_outline = TRUE,
        #     icon = shiny::icon("download")
        #   ),
        #   style="float:right"
        # ),
        
        downloadBttn(
          outputId = "downloadPCA",
          label = "Download PCA Plot (PDF)",
          style = "unite",
          color = "primary",
          size = "sm",
          block = FALSE,
          no_outline = TRUE,
          icon = shiny::icon("download")
        ),
        br(), br(),
        plotOutput(
          outputId = "pcaStatic",
          # height = "80vh",
          # width = "100%",
          inline = T
        )
      ), # close box
    ), # close row
      
    # column(
    #   width = 3,
    #   box(
    #     title = "Brush points",
    #     width = NULL,
    #     solidHeader = TRUE,
    #     status = "primary",
    #     h4("Brushed points"), 
    #     # verbatimTextOutput("brush_info")
    #     tableOutput("brush_info"),
    #     br(), br(),
    #     textInput(
    #       inputId = "removeThesePointsPCA",
    #       label = "Type sample(s) to remove",
    #       value = ""  
    #     ),
    #     tags$b("Remove selected data and redo PCA?"),
    #     actionButton(
    #       inputId = "excludePointsPCA",
    #       label = "Remove selected data",
    #       icon = NULL
    #     )
    #   )
    # ),
    
    # column(
    #   width = 3,
    # ),
    
  fluidRow( 
    column(
      width = 12,
      box(
        title = "Scree Plot",
        width = 4, 
        solidHeader = TRUE,
        status = "primary",
        plotOutput(
          outputId ="pcaScree",
          inline = F,
          width = "100%"
        ),
        # DT::dataTableOutput("pcaVars")
      ),

      box(
        title = "PCA Loadings",
        width = 4, 
        solidHeader = TRUE,
        status = "primary",
        DT::dataTableOutput("pcaLoadings")
      ),

      box(
        title = "Sample removal",
        width = 4,
        solidHeader = TRUE,
        status = "primary",
        
        tags$b(
          "Select samples, then press the button below to redo the PCA without the samples
          displayed below"
        ),
        br(),br(),
        
        # h4("Brushed points"),
        # verbatimTextOutput("brush_info")
        # tableOutput(
        #   outputId = "brushInfo"
        # ),
        
        # DT::dataTableOutput("sampleTablePCA"),
        # checkboxGroupInput(
        #   inputId = "checkboxPCA",
        #   label = "checkboxGroup test",
        #   choices = NULL,
        #   selected = NULL,
        #   inline = FALSE,
        #   width = NULL,
        #   choiceNames = NULL,
        #   choiceValues = NULL
        # ),
        # 
        # verbatimTextOutput("choice"),
        # 
        # uiOutput("checkboxTablePCA"),
        # DT::dataTableOutput('x1'),
        # plotOutput('x2', height = 500),
        
        DT::dataTableOutput('sampleTablePCA'),
        br(),
        tags$head(
          tags$style(HTML("
            /* this will affect only the pre elements under the class myclass */
            .myclass pre {
              color: black;
              background-color: white;
              font-weight: bolder;
            }"))
        ),
        div(class = "myclass",
            verbatimTextOutput("selectedSamplesPCA")
        ),
        # verbatimTextOutput('selectedSamplesPCA'),
        # textOutput('selectedSamplesPCA'),
        br(),
        tags$b("Remove selected data and redo PCA?"),
        br(),br(),
        # actionButton(
        #   inputId = "removeSamplesPCA",
        #   label = "Redo PCA",
        #   class="btn btn-primary",
        #   icon = icon("repeat")
        # ),
        actionBttn(
          inputId = "removeSamplesPCA",
          label = "Redo PCA",
          icon = icon("repeat"),
          style = "simple",
          color = "primary",
          size = "md",
          block = FALSE,
          no_outline = TRUE
        ),
        br(),br(),
        actionButton(
          inputId = "resetSelectionPCA",
          label = "Reset selected samples",
          icon = NULL
        )
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
