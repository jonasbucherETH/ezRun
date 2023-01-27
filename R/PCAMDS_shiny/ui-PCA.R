tabItem(
  tabName = "tab-PCA",
  fluidRow( ### NOTE: 1 row has width = 12
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
          inputId = "sample_labels",
          label = "Display sample labels", 
          value = FALSE
        ),

        ### grouping input (colors & shapes)
        # selectInput("color_by", "Select groups to color by", choices = c("", colnames(grouping_vars))),
        # selectInput("group_by", "Select groups for shapes",  choices = c("", colnames(grouping_vars))),
        # selectInput("color_by", "Select groups to color by", choices = colnames(grouping_vars)),
        # selectInput("group_by", "Select groups for shapes",  choices = colnames(grouping_vars)),

        # varSelectInput("color_group", "Select groups to color by", data = pca_tab[,-PC_indeces]),
        # varSelectInput("shape_group", "Select groups for shapes", data = pca_tab[,-PC_indeces]),
        # selectInput("color_group", "Select groups to color by", choices = list("Sample", "Population")),
        # selectInput("shape_group", "Select groups for shapes", choices = list("Sample", "Population"), selected = "Population"),

        ### no choices, selected = "" as default 
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
        
        ## temporary selectInputs for color & shape
        # selectInput(
        #   inputId = "color_group",
        #   label = "Select groups to color by",
        #   choices = c("Sample", "Population"),
        #   selected = "Population"
        # ),
        # selectInput(
        #   inputId = "shape_group",
        #   label = "Select groups for shapes",
        #   choices = c("", "Sample", "Population"),
        #   selected = ""
        # ),

        ### PC input (axis)
        # selectInput(inputId = "pick_pc_x", label = "Select PC for x-axis", choices = list("PC1", "PC2", "PC3", "PC4")),
        # selectInput(inputId = "pick_pc_y", label = "Select PC for y-axis", choices = list("PC1", "PC2", "PC3", "PC4"), selected = "PC2")

        # selectInput(inputId = "pick_pc_x", label = "Select PC for x-axis", choices = colnames(pca_tab)[PC_indeces]),
        # selectInput(inputId = "pick_pc_y", label = "Select PC for y-axis", choices = colnames(pca_tab)[PC_indeces], selected = colnames(pca_tab)[PC_indeces[2]])

        selectInput(
          inputId = "pick_pc_x",
          label = "Select PC for x-axis",
          choices = "PC1",
          selected = "PC1"
        ),
        selectInput(
          inputId = "pick_pc_y",
          label = "Select PC for y-axis",
          choices = "PC2",
          selected = "PC2"
        ),
        sliderInput("pcaPlotHeight", "Height of plot", min = 100, max = 500, value = 250),
        sliderInput("pcaPlotWidth", "Width of plot", min = 100, max = 500, value = 250),

        # selectInput(
        #   inputId = "pick_pc_x",
        #   label = "Select PC for x-axis",
        #   choices = c("PC1","PC2","PC3","PC4","PC5"),
        #   selected = "PC1"
        # ),
        # selectInput(
        #   inputId = "pick_pc_y",
        #   label = "Select PC for y-axis",
        #   choices = c("PC1","PC2","PC3","PC4","PC5"),
        #   selected = "PC2"
        # ),

        # selectInput("pick_pc_x", "Select PC for x-axis", choices = colnames(pca_tab)[PC_indeces]),
        # selectInput("pick_pc_y", "Select PC for y-axis", choices = colnames(pca_tab)[PC_indeces], selected = "PC2")
      ) # close box
    ), # close column

    column(
      width = 9, # 3 + 9 = 12 to fill row
      box(
        title = "PCA Plots",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        plotOutput(
          outputId = "pca_plot",
          height = "50vh", 
          width = "100%",
          inline = F 
        )
      ), # close box
      
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
      #   title = "Test pick_pc_x",
      #   width = NULL,
      #   solidHeader = TRUE,
      #   status = "primary",
      #   textOutput(outputId = "result")
      # ) # close box

      
      
    ) # close column
  ) # close fluidRow
) # close tabItem
