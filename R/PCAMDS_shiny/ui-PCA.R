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
        collapsed = FALSE,
      
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
          choices = c("PC1","PC2","PC3","PC4","PC5"), 
          selected = "PC1"
        ),
        selectInput(
          inputId = "pick_pc_y",
          label = "Select PC for y-axis",
          choices = c("PC1","PC2","PC3","PC4","PC5"), 
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
      
      #### testing stuff
      ,
      box(
        title = "Test text",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        textOutput(outputId = "test_text")
      ) # close box
      ,
      box(
        title = "Test pick_pc_x",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        textOutput(outputId = "result")
      ) # close box

      
      
    ) # close column
  ) # close fluidRow
) # close tabItem
