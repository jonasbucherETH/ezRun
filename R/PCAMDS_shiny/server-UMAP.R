observe({
  withProgress(message = "Generating UMAP Plot. Please wait...", {
    datasetScaled <- inputDataReactive()$datasetScaled
    groupingVariables <- inputDataReactive()$groupingVariables
    vcfGenind <- inputDataReactive()$vcfGenind
    
    nGrouping <- ncol(groupingVariables)
    gvList <- c(colnames(groupingVariables)[1:nGrouping], "-")
    updateSelectInput(session = session, "colorGroupUMAP", choices = colnames(groupingVariables), selected = colnames(groupingVariables)[2]) # try: selected = colnames(groupingVariables)[2] ([1] is sample name)
    updateSelectInput(session = session, "shapeGroupUMAP", choices = gvList, selected = "-")
    
    # cat("0")
    # updateSelectInput(session = session, "pickFactor1UMAP", choices = pcList, selected = "X1")
    # updateSelectInput(session = session, "pickFactor2UMAP", choices = pcList, selected = "X2")

    observeEvent( # Event number 1
      {
        input$excludeSamplesUMAP
        # input$sampleLabelsUMAP
        # input$pickFactor1UMAP
        # input$pickFactor2UMAP
        # input$colorGroupUMAP
        # input$shapeGroupUMAP
        # input$textSizeUMAP
        # input$pointSizeUMAP
        # input$displayTitleUMAP

      },
      ignoreInit = F, # If TRUE, then, when the eventified object is first created/initialized, don't trigger the action or (compute the value). The default is FALSE.
      ignoreNULL = T, # default = TRUE
      {
        
        umapResult <- umap(
          X = datasetScaled,
          n_neighbors = 15,
          n_components = 2,
          metric = "euclidean"
        ) # returns matrix
        
        UMAPTable <- data.frame(groupingVariables, umapResult, stringsAsFactors = FALSE, row.names = rownames(groupingVariables))
        
        observeEvent( # Event number 2
          {
            input$sampleLabelsUMAP
            # input$pickFactor1UMAP
            # input$pickFactor2UMAP
            input$colorGroupUMAP
            input$shapeGroupUMAP
            input$textSizeUMAP
            input$pointSizeUMAP
            input$displayTitleUMAP
            
          },
            ignoreInit = F, # default = FALSE
            ignoreNULL = T, # default = TRUE
          {
        
            if (input$shapeGroupUMAP == "-") {
              plotUMAP <- UMAPTable %>%
                ggplot(aes(x = X1, y = X2, color = .data[[input$colorGroupUMAP]], fill = .data[[input$colorGroupUMAP]])) +
                geom_point(size = as.numeric(input$pointSizeUMAP))
            } else {
              plotUMAP <- UMAPTable %>%
                ggplot(aes(x = X1, y = X2, color = .data[[input$colorGroupUMAP]], fill = .data[[input$colorGroupUMAP]], shape = .data[[input$shapeGroupUMAP]])) +
                geom_point(size = as.numeric(input$pointSizeUMAP)) +
                scale_shape_manual(values = c(rep(c(21, 22, 23, 24, 25, 8, 3, 4), times = 10))[1:nlevels(as.factor(UMAPTable[[input$shapeGroupUMAP]]))])
              ### Note: 21-25 are simple shapes, but with filling; 8,3,4 are crosses with different amount of lines
            }
            
            if (input$sampleLabelsUMAP) {
              plotUMAP <- plotUMAP +
                geom_text(
                  aes(label = rownames(UMAPTable)),
                  size = (input$textSizeUMAP / 3),
                  hjust = 0.2, vjust = -1.5, check_overlap = T,
                  show.legend = FALSE
                ) 
            }
            
            if (input$displayTitleUMAP) {
              plotUMAP <- plotUMAP + labs(
                title = input$UMAPTitle
              )
            } 
            
            ### themes, axis labels ,legend etc
            plotUMAP <- plotUMAP + labs(
              x = "X1",
              y = "X2",
              color = input$colorGroupUMAP, shape = input$shapeGroupUMAP
            ) +
              theme_bw() +
              theme(
                axis.text.x = element_text(
                  colour = "grey20", size = input$textSizeUMAP, angle = 0, hjust = .5,
                  vjust = .5, face = "plain"
                ),
                axis.text.y = element_text(
                  colour = "grey20", size = input$textSizeUMAP, angle = 0, hjust = 1,
                  vjust = 0.5, face = "plain"
                ),
                axis.title.x = element_text(
                  colour = "grey20", size = input$textSizeUMAP, angle = 0, hjust = .5,
                  vjust = 0, face = "plain"
                ),
                axis.title.y = element_text(
                  colour = "grey20", size = input$textSizeUMAP, angle = 90,
                  hjust = .5, vjust = .5, face = "plain"
                ),
                legend.text = element_text(
                  colour = "grey20", size = input$textSizeUMAP
                ),
                legend.title = element_text(
                  colour = "grey20", size = input$textSizeUMAP
                ),
                title = element_text(colour = "grey20", size = input$textSizeUMAP)
              )
            
            
            # plotUMAP <- plotUMAP +
            #   ylim(c(min(UMAPTable[[input$pickFactor2UMAP]] * 1.1), max(UMAPTable[[input$pickFactor2UMAP]] * 1.2))) +
            #   xlim(c(min(UMAPTable[[input$pickFactor1UMAP]] * 1.1), max(UMAPTable[[input$pickFactor1UMAP]] * 1.2)))
            
            output$UMAPStatic <- renderPlot(
              {
                plotUMAP
              },
              # width = as.numeric(input$TSNEPlotWidth),
              # height = as.numeric(input$TSNEPlotHeight)
              # , res = 96
            )
          } # close Event number 2
        ) # close Event number 2
        
      } # close Event number 1
    ) # close Event number 1
    
  }) # close withProgress
}) # close observe
