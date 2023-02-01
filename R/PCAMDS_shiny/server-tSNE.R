observe({
  withProgress(message = "Generating t-SNE Plot. Please wait...", {
    distanceMatrixTSNE <- inputDataReactive()$distanceMatrixTSNE
    groupingVariables <- inputDataReactive()$groupingVariables
    # colourList <- inputDataReactive()$colourList

    # cat(distanceMatrixTSNE[1:5,1:5])

    maxPerplexity <- floor((nrow(distanceMatrixTSNE) / 3) - 1)
    optimalPerplexity <- sqrt(nrow(distanceMatrixTSNE))
    optimalPerplexity <- round(optimalPerplexity / 5) * 5 # round to nearest 5
    updateNumericInput(
      session = session,
      inputId = "perplexityTSNE",
      value = optimalPerplexity,
      min = max(5, optimalPerplexity - 20),
      max = min(maxPerplexity, optimalPerplexity + 20),
      step = 5
    )

    observeEvent( # Event number 1
      {
        input$excludedSamplesTSNE
        input$goButton
        # input$perplexityTSNE # or precompute?
      },
      ignoreInit = F, # If TRUE, then, when the eventified object is first created/initialized, don't trigger the action or (compute the value). The default is FALSE.
      ignoreNULL = F, # default = TRUE
      {
        # # # TODO: omit selected samples/groups
        # distanceMatrixTSNE <- distanceMatrixTSNE

        # ---------------------------- compute & prepare t-SNE
        nGrouping <- ncol(groupingVariables)

        nDimsTSNE <- 2
        # distanceMatrixTSNE <- as.matrix(distanceMatrixTSNE)
        tsneResult <- Rtsne(
          distanceMatrixTSNE,
          # perplexity=as.numeric(input$perplexityTSNE), # default = 30; good range = 5-50
          perplexity = input$perplexityTSNE, # default = 30; good range = 5-50
          check_duplicates = FALSE,
          is_distance = TRUE,
          max_iter = 5000, # default= 1000
          theta = 0.5, # default = 0.5; exact t-SNE = 0.0
          eta = 200, # default = 200
          dims = nDimsTSNE
        ) # default: dims = 2

        ### perplexity: numeric; Perplexity parameter (should not be bigger than 3 * perplexity < nrow(X) - 1
        ### theta:	numeric; Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact TSNE (default: 0.5)

        tsneTable <- data.frame(groupingVariables, tsneResult$Y, stringsAsFactors = FALSE, row.names = rownames(groupingVariables))
        # colnames(tsneTable)[(nGrouping+1):-1] <- sub("X", "", colnames(tsneTable))

        # all PCs in array for selecting input
        dimList <- colnames(tsneTable)[-(1:nGrouping)]
        # all grouping variables in array for selecting input
        gvList <- c(colnames(tsneTable)[1:nGrouping], "-")


        ##### end of TSNE computation/prep

        ### TODO: put colorGroupTSNE default as 2nd variable in groupingVariables / TSNETable; make it necessary to have (or extend ifelse statement further down)
        # updateSelectInput(session = session, "colorGroupTSNE", choices = gvList, selected = gvList[2]) # try: selected = colnames(groupingVariables)[2] ([1] is sample name)
        # updateSelectInput(session = session, "shapeGroupTSNE", choices = gvList, selected = "")

        # levels(as.factor(TSNETable[[input$shapeGroupTSNE]]))

        # # Can use character(0) to remove all choices
        # if (is.null(x))
        #   x <- character(0)

        updateSelectInput(session = session, "colorGroupTSNE", choices = colnames(groupingVariables), selected = colnames(groupingVariables)[2]) # try: selected = colnames(groupingVariables)[2] ([1] is sample name)
        updateSelectInput(session = session, "shapeGroupTSNE", choices = gvList, selected = "-")

        updateSelectInput(session = session, "pickFactor1TSNE", choices = dimList, selected = "X1")
        updateSelectInput(session = session, "pickFactor2TSNE", choices = dimList, selected = "X2")

        # cat("before bindEvent 2")

        observeEvent(
          { # Event number 2: only need axis inputs here, for the others no need for redrawing plot (?)
            input$sampleLabelsTSNE
            input$pickFactor1TSNE
            input$pickFactor2TSNE
            input$colorGroupTSNE
            input$shapeGroupTSNE
            input$TSNEPlotWidth
            input$TSNEPlotHeight
            input$textSizeTSNE
            input$pointSizeTSNE
            # input$perplexityTSNE
            # input$goButton
          },
          {
            # cat("2")

            if (input$shapeGroupTSNE == "-") {
              plotTSNE <- tsneTable %>%
                ggplot(aes(x = .data[[input$pickFactor1TSNE]], y = .data[[input$pickFactor2TSNE]], color = .data[[input$colorGroupTSNE]], fill = .data[[input$colorGroupTSNE]])) +
                # ggplot(aes(x = .data[[input$pickFactor1TSNE]], y = .data[[input$pickFactor2TSNE]], color = .data[[!!input$colorGroupTSNE]], fill = .data[[!!input$colorGroupTSNE]])) +
                # ggplot(aes(x = .data[[input$pickFactor1TSNE]], y = .data[[input$pickFactor2TSNE]], color = .data[[rownames()]], fill = .data[[rownames()]])) +
                geom_point(size = as.numeric(input$pointSizeTSNE))
            } else {
              plotTSNE <- tsneTable %>%
                ggplot(aes(x = .data[[input$pickFactor1TSNE]], y = .data[[input$pickFactor2TSNE]], color = .data[[input$colorGroupTSNE]], fill = .data[[input$colorGroupTSNE]], shape = .data[[input$shapeGroupTSNE]])) +
                geom_point(size = as.numeric(input$pointSizeTSNE)) +
                scale_shape_manual(values = c(rep(c(21, 22, 23, 24, 25, 8, 3, 4), times = 10))[1:nlevels(as.factor(tsneTable[[input$shapeGroupTSNE]]))])
              ### Note: 21-25 are simple shapes, but with filling; 8,3,4 are crosses with different amount of lines
            }

            # if (input$shapeGroupTSNE == "Population") {
            #   cat("if statement works")
            #   plotTSNE <- TSNETable %>%
            #     ggplot(aes(x = .data[[input$pickFactor1TSNE]], y = .data[[input$pickFactor2TSNE]], color = .data[[input$colorGroupTSNE]], fill = .data[[input$colorGroupTSNE]], shape = .data[[input$shapeGroupTSNE]])) +
            #     geom_point()
            # } else {
            #   plotTSNE <- TSNETable %>%
            #     ggplot(aes(x = .data[[input$pickFactor1TSNE]], y = .data[[input$pickFactor2TSNE]], color = .data[[input$colorGroupTSNE]], fill = .data[[input$colorGroupTSNE]])) +
            #     geom_point()
            #   # scale_shape_manual(values = c(rep(c(21,22,23,24,25,8,3,4), times = 10 ))[1:nlevels(as.factor(datasetTSNE[[input$TSNEFactor2]]))] )
            #   ### Note: 21-25 are simple shapes, but with filling; 8,3,4 are crosses with different amount of lines
            # }
            # plotTSNE <- TSNETable %>%
            #   ggplot(aes(x = .data[[input$pickFactor1TSNE]], y = .data[[input$pickFactor2TSNE]], color = .data[[input$colorGroupTSNE]])) +
            #   geom_point()
            # scale_color_manual(values = colours)




            if (input$sampleLabelsTSNE) {
              plotTSNE <- plotTSNE +
                geom_text(
                  aes(label = rownames(tsneTable)),
                  size = (input$textSizeTSNE / 3),
                  hjust = 0.2, vjust = -1.5, check_overlap = T,
                  show.legend = FALSE
                ) # +
              # ylim(c(min(TSNETable[[input$pickFactor2TSNE]]*1.1), max(TSNETable[[input$pickFactor2TSNE]]*1.2))) +
              # xlim(c(min(TSNETable[[input$pickFactor1TSNE]]*1.1), max(TSNETable[[input$pickFactor1TSNE]]*1.2)))
            }


            ### themes, axis labels ,legend etc
            plotTSNE <- plotTSNE + labs(
              # x = paste0(input$pickFactor1TSNE, " (", format(round(tabVarprop[input$pickFactor1TSNE] * 100, 1), nsmall = 1), "%)"),
              # y = paste0(input$pickFactor2TSNE, " (", format(round(tabVarprop[input$pickFactor2TSNE] * 100, 1), nsmall = 1), "%)"),
              color = input$colorGroupTSNE, shape = input$shapeGroupTSNE
            ) +
              theme_bw() +
              theme(
                axis.text.x = element_text(
                  colour = "grey20", size = input$textSizeTSNE, angle = 0, hjust = .5,
                  vjust = .5, face = "plain"
                ),
                axis.text.y = element_text(
                  colour = "grey20", size = input$textSizeTSNE, angle = 0, hjust = 1,
                  vjust = 0.5, face = "plain"
                ),
                axis.title.x = element_text(
                  colour = "grey20", size = input$textSizeTSNE, angle = 0, hjust = .5,
                  vjust = 0, face = "plain"
                ),
                axis.title.y = element_text(
                  colour = "grey20", size = input$textSizeTSNE, angle = 90,
                  hjust = .5, vjust = .5, face = "plain"
                ),
                legend.text = element_text(
                  colour = "grey20", size = input$textSizeTSNE
                ),
                legend.title = element_text(
                  colour = "grey20", size = input$textSizeTSNE
                ),
                title = element_text(colour = "grey20", size = input$textSizeTSNE)
              )

            plotTSNE <- plotTSNE +
              ylim(c(min(tsneTable[[input$pickFactor2TSNE]] * 1.1), max(tsneTable[[input$pickFactor2TSNE]] * 1.2))) +
              xlim(c(min(tsneTable[[input$pickFactor1TSNE]] * 1.1), max(tsneTable[[input$pickFactor1TSNE]] * 1.2)))

            output$tsneStatic <- renderPlot(
              {
                plotTSNE
              },
              # width = as.numeric(input$TSNEPlotWidth),
              # height = as.numeric(input$TSNEPlotHeight)
              # , res = 96
            )
            
            output$downloadTSNE <- downloadHandler(
              filename = function() {
                paste0(Sys.Date(), "_tSNE.pdf")
              },
              content = function(file) {
                ggsave(
                  filename = file, plot = plotTSNE # ,
                  # width = (as.numeric(input$figWidthPCA) / 3.2),
                  # height = (as.numeric(input$figHeightPCA) / 3.2), limitsize = FALSE,
                  # units = "mm"
                )
              }
            )
            

            # }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")} # close withProgress
            # }) # close tryCatch
          }
        ) # close Event number 2
      }
    ) # close Event number 1
  }) # close withProgress
}) # close observe
