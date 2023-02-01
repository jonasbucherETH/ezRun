observe({
  withProgress(message = "Generating PCA Plots. Please wait...", {
    vcfRaw <- inputDataReactive()$vcfRaw
    vcfGenind <- inputDataReactive()$vcfGenind
    datasetPCA <- inputDataReactive()$datasetPCA
    groupingVariables <- inputDataReactive()$groupingVariables
    colourList <- inputDataReactive()$colourList
    
    ##### -----#####
    # Use observeEvent whenever you want to perform an action in response to an event.
    # (Note that "recalculate a value" does not generally count as performing an action--see eventReactive for that.)
    # The first argument is the event you want to respond to,
    # and the second argument is a function that should be called whenever the event occurs.

    # The observeEvent will only be dependent on the 'event' section in the small piece of code above. (first {...})
    # It will not be dependent on anything in the 'code to run' part. (after ignoreInit)
    ##### -----#####

    observeEvent( # Event number 1
      {
        input$excludedSamplesPCA
        # input$sampleLabelsPCA
        # input$colorGroupPCA
        # input$shapeGroupPCA
        # input$pickFactor1PCA
        # input$pickFactor2PCA
        # input$pcaPlotWidth
        # input$pcaPlotHeight
      },
      ignoreInit = F, # If TRUE, then, when the eventified object is first created/initialized, don't trigger the action or (compute the value). The default is FALSE.
      ignoreNULL = T, # default = TRUE
      {
        # # # TODO: omit selected samples/groups
        # datasetPCA <- datasetPCA

        ##### compute PCA
        # pca <- dudi.pca(X, center = TRUE, scale = TRUE, scan = FALSE, nf = 5)
        pcaResults <- dudi.pca(datasetPCA, center = TRUE, scale = TRUE, scan = FALSE, nf = 5)
        nPC <- pcaResults$nf # number of (> 0) principal components
        # nGrouping <- ncol(groupingVariables)
        nGrouping <- ncol(groupingVariables)
        # eigSum <- sum(pcaResult$eig)
        pcaVarprop <- pcaResults$eig
        # eigSum <- sum(pcaVarprop)
        
        pcaVarprop <- tibble(PC = paste0("PC", factor(1:length(pcaVarprop))), variance = pcaVarprop) %>% 
          mutate(pct = format(variance/sum(variance)*100, digits = 2)) %>%
          # mutate(pct = variance/sum(variance)*100) %>% 
          mutate(pct_cum = cumsum(pct))
        pcaVarprop$PC <- factor(pcaVarprop$PC, levels = pcaVarprop$PC)
        
        # pcaVarprop <- pcaResults$eig / eigSum
        # pcaVarprop <- pcaVarprop[1:nPC]
        pcaTable <- data.frame(groupingVariables, pcaResults$li, stringsAsFactors = FALSE, row.names = rownames(groupingVariables))
        # tabVarprop <- as.data.frame(t(pcaVarprop), stringsAsFactors = FALSE)
        # tabVarprop <- format(round(tabVarprop * 100, 2), nsmall = 2)
        tabVarprop <- pcaVarprop
        # PC_indeces <- seq(1+nGrouping, ncol(pcaTable))
        for (i in 1:nPC) {
          colnames(pcaTable)[i + nGrouping] <- paste0("PC", i)
          # colnames(tabVarprop)[i] <- paste0("PC", i)
        }

        # all PCs in array for selecting input
        pcList <- colnames(pcaTable)[-(1:nGrouping)]
        # all grouping variables in array for selecting input
        gvList <- c(colnames(pcaTable)[1:nGrouping], "-")
        
        # Get the gene loadings (in dudi.pca: $c1)
        pc_loadings <- pcaResults$c1
        colnames(pc_loadings) <- c("PC1", "PC2") 
        # pc_loadings <- pc_loadings %>% 
        #   as_tibble(rownames = "gene")
        pc_loadings$gene <- rownames(pc_loadings)
        top_genes <- pc_loadings %>% 
          select(gene, PC1, PC2) %>%
          pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
          group_by(PC) %>% 
          arrange(desc(abs(loading)))
        

        ##### end of PCA computation/prep

        ### TODO: put colorGroupPCA default as 2nd variable in groupingVariables / pcaTable; make it necessary to have (or extend ifelse statement further down)
        # updateSelectInput(session = session, "colorGroupPCA", choices = gvList, selected = gvList[2]) # try: selected = colnames(groupingVariables)[2] ([1] is sample name)
        # updateSelectInput(session = session, "shapeGroupPCA", choices = gvList, selected = "")

        # levels(as.factor(pcaTable[[input$shapeGroupPCA]]))

        # # Can use character(0) to remove all choices
        # if (is.null(x))
        #   x <- character(0)

        updateSelectInput(session = session, "colorGroupPCA", choices = colnames(groupingVariables), selected = colnames(groupingVariables)[2]) # try: selected = colnames(groupingVariables)[2] ([1] is sample name)
        updateSelectInput(session = session, "shapeGroupPCA", choices = gvList, selected = "-")

        updateSelectInput(session = session, "pickFactor1PCA", choices = pcList, selected = "PC1")
        updateSelectInput(session = session, "pickFactor2PCA", choices = pcList, selected = "PC2")

        # cat("before bindEvent 2")

        observeEvent(
          { # Event number 2: only need axis inputs here, for the others no need for redrawing plot (?)
            input$sampleLabelsPCA

            input$pickFactor1PCA
            input$pickFactor2PCA
            input$colorGroupPCA
            input$shapeGroupPCA
            input$pcaPlotWidth
            input$pcaPlotHeight
            input$textSizePCA
            input$pointSizePCA
            
            input$displayTitlePCA
          },
          {
            # cat("2")

            # colourList <- list()
            # for (i in seq_along(groupingVariables)){
            #   colourList[i] <- paste0(col2hex(input[[paste0("GroupColour", i)]]), "FF")
            # }
            # colours <- NULL
            # for (i in levels(as.factor(pcaTable[[input$colorGroupPCA]]))) {
            #   colours[i] <- colourList[i]
            # }
            # colours <- colours[names(colours) %in% input$colorGroupPCA]
            #
            # print(colourList[1:5])


            # plotPCA <- pcaTable %>%
            #   ggplot(aes(x = .data[[input$pickFactor1PCA]], y = .data[[input$pickFactor2PCA]])) +
            #   geom_point()

            ### TODO: replace deprecated aes_string (above) with .data[[]]
            # if (is.null(input$shapeGroupPCA)) {
            if (input$shapeGroupPCA == "-") {
              # cat("if statement works")
              plotPCA <- pcaTable %>%
                ggplot(aes(x = .data[[input$pickFactor1PCA]], y = .data[[input$pickFactor2PCA]], color = .data[[input$colorGroupPCA]], fill = .data[[input$colorGroupPCA]])) +
                # ggplot(aes(x = .data[[input$pickFactor1PCA]], y = .data[[input$pickFactor2PCA]], color = .data[[!!input$colorGroupPCA]], fill = .data[[!!input$colorGroupPCA]])) +
                # ggplot(aes(x = .data[[input$pickFactor1PCA]], y = .data[[input$pickFactor2PCA]], color = .data[[rownames()]], fill = .data[[rownames()]])) +
                geom_point(size = as.numeric(input$pointSizePCA))
            } else {
              plotPCA <- pcaTable %>%
                ggplot(aes(x = .data[[input$pickFactor1PCA]], y = .data[[input$pickFactor2PCA]], color = .data[[input$colorGroupPCA]], fill = .data[[input$colorGroupPCA]], shape = .data[[input$shapeGroupPCA]])) +
                geom_point(size = as.numeric(input$pointSizePCA)) +
                scale_shape_manual(values = c(rep(c(21, 22, 23, 24, 25, 8, 3, 4), times = 10))[1:nlevels(as.factor(pcaTable[[input$shapeGroupPCA]]))])
              ### Note: 21-25 are simple shapes, but with filling; 8,3,4 are crosses with different amount of lines
            }

            # if (input$shapeGroupPCA == "Population") {
            #   cat("if statement works")
            #   plotPCA <- pcaTable %>%
            #     ggplot(aes(x = .data[[input$pickFactor1PCA]], y = .data[[input$pickFactor2PCA]], color = .data[[input$colorGroupPCA]], fill = .data[[input$colorGroupPCA]], shape = .data[[input$shapeGroupPCA]])) +
            #     geom_point()
            # } else {
            #   plotPCA <- pcaTable %>%
            #     ggplot(aes(x = .data[[input$pickFactor1PCA]], y = .data[[input$pickFactor2PCA]], color = .data[[input$colorGroupPCA]], fill = .data[[input$colorGroupPCA]])) +
            #     geom_point()
            #   # scale_shape_manual(values = c(rep(c(21,22,23,24,25,8,3,4), times = 10 ))[1:nlevels(as.factor(datasetPCA[[input$pcaFactor2]]))] )
            #   ### Note: 21-25 are simple shapes, but with filling; 8,3,4 are crosses with different amount of lines
            # }
            # plotPCA <- pcaTable %>%
            #   ggplot(aes(x = .data[[input$pickFactor1PCA]], y = .data[[input$pickFactor2PCA]], color = .data[[input$colorGroupPCA]])) +
            #   geom_point()
            # scale_color_manual(values = colours)




            if (input$sampleLabelsPCA) {
              plotPCA <- plotPCA +
                geom_text(
                  aes(label = rownames(pcaTable)),
                  size = (input$textSizePCA / 3),
                  hjust = 0.2, vjust = -1.5, check_overlap = T,
                  show.legend = FALSE
                ) # +
              # ylim(c(min(pcaTable[[input$pickFactor2PCA]]*1.1), max(pcaTable[[input$pickFactor2PCA]]*1.2))) +
              # xlim(c(min(pcaTable[[input$pickFactor1PCA]]*1.1), max(pcaTable[[input$pickFactor1PCA]]*1.2)))
            }
            
            if (input$displayTitlePCA) {
              plotPCA <- plotPCA + labs(
                title = input$pcaTitle
              )
            } 


            ### themes, axis labels ,legend etc
            plotPCA <- plotPCA + labs(
              # x = paste0(input$pickFactor1PCA, " (", format(round(tabVarprop[input$pickFactor1PCA] * 100, 1), nsmall = 1), "%)"),
              # y = paste0(input$pickFactor2PCA, " (", format(round(tabVarprop[input$pickFactor2PCA] * 100, 1), nsmall = 1), "%)"),
              x = paste0(input$pickFactor1PCA, " (", tabVarprop$pct[tabVarprop$PC == input$pickFactor1PCA], "% variance)"),
              y = paste0(input$pickFactor2PCA, " (", tabVarprop$pct[tabVarprop$PC == input$pickFactor2PCA], "% variance)"),
              color = input$colorGroupPCA, shape = input$shapeGroupPCA
            ) +
              theme_bw() +
              theme(
                axis.text.x = element_text(
                  colour = "grey20", size = input$textSizePCA, angle = 0, hjust = .5,
                  vjust = .5, face = "plain"
                ),
                axis.text.y = element_text(
                  colour = "grey20", size = input$textSizePCA, angle = 0, hjust = 1,
                  vjust = 0.5, face = "plain"
                ),
                axis.title.x = element_text(
                  colour = "grey20", size = input$textSizePCA, angle = 0, hjust = .5,
                  vjust = 0, face = "plain"
                ),
                axis.title.y = element_text(
                  colour = "grey20", size = input$textSizePCA, angle = 90,
                  hjust = .5, vjust = .5, face = "plain"
                ),
                legend.text = element_text(
                  colour = "grey20", size = input$textSizePCA
                ),
                legend.title = element_text(
                  colour = "grey20", size = input$textSizePCA
                ),
                title = element_text(colour = "grey20", size = input$textSizePCA)
              )


            plotPCA <- plotPCA +
              ylim(c(min(pcaTable[[input$pickFactor2PCA]] * 1.1), max(pcaTable[[input$pickFactor2PCA]] * 1.2))) +
              xlim(c(min(pcaTable[[input$pickFactor1PCA]] * 1.1), max(pcaTable[[input$pickFactor1PCA]] * 1.2)))

            output$pcaStatic <- renderPlot(
              {
                plotPCA
              },
              # width = as.numeric(input$pcaPlotWidth),
              # height = as.numeric(input$pcaPlotHeight)
              # , res = 96
            )
            
            output$pcaScree <- renderPlot({
              tabVarprop2 <- tabVarprop
              # pc_eigenvalues2$PC <- gsub("PC", "", pc_eigenvalues2$PC)
              # pc_eigenvalues2$PC <- factor(pc_eigenvalues2$PC, levels = pc_eigenvalues2$PC)
              if (nrow(tabVarprop2) > 15) {
                tabVarprop2 <- tabVarprop2[1:15,]
              }
              pcaScree <- tabVarprop2 %>%
                  ggplot(aes(x = PC)) +
                  geom_col(aes(y = pct)) +
                  geom_line(aes(y = pct_cum, group = 1)) +
                  geom_point(aes(y = pct_cum)) +
                  labs(x = "Principal component", y = "Fraction variance explained") +
                  theme_classic(base_size = as.numeric(input$textSizePCA))
              pcaScree
            
            }, width = 400, height = 400)
            
            output$downloadPCA <- downloadHandler(
              filename = function() {
                paste0(Sys.Date(), "_PCA.pdf")
              },
              content = function(file) {
                ggsave(
                  filename = file, plot = plotPCA # ,
                  # width = (as.numeric(input$figWidthPCA) / 3.2),
                  # height = (as.numeric(input$figHeightPCA) / 3.2), limitsize = FALSE,
                  # units = "mm"
                )
              }
            )
            
            output$pcaLoadings <- DT::renderDataTable({
              datatable(top_genes, rownames = F) %>% formatRound("loading", digits = 3)
            })
            
            # }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")} # close withProgress
            # }) # close tryCatch
          }
        ) # close Event number 2
      }
    ) # close Event number 1
  }) # close withProgress
}) # close observe
