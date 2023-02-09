observe({
  withProgress(message = "Generating DAPC Plots. Please wait...", {
    vcfRaw <- inputDataReactive()$vcfRaw
    vcfGenind <- inputDataReactive()$vcfGenind
    datasetScaled <- inputDataReactive()$datasetScaled
    groupingVariables <- inputDataReactive()$groupingVariables
    colourList <- inputDataReactive()$colourList

    observeEvent( # Event number 1
      {
        input$excludeSamplesDAPC

      },
      ignoreInit = F, # If TRUE, then, when the eventified object is first created/initialized, don't trigger the action or (compute the value). The default is FALSE.
      ignoreNULL = F, # default = TRUE
      {
        
        ##### compute DAPC
        # 1st input: use as many PCA as possible basically
        # 2nd input: (local) minimum
        groupsDAPC <- find.clusters(
          x = vcfGenind,
          clust = NULL,
          max.n.clust=40
        )
        groupsDAPC$size
        table(pop(vcfGenind), groupsDAPC$grp)
        
        # 1st input retain a few PCs without sacrificing too much information
        # 2nd input: For small number of clusters, all eigenvalues can be
        # retained since all discriminant functions can be examined without
        # difficulty. Whenever more (say, tens of) clusters are analysed,
        # it is likely that the first few dimensions will carry more information
        # than the others, and only those can then be retained and interpreted.
        DAPCResults <- dapc(vcfGenind, groupsDAPC$grp)
        
        scatter(DAPCResults)
        

        ############## ----- to delete
        nPC <- DAPCResults$nf # number of (> 0) principal components
        # nGrouping <- ncol(groupingVariables)
        nGrouping <- ncol(groupingVariables)
        # eigSum <- sum(DAPCResult$eig)
        DAPCVarprop <- DAPCResults$eig
        # eigSum <- sum(DAPCVarprop)
        
        DAPCVarprop <- tibble(PC = paste0("PC", factor(1:length(DAPCVarprop))), variance = DAPCVarprop) %>% 
          mutate(pct = format(variance/sum(variance)*100, digits = 2)) %>%
          # mutate(pct = variance/sum(variance)*100) %>% 
          mutate(pct_cum = cumsum(pct))
        DAPCVarprop$PC <- factor(DAPCVarprop$PC, levels = DAPCVarprop$PC)
        
        # DAPCVarprop <- DAPCResults$eig / eigSum
        # DAPCVarprop <- DAPCVarprop[1:nPC]
        DAPCTable <- data.frame(groupingVariables, DAPCResults$li, stringsAsFactors = FALSE, row.names = rownames(groupingVariables))
        # tabVarprop <- as.data.frame(t(DAPCVarprop), stringsAsFactors = FALSE)
        # tabVarprop <- format(round(tabVarprop * 100, 2), nsmall = 2)
        tabVarprop <- DAPCVarprop
        # PC_indeces <- seq(1+nGrouping, ncol(DAPCTable))
        for (i in 1:nPC) {
          colnames(DAPCTable)[i + nGrouping] <- paste0("PC", i)
          # colnames(tabVarprop)[i] <- paste0("PC", i)
        }
        
        # all PCs in array for selecting input
        pcList <- colnames(DAPCTable)[-(1:nGrouping)]
        # all grouping variables in array for selecting input
        gvList <- c(colnames(DAPCTable)[1:nGrouping], "-")
        
        # Get the gene loadings (in dudi.DAPC: $c1)
        pc_loadings <- DAPCResults$c1
        colnames(pc_loadings) <- c("PC1", "PC2") 
        # pc_loadings <- pc_loadings %>% 
        #   as_tibble(rownames = "gene")
        pc_loadings$gene <- rownames(pc_loadings)
        top_genes <- pc_loadings %>% 
          select(gene, PC1, PC2) %>%
          pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
          group_by(PC) %>% 
          arrange(desc(abs(loading)))
        
        
        ##### end of DAPC computation/prep
        
        ### TODO: put colorGroupDAPC default as 2nd variable in groupingVariables / DAPCTable; make it necessary to have (or extend ifelse statement further down)
        # updateSelectInput(session = session, "colorGroupDAPC", choices = gvList, selected = gvList[2]) # try: selected = colnames(groupingVariables)[2] ([1] is sample name)
        # updateSelectInput(session = session, "shapeGroupDAPC", choices = gvList, selected = "")
        
        # levels(as.factor(DAPCTable[[input$shapeGroupDAPC]]))
        
        # # Can use character(0) to remove all choices
        # if (is.null(x))
        #   x <- character(0)
        
        updateSelectInput(session = session, "colorGroupDAPC", choices = colnames(groupingVariables), selected = colnames(groupingVariables)[2]) # try: selected = colnames(groupingVariables)[2] ([1] is sample name)
        updateSelectInput(session = session, "shapeGroupDAPC", choices = gvList, selected = "-")
        
        updateSelectInput(session = session, "pickFactor1DAPC", choices = pcList, selected = "PC1")
        updateSelectInput(session = session, "pickFactor2DAPC", choices = pcList, selected = "PC2")
        
        # cat("before bindEvent 2")
        
        observeEvent(
          { # Event number 2: only need axis inputs here, for the others no need for redrawing plot (?)
            input$sampleLabelsDAPC
            input$pickFactor1DAPC
            input$pickFactor2DAPC
            input$colorGroupDAPC
            input$shapeGroupDAPC
            input$DAPCPlotWidth
            input$DAPCPlotHeight
            input$textSizeDAPC
            input$pointSizeDAPC
            input$displayTitleDAPC
            # input$DAPCBrush
          },
          ignoreInit = F, # default = FALSE
          ignoreNULL = T, # default = TRUE
          {
           
            ### TODO: replace deprecated aes_string (above) with .data[[]]
            # if (is.null(input$shapeGroupDAPC)) {
            if (input$shapeGroupDAPC == "-") {
              # cat("if statement works")
              plotDAPC <- DAPCTable %>%
                ggplot(aes(x = .data[[input$pickFactor1DAPC]], y = .data[[input$pickFactor2DAPC]], color = .data[[input$colorGroupDAPC]], fill = .data[[input$colorGroupDAPC]])) +
                # ggplot(aes(x = .data[[input$pickFactor1DAPC]], y = .data[[input$pickFactor2DAPC]], color = .data[[!!input$colorGroupDAPC]], fill = .data[[!!input$colorGroupDAPC]])) +
                # ggplot(aes(x = .data[[input$pickFactor1DAPC]], y = .data[[input$pickFactor2DAPC]], color = .data[[rownames()]], fill = .data[[rownames()]])) +
                geom_point(size = as.numeric(input$pointSizeDAPC))
            } else {
              plotDAPC <- DAPCTable %>%
                ggplot(aes(x = .data[[input$pickFactor1DAPC]], y = .data[[input$pickFactor2DAPC]], color = .data[[input$colorGroupDAPC]], fill = .data[[input$colorGroupDAPC]], shape = .data[[input$shapeGroupDAPC]])) +
                geom_point(size = as.numeric(input$pointSizeDAPC)) +
                scale_shape_manual(values = c(rep(c(21, 22, 23, 24, 25, 8, 3, 4), times = 10))[1:nlevels(as.factor(DAPCTable[[input$shapeGroupDAPC]]))])
              ### Note: 21-25 are simple shapes, but with filling; 8,3,4 are crosses with different amount of lines
            }
            
            if (isTRUE(input$DAPCAxesProp)) {
              plotDAPC <- plotDAPC + coord_fixed(ratio = as.numeric(tabVarprop$pct[tabVarprop$PC == input$pickFactor2DAPC])/as.numeric(tabVarprop$pct[tabVarprop$PC == input$pickFactor1DAPC]))
            }
            
            
            if (input$sampleLabelsDAPC) {
              plotDAPC <- plotDAPC +
                geom_text(
                  aes(label = rownames(DAPCTable)),
                  size = (input$textSizeDAPC / 3),
                  hjust = 0.2, vjust = -1.5, check_overlap = T,
                  show.legend = FALSE
                ) # +
              # ylim(c(min(DAPCTable[[input$pickFactor2DAPC]]*1.1), max(DAPCTable[[input$pickFactor2DAPC]]*1.2))) +
              # xlim(c(min(DAPCTable[[input$pickFactor1DAPC]]*1.1), max(DAPCTable[[input$pickFactor1DAPC]]*1.2)))
            }
            
            if (input$displayTitleDAPC) {
              plotDAPC <- plotDAPC + labs(
                title = input$DAPCTitle
              )
            } 
            
            
            ### themes, axis labels ,legend etc
            plotDAPC <- plotDAPC + labs(
              # x = paste0(input$pickFactor1DAPC, " (", format(round(tabVarprop[input$pickFactor1DAPC] * 100, 1), nsmall = 1), "%)"),
              # y = paste0(input$pickFactor2DAPC, " (", format(round(tabVarprop[input$pickFactor2DAPC] * 100, 1), nsmall = 1), "%)"),
              x = paste0(input$pickFactor1DAPC, " (", tabVarprop$pct[tabVarprop$PC == input$pickFactor1DAPC], "% variance)"),
              y = paste0(input$pickFactor2DAPC, " (", tabVarprop$pct[tabVarprop$PC == input$pickFactor2DAPC], "% variance)"),
              color = input$colorGroupDAPC, shape = input$shapeGroupDAPC
            ) +
              theme_bw() +
              theme(
                axis.text.x = element_text(
                  colour = "grey20", size = input$textSizeDAPC, angle = 0, hjust = .5,
                  vjust = .5, face = "plain"
                ),
                axis.text.y = element_text(
                  colour = "grey20", size = input$textSizeDAPC, angle = 0, hjust = 1,
                  vjust = 0.5, face = "plain"
                ),
                axis.title.x = element_text(
                  colour = "grey20", size = input$textSizeDAPC, angle = 0, hjust = .5,
                  vjust = 0, face = "plain"
                ),
                axis.title.y = element_text(
                  colour = "grey20", size = input$textSizeDAPC, angle = 90,
                  hjust = .5, vjust = .5, face = "plain"
                ),
                legend.text = element_text(
                  colour = "grey20", size = input$textSizeDAPC
                ),
                legend.title = element_text(
                  colour = "grey20", size = input$textSizeDAPC
                ),
                title = element_text(colour = "grey20", size = input$textSizeDAPC)
              )
            
            
            plotDAPC <- plotDAPC +
              ylim(c(min(DAPCTable[[input$pickFactor2DAPC]] * 1.1), max(DAPCTable[[input$pickFactor2DAPC]] * 1.2))) +
              xlim(c(min(DAPCTable[[input$pickFactor1DAPC]] * 1.1), max(DAPCTable[[input$pickFactor1DAPC]] * 1.2)))
            
            output$DAPCStatic <- renderPlot(
              {
                plotDAPC
              },
              # width = as.numeric(input$DAPCPlotWidth),
              # height = as.numeric(input$DAPCPlotHeight)
              # , res = 96
            )
            
            # brushedPointsReactiveDAPC <- reactive(brushedPointsDAPC)
            
            # output$brush_info <- renderTable({
            #   brushedPoints(df = DAPCTable, brush = input$DAPCBrush, xvar = input$pickFactor1DAPC, yvar = input$pickFactor2DAPC)[, c(colnames(groupingVariables), input$pickFactor1DAPC, input$pickFactor2DAPC)]
            # })
            
            output$DAPCScree <- renderPlot({
              tabVarprop2 <- tabVarprop
              # pc_eigenvalues2$PC <- gsub("PC", "", pc_eigenvalues2$PC)
              # pc_eigenvalues2$PC <- factor(pc_eigenvalues2$PC, levels = pc_eigenvalues2$PC)
              if (nrow(tabVarprop2) > 15) {
                tabVarprop2 <- tabVarprop2[1:15,]
              }
              DAPCScree <- tabVarprop2 %>%
                ggplot(aes(x = PC)) +
                geom_col(aes(y = pct)) +
                geom_line(aes(y = pct_cum, group = 1)) +
                geom_point(aes(y = pct_cum)) +
                labs(x = "Principal component", y = "Fraction variance explained") +
                theme_classic(base_size = as.numeric(input$textSizeDAPC))
              DAPCScree
              
            }, width = 400, height = 400)
            
            output$downloadDAPC <- downloadHandler(
              filename = function() {
                paste0(Sys.Date(), "_DAPC.pdf")
              },
              content = function(file) {
                ggsave(
                  filename = file, plot = plotDAPC # ,
                  # width = (as.numeric(input$figWidthDAPC) / 3.2),
                  # height = (as.numeric(input$figHeightDAPC) / 3.2), limitsize = FALSE,
                  # units = "mm"
                )
              }
            )
            
            output$DAPCLoadings <- DT::renderDataTable({
              datatable(top_genes, rownames = F) %>% formatRound("loading", digits = 3)
            })
            
            
            # create id for data
            brushedPointsDAPC <- brushedPoints(df = DAPCTable, brush = input$DAPCBrush, xvar = input$pickFactor1DAPC, yvar = input$pickFactor2DAPC, allRows = TRUE)
            # brushedPointsID <- which(brushedPointsDAPC$selected_ == FALSE)
            
            output$brush_info <- renderTable({
              # brushedPoints(df = DAPCTable, brush = input$DAPCBrush, xvar = input$pickFactor1DAPC, yvar = input$pickFactor2DAPC)[, c(colnames(groupingVariables))]
              rownames(brushedPointsDAPC[brushedPointsDAPC$selected_ == TRUE, ])
            })
            
            
            # }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")} # close withProgress
            # }) # close tryCatch
          }
        ) # close Event number 2
      }
    ) # close Event number 1
  }) # close withProgress
}) # close observe
