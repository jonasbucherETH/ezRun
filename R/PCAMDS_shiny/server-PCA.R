observe({
  withProgress(message = "Generating PCA Plots. Please wait...", {
    vcfRaw <- inputDataReactive()$vcfRaw
    vcfGenind <- inputDataReactive()$vcfGenind
    datasetScaled <- inputDataReactive()$datasetScaled
    groupingVariables <- inputDataReactive()$groupingVariables
    colourList <- inputDataReactive()$colourList
    
    # shinyInput <- function(FUN, len, id, ...) {
    #   inputs <- character(len)
    #   for (i in seq_len(len)) {
    #     inputs[i] <- as.character(FUN(paste0(id, i), label = NULL, ...))
    #   }
    #   inputs
    # }

    
    groupingVarsCheckbox <- cbind(groupingVariables, bool = T)
    
    # groupingVarsCheckbox2 <- cbind(
    #   groupingVariables,
    #   Include = shinyInput(checkboxInput, nrow(groupingVarsCheckbox), "checkb")
    # )
    
    # js <- c(
    #   "$('[id^=checkb]').on('click', function(){",
    #   "  var id = this.getAttribute('id');",
    #   "  var i = parseInt(/checkb(\\d+)/.exec(id)[1]);",
    #   "  var value = $(this).prop('checked');",
    #   "  var info = [{row: i, col: 3, value: value}];",
    #   "  Shiny.setInputValue('sampleTablePCA_cell_edit:DT.cellInfo', info);",
    #   "})"
    # )
    
    # output$sampleTablePCA <- DT::renderDataTable({
    #   datatable(groupingVarsCheckbox2,
    #             rownames = F,
    #             escape = FALSE,
    #             editable = list(target = "column", disable = list(columns = 3)),
    #             selection = "none",
    #             callback = JS(js)
    #   )
    # }, server = FALSE
    # )

    # datasetScaledPCA <- as.data.frame(datasetScaled)
    datasetScaledPCA <- datasetScaled
    
    # datasetScaledPCA$selected_ <- rep(FALSE, nrow(datasetScaledPCA))
    # selectedList <- c()
    
    # output$x1 = DT::renderDataTable(groupingVariables, server = FALSE)
    
    # highlight selected rows in the scatterplot
    # output$x2 = renderPlot({
    #   s = input$x1_rows_selected
    #   # par(mar = c(4, 4, 1, .1))
    #   plot(groupingVariables)
    #   # if (length(s)) points(groupingVariables[s, , drop = FALSE], pch = 19, cex = 2)
    # })
    
    # server-side processing
    # groupingVariables = groupingVariables[, 1:2]
    output$sampleTablePCA <- DT::renderDataTable(
      groupingVariables,
      server = TRUE,
      rownames = FALSE,
      filter = "bottom"
    )
    
    # print the selected indices
    output$selectedSamplesPCA <- renderPrint({
      s = input$sampleTablePCA_rows_selected
      if (length(s)) {
        cat('These samples are currently selected:\n\n')
        cat(rownames(datasetScaled)[s], sep = '\n')
      }
    })
    
    # output$selectedSamplesPCA <- renderText({
    #   s = input$sampleTablePCA_rows_selected
    #   # selectedSamplesText <- ""
    #   if (length(s)) {
    #     cat('These samples are currently selected:\n\n')
    #     cat(rownames(datasetScaled)[s], sep = '\n')
    #   }
    # })
    
    # output$cbTablePCA = DT::renderDataTable(groupingVariables, server = TRUE)

    observeEvent( # Event number 1
      {
        input$removeSamplesPCA
        input$resetSelectionPCA
        # input$pcaBrush
        # input$sampleLabelsPCA
        # input$colorGroupPCA
        # input$shapeGroupPCA
        # input$pickFactor1PCA
        # input$pickFactor2PCA
        # input$pcaPlotWidth
        # input$pcaPlotHeight
        # input[["dtable_cell_edit"]]
        # input$sampleTablePCA_cell_edit
        # input$checkboxPCA
        # input$choice
        
      },
      ignoreInit = F, # If TRUE, then, when the eventified object is first created/initialized, don't trigger the action or (compute the value). The default is FALSE.
      ignoreNULL = F, # default = TRUE
      {
        
        # removeRowsPCA <- input$x3_rows_selected
        # if (is.null(removeRowsPCA)) {
        #   removeRowsPCA <- c()
        # }
        # 
        # a <- c(1,3,6)
        # b <- 1:nrow(groupingVariables)
        
        # print(length(input$x3_rows_selected))
        
        removeRowsPCA <- input$sampleTablePCA_rows_selected
        if (length(removeRowsPCA)) {
          # cat('These rows were selected:\n\n')
          s <- c(removeRowsPCA)
          # print(s)
          # print(groupingVariables[-s, ])
          datasetScaledPCA <- datasetScaled[-s, ]
          groupingVariables <- groupingVariables[-s, ]
        }

        # print(removeRowsPCA)


        # datasetScaledPCA <- datasetScaled[-removeRowsPCA, ]
        
        # cbTablePCA <- checkboxGroupTable(
        #   tbl = groupingVariables,
        #   inputId = "checkboxPCA", # for checkboxGroupInput
        #   label = "",
        #   choices = rownames(groupingVariables),
        #   selected = rownames(groupingVariables),
        #   table_label = "",
        #   control_column = 3L,
        #   pixie = . %>% identity(),
        #   # pixie = . %>% 
        #   # sprinkle(bg_pattern_by = "rows") %>%
        #   # sprinkle_table(pad = 7) %>%
        #   # sprinkle_colnames("colnames(groupingVariables)" = "",
        #   #                   control = ""),
        #   display_table = FALSE,
        #   disabled = FALSE,
        #   hidden = FALSE,
        #   disabled_table = FALSE,
        #   hidden_table = FALSE
        # )
        
        # output$checkboxTablePCA <- renderText({
        #   cbTablePCA
        # })
        # 
        # output$choice <- renderText(input$checkboxPCA)
        # 
        # choiceVector <- input$checkboxPCA
        # print(choiceVector[1])
        # print(choiceVector)
        # print(cbTablePCA[1,3])
        # cat(input$checkboxPCA)
        
        # # # TODO: omit selected samples/groups

        # output$sampleTablePCA <- DT::renderDataTable({
        #   datatable(groupingVarsCheckbox2,
        #             rownames = F,
        #             escape = FALSE,
        #             editable = list(target = "cell", disable = list(columns = 3)),
        #             selection = "none",
        #             callback = JS(js)
        #   )
        # }, server = FALSE
        # )
        

        # info <- input[["dtable_cell_edit"]] # this input contains the info of the edit
        # info <- input$sampleTablePCA_cell_edit
        # info <- input[["sampleTablePCA_cell_edit"]]
        
        # groupingVarsCheckbox2[input$sampleTablePCA_cell_edit$row,input$sampleTablePCA_cell_edit$col] <<- input$sampleTablePCA_cell_edit$value
        
        # print(groupingVarsCheckbox2$Include)
        
        # datasetScaledPCA <- editData(datasetScaledPCA, info)
        
        # print(datasetScaledPCA[1,])
        
        # print(datasetScaledPCA[info$row, 1])
        # if (!is.null(nrow(info))) {
        #   datasetScaledPCA <- as.data.frame(datasetScaled)
        #   datasetScaledPCA <- datasetScaledPCA[info$row[info$value==FALSE], ]
        #   groupingVariables <- groupingVariables[info$row[info$value==FALSE], ]
        # 
        # }
        
        # datasetScaledPCA <- as.data.frame(datasetScaled[-info, ])
        
        
        # datasetScaledPCA <- boolDatasetScaledPCA[-info$row, -ncol(boolDatasetScaledPCA)]
        
        # print(info$row)
        # print(info)
        
        # datasetScaled <<- editData(datasetScaled, input$sampleTablePCA_cell_edit, 'sampleTablePCA')
        # datasetScaled <<- datasetScaled[-info$row, ]
        
        # Dat(editData(Dat(), info))

        
        # if (length(selectedList)!=0) {
        #   datasetScaledPCA <- datasetScaledPCA[!(rownames(datasetScaledPCA) %in% selectedList), ]
        #   groupingVariables <- groupingVariables[!(rownames(groupingVariables) %in% selectedList), ]
        # }
          
        # print("up")
        # print(datasetScaledPCA$selected_)
        
        # print(rownames(datasetScaledPCA[!datasetScaledPCA$selected_, ]))
        
        # datasetScaledPCA <- datasetScaledPCA[!datasetScaledPCA$selected_, ]
        # print(selectedList)
        
        # cat("1")
        
        # datasetScaledPCA <- as.data.frame(datasetScaled)
        
        # create id for data
        # if (input$removeThesePointsPCA!="") {
        #   # cat("-------- in if clause -----------")
        #   
        #   # cat(nrow(datasetScaledPCA[rownames(datasetScaledPCA)!=input$removeThesePointsPCA, ]))
        #   # print(input$pcaBrush)
        #   # currentBrush <- brushedPoints(df = pcaTable, brush = input$pcaBrush, xvar = input$pickFactor1PCA, yvar = input$pickFactor2PCA, allRows = T)
        #   # print(currentBrush)
        #   datasetScaledPCA <- datasetScaledPCA[rownames(datasetScaledPCA)!=input$removeThesePointsPCA, ]
        #   groupingVariables <- groupingVariables[rownames(groupingVariables)!=input$removeThesePointsPCA, ]
        # }
        
        # cat(rownames(datasetScaled))
        # cat(nrow(datasetScaled))
        # cat(nrow(datasetScaledPCA))
        # cat(nrow(groupingVariables))
        
        # cat(class(datasetScaled))
        # cat(class(datasetScaledPCA))
        # cat(class(groupingVariables))
        
        # datasetScaled2$id <- 1:nrow(datasetScaled)
        # brushedPointsPCA <- brushedPoints(df = pcaTable, brush = input$pcaBrush, xvar = input$pickFactor1PCA, yvar = input$pickFactor2PCA, allRows = TRUE)
        

        
        # brushedPointsPCA <- brushedPoints(df = pcaTable, brush = input$pcaBrush, xvar = input$pickFactor1PCA, yvar = input$pickFactor2PCA, allRows = TRUE)
        
        # Return the brushed points. See ?shiny::brushedPoints.
        
        # p <- brushedPoints(data, input$brush)
        
        # create vector of ids that match brushed points and data
        # g <- which(brushedPointsID %in% datasetScaled2$id)
        
        # return a subset of the original data without brushed points
        # datasetScaled2 <- datasetScaled2[g,]

        
        # if (exists("brushedPointsPCA")) {
        #   print(brushedPointsPCA)
        # } else {
        #   print("brushed points do not exits/cannot be found --------------------------")
        # }
        # selectedPointsPCA <- rownames(brushedPoints(df = pcaTable, brush = input$pcaBrush, xvar = input$pickFactor1PCA, yvar = input$pickFactor2PCA, allRows = FALSE)) #[, c(colnames(groupingVariables))]
        # brushedPointsPCA <- brushedPoints(df = pcaTable, brush = input$pcaBrush, xvar = input$pickFactor1PCA, yvar = input$pickFactor2PCA, allRows = TRUE)
        # datasetScaled2 <- datasetScaled
        # datasetScaled <- datasetScaled[!(rownames(datasetScaled2) %in% input$removeThesePointsPCA), ]
        
        
        # if (exists("brushedPointsReactivePCA")) {
        #   selectedPointsPCA <- brushedPointsReactivePCA
        #   selectedPointsPCA <- as.data.frame(selectedPointsPCA, row.names = rownames(groupingVariables))
        #   print(selectedPointsPCA)
        # } else {
        #   print("brushed points do not exits/cannot be found --------------------------")
        # }

        # selectedPointsPCA <- selectedPointsPCA[selectedPointsPCA$selected_ == TRUE, ]
        # selectedPointsPCA <- which(selectedPointsPCA$selected_=="TRUE")
        # selectedPointsPCA <- selectedPointsPCA[selectedPointsPCA$selected_ == TRUE, ]
        
        # if (is.null(input$excludePointsPCA)) {
        #   print("nll")
        # }
        
        # observeEvent(input$DeleteSelectedData, {
        #   Var1 <- brushedPoints(rx_de(), input$plot1_brush, allRows = TRUE)
        #   rx_de(Var1[!Var1$selected_, names(Var1) != "selected_", drop = FALSE])
        # })
        # 
        # print(selectedPointsPCA)
        
        # print(selectedPointsPCA)
        # Var1 <- brushedPoints(datasetScaled, input$pcaBrush, allRows = TRUE)
        # datasetScaled(Var1[!Var1$selected_, names(Var1) != "selected_", drop = FALSE])
        # selectedPointsPCA <- brushedPoints(d, input$plot1_brush, allRows = TRUE)
        # datasetScaled <- datasetScaled

        ##### compute PCA
        pcaResults <- dudi.pca(datasetScaledPCA, center = TRUE, scale = TRUE, scan = FALSE, nf = 5)

        
        nPC <- pcaResults$nf # number of (> 0) principal components
        nGrouping <- ncol(groupingVariables)
        pcaVarprop <- pcaResults$eig

        pcaVarprop <- tibble(PC = paste0("PC", factor(1:length(pcaVarprop))), variance = pcaVarprop) %>% 
          mutate(pct = format(variance/sum(variance)*100, digits = 2)) %>%
          # mutate(pct = variance/sum(variance)*100) %>% 
          mutate(pct_cum = cumsum(pct))
        pcaVarprop$PC <- factor(pcaVarprop$PC, levels = pcaVarprop$PC)
        

        pcaTable <- data.frame(groupingVariables, pcaResults$li, stringsAsFactors = FALSE, row.names = rownames(datasetScaledPCA))

        tabVarprop <- pcaVarprop
        for (i in 1:nPC) {
          colnames(pcaTable)[i + nGrouping] <- paste0("PC", i)
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
            input$pcaBrush
            # input$excludeSamplesPCA
          },
            ignoreInit = F, # default = FALSE
            ignoreNULL = F, # default = TRUE
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
            
            if (isTRUE(input$pcaAxesProp)) {
              plotPCA <- plotPCA + coord_fixed(ratio = as.numeric(tabVarprop$pct[tabVarprop$PC == input$pickFactor2PCA])/as.numeric(tabVarprop$pct[tabVarprop$PC == input$pickFactor1PCA]))
            }
            

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
              x = paste0(input$pickFactor1PCA, " (", tabVarprop$pct[tabVarprop$PC == input$pickFactor1PCA], "% variance explained)"),
              y = paste0(input$pickFactor2PCA, " (", tabVarprop$pct[tabVarprop$PC == input$pickFactor2PCA], "% variance explained)"),
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
              width = as.numeric(input$pcaPlotWidth),
              height = as.numeric(input$pcaPlotHeight)
              # , res = 96
            )
            
            # create id for data

            # selectedList <<- append(selectedList, which(brushedPointsPCA$selected_))
            
            # datasetScaledPCA$selected_ <<- brushedPointsPCA$selected_
            
            # print("down")
            # print(datasetScaledPCA$selected_)
            
            # cat(datasetScaledPCA$selected_)
            # brushedPointsID <- which(brushedPointsPCA$selected_ == FALSE)

            # output$brushInfo <- renderTable({
            #   # brushedPoints(df = pcaTable, brush = input$pcaBrush, xvar = input$pickFactor1PCA, yvar = input$pickFactor2PCA)[, c(colnames(groupingVariables))]
            #   rownames(brushedPointsPCA[brushedPointsPCA$selected_ == TRUE, ])
            # },
            # rownames = F,
            # colnames = F
            # )

            # brushedPointsReactivePCA <- reactive(brushedPointsPCA)
            
            # output$brush_info <- renderTable({
            #   brushedPoints(df = pcaTable, brush = input$pcaBrush, xvar = input$pickFactor1PCA, yvar = input$pickFactor2PCA)[, c(colnames(groupingVariables), input$pickFactor1PCA, input$pickFactor2PCA)]
            # })
            
            output$pcaScree <- renderPlot({
              tabVarprop2 <- tabVarprop
              tabVarprop2$PC <- gsub("PC", "", tabVarprop2$PC)
              tabVarprop2$PC <- factor(tabVarprop2$PC, levels = tabVarprop2$PC)
              tabVarprop2$variance <- as.numeric(tabVarprop2$variance)
              tabVarprop2$pct <- as.numeric(tabVarprop2$pct)
              tabVarprop2$pct_cum <- as.numeric(tabVarprop2$pct_cum)
              
              if (nrow(tabVarprop2) > 10) {
                tabVarprop2 <- tabVarprop2[1:10,]
              }
              pcaScree <- tabVarprop2 %>%
                  ggplot(aes(x = PC)) +
                  geom_col(aes(y = pct)) +
                  geom_line(aes(y = pct_cum, group = 1)) +
                  geom_point(aes(y = pct_cum)) +
                  labs(x = "Principal component", y = "Fraction variance explained (%)") +
                  scale_y_continuous(n.breaks = 20) +
                  theme_classic(base_size = as.numeric(input$textSizePCA))
              pcaScree
            
            } 
            # , width = 500, height = 400
            )
            
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
