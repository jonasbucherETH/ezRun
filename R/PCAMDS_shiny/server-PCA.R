observe({
  # withProgress(message = "...", {
    
    vcf <- inputDataReactive()$vcf
    genind <- inputDataReactive()$genind
    X <- inputDataReactive()$X
    grouping_vars <- inputDataReactive()$grouping_vars

    #####-----#####
    # Use observeEvent whenever you want to perform an action in response to an event.
    # (Note that "recalculate a value" does not generally count as performing an action--see eventReactive for that.)
    # The first argument is the event you want to respond to,
    # and the second argument is a function that should be called whenever the event occurs.
    
    # The observeEvent will only be dependent on the 'event' section in the small piece of code above. (first {...})
    # It will not be dependent on anything in the 'code to run' part. (after ignoreInit)
    #####-----#####
    
    # cat("before Event 1")
    
    # updateSelectInput(session = session, "color_group", choices = colnames(grouping_vars), selected = "") # try: selected = colnames(grouping_vars)[2] ([1] is sample name)
    # updateSelectInput(session = session, "shape_group", choices = colnames(grouping_vars), selected = "")
    
    observeEvent( # Event number 1
      {
        input$sample_labels
        input$color_group
        input$shape_group
      }, ignoreInit = T, # If TRUE, then, when the eventified object is first created/initialized, don't trigger the action or (compute the value). The default is FALSE.
      ignoreNULL = F,
      {
 
        ##### compute PCA
        pca <- dudi.pca(X, center = TRUE, scale = TRUE, scan = FALSE, nf = 5)
        n_pcs <- pca$nf # number of (> 0) principal components
        n_grouping <- ncol(grouping_vars)
        eig_sum <- sum(pca$eig)
        pca_varprop <- pca$eig/eig_sum
        pca_varprop <- pca_varprop[1:n_pcs]
        pca_tab <- data.frame(grouping_vars, pca$li, stringsAsFactors = FALSE, row.names = NULL)
        tab_varprop <- as.data.frame(t(pca_varprop), stringsAsFactors = FALSE)
        PC_indeces <- seq(1+n_grouping, ncol(pca_tab))
        for (i in 1:n_pcs){
          colnames(pca_tab)[i+n_grouping] <- paste0("PC", i)
          colnames(tab_varprop)[i] <- paste0("PC", i)}
        
        # all PCs in array for selecting input
        pc_list <- colnames(pca_tab)[-(1:n_grouping)]

        ##### end of PCA computation/prep
        
        ### TODO: put color_group default as 2nd variable in grouping_vars / pca_tab; make it necessary to have (or extend ifelse statement further down)
        updateSelectInput(session = session, "color_group", choices = colnames(grouping_vars), selected = colnames(grouping_vars)[2]) # try: selected = colnames(grouping_vars)[2] ([1] is sample name)
        updateSelectInput(session = session, "shape_group", choices = colnames(grouping_vars), selected = "")
        
        # levels(as.factor(pca_tab[[input$shape_group]]))
        
        # # Can use character(0) to remove all choices
        # if (is.null(x))
        #   x <- character(0)
        
        updateSelectInput(session = session, "pick_pc_x", choices = pc_list, selected = "PC1")
        updateSelectInput(session = session, "pick_pc_y", choices = pc_list, selected = "PC2")
        
        # cat("before bindEvent 2")
        
        observeEvent({ # Event number 2: only need axis inputs here, for the others no need for redrawing plot (?)
          input$pick_pc_x
          input$pick_pc_y
          # input$color_group
          # input$shape_group
        }, {
          
          # pca_ggplot <- pca_tab %>%
          #   ggplot(aes_string(x = input$pick_pc_x, y = input$pick_pc_y)) +
          #   geom_point()
          
          # pca_ggplot <- pca_tab %>%
          #   ggplot(aes(x = .data[[input$pick_pc_x]], y = .data[[input$pick_pc_y]])) +
          #   geom_point()
          
          ## TODO: replace by .data
          # if (input$shape_group == "None") {
          #   pca_ggplot <- pca_tab %>%
          #     ggplot(aes_string(x = input$pick_pc_x, y = input$pick_pc_y, color = input$color_group, fill = input$color_group)) +
          #     geom_point()
          # } else {
          #   pca_ggplot <- pca_tab %>%
          #     ggplot(aes_string(x = input$pick_pc_x, y = input$pick_pc_y, color = input$color_group, fill = input$color_group, shape = input$shape_group)) +
          #     geom_point()
          #     # scale_shape_manual(values = c(rep(c(21,22,23,24,25,8,3,4), times = 10 ))[1:nlevels(as.factor(datasetPCA[[input$pcaFactor2]]))] )
          #     ### Note: 21-25 are simple shapes, but with filling; 8,3,4 are crosses with different amount of lines
          # }
          
          ### TODO: replace deprecated aes_string (above) with .data[[]]
          if (input$shape_group == "None") {
            pca_ggplot <- pca_tab %>%
              ggplot(aes(x = .data[[input$pick_pc_x]], y = .data[[input$pick_pc_y]], color = .data[[input$color_group]], fill = .data[[input$color_group]])) +
              geom_point()
          } else {
            pca_ggplot <- pca_tab %>%
              ggplot(aes(x = .data[[input$pick_pc_x]], y = .data[[input$pick_pc_y]], color = .data[[input$color_group]], fill = .data[[input$color_group]], shape = .data[[input$shape_group]])) +
              geom_point()
              # scale_shape_manual(values = c(rep(c(21,22,23,24,25,8,3,4), times = 10 ))[1:nlevels(as.factor(datasetPCA[[input$pcaFactor2]]))] )
              ### Note: 21-25 are simple shapes, but with filling; 8,3,4 are crosses with different amount of lines
          }
          
          # if (input$sample_labels) {
          #   pca_ggplot <- pca_ggplot + 
          #     geom_text(
          #       aes_string(label = grouping_vars[,1]),
          #       # size = (input$textSizePCA/3),
          #       hjust = 0.2, vjust = -1.5, check_overlap = T,
          #       show.legend  = FALSE
          #     ) +
          #     ylim(c(min(pca_tab[[input$pick_pc_y]]*1.1), max(pca_tab[[input$pick_pc_y]]*1.2))) +
          #     xlim(c(min(pca_tab[[input$pick_pc_x]]*1.1), max(pca_tab[[input$pick_pc_x]]*1.2)))
          #   
          # }
          
          
          if (input$sample_labels) {
            pca_ggplot <- pca_ggplot +
              geom_text(
                aes(label = grouping_vars[,1]),
                # size = (input$textSizePCA/3),
                hjust = 0.2, vjust = -1.5, check_overlap = T,
                show.legend  = FALSE
              )
          }
          
          ### TODO: add this to ui
          # if (input$show_proportion) {
          #   pca_ggplot <- pca_ggplot + 
          #       xlab(paste0("PC1 (", format(round(pca_varprop[1]*100, 1), nsmall = 1), "%)" )) +
          #       ylab(paste0("PC2 (", format(round(pca_varprop[2]*100, 1), nsmall = 1), "%)" )) 
          # }
          
          
          ### themes and stuff
          pca_ggplot <- pca_ggplot +
            theme_classic()

          output$pca_plot <- renderPlot({
            pca_ggplot
          }, 
          # width = as.numeric(input$pcaPlotWidth),
          # height = as.numeric(input$pcaPlotHeight),
          # res = 96
          )

          # output$result <- renderText({
          #   paste("You chose", input$pick_pc_x)
          # })

          # pca_ggplot <- ggplot(data = pca_tab, aes(x = .data[[input$pick_pc_x]], y = .data[[input$pick_pc_y]])) +
          # pca_ggplot <- ggplot(data = pca_tab, aes_string(x = input$pick_pc_x, y = input$pick_pc_y)) +
          # pca_ggplot <- pca_tab %>%
          #   ggplot(aes_string(x = input$pick_pc_x, y = input$pick_pc_y)) +
          #     geom_point() +
          #     theme_classic()
        # output$pca_plot <- renderPlot({
        #       pca_ggplot
        # })

          # test_table <- xtable::xtable(pca_tab)

          # test_table <- matrix(data = NA, nrow = 2, ncol = 2, byrow = FALSE,
          #                      dimnames = NULL)
          #
          # test_table <- xtable::xtable(test_table)
          #
          #
          # output$pca_table <- renderTable(test_table)
          # output$pca_table <- renderTable(pca_tab)


        # }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")} # close withProgress
        }) # close bindEvent number 2



    }) # close bindEvent number 1
  # }) # close withProgress
}) # close observe

    
    

