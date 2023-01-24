observe({
  # withProgress(message = "Generating PCA Plots. Please wait...", {
    # pca <- inputDataReactive()$pca
    # n_pcs <- inputDataReactive()$n_pcs
    # pca_varprop <- inputDataReactive()$pca_varprop
    # pca_tab <- inputDataReactive()$pca_tab
    # tab_varprop <- inputDataReactive()$tab_varprop
    
    vcf <- inputDataReactive()$vcf
    genind <- inputDataReactive()$genind
    grouping_vars <- inputDataReactive()$grouping_vars
    X <- inputDataReactive()$X
    
    cat("abc")
    # expect_match(grouping_vars[1,1], "Eafr_Afr_GE7123")
    # expect_match(grouping_vars[1,1], "Eafr7123")
    
    # observeEvent(input$pick_pc_x, {
    #   print(paste0("Test: ", input$pick_pc_x))
    # })
 
    # n_grouping <- inputDataReactive()$n_grouping
    # pc_list <- inputDataReactive()$pc_list
    
    output$pca_table <- renderTable("pca_tab")
    
    output$test_text <- renderText({n_pcs})
    
    # observeEvent({
    #   updateSelectInput(session = session, inputId = "pick_pc_x", selected = "PC1")
    #   updateSelectInput(session = session, inputId = "pick_pc_y", selected = "PC2")
    # })
    
    # updateVarSelectInput(session = session, "pick_pc_x", data = pca_tab[, PC_indeces], selected = "PC1")
    # updateVarSelectInput(session = session, "pick_pc_x", data = pca_tab[, PC_indeces], selected = "PC2")
    
    # observeEvent({
    #   updateSelectInput(session = session, "pick_pc_x", choices = pc_list, selected = "PC1")
    #   updateSelectInput(session = session, "pick_pc_x", choices = pc_list, selected = "PC2") 
    # })
    
    # updateSelectInput(session = session, "pick_pc_x", choices = pc_list, selected = "PC1")
    # updateSelectInput(session = session, "pick_pc_x", choices = pc_list, selected = "PC2")

    
    # for testing
    # vcf <- read.vcfR("~/sushi_project_JB/data/test_vcf_dataset/ragi_highcov_sa0001_1k.vcf.gz")
    # genind <- vcfR2genind(vcf)
    # grouping_vars <- read.delim("~/sushi_project_JB/data/test_vcf_dataset/populations.txt")
    # # pop(genind) <- populations_txt$Population
    # 
    # X <- scaleGen(genind, NA.method="mean")
    
    observeEvent( # number 1
      {
        # input$pick_pc_x
        # input$pick_pc_y
        # input$color_group
        # input$shape_group
        input$sample_labels
      }, ignoreInit = TRUE,
      {
        vcf <- vcf
        genind <- genind
        grouping_vars <- grouping_vars
        X <- X
        
        # compute PCA
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
          colnames(tab_varprop)[i] <- paste0("PC", i)
        }
        
        # all PCs in array for selecting input
        pc_list <- colnames(pca_tab)[-(1:n_grouping)]
        
        updateSelectInput(session = session, "pick_pc_x", choices = pc_list, selected = "PC1")
        updateSelectInput(session = session, "pick_pc_x", choices = pc_list, selected = "PC2")

      
        observeEvent({ # number 2
          input$pick_pc_x
          input$pick_pc_y
          input$color_group
          input$shape_group
        }, {
      # { # aaa
          
          output$result <- renderText({
            paste("You chose", input$pick_pc_x)
          })

          # pca_ggplot <- ggplot(data = pca_tab, aes(x = .data[[input$pick_pc_x]], y = .data[[input$pick_pc_y]])) +
          # pca_ggplot <- ggplot(data = pca_tab, aes_string(x = input$pick_pc_x, y = input$pick_pc_y)) +
          pca_ggplot <- pca_tab %>%
            ggplot(aes_string(x = input$pick_pc_x, y = input$pick_pc_y)) +
              geom_point() +
              theme_classic()

          output$pca_plot <- renderPlot({
                pca_ggplot
          })
          
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
        }) # close observeEvent number 2
      # } # aaa
        
      
    }) # close observeEvent number 1
  # }) # close withProgress
}) # close observe

    
    

