observe({
  withProgress(message = "Generating PCA Plots. Please wait...", {
    pca <- inputDataReactive()$pca
    n_pcs <- inputDataReactive()$n_pcs
    pca_varprop <- inputDataReactive()$pca_varprop
    pca_tab <- inputDataReactive()$pca_tab
    tab_varprop <- inputDataReactive()$tab_varprop
    
    vcf <- inputDataReactive()$vcf
    genind <- inputDataReactive()$genind
    grouping_vars <- inputDataReactive()$grouping_vars
    X <- inputDataReactive()$X
    
    # for testing
    # vcf <- read.vcfR("~/sushi_project_JB/data/test_vcf_dataset/ragi_highcov_sa0001_1k.vcf.gz")
    # genind <- vcfR2genind(vcf)
    # grouping_vars <- read.delim("~/sushi_project_JB/data/test_vcf_dataset/populations.txt")
    # # pop(genind) <- populations_txt$Population
    # 
    # X <- scaleGen(genind, NA.method="mean")
    
    observeEvent( # number 1
      {
        input$pick_pc_x
        input$pick_pc_y
        input$color_group
        input$shape_group
        input$sample_labels
        
      }, ignoreInit = TRUE,
      {
        
        tryCatch({
          
          # vcf_new <- vcf
          # genind_new <- genind
          # grouping_vars_new <- grouping_vars
          # X_new <- X
          # pca <- dudi.pca(X_new, center = TRUE, scale = TRUE, scan = FALSE, nf = 5)
          # n_pcs <- pca$nf # number of (> 0) principal components
          # eig_sum <- sum(pca$eig)
          # pca_varprop <- pca$eig/eig_sum
          # pca_varprop <- pca_varprop[1:n_pcs]
          # pca_tab <- data.frame(grouping_vars_new, pca$li, stringsAsFactors = FALSE, row.names = NULL)
          # tab_varprop <- as.data.frame(t(pca_varprop), stringsAsFactors = FALSE)
          # PC_indeces <- seq(1+ncol(grouping_vars), ncol(tab))
          # ##
          # for (i in 1:n_pcs){
          #   colnames(pca_tab)[i+ncol(grouping_vars_new)] <- paste0("PC", i)
          #   colnames(tab_varprop)[i] <- paste0("PC", i)
          # }
          # ##
          
          
          # # pca <- dudi.pca(X, center = TRUE, scale = TRUE, scan = FALSE, nf = 5)
          # n_pcs <- pca$nf # number of (> 0) principal components
          # eig_sum <- sum(pca$eig)
          # pca_varprop <- pca$eig/eig_sum
          # pca_varprop <- pca_varprop[1:n_pcs]
          # pca_tab <- data.frame(grouping_vars, pca$li, stringsAsFactors = FALSE, row.names = NULL)
          # tab_varprop <- as.data.frame(t(pca_varprop), stringsAsFactors = FALSE)
          # PC_indeces <- seq(1+ncol(grouping_vars), ncol(tab))
          # ##
          # for (i in 1:n_pcs){
          #   colnames(pca_tab)[i+ncol(grouping_vars)] <- paste0("PC", i)
          #   colnames(tab_varprop)[i] <- paste0("PC", i)
          # }
          # ##
          
          # },
          # 
          # {
          #   ### maybe try updateVarSelectInput (with varSelectInput in ui) 
          #   updateSelectInput(session = "pick_pc_x", inputId = "pick_pc_x", "Select PC for x-axis", choices = colnames(pca_tab)[PC_indeces], selected = colnames(pca_tab)[PC_indeces[1]])
          #   updateSelectInput(session = "pick_pc_y", inputId = "pick_pc_y", "Select PC for x-axis", choices = colnames(pca_tab)[PC_indeces], selected = colnames(pca_tab)[PC_indeces[2]])
          #   
          #   # maybe try PC_indeces[1] and PC_indeces[2] for selected
          #   # updateVarSelectInput(session = session, "pick_pc_x", "Select PC for x-axis", data = pca_tab[, PC_indeces], selected = "PC1")
          #   # updateVarSelectInput(session = session, "pick_pc_x", "Select PC for x-axis", data = pca_tab[, PC_indeces], selected = "PC2")
          #   
          # },
          
          { # AAA
            observeEvent( # number 2
              { 
                input$pick_pc_x
                input$pick_pc_y
                input$color_group
                input$shape_group
              },
              
              {
                # ggplot(data = values$pca, aes(x = .data[[input$pick_pc_x]], y = .data[[input$pick_pc_y]], color=.data[[input$color_group]], shape=.data[[input$shape_group]])) +
                # pca_ggplot <- ggplot(data = pca_tab, aes(x = .data[[input$pick_pc_x]], y = .data[[input$pick_pc_y]], color=!!input$color_group, shape=!!input$shape_group)) +
                #         geom_point() +
                #         # geom_text(hjust=0, vjust=0, label=values$pca$sample.id, label.size=0.1) +
                #         theme_classic() +
                #         xlab(paste0("PC1 (", format(round(tab_varprop[input$pick_pc_x]*100, 1), nsmall = 1), "%)" )) +
                #         ylab(paste0("PC2 (", format(round(tab_varprop[input$pick_pc_y]*100, 1), nsmall = 1), "%)" )) +
                #         labs(color = input$color_by, shape = input$group_by)
                pca_ggplot <- ggplot(data = pca_tab, aes(x = .data[[input$pick_pc_x]], y = .data[[input$pick_pc_y]])) +
                  # pca_ggplot <- ggplot(data = pca_tab, aes(x = !!input$pick_pc_x, y = !!input$pick_pc_y)) +
                  # pca_ggplot <- ggplot(aes_string(x = input$pick_pc_x, y = input$pick_pc_y)) +
                  
                  geom_point() +
                  # geom_text(hjust=0, vjust=0, label=values$pca$sample.id, label.size=0.1) +
                  theme_classic()
                
              },
              
              # see line 141 
              # if (input$show_labels) {
              #   pca_plot
              # }
              
              {
                output$pca_plot <- renderPlot(
                  {
                    pca_ggplot
                  }
                  , height = 400, width = 400
                  # width = as.numeric(input$figWidthPCA),
                  # height = as.numeric(input$figHeightPCA)
                )
              },
              
              {
                test_table <- xtable::xtable(pca_tab)
              },
              
              {
                output$pca_table <- renderTable(test_table)
              }
              
            ) # close observeEvent number 2
          } # AAA
          
          , error=function(e){cat("ERROR :",conditionMessage(e), "\n")} 
        )}# close tryCatch
        
    ) # close observeEvent number 1
      }) # close withProgress
  }) # close observe




