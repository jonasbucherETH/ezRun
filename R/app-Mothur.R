###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMothur = function(input=NA, output=NA, param=NA, 
                                     htmlFile="00index.html"){
    
  require(rmarkdown)
  require(ShortRead)
  require(phyloseq)
  require(plyr)
  require(ape)
  library(scales)
  dataset = input$meta
  sampleNames = input$getNames() 
  isPaired <- param$paired
  mockSample = param$referenceFasta != ""
  ### read fastq files and prepare inputs for Mothur
  ### File preparation: are reads paired? should they be joined? 
  file1PathInDataset <- input$getFullPaths("Read1")
  if(isPaired){
      file2PathInDataset <- input$getFullPaths("Read2")
      
      fastqJoinFun <- function(x,y,z){
          joinedFileName <- paste0(z, ".")
          fastqJoinCmd <- paste("fastq-join", x, y, "-o", joinedFileName)
          ezSystem(fastqJoinCmd)
          joinedFileName <- paste0(joinedFileName,"join")
          joinedFile <- file.path(getwd(),joinedFileName)
          outFileName <- paste0(z,".fastq")
          ezSystem(paste("mv ",joinedFileName,outFileName))
          return(outFileName)
      }
      listOfJoinedFiles <- mapply(fastqJoinFun,file1PathInDataset,file2PathInDataset,
                                  sampleNames)
  }else{
      listOfJoinedFiles <- vector()
      i=0
      for (singleEndFile in file1PathInDataset){
          i=i+1
          listOfJoinedFiles[i] <-  paste0(sampleNames[i],".fastq")
          cpCmd <- paste0("gunzip -c ", singleEndFile, "  > ", listOfJoinedFiles[i])
          ezSystem(cpCmd)
      }
  }
    
    ### create mothur input fasta and group files  
  for (sampleName in sampleNames){
      fastqFileToRead <- readDNAStringSet(paste0(sampleName,".fastq"),format = "fastq")
      readNames <- ldply(strsplit(names(fastqFileToRead)," "), function(x)x[1])$V1
      groupFile <- data.frame(readNames, sampleName, stringsAsFactors = F)
      groupFileName <- "Mothur.groups"
      write.table(groupFile,groupFileName, col.names = F, row.names = F, quote = F,append = TRUE)
      fastaOutName <- "Mothur.fasta"
      fastaFileToWrite <- writeXStringSet(fastqFileToRead, fastaOutName,append=TRUE)
  }
  
  # is there at least a mock sample for the error estimate? The error estimates for the Non-mock samples will be ignored downstream
  mockString = "###seq.error" 
    #  if(mockSample){
    #       copyRefCmd <- paste("cp", param$referenceFasta,"./", sep = " ")
    #       ezSystem(copyRefCmd)
    #       mockString = "seq.error" 
    #       refString = param$referenceFasta
    #      }
    #     write("There is no relevant info in this file",oldErrFile)
    ###
    
    ### cp silva reference locally
  cpSilvaRefCmd <- paste("gunzip -c ",SILVA_DB_MOTHUR_R138, "| tar xvf -")
  ezSystem(cpSilvaRefCmd)
    
    ### update batch file  with parameters and run mothur: step 1, identify region
  updateBatchPart1Cmd1 <- paste0("sed -e s/\"###seq.error\"/", mockString, "/g",
                            " -e s/\"CUTOFF_TAXON\"/", param$cutOffTaxonomy, "/g",
                            " -e s/\"CUTOFF_CLUST\"/", param$cutOffCluster, "/g",
                            " -e s/\"DIFFS\"/", param$diffs, "/g ",
                            file.path(METAGENOMICS_ROOT,UNIFIED_MOTHUR_WORKFLOW_PART1), 
                            " > ",
                            UNIFIED_MOTHUR_WORKFLOW_PART1)
  ezSystem(updateBatchPart1Cmd1)
  updateBatchPart2ACmd1 <- paste0("sed -e s/\"###seq.error\"/", mockString, "/g",
                                 " -e s/\"CUTOFF_TAXON\"/", param$cutOffTaxonomy, "/g",
                                 " -e s/\"CUTOFF_CLUST\"/", param$cutOffCluster, "/g",
                                 " -e s/\"DIFFS\"/", param$diffs, "/g ",
                                 file.path(METAGENOMICS_ROOT,UNIFIED_MOTHUR_WORKFLOW_PART2A), 
                                 " > ",
                                 UNIFIED_MOTHUR_WORKFLOW_PART2A)
  ezSystem(updateBatchPart2ACmd1)
  updateBatchPart2BCmd1 <- paste0("sed -e s/\"###seq.error\"/", mockString, "/g",
                                  " -e s/\"CUTOFF_TAXON\"/", param$cutOffTaxonomy, "/g",
                                  " -e s/\"CUTOFF_CLUST\"/", param$cutOffCluster, "/g",
                                  " -e s/\"DIFFS\"/", param$diffs, "/g ",
                                  file.path(METAGENOMICS_ROOT,UNIFIED_MOTHUR_WORKFLOW_PART2B), 
                                  " > ",
                                  UNIFIED_MOTHUR_WORKFLOW_PART2B)
  ezSystem(updateBatchPart2BCmd1)
  
  cmdMothur1 = paste("mothur",UNIFIED_MOTHUR_WORKFLOW_PART1)
  ezSystem(cmdMothur1)
  
  #check if any chimeras were found and continue the mothur workflow
  
  removeblankfilesCmd1 <- paste("find . -type f -empty -print -delete")
  ezSystem(removeblankfilesCmd1)
  cmdMothur2 = paste("[ -f merged.good.filter.unique.precluster.denovo.vsearch.accnos ] && mothur",UNIFIED_MOTHUR_WORKFLOW_PART2A, "|| mothur", UNIFIED_MOTHUR_WORKFLOW_PART2B)
  ezSystem(cmdMothur2)
  ## create and save QC and chimera summary file 
  groupFile="Mothur.groups"
  ### filter steps
  fastaFiles <- c("Mothur.fasta","Mothur.good.fasta","Mothur.good.unique.fasta",
                  "merged.align",
                  "merged.good.align",
                  "merged.good.filter.unique.fasta")
  
  filterSteps <- c("noFilter","lengthAndHomopPreAlign","deduplication","aligned",
                   "lengthAndHomopPostAlign","dedupPostEndTrimming")
  listOfFilteredSummaries <- mapply(countAndAssignSeqsFromFasta,fastaFiles,
                                    filterSteps,groupFile=groupFile,
                                    SIMPLIFY = FALSE)
  DFforQCPlot <- do.call("rbind",listOfFilteredSummaries)
  ### chimera
  chimFile="merged.good.filter.unique.precluster.denovo.vsearch.chimeras"
  DFforChimeraPlot <- chimeraSummaryTable(chimFile,groupFile)
  QCChimeraObject <- list(chimera=DFforChimeraPlot,filtStep=DFforQCPlot)  
  ### save
  QCChimeraObjectRdata <-  basename(output$getColumn("RObjectQCChimera"))
  saveRDS(QCChimeraObject,QCChimeraObjectRdata)
  
  ## create phyloseq object
  # OTU count
  newOTUsToCountFileName <- basename(output$getColumn("OTUsCountTable"))
  otuMothurFile <- "merged.good.filter.unique.precluster.pick.opti_mcc.shared"
  countOTUs <- phyloSeqOTUFromFile(otuMothurFile)
  write.table(countOTUs,newOTUsToCountFileName,
              row.names = F, col.names = T, quote = F,sep = "\t")
  # taxonomy
  newOTUsToTaxFileName <- basename(output$getColumn("OTUsToTaxonomyFile"))
  taxaMothurFile <- paste("merged.good.filter.unique.precluster.pick.opti_mcc",
                          param$cutOffCluster,"cons.taxonomy",sep = ".")
  taxaOTUs <- phyloSeqTaxaFromFile(taxaMothurFile)
  write.table(taxaOTUs,newOTUsToTaxFileName,
              row.names = F, col.names = T, quote = F,sep = "\t")
  
  ## design Matrix 
  if (param$group){
      factorCols <- grep("Factor",colnames(dataset))
      designMatrix <- data.frame(dataset[,factorCols])
      colnames(designMatrix) <- gsub(" \\[Factor\\]","",colnames(dataset)[factorCols])
      rownames(designMatrix) <- rownames(dataset)
      newdesignMatrixFileName <- basename(output$getColumn("OTUsDesignMatrix"))
      write.table(designMatrix,newdesignMatrixFileName,
                  row.names = F, col.names = T, quote = F,sep = "\t")
  }
    
    ## create phyloseqObject
  phyloseqObjectRdata <-  basename(output$getColumn("RObjectPhyloseq"))
  if (param$group){
      phyloseqObject <- phyloseq(countOTUs,taxaOTUs,sample_data(designMatrix))
  }else{
      phyloseqObject <- phyloseq(countOTUs,taxaOTUs)
  }
  saveRDS(phyloseqObject,phyloseqObjectRdata)
}

##' @template app-template
##' @templateVar method ezMethodMothur()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppMothur <-
    setRefClass("EzAppMothur",
                contains = "EzApp",
                methods = list(
                    initialize = function()
                    {
                        "Initializes the application using its specific defaults."
                        runMethod <<- ezMethodMothur
                        name <<- "EzAppMothur"
                        appDefaults <<- rbind(minLen = ezFrame(Type="integer",  DefaultValue="150",Description="Min length"),     
                                              maxLen= ezFrame(Type="integer",  DefaultValue="450",Description="Max length"),
                                              cutOffTaxonomy= ezFrame(Type="integer",  DefaultValue="80",Description="cut-off taxonomy"),
                                              cutOffCluster= ezFrame(Type="numeric",  DefaultValue="0.03",Description="cut-off cluster"),
                                              diffs= ezFrame(Type="integer",  DefaultValue="0.03",Description="differences")
                        )
                    }
                )
    )
