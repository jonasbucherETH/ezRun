# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodBowtie2 <- function(input = NA, output = NA, param = NA) {
  ref <- getBowtie2Reference(param)
  bamFile <- output$getColumn("BAM")
  sampleName <- sub(".bam", "", basename(bamFile))
  trimmedInput <- ezMethodFastpTrim(input = input, param = param)
  defOpt <- paste("-p", param$cores)
  readGroupOpt <- paste0(
    "--rg-id ", sampleName, " --rg SM:", sampleName,
    " --rg LB:RGLB_", sampleName,
    " --rg PL:illumina", " --rg PU:RGPU_", sampleName
  )
  cmd <- paste(
    "bowtie2", param$cmdOptions, defOpt, readGroupOpt,
    "-x", ref, if (param$paired) "-1", trimmedInput$getColumn("Read1"),
    if (param$paired) paste("-2", trimmedInput$getColumn("Read2")),
    "2>", paste0(sampleName, "_bowtie2.log"), "|",
    "samtools", "view -S -b -", " > bowtie.bam"
  )
  ezSystem(cmd)
  file.remove(trimmedInput$getColumn("Read1"))
  if (param$paired) {
    file.remove(trimmedInput$getColumn("Read2"))
  }

  if (!is.null(param$markDuplicates) && param$markDuplicates) {
    ezSortIndexBam("bowtie.bam", "sorted.bam",
      ram = param$ram, removeBam = TRUE,
      cores = param$cores
    )
    dupBam(inBam = "sorted.bam", outBam = basename(bamFile), operation = "mark")
    file.remove("sorted.bam")
  } else {
    ezSortIndexBam("bowtie.bam", basename(bamFile),
      ram = param$ram,
      removeBam = TRUE, cores = param$cores
    )
  }

  ## write an igv link
  if (param$writeIgvLink) {
    if ("IGV" %in% output$colNames) {
      writeIgvHtml(param, output)
    }
  }
  
  if (param$generateBigWig) {
    bam2bw(file = basename(bamFile), paired = param$paired, method = "Bioconductor", cores = param$cores)
  }
  return("Success")
}

##' @template getref-template
##' @templateVar methodName Bowtie2
##' @param param a list of parameters:
##' \itemize{
##'   \item{ezRef@@refIndex}{ a character specifying the location of the index that is used in the alignment.}
##'   \item{ezRef@@refBuildDir}{ a character specifying the directory of the reference build.}
##'   \item{ezRef@@refFastaFile}{ a character specifying the file path to the fasta file.}
##' }
getBowtie2Reference <- function(param) {
  refBase <- ifelse(param$ezRef["refIndex"] == "",
    file.path(param$ezRef["refBuildDir"], "Sequence/BOWTIE2Index/genome"),
    param$ezRef["refIndex"]
  )
  ## check the ref
  lockFile <- file.path(dirname(refBase), "lock")
  if (!file.exists(dirname(refBase))) {
    ## no lock file and no refFiles, so we build the reference
    dir.create(dirname(refBase))
    ezWrite(Sys.info(), con = lockFile)
    wd <- getwd()
    setwd(dirname(refBase))

    fastaFile <- param$ezRef["refFastaFile"]
    ezSystem(paste("ln -s", fastaFile, "."))
    cmd <- paste("bowtie2-build", "--threads", as.numeric(param$cores), "-f", basename(fastaFile), basename(refBase))
    ezSystem(cmd)
    # ezWriteElapsed(job, "done")
    setwd(wd)
    file.remove(lockFile)
  }
  stopifnot(file.exists(dirname(refBase)))
  i <- 0
  while (file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT) {
    ### somebody else builds and we wait
    Sys.sleep(60)
    i <- i + 1
  }
  if (file.exists(lockFile)) {
    stop(paste("reference building still in progress after", INDEX_BUILD_TIMEOUT, "min"))
  }
  ## there is no lock file
  refFiles <- list.files(dirname(refBase), basename(refBase))
  if (length(refFiles) < 3) {
    ## we assume the index is built and complete
    stop(paste("index not available: ", refBase))
  }
  return(refBase)
}

##' @template app-template
##' @templateVar method ezMethodBowtie2(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
##' @seealso \code{\link{getBowtie2Reference}}
##' @seealso \code{\link{ezMethodFastpTrim}}
EzAppBowtie2 <-
  setRefClass("EzAppBowtie2",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodBowtie2
        name <<- "EzAppBowtie2"
        appDefaults <<- rbind(
          writeIgvSessionLink = ezFrame(Type = "logical", DefaultValue = "TRUE", Description = "should an IGV link be generated"),
          markDuplicates = ezFrame(Type = "logical", DefaultValue = "TRUE", Description = "should duplicates be marked"),
          generateBigWig = ezFrame(Type = "logical", DefaultValue = "FALSE", Description = "should a bigwig file be generated")
        )
      }
    )
  )

ezMethodBowtie <- function(input = NA, output = NA, param = NA) {
  ref <- getBowtieReference(param)
  bamFile <- output$getColumn("BAM")
  trimmedInput <- ezMethodFastpTrim(input = input, param = param)
  defOpt <- paste("--chunkmbs 256", "--sam", "-p", param$cores)
  cmd <- paste(
    "bowtie", param$cmdOptions, defOpt,
    ref, if (param$paired) "-1", trimmedInput$getColumn("Read1"),
    if (param$paired) paste("-2", trimmedInput$getColumn("Read2")),
    "2> bowtie.log", "|", "samtools", "view -S -b -", " > bowtie.bam"
  )
  ezSystem(cmd)
  ezSortIndexBam("bowtie.bam", basename(bamFile),
    ram = param$ram, removeBam = TRUE,
    cores = param$cores
  )

  ## write an igv link
  if (param$writeIgvLink) {
    if ("IGV" %in% output$colNames) {
      writeIgvHtml(param, output)
    }
  }
  return("Success")
}

##' @template getref-template
##' @templateVar methodName Bowtie
##' @inheritParams getBowtie2Reference
getBowtieReference <- function(param) {
  refBase <- ifelse(param$ezRef["refIndex"] == "",
    file.path(param$ezRef["refBuildDir"], "Sequence/BOWTIEIndex/genome"),
    param$ezRef["refIndex"]
  )
  ## check the ref
  lockFile <- file.path(dirname(refBase), "lock")
  if (!file.exists(dirname(refBase))) {
    ## no lock file and no refFiles, so we build the reference
    dir.create(dirname(refBase))
    ezWrite(Sys.info(), con = lockFile)
    wd <- getwd()
    setwd(dirname(refBase))

    fastaFile <- param$ezRef["refFastaFile"]
    # job = ezJobStart("bowtie index")
    ezSystem(paste("ln -s", fastaFile, "."))
    buildOpts <- ""
    if (any(grepl("--threads", system("bowtie-build --help", intern = T)))) {
      buildOpts <- paste("--threads", param$cores)
    }
    cmd <- paste("bowtie-build", buildOpts, "-f", basename(fastaFile), basename(refBase))
    ezSystem(cmd)
    # ezWriteElapsed(job, "done")
    setwd(wd)
    file.remove(lockFile)
  }
  stopifnot(file.exists(dirname(refBase)))
  i <- 0
  while (file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT) {
    ### somebody else builds and we wait
    Sys.sleep(60)
    i <- i + 1
  }
  if (file.exists(lockFile)) {
    stop(paste("reference building still in progress after", INDEX_BUILD_TIMEOUT, "min"))
  }
  ## there is no lock file
  refFiles <- list.files(dirname(refBase), basename(refBase))
  if (length(refFiles) < 3) {
    ## we assume the index is built and complete
    stop(paste("index not available: ", refBase))
  }
  return(refBase)
}

##' @template app-template
##' @templateVar method ezMethodBowtie(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
##' @seealso \code{\link{getBowtieReference}}
##' @seealso \code{\link{ezMethodFastpTrim}}
EzAppBowtie <-
  setRefClass("EzAppBowtie",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodBowtie
        name <<- "EzAppBowtie"
        appDefaults <<- rbind(writeIgvSessionLink = ezFrame(Type = "logical", DefaultValue = "TRUE", Description = "should an IGV link be generated"))
      }
    )
  )

ezMethodSTAR <- function(input = NA, output = NA, param = NA) {
  refDir <- getSTARReference(param)
  bamFile <- output$getColumn("BAM")
  trimmedInput <- ezMethodFastpTrim(input = input, param = param)

  if (!str_detect(param$cmdOptions, "outSAMattributes")) {
    param$cmdOptions <- str_c(param$cmdOptions, "--outSAMattributes All",
      sep = " "
    )
  }
  genomeFn <- param$ezRef@refFastaFile
  
  if (ezIsSpecified(param$controlSeqs)) {
    ## control sequences
    controlSeqsLocalFn <- tempfile(
      pattern = "controlSeqs", tmpdir = getwd(), fileext = ".fa"
    )
    writeXStringSet(getControlSeqs(param$controlSeqs), filepath = controlSeqsLocalFn)
    defer(file.remove(controlSeqsLocalFn))
    
    genomeLocalFn <- tempfile(
      pattern = "genome", tmpdir = getwd(), fileext = ".fa"
    )
    file.copy(from = genomeFn, to = genomeLocalFn)
    writeXStringSet(getControlSeqs(param$controlSeqs),
                    filepath = genomeLocalFn,
                    append = TRUE
    )
    dictFile <- sub(".fa$", ".dict", genomeLocalFn)
    cmd <- str_c(
      prepareJavaTools("picard"), "CreateSequenceDictionary",
      str_c("R=", genomeLocalFn), str_c("O=", dictFile),
      sep = " "
    )
    ezSystem(cmd)
    
    genomeFn <- genomeLocalFn
    defer(file.remove(c(genomeLocalFn, dictFile)))
  }
  

  cmd <- str_c(
    "STAR", "--genomeDir", refDir,
    "--readFilesIn", trimmedInput$getColumn("Read1"),
    if (param$paired) trimmedInput$getColumn("Read2"),
    "--twopassMode", if_else(param$twopassMode, "Basic", "None"),
    "--runThreadN", param$cores, param$cmdOptions,
    ifelse(ezIsSpecified(param$controlSeqs),
           str_c("--genomeFastaFiles", controlSeqsLocalFn, sep = " "), ""),
    "--outStd BAM_Unsorted --outSAMtype BAM Unsorted",
    "--outSAMattrRGline", str_c("ID:", trimmedInput$getNames()), str_c("SM:", trimmedInput$getNames()),
    if_else(str_detect(trimmedInput$getColumn("Read1"), "\\.gz$"), "--readFilesCommand zcat", ""),
    ">  Aligned.out.bam",
    sep = " "
  )
  ezSystem(cmd)

  nSortThreads <- min(param$cores, 8)
  ## if the index is loaded in shared memory we have to use only 10% of the scheduled RAM
  if (str_detect(param$cmdOptions, "--genomeLoad Load")) {
    sortRam <- param$ram / 10
  } else {
    sortRam <- param$ram
  }

  file.rename(
    from = "Log.final.out",
    to = basename(output$getColumn("STARLog"))
  )

  if (!is.null(param$markDuplicates) && param$markDuplicates) {
    ezSortIndexBam("Aligned.out.bam", "sorted.bam",
      ram = sortRam, removeBam = TRUE,
      cores = nSortThreads
    )
    dupBam(inBam = "sorted.bam", outBam = basename(bamFile), operation = "mark")
    file.remove("sorted.bam")
  } else {
    ezSortIndexBam("Aligned.out.bam", basename(bamFile),
      ram = sortRam,
      removeBam = TRUE, cores = nSortThreads
    )
  }

  if (param$getJunctions) {
    file.rename(
      from = "SJ.out.tab",
      to = basename(output$getColumn("Junctions"))
    )
    file.rename(
      from = "Chimeric.out.junction",
      to = basename(output$getColumn("Chimerics"))
    )
  }

  ## check the strandedness
  ezSystem(str_c(
    "infer_experiment.py", "-r", getReferenceFeaturesBed(param),
    "-i", basename(bamFile), "-s 1000000",
    sep = " "
  ),
  stopOnFailure = FALSE
  )

  ## write an igv link
  if (param$writeIgvLink && "IGV" %in% output$colNames) {
    writeIgvHtml(param, output)
  }
  return("Success")
}

##' @template getref-template
##' @templateVar methodName STAR
##' @param param a list of parameters:
##' \itemize{
##'   \item{ezRef@@refIndex}{ a character specifying the location of the index that is used in the alignment.}
##'   \item{ezRef@@refFeatureFile}{ a character specifying the file path to the annotation feature file (.gtf).}
##'   \item{ram}{ an integer specifying how many gigabytes of RAM to use.}
##'   \item{ezRef@@refFastaFile}{ a character specifying the file path to the fasta file.}
##' }
getSTARReference <- function(param) {
  if (ezIsSpecified(param$ezRef["refIndex"])) {
    refDir <- param$ezRef["refIndex"]
  } else {
    if (!ezIsSpecified(param$ezRef["refFeatureFile"])) {
      stop("refFeatureFile not defined")
    }
    refDir <- sub(".gtf$", "_STARIndex", param$ezRef["refFeatureFile"])
    if (ezIsSpecified(param$spikeInSet)) {
      refDir <- paste(refDir, param$spikeInSet, sep = "_")
    }
  }
  ## random sleep to avoid parallel ref building
  Sys.sleep(runif(1, max = 20))

  lockFile <- paste0(refDir, ".lock")
  i <- 0
  while (file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT) {
    ### somebody else builds and we wait
    Sys.sleep(60)
    i <- i + 1
  }
  if (i >= INDEX_BUILD_TIMEOUT) {
    stop(paste("reference building still in progress after", INDEX_BUILD_TIMEOUT, "min"))
  }
  ## there is no lock file
  if (file.exists(file.path(refDir, "SAindex"))) {
    ## we assume the index is built and complete
    return(refDir)
  }

  ## no lock file and no refFiles, so we build the reference
  ezWrite(Sys.info(), con = lockFile)
  dir.create(refDir)

  fai <- ezRead.table(paste0(param$ezRef["refFastaFile"], ".fai"), header = FALSE)
  colnames(fai) <- c("LENGTH", "OFFSET", "LINEBASES", "LINEWDITH")
  if (nrow(fai) > 50) {
    binOpt <- "--genomeChrBinNbits 16"
  } else {
    binOpt <- ""
  }

  genomeLength <- sum(fai$LENGTH)
  readLength <- 150 ## assumption
  indexNBasesOpt <- paste("--genomeSAindexNbases", min(13, floor(log2(genomeLength) / 2 - 1)))
  if (binOpt == "") {
    genomeChrBinNbits <- paste("--genomeChrBinNbits", floor(min(
      18,
      log2(max(genomeLength / nrow(fai), readLength))
    )))
  } else {
    genomeChrBinNbits <- ""
  }

  job <- ezJobStart("STAR genome build")
  if (ezIsSpecified(param$spikeInSet)) {
    spikeInFasta <- paste0(SPIKEINS_ROOT, "/", param$spikeInSet, "/", param$spikeInSet, ".fa")
    spikeInGtf <- paste0(SPIKEINS_ROOT, "/", param$spikeInSet, "/", param$spikeInSet, ".gtf")
    gtfFile <- file.path(refDir, "genesAndSpikes.gtf")
    ezSystem(paste("cp", param$ezRef["refFeatureFile"], gtfFile))
    ezSystem(paste("cat", spikeInGtf, ">>", gtfFile))
    genomeFastaFiles <- paste(param$ezRef["refFastaFile"], spikeInFasta)
  } else {
    gtfFile <- param$ezRef["refFeatureFile"]
    genomeFastaFiles <- param$ezRef["refFastaFile"]
  }
  cmd <- paste(
    "STAR", "--runMode genomeGenerate --genomeDir", refDir, binOpt, indexNBasesOpt, genomeChrBinNbits,
    "--limitGenomeGenerateRAM", format(param$ram * 1e9, scientific = FALSE),
    "--genomeFastaFiles", genomeFastaFiles,
    "--sjdbGTFfile", gtfFile, "--sjdbOverhang 150", "--runThreadN", param$cores, "--genomeSAsparseD 2"
  )
  ezSystem(cmd)
  file.remove(lockFile)
  ezWriteElapsed(job, "done")
  file.remove("Log.out")
  return(refDir)
}

##' @template app-template
##' @templateVar method ezMethodSTAR(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
##' @seealso \code{\link{getSTARReference}}
##' @seealso \code{\link{ezMethodFastpTrim}}
EzAppSTAR <-
  setRefClass("EzAppSTAR",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodSTAR
        name <<- "EzAppSTAR"
        appDefaults <<- rbind(
          getJunctions = ezFrame(Type = "logical", DefaultValue = "FALSE", Description = "should junctions be returned"),
          writeIgvLink = ezFrame(Type = "logical", DefaultValue = "TRUE", Description = "should an IGV link be generated"),
          markDuplicates = ezFrame(Type = "logical", DefaultValue = "TRUE", Description = "should duplicates be marked with picard"),
          twopassMode = ezFrame(Type = "logical", DefaultValue = "TRUE", Description = "1-pass mapping or basic 2-pass mapping")
        )
      }
    )
  )

ezMethodBWA <- function(input = NA, output = NA, param = NA) {
  refIdx <- getBWAReference(param)
  bamFile <- output$getColumn("BAM")
  trimmedInput <- ezMethodFastpTrim(input = input, param = param)
  if (param$algorithm == "aln") {
    cmd <- paste(
      "bwa", param$algorithm, param$cmdOptions, "-t", param$cores,
      refIdx, trimmedInput$getColumn("Read1"), ">", "read1.sai", "2> bwa.log"
    )
    ezSystem(cmd)
    if (param$paired) {
      cmd <- paste(
        "bwa", param$algorithm, param$cmdOptions, "-t", param$cores,
        refIdx, trimmedInput$getColumn("Read2"), ">", "read2.sai", "2> bwa.log"
      )
      ezSystem(cmd)
      cmd <- paste(
        "bwa", "sampe", refIdx, "read1.sai", "read2.sai",
        trimmedInput$getColumn("Read1"), trimmedInput$getColumn("Read2"),
        "2> bwa.log", "|",
        "samtools", "view -S -b -", " > aligned.bam"
      )
      ezSystem(cmd)
    } else {
      cmd <- paste(
        "bwa", "samse", refIdx, "read1.sai",
        trimmedInput$getColumn("Read1"), "2> bwa.log", "|",
        "samtools", "view -S -b -", " > aligned.bam"
      )
      ezSystem(cmd)
    }
  } else {
    if (param$algorithm == "bwasw" && param$paired) {
      stop("paired is not supported for algorithm bwasw")
    }
    cmd <- paste(
      "bwa", param$algorithm, param$cmdOptions, "-t", param$cores,
      refIdx, trimmedInput$getColumn("Read1"),
      if (param$paired) trimmedInput$getColumn("Read2"),
      "2> bwa.log", "|", "samtools", "view -S -b -", " > aligned.bam"
    )
    ezSystem(cmd)
  }
  file.remove(trimmedInput$getColumn("Read1"))
  if (param$paired) {
    file.remove(trimmedInput$getColumn("Read2"))
  }

  if (!is.null(param$markDuplicates) && param$markDuplicates) {
    ezSortIndexBam("aligned.bam", "sorted.bam",
      ram = param$ram, removeBam = TRUE,
      cores = param$cores
    )
    dupBam(inBam = "sorted.bam", outBam = basename(bamFile), operation = "mark")
    file.remove("sorted.bam")
  } else {
    ezSortIndexBam("aligned.bam", basename(bamFile),
      ram = param$ram,
      removeBam = TRUE, cores = param$cores
    )
  }


  ## write an igv link
  if (param$writeIgvLink) {
    if ("IGV" %in% output$colNames) {
      writeIgvHtml(param, output)
    }
  }
  return("Success")
}

ezMethodBWATrimmomatic <- function(input = NA, output = NA, param = NA) { # Perform BWA using fastp for read pre-processing

  refIdx <- getBWAReference(param)
  bamFile <- output$getColumn("BAM")
  trimmedInput <- ezMethodTrim(input = input, param = param)
  if (param$algorithm == "aln") {
    cmd <- paste(
      "bwa", param$algorithm, param$cmdOptions, "-t", param$cores,
      refIdx, trimmedInput$getColumn("Read1"), ">", "read1.sai", "2> bwa.log"
    )
    ezSystem(cmd)
    if (param$paired) {
      cmd <- paste(
        "bwa", param$algorithm, param$cmdOptions, "-t", param$cores,
        refIdx, trimmedInput$getColumn("Read2"), ">", "read2.sai", "2> bwa.log"
      )
      ezSystem(cmd)
      cmd <- paste(
        "bwa", "sampe", refIdx, "read1.sai", "read2.sai",
        trimmedInput$getColumn("Read1"), trimmedInput$getColumn("Read2"),
        "2> bwa.log", "|",
        "samtools", "view -S -b -", " > aligned.bam"
      )
      ezSystem(cmd)
    } else {
      cmd <- paste(
        "bwa", "samse", refIdx, "read1.sai",
        trimmedInput$getColumn("Read1"), "2> bwa.log", "|",
        "samtools", "view -S -b -", " > aligned.bam"
      )
      ezSystem(cmd)
    }
  } else {
    if (param$algorithm == "bwasw" && param$paired) {
      stop("paired is not supported for algorithm bwasw")
    }
    cmd <- paste(
      "bwa", param$algorithm, param$cmdOptions, "-t", param$cores,
      refIdx, trimmedInput$getColumn("Read1"),
      if (param$paired) trimmedInput$getColumn("Read2"),
      "2> bwa.log", "|", "samtools", "view -S -b -", " > aligned.bam"
    )
    ezSystem(cmd)
  }
  file.remove(trimmedInput$getColumn("Read1"))
  if (param$paired) {
    file.remove(trimmedInput$getColumn("Read2"))
  }

  ezSortIndexBam("aligned.bam", basename(bamFile),
    ram = param$ram, removeBam = TRUE,
    cores = param$cores
  )

  ## write an igv link
  if (param$writeIgvLink) {
    if ("IGV" %in% output$colNames) {
      writeIgvHtml(param, output)
    }
  }
  return("Success")
}

##' @template getref-template
##' @templateVar methodName BWA
##' @inheritParams getBowtie2Reference
getBWAReference <- function(param) {
  refPath <- ifelse(param$ezRef["refIndex"] == "",
    file.path(param$ezRef["refBuildDir"], "Sequence/BWAIndex/genome.fa"),
    param$ezRef["refIndex"]
  )
  ## check the ref
  lockFile <- file.path(dirname(refPath), "lock")
  i <- 0
  while (file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT) {
    ### somebody else builds and we wait
    Sys.sleep(60)
    i <- i + 1
  }
  if (file.exists(lockFile)) {
    stop(paste("reference building still in progress after", INDEX_BUILD_TIMEOUT, "min"))
  }
  if (!file.exists(dirname(refPath))) {
    ## no lock file and no refFiles, so we build the reference
    dir.create(dirname(refPath))
    ezWrite(Sys.info(), con = lockFile)
    wd <- getwd()
    setwd(dirname(refPath))

    fastaFile <- param$ezRef["refFastaFile"]
    file.symlink(from = fastaFile, to = ".")
    cmd <- paste("bwa", "index", "-a", "bwtsw", basename(fastaFile))
    ezSystem(cmd)
    setwd(wd)
    file.remove(lockFile)
  }
  stopifnot(file.exists(dirname(refPath)))
  if (!file.exists(paste0(refPath, ".sa"))) {
    stop(paste("sa index not found for:", refPath))
  }
  return(refPath)
}

##' @template app-template
##' @templateVar method ezMethodBWA(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
##' @seealso \code{\link{getBWAReference}}
##' @seealso \code{\link{ezMethodFastpTrim}}
EzAppBWA <-
  setRefClass("EzAppBWA",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodBWA
        name <<- "EzAppBWA"
        appDefaults <<- rbind(
          algorithm = ezFrame(Type = "character", DefaultValue = "mem", Description = "bwa's alignment algorithm. One of aln, bwasw, mem."),
          writeIgvSessionLink = ezFrame(Type = "logical", DefaultValue = "TRUE", Description = "should an IGV link be generated"),
          markDuplicates = ezFrame(Type = "logical", DefaultValue = "TRUE", Description = "should duplicates be marked")
        )
      }
    )
  )

EzAppBWATrimmomatic <-
  setRefClass("EzAppBWATrimmomatic",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodBWATrimmomatic
        name <<- "EzAppBWATrimmomatic"
        appDefaults <<- rbind(
          algorithm = ezFrame(Type = "character", DefaultValue = "mem", Description = "bwa's alignment algorithm. One of aln, bwasw, mem."),
          writeIgvSessionLink = ezFrame(Type = "logical", DefaultValue = "TRUE", Description = "should an IGV link be generated")
        )
      }
    )
  )






ezMethodBismark <- function(input = NA, output = NA, param = NA) {
  
  cat("000")

  ref <- getBismarkReference(param)
  bamFile <- output$getColumn("BAM")
  trimmedInput <- ezMethodFastpTrim(input = input, param = param)
  defOpt <- paste("-p", max(2, param$cores / 2))
  if (param$paired) {
    cmd <- paste(
      "bismark", param$cmdOptions,
      "--path_to_bowtie", paste0("$Bowtie2", "/bin"), defOpt, ref,
      "-1", trimmedInput$getColumn("Read1"),
      if (param$paired) paste("-2", trimmedInput$getColumn("Read2")),
      "2> bismark.log"
    )
    if (param$dirty_harry) {
      cmd <- paste(cmd, "--unmapped")
    }
  } else {
    cmd <- paste(
      "bismark", param$cmdOptions,
      "--path_to_bowtie", paste0("$Bowtie2", "/bin"), defOpt, ref,
      trimmedInput$getColumn("Read1"),
      "2> bismark.log"
    )
  }
  
  ezSystem(cmd)
  bamFileNameBismark <- list.files(".", pattern = "bam$")
  reportFileNameBismark <- list.files(".", pattern = "report.txt$")
  reportFile <- paste0(names(bamFile), ".report.txt")
  ezSystem(paste("mv ", reportFileNameBismark, reportFile))
  if (param$deduplicate) {
    cmd <- paste("deduplicate_bismark", ifelse(param$paired, "-p", "-s"), bamFileNameBismark)
    ezSystem(cmd)
    bamFileNameBismark <- list.files(".", pattern = "deduplicated.bam$")
    deduplicationReportFile <- list.files(".", pattern = "deduplication_report.txt$")
    ezSystem(paste("cat ", deduplicationReportFile, ">>", reportFile))
  }
  
  # cmd <- paste("bismark_methylation_extractor", bamFileNameBismark, ifelse(param$paired, "-p", "-s"), "--bedGraph", "--cytosine_report", "--genome_folder", ref)
  # ezSystem(cmd)
  # 
  # cmd <- paste("samtools", "view -S -b ", bamFileNameBismark, " > bismark.bam")
  # ezSystem(cmd)
  # ezSortIndexBam("bismark.bam", basename(bamFile),
  #                ram = param$ram, removeBam = TRUE,
  #                cores = param$cores
  # )
  
  
  # cmd <- paste("bismark2bedGraph", "-o", bedGraphCpG, "--ample_memory", filesCpG)
  # ezSystem(cmd)
  # cmd <- paste("bismark2bedGraph", "-o", bedGraphCHG, "--CX", "--ample_memory", filesCHG)
  # ezSystem(cmd)
  # cmd <- paste("bismark2bedGraph", "-o", bedGraphCHH, "--CX", "--ample_memory", filesCHH)
  # ezSystem(cmd)
  
  # cmd <- paste("bismark2bedGraph", "-o", bedGraphCpG, "CpG*")
  # ezSystem(cmd)
  # cmd <- paste("bismark2bedGraph", "-o", bedGraphCHG, "--CX", "CHG*")
  # ezSystem(cmd)
  # cmd <- paste("bismark2bedGraph", "-o", bedGraphCHH, "--CX", "CHH*")
  # ezSystem(cmd)
  
  # bedGraphCpG <- paste0("CpG_", names(bamFile), ".bedGraph")
  # bedGraphCHG <- paste0("CHG_", names(bamFile), ".bedGraph")
  # bedGraphCHH <- paste0("CHH_", names(bamFile), ".bedGraph")
  
  # bedGraphCpG <- paste0("CpG_", names(bamFile))
  # bedGraphCHG <- paste0("CHG_", names(bamFile))
  # bedGraphCHH <- paste0("CHH_", names(bamFile))
  
  
  
  cmd <- paste("bismark_methylation_extractor", ifelse(param$paired, "-p", "-s"), "--comprehensive", bamFileNameBismark)
  ezSystem(cmd)
  
  if (param$dirty_harry) {
    unmappedReads1 <- list.files(".", pattern = "*unmapped_reads_1.fq.gz")
    unmappedReads2 <- list.files(".", pattern = "*unmapped_reads_2.fq.gz")
    if (param$directional) {
      cmd <- paste(
        "bismark", "--se", unmappedReads1, 
        "--path_to_bowtie", paste0("$Bowtie2", "/bin"), defOpt, ref,
        "--local" 
      )
      ezSystem(cmd)
      cmd <- paste(
        "bismark", "--se", unmappedReads2, "--pbat",
        "--path_to_bowtie", paste0("$Bowtie2", "/bin"), defOpt, ref,
        "--local"
      )
      ezSystem(cmd)
    } else {
      cmd <- paste(
        "bismark", "--se", unmappedReads1, "--pbat",
        "--path_to_bowtie", paste0("$Bowtie2", "/bin"), defOpt, ref,
        "--local"
      )
      ezSystem(cmd)
      cmd <- paste(
        "bismark", "--se", unmappedReads2, 
        "--path_to_bowtie", paste0("$Bowtie2", "/bin"), defOpt, ref,
        "--local"
      )
      ezSystem(cmd)
    }
    
    # unmappedBam1 <- list.files(".", pattern = paste0(names(bamFile), "*unmapped_reads_1_bismark_bt2.bam"))
    # unmappedBam2 <- list.files(".", pattern = paste0(names(bamFile), "*unmapped_reads_2_bismark_bt2.bam"))
    cat("#################################")
    cat(names(bamFile))
    unmappedBam1 <- list.files(".", pattern = "*unmapped_reads_1_bismark_bt2.bam")
    unmappedBam2 <- list.files(".", pattern = "*unmapped_reads_2_bismark_bt2.bam")
    cmd <- paste("bismark_methylation_extractor", "-s", "--comprehensive", unmappedBam1)
    ezSystem(cmd)
    cmd <- paste("bismark_methylation_extractor", "-s", "--comprehensive", unmappedBam2)
    ezSystem(cmd)
  }
  
  cmd <- paste("samtools", "view -S -b ", bamFileNameBismark, " > bismark.bam")
  ezSystem(cmd)
  ezSortIndexBam("bismark.bam", basename(bamFile),
                 ram = param$ram, removeBam = TRUE,
                 cores = param$cores
  )
  
  if (param$generateBigWig) {
    destination = sub("\\.bam$", "_Cov.bw", basename(bamFile), ignore.case = TRUE)
    bam2bw(file = basename(bamFile), destination = destination, paired = param$paired, method = "Bioconductor", cores = param$cores)
  }
  
  mBiasImages <- list.files(".", pattern = "png$")
  if ((length(mBiasImages)) > 0) {
    for (i in 1:length(mBiasImages)) {
      ezSystem(paste("mv ", mBiasImages[i], paste0(names(bamFile), ".M-bias_R", i, ".png")))
    }
  } else {
    ezSystem(paste("touch", paste0(names(bamFile), ".M-bias_R1.png")))
    ezSystem(paste("touch", paste0(names(bamFile), ".M-bias_R2.png")))
  }
  
  CpGFile <- list.files(".", pattern = "^CpG.*txt$")
  cmd <- paste("bismark2bedGraph --scaffolds", "CpG*", "-o", names(bamFile)) # adds .gz automatically if not already there
  ezSystem(cmd)
  ezSystem(paste("mv ", CpGFile[1], paste0(names(bamFile), ".CpG_context.txt")))
  
  # covFileCpG <- list.files(".", pattern = paste0(names(bamFile), ".gz.bismark.cov.gz$")) 
  # cmd <- paste("coverage2cytosine", "--gzip", "--genome_folder", ref, "-o", paste0(names(bamFile), ".CpG_report.txt"), covFileCpG)
  # ezSystem(cmd)
  
  if (param$allCytosineContexts) {
    CHGFile <- list.files(".", pattern = "^CHG.*txt$")
    CHHFile <- list.files(".", pattern = "^CHH.*txt$")
    cmd <- paste("bismark2bedGraph --scaffolds", "CHG*", "--CX", "-o", paste0(names(bamFile), ".CHG"))
    ezSystem(cmd)
    cmd <- paste("bismark2bedGraph --scaffolds", "CHH*", "--CX", "-o", paste0(names(bamFile), ".CHH"))
    ezSystem(cmd)
    ezSystem(paste("mv ", CHGFile[1], paste0(names(bamFile), ".CHG_context.txt")))
    ezSystem(paste("mv ", CHHFile[1], paste0(names(bamFile), ".CHH_context.txt")))
    
    # covFileCHG <- list.files(".", pattern = "*CHG.gz.bismark.cov.gz$")
    # cmd <- paste("coverage2cytosine", "--CX", "--gzip", "--genome_folder", ref, "-o", paste0(names(bamFile), ".CHG_report.txt"), covFileCHG)
    # ezSystem(cmd)
    #                                                                                
    # covFileCHH <- list.files(".", pattern = "*CHH.gz.bismark.cov.gz$") 
    # cmd <- paste("coverage2cytosine", "--CX", "--gzip", "--genome_folder", ref, "-o", paste0(names(bamFile), ".CHH_report.txt"), covFileCHH)
    # ezSystem(cmd)
  }
  
  # covFiles <- list.files(".", pattern = paste0(names(bamFile), "*.gz.bismark.cov.gz$"))
  # covFiles <- list.files(".", pattern = "*.gz.bismark.cov.gz$")
  
  # cmd <- paste("coverage2cytosine", "--CX", "--gzip", "--genome_folder", ref, "-o", names(bamFile), "*.gz.bismark.cov.gz")
  cmd <- paste("coverage2cytosine", "--CX", "--gzip", "--genome_folder", ref, "-o", names(bamFile), "*.gz.bismark.cov.gz")
  ezSystem(cmd)
  
  if (param$allCytosineContexts) {
    # cmd <- paste("awk -F'\t' '{print > ($6", paste0("\"_report.txt\")}'"), paste0(names(bamFile), "*CX_report*"))
    # cmd <- paste("awk -F '\\t' '{print >", paste0("\"", names(bamFile), ".\"", "($6 ", "\"_report.txt.gz\")}'"), paste0("<(zcat ", names(bamFile), "*CX_report.txt.gz)"))
    # cmd <- paste("awk -F '\\t' '{print >", paste0("\"", names(bamFile), ".\"", "($6 \"_report.txt\")}'"), paste0("<(zcat ", names(bamFile), "*CX_report.txt.gz)"))
    # cmd <- paste("zcat", paste0(names(bamFile), "*CX_report.txt.gz"), "|", "awk -F '\t' '{print >", paste0("\"", names(bamFile), ".\"", "($6 ", "\"_report.txt\")}'"))
    cmd <- paste("zcat", paste0(names(bamFile), "*CX_report.txt.gz"), "> temp_decompressed_file")
    ezSystem(cmd)
    cmd <- paste("awk -F '\\t' '{print >", paste0("\"", names(bamFile), ".\"", "($6 \"_report.txt\")}'"), "temp_decompressed_file")
    ezSystem(cmd)
    ezSystem("rm temp_decompressed_file")
    cmd <- paste("gzip", paste0(names(bamFile), "*_report.txt"))
    ezSystem(cmd)
  }
  
  # bamFile <- "a"
  # names(bamFile) <- "SRR4105496"
  # awk -F'\t' '{print > \"SRR4105496.\"($6 \"_report.txt.gz\")}' <(zcat SRR4105496*CX_report)
  # awk -F'\t' '{print > "SRR4105496."($1 "_report.txt.gz")}' <(zcat /srv/gstore/projects/p1535/Bismark_mm_full3_2023-07-13--11-55-16/SRR4105496.CHG.gz)
  
  # awk -F'\t' '{print > ($6 "_report.A04_BB.txt")}' A04_BB.CX_report.txt
  # awk -F'\t' '{print > ("A04_BB." $6 "_report.txt")}' A04_BB.CX_report.txt
  # awk -F '\t' '{print > "A04_BB."($6 "_report.txt.gz")}' <(zcat A04_BB*CX_report.txt.gz)
  # awk -F '\t' '{print > "A04_BB."($6 "_report.txt")}' <(zcat A04_BB*CX_report.txt.gz) | gzip
  # zcat A04_BB*CX_report.txt.gz | awk -F '\t' '{print $6 | gzip > "A04_BB." $6 "_report.txt.gz")}'
  # zcat A04_BB*CX_report.txt.gz | awk -F '\t' '{print > "A04_BB."($6 "_report.txt.gz")}'
  # echo 1 2 | awk '{print $2 | "gzip > "$1".gz"}'
  # awk -F 't' '{print > "A04_BB."($6 "_report.txt.gz")}' <(zcat A04_BB*CX_report.txt.gz)

  # splittingReportFile <- list.files(".", pattern = "splitting_report.txt$")
  # ezSystem(paste("cat ", splittingReportFile[1], ">>", reportFile))
  
  ## write an igv link
  # if (param$writeIgvSessionLink){
  #  writeIgvSession(genome = getIgvGenome(param), refBuild=param$ezRef["refBuild"], file=basename(output$getColumn("IGV Session")),
  #                  bamUrls = paste(PROJECT_BASE_URL, bamFile, sep="/") )
  #  writeIgvJnlp(jnlpFile=basename(output$getColumn("IGV Starter")), projectId = sub("\\/.*", "", bamFile),
  #               sessionUrl = paste(PROJECT_BASE_URL, output$getColumn("IGV Session"), sep="/"))
  # }
  return("Success")
}


##' @template getref-template
##' @templateVar methodName Bismark
##' @param param a list of parameters:
##' \itemize{
##'   \item{ezRef@@refIndex}{ a character specifying the location of the index that is used in the alignment.}
##'   \item{ezRef@@refBuildDir}{ a character specifying the directory of the reference build.}
##'   \item{ezRef@@refFastaFile}{ a character specifying the file path to the fasta file.}
##' }
getBismarkReference <- function(param) {
  refBase <- ifelse(param$ezRef["refIndex"] == "",
    file.path(param$ezRef["refBuildDir"], "Sequence/WholeGenomeFasta/Bisulfite_Genome/CT_conversion"),
    param$ezRef["refIndex"]
  )
  ## check the ref
  lockFile <- file.path(dirname(refBase), "lock")
  if (!file.exists(dirname(refBase))) {
    ## no lock file and no refFiles, so we build the reference
    dir.create(dirname(refBase))
    ezWrite(Sys.info(), con = lockFile)
    wd <- getwd()
    cmd <- paste(
      "bismark_genome_preparation", dirname(param$ezRef["refFastaFile"]),
      "2> bismarkGenomePrep.log"
    )
    ezSystem(cmd)

    # ezWriteElapsed(job, "done")
    setwd(wd)
    file.remove(lockFile)
  }
  stopifnot(file.exists(dirname(refBase)))
  i <- 0
  while (file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT) {
    ### somebody else builds and we wait
    Sys.sleep(60)
    i <- i + 1
  }
  if (file.exists(lockFile)) {
    stop(paste("reference building still in progress after", INDEX_BUILD_TIMEOUT, "min"))
  }
  ## there is no lock file
  refFiles <- list.files(dirname(refBase), basename(refBase))
  if (length(refFiles) < 1) {
    ## we assume the index is built and complete
    stop(paste("index not available: ", refBase))
  }
  return(dirname(param$ezRef@refFastaFile))
}



##' @template app-template
##' @templateVar method ezMethodBismark(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
EzAppBismark <-
  setRefClass("EzAppBismark",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodBismark
        name <<- "EzAppBismark"
        appDefaults <<- rbind(
            generateBigWig = ezFrame(Type = "logical", DefaultValue = "FALSE", Description = "should a bigwig coverage file be generated")
        )
      }
    )
  )
