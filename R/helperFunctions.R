.printMultipleGraphsPerPage <- function(plots.l, nCol = 1, nRow = 1, pdfFile = NULL, height = NULL, width = NULL, verbose = FALSE) {
    
    checkmate::assertList(plots.l, min.len = 1)
    checkmate::assertInt(nCol, lower = 1)
    checkmate::assertInt(nRow, lower = 1)
    checkmate::assert(checkmate::checkNull(pdfFile), checkmate::checkCharacter(pdfFile))
    if (!checkmate::testNull(pdfFile)) checkmate::assertDirectory(dirname(pdfFile), access = "r")
    checkmate::assert(checkmate::checkNull(height), checkmate::checkNumber(height, lower = 1))
    checkmate::assert(checkmate::checkNull(width), checkmate::checkNumber(width, lower = 1))
    
    if (checkmate::testNull(height)) {
        height = 7
    }
    
    if (checkmate::testNull(width)) {
        width = 7
    }
    
    
    for (i in seq_len(length(plots.l))) {
        #checkmate::assertClass(plots.l[[i]], classes = c("ggplot", "gg", "gtable", "grob"))
    }
    
    
    nPlotsPerPage =  nCol * nRow
    
    plotsNew.l = list()
    
    if (!checkmate::testNull(pdfFile))  {
        .clearOpenDevices()
        grDevices::pdf(pdfFile, height = height, width = width)
    }
    
    index = 0
    
    pb <- progress::progress_bar$new(total = length(plots.l))

    for (indexAll in seq_len(length(plots.l))) {
        
        pb$tick()
        index = index + 1
        plotsNew.l[[index]] = plots.l[[indexAll]]
        
        # Print another page
        if (index %% nPlotsPerPage == 0) { ## print 8 plots on a page
            do.call(gridExtra::grid.arrange,  c(plotsNew.l, list(ncol = nCol, nrow = nRow)))
            plotsNew.l = list() # reset plot 
            index = 0 # reset index
        }
        
    }
    
    # Print the remaining plots in case they don't perfectly fit with the layout
    if (length(plotsNew.l) != 0) { 
        do.call(gridExtra::grid.arrange,  c(plotsNew.l, list(ncol = nCol, nrow = nRow)))
    }
    
    if (!checkmate::testNull(pdfFile)) 
      grDevices::dev.off()
    
    futile.logger::flog.info(paste0("Finished writing plots to file ", pdfFile))
}

.startLogger <- function(logfile, level, removeOldLog = TRUE, appenderName = "consoleAndFile", verbose = TRUE) {
  
  checkmate::assertSubset(level, c("TRACE", "DEBUG", "INFO", "WARN", "ERROR", "FATAL"))
  checkmate::assertFlag(removeOldLog)
  checkmate::assertSubset(appenderName, c("console", "file", "consoleAndFile"))
  checkmate::assertFlag(verbose)
  
  if (appenderName != "console") {
    checkmate::assertDirectory(dirname(logfile), access = "w")
    if (file.exists(logfile)) {
      file.remove(logfile)
    }
  }
  
  # LEVELS: TRACE, DEBUG, INFO, WARN, ERROR, FATAL
  invisible(futile.logger::flog.threshold(level))
  
  
  if (appenderName == "console") {
    invisible(futile.logger::flog.appender(futile.logger::appender.console()))
  } else if (appenderName == "file")  {
    invisible(futile.logger::flog.appender(futile.logger::appender.file(file = logfile)))
  } else {
    invisible(futile.logger::flog.appender(futile.logger::appender.tee(file = logfile)))
  }
  
  
}

.printParametersLog <- function(par.l, verbose = FALSE) {
  
  checkmate::assertList(par.l)
  futile.logger::flog.info(paste0("PARAMETERS:"))
  for (parCur in names(par.l)) {
    
    futile.logger::flog.info(paste0(" ", parCur, "=",  paste0(par.l[[parCur]], collapse = ",")))
    
  }
}

###########################################
# PACKAGE LOADING AND DETACHING FUNCTIONS #
###########################################

.checkOntologyPackageInstallation <- function(pkg) {
    if (!is.installed(pkg)) {
        message = paste0("The package ", pkg, " is not installed, which is however needed for the chosen ontology enrichment. Please install it and re-run this function or change the ontology.")
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
}

.clearOpenDevices <- function() {
  
  while (length(grDevices::dev.list()) > 0) {
    grDevices::dev.off()
  }
}

.createFileList <- function(directory, pattern, recursive = FALSE, ignoreCase = FALSE, verbose = TRUE) {
  
  checkmate::assertCharacter(directory, min.chars = 1, any.missing = FALSE, len = 1)
  checkmate::assertCharacter(pattern, min.chars = 1, any.missing = FALSE, len = 1)
  checkmate::assertFlag(recursive)  
  checkmate::assertFlag(ignoreCase)  
  checkmate::assertFlag(verbose)    
  
  checkmate::assertDirectoryExists(directory)
  
  # Multiple patterns are now supported, integrate over them
  patternAll = strsplit(pattern, ",")[[1]]
  checkmate::assertCharacter(patternAll, min.len = 1)
  
  if(!verbose) futile.logger::flog.threshold(futile.logger::WARN)
  futile.logger::flog.info(paste0("Found ", length(patternAll), " distinct pattern(s) in pattern string."))
  
  nFilesToProcessTotal = 0
  filesToProcess.vec = c()
  
  for (patternCur in patternAll) {
    
    # Replace wildcards by functioning patterns (such as .)
    patternMod = glob2rx(patternCur)
    
    # Remove anchoring at beginning and end
    patternMod = substr(patternMod, 2, nchar(patternMod) - 1)
    
    filesToProcessCur.vec = list.files(path = directory, pattern = patternMod, full.names = TRUE, recursive = recursive, ignore.case = ignoreCase)
    filesToProcess.vec = c(filesToProcess.vec, filesToProcessCur.vec)
    
    futile.logger::flog.info(paste0("Search for files with pattern \"", patternCur, "\" in directory ", directory, " (case insensitive:", ignoreCase, ")"))
    
    nFilesToProcessTotal = nFilesToProcessTotal + length(filesToProcessCur.vec)
  }
  
  
  
  if (nFilesToProcessTotal == 0) {
    message = paste0("No files to process in folder ", directory, " that fulfill the desired criteria (", patternCur, ").")
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
   
  } else {
    
    futile.logger::flog.info(paste0("The following", nFilesToProcessTotal, "files were found:\n", paste0(filesToProcess.vec, collapse = "\n ")))
  }
  
  if(!verbose) futile.logger::flog.threshold(futile.logger::INFO)
  
  filesToProcess.vec
  
}

.testExistanceAndCreateDirectoriesRecursively <- function(directories) {
  
  for (dirname in unique(directories)) {
    
    if (!checkmate::testDirectoryExists(dirname)) {
      dir.create(dirname, recursive = TRUE)
    } else {
      checkmate::assertDirectoryExists(dirname, access = "w")
    }
    
  }
}

.initBiocParallel <- function(nWorkers, verbose = FALSE) {
  
  checkmate::assertInt(nWorkers, lower = 1)    
  checkmate::assertInt(BiocParallel::multicoreWorkers())
  
  if (nWorkers > BiocParallel::multicoreWorkers()) {
    warning("Requested ", nWorkers, " CPUs, but only ", BiocParallel::multicoreWorkers(), " are available and can be used.")
    nWorkers = BiocParallel::multicoreWorkers()
  }
  
  BiocParallel::MulticoreParam(workers = nWorkers, stop.on.error = TRUE)
  
}


.execInParallelGen <- function(nCores, returnAsList = TRUE, listNames = NULL, iteration, abortIfErrorParallel = TRUE, verbose = TRUE, functionName, ...) {
  
  start.time  <-  Sys.time()
  
  checkmate::assertInt(nCores, lower = 1)
  checkmate::assertFlag(returnAsList)
  checkmate::assertFunction(functionName)
  checkmate::assertVector(iteration, any.missing = FALSE, min.len = 1)
  checkmate::assert(checkmate::checkNull(listNames), checkmate::checkCharacter(listNames, len = length(iteration)))
  
  res.l = list()
  
  maxCores = tryCatch(
    {
       out = BiocParallel::multicoreWorkers() / 2
     
    },
    error=function() {

      message = "Could not retrieve the available number of cores. There might be a problem with your installation of BiocParallel. Check whether BiocParallel::multicoreWorkers() returns an integer. Try executing the following line to fix the problem if a re-installation of BiocParallel does not work: library(parallel) and then options(mc.cores=4L), with 4 being the maximum number of cores available for the machine you run the pipeline on For now, the function will just run with 1 core."
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
      return(1)
    }
  )    

  
  if (nCores > maxCores) {
      nCores = max(1, floor(maxCores))
      futile.logger::flog.warn(paste0(" Adjusted nCores down to ", nCores, " (only 50% of the cores are used at max)"))
  }
  
  if (nCores > 1) {
    
    res.l = tryCatch( {
        BiocParallel::bplapply(iteration, functionName, ..., BPPARAM = .initBiocParallel(nCores))
      
    }#, error = function(e) {
    #     warning("An error occured while executing the function with multiple CPUs. Trying again using only only one CPU...")
    #     lapply(iteration, functionName, ...)
    # }
    )
    
    failedTasks = which(!BiocParallel::bpok(res.l))
    if (length(failedTasks) > 0) {
      warning("At least one task failed while executing in parallel, attempting to rerun those that failed: ",res.l[[failedTasks[1]]])
      if (abortIfErrorParallel) stop()
      
      res.l = tryCatch( {
          BiocParallel::bplapply(iteration, functionName, ..., BPPARAM = .initBiocParallel(nCores), BPREDO = res.l)
        
      }, error = function(e) {
        warning("Once again, an error occured while executing the function with multiple CPUs. Trying again using only only one CPU...")
        if (abortIfErrorParallel) stop()
        lapply(iteration, functionName, ...)
      }
      )
    }
    
  } else {
    res.l = lapply(iteration, functionName, ...)
  }
  
  end.time  <-  Sys.time()
  
  
  futile.logger::flog.info(paste0(" Finished execution using ",nCores," cores. TOTAL RUNNING TIME: ", round(end.time - start.time, 1), " ", units(end.time - start.time),"\n"))
  
  
  if (!returnAsList) {
    return(unlist(res.l))
  }
  
  if (!is.null(listNames)) {
    names(res.l) = listNames
  }
  
  res.l
  
}

.printExecutionTime <- function(startTime, prefix = " ") {
  
  endTime  <-  Sys.time()
  futile.logger::flog.info(paste0(prefix, "Finished successfully. Execution time: ", round(endTime - startTime, 1), " ", units(endTime - startTime)))
}

#' @import grDevices
.checkAndLogWarningsAndErrors <- function(object, checkResult, isWarning = FALSE) {
  
  checkmate::assert(checkmate::checkCharacter(checkResult, len = 1), checkmate::checkLogical(checkResult))
  
  if (checkResult != TRUE) {

    objectPart = ""
    if (!is.null(object)) {
      objectname = deparse(substitute(object))
      objectPart = paste0("Assertion on variable \"", objectname, "\" failed: ")
    } 
    
    lastPartError   = "# An error occurred. See details above. If you think this is a bug, please contact us. #\n"
    hashesStrError = paste0(paste0(rep("#", nchar(lastPartError) - 1), collapse = ""), "\n")
    messageError    = paste0(objectPart, checkResult, "\n\n", hashesStrError, lastPartError, hashesStrError)
    
    lastPartWarning = ". \nThis warning may or may not be ignored. Carefully check its significance and whether it may affect the results."
    #hashesStrWarning = paste0(paste0(rep("#", nchar(lastPartWarning) - 1), collapse = ""), "\n")
    messageWarning  = paste0(objectPart, checkResult, lastPartWarning) # , hashesStrWarning)
    
    
    
    if (isWarning) {
      futile.logger::flog.warn(messageWarning)
      warning(messageWarning)
    } else {
      futile.logger::flog.error(messageError)
      # Close all open devices    
      while (grDevices::dev.cur()>1) grDevices::dev.off()
      stop(messageError)
    }
  }
}

.constructGRanges <- function(annotation, seqlengths = NULL, genomeAssembly, zeroBased = TRUE, stopIfError = TRUE, ...) {
  
  # Check if all start and end coordinates are valid
  sharedErrorMessage = paste0("This should not happen and may indicate a wrong genome assembly version. Check the validity of the input data and the parameters of the script. Are these really ", genomeAssembly, " coordinates? ")
    
    
  invalidEnd = which(annotation$end > seqlengths[as.character(annotation$chr)] | annotation$end < 0)
  if (length(invalidEnd) > 0) {
      
      message = paste0(length(invalidEnd), " sequences have ends outside of the chromosome boundaries (e.g., sequences ", 
                        paste0(invalidEnd[seq_len(min(5, length(invalidEnd)))], collapse = ","), ").", sharedErrorMessage)
      
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)    
   
  }
  
  # Check if all start and end coordinates are valid and inside the chromosomes
  invalidStart = which(annotation$start > seqlengths[as.character(annotation$chr)] | annotation$start < 0)
  if (length(invalidStart) > 0) {
      message = paste0(length(invalidStart), " sequences have ends outside of the chromosome boundaries (e.g., sequences ", 
                      paste0(invalidStart[seq_len(min(5, length(invalidStart)))], collapse = ","), ").", sharedErrorMessage)
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)    
  }
  
  gr = GenomicRanges::makeGRangesFromDataFrame(annotation, keep.extra.columns = TRUE, seqinfo = seqlengths, starts.in.df.are.0based = zeroBased, ...)
  
  # Check whether there are out-of-bound sequences, and abort if there are. This should not happen
  seq_outOfBound = which(end(gr) > GenomeInfoDb::seqlengths(gr)[as.character(GenomeInfoDb::seqnames(gr))])
  if (length(seq_outOfBound) > 0) {
    
      index = seq_len(min(length(seq_outOfBound),10))
      annotation_faultyEntries = paste0((GenomeInfoDb::seqnames(gr)[seq_outOfBound] %>% as.character())[index], ":", start(gr)[seq_outOfBound][index] , "-", end(gr)[seq_outOfBound][index], collapse = ",")

      message = paste0(length(seq_outOfBound), " sequences are outside of the chromosome boundaries (e.g., ", annotation_faultyEntries, "). This should not happen and may indicate either a wrong genome assembly version or that coordinates are not 0-based (see the Documentation). Check the validity of the input data and the parameters of the script. Are these really ", genomeAssembly, " coordinates?")
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = !stopIfError)
      
      
  }
  
  if (!is.null(seqlengths)) {
      
      idx <- which(end(gr) > GenomeInfoDb::seqlengths(gr)[as.character(GenomeInfoDb::seqnames(gr))] | start(gr) <= 0)
      if (length(idx) != 0L) {
          futile.logger::flog.warn(paste0(length(idx), " ranges out of ", length(gr), " are outside of the chromosome boundaries and will be removed. Check the validity of the input data."))
          gr = GenomicRanges::trim(gr)
      }
          
  }
  # 
  # nExtraColumns = length(extraColumns)
  # 
  # if (nExtraColumns > 0) {
  #     
  #     values(gr) <-as(annotation[, extraColumns], "DataFrame")
  # }
  # 

  gr
}

.checkAndLoadPackages <- function(packageList, verbose = FALSE) {
    
    for (packageCur in packageList) {
        if (!is.installed(packageCur)) {
            message = paste0("The suggested package \"", packageCur, "\" is not yet installed, but it is required for this function and in this context. Please install it and re-run this function.")
            .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
        }
    }
    
}

.getGenomeObject <- function(genomeAssembly, type = "txbd") {
    
    
    checkmate::assertSubset(type, c("txbd", "BSgenome", "packageName"))
    checkmate::assertSubset(genomeAssembly, c("hg19","hg38", "mm9", "mm10"))
    
    
    if (genomeAssembly == "hg38") {

        if (type == "txbd") {
            .checkAndLoadPackages(c("org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene"), verbose = FALSE) 
            obj <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
        } else if (type == "BSgenome") {
            .checkAndLoadPackages(c("BSgenome.Hsapiens.UCSC.hg38"), verbose = FALSE) 
            obj <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
        } else {
            obj = "org.Hs.eg.db"
        }
        
    } else if (genomeAssembly == "hg19") {
        
        if (type == "txbd") {
            .checkAndLoadPackages(c("org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene"), verbose = FALSE) 
            obj <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
        } else if (type == "BSgenome") {
            .checkAndLoadPackages(c("BSgenome.Hsapiens.UCSC.hg19"), verbose = FALSE) 
            obj <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
        } else {
            obj = "org.Hs.eg.db"
        }
        
        
    } else if (genomeAssembly == "mm10") {
        
        if (type == "txbd") { 
            .checkAndLoadPackages(c("org.Mm.eg.db", "TxDb.Mmusculus.UCSC.mm10.knownGene"), verbose = FALSE) 
            obj <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
        } else if (type == "BSgenome") {
            .checkAndLoadPackages(c("BSgenome.Mmusculus.UCSC.mm10"), verbose = FALSE) 
            obj <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
        } else {
            obj = "org.Mm.eg.db"
        }
        
    } else if (genomeAssembly == "mm9") {
        
        if (type == "txbd") {
            .checkAndLoadPackages(c("org.Mm.eg.db", "TxDb.Mmusculus.UCSC.mm9.knownGene"), verbose = FALSE) 
            obj <- TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene
        } else if (type == "BSgenome") {
            .checkAndLoadPackages(c("BSgenome.Mmusculus.UCSC.mm9"), verbose = FALSE) 
            obj <- BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9
        } else {
            obj = "org.Mm.eg.db"
        }
        
    }
    
    obj
}

.getChrLengths <- function(genomeAssembly) {
  txdb = .getGenomeObject(genomeAssembly)
  GenomeInfoDb::seqlengths(txdb)
}

.checkAndLoadPackagesGenomeAssembly <- function(genomeAssembly) {
    
    checkmate::assertSubset(genomeAssembly, c("hg19","hg38", "mm9", "mm10"))
    
    if (genomeAssembly == "hg38") {
        
        packages = c("org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene")
        
    } else if (genomeAssembly == "hg19") {
        
        packages = c("org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene")
        
    } else if (genomeAssembly == "mm10") {
        
        packages = c("org.Mm.eg.db", "TxDb.Mmusculus.UCSC.mm10.knownGene")
        
    } else if (genomeAssembly == "mm9") {
        
        packages = c("org.Mm.eg.db", "TxDb.Mmusculus.UCSC.mm9.knownGene")
        
    }
    
    packages
}

.shuffleColumns <- function(df, nPermutations, returnUnshuffled = TRUE, returnAsList = TRUE, saveMemory = TRUE) {
    
    start = Sys.time()
   
    df.shuffled.l = list()
    
    futile.logger::flog.info(paste0("Shuffling columns ", nPermutations, " times"))
    
    # Shuffle RNA-Seq data for permutations
    for (permutationCur in 0:nPermutations) {
        
        index = as.character(permutationCur)
        
        if (permutationCur == 0 & returnUnshuffled) {
            df.shuffled.l[[index]]  = df
        } else if (permutationCur > 0) {
          
            sampleOrder = sample(ncol(df))
            
            if (saveMemory) {
              df.shuffled.l[[index]] =  sampleOrder
            } else {
              df.shuffled.l[[index]]  = df[,sampleOrder]
            }
            
        } else {
          # nothing to do
        }
        
    }

    .printExecutionTime(start)
    
    if (nPermutations == 1 & !returnAsList & !returnUnshuffled) {
      df.shuffled.l[["1"]]
    } else {
      df.shuffled.l
    }
    
}

.shuffleRowsPerColumn <- function(df) {
  
  start = Sys.time()
  futile.logger::flog.info(paste0("Shuffling rows per column"))

  
  for (i in seq_len(ncol(df))) {
    df[,i] = (dplyr::select(df, tidyselect::all_of(i)) %>% dplyr::pull()) [sample(nrow(df))]
  }
  
  .printExecutionTime(start)
  
  df
  
}

.checkSelfOverlap <- function(subject) {
    
    overlapsCount = GenomicRanges::countOverlaps(subject, subject, 
                                  minoverlap=1,
                                  type="any",
                                  ignore.strand=TRUE)
    
    length(which(overlapsCount > 1))
    
    
}


.asSparseMatrix <- function(matrix, convertNA_to_zero=TRUE, dimnames = NULL) {
  
  if (convertNA_to_zero) {
    matrix[which(is.na(matrix))]= 0
  }
 
  Matrix::Matrix(matrix, sparse=TRUE, dimnames = dimnames)
}


.asMatrixFromSparse <-function(matrix, convertZero_to_NA=TRUE) {
  
  if (methods::is(matrix,"dgeMatrix") | methods::is(matrix,"dgCMatrix") | methods::is(matrix,"dgRMatrix")) {
    dimNames = dimnames(matrix)
    matrix = as.matrix(matrix)
    
    if (convertZero_to_NA) {
        matrix[which(matrix == 0)]= NA
    }
    
    dimnames(matrix) = dimNames

  }

  matrix
  
}


.prettyNum <- function(number, verbose = TRUE) {
  prettyNum(number, big.mark = ",", scientific = FALSE)
}


.read_tidyverse_wrapper <- function(file, type = "tsv", ncolExpected = NULL, minRows = 0, verbose = TRUE, ...) {
  
  checkmate::assertSubset(type, c("csv", "csv2", "tsv", "delim"))
  
  start = Sys.time()
  if (verbose) futile.logger::flog.info(paste0("Reading file ", file))
  
  
  if (type == "tsv") {
    tbl = readr::read_tsv(file, ...)
  } else if (type == "csv") {
    tbl = readr::read_csv(file, ...)
  } else if (type == "csv2") {
    tbl = readr::read_csv2(file, ...)
  } else if (type == "delim") {
    tbl = readr::read_delim(file, ...)
  }
  
  
  if (nrow(readr::problems(tbl)) > 0) {
    futile.logger::flog.fatal(paste0("Parsing errors: "), readr::problems(tbl), capture = TRUE)
    stop("Error when parsing the file ", file, ", see errors above")
  }
  
  
  if (nrow(tbl) == 0) {
    message = paste0("The file ", file, " is unexpectedly empty.")
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
  if (!is.null(ncolExpected)) {
    if (! ncol(tbl) %in% ncolExpected) {
      message = paste0("The file ", file, " does not have the expected number of ", ncolExpected, " columns, but instead ", ncol(tbl), ".")
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
  }
  
  if (minRows > 0) {
    if (nrow(tbl) < minRows) {
      message = paste0("The file ", file, " does not have the expected minimum number of rows. Expected at least ", minRows, ", but found only ",nrow(tbl), ".")
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
  }
  
  
  if (verbose) .printExecutionTime(start)
  tbl
}


isIntegerMatrix <- function(df) {
  
  res = all.equal(unlist(df), as.integer(unlist(df)), check.attributes = FALSE)
  
  if (is.logical(res)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

.findSuitableColumnsToPlot <- function(metadataTable, remove1LevelFactors = TRUE, verbose = TRUE) {
  
  # Convert character to factor in metadataTable
  metadataTable = dplyr::mutate_if(metadataTable, is.character, as.factor)
  
  # Remove NA only columns
  metadataTable = metadataTable %>%  dplyr::select_if(~sum(!is.na(.)) > 0)
  metadataColumns = colnames(metadataTable) 
  
  # Eliminate 1-factor levels
  if (remove1LevelFactors) {
    columnsToDelete = metadataTable[, sapply(metadataTable, nlevels) == 1, drop = FALSE] %>% colnames()
    if (length(columnsToDelete ) > 0 & verbose) {
      futile.logger::flog.info(paste0("Exclude the following columns because they do not have at least 2 unique values: ", paste0(columnsToDelete, collapse = ",")))
    }
    
    finalColumns = unique(setdiff(metadataColumns, columnsToDelete))
  } else {
    finalColumns = colnames(metadataTable) 
  }
  
  finalColumns
  
}

.getPermStr <- function (permutation) {
  
  dplyr::if_else(permutation == 0, "Real data", "Permuted data")
}


# Helper function to record function arguments
match.call.defaults <- function(asList = TRUE, ...) {
    call <- evalq(match.call(expand.dots = FALSE), parent.frame(2))
    formals <- evalq(formals(), parent.frame(2))
    
    for(i in setdiff(names(formals), names(call)))
        call[i] <- list( formals[[i]] )
    
    
    if (asList) {
        as.list(match.call(sys.function(sys.parent(n = 2)), call)[-1])
    } else {
        match.call(sys.function(sys.parent(n = 2)), call)
    }
}

.firstLetterUppercase <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
}

#' @import utils
is.installed <- function(mypkg){
    is.element(mypkg, installed.packages()[,1])
} 