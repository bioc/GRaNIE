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
    
    futile.logger::flog.info(paste0(" Finished writing plots"))
}

#' @import AnnotationHub
.getAnnotationHub <- function(curAttempt = 1, maxAttempts = 10) {
    
    if (curAttempt < maxAttempts) {
        
        Sys.sleep(curAttempt * 3)
        
        ah <- tryCatch({ 
            AnnotationHub::AnnotationHub()
            
        }, error = function(e) {
            
            futile.logger::flog.warn(paste0("An error occured for AnnotationHub::AnnotationHub(). The error message was: ", e))
            
            # See https://rdrr.io/bioc/AnnotationHub/f/vignettes/TroubleshootingTheCache.Rmd
            futile.logger::flog.info("Trying to fix automatically and re-generate the cache, this may take a while...")
            cache_dir = tools::R_user_dir("AnnotationHub", which = "cache") 
            unlink(cache_dir, recursive = TRUE)
            
            AnnotationHub::getAnnotationHub(curAttempt = curAttempt + 1, maxAttempts = maxAttempts)
            # bfc <- BiocFileCache::BiocFileCache(cache_dir)
            # res <- BiocFileCache::bfcquery(bfc, "annotationhub.index.rds", field = "rname", exact = TRUE)
            # BiocFileCache::bfcremove(bfc, rids = res$rid)
            # AnnotationHub::AnnotationHub()
        })
        
        if (is(ah, "AnnotationHub")) {
            return(ah)
        }
        
        
    }

    

}

.startLogger <- function(logfile, level, removeOldLog = TRUE, appenderName = "consoleAndFile", verbose = TRUE) {
  
  checkmate::assertChoice(level, c("TRACE", "DEBUG", "INFO", "WARN", "ERROR", "FATAL"))
  checkmate::assertFlag(removeOldLog)
  checkmate::assertChoice(appenderName, c("console", "file", "consoleAndFile"))
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

.checkPackageInstallation <- function(pkg, message, isWarning = FALSE, returnLogical = FALSE) {
    
    pkgMissing = c()
    for (packageCur in pkg) {
        if (!is.installed(packageCur)) {
            pkgMissing = c(pkgMissing, packageCur)
        }
    }
    
    if (length(pkgMissing) > 0) {
        message = paste0(message, "\n\nExecute the following line in R to install the missing package(s): `BiocManager::install(c(\"", 
                         paste0(pkgMissing, collapse = "\",\""), 
                         "\"))`")
        
        # Normal behavior
        if (!returnLogical) {
            .checkAndLogWarningsAndErrors(NULL, message, isWarning = isWarning)
        } else {
            
            # Make it just a warning but return FALSE
            .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
            return(FALSE)
        }
        
    } else {
        if (returnLogical) return(TRUE)
    }
   
    
}


.checkAndLoadPackagesGenomeAssembly <- function(genomeAssembly, returnLogical = FALSE) {
    
    baseMessage = paste0("For the chosen genome assembly version and package function, particular genome annotation packages are required but not installed. Please install and re-run this function.")

    if (genomeAssembly == "hg38") {
        res = .checkPackageInstallation(c("org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene", "BSgenome.Hsapiens.UCSC.hg38"), baseMessage, returnLogical = returnLogical)
    } else if (genomeAssembly == "hg19") {
        res = .checkPackageInstallation(c("org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene", "BSgenome.Hsapiens.UCSC.hg19"), baseMessage, returnLogical = returnLogical)
    } else if (genomeAssembly == "mm10") {
        res = .checkPackageInstallation(c("org.Mm.eg.db", "TxDb.Mmusculus.UCSC.mm10.knownGene", "BSgenome.Mmusculus.UCSC.mm10"), baseMessage, returnLogical = returnLogical)
    } else if (genomeAssembly == "mm9") {
        res = .checkPackageInstallation(c("org.Mm.eg.db", "TxDb.Mmusculus.UCSC.mm9.knownGene", "BSgenome.Mmusculus.UCSC.mm9"), baseMessage, returnLogical = returnLogical)
    } else if (genomeAssembly == "rn6") {
        res = .checkPackageInstallation(c("org.Rn.eg.db", "TxDb.Rnorvegicus.UCSC.rn6.refGene", "BSgenome.Rnorvegicus.UCSC.rn6"), baseMessage, returnLogical = returnLogical)
    } else if (genomeAssembly == "rn7") {
        res = .checkPackageInstallation(c("org.Rn.eg.db", "TxDb.Rnorvegicus.UCSC.rn7.refGene", "BSgenome.Rnorvegicus.UCSC.rn7"), baseMessage, returnLogical = returnLogical)
    } else if (genomeAssembly == "dm6") {
        res = .checkPackageInstallation(c("org.Dm.eg.db", "TxDb.Dmelanogaster.UCSC.dm6.ensGene", "BSgenome.Dmelanogaster.UCSC.dm6"), baseMessage, returnLogical = returnLogical)
    } else if (genomeAssembly == "rheMac10") {
        res = .checkPackageInstallation(c("org.Mmu.eg.db", "TxDb.Mmulatta.UCSC.rheMac10.refGene", "BSgenome.Mmulatta.UCSC.rheMac10"), baseMessage, returnLogical = returnLogical)
    } else {
        message = "Genome not listed, this should not happen. Contact the author."
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    if (returnLogical) return(res)

}


.checkPackage_topGO_and_arguments <- function(ontology, algorithm, statistic) {
    
    if (length(intersect(ontology, c("GO_BP", "GO_MF", "GO_CC"))) > 0) {
        
        packageMessage = paste0("The package topGO is not installed but required when selecting any of the three following ontologies: GO_BP, GO_MF, GO_CC. Please install it and re-run this function or change the ontology.")
        .checkPackageInstallation("topGO", packageMessage) 
        
        # This function is calling the topGO:::.onAttach function and needs to be executed once otherwise
        # errors like object 'GOBPTerm' of mode 'environment' was not found are thrown
        suppressMessages(topGO::groupGOTerms()) 
        
        checkmate::assertChoice(algorithm , topGO::whichAlgorithms())
        checkmate::assertChoice(statistic , topGO::whichTests())
        
        # Some statistics do not seem to work properly and cause errors in topGO, we exclude them here
        if (statistic %in% c("sum", "globaltest", "ks.ties")) {
            message = paste0("We stopped supporting the statsitic \"", statistic, "\" because it causes errors in topGO. Please choose a different statistic.")
            .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
        }
        
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
  
  if (!verbose) futile.logger::flog.threshold(futile.logger::WARN)
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
  
  if (!verbose) futile.logger::flog.threshold(futile.logger::INFO)
  
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


.execInParallelGen <- function(nCores, returnAsList = TRUE, listNames = NULL, iteration, abortIfErrorParallel = FALSE, verbose = TRUE, functionName, ...) {
  
  start.time  <-  Sys.time()
  
  checkmate::assertInt(nCores, lower = 1)
  checkmate::assertFlag(returnAsList)
  checkmate::assertFunction(functionName)
  checkmate::assertVector(iteration, any.missing = FALSE, min.len = 1)
  checkmate::assert(checkmate::checkNull(listNames), checkmate::checkCharacter(listNames, len = length(iteration)))
  
  res.l = list()
  

  if (nCores > 1) {
      
    message = paste0("The package BiocParallel is not installed but required if more than 1 core should be used as requested. Please install it and re-run this function to speed-up the computation time or set the nCores parameter to 1 to disable parallel computation.")
    .checkPackageInstallation("BiocParallel", message)
    
    checkFailedTasks = TRUE
    
    res.l = tryCatch( {
        BiocParallel::bplapply(iteration, functionName, ..., BPPARAM = .initBiocParallel(nCores))
      
    }, error = function(e) {
        futile.logger::flog.warn(paste0("The following error occured while executing the function with multiple CPUs using BiocParallel: ", e, ". Trying again using only 1 core..."))
        checkFailedTasks = FALSE
        lapply(iteration, functionName, ...)
     }
    )
    
    failedTasks = which(!BiocParallel::bpok(res.l))
    if (checkFailedTasks && length(failedTasks) > 0) {
        futile.logger::flog.warn(paste0("At least one task failed while executing in parallel, attempting to rerun those that failed: ", res.l[[failedTasks[1]]]))
        
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
.checkAndLogWarningsAndErrors <- function(object, message, isWarning = FALSE) {
  
  checkmate::assert(checkmate::checkCharacter(message, len = 1), checkmate::checkLogical(message))
  
  if (message != TRUE) {

    objectPart = ""
    if (!is.null(object)) {
      objectname = deparse(substitute(object))
      objectPart = paste0("Assertion on variable \"", objectname, "\" failed: ")
    } 
    
    lastPartError   = "# An error occurred. See details above. If you think this is a bug, please contact us. #\n"
    hashesStrError = paste0(paste0(rep("#", nchar(lastPartError) - 1), collapse = ""), "\n")
    messageError    = paste0(objectPart, message, "\n\n", hashesStrError, lastPartError, hashesStrError)
    
    lastPartWarning = ". \nThis warning may or may not be ignored. Carefully check its significance and whether it may affect the results.\n"
    #hashesStrWarning = paste0(paste0(rep("#", nchar(lastPartWarning) - 1), collapse = ""), "\n")
    messageWarning  = paste0(objectPart, message, lastPartWarning) # , hashesStrWarning)
    
    
    
    if (isWarning) {
      futile.logger::flog.warn(messageWarning)
      warning(messageWarning)
    } else {
      futile.logger::flog.error(messageError)
      # Close all open devices    
      while (grDevices::dev.cur() > 1) grDevices::dev.off()
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
  
  # Check if there are any chr names that are not part of seqlengths, which causes an 'seqnames' contains sequence names with no entries in 'seqinfo' error
  missingSeqs = unique(annotation$chr)[which(!unique(annotation$chr) %in% names(seqlengths))]
  
  if (length(missingSeqs) > 0) {
      message = paste0("For ", length(missingSeqs), " chromosomes (", paste0(missingSeqs, collapse = ","), ") and a total of ", length(which(annotation$chr %in% missingSeqs)), 
                       " peaks, their length was not found in biomaRt. These peaks will be discarded")
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)    
      
      # Filter and refactor
      annotation = annotation %>%
          dplyr::filter(!.data$chr %in% missingSeqs) %>%
          dplyr::mutate(chr = factor(.data$chr, levels = unique(.data$chr)))
      
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


# Only needed here: get CG content peaks (can be made optional), and .populatePeakAnnotation (within ChipSeeker)
.getGenomeObject <- function(genomeAssembly, type = "txbd", jasparRelease = 2022) {
    
    checkmate::assertChoice(type, c("txbd", "BSgenome", "packageName", "txID")) #txID: NCBI taxonomy ID
    
    if (type != "txID") {
        checkmate::assertChoice(genomeAssembly, c("hg19","hg38", "mm9", "mm10", "rn6", "rn7", "dm6", "rheMac10"))
    } else {
        availableSpecies.df = rbioapi::rba_jaspar_species(release = jasparRelease)
    }
    
    if (genomeAssembly == "hg38") {
        
        if (type == "txbd") {
            obj <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
        } else if (type == "BSgenome") {
            obj <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
        } else if (type == "txID") {
            obj = availableSpecies.df$tax_id[which(availableSpecies.df$species == "Homo sapiens")]
        } else {
            obj = "org.Hs.eg.db"
        }
        
    } else if (genomeAssembly == "hg19") {
        
        if (type == "txbd") {
            obj <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
        } else if (type == "BSgenome") {
            obj <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
        } else if (type == "txID") {
            obj = availableSpecies.df$tax_id[which(availableSpecies.df$species == "Homo sapiens")]
        } else {
            obj = "org.Hs.eg.db"
        }
        
        
    } else if (genomeAssembly == "mm10") {
        
        if (type == "txbd") { 
            obj <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
        } else if (type == "BSgenome") {
            obj <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
        } else if (type == "txID") {
            obj = availableSpecies.df$tax_id[which(availableSpecies.df$species == "Mus musculus")]
        } else {
            obj = "org.Mm.eg.db"
        }
        
    } else if (genomeAssembly == "mm9") {
        
        if (type == "txbd") {
            obj <- TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene
        } else if (type == "BSgenome") {
            obj <- BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9
        } else if (type == "txID") {
            obj = availableSpecies.df$tax_id[which(availableSpecies.df$species == "Mus musculus")]
        } else {
            obj = "org.Mm.eg.db"
        }
        
    } else if (genomeAssembly == "rn6") {
        
        if (type == "txbd") {
            obj <- TxDb.Rnorvegicus.UCSC.rn6.refGene::TxDb.Rnorvegicus.UCSC.rn6.refGene
        } else if (type == "BSgenome") {
            obj <- BSgenome.Rnorvegicus.UCSC.rn6::BSgenome.Rnorvegicus.UCSC.rn6
        } else if (type == "txID") {
            obj = availableSpecies.df$tax_id[which(availableSpecies.df$species == "Rattus norvegicus")]
        } else {
            obj = "org.Rn.eg.db"
        }
        
    }else if (genomeAssembly == "rn7") {
        
        if (type == "txbd") {
            obj <- TxDb.Rnorvegicus.UCSC.rn7.refGene::TxDb.Rnorvegicus.UCSC.rn7.refGene
        } else if (type == "BSgenome") {
            obj <- BSgenome.Rnorvegicus.UCSC.rn7::BSgenome.Rnorvegicus.UCSC.rn7
        } else if (type == "txID") {
            obj = availableSpecies.df$tax_id[which(availableSpecies.df$species == "Rattus norvegicus")]
        } else {
            obj = "org.Rn.eg.db"
        }
        
    } else if (genomeAssembly == "dm6") {
        
        if (type == "txbd") {
            obj <- TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene
        } else if (type == "BSgenome") {
            obj <- BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6
        } else if (type == "txID") {
            obj = availableSpecies.df$tax_id[which(availableSpecies.df$species == "Drosophila melanogaster")]
        } else {
            obj = "org.Dm.eg.db"
        }
        
    } else if (genomeAssembly == "rheMac10") {
        
        if (type == "txbd") {
            obj <- TxDb.Mmulatta.UCSC.rheMac10.refGene::TxDb.Mmulatta.UCSC.rheMac10.refGene
        } else if (type == "BSgenome") {
            obj <- BSgenome.Mmulatta.UCSC.rheMac10::BSgenome.Mmulatta.UCSC.rheMac10
        } else if (type == "txID") {
            
            # rheMac not available in species list, selecting here the more general Verterbates selection instead then
            obj = availableSpecies.df$tax_id[which(availableSpecies.df$species == "Vertebrata")]
        } else {
            obj = "org.Mmu.eg.db"
        }
        
    }
    
    obj
}



# .getChrLengths <- function(genomeAssembly) {
#   txdb = .getGenomeObject(genomeAssembly)
#   GenomeInfoDb::seqlengths(txdb)
# }

# Try an approach that is lighter and doesnt require large annotation packages

#' @import GenomeInfoDb
.getChrLengths <- function(genomeAssembly) {
    
    chrSizes = GenomeInfoDb::getChromInfoFromUCSC(genomeAssembly)

    chrSizes.vec = chrSizes$size
    names(chrSizes.vec) = chrSizes$chrom
    
    chrSizes.vec
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
                                  minoverlap = 1,
                                  type = "any",
                                  ignore.strand = TRUE)
    
    length(which(overlapsCount > 1))
    
    
}

#' @importFrom Matrix Matrix
.asSparseMatrix <- function(matrix, convertNA_to_zero=TRUE, dimnames = NULL) {
  
  if (convertNA_to_zero) {
    matrix[which(is.na(matrix))] = 0
  }
 
    # TODO: as("lgCMatrix") may be a bit more space-saing BUT then we cannot add the isFiltered column to it that the code currently relies on
  Matrix::Matrix(matrix, sparse = TRUE, dimnames = dimnames) %>% methods::as("dMatrix")
}

#' @importFrom Matrix Matrix
.asMatrixFromSparse <- function(matrix, convertZero_to_NA=TRUE) {
  
  if (methods::is(matrix,"dgeMatrix") | methods::is(matrix,"dgCMatrix") | methods::is(matrix,"dgRMatrix") | methods::is(matrix,"lgCMatrix")) {
    dimNames = dimnames(matrix)
    matrix = as.matrix(matrix)
    
    if (convertZero_to_NA) {
        matrix[which(matrix == 0)] = NA
    }
    
    dimnames(matrix) = dimNames

  } else {
      # 
      # if (!is.matrix(matrix) & !is.data.frame(matrix)) {
      #     message = paste0("Unknown sparse matrix type \"", class(matrix), "\", this should not happen. Contact the authors.")
      #     .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
      # }
     
  }

  matrix
  
}


.prettyNum <- function(number, verbose = TRUE) {
  prettyNum(number, big.mark = ",", scientific = FALSE)
}


.read_tidyverse_wrapper <- function(file, type = "tsv", ncolExpected = NULL, minRows = 0, verbose = TRUE, ...) {
  
  checkmate::assertChoice(type, c("csv", "csv2", "tsv", "delim"))
  
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
    if (!ncol(tbl) %in% ncolExpected) {
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
    
  # Quick test first: test only first row.
  res = all.equal(unlist(df[1,]), as.integer(unlist(df[1,])), check.attributes = FALSE)  
  if (!is.logical(res)) {
      return(FALSE)
  } else {
      res = all.equal(unlist(df), as.integer(unlist(df)), check.attributes = FALSE)
      if (is.logical(res)) {
          return(TRUE)
      } else {
          return(FALSE)
      }
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

.getPermStr <- function(permutation) {
  
  dplyr::if_else(permutation == 0, "real data", "permuted data")
}


# Helper function to record function arguments
match.call.defaults <- function(asList = TRUE, ...) {
    call <- evalq(match.call(expand.dots = FALSE), parent.frame(2))
    formals <- evalq(formals(), parent.frame(2))
    
    for (i in setdiff(names(formals), names(call)))
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

# Taken from https://www.bioconductor.org/packages/release/bioc/vignettes/BiocFileCache/inst/doc/BiocFileCache.html
.get_cache <- function() {
    # cache <- tools::R_user_dir(utils::packageName(), which="cache")
    cache <- tools::R_user_dir("GRaNIE", which = "cache")
    BiocFileCache::BiocFileCache(cache, ask = FALSE)
}

.checkOutputFile <- function(fileCur) {
  
  if (!checkmate::testDirectory(dirname(fileCur), access = "w")) {
    message = paste0("The specified output file or directory (", fileCur, ") is not writable or does not exist. Please adjust the output folder and rerun the function or change the output folder globally using the function changeOutputDirectory.")
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
}

.getBiomartParameters <- function(genomeAssembly, suffix = "") {
    
    host = "https://www.ensembl.org"
    
    if (grepl(x = genomeAssembly, pattern = "^hg\\d+")) {
        dataset = "hsapiens"
        if (genomeAssembly == "hg38") {
        } else if (genomeAssembly == "hg19") {
            host = "https://grch37.ensembl.org"
        }
    } else if (grepl(x = genomeAssembly, pattern = "^mm\\d+")) {
        dataset = "mmusculus"
    } else if (grepl(x = genomeAssembly, pattern = "^rn\\d+")) {
        dataset = "rnorvegicus"
    } else if (grepl(x = genomeAssembly, pattern = "^dm\\d+")) {
        dataset = "dmelanogaster"
    }
    
    
    
    list(dataset = paste0(dataset, suffix), host = host)
}


.biomart_getEnsembl <- function(biomart, version, host, dataset, maxAttempts = 40) {
    
    ensembl <- NULL
    mirrorIndex <- 0
    attemptsDone = 0
    
    mirrors = c('www', 'uswest', 'useast', 'asia')
    while (!"Mart" %in% class(ensembl) && attemptsDone <= maxAttempts ) {
        mirrorIndex <- (mirrorIndex %% 4) + 1
        attemptsDone = attemptsDone + 1
        
        ensembl = tryCatch({ 
            biomaRt::useEnsembl(biomart = biomart , version = version, host = host,  
                                dataset = dataset, mirror = mirrors[mirrorIndex])
            
            
        }, error = function(e) {
            futile.logger::flog.warn(paste0("Attempt ", attemptsDone, " out of ", maxAttempts, " failed. Another automatic attempt will be performed using a different mirror. The error message was: ", e))
        })
    } 
    
    if (!"Mart" %in% class(ensembl)) {
        
        error_Biomart = paste0("A temporary error occured with biomaRt::getBM or biomaRt::useEnsembl despite trying ", 
                               maxAttempts, 
                               " times via different mirrors. This is often caused by an unresponsive Ensembl site.", 
                               " Try again at a later time. Note that this error is not caused by GRaNIE but external services.")
        .checkAndLogWarningsAndErrors(NULL, error_Biomart, isWarning = FALSE)
        return(NULL)
        
    } 
    
    futile.logger::flog.info(paste0(" Retrieving BioMart database succeeded"))
    
    ensembl
}

.callBiomart <- function(mart, attributes = NULL, values = "", filters = "", maxAttempts = 40) {
    
    result.df <- NULL
    attemptsDone = 0

    while (!is.data.frame(result.df) && attemptsDone <= maxAttempts ) {
        attemptsDone = attemptsDone + 1
        
        result.df = tryCatch({ 
            biomaRt::getBM(attributes = attributes,
                           filters = filters,
                           values = values,
                           mart = mart) 
            
        }, error = function(e) {
            futile.logger::flog.warn(paste0("Attempt ", attemptsDone, " out of ", maxAttempts, " failed. Another automatic attempt will be performed using a different mirror. The error message was: ", e))
        })
    } 
    
    if (!is.data.frame(result.df)) {
        
        error_Biomart = paste0("A temporary error occured with biomaRt::getBM or biomaRt::useEnsembl despite trying ", 
                               maxAttempts, 
                               " times via different mirrors. This is often caused by an unresponsive Ensembl site.", 
                               " Try again at a later time. Note that this error is not caused by GRaNIE but external services.")
        .checkAndLogWarningsAndErrors(NULL, error_Biomart, isWarning = FALSE)
        return(NULL)
        
    } 
    
    futile.logger::flog.info(paste0(" Retrieving genome annotation succeeded"))
    
    result.df
    
    
}


.correlateData <- function(x, y, corMethod) {
    
    if (corMethod %in% c("pearson", "spearman")) {
        
        cor.res = suppressWarnings(cor(x,y, method = corMethod))
        
    } else if (corMethod == "bicor") {
        
        cor.res = suppressWarnings(WGCNA::bicor(x, y))
    }
    
    cor.res
}