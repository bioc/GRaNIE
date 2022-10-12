######## Init GRN ########

#' Create and initialize a \code{\linkS4class{GRN}} object.
#' 
#' Executing this function is the very first step in the *GRaNIE* workflow. After its execution, data can be added to the object. 
#' \strong{Depending on the genome assembly version, additional genome annotation packages are required, as follows:} 
#' For \code{hg19} and \code{hg38},
#' the packages \code{org.Hs.eg.db} as well as \code{BSgenome.Hsapiens.UCSC.hg19}+\code{TxDb.Hsapiens.UCSC.hg19.knownGene} or 
#' \code{BSgenome.Hsapiens.UCSC.hg38}+\code{TxDb.Hsapiens.UCSC.hg38.knownGene} are required, respectively. 
#' For \code{mm9} and \code{mm10},
#' the packages \code{org.Mm.eg.db} as well as \code{BSgenome.Mmusculus.UCSC.mm9}+\code{TxDb.Mmusculus.UCSC.mm9.knownGene} or 
#' \code{Mmusculus.UCSC.mm10}+\code{TxDb.Mmusculus.UCSC.mm10.knownGene} are required, respectively.
#' For more information, see the error message if any of these packages is missing or the 
#' \href{https://grp-zaugg.embl-community.io/GRaNIE/articles/GRaNIE_packageDetails.html#installation}{Package Details Vignette}.
#' 
#' @export
#' @param objectMetadata List. Default \code{list()}. Optional (named) list with an arbitrary number of elements, all of which 
#' capture metadata for the object. This is mainly used to distinguish GRN objects from one another by storing object-specific metadata along with the data.
#' @param outputFolder Output folder, either absolute or relative to the current working directory. Default \code{"."}. 
#' Default output folder where all pipeline output will be put unless specified otherwise. We recommend specifying an absolute path. 
#' Note that for Windows-based systems, the path must be correctly specified with "/" as path separator.
#' @param genomeAssembly Character. No default. The genome assembly of all data that to be used within this object. 
#' Currently, supported genomes are: \code{hg19}, \code{hg38}, \code{mm9} and \code{mm10}. See function description for further information and notes.
#' @return Empty \code{\linkS4class{GRN}} object
#' @examples 
#' meta.l = list(name = "exampleName", date = "01.03.22")
#' GRN = initializeGRN(objectMetadata = meta.l, outputFolder = "output", genomeAssembly = "hg38")
#' @export
#' @importFrom stats sd median cor cor.test quantile
initializeGRN <- function(objectMetadata = list(),
                          outputFolder = ".", 
                          genomeAssembly) {
  
  start = Sys.time()   
    
  checkmate::assert(checkmate::checkNull(objectMetadata), checkmate::checkList(objectMetadata))
  checkmate::assertChoice(genomeAssembly, c("hg19","hg38", "mm9", "mm10"))

  .checkAndLoadPackagesGenomeAssembly(genomeAssembly)
  
  # Create the folder first if not yet existing
  checkmate::assertCharacter(outputFolder, min.chars = 1, len = 1)
  if (!dir.exists(outputFolder)) {
    res = dir.create(outputFolder)
    if (!res) {
        message = paste0("Could not create the output directory ", outputFolder, ". Check the path /and/or access rights.")
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    checkmate::assertDirectoryExists(outputFolder, access = "w")
  }
  
  # Create an absolute path out of the given outputFolder now that it exists
  outputFolder = tools::file_path_as_absolute(outputFolder)
  checkmate::assertDirectory(outputFolder, access = "w")
  
  if (!endsWith(outputFolder, .Platform$file.sep)) {
    outputFolder = paste0(outputFolder, .Platform$file.sep)
  }
  
  
  dir_output_plots = paste0(outputFolder, "plots", .Platform$file.sep)
  if (!dir.exists(dir_output_plots)) {
    dir.create(dir_output_plots)
  }
  
  GRN = .createGRNObject()
  GRN@config$functionParameters = list()
  
  GRN = .addFunctionLogToObject(GRN)
  
  GRN@config$isFiltered = FALSE
  
  par.l = list()
  
  packageName = utils::packageName()
  par.l$packageVersion = ifelse(is.null(packageName), NA, paste0(utils::packageVersion(packageName)[[1]], collapse = "."))
  par.l$genomeAssembly = genomeAssembly

  
  # Make an internal subslot
  par.l$internal = list()
  
  # Recommended to leave at 1, more permutations are currently not necessary
  par.l$internal$nPermutations = 1 
  
  # Step size for the TF-peak FDR calculation
  par.l$internal$stepsFDR = round(seq(from = -1, to = 1, by = 0.05),2)
  
  # Stringencies for AR classification
  par.l$internal$allClassificationThresholds = c(0.1, 0.05, 0.01, 0.001)
  
  # Minimum number of TFBS to include a TF in the heatmap
  par.l$internal$threshold_minNoTFBS_heatmap = 100
  
  # Colors for the different classifications
  par.l$internal$colorCategories = c("activator" = "#4daf4a", "undetermined" = "black", "repressor" = "#e41a1c", "not-expressed" = "Snow3") # diverging, modified
  
  
  GRN@config$parameters = par.l
  GRN@config$metadata = objectMetadata
  
  
  # OUTPUT
  GRN@config$directories$outputRoot         = outputFolder
  GRN@config$directories$output_plots       = dir_output_plots 
  GRN@config$files$output_log               = paste0(outputFolder, "GRN.log")
  
  .testExistanceAndCreateDirectoriesRecursively(c(outputFolder, dir_output_plots))
  
  checkmate::assertDirectory(outputFolder, access = "w")
  
  
  .startLogger(GRN@config$files$output_log , "INFO",  removeOldLog = FALSE)
  #.printParametersLog(par.l)
  
  futile.logger::flog.info(paste0("Empty GRN object created successfully. Type the object name (e.g., GRN) to retrieve summary information about it at any time."))
  
  futile.logger::flog.info(paste0(" Default output folder: ", GRN@config$directories$outputRoot))
  futile.logger::flog.info(paste0(" Genome assembly: ", genomeAssembly))
  
  .printExecutionTime(start, prefix = "")
  
  GRN
}

######## Add and filter data ########

#' Add data to a \code{\linkS4class{GRN}} object.
#' 
#' This function adds both RNA and peak data to a \code{\linkS4class{GRN}} object, along with data normalization.
#' In addition, and highly recommended, sample metadata can be optionally provided.
#' 
#' If the \code{ChIPseeker} package is installed, additional peak annotation is provided in the annotation slot and a peak annotation QC plot is produced as part of peak-gene QC.
#' This is fully optional, however, and has no consequences for downstream functions.
#' 
#' @export
#' @template GRN 
#' @param counts_peaks Data frame. No default. Counts for the peaks, with raw or normalized counts for each peak (rows) across all samples (columns). 
#' In addition to the count data, it must also contain one ID column with a particular format, see the argument \code{idColumn_peaks} below. 
#' Row names are ignored, column names must be set to the sample names and must match those from the RNA counts and the sample metadata table.
#' @param normalization_peaks Character. Default \code{DESeq_sizeFactor}. 
#' Normalization procedure for peak data. Must be one of \code{DESeq_sizeFactor}, \code{none}, or \code{quantile}.
#' @param idColumn_peaks Character. Default \code{peakID}. Name of the column in the counts_peaks data frame that contains peak IDs. 
#' The required format must be {chr}:{start}-{end}", with {chr} denoting the abbreviated chromosome name, and {start} and {end} the begin and end 
#' of the peak coordinates, respectively. End must be bigger than start. Examples for valid peak IDs are \code{chr1:400-800} or \code{chrX:20-25}.
#' @param counts_rna Data frame. No default. Counts for the RNA-seq data, with raw or normalized counts for each gene (rows) across all samples (columns). 
#' In addition to the count data, it must also contain one ID column with a particular format, see the argument \code{idColumn_rna} below. 
#' Row names are ignored, column names must be set to the sample names and must match those from the RNA counts and the sample metadata table.
#' @param normalization_rna Character. Default \code{quantile}. Normalization procedure for peak data. 
#' Must be one of "DESeq_sizeFactor", "none", or "quantile". For "quantile", \code{limma::normalizeQuantiles} is used for normalization.
#' @param idColumn_RNA Character. Default \code{ENSEMBL}. Name of the column in the \code{counts_rna} data frame that contains Ensembl IDs.
#' @param sampleMetadata Data frame. Default \code{NULL}. Optional, additional metadata for the samples, such as age, sex, gender etc. 
#' If provided, the @seealso [plotPCA_all()] function can then incorporate and plot it. Sample names must match with those from both peak and RNA-Seq data. The first column is expected to contain the sample IDs, the actual column name is irrelevant.
#' @param allowOverlappingPeaks \code{TRUE} or \code{FALSE}. Default \code{FALSE}. Should overlapping peaks be allowed (then only a warning is issued 
#' when overlapping peaks are found) or (the default) should an error be raised?
#' @template forceRerun
#' @return An updated \code{\linkS4class{GRN}} object, with added data from this function (e.g., slots \code{GRN@data$peaks} and \code{GRN@data$RNA})
#' @seealso \code{\link{plotPCA_all}}
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' # library(readr)
#' # rna.df   = read_tsv("https://www.embl.de/download/zaugg/GRaNIE/rna.tsv.gz")
#' # peaks.df = read_tsv("https://www.embl.de/download/zaugg/GRaNIE/peaks.tsv.gz")
#' # meta.df  = read_tsv("https://www.embl.de/download/zaugg/GRaNIE/sampleMetadata.tsv.gz")
#' # GRN = loadExampleObject()
#' # We omit sampleMetadata = meta.df in the following line, becomes too long otherwise
#' # GRN = addData(GRN, counts_peaks = peaks.df, counts_rna = rna.df, forceRerun = FALSE)

# TODO: add a isSingleCell argument or something similar

addData <- function(GRN, counts_peaks, normalization_peaks = "DESeq_sizeFactor", idColumn_peaks = "peakID", 
                    counts_rna, normalization_rna = "quantile", idColumn_RNA = "ENSEMBL", sampleMetadata = NULL,
                    allowOverlappingPeaks= FALSE,
                    forceRerun = FALSE) {
  
  start = Sys.time()
  
  checkmate::assertClass(GRN, "GRN")
  GRN = .addFunctionLogToObject(GRN)    

  checkmate::assertDataFrame(counts_peaks, min.rows = 1, min.cols = 2)
  checkmate::assertDataFrame(counts_rna, min.rows = 1, min.cols = 2)
  checkmate::assertCharacter(idColumn_peaks, min.chars = 1, len = 1)
  checkmate::assertCharacter(idColumn_RNA, min.chars = 1, len = 1)
  checkmate::assertChoice(normalization_peaks, c("none", "DESeq_sizeFactor", "quantile"))
  checkmate::assertChoice(normalization_rna, c("none", "DESeq_sizeFactor", "quantile"))  
  checkmate::assertFlag(allowOverlappingPeaks)
  checkmate::assertFlag(forceRerun)
  
  if (is.null(GRN@data$peaks$counts) |
      is.null(GRN@data$peaks$counts_metadata) | 
      is.null(GRN@data$RNA$counts) |
      is.null(GRN@data$RNA$counts_metadata) |
      forceRerun) {
  
    # Normalize ID column names
    if (idColumn_peaks != "peakID") {
      counts_peaks = dplyr::rename(counts_peaks, peakID = !!(idColumn_peaks))
      idColumn_peaks ="peakID"
    }
    if (idColumn_RNA != "ENSEMBL") {
      counts_rna = dplyr::rename(counts_rna, ENSEMBL = !!(idColumn_RNA))
      idColumn_RNA = "ENSEMBL"
    }
      
    # Check existence of correct ID column now
    checkmate::assertSubset(idColumn_peaks, colnames(counts_peaks))
    checkmate::assertSubset(idColumn_RNA, colnames(counts_rna))

    # Check ID columns for missing values and remove
    rna_missing_ID =  which(is.na(counts_rna$ENSEMBL))
    if (length(rna_missing_ID) > 0) {
      message = paste0(" Found ", length(rna_missing_ID), " missing IDs in the ID column of the RNA counts. These rows will be removed.")
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
      counts_rna = dplyr::slice(counts_rna, -rna_missing_ID)
    }
    
    peaks_missing_ID =  which(is.na(counts_peaks$peakID))
    if (length(peaks_missing_ID) > 0) {
      message = paste0(" Found ", length(peaks_missing_ID), " missing IDs in the ID column of the peaks counts. These rows will be removed.")
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
      counts_peaks = dplyr::slice(counts_peaks, -peaks_missing_ID)
    }
    
    
    # Remove potential scientific notation from peak IDs
    peaks_eNotation = which(grepl("e+", counts_peaks$peakID))
    if (length(peaks_eNotation) > 0) {
      message = paste0("Found at least one peak (", paste0(counts_peaks$peakID[peaks_eNotation], collapse = ",") , ") for which the position contains the scientific notation, attempting to fix.")
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
      counts_peaks$peakID[peaks_eNotation] = .removeScientificNotation_positions(counts_peaks$peakID[peaks_eNotation])
      
      
    }
    
    # Clean Ensembl IDs
    counts_rna$ENSEMBL = gsub("\\..+", "", counts_rna$ENSEMBL, perl = TRUE)
    
    # Check uniqueness of IDs
    nDuplicateRows = nrow(counts_rna) - length(unique(counts_rna$ENSEMBL))
    if (nDuplicateRows > 0) {
      message = paste0(" Found ", nDuplicateRows, " duplicate rows in RNA-Seq data, consolidating them by summing them up.")
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
        
      counts_rna = counts_rna %>%
        dplyr::group_by(.data$ENSEMBL) %>%
        dplyr::summarise_if(is.numeric, sum) 
      # dplyr::summarise_if(is.numeric, sum, .groups = 'drop') # the .drop caused an error with dplyr 1.0.5
    }
    
    
    # Normalize counts
    countsPeaks.norm.df  = .normalizeCounts(counts_peaks, method = normalization_peaks, idColumn = "peakID")
    countsRNA.norm.df  = .normalizeCounts(counts_rna, method = normalization_rna, idColumn = "ENSEMBL")
    
    GRN@config$parameters$normalization_rna = normalization_rna
    GRN@config$parameters$normalization_peaks = normalization_peaks
    
    # We have our first connection type, the default one; more may be added later
    GRN@config$TF_peak_connectionTypes = "expression"
    
    # Make sure ENSEMBL is the first column
    countsRNA.norm.df = dplyr::select(countsRNA.norm.df, .data$ENSEMBL, tidyselect::everything())
    countsPeaks.norm.df = dplyr::select(countsPeaks.norm.df, .data$peakID, tidyselect::everything())
    
    samples_rna  = colnames(countsRNA.norm.df)
    samples_peaks =  colnames(countsPeaks.norm.df)
    allSamples =  unique(c(samples_rna, samples_peaks)) %>% setdiff(c("ENSEMBL", "isFiltered", "peakID"))
    
    # Subset data to retain only samples that appear in both RNA and peaks
    data.l = .intersectData(countsRNA.norm.df, countsPeaks.norm.df)
    
    # Generate metadata first to determine the nmberof shared samples etc
    if (!is.null(sampleMetadata)) {
      
      futile.logger::flog.info("Parsing provided metadata...")
      GRN@data$metadata = sampleMetadata %>% dplyr::distinct() %>% tibble::tibble(.name_repair = "universal")
      
      # Force the first column to be the ID column
      if ("sampleID" %in% colnames(GRN@data$metadata)) {
        message = paste0("Renaming the first column to sampleID. However, this column already exists, it will be renamed accordingly.")
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)  
          
        colnames(GRN@data$metadata)[which(colnames(GRN@data$metadata) == "sampleID")] = "sampleID_original"
        
      } 
      colnames(GRN@data$metadata)[1] = "sampleID"
      
      # Assume the ID is in column 1, has to be unique
      if (nrow(GRN@data$metadata) > length(unique(GRN@data$metadata$sampleID))) {
        message = paste0("The first column in the sample metadata table must contain only unique values, as it is used as sample ID. Make sure the values are unique.")
        tbl_ids = table(GRN@data$metadata$sampleID)
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
      }
      
      missingIDs = which(! allSamples %in% GRN@data$metadata$sampleID)
      if (length(missingIDs) > 0) {
        GRN@data$metadata = tibble::add_row(GRN@data$metadata, sampleID = allSamples[ missingIDs])
      }
    } else {
      GRN@data$metadata = tibble::tibble(sampleID = allSamples)
    }
    
    GRN@data$metadata =  GRN@data$metadata %>%
      dplyr::mutate(has_RNA = .data$sampleID  %in% samples_rna,
                    has_peaks = .data$sampleID %in% samples_peaks,
                    has_both = .data$has_RNA & .data$has_peaks
      )
    
    GRN@config$sharedSamples = dplyr::filter(GRN@data$metadata, .data$has_both) %>% dplyr::pull(.data$sampleID) %>% as.character()
    
    ### COUNT MATRICES
    # Store the matrices either as normal or sparse matrix
    
    GRN@data$peaks$counts = .storeAsMatrixOrSparseMatrix(GRN, df = data.l[["peaks"]], ID_column = "peakID", slotName = "GRN@data$peaks$counts")

    # includeFiltered = TRUE here as it doesnt make a difference and because getCounts requires counts_metadata$isFiltered to be set already
    GRN@data$peaks$counts_metadata = .createConsensusPeaksDF(getCounts(GRN, type = "peaks", permuted = FALSE, asMatrix = FALSE, includeFiltered = TRUE)) 
    stopifnot(c("chr", "start", "end", "peakID", "isFiltered") %in% colnames(GRN@data$peaks$counts_metadata))
    
    GRN@data$RNA$counts   = .storeAsMatrixOrSparseMatrix(GRN, df =  data.l[["RNA"]], ID_column = "ENSEMBL", slotName = "GRN@data$RNA$counts")
      
    GRN@data$RNA$counts_metadata = tibble::tibble(ID = data.l[["RNA"]]$ENSEMBL, isFiltered = FALSE)
    
    GRN@data$RNA$counts_permuted_index = .shuffleColumns(countsRNA.norm.df %>% dplyr::select(-.data$ENSEMBL), .getMaxPermutation(GRN), returnUnshuffled = FALSE, returnAsList = FALSE)
    
    futile.logger::flog.info(paste0( "Final dimensions of data:"))
    futile.logger::flog.info(paste0( " RNA  : ", nrow(countsRNA.norm.df)  , " x ", ncol(countsRNA.norm.df)   - 2, " (rows x columns)"))
    futile.logger::flog.info(paste0( " peaks: ", nrow(countsPeaks.norm.df), " x ", ncol(countsPeaks.norm.df) - 1, " (rows x columns)"))
    # Create permutations for RNA
    futile.logger::flog.info(paste0( "Generate ", .getMaxPermutation(GRN), " permutations of RNA-counts"))
    
    
    futile.logger::flog.info(paste0("Creating consensus peaks and check for overlapping peaks..."))
    
    # Consensus peaks

    # Assume 0-based exclusive format, see https://arnaudceol.wordpress.com/2014/09/18/chromosome-coordinate-systems-0-based-1-based/ and http://genome.ucsc.edu/FAQ/FAQformat.html#format1 for details
    consensus.gr   = .constructGRanges(GRN@data$peaks$counts_metadata, seqlengths = .getChrLengths(GRN@config$parameters$genomeAssembly), GRN@config$parameters$genomeAssembly, zeroBased = TRUE)

    overlappingPeaks = which(GenomicRanges::countOverlaps(consensus.gr ,consensus.gr) >1)
    
    if (length(overlappingPeaks) > 0) {
      
      ids = (consensus.gr[overlappingPeaks] %>% as.data.frame())$peakID
      
      messageAll = paste0(" ", length(overlappingPeaks), 
                          " overlapping peaks have been identified. The first ten are: ", paste0(ids[seq_len(min(10, length(ids)))], collapse = ","),
                          ". This may not be what you want, since overlapping peaks may have a heigher weight in the network. "
      )
      
      
      if (allowOverlappingPeaks) {
        
        message = paste0(messageAll, "As allowOverlappingPeaks has been set to TRUE, this is only a warning and not an error.")
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
      } else {
        message = paste0(messageAll, "As allowOverlappingPeaks = FALSE (the default), this is an error and not a warning. You may want to regenerate the peak file, eliminate peak overlaps, and rerun this function")
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
      }
      
    }
    
  }
  futile.logger::flog.info(paste0("Adding peak and gene annotation..."))
  
  # Add peak annotation once
  GRN = .populatePeakAnnotation(GRN)
  
  # Add gene annotation once
  GRN = .populateGeneAnnotation(GRN)
  
  .printExecutionTime(start, prefix = "")
  
  GRN
  
}

.storeAsMatrixOrSparseMatrix <- function (GRN, df, ID_column, slotName, threshold = 0.1) {
    
    stopifnot(identical(GRN@config$sharedSamples, colnames(df)[-1]))
    
    # Store as sparse matrix if enough 0s
    checkmate::assertIntegerish(length(GRN@config$sharedSamples), lower = 1)

    df.m = df %>% 
        dplyr::select(tidyselect::one_of(ID_column, GRN@config$sharedSamples)) %>% 
        tibble::column_to_rownames(ID_column) %>%
        as.matrix()
    
    # Determine sparsity
    fractionZero = (length(df.m) - Matrix::nnzero(df.m)) / length(df.m)
    
    
    if (fractionZero > threshold) {
        futile.logger::flog.info(paste0("Storing ", slotName, " matrix as sparse matrix because fraction of 0s is > ", threshold, " (", fractionZero, ")"))
        df.m = .asSparseMatrix(df.m, convertNA_to_zero = FALSE, 
                               dimnames = list(df[, ID_column] %>% unlist(use.names = FALSE), GRN@config$sharedSamples))
    } 
    
    df.m
    
}

.createConsensusPeaksDF <- function(countsPeaks, idColumn = "peakID") {
  
  checkmate::assertChoice(idColumn, colnames(countsPeaks))
  
  ids.split = strsplit(countsPeaks %>% dplyr::pull(!!(idColumn)), split = "[:-]+")
  ids.split.length = sapply(ids.split, length)
  if (!all(ids.split.length == 3)) {
    message = paste0(" At least one of the IDs in the peaks data has an unsupported format. Make sure all peakIDs are in the format \"chr:start-end\"")
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
  
  consensus.df = tibble::tibble(chr   = as.factor(sapply(ids.split, "[[", 1)),
                                start = as.numeric(sapply(ids.split, "[[", 2)), 
                                end   = as.numeric(sapply(ids.split, "[[", 3)),
                                peakID = paste0(.data$chr, ":", .data$start, "-", .data$end),
                                isFiltered = FALSE)
  consensus.df
}

.removeScientificNotation_positions <- function(peakIDs.vec) {
  ids = strsplit(peakIDs.vec, split = ":", fixed = TRUE)
  ids_chr = sapply(ids, "[[", 1)
  ids_pos = sapply(ids, "[[", 2)
  ids_pos = strsplit(ids_pos, split = "-", fixed = TRUE)
  start = sapply(ids_pos, "[[", 1)
  end   = sapply(ids_pos, "[[", 2)
  
  paste0(ids_chr, ":", format(as.integer(start), scientific = FALSE), "-", format(as.integer(end), scientific= FALSE))
}

#' @importFrom biomaRt useEnsembl getBM
.retrieveAnnotationData <- function(genomeAssembly) {
    
    futile.logger::flog.info(paste0("Retrieving genome annotation data from biomaRt for ", genomeAssembly, "... This may take a while"))
    
    params.l = .getBiomartParameters(genomeAssembly)
    
    columnsToRetrieve = c("chromosome_name", "start_position", "end_position",
                          "strand", "ensembl_gene_id", "gene_biotype", "external_gene_name")
    
    geneAnnotation <- NULL
    attempt <- 0
    mirrors = c('www', 'uswest', 'useast', 'asia')
    while(!is.data.frame(geneAnnotation) && attempt <= 3 ) {
        attempt <- attempt + 1
        geneAnnotation = tryCatch({ 
            ensembl = biomaRt::useEnsembl(biomart = "genes", host = params.l[["host"]], dataset = params.l[["dataset"]], mirror = mirrors[attempt])
            biomaRt::getBM(attributes = columnsToRetrieve, mart = ensembl)
            
        }
        )
    } 
    
    
    if (!is.data.frame(geneAnnotation)) {
        
        error_Biomart = "A temporary error occured with biomaRt::getBM or biomaRt::useEnsembl. This is often caused by an unresponsive Ensembl site. Try again at a later time. Note that this error is not caused by GRaNIE but external services."
        .checkAndLogWarningsAndErrors(NULL, error_Biomart, isWarning = FALSE)
        return(NULL)
        
    }
    
    geneAnnotation %>%
        tibble::as_tibble() %>%
        dplyr::filter(stringr::str_length(.data$chromosome_name) <= 5) %>%
        dplyr::mutate(chromosome_name = paste0("chr", .data$chromosome_name)) %>%
        dplyr::rename(gene.chr = .data$chromosome_name, gene.start = .data$start_position, gene.end = .data$end_position, 
                      gene.strand = .data$strand, gene.ENSEMBL = .data$ensembl_gene_id, gene.type = .data$gene_biotype, gene.name = .data$external_gene_name) %>%
        tidyr::replace_na(list(gene.type = "unknown")) %>%
        dplyr::mutate_if(is.character, as.factor) %>%
        dplyr::mutate(gene.type = dplyr::recode_factor(.data$gene.type, lncRNA = "lincRNA")) %>%  # there seems to be a name change from lincRNA -> lncRNA, lets change it here 
        dplyr::mutate(gene.strand = factor(.data$gene.strand, levels = c(1,-1), labels = c("+", "-")))
    
}


.populatePeakAnnotation <- function (GRN) {
  
  countsPeaks.clean = getCounts(GRN, type = "peaks", permuted = FALSE, asMatrix = TRUE, includeFiltered = TRUE)
  
  futile.logger::flog.info(paste0(" Calculate statistics for each peak (mean and CV)"))
  
  rowMeans_peaks   = rowMeans(countsPeaks.clean)
  rowMedians_peaks = matrixStats::rowMedians(countsPeaks.clean)
  CV_peaks = matrixStats::rowSds(countsPeaks.clean) /  rowMeans_peaks
  
  metadata_peaks = tibble::tibble(peak.ID = rownames(countsPeaks.clean), 
                                  peak.mean = rowMeans_peaks, 
                                  peak.median = rowMedians_peaks, 
                                  peak.CV = CV_peaks)
  
  GRN@annotation$peaks = metadata_peaks

  if (!is.installed("ChIPseeker")) {
      packageMessage = paste0("The package ChIPseeker is currently not installed, which is needed for additional peak annotation that can be useful for further downstream analyses. ", 
                              " You may want to install it and re-run this function. However, this is optional and except for some missing additional annotation columns, there is no limitation.")
      .checkPackageInstallation("ChIPseeker", packageMessage, isWarning = TRUE)
  } else {
    
    futile.logger::flog.info(paste0(" Retrieve peak annotation using ChipSeeker. This may take a while"))
    genomeAssembly = GRN@config$parameters$genomeAssembly
    consensusPeaks     = GRN@data$peaks$counts_metadata %>% dplyr::filter(!.data$isFiltered)
    consensusPeaks.gr  = .constructGRanges(GRN@data$peaks$counts_metadata, seqlengths = .getChrLengths(genomeAssembly), genomeAssembly)
    
    # Add ChIPSeeker anotation
    peaks.annotated = suppressMessages(ChIPseeker::annotatePeak(
      consensusPeaks.gr,
      tssRegion = c(-5000, 5000), # extended from 3kb to 5
      TxDb = .getGenomeObject(genomeAssembly, type = "txbd"),
      level = "gene", 
      assignGenomicAnnotation = TRUE,  # the default
      genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                    "Downstream", "Intergenic"),  # the default
      annoDb = .getGenomeObject(genomeAssembly, type = "packageName"), # optional, if provided, extra columns including SYMBOL, GENENAME, ENSEMBL/ENTREZID will be added
      sameStrand = FALSE, # the default
      ignoreOverlap = FALSE, # the default
      ignoreUpstream = FALSE, # the default
      ignoreDownstream = FALSE, # the default
      overlap = "TSS", # the default
      verbose = TRUE # the default
    ))
    
    GRN@annotation$peaks_obj = peaks.annotated
    
    peaks.annotated.df = as.data.frame(peaks.annotated)
    peaks.annotated.df$annotation[grepl("Exon", peaks.annotated.df$annotation)] = "Exon"
    peaks.annotated.df$annotation[grepl("Intron", peaks.annotated.df$annotation)] = "Intron"
    
    GRN@annotation$peaks = dplyr::left_join(GRN@annotation$peaks, 
                                                     peaks.annotated.df  %>% 
                                                       dplyr::select(.data$peakID, .data$annotation, tidyselect::starts_with("gene"), -.data$geneId, .data$distanceToTSS, .data$ENSEMBL, .data$SYMBOL, .data$GENENAME) %>%
                                                       dplyr::mutate(annotation  = as.factor(.data$annotation), 
                                                                     ENSEMBL = as.factor(.data$ENSEMBL), 
                                                                     GENENAME = as.factor(.data$GENENAME),
                                                                     SYMBOL = as.factor(.data$SYMBOL)),
                                                     by = c("peak.ID" = "peakID")) %>%
      dplyr::rename(peak.nearestGene.chr = .data$geneChr,
                    peak.nearestGene.start = .data$geneStart, 
                    peak.nearestGene.end = .data$geneEnd, 
                    peak.nearestGene.length = .data$geneLength, 
                    peak.nearestGene.strand = .data$geneStrand, 
                    peak.nearestGene.name = .data$GENENAME,
                    peak.nearestGene.distanceToTSS = .data$distanceToTSS,
                    peak.nearestGene.ENSEMBL = .data$ENSEMBL,
                    peak.nearestGene.symbol = .data$SYMBOL,
                    peak.annotation = .data$annotation
      )
    
    
  }
  
  
  
  
  # Also add GC content as annotation columns
  GRN = .calcGCContentPeaks(GRN)
  
  GRN
  
}

.populateGeneAnnotation <- function (GRN) {


  countsRNA.m  = getCounts(GRN, type = "rna", permuted = FALSE, asMatrix = TRUE, includeFiltered = TRUE)
  
  futile.logger::flog.info(paste0(" Calculate statistics for each of the ", nrow(countsRNA.m), " genes that were provided with the RNA-seq data (mean and CV)"))
  
  
  rowMeans_rna = rowMeans(countsRNA.m)
  rowMedians_rna = matrixStats::rowMedians(countsRNA.m)
  CV_rna = matrixStats::rowSds(countsRNA.m) /  rowMeans_rna
  
  genomeAnnotation.df = .retrieveAnnotationData(GRN@config$parameters$genomeAssembly)
  
  metadata_rna = tibble::tibble(gene.ENSEMBL = rownames(countsRNA.m), 
                                gene.mean = rowMeans_rna, 
                                gene.median = rowMedians_rna, 
                                gene.CV = CV_rna) %>%
    dplyr::left_join(genomeAnnotation.df, by = c("gene.ENSEMBL")) %>%
    dplyr::mutate(gene.type = forcats::fct_explicit_na(.data$gene.type, na_level = "unknown/missing"))
  
  GRN@annotation$genes = metadata_rna
  
  GRN
  
}

.populateGOAnnotation <- function(GRN, results.tbl, ontology){
  
  GRN@annotation$GO[[ontology]] = results.tbl[,c("ID", "Term")]
  GRN
  
}

#' @importFrom rlang .data `:=`
.calcGCContentPeaks <- function(GRN) {
  
  futile.logger::flog.info(paste0("Calculate GC-content for peaks. This may take a while"))
  start = Sys.time()
  genomeAssembly = GRN@config$parameters$genomeAssembly

  genome = .getGenomeObject(genomeAssembly, type = "BSgenome")
  
  # Get peaks as GRanges object
  query   = .constructGRanges(GRN@data$peaks$counts_metadata, 
                              seqlengths = .getChrLengths(genomeAssembly), 
                              genomeAssembly)
  
  # Get DNAStringSet object
  seqs_peaks = Biostrings::getSeq(genome, query)
  
  GC_content.df = Biostrings::letterFrequency(seqs_peaks, "GC") / Biostrings::letterFrequency(seqs_peaks, "ACGT")
  
  GC_content.df = GC_content.df %>%
    tibble::as_tibble() %>%
    dplyr::mutate(length = GenomicRanges::width(query),
                  peak.ID = query$peakID,
                  GC_class = cut(.data$`G|C`, breaks = seq(0,1,0.1), include.lowest = TRUE, ordered_result = TRUE))
  
  GC_classes.df = GC_content.df %>%
    dplyr::group_by(.data$GC_class) %>%
    dplyr::summarise(n= dplyr::n(), peak_width_mean = mean(length), peak_width_sd = sd(length)) %>%
    dplyr::ungroup() %>% 
    tidyr::complete(.data$GC_class, fill = list(n = 0)) %>%
    dplyr::mutate(n_rel = .data$n / length(query))
  
  # TODO: Put where
  #ggplot2::ggplot(GC_content.df, ggplot2::aes(GC.class)) + geom_histogram(stat = "count") + ggplot2::theme_bw()
  
  #ggplot2::ggplot(GC_classes.df , ggplot2::aes(GC.class, n_rel)) + geom_bar(stat = "identity") + ggplot2::theme_bw()
  
  GRN@annotation$peaks = dplyr::left_join(GRN@annotation$peaks, GC_content.df, by = "peak.ID") %>%
    dplyr::rename( peak.GC.perc    = .data$`G|C`,
                   peak.width      = .data$length,
                   peak.GC.class   = .data$GC_class)
  
  GRN@stats$peaks = list()
  GRN@stats$peaks$GC = GC_classes.df
  
  .printExecutionTime(start)
  
  GRN
}

#' Filter RNA-seq and/or peak data from a \code{\linkS4class{GRN}} object
#' 
#' This function marks genes and/or peaks as \code{filtered} depending on the chosen filtering criteria. Filtered genes / peaks will then be 
#' disregarded when adding connections in subsequent steps via \code{\link{addConnections_TF_peak}} and  \code{\link{addConnections_peak_gene}}. \strong{This function does NOT (re)filter existing connections when the \code{\linkS4class{GRN}} object already contains connections. Thus, upon re-execution of this function with different filtering criteria, all downstream steps have to be re-run.}
#' 
#' All this function does is setting (or modifying) the filtering flag in \code{GRN@data$peaks$counts_metadata} and \code{GRN@data$RNA$counts_metadata}, respectively.
#' 
#' @template GRN 
#' @param minNormalizedMean_peaks Numeric[0,] or \code{NULL}. Default 5. Minimum mean across all samples for a peak to be retained for the normalized counts table. Set to \code{NULL} for not applying the filter.
#' @param maxNormalizedMean_peaks Numeric[0,] or \code{NULL}. Default \code{NULL}. Maximum mean across all samples for a peak to be retained for the normalized counts table. Set to \code{NULL} for not applying the filter.
#' @param minNormalizedMeanRNA Numeric[0,] or \code{NULL}. Default 5. Minimum mean across all samples for a gene to be retained for the normalized counts table. Set to \code{NULL} for not applying the filter.
#' @param maxNormalizedMeanRNA Numeric[0,] or \code{NULL}. Default \code{NULL}. Maximum mean across all samples for a gene to be retained for the normalized counts table. Set to \code{NULL} for not applying the filter.
#' @param chrToKeep_peaks Character vector or \code{NULL}. Default \code{NULL}. Vector of chromosomes that peaks are allowed to come from. This filter can be used to filter sex chromosomes from the peaks, for example (e.g, \code{c(paste0("chr", 1:22), "chrX", "chrY")})
#' @param minSize_peaks Integer[1,] or \code{NULL}. Default \code{NULL}. Minimum peak size (width, end - start) for a peak to be retained. Set to \code{NULL} for not applying the filter.
#' @param maxSize_peaks Integer[1,] or \code{NULL}. Default 10000. Maximum peak size (width, end - start) for a peak to be retained. Set to \code{NULL} for not applying the filter.
#' @param minCV_peaks Numeric[0,] or \code{NULL}. Default \code{NULL}. Minimum CV (coefficient of variation, a unitless measure of variation) for a peak to be retained. Set to \code{NULL} for not applying the filter.
#' @param maxCV_peaks Numeric[0,] or \code{NULL}. Default \code{NULL}. Maximum CV (coefficient of variation, a unitless measure of variation) for a peak to be retained. Set to \code{NULL} for not applying the filter.
#' @param minCV_genes Numeric[0,] or \code{NULL}. Default \code{NULL}. Minimum CV (coefficient of variation, a unitless measure of variation) for a gene to be retained. Set to \code{NULL} for not applying the filter.
#' @param maxCV_genes Numeric[0,] or \code{NULL}. Default \code{NULL}. Maximum CV (coefficient of variation, a unitless measure of variation) for a gene to be retained. Set to \code{NULL} for not applying the filter.
#' @template forceRerun
#' @return An updated \code{\linkS4class{GRN}} object, with added data from this function.
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' GRN = filterData(GRN, forceRerun = FALSE)
#' @export
filterData <- function (GRN, 
                        minNormalizedMean_peaks = 5, maxNormalizedMean_peaks = NULL, 
                        minNormalizedMeanRNA = 1,  maxNormalizedMeanRNA = NULL,
                        chrToKeep_peaks = NULL,
                        minSize_peaks = NULL, maxSize_peaks = 10000,
                        minCV_peaks = NULL, maxCV_peaks = NULL,
                        minCV_genes = NULL, maxCV_genes = NULL,
                        forceRerun = FALSE) {
  
  start = Sys.time()
    
  checkmate::assertClass(GRN, "GRN")
  GRN = .addFunctionLogToObject(GRN) 
  
  GRN = .makeObjectCompatible(GRN)

  checkmate::assertNumber(minNormalizedMean_peaks, lower = 0, null.ok = TRUE)
  checkmate::assertNumber(minNormalizedMeanRNA, lower = 0, null.ok = TRUE)
  checkmate::assertNumber(maxNormalizedMean_peaks, lower = minNormalizedMean_peaks , null.ok = TRUE)
  checkmate::assertNumber(maxNormalizedMeanRNA, lower = minNormalizedMeanRNA, null.ok = TRUE)
  checkmate::assertCharacter(chrToKeep_peaks, min.len = 1, any.missing = FALSE, null.ok = TRUE)
  checkmate::assertSubset(chrToKeep_peaks, GRN@data$peaks$counts_metadata %>% dplyr::pull(.data$chr) %>% unique() %>% as.character())
  
  checkmate::assertIntegerish(minSize_peaks, lower = 1, null.ok = TRUE)
  checkmate::assertIntegerish(maxSize_peaks, lower = dplyr::if_else(is.null(minSize_peaks), 1, minSize_peaks), null.ok = TRUE)
  checkmate::assertNumber(minCV_peaks, lower = 0, null.ok = TRUE)
  checkmate::assertNumber(maxCV_peaks, lower = dplyr::if_else(is.null(minCV_peaks), 0, minCV_peaks), null.ok = TRUE)
  checkmate::assertNumber(minCV_genes, lower = 0, null.ok = TRUE)
  checkmate::assertNumber(maxCV_genes, lower = dplyr::if_else(is.null(minCV_genes), 0, minCV_genes), null.ok = TRUE)
  checkmate::assertFlag(forceRerun)
  
  if (!GRN@config$isFiltered | forceRerun) {
    
    GRN@data$peaks$counts_metadata$isFiltered = FALSE
    
    if(!is.null(GRN@data$TFs$TF_peak_overlap)) {
      GRN@data$TFs$TF_peak_overlap[, "isFiltered"] = 0
    }
    
    
    # Filter peaks
    futile.logger::flog.info("FILTER PEAKS")
    peakIDs.CV = .filterPeaksByMeanCV(GRN, 
                                      minMean = minNormalizedMean_peaks, maxMean = maxNormalizedMean_peaks, 
                                      minCV = minCV_peaks, maxCV = maxCV_peaks) 
    
    # Clean peaks from alternative contigs etc 
    GRN@config$parameters$chrToKeep =  chrToKeep_peaks
    peakIDs.chr = .filterPeaksByChromosomeAndSize(GRN, 
                                                  chrToKeep_peaks, 
                                                  minSize_peaks = minSize_peaks, maxSize_peaks = maxSize_peaks)
    
    nPeaksBefore = nrow(GRN@data$peaks$counts_metadata)
    peaks_toKeep = intersect(peakIDs.chr, peakIDs.CV)
    futile.logger::flog.info(paste0("Collectively, filter ", nPeaksBefore -length(peaks_toKeep), " out of ", nPeaksBefore, " peaks."))
    futile.logger::flog.info(paste0("Number of remaining peaks: ", length(peaks_toKeep)))
    
    if (length(peaks_toKeep) < 1000) {
        message = paste0("Too few peaks (", length(peaks_toKeep), ") remain after filtering. At least 1000 peaks must remain. Adjust the filtering settings.")
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    GRN@data$peaks$counts_metadata$isFiltered  = ! GRN@data$peaks$counts_metadata$peakID  %in% peaks_toKeep
    #GRN@data$peaks$counts_raw$isFiltered = ! GRN@data$peaks$counts_raw$peakID  %in% peaks_toKeep
    GRN@data$peaks$counts_metadata$isFiltered = ! GRN@data$peaks$counts_metadata$peakID  %in% peaks_toKeep
    
    
    if(!is.null(GRN@data$TFs$TF_peak_overlap)) {
      GRN@data$TFs$TF_peak_overlap[, "isFiltered"] = as.integer (! rownames(GRN@data$TFs$TF_peak_overlap) %in% peaks_toKeep)
    }
    
    
    # Remove genes with small rowMeans
    #Only for real data, not for permuted (rowmeans is equal anyway)
    # Filter peaks
    futile.logger::flog.info("FILTER RNA-seq")
    genes.CV = .filterGenesByMeanCV(GRN, 
                                    minMean = minNormalizedMeanRNA, maxMean = maxNormalizedMeanRNA, 
                                    minCV = minCV_genes, maxCV = maxCV_genes) 
    
    
    if (length(genes.CV) < 100) {
        message = paste0("Too few genes (", length(genes.CV), ") remain after filtering. At least 100 genes must remain. Adjust the filtering settings.")
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    rowMeans = rowMeans(getCounts(GRN, type = "rna", asMatrix = TRUE, includeFiltered = TRUE))
    GRN@data$RNA$counts_metadata$isFiltered = rowMeans < minNormalizedMeanRNA 

    nRowsFlagged = length(which(GRN@data$RNA$counts_metadata$isFiltered))
    
    # Raw counts are left untouched and filtered where needed only
    futile.logger::flog.info(paste0(" Flagged ", nRowsFlagged, " rows because the row mean was smaller than ", minNormalizedMeanRNA))
    
    GRN@config$isFiltered = TRUE
  } 
  
  .printExecutionTime(start, prefix = "")
  
  GRN
}


.filterPeaksByChromosomeAndSize <- function(GRN, chrToKeep, minSize_peaks = NULL, maxSize_peaks, idColumn = "peakID") {
  
  startTime = Sys.time()
  
  if (is.null(minSize_peaks)) {
    minSize_peaks = 1
  }
  
  if (is.null(chrToKeep)) {
    chrToKeep = GRN@data$peaks$counts_metadata %>% dplyr::pull(.data$chr) %>% unique()
  } else {
    futile.logger::flog.info(paste0("Filter and sort peaks and remain only those on the following chromosomes: ", paste0(chrToKeep, collapse = ",")))
  }
  

  futile.logger::flog.info(paste0("Filter and sort peaks by size and remain only those smaller than : ", maxSize_peaks))
  futile.logger::flog.info(paste0(" Number of peaks before filtering: ", nrow(GRN@data$peaks$counts_metadata)))
  
  countsPeaks.clean = GRN@data$peaks$counts_metadata %>%
    dplyr::mutate(size = end-start) %>%
    dplyr::filter(.data$chr %in% chrToKeep, .data$size <= maxSize_peaks, .data$size >= minSize_peaks) %>%
    # arrange(chr, start) %>%
    dplyr::rename(peakID = !!(idColumn)) %>%
    dplyr::select(-.data$chr,-.data$start, -.data$end, -.data$size) %>%
    dplyr::select(.data$peakID, tidyselect::everything())
  
  futile.logger::flog.info(paste0(" Number of peaks after filtering : ", nrow(countsPeaks.clean)))
  
  .printExecutionTime(startTime)
  countsPeaks.clean$peakID
}

.filterPeaksByCV <- function (GRN, minCV = 0, maxCV = 1e7) {
  
  startTime = Sys.time()
  
  if (is.null(maxCV)) {
    futile.logger::flog.info(paste0("Filter peaks by CV: Min CV = ", minCV))
    
    peaksFiltered = dplyr::filter(GRN@annotation$peaks, .data$peak.CV >= minCV)
    
  } else {
    futile.logger::flog.info(paste0("Filter peaks by CV: Min CV = ", minCV, ", max CV = ", maxCV))
    
    peaksFiltered = dplyr::filter(GRN@annotation$peaks, .data$peak.CV >= minCV, .data$peak.CV <= maxCV)
  }
  
  futile.logger::flog.info(paste0(" Number of peaks after filtering : ", nrow(peaksFiltered)))
  
  .printExecutionTime(startTime)
  
  peaksFiltered$peak.ID
}

.filterPeaksByMeanCV <- function (GRN, minMean = 0, maxMean = NULL, minCV = 0, maxCV = NULL) {
  
  startTime = Sys.time()
  
  futile.logger::flog.info(paste0(" Number of peaks before filtering : ", nrow(GRN@annotation$peaks)))
  
  if (is.null(minCV)) {
    minCV = 0
  }
  
  if (is.null(maxCV)) {
    futile.logger::flog.info(paste0("  Filter peaks by CV: Min = ", minCV))
    maxCV = 9e+99
    
  } else {
    futile.logger::flog.info(paste0("  Filter peaks by CV: Min = ", minCV, ", Max = ", maxCV))
  }
  
  
  if (is.null(minMean)) {
    
    # As data can be pre-normalized, set the minimum to a very small value so the filter is effectively off
    minMean = -9e+99
  }
  
  if (is.null(maxMean)) {
    futile.logger::flog.info(paste0("  Filter peaks by mean: Min = ", minMean))
    maxMean = 9e+99
  } else {
    futile.logger::flog.info(paste0("  Filter peaks by mean: Min = ", minMean, ", Max = ", maxMean))  
  }   
  
  
  peaksFiltered = dplyr::filter(GRN@annotation$peaks, 
                                .data$peak.CV >= minCV, .data$peak.CV <= maxCV, 
                                .data$peak.mean >= minMean, .data$peak.mean <= maxMean)
  
  futile.logger::flog.info(paste0(" Number of peaks after filtering : ", nrow(peaksFiltered)))
  
  .printExecutionTime(startTime)
  
  peaksFiltered$peak.ID
}

.filterGenesByMeanCV <- function (GRN, minMean = 0, maxMean = NULL, minCV = 0, maxCV = NULL) {
  
  startTime = Sys.time()
  
  futile.logger::flog.info(paste0(" Number of genes before filtering : ", nrow(GRN@annotation$genes)))
  
  if (is.null(minCV)) {
    minCV = 0
  }
  
  if (is.null(maxCV)) {
    futile.logger::flog.info(paste0("  Filter genes by CV: Min = ", minCV))
    maxCV = 9e+99
    
  } else {
    futile.logger::flog.info(paste0("  Filter genes by CV: Min = ", minCV, ", Max = ", maxCV))
  }
  
  messageMean = paste0("  Filter genes by mean:")
  
  if (is.null(minMean)) {
    minMean = -9e+99
  } else {
    messageMean = paste0(messageMean, " Min = ", minMean)
  }
  

  if (is.null(maxMean)) {
    maxMean = 9e+99
  } else {
    messageMean = paste0(messageMean, " Max = ", maxMean)
  }   
  
  futile.logger::flog.info(messageMean)
  
  
  genesFiltered = dplyr::filter(GRN@annotation$genes, 
                                .data$gene.CV >= minCV, .data$gene.CV <= maxCV, 
                                .data$gene.mean >= minMean, .data$gene.mean <= maxMean)
  
  
  futile.logger::flog.info(paste0(" Number of genes after filtering : ", nrow(genesFiltered)))
  
  .printExecutionTime(startTime)
  
  genesFiltered$gene.ENSEMBL 
}



######## TFBS ########

#' Add TFBS to a \code{\linkS4class{GRN}} object
#' 
#' @template GRN 
#' @param motifFolder Character. No default. Path to the folder that contains the TFBS predictions. The files must be in BED format, 6 columns, one file per TF. See the other parameters for more details.
#' @param TFs Character vector. Default \code{all}. Vector of TF names to include. The special keyword \code{all} can be used to include all TF found in the folder as specified by \code{motifFolder}. If \code{all} is specified anywhere, all TFs will be included. TF names must otherwise match the file names that are found in the folder, without the file suffix.
#' @param filesTFBSPattern Character. Default \code{"_TFBS"}. Suffix for the file names in the TFBS folder that is not part of the TF name. Can be empty. For example, for the TF CTCF, if the file is called \code{CTCF.all.TFBS.bed}, set this parameter to \code{".all.TFBS"}.
#' @param fileEnding Character. Default \code{".bed"}. File ending for the files from the motif folder.
#' @param nTFMax \code{NULL} or Integer[1,]. Default \code{NULL}. Maximal number of TFs to import. Can be used for testing purposes, e.g., setting to 5 only imports 5 TFs even though the whole \code{motifFolder} has many more TFs defined.
#' @template forceRerun
#' @return An updated \code{\linkS4class{GRN}} object, with additional information added from this function (\code{GRN@annotation$TFs} in particular)
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' @export
addTFBS <- function(GRN, motifFolder, TFs = "all", nTFMax = NULL, filesTFBSPattern = "_TFBS", fileEnding = ".bed", forceRerun = FALSE) {
  
  start = Sys.time()
  checkmate::assertClass(GRN, "GRN")
  GRN = .addFunctionLogToObject(GRN)
  
  GRN = .makeObjectCompatible(GRN)

  checkmate::assertDirectoryExists(motifFolder)
  checkmate::assertCharacter(TFs, min.len = 1)
  checkmate::assert(checkmate::testNull(nTFMax), checkmate::testIntegerish(nTFMax, lower = 1))
  checkmate::assertCharacter(filesTFBSPattern, len = 1, min.chars = 0)
  checkmate::assertCharacter(fileEnding, len = 1, min.chars = 1)
  checkmate::assertFlag(forceRerun)
  
  if (is.null(GRN@annotation$TFs) | is.null(GRN@annotation$TFs) | is.null(GRN@config$allTF)  | is.null(GRN@config$directories$motifFolder) | forceRerun) {
    
    GRN@config$TFBS_fileEnding  = fileEnding
    GRN@config$TFBS_filePattern = filesTFBSPattern
    GRN@annotation$TFs = .getFinalListOfTFs(motifFolder, filesTFBSPattern, fileEnding, TFs, nTFMax, getCounts(GRN, type = "rna", permuted = FALSE))
    
    
    GRN@annotation$TFs = GRN@annotation$TFs %>%
      dplyr::rename(TF.ENSEMBL = .data$ENSEMBL, TF.HOCOID = .data$HOCOID)  %>% 
      dplyr::mutate(TF.name = .data$TF.HOCOID)  %>%
      dplyr::select(.data$TF.name, .data$TF.ENSEMBL, .data$TF.HOCOID)
    
    # TODO: Change here and make it more logical what to put where
    GRN@config$allTF = GRN@annotation$TFs$TF.name
    
    #Store all data-dependent TF information
    # GRN@config$TF_list = list()
    # GRN@config$TF_list[["all_TFBS"]] =GRN@config$allTF
    GRN@config$directories$motifFolder = motifFolder
    
    
  } 
  
  .printExecutionTime(start, prefix = "")
  
  GRN
  
}


.getFinalListOfTFs <- function(folder_input_TFBS, filesTFBSPattern, fileEnding, TFs, nTFMax, countsRNA) {
  
  futile.logger::flog.info(paste0("Checking database folder for matching files: ", folder_input_TFBS))
  files = .createFileList(folder_input_TFBS, "*.bed*", recursive = FALSE, ignoreCase = FALSE, verbose = FALSE)
  TFsWithTFBSPredictions = gsub(pattern = filesTFBSPattern, "", tools::file_path_sans_ext(basename(files), compression = TRUE))
  TFsWithTFBSPredictions = gsub(pattern = fileEnding, "", TFsWithTFBSPredictions)
  
  futile.logger::flog.info(paste0("Found ", length(TFsWithTFBSPredictions), " matching TFs: ", paste0(TFsWithTFBSPredictions, collapse = ", ")))
  
  
  # Filter TFs
  if (length(TFs) == 1 && TFs == "all") {
    
    futile.logger::flog.info(paste0("Use all TF from the database folder ", folder_input_TFBS))
    
  } else {
    
    futile.logger::flog.info(paste0("Subset TFs to user-specified list: ", paste0(TFs, collapse = ", ")))
    TFsWithTFBSPredictions = TFsWithTFBSPredictions[TFsWithTFBSPredictions %in% TFs]
    
    if (length(TFsWithTFBSPredictions) == 0) {
      message = paste0("No TFs are left after subsetting. Make sure the TF names are identical to the names in the database folder.")
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    futile.logger::flog.info(paste0("List of TFs: ", paste0(TFs, collapse = ", ")))
    
  }
  
  file_input_HOCOMOCO = paste0(folder_input_TFBS, .Platform$file.sep, "translationTable.csv")
  HOCOMOCO_mapping.df = .readHOCOMOCOTable(file_input_HOCOMOCO, delim = " ")
  
  TF_notExpressed = sort(dplyr::filter(HOCOMOCO_mapping.df, ! .data$ENSEMBL %in% countsRNA$ENSEMBL, .data$HOCOID %in% TFsWithTFBSPredictions) %>% dplyr::pull(.data$HOCOID))
  
  if (length(TF_notExpressed) > 0) {
    futile.logger::flog.info(paste0("Filtering the following ", length(TF_notExpressed), " TFs as they are not present in the RNA-Seq data: ", paste0(TF_notExpressed, collapse = ",")))
    
  }
  
  allTF = sort(dplyr::filter(HOCOMOCO_mapping.df, 
                             .data$ENSEMBL %in% countsRNA$ENSEMBL, 
                             .data$HOCOID %in% TFsWithTFBSPredictions) %>% dplyr::pull(.data$HOCOID))
  
  nTF = length(allTF)
  if (nTF == 0) {
    message = paste0("No shared Tfs.")
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
  if (!is.null(nTFMax)) {
    
    if (!is.null(nTFMax) && nTFMax < nTF) {
      futile.logger::flog.info(paste0("Use only the first ", nTFMax, " TFs because nTFMax has been set."))
      allTF = allTF[seq_len(nTFMax)]
      futile.logger::flog.info(paste0("Updated list of TFs: ", paste0(allTF, collapse = ", ")))
    } 
    
  }
  
  futile.logger::flog.info(paste0("Running the pipeline for ", nTF, " TF in total."))
  
  HOCOMOCO_mapping.df.exp = dplyr::filter(HOCOMOCO_mapping.df, .data$HOCOID %in% allTF)
  if (nrow(HOCOMOCO_mapping.df.exp) == 0) {
    message = paste0("Number of rows of HOCOMOCO_mapping.df.exp is 0. Something is wrong with the mapping table or the filtering")
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
  
  HOCOMOCO_mapping.df.exp
}

#' Overlap peaks and TFBS for a \code{\linkS4class{GRN}} object
#' 
#' @template GRN
#' @template nCores
#' @template forceRerun
#' @return An updated \code{\linkS4class{GRN}} object, with added data from this function (\code{GRN@data$TFs$TF_peak_overlap} in particular)
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' GRN = overlapPeaksAndTFBS(GRN, nCores = 2, forceRerun = FALSE)
#' @export
overlapPeaksAndTFBS <- function(GRN, nCores = 2, forceRerun = FALSE) {
  
  start = Sys.time()
    
  checkmate::assertClass(GRN, "GRN")
  GRN = .addFunctionLogToObject(GRN)
  
  GRN = .makeObjectCompatible(GRN)
  
  checkmate::assertIntegerish(nCores, lower = 1)
  checkmate::assertFlag(forceRerun)
  
  if (is.null(GRN@data$TFs$TF_peak_overlap) | forceRerun) {
    
    
    futile.logger::flog.info(paste0("Overlap peaks and TFBS using ", nCores, " cores. This may take a while, particularly if the number of samples is large..."))
    
    genomeAssembly = GRN@config$parameters$genomeAssembly
    seqlengths = .getChrLengths(genomeAssembly)
    
    if (!is.null(GRN@config$TFBS_filePattern)) {
      filesTFBSPattern = GRN@config$TFBS_filePattern
    } else {
      message = "Could not retrieve value from GRN@config$TFBS_filePattern. Please rerun the function addTFBS, as this was added in a recent version of the package."
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    
    # Check whether we have peaks on chromosomes not part of the sequence length reference. If yes, discard them
    annotation_discared = dplyr::filter(GRN@data$peaks$counts_metadata, ! .data$chr %in% names(seqlengths))
    
    if (nrow(annotation_discared) > 0) {
      
      tbl_discarded = table(annotation_discared$chr)
      tbl_discarded = tbl_discarded[which(tbl_discarded > 0)]
      
      message = paste0("Found ", sum(tbl_discarded), " regions from chromosomes without a reference length. ", 
                       "Typically, these are random fragments from known or unknown chromosomes. The following regions will be discarded: \n",
                       paste0(names(tbl_discarded), " (", tbl_discarded, ")", collapse = ","))
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)  
      
      GRN@data$peaks$counts_metadata = dplyr::filter(GRN@data$peaks$counts_metadata, .data$chr %in% names(seqlengths))
    }
    
    # Construct GRanges
    consensus.gr   = .constructGRanges(GRN@data$peaks$counts_metadata, seqlengths = seqlengths, genomeAssembly)
    
    res.l = .execInParallelGen(nCores, returnAsList = TRUE, listNames = GRN@config$allTF, 
                               iteration = seq_len(length(GRN@config$allTF)), 
                               verbose = FALSE, 
                               functionName = .intersectTFBSPeaks, GRN = GRN, consensusPeaks = consensus.gr, filesTFBSPattern = filesTFBSPattern)
    
    # Sanity check
    
    TFBS_bindingMatrix.df = tibble::as_tibble(res.l)
    
    if (!all(colnames(TFBS_bindingMatrix.df) %in% GRN@config$allTF)) {
      
      message = paste0("Internal mismatch detected between the TF names and the TF names derived from the translation file (see log, column HOCOID).", 
                       "This may happen if the genome assembly version has been changed, but intermediate files have not been properly recreated. ",
                       "Set the parameter forceRerun to TRUE and rerun the script.")
      
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    # new 
    filteredPeaks = dplyr::filter(GRN@data$peaks$counts_metadata, .data$isFiltered) %>% dplyr::pull(.data$peakID)
    # Collect binary 0/1 binding matrix from all TF and concatenate
    GRN@data$TFs$TF_peak_overlap = TFBS_bindingMatrix.df %>%
      dplyr::mutate(peakID = GRN@data$peaks$counts_metadata$peakID,
                    isFiltered = .data$peakID %in% filteredPeaks) %>%
      dplyr::mutate_if(is.logical, as.numeric) %>%
      dplyr::select(tidyselect::all_of(sort(GRN@config$allTF)), .data$isFiltered)
    
    GRN@data$TFs$TF_peak_overlap = .asSparseMatrix(as.matrix(GRN@data$TFs$TF_peak_overlap), 
                                                   convertNA_to_zero = FALSE, 
                                                   dimnames = list(GRN@data$peaks$counts_metadata$peakID, colnames(GRN@data$TFs$TF_peak_overlap)))
    
    # The order of rows is here the sorted version as it originates from the sorted consensus peak file
    # We resort it to match the countsPeaks.norm
    # TODO: Here could be an error due to differences in sorting. Also, whether or not all or the filtered version shall be used
    stopifnot(identical(rownames(GRN@data$TFs$TF_peak_overlap), GRN@data$peaks$counts_metadata$peakID))
    
  } 
  
  .printExecutionTime(start, prefix = "")
  
  GRN
}


#' @import GenomicRanges
.intersectTFBSPeaks <- function(GRN, TFIndex, consensusPeaks, filesTFBSPattern, verbose = FALSE) {
  
  TFCur = GRN@config$allTF[TFIndex]
  
  file_tfbs_in  = paste0(GRN@config$directories$motifFolder, .Platform$file.sep, TFCur, filesTFBSPattern, GRN@config$TFBS_fileEnding)
  
  # Intersect consensusPeaks GR with bed file GR 
  TFBS.df = .readTFBSFile(file_tfbs_in) 
  
  subject.gr  = .constructGRanges(TFBS.df, seqlengths = .getChrLengths(GRN@config$parameters$genomeAssembly), GRN@config$parameters$genomeAssembly)
  
  intersect.gr = GenomicRanges::intersect(subject.gr, consensusPeaks, ignore.strand=TRUE)
  
  query.gr = consensusPeaks
  
  overlapsAll = GenomicRanges::findOverlaps(query.gr, subject.gr, 
                                            minoverlap=1,
                                            type="any",
                                            select="all",
                                            ignore.strand=TRUE)
  
  query_row_ids  = S4Vectors::queryHits(overlapsAll)
  subject_rowids = S4Vectors::subjectHits(overlapsAll)
  
  subject_overlap_df = as.data.frame(S4Vectors::elementMetadata(subject.gr)[subject_rowids, ])
  subject_overlap_df$tfbs_chr = as.character(GenomeInfoDb::seqnames(subject.gr))[subject_rowids]
  subject_overlap_df$tfbs_start = start(subject.gr)[subject_rowids]
  subject_overlap_df$tfbs_end   = end(subject.gr)  [subject_rowids]
  
  query_overlap_df = as.data.frame(S4Vectors::elementMetadata(query.gr)  [query_row_ids, "peakID", drop = FALSE])
  query_overlap_df$peak_chr    = as.character(GenomeInfoDb::seqnames(query.gr))[query_row_ids]
  query_overlap_df$peak_start  = start(query.gr)[query_row_ids]
  query_overlap_df$peak_end    = end(query.gr)  [query_row_ids]
  
  final.df = cbind.data.frame(query_overlap_df, subject_overlap_df) %>%
    dplyr::select(-.data$score, -.data$annotation) %>%
    dplyr::mutate(tfbsID = paste0(.data$tfbs_chr, ":", .data$tfbs_start, "-", .data$tfbs_end),
                  coordCentTfbs = round((.data$tfbs_start + .data$tfbs_end)/2, 0),
                  coordSummit   = round((.data$peak_start + .data$peak_end)/2, 0),
                  distance = abs(.data$coordSummit - .data$coordCentTfbs))  %>%
    dplyr::group_by(.data$peakID) %>%
    dplyr::slice(which.min(.data$distance)) %>%
    #arrange(distance, .by_group = TRUE) %>%
    # top_n(n = 2, dplyr::desc(distance)) %>%
    dplyr::ungroup()
  
  futile.logger::flog.info(paste0(" Calculating intersection for TF ", TFCur, " finished. Number of overlapping TFBS after filtering: ", nrow(final.df)))
  
  
  return(GRN@data$peaks$counts_metadata$peakID %in% final.df$peakID)
  
}

.readTFBSFile <- function(file_tfbs_in) {
  
  TFBS.df = suppressMessages(.read_tidyverse_wrapper(file_tfbs_in, type = "tsv", col_names = FALSE, ncolExpected = 3:11, verbose = FALSE))
  if (ncol(TFBS.df) == 3) {
    colnames(TFBS.df) = c("chr", "start", "end")
  } else if (ncol(TFBS.df) == 4) {
    colnames(TFBS.df) = c("chr", "start", "end", "annotation")
  } else if (ncol(TFBS.df) == 5) {
    colnames(TFBS.df) = c("chr", "start", "end", "annotation", "strand")
  } else if (ncol(TFBS.df) >= 6) {
    
    if (ncol(TFBS.df) > 6) {
      message = paste0(" File ", file_tfbs_in, " had more than 6 columns, only the first 6 will be used.")
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)  
        
      TFBS.df = TFBS.df[,seq_len(6)]
    }
    colnames(TFBS.df) = c("chr", "start", "end", "annotation", "score", "strand")
  }
  
  TFBS.df
}



# TODO: Add columns for TF availability here also
# GRN@config$TF_list[["all_TFBS"]] =GRN@config$allTF
.correlateMatrices <- function(matrix1, matrix_peaks, HOCOMOCO_mapping, corMethod = "pearson", whitespacePrefix = " ") {
  
  start = Sys.time()
  
  # Set the column name to just ENSEMBL to avoid column mismatch issues
  HOCOMOCO_mapping$ENSEMBL = HOCOMOCO_mapping$TF.ENSEMBL
  
  # Filter to only the TFs
  # In addition, the no of TF because multiple TFs can map to the same gene/ ENSEMBL ID
  # Also filter 0 count genes because they otherwise cause errors downstream
  rowSums = rowSums(dplyr::select(matrix1, -.data$ENSEMBL))
  
  # Keep only Ensembl IDs from TFs we have data from
  matrix1.norm.TFs.df = dplyr::filter(matrix1, .data$ENSEMBL %in% HOCOMOCO_mapping$TF.ENSEMBL, rowSums != 0)
  
  nFiltered1 = dplyr::filter(matrix1, ! .data$ENSEMBL %in% HOCOMOCO_mapping$TF.ENSEMBL) %>% nrow()
  nFiltered2 = dplyr::filter(matrix1, rowSums == 0) %>% nrow()
  
  diff = nrow(matrix1) - nrow(matrix1.norm.TFs.df)
  if (diff > 0) {
    message = paste0(whitespacePrefix, "Retain ", nrow(matrix1.norm.TFs.df), " unique genes from TF/gene data out of ", nrow(matrix1), 
                     " (filter ",  nFiltered1, " non-TF genes and ", nFiltered2, 
                     " TF genes with 0 counts throughout).")
    futile.logger::flog.info(message)
  }
  
  if (nrow(matrix1.norm.TFs.df) == 0) {
    message = " No rows remaining from TF/gene data after filtering against ENSEMBL IDs from HOCOMOCO. Check your ENSEMBL IDs for overlap with the HOCOMOCO translation table."
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
  HOCOMOCO_mapping.exp = dplyr::filter(HOCOMOCO_mapping, .data$TF.ENSEMBL %in% matrix1.norm.TFs.df$ENSEMBL)
  futile.logger::flog.info(paste0(whitespacePrefix, "Correlate TF/gene data for ", nrow(matrix1.norm.TFs.df), " unique Ensembl IDs (TFs) and peak counts for ", nrow(matrix_peaks), " peaks."))
  futile.logger::flog.info(paste0(whitespacePrefix, "Note: For subsequent steps, the same gene may be associated with multiple TF, depending on the translation table."))
  # Correlate TF gene counts with peak counts 
  # matrix1:  rows: all TF genes, columns: all samples
  # matrix_peaks: rows: peak IDs, columns: samples
  # Transverse both for the cor function then
  
  # counts for peaks may be 0 throughout, then a warning is thrown
  
  #If the sd is zero, a warning is issued. We suppress it here to not confuse users as this is being dealt with later by ignoring the NA entries
  cor.m = suppressWarnings(t(cor(t(dplyr::select(matrix1.norm.TFs.df, -.data$ENSEMBL)), t(dplyr::select(matrix_peaks, -.data$peakID)), method = corMethod)))
  
  colnames(cor.m) = matrix1.norm.TFs.df$ENSEMBL
  rownames(cor.m) = matrix_peaks$peakID
  
  # Some entries in the HOCOMOCO mapping can be repeated (i.e., the same ID for two different TFs, such as ZBTB4.S and ZBTB4.D)
  # Originally, we deleted these rows from the mapping and took the first entry only
  # However, since TFs with the same ENSEMBL ID can still be different with respect to their TFBS, we now duplicate such genes also in the correlation table
  #HOCOMOCO_mapping.exp = HOCOMOCO_mapping.exp[!duplicated(HOCOMOCO_mapping.exp[, c("ENSEMBL")]),]
  #checkmate::assertSubset(as.character(HOCOMOCO_mapping.exp$ENSEMBL), colnames(sort.cor.m))
  
  # If a peak has identical counts across all samples,
  sort.cor.m = cor.m[,names(sort(colMeans(cor.m, na.rm = TRUE)))] 
  # Change the column names from ENSEMBL ID to TF names. 
  # Reorder to make sure the order is the same. Due to the duplication ID issue, the number of columns may increase after the column selection
  
  # Some columns may be removed here due to zero standard deviation
  HOCOMOCO_mapping.exp.filt = HOCOMOCO_mapping.exp %>% dplyr::filter(.data$TF.ENSEMBL %in% colnames(sort.cor.m))
  
  sort.cor.m = sort.cor.m[,as.character(HOCOMOCO_mapping.exp.filt$ENSEMBL)] 
  colnames(sort.cor.m) = as.character(HOCOMOCO_mapping.exp.filt$TF.HOCOID)
  
  .printExecutionTime(start, prefix = whitespacePrefix)
  sort.cor.m
}

.filterSortAndShuffle_peakTF_overlapTable <- function(GRN, perm, TF_peak_cor = NULL, shuffle = TRUE) {
  
  peak_TF_overlapCur.df = .asMatrixFromSparse(GRN@data$TFs$TF_peak_overlap, convertZero_to_NA = FALSE) %>% 
    tibble::as_tibble() %>%
    dplyr::filter(!.data$isFiltered) %>%  # Works because 1 / 0 is interpreted here as logical and not 1/0
    dplyr::select(-.data$isFiltered) 
  
  if (perm > 0 & shuffle) {
    peak_TF_overlapCur.df = .shuffleRowsPerColumn(peak_TF_overlapCur.df)
  }
  
  peak_TF_overlapCur.df
  
  
}



###### TF Activity functions and other data types ######

#' Add TF activity data to GRN object using a simplified procedure for estimating it. EXPERIMENTAL.
#' 
#' We do not yet provide full support for this function. It is currently being tested. Use at our own risk.
#' 
#' @template GRN
#' @template normalization_TFActivity
#' @template name_TFActivity
#' @template forceRerun
#' @return An updated \code{\linkS4class{GRN}} object, with added data from this function
#' (\code{GRN@data$TFs[[name]]} in particular, with \code{name} referring to the value of tje \code{name} parameter) 

addData_TFActivity <- function(GRN, normalization = "cyclicLoess", name = "TF_activity", forceRerun = FALSE) {
  
  checkmate::assertClass(GRN, "GRN")
  GRN = .addFunctionLogToObject(GRN)
  start = Sys.time()
  
  
  checkmate::assertChoice(normalization, c("cyclicLoess", "sizeFactors", "quantile", "none"))
  checkmate::assertCharacter(name, min.chars = 1, len = 1)
  checkmate::assertFlag(forceRerun)
  
  message = paste0("This function is currently under development and testing and not fully functional yet.")
  .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  
  forbiddenNames = "expression"
  .checkForbiddenNames(name, forbiddenNames)
  
  if (is.null(GRN@data$TFs[[name]]) | forceRerun) {
    
    futile.logger::flog.info(paste0("Calculate sample-specific TF activity from peaks data. This may take a while."))
    
    
    # Input: Raw peak counts per TF and TF-binding matrix
    # TODO: How to add raw counts here
    counts.df = getCounts(GRN, type = "peaks", permuted = FALSE) %>%
      tibble::as_tibble()
    
    # TODO replace by getter
    countsPeaks = .normalizeNewDataWrapper(GRN@data$peaks$counts_orig, normalization = normalization)
    
    stopifnot(identical(nrow(countsPeaks), nrow(GRN@data$TFs$TF_peak_overlap)))
    
    #Select a maximum set of TFs to run this for
    allTF = GRN@annotation$TFs$TF.name
    
    # rownamesTFs = GRN@annotation$TFs$ENSEMBL[match(allTF, GRN@annotation$TFs$HOCOID)] 
    
    # Calculating TF activity is done for all TF that are available
    TF.activity.m = matrix(NA, nrow = length(allTF), ncol = length(GRN@config$sharedSamples), 
                           dimnames = list(allTF, GRN@config$sharedSamples))
    
    
    pb <- progress::progress_bar$new(total = length(allTF))
    
    for (TFCur in allTF) {
      
      pb$tick()
      
      # Filter count matrix to those peaks with TFBS
      # Remove names from vector
      rows1 = as.vector(which(GRN@data$TFs$TF_peak_overlap[, TFCur] == 1))
      
      # Derive normalized counts for all peaks from the foreground (i.e., peaks with a predicted TFBS)
      Peaks_raw.cur.fg = countsPeaks[rows1,]
      
      # Derive z-scores
      scaled = t(scale(t(Peaks_raw.cur.fg)))
      
      colmeansCur = colMeans(scaled)
      TF.activity.m[TFCur, ] = colmeansCur
      stopifnot(identical(names(colmeansCur), GRN@config$sharedSamples))
    }
    
    # Store as data frame with both TF names and Ensembl IDs, in analogy to the other types of TF data that can be imported
    GRN@data$TFs[[name]] = TF.activity.m %>%
      as.data.frame() %>% 
      tibble::rownames_to_column("TF.name") %>%
      tibble::as_tibble() %>%
      dplyr::left_join(GRN@annotation$TFs, by = "TF.name") %>%
      dplyr::select(.data$ENSEMBL, .data$TF.name, tidyselect::all_of(GRN@config$sharedSamples))
    
    # Update available connection types
    GRN@config$TF_peak_connectionTypes = unique(c(GRN@config$TF_peak_connectionTypes, name))
    futile.logger::flog.info(paste0("TF activity successfully calculated. Data has been stored under the name ", name))
    
  } else {
    
    futile.logger::flog.info(paste0("Data already exists in object (slot ", name, "), nothing to do. Set forceRerun = TRUE to regenerate and overwrite."))
  }
  
  
  .printExecutionTime(start)
  
  GRN
  
}

.normalizeNewDataWrapper <- function(data, normalization, idColumn = "peakID") {
  
  if (checkmate::testClass(data, "DESeqDataSet")) {
    
    counts_raw = DESeq2::counts(data, normalized = FALSE)
    
  } else {
    
    # Capture incompatible cases
    if (normalization == "cyclicLoess" | normalization == "sizeFactors") {
      message = paste0("Selected normalization method for TF activity (", normalization, ") cannot be performed because the provided counts are not integer only. Select either \"quantile\" or \"none\" as normalization method for TF activity.")
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    counts_raw = data
  }
  
  
  if (normalization == "cyclicLoess") {
    
    futile.logger::flog.info(paste0(" Normalizing data using cyclic LOESS"))
      
      
    packageMessage = paste0("The package csaw is not installed but required for the cyclic LOESS normalization. Please install it and re-run this function or change the normalization type (if possible).")
    .checkPackageInstallation("csaw", packageMessage)   
      
    # Perform a cyclic loess normalization
    # We use a slighlty more complicated setup to derive size factors for library normalization
    # Instead of just determining the size factors in DeSeq2 via cirtual samples, we use 
    # a normalization from the csaw package (see https://www.rdocumentation.org/packages/csaw/versions/1.6.1/topics/normOffsets)
    # and apply a non-linear normalization. 
    # For each sample, a lowess curve is fitted to the log-counts against the log-average count. 
    # The fitted value for each bin pair is used as the generalized linear model offset for that sample. 
    # The use of the average count provides more stability than the average log-count when low counts are present for differentially bound regions.
    
    # since counts returns,by default, non-normalized counts, the following code should be fine and there is no need to also run estimateSizeFactors beforehand
    
    if (packageVersion("csaw") <= "1.14.1") {
      normFacs = exp((csaw::normOffsets(data, lib.sizes = colSums(data), type = "loess")))
    } else {
      object = SummarizedExperiment::SummarizedExperiment(list(counts=counts_raw))
      object$totals = colSums(counts_raw)
      normFacs  = exp(csaw::normOffsets(object, se.out = FALSE))
    }
    
    rownames(normFacs) = rownames(data)
    colnames(normFacs) = colnames(data)
    
    # We now provide gene-specific normalization factors for each sample as a matrix, which will preempt sizeFactors
    DESeq2::normalizationFactors(data) <- normFacs
    dataNorm = DESeq2::counts(data, normalized=TRUE)
    
  } else if (normalization == "sizeFactors") {
    
    futile.logger::flog.info(paste0(" Normalizing data using DESeq size factors"))
    data = DESeq2::estimateSizeFactors(data)
    dataNorm = DESeq2::counts(data, normalized=TRUE)
    
  } else if (normalization == "quantile") {
    
    futile.logger::flog.info(paste0(" Normalizing data using quantile normalization"))
    dataNorm = .normalizeCounts(data, method = "quantile",  idColumn = idColumn)
    
  } else if (normalization == "none") {
    dataNorm = data
    futile.logger::flog.info(paste0(" Skip normalization."))
    # Nothing to do, leave countsPeaks as they are
  }
  
  dataNorm
}

#' Import externally derived TF Activity data. EXPERIMENTAL.
#' 
#' We do not yet provide full support for this function. It is currently being tested. Use at our own risk.
#' 
#' @template GRN
#' @param data Data frame. No default. Data with TF data.
#' @template name_TFActivity
#' @param idColumn Character. Default \code{ENSEMBL}. Name of the ID column. Must not be unique as some TFs may correspond to the same ID.
#' @param nameColumn Character. Default \code{TF.name}. Must be unique for each TF / row.
#' @template normalization_TFActivity
#' @template forceRerun
#' @return An updated \code{\linkS4class{GRN}} object, with added data from this function.  
importTFData <- function(GRN, data, name, idColumn = "ENSEMBL", nameColumn = "TF.name", normalization = "none", forceRerun = FALSE) {
  
  checkmate::assertClass(GRN, "GRN")
  GRN = .addFunctionLogToObject(GRN)
  
  start = Sys.time()
  
  
  checkmate::assertDataFrame(data, min.cols = 2, min.rows = 1)
  checkmate::assertChoice(idColumn, colnames(data))
  checkmate::assertChoice(nameColumn, colnames(data))
  checkmate::assertSubset(GRN@config$sharedSamples, colnames(data), empty.ok = FALSE)
  checkmate::assertCharacter(name, min.chars = 1, any.missing = FALSE, len = 1)
  checkmate::assertChoice(normalization, c("cyclicLoess", "sizeFactors", "quantile", "none"))
  
  
  if (is.null(GRN@data$TFs[[name]]) | forceRerun) {
    
    futile.logger::flog.info(paste0("Importing external TF data under the name ", name)) 
    
    # Check whether TF have been added already
    if (is.null(GRN@annotation$TFs)) {
      message = paste0("No TFBS afound in the object. Make sure to run the function addTFBS first.")
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    # Replace spaces
    name = stringr::str_replace(name, "\\s", "_")
    
    # Rename idColumn and name column to "default"
    data = dplyr::rename(data, ENSEMBL = !!(idColumn))
    idColumn = "ENSEMBL"
    
    forbiddenNames = c("TF_activity", "expression")
    .checkForbiddenNames(name, forbiddenNames)
    
    idColumns = idColumn
    data = dplyr::rename(data, TF.name = !!(nameColumn))
    idColumns = c(idColumns, "TF.name")
    
    # Check uniqueness of TF names
    if (length(unique(data$TF.name)) < nrow(data)) {
      message = "TF names must be unique, but at least 2 TFs have the same TF name or TF names are missing."
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    # Make sure the column is reset
    GRN@annotation$TFs[[paste0("TF.name.", name)]] = NULL
    
    # TODO: Repeated execution results in more and more rows
    GRN@annotation$TFs = dplyr::left_join(GRN@annotation$TFs, data[, idColumns], 
                                                     by = "TF.name", suffix = c("", paste0(".", name)))
    
    # data = dplyr::select(data, -tidyselect::one_of(nameColumn))
    
    # Only TF.names are unique
    countsNorm = .normalizeNewDataWrapper(data %>% dplyr::select(-.data$ENSEMBL), normalization = normalization, idColumn = "TF.name")
    
    # Check overlap of ENSEMBL IDs
    countsNorm$ENSEMBL = data$ENSEMBL
    
    nRowBefore = nrow(countsNorm)
    countsNorm.subset = dplyr::filter(countsNorm, .data$ENSEMBL %in% GRN@annotation$TFs$TF.ENSEMBL)
    nRowAfter = nrow(countsNorm.subset)
    if (nRowAfter == 0) {
      message = "No rows overlapping with translation table, check ENSEMBL IDs."
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    } else if (nRowAfter < nRowBefore) {
      message = paste0("Retain ", nRowAfter, " from ", nRowBefore, " rows after filtering for overlapping ENSEMBL IDs from the translation table. This will raise a warning but this is usually expected that some rows are filtered")
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
    }
    
    
    # Check overlap of sample names
    nColBefore = ncol(countsNorm.subset)
    countsNorm.subset = dplyr::select(countsNorm.subset, tidyselect::all_of(idColumns), tidyselect::all_of(GRN@config$sharedSamples))
    nColAfter = ncol(countsNorm.subset)
    if (nColBefore > nColAfter) {
      
      if (nColAfter == length(idColumns)) {
        message = "No samples overlapping."
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
      } else {
        message = "Not all samples overlap"
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
      }
    }
    
    # Ensembl IDs may not be unique as different TFs can have the same Ensembl ID. 
    # Therefore, use TF names as row names, same as with the TF Activity matrix
    
    GRN@data$TFs[[name]] = countsNorm.subset %>%
      dplyr::select(.data$ENSEMBL, .data$TF.name, tidyselect::all_of(GRN@config$sharedSamples)) %>%
      tibble::as_tibble()
    
    # Update available connection types
    GRN@config$TF_peak_connectionTypes = unique(c(GRN@config$TF_peak_connectionTypes, name))
    
    
  } else {
    
    futile.logger::flog.info(paste0("Data already exists in object, nothing to do. Set forceRerun = TRUE to regenerate and overwrite."))
    
  }
  
  .printExecutionTime(start)
  
  GRN
}

######## AR classification ########

#' Run the activator-repressor classification for the TFs for a \code{\linkS4class{GRN}} object
#' 
#' @template GRN
#' @param significanceThreshold_Wilcoxon Numeric[0,1]. Default 0.05. Significance threshold for Wilcoxon test that is run in the end for the final classification. See the Vignette and *diffTF* paper for details.
#' @param plot_minNoTFBS_heatmap Integer[1,]. Default 100. Minimum number of TFBS for a TF to be included in the heatmap that is part of the output of this function.
#' @param deleteIntermediateData \code{TRUE} or \code{FALSE}.  Default \code{TRUE}. Should intermediate data be deleted before returning the object after a successful run? Due to the size of the produced intermediate data, we recommend setting this to \code{TRUE}, but if memory or object size are not an issue, the information can also be kept.
#' @template plotDiagnosticPlots
#' @template outputFolder
#' @template corMethod
#' @template forceRerun
#' @return An updated \code{\linkS4class{GRN}} object, with additional information added from this function. 
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' # GRN = loadExampleObject()
#' # GRN = AR_classification_wrapper(GRN, outputFolder = ".", forceRerun = FALSE)
#' @export
AR_classification_wrapper<- function (GRN, significanceThreshold_Wilcoxon = 0.05, 
                                      plot_minNoTFBS_heatmap = 100, deleteIntermediateData = TRUE,
                                      plotDiagnosticPlots = TRUE, outputFolder= NULL,
                                      corMethod = "pearson",
                                      forceRerun = FALSE) {
  
  start = Sys.time()
    
  checkmate::assertClass(GRN, "GRN")
  GRN = .addFunctionLogToObject(GRN)
  
  GRN = .makeObjectCompatible(GRN)
  
  checkmate::assertNumber(significanceThreshold_Wilcoxon, lower = 0, upper = 1)
  checkmate::assertNumber(plot_minNoTFBS_heatmap, lower = 1)
  checkmate::assertFlag(deleteIntermediateData)
  checkmate::assertFlag(plotDiagnosticPlots)
  checkmate::assertChoice(corMethod, c("pearson", "spearman"))
  checkmate::assertFlag(forceRerun)
  
  outputFolder = .checkOutputFolder(GRN, outputFolder)
  
  GRN@data$TFs$classification$TF.translation.orig = GRN@annotation$TFs %>%
    dplyr::mutate(TF.name = .data$TF.HOCOID)
  
  if (is.null(GRN@data$TFs$TF_peak_overlap)) {
    message = paste0("Could not find peak - TF matrix. Run the function overlapPeaksAndTFBS first / again")
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
  GRN@config$parameters$internal$plot_minNoTFBS_heatmap = plot_minNoTFBS_heatmap
  
  allPermutations = 0:.getMaxPermutation(GRN)
  
  connectionTypes = as.character(unique(GRN@connections$TF_peaks[["0"]]$main$TF_peak.connectionType))
  
  for (connectionTypeCur in connectionTypes) {
    
    futile.logger::flog.info(paste0(" Connection type ", connectionTypeCur, "\n"))
    
    for (permutationCur in allPermutations) {
      
      futile.logger::flog.info(paste0(" ", .getPermStr(permutationCur), "\n"))
      permIndex = as.character(permutationCur)
      permSuffix = ifelse(permutationCur == 0, "", ".permuted")
      
      if (is.null(GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]])) {
        if (is.null(GRN@data$TFs$classification[[permIndex]])) {
          GRN@data$TFs$classification[[permIndex]] = list()
        }
        GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]] = list()
      }
      
      if (is.null(GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_cor_median_foreground) |
          is.null(GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_cor_median_background) |
          is.null(GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_peak_cor_foreground) |
          is.null(GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_peak_cor_background) |
          is.null(GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_peak_cor) |
          forceRerun
      ) {
        
        GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]] = list()
        
        if (connectionTypeCur == "expression") {
          counts1 = getCounts(GRN, type = "rna", permuted = as.logical(permutationCur))
          
        } else {
          
          # TF activity data
          counts1 = GRN@data$TFs[[connectionTypeCur]] %>% 
            dplyr::select(-.data$TF.name)
          
        } 
        
        futile.logger::flog.info(paste0(" Correlate ", connectionTypeCur, " and peak counts"))
        
        counts_peaks = getCounts(GRN, type = "peaks", permuted = FALSE)
        
        TF_peak_cor = .correlateMatrices(matrix1      = counts1, 
                                         matrix_peaks = counts_peaks, 
                                         GRN@annotation$TFs, corMethod)
        
        peak_TF_overlapCur.df = .filterSortAndShuffle_peakTF_overlapTable(GRN, permutationCur, TF_peak_cor)
        res.l = .computeForegroundAndBackgroundMatrices(peak_TF_overlapCur.df, TF_peak_cor)
        GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_cor_median_foreground = res.l[["median_foreground"]]
        GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_cor_median_background = res.l[["median_background"]]
        GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_peak_cor_foreground   = res.l[["foreground"]]
        GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_peak_cor_background   = res.l[["background"]]
        
        GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_peak_cor = TF_peak_cor
      }
      
      # Final classification: Calculate thresholds by calculating the quantiles of the background and compare the real values to the background
      # TODO: Clarify whether the default convertZero_to_NA=TRUE is really needed here
      if (is.null(GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$act.rep.thres.l) | forceRerun) {
        GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$act.rep.thres.l = 
          .calculate_classificationThresholds(.asMatrixFromSparse(GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_peak_cor_background), 
                                              GRN@config$parameters)
      }
      
      if (is.null(GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF.classification) | forceRerun) {
        
        GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF.classification = 
          .finalizeClassificationAndAppend(
            output.global.TFs = GRN@data$TFs$classification$TF.translation.orig %>% dplyr::mutate(TF = .data$TF.name), 
            median.cor.tfs = GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_cor_median_foreground, 
            act.rep.thres.l = GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$act.rep.thres.l, 
            par.l = GRN@config$parameters, 
            t.cor.sel.matrix = .asMatrixFromSparse(GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_peak_cor_foreground), 
            t.cor.sel.matrix.non = .asMatrixFromSparse(GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_peak_cor_background), 
            significanceThreshold_Wilcoxon = significanceThreshold_Wilcoxon)
      }
      
      
      # PLOTS FOR THE RNA-SEQ CLASSIFICATION
      
      if (plotDiagnosticPlots) {
        
        outputFolder = .checkOutputFolder(GRN, outputFolder)
        
        suffixFile = .getPermutationSuffixStr(permutationCur)
        
        
        fileCur = paste0(outputFolder, .getOutputFileName("plot_class_density"), "_", connectionTypeCur, suffixFile, ".pdf")
        if (!file.exists(fileCur) | forceRerun) {
          .plot_density(.asMatrixFromSparse(GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_peak_cor_foreground),
                        .asMatrixFromSparse(GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_peak_cor_background), 
                        corMethod,
                        fileCur, width = 5, height = 5)
        } else {
          futile.logger::flog.info(paste0("  File ", fileCur, " already exists, not overwriting since forceRerun = FALSE"))
        }
        
        fileCur = paste0(outputFolder, .getOutputFileName("plot_class_medianClass"), "_", connectionTypeCur, suffixFile, ".pdf")
        if (!file.exists(fileCur) | forceRerun) {
          .plot_AR_thresholds(
            median.cor.tfs = .asMatrixFromSparse(GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_cor_median_foreground), 
            median.cor.tfs.non = .asMatrixFromSparse(GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_cor_median_background), 
            par.l = GRN@config$parameters, 
            act.rep.thres.l = GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$act.rep.thres.l, 
            corMethod = corMethod,
            file = fileCur,  width = 4, height = 8)
        } else {
          futile.logger::flog.info(paste0("  File ", fileCur, " already exists, not overwriting since forceRerun = FALSE"))
        }
        
        fileCur = paste0(outputFolder, .getOutputFileName("plot_class_densityClass"), "_", connectionTypeCur, suffixFile, ".pdf")
        if (!file.exists(fileCur) | forceRerun) {
          
          TF_peak_cor = GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_peak_cor
          peak_TF_overlapCur.df = .filterSortAndShuffle_peakTF_overlapTable(GRN, permutationCur, TF_peak_cor)
          .plot_heatmapAR(TF.peakMatrix.df = peak_TF_overlapCur.df, 
                          HOCOMOCO_mapping.df.exp = GRN@annotation$TFs %>% dplyr::mutate(TF = .data$TF.name), 
                          sort.cor.m = TF_peak_cor, 
                          par.l = GRN@config$parameters, 
                          corMethod = corMethod,
                          median.cor.tfs = .asMatrixFromSparse(GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_cor_median_foreground), 
                          median.cor.tfs.non = .asMatrixFromSparse(GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_cor_median_background), 
                          act.rep.thres.l = GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$act.rep.thres.l, 
                          finalClassification = GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF.classification,
                          file = fileCur, width = 8, height = 15)
        } else {
          futile.logger::flog.info(paste0("  File ", fileCur, " already exists, not overwriting since forceRerun = FALSE"))
        }
      }
      
      
      if (deleteIntermediateData) {
        
        GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_cor_median_foreground = NULL
        GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_cor_median_background = NULL
        GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_peak_cor_foreground = NULL
        GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_peak_cor_background = NULL
        GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$act.rep.thres.l = NULL
        
      } else {
        # Save as sparse matrices
        GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_peak_cor_foreground = 
          .asSparseMatrix(as.matrix(GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_peak_cor_foreground), convertNA_to_zero = TRUE)
        GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_peak_cor_background = 
          .asSparseMatrix(as.matrix(GRN@data$TFs$classification[[permIndex]] [[connectionTypeCur]]$TF_peak_cor_background), convertNA_to_zero = TRUE)
      }
      
    } # end for all permutations
    
  } # end of for each connection type
  
  .printExecutionTime(start, prefix = "")
  
  GRN
  
}

.checkAndUpdateConnectionTypes <- function(GRN) {
  
  if (is.null(GRN@config$TF_peak_connectionTypes)) {
    GRN@config$TF_peak_connectionTypes = "expression"
    if (!is.null(GRN@data$TFs$TF_activity)) {
      GRN@config$TF_peak_connectionTypes = c(GRN@config$TF_peak_connectionTypes, "TF_activity")
    }
  }
  GRN
}


######## Connections ########

#' Add TF-peak connections to a \code{\linkS4class{GRN}} object
#' 
#' After the execution of this function, QC plots can be plotted with the function \code{\link{plotDiagnosticPlots_TFPeaks}} unless this has already been done by default due to \code{plotDiagnosticPlots = TRUE}
#' 
#' @template GRN 
#' @template plotDiagnosticPlots
#' @template plotDetails
#' @template outputFolder
#' @template corMethod
#' @param connectionTypes Character vector. Default \code{expression}. Vector of connection types to include for the TF-peak connections. If an additional connection type is specified here, it has to be available already within the object (EXPERIMENTAL). See the function \code{\link{addData_TFActivity}} for details.
#' @param removeNegativeCorrelation  Vector of \code{TRUE} or \code{FALSE}. Default \code{FALSE}. EXPERIMENTAL. Must be a logical vector of the same length as the parameter \code{connectionType}. Should negatively correlated TF-peak connections be removed for the specific connection type? For connection type expression, the default is \code{FALSE}, while for any TF Activity related connection type, we recommend setting this to \code{TRUE}.  
#' @param maxFDRToStore Numeric[0,1]. Default 0.3. Maximum TF-peak FDR value to permanently store a particular TF-peak connection in the object? This parameter has a large influence on the overall memory size of the object, and we recommend not storing connections with a high FDR due to their sheer number.
#' @param addForPermuted \code{TRUE} or \code{FALSE}.  Default \code{FALSE}. Add connections also for permuted data. Leave at \code{TRUE} unless you know what you are doing.
#' @param useGCCorrection \code{TRUE} or \code{FALSE}.  Default \code{FALSE}. EXPERIMENTAL. Should a GC-matched background be used when calculating FDRs?
#' @param percBackground_size Numeric[0,100]. Default 75. EXPERIMENTAL. Description will follow. Only relevant if \code{useGCCorrection} is set to \code{TRUE}, ignored otherwise.
#' @param percBackground_resample \code{TRUE} or \code{FALSE}.  Default \code{TRUE}. EXPERIMENTAL. Should resampling be enabled for those GC bins for which not enough background peaks are available?. Only relevant if \code{useGCCorrection} is set to \code{TRUE}, ignored otherwise.
#' @template forceRerun
#' @seealso \code{\link{plotDiagnosticPlots_TFPeaks}}
#' @return An updated \code{\linkS4class{GRN}} object, with additional information added from this function. 
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' GRN = addConnections_TF_peak(GRN, plotDiagnosticPlots = FALSE, forceRerun = FALSE)
#' @export
addConnections_TF_peak <- function (GRN, plotDiagnosticPlots = TRUE, plotDetails = FALSE, outputFolder = NULL, 
                                    corMethod = "pearson", 
                                    connectionTypes = c("expression"),
                                    removeNegativeCorrelation = c(FALSE),
                                    maxFDRToStore = 0.3, 
                                    addForPermuted = TRUE,
                                    useGCCorrection = FALSE, percBackground_size = 75, percBackground_resample = TRUE,
                                    forceRerun = FALSE) {
  
  start = Sys.time()

  checkmate::assertClass(GRN, "GRN")
  GRN = .addFunctionLogToObject(GRN)
  
  GRN = .makeObjectCompatible(GRN)

  checkmate::assertFlag(plotDiagnosticPlots)
  checkmate::assertFlag(plotDetails)
  checkmate::assertChoice(corMethod, c("pearson", "spearman"))
  checkmate::assertFlag(addForPermuted)
  
  GRN = .checkAndUpdateConnectionTypes(GRN) # For compatibility with older versions
  checkmate::assertSubset(connectionTypes, GRN@config$TF_peak_connectionTypes, empty.ok = FALSE)
  
  checkmate::assertLogical(removeNegativeCorrelation, any.missing = FALSE, len = length(connectionTypes))
  
  
  #checkmate::assert(checkmate::checkSubset(add_TFActivity, c("none", "calculate", names(slot(GRN, "data")[["TFs"]]))))
  
  checkmate::assertNumber(maxFDRToStore, lower = 0, upper = 1)
  checkmate::assertFlag(useGCCorrection)
  checkmate::assertNumber(percBackground_size, lower = 0, upper = 100)
  checkmate::assertFlag(percBackground_resample)
  checkmate::assertFlag(forceRerun)
  
  
  if (is.null(GRN@connections$TF_peaks) | forceRerun) {
    
    GRN@connections$TF_peaks = list()
    
    GRN@config$parameters$corMethod_TF_Peak = corMethod
    GRN@config$parameters$useGCCorrection = useGCCorrection 
    
    if (is.null(GRN@data$TFs$TF_peak_overlap)) {
      message = paste0("Could not find peak - TF matrix. Run the function overlapPeaksAndTFBS first / again")
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    for (permutationCur in 0:.getMaxPermutation(GRN)) {
        
      if(!addForPermuted & permutationCur != 0) {
          next
      }
      
      futile.logger::flog.info(paste0("\n", .getPermStr(permutationCur), "\n"))
      permIndex = as.character(permutationCur)
      
      resFDR.l  = .computeTF_peak.fdr(GRN, perm = permutationCur, connectionTypes = connectionTypes, corMethod = corMethod, 
                                      removeNegativeCorrelation = removeNegativeCorrelation, 
                                      maxFDRToStore = maxFDRToStore, useGCCorrection = useGCCorrection,
                                      percBackground_size = percBackground_size, 
                                      percBackground_resample = percBackground_resample, plotDetails = plotDetails)
      
      # TODO: remove extra columns again
      GRN@connections$TF_peaks[[permIndex]]$main            = .optimizeSpaceGRN(stats::na.omit(resFDR.l[["main"]]))
      GRN@connections$TF_peaks[[permIndex]]$connectionStats = resFDR.l[["connectionStats"]] 
      
      futile.logger::flog.info(paste0("Finished. Stored ", nrow(GRN@connections$TF_peaks[[permIndex]]$main), " connections with an FDR <= ", maxFDRToStore))
      
      
      # GC plots, empty when no GC correction should be done
      GRN@stats$plots_GC = resFDR.l[["plots_GC"]]
      rm(resFDR.l)
      
      
      
    } #end for each permutation
    
    
    if (plotDiagnosticPlots) {
      
      plotDiagnosticPlots_TFPeaks(GRN, outputFolder = outputFolder, plotDetails = FALSE, forceRerun = forceRerun)
      
    }
    
  }
  
  .printExecutionTime(start, prefix = "")
  
  GRN
}



.computeTF_peak.fdr <- function(GRN, perm, connectionTypes, corMethod = "pearson", useGCCorrection = FALSE, 
                                removeNegativeCorrelation, maxFDRToStore = 0.3 , percBackground_size = 75, percBackground_resample = TRUE, plotDetails = FALSE, backgroundSize_min = 1000) {
  
  start = Sys.time()
  checkmate::assertIntegerish(backgroundSize_min, lower = 100)
  
  if (plotDetails) {
    message = "Plotting details is not supported currently. Set plotDetails = FALSE."
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
  peak_TF_overlap.df= .filterSortAndShuffle_peakTF_overlapTable(GRN, perm)
  
  plots_GC.l = list()
  
  # Lists that contain all the data
  connections_all.l = list()
  connectionStats_all.l = list()
  
  # List of connection types for which r < 0 should be filtered
  connectionTypes_removeNegCor = connectionTypes[removeNegativeCorrelation]
  
  for (connectionTypeCur in connectionTypes) {
    
    futile.logger::flog.info(paste0("Calculate TF-peak links for connection type ", connectionTypeCur))
    start2 = Sys.time()
    
    if (connectionTypeCur == "expression") {
      
      counts_connectionTypeCur = getCounts(GRN, type = "rna", permuted = as.logical(perm))
      
    } else {
      
      # Keep only Ensembl ID here
      counts_connectionTypeCur = GRN@data$TFs[[connectionTypeCur]] %>% 
        dplyr::select(-.data$TF.name)
      
    } 
    
    futile.logger::flog.info(paste0(" Correlate ", connectionTypeCur, " and peak counts"))
    
    counts_peaks = getCounts(GRN, type = "peaks", permuted = FALSE)
    
    # Filtering of the matrices happens automatically within the next function
    peaksCor.m = .correlateMatrices( matrix1= counts_connectionTypeCur, 
                                     matrix_peaks = counts_peaks, 
                                     GRN@annotation$TFs, 
                                     corMethod,
                                     whitespacePrefix = "  ")
    

    allTF = intersect(colnames(peak_TF_overlap.df), colnames(peaksCor.m))
    checkmate::assertIntegerish(length(allTF), lower = 1)
    
    futile.logger::flog.info(paste0(" Run FDR calculations for ", length(allTF), " TFs for which TFBS predictions and ",
                                    connectionTypeCur, " data for the corresponding gene are available."))
    if (length(allTF) < ncol(peak_TF_overlap.df) | length(allTF) < ncol(peaksCor.m)) {
      
      TF_missing = setdiff(colnames(peak_TF_overlap.df), colnames(peaksCor.m))
      if (length(TF_missing) > 0) futile.logger::flog.info(paste0("  Skip the following ", length(TF_missing), " TF due to missing data or because they are marked as filtered: ", paste0(TF_missing, collapse = ",")))
    }
    
    stopifnot(nrow(peaksCor.m) == nrow(peak_TF_overlap.df))
    
    peak_TF_overlap.df <- peak_TF_overlap.df[,allTF]
    
    sort.cor.m.sort<-peaksCor.m[,colnames(peak_TF_overlap.df)]
    
    stopifnot(identical(colnames(sort.cor.m.sort), colnames(peak_TF_overlap.df)))
    
    futile.logger::flog.info(paste0("  Compute FDR for each TF. This may take a while..."))
    
    pb <- progress::progress_bar$new(total = length(allTF))
    
    peaksFiltered = GRN@data$peaks$counts_metadata %>% dplyr::filter(!.data$isFiltered) 
    
    if (!useGCCorrection) {
      minPerc = 100
      # Since we do not control for this, we set it to NA
      background_match_success = NA
    }
    
    for (TFCur in allTF) {
      
      pb$tick()
      #start = Sys.time()
      #for(TFCur in colnames(corr_TF.sort)){
      overlapYes = which(peak_TF_overlap.df[,TFCur]==1)
      overlapNo  = which(peak_TF_overlap.df[,TFCur]==0)
      
      tp <- sort.cor.m.sort[overlapYes, TFCur]
      n_tp       = length(tp)
      
      peaksForeground =                     peaksFiltered %>% dplyr::slice(overlapYes)
      peaksBackground = peaksBackgroundGC = peaksFiltered %>% dplyr::slice(overlapNo)
      
      nPeaksForeground = nrow(peaksForeground)
      nPeaksBackground = nPeaksBackgroundGC = nrow(peaksBackground)
      
      #.printExecutionTime(start, "Interval 1: ")
      #start = Sys.time()
      # GC-adjust background and select background regions according to foreground
      if (useGCCorrection) {
        
        fp_orig <- sort.cor.m.sort[overlapNo, TFCur]
        n_fp_orig  = length(fp_orig)
        
        # Get GC info from those peaks from the foreground
        GC_classes_foreground.df = peaksForeground %>%
          dplyr::group_by(.data$GC_class) %>%
          dplyr::summarise(n= dplyr::n(), peak_width_mean = mean(.data$peak_width), peak_width_sd = sd(.data$peak_width)) %>%
          dplyr::ungroup() %>% 
          tidyr::complete(.data$GC_class, fill = list(n = 0)) %>%
          dplyr::mutate(n_rel = .data$n / nPeaksForeground, type = "foreground") %>%
          dplyr::arrange(dplyr::desc(.data$n_rel))
        
        GC_classes_background.df = peaksBackground %>%
          dplyr::group_by(.data$GC_class) %>%
          dplyr::summarise(n= dplyr::n(), peak_width_mean = mean(.data$peak_width), peak_width_sd = sd(.data$peak_width)) %>%
          dplyr::ungroup() %>% 
          tidyr::complete(.data$GC_class, fill = list(n = 0)) %>%
          dplyr::mutate(n_rel = .data$n / nPeaksBackground, type = "background_orig")
        
        
        
        background_match_success = TRUE
        
        if (!is.null(percBackground_size)) {
          minPerc = percBackground_size
        } else {
          
          threshold_percentage = 0.05
          minPerc = .findMaxBackgroundSize(GC_classes_foreground.df, GC_classes_background.df, peaksBackground, threshold_percentage =  threshold_percentage)
          
          if (minPerc == 0) {
            #futile.logger::flog.warn(paste0(" Mimicking the foreground failed for TF ", TFCur, ". The background will only be approximated as good as possible using 5% of the peaks."))
            minPerc = 5
            background_match_success = FALSE
          }
        }
        
        targetNoPeaks = minPerc/100 * nPeaksBackground
        
        # Ensure a minimum no of points in background, even if this sacrifices the mimicking of the distributions
        if (targetNoPeaks < backgroundSize_min) {
          
          if (!is.null(percBackground_size)) {
            #futile.logger::flog.warn(paste0("Number of peaks in background is smaller than 1000 for TF ", TFCur, ". Increasing in steps of 5% until a value > 1000 is found.This warning results from a too low value for the parameter percBackground_size"))
          } else {
            background_match_success = FALSE
          }
          
          targetNoPeaksNew = targetNoPeaks
          while(targetNoPeaksNew < 1000) {
            minPerc = minPerc + 5
            targetNoPeaksNew = minPerc/100 * nPeaksBackground
          }
          
          
        }
        
        # Add sel. minPerc to table and calculate required frequencies
        
        GC_classes_all.df = dplyr::full_join(GC_classes_foreground.df, GC_classes_background.df, suffix = c(".fg",".bg"), by = "GC_class") %>%
          dplyr::mutate(maxSizeBackground = .data$n.bg / .data$n_rel.fg,
                        n.bg.needed = floor(.data$n_rel.fg * targetNoPeaks), 
                        n.bg.needed.perc = .data$n.bg / .data$n.bg.needed) 
        
        
        #futile.logger::flog.info(paste0( " GC-adjustment: Randomly select a total of ", round(targetNoPeaks,0), 
        #                  " peaks (", minPerc, " %) from the background (out of ", nrow(peaksBackground), 
        #                  " overall) in a GC-binwise fashion to mimick the foreground distribution"))
        
        # Now we know the percentage, lets select an appropriate background
        # Sample peaks from background for each GC-bin specifically
        peakIDsSel = c()
        for (i in seq_len(nrow(GC_classes_foreground.df))) {
          
          peaksBackgroundGCCur =  peaksBackground %>% dplyr::filter(.data$GC_class == GC_classes_foreground.df$GC_class[i])
          
          if (nrow( peaksBackgroundGCCur) == 0) {
            next
          }
          
          #Select the minimum, which for low % classes is smaller than the required number to mimic the foreground 100%
          if (percBackground_resample) {
            nPeaksCur = GC_classes_all.df$n.bg.needed[i]    
          } else {
            nPeaksCur = min(GC_classes_all.df$n.bg.needed[i], nrow(peaksBackgroundGCCur))
          }
          
          if (GC_classes_all.df$n.bg.needed[i] > nrow(peaksBackgroundGCCur)) {
            peakIDsSel = c(peakIDsSel, peaksBackgroundGCCur %>% dplyr::sample_n(nPeaksCur, replace = percBackground_resample) %>% dplyr::pull(.data$peakID))  
          } else {
            peakIDsSel = c(peakIDsSel, peaksBackgroundGCCur %>% dplyr::sample_n(nPeaksCur, replace = FALSE) %>% dplyr::pull(.data$peakID))
          }
          # Take a sample from the background, and record the row numbers
          
        }
        
        # We cannot simply select now the peakIDs, as some peaks may be present multiple times
        #peaksBackgroundGC = peaksBackground %>% dplyr::filter(peakID %in% peakIDsSel) 
        peaksBackgroundGC = peaksBackground[match(peakIDsSel, peaksBackground$peakID),] 
        
        
        if (is.null(plots_GC.l[[TFCur]])) {
          plots_GC.l[[TFCur]] = .generateTF_GC_diagnosticPlots(TFCur, GC_classes_foreground.df, GC_classes_background.df, GC_classes_all.df, peaksForeground, peaksBackground, peaksBackgroundGC)
        }
        
        
        # Select the rows by their peak IDs
        fp    = sort.cor.m.sort[peakIDsSel, TFCur]
        
        fp_orig = sort.cor.m.sort[overlapNo, TFCur]
        n_fp_orig  = length(fp_orig)
        
      } else {
        # TODO: Redudnant so far in this case
        fp    =  fp_orig = sort.cor.m.sort[overlapNo, TFCur]
        n_fp_orig  = length(fp_orig)
      }
      
      n_fp = length(fp)
      
      cor.peak.tf = tibble::tibble(peak.ID    = rownames(sort.cor.m.sort)[overlapYes], 
                                   TF_peak.r  = sort.cor.m.sort[overlapYes, TFCur],
                                   TF.name    = as.factor(TFCur),
                                   TF_peak.connectionType = as.factor(connectionTypeCur))
      
      #val.sign      = (median(tp) - median(fp))
      # val.sign_orig = (median(tp) - median(fp_orig))
      
      #.printExecutionTime(start, "Interval 2: ")
      
      # Determine unique levels so plotting is identical for all
      #seq_pos<-unique(as.character(cut(GRN@config$parameters$internal$stepsFDR, breaks = GRN@config$parameters$internal$stepsFDR,      right = FALSE, include.lowest = TRUE )))
      #seq_neg<-unique(as.character(cut(GRN@config$parameters$internal$stepsFDR, breaks = rev(GRN@config$parameters$internal$stepsFDR), right = TRUE,  include.lowest = TRUE )))
      
      for (directionCur in c("pos","neg")) {
        
        indexStr = paste0(connectionTypeCur, "_", TFCur, "_", directionCur)
        
        #start = Sys.time()
        n_tp2.vec = n_fp2.vec =  n_fp2_orig.vec = rep(NA_real_, length(GRN@config$parameters$internal$stepsFDR)) 
        
        if (directionCur == "pos") {
          
          stepsCur = GRN@config$parameters$internal$stepsFDR
          rightOpen = FALSE
        } else { 
          
          stepsCur = rev(GRN@config$parameters$internal$stepsFDR)  
          rightOpen = TRUE
        }
        
        # Unique necessary to eliminate a duplication for one bin [-1,-0.95]
        levelsBins = unique(as.character(cut(stepsCur, breaks = stepsCur, right = rightOpen, include.lowest = TRUE)))
        
        cor.peak.tf$TF_peak.r_bin <- as.character(cut(cor.peak.tf$TF_peak.r, breaks = stepsCur, right = rightOpen, include.lowest = TRUE))
        #.printExecutionTime(start, "Interval 3a: ")
        #start = Sys.time()
        
        
        i = 0
        for (thres in stepsCur) {
          i = i + 1
          
          # na.rm = TRUE for all sums here to make sure NAs will not cause a problem
          if (directionCur == "pos") {
            
            n_tp2.vec[i]      = sum(tp>=thres, na.rm = TRUE)
            n_fp2.vec[i]      = sum(fp>=thres, na.rm = TRUE)
            n_fp2_orig.vec[i] = sum(fp_orig>=thres, na.rm = TRUE)
            
          } else {
            
            n_tp2.vec[i]      = sum(tp<thres, na.rm = TRUE)
            n_fp2.vec[i]      = sum(fp<thres, na.rm = TRUE)
            n_fp2_orig.vec[i] = sum(fp_orig<thres, na.rm = TRUE)
          }
          
        }
        
        #.printExecutionTime(start, "Interval 3b_new: ")
        #start = Sys.time()
        
        # Normalize the false positives to make them comparable to the true positives by dividing by the ratio
        # The maximum number is then identical to the maximum for the true positives
        n_fp2_norm.vec   = (n_fp2.vec/(n_fp/n_tp))
        n_fp2_orig_norm.vec = (n_fp2_orig.vec/(n_fp_orig/n_tp))
        
        #TODO: Decide for a variant.  +1 for raw or unnormalized values?
        
        
        fdr.curve = tibble::tibble( 
          TF_peak.r_bin2  = stepsCur,
          tpvalue = n_tp2.vec, 
          
          fpvalue = n_fp2.vec,
          fpvalue_orig = n_fp2_orig.vec,
          
          fpvalue_norm = n_fp2_norm.vec,
          fpvalue_norm_orig = n_fp2_orig_norm.vec,
          
          TF_peak.fdr_orig   = (n_fp2_orig_norm.vec) / (n_fp2_orig_norm.vec + n_tp2.vec),
          TF_peak.fdr     = (n_fp2_norm.vec) / (n_fp2_norm.vec + n_tp2.vec),
          # TF_peak.fdr_orig   = (n_fp2_orig_norm.vec + 1) / (n_fp2_orig_norm.vec + n_tp2.vec + 1),
          # TF_peak.fdr     = (n_fp2_norm.vec +1) / (n_fp2_norm.vec + n_tp2.vec + 1),
          # 
          TF_peak.fdr_direction     = directionCur,
          TF_peak.r_bin      = 
            as.character(cut(.data$TF_peak.r_bin2, breaks = stepsCur, right = rightOpen, include.lowest = TRUE))
        ) 
        
        # Derive connection summaries for all TF for both directions
        # Get no. of connections per bin, here make sure to also include that have n = 0
        
        
        # Remove negatively correlated connections for the specific connection type for which this was asked for
        if (connectionTypeCur %in% connectionTypes_removeNegCor ) {
          
          # futile.logger::flog.info(paste0(" Remove negatively correlated TF-peak pairs for connection type ", connectionTypeCur))
          cor.peak.tf = dplyr::filter(cor.peak.tf, .data$TF_peak.r >= 0)
          
        }
        
        connectionStats_all.l[[indexStr]] =  cor.peak.tf %>%
          dplyr::group_by(.data$TF_peak.r_bin) %>%
          dplyr::summarise(n = dplyr::n()) %>%
          dplyr::ungroup() %>%
          dplyr::right_join(fdr.curve, by = "TF_peak.r_bin") %>%
          dplyr::mutate(n = tidyr::replace_na(.data$n, replace = 0), 
                        TF.name = as.factor(TFCur), 
                        TF_peak.connectionType = factor(connectionTypeCur, levels = connectionTypes),
                        TF_peak.fdr_direction  = factor(directionCur, levels = c("pos", "neg")),
                        TF_peak.r_bin = factor(.data$TF_peak.r_bin, levels = levelsBins),
                        
                        # Collect extra information, currently however a bit repetitively stored
                        nForeground              = nPeaksForeground,
                        nBackground              = nrow(peaksBackgroundGC),
                        nBackground_orig         = nPeaksBackground,
                        percBackgroundUsed       = minPerc,
                        background_match_success = background_match_success) %>%
          dplyr::select(.data$TF.name, .data$TF_peak.r_bin, 
                        .data$n, .data$tpvalue, .data$fpvalue, .data$fpvalue_norm, 
                        .data$TF_peak.fdr, 
                        .data$TF_peak.fdr_orig, .data$TF_peak.fdr_direction, 
                        .data$TF_peak.connectionType,
                        tidyselect::contains("ground")) %>%
          dplyr::rename(n_tp = .data$tpvalue, n_fp = .data$fpvalue, n_fp_norm = .data$fpvalue_norm) %>%
          dplyr::distinct() %>%
          dplyr::arrange(.data$TF_peak.r_bin)
        
        
        # Collect data for additional QC plots before they are filtered
        
        # Filter now high FDR connections to save space and time
        # DISCARD other rows altogether
        # Left join here is what we want, as we need this df only for "real" data
        tblFilt.df = dplyr::left_join(cor.peak.tf, fdr.curve, by = "TF_peak.r_bin") %>%
          dplyr::filter(.data$TF_peak.fdr <= maxFDRToStore) %>%
          dplyr::select(.data$TF.name, .data$TF_peak.r_bin, .data$TF_peak.r, .data$TF_peak.fdr, 
                        .data$TF_peak.fdr_orig, .data$peak.ID, .data$TF_peak.fdr_direction, 
                        .data$TF_peak.connectionType, tidyselect::contains("value"))
        
        
        if (!plotDetails) {
          tblFilt.df = dplyr::select(tblFilt.df, -tidyselect::contains("value"))
        }
        
        connections_all.l[[indexStr]] = tblFilt.df
        
        #.printExecutionTime(start, "Interval 4: ")
        #start = Sys.time()
        
      } # end for directionCur in c("pos", "neg")
      
    } # end for each TF
    
    
    .printExecutionTime(start2, prefix = "  ")
    
  } # end for each connectionType
  
  
  .printExecutionTime(start, prefix = "")
  list(main            = data.table::rbindlist(connections_all.l), 
       connectionStats = data.table::rbindlist(connectionStats_all.l), 
       plots_GC        = plots_GC.l
  )
}

.findMaxBackgroundSize <- function (GC_classes_foreground.df, GC_classes_background.df, peaksBackground, threshold_percentage = 0.05) {
  
  # Iterate over different background sizes
  minPerc = 0
  for (percCur in c(seq(100,10,-5),5)) {
    
    if (minPerc > 0) next
    targetNoPeaks = percCur/100 * nrow(peaksBackground)
    
    #futile.logger::flog.info(paste0("Percentage of background: ", percCur))
    
    #Check for each GC bin in the foreground, starting with the most abundant, whether we have enough background peaks to sample from
    
    # threshold_percentage: From which percentage of GC bin frequency from the foreground should the mimicking fail?
    #The motivation is that very small bins that have no weight in the foreground will not cause a failure of the mimicking
    
    
    for (i in seq_len(nrow(GC_classes_foreground.df))) {
      
      n_rel    = GC_classes_foreground.df$n_rel[i]
      GC_class.cur = GC_classes_foreground.df$GC_class[i]
      
      requiredNoPeaks = round(n_rel * targetNoPeaks, 0)
      # Check in background
      availableNoPeaks = GC_classes_background.df %>% 
        dplyr::filter(.data$GC_class == GC_class.cur) %>%
        dplyr::pull(.data$n)
      #futile.logger::flog.info(paste0(" GC.class ", GC.class.cur, ": Required: ", requiredNoPeaks, ", available: ", availableNoPeaks))
      if ( availableNoPeaks < requiredNoPeaks) {
        #futile.logger::flog.info(paste0("  Not enough"))
      }
      if (availableNoPeaks < requiredNoPeaks & n_rel > threshold_percentage) {
        #futile.logger::flog.info(paste0(" Mimicking distribution FAILED (GC class ", GC.class.cur, " could not be mimicked"))
        break
      }
      
      if (i == nrow(GC_classes_foreground.df)) {
        minPerc = percCur
        #futile.logger::flog.info(paste0("Found max. percentage of background that is able to mimick the foreground: ", percCur))
        
      }
      
    }  # end of  for (i in 1:nrow(GC_classes_foreground.df)) {
    
    
  } # end for all percentages
  
  minPerc
}






#' Add peak-gene connections to a \code{\linkS4class{GRN}} object
#' 
#' After the execution of this function, QC plots can be plotted with the function \code{\link{plotDiagnosticPlots_peakGene}} unless this has already been done by default due to \code{plotDiagnosticPlots = TRUE}
#' 
#' @export
#' @template GRN
#' @param  overlapTypeGene Character. \code{"TSS"} or \code{"full"}. Default \code{"TSS"}. If set to \code{"TSS"}, only the TSS of the gene is used as reference for finding genes in the neighborhood of a peak. If set to \code{"full"}, the whole annotated gene (including all exons and introns) is used instead. 
#' @template corMethod
#' @param  promoterRange Integer >=0. Default 250000. The size of the neighborhood in bp to correlate peaks and genes in vicinity. Only peak-gene pairs will be correlated if they are within the specified range. Increasing this value leads to higher running times and more peak-gene pairs to be associated, while decreasing results in the opposite.
#' @param TADs Data frame with TAD domains. Default \code{NULL}. If provided, the neighborhood of a peak is defined by the TAD domain the peak is in rather than a fixed-sized neighborhood. The expected format is a BED-like data frame with at least 3 columns in this particular order: chromosome, start, end, the 4th column is optional and will be taken as ID column. All additional columns as well as column names are ignored. For the first 3 columns, the type is checked as part of a data integrity check.
#' @template nCores
#' @template plotDiagnosticPlots
#' @param plotGeneTypes List of character vectors. Default \code{list(c("all"), c("protein_coding"))}. Each list element may consist of one or multiple gene types that are plotted collectively in one PDF. The special keyword \code{"all"} denotes all gene types that are found (be aware: this typically contains 20+ gene types, see \url{https://www.gencodegenes.org/pages/biotypes.html} for details).
#' @template outputFolder
#' @template addRobustRegression
#' @template forceRerun
#' @seealso \code{\link{plotDiagnosticPlots_peakGene}}
#' @return An updated \code{\linkS4class{GRN}} object, with additional information added from this function. 
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' GRN = addConnections_peak_gene(GRN, promoterRange=10000, plotDiagnosticPlots = FALSE)
addConnections_peak_gene <- function(GRN, overlapTypeGene = "TSS", corMethod = "pearson",
                                     promoterRange = 250000, TADs = NULL,
                                     nCores = 4, 
                                     plotDiagnosticPlots = TRUE, 
                                     plotGeneTypes = list(c("all"), c("protein_coding")), 
                                     outputFolder = NULL,
                                     addRobustRegression = FALSE,
                                     forceRerun = FALSE) {
  
  start = Sys.time() 
    
  checkmate::assertClass(GRN, "GRN")
  GRN = .addFunctionLogToObject(GRN)
  
  GRN = .makeObjectCompatible(GRN)
  
  checkmate::assertChoice(overlapTypeGene, c("TSS", "full"))
  checkmate::assertChoice(corMethod, c("pearson", "spearman"))
  checkmate::assertIntegerish(promoterRange, lower = 0)
  checkmate::assert(checkmate::testNull(TADs), checkmate::testDataFrame(TADs))
  checkmate::assertIntegerish(nCores, lower = 1)
  checkmate::assertFlag(plotDiagnosticPlots) 
  for (elemCur in plotGeneTypes) {
      checkmate::assertSubset(elemCur, c("all", unique(as.character(GRN@annotation$genes$gene.type))) %>% stats::na.omit(), empty.ok = FALSE)
  }
  
  checkmate::assert(checkmate::testNull(outputFolder), checkmate::testDirectoryExists(outputFolder))
  checkmate::assertFlag(addRobustRegression)
  checkmate::assertFlag(forceRerun)
  
  if (addRobustRegression) {
      packageMessage = paste0("The package robust is not installed, but needed here due to addRobustRegression = TRUE. Please install it and re-run this function or change addRobustRegression to FALSE.")
      .checkPackageInstallation("robust", packageMessage)  
  }
  
  # As this is independent of the underlying GRN, it has to be done only once
  
  if (is.null(GRN@connections$peak_genes[["0"]]) | forceRerun) {
    
    GRN@config$parameters$promoterRange = promoterRange
    GRN@config$parameters$corMethod_peak_gene = corMethod
    
    
    GRN@connections$peak_genes = list()
    
    
    # Check which gene types are available for the particular genome annotation
    # Use all of them to collect statistics. Filtering can be done later
    # Just remove NA
    gene.types = unique(GRN@annotation$genes$gene.type) %>% stats::na.omit()
    
    
    for (permutationCur in 0:.getMaxPermutation(GRN)) {
      
      futile.logger::flog.info(paste0("\n", .getPermStr(permutationCur), "\n"))
      
      if (permutationCur == 0) {
        futile.logger::flog.info(paste0("Calculate peak-gene correlations for neighborhood size ", promoterRange))
        randomizePeakGeneConnections = FALSE
      } else {
        futile.logger::flog.info(paste0("Calculate random peak-gene correlations for neighborhood size ", promoterRange))
        randomizePeakGeneConnections = TRUE
      }
      
      GRN@connections$peak_genes[[as.character(permutationCur)]] = 
        .calculatePeakGeneCorrelations(GRN, permutationCur,
                                       TADs = TADs,
                                       neighborhoodSize = promoterRange,
                                       gene.types = as.character(gene.types),
                                       corMethod = corMethod,
                                       randomizePeakGeneConnections = randomizePeakGeneConnections,
                                       overlapTypeGene,
                                       nCores = nCores,
                                       debugMode_nPlots = 0,
                                       addRobustRegression = addRobustRegression
        )
      
    }
    
  } 
  
  if (plotDiagnosticPlots) {
    
    plotDiagnosticPlots_peakGene(GRN, outputFolder, gene.types = plotGeneTypes, useFiltered = FALSE, forceRerun= forceRerun)
    
  }
  
  .printExecutionTime(start, prefix = "")
  
  GRN

}

.calculatePeakGeneOverlaps <- function(GRN, allPeaks, peak_TAD_mapping = NULL, neighborhoodSize = 250000, genomeAssembly, 
                                       gene.types, overlapTypeGene = "TSS", removeEnsemblNA = TRUE) {
  
  start = Sys.time()
  futile.logger::flog.info(paste0("Calculate peak gene overlaps..."))
  
  # EXTEND PEAKS #
  
  
  # Add Hi-C domain data to query metadata if available
  if (!is.null(peak_TAD_mapping)) {
    
    futile.logger::flog.info(paste0("Integrate Hi-C data to extend peaks"))
    
    futile.logger::flog.info(paste0(" For peaks overlapping multiple TADs, use the union of all to define the neighborhood"))
    peak_TAD_mapping = peak_TAD_mapping %>%
      dplyr::group_by(.data$peakID) %>%
      dplyr::mutate(tadStart2 = min(.data$tadStart), # If a peak overlaps multiple TADs, merge them
                    tadEnd2   = max(.data$tadEnd),  # If a peak overlaps multiple TADs, merge them
                    tad.ID_all = paste0(.data$tad.ID, collapse = "|")) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
    
    query   = .constructGRanges(peak_TAD_mapping, seqlengths = .getChrLengths(genomeAssembly), genomeAssembly)
    
    # Remove rows with NA for the TAD 
    peaksNATAD = which(is.na(query$tad.ID))
    if (length(peaksNATAD) > 0) {
      futile.logger::flog.info(paste0(" ", length(peaksNATAD), " out of ", length(query), " peaks will not be tested for gene associations because they had no associated TAD"))
      query = query[!is.na(query$tad.ID)]
    }
    
    # Store original start and end positions before modifying them
    query$orig_start = start(query)
    query$orig_end   = end(query)
    
    # Extend GRanges by integrating Hi-C data. Use the newly defined TAD coordinates
    start(query) = suppressMessages(query$tadStart2)
    end(query)   = suppressMessages(query$tadEnd2)
    
    
  } else {
    
    # Without Hi-C data, we simply extend the ranges by a user-defined amount of bp, 250 kb being the default
    futile.logger::flog.info(paste0("Extend peaks based on user-defined extension size of ", neighborhoodSize, " up- and downstream."))
    
    query   = .constructGRanges(allPeaks, seqlengths = .getChrLengths(genomeAssembly), genomeAssembly)
    
    # Store original start and end positions before modifying them
    query$orig_start = start(query)
    query$orig_end   = end(query)
    
    suppressWarnings({start(query) = start(query) - neighborhoodSize})
    suppressWarnings({end(query)   = end(query) + neighborhoodSize})
    
    # Correct negative 
  }
  
  # correct ranges if within the promoterRange from the chr. starts and ends
  query = GenomicRanges::trim(query)

  subject = GRN@annotation$genes %>%
      dplyr::select(-.data$gene.mean, -.data$gene.median, -.data$gene.CV) %>%
      dplyr::filter(!is.na(.data$gene.start), !is.na(.data$gene.end))
  
  if (!is.null(gene.types)) {
      if (! "all" %in% gene.types) {
          subject = dplyr::filter(subject, .data$gene.type %in% gene.types)
      }
  }
      
  
  requiredColnames = c("gene.ENSEMBL","gene.type", "gene.name")
  checkmate::assertSubset(requiredColnames, colnames(subject), empty.ok = FALSE)
  
  subject.gr = GenomicRanges::makeGRangesFromDataFrame(subject, keep.extra.columns = TRUE)

  if (overlapTypeGene == "TSS") {
    # Take only the 5' end of the gene (start site and NOT the full gene length)
    end(subject.gr) = start(subject.gr)
  } 
  
  overlapsAll = suppressWarnings(GenomicRanges::findOverlaps(query, subject.gr, 
                                            minoverlap=1,
                                            type="any",
                                            select="all",
                                            ignore.strand=TRUE))
  
  query_row_ids  = S4Vectors::queryHits(overlapsAll)
  subject_rowids = S4Vectors::subjectHits(overlapsAll)
  
  
  
  #subject_overlap_df = as.data.frame(S4Vectors::elementMetadata(subject)[subject_rowids, c("ENSEMBL","ENTREZID", "SYMBOL")])
  subject_overlap_df = as.data.frame(S4Vectors::elementMetadata(subject.gr)[subject_rowids, requiredColnames])
  subject_overlap_df$gene.chr = as.character(GenomeInfoDb::seqnames(subject.gr))[subject_rowids]
  subject_overlap_df$gene.start = start(subject.gr)[subject_rowids]
  subject_overlap_df$gene.end   = end(subject.gr)  [subject_rowids]
  
  # Some entries in here will have only NAs
  
  query_overlap_df = as.data.frame(S4Vectors::elementMetadata(query)  [query_row_ids, "peakID"])
  
  query_overlap_df$ext_peak_chr    = as.character(GenomeInfoDb::seqnames(query))[query_row_ids]
  query_overlap_df$ext_peak_start  = start(query)[query_row_ids]
  query_overlap_df$ext_peak_end    = end(query)  [query_row_ids]
  query_overlap_df$orig_peak_start = query$orig_start[query_row_ids]
  query_overlap_df$orig_peak_end   = query$orig_end  [query_row_ids]    
  
  overlaps.df = cbind.data.frame(query_overlap_df, subject_overlap_df)
  colnames(overlaps.df)[1] = c("peak.ID")
  
  # Always compute distance to 5' of the gene:gene.start
  
  overlaps.sub.df = overlaps.df %>%
    dplyr::distinct() %>%
    dplyr::mutate(peak_gene.distance = dplyr::case_when(gene.start >= orig_peak_start & gene.start <= orig_peak_end ~ 0L,
                                                        TRUE ~ pmin(abs(orig_peak_start - gene.start), abs(orig_peak_end - gene.start))))
  
  if (removeEnsemblNA) {
    overlaps.sub.df = dplyr::filter(overlaps.sub.df, !is.na(.data$gene.ENSEMBL))
  }
  
  .printExecutionTime(start)
  
  overlaps.sub.df
  
}



.calculatePeakGeneCorrelations <- function(GRN, perm,
                                           TADs = NULL, 
                                           mergeOverlappingTADs = FALSE, 
                                           neighborhoodSize = 250000,
                                           gene.types = c("protein_coding"),
                                           overlapTypeGene = "TSS",
                                           corMethod = "pearson",
                                           randomizePeakGeneConnections = FALSE, 
                                           nCores = 1,
                                           chunksize = 50000,
                                           addRobustRegression = TRUE,
                                           debugMode_nPlots = 0) {
  
  start.all = Sys.time()
  
  genomeAssembly = GRN@config$parameters$genomeAssembly
  

  if (is.null(nCores)) {
    nCores = 1
  }
  
  consensusPeaks = GRN@data$peaks$counts_metadata %>% dplyr::filter(!.data$isFiltered)
  # Preprocess TAD boundaries
  if (!is.null(TADs)) {
    
    futile.logger::flog.info(paste0("Integrate Hi-C data and overlap peaks and HI-C domains"))  
    
    # Check format
    checkmate::assertCharacter(TADs$X1)
    checkmate::assertIntegerish(TADs$X2, lower = 1)
    checkmate::assertIntegerish(TADs$X3, lower = 1)
    
    colnames(TADs)[seq_len(3)] = c("chr", "start", "end")
    
    if (ncol(TADs) < 4) {
      TADs = dplyr::mutate(TADs, ID = paste0(.data$chr, ":", .data$start, "-", .data$end)) 
    } else {
      colnames(TADs)[4] = "ID"
    }
    
    
    # Construct GRanges
    query   = .constructGRanges(consensusPeaks, seqlengths = .getChrLengths(genomeAssembly), genomeAssembly)
    subject = .constructGRanges(TADs, seqlengths = .getChrLengths(genomeAssembly), genomeAssembly)
    
    TADOverlaps = GenomicRanges::countOverlaps(subject, subject)
    TADOverlaps_min2 = length(which(TADOverlaps > 1))
    
    futile.logger::flog.info(paste0(TADOverlaps_min2, "  TADs overlap each other"))
    
    # Merge overlapping TADs. min.gapwidth is set to 0 to prevent that directly adjacent TADs are merged
    if (mergeOverlappingTADs & TADOverlaps_min2 > 0) {
      futile.logger::flog.info(paste0("Merge overlapping TAD domains to one domain"))  
      subject = GenomicRanges::reduce(subject, min.gapwidth=0L)
      # Metadata has been lost, redefine it with the new boundaries
      subject$ID = paste0(as.character(GenomeInfoDb::seqnames(subject)), ":", start(subject), "-", end(subject))
    } else {
      futile.logger::flog.info(paste0("Overlapping TADs will not be merged"))  
    }
    
    
    # Check whether TAD boundaries overlap and print a warning if so
    nMultipleOverlaps = .checkSelfOverlap(subject)
    if (nMultipleOverlaps > 0) {
        message =paste0(nMultipleOverlaps, " out of ", length(subject), " TADs overlap with at least one other TAD. Please verify whether this is intended or a mistake. Particularly 1bp overlaps may not resembl the truth.")
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
    }
    
    
    # Finally, do the actual overlaps
    overlapsAll = suppressWarnings(GenomicRanges::findOverlaps(query, subject, 
                                              minoverlap = 1,
                                              type = "any",
                                              select = "all",
                                              ignore.strand = TRUE))
    
    
    query_row_ids  = S4Vectors::queryHits(overlapsAll)
    subject_rowids = S4Vectors::subjectHits(overlapsAll)
    
    subject_overlap_df = as.data.frame(S4Vectors::elementMetadata(subject)[subject_rowids, c("ID")]) %>%
      dplyr::mutate(tadChr = as.character(GenomeInfoDb::seqnames(subject))[subject_rowids],
                    tadStart = start(subject)[subject_rowids],
                    tadEnd = end(subject)[subject_rowids])
    # Some entries in here will have only NAs
    
    query_overlap_df   = as.data.frame(S4Vectors::elementMetadata(query)  [query_row_ids, "peakID"])
    
    overlaps.df = cbind.data.frame(query_overlap_df,subject_overlap_df)
    colnames(overlaps.df)[seq_len(2)] = c("peakID","tad.ID")
    
    peak.TADs.df = suppressWarnings(dplyr::left_join(consensusPeaks, overlaps.df, by = "peakID") )
    
    
    nPeaks = nrow(consensusPeaks)
    nPeaksWithOutTAD = length(which(is.na(peak.TADs.df$tad.ID)))
    futile.logger::flog.info(paste0(" Out of the ", nPeaks, " peaks, ", nPeaksWithOutTAD, " peaks are not within a TAD domain. These will be ignored for subsequent overlaps"))   
    
    nPeaksWithMultipleTADs = peak.TADs.df %>% dplyr::group_by(.data$peakID) %>% dplyr::summarize(n = dplyr::n()) %>% dplyr::filter(.data$n > 1) %>% nrow()
    
    if (nPeaksWithMultipleTADs > 0) {
      futile.logger::flog.info(paste0(" Out of the ", nPeaks, " peaks, ", nPeaksWithMultipleTADs, " overlap with more than one TAD. This usually means they are crossing TAD borders.")) 
    }
    
    
  } else {
    
    peak.TADs.df = NULL
  }
  # if a peak overlaps with a gene, should the same gene be reported as the connection?
  
  # OVERLAP OF PEAKS AND EXTENDED GENES
  overlaps.sub.df = .calculatePeakGeneOverlaps(GRN, allPeaks = consensusPeaks, peak.TADs.df, 
                                               neighborhoodSize = neighborhoodSize, 
                                               genomeAssembly = genomeAssembly, 
                                               gene.types = gene.types, overlapTypeGene = overlapTypeGene) 
  
  
  overlaps.sub.filt.df = overlaps.sub.df %>%
    dplyr::mutate(gene.ENSEMBL = gsub("\\..+", "", .data$gene.ENSEMBL, perl = TRUE)) # Clean ENSEMBL IDs
  
  # Only now we shuffle to make sure the background of possible connections is the same as in the foreground, as opposed to completely random
  # which would also include peak-gene connections that are not in the foreground at all
  if (randomizePeakGeneConnections) {
    futile.logger::flog.info(paste0(" Randomize gene-peak links by shuffling the peak IDs."))
    overlaps.sub.filt.df$peak.ID = sample(overlaps.sub.filt.df$peak.ID, replace = FALSE)
    overlaps.sub.filt.df$peak_gene.distance = NA
    
  }
  
  
  # Set to empty df to simplify the code below
  if (is.null(peak.TADs.df)) {
    peak.TADs.df = tibble::tibble(peak.ID = "", tad.ID = "")
  }
  
  # ITERATE THROUGH ALL PEAK_GENE PAIRS
  
  permIndex = as.character(perm)
  
  countsPeaks.clean = getCounts(GRN, type = "peaks",  permuted = FALSE, includeIDColumn = FALSE)
  countsRNA.clean   = getCounts(GRN, type = "rna", permuted = as.logical(perm), includeIDColumn = FALSE)
  
  # Cleverly construct the count matrices so we do the correlations in one go
  map_peaks = match(overlaps.sub.filt.df$peak.ID,  getCounts(GRN, type = "peaks", permuted = FALSE)$peakID)
  map_rna  = match(overlaps.sub.filt.df$gene.ENSEMBL, getCounts(GRN, type = "rna", permuted = as.logical(perm))$ENSEMBL) # may contain NA values because the gene is not actually in the RNA-seq counts
  
  # There should not b any NA because it is about the peaks
  stopifnot(all(!is.na(map_peaks)))
  # Some NAs might be expected, given our annotation contains all known genes
  stopifnot(!all(is.na(map_rna)))
  
  #res.m = matrix(NA, ncol = 2, nrow = nrow(overlaps.sub.filt.df), dimnames = list(NULL, c("p.raw", "peak_gene.r")))
  
  futile.logger::flog.info(paste0(" Iterate through ", nrow(overlaps.sub.filt.df), " peak-gene combinations and (if possible) calculate correlations using ", nCores, " cores. This may take a few minutes."))
  
  # parallel version of computing peak-gene correlations
  maxRow = nrow(overlaps.sub.filt.df)
  startIndexMax = ceiling(maxRow / chunksize) - 1 # -1 because we count from 0 onwards
  
  
  if (debugMode_nPlots > 0) {
    nCores = 1
  }
  
  
  res.l = .execInParallelGen(nCores, returnAsList = TRUE, listNames = NULL, iteration = 0:startIndexMax, verbose = FALSE, functionName = .correlateData, 
                             chunksize = chunksize, maxRow = maxRow, 
                             counts1 = countsPeaks.clean, counts2 = countsRNA.clean, map1 = map_peaks, map2 = map_rna, 
                             corMethod = corMethod, debugMode_nPlots = debugMode_nPlots, addRobustRegression = addRobustRegression)
  
  res.m  = do.call(rbind, res.l)
  
  futile.logger::flog.info(paste0(" Finished with calculating correlations, creating final data frame and filter NA rows due to missing RNA-seq data"))
  futile.logger::flog.info(paste0(" Initial number of rows: ", nrow(res.m)))
  # Neighborhood size not relevant for TADs
  if (!is.null(TADs)) {
    neighborhoodSize = -1
  }
  
  selectColumns = c("peak.ID", "gene.ENSEMBL", "peak_gene.distance", "tad.ID", "r", "p.raw")
  if (addRobustRegression) {
    selectColumns = c(selectColumns, "p_raw.robust", "r_robust", "bias_M_p.raw", "bias_LS_p.raw")
  }
  
  
  # Make data frame and adjust p-values
  res.df = suppressMessages(tibble::as_tibble(res.m) %>%
                              dplyr::mutate(peak.ID = getCounts(GRN, type = "peaks", permuted = FALSE)$peakID[map_peaks],
                                            gene.ENSEMBL = getCounts(GRN, type = "rna", permuted = as.logical(perm))$ENSEMBL[map_rna]) %>%
                              dplyr::filter(!is.na(.data$gene.ENSEMBL)) %>%  # For some peak-gene combinations, no RNA-Seq data was available, these NAs are filtered
                              # Add gene annotation and distance
                              dplyr::left_join(overlaps.sub.filt.df, by = c("gene.ENSEMBL", "peak.ID")) %>%
                              # Integrate TAD IDs also
                              dplyr::left_join(dplyr::select(peak.TADs.df, .data$peak.ID, .data$tad.ID), by = c("peak.ID")) %>%
                              
                              dplyr::select(tidyselect::all_of(selectColumns))) %>%
    dplyr::mutate(peak.ID = as.factor(.data$peak.ID),
                  gene.ENSEMBL = as.factor(.data$gene.ENSEMBL), 
                  tad.ID = as.factor(.data$tad.ID)) %>%
    dplyr::rename(peak_gene.r = .data$r, 
                  peak_gene.p_raw = .data$p.raw)
  
  if (addRobustRegression) {
    res.df = dplyr::rename(res.df, 
                           peak_gene.p_raw.robust = .data$p_raw.robust, 
                           peak_gene.r_robust = .data$r_robust,
                           peak_gene.bias_M_p.raw = .data$bias_M_p.raw,
                           peak_gene.bias_LS_p.raw = .data$bias_LS_p.raw)
  }
  
  if (is.null(TADs)) {
    res.df = dplyr::select(res.df, -.data$tad.ID)
  }
  
  futile.logger::flog.info(paste0(" Finished. Final number of rows: ", nrow(res.df)))
  
  
  .printExecutionTime(start.all)
  res.df
}

#' @import ggplot2
.correlateData <- function(startIndex, chunksize, maxRow, counts1, counts2, map1, map2, corMethod, debugMode_nPlots = 0, addRobustRegression = TRUE) {
  
  start = chunksize * startIndex + 1
  end = min(start +  chunksize - 1, maxRow)
  
  if (addRobustRegression) {
    res.m = matrix(NA, ncol = 6, nrow = end - start + 1, 
                   dimnames = list(NULL, c("p.raw", "r", "p_raw.robust", "r_robust", "bias_M_p.raw", "bias_LS_p.raw")))
  } else {
    res.m = matrix(NA, ncol = 2, nrow = end - start + 1, dimnames = list(NULL, c("p.raw", "r")))
  }
  
  
  rowCur = 0
  nPlotted = 0
  
  for (i in start:end) {
    
    rowCur = rowCur + 1
    
    if (is.na(map1[i]) | is.na(map2[i])) {
      next
    }
    
    data1   = unlist(counts1 [map1 [i],])
    data2   = unlist(counts2 [map2 [i],])
    
    res =  suppressWarnings(stats::cor.test(data1, data2, method = corMethod))
    
    res.m[rowCur, "p.raw"] = res$p.value
    res.m[rowCur, "r"]     = res$estimate
    
    if (addRobustRegression) {
        
      lmRob.sum = tryCatch( {
        summary(robust::lmRob(data1 ~data2))
        
      }, error = function(e) {
        # warning("Error running lmRob")
      }, warning = function(w) {
        # dont print anything
      }
      )
      
      if (methods::is(lmRob.sum, "summary.lmRob")) {
        res.m[rowCur, "p_raw.robust"]       = lmRob.sum$coefficients[2,4]
        # TODO: This is no r
        res.m[rowCur, "r_robust"]           = lmRob.sum$coefficients[2,1]
        res.m[rowCur, "bias_M_p.raw"]       = lmRob.sum$biasTest[1,2]
        res.m[rowCur, "bias_LS_p.raw"]      = lmRob.sum$biasTest[2,2]
      }
    }
    
    # https://stats.stackexchange.com/questions/205614/p-values-and-significance-in-rlm-mass-package-r
    
    
    if (nPlotted < debugMode_nPlots & res$p.value < 0.02) {
      dataCur.df = tibble::tibble(data_Peaks =unlist(counts1 [map1[i],]), 
                                  data_RNA = unlist(counts2 [map2[i],]) )
      
      # TODO: perm and GRN not known here
      stop("Not implemeneted yet")
      # g = ggplot2::ggplot(dataCur.df, ggplot2::aes(data1, data2)) + ggplot2::geom_point() + ggplot2::geom_smooth(method ="lm") + 
      #   ggplot2::ggtitle(paste0(getCounts(GRN, type = "peaks", permuted = FALSE)$peakID[map1[i]], 
      #                  ", ", 
      #                  getCounts(GRN, type = "rna", permuted = as.logical(perm))$ENSEMBL[map2 [i]], 
      #                  ", p = ", round(res$p.value, 3))) + 
      #   ggplot2::theme_bw()
      # plot(g)
      # 
      # nPlotted = nPlotted + 1
    }
    
    if (debugMode_nPlots > 0 & nPlotted >=  debugMode_nPlots){
      return(res.m)
    }
    
  }
  
  
  res.m
  
} # end function



#' Filter TF-peaks and peak-gene connections and combine them to TF-peak-gene connections to construct an eGRN.
#' 
#' This is one of the main integrative functions of the \code{GRaNIE} package. It has two main functions: 
#' First, filtering both TF-peak and peak-gene connections according to different criteria such as FDR and other properties 
#' Second, joining the three major elements that an eGRN consist of (TFs, peaks, genes) into one data frame, with one row per unique TF-peak-gene connection.
#' \strong{After successful execution, the connections (along with additional feature metadata) can be retrieved with the function \code{\link{getGRNConnections}}.}
#' \strong{Note that a previously stored eGRN graph is reset upon successful execution of this function along with printing a descriptive warning,
#'  and re-running the function \code{\link{build_eGRN_graph}} is necessary when any of the network functions of the package shall be executed. 
#' If the filtered connections changed, all network related enrichment functions also have to be rerun.}

#' Internally, before joining them, both TF-peak links and peak-gene connections are filtered separately for reasons of memory and computational efficacy:
#' First filtering out unwanted links dramatically reduces the memory needed for the full eGRN. Peak-gene p-value adjustment is only done after all filtering steps on the remaining set of
#' connections to lower the statistical burden of multiple-testing adjustment; therefore, this may lead to initially counter-intuitive effects such as a particular connections not being included anymore as compared to a 
#' filtering based on different thresholds, or the FDR being different for the same reason.
#' @template GRN
#' @param TF_peak.fdr.threshold Numeric[0,1]. Default 0.2. Maximum FDR for the TF-peak links. Set to 1 or NULL to disable this filter.
#' @template TF_peak.connectionTypes
#' @param peak_gene.p_raw.threshold Numeric[0,1]. Default NULL. Threshold for the peak-gene connections, based on the raw p-value. All peak-gene connections with a larger raw p-value will be filtered out.
#' @param peak_gene.fdr.threshold Numeric[0,1]. Default 0.2. Threshold for the peak-gene connections, based on the FDR. All peak-gene connections with a larger FDR will be filtered out.
#' @param peak_gene.fdr.method Character. Default "BH". One of: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none", "IHW". 
#' Method for adjusting p-values for multiple testing. 
#' If set to "IHW", the package \code{IHW} is required (as it is listed under \code{Suggests}, it may not be installed), 
#' and independent hypothesis weighting will be performed, and a suitable covariate has to be specified for the parameter \code{peak_gene.IHW.covariate}.
#' @param peak_gene.IHW.covariate Character. Default NULL. Name of the covariate to use for IHW (column name from the table thatis returned with the function \code{getGRNConnections}. Only relevant if \code{peak_gene.fdr.method} is set to "IHW". You have to make sure the specified covariate is suitable or IHW, see the diagnostic plots that are generated in this function for this. For many datasets, the peak-gene distance (called \code{peak_gene.distance} in the object) seems suitable.
#' @param peak_gene.IHW.nbins Integer or "auto". Default "auto". Number of bins for IHW. Only relevant if \code{peak_gene.fdr.method} is set to "IHW".
#' @template gene.types
#' @param allowMissingTFs \code{TRUE} or \code{FALSE}.  Default \code{FALSE}. Should connections be returned for which the TF is NA (i.e., connections consisting only of peak-gene links?). If set to \code{TRUE}, this generally greatly increases the number of connections but it may not be what you aim for.
#' @param allowMissingGenes \code{TRUE} or \code{FALSE}.  Default \code{TRUE}. Should connections be returned for which the gene is NA (i.e., connections consisting only of TF-peak links?). If set to \code{TRUE}, this generally increases the number of connections.
#' @param peak_gene.r_range Numeric(2). Default \code{c(0,1)}. Filter for lower and upper limit for the peak-gene links. Only links will be retained if the correlation coefficient is within the specified interval. This filter is usually used to filter out negatively correlated peak-gene links.
#' @param peak_gene.selection \code{"all"} or \code{"closest"}. Default \code{"all"}. Filter for the selection of genes for each peak. If set to \code{"all"}, all previously identified peak-gene are used, while \code{"closest"} only retains the closest gene for each peak that is retained until the point the filter is applied.
#' @param peak_gene.maxDistance Integer >0. Default \code{NULL}. Maximum peak-gene distance to retain a peak-gene connection.
#' @param filterTFs Character vector. Default \code{NULL}. Vector of TFs (as named in the GRN object) to retain. All TFs not listed will be filtered out.
#' @param filterGenes Character vector. Default \code{NULL}. Vector of gene IDs (as named in the GRN object) to retain. All genes not listed will be filtered out.
#' @param filterPeaks Character vector. Default \code{NULL}. Vector of peak IDs (as named in the GRN object) to retain. All peaks not listed will be filtered out.
#' @param TF_peak_FDR_selectViaCorBins \code{TRUE} or \code{FALSE}.  Default \code{FALSE}. Use a modified procedure for selecting TF-peak links that is based on the user-specified FDR but that retains also links that may have a higher FDR but a more extreme correlation.
#' @param silent \code{TRUE} or \code{FALSE}.  Default \code{FALSE}. Print progress messages and filter statistics.
#' @param resetGraphAndStoreInternally \code{TRUE} or \code{FALSE}.  Default \code{TRUE}. If set to \code{TRUE}, the stored eGRN graph slot graph) is reset due to the potentially changed connections that
#' would otherwise cause conflicts in the information stored in the object. Also, a GRN object is returned. If set to \code{FALSE}, only the new filtered connections are returned and the object is not altered.
#' @param filterLoops  \code{TRUE} or \code{FALSE}. Default \code{TRUE}. If a TF regulates itself (i.e., the TF and the gene are the same entity), should such loops be filtered from the GRN?
#' @template outputFolder
#' @return An updated \code{\linkS4class{GRN}} object, with additional information added from this function. 
#' The filtered and merged TF-peak and peak-gene connections in the slot \code{GRN@connections$all.filtered} and can be retrieved (along with other feature metadata) using the function \code{\link{getGRNConnections}}.
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' GRN = filterGRNAndConnectGenes(GRN)
#' @seealso \code{\link{visualizeGRN}}
#' @seealso \code{\link{addConnections_TF_peak}} 
#' @seealso \code{\link{addConnections_peak_gene}} 
#' @seealso \code{\link{build_eGRN_graph}} 
#' @seealso \code{\link{getGRNConnections}} 
#' @importFrom rlang .data
#' @importFrom magrittr `%>%`
#' @export
filterGRNAndConnectGenes <- function(GRN,
                                     TF_peak.fdr.threshold = 0.2, 
                                     TF_peak.connectionTypes = "all",
                                     peak_gene.p_raw.threshold = NULL, 
                                     peak_gene.fdr.threshold= 0.2,
                                     peak_gene.fdr.method = "BH",
                                     peak_gene.IHW.covariate = NULL,
                                     peak_gene.IHW.nbins = "auto",
                                     gene.types = c("protein_coding"), 
                                     allowMissingTFs = FALSE, allowMissingGenes = TRUE,
                                     peak_gene.r_range = c(0,1), 
                                     peak_gene.selection = "all",
                                     peak_gene.maxDistance = NULL,
                                     filterTFs = NULL, filterGenes = NULL, filterPeaks = NULL, 
                                     TF_peak_FDR_selectViaCorBins = FALSE,
                                     filterLoops = TRUE,
                                     outputFolder = NULL,
                                     resetGraphAndStoreInternally = TRUE,
                                     silent = FALSE) {
  
  start = Sys.time()  

  checkmate::assertClass(GRN, "GRN")
  GRN = .addFunctionLogToObject(GRN)
  
  GRN = .makeObjectCompatible(GRN)
  
  
  checkmate::assertCharacter(TF_peak.connectionTypes, min.len = 1, any.missing = FALSE)
  checkmate::assert(checkmate::checkNull(peak_gene.p_raw.threshold), checkmate::checkNumber(peak_gene.p_raw.threshold, lower = 0, upper = 1))
  checkmate::assertNumeric(peak_gene.r_range, lower = -1, upper = 1, len = 2)
  checkmate::assertSubset(gene.types, c("all", unique(as.character(GRN@annotation$genes$gene.type))) %>% stats::na.omit(), empty.ok = FALSE)
  checkmate::assertNumber(TF_peak.fdr.threshold, lower = 0, upper = 1)
  checkmate::assert(checkmate::checkNull(peak_gene.fdr.threshold), checkmate::checkNumber(peak_gene.fdr.threshold, lower = 0, upper = 1))
  
  checkmate::assertChoice(peak_gene.fdr.method, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none", "IHW"))
  checkmate::assert(checkmate::checkNull(peak_gene.IHW.covariate), checkmate::checkCharacter(peak_gene.IHW.covariate, min.chars = 1, len = 1))
  
  checkmate::assert(checkmate::checkIntegerish(peak_gene.IHW.nbins, lower = 1), checkmate::checkSubset(peak_gene.IHW.nbins, "auto"))
  
  checkmate::assertFlag(resetGraphAndStoreInternally)
  checkmate::assertFlag(silent)
  checkmate::assertChoice(peak_gene.selection, c("all", "closest"))
  checkmate::assertFlag(allowMissingTFs)
  checkmate::assertFlag(allowMissingGenes)
  
  checkmate::assertFlag(TF_peak_FDR_selectViaCorBins)
  checkmate::assertFlag(filterLoops)
  
  if (peak_gene.fdr.method == "IHW" & !is.installed("IHW")) {
    packageMessage = "IHW has been selected for p-value adjustment, but IHW is currently not installed. Please install it and re-run the function or choose a different method."
    .checkPackageInstallation("IHW", packageMessage)  

  }
  
  if (!is.null(peak_gene.p_raw.threshold) & !is.null(peak_gene.fdr.threshold)) {
    message = "Both parameters peak_gene.p_raw.threshold and peak_gene.fdr.threshold have been specified, choose only either of them."
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
  if (silent) {
    futile.logger::flog.threshold(futile.logger::WARN)
  }
  
  futile.logger::flog.info(paste0("Filter GRN network"))
  
  # Reset TF-gene links as these are recomputed with the filtered set and this cannot be done beforehand
  GRN@connections$TF_genes.filtered = NULL
  GRN@connections$all.filtered = list()
  
  for (permutationCur in 0:.getMaxPermutation(GRN)){
    
    futile.logger::flog.info(paste0("\n\n", .getPermStr(permutationCur)))
    permIndex = as.character(permutationCur)
    
    if (is.null(GRN@connections$peak_genes[[permIndex]])) {
      message = "No peak-gene connections found. Run the function addConnections_peak_gene first"
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    # Only select the absolute necessary here, no additional metadata
    ann.gene.red = GRN@annotation$genes %>%
        dplyr::mutate(gene.ENSEMBL = as.character(.data$gene.ENSEMBL)) %>%
        dplyr::select(.data$gene.ENSEMBL, .data$gene.name, .data$gene.type)
    
    peakGeneCorrelations = GRN@connections$peak_genes[[permIndex]] %>%
      dplyr::mutate(gene.ENSEMBL = as.character(.data$gene.ENSEMBL)) %>%
      dplyr::left_join(ann.gene.red, by = "gene.ENSEMBL")


    
    # Add TF Ensembl IDs
    grn.filt = GRN@connections$TF_peaks[[permIndex]]$main  %>% 
      tibble::as_tibble() %>%
      dplyr::left_join(GRN@annotation$TFs %>% dplyr::select(.data$TF.name, .data$TF.ENSEMBL), by = c("TF.name")) %>%
      dplyr::select(-.data$TF_peak.fdr_orig) %>%
      dplyr::mutate(TF.ENSEMBL = as.factor(.data$TF.ENSEMBL))
    
    if (is.null(grn.filt)) {
      message = "No TF-peak connections found. Run the function addConnections_TF_peak first"
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    # Filter network #
    futile.logger::flog.info(paste0("Inital number of rows left before all filtering steps: ", nrow(grn.filt)))
    
    if (! "all" %in% TF_peak.connectionTypes) {
      checkmate::assertSubset(TF_peak.connectionTypes, .getAll_TF_peak_connectionTypes(GRN), empty.ok = FALSE)
      futile.logger::flog.info(paste0(" Filter network and retain only rows with one of the following TF-peak connection types: ", paste0(TF_peak.connectionTypes, collapse = ", ")))
      futile.logger::flog.info(paste0("  Number of TF-peak rows before filtering connection types: ", nrow(grn.filt)))
      grn.filt = dplyr::filter(grn.filt, .data$TF_peak.connectionType %in% TF_peak.connectionTypes)
      futile.logger::flog.info(paste0("  Number of TF-peak rows after filtering connection types: ", nrow(grn.filt)))
    }
    
    if (!is.null(TF_peak.fdr.threshold)) {
      futile.logger::flog.info(paste0(" Filter network and retain only rows with TF-peak connections with an FDR < ", TF_peak.fdr.threshold))
      futile.logger::flog.info(paste0("  Number of TF-peak rows before filtering TFs: ", nrow(grn.filt)))
      
      if (!TF_peak_FDR_selectViaCorBins) {
        grn.filt = dplyr::filter(grn.filt, .data$TF_peak.fdr < TF_peak.fdr.threshold)
        futile.logger::flog.info(paste0("  Number of TF-peak rows after filtering TFs: ", nrow(grn.filt)))
      } else {
        
        # Add a new ID column
        grn.filt$row.ID = seq_len(nrow(grn.filt))
        
        # For each TF, identify those TF-peak correlation bins that are more extreme than the first correlation bin that is beyond the user-specified bin
        # Add one additional column to the table, and filter later by this column
        idsRowsKeep = c()
        for (TFCur in unique( grn.filt$TF.name)) {
          
          for (connectionTypeCur in TF_peak.connectionTypes) {
            
            grn.filt.TF = dplyr::filter(grn.filt, .data$TF.name == TFCur, .data$TF_peak.connectionType == connectionTypeCur)
            
            for (dirCur in c("pos", "neg")) {
              
              if (dirCur == "pos") {
                
                grn.filt.TF.dir = dplyr::filter(grn.filt.TF, .data$TF_peak.fdr_direction == "pos")
                stepsCur = GRN@config$parameters$internal$stepsFDR
                rightOpen = FALSE
                grn.filt.TF.dir$TF_peak.r_bin = cut(grn.filt.TF.dir$TF_peak.r, breaks = stepsCur, 
                                                    right = rightOpen, include.lowest = TRUE, ordered_result = TRUE)
                
                binThresholdNeg = grn.filt.TF.dir  %>% 
                  dplyr::select(.data$TF_peak.r_bin, .data$TF_peak.fdr) %>% 
                  dplyr::distinct() %>% 
                  dplyr::arrange(.data$TF_peak.r_bin)
                
                relBins =which(binThresholdNeg$TF_peak.fdr < TF_peak.fdr.threshold)
                if (length(relBins) > 0) {
                  relBins = binThresholdNeg$TF_peak.r_bin[min(relBins):nrow(binThresholdNeg)]
                  idsRowsKeep = c(idsRowsKeep, grn.filt.TF.dir %>% 
                                    dplyr::filter(.data$TF_peak.r_bin %in% relBins) %>% 
                                    dplyr::pull(.data$row.ID))
                }
                
                
              } else {
                
                grn.filt.TF.dir = dplyr::filter(grn.filt.TF, .data$TF_peak.fdr_direction == "neg")
                
                stepsCur = rev(GRN@config$parameters$internal$stepsFDR)  
                rightOpen = TRUE
                
                grn.filt.TF.dir$TF_peak.r_bin = cut(grn.filt.TF.dir$TF_peak.r, breaks = stepsCur, 
                                                    right = rightOpen, include.lowest = TRUE, ordered_result = TRUE)
                
                binThresholdNeg = grn.filt.TF.dir  %>% 
                  dplyr::select(.data$TF_peak.r_bin, .data$TF_peak.fdr) %>% 
                  dplyr::distinct() %>% 
                  dplyr::arrange(.data$TF_peak.r_bin)
                
                relBins =which(binThresholdNeg$TF_peak.fdr < TF_peak.fdr.threshold)
                if (length(relBins) > 0) {
                  relBins = binThresholdNeg$TF_peak.r_bin[seq_len(max(relBins))]
                  idsRowsKeep = c(idsRowsKeep, grn.filt.TF.dir %>% 
                                    dplyr::filter(.data$TF_peak.r_bin %in% relBins) %>% 
                                    dplyr::pull(.data$row.ID))
                }
              }
              
            } # end for both directions
          } # end for all TF-peak link types
        } # end for all TF
        
        
        grn.filt = grn.filt %>%
          dplyr::filter(.data$row.ID %in% idsRowsKeep) %>%
          dplyr::select(-.data$row.ID)
        
        futile.logger::flog.info(paste0("  Number of TF-peak rows after filtering TFs: ", nrow(grn.filt)))
        
      }
      
    }
    
    
    if (!is.null(filterTFs)) {
      futile.logger::flog.info(paste0(" Filter network to the following TF: ", paste0(filterTFs, collapse = ",")))
      futile.logger::flog.info(paste0("  Number of TF-peak rows before filtering TFs: ", nrow(grn.filt)))
      grn.filt = dplyr::filter(grn.filt, .data$TF.name %in% filterTFs)
      futile.logger::flog.info(paste0("  Number of TF-peak rows after filtering TFs: ", nrow(grn.filt)))
    }
    
    if (!is.null(filterPeaks)) {
      futile.logger::flog.info(paste0(" Filter network to the following peaks: ", paste0(filterPeaks, collapse = ",")))
      futile.logger::flog.info(paste0("  Number of TF-peak rows before filtering peaks: ", nrow(grn.filt)))
      grn.filt = dplyr::filter(grn.filt, .data$peak.ID %in% filterPeaks)
      futile.logger::flog.info(paste0("  Number of TF-peak rows after filtering peaks: ", nrow(grn.filt)))
    }
    
    # Filters on peak-genes
    
    futile.logger::flog.info("2. Filter peak-gene connections")
    
    if (!is.null(filterGenes)) {
        futile.logger::flog.info(paste0(" Filter peak-gene connections for the following genes: ", paste0(filterGenes, collapse = ",")))
        futile.logger::flog.info(paste0("  Number of rows before filtering genes: ", nrow(peakGeneCorrelations)))
        peakGeneCorrelations = dplyr::filter(peakGeneCorrelations, .data$gene.ENSEMBL %in% filterGenes)
        futile.logger::flog.info(paste0("  Number of rows after filtering genes: ", nrow(peakGeneCorrelations)))
    }
    
    if (!is.null(peak_gene.maxDistance)) {
      checkmate::assertIntegerish(peak_gene.maxDistance, lower = 0)
      futile.logger::flog.info(paste0(" Filter peak-gene connections for their distance and keep only connections with a maximum distance of  ", peak_gene.maxDistance))
      futile.logger::flog.info(paste0("  Number of peak-gene rows before filtering connection types: ", nrow(peakGeneCorrelations)))
      peakGeneCorrelations = dplyr::filter(peakGeneCorrelations, .data$peak_gene.distance < peak_gene.maxDistance)
      futile.logger::flog.info(paste0("  Number of peak-gene rows after filtering connection types: ", nrow(peakGeneCorrelations)))
    }
    
    if (! "all" %in% gene.types) {
      futile.logger::flog.info(paste0(" Filter genes by gene type, keep only the following gene types: ", paste0(gene.types, collapse = ", ")))
      futile.logger::flog.info(paste0("  Number of peak-gene rows before filtering by gene type: ", nrow(peakGeneCorrelations)))
      peakGeneCorrelations = dplyr::filter(peakGeneCorrelations, .data$gene.type %in% gene.types)
      futile.logger::flog.info(paste0("  Number of peak-gene rows after filtering by gene type: ", nrow(peakGeneCorrelations)))
    }
    
   
    
    
    futile.logger::flog.info(paste0("3. Merging TF-peak with peak-gene connections and filter the combined table..."))
    # Now we need the connected genes. All fitters that are independent of that have been done
    # Don't warn about the coercing of factors etc
    
    if (allowMissingTFs) {
      grn.filt = suppressWarnings(dplyr::full_join(grn.filt, peakGeneCorrelations, by = "peak.ID"))
    } else {
      grn.filt = suppressWarnings(dplyr::left_join(grn.filt, peakGeneCorrelations, by = "peak.ID"))
    }
    
    
    futile.logger::flog.info(paste0("Inital number of rows left before filtering steps: ", nrow(grn.filt)))
    
    if (filterLoops) {
      futile.logger::flog.info(paste0(" Filter TF-TF self-loops"))
      futile.logger::flog.info(paste0("  Number of rows before filtering genes: ", nrow(grn.filt)))
      
      # Be aware of NA values here in the selection, depending on allowMissingTFs
      grn.filt = dplyr::filter(grn.filt, 
                               is.na(.data$TF.ENSEMBL) | 
                               (!is.na(.data$TF.ENSEMBL) & (as.character(.data$gene.ENSEMBL) != as.character(.data$TF.ENSEMBL))))
      
      futile.logger::flog.info(paste0("  Number of rows after filtering genes: ", nrow(grn.filt)))
    }
    
    
    
    if (allowMissingGenes) {
      # Nothing to do here
    } else {
      
      futile.logger::flog.info(paste0(" Filter rows with missing ENSEMBL IDs"))
      futile.logger::flog.info(paste0("  Number of rows before filtering: ", nrow(grn.filt)))
      grn.filt = dplyr::filter(grn.filt, !is.na(.data$gene.ENSEMBL))
      futile.logger::flog.info(paste0("  Number of rows after filtering: ", nrow(grn.filt)))
      
    }
    
    # TODO: Make order more logical
    
    if (!is.null(peak_gene.r_range)) {
      
      futile.logger::flog.info(paste0(" Filter network and retain only rows with peak_gene.r in the following interval: (", 
                                      peak_gene.r_range[1], " - ", peak_gene.r_range[2], "]"))
      
      futile.logger::flog.info(paste0("  Number of rows before filtering: ", nrow(grn.filt)))
      grn.filt = dplyr::filter(grn.filt, 
                               is.na(.data$peak_gene.r) | .data$peak_gene.r  > peak_gene.r_range[1], 
                               is.na(.data$peak_gene.r) | .data$peak_gene.r <= peak_gene.r_range[2])
      futile.logger::flog.info(paste0("  Number of rows after filtering: ", nrow(grn.filt)))
      
    }
    
    
    if (!is.null(peak_gene.p_raw.threshold)) {
      
      futile.logger::flog.info(paste0(" Filter network and retain only rows with peak-gene connections with p.raw < ", peak_gene.p_raw.threshold))
      futile.logger::flog.info(paste0("  Number of rows before filtering TFs: ", nrow(grn.filt)))
      grn.filt = dplyr::filter(grn.filt, is.na(.data$peak_gene.p_raw) | .data$peak_gene.p_raw < peak_gene.p_raw.threshold)
      futile.logger::flog.info(paste0("  Number of rows after filtering TFs: ", nrow(grn.filt)))
      
    }
    
    if (!is.null(peak_gene.fdr.threshold)) {
      
      futile.logger::flog.info(paste0(" Calculate FDR based on remaining rows, filter network and retain only rows with peak-gene connections with an FDR < ",  peak_gene.fdr.threshold))
      
      futile.logger::flog.info(paste0("  Number of rows before filtering genes (including/excluding NA): ", nrow(grn.filt), "/", nrow(grn.filt %>% dplyr::filter(!is.na(.data$peak_gene.p_raw)))))
      
      
      # Adjusted p-value is calculated dynamically here and therefore, for different filters, the numbers may vary
      # After a discussion in the group, this procedure was agreed upon even though it can sometimes yield confusing results when compared among each other
      
      if (peak_gene.fdr.method == "IHW") {
        
        # Identify those entries for which both p-value and covariate are not NA
        
        covariate_val = grn.filt %>% dplyr::pull(!!(peak_gene.IHW.covariate))
        indexes = which(!is.na(grn.filt$peak_gene.p_raw) & !is.na(covariate_val))
        
        if (length(indexes) < nrow(grn.filt)) {
          message = paste0("Could only take ", length(indexes), " rows out of ", nrow(grn.filt), " because some entries for either p-value or covariate were NA. The remaining ", nrow(grn.filt) - length(indexes), " rows will be ignored for p-value adjustment and set to NA.")
          .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
        }
        
        
        if (peak_gene.fdr.method == "IHW" && length(indexes) < 1000) {
          message = paste0("IHW should only be performed with at least 1000 p-values, but only ", length(indexes), " are available. Switching to BH adjustment as fallback. This is to be expected for the permuted data but not for the non-permuted one.")
          
          if (permutationCur == 0) {
            .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
          } else {
            futile.logger::flog.info(message)
          }
          
          peak_gene.fdr.method = "BH"
        }
        
      }
      
      if (peak_gene.fdr.method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")) {
        
        grn.filt = dplyr::mutate(grn.filt, peak_gene.p_adj = stats::p.adjust(.data$peak_gene.p_raw, method = peak_gene.fdr.method))
        
      } else { # Do IHW
        
        suffixFile = .getPermutationSuffixStr(permutationCur)
        
        outputFolder = .checkOutputFolder(GRN, outputFolder)
        outputFile = paste0(outputFolder, .getOutputFileName("plot_peakGene_IHW_diag"), suffixFile, ".pdf")
        
        IHW.res = .performIHW(grn.filt$peak_gene.p_raw[indexes], 
                              covariate_val[indexes] %>% unlist() %>% unname(), 
                              alpha = peak_gene.fdr.threshold, nbins = peak_gene.IHW.nbins, pdfFile = outputFile)
        
        
        grn.filt$peak_gene.p_adj  = NA
        grn.filt$peak_gene.p_adj[indexes] = IHW::adj_pvalues(IHW.res$ihwResults)
        
      }
      
      if (allowMissingGenes) {
        grn.filt = dplyr::filter(grn.filt, is.na(.data$peak_gene.p_adj) | .data$peak_gene.p_adj <  peak_gene.fdr.threshold) # keep NA here due to completeCases variable 
      } else {
        grn.filt = dplyr::filter(grn.filt, .data$peak_gene.p_adj <  peak_gene.fdr.threshold) 
      }
      
      futile.logger::flog.info(paste0("  Number of rows after filtering genes (including/excluding NA): ", nrow(grn.filt), "/", nrow(grn.filt %>% dplyr::filter(!is.na(.data$peak_gene.p_adj)))))
    }
    
    
    if (peak_gene.selection == "closest") {
      
      # Select only the closest gene for each peak
      # Currently, this filter is applied BEFORE any of the other peak-gene filters
      futile.logger::flog.info(paste0(" Filter network and retain only the closest genes for each peak. Note that previous filters may already have eliminated the overall closest gene for a particular peak. To make sure to always use the closest gene in the network, set the other peak_gene filters to NULL."))
      
      # NA distances should be kept, only genes with non NA-values should be filtered
      
      grn.filt = grn.filt %>%
        dplyr::filter(!is.na(.data$peak_gene.distance)) %>%
        dplyr::group_by(.data$peak.ID) %>%
        dplyr::slice(which.min(.data$peak_gene.distance)) %>%
        dplyr::ungroup() %>%
        rbind(dplyr::filter(grn.filt, is.na(.data$peak_gene.distance))) # rbind the na rows separately here
      
      futile.logger::flog.info(paste0("  Number of rows after filtering: ", nrow(grn.filt)))
      
    }
    
    
    grn.filt = grn.filt %>%
      dplyr::select(tidyselect::starts_with("TF."), 
                    tidyselect::starts_with("TF_peak."), 
                    tidyselect::starts_with("peak."), 
                    tidyselect::starts_with("peak_gene."),
                    tidyselect::starts_with("gene."),
                    tidyselect::everything()) %>%
      dplyr::mutate(peak.ID      = as.factor(.data$peak.ID),
                    gene.ENSEMBL = as.factor(.data$gene.ENSEMBL),
                    TF.name      = as.factor(.data$TF.name))
    
    
    GRN@connections$all.filtered[[permIndex]] = grn.filt
    
    futile.logger::flog.info(paste0("Final number of rows left after all filtering steps: ", nrow(grn.filt)))
    
    if (nrow(grn.filt) == 0 & permutationCur == 0 & !silent) {
      message = "No connections passed the filter steps. Rerun the function and be less stringent."
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
    }
    
    
    
    
  } #end for all permutations
  
  
  if (!is.null(GRN@graph) & resetGraphAndStoreInternally) {
    message = "To avoid object inconsistencies and unexpected/non-reproducible results, the graph slot in the object has been reset. For all network-related functions as well as eGRN visualization, rerun the method build_eGRN_graph and all other network-related ans enrichment functions to update to the new set of filtered connections"
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
    GRN@graph = list()
  }
 
  
  
  if (silent) futile.logger::flog.threshold(futile.logger::INFO)
  
  if (!silent) .printExecutionTime(start)
  
  if (!resetGraphAndStoreInternally) {
    return(GRN@connections$all.filtered)
  } else {
    return(GRN)
  }
 
}




.performIHW <- function(pvalues, covariates, alpha = 0.1, covariate_type = "ordinal",
                        nbins = "auto", m_groups = NULL, quiet = TRUE, nfolds = 5L,
                        nfolds_internal = 5L, nsplits_internal = 1L, lambdas = "auto",
                        seed = 1L, distrib_estimator = "grenander", lp_solver = "lpsymphony",
                        adjustment_type = "BH", return_internal = FALSE, doDiagnostics = TRUE, pdfFile = NULL, verbose = TRUE, ...) {
  
  start = Sys.time()
  
  checkmate::assertNumeric(pvalues, lower = 0, upper = 1, any.missing = TRUE)
  checkmate::assert(checkmate::checkNumeric(covariates), checkmate::checkFactor(covariates))
  stopifnot(length(pvalues) == length(covariates))
  checkmate::assertNumber(alpha, lower = 0, upper = 1)
  checkmate::assertChoice(covariate_type, c("ordinal", "nominal"))
  checkmate::assert(checkmate::checkInt(nbins, lower = 1), checkmate::checkSubset(nbins, c("auto")))
  checkmate::assert(checkmate::checkNull(m_groups), checkmate::checkInteger(m_groups, len = length(covariates)))
  checkmate::assertLogical(quiet)
  checkmate::assertInt(nfolds, lower = 1)
  checkmate::assertInt(nfolds_internal, lower = 1)
  checkmate::assertInt(nsplits_internal, lower = 1)
  checkmate::assert(checkmate::checkNumeric(lambdas), checkmate::checkSubset(lambdas, c("auto")))
  checkmate::assert(checkmate::checkNull(seed), checkmate::checkInteger(seed)) 
  checkmate::assertChoice(distrib_estimator, c("grenander", "ECDF"))
  checkmate::assertChoice(lp_solver, c("lpsymphony", "gurobi"))
  checkmate::assertChoice(adjustment_type, c("BH", "bonferroni"))
  checkmate::assertLogical(return_internal)
  checkmate::assertLogical(doDiagnostics)
  checkmate::assert(checkmate::checkNull(pdfFile), checkmate::checkDirectoryExists(dirname(pdfFile), access = "w"))
  checkmate::assertLogical(verbose)
  
  
  res.l = list()
  
  
  futile.logger::flog.info(paste0("Perform IHW based on ", length(pvalues), " p-values"))
  
  ihw_res = IHW::ihw(pvalues, covariates, alpha, covariate_type = covariate_type,
                     nbins = nbins, m_groups = m_groups, quiet = quiet, nfolds = nfolds,
                     nfolds_internal = nfolds_internal, nsplits_internal = nsplits_internal, lambdas = lambdas,
                     seed = seed, distrib_estimator = distrib_estimator, lp_solver = lp_solver,
                     adjustment_type = adjustment_type, return_internal = return_internal, ...)
  
  res.l$ihwResults = ihw_res
  
  if (IHW::nbins(ihw_res) == 1) {
    message = "Only 1 bin, IHW reduces to Benjamini Hochberg (uniform weights). Skipping diagnostic plots"
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
    return(res.l)
  } else {
    futile.logger::flog.info(paste0("  Number of chosen bins (should be >1): ", IHW::nbins(ihw_res)))
  }
  
  
  if (doDiagnostics) {
    
    futile.logger::flog.info("  Generate diagnostic plots for IHW results and data...")
    
    # We can compare this to the result of applying the method of Benjamini and Hochberg to the p-values only:
    
    padj_bh = stats::p.adjust(pvalues, method = "BH")
    rejectionsBH = sum(padj_bh <= alpha, na.rm = TRUE)
    rejectionsIHW = IHW::rejections(ihw_res)
    
    futile.logger::flog.info(paste0("  Number of rejections for IHW: ", 
                                    rejectionsIHW, " (for comparison, number of rejections for BH: ", 
                                    rejectionsBH, " [the latter should be lower])"))
    
    if (rejectionsIHW < rejectionsBH) {
      message = "  Number of rejections for IHW smaller than for BH, something might be wrong. The covariate chosen might not be appropriate."
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
    }
    
    res.l$ihwPlots = list()
    
    ##                                                                           ##
    ## 1. Build in diagnostic functions: Estimated weights and decision boundary ##
    ##                                                                           ##
    
    # Plot No. 1
    res.l$ihwPlots$estimatedWeights  = 
      IHW::plot(ihw_res, what = "weights") + 
      ggplot2::ggtitle("Estimated weights")
    
    # Plot No. 2
    res.l$ihwPlots$decisionBoundary  = 
      IHW::plot(ihw_res, what = "decisionboundary") + 
      ggplot2::ggtitle("Decision boundary")
    
    # Plot No. 3
    res.l$ihwPlots$rawVsAdjPVal_all <- 
      as.data.frame(ihw_res) %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$pvalue, y = .data$adj_pvalue, col = .data$group)) + 
      ggplot2::geom_point(size = 0.25) + ggplot2::scale_colour_hue(l = 70, c = 150, drop = FALSE) + 
      ggplot2::ggtitle("Raw versus adjusted p-values")
    
    # TODO why hard-coded
    pValThreshold = 0.2
    
    # Plot No. 4
    res.l$ihwPlots$rawVsAdjPVal_subset = 
      res.l$ihwPlots$rawVsAdjPVal_all %+% subset(IHW::as.data.frame(ihw_res), .data$adj_pvalue <= pValThreshold) + 
      ggplot2::ggtitle("raw versus adjusted p-values (zoom)")
    
    ##                       ##
    ## 2. p-value histograms ##
    ##                       ##
    
    data.df = data.frame(pValues = pvalues, covariate = covariates)
    
    # One of the most useful diagnostic plots is the p-value histogram (before applying any multiple testing procedure)
    
    # Plot No. 5
    
    res.l$ihwPlots$pValHistogram = 
      ggplot2::ggplot(data.df, ggplot2::aes(x = .data$pValues)) + ggplot2::geom_histogram(binwidth = 0.025, boundary = 0) + 
      ggplot2::ggtitle("p-Value histogram independent of the covariate")
    
    # Stratified p-value histograms by covariate
    # Plot No. 6
    
    nGroups = ihw_res@nbins
    
    data.df$covariate_group <- IHW::groups_by_filter(data.df$covariate, nGroups)
    res.l$ihwPlots$pValHistogramGroupedByCovariate = 
      ggplot2::ggplot(data.df, ggplot2::aes(x=.data$pValues)) + ggplot2::geom_histogram(binwidth = 0.025, boundary = 0) + 
        ggplot2::facet_wrap( ~covariate_group, nrow = 2) + 
      ggplot2::ggtitle("p-Value histogram grouped by the covariate")
    
    # Plot No. 7
    res.l$ihwPlots$pValHistogramGroupedByCovariateECDF = 
      ggplot2::ggplot(data.df, ggplot2::aes(x = .data$pValues, col = .data$covariate_group)) + ggplot2::stat_ecdf(geom = "step")  + 
      ggplot2::ggtitle("p-Value histogram grouped by the covariate (ECDF)")
    
    
    # Check whether the covariate is informative about power under the alternative (property 1), 
    # plot the −log10(p-values) against the ranks of the covariate:
    data.df <- stats::na.omit(data.df)
    data.df$covariateRank = rank(data.df$covariate)/nrow(data.df)
    
    # Plot No. 8
    res.l$ihwPlots$pValAgainstRankCovariate = 
      ggplot2::ggplot(data.df, ggplot2::aes(x= .data$covariateRank, y = -log10(.data$pValues))) + ggplot2::geom_hex(bins = 100) + 
      ggplot2::ggtitle("p-Value against rank of the covariate")
    
    
    if (!is.null(pdfFile)) {
      
      .printMultipleGraphsPerPage(res.l$ihwPlots, nCol = 1, nRow = 1, pdfFile = pdfFile)
      futile.logger::flog.info(paste0("Diagnostic plots written to file ", pdfFile))
    }
    
  }
  
  .printExecutionTime(start)
  return(res.l)
  
  
}

#' Add TF-gene correlations to a \code{\linkS4class{GRN}} object. The information is currently stored in \code{GRN@connections$TF_genes.filtered}. Note that raw p-values are not adjusted.
#' 
#' @export
#' @template GRN
#' @template corMethod
#' @template addRobustRegression
#' @template nCores
#' @template forceRerun
#' @return An updated \code{\linkS4class{GRN}} object, with additional information added from this function.
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' GRN = add_TF_gene_correlation(GRN, forceRerun = FALSE)
add_TF_gene_correlation <- function(GRN, corMethod = "pearson", addRobustRegression = FALSE, nCores = 1, forceRerun = FALSE) {
  
  start = Sys.time() 
  
  checkmate::assertClass(GRN, "GRN")  
  GRN = .addFunctionLogToObject(GRN)    
  
  GRN = .makeObjectCompatible(GRN)
  
  checkmate::assertChoice(corMethod, c("pearson", "spearman"))
  checkmate::assertFlag(addRobustRegression)
  checkmate::assertIntegerish(nCores, lower = 1)
  checkmate::assertFlag(forceRerun)
  
  if (addRobustRegression) {
      packageMessage = paste0("The package robust is not installed, but needed here due to addRobustRegression = TRUE. Please install it and re-run this function or change addRobustRegression to FALSE.")
      .checkPackageInstallation("robust", packageMessage)  
  }
  
  if (is.null(GRN@connections$TF_genes.filtered) | forceRerun) {
    
    GRN@connections$TF_genes.filtered = list()
    
    .checkExistanceFilteredConnections(GRN)
    
    futile.logger::flog.info(paste0("Calculate correlations for TF and genes from the filtered set of connections"))
    
    for (permutationCur in 0:.getMaxPermutation(GRN)) {
      
      futile.logger::flog.info(paste0(" ", .getPermStr(permutationCur)))
      # Get all TF peak pairs to check
      
      permIndex = as.character(permutationCur)
      TF_genePairs = GRN@connections$all.filtered[[permIndex]] %>%
        dplyr::filter(!is.na(.data$gene.ENSEMBL)) %>%
        dplyr::select(.data$TF.name, .data$gene.ENSEMBL) %>%
        dplyr::distinct() %>%
        dplyr::left_join(GRN@annotation$TFs, by = c("TF.name"), suffix = c("", ".transl")) # %>%
      # dplyr::distinct(ENSEMBL, ENSEMBL.transl)
      # TODO: Improve: Only loop over distinct ENSMBL_TF and ENSEMBL_gene pairs
      
      maxRow = nrow(TF_genePairs)
      if ( maxRow > 0) {
        futile.logger::flog.info(paste0("  Iterate through ", maxRow, " TF-gene combinations and (if possible) calculate correlations using ", nCores, " cores. This may take a few minutes."))
        
        
        countsRNA.clean  = getCounts(GRN, type = "rna",  permuted = as.logical(permutationCur), includeIDColumn = FALSE)
        
        map_TF =   match(TF_genePairs$TF.ENSEMBL, getCounts(GRN, type = "rna", permuted = as.logical(permutationCur))$ENSEMBL)
        map_gene = match(TF_genePairs$gene.ENSEMBL, getCounts(GRN, type = "rna",permuted = as.logical(permutationCur))$ENSEMBL)
        
        # Some NAs might be expected, given our annotation contains all known genes
        stopifnot(!all(is.na(map_TF)))
        stopifnot(!all(is.na(map_gene)))
        
        chunksize = 10000
        startIndexMax = ceiling(maxRow / chunksize) - 1 # -1 because we count from 0 onwards
        
        debugMode_nPlots = 0
        if (debugMode_nPlots > 0) {
          nCores = 1
        }
        
        res.l = .execInParallelGen(nCores, returnAsList = TRUE, listNames = NULL, iteration = 0:startIndexMax, verbose = FALSE, functionName = .correlateData, 
                                   chunksize = chunksize, maxRow = maxRow, counts1 = countsRNA.clean, counts2 = countsRNA.clean, 
                                   map1 = map_TF, map2 = map_gene, corMethod = corMethod, debugMode_nPlots = debugMode_nPlots, addRobustRegression = addRobustRegression)
        
        res.m  = do.call(rbind, res.l)
        
        
        selectColumns = c("gene.ENSEMBL", "TF.ENSEMBL", "r", "p.raw", "TF.name")
        if (addRobustRegression) {
          selectColumns = c(selectColumns, "p_raw.robust", "r_robust", "bias_M_p.raw", "bias_LS_p.raw")
        }
        
        
        futile.logger::flog.info(paste0("  Done. Construct the final table, this may result in an increased number of TF-gene pairs due to different TF names linked to the same Ensembl ID."))
        
        
        # Make data frame and adjust p-values
        res.df = suppressMessages(tibble::as_tibble(res.m) %>%
                                    dplyr::mutate(TF.ENSEMBL   = getCounts(GRN, type = "rna", permuted = as.logical(permutationCur))$ENSEMBL[map_TF],
                                                  gene.ENSEMBL = getCounts(GRN, type = "rna", permuted = as.logical(permutationCur))$ENSEMBL[map_gene]) %>%
                                    dplyr::filter(!is.na(.data$gene.ENSEMBL), !is.na(.data$TF.ENSEMBL)) %>%  # For some peak-gene combinations, no RNA-Seq data was available, these NAs are filtered
                                    dplyr::left_join(GRN@annotation$TFs, by = c("TF.ENSEMBL")) %>%
                                    dplyr::select(tidyselect::all_of(selectColumns))) %>%
          dplyr::mutate(gene.ENSEMBL = as.factor(.data$gene.ENSEMBL), 
                        TF.ENSEMBL   = as.factor(.data$TF.ENSEMBL),
                        TF.name           = as.factor(.data$TF.name)) %>%
          dplyr::rename(TF_gene.r     = .data$r, 
                        TF_gene.p_raw = .data$p.raw) %>%
          dplyr::select(.data$TF.name, .data$TF.ENSEMBL, .data$gene.ENSEMBL, tidyselect::everything())
        
        
        if (addRobustRegression) {
          res.df = dplyr::rename(res.df, 
                                 TF_gene.p_raw.robust = .data$p_raw.robust, 
                                 TF_gene.r_robust = .data$r_robust,
                                 TF_gene.bias_M_p.raw = .data$bias_M_p.raw,
                                 TF_gene.bias_LS_p.raw = .data$bias_LS_p.raw)
        }
        
      } else {
        futile.logger::flog.info(paste0(" Nothing to do, skip."))
       
          
        if (addRobustRegression) {
          res.df = tibble::tribble(~TF.name, ~TF.ENSEMBL, ~gene.ENSEMBL, ~TF_gene.r, ~TF_gene.p_raw, ~TF_gene.p_raw.robust, 
                                   ~TF_gene.r_robust, ~TF_gene.bias_M_p.raw, ~TF_gene.bias_LS_p.raw)
        } else {
          res.df = tibble::tribble(~TF.name, ~TF.ENSEMBL, ~gene.ENSEMBL, ~TF_gene.r, ~TF_gene.p_raw)
        }
        
      }
      
      
      
      GRN@connections$TF_genes.filtered[[permIndex]] = res.df
      
    } # end for each permutation
  }
  
  .printExecutionTime(start, prefix = "")
  
  GRN
}




######### SNP functions ################


# Not ready and not tested
addSNPOverlap <- function(grn, SNPData, col_chr = "chr", col_pos = "pos", col_peakID = "peakID", 
                          genomeAssembly_peaks, genomeAssembly_SNP, addAllColumns = TRUE) {
  
  start = Sys.time()
  futile.logger::flog.info(paste0("Adding SNP overlap"))
  
  futile.logger::flog.info(paste0(" Checking validity of input"))
  checkmate::assertDataFrame(SNPData, min.rows = 1)
  checkmate::assertSubset(c(col_chr, col_pos), colnames(SNPData), empty.ok = FALSE)
  checkmate::assertSubset(c(col_peakID), colnames(grn), empty.ok = FALSE)
  
  if (genomeAssembly_peaks != genomeAssembly_SNP) {
    
    message = "Genome assemblies of peaks and SNPs must be identical"
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
  if (!all(grepl("chr", unlist(SNPData[, col_chr])))) {
    SNPData[, col_chr] = paste0("chr", SNPData[, col_chr])
  }
  
  # Check for chr23 = chrX
  if (genomeAssembly_SNP %in% c("hg19", "hg38")) {
    index_chr23 = which(grepl("chr23", SNPData[, col_chr]))
    if (length(index_chr23) > 0) {
      message =paste0(" Chr23 found in SNP data. Renaming to chrX")
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
      SNPData[index_chr23, col_chr] = stringr::str_replace(SNPData[index_chr23, col_chr], "chr23", "chrX")
    }
  }
  
  futile.logger::flog.info(paste0(" Overlapping peaks and SNPs"))
  txdb = .getGenomeObject(genomeAssembly_SNP)
  # GRanges for peaks and SNP
  
  SNPData.mod = SNPData
  SNPData.mod$chr = as.character(SNPData.mod[, col_chr])
  SNPData.mod$end = as.numeric(SNPData.mod[, col_pos])
  SNPData.mod$start = as.numeric(SNPData.mod[, col_pos])
  
  SNPs.gr = .constructGRanges(dplyr::select(SNPData.mod, .data$chr, {{col_pos}}, .data$end), seqlengths = GenomeInfoDb::seqlengths(txdb), genomeAssembly_SNP,  
                              start.field = col_pos, seqnames.field = col_chr)
  
  
  
  peaks.df = .createDF_fromCoordinates(as.character(unlist(unique(grn[,col_peakID]))), col_peakID)
  
  txdb = .getGenomeObject(genomeAssembly_peaks)
  peaks.gr = .constructGRanges(peaks.df, seqlengths = GenomeInfoDb::seqlengths(txdb), genomeAssembly_peaks)
  
  # Overlap
  overlapsAll = GenomicRanges::findOverlaps(peaks.gr, SNPs.gr, 
                                            minoverlap = 1,
                                            type = "any",
                                            select = "all",
                                            ignore.strand = TRUE)
  
  futile.logger::flog.info(paste0(" Summarizing overlap"))
  
  query_row_ids   = S4Vectors::queryHits(overlapsAll)
  subject_row_ids = S4Vectors::subjectHits(overlapsAll)
  
  query_overlap_df     = as.data.frame(S4Vectors::elementMetadata(peaks.gr)[query_row_ids, col_peakID])
  subject_overlap_df   = as.data.frame( SNPs.gr)[subject_row_ids, c("seqnames", "start")]
  
  overlaps.df = cbind.data.frame(query_overlap_df,subject_overlap_df) %>% dplyr::mutate(seqnames = as.character(.data$seqnames))
  colnames(overlaps.df) = c("peakID", "SNP_chr", "SNP_start")
  
  if (addAllColumns) {
    
    overlaps.df = dplyr::left_join(overlaps.df, SNPData.mod, by = c("SNP_chr" = "chr", "SNP_start" = "start")) %>% dplyr::select(-end)
  }
  
  
  myfun <- function(x) {
    paste0(x,collapse = "|")
  }
  
  colnamesNew = setdiff(colnames(overlaps.df), col_peakID)
  overlaps.summarized.df =  overlaps.df %>%
    dplyr::group_by(.data$peak) %>%
    dplyr::summarise_at(colnamesNew, myfun) %>%
    dplyr::ungroup()
  
  # Add no of SNPs also as column
  overlaps.summarized.df = overlaps.summarized.df %>%
    dplyr::mutate(SNP_nOverlap = stringr::str_count(.data$SNP_chr, stringr::fixed("|")) + 1)
  
  grn.SNPs = dplyr::left_join(grn, overlaps.summarized.df, by = col_peakID)
  
  .printExecutionTime(start)
  
  grn.SNPs
  
}

# TODO: Only called for SNPs
.createDF_fromCoordinates <- function (IDs, colnameID = "peakID") {
  
  splits = strsplit(as.character(IDs), split=c(":|-"), fixed = FALSE,  perl = TRUE)
  df = tibble::tibble(      chr   = sapply(splits,"[[", 1),
                            start = as.numeric(sapply(splits,"[[", 2)),
                            end   = as.numeric(sapply(splits,"[[", 3))) %>%
    dplyr::mutate({{colnameID}} := paste0(.data$chr, ":", .data$start, "-",.data$end))
  
  df
}


####### STATS #########


#' Generate a summary for the number of connections for different filtering criteria for a \code{\linkS4class{GRN}} object. 
#' 
#' This functions calls \code{\link{filterGRNAndConnectGenes}} repeatedly and stores the total number of connections and other statistics each time to summarize them afterwards. 
#' All arguments are identical to the ones in \code{\link{filterGRNAndConnectGenes}}, see the help for this function for details.
#' The function \code{\link{plot_stats_connectionSummary}} can be used afterwards for plotting.
#' 
#' @export
#' @template GRN 
#' @param TF_peak.fdr Numeric vector[0,1]. Default \code{c(0.001, 0.01, 0.05, 0.1, 0.2)}. TF-peak FDR values to iterate over.
#' @template TF_peak.connectionTypes
#' @param peak_gene.fdr Numeric vector[0,1]. Default \code{c(0.001, 0.01, 0.05, 0.1, 0.2)}. Peak-gene FDR values to iterate over.
# #' @param peak_gene.p_raw  Numeric vector[0,1]. Default \code{NULL}. Peak-gene raw p-value values to iterate over. Skipped if set to NULL.
#' @param peak_gene.r_range Numeric vector of length 2[-1,1]. Default \code{c(0,1)}. The correlation range of peak-gene connections to keep.
#' @template gene.types
#' @param allowMissingGenes Logical vector of length 1 or 2. Default \code{c(FALSE, TRUE)}. Allow genes to be missing for peak-gene connections? If both \code{FALSE} and \code{TRUE} are given, the code loops over both
#' @param allowMissingTFs Logical vector of length 1 or 2. Default \code{c(FALSE)}. Allow TFs to be missing for TF-peak connections?  If both \code{FALSE} and \code{TRUE} are given, the code loops over both
#' @template forceRerun
#' @seealso \code{\link{plot_stats_connectionSummary}}
#' @seealso \code{\link{filterGRNAndConnectGenes}}
#' @return An updated \code{\linkS4class{GRN}} object, with additional information added from this function.
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' GRN = generateStatsSummary(GRN, TF_peak.fdr = c(0.01, 0.1), peak_gene.fdr = c(0.01, 0.1))
#' 
generateStatsSummary <- function(GRN, 
                                 TF_peak.fdr = c(0.001, 0.01, 0.05, 0.1, 0.2),
                                 TF_peak.connectionTypes = "all",
                                 peak_gene.fdr = c(0.001, 0.01, 0.05, 0.1, 0.2),
#                                 peak_gene.p_raw = NULL,
                                 peak_gene.r_range = c(0,1),
                                 gene.types = c("protein_coding"),
                                 allowMissingGenes = c(FALSE, TRUE),
                                 allowMissingTFs = c(FALSE),
                                 forceRerun = FALSE) {
  
  start = Sys.time()   
  
  checkmate::assertClass(GRN, "GRN")
  GRN = .addFunctionLogToObject(GRN)
  
  GRN = .makeObjectCompatible(GRN)

  checkmate::assertNumeric(TF_peak.fdr, lower = 0, upper = 1, min.len = 1)
  checkmate::assertSubset(TF_peak.connectionTypes, c("all", GRN@config$TF_peak_connectionTypes), empty.ok = FALSE)
  checkmate::assertNumeric(peak_gene.fdr, lower = 0, upper = 1, min.len = 1)
  #   checkmate::assert(checkmate::checkNull(peak_gene.p_raw), checkmate::checkNumeric(peak_gene.p_raw, lower = 0, upper = 1, min.len = 1))
  checkmate::assertNumeric(peak_gene.r_range, lower = -1, upper = 1, len = 2)
  checkmate::assertSubset(gene.types, c("all", unique(as.character(GRN@annotation$genes$gene.type))) %>% stats::na.omit(), empty.ok = FALSE)
  checkmate::assertSubset(allowMissingGenes, c(TRUE, FALSE))
  checkmate::assertSubset(allowMissingTFs, c(TRUE, FALSE))
  checkmate::assertFlag(forceRerun)
  
  if (is.null(GRN@stats$connections) | is.null(GRN@stats$connectionDetails.l) | forceRerun) {
    
    GRN@stats$connections = .initializeStatsDF()
    
    if (TF_peak.connectionTypes == "all") {
      TF_peak.connectionTypesAllComb =.getAll_TF_peak_connectionTypes(GRN)
    } else {
      TF_peak.connectionTypesAllComb = unique(TF_peak.connectionTypes)
    }
    
    
    futile.logger::flog.info(paste0("Generating summary. This may take a while..."))
    
    
    for (permutationCur in 0:.getMaxPermutation(GRN)) {
      
      futile.logger::flog.info(paste0("\n", .getPermStr(permutationCur), "...\n"))
      
      permIndex = as.character(permutationCur)
      GRN@stats$connectionDetails.l[[permIndex]] = list()
      
      # Iterate over different stringency thresholds and collect statistics
      
      for (TF_peak.fdr_cur in TF_peak.fdr) {
        
        TF_peak.fdr_cur.str = as.character(TF_peak.fdr_cur)
        
        futile.logger::flog.info(paste0("Calculate network stats for TF-peak FDR of ", TF_peak.fdr_cur))
        GRN@stats$connectionDetails.l[[permIndex]] [[TF_peak.fdr_cur.str]] =list()
        
        futile.logger::flog.debug(paste0("Iterating over different peak-gene FDR thresholds..."))
        for (peak_gene.fdr_cur in peak_gene.fdr) {
          
          futile.logger::flog.debug(paste0("Peak-gene FDR = ", peak_gene.fdr_cur))
          peak_gene.fdr_cur.str = as.character(peak_gene.fdr_cur)
          GRN@stats$connectionDetails.l[[permIndex]] [[TF_peak.fdr_cur.str]] [[peak_gene.fdr_cur.str]] =list()
          
          for (allowMissingTFsCur in allowMissingTFs) {
            
            futile.logger::flog.debug(paste0("  allowMissingTFs = ", allowMissingTFsCur))
            allowMissingTFsCur.str = as.character(allowMissingTFsCur)
            GRN@stats$connectionDetails.l[[permIndex]] [[TF_peak.fdr_cur.str]] [[peak_gene.fdr_cur.str]] [[allowMissingTFsCur.str]] = list()
            
            for (allowMissingGenesCur in allowMissingGenes) {
              
              futile.logger::flog.debug(paste0("  allowMissingGenes = ", allowMissingGenesCur))
              allowMissingGenesCur.str = as.character(allowMissingGenesCur)
              GRN@stats$connectionDetails.l[[permIndex]] [[TF_peak.fdr_cur.str]] [[peak_gene.fdr_cur.str]] [[allowMissingTFsCur.str]] [[allowMissingGenesCur.str]] = list()
              
              for (TF_peak.connectionTypeCur in TF_peak.connectionTypesAllComb) {
                
                GRN@stats$connectionDetails.l[[permIndex]] [[TF_peak.fdr_cur.str]] [[peak_gene.fdr_cur.str]] [[allowMissingTFsCur.str]] [[allowMissingGenesCur.str]] [[TF_peak.connectionTypeCur]] = list()
                
                futile.logger::flog.debug(paste0("    TF_peak.connectionType = ", TF_peak.connectionTypeCur))
                
                futile.logger::flog.threshold(futile.logger::WARN)
                con.filt.l = filterGRNAndConnectGenes(GRN, 
                                               TF_peak.fdr.threshold = TF_peak.fdr_cur, 
                                               TF_peak.connectionTypes = TF_peak.connectionTypeCur, 
                                               peak_gene.p_raw.threshold = NULL, 
                                               peak_gene.fdr.threshold= peak_gene.fdr_cur,
                                               gene.types = gene.types, 
                                               allowMissingGenes = allowMissingGenesCur, 
                                               allowMissingTFs = allowMissingTFsCur,
                                               peak_gene.r_range = peak_gene.r_range,
                                               filterTFs = NULL, filterGenes = NULL, filterPeaks = NULL,
                                               resetGraphAndStoreInternally = FALSE,
                                               silent = TRUE)
                
                futile.logger::flog.threshold(futile.logger::INFO)
                
                results.l = .addStats(GRN@stats$connections, con.filt.l[[permIndex]], 
                                      perm = permutationCur, 
                                      TF_peak.fdr = TF_peak.fdr_cur, TF_peak.connectionType = TF_peak.connectionTypeCur,
                                      peak_gene.p_raw = NA,
                                      peak_gene.fdr = peak_gene.fdr_cur, 
                                      peak_gene.r_range = paste0(peak_gene.r_range, collapse = ","),
                                      gene.types = paste0(gene.types, collapse = ","),
                                      allowMissingGenes = allowMissingGenesCur, 
                                      allowMissingTFs   = allowMissingTFsCur)
                
                GRN@stats$connections  = results.l[["summary"]]
                
                GRN@stats$connectionDetails.l[[permIndex]] [[TF_peak.fdr_cur.str]] [[peak_gene.fdr_cur.str]] [[allowMissingTFsCur.str]] [[allowMissingGenesCur.str]] [[TF_peak.connectionTypeCur]] = 
                  results.l[["details"]]
                
              } # end of  for (TF_peak.connectionTypeCur in TF_peak.connectionTypesAllComb)
              
            } # end of  for (allowMissingGenesCur in allowMissingGenes) 
            
          } # end of  for (allowMissingTFsCur in allowMissingTFs)
          
        } # end of for (peak_gene.fdr_cur in peak_gene.fdr)
        
        
        # REMOVED FOR NOW, WAS INCOMPLETE AND NOT NEEDED ANYWAY
        # if (!is.null(peak_gene.p_raw)) {
        #   
        #   futile.logger::flog.info(paste0(" Iterating over different peak-gene raw p-value thresholds..."))
        #   for (peak_gene.p_raw_cur in peak_gene.p_raw) {
        #     
        #     futile.logger::flog.info(paste0("  Peak-gene raw p-value = ", peak_gene.p_raw_cur))
        #     
        #     # TODO: Add the other ones here also
        #     for (allowMissingGenesCur in allowMissingGenes) {
        #       
        #       GRN = filterGRNAndConnectGenes(GRN, 
        #                                      TF_peak.fdr.threshold = TF_peak.fdr_cur, peak_gene.p_raw.threshold = peak_gene.p_raw_cur, 
        #                                      peak_gene.fdr.threshold= NULL,
        #                                      gene.types = gene.types, 
        #                                      allowMissingGenes = allowMissingGenesCur, peak_gene.r_range = peak_gene.r_range,
        #                                      filterTFs = NULL, filterGenes = NULL, filterPeaks = NULL,
        #                                      silent = TRUE)
        #       
        #       
        #       GRN@stats$connections = .addStats(GRN@stats$connections, GRN@connections$all.filtered[[permIndex]], 
        #                                         perm = permutationCur, 
        #                                         TF_peak.fdr = TF_peak.fdr_cur, TF_peak.connectionType = TF_peak.connectionTypeCur,
        #                                         peak_gene.fdr = NA, 
        #                                         peak_gene.p_raw = peak_gene.p_raw_cur,
        #                                         peak_gene.r_range = paste0(peak_gene.r_range, collapse = ","),
        #                                         gene.types = paste0(gene.types, collapse = ","),
        #                                         allowMissingGenes = allowMissingGenesCur,
        #                                         allowMissingTFs   = allowMissingTFsCur)
        #       
        #     } # end of for (allowMissingGenesCur in allowMissingGenes)
        #     
        #   } # end of   for (peak_gene.p_raw_cur in peak_gene.p_raw) 
        # } # end of  if (!is.null(peak_gene.p_raw))
        
        
      }
    } # end for each permutation
    
  } else {
    futile.logger::flog.info(paste0("Data already exists in object. Set forceRerun = TRUE to regenerate and overwrite."))
  }
  
  .printExecutionTime(start, prefix = "")
  
  GRN
  
}


.addStats <- function(stats.df, connections.df, perm, 
                      TF_peak.fdr, TF_peak.connectionType, 
                      peak_gene.p_raw, peak_gene.fdr, peak_gene.r_range, 
                      gene.types,
                      allowMissingGenes, allowMissingTFs) {
  
  TF.stats   = dplyr::select(connections.df, .data$TF.name, .data$peak.ID)   %>% dplyr::filter(!is.na(.data$peak.ID)) %>% dplyr::pull(.data$TF.name)   %>% as.character() %>% table() 
  gene.stats = dplyr::select(connections.df, .data$peak.ID, .data$gene.ENSEMBL) %>% dplyr::filter(!is.na(.data$gene.ENSEMBL)) %>% dplyr::pull(.data$gene.ENSEMBL) %>% as.character() %>% table() 
  
  peak_gene.stats = dplyr::select(connections.df, .data$peak.ID, .data$gene.ENSEMBL) %>% dplyr::filter(!is.na(.data$gene.ENSEMBL),!is.na(.data$peak.ID)) %>% dplyr::pull(.data$peak.ID) %>% as.character() %>% table() 
  peak.TF.stats   = dplyr::select(connections.df, .data$peak.ID, .data$TF.name)      %>% dplyr::filter(!is.na(.data$TF.name),     !is.na(.data$peak.ID)) %>% dplyr::pull(.data$peak.ID) %>% as.character() %>% table() 
  
  if (length(TF.stats) > 0){
    TF.connections = c(min(TF.stats, na.rm = TRUE), 
                       mean(TF.stats, na.rm = TRUE), 
                       median(TF.stats, na.rm = TRUE), 
                       max(TF.stats, na.rm = TRUE))
  } else {
    TF.connections = rep(0,4)
  }
  
  if (length(gene.stats) > 0){
    gene.connections = c(min(gene.stats, na.rm = TRUE), 
                         mean(gene.stats, na.rm = TRUE), 
                         median(gene.stats, na.rm = TRUE), 
                         max(gene.stats, na.rm = TRUE))
  } else {
    gene.connections = rep(0,4)
  }
  
  if (length(peak_gene.stats) > 0){
    peak_gene.connections = c(min(peak_gene.stats, na.rm = TRUE), 
                              mean(peak_gene.stats, na.rm = TRUE), 
                              median(peak_gene.stats, na.rm = TRUE), 
                              max(peak_gene.stats, na.rm = TRUE))
  } else {
    peak_gene.connections = rep(0,4)
  }
  
  if (length(peak.TF.stats) > 0){
    peak_TF.connections = c(min(peak.TF.stats, na.rm = TRUE), 
                            mean(peak.TF.stats, na.rm = TRUE), 
                            median(peak.TF.stats, na.rm = TRUE), 
                            max(peak.TF.stats, na.rm = TRUE))
  } else {
    peak_TF.connections = rep(0,4)
  }
  
  stats.df = tibble::add_row(stats.df, 
                             perm = perm, 
                             TF_peak.fdr = TF_peak.fdr,
                             TF_peak.connectionType = TF_peak.connectionType,
                             peak_gene.p_raw = peak_gene.p_raw,
                             peak_gene.fdr = peak_gene.fdr,
                             peak_gene.r_range = peak_gene.r_range,
                             gene.types = gene.types,
                             allowMissingGenes = allowMissingGenes,
                             allowMissingTFs = allowMissingTFs,
                             
                             nGenes = length(unique(connections.df$gene.ENSEMBL)),
                             nPeaks = length(unique(connections.df$peak.ID)),
                             nTFs   = length(unique(connections.df$TF.name)),
                             
                             TF.connections_min           = TF.connections[1],
                             TF.connections_mean          = TF.connections[2],
                             TF.connections_median        = TF.connections[3],
                             TF.connections_max           = TF.connections[4],
                             
                             peak_TF.connections_min      = peak_TF.connections[1],
                             peak_TF.connections_mean     = peak_TF.connections[2],
                             peak_TF.connections_median   = peak_TF.connections[3],
                             peak_TF.connections_max      = peak_TF.connections[4],
                             
                             peak_gene.connections_min    = peak_gene.connections[1],
                             peak_gene.connections_mean   = peak_gene.connections[2],
                             peak_gene.connections_median = peak_gene.connections[3],
                             peak_gene.connections_max    = peak_gene.connections[4],
                             
                             gene.connections_min         = gene.connections[1],
                             gene.connections_mean        = gene.connections[2],
                             gene.connections_median      = gene.connections[3],
                             gene.connections_max         = gene.connections[4],
  )
  
  list(summary = stats.df, details = list(TF = TF.stats, gene = gene.stats, peak_gene = peak_gene.stats, peak.TF = peak.TF.stats))
}


.initializeStatsDF <- function(){
  
  tibble::tribble(~perm, 
                  ~TF_peak.fdr, ~TF_peak.connectionType,
                  ~peak_gene.p_raw, ~peak_gene.fdr, ~peak_gene.r_range, 
                  ~gene.types,
                  ~allowMissingGenes, ~allowMissingTFs,
                  ~nGenes, ~nPeaks, ~nTFs,
                  ~TF.connections_min, ~TF.connections_mean, ~TF.connections_median, ~TF.connections_max,
                  ~peak_TF.connections_min, ~peak_TF.connections_mean, ~peak_TF.connections_median,  ~peak_TF.connections_max,
                  ~peak_gene.connections_min, ~peak_gene.connections_mean, ~peak_gene.connections_median, ~peak_gene.connections_max,
                  ~gene.connections_min, ~gene.connections_mean, ~gene.connections_median,  ~gene.connections_max)
  
}



####### Misc functions #########

#' Load example GRN dataset
#' 
#' Loads an example GRN object with 6 TFs, ~61.000 peaks, ~19.000 genes, 259 filtered connections and pre-calculated enrichments from the internet. 
#' This function uses \code{BiocFileCache} if installed to cache the example object, which is 
#' considerably faster than re-downloading the file anew every time the function is executed.
#' If not, the file is re-downloaded every time anew. Thus, to enable caching, you may install the package \code{BiocFileCache}.
#' 
#' @export
#' @param forceDownload \code{TRUE} or \code{FALSE}. Default \code{FALSE}. Should the download be enforced even if the local cached file is already present?
#' @param fileURL Character. Default \url{https://www.embl.de/download/zaugg/GRaNIE/GRN.rds}. URL to the GRN example object in rds format.
#' @examples 
#' GRN = loadExampleObject()
#' @return An small example \code{\linkS4class{GRN}} object
loadExampleObject <- function(forceDownload = FALSE, fileURL = "https://www.embl.de/download/zaugg/GRaNIE/GRN.rds") {
  
  checkmate::assertFlag(forceDownload)
  options(timeout=200)
    
  if (!is.installed("BiocFileCache")) {
      
      message = paste0("The package BiocFileCache is not installed, but recommended if you want to speed-up the retrieval of the GRN example object ",
      "via this function when using it multiple times. If not installed, the example object has to be downloaded anew every time you use this function.")
      .checkPackageInstallation("BiocFileCache", message, isWarning = TRUE)
      
      GRN = readRDS(url(fileURL))

  } else {
      
      
      # Taken and modified from https://www.bioconductor.org/packages/release/bioc/vignettes/BiocFileCache/inst/doc/BiocFileCache.html
      
      bfc <- .get_cache()
      
      rid <- BiocFileCache::bfcquery(bfc, "GRaNIE_object_example")$rid
      if (!length(rid)) {
          rid <- names(BiocFileCache::bfcadd(bfc, "GRaNIE_object_example", fileURL))
      }
      if (!isFALSE(BiocFileCache::bfcneedsupdate(bfc, rid)) | forceDownload) {
          messageStr = paste0("Downloading GRaNIE example object from ", fileURL)
          message(messageStr)
          filePath = BiocFileCache::bfcdownload(bfc, rid, ask = FALSE)
      }
      
      
      filePath = BiocFileCache::bfcrpath(bfc, rids = rid)
      
      # Now we can read in the locally stored file
      GRN = readRDS(filePath)
  }

  
  # Change the default path to the current working directory
  GRN@config$directories$outputRoot = "."
  GRN@config$directories$output_plots = "."
  GRN@config$directories$motifFolder = "."
  GRN@config$files$output_log = "GRN.log"
  
  GRN = .makeObjectCompatible(GRN)
  
  # TODO remove once done already in the object
  GRN = add_TF_gene_correlation(GRN)
  
  # GRN = add_featureVariation(GRN, formula = "auto", metadata = "mt_frac")
  
  GRN
  
}


#' Get counts for the various data defined in a \code{\linkS4class{GRN}} object
#' 
#' Get counts for the various data defined in a \code{\linkS4class{GRN}} object.
#' \strong{Note: This function, as all \code{get} functions from this package, does NOT return a \code{\linkS4class{GRN}} object.}
#' 
#' @template GRN 
#' @param type Character. Either \code{peaks} or \code{rna}. \code{peaks} corresponds to the counts for the open chromatin data, while \code{rna} refers to th RNA-seq counts. If set to \code{rna}, both permuted and non-permuted data can be retrieved, while for \code{peaks}, only the non-permuted one (i.e., 0) can be retrieved.
#' @param asMatrix Logical. \code{TRUE} or \code{FALSE}. Default \code{FALSE}. If set to \code{FALSE}, counts are returned as a data frame with or without an ID column (see \code{includeIDColumn}). If set to \code{TRUE}, counts are returned as a matrix with the ID column as row names.
#' @param includeIDColumn Logical. \code{TRUE} or \code{FALSE}. Default \code{TRUE}. Only relevant if \code{asMatrix = FALSE}. If set to \code{TRUE}, an explicit ID column is returned (no row names). If set to \code{FALSE}, the IDs are in the row names instead.
#' @param includeFiltered  Logical. \code{TRUE} or \code{FALSE}. Default \code{FALSE}. If set to \code{FALSE}, genes or peaks marked as filtered (after running the function \code{filterData}) will not be returned. If set to \code{TRUE}, all elements are returned regardless of the currently active filter status.
#' @template permuted
#' @export
#' @import tibble
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' counts.df = getCounts(GRN, type = "peaks", permuted = FALSE)
#' @return Data frame of counts, with the type as indicated by the function parameters. This function does **NOT** return a \code{\linkS4class{GRN}} object.
getCounts <- function(GRN, type,  permuted = FALSE, asMatrix = FALSE, includeIDColumn = TRUE, includeFiltered = FALSE) {
  
  checkmate::assertClass(GRN, "GRN")
  GRN = .addFunctionLogToObject(GRN)     
  
  GRN = .makeObjectCompatible(GRN)
  
  checkmate::assertChoice(type, c("peaks", "rna"))
  checkmate::assertFlag(permuted)
  
  permIndex = dplyr::if_else(permuted, "1", "0")
  
  if (type == "peaks") {
    
    if (permuted) {
      message = "Could not find permuted peak counts in GRN object. Peaks are not stored as permuted, set permuted = FALSE for type = \"peaks\""
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    classPeaks = class(GRN@data$peaks$counts)
    if ("matrix" %in% classPeaks) {
      result = GRN@data$peaks$counts
    } else if ("dgCMatrix" %in% classPeaks) {
      result = .asMatrixFromSparse(GRN@data$peaks$counts, convertZero_to_NA = FALSE)
    } else {
      message = paste0("Unsupported class for GRN@data$peaks$counts. Contact the authors.")
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    
  } else if (type == "rna") {
    
    classRNA = class(GRN@data$RNA$counts)
    if ("matrix" %in% classRNA) {
      result = GRN@data$RNA$counts
    } else if ("dgCMatrix" %in% classRNA) {
      result = .asMatrixFromSparse(GRN@data$RNA$counts, convertZero_to_NA = FALSE)
    } else {
      message = paste0("Unsupported class for GRN@data$RNA$counts. Contact the authors.")
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
  }
  
  
  if (permuted) {
    
    if (type == "rna") {
      
      # Columns are shuffled so that non-matching samples are compared throughout the pipeline for all correlation-based analyses
      colnames(result) = colnames(result)[GRN@data$RNA$counts_permuted_index]
    }
  } 
  
  if (!includeFiltered) {
    
    if (type == "rna") {
      nonFiltered = which(! GRN@data$RNA$counts_metadata$isFiltered)
    } else {
      nonFiltered = which(! GRN@data$peaks$counts_metadata$isFiltered)
    }
    
    result = result[nonFiltered,]
  }

  
  
  if (!asMatrix) {
 
    result.df =  result %>%
      as.data.frame()
    
    if (includeIDColumn)  {
      
      ID_column = dplyr::if_else(type == "rna", "ENSEMBL", "peakID")
      
      result.df = result.df %>%
        tibble::rownames_to_column(ID_column) 
    } 
    
    return(result.df)

  }
  
  result
}


#' Extract connections or links from a \code{\linkS4class{GRN}} object as da data frame.
#' 
#' Returns stored connections/links (either TF-peak, peak-genes or the filtered set of connections as produced by \code{\link{filterGRNAndConnectGenes}}). 
#' \strong{Note: This function, as all \code{get} functions from this package, does NOT return a \code{\linkS4class{GRN}} object.}
#' 
#' @export
#' @template GRN 
#' @template permuted
#' @param type Character. One of \code{TF_peaks}, \code{peak_genes}, \code{TF_genes} or \code{all.filtered}. Default \code{all.filtered}. The type of connections to retrieve.
#' @param include_TF_gene_correlations Logical. \code{TRUE} or \code{FALSE}. Default \code{FALSE}. Should TFs and gene correlations be returned as well? If set to \code{TRUE}, they must have been computed beforehand with \code{\link{add_TF_gene_correlation}}. Only relevant for type = "all.filtered"
#' @param include_TFMetadata Logical. \code{TRUE} or \code{FALSE}. Default \code{FALSE}. Should TF metadata be returned as well? Only relevant for type = "all.filtered"
#' @param include_peakMetadata Logical. \code{TRUE} or \code{FALSE}. Default \code{FALSE}. Should peak metadata be returned as well?  Only relevant for type = "all.filtered"
#' @param include_geneMetadata Logical. \code{TRUE} or \code{FALSE}. Default \code{FALSE}. Should gene metadata be returned as well?  Only relevant for type = "all.filtered"
#' @param include_variancePartitionResults Logical. \code{TRUE} or \code{FALSE}. Default \code{FALSE}. 
#' Should the results from the function \code{\link{add_featureVariation}} be included? 
#' If set to \code{TRUE}, they must have been computed beforehand with \code{\link{add_featureVariation}}; otherwise, an error is thrown.
#' Only relevant for type = "all.filtered"
#' @return A data frame with the requested connections. This function does **NOT** return a \code{\linkS4class{GRN}} object.
#' @seealso \code{\link{filterGRNAndConnectGenes}}
#' @seealso \code{\link{add_featureVariation}}
#' @seealso \code{\link{add_TF_gene_correlation}}
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' GRN_con.all.df = getGRNConnections(GRN)
getGRNConnections <- function(GRN, type = "all.filtered",  permuted = FALSE, 
                              include_TF_gene_correlations = FALSE, 
                              include_TFMetadata = FALSE,
                              include_peakMetadata = FALSE,
                              include_geneMetadata = FALSE,
                              include_variancePartitionResults = FALSE) {
  
  checkmate::assertClass(GRN, "GRN")  
  GRN = .addFunctionLogToObject(GRN)
  
  GRN = .makeObjectCompatible(GRN)
   
  checkmate::assertChoice(type, c("TF_peaks", "peak_genes", "TF_genes", "all.filtered"))
  checkmate::assertFlag(permuted)
  #checkmate::assertIntegerish(permutation, lower = 0, upper = .getMaxPermutation(GRN))
  checkmate::assertFlag(include_TF_gene_correlations)
  checkmate::assertFlag(include_variancePartitionResults)
  checkmate::assertFlag(include_TFMetadata)
  checkmate::assertFlag(include_peakMetadata)
  checkmate::assertFlag(include_geneMetadata)
  
  permIndex = dplyr::if_else(permuted, "1", "0")
  
  if (type == "all.filtered") {
      
    merged.df = GRN@connections$all.filtered[[permIndex]]

    if (include_TF_gene_correlations) {
        
        if (is.null(GRN@connections$TF_genes.filtered)) {
            message = "Please run the function add_TF_gene_correlation first. "
            .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
        }
        
        # Merge with TF-gene table
        merged.df = merged.df %>%
            dplyr::left_join(GRN@connections$TF_genes.filtered[[permIndex]], 
                             by = c("TF.name", "TF.ENSEMBL", "gene.ENSEMBL")) 
    }
    
    if (include_variancePartitionResults) {
        
        if (ncol(GRN@annotation$genes %>% dplyr::select(tidyselect::starts_with("variancePartition"))) == 0 |
            ncol(GRN@annotation$peaks %>% dplyr::select(tidyselect::starts_with("variancePartition"))) == 0) {
            message = "Please run the function add_featureVariation first. "
            .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
        }
        
        merged.df = merged.df %>%
            dplyr::left_join(GRN@annotation$TFs %>% 
                                 dplyr::select(.data$TF.ENSEMBL, tidyselect::starts_with("TF.variancePartition")), 
                             by = "TF.ENSEMBL")  %>%
            dplyr::left_join(GRN@annotation$genes %>% 
                                 dplyr::select(.data$gene.ENSEMBL, tidyselect::starts_with("gene.variancePartition")), 
                             by = "gene.ENSEMBL")  %>%
            dplyr::left_join(GRN@annotation$peaks %>% 
                                 dplyr::select(.data$peak.ID, tidyselect::starts_with("peak.variancePartition")), 
                             by = "peak.ID")
        
    }
    
    if (include_TFMetadata) {
        
        colsMissing = setdiff(colnames(GRN@annotation$TFs), colnames(merged.df))
        if (length(colsMissing) > 0) {
            merged.df = merged.df %>%
                dplyr::left_join(GRN@annotation$TFs, by = c("TF.name", "TF.ENSEMBL"))
        }
    }
    
    if (include_peakMetadata) {
        
        colsMissing = setdiff(colnames(GRN@annotation$peaks), colnames(merged.df))
        if (length(colsMissing) > 0) {
            merged.df = merged.df %>%
                dplyr::left_join(GRN@annotation$peaks, by = "peak.ID")
        }
    }
    
    if (include_geneMetadata) {
        
        colsMissing = setdiff(colnames(GRN@annotation$genes), colnames(merged.df))
        if (length(colsMissing) > 0) {
            merged.df = merged.df %>%
                dplyr::left_join(GRN@annotation$genes, by = "gene.ENSEMBL")
        }
    }

    
    merged.df = merged.df %>%
        dplyr::select(tidyselect::starts_with("TF."), 
                      tidyselect::starts_with("peak."), 
                      tidyselect::starts_with("TF_peak."), 
                      tidyselect::starts_with("gene."), 
                      tidyselect::starts_with("peak_gene."), 
                      tidyselect::starts_with("TF_gene."), 
                      tidyselect::everything()) %>%
        tibble::as_tibble()

    return(merged.df)
    
  } else if (type == "TF_peaks") {
    
    return(tibble::as_tibble(GRN@connections$TF_peaks[[permIndex]]$main))
    
  } else if (type == "peak_genes") {
    
    return(tibble::as_tibble(GRN@connections$peak_genes[[permIndex]]))
    
  } else if (type == "TF_genes") {
      
    if (is.null(GRN@connections$TF_genes.filtered)) {
      message = "Please run the function add_TF_gene_correlation first. "
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
      
    return(tibble::as_tibble(GRN@connections$TF_genes.filtered[[permIndex]]))
      
  } 
  
}

.getAll_TF_peak_connectionTypes <- function (GRN) {
  
  TF_peak.connectionTypesAll = unique(as.character(GRN@connections$TF_peaks[["0"]]$main$TF_peak.connectionType))
  if (length(TF_peak.connectionTypesAll) > 1) {
    # Dont use all, always separate between the different connection types
    # TF_peak.connectionTypesAll = c(TF_peak.connectionTypesAll, "all")
  }
  TF_peak.connectionTypesAll 
}


.checkOutputFolder <- function (GRN, outputFolder) {
  
  if (!is.null(outputFolder)) {
    
    checkmate::assertDirectory(outputFolder, access = "w")
    
    if (!endsWith(outputFolder, .Platform$file.sep)) {
      outputFolder = paste0(outputFolder, .Platform$file.sep)
    }
    
    if (!dir.exists(outputFolder)) {
      dir.create(outputFolder)
    }
  } else {
    
    #  Re-create the output folder here and adjust to the OS-specific path separator, do not rely on what is stored in the object
    if (.Platform$OS.type == "windows") {
      GRN@config$directories$output_plots = gsub('/', ('\\'), GRN@config$directories$output_plots, fixed = TRUE)
    } else {
      GRN@config$directories$output_plots = gsub("\\", "/", GRN@config$directories$output_plots, fixed = TRUE)
    }
      
    if (!endsWith(GRN@config$directories$output_plots, .Platform$file.sep)) {
        GRN@config$directories$output_plots = paste0(GRN@config$directories$output_plots, .Platform$file.sep)
    }
    
    if (!dir.exists(GRN@config$directories$output_plots)) {
      dir.create(GRN@config$directories$output_plots, recursive = TRUE)
    }
  }
  
  dplyr::if_else(is.null(outputFolder), GRN@config$directories$output_plots, outputFolder)
  
  
  
}


.optimizeSpaceGRN <- function(df) {
  
  if (is.null(df)) {
    message = "Data frame for optimization not found"
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  } 
  
  df  %>%
    dplyr::mutate(TF.name = as.factor(.data$TF.name ),
                  TF_peak.r_bin = as.factor(.data$TF_peak.r_bin),
                  peak.ID = as.factor(.data$peak.ID),
                  TF_peak.fdr_direction = as.factor(.data$TF_peak.fdr_direction),
                  TF_peak.connectionType = as.factor(.data$TF_peak.connectionType))
  
}

.getOutputFileName <- function (name) {
  
  # List of all output file names
  allNames = list(
    
    "plot_pca"                       = "PCA_sharedSamples",
    "plot_TFPeak_fdr"                = "TF_peak.fdrCurves",
    "plot_TFPeak_fdr_GC"             = "TF_peak.GCCorrection",
    "plot_TFPeak_TFActivity_QC"      = "TF_peak.TFActivity_QC",
    "plot_class_density"             = "TF_classification_densityPlots",
    "plot_class_medianClass"         = "TF_classification_stringencyThresholds",
    "plot_class_densityClass"        = "TF_classification_summaryHeatmap",
    "plot_peakGene_diag"             = "peakGene_diagnosticPlots",
    "plot_peakGene_IHW_diag"         = "peakGene_IHW.diagnosticPlots",
    "plot_connectionSummary_heatmap" = "GRN.connectionSummary_heatmap",
    "plot_connectionSummary_boxplot" = "GRN.connectionSummary_boxplot",
    "plot_generalEnrichment"         = "GRN.overall_enrichment",
    "plot_communityStats"            = "GRN.community_stats",
    "plot_communityEnrichment"       = "GRN.community_enrichment",
    "plot_generalNetworkStats"       = "GRN.overall_stats",
    "plot_TFEnrichment"              = "GRN.TF_enrichment",
    "plot_network"                   = "GRN.network_visualisation"
    
  )
  
  if (is.null(allNames[[name]])) {
    message = paste0("Name ", name, " not defined in list allNames.")
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  } 
  
  return(allNames[[name]])
  
  
}

.getPermutationSuffixStr <- function(permutation) {
  
  if (permutation == 0) {
    suffixFile = "_original"
  } else {
    suffixFile = "_permuted"
  }
  
  suffixFile
}


# Helper function to automatically record the function calls to keep a better record
.addFunctionLogToObject <- function(GRN) {
  
  #listName = gsub("\\(|\\)", "", match.call()[1], perl = TRUE)
  functionName = evalq(match.call(), parent.frame(1))[1]
  listName = gsub("\\(|\\)", "", functionName, perl = TRUE)
  listName = gsub("GRNdev::", "", listName, fixed = TRUE)
  listName = gsub("GRN::", "", listName, fixed = TRUE)
  
  # Compatibility with old objects
  if (is.null(GRN@config$functionParameters)) {
    GRN@config$functionParameters = list()
  }
  
  GRN@config$functionParameters[[listName]] = list()
  GRN@config$functionParameters[[listName]]$date = Sys.time()
  
  GRN@config$functionParameters[[listName]]$call = match.call.defaults(asList = FALSE)
  GRN@config$functionParameters[[listName]]$parameters = match.call.defaults()
  
  GRN
}


#' Retrieve parameters for previously used function calls and general parameters for a \code{\linkS4class{GRN}} object. 
#' 
#' \strong{Note: This function, as all \code{get} functions from this package, does NOT return a \code{\linkS4class{GRN}} object.}
#' 
#' @export
#' @template GRN 
#' @param name Character. Default \code{all}. Name of parameter or function name to retrieve. Set to the special keyword \code{all} to retrieve all parameters.
#' @param type Character. Either \code{function} or \code{parameter}. Default \code{parameter}. When set to \code{function}, a valid \code{GRaNIE} function name must be given that has been run before. When set to \code{parameter}, in combination with \code{name}, returns a specific parameter (as specified in \code{GRN@config})).
#' @return The requested parameters. This function does **NOT** return a \code{\linkS4class{GRN}} object.
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' params.l = getParameters(GRN, type = "parameter", name = "all")
getParameters <- function (GRN, type = "parameter", name = "all") {
  
  checkmate::assertClass(GRN, "GRN")
  GRN = .addFunctionLogToObject(GRN)
    
  GRN = .makeObjectCompatible(GRN)
    
  checkmate::assertChoice(type, c("function", "parameter"))
  
  if (type == "function") {
    
    checkmate::assertCharacter(name, any.missing = FALSE, len = 1)
    functionParameters = GRN@config$functionParameters[[name]]
    if (is.null(functionParameters)) {
      checkmate::assertChoice(name, ls(paste0("package:", utils::packageName())))
    } 
    
    return(functionParameters)
    
  } else {
    
    if (name == "all") {
      return(GRN@config)
    } else {
      parameters = GRN@config[[name]]
      if (is.null(parameters)) {
        checkmate::assertChoice(name, names(GRN@config$parameters))
      } 
      
      return(parameters)
    }
    
  }
  
}

.getMaxPermutation <- function (GRN) {
  
  if (!is.null(GRN@config$parameters$internal$nPermutations)) {
    return(GRN@config$parameters$internal$nPermutations)
  } else {
    
    # Compatibility mode for previous object versions
    return(GRN@config$parameters$nPermutations)
  }
}


#' Optional convenience function to delete intermediate data from the function \code{\link{AR_classification_wrapper}} and summary statistics that may occupy a lot of space
#' @export
#' @template GRN
#' @return An updated \code{\linkS4class{GRN}} object, with some slots being deleted (\code{GRN@data$TFs$classification} as well as \code{GRN@stats$connectionDetails.l})
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' GRN = deleteIntermediateData(GRN)
deleteIntermediateData <- function(GRN) {
  
  checkmate::assertClass(GRN, "GRN")
  GRN = .addFunctionLogToObject(GRN)
  
  GRN = .makeObjectCompatible(GRN)
  
  for (permutationCur in 0:.getMaxPermutation(GRN)) {
    
    permIndex = as.character(permutationCur)
    GRN@data$TFs$classification[[permIndex]]$TF_cor_median_foreground = NULL
    GRN@data$TFs$classification[[permIndex]]$TF_cor_median_background = NULL
    GRN@data$TFs$classification[[permIndex]]$TF_peak_cor_foreground = NULL
    GRN@data$TFs$classification[[permIndex]]$TF_peak_cor_background = NULL
    GRN@data$TFs$classification[[permIndex]]$act.rep.thres.l = NULL
  }
  
  GRN@stats$connectionDetails.l = NULL
  
  GRN
  
  
}

.checkForbiddenNames <- function(name, forbiddenNames) {
  
  if (name %in% forbiddenNames) {
    message = paste0("Name must not be one of the following reserved names: ", paste0(forbiddenNames, collapse = ","), ". Please choose a different one.")
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
}


getBasic_metadata_visualization <- function(GRN, forceRerun = FALSE) {
  
  checkmate::assertClass(GRN, "GRN")
  GRN = .addFunctionLogToObject(GRN) 
  checkmate::assertFlag(forceRerun)
  
  # Get base mean expression for genes and TFs and mean accessibility from peaks
  
  
  # TODO: Do we need this for the shuffled one?
  
  if (is.null(GRN@visualization$metadata) | forceRerun) {
    
    expMeans.m = getCounts(GRN, type = "rna", permuted = FALSE, asMatrix = TRUE)
    
    baseMean = rowMeans(expMeans.m)
    expression.df = tibble::tibble(ENSEMBL_ID = getCounts(GRN, type = "rna", permuted = FALSE, includeIDColumn = TRUE)$ENSEMBL, baseMean = baseMean) %>%
      dplyr::mutate(ENSEMBL_ID = gsub("\\..+", "", .data$ENSEMBL_ID, perl = TRUE),
                    baseMean_log = log2(baseMean + 0.01))
    
    expression_TF.df = dplyr::filter(expression.df, .data$ENSEMBL_ID %in% GRN@annotation$TFs$TF.ENSEMBL) %>%
      dplyr::left_join(GRN@annotation$TFs, by = c("ENSEMBL_ID" = "TF.ENSEMBL"))
    
    meanPeaks.df = tibble::tibble(peakID = getCounts(GRN, type = "peaks", permuted = FALSE)$peakID, 
                                  mean = rowMeans(getCounts(GRN, type = "peaks", permuted = FALSE, asMatrix = TRUE))) %>%
      dplyr::mutate(mean_log = log2(mean + 0.01))
    
    GRN@visualization$metadata = list("RNA_expression_genes" = expression.df,
                    "RNA_expression_TF"    = expression_TF.df,
                    "Peaks_accessibility"   = meanPeaks.df)
    
  } 
  
  GRN
  
}

#' Change the output directory of a GRN object
#' 
#' @export
#' @template GRN
#' @param outputDirectory Character. Default \code{.}. New output directory for all output files unless overwritten via the parameter \code{outputFolder}.
#' @return An updated \code{\linkS4class{GRN}} object, with the output directory being adjusted accordingly
#' @examples 
#' GRN = loadExampleObject()
#' GRN = changeOutputDirectory(GRN, outputDirectory = ".")
changeOutputDirectory <- function(GRN, outputDirectory = ".") {
  
    start = Sys.time()  
    
    checkmate::assertClass(GRN, "GRN")
    GRN = .addFunctionLogToObject(GRN)
    
    GRN = .makeObjectCompatible(GRN)
    
    checkmate::assertCharacter(outputDirectory, len = 1, min.chars = 1)
    
    GRN@config$directories$outputRoot   = outputDirectory
    GRN@config$directories$output_plots = paste0(outputDirectory, .Platform$file.sep, "plots", .Platform$file.sep)
    GRN@config$files$output_log         = paste0(outputDirectory, .Platform$file.sep, "GRN.log")
    
    futile.logger::flog.info(paste0("Output directory changed in the object to " , outputDirectory))
    
    
    GRN
  
}


.createTables_peakGeneQC <- function(peakGeneCorrelations.all.cur, networkType_details, colors_vec, range) {
    
    d = peakGeneCorrelations.all.cur %>% 
        dplyr::group_by(.data$r_positive, class, .data$peak_gene.p.raw.class) %>%
        dplyr::summarise(n = dplyr::n()) %>%
        dplyr::ungroup()
    
    # Some classes might be missing, add them here with explicit zeros
    for (r_pos in c(TRUE, FALSE)) {
        for (classCur in networkType_details){
            for (pclassCur in levels(peakGeneCorrelations.all.cur$peak_gene.p.raw.class)) {
                
                row = which(d$r_positive == r_pos & d$class == classCur & as.character(d$peak_gene.p.raw.class) == as.character(pclassCur))
                if (length(row) == 0) {
                    d = tibble::add_row(d, r_positive = r_pos, class = classCur, peak_gene.p.raw.class = pclassCur, n = 0)
                }
            }
        }
    }
    
    # Restore the "ordered" factor for class
    d$peak_gene.p.raw.class = factor(d$peak_gene.p.raw.class, ordered = TRUE, levels =  levels(peakGeneCorrelations.all.cur$peak_gene.p.raw.class))
    
    
    # Normalization factors
    dsum = d %>%
        dplyr::group_by(.data$r_positive, .data$class) %>%
        dplyr::summarise(sum_n = sum(.data$n))
    
    
    # Summarize per bin
    d2 = d %>%
        dplyr::group_by(class, .data$peak_gene.p.raw.class)%>%
        dplyr::summarise(sum_pos = .data$n[.data$r_positive],
                         sum_neg = .data$n[!.data$r_positive]) %>%
        dplyr::mutate(ratio_pos_raw = .data$sum_pos / .data$sum_neg,
                      enrichment_pos = .data$sum_pos / (.data$sum_pos + .data$sum_neg),
                      ratio_neg = 1 - .data$enrichment_pos) %>%
        dplyr::ungroup()
    
    # Compare between real and permuted
    normFactor_real = dplyr::filter(dsum, class ==  !! (networkType_details[1])) %>%  dplyr::pull(.data$sum_n) %>% sum() /
        dplyr::filter(dsum, class ==  !! (networkType_details[2])) %>%  dplyr::pull(.data$sum_n) %>% sum()
    
    # ratio_norm not used currently, no normalization necessary here or not even useful because we dont want to normalize the r_pos and r_neg ratios: These are signal in a way. Only when comparing between real and permuted, we have to account for sample size for corrections
    d3 = d %>%
        dplyr::group_by(.data$peak_gene.p.raw.class, .data$r_positive) %>%
        dplyr::summarise(n_real     = .data$n[class == !! (names(colors_vec)[1]) ],
                         n_permuted = .data$n[class == !! (names(colors_vec)[2]) ]) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(ratio_real_raw = .data$n_real / .data$n_permuted,
                      ratio_real_norm = .data$ratio_real_raw / normFactor_real,
                      enrichment_real      = .data$n_real / (.data$n_real + .data$n_permuted),
                      enrichment_real_norm = (.data$n_real / normFactor_real) / ((.data$n_real / normFactor_real) + .data$n_permuted)) 
    
    
    stopifnot(identical(levels(d2$peak_gene.p.raw.class), levels(d3$peak_gene.p.raw.class)))
    # 2 enrichment bar plots but combined using facet_wrap
    d2$set = "r+ / r-"; d3$set = "real / permuted" 
    d_merged <- tibble::tibble(peak_gene.p.raw.class = c(as.character(d2$peak_gene.p.raw.class), 
                                                         as.character(d3$peak_gene.p.raw.class)),
                               ratio = c(d2$ratio_pos_raw, d3$ratio_real_norm),
                               classAll = c(as.character(d2$class), d3$r_positive),
                               set = c(d2$set, d3$set)) %>%
        dplyr::mutate(classAll = factor(.data$classAll, levels=c(paste0("real_",range), paste0("random_",range), "TRUE", "FALSE")),
                      peak_gene.p.raw.class = factor(.data$peak_gene.p.raw.class, levels = levels(d2$peak_gene.p.raw.class)))
    
    d4 = tibble::tibble(peak_gene.p.raw.class = unique(d$peak_gene.p.raw.class), 
                        n_rpos_real = NA_integer_, n_rpos_random = NA_integer_,
                        n_rneg_real = NA_integer_, n_rneg_random = NA_integer_,
                        ratio_permuted_real_rpos_norm = NA_real_,
                        ratio_permuted_real_rneg_norm = NA_real_)
    
    for (i in seq_len(nrow(d4))) {
        row_d2 = which(d2$class == networkType_details[1] & d2$peak_gene.p.raw.class == d4$peak_gene.p.raw.class[i])
        stopifnot(length(row_d2) == 1)
        d4[i, "n_rpos_real"] = d2[row_d2, "sum_pos"] %>% unlist()
        d4[i, "n_rneg_real"] = d2[row_d2, "sum_neg"] %>% unlist()
        row_d2 = which(d2$class == paste0("random_",range) & d2$peak_gene.p.raw.class == d4$peak_gene.p.raw.class[i])
        d4[i, "n_rpos_random"] = d2[row_d2, "sum_pos"] %>% unlist()
        d4[i, "n_rneg_random"] = d2[row_d2, "sum_neg"] %>% unlist()
        
        row_d3 = which(d3$r_positive == TRUE & d3$peak_gene.p.raw.class == d4$peak_gene.p.raw.class[i])
        d4[i, "ratio_permuted_real_rpos_norm"] = 1- d3[row_d3, "ratio_real_norm"] %>% unlist()
        row_d3 = which(d3$r_positive == FALSE & d3$peak_gene.p.raw.class == d4$peak_gene.p.raw.class[i])
        d4[i, "ratio_permuted_real_rneg_norm"] = 1- d3[row_d3, "ratio_real_norm"] %>% unlist()
    }
    
    d4 = d4 %>%
        dplyr::mutate(ratio_rneg_rpos_real = .data$n_rneg_real / (.data$n_rneg_real + .data$n_rpos_real),
                      ratio_rneg_rpos_random = .data$n_rneg_random / (.data$n_rneg_random + .data$n_rpos_random),
                      peak_gene.p.raw.class.bin = as.numeric(.data$peak_gene.p.raw.class)) %>%
        dplyr::arrange(.data$peak_gene.p.raw.class.bin)
    
    d4_melt = reshape2::melt(d4, id  = c("peak_gene.p.raw.class.bin", "peak_gene.p.raw.class")) %>%
        dplyr::filter(grepl("ratio", .data$variable))
    
    
    list(d = d, d2 = d2, d3 = d3, d4 = d4, d4_melt = d4_melt, d_merged = d_merged)
    
}

.classFreq_label <- function(tbl_freq) {
    paste0(names(tbl_freq), " (", tbl_freq, ")")
}

######## Quantify source of variation ########

#' Quantify and interpret multiple sources of biological and technical variation for features (TFs, peaks, and genes) in a \code{\linkS4class{GRN}} object
#' 
#' Runs the main function \code{fitExtractVarPartModel} of the package \code{variancePartition}: Fits a linear (mixed) model to estimate contribution of multiple sources of variation while simultaneously correcting for all other variables for the features in a GRN object (TFs, peaks, and genes) given particular metadata. The function reports the fraction of variance attributable to each metadata variable.
#' \strong{Note: The results are not added to \code{GRN@connections$all.filtered}, rerun the function \code{\link{getGRNConnections}} and set \code{include_variancePartitionResults} to \code{TRUE} to do so}.
#' The results object is stored in \code{GRN@stats$variancePartition} and can be used for the various diagnostic and plotting functions from \code{variancePartition}.
#' 
#' The normalized count matrices are used as input for \code{fitExtractVarPartModel}. 
#' 
#' @template GRN 
#' @param formula Character(1). Either \code{auto} or a manually defined formula to be used for the model fitting. Default \code{auto}. Must include only terms that are part of the metadata as specified with the \code{metadata} parameter. If set to \code{auto}, the formula will be build automatically based on all metadata variables as specified with the \code{metadata} parameter. By default, numerical variables will be modeled as fixed effects, while variables that are defined as factors or can be converted to factors (characters and logical variables) are modeled as random effects as recommended by the \code{variancePartition} package.
#' @param metadata Character vector. Default \code{all}. Vector of column names from the metadata data frame that was provided when using the function 
#' \code{\link{addData}}. Must either contain the special keyword \code{all} to denote that all (!) metadata columns from \code{GRN@data$metadata} are taken
#' or if not, a subset of the column names from \code{GRN@data$metadata}to include in the model fitting for \code{fitExtractVarPartModel}..
#' @param features Character(1). Either \code{all_filtered} or \code{all}. Default \code{all_filtered}. Should \code{variancePartition} only be run for the features (TFs, peaks and genes) from the filtered set of connections (the result of \code{\link{filterGRNAndConnectGenes}}) or for all genes that are defined in the object? If set to \code{all}, the running time is greatly increased.
#' @template nCores
#' @template forceRerun
#' @param ... Additional parameters passed on to \code{variancePartition::fitExtractVarPartModel} beyond \code{exprObj}, \code{formula} and \code{data}. See the function help for more information
# #' @seealso \code{\link{plotDiagnosticPlots_featureVariation}}
#' @seealso \code{\link{addData}}
#' @seealso \code{\link{getGRNConnections}}
#' @return An updated \code{\linkS4class{GRN}} object, with additional information added from this function to \code{GRN@stats$variancePartition} as well as the elements \code{genes}, \code{consensusPeaks} and \code{TFs} within \code{GRN@annotation}. 
#' As noted above, the results are not added to \code{GRN@connections$all.filtered}; rerun the function \code{\link{getGRNConnections}} and set \code{include_variancePartitionResults} to \code{TRUE} to include the results in the eGRN output table.
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' GRN = add_featureVariation(GRN, metadata = c("mt_frac"), forceRerun = TRUE)
#' @export
add_featureVariation <- function (GRN, 
                                    formula = "auto", 
                                    metadata = c("all"),
                                    features = "all_filtered", 
                                    nCores = 1,
                                    forceRerun = FALSE, 
                                    ...) {
    
    start = Sys.time()  
    
    checkmate::assertClass(GRN, "GRN")
    GRN = .addFunctionLogToObject(GRN)
    
    GRN = .makeObjectCompatible(GRN)
    
    checkmate::assertCharacter(formula, min.chars = 2, len = 1)
    
    packageMessage = paste0("The package variancePartition is not installed but required for this function.")
    .checkPackageInstallation("variancePartition", packageMessage, isWarning = FALSE)

    if (is.null(GRN@data$metadata)) {
        message = paste0("No metadata was found (GRN@data$metadata is NULL), cannot run variancePartition without metadata. Reren the addData function and provide metadata.")
        
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)    
    }
    checkmate::assert(checkmate::checkChoice(metadata, "all"), checkmate::checkSubset(metadata, colnames(GRN@data$metadata)))

    checkmate::assertChoice(features, c("all_filtered", "all"))
    checkmate::assertIntegerish(nCores, lower = 1)
    checkmate::assertFlag(forceRerun)
    
    if (is.null(GRN@stats$variancePartition$RNA) | is.null(GRN@stats$variancePartition$peaks) | forceRerun) {
        
        if (features == "all_filtered") {
            .checkExistanceFilteredConnections(GRN)
        }
        
        # Prepare the metadata
        if ("all" %in% metadata) {
            columnsToSelect = colnames(GRN@data$metadata)
        } else {
            columnsToSelect = unique(c(metadata, "has_both"))
        }
        meta <- GRN@data$metadata %>% 
            dplyr::filter(.data$has_both == TRUE) %>%
            dplyr::select(tidyselect::one_of(columnsToSelect)) %>%
            dplyr::mutate_if(is.character, as.factor) %>%
            dplyr::mutate_if(is.logical, as.factor) %>%
            dplyr::select(-.data$has_both)
        
        # Remove factors with only one level
        coltypes = meta %>% dplyr::summarise_all(class)
        factorVariables  = which(coltypes[1,] == "factor")
        numericVariables = which(coltypes[1,] == "numeric")
        factorVariablesNames = colnames(coltypes)[factorVariables]
        numericVariablesNames = colnames(coltypes)[numericVariables]

        nLevels = sapply(factorVariables, function(x) {nlevels(dplyr::pull(meta, colnames(coltypes)[x]))})
        nLevelOne = which(nLevels == 1)
        if (length(nLevelOne) > 0) {
            meta = dplyr::select(meta, -tidyselect::one_of(factorVariablesNames[nLevelOne]))
            factorVariablesNames = intersect(colnames(meta), factorVariablesNames)
        }
        
        if (formula == "auto") {
            
            futile.logger::flog.info(paste0("Due to formula = \"auto\", all factors will be modeled as random effects and all numerical variables as fixed effects, as recommended in the variancePartition vignette."))
            
            # See https://bioconductor.org/packages/release/bioc/vignettes/variancePartition/inst/doc/variancePartition.pdf for details
            # Factor variables are modeled as random effect
            if (length(factorVariablesNames) > 0) {
                formulaElems = paste0("(1|", factorVariablesNames, ")")
            } else {
                formulaElems  =""
            }
            
            if (length(numericVariablesNames) > 0) {
                if (length(factorVariablesNames) > 0) {
                    separator = " + " 
                } else {
                    separator = ""
                }
                # Numeric variables are modeled as fixed effect
                formulaElems2 = paste0(numericVariablesNames)
                
                message = paste0("Make sure that all variables from the metadata that are numeric are indeed really numeric and not (hidden) factors. The following variables will be treated as numeric: ", paste0(numericVariablesNames, collapse = ",", ". If one of these is in fact a factor, change the type in GRN@data$metadata and re-run, or provide the design formula manually"))
                .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
                
            } else {
                separator = ""
                formulaElems2 = c()
            }
            
            designFormula <- stats::as.formula(paste0("~ ", paste0(formulaElems, collapse = " + "), separator, paste0(formulaElems2, collapse = " + "))) 
            
        } else {
            
            message = paste0("A custom formula has been provides. This is currently not being checked for correctness and therefore may reuslt in errors when running variancePartition. In that case, make sure the provided formula is correct.")
            .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
            designFormula = formula
        }
        
        futile.logger::flog.info(paste0("The following formula will be used for model fitting: \"", deparse1(designFormula), "\". Make sure this is appropriate."))
        
        
        # RNA and ATAC norm counts
        if (features == "all_filtered") {
            genesToInclude = unique(c(GRN@connections$all.filtered[["0"]]$TF.ENSEMBL, GRN@connections$all.filtered[["0"]]$gene.ENSEMBL))
            peaksToInclude = unique(GRN@connections$all.filtered[["0"]]$peak.ID)
            
            futile.logger::flog.info(paste0("Using ", length(peaksToInclude), " peaks and ", length(genesToInclude), " TF/genes from the filtered connections."))
            
        } else if (features == "all") {
            genesToInclude = GRN@data$RNA$counts_metadata$ENSEMBL
            peaksToInclude = GRN@data$peaks$counts_metadata$peakID
            
            futile.logger::flog.info(paste0("Using all ", length(peaksToInclude), " peaks and ", length(genesToInclude), " TF/genes. This may take a long time. Consider setting features = \"all_filtered\"."))
        }
        
        data.l = list()
        #TODO
        data.l[["RNA"]] <- getCounts(GRN, type = "rna", permuted = FALSE, includeFiltered = FALSE, includeIDColumn = FALSE)
          
        data.l[["peaks"]] <-  getCounts(GRN, type = "peaks", permuted = FALSE, includeFiltered = FALSE, includeIDColumn = FALSE)
        
        
        for (dataType in c("RNA", "peaks")) {
            
            # row_order <- matrixStats::rowVars(as.matrix(data.l[[dataType]]) ) %>% order(decreasing = T)
            # data_set <- data.l[[dataType]][row_order,]
            
            # fit model
            futile.logger::flog.info(paste0("Running variancePartition and fit models for data type ", dataType, " using ", nCores, " core(s). This may take a while."))
            varPart <- variancePartition::fitExtractVarPartModel(exprObj = data.l[[dataType]], 
                                                                 formula = designFormula, 
                                                                 data = meta, 
                                                                 BPPARAM = .initBiocParallel(nCores),
                                                                 ...)
            GRN@stats$variancePartition[[dataType]]  = varPart
            
            res = as.data.frame(GRN@stats$variancePartition[[dataType]]) %>%
                dplyr::rename_all( .funs = ~ paste0("variancePartition_", .x)) %>%
                tibble::as_tibble(rownames = "ID")
            
            # Update annotation slots
            if (dataType == "RNA") {
                
                colnames(res)[2:ncol(res)] = paste0("gene.", colnames(res)[2:ncol(res)])
                
                GRN@annotation$genes = GRN@annotation$genes %>%
                    dplyr::select(-tidyselect::starts_with("gene.variancePart")) %>%
                    dplyr::left_join(res, by = c("gene.ENSEMBL" = "ID"))
                
                colnames(res) = gsub("gene.", "TF.", colnames(res))
                
                GRN@annotation$TFs = GRN@annotation$TFs %>%
                    dplyr::select(-tidyselect::starts_with("TF.variancePart")) %>%
                    dplyr::left_join(res, by = c("TF.ENSEMBL" = "ID"))
                
                
            } else {
                
                colnames(res)[2:ncol(res)] = paste0("peak.", colnames(res)[2:ncol(res)])
                
                GRN@annotation$peaks = GRN@annotation$peaks %>%
                    dplyr::select(-tidyselect::starts_with("variancePart")) %>%
                    dplyr::left_join(res, by = c("peak.ID" = "ID"))
            }
            
        }
        
        futile.logger::flog.info(paste0("The result objects have been stored in GRN@stats$variancePartition for both RNA and peaks."))
        
    } else {
        futile.logger::flog.info(paste0("Data already exists in object, nothing to do due to forceRerun = FALSE"))
    }
    
    
    .printExecutionTime(start, prefix = "")
    
    GRN
    
}


.checkExistanceFilteredConnections <- function(GRN) {
    
    if (is.null(GRN@connections$all.filtered[["0"]])) {
        message = "Could not find filtered connections (the slot GRN@connections$all.filtered[[\"0\"]] is NULL). Run the function filterGRNAndConnectGenes."
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
}

# Converts from an older GRN object format to the most current one due to internal optimizations
.makeObjectCompatible <- function(GRN) {
    
    if (is.null(GRN@annotation$TFs) & !is.null(GRN@data$TFs$translationTable)) {
        GRN@annotation$TFs = GRN@data$TFs$translationTable
        GRN@data$TFs$translationTable = NULL
    }
    
    if (!is.null(GRN@annotation$TFs)) {
    
        if (!"TF.ENSEMBL" %in% colnames(GRN@annotation$TFs)) {
          GRN@annotation$TFs = dplyr::rename(GRN@annotation$TFs, TF.ENSEMBL = .data$ENSEMBL)
        }
        if (!"TF.HOCOID" %in% colnames(GRN@annotation$TFs)) {
          GRN@annotation$TFs = dplyr::rename(GRN@annotation$TFs, TF.HOCOID = .data$HOCOID)
        }
        
        if ("ENSEMBL" %in% colnames(GRN@annotation$TFs)) {
          GRN@annotation$TFs = dplyr::select(GRN@annotation$TFs, -.data$ENSEMBL)
        }
    }

    
    # Temporary fix: Replace lincRNA with lncRNA due to a recent change in biomaRt until we update the object directly
    if ("lncRNA" %in% levels(GRN@annotation$genes$gene.type)) {
        GRN@annotation$genes = dplyr::mutate(GRN@annotation$genes, gene.type = dplyr::recode_factor(.data$gene.type, lncRNA = "lincRNA")) 
    }
    
    if (is.null(GRN@annotation$peaks) & !is.null(GRN@annotation$consensusPeaks)) {
        GRN@annotation[["peaks"]] = GRN@annotation[["consensusPeaks"]]
        GRN@annotation[["consensusPeaks"]] = NULL
        GRN@annotation[["peaks_obj"]] = GRN@annotation[["consensusPeaks_obj"]]
        GRN@annotation[["consensusPeaks_obj"]] = NULL
        
    }
    
    # Renamed count slots and their structure
    # 1. peaks
    if (!is.null(GRN@data$peaks[["counts_orig"]])) {
        GRN@data$peaks[["counts_orig"]] = NULL
    }
    if (is.null(GRN@data$peaks[["counts"]])) {
        GRN@data$peaks[["counts"]] = .storeAsMatrixOrSparseMatrix(GRN, df = GRN@data$peaks$counts_norm %>% dplyr::select(-.data$isFiltered), 
                                                                  ID_column = "peakID", slotName = "GRN@data$peaks$counts")
    }
    if (!is.null(GRN@data$peaks[["counts_norm"]])) {
        GRN@data$peaks[["counts_norm"]] = NULL
    }
    
    if (is.null(GRN@data$peaks[["counts_metadata"]])) {
        GRN@data$peaks[["counts_metadata"]] = .createConsensusPeaksDF(GRN@data$peaks[["counts"]] %>% as.data.frame() %>% tibble::rownames_to_column("peakID")) 
        stopifnot(c("chr", "start", "end", "peakID", "isFiltered") %in% colnames(GRN@data$peaks$counts_metadata))
    }
    if (!is.null(GRN@data$peaks[["consensusPeaks"]])) {
        GRN@data$peaks[["consensusPeaks"]] = NULL
    }
    
    # 2. RNA
    if (!is.null(GRN@data$RNA[["counts_orig"]])) {
        GRN@data$RNA[["counts_orig"]] = NULL
    }
    if (is.null(GRN@data$RNA[["counts"]]) & !is.null(GRN@data$RNA$counts_norm.l[["0"]])) {
        GRN@data$RNA[["counts"]] = .storeAsMatrixOrSparseMatrix(GRN, df = GRN@data$RNA$counts_norm.l[["0"]] %>% dplyr::select(-.data$isFiltered), 
                                                                ID_column = "ENSEMBL", slotName = "GRN@data$RNA$counts")
    }
    if (!is.null(GRN@data$RNA[["counts_norm.l"]])) {
        GRN@data$RNA[["counts_norm.l"]] = NULL
    }
    
    if (is.null(GRN@data$RNA[["counts_metadata"]]) & !is.null(GRN@data$RNA$counts)) {
        GRN@data$RNA[["counts_metadata"]] = tibble::tibble(ID = GRN@data$RNA$counts %>% as.data.frame() %>% tibble::rownames_to_column("ENSEMBL"), isFiltered = FALSE)
    }
    
    if (is.null(GRN@data$RNA[["counts_permuted_index"]]) & !is.null(GRN@data$RNA$counts)) {
        GRN@data$RNA[["counts_permuted_index"]] = .shuffleColumns(as.data.frame(GRN@data$RNA$counts), 
                                                             .getMaxPermutation(GRN), returnUnshuffled = FALSE, returnAsList = FALSE)
    }

    
    GRN
}



  