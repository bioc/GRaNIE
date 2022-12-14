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
        
    } else {
        .printDataAlreadyExistsMessage()
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
    
    futile.logger::flog.info(paste0(" Iterate through ", nrow(overlaps.sub.filt.df), " peak-gene combinations and calculate correlations using ", nCores, " cores. This may take a few minutes."))
    
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