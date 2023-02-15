###################
# DATA PROCESSING #
###################

# Needed
.filterPeaksByRowMeans <- function(peakCounts, TF.peakMatrix = NULL, minMean = 1, idColumn = "peakID") {
  
  start = Sys.time()
  futile.logger::flog.info(paste0("Filter peaks with a mean across samples of smaller than ", minMean))
  
  index_lastColumn = which(colnames(peakCounts) == idColumn)
  stopifnot(length(index_lastColumn) == 1)
  
  rowMeans2 = rowMeans(peakCounts[,-index_lastColumn])
  rowsToDelete = which(rowMeans2 < minMean)
  if (length(rowsToDelete) > 0) {
    futile.logger::flog.info(paste0("Removed ", length(rowsToDelete), " peaks out of ", nrow(peakCounts), 
                                    " because they had a row mean smaller than ", minMean, "."))
    
    if (!is.null(TF.peakMatrix)) {
      # Filter these peaks also from the peakCount matrix  
      stopifnot(nrow(TF.peakMatrix) == nrow(peakCounts))
      TF.peakMatrix = TF.peakMatrix[-rowsToDelete,]
    }
    
    peakCounts = peakCounts[-rowsToDelete,]
  }
  
  .printExecutionTime(start)
  
  if (!is.null(TF.peakMatrix)) { 
    list(peakCounts = peakCounts, bindingMatrix = TF.peakMatrix)
  } else {
    peakCounts
  }
  
}


# Needed
.intersectData <- function(countsRNA, countsPeaks, idColumn_RNA = "ENSEMBL", idColumn_peaks = "peakID") {
  
  start = Sys.time()
  futile.logger::flog.info(paste0("Subset RNA and peaks and keep only shared samples"))
  
  stopifnot(idColumn_RNA %in% colnames(countsRNA))
  stopifnot(idColumn_peaks %in% colnames(countsPeaks))
  
  futile.logger::flog.info(paste0(" Number of samples for RNA before filtering: " , ncol(countsRNA) - 1))
  futile.logger::flog.info(paste0(" Number of samples for peaks before filtering: ", ncol(countsPeaks) - 1))
  
  # Subset peaks and RNA to the same set of samples
  sharedColumns = intersect(colnames(countsRNA), colnames(countsPeaks))
  
  if (length(sharedColumns) == 0) {
    message = "RNA and peaks counts have no shared samples. Verify that the colum names in the RNA-seq counts file are identical to the names in the sample table."
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  } else {
    
    countsRNA.df  = dplyr::select(countsRNA,  tidyselect::all_of(c(idColumn_RNA,  sharedColumns)))
    countsPeaks.df = dplyr::select(countsPeaks, tidyselect::all_of(c(idColumn_peaks, sharedColumns)))
    
    futile.logger::flog.info(paste0(" ", length(sharedColumns), " samples (", 
                                    dplyr::if_else(length(sharedColumns) > 100, "names omitted here", paste0(sharedColumns, collapse = ",")), 
                                    ") are shared between the peaks and RNA-Seq data"))
    
    notIntersecting_peaks = setdiff(colnames(countsPeaks),c(sharedColumns,idColumn_peaks))
    notIntersecting_RNA  = setdiff(colnames(countsRNA),c(sharedColumns,idColumn_RNA))
    if (length(notIntersecting_peaks) > 0 ) {
      futile.logger::flog.warn(paste0("The following samples from the peaks will be ignored for the classification due to missing overlap with RNA-Seq: ", paste0(notIntersecting_peaks, collapse = ",")))
    }
    if (length(notIntersecting_RNA) > 0) {
      futile.logger::flog.warn(paste0("The following samples from the RNA-Seq data will be ignored for the classification due to missing overlap with peaks data: ", paste0(notIntersecting_RNA, collapse = ",")))
    }
  }
  
  futile.logger::flog.info(paste0(" Number of samples for RNA after filtering: " , ncol(countsRNA.df) - 1))
  futile.logger::flog.info(paste0(" Number of samples for peaks data after filtering: ", ncol(countsPeaks.df) - 1))
  .printExecutionTime(start)
  
  list(RNA = countsRNA.df, peaks = countsPeaks.df)
}

#####################
# AR CLASSIFICATION #
#####################

.my.median = function(x) median(x, na.rm = TRUE)
.my.mean   = function(x) mean(x, na.rm = TRUE)


# Needed
.readTranslationTable <- function(file, delim = " ") {

  mapping.df = .read_tidyverse_wrapper(file, type = "delim", delim = delim, col_types = readr::cols()) 
  
  # Make compatible to old translation tables
  if ("HOCOID" %in% colnames(mapping.df)) {
      mapping.df = dplyr::rename(mapping.df, ID = "HOCOID")
  }
  
  checkmate::assertSubset(c("ENSEMBL", "ID"), colnames(mapping.df))
  
  mapping.df = dplyr::mutate(mapping.df , ENSEMBL = gsub("\\..+", "", .data$ENSEMBL, perl = TRUE)) # Clean ENSEMBL IDs
  

  
  mapping.df
}


.computeForegroundAndBackgroundMatrices <- function(peakMatrix, sort.cor.m) {
    
    start = Sys.time()
    futile.logger::flog.info(paste0("Compute foreground and background as well as their median values per TF"))
    # TODO: Extend with a few steps before even
    
    # This binary matrix has peaks as rows and TFs as columns of whether or not a particular peak has a TFBS from this TF or not
    # The unique prevents colnames from changeing, with ".1" being added to it automatically in case of duplicate column names
    sel.TF.peakMatrix.df = peakMatrix[,unique(colnames(sort.cor.m))]
    
    
    # 1. Focus on peaks with TFBS overlaps
    t.cor.sel.matrix = sort.cor.m
    # Transform 0 values into NA for the matrix to speed up subsequent analysis
    t.cor.sel.matrix[(sel.TF.peakMatrix.df == 0)] = NA
    # Goal: Eliminate all correlation values for cases in which a peak has no TFBS for the particular TF
    # Same dimensions as the two matrices used for input.
    # Matrix multiplication here essentially means we only multiply each individual entry with either 1 or NA
    # The result is a matrix full of NAs and the remaining entries are the correlation values for peaks that have a TFBS for the particular TF
    t.cor.sel.matrix = sel.TF.peakMatrix.df * t.cor.sel.matrix
    # Gives one value per TF, designating the median correlation per TF across all peaks
    median.cor.tfs = sort(apply(t.cor.sel.matrix, MARGIN = 2, FUN = .my.median))
    
    # 2. Background
    # Start with the same correlation matrix
    t.cor.sel.matrix.non = sort.cor.m
    # Transform 1 values into NA for the matrix to speed up subsequent analysis
    t.cor.sel.matrix.non[(sel.TF.peakMatrix.df == 1)] = NA
    t.cor.sel.matrix.non = sel.TF.peakMatrix.df + t.cor.sel.matrix.non
    # Gives one value per TF, designating the median correlation per TF
    median.cor.tfs.non <- sort(apply(t.cor.sel.matrix.non, MARGIN=2, FUN = .my.median))
    
    # Not used thereafter
    median.cor.tfs.rest <- sort(median.cor.tfs - median.cor.tfs.non[names(median.cor.tfs)])
    
    .printExecutionTime(start)
    
    list(median_foreground = median.cor.tfs, 
         median_background = median.cor.tfs.non[names(median.cor.tfs)],
         foreground = t.cor.sel.matrix,
         background = t.cor.sel.matrix.non
    )
}


.calculate_classificationThresholds <- function(background, par.l) {
    
    start = Sys.time()
    futile.logger::flog.info(paste0("Calculate classification thresholds for repressors / activators"))
    act.rep.thres.l = list()
    
    if (is.null(par.l$internal$allClassificationThresholds)) {
        message = paste0("GRN object needs to be updated. Slot internal$allClassificationThresholds is NULL.")
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    for (thresholdCur in par.l$internal$allClassificationThresholds) {
        
        act.rep.cur = quantile(sort(apply(background, MARGIN = 2, FUN = .my.median)), 
                               probs = c(thresholdCur, 1-thresholdCur))
        # Enforce the thresholds to be at least 0, so we never have an activator despite a negative median correlation 
        # and a repressor despite a positive one
        act.rep.thres.l[[as.character(thresholdCur)]][1] = min(0, act.rep.cur[1])
        act.rep.thres.l[[as.character(thresholdCur)]][2] = max(0, act.rep.cur[2])
        
        futile.logger::flog.info(paste0(" Stringency ", thresholdCur, ": ", 
                         round(act.rep.thres.l[[as.character(thresholdCur)]][1],4), 
                         " / ", 
                         round(act.rep.thres.l[[as.character(thresholdCur)]][2],4)))
        
    }
    
    .printExecutionTime(start)
    act.rep.thres.l
}


.finalizeClassificationAndAppend <- function(output.global.TFs, median.cor.tfs, act.rep.thres.l, par.l, t.cor.sel.matrix, t.cor.sel.matrix.non, significanceThreshold_Wilcoxon = 0.05) {
    
    start = Sys.time()
    futile.logger::flog.info(paste0("Finalize classification"))
    colnameMedianCor         = paste0("median.cor.tfs")
    colnameClassificationPVal= paste0("classification_distr_rawP")
    
    AR.data = as.data.frame(median.cor.tfs)
    AR.data$TF = rownames(AR.data)
    colnames(AR.data)[1] = colnameMedianCor
    
    output.global.TFs[,colnameMedianCor] = NULL
    
    if ("TF" %in% colnames(output.global.TFs)) {
      output.global.TFs = merge(output.global.TFs, AR.data, by = "TF",all.x = TRUE)
    } else if ("TF.name" %in% colnames(output.global.TFs)) {
      output.global.TFs = merge(output.global.TFs, AR.data, by = "TF.name",all.x = TRUE)
    } else {
      message = paste0("Could npt find column for merging.")
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    
    
    # Define classes, tidyverse style
    for (thresCur in names(act.rep.thres.l)) {
        thresCur.v = act.rep.thres.l[[thresCur]]
        colnameClassificationCur = paste0("classification_q", thresCur)
        output.global.TFs[, colnameClassificationCur] = dplyr::case_when( is.na(output.global.TFs[,colnameMedianCor]) ~ "not-expressed",
                                                                   output.global.TFs[,colnameMedianCor] < thresCur.v[1] ~ "repressor",
                                                                   output.global.TFs[,colnameMedianCor] > thresCur.v[2] ~ "activator",
                                                                   TRUE ~ "undetermined")
        output.global.TFs[, colnameClassificationCur] = factor(output.global.TFs[, colnameClassificationCur], levels = names(par.l$internal$colorCategories))
        
    }
    
    if (!is.null(significanceThreshold_Wilcoxon)) {
        
        checkmate::assertNumber(significanceThreshold_Wilcoxon, lower = 0, upper = 1)
        futile.logger::flog.info(paste0(" Perform Wilcoxon test for each TF. This may take a few minutes."))
        output.global.TFs[,colnameClassificationPVal] = NULL
        # Do a Wilcoxon test for each TF as a 2nd filtering criterion
        
        pb <- progress::progress_bar$new(total = ncol(t.cor.sel.matrix))
   
        for (TFCur in colnames(t.cor.sel.matrix)) {
            
            pb$tick()
            rowNo = which(output.global.TFs$TF == TFCur)
            
            # Should normally not happen
            if (length(rowNo) != 1) {

                message = paste0("Mismatch detected between TF names in the correlation matrix and the output table. Error occured for the TF ", TFCur, ". This should not happen. Contact the authors.")
                .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
                
            }
            
            # Removing NAs actually makes a difference, as these are "artifical" anyway here due to the two matrices let's remove them
            dataMotif      = stats::na.omit(t.cor.sel.matrix[,TFCur])
            dataBackground = stats::na.omit(t.cor.sel.matrix.non[,TFCur])
            
            # Skip if NA for median correlation
            if (is.na(output.global.TFs[rowNo, colnameMedianCor])) {
                output.global.TFs[rowNo,colnameClassificationPVal] = NA
                next
            }
            
            # Test the distributions
            if (output.global.TFs[rowNo, colnameMedianCor] > 0) {
                alternativeTest = "greater"
            } else {
                alternativeTest = "less"
            }
            
            testResults = stats::wilcox.test(dataMotif, dataBackground, alternative = alternativeTest)
            
            stopifnot(length(rowNo) == 1)
            output.global.TFs[rowNo,colnameClassificationPVal] = testResults$p.value
            
        }
        
        
        ################################################
        # POST-FILTER: CHANGE SOME TFs TO UNDETERMINED #
        ################################################
        
        # Change the classification with the p-value from the distribution test
        

        
        for (thresholdCur in par.l$internal$allClassificationThresholds) {
            
            futile.logger::flog.info(paste0("  Stringency ", thresholdCur))
            
            colnameClassification      = paste0("classification_q", thresholdCur)
            colnameClassificationFinal = paste0("classification_q", thresholdCur, "_final")
            
            output.global.TFs[,colnameClassificationFinal] = output.global.TFs[,colnameClassification]
            
            TFs_to_change = dplyr::filter(output.global.TFs, (!!as.name(colnameClassification) == "activator" | 
                                                                  !!as.name(colnameClassification) == "repressor") & 
                                              !!as.name( colnameClassificationPVal) > !!significanceThreshold_Wilcoxon) %>%
                            dplyr::pull(.data$TF)
            
            # Filter some TFs to be undetermined
            if (length(TFs_to_change) > 0) {
                futile.logger::flog.info(paste0("   Change the following TFs to 'undetermined' as they were classified as activator/repressor before but the Wilcoxon test was not significant: ", paste0(TFs_to_change, collapse = ",")))
                
                output.global.TFs[which(output.global.TFs$TF  %in% TFs_to_change), colnameClassificationFinal] = "undetermined"
            }
            
            
        } # enhd for each threshold
        
    } # end if doWilcoxon
    
    # Print a summary of the classification
    futile.logger::flog.info(" Summary of classification:")
    colnamesIndex = which(grepl("final", colnames(output.global.TFs)))
    for (colnameCur in colnamesIndex) {
        futile.logger::flog.info(paste0("  Column ", colnames(output.global.TFs)[colnameCur]))
        tbl = table(output.global.TFs %>% dplyr::pull(colnameCur))
        futile.logger::flog.info(paste0("   ", paste0(names(tbl), ": ", tbl), collapse = ", "))
    }
    
    .printExecutionTime(start)
    output.global.TFs
    
}


##################
# PLOT FUNCTIONS #
##################

# Density plots for TFs
#' @import graphics
.plot_density <- function(foreground.m, background.m, corMethod, file = NULL, ...) {
    
    start = Sys.time()
    futile.logger::flog.info(paste0("Plotting density plots with foreground and background for each TF", ifelse(is.null(file), "", paste0(" to file ", file))))
    stopifnot(identical(colnames(foreground.m), colnames(background.m)))
    
    # 1. Determine maximum y-values across all TFs
    yMax = 2
    for (colCur in seq_len(ncol(foreground.m))) {
      
       n_notNA = length(which(!is.na(foreground.m[,colCur])))
       if (n_notNA > 1){
         yMaxCur = max(c(stats::density(foreground.m[,colCur], na.rm = TRUE)$y, stats::density(background.m[,colCur], na.rm = TRUE)$y))
         
         if (yMaxCur > yMax) {
           yMax = yMaxCur + 0.1
         }
       }
       
    }
    
    if (!is.null(file)) {
         grDevices::pdf(file, ...)
    }
    
    for (colCur in seq_len(ncol(foreground.m))) {
        
        TFCur = colnames(foreground.m)[colCur]
        dataMotif      = foreground.m[,colCur]
        dataBackground = background.m[,colCur]
        mainLabel = paste0(TFCur," (#TFBS = ",length(which(!is.na(dataMotif)))," )")
        
        nForeground = length(which(!is.na(dataMotif)))
        nBackground = length(which(!is.na(dataBackground)))
        
        if (nForeground > 1 & nBackground > 1 ){
          
          # Show n for both foreground and background
          legendAnno = c(paste0("Motif (n=", .prettyNum(nForeground), ")"), paste0("Non-motif (n=", .prettyNum(nBackground), ")"))
          
          plot(stats::density(dataMotif, na.rm = TRUE), xlim = c(-1,1), ylim = c(0,yMax),
               main= mainLabel, lwd=2.5, col="red", axes = FALSE, xlab = paste0(.firstLetterUppercase(corMethod), " correlation"))
          abline(v=0, col="black", lty=2)
          legend("topleft",box.col = grDevices::adjustcolor("white",alpha.f = 0), legend = legendAnno, lwd = c(2,2),cex = 0.8, col = c("red","darkgrey"), lty = c(1,1) )
          axis(side = 1, lwd = 1, line = 0)
          axis(side = 2, lwd = 1, line = 0, las = 1)
          
          lines(stats::density(dataBackground, na.rm = TRUE), lwd = 2.5, col = "darkgrey")
          
        } else {
          futile.logger::flog.warn(paste0(" Not enough data for estimating densities for TF ", TFCur, ", skip for plotting."))
        }
        
        
    } 
    
    if (!is.null(file)) {
         grDevices::dev.off()
    }
    
    .printExecutionTime(start)
    
}


#' @import graphics
.plot_AR_thresholds  <- function(median.cor.tfs, median.cor.tfs.non, par.l, act.rep.thres.l, corMethod, file = NULL, ...) {
    
    start = Sys.time()
    futile.logger::flog.info(paste0("Plotting AR summary plot", ifelse(is.null(file), "", paste0(" to file ", file))))
    
    xlab=paste0("Median ", .firstLetterUppercase(corMethod), " correlation (r)")
    ylab=""
    
    xlimMax = max(abs(c(range(median.cor.tfs), range(median.cor.tfs.non))))
    xlim =  c(-(xlimMax * 1.1), (xlimMax * 1.1))
    ylim = c(1,length(median.cor.tfs.non))
    
    nTF = length(median.cor.tfs)
    
    stopifnot(identical(names(median.cor.tfs), names(median.cor.tfs.non)))
    
    
    if (!is.null(file)) {
         grDevices::pdf(file,  ...)
    }
    
    par(mfrow = c(1,1))
    
    for (thresCur in names(act.rep.thres.l)) {
        thresCur.v = act.rep.thres.l[[thresCur]]
        
        thresCur_upper = (1 - as.numeric(thresCur)) * 100
        thresCur_lower = as.numeric(thresCur) * 100
        
        mainCur = paste0("Stringency: ", thresCur)
        
        plot(median.cor.tfs.non, seq_len(nTF),
             xlim = xlim, ylim = ylim, main=mainCur, xlab=xlab, ylab=ylab,
             col = grDevices::adjustcolor("darkgrey",alpha.f=1), pch = 16, cex = 0.5, axes = FALSE)
        points(median.cor.tfs, seq_len(nTF),
               pch=16,  cex = 0.5, 
               col= dplyr::case_when(median.cor.tfs > thresCur.v[2] ~ par.l$internal$colorCategories["activator"],
                              median.cor.tfs < thresCur.v[1] ~ par.l$internal$colorCategories["repressor"],
                              TRUE ~ par.l$internal$colorCategories["undetermined"])
        ) 
        

        
        text(x  = c((thresCur.v[1]),(thresCur.v[2])),
             y = c(nTF, nTF), pos = c(2,4),
             labels  = c(paste0(thresCur_lower, " percentile\n", "(", round(thresCur.v[1],5), ")"), 
                       paste0(thresCur_upper, " percentile\n", "(", round(thresCur.v[2],5), ")")),
             cex = 0.7, col = c("black","black"))
        abline(v = thresCur.v[1], col= par.l$internal$colorCategories["repressor"])
        abline(v = thresCur.v[2], col= par.l$internal$colorCategories["activator"])
        
        dataPoints = c(median.cor.tfs.non, median.cor.tfs)
        xAxisLimits = c(min(dataPoints) * 1.1, max(dataPoints) * 1.1)
        xAxisLimits[1] = max(-1, xAxisLimits[1])
        xAxisLimits[2] = min(1, xAxisLimits[2])
        
        # Set the limits dynamically
        if (xAxisLimits[1] > -0.1 | xAxisLimits[2] < 0.1) {
            defaultLimits = seq(-1,1,0.02)
        } else {
            defaultLimits = seq(-1,1,0.1)
        }
        
        defaultLimits = defaultLimits[-c(which(defaultLimits < 0 & defaultLimits < xAxisLimits[1] - 0.2), 
                                         which(defaultLimits > 0 & defaultLimits > xAxisLimits[2] + 0.2))]
        
  
        axis(side = 1, at = defaultLimits, lwd = 1, line = 0, cex = 1)
        
    }
    
    if (!is.null(file)) {
         grDevices::dev.off()
    }
    
    .printExecutionTime(start)
    
}


# Code from Armando Reyes
.plot_heatmapAR <- function(TF.peakMatrix.df, TF_mapping.df.exp, sort.cor.m, par.l, corMethod,
                           median.cor.tfs, median.cor.tfs.non, act.rep.thres.l, finalClassification = NULL,  file = NULL, ...) {
    
    start = Sys.time()
    futile.logger::flog.info(paste0("Plotting AR heatmap", ifelse(is.null(file), "", paste0(" to file ", file))))
    
    
    missingGenes = which(!TF_mapping.df.exp$TF.ID %in% colnames(sort.cor.m))
    if (length(missingGenes) > 0) {
        TF_mapping.df.exp = dplyr::filter(TF_mapping.df.exp, .data$TF.ID %in% colnames(sort.cor.m))
    }
    
    if (!is.null(file)) {
         grDevices::pdf(file, ...)
    }
    
    cor.r.filt.m <- sort.cor.m[,as.character(TF_mapping.df.exp$TF.ID)]
    
    stopifnot(identical(colnames( cor.r.filt.m), as.character(TF_mapping.df.exp$TF.ID)))
    
    BREAKS = seq(-1,1,0.05)
    diffDensityMat = matrix(NA, nrow = ncol( cor.r.filt.m), ncol = length(BREAKS) - 1)
    rownames(diffDensityMat) = TF_mapping.df.exp$TF.ID
    
    TF_Peak_all.m <- TF.peakMatrix.df
    TF_Peak.m <- TF_Peak_all.m
    
    for (i in seq_len(ncol(cor.r.filt.m))) {
        TF = colnames( cor.r.filt.m)[i]
        ## for the background, use all peaks
        h_noMotif = hist( cor.r.filt.m[,TF][TF_Peak_all.m[,TF] == 0], breaks = BREAKS, plot = FALSE)
        ## for the foreground use only peaks with less than min_mot_n different TF motifs
        h_Motif   = hist( cor.r.filt.m[,TF][TF_Peak.m[,TF]     != 0], breaks = BREAKS, plot = FALSE)
        diff_density = h_Motif$density - h_noMotif$density
        diffDensityMat[rownames(diffDensityMat) == TF, ] <- diff_density
    }
    diffDensityMat = diffDensityMat[!is.na(diffDensityMat[,1]),]
    colnames(diffDensityMat) = signif(h_Motif$mids,2)
    # quantile(diffDensityMat)
    
    ## check to what extent the number of TF motifs affects the density values
    n_min = dplyr::if_else(colSums(TF_Peak.m) < nrow(TF_Peak.m),colSums(TF_Peak.m), nrow(TF_Peak.m) - colSums(TF_Peak.m))
    names(n_min) = TF_mapping.df.exp$TF.ID#[match(names(n_min), as.character(tf2ensg$ENSEMBL))]
    n_min <- sapply(split(n_min,names(n_min)),sum)
    
    # Make sure n_min and diffDenityMat are compatible because some NA rows may have been filtered out for diffDensityMat
    n_min <- n_min[rownames(diffDensityMat)]
    #quantile(n_min)
    remove_smallN = which(n_min < par.l$internal$plot_minNoTFBS_heatmap)
    cor(n_min[-remove_smallN], matrixStats::rowMaxs(diffDensityMat)[-remove_smallN], method = corMethod)
    
    factorClassificationPlot <- sort(median.cor.tfs, decreasing = TRUE)
    diffDensityMat_Plot = diffDensityMat[match(names(factorClassificationPlot), rownames(diffDensityMat)), ]
    diffDensityMat_Plot = diffDensityMat_Plot[!is.na(rownames(diffDensityMat_Plot)),]
    annotation_rowDF = data.frame(median_diff = factorClassificationPlot[match(rownames(diffDensityMat_Plot), 
                                  names(factorClassificationPlot))])
    
    # Define the annotation row data frame with one column per threshold, with each TF colored according to its classification status
    anno_rowDF = data.frame(matrix(NA, nrow = nrow(diffDensityMat_Plot), ncol = 0))
    rownames(anno_rowDF) = rownames(diffDensityMat_Plot)
    annotation_colors = list()
    for (thresCur in names(act.rep.thres.l)) {

        nameCur = paste0(as.numeric(thresCur)*100, " / ", (1 - as.numeric(thresCur))*100, " %")
        colBreaks = unique(c((-1),
                             act.rep.thres.l[[thresCur]][1], 
                             act.rep.thres.l[[thresCur]][2],
                             1))

        anno_rowDF[,nameCur] = cut(annotation_rowDF$median_diff, breaks = colBreaks, labels = c("repressor", "undetermined", "activator"))
        
        if (is.null(finalClassification)) {
            # Plot original colors
            colors = c(par.l$internal$colorCategories["repressor"],par.l$internal$colorCategories["not-expressed"], par.l$internal$colorCategories["activator"])
            names(colors) = levels(anno_rowDF[,nameCur])
            annotation_colors[[nameCur]] = colors
        } else {
            # Lighter colors because not final classification
            colors = c("#f18384", par.l$internal$colorCategories["not-expressed"], "#aadaa8")
            names(colors) = levels(anno_rowDF[,nameCur])
            annotation_colors[[nameCur]] = colors
        }
        
        
        # Incorporate the provided final classification in here
        if (!is.null(finalClassification)) {
            
            colnameTable = paste0("classification_q", thresCur, "_final")
            
            checkmate::assertDataFrame(finalClassification)
            checkmate::assertSubset(c(colnameTable, "TF"), colnames(finalClassification))
            
            colnameAnno = paste0(nameCur, " final")
            matchTables = match(rownames(anno_rowDF), finalClassification$TF)
            anno_rowDF[,colnameAnno] = as.character(finalClassification[matchTables, colnameTable])
            
            # QC: Check if there are any impossible transitions
            nRows = nrow(anno_rowDF[which(anno_rowDF[,nameCur] == "undetermined" & anno_rowDF[,colnameAnno] == "activator"),])
            if (nRows > 0) {
                stop("Inconsistency deteced")
            }
            colors = c(par.l$internal$colorCategories["repressor"],par.l$internal$colorCategories["not-expressed"], par.l$internal$colorCategories["activator"])
            
            names(colors) = levels(anno_rowDF[,nameCur])
            annotation_colors[[colnameAnno]] = colors
            
        }

    }

    labelMain = "Summary density heatmap (foreground - background, sorted)\nfor each TF and classifications across stringencies"

    col.list = list()
    for (colCur in colnames(anno_rowDF)) {
        col.list[[colCur]] = c("repressor" = "#e41a1c", "undetermined" = "Snow3", "activator" = "#4daf4a")
    }
    
    anno_rowDF= dplyr::select(anno_rowDF, dplyr::contains("final"))
    colnames(anno_rowDF) = gsub(" final", "", colnames(anno_rowDF))
    
    left_annotation =  ComplexHeatmap::rowAnnotation(df = anno_rowDF, col = col.list, show_legend = FALSE)
    
    fontsize_row = 2
    heatmapCur = ComplexHeatmap::Heatmap(
        diffDensityMat_Plot,
        name = "Density\n(foreground -\nbackground)",
        col = NULL,
        cluster_columns = FALSE, cluster_rows = FALSE,
        row_names_side = "left", row_names_gp = grid::gpar(fontsize = fontsize_row), 
        column_title = labelMain,
        column_names_gp = grid::gpar(fontsize = 10),
        left_annotation = left_annotation,
        row_names_max_width = ComplexHeatmap::max_text_width(
            rownames(diffDensityMat_Plot), 
            gp = grid::gpar(fontsize = fontsize_row)
        ),
        row_title = "TF"
    )
    
    lgd_list = list(
        ComplexHeatmap::Legend(labels = c("activator", "undetermined", "repressor"), title = "Classification\n(shared for all\nstringencies)",
                               legend_gp = grid::gpar(fill = c("#4daf4a", "Snow3", "#e41a1c")))
    )

    
    ComplexHeatmap::draw(heatmapCur, annotation_legend_list = lgd_list, merge_legend = TRUE)
    

    
    if (!is.null(file)) {
         grDevices::dev.off()
    }
    
    .printExecutionTime(start)
    
    
} # end function

