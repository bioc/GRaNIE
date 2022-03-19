
############### PCA ###############


#' Produce a PCA plot of the data from a \code{\linkS4class{GRN}} object
#'
#' @template GRN 
#' @template outputFolder
#' @param type Character. Either "peaks" or "rna" or "all". Type of data to plot a PCA for. "peaks" corresponds to the the open chromatin data, while "rna" refers to the RNA-seq counts. If set to "all", PCA will be done for both data modalities. In any case, PCA will be based on the original provided data before any additional normalization has been run (i.e., usually the raw data).
#' @param topn Integer vector. Default c(500,1000,5000). Number of top variable features to do PCA for. Can be a vector of different numbers (see default).
#' @template forceRerun
#' @template basenameOutput
#' @template plotAsPDF
#' @template pdf_width
#' @template pdf_height
#' @return The same \code{\linkS4class{GRN}} object, without modifications. In addition, for each specified \code{type}, a PDF file is produced with a PCA. We refer to the Vignettes for details and further explanations.
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' GRN = plotPCA_all(GRN, topn = 500, outputFolder = ".", type = "rna", plotAsPDF=FALSE)
#' @export
plotPCA_all <- function(GRN, outputFolder = NULL, basenameOutput = NULL, 
                        type = c("rna", "peaks"), topn = c(500,1000,5000), 
                        plotAsPDF = TRUE, pdf_width = 12, pdf_height = 12,
                        forceRerun = FALSE) {
  
  GRN = .addFunctionLogToObject(GRN)
  
  checkmate::assertClass(GRN, "GRN")
  checkmate::assertSubset(type , c("rna", "peaks", "all"))
  checkmate::assertFlag(plotAsPDF)
  checkmate::assertNumeric(pdf_width, lower = 5, upper = 99)
  checkmate::assertNumeric(pdf_height, lower = 5, upper = 99)
  checkmate::assertFlag(forceRerun)
  checkmate::assert(checkmate::checkNull(basenameOutput), checkmate::checkCharacter(basenameOutput, len = 1, min.chars = 1, any.missing = FALSE))
  
  outputFolder = .checkOutputFolder(GRN, outputFolder)
  
  if (is.null(GRN@data$metadata)) {
    message = "Slot GRN@data$metadata is NULL"
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
  uniqueSamples = GRN@config$sharedSamples
  
  if ("rna" %in% type | "all" %in% type) {
    
    # Check if we have raw integer counts. Ust vst transform if, otherwise "regular log2
    transformation = dplyr::if_else(checkmate::testClass(GRN@data$RNA$counts_orig, "DESeqDataSet"), "vst", "log2")
    
    normVector = c("raw", "normalized")
    # If data is pre-normalized, dont do any transformation, as "raw" and "normalized" are supposed to be identical
    if (GRN@config$parameters$normalization_rna == "none") {
      normVector = c("normalized")
    }
    
    for (norm in normVector) {
      
      # Only do a transformation if this is raw counts, already normalized counts dont need another transformation
      transformationCur = dplyr::if_else(norm == "raw", transformation, "none")
      
      # Changed on 08.02.22: We now never do a(nother) log transform, as the input for the density plot is already transformed counts. 
      # logTransformDensity = dplyr::if_else(norm == "raw", TRUE, FALSE)
      logTransformDensity = FALSE
      
      # Get filtered rows from normalized version always, this is correct, as isFiltered is missing in raw counts
      nonFilteredRows = which(! GRN@data$RNA$counts_norm.l[["0"]]$isFiltered)
      
      if (norm == "raw") {
        
        matrixCur = GRN@data$RNA$counts_orig[nonFilteredRows, uniqueSamples]
        
      } else {
        
        matrixCur = GRN@data$RNA$counts_norm.l[["0"]][nonFilteredRows, uniqueSamples] %>% as.matrix()
        rownames(matrixCur) = GRN@data$RNA$counts_norm.l[["0"]]$ENSEMBL[nonFilteredRows]
      }
      
      fileCur = paste0(outputFolder, 
                       dplyr::if_else(is.null(basenameOutput), .getOutputFileName("plot_pca"), basenameOutput),
                       "_RNA.", norm, ".pdf")
      if (!file.exists(fileCur) | forceRerun) {
        
        futile.logger::flog.info(paste0("\nPlotting PCA and metadata correlation of ", norm, 
                                        " RNA data for all shared samples to file ", fileCur , 
                                        "... This may take a few minutes"))
        
        .plot_PCA_wrapper(GRN, matrixCur, transformation = transformationCur, logTransformDensity = logTransformDensity, file = fileCur, topn = topn, pdf_width = pdf_width, pdf_height = pdf_height, plotAsPDF = plotAsPDF)
        
      } else {
        futile.logger::flog.info(paste0("File ", fileCur, " already exists, not overwriting due to forceRerun = FALSE"))
        
      }
    }
    
  } # end if plot RNA PCA
  
  if ("peaks" %in% type | "all" %in% type) {
    
    # Check if we have raw integer counts
    transformation = dplyr::if_else(checkmate::testClass(GRN@data$peaks$counts_orig, "DESeqDataSet"), "vst", "log2")
    
    normVector = c("raw", "normalized")
    # If data is pre-normalized, dont do any transformation, as "raw" and "normalized" are supposed to be identical
    if (GRN@config$parameters$normalization_peaks == "none") {
      normVector = c("normalized")
    }
    
    for (norm in normVector) {
      
      transformationCur = dplyr::if_else(norm == "raw", transformation, "none")
      logTransformDensity = dplyr::if_else(norm == "raw", TRUE, FALSE)
      
      nonFilteredRows = which(! GRN@data$peaks$counts_norm$isFiltered)
      
      if (norm == "raw") {
        
        matrixCur = GRN@data$peaks$counts_orig[nonFilteredRows, uniqueSamples]
        
      } else {
        
        matrixCur = GRN@data$peaks$counts_norm[nonFilteredRows, uniqueSamples] %>% as.matrix()
        rownames(matrixCur) = GRN@data$peaks$counts_norm$peakID[nonFilteredRows]
      }
      
      fileCur = paste0(outputFolder, dplyr::if_else(is.null(basenameOutput), .getOutputFileName("plot_pca"), basenameOutput), "_peaks.", norm, ".pdf")
      
      if (!file.exists(fileCur) | forceRerun) {
        
        futile.logger::flog.info(paste0("Plotting PCA and metadata correlation of ", norm, 
                                        " peaks data for all shared samples to file ", fileCur , 
                                        "... This may take a few minutes"))
        
        .plot_PCA_wrapper(GRN, matrixCur, transformation = transformationCur, logTransformDensity = logTransformDensity, file = fileCur, topn = topn, pdf_width = pdf_width, pdf_height = pdf_height, plotAsPDF = plotAsPDF)
        
      } else {
        futile.logger::flog.info(paste0("File ", fileCur, " already exists, not overwriting due to forceRerun = FALSE"))
        
      }
      
    } # end for raw and norm
    
  }
  
  
  GRN
  
  
}


.plot_PCA_wrapper <- function(GRN, counts, transformation = "vst", scale = TRUE, logTransformDensity, file, topn, pdf_width, pdf_height, plotAsPDF) {
  
  checkmate::assertChoice(transformation, c("vst", "log2", "none"))
  checkmate::assertIntegerish(topn, lower = 50, any.missing = FALSE,min.len = 1)
  
  if (transformation == "vst") {
    checkmate::assertClass(counts, "DESeqDataSet")
    # Use ALL genes for variance stabilization, and not only the top X
    dd.vst = DESeq2::varianceStabilizingTransformation(counts)
    counts.transf = SummarizedExperiment::assay(dd.vst)
    metadata = SummarizedExperiment::colData(dd.vst) %>% as.data.frame()
  } else if (transformation == "log2") {
    checkmate::assertMatrix(counts)
    counts.transf = log2(counts + 1)
    metadata = GRN@data$metadata %>% dplyr::filter(sampleID %in% colnames(counts.transf)) %>% tibble::column_to_rownames("sampleID") %>% as.data.frame()
  } else {
    checkmate::assertMatrix(counts)
    counts.transf = counts
    metadata = GRN@data$metadata %>% dplyr::filter(sampleID %in% colnames(counts.transf)) %>% tibble::column_to_rownames("sampleID") %>% as.data.frame()
  }
  
  # Add sampleID as explicit column so it is alwys included in the PCA plot to identify outliers
  metadata$sampleID = rownames(metadata)
  
  futile.logger::flog.info(paste0("Prepare PCA. Count transformation: ", transformation))
  
  
  if (plotAsPDF) {
    .checkOutputFile(file)
    grDevices::pdf(file, width = pdf_width, height = pdf_height)
    futile.logger::flog.info(paste0("Plotting to file ", file))
  } else {
    futile.logger::flog.info(paste0("Plotting directly"))
  }
  
  # calculate the variance for each gene
  rv <- matrixStats::rowVars(counts.transf)
  
  # 0. Density plot for ALL features
  .plotDensity(counts.transf, logTransform = logTransformDensity, legend = dplyr::if_else(ncol(counts.transf) > 20, FALSE, TRUE))
  
  for (topn in topn) {
    
    # select the ntop genes by variance
    select <- order(rv, decreasing=TRUE)[seq_len(min(topn, length(rv)))]
    
    # perform a PCA on the data in assay(x) for the selected genes
    # For already normalized data, scale is not necessary
    pca <- stats::prcomp(t(counts.transf[select,]), scale = scale)
    
    nPCS = min(10, ncol(pca))
    
    # the contribution to the total variance for each component
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    
    screeplot.df = tibble::tibble(PCs = factor(paste0("PC ", seq_len(nPCS)), levels = paste0("PC ", seq_len(nPCS))), variation = percentVar[seq_len(nPCS)] * 100) %>%
      dplyr::mutate(variation_sum = cumsum(variation),
                    pos = seq_len(nPCS))
    
    # 1. Scree plot
    g = ggplot(screeplot.df, aes(PCs, variation)) + geom_bar(stat = "identity") + 
      xlab("Principal component (PC)") + ylab("Explained variation / Cumultative sum (in %)") + 
      geom_point(data = screeplot.df, aes(pos, variation_sum), color = "red") + 
      geom_line(data = screeplot.df, aes(pos, variation_sum), color = "red") +
      ggtitle(paste0("Screeplot for top ", topn, " variable features")) + 
      theme_bw()
    
    plot(g)
    
    
    # 2. Metadata correlation with PCs
    metadata_columns = .findSuitableColumnsToPlot(metadata, remove1LevelFactors = TRUE, verbose = FALSE)
    .plotPCA_QC(pcaResult = pca, dataCols = colnames(counts.transf), metadataTable = metadata, metadataColumns = metadata_columns, 
                varn = topn, metadataColumns_fill = metadata_columns, 
                file = NULL, logTransform = FALSE)
    
    
    
    # 3. Correlation plots colored by metadata
    metadata_columns = .findSuitableColumnsToPlot(metadata, remove1LevelFactors = FALSE, verbose = FALSE)
    for (varCur in metadata_columns) {
      
      plotTitle = paste0("Top ", topn, " variable features, colored by \n", varCur)
      
      skipLegend = FALSE
      vecCur = metadata[, varCur] %>%unlist()
      if (is.character(vecCur) & length(unique(vecCur)) > 30) {
        skipLegend = TRUE
        plotTitle = paste0(plotTitle, " (legend skipped)")
      }
      
      result = tryCatch({
        
        # The following code is taken from plotPCA from DESeq2 and adapted here to be more flexible
        intgroup.df <- as.data.frame(metadata[, varCur, drop=FALSE])
        
        # add the varCur factors together to create a new grouping factor
        group = metadata[[varCur]]
        
        # assembly the data for the plot
        d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(counts.transf))
        
        d = d %>%
          dplyr::mutate_if(is.character, as.factor) %>%
          dplyr::mutate_if(is.logical, as.factor)
        
        
        g = g = ggplot(data=d, aes(x=PC1, y=PC2, color=group)) + geom_point(size=3) + 
          xlab(paste0("PC1: ",round(percentVar[1] * 100, 1),"% variance")) +
          ylab(paste0("PC2: ",round(percentVar[2] * 100, 1),"% variance")) +
          ggtitle(plotTitle) + coord_fixed()
        
        
        
        if (is.factor(d$group)) {
          g = g + scale_color_viridis_d()
          
        } else {
          
          is.date <- function(x) inherits(x, 'Date')
          
          if (!is.date(d$group)) {
            g = g + scale_color_viridis_c()
          } else {
            g = g + scale_color_viridis_c(trans = "date")
          }
          
          
        }
        
        
        if (skipLegend) {
          g = g + theme(legend.position = "none")
        } 
        
        plot(g)
        
      }, error = function(e) {
        futile.logger::flog.warn(paste0(" Could not plot PCA with variable ", varCur))
        plot(c(0, 1), c(0, 1), ann = FALSE, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        message = paste0(" Could not plot PCA with variable ", varCur)
        text(x = 0.5, y = 0.5, message, cex = 1.6, col = "red")
      })
      
      
    }
    
  }
  
  if (plotAsPDF) grDevices::dev.off()
  
}


.plotDensity <- function(data, logTransform = TRUE, legend = FALSE) {
  
  if (logTransform) {
    data = log2(data+1)
  } 
  
  data.df = reshape2::melt(tibble::as_tibble(data), measure.vars = colnames(data)) %>% dplyr::rename(sample = variable)
  g = ggplot(data.df, aes(value, color = sample)) + geom_density() + theme_bw()
  
  title = paste0("Density of values per sample for all features")
  if (! legend) {
    g = g + theme(legend.position = "none")
    title = paste0(title, " (Legend omitted, too many samples)")
  }
  g = g + ggtitle(title)
  
  plot(g)
  
}

# data: Matrix of your data, with column names. May or may not be vsd transformed, log2 transformed etc
# metadata: The metadata data frame. Requires row.names that are identical to the colnames in data. See examples below
# varn: Number of most variable rows (features) to take.
# metadataColumns: Covariates you want to compare against PCs in character form. Must be a subset of the column names of metadata
# scale: prcomp is used for PCA, should your data be scaled? Default TRUE, scaling is advisable
# pcCount: How many PCs do you want to compare against?
# metadataColumns_fill: One of the metadataColumns variables or multiple ones
# logTransform: log2 the data for the density plot?
# ... Used for the pdf function, e.g., for specifying the pdf page sizes

# Typical use case:
# library(SummarizedExperiment)
# plotPCA_QC(assay(vsd), metadata = coldataCur, metadataColumns = PCAVariables_red, varn  = 500, metadataColumns_fill = PCAVariables_red, scale = TRUE, pcCount = 5, logTransform = FALSE, file = "~/PCA.pdf", height = 10, width = 10)

# EXAMPLE
# data = matrix(rnorm(1000 * 10), nrow = 1000, dimnames = list(c(1:1000), paste0("sample", 1:10)))
# p = rep(0.5,10)
# metadata = data.frame(sampleID = paste0("sample", 1:10), 
#                       gender   = rep(c("male","female"), 5), 
#                       age      = sample.int(100, 10),
#                       obese    = runif(length(p)) < p,
#                      row.names = paste0("sample", 1:10) )

# metadataColumns = colnames(metadata)
# metadataColumns_fill = c("gender", "age", "obese") # Could be a subset of all metadata
# varn = 500  # Number of most variable rows in the data for calculating the associations, 
# pcCount = 5 # Number of PCs to correlate metadata with

# plotPCA_QC(data = data, metadata = metadata, metadataColumns = metadataColumns, varn = varn, metadataColumns_fill = metadataColumns_fill, 
#            scale = TRUE, pcCount = pcCount, file = "test.pdf")

.plotPCA_QC <- function(pcaResult, dataCols, metadataTable, metadataColumns, varn = 500, metadataColumns_fill, pcCount = 10, file = NULL, logTransform = FALSE,  pdf_width = 7, pdf_height = 7) {
  
  
  checkmate::assertDataFrame(metadataTable)
  metadataTable = as.data.frame(metadataTable)
  checkmate::assertIntegerish(varn, lower = 100)
  checkmate::assertFlag(logTransform)
  
  checkmate::assertSubset(metadataColumns, colnames(metadataTable))
  checkmate::assertSubset(metadataColumns_fill, metadataColumns)
  
  # Enforce checking that they are identical
  checkmate::assertSetEqual(rownames(metadataTable), dataCols)
  
  futile.logger::flog.info(paste0("  Performing and summarizing PCs across metadata for top ", varn, " features", dplyr::if_else(is.null(file), "", paste0(" to file ", file))))
  
  # Helper functions
  r2_fun <- function(x, y) {
    summary(stats::lm(y ~ x))$r.squared
  } #function to get adj. R2 values
  
  
  if (!is.null(file)) {
    .checkOutputFile(file)
    grDevices::pdf(file, width = pdf_width, height = pdf_height)
  }
  
  pc = pcaResult
  
  # The number of PCs may be smaller than the asked number, depending on the tol parameter. Check here
  if (pcCount > ncol(pc$x)) {
    message = paste0("Only ", ncol(pc$x), " PCs can be used instead of ",pcCount )
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
    pcCount = ncol(pc$x)
  }
  
  pcs_cv <- pc$x[,seq_len(pcCount)]   #Top PCs from top 1000 variable genes
  
  eig.val <- 100*summary(pc)$importance[2,seq_len(pcCount)]
  
  metadataTable = dplyr::mutate_if(metadataTable, is.character, as.factor)
  
  covariates <- data.frame(metadataTable[,metadataColumns], row.names=dataCols)
  
  r2_cv <- matrix(NA, nrow = ncol(covariates), ncol = ncol(pcs_cv), dimnames = list(colnames(covariates), colnames(pcs_cv)))
  
  for (cov in colnames(covariates)) {
    for (pc in colnames(pcs_cv)) {
      r2_cv[cov, pc] <- r2_fun(covariates[,cov], pcs_cv[,pc])
    }
  }  #get adj. R2 value for covariates ~ PC 1-6
  
  colnames(r2_cv) = paste0(colnames(pcs_cv), " (", round(eig.val, 1),  " %)")
  
  breaks = seq(0,1,0.02)
  colors = grDevices::colorRampPalette(c("white", "red"))(length(breaks))
  #colors = colorRampPalette((brewer.pal(n = 7, name = "Reds")))(length(breaks))
  
  
  ComplexHeatmap::Heatmap(
    r2_cv,
    name = "Correlation",
    col = colors,
    cluster_columns = FALSE, cluster_rows = TRUE,
    row_names_side = "right", row_names_gp = grid::gpar(fontsize = 10), 
    column_title = paste0("Metadata correlation with PCs for top ", varn, " features"),
    column_names_gp = grid::gpar(fontsize = 10),
    heatmap_legend_param = list(title = "Correlation\nwith PC",  legend_direction = "vertical"),
    row_title = "Metadata",
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.text(sprintf("%.2f", r2_cv[i, j]), x, y, gp = grid::gpar(fontsize = 10))
    }
  )
  
  
  
  if (!is.null(file)) {
    grDevices::dev.off()
  }
}



############### TF-peak QC ###############




#' Plot diagnostic plots for TF-peak connections for a \code{\linkS4class{GRN}} object
#' 
#' Due to the number of plots that this functions produces, we currently provide only the option to plot as PDF. This may change in the future.
#'
#' @template GRN 
#' @template outputFolder
#' @template basenameOutput
#' @template plotDetails
#' @param plotPermuted \code{TRUE} or \code{FALSE}. Default  \code{TRUE}. Also produce the diagnostic plots for permuted data?
#' @param nTFMax \code{NULL} or Integer. Default \code{NULL}. Maximum number of TFs to process. Can be used for testing purposes by setting this to a small number i(.e., 10)
#' @template forceRerun
#' @return The same \code{\linkS4class{GRN}} object, with added data from this function.
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' GRN = plotDiagnosticPlots_TFPeaks(GRN, outputFolder = ".", plotPermuted = FALSE, nTFMax = 2)
#' @export
plotDiagnosticPlots_TFPeaks <- function(GRN, 
                                        outputFolder = NULL, 
                                        basenameOutput = NULL, 
                                        plotDetails = FALSE,
                                        plotPermuted = TRUE,
                                        nTFMax = NULL,
                                        forceRerun = FALSE) {
  
  GRN = .addFunctionLogToObject(GRN)
  
  checkmate::assertClass(GRN, "GRN")
  checkmate::assert(checkmate::checkNull(outputFolder), checkmate::checkCharacter(outputFolder, min.chars = 1))
  checkmate::assert(checkmate::checkNull(basenameOutput), checkmate::checkCharacter(basenameOutput, len = 1, min.chars = 1, any.missing = FALSE))
  checkmate::assertFlag(plotDetails)
  checkmate::assertFlag(plotPermuted)
  checkmate::assert(checkmate::checkNull(nTFMax), checkmate::checkIntegerish(nTFMax))
  checkmate::assertFlag(forceRerun)
  
  useGCCorrection = GRN@config$parameters$useGCCorrection
  checkmate::assertFlag(useGCCorrection, na.ok = FALSE)
  
  outputFolder = .checkOutputFolder(GRN, outputFolder)
  
  for (permutationCur in 0:.getMaxPermutation(GRN)) {
    
    if (!plotPermuted & permutationCur != 0) {
      next
    }
    
    suffixFile = .getPermutationSuffixStr(permutationCur)
    
    fileCur = paste0(outputFolder, 
                     dplyr::if_else(is.null(basenameOutput), .getOutputFileName("plot_TFPeak_fdr"), basenameOutput), 
                     suffixFile, ".pdf")
    
    if (!file.exists(fileCur) | forceRerun) {
      
      heightCur = 8* length(GRN@config$TF_peak_connectionTypes)
      .plot_TF_peak_fdr(GRN, perm = permutationCur, useGCCorrection = useGCCorrection, 
                        plotDetails = plotDetails, fileCur, width = 7, height = heightCur,
                        nPagesMax = nTFMax) 
    }
    
    fileCur = paste0(outputFolder, .getOutputFileName("plot_TFPeak_fdr_GC"), suffixFile)
    if (useGCCorrection & (!file.exists(fileCur) | forceRerun)) {
      .plotTF_peak_GC_diagnosticPlots(GRN, perm = permutationCur, fileCur, width = 7, height = 8) 
    }
    
    # TODO: handle multiple activities and actually write this function
    # fileCur = paste0(outputFolder, .getOutputFileName("plot_TFPeak_TFActivity_QC"), suffixFile)
    # if ("TFActivity" %in% GRN@config$TF_peak_connectionTypes & (!file.exists(fileCur) | forceRerun)) {
    #   .plotTF_peak_TFActivity_QC(GRN, perm = permutationCur, fileCur, width = 7, height = 8) 
    # }
  }
  
  GRN
  
  
}


.generateTF_GC_diagnosticPlots <- function(TFCur, GC_classes_foreground.df, GC_classes_background.df, GC_classes_all.df, peaksForeground, peaksBackground, peaksBackgroundGC) {
  
  GC_classes_background_GC.df = peaksBackgroundGC %>%
    dplyr::group_by(GC_class) %>%
    dplyr::summarise(n=  dplyr::n(), peak_width_mean = mean(peak_width), peak_width_sd = sd(peak_width)) %>%
    dplyr::ungroup() %>% 
    tidyr::complete(GC_class, fill = list(n = 0)) %>%
    dplyr::mutate(n_rel = .data$n / nrow(peaksBackgroundGC), type = "background_GC")
  
  # TODO
  GC_classes_all2.df = rbind(GC_classes_foreground.df, GC_classes_background.df, GC_classes_background_GC.df) %>%
    dplyr::mutate(type = forcats::fct_relevel(.data$type, "foreground", "background_GC", "background_orig")) %>%
    dplyr::left_join(GC_classes_all.df, by = "GC_class") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(n.bg.needed.relFreq = dplyr::if_else(round(n.bg.needed.perc * 100, 0) > 100, "", paste0(as.character(round(n.bg.needed.perc * 100, 0)),"%")))
  
  # Current hack: Put those numbers where the no. of peaks in the GC background < the required one (i.e., where resampling occured)
  GC_classes_all2.df$n.bg.needed.relFreq[GC_classes_all2.df$type != "background_GC" | GC_classes_all2.df$n == 0] = ""
  
  label_foreground = paste0("Foreground\n(", .prettyNum(nrow(peaksForeground)), " peaks)\n")
  label_background_orig = paste0("Background (raw)\n(", .prettyNum(nrow(peaksBackground)), " peaks)\n")
  label_background_GC = paste0("Background (GC-adj.)\n(", .prettyNum(nrow(peaksBackgroundGC)), " peaks)\n")
  
  labelVec = c("background_orig" = label_background_orig, "background_GC" = label_background_GC, "foreground" = label_foreground)
  colorVec = c("background_orig" = "lightgray", "background_GC" = "darkgray", "foreground" = "black")
  
  g1 = ggplot(GC_classes_all2.df , aes(GC_class, n_rel, group = .data$type, fill = .data$type, label = n.bg.needed.relFreq)) + 
    geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.5) + 
    geom_text(vjust=-0.5, size = 3)  + 
    ylim(0,0.75) + 
    #geom_line(aes(color= type), size = 2) +
    scale_fill_manual(name = "Peakset origin", labels  = labelVec, values = colorVec) + 
    #scale_color_manual(name = "Signal", labels = labelVec, values = colorVec) + 
    theme_bw() + ggtitle(paste0(TFCur, ": Before and after GC-adjustment")) + 
    xlab("GC class from peaks") +
    ylab("Relative frequencies") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) 
  
  
  g2 = ggplot(GC_classes_all2.df , aes(GC_class, log10(.data$n+1), group = .data$type, fill = .data$type)) + 
    geom_bar(stat = "identity", position = position_dodge(width = 0.5), alpha = 1, width = 0.5) + 
    #geom_line(aes(color= type), size = 2) +
    scale_fill_manual(name = "Peakset origin", labels  = labelVec, values = colorVec) + 
    #scale_color_manual(name = "Signal", labels = labelVec, values = colorVec) + 
    theme_bw() + ggtitle(paste0(TFCur, ": Before and after GC-adjustment")) + 
    xlab("GC class from peaks") +
    ylab("Absolute frequencies (log10 + 1)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) 
  
  list(g1,g2)
  
}

.plotTF_peak_GC_diagnosticPlots <- function(GRN, perm = 0, file, height = 7, width = 7) {
  
  
  start = Sys.time()
  futile.logger::flog.info(paste0("Plotting GC diagnostic plots (for each TF)summary + TF-specific", dplyr::if_else(is.null(file), "", paste0(" to file ", file))))
  
  permIndex = as.character(perm)
  summary.df = GRN@connections$TF_peaks[[permIndex]]$summary %>% dplyr::filter(TF_peak.connectionType == "expression")
  
  g1 = ggplot(summary.df %>% dplyr::filter(direction == "pos"), aes(nForeground)) + geom_histogram(bins = 100) +
    xlab("Number peaks in foreground") + 
    theme_bw()
  
  g2 = ggplot(summary.df %>% dplyr::filter(direction == "pos"), aes(nBackground_orig)) + geom_histogram(bins = 100) +
    xlab("Number peaks in background (raw)") + 
    theme_bw()
  
  g3 = ggplot(summary.df %>% dplyr::filter(direction == "pos"), aes(nBackground)) + geom_histogram(bins = 100) +
    xlab("Number peaks in background (after GC-adjustment)") + 
    theme_bw()
  
  g4 = ggplot(summary.df %>% dplyr::filter(direction == "pos"), aes(ratio_fg_bg_orig)) + geom_histogram(bins = 100) +
    xlab("Ratio of foreground peaks / background peaks (raw)") + 
    theme_bw()
  
  g5 = ggplot(summary.df %>% dplyr::filter(direction == "pos"), aes(ratio_fg_bg)) + geom_histogram(bins = 100) +
    xlab("Ratio of foreground peaks / background peaks (after GC-adjustment)") + 
    theme_bw()
  
  g6 = ggplot(summary.df %>% dplyr::filter(direction == "pos"), aes(background_match_success)) + geom_histogram(stat = "count") +
    xlab("GC-background matching successful") + 
    theme_bw()
  
  g7 = ggplot(summary.df %>% dplyr::filter(direction == "pos"), aes(percBackgroundUsed)) + geom_histogram(stat = "count") +
    xlab("Percentage of used background for GC matching") +
    theme_bw()
  
  g8 = ggplot()
  
  # plots.l is nested, unnest 
  plots.mod.l = list()
  for (nameCur in names(GRN@stats$plots_GC)) {
    plots.mod.l[[paste0(nameCur,".rel")]] = GRN@stats$plots_GC[[nameCur]][[1]]
    plots.mod.l[[paste0(nameCur,".abs")]] = GRN@stats$plots_GC[[nameCur]][[2]]
  }
  
  
  plots_all.l = c(list(g1,g2,g3,g4,g5,g6,g7,g8), plots.mod.l)
  
  .printMultipleGraphsPerPage(plots_all.l, nCol = 1, nRow = 2, pdfFile = file, height = height, width = width)
  
}


.plot_TF_peak_fdr <- function(GRN, perm, useGCCorrection, plotDetails = FALSE, file = NULL, width = 7, height = 7, nPagesMax = NULL) {
  
  start = Sys.time()
  futile.logger::flog.info(paste0("Plotting FDR curves for each TF", dplyr::if_else(is.null(file), "", paste0(" to file ", file))))
  
  if (!is.null(file)) {
    plots.l = list()
    index = 1
  } 
  
  # Dont take all TF, some might be missing.
  connections_TF_peak = GRN@connections$TF_peaks[[as.character(perm)]]$connectionStats
  allTF = unique(connections_TF_peak$TF.name)
  nTF = ifelse(is.null(nPagesMax), length(allTF), nPagesMax)
  futile.logger::flog.info(paste0(" Including a total of ", nTF,  " TF. Preparing plots..."))
  
  
  # TODO: Check difference between TFActivity TFs and expression TFs
  

  pb <- progress::progress_bar$new(total = nTF)
  
  levels_pos<-unique(as.character(cut(GRN@config$parameters$internal$stepsFDR, breaks = GRN@config$parameters$internal$stepsFDR, right = FALSE, include.lowest = TRUE )))
  levels_neg<-unique(as.character(cut(GRN@config$parameters$internal$stepsFDR, breaks = rev(GRN@config$parameters$internal$stepsFDR), right = TRUE,  include.lowest = TRUE )))
  

  
  for (i in seq_len(nTF)) {
    pb$tick()
    TFCur = allTF[i]
    
    connections_TF_peakCur = dplyr::filter(connections_TF_peak, TF.name == TFCur)
    
    # Produce two FDR plots, coming from both directions
    for (typeCur in c("pos","neg")) {
      
      if (typeCur == "pos") {
        levelsCur = levels_pos
      } else {
        levelsCur = levels_neg
      }
      
      
      connections_TF_peak.filt = dplyr::filter(connections_TF_peakCur, TF_peak.fdr_direction == typeCur ) %>%
        dplyr::select(-TF.name, -TF_peak.fdr_direction) %>%
        dplyr::mutate(TF_peak.r_bin = factor(TF_peak.r_bin, levels = levelsCur))  %>%
        reshape2::melt(id = c("TF_peak.r_bin", "TF_peak.connectionType", "n")) 
      
      if (!useGCCorrection) {
        
        g = suppressWarnings(connections_TF_peak.filt %>% 
                               dplyr::filter(variable %in% c("TF_peak.fdr")) %>%
                               ggplot(aes(TF_peak.r_bin, value, color = .data$n, shape = variable)) +
                               geom_point(na.rm = TRUE) +
                               facet_wrap(~TF_peak.connectionType, ncol = 1) + 
                               labs(x="Correlation bin", y="FDR TF-peak", title= paste0(TFCur, ": Direction ", typeCur)) +
                               theme_bw() +
                               theme(axis.text.x=element_text(angle=60, hjust=1, size= 8)) +
                               scale_x_discrete(drop=FALSE) + 
                               scale_y_continuous(limits = c(0,1),minor_breaks = seq(0 , 1, 0.1), breaks = seq(0, 1, 0.1)) +
                               scale_shape_manual("Background", 
                                                  label = c("TF_peak.fdr" = "original"),
                                                  values = c("TF_peak.fdr" = 1)) +
                               scale_color_viridis_c("No. of connections") 
        )
        
        if (plotDetails) {
          
          g2 = connections_TF_peak.filt %>% dplyr::filter(grepl("value", variable), !grepl("orig", variable)) %>%
            ggplot(aes(TF_peak.r_bin, log10(value+1), color = variable)) + 
            geom_point(size = 1) +
            facet_wrap(~TF_peak.connectionType, ncol = 1) + 
            labs(x="Correlation bin", y="log10(Observed frequencies +1)", title= paste0(TFCur, ": Direction ", typeCur)) +
            theme_bw() +
            theme(axis.text.x=element_text(angle=60, hjust=1, size= 8)) +
            scale_x_discrete(drop=FALSE) + 
            geom_hline(yintercept = log10(1), linetype = "dotted", color = "gray") +
            scale_color_manual("Observed data for\nFDR calculation", 
                               values = c("tpvalue" = "black", 
                                          "fpvalue" = "brown1", "fpvalue_norm" = "darkred"),
                               labels = c("tpvalue" = "True positives\n(foreground, scaling ref.)\n", 
                                          "fpvalue" = "False positives\n(background, unscaled)\n", 
                                          "fpvalue_norm" = "False positives\n(background, scaled to ref.)"))
        }
        
        
        
      } else {
        
        g = suppressWarnings(connections_TF_peak.filt %>% 
                               dplyr::filter(stringr::str_starts(variable, "TF_peak.fdr")) %>% 
                               ggplot(aes(TF_peak.r_bin, value, color = .data$n, shape = variable, alpha = variable)) + 
                               geom_point(na.rm = TRUE) +
                               facet_wrap(~TF_peak.connectionType, ncol = 1) + 
                               labs(x="Correlation bin", y="FDR TF-peak", title= paste0(TFCur, ": Direction ", typeCur)) +
                               theme_bw() +
                               theme(axis.text.x=element_text(angle=60, hjust=1, size= 8)) +
                               scale_x_discrete(drop=FALSE) + 
                               scale_y_continuous(limits = c(0,1),minor_breaks = seq(0 , 1, 0.1), breaks = seq(0, 1, 0.1)) +
                               scale_shape_manual("Background", 
                                                  label = c("TF_peak.fdr" = "GC-adjusted","TF_peak.fdr_orig" = "original"),
                                                  values = c("TF_peak.fdr" = 3,"TF_peak.fdr_orig" = 1)) +
                               scale_color_viridis_c("No. of connections") + 
                               scale_alpha_discrete(range = c("TF_peak.fdr" = 0.5,"TF_peak.fdr_orig" = 1))  +
                               guides(alpha = FALSE))
        
        if (plotDetails) {
          g2 = ggplot(connections_TF_peak.filt %>% dplyr::filter(grepl("value", variable)), aes(TF_peak.r_bin, log10(value+1), color = variable)) + 
            geom_point(size = 1) +
            facet_wrap(~TF_peak.connectionType, ncol = 1) + 
            labs(x="Correlation bin", y="log10(Observed frequencies +1)", title= paste0(TFCur, ": Direction ", typeCur)) +
            theme_bw() +
            theme(axis.text.x=element_text(angle=60, hjust=1, size= 8)) +
            scale_x_discrete(drop=FALSE) + 
            geom_hline(yintercept = log10(1), linetype = "dotted", color = "gray") +
            scale_color_manual("Observed data for\nFDR calculation", 
                               values = c("tpvalue" = "black", 
                                          "fpvalue_orig" = "lightgray",  "fpvalue_norm_orig" = "darkgray",
                                          "fpvalue" = "brown1", "fpvalue_norm" = "darkred"),
                               labels = c("tpvalue" = "True positives\n(foreground, scaling ref.)\n", 
                                          "fpvalue_orig" = "False positives before GC\n(background, unscaled)\n",  
                                          "fpvalue_norm_orig" = "False positives without GC\n(background, scaled to ref.)\n",
                                          "fpvalue" = "False positives with GC\n(background, unscaled)\n", 
                                          "fpvalue_norm" = "False positives with GC\n(background, scaled to ref.)"))
        }
        
      }
      
      
      
      if (is.null(file)) {
        plot(g)
      } else {
        plots.l[[index]] = g
        index = index + 1
        
        if (plotDetails) {
          plots.l[[index]] = g2
          index = index + 1
        }
        
      }
      
      
    }
  } # end allTF
  
  if (!is.null(file)) {
    futile.logger::flog.info(paste0(" Finished generating plots, start plotting to file ",file, ". This may take a few minutes."))
    .printMultipleGraphsPerPage(plots.l, nCol = 1, nRow = 2, pdfFile = file, height = height, width = width)
    
  }
  
  .printExecutionTime(start)
}

.plotTF_peak_TFActivity_QC <- function(GRN, perm = permIndex, fileCur, width = 7, height = 8) {
  
  # TODO
}


######### peak-gene QC ##########

#' Plot diagnostic plots for peak-gene connections for a \code{\linkS4class{GRN}} object
#'
#' @template GRN 
#' @template outputFolder
#' @template basenameOutput
#' @param gene.types List of character vectors. Default list(c("all"), c("protein_coding"), c("protein_coding", "lincRNA")). Vectors of gene types to consider for the diagnostic plots. Multiple distinct combinations of gene types can be specified. For example, if set to \code{list(c("protein_coding", "lincRNA"), c("protein_coding"), c("all"))}, 3 distinct PDFs will be produced, one for each element of the list. The first file would only consider protein-coding and lincRNA genes, while the second plot only considers protein-coding ones. The special keyword "all" denotes all gene types found (usually, there are many gene types present, also more exotic and rare ones).
#' @param useFiltered Logical. TRUE or FALSE. Default FALSE. If set to \code{FALSE}, the diagnostic plots will be produced based on all peak-gene connections. This is the default and will usually be best to judge whether the background behaves as expected. If set to TRUE, the diagnostic plots will be produced based on the filtered set of connections. For this, the function \code{link{filterGRNAndConnectGenes}} must have been run before.
#' @template plotDetails
#' @param plotPerTF Logical. TRUE or FALSE. Default \code{FALSE}. If set to \code{FALSE}, the diagnostic plots will be done across all TF (the default), while setting it to \code{TRUE} will generate the QC plots TF-specifically, including "all" TF, sorted by the number of connections.
#' @template plotAsPDF
#' @template pdf_width
#' @template pdf_height
#' @template forceRerun
#' @return The same \code{\linkS4class{GRN}} object, with added data from this function.
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' # GRN = loadExampleObject()
#' # types = list(c("protein_coding"))
#' # GRN = plotDiagnosticPlots_peakGene(GRN, outputFolder=".", gene.types=types, plotAsPDF=FALSE)
#' @export
# TODO: implement forceRerun correctly
plotDiagnosticPlots_peakGene <- function(GRN, 
                                         outputFolder = NULL, 
                                         basenameOutput = NULL, 
                                         gene.types = list(c("all"), c("protein_coding"), c("protein_coding", "lincRNA")), 
                                         useFiltered = FALSE, 
                                         plotDetails = FALSE,
                                         plotPerTF = FALSE,
                                         plotAsPDF = TRUE, pdf_width = 12, pdf_height = 12,
                                         forceRerun= FALSE) {
  
  GRN = .addFunctionLogToObject(GRN)
  
  checkmate::assertClass(GRN, "GRN")
  checkmate::assert(checkmate::checkNull(outputFolder), checkmate::checkCharacter(outputFolder, min.chars = 1))
  checkmate::assert(checkmate::checkNull(basenameOutput), checkmate::checkCharacter(basenameOutput, len = 1, min.chars = 1, any.missing = FALSE))
  checkmate::assertList(gene.types, any.missing = FALSE, min.len = 1, types = "character")
  checkmate::assertFlag(useFiltered)
  checkmate::assertFlag(plotDetails)
  checkmate::assertFlag(plotPerTF)
  checkmate::assertFlag(plotAsPDF)
  checkmate::assertNumeric(pdf_width, lower = 5, upper = 99)
  checkmate::assertNumeric(pdf_height, lower = 5, upper = 99)
  checkmate::assertFlag(forceRerun)
  
  
  # For compatibility reasons, re-create the genes and peaks annotation if not present in the object
  if (! "peak.annotation" %in% colnames(GRN@annotation$consensusPeaks)) {
    GRN = .populatePeakAnnotation(GRN)
  }
  if (! "gene.CV" %in% colnames(GRN@annotation$genes)) {
    GRN = .populateGeneAnnotation(GRN)
  }
  
  
  outputFolder = .checkOutputFolder(GRN, outputFolder)
  
  if (is.null(GRN@connections$peak_genes[["0"]])) {
    message = paste0("Could not find peak-gene connections in GRN object. Run the function addConnections_peak_gene first")
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
  if (!plotAsPDF) {
    filenameCurBase = NULL
  } else {
    filenameCurBase = paste0(outputFolder, dplyr::if_else(is.null(basenameOutput), .getOutputFileName("plot_peakGene_diag"), basenameOutput), "_")
  }
  
  .plotDiagnosticPlots_peakGene_all(GRN, gene.types = gene.types, fileBase = filenameCurBase, 
                                    useFiltered = useFiltered, 
                                    plotPerTF = plotPerTF,
                                    pdf_width = pdf_width, pdf_height = pdf_height,
                                    plotDetails = plotDetails,
                                    forceRerun = forceRerun)
  
  # Summarize all data, both real and random
  # filteredStr = dplyr::if_else(useFiltered, "_filtered", "")
  # for (geneTypesSelected in gene.types) {
  #   
  #   filenameCur = paste0(outputFolder, .getOutputFileName("plot_peakGene_diag"), "_", paste0(geneTypesSelected, collapse = "+"), filteredStr, ".pdf")
  #   if (!file.exists(filenameCur) | forceRerun) {
  #    
  #   } else {
  #     futile.logger::flog.info(paste0("File ", filenameCur, " already exists. As forceRerun = FALSE, not overwriting it."))
  #   }
  #   
  # }
  
  GRN
  
}
#' @import patchwork
#' @importFrom rlang .data
.plotDiagnosticPlots_peakGene_all <- function(GRN, gene.types,
                                              useFiltered = FALSE,
                                              plotDetails = FALSE,
                                              plotPerTF = FALSE,
                                              fileBase = NULL,
                                              pdf_width = 12,
                                              pdf_height = 12,
                                              forceRerun = FALSE) {
  
  start = Sys.time()
  
  # get list of all filenames that are going to be created
  filenames.all = c()
  filteredStr = dplyr::if_else(useFiltered, "_filtered", "")
  for (geneTypesSelected in gene.types) {
    filenameCur = paste0(fileBase, paste0(geneTypesSelected, collapse = "+"), filteredStr, ".pdf")
    filenames.all = c(filenames.all, filenameCur)
  }
  
  if (!all(file.exists(filenames.all)) | forceRerun) {
    
    futile.logger::flog.info(paste0("Plotting diagnostic plots for peak-gene correlations", dplyr::if_else(is.null(fileBase), "", paste0(" to file(s) with basename ", fileBase))))
    
    cols_keep = c("peak.ID", "gene.ENSEMBL", "peak_gene.r", "peak_gene.p_raw", "peak_gene.distance")
    
    # Change from 10 to 5 
    nCategoriesBinning = 5
    probs = seq(0,1,1/nCategoriesBinning)
    
    
    range = GRN@config$parameters$promoterRange
    
    networkType_details = c(paste0("real_",range), paste0("random_",range))
    
    colors_vec = c("black", "darkgray")
    networkType_vec = c("real", "permuted")
    names(colors_vec) = names(networkType_vec) = networkType_details
    
    options(dplyr.summarise.inform = FALSE) 
    
    includeRobustCols = FALSE
    #Either take all or the filtered set of connections
    if (useFiltered) {
      
      # TODO: Robust columns filter
      peakGeneCorrelations.all = 
        rbind(dplyr::select(GRN@connections$all.filtered[["0"]], dplyr::everything()) %>%
                dplyr::mutate(class = paste0("real_",range))) %>%
        rbind(dplyr::select(GRN@connections$all.filtered[["1"]],  dplyr::everything()) %>% 
                dplyr::mutate(class = paste0("random_",range)))
      
    } else {
      
      robustColumns = c("peak_gene.p_raw.robust", "peak_gene.bias_M_p.raw", "peak_gene.bias_LS_p.raw", "peak_gene.r_robust")
      if (all(robustColumns %in% colnames(GRN@connections$peak_genes[["0"]]))) {
        includeRobustCols = TRUE
        cols_keep = c(cols_keep, robustColumns)
      }
      
      class_levels = c(paste0("real_",range), paste0("random_",range))
      
      peakGeneCorrelations.all = 
        rbind(
          dplyr::select(GRN@connections$peak_genes[["0"]], tidyselect::all_of(cols_keep)) %>%
            dplyr::mutate(class = factor(paste0("real_",range), levels = class_levels)),
          dplyr::select(GRN@connections$peak_genes[["1"]],  tidyselect::all_of(cols_keep)) %>% 
            dplyr::mutate(class = factor(paste0("random_",range), levels = class_levels))) %>%
        dplyr::left_join(dplyr::select(GRN@annotation$genes, gene.ENSEMBL, gene.type, gene.mean, gene.median, gene.CV), by = "gene.ENSEMBL") %>%
        dplyr::left_join(GRN@annotation$consensusPeaks %>% 
                           dplyr::select(-dplyr::starts_with("peak.gene."), -peak.GC.perc), by = "peak.ID") %>%
        dplyr::select(-gene.ENSEMBL)
      
    }
    
    if (!plotPerTF) {
      peakGeneCorrelations.all = dplyr::select(peakGeneCorrelations.all, -peak.ID)
    }
    
    nClasses_distance  = 10
    
    peakGeneCorrelations.all = peakGeneCorrelations.all %>%
      dplyr::mutate(r_positive = peak_gene.r > 0,
                    peak_gene.distance_class = 
                      forcats::fct_explicit_na(addNA(cut(peak_gene.distance, breaks = nClasses_distance, include.lowest = TRUE)), "random"),
                    peak_gene.distance_class_abs = forcats::fct_explicit_na(addNA(cut(abs(peak_gene.distance), 
                                                                                      breaks = nClasses_distance, include.lowest = TRUE, ordered_result = TRUE)), "random"),
                    peak_gene.p.raw.class = cut(peak_gene.p_raw, breaks = seq(0,1,0.05), include.lowest = TRUE, ordered_result = TRUE),
                    peak_gene.r.class = cut(peak_gene.r, breaks = seq(-1,1,0.05), include.lowest = TRUE, ordered_result = TRUE)) %>%
      dplyr::filter(!is.na(peak_gene.r)) # Eliminate rows with NA for peak_gene.r. This can happen if the normalized gene counts are identical across ALL samples, and due to the lack of any variation, the correlation cannot be computed
    
    
    # Oddity of cut: When breaks is specified as a single number, the range of the data is divided into breaks pieces of equal length, and then the outer limits are moved away by 0.1% of the range to ensure that the extreme values both fall within the break intervals. 
    levels(peakGeneCorrelations.all$peak_gene.distance_class_abs)[1] = 
      gsub("(-\\d+)", "0", levels(peakGeneCorrelations.all$peak_gene.distance_class_abs)[1], perl = TRUE)
    
    
    if (includeRobustCols) {
      peakGeneCorrelations.all = peakGeneCorrelations.all %>%
        dplyr::mutate(peak_gene.p_raw.robust.class = 
                        cut(peak_gene.p_raw.robust, breaks = seq(0,1,0.05), include.lowest = TRUE, ordered_result = TRUE))
    }
    
    
    # Prepare plots #
    
    colors_class = c("black", "black")
    names(colors_class)= unique(peakGeneCorrelations.all$class)
    colors_class[which(grepl("random", names(colors_class)))] = "darkgray"
    
    r_pos_class = c("black", "darkgray")
    names(r_pos_class) =c("TRUE", "FALSE")
    
    dist_class = c("dark red", "#fc9c9c")
    names(dist_class) = class_levels
    
    freqs= table(peakGeneCorrelations.all$class)
    freq_class = paste0(gsub(names(freqs), pattern = "(.+)(_.*)", replacement = "\\1"), " (n=", .prettyNum(freqs) , ")")
    # Change upstream and go with "permuted" everywhere
    freq_class = gsub(freq_class, pattern = "random", replacement = "permuted")
    names(freq_class) <- names(freqs)
    
    
    
    xlabels_peakGene_r.class = levels(peakGeneCorrelations.all$peak_gene.r.class)
    nCur = length(xlabels_peakGene_r.class)
    xlabels_peakGene_r.class[setdiff(seq_len(nCur), c(1, floor(nCur/2), nCur))] <- ""
    
    # For the last plot, which is wider, we label a few more
    xlabels_peakGene_r.class2 = levels(peakGeneCorrelations.all$peak_gene.r.class)
    nCur = length(xlabels_peakGene_r.class2)
    xlabels_peakGene_r.class2[setdiff(seq_len(nCur), c(1, floor(nCur/4), floor(nCur/2), floor(nCur/4*3), nCur))] <- ""
    
    xlabels_peakGene_praw.class = levels(peakGeneCorrelations.all$peak_gene.p.raw.class)
    nCur = length(xlabels_peakGene_praw.class)
    xlabels_peakGene_praw.class[setdiff(seq_len(nCur), c(1, floor(nCur/2), nCur))] <- ""
    
    
    for (geneTypesSelected in gene.types) {
      
      futile.logger::flog.info(paste0(" Gene type ", paste0(geneTypesSelected, collapse = "+")))
      
      if ("all" %in% geneTypesSelected) {
        indexCur = seq_len(nrow(peakGeneCorrelations.all))
      } else {
        indexCur = which(peakGeneCorrelations.all$gene.type %in% geneTypesSelected)
      }
      
      
      # START PLOTTING #
      
      if (!is.null(fileBase)) {
        filenameCur = paste0(fileBase, paste0(geneTypesSelected, collapse = "+"), filteredStr, ".pdf")
        .checkOutputFile(filenameCur)
        grDevices::pdf(file = filenameCur, width = pdf_width, height = pdf_height)
      }
      
      if (plotPerTF) {
        
        TF.nRows = rep(-1, length(GRN@config$allTF))
        #TF.nRows = rep(-1, 10)
        TF.peaks = list()
        names(TF.nRows) = GRN@config$allTF
        for (TFCur in GRN@config$allTF) {
          TF.peaks[[TFCur]] = names(which(GRN@data$TFs$TF_peak_overlap[,TFCur] == 1))
          TF.nRows[TFCur] = peakGeneCorrelations.all[indexCur,] %>% dplyr::filter(peak.ID %in% TF.peaks[[TFCur]]) %>% nrow()
        }
        
        TFs_sorted = names(sort(TF.nRows, decreasing = TRUE))
        
        allTF = c("all", TFs_sorted)
        
      } else {
        allTF = "all"
      }
      
      counter = 0
      for (TFCur in allTF) {
        
        counter = counter + 1
        
        if (length(allTF) > 1) {
          futile.logger::flog.info(paste0(" QC plots for TF ", TFCur, " (", counter, " of ", length(allTF), ")"))
        }
        
        
        if ("all" %in% geneTypesSelected) {
          indexCur = seq_len(nrow(peakGeneCorrelations.all))
        } else {
          indexCur = which(peakGeneCorrelations.all$gene.type %in% geneTypesSelected)
        }
        
        
        if (TFCur != "all") {
          indexCur = intersect(indexCur, which(peakGeneCorrelations.all$peak.ID %in% TF.peaks[[TFCur]]))
        }
        
        # Get subset also for just the real data
        indexCurReal = intersect(indexCur, which(peakGeneCorrelations.all$class == names(dist_class)[1]))
        
        
        xlabel = paste0("Correlation raw p-value")
        
        # DENSITY PLOTS P VALUE
        
        # TODO: Densities as ratio
        # https://stackoverflow.com/questions/58629959/how-can-i-extract-data-from-a-kernel-density-function-in-r-for-many-samples-at-o#:~:text=To%20compute%20the%20density%20you,can%20use%20the%20package%20spatstat%20.
        
        # Produce the labels for the class-specific subtitles
        customLabel_class = .customLabeler(table(peakGeneCorrelations.all[indexCur,]$class))
        
        r_pos_freq = table(peakGeneCorrelations.all[indexCur,]$r_positive)
        labeler_r_pos = labeller(r_positive = c("TRUE"  = paste0("r positive (", r_pos_freq["TRUE"], ")"), 
                                                "FALSE" = paste0("r negative (", r_pos_freq["FALSE"], ")")) )
        theme_main = theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = rel(0.8)),
                           axis.text.y = element_text(size=rel(0.8)),
                           panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        
        
        
        ## p-val density curves stratified by real/permuted ##
        
        gA2 = ggplot(peakGeneCorrelations.all[indexCur,], aes(peak_gene.p_raw, color = r_positive)) + geom_density()  +
          facet_wrap(~ class, labeller = labeller(class=freq_class) ) +
          xlab(xlabel) + ylab("Density") +  theme_bw() +
          scale_color_manual(labels = names(r_pos_class), values = r_pos_class) +
          theme(legend.position = "none", axis.text=element_text(size=rel(0.6)), strip.text.x = element_text(size = rel(0.8)))
        
        # Helper function to retrieve all tables and data aggregation steps for subsequent visualization
        tbl.l = .createTables_peakGeneQC(peakGeneCorrelations.all[indexCur,], networkType_details, colors_vec, range)
        
        
        xlabel = paste0("Correlation raw\np-value (binned)")
        
        
        xlabels = levels(tbl.l$d_merged$peak_gene.p.raw.class)
        xlabels[setdiff(seq_len(length(xlabels)), c(1, floor(length(xlabels)/2), length(xlabels)))] <- ""
        
        gB3 = ggplot(tbl.l$d_merged, aes(peak_gene.p.raw.class, ratio, fill = classAll)) + 
          geom_bar(stat = "identity", position="dodge", na.rm = TRUE, width = 0.5) + 
          geom_hline(yintercept = 1, linetype = "dotted") + 
          xlab(xlabel) + ylab("Ratio") +
          scale_fill_manual("Class", values = c(dist_class, r_pos_class), 
                            labels = c("real", "permuted", "r+ (r>0)", "r- (r<=0)"), 
          ) + # labels vector can be kind of manually specified here because the levels were previosly defined in a certain order
          scale_x_discrete(labels = xlabels_peakGene_praw.class) +
          theme_bw() +  
          #theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), strip.background = element_blank(), strip.placement = "outside", axis.title.y = element_blank()) +
          # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8) , axis.title.y = element_blank()) +
          theme_main +
          facet_wrap(~ factor(set), nrow = 2, scales = "free_y", strip.position = "left") 
        
        #  plot two FDR plots as well: (fraction of negative / negative+positive) and (fraction of permuted / permuted + real)
        # that could give an indication of whether an FDR based on the permuted or based on the negative correlations would be more stringent
        
        # R PEAK GENE #
        
        xlabel = paste0("Correlation coefficient r")
        
        sum_real = table(peakGeneCorrelations.all[indexCur,]$class)[names(dist_class)[1]]
        sum_rnd  = table(peakGeneCorrelations.all[indexCur,]$class)[names(dist_class)[2]]
        binData.r = peakGeneCorrelations.all[indexCur,] %>%
          dplyr::group_by(class) %>%
          dplyr::count(peak_gene.r.class) %>%
          dplyr::mutate(nnorm = dplyr::case_when(class == !! (names(dist_class)[1]) ~ .data$n / (sum_real / sum_rnd), 
                                                 TRUE ~ as.double(.data$n)))
        
        xlabel = paste0("Correlation coefficient r (binned)")
        
        gD = ggplot(binData.r, aes(peak_gene.r.class, nnorm, group = class, fill = class)) + 
          geom_bar(stat = "identity", position = position_dodge(preserve = "single"), na.rm = FALSE, width = 0.5) +
          geom_line(aes(peak_gene.r.class, nnorm, group = class, color= class), stat = "identity") +
          scale_fill_manual("Group", labels = names(dist_class), values = dist_class) +
          scale_color_manual("Group", labels = names(dist_class), values = dist_class) +
          scale_x_discrete(labels = xlabels_peakGene_r.class2, drop = FALSE) +
          theme_bw() + theme(legend.position = "none") +
          xlab(xlabel) + ylab("Abundance") +
          theme_main  +
          scale_y_continuous(labels = scales::scientific)
        
        
        mainTitle = paste0("Summary QC (TF: ", TFCur, ", gene type: ", paste0(geneTypesSelected, collapse = "+"), ",\n", .prettyNum(range), " bp promoter range)")
        
        plots_all = ( ((gA2 | gB3 ) + 
                         patchwork::plot_layout(widths = c(2.5,1.5))) / ((gD) + 
                                                                           patchwork::plot_layout(widths = c(4))) ) + 
          patchwork::plot_layout(heights = c(2,1), guides = 'collect') +
          patchwork::plot_annotation(title = mainTitle, theme = theme(plot.title = element_text(hjust = 0.5)))
        
        plot(plots_all)
        
        
      } # end for each TF
      
      
      
      if (!is.null(GRN@annotation$consensusPeaks_obj) & is.installed("ChIPseeker")) {
        #plot(ChIPseeker::plotAnnoBar(GRN@annotation$consensusPeaks_obj))
        
        # no plot, as this is somehow just a list and no ggplot object
        ChIPseeker::plotAnnoPie(GRN@annotation$consensusPeaks_obj, 
                                main = paste0("\nPeak annotation BEFORE filtering (n = ", nrow(GRN@annotation$consensusPeaks), ")"))
        plot(ChIPseeker::plotDistToTSS(GRN@annotation$consensusPeaks_obj))
      }
      
      # PEAK AND GENE PROPERTIES #
      
      xlabel = paste0("Correlation raw\np-value (binned)")
      mainTitle = paste0("Summary QC (TF: ", TFCur, ", gene type: ", paste0(geneTypesSelected, collapse = "+"), ", ", .prettyNum(range), " bp promoter range)")
      
      allVars = c("peak.GC.class",  "peak.width", "peak.mean","peak.median",
                  "peak.CV", "gene.median",  "gene.mean", "gene.CV", "peak.gene.combined.CV")
      
      # If ChIPseeker is not installed, remove peak.annotation from it
      if ("peak.annotation" %in% colnames(peakGeneCorrelations.all)) {
        allVars = c(allVars, "peak.annotation")
      }
      
      
      for (varCur in allVars) {
        #next
        # Save memory and prune the table and add only the variable we need here
        if (varCur != "peak.gene.combined.CV") {
          dataCur = peakGeneCorrelations.all[indexCurReal,] %>%
            dplyr::select(peak_gene.p_raw, tidyselect::all_of(varCur), class, r_positive, peak_gene.p.raw.class, peak_gene.distance) 
        } else {
          dataCur = peakGeneCorrelations.all[indexCurReal,] %>%
            dplyr::select(peak_gene.p_raw, class, gene.CV, peak.CV, r_positive, peak_gene.p.raw.class, peak_gene.distance)
        }
        
        
        # Choose colors depending on type of variable: gradient or not?
        
        if (varCur %in% c("peak.annotation","peak.GC.class")) {
          
          newColName = varCur
          
        } else {
          
          newColName = paste0(varCur, ".class")
          
          if (varCur == "peak.gene.combined.CV") {
            
            dataCur = dataCur %>%
              dplyr::mutate(!!(newColName) := dplyr::case_when(
                gene.CV < 0.5                & peak.CV < 0.5                ~ "gene.CV+peak.CV<0.5", ##
                # gene.CV >= 0.5 & gene.CV < 1 & peak.CV < 0.5                ~ "gene.CV<1_peak.CV<0.5",
                # gene.CV >= 1                 & peak.CV < 0.5                ~ "gene.CV>1_peak.CV<0.5",
                # 
                # gene.CV < 0.5                & peak.CV >= 0.5 & peak.CV < 1 ~ "gene.CV<0.5_peak.CV<1",
                # gene.CV >= 0.5 & gene.CV < 1 & peak.CV >= 0.5 & peak.CV < 1 ~ "gene.CV<1_peak.CV<1",
                # gene.CV >= 1                 & peak.CV >= 0.5 & peak.CV < 1 ~ "gene.CV>1_peak.CV<1",
                # 
                # gene.CV < 0.5                & peak.CV >= 1                 ~ "gene.CV<0.5_peak.CV>1",
                # gene.CV >= 0.5 & gene.CV < 1 & peak.CV >= 1                 ~ "gene.CV<1_peak.CV>1",
                gene.CV >= 1                 & peak.CV >= 1                 ~ "gene.CV+peak.CV>1", ##
                TRUE ~ "other"
              )) %>%
              dplyr::select(-gene.CV, -peak.CV)
            
          } else {
            
            dataCur = dataCur %>%
              dplyr::mutate(!!(newColName) := cut(.data[[varCur]], breaks = unique(quantile(.data[[varCur]], probs = probs, na.rm = TRUE)), 
                                                  include.lowest = TRUE, ordered_result = TRUE)) %>%
              dplyr::select(-tidyselect::all_of(varCur))
          }
          
          
        }
        
        # Filter groups with fewer than 100 observations
        nGroupsMin = 100
        dataCur = dataCur %>%
          dplyr::group_by(.data[[newColName]]) %>%  
          dplyr::filter(dplyr::n() >= nGroupsMin) %>%
          dplyr::ungroup()
        
        var.label = .prettyNum(table(dataCur[, newColName]))
        var.label = paste0(names(var.label), "\n(n=", var.label , ")\n")
        
        # Set colors
        if (varCur != "peak.annotation") {
          mycolors <- viridis::viridis(length(var.label))
        } else {
          # Only here, we want to have colors that are not a gradient
          mycolors = var.label
          downstream  = which(grepl("downstream", var.label, ignore.case = TRUE))
          promoter    = which(grepl("promoter", var.label, ignore.case = TRUE))
          
          # Remove the first, white-like color from the 2 palettes
          mycolors[downstream] = RColorBrewer::brewer.pal(length(downstream) + 1, "Greens")[-1]
          mycolors[promoter]   = RColorBrewer::brewer.pal(length(promoter) + 1, "Purples")[-1]
          mycolors[which(grepl("3' UTR", mycolors))] = "yellow"
          mycolors[which(grepl("5' UTR", mycolors))] = "orange"
          mycolors[which(grepl("Distal Intergenic", mycolors))] = "red"
          mycolors[which(grepl("Exon", mycolors))] = "maroon"
          mycolors[which(grepl("Intron", mycolors))] = "lightblue"
          
        }
        
        r_pos_tbl = dataCur %>% dplyr::group_by(r_positive) %>% dplyr::pull(r_positive) %>% table()
        r_positive_label = c("TRUE"   = paste0("r+(r>0, n=", .prettyNum(r_pos_tbl[["TRUE"]]), ")"), 
                             "FALSE" = paste0("r-(r<=0, n=" ,.prettyNum(r_pos_tbl[["FALSE"]]), ")"))
        
        # Class in the facet_wrap has been removed, as this is only for real data here
        gA3 = ggplot(dataCur, aes(peak_gene.p_raw, color = .data[[newColName]])) + geom_density(size = 0.2)  +
          facet_wrap(~r_positive, labeller = labeller(r_positive = r_positive_label), nrow = 2) +
          geom_density(aes(color = classNew), color = "black",  linetype = "dotted", alpha = 1) + 
          xlab(xlabel) + ylab("Density (for real data only)") +  theme_bw() +
          scale_color_manual(newColName, values = mycolors, labels = var.label, drop = FALSE ) +
          theme(axis.text=element_text(size=rel(0.6)), strip.text.x = element_text(size = rel(0.8)), 
                legend.text=element_text(size=rel(0.6)), legend.position = "none")
        
        # The following causes an error with patchwork
        # gA3 = ggplot(dataCur, aes(peak_gene.p_raw, color = .data[[newColName]])) + geom_density(size = 0.2)  +
        #   facet_wrap(~ class + r_positive, labeller = labeller(class=freq_class, r_positive = r_positive_label), nrow = 2) +
        #   geom_density(aes(color = classNew), color = "black",  linetype = "dotted", alpha = 1) + 
        #   xlab(xlabel) + ylab("Density") +  theme_bw() +
        #   scale_color_manual(newColName, values = mycolors, labels = var.label, drop = FALSE ) +
        #   theme(axis.text=element_text(size=rel(0.6)), strip.text.x = element_text(size = rel(0.8)), 
        #         legend.text=element_text(size=rel(0.6)), legend.position = "none")
        
        
        # Ratios for r+ / r-
        freq =  dataCur %>%
          dplyr::group_by(class, .data[[newColName]], peak_gene.p.raw.class, r_positive) %>%
          dplyr::summarise(n = dplyr::n()) %>%
          dplyr::ungroup() %>%
          tidyr::complete(class, .data[[newColName]], peak_gene.p.raw.class, r_positive, fill = list(n = 0)) %>% # SOme cases might be missing
          dplyr::group_by(class, .data[[newColName]], peak_gene.p.raw.class) %>% # dont group by r_positive because we want to calculate the ratio within each group
          dplyr::mutate(  
            n = .data$n + 1, # to allow ratios to be computed even for 0 counts
            ratio_pos_raw = .data$n[r_positive] / .data$n[!r_positive]) %>%
          dplyr::filter(r_positive, class == names(dist_class)[1])# Keep only one r_positive row per grouping as we operate via the ratio and this data is duplicated otherwise. Remove random data also because these have been filtered out before are only back due to the complete invokation.
        
        # Cap ratios > 10 at 10 to avoid visual artefacts
        freq$ratio_pos_raw[which(freq$ratio_pos_raw > 10)] = 10
        
        # Without proper colors for now, this will be added after the next plot
        gB3 = ggplot(freq, aes(peak_gene.p.raw.class, ratio_pos_raw, fill = .data[[newColName]])) + 
          geom_bar(stat = "identity", position="dodge", na.rm = TRUE, width = 0.8) + 
          geom_hline(yintercept = 1, linetype = "dotted") + 
          xlab(xlabel) + ylab("Ratio r+ / r- (capped at 10)") +
          scale_x_discrete(labels = xlabels_peakGene_praw.class) +
          scale_fill_manual (varCur, values = mycolors, labels = var.label, drop = FALSE )  +
          theme_bw() +  
          #theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), strip.background = element_blank(), strip.placement = "outside", axis.title.y = element_blank()) +
          # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8) , axis.title.y = element_blank()) +
          theme_main +
          facet_wrap(~ factor(class), nrow = 2, scales = "free_y", strip.position = "left", labeller = labeller(class=freq_class)) 
        
        
        plots_all = ( ((gA3 | gB3 ) + 
                         patchwork::plot_layout(widths = c(2.5,1.5))) ) + 
          patchwork::plot_layout(guides = 'collect') +
          patchwork::plot_annotation(title = mainTitle, theme = theme(plot.title = element_text(hjust = 0.5)))
        
        plot(plots_all)
        
        
        # VERSION JUDITH: simplified peak.gene distance + another variable as histogram, no pemuted data
        
        datasetName = ""
        if (!is.null(GRN@config$metadata$name)) {
          datasetName = GRN@config$metadata$name
        }
        
        mainTitle2 = paste0(datasetName, "\nSummary QC (TF: ", TFCur, ", gene type: ", paste0(geneTypesSelected, collapse = "+"), ", ",
                            .prettyNum(range), " bp promoter range, stratified by distance + ", varCur, ")")
        
        # r+ and r-
        binwidth = 0.1
        mycolors <- viridis::viridis(2)
        xlabel = paste0("Correlation raw p-value (", binwidth, " bins)")
        
        dataCur = dataCur %>%
          dplyr::mutate(peak_gene.distance.class250k = factor(dplyr::if_else(peak_gene.distance <= 250000, "<=250k", ">250k"))) %>%
          dplyr::select(-peak_gene.distance)
        
        # closed = "left", boundary = 0, ensure correct numbers. See https://github.com/tidyverse/ggplot2/issues/1739
        # Without boundary = 0, counts are actually wrong
        
        
        nrows_plot = 2
        if (length(unique(dataCur[[newColName]])) > 9) {
          nrows_plot = 3
        }
        
        gA5 = ggplot(dataCur, aes(peak_gene.p_raw, fill = r_positive)) + 
          geom_histogram(binwidth = binwidth, position="dodge", closed = "left", boundary = 0)  +
          facet_wrap(~ peak_gene.distance.class250k  + .data[[newColName]], nrow = nrows_plot, scales = "free_y") +
          xlab(xlabel) + ylab(paste0("Abundance for classes with n>=", nGroupsMin)) +  theme_bw() +
          ggtitle(mainTitle2) + 
          scale_fill_manual("Class for r", values = mycolors, labels = r_positive_label, drop = FALSE ) +
          theme(axis.text=element_text(size=rel(0.6)), strip.text.x = element_text(size = rel(0.6)), 
                legend.text=element_text(size=rel(0.7)))
        
        
        plot(gA5)
        
        
        
      } #end for each variable
      
      
      # DISTANCE FOCUSED #
      
      
      # Here, we focus on distance and exclude distance classes with too few points and create a new subset of the data
      
      # Filter distance classes with too few points
      distance_class_abund = table(peakGeneCorrelations.all[indexCur,]$peak_gene.distance_class_abs)
      indexFilt = which(peakGeneCorrelations.all$peak_gene.distance_class_abs %in% 
                          names(distance_class_abund)[which(distance_class_abund > 50)])
      indexFilt = intersect(indexFilt, indexCur)
      
      if (length(indexFilt) > 0) {
        g = ggplot(peakGeneCorrelations.all[indexFilt,], aes(peak_gene.p_raw, color = peak_gene.distance_class_abs)) + geom_density() + 
          ggtitle(paste0("Density of the raw p-value distributions")) + 
          facet_wrap(~ r_positive, ncol = 2, labeller = labeler_r_pos) + 
          scale_color_viridis_d(labels = .classFreq_label(table(peakGeneCorrelations.all[indexFilt,]$peak_gene.distance_class_abs))) +
          theme_bw()
        plot(g)
        
        
        if (includeRobustCols) {
          g = ggplot(peakGeneCorrelations.all[indexFilt,], aes(peak_gene.p_raw.robust, color = peak_gene.distance_class_abs)) + geom_density() + 
            ggtitle(paste0("Density of the raw p-value distributions\n(stratified by whether r is positive)")) + 
            facet_wrap(~ r_positive, ncol = 2, labeller = labeler_r_pos) + 
            scale_color_viridis_d(labels = .classFreq_label(table(peakGeneCorrelations.all[indexFilt,]$peak_gene.distance_class_abs))) +
            theme_bw()
          plot(g)
          
          
        }
        
      } # end if (length(indexFilt) > 0)
      
      
      
      
      if (length(indexFilt) > 0) {
        g = ggplot(peakGeneCorrelations.all[indexFilt,], aes(peak_gene.p_raw, color = r_positive)) + 
          geom_density() + 
          ggtitle(paste0("Density of the raw p-value distributions")) + 
          facet_wrap(~ peak_gene.distance_class_abs,  ncol = 2, labeller = .customLabeler(table(peakGeneCorrelations.all$peak_gene.distance_class_abs))) +
          scale_color_manual(labels = .classFreq_label(table(peakGeneCorrelations.all[indexFilt,]$r_positive)), values = r_pos_class) +
          theme_bw()
        plot(g)
        
        
      }
      
      
      #
      # Focus on peak_gene.r
      #
      
      g = ggplot(peakGeneCorrelations.all[indexCur,], aes(peak_gene.r, color = peak_gene.distance_class_abs)) + geom_density() + 
        geom_density(data = peakGeneCorrelations.all[indexCur,], aes(peak_gene.r), color = "black") + 
        ggtitle(paste0("Density of the correlation coefficients")) + 
        scale_color_viridis_d(labels = .classFreq_label(table(peakGeneCorrelations.all[indexCur,]$peak_gene.distance_class_abs))) +
        theme_bw()
      plot(g)
      
      
      if (!is.null(fileBase)) {
        grDevices::dev.off()
      }
      
      
    }
    
  } else {
    futile.logger::flog.info(paste0("All output files already exist. Set forceRerun = TRUE to regenerate and overwrite."))
    
  }
  
  
  
  .printExecutionTime(start)
  
}

.createTables_peakGeneQC <- function(peakGeneCorrelations.all.cur, networkType_details, colors_vec, range) {
  
  d = peakGeneCorrelations.all.cur %>% 
    dplyr::group_by(r_positive, class, peak_gene.p.raw.class) %>%
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
    dplyr::group_by(r_positive, class) %>%
    dplyr::summarise(sum_n = sum(.data$n))
  
  
  # Summarize per bin
  d2 = d %>%
    dplyr::group_by(class, peak_gene.p.raw.class)%>%
    dplyr::summarise(sum_pos = .data$n[r_positive],
                     sum_neg = .data$n[!r_positive]) %>%
    dplyr::mutate(ratio_pos_raw = sum_pos / sum_neg,
                  enrichment_pos = sum_pos / (sum_pos + sum_neg),
                  ratio_neg = 1 - enrichment_pos) %>%
    dplyr::ungroup()
  
  # Compare between real and permuted
  normFactor_real = dplyr::filter(dsum, class ==  !! (networkType_details[1])) %>%  dplyr::pull(sum_n) %>% sum() /
    dplyr::filter(dsum, class ==  !! (networkType_details[2])) %>%  dplyr::pull(sum_n) %>% sum()
  
  # ratio_norm not used currently, no normalization necessary here or not even useful because we dont want to normalize the r_pos and r_neg ratios: These are signal in a way. Only when comparing between real and permuted, we have to account for sample size for corrections
  d3 = d %>%
    dplyr::group_by(peak_gene.p.raw.class, r_positive) %>%
    dplyr::summarise(n_real     = .data$n[class == !! (names(colors_vec)[1]) ],
                     n_permuted = .data$n[class == !! (names(colors_vec)[2]) ]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(ratio_real_raw = n_real / n_permuted,
                  ratio_real_norm = ratio_real_raw / normFactor_real,
                  enrichment_real      = n_real / (n_real + n_permuted),
                  enrichment_real_norm = (n_real / normFactor_real) / ((n_real / normFactor_real) + n_permuted)) 
  
  
  stopifnot(identical(levels(d2$peak_gene.p.raw.class), levels(d3$peak_gene.p.raw.class)))
  # 2 enrichment bar plots but combined using facet_wrap
  d2$set = "r+ / r-"; d3$set = "real / permuted" 
  d_merged <- tibble::tibble(peak_gene.p.raw.class = c(as.character(d2$peak_gene.p.raw.class), 
                                                       as.character(d3$peak_gene.p.raw.class)),
                             ratio = c(d2$ratio_pos_raw, d3$ratio_real_norm),
                             classAll = c(as.character(d2$class), d3$r_positive),
                             set = c(d2$set, d3$set)) %>%
    dplyr::mutate(classAll = factor(classAll, levels=c(paste0("real_",range), paste0("random_",range), "TRUE", "FALSE")),
                  peak_gene.p.raw.class = factor(peak_gene.p.raw.class, levels = levels(d2$peak_gene.p.raw.class)))
  
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
    dplyr::mutate(ratio_rneg_rpos_real = n_rneg_real / (n_rneg_real + n_rpos_real),
                  ratio_rneg_rpos_random = n_rneg_random / (n_rneg_random + n_rpos_random),
                  peak_gene.p.raw.class.bin = as.numeric(peak_gene.p.raw.class)) %>%
    dplyr::arrange(peak_gene.p.raw.class.bin)
  
  d4_melt = reshape2::melt(d4, id  = c("peak_gene.p.raw.class.bin", "peak_gene.p.raw.class")) %>%
    dplyr::filter(grepl("ratio", variable))
  
  
  list(d = d, d2 = d2, d3 = d3, d4 = d4, d4_melt = d4_melt, d_merged = d_merged)
  
}



.customLabeler <- function(tbl_freq) {
  tbl_freq_label = paste0(names(tbl_freq), " (", tbl_freq, ")")
  names(tbl_freq_label) = names(tbl_freq)
  ggplot2::as_labeller(tbl_freq_label)
}

.classFreq_label <- function(tbl_freq) {
  paste0(names(tbl_freq), " (", tbl_freq, ")")
}


###### Connection summaries ########

#' Plot various network connectivity summaries for a \code{\linkS4class{GRN}} object
#'
#' @template GRN 
#' @param type Character. Either \code{"heatmap"} or \code{"boxplot"}. Default \code{"heatmap"}. Which plot type to produce?
#' @template outputFolder
#' @template basenameOutput
#' @template plotAsPDF
#' @template pdf_width
#' @template pdf_height
#' @template forceRerun
#' @return The same \code{\linkS4class{GRN}} object, without modifications. In addition, for the specified \code{type}, a PDF file (default filename is \code{GRN.connectionSummary_{type}.pdf}) is produced with a connection summary.
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' GRN = plot_stats_connectionSummary(GRN, outputFolder = ".", forceRerun = FALSE, plotAsPDF=FALSE)
#' @export
#' @importFrom circlize colorRamp2
plot_stats_connectionSummary <- function(GRN, type = "heatmap", 
                                         outputFolder = NULL, basenameOutput = NULL, 
                                         plotAsPDF = TRUE, pdf_width = 12, pdf_height = 12, 
                                         forceRerun = FALSE) {
  
  GRN = .addFunctionLogToObject(GRN)
  
  checkmate::assertClass(GRN, "GRN")
  checkmate::assertChoice(type, c("heatmap", "boxplot", "density"))
  checkmate::assertFlag(plotAsPDF)
  checkmate::assertNumeric(pdf_width, lower = 5, upper = 99)
  checkmate::assertNumeric(pdf_height, lower = 5, upper = 99)
  checkmate::assertFlag(forceRerun)
  checkmate::assert(checkmate::checkNull(outputFolder), checkmate::checkCharacter(outputFolder, min.chars = 1))
  checkmate::assert(checkmate::checkNull(basenameOutput), checkmate::checkCharacter(basenameOutput, len = 1, min.chars = 1, any.missing = FALSE))
  
  outputFolder = .checkOutputFolder(GRN, outputFolder)
  
  if (type ==  "heatmap") {
    if (plotAsPDF) {
      file = paste0(outputFolder, dplyr::if_else(is.null(basenameOutput), .getOutputFileName("plot_connectionSummary_heatmap"), basenameOutput), ".pdf")
    } else {
      file = NULL
    }
    
    .plot_stats_connectionSummaryHeatmap(GRN, file = file, pdf_width = pdf_width, pdf_height = pdf_height, forceRerun = forceRerun) 
  
  } else if (type ==  "boxplot") {
    
    if (plotAsPDF) {
      file = paste0(outputFolder, dplyr::if_else(is.null(basenameOutput), .getOutputFileName("plot_connectionSummary_boxplot"), basenameOutput), ".pdf")
    } else {
      file = NULL
    }
    
    .plot_stats_connectionSummaryBoxplot(GRN, file = file, pdf_width = pdf_width, pdf_height = pdf_height, forceRerun = forceRerun) 
    
  } else if (type ==  "density") {
    
    stop("Not yet implemented")
    if (plotAsPDF) {
      file = "TODO"
    } else {
      file = NULL
    }
    .plot_stats_connectionSummaryDensity(GRN, file = file, pdf_width = pdf_width, pdf_height = pdf_height, forceRerun = forceRerun) 
  }
  
  GRN
  
}


.plot_stats_connectionSummaryHeatmap <- function(GRN, file = NULL, pdf_width = 12, pdf_height = 12, forceRerun = FALSE) {
  
  start = Sys.time()
  
  futile.logger::flog.info(paste0("Plotting connection summary", dplyr::if_else(is.null(file), "", paste0(" to file ", file))))
  
  if ((!is.null(file) && !file.exists(file))| forceRerun) {
    
    if (nrow(GRN@stats$connections) == 0) {
      message = paste0("Statistics summary missing from object, please run the function generateStatsSummary first")
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    index = 0
    
    
    if (!is.null(file)) {
      .checkOutputFile(file)
      grDevices::pdf(file = file, width = pdf_width, height = pdf_height)
    }
    
    
    for (allowMissingTFsCur in unique(GRN@stats$connections$allowMissingTFs)) {
      
      for (allowMissingGenesCur in unique(GRN@stats$connections$allowMissingGenes)) {
        
        for (TF_peak.connectionTypeCur in .getAll_TF_peak_connectionTypes(GRN)) {
          
          # Plot a summary page with the parameters used for the next page
          # 
          #index = index + 1
          
          for (elemCur in c("nTFs", "nPeaks", "nGenes")) {
            
            plotData.l = list()
            for (permCur in 0:.getMaxPermutation(GRN)) {
              
              permIndex = as.character(permCur)
              
              stats_filtered.df = GRN@stats$connections %>%
                dplyr::filter(allowMissingGenes == allowMissingGenesCur, 
                              allowMissingTFs == allowMissingTFsCur,
                              TF_peak.connectionType == TF_peak.connectionTypeCur, 
                              is.na(peak_gene.p_raw), 
                              perm == permCur) %>%
                dplyr::select(TF_peak.fdr, peak_gene.fdr, nGenes, nTFs, nPeaks)
              
              # Stratify by TF, peak or gene and produce a simple matrix
              if (nrow(stats_filtered.df) > 0) {
                plotData.l[[permIndex]] = as.matrix(reshape2::dcast(stats_filtered.df, TF_peak.fdr ~ peak_gene.fdr, value.var = elemCur)[,-1])
                
                colnames(plotData.l[[permIndex]]) = paste0("peak_gene\n.fdr_", colnames(plotData.l[[permIndex]]))
                rownames(plotData.l[[permIndex]]) = paste0("TF_peak\n.fdr_", sort(unique(stats_filtered.df$TF_peak.fdr)))
              }
              
              
            }
            
            # Now we have the data for both permutations, let's plot
            
            maxDataRange = max(unlist(lapply(plotData.l,FUN=max))) *1.1
            
            for (permCur in 0:.getMaxPermutation(GRN)) {
              
              index = index + 1
              permIndex = as.character(permCur)
              permSuffix = paste0(dplyr::if_else(permCur == 0, " (real)", " (permuted)"))
              
              titlePlot =  paste0("Number of unique ", gsub("^n", "", elemCur), permSuffix,
                                  "\nallowMissingTFs: ", allowMissingTFsCur, 
                                  ", allowMissingGenes: ", allowMissingGenesCur, 
                                  ",\n TF_peak.connectionType: ", TF_peak.connectionTypeCur)
              
              
              colors = circlize::colorRamp2(c(0, maxDataRange), c("white", "red"))
              
              ComplexHeatmap::Heatmap(
                plotData.l[[permIndex]],
                name = "Number of\nconnections",
                col = colors,
                cluster_columns = FALSE, cluster_rows = FALSE,
                row_names_side = "right", row_names_gp = grid::gpar(fontsize = 10), 
                column_title = titlePlot,
                column_names_gp = grid::gpar(fontsize = 10),
                row_title = "TF-peak FDR",
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid::grid.text(sprintf("%.0f", plotData.l[[permIndex]][i, j]), x, y, gp = grid::gpar(fontsize = 10))
                }
              ) %>% plot()
              
              
              
              
            }
            
            
          } # end for (elemCur in c("nTFs", "nPeaks", "nGenes"))
          
        } # end for (TF_peak.connectionTypeCur in TF_peak.connectionTypesAll)
      }
    }
    
    if (!is.null(file)) {
      grDevices::dev.off()
    }
    
    
  }
  
  .printExecutionTime(start)
  
}

.plot_stats_connectionSummaryDensity <- function(GRN, file = file, TF_peak.connectionType, allowMissingTFs = FALSE, allowMissingGenes = FALSE, pdf_width = 12, pdf_height = 12, forceRerun = FALSE) {
  
  # TODO
  # Same for both permutations, in one plot
  main.df = GRN@connections$TF_peaks[["0"]]$main
  main.filt.df = dplyr::filter(main.df, TF_peak.fdr < 0.3)
  ggplot(main.df, aes(TF_peak.fdr)) + geom_density() + theme_bw()
  ggplot(main.df, aes(TF_peak.fdr)) + geom_histogram(binwidth = 0.01) +  
    geom_vline(xintercept = c(0.001, 0.01, 0.05, 0.1, 0.2, 0.3), linetype = "dotted", color = "darkgray") + theme_bw()
  ggplot(main.filt.df, aes(TF_peak.fdr)) + geom_histogram(binwidth = 0.01) +
    geom_vline(xintercept = c(0.001, 0.01, 0.05, 0.1, 0.2, 0.3), linetype = "dotted", color = "darkgray") + theme_bw()
  
  
  # Connectivity of TF-peaks in dependency of TF_peak.fdr
  
  
}


# Compare real and permuted one, replicate Darias plots
.plot_stats_connectionSummaryBoxplot <- function(GRN, file = NULL, pdf_width = 12, pdf_height = 12, forceRerun = FALSE) {
  
  start = Sys.time()
  
  if ((!is.null(file) && !file.exists(file))| forceRerun) {
    
    
    futile.logger::flog.info(paste0("Plotting diagnostic plots for network connections", dplyr::if_else(is.null(file), "", paste0(" to file ", file))))
    
    
    if (!is.null(file)) {
      .checkOutputFile(file)
      grDevices::pdf(file = file, width = pdf_width, height = pdf_height)
    }
    
    stats_details.l = GRN@stats$connectionDetails.l
    
    TF_peak.fdrs      = names(GRN@stats$connectionDetails.l[["0"]])
    peak_gene.fdrs    = names(GRN@stats$connectionDetails.l[["0"]][[TF_peak.fdrs[1]]])
    allowMissingTFs   = names(GRN@stats$connectionDetails.l[["0"]][[TF_peak.fdrs[1]]][[peak_gene.fdrs[1]]])
    allowMissingGenes = names(GRN@stats$connectionDetails.l[["0"]][[TF_peak.fdrs[1]]][[peak_gene.fdrs[1]]][[allowMissingTFs[1]]])
    
    
    nElemsTotal = length(TF_peak.fdrs) * length(peak_gene.fdrs) * length(allowMissingTFs) * length(allowMissingGenes) * length(.getAll_TF_peak_connectionTypes(GRN))
    
    pb <- progress::progress_bar$new(total = nElemsTotal)
    
    
    
    for (TF_peak.fdr_cur in rev(TF_peak.fdrs)) {
      
      for (peak_gene.fdr_cur in peak_gene.fdrs) {
        
        for (allowMissingTFsCur in allowMissingTFs) {
          
          for (allowMissingGenesCur in allowMissingGenes) {
            
            for (TF_peak.connectionTypeCur in .getAll_TF_peak_connectionTypes(GRN)) {
              
              pb$tick()
              elemCur = "TF"
              dataCur0.df = as.data.frame(stats_details.l[["0"]][[TF_peak.fdr_cur]][[peak_gene.fdr_cur]] [[allowMissingTFsCur]] [[allowMissingGenesCur]] [[TF_peak.connectionTypeCur]] [[elemCur]])
              dataCur1.df = as.data.frame(stats_details.l[["1"]][[TF_peak.fdr_cur]][[peak_gene.fdr_cur]] [[allowMissingTFsCur]] [[allowMissingGenesCur]] [[TF_peak.connectionTypeCur]] [[elemCur]])
              
              if (nrow(dataCur0.df) == 0) {
                dataCur0.df = data.frame(. = "MISSING", Freq = 0)
              }
              if (nrow(dataCur1.df) == 0) {
                dataCur1.df = data.frame(. = "MISSING", Freq = 0)
              }
              
              dataCur.df = rbind(dplyr::mutate(dataCur0.df, networkType = "real", connectionType = elemCur), 
                                 dplyr::mutate(dataCur1.df, networkType = "permuted", connectionType = elemCur)) %>% 
                tibble::as_tibble()
              
              for (elemCur in c("gene","peak.gene","peak.TF")) {
                
                dataCur0.df = as.data.frame(stats_details.l[["0"]][[TF_peak.fdr_cur]][[peak_gene.fdr_cur]] [[allowMissingTFsCur]] [[allowMissingGenesCur]] [[TF_peak.connectionTypeCur]] [[elemCur]])
                dataCur1.df = as.data.frame(stats_details.l[["1"]][[TF_peak.fdr_cur]][[peak_gene.fdr_cur]] [[allowMissingTFsCur]] [[allowMissingGenesCur]] [[TF_peak.connectionTypeCur]] [[elemCur]])
                
                if (nrow(dataCur0.df) == 0) {
                  dataCur0.df = data.frame(. = "MISSING", Freq = 0)
                }
                if (nrow(dataCur1.df) == 0) {
                  dataCur1.df = data.frame(. = "MISSING", Freq = 0)
                }
                
                dataCur.df = rbind(dataCur.df,
                                   dplyr::mutate(dataCur0.df, networkType = "real", connectionType = elemCur), 
                                   dplyr::mutate(dataCur1.df, networkType = "permuted", connectionType = elemCur))
                
              }
              
              dataCur.df$networkType =as.factor(dataCur.df$networkType)
              
              titlePlot =  paste0("TF_peak.fdr: ", TF_peak.fdr_cur, 
                                  ", peak_gene.fdr: ", peak_gene.fdr_cur, 
                                  "\nallowMissingTFs: ", allowMissingTFsCur, 
                                  ", allowMissingGenes: ", allowMissingGenesCur, 
                                  ",\n TF_peak.connectionType: ", TF_peak.connectionTypeCur)
              
              if (max(dataCur.df$Freq)> 0) {
                
                g = ggplot(dataCur.df, aes(networkType, log10(Freq), fill = networkType)) + geom_boxplot() +theme_bw() + 
                  xlab("Network type")  +ylab(paste0("log10(Number of connections per category",  ", empty = 0)")) +
                  ggtitle(titlePlot) + 
                  facet_wrap(~connectionType, scales = "free") + 
                  scale_fill_manual("Network type", values = c("real" =  "red", "permuted" = "gray"), 
                                    labels = c("real" =  "real", "permuted" = "permuted"), drop = FALSE)
                
                suppressWarnings(plot(g))
                
              } else {
                
                plot(c(0, 1), c(0, 1), ann = FALSE, bty = 'n', type = 'n', axes=FALSE, main = titlePlot)
                message = paste0(titlePlot, "\nNo data to show")
                text(x = 0.5, y = 0.5, message, cex = 1.2, col = "red")
              }
              
            } #end  for (TF_peak.connectionTypeCur in TF_peak.connectionTypesAll)
            
          } # end  for (allowMissingGenesCur in c(FALSE, TRUE))
          
        } # end for (allowMissingTFsCur in c(FALSE, TRUE))
        
      } # end for each peak_gene.fdr_cur
      
    } # end for TF_peak.fdr_cur
    
    if (!is.null(file)) {
      grDevices::dev.off()
    }
  }
  
  
  .printExecutionTime(start)
  
}





######## General & communities Stats Functions ########


#' Plot general structure and connectivity statistics for a filtered \code{\linkS4class{GRN}}
#' 
#' This function generates graphical summaries about the structure and connectivity of the TF-peak-gene and TF-gene graphs. These include, distribution of vertex types (TF, peak, gene) and edge types (tf-peak, peak-gene), the distribution of vertex degrees, and the most "important" vertices according to degree centrality and eigenvector centrality scores.
#' @template GRN
#' @template outputFolder
#' @template basenameOutput
#' @template forceRerun
#' @template plotAsPDF
#' @template pdf_width
#' @template pdf_height
#' @return The same \code{\linkS4class{GRN}} object with no changes. The results are output to a file.
#' @seealso \code{\link{plotGeneralEnrichment}}
#' @seealso \code{\link{plotCommunitiesStats}}
#' @seealso \code{\link{plotCommunitiesEnrichment}}
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' GRN = plotGeneralGraphStats(GRN, outputFolder = ".", plotAsPDF=FALSE)
#' @export
plotGeneralGraphStats <- function(GRN, outputFolder = NULL, basenameOutput = NULL, 
                                  plotAsPDF = TRUE, pdf_width = 12, pdf_height = 12, 
                                  forceRerun = FALSE) {
  
  start = Sys.time()
  GRN = .addFunctionLogToObject(GRN)
  
  checkmate::assertClass(GRN, "GRN")
  checkmate::assertFlag(plotAsPDF)
  checkmate::assertNumeric(pdf_width, lower = 5, upper = 99)
  checkmate::assertNumeric(pdf_height, lower = 5, upper = 99)
  checkmate::assertFlag(forceRerun)
  checkmate::assert(checkmate::checkNull(outputFolder), checkmate::checkCharacter(outputFolder, min.chars = 1))
  checkmate::assert(checkmate::checkNull(basenameOutput), checkmate::checkCharacter(basenameOutput, len = 1, min.chars = 1, any.missing = FALSE))
  
  outputFolder = .checkOutputFolder(GRN, outputFolder)
  
  fileCur = paste0(outputFolder, dplyr::if_else(is.null(basenameOutput), .getOutputFileName("plot_generalNetworkStats"), basenameOutput), ".pdf")
  if (!file.exists(fileCur) | forceRerun | !plotAsPDF) {
    
    if (plotAsPDF) {
      .checkOutputFile(fileCur)
      grDevices::pdf(fileCur, width = pdf_width, height = pdf_height)
      futile.logger::flog.info(paste0("Plotting to file ", fileCur))
    } else {
      futile.logger::flog.info(paste0("Plotting directly"))   
    }
    
    TF_peak_gene.df = GRN@graph$TF_peak_gene$table
    TF_gene.df = GRN@graph$TF_gene$table
    
    # pie charts
    totalVerteces = data.frame(Class = c("TF", "Peak", "Gene"),
                               Count = c(length(unique(stats::na.omit(GRN@connections$all.filtered$`0`$TF.name))),
                                         length(unique(stats::na.omit(GRN@connections$all.filtered$`0`$peak.ID))),
                                         length(unique(stats::na.omit(GRN@connections$all.filtered$`0`$gene.ENSEMBL))) ))
    
    theme_plots = theme(axis.ticks = element_blank(), axis.title = element_blank(), axis.text = element_blank(), 
                        legend.title = element_blank(), panel.background = element_rect(fill = "white"),
                        plot.title = element_text(hjust = 0.5))
    
    
    gVertexDist = ggplot(totalVerteces, aes(x="", y=Count, fill=Class)) + geom_bar(stat="identity") +
      coord_polar("y", start=0) +
      scale_fill_manual(values = c("#3B9AB2", "#F21A00", "#E1AF00")) +
      geom_text(aes(label = paste0(round(Count/sum(Count) *100, digits = 1), "%")),
                position = position_stack(vjust = 0.5)) +
      theme_plots +
      ggtitle(paste0("Vertices (n=", sum(totalVerteces$Count), ")"))
    
    geneDist = as.data.frame(table(droplevels(GRN@connections$all.filtered$`0`$gene.type)))
    if (nrow(geneDist) == 0) {
      gGeneDist = 
        ggplot() + 
        annotate("text", x = 4, y = 25, size=5, color = "red", label = "Nothing to plot:\nNo genes found") + 
        theme_void()
      
      
      message = "No genes found in the GRN object. Make sure the filtered connections contain also genes."
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
      
    } else {
      gGeneDist = ggplot(geneDist, aes(x= "", y = Freq, fill = Var1)) + geom_bar(stat="identity") +
        coord_polar("y") +
        geom_text(aes(label = paste0(round(Freq/sum(Freq) *100, digits = 1 ), "%")),
                  position = position_stack(vjust = 0.5)) +
        theme_plots +
        ggtitle(paste0("Genes (n=", sum(geneDist$Freq), ")"))
    }
    
    
    totalEdges = as.data.frame(table(TF_peak_gene.df$connectionType))
    if (nrow(totalEdges) == 0) {
      gEdgeDist= ggplot() + 
        annotate("text", x = 4, y = 25, size=5, color = "red", label = "Nothing to plot:\nNo genes found") + 
        theme_void()
    } else {
      gEdgeDist = ggplot(totalEdges, aes(x="", y=Freq, fill=Var1)) + geom_bar(stat="identity") +
        coord_polar("y", start=0) +
        scale_fill_manual(values = c("#BDC367", "#6BB1C1")) +
        geom_text(aes(label = paste0(round(Freq/sum(Freq) *100, digits = 1 ), "%")),
                  position = position_stack(vjust = 0.5)) +
        theme_plots +
        ggtitle(paste0("Edges (n=", sum(totalEdges$Freq), ")"))
    }
    
    # Page 1: Pie charts of the connections
    print((gVertexDist + gEdgeDist)/gGeneDist + patchwork::plot_layout(widths = c(2,1), guides = "collect"))
    
    # First, we focus on the TF-peak-gene graph
    # Get degree stats and central vertexes in the TF-peak-gene graph
    TF_peak_gene.degree.stats = .getDegreeStats(GRN, TF_peak_gene.df)
    suffix = paste0(" in the filtered TF-peak-gene eGRN")
    if (nrow(TF_peak_gene.df) > 0) {
      #degreeDist.tbl = TF_peak_gene.degree.stats$tbl$degrees
      gDegrees  = TF_peak_gene.degree.stats$figures$degreeDist + ggtitle(paste0("Distribution of vertex degrees", suffix))
      gTopGenes = TF_peak_gene.degree.stats$figures$topGenes   + ggtitle(paste0("Top degree-central genes", suffix))
      gTopTFs   = TF_peak_gene.degree.stats$figures$topTFs     + ggtitle(paste0("Top degree-central TFs", suffix))
      
      eigen_stats = .getEigenCentralVertices(GRN, "TF_peak_gene")
      gTopEigenGenes = eigen_stats$topGenes + ggtitle(paste0("Top eigenvector-central genes", suffix))
      gTopEigenTFs   = eigen_stats$topTFs   + ggtitle(paste0("Top eigenvector-central TFs", suffix))
      
      # Pages 2 to 4: TF-peak-gene network: Degree distribution, and top genes for 2 different measures
      print(gDegrees)
      print(gTopGenes/gTopTFs)
      print(gTopEigenGenes/gTopEigenTFs)
      
      ##
      # Now we focus on the TF-gene graph only, not the TF-peak-gene one
      ###
      suffix = paste0(" in the filtered TF-gene eGRN")
      # Get degree stats and central vertexes in the filtered GRN
      TF_gene.degree.stats = .getDegreeStats(GRN, TF_gene.df)
      degreeDist.tbl = TF_gene.degree.stats$tbl$degrees
      gDegrees = TF_gene.degree.stats$figures$degreeDist + ggtitle(paste0("Distribution of vertex degrees", suffix))
      gTopGenes = TF_gene.degree.stats$figures$topGenes + ggtitle(paste0("Top degree-central genes", suffix))
      gTopTFs = TF_gene.degree.stats$figures$topTFs + ggtitle(paste0("Top degree-central TFs", suffix))
      
      # get top eigenvalue central TFs and genes
      eigen_stats = .getEigenCentralVertices(GRN, "TF_gene")
      gTopEigenGenes = eigen_stats$topGenes + ggtitle(paste0("Top eigenvector-central genes", suffix))
      gTopEigenTFs = eigen_stats$topTFs + ggtitle(paste0("Top eigenvector-central TFs", suffix))
      
      # Pages 5 to 7: TF-gene network: Degree distribution, and top genes for 2 different measures
      print(gDegrees)
      print(gTopGenes/gTopTFs)
      print(gTopEigenGenes/gTopEigenTFs)
    }
    
    
    if (plotAsPDF) {grDevices::dev.off()}
    
  } else {
    futile.logger::flog.info(paste0("File ", fileCur, " already exists, not overwriting due to forceRerun = FALSE"))
  }
  
  .printExecutionTime(start)
  
  GRN
  
}



#' Plot the general enrichement results
#' 
#' This function plots the results of the general enrichment analysis for every specified ontology.
#' 
#' @template GRN
#' @template outputFolder
#' @template basenameOutput
#' @template plotAsPDF
#' @template pdf_width
#' @template pdf_height
#' @template forceRerun
#' @param ontology Character. \code{NULL} or vector of ontology names. Default \code{NULL}. Vector of ontologies to plot. The results must have been previously calculated otherwise an error is thrown.
#' @param p Numeric. Default 0.05. p-value threshold to determine significance.
#' @param topn_pvalue Numeric. Default 30. Maximum number of ontology terms that meet the p-value significance threshold to display in the enrichment dot plot
#' @param display_pAdj \code{TRUE} or \code{FALSE}. Default \code{FALSE}. Is the p-value being displayed in the plots the adjusted p-value? This parameter is relevant for KEGG, Disease Ontology, and Reactome enrichments, and does not affect GO enrichments.
#' @return The same \code{\linkS4class{GRN}} object, without modifications. A single PDF file is produced with the results.
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' GRN = plotGeneralEnrichment(GRN, outputFolder = ".", plotAsPDF=FALSE)
#' @export
plotGeneralEnrichment <- function(GRN, outputFolder = NULL, basenameOutput = NULL, 
                                  ontology = NULL, topn_pvalue = 30, p = 0.05, 
                                  display_pAdj = FALSE, 
                                  plotAsPDF = TRUE, pdf_width = 12, pdf_height = 12, 
                                  forceRerun = FALSE) {
  
  start = Sys.time()
  GRN = .addFunctionLogToObject(GRN)
  
  checkmate::assertClass(GRN, "GRN")
  
  ontologiesFound = names(GRN@stats$Enrichment$general)
  
  checkmate::assert(checkmate::checkNull(outputFolder), checkmate::checkCharacter(outputFolder, min.chars = 1))
  checkmate::assert(checkmate::checkNull(basenameOutput), checkmate::checkCharacter(basenameOutput, len = 1, min.chars = 1, any.missing = FALSE))
  
  checkmate::assert(checkmate::checkNull(ontology), checkmate::checkSubset(ontology, ontologiesFound))
  checkmate::assertNumeric(p, lower = 0, upper = 1)
  checkmate::assertIntegerish(topn_pvalue, lower = 1)
  checkmate::assertFlag(display_pAdj)
  checkmate::assertFlag(plotAsPDF)
  checkmate::assertNumeric(pdf_width, lower = 5, upper = 99)
  checkmate::assertNumeric(pdf_height, lower = 5, upper = 99)
  checkmate::assertFlag(forceRerun)
  
  outputFolder = .checkOutputFolder(GRN, outputFolder)
  
  if (length(ontologiesFound) == 0) {
    message = "No ontologies found in GRN object for general enrichment. Run the function calculateGeneralEnrichment first"
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
  fileCur = paste0(outputFolder, dplyr::if_else(is.null(basenameOutput), .getOutputFileName("plot_generalEnrichment"), basenameOutput), ".pdf")
  if (!file.exists(fileCur) | forceRerun | !plotAsPDF) {
    
    if (plotAsPDF) {
      .checkOutputFile(fileCur)
      grDevices::pdf(fileCur, width = pdf_width, height = pdf_height)
      futile.logger::flog.info(paste0("Plotting to file ", fileCur))
    } else {
      futile.logger::flog.info(paste0("Plotting directly"))
    }
    
    ontologies =  ontologiesFound
    if (!is.null(ontology)) {
      ontologies = ontology
    } 
    
    futile.logger::flog.info(paste0("Found the following ontology results for plotting: ", paste0(ontologiesFound, collapse = ",")))
    futile.logger::flog.info(paste0("Plotting for the following user-selected ontologies: ", paste0(ontologies, collapse = ",")))
    
    for (ontologyCur in ontologies) {
      
      futile.logger::flog.info(paste0(" Ontology ", ontologyCur))   
      .plotEnrichmentGeneral(data = GRN@stats$Enrichment$general[[ontologyCur]], 
                             type = ontologyCur, 
                             prefixTitle = "General Enrichment Analysis",
                             topn_pvalue, p = p, maxWidth_nchar_plot = 100, display_pAdj = display_pAdj)
    }
    
    if (plotAsPDF) grDevices::dev.off()
    
  } else {
    futile.logger::flog.info(paste0("File ", fileCur, " already exists, not overwriting due to forceRerun = FALSE"))
  }
  
  .printExecutionTime(start)
  
  GRN
}


.plotEnrichmentGeneral <- function(data, type, prefixTitle, topn_pvalue = 30, maxWidth_nchar_plot = 100, p = 0.05, display_pAdj = FALSE) {
  
  dataCur   = data[["results"]]
  paramsCur = data[["parameters"]]
  
  parameterStr = paste(names(paramsCur), paramsCur, sep = ":", collapse = ", ")
  # Add extra line space
  parameterStr = gsub(", nBackground:", ",\nnBackground:", parameterStr, fixed = TRUE)
  titleCur = paste0(prefixTitle, " for ontology ", type, "\n(", parameterStr, ")") 
  
  pValuePrefix = "Raw "
  if (display_pAdj & !stringr::str_starts(type, "GO_")) {
    dataCur = dataCur %>%
      dplyr::mutate(pval = .data$p.adjust)
    
    pValuePrefix = "Adj."
  }
  
  # if (!is.data.frame(dataCur)) {
  #     
  #     plot(c(0, 1), c(0, 1), ann = FALSE, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n', main = title)
  #     
  #     message = paste0(titleCur, "\n\nThe enrichment calculation previously returned an error.")  
  #     message = paste0(message, "\n\nYou may run it again in case the problem was temporary\nor enrichment calculation is not possible to calculate here.")
  #     
  #     text(x = 0.5, y = 0.5, message, cex = 1.6, col = "red")
  #     
  #     return()
  # } 
  # 
  
  dataCur = dataCur %>%
    dplyr::mutate(dataCur, Term.mod = stringr::str_trunc(as.character(Term), width = maxWidth_nchar_plot, side = "right")) %>%
    dplyr::mutate(pval = as.numeric(gsub(">|<", "", pval))) %>%
    dplyr::arrange(pval) %>%
    dplyr::filter(pval <= p) # filter out non significant results and of the remaining take the top n
  
  legendsTitle = paste0(pValuePrefix, " p-value")
  labelAxis = paste0("(Truncated) ", type, " Term")
  if (nrow(dataCur) > topn_pvalue) {
    labelAxis = paste0(labelAxis, " (top ", topn_pvalue, " terms only)")
  }  
  
  dataCur = dataCur %>%
    dplyr::slice(seq_len(topn_pvalue)) %>%
    dplyr::arrange(GeneRatio)
  #dplyr::top_n(n = topn_pvalue, wt = -pval)
  
  
  
  if (nrow(dataCur) > 0) {
    averagePvalue = mean(dataCur$pval, na.rm = TRUE)
    
    g = ggplot(dataCur , aes(x = stats::reorder(Term, GeneRatio), y = GeneRatio, colour = pval, size = Found) ) +
      geom_point(na.rm = TRUE) +
      scale_x_discrete(labelAxis, labels= dataCur$Term.mod) +
      scale_color_gradient2(legendsTitle, midpoint = averagePvalue, low = "#F21A00" , mid = "#EBCC2A", high = "#3B9AB2") +
      scale_size_continuous("Number of\ngenes found (n)") + 
      ggtitle(titleCur) +
      ylab(paste0("Gene ratio [n / total foreground (", paramsCur$nForeground, ")]")) +
      coord_flip() +
      theme_bw() +
      theme(plot.title = element_text(size=rel(0.9)))
    
    plot(g)
  } else {
    
    plot(c(0, 1), c(0, 1), ann = FALSE, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n', main = title)
    
    message = paste0(titleCur, "\n\nNo enrichment found.")  
    if (pValuePrefix == "Adj.") {
      message = paste0(message, "\n\nYou may run again based on\nraw p-value and not adjusted one.")
    }  
    
    text(x = 0.5, y = 0.5, message, cex = 1.6, col = "red")
    
  }
  
  
  
}


#' Plot general structure & connectivity statistics for each community in a filtered \code{\linkS4class{GRN}}
#' 
#' Similarly to the statistics produced by \code{\link{plotGeneralGraphStats}}, summaries regarding the vertex degrees and the most important vertices per community are generated. Note that the communities need to first be calculated using the \code{\link{calculateCommunitiesStats}} function
#'
#' @template GRN
#' @template outputFolder
#' @template plotAsPDF
#' @template pdf_width
#' @template pdf_height
#' @template basenameOutput
#' @inheritParams plotCommunitiesEnrichment
#' @param communities Numeric vector. Default \code{seq_len(10)}. Depending on what was specified in the \code{display} parameter, this parameter would indicate either the rank or the label of the communities to be plotted. i.e. for \code{communities = c(1,4)}, if \code{display = "byRank"} the results for the first and fourth largest communities will be plotted. if \code{display = "byLabel"}, the results for the communities labeled \code{"1"}, and \code{"4"} will be plotted. If set to \code{NULL}, all communities will be plotted
#' @template forceRerun
#' @param topnGenes Integer. Default 20. Number of genes to plot, sorted by their rank or label.
#' @param topnTFs Integer. Default 20. Number of TFs to plot, sorted by their rank or label.
#' @return The same \code{\linkS4class{GRN}} object, without modifications. A single PDF file is produced with the statistics.
#' @seealso \code{\link{plotGeneralGraphStats}}
#' @seealso \code{\link{calculateCommunitiesStats}}
#' @seealso \code{\link{calculateCommunitiesEnrichment}}
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' GRN = plotCommunitiesStats(GRN, outputFolder = ".", plotAsPDF=FALSE)
#' @export
plotCommunitiesStats <- function(GRN, outputFolder = NULL, basenameOutput = NULL, 
                                 display = "byRank", communities = seq_len(10), 
                                 topnGenes = 20, topnTFs = 20, 
                                 plotAsPDF = TRUE, pdf_width = 12, pdf_height = 12, 
                                 forceRerun = FALSE){
  
  start = Sys.time()
  GRN = .addFunctionLogToObject(GRN)
  checkmate::assertClass(GRN, "GRN")
  
  checkmate::assert(checkmate::checkNull(outputFolder), checkmate::checkCharacter(outputFolder, min.chars = 1))
  checkmate::assert(checkmate::checkNull(basenameOutput), checkmate::checkCharacter(basenameOutput, len = 1, min.chars = 1, any.missing = FALSE))
  
  checkmate::assertIntegerish(topnGenes)
  checkmate::assertIntegerish(topnTFs)
  checkmate::assertSubset(display , c("byRank", "byLabel"))
  checkmate::assert(checkmate::checkNull(communities), checkmate::checkNumeric(communities, lower = 1, any.missing = FALSE, min.len = 1))
  checkmate::assertFlag(plotAsPDF)
  checkmate::assertNumeric(pdf_width, lower = 5, upper = 99)
  checkmate::assertNumeric(pdf_height, lower = 5, upper = 99)
  checkmate::assertFlag(forceRerun)
  
  outputFolder = .checkOutputFolder(GRN, outputFolder)
  fileCur = paste0(outputFolder, dplyr::if_else(is.null(basenameOutput), .getOutputFileName("plot_communityStats"), basenameOutput), ".pdf")
  
  if (!file.exists(fileCur) | forceRerun | !plotAsPDF) {
    
    if (plotAsPDF) {
      .checkOutputFile(fileCur)
      grDevices::pdf(fileCur, width = pdf_width, height = pdf_height)
      futile.logger::flog.info(paste0("Plotting to file ", fileCur))
    } else {
      futile.logger::flog.info(paste0("Plotting directly"))
    }
    
    
    vertexMetadata = as.data.frame(igraph::vertex.attributes(GRN@graph$TF_gene$graph))
    
    
    if (is.null(vertexMetadata$community)) {
      
      g = ggplot() + annotate("text", x = 4, y = 25, size=5, color = "red", 
                              label = "Nothing to plot:\nNo communities found") + theme_void()
      plot(g)
      
      message = "No communities found in GRN object. Make sure the filtered connections contain also genes."
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
      
    } else {
      if (!is.null(communities)) {
        if (display == "byRank"){
          # Only display communities we have data for, in a reasonable order
          communitiesDisplay = .selectCommunitesByRank(GRN, communities)
        } else{ # byLabel
          communitiesDisplay = as.character(communities)
        }
      } else {
        communitiesDisplay = unique(vertexMetadata$community)
      }
      
      # TODO here
      # communityVertices = GRN@stats$communityVertices %>%
      #     dplyr::filter(community %in% communitiesDisplay)
      # 
      communityVertices = vertexMetadata %>%
        dplyr::filter(community %in% communitiesDisplay) %>%
        dplyr::mutate(Class = dplyr::if_else(isTF, "TF", "gene"))
      
      # Class: TF or gene
      gCommunityVertices = ggplot(communityVertices, aes(x = community, fill = Class)) +
        geom_bar(position = "stack") +
        ggtitle("Vertices per community") +
        xlab("Community") +
        ylab("Vertex count") +
        scale_fill_manual(values = c("#3B9AB2", "#E1AF00")) + 
        theme_bw()
      
      plot(gCommunityVertices)
      
      for (communityCur in as.character(communitiesDisplay)) { # change this to select communities
        
        communityVerticesCur = communityVertices %>%
          dplyr::filter(community == communityCur) %>%
          dplyr::pull(name)
        
        community_subgraph.df = 
          igraph::induced_subgraph(graph = GRN@graph$TF_gene$graph, 
                                   vids = communityVerticesCur) %>%
          igraph::as_long_data_frame() %>% 
          dplyr::select(from_name, to_name, connectionType, from_names_TF_all, to_names_gene) %>%
          dplyr::rename(V1 = from_name, V2 = to_name, V1_name = from_names_TF_all, V2_name = to_names_gene) %>%
          dplyr::mutate(community = communityCur, isTFTFinteraction = FALSE)
        
        # Handle cases where TF-TF interactions occur and the V2_name is NA because the gene name is NA
        allTFs = GRN@graph$TF_gene$table$V1
        index_V2_isTF = which(community_subgraph.df$V2 %in% allTFs)
        community_subgraph.df$isTFTFinteraction[index_V2_isTF] = TRUE
        match_TF_name = match(community_subgraph.df$V2[index_V2_isTF], vertexMetadata$name)
        # Change only those V2_names for which V2 is a TF, put the TF name in this case
        community_subgraph.df$V2_name[index_V2_isTF] = vertexMetadata$names_TF_all[match_TF_name]
        
        
        community.degreeStats = .getDegreeStats(GRN, community_subgraph.df, 
                                                nCentralGenes = topnGenes,
                                                nCentralTFs = topnTFs)
        gDegreeDist = community.degreeStats$figures$degreeDist + ggtitle("Degree Distribution")
        
        
        gTopGenes   = community.degreeStats$figures$topGenes + ggtitle(paste0("Top ", topnGenes," degree-central genes"))
        gTopTFs     = community.degreeStats$figures$topTFs   + ggtitle(paste0("Top ", topnTFs  ," degree-central TFs"))
        
        # Only needed here temporarily
        # futile.logger::flog.info(paste0("Generating community graph..."))
        
        GRN@graph$communitySubgraph = list() 
        GRN@graph$communitySubgraph$table = community_subgraph.df
        GRN@graph$communitySubgraph$graph = .buildGraph(GRN@graph$communitySubgraph$table, 
                                                        directed = GRN@graph$parameters$directed, 
                                                        allowLoops = GRN@graph$parameters$allowLoops, 
                                                        removeMultiple = GRN@graph$parameters$removeMultiple,
                                                        silent = TRUE
        )
        
        
        community.eigenStats = .getEigenCentralVertices(GRN, graphType = "communitySubgraph", 
                                                        nCentralGenes = topnGenes,
                                                        nCentralTFs = topnTFs
        )
        GRN@graph$communitySubgraph = NULL
        
        gTopEigenGenes = community.eigenStats$topGenes + ggtitle(paste0("Top ", topnGenes," eigenvector-central genes"))
        gTopEigenTFs   = community.eigenStats$topTFs   + ggtitle(paste0("Top ", topnTFs,  " eigenvector-central TFs"))
        
        plot(gDegreeDist + (gTopGenes/gTopTFs) + 
               patchwork::plot_layout(widths = c(2,2)) + 
               patchwork::plot_annotation(title = paste0("Community ", communityCur)) )
        
        plot(gTopEigenGenes/gTopEigenTFs)
      }
      
    }
    
    
    if (plotAsPDF) grDevices::dev.off()
    
  } else {
    futile.logger::flog.info(paste0("File ", fileCur, " already exists and forceRerun has been set to FALSE. Do nothing."))
  }
  
  .printExecutionTime(start)
  
  GRN
}

#' Plot community-based enrichment results
#' 
#' Similarly to \code{\link{plotGeneralEnrichment}}, the results of the community-based enrichment analysis are plotted.. By default, the results for the 10 largest communities are displayed. Additionally, if a general enrichment analysis was previously generated, this function plots an additional heatmap to compare the general enrichment with the community based enrichment. A reduced version of this heatmap is also produced where terms are filtered out to improve visibility and display and highlight the most significant terms.
#' 
#' @inheritParams plotGeneralEnrichment
#' @param display Character. Default \code{"byRank"}. One of: \code{"byRank"}, \code{"byLabel"}. Specify whether the communities will by displayed based on their rank, where the largest community (with most vertices) would have a rank of 1, or by their label. Note that the label is independent of the rank.
#' @param communities \code{NULL} or numeric vector. Default \code{NULL}. If set to \code{NULL}, the default, all communities enrichments that have been calculated before are plotted. If a numeric vector is specified: Depending on what was specified in the \code{display} parameter, this parameter indicates either the rank or the label of the communities to be plotted. i.e. for \code{communities = c(1,4)}, if \code{display = "byRank"} the results for the first and fourth largest communities are plotted. if \code{display = "byLabel"}, the results for the communities labeled \code{"1"}, and \code{"4"} are plotted. 
#' @param nSignificant Numeric. Default 3. Threshold to filter out an ontology term with less than \code{nSignificant} overlapping genes. 
#' @param nID Numeric. Default 10. For the reduced heatmap, number of top terms to select per community.
#' @param maxWidth_nchar_plot Integer (>=10). Default 100. Maximum number of characters for a term before it is truncated.
#' @return  The same \code{\linkS4class{GRN}} object, without modifications. A single PDF file is produced with the results.
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' GRN = plotCommunitiesEnrichment(GRN, outputFolder = ".", plotAsPDF=FALSE)
#' @export
#' @import ggplot2
#' @importFrom grid gpar
plotCommunitiesEnrichment <- function(GRN, outputFolder = NULL, basenameOutput = NULL, 
                                      display = "byRank", communities = NULL,
                                      topn_pvalue = 30, p = 0.05, nSignificant = 2, nID = 10, maxWidth_nchar_plot = 100,
                                      display_pAdj = FALSE,
                                      plotAsPDF = TRUE, pdf_width = 12, pdf_height = 12,
                                      forceRerun = FALSE) {
  
  start = Sys.time()
  GRN = .addFunctionLogToObject(GRN)
  checkmate::assertClass(GRN, "GRN")
  
  checkmate::assertSubset(display , c("byRank", "byLabel"))
  checkmate::assert(checkmate::checkNull(communities), checkmate::checkNumeric(communities, lower = 1, any.missing = FALSE, min.len = 1))
  checkmate::assertNumeric(p, lower = 0, upper = 1)
  checkmate::assertIntegerish(topn_pvalue, lower = 1)
  checkmate::assertNumeric(nSignificant, lower = 1)
  checkmate::assertNumeric(nID, lower = 1)
  checkmate::assertFlag(display_pAdj)
  checkmate::assertFlag(plotAsPDF)
  checkmate::assertIntegerish(maxWidth_nchar_plot, lower = 10, upper = 500)
  checkmate::assertNumeric(pdf_width, lower = 5, upper = 99)
  checkmate::assertNumeric(pdf_height, lower = 5, upper = 99)
  checkmate::assertFlag(forceRerun)
  checkmate::assert(checkmate::checkNull(outputFolder), checkmate::checkCharacter(outputFolder, min.chars = 1))
  checkmate::assert(checkmate::checkNull(basenameOutput), checkmate::checkCharacter(basenameOutput, len = 1, min.chars = 1, any.missing = FALSE))
  
  if (is.null(GRN@stats$Enrichment$byCommunity)) {
    message = "No communities found, cannot calculate enrichment. Run the functioncalculateCommunitiesStats first. If you did already, it looks like no communities could be identified before"
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
  outputFolder = .checkOutputFolder(GRN, outputFolder)
  fileCur = paste0(outputFolder, dplyr::if_else(is.null(basenameOutput), .getOutputFileName("plot_communityEnrichment"), basenameOutput), ".pdf")
  
  futile.logger::flog.info(paste0("Including terms only if overlap is at least ", nSignificant, " genes."))
  
  if (!file.exists(fileCur) | forceRerun | !plotAsPDF) {
    
    if (plotAsPDF) {
      .checkOutputFile(fileCur)
      grDevices::pdf(fileCur, width = pdf_width, height = pdf_height)
      futile.logger::flog.info(paste0("Plotting to file ", fileCur))
    } else {
      futile.logger::flog.info(paste0("Plotting directly"))
    }
    
    
    if (is.null(GRN@stats$Enrichment$general)){
      message = paste0("Could not find general enrichment analysis. Please run the function calculateGeneralEnrichment first.")
      .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    } 
    
    GRN@stats$Enrichment$byCommunity$combined = NULL
    allCalculatedCommunities = names(GRN@stats$Enrichment$byCommunity)
    
    if (display == "byLabel"){
      
      if (is.null(communities)) {
        message = paste("If display = \"byLabel\", the parameter \"communities\" cannot be NULL.")
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
      }
      
      communitiesDisplay = as.character(communities)
      # issue a warning if the community label does not exist
      diff.communities = setdiff(communitiesDisplay, allCalculatedCommunities)
      if (length(diff.communities)>0){
        message = paste("The following communities do not exist and will not be in the analysis: ", paste0(diff.communities, collapse = " + "))
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
        communitiesDisplay = setdiff(communitiesDisplay, diff.communities)
      }
    } else { #byRank
      
      if (is.null(communities)) {
        communitiesDisplay = allCalculatedCommunities
      } else {
        communitiesDisplay = .selectCommunitesByRank(GRN, communities = communities)
      }
      
    }
    
    futile.logger::flog.info(paste0(" Plotting the enrichment for the following communities: ", paste0(communitiesDisplay, collapse = ",")))
    
    # Preparing the heatmap, has to be calculated only once
    vertexMetadata = as.data.frame(igraph::vertex.attributes(GRN@graph$TF_gene$graph))
    # Get the number of vertexes per community as additional annotation column for the heatmap
    geneCounts = vertexMetadata %>%
      dplyr::select(.data$name, .data$community) %>%
      dplyr::distinct() %>%
      dplyr::count(.data$community)
    
    
    allOntologies = .checkEnrichmentCongruence_general_community(GRN, type = "community")
    
    # Reset, as this is created anew and otherwise raises a warning due to the missing "result" slot
    GRN@stats$Enrichment$byCommunity[["combined"]] = list()
    
    for (ontologyCur in allOntologies) {
      
      futile.logger::flog.info(paste0(" Plotting results for ontology ", ontologyCur))
      
      pValPrefix = "raw "
      # p-adjust only available for non-GO ontologies
      if (display_pAdj && !stringr::str_starts("GO_", ontologyCur)) {
        pValPrefix = "adj. "
      }
      
      for (communityCur in as.character(communitiesDisplay)) {
        
        if (is.null(GRN@stats$Enrichment$byCommunity[[communityCur]][[ontologyCur]])){
          # GRN = calculateCommunitiesEnrichment(GRN = GRN, selection = "byLabel", communities = communityCur)
          message = paste0("Could not find community enrichment results for ", ontologyCur, " for community ", communityCur, ". Please run calculateCommunitiesEnrichment first or change the communities to plot")
          .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
        }
        
        
        # dataCur this contains all ontologies, and within each ontology, 2 slots: results, parameters
        dataCur = GRN@stats$Enrichment$byCommunity[[communityCur]] [[ontologyCur]]
        titleCur = paste0("Enrichment Analysis - Community ", communityCur)
        
        if (is.null(dataCur)) {
          message = paste0("Could not find enrichment results for ontology ", ontologyCur, " and community ", communityCur, ". Rerun the function calculateCommunitiesEnrichment.")
          .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
        }
        .plotEnrichmentGeneral(dataCur, ontologyCur, titleCur, 
                               topn_pvalue = topn_pvalue, p = p, display_pAdj = display_pAdj)
        
      } # end for all communities 
      
      # Summary heatmap to compare terms enriched in the general network to their enrichment in the communities:
      # Currently assumes that the general enrichment has been run
      
      # Take one community as example for which ontologies have been produced
      GRN@stats$Enrichment$byCommunity[["combined"]][[ontologyCur]] = 
        .combineEnrichmentResults(GRN, type = "byCommunity", ontologyCur, 
                                  p = p, nSignificant = nSignificant, display_pAdj) %>%
        dplyr::filter(community %in% c("all", as.character(communitiesDisplay)))
      
      # It can happen that some TFs are filtered out here because all terms are not significant. Thus, even though
      # 5 TFs have been requested only 4 actually show any enriched terms
      
      if (nrow(GRN@stats$Enrichment$byCommunity[["combined"]][[ontologyCur]]) == 0) {
        message = paste0("  No enrichment terms passed the filters when creating the across-community summary plot for ontology ", ontologyCur, ". Skipping. You may adjust the parameter nSignificant to a lower value")
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
        next
      }
      
      # Convert to wide format and filter those terms that are significant at least once
      all.df.wide = GRN@stats$Enrichment$byCommunity[["combined"]][[ontologyCur]] %>% 
        dplyr::select(community, ID, pval) %>%
        tidyr::pivot_wider(names_from = community, values_from = pval) %>%
        dplyr::mutate_at(dplyr::vars(!dplyr::contains("ID")), as.numeric) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(nSig = sum(dplyr::c_across(where(is.numeric)) <= p, na.rm = TRUE)) %>% # In how many columns significant?
        dplyr::ungroup() %>%
        dplyr::filter(nSig > 0)
      
      
      # sort by community size
      communities.order = c("all", intersect(.selectCommunitesByRank(GRN , communities = NULL), colnames(all.df.wide)))
      
      
      matrix.m = all.df.wide %>%
        dplyr::mutate(Term = GRN@stats$Enrichment$byCommunity[["combined"]][[ontologyCur]]$Term[match(ID, GRN@stats$Enrichment$byCommunity[["combined"]][[ontologyCur]]$ID)]) %>%
        dplyr::filter(!is.na(Term)) %>%
        dplyr::select(Term, dplyr::any_of(communities.order)) %>% # reorder the table based on the previously generated custom order
        dplyr::mutate_at(dplyr::vars(!dplyr::contains("Term")), function(x){return(-log10(x))}) %>%
        dplyr::mutate(Term = stringr::str_trunc(as.character(Term), width = maxWidth_nchar_plot, side = "right")) %>%
        tibble::column_to_rownames("Term") %>%
        as.matrix()
      
 
      geneCounts_communities = geneCounts %>%
        dplyr::filter(community %in% as.character(colnames(matrix.m)),
                      community %in% geneCounts$community) %>%
        dplyr::arrange(dplyr::desc(.data$n))
      
      # Sanity check
      stopifnot(identical(as.character(geneCounts_communities$community), colnames(matrix.m)[-1]))
        
      # Common heatmap parameters for both p1 and p2
      top_annotation = ComplexHeatmap::HeatmapAnnotation(
        nGenes = ComplexHeatmap::anno_barplot(
          x = c(sum(geneCounts$n), geneCounts_communities$n), 
          border = FALSE,  bar_width = 0.8,  gp = grid::gpar(fill = "#046C9A")),
        annotation_name_gp = grid::gpar(fontsize=9), annotation_name_side = "left", annotation_name_rot = 90)
      
      axixStr = paste0("-log10\n(", pValPrefix, "p-value)")
      colors_pvalue = viridis::viridis(n=30, direction = -1)
      
      futile.logger::flog.info(paste0("  Including ", nrow(matrix.m), " terms in the full summary heatmap and " , ncol(matrix.m), " columns"))
      
      p1 = .drawCombinedHeatmap(matrix = matrix.m, pdf_width = pdf_width, pdf_height = pdf_height,
                                name = axixStr, 
                                colors_pvalue = colors_pvalue, 
                                column_title = paste0("Summary of all significantly enriched terms\nacross all communities and overall (Ontology: ", ontologyCur,")"),
                                top_annotation = top_annotation)
      
      print(p1)
      
      
      # Now focus on the top X only per community
      ID_subset =  GRN@stats$Enrichment$byCommunity[["combined"]][[ontologyCur]] %>% 
        dplyr::group_by(community) %>% 
        dplyr::arrange(.data$pval) %>% 
        dplyr::slice(seq_len(nID)) %>%
        dplyr::pull(ID) %>% as.character()
      
      
      matrix.m = all.df.wide %>%
        dplyr::filter(ID %in% ID_subset) %>%
        dplyr::mutate(Term = GRN@stats$Enrichment$byCommunity[["combined"]][[ontologyCur]]$Term[match(ID, GRN@stats$Enrichment$byCommunity[["combined"]][[ontologyCur]]$ID)]) %>%
        dplyr::filter(!is.na(Term)) %>%
        dplyr::select(-nSig, -ID) %>%
        tibble::column_to_rownames("Term") %>%
        dplyr::mutate_at(dplyr::vars(!dplyr::contains("ID")), function(x){return(-log10(x))}) %>%
        dplyr::select(dplyr::any_of(communities.order)) %>% # reorder based on the previously generated custom order
        as.matrix()
      
      futile.logger::flog.info(paste0("  Including ", nrow(matrix.m), " terms in the reduced summary heatmap and " , ncol(matrix.m), " columns"))
      
      p2 = .drawCombinedHeatmap(matrix = matrix.m, pdf_width = pdf_width, pdf_height = pdf_height, 
                                name = axixStr, 
                                colors_pvalue = colors_pvalue, 
                                column_title = paste0("Summary of top ", nID, " enriched terms\nper community and overall (Ontology: ", ontologyCur,")"),
                                top_annotation = top_annotation)
      
      print(p2)
      
    }
    
    
    
    if (plotAsPDF) grDevices::dev.off()
    
  } else {
    futile.logger::flog.info(paste0("File ", fileCur, " already exists, not overwriting due to forceRerun = FALSE"))
  }
  
  .printExecutionTime(start)
  
  GRN
  
}


.drawCombinedHeatmap <- function(matrix, pdf_width, pdf_height, name, colors_pvalue, column_title, top_annotation) {
  
  nTerms = nrow(matrix)
  
  # TODO: very crude so far, and not dependent on the pdf_height
  fontsize_row_names_gp = dplyr::case_when(
    nTerms < 20  ~ 12,
    nTerms < 40  ~ 10,
    nTerms < 50  ~ 8,
    nTerms < 80  ~ 6,
    nTerms < 110  ~ 5,
    nTerms < 140  ~ 4,
    nTerms < 200  ~ 2,
    TRUE ~ 1
  )
  
  # TODO:
  fontsize_column_names_gp = dplyr::case_when(
    TRUE ~ 8
  )
  
  # Prevent :Error: You should have at least two distinct break values.
  if (unique(matrix) %>% as.vector() %>% length() == 1) {
    colors_pvalue = colors_pvalue[1]
  }
  
  ComplexHeatmap::Heatmap(
    matrix,
    # Do NOT set the width here, this messes it up somehow, seems to be much better when omitting it
    #  heatmap_width  = unit(pdf_width, "cm"),
    #  heatmap_height = unit(pdf_height, "cm"),
    name = name,
    col = colors_pvalue,
    cluster_columns = FALSE, cluster_rows = FALSE,
    row_names_side = "left", row_names_gp = grid::gpar(fontsize = fontsize_row_names_gp), 
    column_title = column_title,
    column_names_gp = grid::gpar(fontsize = fontsize_column_names_gp),
    # heatmap annotation: bar plot displaying the number of genes in each subgroup
    top_annotation = top_annotation,
    row_names_max_width = ComplexHeatmap::max_text_width(
      rownames(matrix), 
      gp = grid::gpar(fontsize = fontsize_row_names_gp)
    )
  )
  
  
  
}



#' Plot TF-based GO enrichment results
#' 
#' This function plots the enrichment results. The result consist of a dot plot per specified TF, as well as two comparative heatmaps. The first heatmap displays the p value for each GO term across the TFs. Terms that The second heatmap is a subset of the first, where select terms are kept or filtered out for better visibility and display.
#' 
#' @inheritParams plotGeneralEnrichment
#' @inheritParams plotCommunitiesEnrichment
#' @param TF.names \code{NULL} or character vector. Default \code{NULL}. For \code{rankType="custom"} the names of the TFs to plot. Ignored otherwise.
#' @param rankType Character. One of: "degree", "EV", "custom". This parameter will determine the criterion to be used to identify the "top" nodes. If set to "degree", the function will select top nodes based on the number of connections they have, i.e. based on their degree-centrality. If set to "EV" it will select the top nodes based on their eigenvector-centrality score in the network.
#' @param n NULL or numeric. Default NULL. If set to NULL, all previously calculated TF enrichments will be plotted. If set to a value between (0,1), it is treated as a percentage of top nodes. If the value is passed as an integer it will be treated as the number of top nodes. This parameter is not relevant if rankType = "custom".
#' @return The same \code{\linkS4class{GRN}} object, without modifications. A single PDF file is produced with the results.
#' @seealso \code{\link{calculateTFEnrichment}}
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' GRN = plotTFEnrichment(GRN, n = 5, outputFolder = ".", plotAsPDF=FALSE)
#' @export
#' @importFrom grid gpar
plotTFEnrichment <- function(GRN, rankType = "degree", n = NULL, TF.names = NULL,
                             topn_pvalue = 30, p = 0.05, 
                             nSignificant = 2, nID = 10,
                             display_pAdj = FALSE,
                             maxWidth_nchar_plot = 100,
                             outputFolder = NULL, 
                             basenameOutput = NULL, 
                             plotAsPDF = TRUE, pdf_width = 12, pdf_height = 12,
                             forceRerun = FALSE) {
  
  start = Sys.time()
  GRN = .addFunctionLogToObject(GRN)
  
  checkmate::assertClass(GRN, "GRN")
  checkmate::assertSubset(rankType, c("degree", "EV", "custom"))
  checkmate::assert(checkmate::checkNull(n), checkmate::checkNumeric(n))
  checkmate::assertNumeric(p, lower = 0, upper = 1)
  checkmate::assertNumeric(topn_pvalue, lower = 1)
  checkmate::assertNumeric(nSignificant, lower = 0)
  checkmate::assertNumeric(nID, lower = 0)
  checkmate::assertFlag(plotAsPDF)
  checkmate::assertNumeric(pdf_width, lower = 5, upper = 99)
  checkmate::assertNumeric(pdf_height, lower = 5, upper = 99)
  checkmate::assertFlag(display_pAdj)
  checkmate::assert(checkmate::checkNull(outputFolder), checkmate::checkCharacter(outputFolder, min.chars = 1))
  checkmate::assert(checkmate::checkNull(basenameOutput), checkmate::checkCharacter(basenameOutput, len = 1, min.chars = 1, any.missing = FALSE))
  checkmate::assertFlag(forceRerun)
  
  
  
  outputFolder = .checkOutputFolder(GRN, outputFolder)
  
  
  if (rankType == "custom"){
    if(is.null(TF.names)){
      futile.logger::flog.error("To plot the GO enrichment for a custom set of TFs, you must provide the TF names in the 'TF.names' parameter.")
    }
    wrongTFs = setdiff(TF.names, unique(GRN@connections$all.filtered$`0`$TF.name))
    if (length(wrongTFs)>0){
      futile.logger::flog.warn(paste0("The TFs ",  paste0(wrongTFs, collapse = " + "), " are not in the filtered GRN. They will be ommited from the results."))
    }
    TFset = setdiff(TF.names, wrongTFs) 
  } else{ #rankType = "degree"
    
    if (is.null(n)) {
      TFset = setdiff(names(GRN@stats$Enrichment$byTF), "combined")
    } else {
      TFset = getTopNodes(GRN, nodeType = "TF", rankType = rankType, n = n) %>% dplyr::pull(TF.name)
    }
    
  }
  
  fileCur = paste0(outputFolder, 
                   dplyr::if_else(is.null(basenameOutput), .getOutputFileName("plot_TFEnrichment"), basenameOutput), 
                   ".pdf")
  
  if (!file.exists(fileCur) | forceRerun | !plotAsPDF) {
    
    if (plotAsPDF) {
      .checkOutputFile(fileCur)
      grDevices::pdf(fileCur, width = pdf_width, height = pdf_height)
      futile.logger::flog.info(paste0("Plotting to file ", fileCur))
    } else {
      futile.logger::flog.info(paste0("Plotting directly"))
    }
    
    
    # Reset, as this is created anew and otherwise raises a warning due to the missing "result" slot
    GRN@stats$Enrichment$byCommunity[["combined"]] = NULL
    
    
    # Get the number of vertexes per TF as additional annotation column for the heatmap
    vertexMetadata = as.data.frame(igraph::vertex.attributes(GRN@graph$TF_gene$graph))
    nodeDegree = igraph::degree(GRN@graph$TF_gene$graph)
    
    # nodeDegree_TFset = sort(nodeDegree[GRN@data$TFs$translationTable %>% dplyr::filter(TF.name %in% as.character(TFset)) %>% 
    #                                      dplyr::pull(TF.ENSEMBL)], decreasing = TRUE)
    
    nodeDegree_TFset = dplyr::left_join(GRN@data$TFs$translationTable, 
                                        as.data.frame(nodeDegree) %>% tibble::rownames_to_column("ENSEMBL"), by = "ENSEMBL") %>%
                        dplyr::filter(TF.name %in% as.character(TFset))
    
    
    allOntologies = .checkEnrichmentCongruence_general_community(GRN, type = "TF")
    for (ontologyCur in allOntologies) {
      
      futile.logger::flog.info(paste0(" Plotting results for ontology ", ontologyCur))
      
      pValPrefix = "raw "
      # p-adjust only available for non-GO ontologies
      if (display_pAdj && !stringr::str_starts("GO_", ontologyCur)) {
        pValPrefix = "adj. "
      }
      
      for (TFCur in as.character(TFset)){
        
        TF.ENSEMBL = GRN@data$TFs$translationTable %>% dplyr::filter(TF.name == TFCur) %>% dplyr::pull(TF.ENSEMBL)
        
        TF.name.full = paste0(TFCur, " (", TF.ENSEMBL, ")")
        
        futile.logger::flog.info(paste0("  TF ", TF.name.full))
        
        # if the enrichment slot for the TFs is empty, calculate the enrichment
        if (is.null(GRN@stats$Enrichment[["byTF"]][[TFCur]])){ 
          message = paste0("Could not find TF enrichment results. Run the function calculateTFEnrichment first or make sure the parameter n that was used for calculateTFEnrichmentn was larger or equal with the current value for n for this function.")
          .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
        }
        
        dataCur = GRN@stats$Enrichment[["byTF"]][[TFCur]][[ontologyCur]]
        if (is.null(dataCur)) {
          message = paste0("Could not find enrichment results for ontology ", ontologyCur, " and TF ", TFCur, ". Rerun the function calculateTFEnrichment.")
          .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
        }
        titleCur = paste0("Enrichment Analysis - TF: ", TF.name.full)
        .plotEnrichmentGeneral(dataCur, ontologyCur, titleCur, topn_pvalue = topn_pvalue,  p = p, display_pAdj = display_pAdj)
        
      }
      
      # Summary heatmap to compare terms enriched in the general network to their enrichment in the communities:
      # Currently assumes that the general enrichment has been run
      
      # Take one community as example for which ontologies have been produced
      GRN@stats$Enrichment$byTF[["combined"]][[ontologyCur]] = 
        .combineEnrichmentResults(GRN, type = "byTF", ontologyCur, 
                                  p = p, nSignificant = nSignificant, display_pAdj) %>%
        dplyr::filter(TF.name %in% c("all", as.character(TFset)))
      
      if (nrow(GRN@stats$Enrichment$byTF[["combined"]][[ontologyCur]]) == 0) {
        message = paste0("  No enrichment terms passed the filters when creating the across-community summary plot for ontology ", ontologyCur, ". Skipping. You may adjust the parameter nSignificant to a lower value")
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
        next
      }
      
      # Convert to wide format and filter those terms that are significant at least once
      all.df.wide = GRN@stats$Enrichment$byTF[["combined"]][[ontologyCur]] %>% 
        dplyr::select(TF.name, ID, pval) %>%
        tidyr::pivot_wider(names_from = TF.name, values_from = pval) %>%
        dplyr::mutate_at(dplyr::vars(!dplyr::contains("ID")), as.numeric) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(nSig = sum(dplyr::c_across(where(is.numeric)) <= p, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::filter(nSig > 0)
      
      
      TF.order = c("all", setdiff(unique(GRN@stats$Enrichment$byTF[["combined"]][[ontologyCur]]$TF.name), "all"))
      
      
      matrix.m = all.df.wide %>%
        dplyr::mutate(Term = GRN@stats$Enrichment$byTF[["combined"]][[ontologyCur]]$Term[match(ID, GRN@stats$Enrichment$byTF[["combined"]][[ontologyCur]]$ID)]) %>%
        dplyr::filter(!is.na(Term)) %>%
        dplyr::select(Term, dplyr::all_of(TF.order)) %>% # reorder the table based on the previously generated custom order
        dplyr::mutate_at(dplyr::vars(!dplyr::contains("Term")), function(x){return(-log10(x))}) %>%
        dplyr::mutate(Term = stringr::str_trunc(as.character(Term), width = maxWidth_nchar_plot, side = "right")) %>%
        tibble::column_to_rownames("Term") %>%
        as.matrix()
      
      
      # Make sure the top annotation has the same dimensionality as the resulting matrix
      nodeDegree_TFset_numbers =  nodeDegree_TFset %>%
        dplyr::filter(TF.name %in% colnames(matrix.m)) %>%
        dplyr::arrange(dplyr::desc(nodeDegree)) %>%
        dplyr::pull(nodeDegree)
      
      top_annotation = ComplexHeatmap::HeatmapAnnotation(
        nodeDegree = ComplexHeatmap::anno_barplot(
          x = c("all" = nrow(vertexMetadata),  nodeDegree_TFset_numbers ), 
          border = FALSE,  bar_width = 0.8,  gp = grid::gpar(fill = "#046C9A")),
        annotation_name_gp = grid::gpar(fontsize=5), annotation_name_side = "left", annotation_name_rot = 90)
      
      axixStr = paste0("-log10\n(", pValPrefix, "p-value)")
      
      colors_pvalue = viridis::viridis(n=30, direction = -1)
      # Why colors here like this
      #colors_pvalue = c("#3B9AB2", "#9EBE91", "#E4B80E", "#F21A00")
      futile.logger::flog.info(paste0("  Including ", nrow(matrix.m), " terms in the full heatmap and " , ncol(matrix.m), " columns"))
      
      p1 = .drawCombinedHeatmap(matrix = matrix.m, pdf_width = pdf_width, pdf_height = pdf_height,
                                name = axixStr, 
                                colors_pvalue = colors_pvalue, 
                                column_title = paste0("Summary of all significantly enriched terms\nacross all TFs and overall (Ontology: ", ontologyCur,")"),
                                top_annotation = top_annotation)
      
      print(p1)
      
      
      # Now focus on the top X only per community
      ID_subset =  GRN@stats$Enrichment$byTF[["combined"]][[ontologyCur]] %>% 
        dplyr::group_by(TF.name) %>% 
        dplyr::arrange(pval) %>% 
        dplyr::slice(seq_len(nID)) %>%
        dplyr::pull(ID) %>% as.character()
      
      
      matrix.m = all.df.wide %>%
        dplyr::filter(ID %in% ID_subset) %>%
        dplyr::mutate(Term = GRN@stats$Enrichment$byTF[["combined"]][[ontologyCur]]$Term[match(ID, GRN@stats$Enrichment$byTF[["combined"]][[ontologyCur]]$ID)]) %>%
        dplyr::filter(!is.na(Term)) %>%
        dplyr::select(-nSig, -ID) %>%
        tibble::column_to_rownames("Term") %>%
        dplyr::mutate_at(dplyr::vars(!dplyr::contains("ID")), function(x){return(-log10(x))}) %>%
        dplyr::select(dplyr::all_of(TF.order)) %>% # reorder based on the previously generated custom order
        as.matrix()
      
      futile.logger::flog.info(paste0("  Including ", nrow(matrix.m), " terms in the reduced summary heatmap and " , ncol(matrix.m), " columns"))
      
      p2 = .drawCombinedHeatmap(matrix = matrix.m, pdf_width = pdf_width, pdf_height = pdf_height, 
                                name = axixStr, 
                                colors_pvalue = colors_pvalue, 
                                column_title = paste0("Summary of top ", nID, " enriched terms\n per TF and overall (Ontology: ", ontologyCur,")"),
                                top_annotation = top_annotation)
      
      print(p2)
      
      
    }
    #  old code, restructured
    # # plot the comparative heatmap:
    # for (ontologyCur in allOntologies){
    #   
    #   enrichmentData = purrr::transpose(purrr::transpose(GRN@stats$Enrichment[["byTF"]])$results)[[ontologyCur]][TFset]
    #   
    #   enrichmentData.all <- enrichmentData %>%
    #     dplyr::bind_rows(.id = "TF") %>%
    #     dplyr::filter(!is.na(Term))
    #     
    #   if (!stringr::str_starts(ontologyCur, "GO_") && display_pAdj){
    #     enrichmentData.all$pval = enrichmentData.all$p.adjust
    #   }
    #    
    #   enrichmentData.all = enrichmentData.all %>% 
    #     dplyr::filter(Significant >= nSignificant) %>% # filter out terms with less than a certain number of significant terms
    #     dplyr::mutate(pval = as.numeric(gsub(">|<", "", pval))) %>%
    #     dplyr::filter(pval <= p ) # filter out terms below a certain p value threshold
    #   
    #   # convert to wide format and only keep GO terms/ rows that are significant in at least 1 TF
    #   dataCur = enrichmentData.all %>% 
    #     dplyr::select(TF, ID,  pval) %>%
    #     tidyr::pivot_wider(names_from = TF, values_from = pval) %>%
    #     dplyr::rowwise() %>%
    #     dplyr::mutate(nSig = sum(dplyr::c_across(tidyselect::where(is.numeric)) <= p, na.rm = TRUE)) %>%
    #     dplyr::ungroup() %>%
    #     dplyr::filter(nSig >0) 
    #   
    #   # tweaks to the matrix passed to the heatmap
    #   matrix.m = dataCur %>%
    #     dplyr::select(-nSig, -ID) %>%
    #     dplyr::mutate_at(dplyr::vars(!dplyr::contains("GO")), function(x){return(-log10(x))}) %>% # modify pval to -log10 pval
    #     dplyr::select(sort(colnames(.))) %>% # to alphabetically sort the TFs; similar TF names will be next to each other
    #     dplyr::arrange(dplyr::desc(dplyr::across())) %>% # to order the GO terms, the terms will be grouped 
    #     as.matrix()
    #   
    #   # dplyr::mutate(Term = stringr::str_trunc(as.character(Term), width = maxWidth_nchar_plot, side = "right")) %>%
    #   #     tibble::column_to_rownames("Term") %>%
    #   
    #   p1 = ComplexHeatmap::Heatmap(matrix.m,
    #                                cluster_rows = FALSE, cluster_columns = FALSE,
    #                                name = paste0("-log10\n(", pValPrefix, "p-value)"),
    #                                column_title = paste0("TF Enrichment Analysis (", ontologyCur, ")"),
    #                                row_names_side = "left", row_names_gp = grid::gpar(fontsize = 8),
    #                                #col = (wes_palette(name = "Zissou1", n=4 ,type = "continuous")),
    #                                col = c("#3B9AB2", "#9EBE91", "#E4B80E", "#F21A00")
    #   )
    #   print(p1)
    #   
    #   ##### reduced heatmap ####
    #   
    #   # from the previously generated matrix/heatmap, select the top n GO terms for each TF to have a reduced and more legible heatmap.
    #   GO_subset = as.vector((enrichmentData.all %>% dplyr::group_by(TF) %>% dplyr::arrange(pval) %>% dplyr::group_modify(~ head(.x, n= nID)))$ID)
    #   
    #   matrix.m = dataCur %>%
    #     dplyr::filter(ID %in% GO_subset) %>%
    #     dplyr::mutate(GO.term = enrichmentData.all$Term[match(ID, enrichmentData.all$ID)]) %>%
    #     tibble::column_to_rownames("GO.term") %>%
    #     dplyr::select(-nSig, -ID) %>%
    #     dplyr::mutate_at(dplyr::vars(!dplyr::contains("GO")), function(x){return(-log10(x))}) %>%
    #     dplyr::select(sort(colnames(.))) %>%
    #     dplyr::arrange(dplyr::desc(dplyr::across())) %>%
    #     as.matrix()
    #   
    #   p2 = ComplexHeatmap::Heatmap(matrix.m,
    #                                name = paste0("-log10\n(", pValPrefix, "p-value)"),
    #                                heatmap_width = unit(10,"cm"),
    #                                cluster_rows = FALSE, cluster_columns = FALSE,
    #                                column_title = paste0("TF Enrichment Analysis (", ontologyCur, ")"),
    #                                row_names_side = "left", row_names_gp = grid::gpar(fontsize = 5),
    #                                #col = wes_palette(name = "Zissou1", n=4 ,type = "continuous")
    #                                col = c("#3B9AB2", "#9EBE91", "#E4B80E", "#F21A00")
    #   )
    #   print(p2)
    # }
    
    if (plotAsPDF) grDevices::dev.off()
    
  } else {
    futile.logger::flog.info(paste0("File ", fileCur, " already exists, not overwriting due to forceRerun = FALSE"))
  }
  
  .printExecutionTime(start)
  
  GRN
  
}

.getDegreeStats <- function(GRN, df, nCentralGenes = 20, nCentralTFs = 20){
  
  # This function takes as input the GRN and the graph in dataframe format 
  # returns histogram of degree distribution and plots of top degree-central TFs and genes
  
  TF.degrees = df %>%
    dplyr::filter(stringr::str_detect(connectionType, "^tf")) %>% 
    dplyr::mutate(name_plot = paste0(V1_name, "\n(", V1, ")")) %>%
    dplyr::count(name_plot, V1) %>% 
    dplyr::rename(ID = V1, ID_all = name_plot, Degree = .data$n) 
  
  # TODO modify
  peak_TFend.degrees = df %>%
    dplyr::filter(stringr::str_detect(connectionType, "peak$")) %>%
    dplyr::count(V2) %>% dplyr::rename(ID = V2, Degree = .data$n)
  
  peak_geneend.degrees = df %>%
    dplyr::filter(stringr::str_detect(connectionType, "^peak")) %>%
    dplyr::count(V1) %>% dplyr::rename(ID = V1, Degree = .data$n)
  
  gene.degrees = df %>%
    dplyr::filter(stringr::str_detect(connectionType, "gene$")) %>%
    dplyr::mutate(name_plot = paste0(V2_name, "\n(", V2, ")")) %>%  
    dplyr::count(name_plot, V2) %>% 
    dplyr::rename(ID = V2, ID_all = name_plot, Degree = .data$n)
  
  degrees.table = dplyr::bind_rows( TF = TF.degrees, 
                                    `peak (TF end)` = peak_TFend.degrees, 
                                    `peak (gene end)` = peak_geneend.degrees, 
                                    gene = gene.degrees,
                                    .id = "class") %>%
    dplyr::mutate(class = droplevels(factor(class, levels = c("TF", "peak (TF end)", "peak (gene end)", "gene")))) %>% # in case it was the TF-gene df that was passed, the peak levels will be dropped
    dplyr::arrange(class, dplyr::desc(Degree))
  
  colours = c("TF" = "#E1AF00", "peak (TF end)" = "#F21A00", "peak (gene end)" = "#F21A00", "gene" = "#3B9AB2")
  
  ## General Degree Distribution ##
  gDegrees = ggplot(degrees.table, aes(x=Degree, fill = class)) +
    geom_histogram(bins = 50) +
    scale_fill_manual(values= colours) +
    theme(legend.position = "none",  axis.text.x = element_text(angle = 45, vjust = 1, hjust =1), 
          strip.text = element_text(size=8)) +
    ylab("Abundance") +
    facet_wrap(~class, labeller = labeller(class=.facetLabel(degrees.table)), scales = "free", ncol = 2) +
    theme_bw()
  
  ## Top degree-central TFs and Genes ##
  top.n.gene = degrees.table %>% 
    dplyr::filter(class =="gene") %>%
    dplyr::arrange(dplyr::desc(Degree)) %>% 
    dplyr::slice(seq_len(nCentralGenes)) 
  
  top.n.tf = degrees.table %>% 
    dplyr::filter(class =="TF") %>%
    dplyr::rename(TF.name = ID) %>%
    dplyr::arrange(dplyr::desc(Degree)) %>% 
    dplyr::slice(seq_len(nCentralTFs))
  
  theme_all = theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =1, size = 8))
  
  gTopGenes = ggplot(data=top.n.gene, aes(x=stats::reorder(ID_all, -Degree), y = Degree)) +
    geom_bar(stat = "identity", fill = "#3B9AB2") +
    geom_text(aes(label = Degree), size = 3) + 
    xlab("Gene name") +
    theme_all
  
  gTopTFs = ggplot(data=top.n.tf, aes(x=stats::reorder(ID_all, -Degree), y = Degree)) +
    geom_bar(stat = "identity", fill = "#E1AF00") +
    geom_text(aes(label = Degree), size = 3) + 
    xlab("TF name") +
    theme_all
  
  results = list()
  results[["tbl"]][["degrees"]]        = degrees.table
  results[["figures"]][["degreeDist"]] = gDegrees
  results[["figures"]][["topGenes"]]   = gTopGenes
  results[["figures"]][["topTFs"]]     = gTopTFs
  
  return(results)
}

.getEigenCentralVertices <- function(GRN, graphType, nCentralGenes = 20, nCentralTFs = 20){
  
  # returns a list of figures for the top eigenvector central Genes and TFs
  if (graphType == "TF_peak_gene") {
    geneSet = GRN@graph$TF_peak_gene$table %>%
      dplyr::filter(connectionType == "peak-gene") %>%
      dplyr::pull(V2) %>%
      unique()
    TFSet   = GRN@graph$TF_peak_gene$table %>%
      dplyr::filter(connectionType == "tf-peak") %>%
      dplyr::pull(V1) %>%
      unique()
    
  } else if (graphType == "TF_gene") {
    
    geneSet = unique(GRN@graph$TF_peak_gene$table$V2)
    TFSet   = unique(GRN@graph$TF_peak_gene$table$V1)
    
  } else if (graphType == "communitySubgraph") {
    # Also a TF-gene graph
    geneSet = unique(GRN@graph$communitySubgraph$table$V2)
    TFSet   = unique(GRN@graph$communitySubgraph$table$V1)
    
  } else {
    message = "Unknown graphType in.getEigenCentralVertices()"
    .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
  # get the centrality score of each vertex and extract top n genes and TFs
  eigen_ctr = igraph::eigen_centrality(GRN@graph[[graphType]]$graph)$vector
  
  
  nCentralGenes = seq_len(min(nCentralGenes, length(eigen_ctr[geneSet])))
  nCentralTFs   = seq_len(min(nCentralTFs  , length(eigen_ctr[TFSet])))
  
  top.n.genes = tibble::tibble(
    gene.ENSEMBL = names(sort(eigen_ctr[geneSet], decreasing = TRUE)[nCentralGenes]),
    Score = sort(eigen_ctr[geneSet], decreasing = TRUE)[nCentralGenes]) %>% 
    stats::na.omit() %>% # Remove NA rows in case there are not enough rows without NA
    dplyr::distinct() %>% 
    dplyr::left_join(GRN@graph[[graphType]]$table %>% dplyr::select(V2, V2_name) %>% dplyr::distinct(), 
                     by = c("gene.ENSEMBL" = "V2") ) %>%
    #dplyr::mutate(gene.name = GRN@connections$all.filtered$`0`$gene.name[match(gene.ENSEMBL, GRN@connections$all.filtered$`0`$gene.ENSEMBL)])
    dplyr::rename(gene.name = V2_name) %>%
    dplyr::mutate(name_plot = paste0(gene.name, "\n(", gene.ENSEMBL, ")")) %>%
    dplyr::filter(Score > 0) # Require the score to be > 0 because we dont want top nodes to have a score of 0
  
  top.n.tf =  tibble::tibble(TF.ENSEMBL = names(sort(eigen_ctr[TFSet], decreasing = TRUE)[nCentralTFs]),
                             Score   =    sort(eigen_ctr[TFSet], decreasing = TRUE)[nCentralTFs]) %>% 
    stats::na.omit() %>% # Remove NA rows in case there are not enough rows without NA
    dplyr::distinct() %>%
    dplyr::left_join(GRN@graph[[graphType]]$table %>% dplyr::select(V1, V1_name) %>% dplyr::distinct(), 
                     by = c("TF.ENSEMBL" = "V1") ) %>%
    dplyr::rename(TF.name = V1_name) %>%
    dplyr::mutate(name_plot = paste0(TF.name, "\n(", TF.ENSEMBL, ")")) %>%
    dplyr::filter(Score > 0) # Require the score to be > 0 because we dont want top nodes to have a score of 0
  
  # dplyr::mutate(TF.name = GRN@connections$all.filtered$`0`$gene.name[match(gene.ENSEMBL, GRN@connections$all.filtered$`0`$gene.ENSEMBL)])
  
  theme_all = theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =1, size = 8))
  
  
  gTopEigenGenes = ggplot(data=top.n.genes, aes(x=stats::reorder(name_plot, -Score), y = Score)) +
    geom_bar(stat = "identity", fill = "#3B9AB2") +
    xlab("Gene name") +
    theme_all
  
  gTopEigenTFs = ggplot(data=top.n.tf, aes(x=stats::reorder(name_plot, -Score), y = Score)) +
    geom_bar(stat = "identity", fill = "#E1AF00") +
    xlab("TF name") +
    theme_all
  
  return (list(topGenes = gTopEigenGenes, topTFs = gTopEigenTFs))
  
}

.facetLabel = function(degrees.table){
  
  summary = degrees.table %>%
    dplyr::group_by(class) %>%
    dplyr::summarise(mean = mean(Degree), median = median(Degree), max = max(Degree), sd = sd(Degree)) %>%
    dplyr::mutate(label = paste0(class, "\n(mean: ", round(mean, 1),
                                 ", median: ", round(median, 1), 
                                 ", max: ", round(max, 1),
                                 ", sd:", round(sd, 1), ")"))
  
  labels <- summary$label
  names(labels) <- summary$class
  
  return (labels)
}


############### GRN visualization ###############

#' Visualize a filtered GRN (in development)
#'
#' @template GRN 
#' @template permuted
#' @param file TODO
#' @param title TODO
#' @param maxRowsToPlot TODO
#' @param useDefaultMetadata TODO
#' @param vertice_color_TFs TODO
#' @param vertice_color_peaks TODO
#' @param vertice_color_genes TODO
#' @template plotAsPDF
#' @template pdf_width
#' @template pdf_height
#' @template forceRerun

#' @return TODO
# #' @export
visualizeGRN <- function(GRN, permuted, file, title = NULL, maxRowsToPlot = 500, useDefaultMetadata = TRUE,
                         vertice_color_TFs = NULL, vertice_color_peaks = NULL, vertice_color_genes = NULL,
                         plotAsPDF = TRUE, pdf_width = 12, pdf_height = 12,
                         forceRerun = FALSE
) {
  
  start = Sys.time()
  GRN = .addFunctionLogToObject(GRN)
  
  checkmate::assertFlag(plotAsPDF)
  checkmate::assertNumeric(pdf_width, lower = 5, upper = 99)
  checkmate::assertNumeric(pdf_height, lower = 5, upper = 99)
  
  futile.logger::flog.info(paste0("Plotting GRN network", dplyr::if_else(is.null(file), "", paste0(" to file ", file))))
  
  if (useDefaultMetadata) {
    metadata_visualization.l = getBasic_metadata_visualization(GRN)
    vertice_color_TFs   = list(metadata_visualization.l[["RNA_expression_TF"]],    "HOCOID",     "baseMean_log")
    vertice_color_genes = list(metadata_visualization.l[["RNA_expression_genes"]], "ENSEMBL_ID", "baseMean_log")
    vertice_color_peaks = list(metadata_visualization.l[["peaks_accessibility"]],   "peakID",     "mean_log")
  }
  
  
  
  grn.merged = getGRNConnections(GRN, permuted = permuted, type = "all.filtered")
  nRows = nrow(grn.merged)
  
  futile.logger::flog.info(paste0("Number of rows: ",nRows))
  if (maxRowsToPlot > 500 & nRows > 500) {
    futile.logger::flog.info(paste0("Plotting many connections takes a lot of time and memory"))
  }
  
  
  
  if (plotAsPDF) {
    .checkOutputFile(file)
    grDevices::pdf(file = file, width = pdf_width, height = pdf_height)
    futile.logger::flog.info(paste0("Plotting to file ", fileCur))
  } else {
    futile.logger::flog.info(paste0("Plotting directly"))
  }
  
  if (nRows > maxRowsToPlot) {
    futile.logger::flog.info(paste0("Number of rows to plot (", nRows, ") exceeds limit of the maxRowsToPlot parameter. Plotting only empty page"))
    plot(c(0, 1), c(0, 1), ann = FALSE, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n', main = title)
    message = paste0(title, "\n\nPlotting omitted.\n\nThe number of rows in the GRN (", nRows, ") exceeds the maximum of ", maxRowsToPlot, ".\nSee the maxRowsToPlot parameter to increase the limit")
    text(x = 0.5, y = 0.5, message, cex = 1.6, col = "red")
    
    if (!is.null(file)) {
      grDevices::dev.off()
    }
    
    .printExecutionTime(start)
    return()
  }
  
  if (nrow(grn.merged) == 0) {
    
    futile.logger::flog.warn(paste0("No rows left in the GRN. Creating empty plot."))
    
    plot(c(0, 1), c(0, 1), ann = FALSE, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n', main = title)
    message = paste0(title, "\n\nThe GRN has no edges that pass the filter criteria.")
    text(x = 0.5, y = 0.5, message, cex = 1.6, col = "red")
    
  } else {
    
    # Fix with length
    colors_categories.l = list(
      "TF"   = RColorBrewer::brewer.pal(3,"Set1")[1], 
      "PEAK" = RColorBrewer::brewer.pal(3,"Set1")[2], 
      "GENE" = RColorBrewer::brewer.pal(3,"Set1")[3])
    
    symbols_categories.l = list(
      "TF"   = 15, # square
      "PEAK" = 21, # circle
      "GENE" = 21 # circle
    )
    
    
    nBins_orig = 100
    nBins_discard = 25
    nBins_real = nBins_orig - nBins_discard
    
    if (!is.null(vertice_color_TFs)) {
      color_gradient = rev(colorspace::sequential_hcl(nBins_orig, h = 10, c = 85, l = c(25, 95)))[(nBins_discard + 1):nBins_orig] # red
      colors_categories.l[["TF"]]  = c(color_gradient[1], color_gradient[nBins_real]) 
      colors_categories.l[["TF"]]  = color_gradient 
      symbols_categories.l[["TF"]] = c(15,NA,15)
    }
    
    
    if (!is.null(vertice_color_peaks)) {
      color_gradient = rev(colorspace::sequential_hcl(nBins_orig, h = 135, c = 45, l = c(35, 95)))[(nBins_discard + 1):nBins_orig]  # green
      colors_categories.l[["PEAK"]] = c(color_gradient[1], color_gradient[nBins_real])
      colors_categories.l[["PEAK"]] = color_gradient
      symbols_categories.l[["PEAK"]] = c(21,NA,21)
    }
    
    
    if (!is.null(vertice_color_genes)) {
      color_gradient = rev(colorspace::sequential_hcl(nBins_orig, h = 260, c = 80, l = c(30, 90)))[(nBins_discard + 1):nBins_orig] # blue
      colors_categories.l[["GENE"]] = c(color_gradient[1], color_gradient[nBins_real]) 
      colors_categories.l[["GENE"]] = color_gradient
      symbols_categories.l[["GENE"]] = c(21,NA,21)
    }
    
    
    #plot(NULL, xlim=c(0,length(color_gradient)), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n")
    #rect(0:(length(color_gradient)-1), 0, 1:length(color_gradient), 1, col=color_gradient)
    
    
    
    
    ## EDGES ##
    
    # Construct items for igraph
    unique_TF_peak.con = grn.merged  %>%
      dplyr::select(TF, peak, TF_peak.fdr, TF_peak.r) %>%
      dplyr::mutate_if(is.factor, as.character) %>%
      dplyr::distinct()
    
    
    unique_peak_gene.con = grn.merged  %>%
      dplyr::select(peak, ENSEMBL, peak_gene.r, peak_gene.distance) %>%
      dplyr::filter(!is.na(ENSEMBL)) %>%
      dplyr::mutate_if(is.factor, as.character) %>%
      dplyr::distinct()  
    
    # Cannot have 0 rows
    edges_final = tibble::tibble(
      from   = unique_TF_peak.con$TF,
      to     = unique_TF_peak.con$peak,
      weight = 1- unique_TF_peak.con$TF_peak.fdr,
      R      = unique_TF_peak.con$TF_peak.r,
      linetype = "solid"
    )
    
    if (nrow(unique_peak_gene.con) > 0) {
      
      edges_final =  edges_final %>%
        tibble::add_row(
          from   = unique_peak_gene.con$peak,
          to     = unique_peak_gene.con$ENSEMBL,
          weight = unique_peak_gene.con$peak_gene.r,
          R      = unique_peak_gene.con$peak_gene.r,
          linetype = "solid"
        )
    } 
    
    edges_final =  edges_final %>%
      dplyr::mutate(weight_transformed = dplyr::case_when(weight < 0.2 ~ 1,
                                                          weight < 0.4 ~ 1.5,
                                                          weight < 0.6 ~ 2,
                                                          weight < 0.8 ~ 2.5,
                                                          TRUE ~ 3),
                    R_direction = dplyr::case_when(R < 0 ~ "neg", TRUE ~ "pos"),
                    color       = dplyr::case_when(R < 0 ~ "blue", TRUE ~ "red")
      )
    
    
    ## VERTICES ##
    
    shape_vertex = c("square","circle", "circle")
    names(shape_vertex) = names(colors_categories.l)
    
    
    vertices = tibble::tribble(~id,
                               ~type,
                               ~label,
                               ~color_raw,
                               ~color_bin,
                               ~color_final)
    
    ## 1. TFs ##
    
    # Make the vertices unique, so that the same peak has only one vertice 
    vertices_TFs = unique_TF_peak.con %>%
      dplyr::group_by(TF) %>%
      dplyr::summarize(label = unique(TF)) %>%
      dplyr::ungroup()
    
    if (nrow(vertices_TFs) > 0) {
      
      if (!is.null(vertice_color_TFs)) {
        
        .verifyArgument_verticeType(vertice_color_TFs)
        
        vertices_TFs = vertices_TFs %>%
          dplyr::left_join(vertice_color_TFs[[1]], by = c("TF" = vertice_color_TFs[[2]])) %>%
          dplyr::rename(color_raw = !!(vertice_color_TFs[[3]])) %>%
          dplyr::mutate(color_bin = as.character(cut(color_raw, nBins_real, labels = colors_categories.l[["TF"]], ordered_result = TRUE)))  # Transform the colors for the vertices
        
        
      } else {
        vertices_TFs = dplyr::mutate(vertices_TFs, color_raw = NA, color_bin = colors_categories.l[["TF"]])
      } 
      
      vertices = tibble::add_row(vertices, 
                                 id = vertices_TFs$TF, 
                                 type = "TF", 
                                 label = vertices_TFs$label, 
                                 color_raw = vertices_TFs$color_raw,
                                 color_bin = vertices_TFs$color_bin) 
    }
    
    
    ## 2. PEAKS ##
    
    # Make the vertices unique, so that the same peak has only one vertice 
    vertices_peaks = tibble::tibble(peak = unique(c(unique_peak_gene.con$peak, unique_TF_peak.con$peak)), label = NA)
    
    if (nrow(vertices_peaks) > 0) {
      
      if (!is.null(vertice_color_peaks)) {
        
        .verifyArgument_verticeType(vertice_color_peaks)
        
        vertices_peaks = vertices_peaks %>%
          dplyr::left_join(vertice_color_peaks[[1]], by = c("peak" = vertice_color_peaks[[2]])) %>%
          dplyr::rename(color_raw = !!(vertice_color_peaks[[3]])) %>%
          dplyr::mutate(color_bin = as.character(cut(color_raw, nBins_real, labels = colors_categories.l[["PEAK"]], ordered_result = TRUE)))  # Transform the colors for the vertices
        
      } else {
        vertices_peaks = dplyr::mutate(vertices_peaks, color_raw = NA, color_bin = colors_categories.l[["PEAK"]])
      } 
      
      vertices = tibble::add_row(vertices, 
                                 id = vertices_peaks$peak, 
                                 type = "PEAK", 
                                 label = vertices_peaks$label, 
                                 color_raw = vertices_peaks$color_raw,
                                 color_bin = vertices_peaks$color_bin) 
    }
    
    
    ## 3. GENES ##
    
    # Make the vertices unique, so that the same gene has only one vertice 
    vertices_genes = unique_peak_gene.con %>%
      dplyr::group_by(ENSEMBL) %>%
      dplyr::summarize(label = NA) %>% #, id2 = paste0(SYMBOL, collapse = ",")) %>%
      dplyr::ungroup()
    
    if (nrow(vertices_genes) > 0) {
      
      if (!is.null(vertice_color_genes)) {
        
        .verifyArgument_verticeType(vertice_color_genes)
        
        vertices_genes = vertices_genes %>%
          dplyr::left_join(vertice_color_genes[[1]], by = c("ENSEMBL" = vertice_color_genes[[2]])) %>%
          dplyr::rename(color_raw = !!(vertice_color_genes[[3]])) %>%
          dplyr::mutate(color_bin = as.character(cut(color_raw, nBins_real, labels = colors_categories.l[["GENE"]], ordered_result = TRUE)))  # Transform the colors for the vertices
        
      } else {
        vertices_genes = dplyr::mutate(vertices_genes, color_raw = NA, color_bin = colors_categories.l[["GENE"]])
      } 
      
      
      vertices = tibble::add_row(vertices, 
                                 id = vertices_genes$ENSEMBL, 
                                 type = "GENE", 
                                 label = vertices_genes$label, 
                                 color_raw = vertices_genes$color_raw,
                                 color_bin = vertices_genes$color_bin) 
      
    }
    
    
    
    vertices = vertices %>%
      dplyr::mutate(size = dplyr::case_when(type == "TF" ~ 6,
                                            type == "PEAK" ~ 3,
                                            TRUE ~ 4),
                    size_transformed = NA) 
    
    
    vertices_colorRanges = vertices %>% dplyr::group_by(.data$type) %>% dplyr::summarize(min = min(color_raw, na.rm = TRUE), max = max(color_raw, na.rm = TRUE))
    
    if (nrow(dplyr::filter(vertices_colorRanges, .data$type == "GENE")) == 0) {
      vertices_colorRanges = tibble::add_row(vertices_colorRanges, type = "GENE", min = NA, max = NA)
    }
    
    text_categories.l = list(
      "TF"   = "TF",
      "PEAK" = "PEAK",
      "GENE" = "GENE"
    )
    
    if (!is.null(vertice_color_TFs)) {
      subsetCur = dplyr::filter(vertices_colorRanges, .data$type == "TF") 
      text_categories.l[["TF"]] = c(signif(dplyr::pull(subsetCur, min),2), 
                                    paste0("TF expression (", vertice_color_TFs[[3]], ")"),
                                    signif(dplyr::pull(subsetCur, max),2)
      )
    }
    
    
    if (!is.null(vertice_color_peaks)) {
      subsetCur = dplyr::filter(vertices_colorRanges, .data$type == "PEAK") 
      text_categories.l[["PEAK"]] = c(signif(dplyr::pull(subsetCur, min),2), 
                                      paste0("Peak accessibility (", vertice_color_peaks[[3]], ")"),
                                      signif(dplyr::pull(subsetCur, max),2)
      )
    }
    
    
    if (!is.null(vertice_color_genes)) {
      subsetCur = dplyr::filter(vertices_colorRanges, .data$type == "GENE") 
      text_categories.l[["GENE"]] = c(signif(dplyr::pull(subsetCur, min),2), 
                                      paste0("Gene expression (", vertice_color_genes[[3]], ")"),
                                      signif(dplyr::pull(subsetCur, max),2)
      )
    }
    
    
    
    # TODO take already existing graph structure from GRN@graph and do not build anew here
    net <- igraph::graph_from_data_frame(d=edges_final, vertices = vertices, directed=TRUE) 
    
    # TODO: Integrate network stats: https://kateto.net/networks-r-igraph
    # Make a separate df_to_igraph function for the entwork stats
    
    
    ########### Color and Shape parameters 
    
    # TODO: https://stackoverflow.com/questions/48490378/order-vertices-within-layers-on-tripartite-igraph
    
    net <- igraph::simplify(net, remove.multiple = FALSE, remove.loops = TRUE)
    deg <- igraph::degree(net, mode="all")
    # V(net)$size <- deg*2
    igraph::V(net)$vertex_degree <-  deg*4
    igraph::V(net)$vertex.color = vertices$color_bin
    igraph::V(net)$vertex.size = vertices$size
    # https://rstudio-pubs-static.s3.amazonaws.com/337696_c6b008e0766e46bebf1401bea67f7b10.html
    # TODO: E(net)$weight <- edges_final$weight_transformed
    igraph::E(net)$color = edges_final$color
    #E(net)$width <- E(net)$weight
    #change arrow size and edge color:
    igraph::E(net)$arrow.size <- .1
    igraph::E(net)$edge.color <- edges_final$color
    # TODO: E(net)$lty = edges_final$linetype
    # TODO: E(net)$width <- 1+E(net)$weight/12
    
    l <- igraph::layout_with_fr(net)
    #test.layout <- layout_(net,with_dh(weight.edge.lengths = edge_density(net)/1000))
    
    # MyLO = matrix(0, nrow=vcount(net), ncol=2)
    # 
    # ## Horizontal position is determined by layer
    # layer <- rep(NA, length(V(net)$name))
    # layer[vertices$type == "TF"]   = 1
    # layer[vertices$type == "PEAK"] = 2
    # layer[vertices$type == "GENE"] = 3
    # MyLO[,1] = layer
    # 
    # ## Vertical position is determined by sum of sorted vertex_degree
    # for(i in 1:3) {
    #     L  = which(layer ==i)
    #     OL = order(V(net)$vertex_degree[L], decreasing=TRUE)
    #     MyLO[L[OL],2] = cumsum(V(net)$vertex_degree[L][OL])
    # }
    # 
    # layout = layout_with_sugiyama(net, layers=layer)
    # plot(net,
    #      layout=cbind(layer,layout$layout[,1]),edge.curved=0,
    #      vertex.shape=c("square","circle","square")[layer],
    #      vertex.frame.color = c("darkolivegreen","darkgoldenrod","orange3")[layer],
    #      vertex.color=c("olivedrab","goldenrod1","orange1")[layer],
    #      vertex.label.color="white",
    #      vertex.label.font=1,
    #      vertex.size=V(net)$vertex_degree,
    #      vertex.label.dist=c(0,0,0)[layer],
    #      vertex.label.degree=0)
    
    # 
    # vertex.color	 Node color
    # vertex.frame.color	 Node border color
    # vertex.shape	 One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
    # vertex.size	 Size of the node (default is 15)
    # vertex.size2	 The second size of the node (e.g. for a rectangle)
    # vertex.label	 Character vector used to label the nodes
    # vertex.label.family	 Font family of the label (e.g.“Times”, “Helvetica”)
    # vertex.label.font	 Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
    # vertex.label.cex	 Font size (multiplication factor, device-dependent)
    # vertex.label.dist	 Distance between the label and the vertex
    # vertex.label.degree	 The position of the label in relation to the vertex, where 0 right, “pi” is left, “pi/2” is below, and “-pi/2” is above
    # EDGES	 
    # edge.color	 Edge color
    # edge.width	 Edge width, defaults to 1
    # edge.arrow.size	 Arrow size, defaults to 1
    # edge.arrow.width	 Arrow width, defaults to 1
    # edge.lty	 Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
    # edge.label	 Character vector used to label edges
    # edge.label.family	 Font family of the label (e.g.“Times”, “Helvetica”)
    # edge.label.font	 Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
    # edge.label.cex	 Font size for edge labels
    # edge.curved	 Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
    # arrow.mode	 Vector specifying whether edges should have arrows,
    # possible values: 0 no arrow, 1 back, 2 forward, 3 both
    
    #par(mar=c(5, 4, 4, 2) + 0.1)
    par(mar=c(7,0,2,0) + 0.1)
    #par(xpd = TRUE, mar = par()$mar + c(-5,7,0,0))
    
    # Calling plot.new() might be necessary here
    plot(net, layout=l,
         #edge.arrow.size= 0.4, 
         # TODO: edge.arrow.width= E(net)$weight, 
         # TODO: edge.width = E(net)$weight, 
         edge.font= 2,
         # TODO: edge.lty = E(net)$lty,
         vertex.size= igraph::V(net)$vertex.size,
         edge.color = igraph::E(net)$color,
         vertex.color=igraph::V(net)$vertex.color,
         vertex.label.font=1, 
         vertex.label.cex = 0.3, 
         vertex.label.family="Helvetica", 
         vertex.label.color = "black",
         #vertex.label.color = col_vector[V(net)$type],
         vertex.label.degree= -pi/2,
         vertex.label=igraph::V(net)$label,
         vertex.label.dist= 0,
         vertex.shape = shape_vertex[igraph::V(net)$type],
         main = title
    )
    
    text_final    = c(c(paste0(sapply(text_categories.l, '[[', 1), " (sel. min.)"), "negative (fixed color)"),   
                      c(sapply(text_categories.l, '[[', 2), "Correlation between vertices"),    
                      c(paste0(sapply(text_categories.l, '[[', 3), " (sel. max.)"), "positive (fixed color)")#, 
                      #c("Negative correlation", "Positive correlation", "bla")
    )
    symbols_final = c(c(sapply(symbols_categories.l, '[[', 1), 20),
                      c(sapply(symbols_categories.l, '[[', 2), NA), 
                      c(sapply(symbols_categories.l, '[[', 3), 20)#,
                      #c(20,20,20)
    )
    colors_final  = c(c(sapply(colors_categories.l , '[[', 1), "blue"), 
                      c(rep(NA,3), NA), 
                      c(sapply(colors_categories.l , '[[', nBins_real), "red")#,
                      #c("red","blue","green")
    )
    
    # https://stackoverflow.com/questions/24933703/adjusting-base-graphics-legend-label-width
    #labels = c("6.4", "blaaaaaaaaaaaaaaaaaaaaaaaa", "6.4")
    
    # xcoords <- c(0, 10, 30, 60)
    # secondvector <- (1:length(legtext))-1
    # textwidths <- xcoords/secondvector # this works for all but the first element
    # textwidths[1] <- 0
    
    #par(xpd = TRUE, mar = c(6,0,2,0))
    legend(x= "bottom", text_final, 
           pch=symbols_final,
           col=colors_final, 
           pt.bg=colors_final, 
           pt.cex=1, cex=.8, bty="n", xpd = TRUE, ncol=3, xjust = 0.5, yjust = 0.5, 
           inset=c(0,-0.1) #,
           # text.width=c(0, 0.5,0.8,0.5)
    ) # , ) #  y.intersp=0.8,x.intersp=0.8, text.width=strwidth(labels, units="user")
    par(mar=c(5, 4, 4, 2) + 0.1)
    #bottom, left, top, and right.
  }
  
  if (!is.null(file)) {
    grDevices::dev.off()
  }
  
  .printExecutionTime(start)
  
}



.verifyArgument_verticeType <- function(vertice_color_list) {
  
  checkmate::assertList(vertice_color_list, len = 3)
  checkmate::assertDataFrame(vertice_color_list[[1]])
  checkmate::assertSubset(c(vertice_color_list[[2]], vertice_color_list[[3]]), colnames(vertice_color_list[[1]]))
}
