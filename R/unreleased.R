
plotCorrelations <- function(GRN, TF_peak.fdr.threshold = 0.2, TF_peak.r.abs.threshold = 0.3, randomized = TRUE,
                             topn = 100, maxAll = 30000, TF.names = NULL, peak.IDs = NULL,
                             corMethod = "pearson", filePDF) {
    
    checkmate::assertChoice(sortby, colnames(GRN@connections$TF_peaks$`0`$main))
    
    con.filt = GRN@connections$TF_peaks$`0`$main %>%
        dplyr::mutate(TF_peak.r_abs = abs(TF_peak.r)) %>%
        dplyr::filter(.data$TF_peak.fdr <= TF_peak.fdr.threshold,
                      .data$TF_peak.r_abs >= TF_peak.r.abs.threshold)
    
    
    if (!is.null(TF.names)) {
        con.filt = dplyr::filter(con.filt, .data$Tf.name %in% TF.names)
    }
    
    if (!is.null(peak.IDs)) {
        con.filt = dplyr::filter(con.filt, .data$peak.ID %in% peak.IDs)
    }
    
    if (randomized) {
        con.filt = dplyr::slice_sample(con.filt, n = topn) 
    } else {
        con.filt = con.filt %>%
            dplyr::arrange(desc(TF_peak.r_abs)) %>%
            dplyr::slice_head(con.filt, n = topn)
    }

    pdf(filePDF, height = 10)
    
    for (i in seq_len(topn)) {
        
        TF.name = con.filt$TF.name[i]
        TF.ENSEMBL = GRN@annotation$TFs$TF.ENSEMBL[which(GRN@annotation$TFs$TF.name == TF.name)]
        peak.peakID = con.filt$peak.ID[i]
        corCalc = con.filt$TF_peak.r[i]
        fdrCur = con.filt$TF_peak.fdr[i]
        
        cor_cur = cor(GRN@data$RNA$counts[TF.ENSEMBL, ], GRN@data$peaks$counts[peak.peakID,], method = corMethod)
        
        data.df = tibble(TF.exp.norm = GRN@data$RNA$counts[TF.ENSEMBL, ], peak.acc.norm = GRN@data$peaks$counts[peak.peakID,])
        
        # For testing only
        stopifnot(abs(cor_cur - corCalc) < 0.01)
        
        colorscale = scale_fill_gradientn(
            colors = RColorBrewer::brewer.pal(9, "YlGnBu"),
            values = c(0, exp(seq(-5, 0, length.out = 100))))
        
        g1 = ggplot(data.df, aes(.data$TF.exp.norm, .data$peak.acc.norm)) + 
            geom_smooth(method = "lm", formula = "y ~ x", col = "black") + 
            geom_hex(bins = 50) + colorscale  +
            xlim(-0.01,NA) + ylim(-0.01,NA) + 
            ggtitle(paste0("n = ", nrow(data.df), " (all), Cor = ", round(corCalc,2)))
        
        # Remove entries with pure zeroes
        data.filt.df = dplyr::filter(data.df, .data$TF.exp.norm > 0, .data$peak.acc.norm > 0)
        
        cor_cur = cor(dplyr::pull(data.filt.df, .data$TF.exp.norm), dplyr::pull(data.filt.df, .data$peak.acc.norm), method = corMethod)
        
        g2 = ggplot(data.filt.df, aes(TF.exp.norm, peak.acc.norm)) + 
            geom_smooth(method = "lm", formula = "y ~ x", col = "black") + 
            geom_hex(bins = 50) + colorscale  +
            xlim(-0.01,NA) + ylim(-0.01,NA) + 
            ggtitle(paste0("n = ", nrow(data.filt.df), " (only >0 for both TF expr. & peak acc.), Cor = ", round(cor_cur,2)))    
        
        data.filt2.df = dplyr::filter(data.df, TF.exp.norm > 0 | peak.acc.norm > 0)
        
        cor_cur = cor(dplyr::pull(data.filt2.df, .data$TF.exp.norm), dplyr::pull(data.filt2.df, .data$peak.acc.norm), method = corMethod)
        
        g3 = ggplot(data.filt2.df, aes(TF.exp.norm, peak.acc.norm)) + 
            geom_smooth(method = "lm", formula = "y ~ x", col = "black") + 
            geom_hex(bins = 50) + colorscale  +
            xlim(-0.01,NA) + ylim(-0.01,NA) + 
            ggtitle(paste0("n = ", nrow(data.filt2.df), " (only excluding (0,0) for both TF expr. & peak acc.), Cor = ", round(cor_cur,2)))    
        
        
        mainTitle = paste0("TF: ", TF.name, "(", TF.ENSEMBL, "), peak: ", peak.peakID, " (FDR: ",  round(fdrCur, 2), ")")
        
        plots_all =  g1 / g2 / g3 + 
            patchwork::plot_annotation(title = mainTitle, theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
        plot(plots_all)
        # 
        # r1 = rnorm(1000)
        # r2 = rnorm(1000, sd = 1)
        # plot(r1,r2, main = paste0("Cor = ", round(cor(r1,r2, method = corMethod), 2)))
        # 
        # r1 = c(r1, rep(0,5000))
        # r2 = c(r2, rep(0,5000))
        # plot(r1,r2, main = paste0("Cor = ", round(cor(r1,r2, method = corMethod), 2)))
    }

    if (maxAll > 0) {
        
        con.filt = GRN@connections$TF_peaks$`0`$main
        nRowsReal = min(maxAll, nrow(con.filt))
        
        randomRows = sample.int(nrow(con.filt), nRowsReal)
        
        res.df = data.frame(matrix(NA, nrow = nRowsReal, ncol = 5))
        colnames(res.df) = c("TF.name", "peak.peakID", "cor_all", "cor_nonZeroAll", "cor_nonZero")
        #res.df = tribble(~TF.name, ~peak.peakID, ~cor_all, ~cor_nonZeroAll, ~cor_nonZero)
        for (i in 1:nRowsReal) {
            
            if (i %% 1000 == 1) futile.logger::flog.info(paste0("Row ", i, " out of ", nRowsReal))
            index =randomRows[i]
            TF.name = con.filt$TF.name[index]
            TF.ENSEMBL = GRN@annotation$TFs$TF.ENSEMBL[which(GRN@annotation$TFs$TF.name == TF.name)]
            peak.peakID = con.filt$peak.ID[index]
            corCalc = con.filt$TF_peak.r[index]
            
            
            #cor_cur1 = cor(GRN@data$RNA$counts[TF.ENSEMBL, ], GRN@data$peaks$counts[peak.peakID,], method = corMethod)
            x_orig = GRN@data$RNA$counts[TF.ENSEMBL, ]
            y_orig = GRN@data$peaks$counts[peak.peakID,]
            keep = which(x_orig > 0 & y_orig > 0)
            x = x_orig[keep]
            y = y_orig[keep]
            
            #data.df = tibble(TF.exp.norm = GRN@data$RNA$counts[TF.ENSEMBL, ], peak.acc.norm = GRN@data$peaks$counts[peak.peakID,])
            
            #data.filt.df = dplyr::filter(data.df, .data$TF.exp.norm > 0, .data$peak.acc.norm > 0)
            
            #cor_cur2 = cor(dplyr::pull(data.filt.df, .data$TF.exp.norm), dplyr::pull(data.filt.df, .data$peak.acc.norm), method = corMethod)
            cor_cur2 = cor(x, y, method = corMethod)
            
            keep = which(x_orig > 0 | y_orig > 0)
            x = x_orig[keep]
            y = y_orig[keep]
            
            # data.filt2.df = dplyr::filter(data.df, .data$TF.exp.norm > 0 | peak.acc.norm > 0)
            cor_cur3 = cor(x, y, method = corMethod)
            # cor_cur3 = cor(dplyr::pull(data.filt2.df, .data$TF.exp.norm), dplyr::pull(data.filt2.df, .data$peak.acc.norm), method = corMethod)
            res.df[i,] = c(as.character(TF.name), as.character(peak.peakID), as.numeric(corCalc), as.numeric(cor_cur2), as.numeric(cor_cur3))
            
            # res.df = add_row(res.df, 
            #                  TF.name = TF.name, peak.peakID = peak.peakID, 
            #                  cor_all = corCalc, cor_nonZeroAll = cor_cur2, cor_nonZero = cor_cur3)
        }
        
        res.df$cor_all= as.numeric(res.df$cor_all)
        res.df$cor_nonZeroAll= as.numeric(res.df$cor_nonZeroAll)
        res.df$cor_nonZero= as.numeric(res.df$cor_nonZero)
        
        cor(res.df$cor_all, res.df$cor_nonZeroAll)
        cor(res.df$cor_all, res.df$cor_nonZero)
        
        res.df = res.df %>%
            dplyr::mutate(diff_all_nonZeroAll = cor_all - cor_nonZeroAll,
                          diff_all_nonZero = cor_all - cor_nonZero)
        
        plot(ggplot(res.df, aes(.data$diff_all_nonZeroAll)) + 
                 geom_histogram(binwidth = 0.05) + 
                 ggtitle(paste0("No. rows: ", nRowsReal)) + 
                 xlim(-1,1) + xlab("Only pairs for which both values are > 0")
        )
        
        
        plot(ggplot(res.df, aes(.data$diff_all_nonZero)) + 
                 geom_histogram(binwidth = 0.05) + 
                 ggtitle(paste0("No. rows: ", nRowsReal))  + 
                 xlim(-1,1) + xlab("Only pairs for which at least one point is > 0")
        )
        
    }
    
    dev.off()
    
    GRN
    
}


plot_peaks <- function (GRN) {
    
    genomeAssembly = GRN@config$parameters$genomeAssembly
    query   = .constructGRanges(GRN@data$peaks$counts_metadata, 
                                seqlengths = .getChrLengths(genomeAssembly), 
                                genomeAssembly)
    
    x = GenomicDistributions::calcChromBinsRef(query, genomeAssembly)
    plotChromBins(x)
}