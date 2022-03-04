######## GRN graph functions ########
#' Builds a graph out of a set of connections
#'
#' @template GRN
#' @param model_TF_gene_nodes_separately \code{TRUE} or \code{FALSE}.  Default \code{FALSE}. TODO description follows 
#' @param allowLoops \code{TRUE} or \code{FALSE}.  Default \code{FALSE}. Allow loops in the network (i.e., a TF that regulates itself)
#' @param removeMultiple \code{TRUE} or \code{FALSE}.  Default \code{FALSE}. Remove loops with the same start and end point? This can happen if multiple TF originate from the same gene, for example.
#' @param directed \code{TRUE} or \code{FALSE}.  Default \code{FALSE}. Should the network be directed?
#' @template forceRerun
#' @export
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' GRN = build_eGRN_graph(GRN, forceRerun = FALSE)
#' @return The same \code{\linkS4class{GRN}} object.
build_eGRN_graph <- function(GRN, model_TF_gene_nodes_separately = FALSE, 
                             allowLoops = FALSE, removeMultiple = FALSE, directed = FALSE, forceRerun = FALSE) {
    
    checkmate::assertClass(GRN, "GRN")
    checkmate::assertFlag( model_TF_gene_nodes_separately)
    checkmate::assertFlag(allowLoops)
    checkmate::assertFlag(removeMultiple)
    checkmate::assertFlag(directed)
    checkmate::assertFlag(forceRerun)
    
    # This function returns the tf-peak-gene and the tf-gene graphs in dataframe format
    # tf-peak-gene graph is weighted (r), tf-gene graph is unweighted
    
    if (is.null(GRN@graph$TF_gene) | is.null(GRN@graph$TF_peak_gene) | forceRerun) {
        
        # Should the TF nodes and gene nodes represent the same or different nodes? 
        # If set to TRUE, the new default, self-loops can happen and the graph is not strictly tripartite anymore
        if (model_TF_gene_nodes_separately) {
            TF_peak.df = GRN@connections$all.filtered[["0"]] %>%
                dplyr::filter(!is.na(gene.ENSEMBL)) %>% 
                dplyr::select(c("TF.name", "peak.ID", "TF.ENSEMBL", "TF_peak.r")) %>%
                stats::na.omit() %>% 
                dplyr::mutate(V2_name = NA) %>%
                unique() %>%
                dplyr::rename(V1 = TF.name, V2 = peak.ID, V1_name = TF.ENSEMBL, r = TF_peak.r) %>%
                dplyr::mutate_at(c("V1","V2"), as.vector)
        } else {
            # Get Ensembl ID for TFs here to make a clean join and force a TF that is regulated by a peak to be the same node
            TF_peak.df = GRN@connections$all.filtered[["0"]] %>%
                dplyr::filter(!is.na(gene.ENSEMBL)) %>% 
                dplyr::select(c("TF.ENSEMBL", "peak.ID", "TF.name", "TF_peak.r")) %>%
                stats::na.omit() %>% 
                dplyr::mutate(V2_name = NA) %>%
                unique() %>%
                dplyr::rename(V1 = TF.ENSEMBL, V2 = peak.ID, V1_name = TF.name, r = TF_peak.r) %>%
                dplyr::mutate_at(c("V1","V2"), as.vector)
        }
        
        
        peak_gene.df = GRN@connections$all.filtered[["0"]][,c("peak.ID", "gene.ENSEMBL", "gene.name", "peak_gene.r")] %>% 
            stats::na.omit() %>% 
            dplyr::mutate(V1_name = NA) %>%
            unique() %>%
            dplyr::rename(V1 = peak.ID, V2 = gene.ENSEMBL, V2_name = gene.name, r = peak_gene.r) %>%
            dplyr::mutate_at(c("V1","V2"), as.vector)
        
        TF_peak_gene.df = dplyr::bind_rows(list(`tf-peak`=TF_peak.df, `peak-gene` = peak_gene.df), .id = "connectionType") %>%
            dplyr::select(V1, V2, V1_name, V2_name, r, connectionType)
        
        TF_gene.df = dplyr::inner_join(TF_peak.df, peak_gene.df, by = c("V2"="V1"), suffix = c(".TF_peak", ".peak_gene")) %>% 
            dplyr::select(V1, V2.peak_gene, V1_name.TF_peak, V2_name.peak_gene) %>%
            dplyr::rename(V1_name = V1_name.TF_peak, V2 = V2.peak_gene, V2_name = V2_name.peak_gene) %>%
            dplyr::distinct() %>%
            dplyr::mutate(connectionType = "tf-gene") 
        
        # If the graph is NOT directed, retrieving the graph structure as data frame may result in V1 and V2 switched
        # This happens when TF-TF interactions occur. The order (V1, V2) is irrelevant for an undirected graph anyway
        futile.logger::flog.info(paste0("Building TF-peak-gene graph..."))
        GRN@graph$TF_peak_gene = list(table = TF_peak_gene.df,
                                      graph = .buildGraph(TF_peak_gene.df, 
                                                          directed = directed, 
                                                          allowLoops = allowLoops, 
                                                          removeMultiple = removeMultiple))
        
        futile.logger::flog.info(paste0("Building TF-gene graph..."))
        GRN@graph$TF_gene      = list(table = TF_gene.df,
                                      graph = .buildGraph(TF_gene.df, 
                                                          directed = directed, 
                                                          allowLoops = allowLoops, 
                                                          removeMultiple = removeMultiple))
        
        
        GRN@graph$parameters = list()
        GRN@graph$parameters$directed       = directed
        GRN@graph$parameters$allowLoops     = allowLoops
        GRN@graph$parameters$removeMultiple = removeMultiple
        
    }
    
    
    GRN
    
}

.buildGraph <- function(df, directed, allowLoops, removeMultiple = FALSE, silent = FALSE) {
    
    # Remove V1_name and V2_name as igraph treats additional columns as edge attribute, which can become messed up as it here refers to vertice attribute
    df_mod = df %>% dplyr::select(-V1_name, -V2_name)
    
    TF_vertices = df %>%
        dplyr::select(V1, V1_name) %>% 
        dplyr::rename(nodeID = V1) %>%
        dplyr::distinct() %>%
        dplyr::group_by(nodeID) %>%
        dplyr::summarise(names_TF_all = paste0(V1_name, collapse="|"),
                         nTF = dplyr::n(),
                         isTF = TRUE, .groups = "keep") %>%
        dplyr::ungroup()
    
    gene_vertices = df %>%
        dplyr::select(V2, V2_name) %>% 
        dplyr::distinct() %>%
        dplyr::mutate(isGene = TRUE) %>%
        dplyr::rename(names_gene = V2_name, nodeID = V2) %>%
        dplyr::ungroup()
    
    # Combine vertex metadata
    vertexMetadata = dplyr::full_join(TF_vertices, gene_vertices, by = "nodeID")
    
    # Fix the isTF column
    vertexMetadata$isTF[is.na(vertexMetadata$isTF)] = FALSE
    
    graph = igraph::graph_from_data_frame(d=df_mod, directed = directed, vertices = vertexMetadata)
    
    if (!igraph::is_simple(graph)) {
        if (!silent) futile.logger::flog.info(paste0(" Graph contains either loops and/or multiple edges. A simplification is possible."))
        
        .printLoopsGraph(df, silent = silent)
        .printMultipleEdges(df, silent = silent)
        
        if (removeMultiple | !allowLoops) {
            if (!silent) futile.logger::flog.info(paste0(" Simplify graph..."))
            graph <- igraph::simplify(graph, remove.multiple = removeMultiple, remove.loops = !allowLoops)
        } else {
            if (!silent) futile.logger::flog.info(paste0(" Not doing any graph simplification, see the parameters removeMultiple and allowLoops to change it."))
        }
        
    }
    
    .printGraphSummary(GRN, graph, silent = silent)

    if (!silent) futile.logger::flog.info(paste0(" Done. Graphs are saved in GRN@graph"))
    
    graph
}

.printLoopsGraph <- function(graph_table, silent = FALSE) {
    
    loop_vertices = graph_table %>%
        dplyr::filter(V1 == V2) %>%
        dplyr::mutate(V1_name_combined = paste0(V1, " (", V1_name, ")")) %>%
        dplyr::pull(V1_name_combined)
    
    if (length(loop_vertices) > 0) {
        if (!silent) futile.logger::flog.info(paste0(" The following nodes / vertices have loop edges (TF regulating itself):\n", paste0(loop_vertices, collapse = ", ")))
    }
    
}

.printMultipleEdges <- function(graph_table, silent = FALSE) {
    
    
    multipleEdges = graph_table %>%
        dplyr::group_by(V1,V2) %>% 
        dplyr::summarize(n = dplyr::n(), .groups = "keep") %>% 
        dplyr::ungroup() %>%
        dplyr::filter(.data$n>1)
    
    if (nrow(multipleEdges) > 0) {
        if (!silent) futile.logger::flog.info(paste0(" ", nrow(multipleEdges), " edges have the same vertices. This is often caused by multiple TF belonging to the same gene ID."))
    }
    
}


.printGraphSummary <- function(GRN, graph, silent = FALSE) {
    if (!silent) futile.logger::flog.info(paste0(" Graph summary:"))
    nVertices = length(igraph::V(graph))
    nEdges = length(igraph::E(graph))
    
    if (!silent) futile.logger::flog.info(paste0("  Nodes (vertices): ", nVertices))
    if (!silent) futile.logger::flog.info(paste0("  Edges: ", nEdges))
    
}


#' Perform all network-related statistical and descriptive analyses, including community and enrichment analyses. 
#' 
#' A convenience function that calls all network-related functions in one-go, using selected default parameters and a set of adjustable ones also. For full adjustment, run the individual functions separately.
#'
#' @inheritParams calculateGeneralEnrichment
#' @inheritParams plotCommunitiesStats
#' @inheritParams plotCommunitiesEnrichment
#' @inheritParams calculateCommunitiesStats
#' @export
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' GRN = performAllNetworkAnalyses(GRN, forceRerun = FALSE)
#' @return The same \code{\linkS4class{GRN}} object, with added data from this function.
performAllNetworkAnalyses <- function(GRN, ontology = c("GO_BP", "GO_MF"), 
                                      algorithm = "weight01", statistic = "fisher",
                                      background = "neighborhood", 
                                      clustering = "louvain",
                                      communities = seq_len(10), display = "byRank",
                                      topnGenes = 20, topnTFs = 20,
                                      display_pAdj = FALSE,
                                      outputFolder = NULL,
                                      forceRerun = FALSE) {
    
    start = Sys.time()
    GRN = .addFunctionLogToObject(GRN)
    
    GRN = build_eGRN_graph(GRN, model_TF_gene_nodes_separately = FALSE, allowLoops = FALSE, directed = FALSE, removeMultiple = FALSE)
    
    GRN = plotGeneralGraphStats(GRN, outputFolder = outputFolder, forceRerun = forceRerun) 
    
    GRN = calculateGeneralEnrichment(GRN, ontology = ontology, algorithm = algorithm, statistic = statistic, 
                                     background = background, forceRerun = forceRerun)
    GRN = plotGeneralEnrichment(GRN, outputFolder = outputFolder, display_pAdj = display_pAdj, forceRerun = forceRerun) 
    
    
    GRN = calculateCommunitiesStats(GRN, clustering = clustering, forceRerun = forceRerun)
    
    GRN = plotCommunitiesStats     (GRN, outputFolder = outputFolder, display = display, communities = communities, 
                                    forceRerun = forceRerun, topnGenes = topnGenes, topnTFs = topnTFs)
    
    GRN = calculateCommunitiesEnrichment(GRN, ontology = ontology, algorithm = algorithm, statistic = statistic, 
                                         background = background, forceRerun = forceRerun)
    
    GRN = plotCommunitiesEnrichment(GRN, outputFolder = outputFolder, display = display, communities = communities, 
                                    display_pAdj = display_pAdj, forceRerun = forceRerun)
    
    GRN = calculateTFEnrichment(GRN, ontology = ontology, algorithm = algorithm, statistic = statistic,
                                background = background, pAdjustMethod = "BH",
                                forceRerun = forceRerun)

    GRN = plotTFEnrichment(GRN, display_pAdj = display_pAdj, forceRerun = forceRerun)
    
    
    .printExecutionTime(start)
    
    GRN
    
}



# Retrieve set of background genes (as vector) used for enrichment analyses from a GRN object
.getBackgroundGenes <-function(GRN, type = "neighborhood") {
    
    checkmate::assertSubset(type, c("all_annotated", "all_RNA", "neighborhood"))
    
    if (type == "all_annotated") {
        
        backgroundGenes = geneAnnotation[[GRN@config$parameters$genomeAssembly]]$gene.ENSEMBL
        
    } else if (type == "all_RNA") {
        
        if (checkmate::testClass(GRN@data$RNA$counts_orig, "DESeqDataSet")) {
            backgroundGenes = rownames(GRN@data$RNA$counts_orig)
        } else {
            backgroundGenes = GRN@data$RNA$counts_orig$ENSEMBL
        }
        
        
    } else if (type == "neighborhood") {
        
        # Retrieve only those who are in the neighborhood of genes
        backgroundGenes = levels(GRN@connections$peak_genes[["0"]]$gene.ENSEMBL)
    }
    
    
    backgroundGenes 
    
}


#' Run an enrichment analysis for the genes in the filtered \code{\linkS4class{GRN}}
#' 
#' This function runs an enrichment analysis for the genes in the filtered network.
#' @template GRN
#' @param ontology Character vector of ontologies. Default c("GO_BP", "GO_MF"). Valid values are "GO_BP", "GO_MF", "GO_CC", "KEGG", "DO", and "Reactome", referring to GO Biological Process, GO Molecular Function, GO Cellular Component, KEGG, Disease Ontology, and Reactome Pathways. The GO enrichments 
#' @param algorithm Character. Default "weight01". One of: "classic", "elim", "weight", "weight01", "lea", "parentchild." Only relevant if ontology is GO related (GO_BP, GO_MF, GO_CC), ignored otherwise. Name of the algorithm that handles the GO graph structures. Valid inputs are those supported by the topGO library.
#' @param statistic Character. Default "fisher". One of: "fisher", "ks", "t", "globaltest", "sum", "ks.ties". Statistical test to be used. Only relevant if ontology is GO related (GO_BP, GO_MF, GO_CC), and valid inputs are those supported by the topGO library, ignored otherwise. For the other ontologies the test statistic is always Fisher. 
#' @param background Character. Default "neighborhood". One of: "all_annotated", "all_RNA", "neighborhood". Set of genes to be used to construct the background for the enrichment analysis. This can either be all annotated genes in the reference genome (all_annotated), all differentially expressed genes (all_RNA), or all the genes that are within the neighbourhood of a peak in the GRN (neighbourhood)
#' @param pAdjustMethod Character. Default "BH". One of: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr". This parameter is only relevant for the following ontologies: KEGG, DO, Reactome. For the other ontologies, the algorithm serves as an adjustment.
#' @template forceRerun
#' @return The same \code{\linkS4class{GRN}} object, with the enrichment results stored in the \code{stats$Enrichment$general} slot.
#' @seealso \code{\link{plotGeneralEnrichment}}
#' @seealso \code{\link{calculateCommunitiesEnrichment}}
#' @seealso \code{\link{plotCommunitiesEnrichment}}
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN =  loadExampleObject()
#' GRN =  calculateGeneralEnrichment(GRN, ontology = "GO_BP", forceRerun = FALSE)
#' @export
#' @import topGO
#' @import BiocManager
calculateGeneralEnrichment <- function(GRN, ontology = c("GO_BP", "GO_MF"), 
                                       algorithm = "weight01", statistic = "fisher",
                                       background = "neighborhood",  pAdjustMethod = "BH", forceRerun = FALSE) {
    
    start = Sys.time()
    
    checkmate::assertClass(GRN, "GRN")
    GRN = .addFunctionLogToObject(GRN)
    
    checkmate::assertSubset(ontology , c("GO_BP", "GO_MF", "GO_CC", "KEGG", "DO", "Reactome"))
    checkmate::assertSubset(algorithm , topGO::whichAlgorithms())
    checkmate::assertSubset(statistic , topGO::whichTests())
    checkmate::assertSubset(background, c("all_annotated", "all_RNA", "neighborhood"))
    checkmate::assertSubset(pAdjustMethod, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"))
    checkmate::assertFlag(forceRerun)
    
    
    
    
    mapping = .getGenomeObject(GRN@config$parameters$genomeAssembly, type = "packageName")
    backgroundGenes = .getBackgroundGenes(GRN, type = background)
    
    futile.logger::flog.info(paste0("Calculating general enrichment statistics. This may take a while"))
    
    if (is.null(GRN@stats[["Enrichment"]][["general"]])) {
        GRN@stats[["Enrichment"]][["general"]] = list()
    }
    
    for (ontologyCur in ontology) {
        
        if (is.null(GRN@stats[["Enrichment"]][["general"]][[ontologyCur]]) | forceRerun) {
            
            # run general enrichment analysis and store tabulated results in GRN object
            # Only use the "targets", i.e., genes as foreground because it would artificially enrich for TF terms, such as "DNA-binding" "transcription activation" type terms.
            GRN@stats[["Enrichment"]][["general"]][[ontologyCur]] =
                .runEnrichment(foreground = GRN@graph$TF_gene$table$V2, 
                               background = backgroundGenes,
                               backgroundStr = background, 
                               ontology = ontologyCur, algorithm = algorithm, statistic = statistic,
                               mapping = mapping,
                               pAdjustMethod =  pAdjustMethod)
            
            futile.logger::flog.info(paste0("Result stored in GRN@stats$Enrichment$general$", ontologyCur, "$results"))
            
        } else {
            
            futile.logger::flog.info(paste0("Results already found / previously calculated. Not re-running as forceRerun = FALSE"))
        }
        
    }
    
    
    
    
    
    .printExecutionTime(start, prefix = "")
    
    GRN
}


.getBackgroundGenes <-function(GRN, type = "neighborhood") {
    # Retrieve set of background genes (as vector) used for enrichment analyses from a GRN object
    checkmate::assertSubset(type, c("all_annotated", "all_RNA", "neighborhood"))
    if (type == "all_annotated") {
        backgroundGenes = geneAnnotation[[GRN@config$parameters$genomeAssembly]]$gene.ENSEMBL
    } else if (type == "all_RNA") {
        if (checkmate::testClass(GRN@data$RNA$counts_orig, "DESeqDataSet")) {
            backgroundGenes = rownames(GRN@data$RNA$counts_orig)
        } else {
            backgroundGenes = GRN@data$RNA$counts_orig$ENSEMBL
        }
    } else if (type == "neighborhood") {
        # Retrieve only those who are in the neighborhood of genes
        backgroundGenes = levels(GRN@connections$peak_genes[["0"]]$gene.ENSEMBL)
    }
    backgroundGenes 
}


.combineEnrichmentResults <- function(GRN, type, ontology, p, nSignificant, display_pAdj) {
    
    if (type == "byCommunity") {
        idMerge = "community"
    } else if (type == "byTF") {
        idMerge = "TF.name"
    }
    
    # Merge all community-specific results to one data frame
    resultsCombined.df = suppressWarnings(GRN@stats[["Enrichment"]][[type]] %>%
        lapply(function(x) {x[[ontology]]$results}) %>%
        dplyr::bind_rows(.id = idMerge) %>%
        dplyr::select(-tidyselect::starts_with("topG")) %>%
        tibble::as_tibble())
    
    # p-adjust only available for non-GO ontologies
    if (display_pAdj && !stringr::str_starts("GO_", ontology)) {
        resultsCombined.df$pval = resultsCombined.df$p.adjust
    }
    
    # Add general enrichment
    
    if (is.null(GRN@stats$Enrichment$general[[ontology]]$results)) {
        message = paste0("Could not find enrichment results for general enrichment for ontology ", ontology, ".. Please (re)run the function calculateGeneralEnrichment for the ontology ", ontology)
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    enrichmentGeneral = GRN@stats$Enrichment$general[[ontology]]$results %>%
        dplyr::mutate({{idMerge}} := "all") %>%
        dplyr::select(colnames(resultsCombined.df))
    
    # get enriched terms from general enrichment and make sure it is kept for the community enrichment
    enrichedTermsGeneral = enrichmentGeneral %>% 
        dplyr::filter(pval <= p, Found >= nSignificant) %>% 
        dplyr::pull(ID)
    
    enrichedTermsGrouped = resultsCombined.df %>% 
        dplyr::filter(pval <= p, Found >= nSignificant) %>% 
        dplyr::pull(ID)
    
    all.df = resultsCombined.df %>%
        rbind(enrichmentGeneral) %>%
        dplyr::mutate(ID = as.factor(ID),
                      pval =as.numeric(gsub(">|<", "", pval))) %>%
        dplyr::filter(pval <= p & (Found >= nSignificant | ID %in% c(enrichedTermsGeneral, enrichedTermsGrouped)))

    all.df[, idMerge] = as.factor(all.df[, idMerge, drop = TRUE])
    
    all.df
        
}

.checkEnrichmentCongruence_general_community <-function(GRN, type = "community") {
    
    allOntologiesGeneral = sort(names(GRN@stats$Enrichment$general))
    
    if (type == "community") {
        allOntologiesGroup1 = sort(names(GRN@stats$Enrichment$byCommunity[[1]]))
    } else if (type == "TF") {
        allOntologiesGroup1 = sort(names(GRN@stats$Enrichment$byTF[[1]]))
    }
   
    if (!identical(allOntologiesGeneral, allOntologiesGroup1)) {
        message = paste0("General enrichment and ", type, " enrichment do not have the same ontologies precalculated (",
                         paste0(allOntologiesGeneral, collapse = " & "), " vs. ", paste0(allOntologiesGroup1, collapse = "&"), ")")
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    allOntologiesGeneral
        
}

.runEnrichment <- function(foreground, background, backgroundStr, database = "GO", ontology, 
                           description = "Enrichment Analysis",
                           algorithm="weight01", statistic = "fisher", mapping, pAdjustMethod = "BH", minGSSize = 0, maxGSSize = 5000){
    
    
    result.list = list()
    # Implementation change: Allow only one ontology term here, and force the calling function to handle multiple ontologies.
    # Advantage: Prevents recalculation if enrichment for ontology has already been calculated
    checkmate::assertCharacter(ontology, len = 1)
    
    
    foreground = as.character(foreground) %>% unique()
    background = as.character(background) %>% unique()
    
    nForeground = length(foreground)
    nBackground = length(background)
    

    if (ontology %in% c("KEGG", "DO", "Reactome")) {
        # the ENSEMBL IDs will need to be mapped to Enntrez IDs for these ontologies
        
        if (statistic != "fisher") {
            statistic = "fisher"
            message = paste0("For KEGG, DO and Reacome enrichment, the parameter statistic can only be \"fisher\". it has been changed acoordingly.")
            .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
        }
        
        .checkOntologyPackageInstallation("biomaRt")
        
        if (grep(x = GRN@config$parameters$genomeAssembly, pattern = "^hg\\d\\d" )){
            dataset = "hsapiens_gene_ensembl"
        } else if (grep(x = GRN@config$parameters$genomeAssembly, pattern = "^mm\\d\\d")){
            dataset = "mmusculus_gene_ensembl"
        }
        
        errorOcured = FALSE
        
        ensembl = tryCatch({ 
            biomaRt::useEnsembl(biomart = "genes", dataset = dataset)
            
        }, error = function(e) {
            errorOcured = TRUE 
        }
        )
        
    
        foreground_entrez = tryCatch({ 
            biomaRt::getBM(mart = ensembl,
                           attributes =  "entrezgene_id",
                           filters = "external_gene_name",
                           values = GRN@annotation$genes$gene.name[
                               match(foreground, GRN@annotation$genes$gene.ENSEMBL)])[,1] %>%
                stats::na.omit() %>% as.character()
            
        }, error = function(e) {
            errorOcured = TRUE 
        }
        )
        
        background_entrez = tryCatch({ 
            biomaRt::getBM(mart = ensembl,
                           attributes = "entrezgene_id",
                           filters = "external_gene_name",
                           values = GRN@annotation$genes$gene.name[match(background, GRN@annotation$genes$gene.ENSEMBL)])[,1] %>%
                stats::na.omit() %>% as.character()
            
        }, error = function(e) {
            errorOcured = TRUE 
        }
        )
        
        if (errorOcured) {

            error_Biomart = "A temporary error occured with biomaRt::getBM or biomaRt::useEnsembl. This is often caused by an unresponsive Ensembl site and may be caused by the ontology type (e.g, it may work for the GO ontologies but not for KEGG). Try again at a later time or change ontologies. For now, this ontology has been skipped. Note that this error is not caused by GRaNIE but external services."
            .checkAndLogWarningsAndErrors(NULL, error_Biomart, isWarning = TRUE)
            return(NULL)
            
        }
        
        
    }
    
    
    geneList = factor(as.integer(unique(background) %in% unique(foreground)))
    names(geneList) = unique(background)
    
    futile.logger::flog.info(paste0("   Running enrichment analysis for ontology ", ontology, " using ", nForeground, " and ", nBackground, " genes as foreground and background (", backgroundStr, "), respectively. This may take a while."))
    
    
    if (ontology %in% c("GO_BP","GO_MF","GO_CC")){
        
        go_enrichment = suppressMessages(new("topGOdata",
                                             ontology = gsub("GO_", "", ontology),
                                             allGenes = geneList,
                                             description = description,
                                             nodeSize = 5,
                                             annot = topGO::annFUN.org,
                                             mapping = mapping, 
                                             ID = "ensembl"))
        
        result = suppressMessages(topGO::runTest(go_enrichment, algorithm = algorithm, statistic = statistic))
        # Dont trim GO terms here, happens later when plotting
        result.tbl = unique(topGO::GenTable(go_enrichment, pval = result, orderBy = "pval", numChar = 1000, 
                                            topNodes = length(topGO::score(result))) ) %>%
            dplyr::rename(ID = GO.ID, Found = Significant)  %>%      # make it more clear what Significant refers to here
            dplyr::mutate(GeneRatio = Found / nForeground)
        
        
        result.list[["results"]] = result.tbl
        
    }
    
    # Shared error message for different ontologies
    enrichmentErrorMessage = "Could not calculate enrichment, the server returned an error. This may happen for multiple reasons, for example if no gene can be mapped. The results will be set to NA."
    
    if (ontology == "KEGG"){
        
        .checkOntologyPackageInstallation("clusterProfiler")
        
        if (grep(x = GRN@config$parameters$genomeAssembly, pattern = "^hg\\d\\d" )){
            org = "hsa"
        } else if (grep(x = GRN@config$parameters$genomeAssembly, pattern = "^mm\\d\\d")){
            org = "mmu"
        }
        
        kegg_enrichment = tryCatch({ 
            clusterProfiler::enrichKEGG(
                gene = foreground_entrez,
                universe = background_entrez,
                keyType = "ncbi-geneid",
                organism = org,
                pvalueCutoff = 1,
                qvalueCutoff = 1,
                minGSSize = minGSSize,
                maxGSSize = maxGSSize,
                pAdjustMethod = pAdjustMethod)
            
        }, error = function(e) {
            message = enrichmentErrorMessage
            .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
        }, warning = function(w) {
            message = enrichmentErrorMessage
            .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
        }
        )
        
        result.list[["results"]] = .createEnichmentTable(kegg_enrichment)
        
    }
    
    if (ontology == "Reactome"){
        
        .checkOntologyPackageInstallation("ReactomePA")
        
        reactome_enrichment = tryCatch({ 
            ReactomePA::enrichPathway(
                gene=foreground_entrez,
                universe = background_entrez,
                pvalueCutoff = 1,
                qvalueCutoff = 1,
                minGSSize = minGSSize,
                maxGSSize = maxGSSize,
                pAdjustMethod = pAdjustMethod)
            
        }, error = function(e) {
            message = enrichmentErrorMessage
            .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
        }, warning = function(w) {
            # dont print anything
        }
        )
        
        result.list[["results"]] = .createEnichmentTable(reactome_enrichment)
        
    }
    
    if (ontology == "DO"){
        
        .checkOntologyPackageInstallation("DOSE")
        
        DO_enrichment = tryCatch({
            
            DOSE::enrichDO(gene          = foreground_entrez,
                           universe      = background_entrez,
                           ont           = "DO",
                           pAdjustMethod = pAdjustMethod,
                           pvalueCutoff  = 1,
                           qvalueCutoff  = 1,
                           minGSSize     = minGSSize,
                           maxGSSize     = maxGSSize)
            
        }, error = function(e) {
            message = enrichmentErrorMessage
            .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
        }, warning = function(w) {
            # dont print anything
        }
        )
        
        result.list[["results"]] = .createEnichmentTable(DO_enrichment)
        
    }
    

    # Common parameter list
    result.list[["parameters"]] = list(
        "statistic"   = statistic,
        "background"  = backgroundStr,
        "nBackground" = nBackground,
        "nForeground" = nForeground)
    
    # Non-shared parameters
    if (ontology %in% c("KEGG", "DO", "Reactome")) {
        result.list[["parameters"]]$pAdjustMethod = pAdjustMethod
    } else {
        result.list[["parameters"]]$algorithm = algorithm
    }
    
    
    
    
    return(result.list)
}


.createEnichmentTable <- function (enrichmentObj) {
    
    if (!is.null(enrichmentObj)) {
        
        result.tbl = enrichmentObj@result %>%
            dplyr::rename(Term = Description, pval = .data$pvalue, Found = Count) %>%
            dplyr::mutate(GeneRatio = sapply(parse (text = enrichmentObj@result$GeneRatio), eval))
        
    } else {
        
        # Set an empty data frame so downstream aggregation functions dont stumble upon this
        result.tbl = tibble::tribble(~ID, ~Term, ~GeneRatio, ~BgRatio, ~pval, ~p.adjust, ~qvalue, ~geneID, ~Found)
    }
    
    result.tbl
    
}

# getEnrichmentResults <- function(GRN, enrichmentGroup, ontology, filePath = NULL){
#   
#   start = Sys.time()
#   GRN = .addFunctionLogToObject(GRN)
#   
#   checkmate::assertClass(GRN, "GRN")
#   checkmate::assertSubset(enrichmentType, c("general", "byCommunity", "byTF"))
#   checkmate::assertSubset(ontology, c("GO_BP", "GO_MF", "GO_CC"))
#   
#   if (enrichmentGroup == "general"){
#     
#   }
#   if (enrichmentGroup == "byCommunity")
#   
#   bind_rows(GRN@stats$Enrichment$general[c("GO_BP")], .id = "enrichmentGroup")
#   transpose(GRN@stats$Enrichment$byCommunity)
#   
# }


#' Generate graph communities and their summarizing statistics
#' 
#' This function generates the TF-gene graph from the filtered GRN object, and clusters its vertices into communities using established community detection algorithms.
#' @template GRN
#' @param clustering Character. Default \code{louvain}. One of: \code{louvain}, \code{leiden}, \code{leading_eigen}, \code{fast_greedy}, \code{optimal}, \code{walktrap}. The community detection algorithm to be used. Please bear in mind the robustness and time consumption of the algorithms when opting for an alternative to the default. 
#' @param ... Additional parameters for the used clustering method, see the \code{igraph::cluster_*} methods for details on the specific parameters and what they do. For \code{leiden} clustering, for example, you may add a \code{resolution_parameter} to control the granularity of the community detection or \code{n_iterations} to modify the number of iterations.
#' @template forceRerun
#' @return The same \code{\linkS4class{GRN}} object, with a table that consists of the connections clustered into communities stored in the \code{stats$communities} slot.
#' @import patchwork
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' GRN = calculateCommunitiesStats(GRN, forceRerun = FALSE)
#' @export
calculateCommunitiesStats <- function(GRN, clustering = "louvain", forceRerun = FALSE, ...){
    
    start = Sys.time()
    GRN = .addFunctionLogToObject(GRN)
    
    checkmate::assertClass(GRN, "GRN")
    checkmate::assertSubset(clustering, c("louvain", "leading_eigen", "fast_greedy", "optimal", "walktrap", "leiden"))
    checkmate::assertFlag(forceRerun)
    
    if (is.null(igraph::vertex.attributes(GRN@graph$TF_gene$graph)$community) | forceRerun) {
        
        futile.logger::flog.info(paste0("Calculating communities for clustering type ", clustering, "..."))
        
        if (clustering == "louvain") {
            
            communities_cluster = igraph::cluster_louvain(GRN@graph$TF_gene$graph, weights = NA, ...)
            
        } else if (clustering == "leading_eigen") {
            
            communities_cluster = igraph::cluster_leading_eigen(GRN@graph$TF_gene$graph, ...)
            
        } else if (clustering == "fast_greedy") {
            
            communities_cluster = igraph::cluster_fast_greedy(GRN@graph$TF_gene$graph, ...)
            
        } else if (clustering == "optimal") {
            
            communities_cluster = igraph::cluster_optimal(GRN@graph$TF_gene$graph, ...)
            
        } else if (clustering == "walktrap") {
            
            communities_cluster = igraph::cluster_walktrap(GRN@graph$TF_gene$graph, ...)
            
        } else if (clustering == "leiden") {
            
            # The default, see https://www.nature.com/articles/s41598-019-41695-z for a reasoning
            communities_cluster = igraph::cluster_leiden(GRN@graph$TF_gene$graph, ...)
   
        }
        
       
        
        # TODO: How redundant is it to store this separately?
        GRN@graph$TF_gene$clusterGraph = communities_cluster 
        
        # Add the community to the vertex metadata. First, sort them according to their size
        communities_count = sort(table(communities_cluster$membership), decreasing = TRUE)
        stopifnot(identical(igraph::vertex.attributes(GRN@graph$TF_gene$graph)$name, communities_cluster$names))
        igraph::vertex.attributes(GRN@graph$TF_gene$graph)$community = factor(communities_cluster$membership, levels = names(communities_count))
        
        nClustersMax = min(length(communities_count), 10)
        futile.logger::flog.info(paste0("Community summary for largest ", nClustersMax, " communities (Number of nodes per community, sorted by community size):"))
        for (clusterCur in seq_len(nClustersMax)) {
            futile.logger::flog.info(paste0(" Community ", names(communities_count)[clusterCur], ": ", communities_count[clusterCur], " nodes"))
        }
   

        
    } else {
        
        futile.logger::flog.info(paste0("Data already exists in object, nothing to do"))
    }  
    
    .printExecutionTime(start)
    
    GRN
}

#' Enrichment analysis for the genes in each community in the filtered \code{\linkS4class{GRN}}
#' 
#' After the vertices of the filtered GRN are clustered into communities using \code{\link{calculateCommunitiesStats}}, this function will run a per-community enrichment analysis.
#' @inheritParams calculateGeneralEnrichment
#' @param selection Character. Default "byRank". One of: "byRank", "byLabel". Specify whether the communities enrichment will by calculated based on their rank, where the largest community (with most vertices) would have a rank of 1, or by their label. Note that the label is independent of the rank.
#' @param communities Numeric vector. Default c(1:10). Depending on what was specified in the \code{display} parameter, this parameter would indicate either the rank or the label of the communities to be plotted. i.e. for \code{communities} = c(1,4), if \code{display} = "byRank" the GO enrichment for the first and fourth largest communities will be calculated if \code{display} = "byLabel", the results for the communities labeled "1", and "4" will be plotted.
#' @return The same \code{\linkS4class{GRN}} object, with the enrichment results stored in the \code{stats$Enrichment$byCommunity} slot.
#' @seealso \code{\link{plotCommunitiesEnrichment}}
#' @seealso \code{\link{plotGeneralEnrichment}}
#' @seealso \code{\link{calculateGeneralEnrichment}}
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' GRN = calculateCommunitiesEnrichment(GRN, ontology = c("GO_BP"), forceRerun = FALSE)
#' @export
calculateCommunitiesEnrichment <- function(GRN, 
                                           ontology = c("GO_BP", "GO_MF"), algorithm = "weight01", 
                                           statistic = "fisher", background = "neighborhood", 
                                           selection = "byRank", communities = seq_len(10),
                                           pAdjustMethod = "BH",
                                           forceRerun = FALSE) {
    
    start = Sys.time()
    GRN = .addFunctionLogToObject(GRN)
    
    checkmate::assertClass(GRN, "GRN")
    checkmate::assertSubset(ontology , c("GO_BP", "GO_MF", "GO_CC", "KEGG", "DO", "Reactome"))
    checkmate::assertSubset(algorithm , topGO::whichAlgorithms())
    checkmate::assertSubset(statistic , topGO::whichTests())
    checkmate::assertSubset(background, c("all_annotated", "all_RNA", "neighborhood"))
    checkmate::assertSubset(selection, c("byRank", "byLabel"))
    checkmate::assertSubset(pAdjustMethod, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"))
    checkmate::assertFlag(forceRerun)
    
    
    vertexMetadata = as.data.frame(igraph::vertex.attributes(GRN@graph$TF_gene$graph))
    
    foundCommunities = as.character(unique(vertexMetadata)$community)
    if (length(foundCommunities) == 0) {
        message = "No communities found, cannot calculate enrichment. Run the function calculateCommunitiesStats first. If you did already, it looks like no communities could be identified before"
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    if (selection == "byLabel"){
        communitiesDisplay = as.character(communities)
        # issue a warning if the community label does not exist
        diff.communities = setdiff(communitiesDisplay, foundCommunities)
        if (length(diff.communities)>0){
            message = paste("The following communities do not exist and will not be in the analysis: ", paste0(diff.communities, collapse = " + "))
            .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
            communitiesDisplay = setdiff(communitiesDisplay, diff.communities)
        }
    } else { # byRank
        
        communitiesDisplay = .selectCommunitesByRank(GRN, communities, graph = "TF_gene")
        
        

    }
    
    futile.logger::flog.info(paste0("Running enrichment analysis for selected ", length(communitiesDisplay), " communities. This may take a while..."))
    
    mapping = .getGenomeObject(GRN@config$parameters$genomeAssembly, type = "packageName")
    backgroundGenes = .getBackgroundGenes(GRN, type = background)
    
    if (is.null(GRN@stats[["Enrichment"]][["byCommunity"]])) {
        GRN@stats[["Enrichment"]][["byCommunity"]] = list()
    }
    
    for (communityCur in communitiesDisplay){
        
        futile.logger::flog.info(paste0(" Community ", communityCur))
        
        if (is.null(GRN@stats[["Enrichment"]][["byCommunity"]][[communityCur]])) {
            GRN@stats[["Enrichment"]][["byCommunity"]][[communityCur]] = list()
        }
        
       
        foregroundCur = vertexMetadata %>%
            dplyr::filter(community == communityCur) %>%
            dplyr::pull(name)
            
        for (ontologyCur in ontology) {
            
            if (is.null(GRN@stats[["Enrichment"]][["byCommunity"]][[communityCur]][[ontologyCur]]) | forceRerun) {
                
                GRN@stats[["Enrichment"]][["byCommunity"]][[communityCur]][[ontologyCur]] = 
                    .runEnrichment(foreground = foregroundCur, 
                                   background = backgroundGenes, 
                                   backgroundStr = background,
                                   ontology = ontologyCur, 
                                   algorithm = algorithm, 
                                   statistic = statistic,
                                   mapping = mapping,
                                   pAdjustMethod =  pAdjustMethod)
                
                futile.logger::flog.info(paste0("Result stored in GRN@stats$Enrichment$byCommunity[[\"", communityCur,  "\"]]$", ontologyCur, "$results"))
                
            } else {
                futile.logger::flog.info(paste0("Results already found / previously calculated. Not re-running as forceRerun = FALSE"))
            }
            
            
            
        }
    }
    
    .printExecutionTime(start)
    
    GRN
}

.selectCommunitesByRank <- function(GRN, communities, graph = "TF_gene") {
    
    df = igraph::vertex.attributes(GRN@graph[[graph]]$graph) %>%
        as.data.frame() %>% 
        dplyr::count(community)
    
    
    if (is.null(communities)) {
        communities = seq_len(nrow(df))
    }
    df %>% 
        dplyr::arrange(dplyr::desc(.data$n)) %>%
        dplyr::slice(communities) %>%
        dplyr::pull(community) %>%
        as.character()
}

#' Retrieve top Nodes in the filtered \code{\linkS4class{GRN}}
#' 
#' @template GRN
#' @param nodeType Character. One of: "gene", "TF". 
#' @param rankType Character. One of: "degree", "EV". This parameter will determine the criterion to be used to identify the "top" nodes. If set to "degree", the function will select top nodes based on the number of connections they have, i.e. based on their degree-centrality. If set to "EV" it will select the top nodes based on their eigenvector-centrality score in the network.
#' @param n Numeric. Default 0.1. If this parameter is passed as a value between (0,1), it is treated as a percentage of top nodes. If the value is passed as an integer it will be treated as the number of top nodes.
#' @param use_TF_gene_network \code{TRUE} or \code{FALSE}. Default \code{TRUE}. Should the TF-gene network be used (\code{TRUE}) or the TF-peak-gene network (\code{FALSE})?
#' @return A dataframe with the node names and the corresponding scores used to rank them
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' topGenes = getTopNodes(GRN, nodeType = "gene", rankType = "degree", n = 3)
#' topTFs = getTopNodes(GRN, nodeType = "TF", rankType = "EV", n = 5)
#' @export
getTopNodes <- function(GRN, nodeType, rankType, n = 0.1, use_TF_gene_network = TRUE) { # },
    #        TFConnectionType = "tf-gene", geneConnectionType = "peak-gene"){
    
    start = Sys.time()
    
    checkmate::assertClass(GRN, "GRN")
    checkmate::assertSubset(nodeType, c("gene", "TF"))
    checkmate::assertSubset(rankType, c("degree", "EV"))
    checkmate::assertFlag(use_TF_gene_network)
    #checkmate::assertSubset(TFConnectionType, c("tf-gene", "tf-peak"))
    #checkmate::assertSubset(geneConnectionType, c("peak-gene", "tf-gene"))
    checkmate::assertNumeric(n)
    
    if (nodeType == "gene") {
        slot = "gene.ENSEMBL"
        link = dplyr::if_else(use_TF_gene_network, "tf-gene", "peak-gene")
    } else {
        slot = "TF.name"
        slot = "TF.ENSEMBL"
        link = dplyr::if_else(use_TF_gene_network, "tf-gene", "tf-peak")
    } 
    
    graphType = dplyr::if_else(use_TF_gene_network, "TF_gene", "TF_peak_gene")
    
    
    if(n<1){
        # Get the total number of distinct nodes and calculate a percentage of that irrespective of ndoe degree
        top.n =  (GRN@connections$all.filtered$`0`[[slot]] %>% 
                      unique() %>% 
                      stats::na.omit() %>% 
                      length() * n) %>% round()
    }else{
        top.n = n
    }
    
    futile.logger::flog.info(paste0("n = ", n, " equals finding the top ", top.n, " ", rankType , "-central ", nodeType, "s in the network"))
    
    graph.df = GRN@graph[[graphType]]$table
    
    if (rankType == "degree"){
        col = dplyr::if_else(nodeType == "gene", "V2", "V1")
        topNodes = graph.df %>%
            dplyr::filter(connectionType == link) %>%
            dplyr::count(!!as.name(col), sort = TRUE) %>%
            # dplyr::rename(!!slot := V1, Connections = n) %>%
            dplyr::rename(Connections = n) %>%
            dplyr::arrange(dplyr::desc(Connections)) %>%
            dplyr::slice(seq_len(top.n)) 
        
        # TODO: change column names
        if (nodeType == "gene") {
            topNodes = topNodes  %>%
                dplyr::left_join(graph.df %>% dplyr::select(V2, V2_name) %>% dplyr::distinct(), by = "V2") %>%
                dplyr::rename(gene.ENSEMBL = V2, gene.name = V2_name)
        } else {
            topNodes = topNodes  %>%
                dplyr::left_join(graph.df %>% dplyr::select(V1, V1_name) %>% dplyr::distinct(), by = "V1") %>%
                dplyr::rename(TF.ENSEMBL = V1, TF.name = V1_name)
        }
        
        
        
    } else{ # if EV
        slot2 = dplyr::if_else(nodeType == "gene", "topGenes", "topTFs")
        topNodes = .getEigenCentralVertices(GRN, graphType = graphType, nCentralGenes = top.n, nCentralTFs = top.n)[[slot2]][["data"]]
    }
    return (topNodes)
}



#' Calculate TF-based GO enrichment
#' 
#' This function calculates the GO enrichment per TF, i.e. for the set of genes a given TF is connected to in the filtered \code{\linkS4class{GRN}}. 
#' 
#' @inheritParams calculateGeneralEnrichment
#' @param rankType Character. One of: "degree", "EV", "custom". This parameter will determine the criterion to be used to identify the "top" TFs. If set to "degree", the function will select top TFs based on the number of connections to genes they have, i.e. based on their degree-centrality. If set to "EV" it will select the top TFs based on their eigenvector-centrality score in the network. If set to custom, a set of TF names will have to be passed to the "TF.names" parameter.
#' @param n Numeric. Default 0.1. If this parameter is passed as a value between (0,1), it is treated as a percentage of top nodes. If the value is passed as an integer it will be treated as the number of top nodes. This parameter is not relevant if rankType = "custom".
#' @param TF.names Character vector. If the rank type is set to "custom", a vector of TF names for which the GO enrichment should be calculated should be passed to this parameter.
#' @return The same \code{\linkS4class{GRN}} object, with the enrichment results stored in the \code{stats$Enrichment$byTF} slot.
#' @examples 
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN =  loadExampleObject()
#' GRN =  calculateTFEnrichment(GRN, n = 5, ontology = "GO_BP", forceRerun = FALSE)
#' GRN =  calculateTFEnrichment(GRN, n = 5, ontology = "GO_BP", forceRerun = FALSE)
#' @export
calculateTFEnrichment <- function(GRN, rankType = "degree", n = 0.1, TF.names = NULL,
                                  ontology = c("GO_BP", "GO_MF"), algorithm = "weight01", 
                                  statistic = "fisher", background = "neighborhood",
                                  pAdjustMethod = "BH",
                                  forceRerun = FALSE){
    
    start = Sys.time()
    GRN = .addFunctionLogToObject(GRN)
    
    checkmate::assertClass(GRN, "GRN")
    checkmate::assertSubset(rankType, c("degree", "EV", "custom"))
    checkmate::assertNumeric(n)
    checkmate::assertSubset(ontology , c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "DO"))
    checkmate::assertSubset(algorithm , topGO::whichAlgorithms())
    checkmate::assertSubset(statistic , topGO::whichTests())
    checkmate::assertSubset(background, c("all_annotated", "all_RNA", "neighborhood"))
    checkmate::assertSubset(pAdjustMethod, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"))
    checkmate::assertFlag(forceRerun)
    
    if (rankType == "custom"){
        if(is.null(TF.names)){
            futile.logger::flog.error("To calculate the GO enrichment for a custom set of TFs, you must provide the TF names in the 'TF.names' parameter.")
        }
        wrongTFs = setdiff(TF.names, unique(GRN@connections$all.filtered$`0`$TF.name))
        if (length(wrongTFs)>0){
            futile.logger::flog.warn(paste0("The following TFs are not in the filtered GRN and will be ommited from the analysis: ",  paste0(wrongTFs, collapse = ", ")))
        }
        TFset = setdiff(TF.names, wrongTFs) 
    } else{
        
        # TF.name is always there, irrespective of whether ENSEMBL ID or TF name is used as primary ID type
        TFset = getTopNodes(GRN, nodeType = "TF", rankType = rankType, n, use_TF_gene_network = TRUE) %>% dplyr::pull(TF.name)
    }
    
    # TODO: Continue working on the TF.name level or switch to Ensembl? Should be in concordance with the graph!
    
    mapping = .getGenomeObject(GRN@config$parameters$genomeAssembly, type = "packageName")
    
    if (is.null(GRN@stats[["Enrichment"]][["byTF"]])) {
        
        GRN@stats[["Enrichment"]][["byTF"]] = list()
    }
    
    
    if (length(TFset) > 0) {
        futile.logger::flog.info(paste0("Running enrichment analysis for the following TFs: ", paste0(TFset, collapse = ", ")))
        
    } else {
        message = paste0("No TF fulfills the chosen criteria. Try increasing the value of the parameter n")
        .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
    }
    
    for (TF in as.character(TFset)){
        
        futile.logger::flog.info(paste0(" Running enrichment analysis for genes connected to the TF ", TF))
        
        # get the genes associated with current top TF
        curGenes = GRN@connections$all.filtered$`0` %>% 
            dplyr::filter(TF.name == TF) %>% 
            dplyr::pull(gene.ENSEMBL) %>%
            unique()
        
        
        backgroundGenes = .getBackgroundGenes(GRN, type = background)
        
        for (ontologyCur in ontology) {
            
            futile.logger::flog.info(paste0("  Ontology ", ontologyCur))
            
            if (is.null(GRN@stats[["Enrichment"]][["byTF"]][[TF]][[ontologyCur]]) | forceRerun) {
                
                if (is.null(GRN@stats[["Enrichment"]][["byTF"]][[TF]])) {
                    GRN@stats[["Enrichment"]][["byTF"]][[TF]] = list()
                }
                
                
                GRN@stats[["Enrichment"]][["byTF"]][[TF]][[ontologyCur]] =  
                    .runEnrichment(foreground = curGenes,
                                   background = backgroundGenes, 
                                   backgroundStr = background,
                                   ontology = ontologyCur, 
                                   algorithm = algorithm, 
                                   statistic = statistic,
                                   mapping = mapping,
                                   pAdjustMethod =  pAdjustMethod)
                
                futile.logger::flog.info(paste0("   Results stored in GRN@stats$Enrichment$byTF[[\"", TF, "\"]]$", ontologyCur, "$results"))
                
            }else {
                futile.logger::flog.info(paste0("   Results already found / previously calculated for TF ", TF, ". Not re-running as forceRerun = FALSE"))
            }
        }
        
    }
    
    .printExecutionTime(start)
    GRN
}
