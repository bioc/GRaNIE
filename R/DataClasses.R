
#### S4 class definition
#' Create, represent, investigate, quantify and visualize enhancer-mediated gene regulatory networks (\strong{eGRNs})
#' 
#' The class \code{\linkS4class{GRN}} stores data and information related to our \code{eGRN} approach to construct enhancer-mediated gene regulatory networks out of open chromatin and RNA-Seq data. See the description below for more details, and \strong{visit our project website at https://grp-zaugg.embl-community.io/GRaNIE and have a look at the various Vignettes}.
#' @section Constructors:
#' Currently, a \code{\linkS4class{GRN}} object is created by executing the function \code{\link{initializeGRN}}.
#' 
#' @section Accessors:
#' In the following code snippets, \code{GRN} is a \code{\linkS4class{GRN}} object.
#' 
#'# Get general annotation of a GRaNIE object
#' 
#' \code{nPeaks(GRN))} and  \code{nGenes(GRN))}: Retrieve the number of peaks and genes, respectively, that have been added to the object (both before and after filtering)
#'
#' @slot data Currently stores 4 different types of data:\cr 
#' \itemize{
#' \item \code{peaks}:
#'  \itemize{
#'  \item \code{counts_norm}:
#'  \item \code{counts_orig}:
#'  \item \code{consensusPeaks}:
#'  }
#' \item \code{RNA}: 
#'  \itemize{
#'  \item \code{counts_norm.l}:
#'  \item \code{counts_orig}:
#'  }
#' \item \code{metadata}: 
#' \item \code{TFs}: 
#'  \itemize{
#'  \item \code{translationTable}:
#'  \item \code{TF_peak_overlap}:
#'  \item \code{classification}:
#'  }
#' }

#' @slot config Contains general configuration data and parameters such as parameters, files, directories, flags, and recorded function parameters.
#' 
#' @slot connections  Stores various types of connections
#' 
#' @slot annotation Stores annotation data for peaks and genes
#' 
#' @slot stats Stores statistical and summary information for a \code{GRN} network. Currently, connection details are stored here.
#' 
#' @slot visualization Stores visualization results, currently always empty. Feature in development.
#' @slot isDev Flag whether this is an object from the development version of the package

#' @keywords GRN-class, GRN
#' @name GRN-class
#' @aliases GRN-class
#
setClass("GRN", representation = representation(
  annotation    = "list",
  config        = "list", 
  connections   = "list",
  data          = "list",
  stats         = "list",
  visualization = "list",
  graph         = "list",
  isDev         = "logical"
  )
)

# setMethod("initialize", "GRN", function(.Object, x) {
#   # class(.Object) ="GRN"
#   .Object
# })

#' @import checkmate
#' @import GenomicRanges
#' @importFrom GenomicRanges mcols
.validGRNObj <- function(object) {
  
  valid = TRUE
  msg   = NULL
  
  # 
  # if (object@internal$countType == "enrichment" & !parameters(object)$normByInput) {
  #     valid = FALSE
  #     msg = c(msg, "The element countType in the slot \"internal\" must match with the value of the parameter normByInput.\n")
  # }
  # 
  if (valid) TRUE else msg
  
}

setValidity("GRN", .validGRNObj)        


#' @importFrom methods new 
.createGRNObject <- function() {
  
  
  new(methods::getClassDef("GRN", package = utils::packageName(), inherits = FALSE))
}

# # OTHER FUNCTIONS #

#' @importFrom utils tail
setMethod("show",
          "GRN",
          function(object) {
            
            GRN = object
            # .checkObjectValidity(GRN)
            
            packageName = utils::packageName()
            packageVersion = ifelse(is.null(packageName), NA, utils::packageVersion(packageName))
            
            cat("Object of class:", packageName," ( version", paste0(packageVersion[[1]], collapse = "."), ")\n")
            
            
            
            
            #nGenes       = nGenes(GRN)
            
            cat("Data summary:\n")
            if (!is.null(GRN@data$peaks$consensusPeaks)) {
              nPeaks_filt  = nPeaks(GRN, filter = TRUE)
              nPeaks_all   = nPeaks(GRN, filter = FALSE)
              cat(" Number of peaks (filtered, all): ", nPeaks_filt, ", ",  nPeaks_all, "\n", sep = "")
            } else {
              cat (" Number of peaks: No peak data found.\n")
            }
            
            if (!is.null(GRN@data$RNA$counts_orig)) {
               
              nGenes_filt = nGenes(GRN, filter = TRUE)
              nGenes_all  = nGenes(GRN, filter = FALSE)
              
              if (checkmate::testClass(GRN@data$RNA$counts_orig, "DESeqDataSet")) {
                
              } else {
                
              }
             
              cat(" Number of genes (filtered, all): ", nGenes_filt, ", ",  nGenes_all, "\n", sep = "")
            } else {
              cat (" Number of genes: No RNA-seq data found.\n")
            }
            
            cat("Parameters:\n")
            if (!is.null(GRN@config$metadata$genomeAsembly)) {
              cat(" Genome assembly: ", GRN@config$parameters$genomeAssembly, "\n")
            }
            
            cat("Provided metadata:\n")
            
            for (nameCur in names(GRN@config$metadata)) {
              cat(" ", nameCur, ": ", unlist(GRN@config$metadata[nameCur]), "\n")
 
            }
            
            
            cat("Connections:\n")
            if (!is.null(GRN@connections)) {

                if (!is.null(GRN@connections$TF_peaks$`0`$main)) {
                    maxFDR = round(max(GRN@connections$TF_peaks$`0`$main$TF_peak.fdr), 1)
                    cat(" TF-peak links (", nrow(GRN@connections$TF_peaks$`0`$main), " with FDR < ", maxFDR, ")\n", sep = "")
                } else {
                    cat(" TF-peak links: none found\n")
                }
                if (!is.null(GRN@connections$peak_genes$`0`)) {
                    cat(" peak-gene links (", nrow(GRN@connections$peak_genes$`0`)," with ", GRN@config$parameters$promoterRange, " bp promoter range)\n", sep = "")
                } else {
                    cat(" peak-gene links: none found\n")
                }
                
                
                if (!is.null(GRN@connections$all.filtered$`0`)) {
                    nRows = nrow(GRN@connections$all.filtered$`0`)
                    max_TF_peak_FDR = round(max(GRN@connections$all.filtered$`0`$TF_peak.fdr), 1)
                    max_peak_gene_FDR = round(max(GRN@connections$all.filtered$`0`$peak_gene.p_adj), 1)
                    cat(" TF-peak-gene links (", nRows," with TF-peak FDR ", max_TF_peak_FDR, " and peak-gene FDR ",  max_peak_gene_FDR , ")\n", sep = "")
                    
                } else {
                    cat(" TF-peak-gene links (filtered): none found\n")
                }
            

              
            } else {
              cat (" Connections: none found")
            }
            
            cat("Network-related:\n")
            
            if (!is.null(GRN@graph$TF_gene)) {
                
                cat(" eGRN network:\n")
                nNodes = length(igraph::V(GRN@graph$TF_gene$graph))
                nEdges = length(igraph::E(GRN@graph$TF_gene$graph))
                cat("   TF-gene eGRN: ", nNodes, " nodes and ", nEdges, " edges\n",sep = "")
 
                nNodes = length(igraph::V(GRN@graph$TF_peak_gene$graph))
                nEdges = length(igraph::E(GRN@graph$TF_peak_gene$graph))
                cat("   TF-peak-gene eGRN: ", nNodes, " nodes and ", nEdges, " edges\n",sep = "")              
                
                # Community identification (no, yes and how many and how many nodes each)
                cat(" Communities (TF-gene):\n")
                df = igraph::vertex.attributes(GRN@graph[["TF_gene"]]$graph)
                if (!is.null(df) & "community" %in% names(df)) {
                    communities = df %>% as.data.frame() %>% dplyr::count(.data$community) %>% dplyr::arrange(dplyr::desc(.data$n))
                    cat("  Communities, sorted by size (n = Number of nodes): ", paste0(communities$community, " (n=", communities$n, collapse = "), "), ")\n", sep = "")
                } else {
                    cat("  None found\n")
                }
                
                cat(" Enrichment (TF-gene):\n")
                if (!is.null(GRN@stats$Enrichment$byTF)) {
                    ontologies = names(GRN@stats$Enrichment$general)
                    cat("  Overall network: ", paste0(ontologies, collapse = " & "), "\n",sep = "")
                } else {
                    cat("  Overall network: no enrichment results found\n")
                }
                
                if (!is.null(GRN@stats$Enrichment$byTF)) {
                    ontologies = names(GRN@stats$Enrichment$byCommunity[[1]])
                    cat("  Communities: ", paste0(ontologies, collapse = " & "), "\n",sep = "")
                } else {
                    cat("  Communities: no enrichment results found\n")
                }
                
                if (!is.null(GRN@stats$Enrichment$byTF)) {
                    ontologies = names(GRN@stats$Enrichment$byTF[[1]])
                    cat("  TFs: ", paste0(ontologies, collapse = " & "), "\n",sep = "")
                } else {
                    cat("  TFs: no enrichment results found\n")
                }
                
            } else {
                cat("  eGRN network: not found\n")
            }
            
          }
)



#' Get the number of peaks for a \code{\linkS4class{GRN}} object.
#' 
#' Returns the number of peaks (all or only non-filtered ones) from the provided peak datain the \code{\linkS4class{GRN}} object.
#'
#' @template GRN
#' @param filter TRUE or FALSE. Default TRUE. Should peaks marked as filtered be included in the count?
#' @return Integer. Number of peaks that are defined in the \code{\linkS4class{GRN}} object, either by excluding (filter = TRUE) or including (filter = FALSE) peaks that are currently marked as \emph{filtered}.
#' @seealso \code{\link{nTFs}}
#' @seealso \code{\link{nGenes}}
#' @examples
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' nPeaks(GRN, filter = TRUE)
#' nPeaks(GRN, filter = FALSE)
#' @export
#' @aliases peaks
nPeaks <- function(GRN, filter = TRUE) {
  
  checkmate::assertClass(GRN, "GRN")
  checkmate::assertFlag(filter)
  
  if (is.null(GRN@data$peaks$consensusPeaks)) {
    return(NA)
  }
  
  if (filter) {
    nPeaks = GRN@data$peaks$consensusPeaks %>% dplyr::filter(!.data$isFiltered) %>% nrow()
  } else {
    nPeaks = GRN@data$peaks$consensusPeaks %>% nrow()
  }
  
  nPeaks
  
}


#' Get the number of genes for a \code{\linkS4class{GRN}} object.
#' 
#' Returns the number of genes (all or only non-filtered ones) from the provided RNA-seq data in the \code{\linkS4class{GRN}} object.
#'
#' @template GRN
#' @param filter TRUE or FALSE. Default TRUE. Should genes marked as filtered be included in the count?
#' @return Integer. Number of genes that are defined in the \code{\linkS4class{GRN}} object, either by excluding (filter = TRUE) or including (filter = FALSE) genes that are currently marked as \emph{filtered}.
#' @seealso \code{\link{nTFs}}
#' @seealso \code{\link{nPeaks}}
#' @examples
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' nGenes(GRN, filter = TRUE)
#' nGenes(GRN, filter = FALSE)
#' @export
#' @aliases genes
nGenes <- function(GRN, filter = TRUE) {
  
  checkmate::assertClass(GRN, "GRN")
  checkmate::assertFlag(filter)
  
  if (is.null(GRN@data$RNA$counts_norm.l)) {
    return(NA)
  }
  
  if (filter) {
    nGenes = GRN@data$RNA$counts_norm.l[["0"]] %>% dplyr::filter(!.data$isFiltered) %>% nrow()
  } else {
    nGenes = GRN@data$RNA$counts_norm.l[["0"]] %>% nrow()
  }
  
  nGenes
  
}



#' Get the number of TFs for a \code{\linkS4class{GRN}} object.
#' 
#' Returns the number of TFs  from the provided TFBS data in the \code{\linkS4class{GRN}} object.
#'
#' @template GRN
#' @return Integer. Number of TFs that are defined in the \code{\linkS4class{GRN}} object.
#' @seealso \code{\link{nGenes}}
#' @seealso \code{\link{nPeaks}}
#' @examples
#' # See the Workflow vignette on the GRaNIE website for examples
#' GRN = loadExampleObject()
#' nTFs(GRN)
#' @export
#' @aliases TFs
nTFs <- function(GRN) {
    
    checkmate::assertClass(GRN, "GRN")
    
    if (is.null(GRN@data$TFs$translationTable)) {
        return(NA)
    }
    
    nrow(GRN@data$TFs$translationTable)
    
}