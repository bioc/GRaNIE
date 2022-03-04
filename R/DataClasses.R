
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
#'# Get general annotation of a GRN object
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
  
  #In the dev version, we add an additional slot to differentiate this from a GRN object so that the class cache does not load the GRN package instead automatically
  # GRN  <-  new("GRN", 
  #              annotation    = list(),
  #              config        = list(),
  #              connections   = list(),
  #              data          = list(),
  #              stats         = list(),
  #              visualization = list(),
  #              isDev = TRUE
  # )
  
  new(methods::getClassDef("GRN", package = utils::packageName(), inherits = FALSE))
}

# # 
# # 
# # # SLOT config #
# # 
# #' Retrieve the parameters of an object.
# #'
# #' @param object An object containing parameters with which it was created.
# #' @param ... Additional arguments, for use in specific methods.
# #' @docType methods
# #' @rdname parameters-methods
# #' @aliases parameters
# #' @export
# setGeneric("parameters", function(object, ...) standardGeneric("parameters"))
# 
# #' Retrieve the parameters of a \code{SNPhood} object.
# #'
# #' @docType methods
# #' @rdname parameters-methods
# #' @aliases parameters
# #' @return A named list with all parameters and its current values of the \code{GRN} object.
# #' @examples
# # #' data(GRN, package="SNPhood")
# # #' parameters(GRN)
# #' @export
# setMethod("parameters", "GRN", function(object, ...) {object@config})

# # # SLOT config #
# # 
# #' PCA plot for an object
# #'
# #' @param object An object containing parameters with which it was created.
# #' @param ... Additional arguments, for use in specific methods.
# #' @docType methods
# #' @rdname plotPCA-methods
# #' @aliases plotPCA
# #' @export
# setGeneric("plotPCA", function(object, ...) standardGeneric("plotPCA"))


# #' PCA plots for ATAC and RNA data from a \code{GRN} object.
# #'
# #' @export
# setMethod("plotPCA", "GRN", function(object, ...) {.plotPCA_generic(object, ...)})
# 
# .plotPCA_generic <- function (GRN, outputFolder, type = c("RNA", "ATAC"), forceRerun = FALSE) {
# 
#   plotPCA(GRN, outputFolder, type, forceRerun) 
# }

# setMethod("enrichment", "SNPhood", function(object, readGroup = NULL, dataset = NULL, ...) {.getEnrichment(object, readGroup, dataset)})
# 
# .getEnrichment <- function(GRN, readGroup = NULL, dataset = NULL) {
# 
#   .getCounts(GRN, type = "enrichmentBinned", readGroup, dataset)
# }

# 
# #' Extract count data from a \code{\link{SNPhood}} object.
# #' 
# #' \code{counts} extracts count data from a \code{\link{SNPhood}} object. The full count data or only a subset can be extracted 
# #' by settings the parameters \code{type}, \code{readGroup} and \code{dataset} accordingly. Either the count data
# #' for the unbinned or binned SNP regions can be extracted.
# #' @template object
# #' @param type Character(1). Default "binned". Either "binned" or "unbinned" to extract counts after or before binning the SNP regions, respectively.  
# #' @template readGroup
# #' @template dataset
# #' @param ... not used
# #' @return A named nested list with the requested count data, organized after read group and dataset.
# #' @aliases counts
# #' @docType methods
# #' @rdname counts-method
# #' @examples
# #' data(GRN, package="SNPhood")
# #' str(counts(GRN))
# #' str(counts(GRN, readGroup = "paternal", dataset = 1))
# #' str(counts(GRN, readGroup = c("maternal", "paternal"), dataset = 1))
# #' @seealso \code{\link{SNPhood}}, \code{\link{enrichment}}
# #' @export
# #' @importFrom BiocGenerics counts
# setMethod("counts", "SNPhood", function(object, type = "binned", readGroup = NULL, dataset = NULL, ...) {.getCounts(object, type, readGroup, dataset)})
# 
# #' @import checkmate
# .getCounts <- function(GRN, type, readGroup = NULL, dataset = NULL, ...) {
# 
#   .checkObjectValidity(GRN)
# 
#   # TODO: extend and enable datasets and not just one and readGroups and not just one
# 
#   assertChoice(type, c("binned", "unbinned", "enrichmentBinned"))
#   assertSubset(readGroup, GRN@annotation$readGroups)
#   assert(checkNull(dataset), 
#          checkIntegerish(dataset, lower = 1, upper = length(GRN@readCountsUnbinned[[1]]), any.missing = FALSE, min.len = 1, unique = TRUE), 
#          checkSubset(dataset, names(GRN@readCountsUnbinned[[1]])))
# 
# 
# 
#   # get the read groups the user does not want to have
#   discardReadGroups = NULL
#   if (!testNull(readGroup)) {
#     discardReadGroups = setdiff(GRN@annotation$readGroups, readGroup)       
#   }
# 
#   discardDatasets = NULL
#   if (!testNull(dataset)) {
# 
#     if (testSubset(dataset, names(GRN@readCountsUnbinned[[1]]))) {
#       discardDatasets = setdiff(names(GRN@readCountsUnbinned[[1]]), dataset)       
# 
#     } else {
#       discardDatasets = setdiff(seq_len(length(GRN@readCountsUnbinned[[1]])), dataset)
#     }
#   }
# 
#   if (type == "unbinned") {
#     counts.l = GRN@readCountsUnbinned
#   } else if (type == "binned") {
#     counts.l = GRN@readCountsBinned
#   } else {
#     counts.l = GRN@enrichmentBinned
#   }
# 
# 
# 
#   if (!(testNull(readGroup) & testNull(dataset))) {
# 
#     # Filter read groups
#     for (readGroupCur in discardReadGroups) {
#       counts.l[[readGroupCur]] = NULL
#     }
# 
#     # Filter datasets
#     for (readGroupCur in names(counts.l)) {
# 
#       # Delete datasets, but presort so the indexes are still correct
#       if (!is.null(discardDatasets)) {
#         for (datasetCur in sort(discardDatasets, decreasing = TRUE)) {
#           counts.l[[readGroupCur]][[datasetCur]] = NULL
#         }
#       }
# 
# 
#     } 
# 
#   }
# 
# 
#   for (readGroupCur in names(counts.l)) {
# 
#     for (datasetCur in names(counts.l[[readGroupCur]])) {
# 
#       if (type == "binned" | type == "enrichmentBinned") {
#         rownames(counts.l[[readGroupCur]] [[datasetCur]]) = annotationRegions(GRN)
#         colnames(counts.l[[readGroupCur]] [[datasetCur]]) = GRN@annotation$bins
#       } else {
#         names(counts.l[[readGroupCur]] [[datasetCur]]) = annotationRegions(GRN)
# 
#       }
# 
#     }
# 
#   } 
# 
#   # Reduce the complexity of the reuslting count object
#   # First for the read group
#   if (length(counts.l) == 1) {
#     counts.l = counts.l[[1]]
# 
#     # Repeat it for the dataset
#     if (length(counts.l) == 1) {
#       counts.l = counts.l[[1]]
#     }
# 
#   } else { # If multiple read groups are defined, but only one dataset
# 
#     for (i in seq_len(length(counts.l))) {
#       if (length(counts.l[[i]]) == 1) {
#         counts.l[[i]] = counts.l[[i]][[1]]
#       }
#     }
# 
#   }
# 
#   if (length(counts.l) == 0) {
# 
#     if (GRN@config$onlyPrepareForDatasetCorrelation) {
#       warning("Returning an empty list, could not find the requested data in the object. The parameter \"onlyPrepareForDatasetCorrelation\" has been set to TRUE. Rrun the function analyzeSNPhood again and set the parameter to FALSE. See also the help pages.")
# 
#     } else {
#       warning("Returning an empty list, could not find the requested data in the object. Did you ask for the correct type of data? See also the help pages.")
# 
#     }
#   }
# 
#   return(counts.l)
# 
# }
# 
# 
# #' Extract enrichment data from an object.
# #' 
# #' @param object An object containing enrichment information.
# #' @param ... Additional arguments, for use in specific methods.
# #' @return Enrichment of the object or the objects components.
# #' @export
# #' @docType methods
# #' @rdname enrichment-methods
# #' @aliases enrichment
# setGeneric("enrichment", function(object, ...) standardGeneric("enrichment")) 
# 
# #' Extract enrichment data from a \code{\link{SNPhood}} object.
# #' 
# #' \code{enrichment} extracts enrichment data from a \code{\link{SNPhood}} object. The full count data or only a subset can be extracted 
# #' by settings the parameters \code{type}, \code{readGroup} and \code{dataset} accordingly. Principally, either the count data
# #' for the unbinned or binned SNP regions can be extracted.
# #' @template readGroup
# #' @template dataset
# #' @return Named list with the requested enrichment matrices from the \code{\link{SNPhood}} object, organized by read group and dataset
# #' @seealso \code{\link{counts}}
# #' 
# #' @export
# #' @docType methods
# #' @rdname enrichment-methods
# #' @aliases enrichment
# #' @examples
# #' data(GRN, package="SNPhood")
# #' str(enrichment(GRN), list.len=5)
# setMethod("enrichment", "SNPhood", function(object, readGroup = NULL, dataset = NULL, ...) {.getEnrichment(object, readGroup, dataset)})
# 
# .getEnrichment <- function(GRN, readGroup = NULL, dataset = NULL) {
# 
#   .getCounts(GRN, type = "enrichmentBinned", readGroup, dataset)
# }
# 
# 
# 
# # SLOT ANNOTATION #
# 
# #' Retrieve the annotation of a \code{SNPhood} object.
# #' 
# #' Specific elements within the annotation slot may also be extracted by using the \code{elements} parameter.
# #' @template object
# #' @param elements Character. The name of the element(s) in the annotation slot to be extracted. 
# #' If set to \code{NULL}, the full annotation slot is returned.
# #' @param ... not supported
# #' @export
# #' @return If only a single value for \code{elements} is provided, the element is returned directly. 
# #' If multiple values are provided, a named list with the requested elements is returned.
# #' @docType methods
# #' @rdname annotation-methods
# #' @aliases annotation
# #' @importFrom BiocGenerics annotation
# #' @examples 
# #' data(GRN, package="SNPhood")
# #' annotation(GRN)
# #' annotation(GRN, elements = "regions")
# #' annotation(GRN, elements = c("regions", "bins"))
# setMethod("annotation", "SNPhood", function(object, elements = NULL, ...) {.getAnnotation(object, elements)})
# 
# .getAnnotation <- function(GRN, elements) {
# 
#   .checkObjectValidity(GRN)
# 
#   result = NULL
# 
#   assertSubset(elements, names(GRN@annotation))
# 
#   if (!testNull(elements)) {
# 
#     result = GRN@annotation[elements]
# 
#     if (length(elements) == 1) {
#       result = result[[elements]]
#     } 
# 
#   } else {
# 
#     result = GRN@annotation
#   }
# 
# 
#   result
# 
# }
# 
# 
# 
# #' Get results of various analyses performed with a \code{SNPhood} object.
# #' 
# #' Return the results of a particular analyis that is stored in the \code{SNPhood} object.
# #'
# #' @template SNPhood
# #' @param type Character(1). Name of analyses one wants to retrieve the results for. 
# #' Currently supported are "allelicBias", "clustering", "genotype" and "samplesCorrelation".
# #' @param elements Character. Default NULL. Which elements of the resulting list structure should be returned? 
# #' If set to NULL, all elements will be returned. Otherwise, if names are provided, only the requested subset elements will be returned.
# #' If type equals "allelicBias", valid values are "pValue", "confIntervalMin", "confIntervalMax", "fractionEstimate", "background", "FDR_results", and "parameters".
# #' If type equals "clustering", valid values are the defined read groups in the object.
# #' If type equals "genotype", valid values are "strongGenotypes", "weakGenotypes", and "invariantGenotypes".
# #' If type equals "samplesCorrelation", valid values are "corTable", and "transl".
# 
# #' @return A list with the results of the requested analysis and elements within.
# #' @examples
# #' data(GRN, package="SNPhood")
# #' head(results(GRN, type="allelicBias", elements = "parameters"))
# #' head(results(GRN, type="allelicBias"))
# #' @export
# #' @import checkmate
# #' @importFrom utils head
# results <- function(GRN, type, elements = NULL) {
# 
#   # Check types and validity of arguments 
#   .checkObjectValidity(GRN)
#   assertChoice(type, GRN@internal$addResultsElementsAdded)
# 
#   assert(checkNull(elements), 
#          checkIntegerish(elements, any.missing = FALSE, lower = 1, upper = length(GRN@additionalResults[[type]])), 
#          checkSubset(elements, names(GRN@additionalResults[[type]]))
#   )
# 
#   if (type == "allelicBias") {
# 
#     if (testNull(GRN@additionalResults$allelicBias) | testNull(names(GRN@additionalResults$allelicBias))) {
#       stop("Could not find results of allelic bias test, did you call testForAllelicBiases before?")
#     }
# 
#     results = GRN@additionalResults$allelicBias
# 
#     # Add row and coumn names for each matrix. Only do it now to save space
#     elemsToAnnotate = c("pValue", "confIntervalMin", "confIntervalMax", "fractionEstimate")
# 
# 
#     for (elemCur in elemsToAnnotate) {
#       for (datasetCur in names(results[[elemCur]])) {
#         dimnames(results[[elemCur]][[datasetCur]]) = list(annotationRegions(GRN), GRN@annotation$bins)
#       }
# 
#     }
# 
#     if (results$parameters$calculateBackgroundDistribution) {
#       for (datasetCur in names(results[["background"]])) {
#         dimnames(results[[elemCur]][[datasetCur]]) = list(annotationRegions(GRN), GRN@annotation$bins)
#       }
#     }
# 
# 
#     if (!testNull(elements)) {
#       for (elementCur in names(results)) {
#         if (!elementCur %in% elements) results[[elementCur]] = NULL
#       }
# 
#       results = results[[elements]]
#     }
# 
# 
#     return(results)
# 
#   } else if (type == "genotype") {
# 
#     if (testNull(GRN@additionalResults[[type]]) | testNull(names(GRN@additionalResults[[type]]))) {
#       stop("Could not find results of genotype analysis. Did you run the function calculateWeakAndStrongGenotype before?")
#     }
# 
# 
#   } else if (type == "samplesCorrelation") {
# 
#     if (testNull(GRN@additionalResults[[type]]) | testNull(names(GRN@additionalResults[[type]]))) {
#       stop("Could not find results of samplesCorrelation analysis. Did you run the function calcCorrelationSamples before?")
#     }
# 
# 
# 
#   } else if (type == "clustering") {
# 
#     if (testNull(GRN@additionalResults[[type]]) | testNull(names(GRN@additionalResults[[type]]))) {
#       stop("Could not find results of samplesCorrelation analysis. Did you run the function clusterCountMatrix before?")
#     }
# 
#   }
# 
#   if (length(elements) == 0) {
#     results = GRN@additionalResults[[type]]
#   } else if (length(elements) == 1) {
#     results = GRN@additionalResults[[type]][[elements]]
#   } else {
#     results = GRN@additionalResults[[type]][elements]
#   }
# 
#   return(results)
# 
# }
# 
# 
# 
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
                if (!is.null(df)) {
                    communities = df %>% as.data.frame() %>% dplyr::count(community) %>% dplyr::arrange(desc(n))
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


# # SNPhood OBJECT VALIDITY #
# 
# 
# #' @import checkmate
# #' @importFrom methods validObject
# .checkObjectValidity <- function(object, verbose = FALSE) {
# 
#   assertClass(object, "SNPhood")
# 
#   if (!is.null(object@internal$disableObjectIntegrityChecking)) {
#     if (object@internal$disableObjectIntegrityChecking == FALSE) {
#       # if (verbose) message("Check object integrity and validity. For large objects, this may take some time. Use the function changeObjectIntegrityChecking to disable this check for the object.")
# 
#     }
#   }
# 
#   res = validObject(object, test = TRUE, complete = TRUE)
#   if (testFlag(res)) {
#     if (!res) {
#       stop("The SNPhood object is not valid for the following reason(s): \n\n", res)
#     }
#   } else {
#     stop("The SNPhood object is not valid for the following reason(s): \n\n", res)
#   }   
# 
# }
# 
# #' @import checkmate
# #' @import GenomicRanges
# # @importFrom GenomicRanges mcols
# .validSNPhoodObj <- function(object) {
# 
#   valid = TRUE
#   msg   = NULL
# 
#   #saveRDS(object,file = "/home/christian/SNPhood.rds")
# 
#   # Skip validity check if object has not yet been fully constructed
#   if (object@internal$disableObjectIntegrityChecking) return(TRUE) 
# 
# 
#   if (object@config$onlyPrepareForDatasetCorrelation) return(TRUE) 
# 
# 
#   # SLOT ANNOTATION #
# 
# 
#   nRows = length(object@annotation$regions)
#   nCols = length(object@annotation$bins)
#   nDatasets = length(object@annotation$files)
# 
#   if (testNull(object@annotation) | testNull(names(object@annotation))) {
# 
#     valid = FALSE
#     msg = c(msg, "Slot \"annotation\" or its names must not be NULL.\n")
# 
#   } else {
# 
#     validNames = c("regions", "genotype", "files", "readGroups", "bins")
#     if (!testSubset(validNames, names(object@annotation), empty.ok = FALSE)) {
#       valid = FALSE
#       msg = c(msg, "Slot \"annotation\" must contain the elements ",paste0(validNames, collapse = ","),".\n")
#     }
# 
#     ## Read groups ##
#     if (!testCharacter(object@annotation$readGroups, min.chars = 1, any.missing = FALSE, min.len = 1)) {
#       valid = FALSE
#       msg = c(msg, "Element \"readGroups\" in slot \"annotation\" must be a character vector with at least one element.\n")
#     }  
# 
#     ## Files ##
#     if (!testList(object@annotation$files, any.missing = FALSE, min.len = 1, types = c("list"))) {
#       valid = FALSE
#       msg = c(msg, "Element \"files\" in slot \"annotation\" must be a named list of lists with at least one element.\n")
#     } 
# 
#     for (i in seq_len(length(object@annotation$files))) {
#       if (!testList(object@annotation$files[[i]], any.missing = FALSE, min.len = 4, types = c("character", "logical", "list"))) {
#         valid = FALSE
#         msg = c(msg, "Each element in \"files\" in slot \"annotation\" must be a list with four elements of type \"character\", \"logical\" or \"list\".\n")
#       } 
#     }
# 
#     ## Bins ##
#     if (!object@config$onlyPrepareForDatasetCorrelation) {
# 
# 
#       if (!object@config$normByInput) {
#         lenMatrix = ncol(object@readCountsBinned[[1]][[1]])
#       } else {
#         lenMatrix = ncol(object@enrichmentBinned[[1]][[1]])
#       }
# 
#       if (!testCharacter(object@annotation$bins, min.chars = 1, any.missing = FALSE, len = lenMatrix)) {
#         valid = FALSE
#         msg = c(msg, "The elements bins in the slot \"annotation\" must be of type character with a minimal length of 1 and no missing elements.\n")
#       }  
# 
#     }
# 
# 
#     ## Genotype ##
#     validNames = c("readsDerived", "external", "heterozygosity", "mismatches")
#     if (!testSubset(names(object@annotation$genotype), validNames)) {
# 
#       valid = FALSE
#       msg = c(msg, "Element \"genotype\" in slot annotation must must only contain the elements ",paste0(validNames, collapse = ","),".\n")
# 
#     }
# 
#     # 1. readsDerived
#     for (i in seq_len(length(object@annotation$readGroups))) {
# 
#       if (!testSubset(names(object@annotation$genotype$readsDerived)[i], object@annotation$readGroups[i])) {
#         valid = FALSE
#         msg = c(msg, "The names of the read groups are incorrect for the element readsDerived in the slot \"annotation$genotype\".\n")
#       }
# 
#       for (j in seq_len(length(object@annotation$files))) {
# 
#         if (!testSubset(names(object@annotation$genotype$readsDerived[[i]])[j], names(object@annotation$files)[j])) {
#           valid = FALSE
#           msg = c(msg, "The names of the datasets are incorrect for the element readsDerived in the slot \"annotation$genotype\".\n")
#         }
# 
#         if (!testMatrix(object@annotation$genotype$readsDerived[[i]][[j]], any.missing = FALSE, nrows = 4, ncols = length(object@annotation$regions))) {
#           valid = FALSE
#           msg = c(msg, "The genotype distribution element is invalid in the slot \"annotation$genotype$readsDerived\".\n")
#         }
# 
#         if (!testSubset(rownames(object@annotation$genotype$readsDerived[[i]][[j]]), c("A","C","G","T"))) {
#           valid = FALSE
#           msg = c(msg, "The names of the datasets are incorrect for the element readsDerived in the slot \"annotation$genotype\".\n")
#         }
# 
#         if (!testIntegerish(object@annotation$genotype$readsDerived[[i]][[j]], lower = 0, any.missing = FALSE)) {
#           valid = FALSE
#           msg = c(msg, "The genotype distribution element is invalid in the slot \"annotation$genotype$readsDerived\".\n")
#         }
# 
#       }
#     }
# 
#     # 2. External
# 
#     if (!testNull(object@annotation$genotype$external)) {
# 
#       if (!testDataFrame(object@annotation$genotype$external, min.cols = 4, all.missing = FALSE, types = c("character", "logical"))) {
#         valid = FALSE
#         msg = c(msg, "The element genotype$external in the slot \"annotation\" must be data frame.\n")
#       }
# 
#       validNames = c("alleleRef", "alleleAlt", "genotypeMismatch")
#       if (!testSubset(colnames(object@annotation$genotype$external)[1:3],  validNames)) {
#         valid = FALSE
#         msg = c(msg, "The element \"external\" in the slot \"annotation\" must contain only contain the following elements: ",paste0(validNames, collapse = ","),".\n")
#       }
# 
#       validNames = c()
#       for (i in seq_len(length(object@annotation$files))) {
#         validNames = c(validNames, paste0(object@annotation$files[[i]]$genotypeFile,":",object@annotation$files[[i]]$genotypeFileSampleName))
#       }
# 
#       nColsGenotype = ncol(object@annotation$genotype$external)
#       if (nColsGenotype > 3) {
#         testSubset(colnames(object@annotation$genotype$external)[4:nColsGenotype],  validNames)
#       }
# 
#       if (!testLogical(object@annotation$genotype$external$genotypeMismatch)) {
#         valid = FALSE
#         msg = c(msg, "The column genotypeMismatch in genotype$external in the slot \"annotation\" must be of type logical.\n")
#       }
# 
#       validValues = c("A","C","G","T", NA)
#       if (!testSubset(object@annotation$genotype$external$alleleRef, validValues)) {
#         valid = FALSE
#         msg = c(msg, "The column alleleRef in genotype$external in the slot \"annotation\" must only contain the values ",paste0(validValues, collapse = ","),".\n")
#       }
# 
#     }
# 
#     #3. heterozygosity
#     if (!testMatrix(object@annotation$genotype$heterozygosity, mode = "logical", nrows = nRows, ncols = nDatasets + 1, any.missing = TRUE)) {
#       valid = FALSE
#       msg = c(msg, "The element \"heterozygosity\" in genotype in the slot \"annotation\" must be a matrix containing only logical or NA values.\n")
#     }
# 
# 
# 
#     ## Regions ##  
#     if (!testClass(object@annotation$regions, "GRanges")) {
#       valid = FALSE
#       msg = c(msg, "Element \"regions\" in slot \"annotation\" must be an object of class \"GRanges\".\n")
#     } 
# 
# 
#     if (!assert(checkNull(object@annotation$regions), checkSubset(c("annotation", "SNPPos"), names(mcols(object@annotation$regions)), empty.ok = FALSE))) {
#       valid = FALSE
#       msg = c(msg, "Element \"regions\" in slot \"annotation\" must contain at least the metadata \"annotation\" and \"SNPPos\"\n")
#     }
# 
#     if (length(object@annotation$regions) < 1) {
#       valid = FALSE
#       msg = c(msg, "Element \"regions\" in slot \"annotation\" must contain at least one SNP region\n")
#     }  
#   }
# 
# 
#   # Slot config #
# 
#   if (testNull(object@config)) {
#     valid = FALSE
#     msg = c(msg, "Slot \"config\" must not be NULL\n")
#   }
# 
# 
#   if (!testList(object@config, min.len = 1)) {
#     valid = FALSE
#     msg = c(msg, "Slot \"config\" must be a list\n")
# 
#   } else {
# 
#     if (!testSubset(c("onlyPrepareForDatasetCorrelation", "input", names(.getListOfSupportedParameters())), 
#                     names(object@config))) {
#       valid = FALSE
#       msg = c(msg, "Slot \"config\" must contain all required parameters as elements, but at least one is missing.\n")
#     }
# 
#     if (!testDataFrame(object@config$input, types = c("logical","character"), min.cols = 4, min.rows = 1, col.names = "named") |
#         !testSubset(c("signal", "input", "individual", "genotype"), colnames(object@config$input)) ) {
# 
#       valid = FALSE
#       msg = c(msg, "Element \"input\" in \"config\" must contain the columns signal, input, genotype and individual.\n")
#     }
# 
# 
#   }
# 
#   # Slot INTERNAL #
# 
#   if (testNull(object@internal)) {
#     valid = FALSE
#     msg = c(msg, "Slot \"internal\" must not be NULL\n")
# 
#   } else {
# 
#     requiredElems = c("isAllelicRatio", 
#                       "sizeFactors", 
#                       "mergedReadGroups", 
#                       "plot_origBinSNPPosition",
#                       "plot_labelBins",
#                       "countType",
#                       "readWidth",
#                       "readStartPos",
#                       "globalBackground",
#                       "disableObjectIntegrityChecking",
#                       "addResultsElementsAdded"
#     )
# 
# 
#     if (!testSubset(requiredElems, names(object@internal))) {
#       valid = FALSE
#       msg = c(msg, "Slot \"internal\" must contain the elements ",paste0(requiredElems,collapse = ","),".\n")
#     }
# 
#     if (!testFlag(object@internal$isAllelicRatio)) {
#       valid = FALSE
#       msg = c(msg, "The element isAllelicRatio in the slot \"internal\" must be of type logical.\n")
#     }
# 
#     if (!testList(object@internal$sizeFactors, any.missing = FALSE, len = length(object@annotation$readGroups), types = "list") |
#         !testSetEqual(object@annotation$readGroups, names(object@internal$sizeFactors))) {
#       valid = FALSE
#       msg = c(msg, "Elements sizeFactors in the slot \"internal\" must be a list of lists with the number of elements identical to the element readGroups in slot annotation.\n")
#     }  
# 
#     if (!testFlag(object@internal$mergedReadGroups)) {
#       valid = FALSE
#       msg = c(msg, "The element mergedReadGroups in the slot \"internal\" must be of type logical.\n")
#     }
# 
# 
#     if (!testFlag(object@internal$disableObjectIntegrityChecking)) {
#       valid = FALSE
#       msg = c(msg, "The element disableObjectIntegrityChecking in the slot \"internal\" must be of type logical.\n")
#     }
# 
#     if (!testCharacter(object@internal$addResultsElementsAdded, min.chars = 1, any.missing = FALSE)) {
#       valid = FALSE
#       msg = c(msg, "The element addResultsElementsAdded in the slot \"internal\" must be of type character.\n")
#     }
# 
# 
#     if (!testNumeric(object@internal$plot_labelBins, any.missing = FALSE, len = length(object@annotation$bins), unique = TRUE)) {
#       valid = FALSE
#       msg = c(msg, "The element plot_labelBins in the slot \"internal\" must be of type character.\n")
#     }
# 
#     if (!testSubset(as.numeric(names(object@internal$plot_labelBins)), seq_len(length(object@annotation$bins)))) {
#       valid = FALSE
#       msg = c(msg, "The names of the element plot_labelBins in the slot \"internal\" must be from 1:", length(object@annotation$bins),".\n")
#     }
# 
#     if (!testNumber(object@internal$plot_origBinSNPPosition, lower = 1)) {
#       valid = FALSE
#       msg = c(msg, "The element plot_origBinSNPPosition in the slot \"internal\" must be a single number.\n")
#     }
# 
#     validValues = c("readCountsRaw", "readCountsNormalized", "enrichment")
#     if (!testSubset(object@internal$countType, validValues)) {
#       valid = FALSE
#       msg = c(msg, "The element countType in the slot \"internal\" must be one of ",paste0(validValues,collapse = ","), ".\n")
#     }
# 
#     if (object@internal$countType == "enrichment" & !parameters(object)$normByInput) {
#       valid = FALSE
#       msg = c(msg, "The element countType in the slot \"internal\" must match with the value of the parameter normByInput.\n")
#     }
# 
# 
#     if (length(object@annotation$readGroups) > 1) {
# 
#       for (i in seq_len(length(object@annotation$readGroups))) {
# 
#         if (!testSubset(names(object@internal$readStartPos)[i], object@annotation$readGroups[i])) {
#           valid = FALSE
#           msg = c(msg, "The names of the read groups are incorrect for the element readStartPos in the slot \"internal\":",paste0(names(object@internal$readStartPos),collapse = ",")," but expected ",paste0(object@annotation$readGroups,collapse = ","),".\n")
#         }
# 
#         for (j in seq_len(length(object@annotation$files))) {
# 
#           if (!testSubset(names(object@internal$readStartPos[[i]])[j], names(object@annotation$files)[j])) {
#             valid = FALSE
#             msg = c(msg, "The names of the datasets are incorrect for the element readStartPos in the slot \"internal\".\n")
#           }
# 
#           # Check length, should equal nRegions
#           if (length(object@internal$readStartPos[[i]][[j]]) != length(object@annotation$regions)) {
#             valid = FALSE
#             msg = c(msg, "The dimension of the element readStartPos are incorrect in the slot \"internal\".\n")
#           }
# 
# 
#           for (l in seq_len(length(object@annotation$regions))) {
#             if (!testInteger(object@internal$readStartPos[[i]][[j]][[l]], lower = 1, any.missing = FALSE)) {
#               valid = FALSE
#               msg = c(msg, paste0("At least one element in readStartPos in the slot \"internal\" is invalid and does not contain a vector of type Integer for read group ", i, ", and sample ", j, ".\n"))
#             } 
#           }
# 
#         }
#       }
#     }
#   }
# 
# 
#   # Slot additionalResults #
# 
#   #assertSubset(elemsToAnnotate, names(results), empty.ok = FALSE)
#   #assertLogical(results$parameters$calculateBackgroundDistribution, len = 1, any.missing = FALSE)
# 
#   #TODO
# 
# 
# 
#   #####################################################################
#   # Slot readCountsUnbinned and readCountsBinned and enrichmentBinned #
#   #####################################################################
# 
#   if (testNull(object@readCountsUnbinned) | testNull(object@readCountsBinned) | testNull(object@enrichmentBinned)) {
#     valid = FALSE
#     msg = c(msg, "Slots \"readCountsUnbinned\", \"readCountsBinned\" , and \"enrichmentBinned\"  must not be NULL\n")
# 
#   } else {
# 
#     # Test if names are correct
#     for (readGroup in object@annotation$readGroups) {
# 
#       testSetEqual(names(object@annotation$files), names(object@readCountsUnbinned[[readGroup]]))
# 
#       if (!object@config$onlyPrepareForDatasetCorrelation)  {
# 
#         testSetEqual(names(object@annotation$files), names(object@readCountsBinned[[readGroup]])) 
# 
#         # Test enrichment slot also
#         if (length(object@enrichmentBinned[[1]]) > 0) {
#           testSetEqual(names(object@annotation$files), names(object@enrichmentBinned[[readGroup]])) 
#         }
#       } 
# 
# 
# 
#     }
# 
#     if (!object@config$onlyPrepareForDatasetCorrelation)  {
# 
#       if (!testList(object@readCountsUnbinned, any.missing = FALSE, len = length(object@annotation$readGroups), types = "list") |
#           !testList(object@readCountsBinned,   any.missing = FALSE, len = length(object@annotation$readGroups), types = "list") |
#           !testList(object@enrichmentBinned,       any.missing = FALSE, len = length(object@annotation$readGroups), types = "list") |
#           !testSetEqual(object@annotation$readGroups, names(object@readCountsUnbinned)) | 
#           !testSetEqual(object@annotation$readGroups, names(object@readCountsBinned))    |
#           !testSetEqual(object@annotation$readGroups, names(object@enrichmentBinned))     
#       ) {
#         valid = FALSE
#         msg = c(msg, "Slots \"readCountsUnbinned\", \"readCountsBinned\" and \"enrichmentBinned\" must be named lists of lists with the number of elements identical to the element readGroups in slot annotation\n")
#       }
# 
#       # test validity and equality of names within counts of regions and bins
#       if (length(object@readCountsBinned[[1]]) > 0) {
# 
#         for (i in seq_len(length(object@annotation$readGroups))) {
#           if (!testSetEqual(names(object@readCountsUnbinned[[i]]), names(object@readCountsBinned[[i]]))) {
#             valid = FALSE
#             msg = c(msg, "Slots \"readCountsUnbinned\" and \"readCountsBinned\" must contain identical names.\n")
# 
#           }
#         }
# 
#       }  
# 
#       # Integrity checks for the enrichmentBinned slot
#       if (length(object@enrichmentBinned[[1]]) > 0) {
# 
#         for (i in seq_len(length(object@annotation$readGroups))) {
#           for (j in names(object@annotation$files)) {
# 
#             # Skip input files, enrichment has only been calculated for the signal files
#             if (object@annotation$files[[j]]$type == "input") {
#               next
#             }
# 
#             if (!j %in% names(object@enrichmentBinned[[i]])) {
#               valid = FALSE
#               msg = c(msg, "Could not find element ", j, " in slot \"enrichmentBinned\".\n")
# 
#             }
# 
# 
#           } # end all files
# 
#         } # end all read groups
# 
#       }  # end if enrichment slot not empty
#     }
# 
# 
# 
#     anyMatrixValueMissingAllowed = FALSE
#     if (object@internal$isAllelicRatio | object@config$normByInput) {
#       anyMatrixValueMissingAllowed = TRUE
#     } 
# 
#     # Test dimensions and matrix of all elements in readCountsUnbinned and readCountsBinned
#     for (i in seq_len(length(object@annotation$readGroups))) {
# 
#       for (j in names(object@annotation$files)) {
# 
#         # Skip check under certain conditions due to memory saving
#         if (!object@config$keepAllReadCounts & object@config$normByInput) {
# 
#           next
#         }
# 
#         # Test slot readCountsUnbinned
#         if (!testVector(object@readCountsUnbinned[[i]][[j]], strict = TRUE, any.missing = anyMatrixValueMissingAllowed, len = nRows)) {
#           valid = FALSE
#           msg = c(msg, paste0("Slot \"readCountsUnbinned\" must contain only numerical vectors of read counts for each element. However, at position [[", i,"]] [[", j,"]], a violation was found.\n") )              
#         }
# 
#         if (!object@internal$isAllelicRatio) {
#           if (!testNumeric(object@readCountsUnbinned[[i]][[j]], lower = 0, any.missing = FALSE, len = nRows)) {
#             valid = FALSE
#             msg = c(msg, paste0("Slot \"readCountsUnbinned\" must contain only numerical vectors of read counts for each element. However, at position [[", i,"]] [[", j,"]], a violation was found.\n"))               
#           }
#         } else {
#           if (!testNumeric(object@readCountsUnbinned[[i]][[j]], lower = 0, upper = 1, any.missing = TRUE, len = nRows)) {
#             valid = FALSE
#             msg = c(msg, paste0("Slot \"readCountsUnbinned\" must contain only numerical vectors of read counts for each element. However, at position [[", i,"]] [[", j,"]], a violation was found.\n"))            
#           }
#         }
# 
#         if (!object@config$onlyPrepareForDatasetCorrelation)  {
#           # Test slot readCountsBinned         
#           if (nCols > 1) {
#             if (!testMatrix(object@readCountsBinned[[i]][[j]], any.missing = anyMatrixValueMissingAllowed, nrows = nRows, ncols = nCols)) {
#               valid = FALSE
#               msg = c(msg, paste0("Slot \"readCountsBinned\" 1 must contain only matrices of read counts for each SNP and bin. However, at position [[", i,"]] [[", j,"]], a violation was found.\n") )             
# 
#             }
#           } else {
#             if (!testNumeric(object@readCountsBinned[[i]][[j]], any.missing = anyMatrixValueMissingAllowed, len = nRows)) {
#               valid = FALSE
#               msg = c(msg, paste0("Slot \"readCountsBinned\" 2 must contain only matrices of read counts for each SNP and bin. However, at position [[", i,"]] [[", j,"]], a violation was found.\n") )             
#             }
#           }
# 
# 
# 
#           if (!object@internal$isAllelicRatio) {
# 
#             if (object@internal$countType == "readCountsRaw") {
#               if (!testIntegerish(object@readCountsBinned[[i]][[j]], lower = 0, any.missing = FALSE)) {
#                 valid = FALSE
#                 msg = c(msg, paste0("Slot \"readCountsBinned\" 3 must contain only numerical vectors of read counts for each element. However, at position [[", i,"]] [[", j,"]], a violation was found.\n"))              
#               }
#             } else {
#               if (!testNumeric(object@readCountsBinned[[i]][[j]], lower = 0, any.missing = TRUE)) {
#                 valid = FALSE
#                 msg = c(msg, paste0("Slot \"readCountsBinned\" 4 must contain only numerical vectors of read counts for each element. However, at position [[", i,"]] [[", j,"]], a violation was found.\n"))             
#               }
#             }
# 
#           } else {
# 
#             if (!testNumeric(object@readCountsBinned[[i]][[j]], lower = 0, upper = 1, any.missing = TRUE)) {
#               valid = FALSE
#               msg = c(msg, paste0("Slot \"readCountsBinned\" 5 must contain only numerical vectors of read counts for each element. However, at position [[", i,"]] [[", j,"]], a violation was found.\n"))             
#             }
#           }
# 
# 
#         } # end if prepare only for correlation
# 
#       } # end all files
# 
# 
#       # Test slot enrichmentBinned
#       if (object@config$normByInput) {
# 
#         for (j in names(object@annotation$files)) {
# 
#           # Skip input files, enrichment has only been calculated for the signal files
#           if (object@annotation$files[[j]]$type == "input") {
#             next
#           }
# 
#           if (!testMatrix(object@enrichmentBinned[[i]][[j]], any.missing =  anyMatrixValueMissingAllowed, nrows = nRows, ncols = nCols)) {
#             valid = FALSE
#             msg = c(msg, paste0("Slot \"enrichmentBinned\" must contain only matrices of read counts for each SNP and bin. However, at position [[", i,"]] [[", j,"]], a violation was found.\n"))             
#           }
# 
#         }
#       }
# 
# 
# 
#     } # end all read groups
# 
#   }
# 
# 
#   if (valid) TRUE else msg
# 
# }
# 
# setValidity("SNPhood", .validSNPhoodObj)        
# 
# 

#' Get the number of peaks for a \code{\linkS4class{GRN}} object.
#' 
#' Return the number of peaks (all or only non-filtered ones) that are defined in the \code{\linkS4class{GRN}} object.
#'
#' @template GRN
#' @param filter TRUE or FALSE. Default TRUE. Should peaks marked as filtered be included in the count?
#' @return Integer. Number of peaks hat are defined in the \code{\linkS4class{GRN}} object, either by excluding (filter = TRUE) or including (filter = FALSE) peaks that are currently marked as \emph{filtered} (see method TODO)
#' @examples
# #' data(GRN, package="GRN")
#' nPeaks(GRN, filter = TRUE)
# #' nPeaks(GRN, filter = FALSE)
#' @export
#' @aliases peaks
#' @rdname peaks-methods
nPeaks <- function(GRN, filter = TRUE) {
  
  checkmate::assertClass(GRN, "GRN")
  checkmate::assertFlag(filter)
  
  if (is.null(GRN@data$peaks$consensusPeaks)) {
    return(NA)
  }
  
  if (filter) {
    nPeaks = GRN@data$peaks$consensusPeaks %>% dplyr::filter(!isFiltered) %>% nrow()
  } else {
    nPeaks = GRN@data$peaks$consensusPeaks %>% nrow()
  }
  
  nPeaks
  
}


#' Get the number of genes for a \code{\linkS4class{GRN}} object.
#' 
#' Return the number of genes (all or only non-filtered ones) that are defined in the \code{\linkS4class{GRN}} object.
#'
#' @template GRN
#' @param filter TRUE or FALSE. Default TRUE. Should genes marked as filtered be included in the count?
#' @return Integer. Number of genes hat are defined in the \code{\linkS4class{GRN}} object, either by excluding (filter = TRUE) or including (filter = FALSE) genes that are currently marked as \emph{filtered} (see method TODO)
#' @examples
# #' data(GRN, package="GRN")
#' nGenes(GRN, filter = TRUE)
# #' nGenes(GRN, filter = FALSE)
#' @export
#' @aliases genes
#' @rdname genes-methods
nGenes <- function(GRN, filter = TRUE) {
  
  checkmate::assertClass(GRN, "GRN")
  checkmate::assertFlag(filter)
  
  if (is.null(GRN@data$RNA$counts_norm.l)) {
    return(NA)
  }
  
  if (filter) {
    nGenes = GRN@data$RNA$counts_norm.l[["0"]] %>% dplyr::filter(!isFiltered) %>% nrow()
  } else {
    nGenes = GRN@data$RNA$counts_norm.l[["0"]] %>% nrow()
  }
  
  nGenes
  
}