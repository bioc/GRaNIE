% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network.R
\name{calculateGeneralEnrichment}
\alias{calculateGeneralEnrichment}
\title{Run an enrichment analysis for the genes in the whole network in the filtered \code{\linkS4class{GRN}} object}
\usage{
calculateGeneralEnrichment(
  GRN,
  ontology = c("GO_BP", "GO_MF"),
  algorithm = "weight01",
  statistic = "fisher",
  background = "neighborhood",
  background_geneTypes = "all",
  pAdjustMethod = "BH",
  forceRerun = FALSE
)
}
\arguments{
\item{GRN}{Object of class \code{\linkS4class{GRN}}}

\item{ontology}{Character vector of ontologies. Default \code{c("GO_BP", "GO_MF")}. 
Valid values are \code{"GO_BP"}, \code{"GO_MF"}, \code{"GO_CC"}, \code{"KEGG"}, \code{"DO"}, and \code{"Reactome"}, 
referring to \emph{GO Biological Process}, \emph{GO Molecular Function}, \emph{GO Cellular Component}, \emph{KEGG}, \emph{Disease Ontology}, 
and \emph{Reactome Pathways}, respectively. \code{GO} ontologies require the \code{topGO}, 
\code{"KEGG"} the \code{clusterProfiler}, \code{"DO"} the \code{DOSE}, and \code{"Reactome"} the \code{ReactomePA} packages, respectively.
As they are listed under \code{Suggests}, they may not yet be installed, and the function will throw an error if they are missing.}

\item{algorithm}{Character. Default \code{"weight01"}. One of: \code{"classic"}, \code{"elim"}, \code{"weight"}, \code{"weight01"}, \code{"lea"}, \code{"parentchild"}. Only relevant if ontology is GO related (GO_BP, GO_MF, GO_CC), ignored otherwise. Name of the algorithm that handles the GO graph structures. Valid inputs are those supported by the \code{topGO} library. 
For general information about the algorithms, see \url{https://academic.oup.com/bioinformatics/article/22/13/1600/193669}. \code{weight01} is a mixture between the \code{elim} and the \code{weight} algorithms.}

\item{statistic}{Character. Default \code{"fisher"}. One of: \code{"fisher"}, \code{"ks"}, \code{"t"}. Statistical test to be used. Only relevant if ontology is GO related (\code{GO_BP}, \code{GO_MF}, \code{GO_CC}), and valid inputs are a subset of those supported by the \code{topGO} library (we had to remove some as they do not seem to work properly in \code{topGO} either), ignored otherwise. For the other ontologies the test statistic is always Fisher.}

\item{background}{Character. Default \code{"neighborhood"}. One of: \code{"all_annotated"}, \code{"all_RNA"}, \code{"all_RNA_filtered"}, \code{"neighborhood"}. Set of genes to be used to construct the background for the enrichment analysis. This can either be all annotated genes in the reference genome (\code{all_annotated}), all genes from the provided RNA data (\code{all_RNA}), all genes from the provided RNA data excluding those marked as filtered after executing \code{filterData} (\code{all_RNA_filtered}), or all the genes that are within the neighborhood of any peak (before applying any filters except for the user-defined \code{promoterRange} value in \code{addConnections_peak_gene}) (\code{neighborhood}).}

\item{background_geneTypes}{Character vector of gene types that should be considered for the background. Default \code{"all"}. 
Only gene types as defined in the \code{\linkS4class{GRN}} object, slot \code{GRN@annotation$genes$gene.type} are allowed. 
The special keyword \code{"all"} means no filter on gene type.}

\item{pAdjustMethod}{Character. Default \code{"BH"}. One of: \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"}. This parameter is only relevant for the following ontologies: KEGG, DO, Reactome. For the other ontologies, the algorithm serves as an adjustment.}

\item{forceRerun}{\code{TRUE} or \code{FALSE}. Default \code{FALSE}. Force execution, even if the GRN object already contains the result. Overwrites the old results.}
}
\value{
An updated \code{\linkS4class{GRN}} object, with the enrichment results stored in the \code{stats$Enrichment$general} slot.
}
\description{
The enrichment analysis is based on the whole network, see \code{\link{calculateCommunitiesEnrichment}} and \code{\link{calculateTFEnrichment}} for 
community- and TF-specific enrichment, respectively.
This function requires the existence of the eGRN graph in the \code{\linkS4class{GRN}} object as produced by \code{\link{build_eGRN_graph}}.
Results can subsequently be visualized with the function \code{\link{plotGeneralEnrichment}}.
}
\details{
All enrichment functions use the TF-gene graph as defined in the `GRN` object. See the `ontology` argument for currently supported ontologies.
Also note that some parameter combinations for `algorithm` and `statistic` are incompatible, an error message will be thrown in such a case.
}
\examples{
# See the Workflow vignette on the GRaNIE website for examples
GRN =  loadExampleObject()
GRN =  calculateGeneralEnrichment(GRN, ontology = "GO_BP", forceRerun = FALSE)
}
\seealso{
\code{\link{plotGeneralEnrichment}}

\code{\link{calculateCommunitiesEnrichment}}

\code{\link{calculateTFEnrichment}}

\code{\link{plotCommunitiesEnrichment}}
}
