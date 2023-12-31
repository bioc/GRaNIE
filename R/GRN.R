#' \strong{GRaNIE} (\strong{G}ene \strong{R}egul\strong{a}tory \strong{N}etwork \strong{I}nference including \strong{E}nhancers): Reconstruction and evaluation of data-driven, cell type specific gene regulatory networks including enhancers using chromatin accessibility and RNAseq data (general package information)
#'
#' Genetic variants associated with diseases often affect non-coding regions, thus likely having a regulatory role. To understand the effects of genetic variants in these regulatory regions, identifying genes that are modulated by specific regulatory elements (REs) is crucial. The effect of gene regulatory elements, such as enhancers, is often cell-type specific, likely because the combinations of transcription factors (TFs) that are regulating a given enhancer have celltype specific activity. This TF activity can be quantified with existing tools such as \code{diffTF} and captures differences in binding of a TF in open chromatin regions. Collectively, this forms a gene regulatory network (\code{eGRN}) with cell-type and data-specific TF-RE and RE-gene links. Here, we reconstruct such a \code{eGRN} using bulk RNAseq and open chromatin (e.g., using ATACseq or ChIPseq for open chromatin marks) and optionally TF activity data. Our network contains different types of links, connecting TFs to regulatory elements, the latter of which is connected to genes in the vicinity or within the same chromatin domain (TAD). We use a statistical framework to assign empirical FDRs and weights to all links using a permutation-based approach.
#' 
#' @section Package functions:
#' See the Vignettes for a workflow example and more generally \url{https://grp-zaugg.embl-community.io/GRaNIE} for all project-related information.
#' 
#' @section GRN object:
#' The \code{GRaNIE} package works with \code{GRN} objects. See \code{\linkS4class{GRN}} for details.
#' 
#' @section Contact Information: 
#' Please check out \url{https://grp-zaugg.embl-community.io/GRaNIE} for how to get in contact with us.
#' 
#' @docType package
#' @keywords GRaNIE, GRaNIE-package
#' @name GRaNIE

NULL

