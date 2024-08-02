#' @param gene.types  Character vector of supported gene types. Default \code{"all"}. 
#' Filter for gene types to retain, genes with gene types not listed here are filtered. The special keyword "all" indicates no filter and retains all gene types.
#' To filter for only protein-coding genes, use \code{"protein_coding"}.
#' The specified names must match the names as stored in the \code{\linkS4class{GRN}} object (\code{see GRN@annotation$genes$gene.type}) and
#' correspond 1:1 to the gene type names as provided by \code{biomaRt}, with the exception of \code{lncRNAs}, 
#' which is internally renamed to \code{lincRNAs} when first fetching all gene types. This is done due to a recent change in \code{biomaRt} and aims at 
#' keeping backwards compatibility with \code{\linkS4class{GRN}} objects.
