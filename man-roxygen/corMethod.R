#' @param corMethod Character. One of \code{pearson}, \code{spearman} or \code{bicor}. Default \code{pearson}. Method for calculating the correlation coefficient. 
#' For \code{pearson} and \code{spearman} , see \link{cor} for details. \code{bicor} denotes the *biweight midcorrelation*, a correlation measure based on medians as
#' calculated by \code{WGCNA::bicorAndPvalue}. Both \code{spearman} and \code{bicor} are considered more robust measures that are less prone to be affected by outliers.
