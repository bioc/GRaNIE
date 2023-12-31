% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{add_featureVariation}
\alias{add_featureVariation}
\title{Quantify and interpret multiple sources of biological and technical variation for features (TFs, peaks, and genes) in a \code{\linkS4class{GRN}} object}
\usage{
add_featureVariation(
  GRN,
  formula = "auto",
  metadata = c("all"),
  features = "all_filtered",
  nCores = 1,
  forceRerun = FALSE,
  ...
)
}
\arguments{
\item{GRN}{Object of class \code{\linkS4class{GRN}}}

\item{formula}{Character(1). Either \code{auto} or a manually defined formula to be used for the model fitting. Default \code{auto}. Must include only terms that are part of the metadata as specified with the \code{metadata} parameter. If set to \code{auto}, the formula will be build automatically based on all metadata variables as specified with the \code{metadata} parameter. By default, numerical variables will be modeled as fixed effects, while variables that are defined as factors or can be converted to factors (characters and logical variables) are modeled as random effects as recommended by the \code{variancePartition} package.}

\item{metadata}{Character vector. Default \code{all}. Vector of column names from the metadata data frame that was provided when using the function 
\code{\link{addData}}. Must either contain the special keyword \code{all} to denote that all (!) metadata columns from \code{GRN@data$metadata} are taken
or if not, a subset of the column names from \code{GRN@data$metadata}to include in the model fitting for \code{fitExtractVarPartModel}..}

\item{features}{Character(1). Either \code{all_filtered} or \code{all}. Default \code{all_filtered}. Should \code{variancePartition} only be run for the features (TFs, peaks and genes) from the filtered set of connections (the result of \code{\link{filterGRNAndConnectGenes}}) or for all genes that are defined in the object? If set to \code{all}, the running time is greatly increased.}

\item{nCores}{Integer >0. Default 1. Number of cores to use. 
A value >1 requires the \code{BiocParallel} package (as it is listed under \code{Suggests}, it may not be installed yet).}

\item{forceRerun}{\code{TRUE} or \code{FALSE}. Default \code{FALSE}. Force execution, even if the GRN object already contains the result. Overwrites the old results.}

\item{...}{Additional parameters passed on to \code{variancePartition::fitExtractVarPartModel} beyond \code{exprObj}, \code{formula} and \code{data}. See the function help for more information}
}
\value{
An updated \code{\linkS4class{GRN}} object, with additional information added from this function to \code{GRN@stats$variancePartition} as well as the elements \code{genes}, \code{consensusPeaks} and \code{TFs} within \code{GRN@annotation}. 
As noted above, the results are not added to \code{GRN@connections$all.filtered}; rerun the function \code{\link{getGRNConnections}} and set \code{include_variancePartitionResults} to \code{TRUE} to include the results in the eGRN output table.
}
\description{
Runs the main function \code{fitExtractVarPartModel} of the package \code{variancePartition}: Fits a linear (mixed) model to estimate contribution of multiple sources of variation while simultaneously correcting for all other variables for the features in a GRN object (TFs, peaks, and genes) given particular metadata. The function reports the fraction of variance attributable to each metadata variable.
\strong{Note: The results are not added to \code{GRN@connections$all.filtered}, rerun the function \code{\link{getGRNConnections}} and set \code{include_variancePartitionResults} to \code{TRUE} to do so}.
The results object is stored in \code{GRN@stats$variancePartition} and can be used for the various diagnostic and plotting functions from \code{variancePartition}.
}
\details{
The normalized count matrices are used as input for \code{fitExtractVarPartModel}.
}
\examples{
# See the Workflow vignette on the GRaNIE website for examples
# GRN = loadExampleObject()
# GRN = add_featureVariation(GRN, metadata = c("mt_frac"), forceRerun = TRUE)
}
\seealso{
\code{\link{addData}}

\code{\link{getGRNConnections}}
}
