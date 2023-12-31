% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{deleteIntermediateData}
\alias{deleteIntermediateData}
\title{Optional convenience function to delete intermediate data from the function \code{\link{AR_classification_wrapper}} and summary statistics that may occupy a lot of space}
\usage{
deleteIntermediateData(GRN)
}
\arguments{
\item{GRN}{Object of class \code{\linkS4class{GRN}}}
}
\value{
An updated \code{\linkS4class{GRN}} object, with some slots being deleted (\code{GRN@data$TFs$classification} as well as \code{GRN@stats$connectionDetails.l})
}
\description{
Optional convenience function to delete intermediate data from the function \code{\link{AR_classification_wrapper}} and summary statistics that may occupy a lot of space
}
\examples{
# See the Workflow vignette on the GRaNIE website for examples
GRN = loadExampleObject()
GRN = deleteIntermediateData(GRN)
}
