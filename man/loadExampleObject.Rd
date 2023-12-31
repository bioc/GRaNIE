% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{loadExampleObject}
\alias{loadExampleObject}
\title{Load example GRN dataset}
\usage{
loadExampleObject(
  forceDownload = FALSE,
  fileURL = "https://git.embl.de/grp-zaugg/GRaNIE/-/raw/master/data/GRN.rds"
)
}
\arguments{
\item{forceDownload}{\code{TRUE} or \code{FALSE}. Default \code{FALSE}. Should the download be enforced even if the local cached file is already present?}

\item{fileURL}{Character. Default \url{https://git.embl.de/grp-zaugg/GRaNIE/-/raw/master/data/GRN.rds}. URL to the GRN example object in rds format.}
}
\value{
An small example \code{\linkS4class{GRN}} object
}
\description{
Loads an example GRN object with 6 TFs, ~61.000 peaks, ~19.000 genes, 259 filtered connections and pre-calculated enrichments. 
This function uses \code{BiocFileCache} if installed to cache the example object, which is 
considerably faster than re-downloading the file anew every time the function is executed.
If not, the file is re-downloaded every time anew. Thus, to enable caching, you may install the package \code{BiocFileCache}.
}
\examples{
GRN = loadExampleObject()
}
