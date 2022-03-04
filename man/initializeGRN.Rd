% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{initializeGRN}
\alias{initializeGRN}
\title{Initialize a \code{\linkS4class{GRN}} object}
\usage{
initializeGRN(objectMetadata = list(), outputFolder, genomeAssembly)
}
\arguments{
\item{objectMetadata}{List. Default list(). Optional (named) list with an arbitrary number of elements, all of which capture metadata for the object. This is mainly used to distinguish GRN objects from one another by storing object-specific metadata along with the data.}

\item{outputFolder}{Absolute file path. No default. Default output folder where all pipeline output will be put unless specified otherwise. We recommend specifying an absolute path. Note that for Windows-based systems, the path must be correctly specified with "/" as path separator.}

\item{genomeAssembly}{Character. No default. The genome assembly of all data that to be used within this object. Currently, supported genomes are: \code{hg19}, \code{hg38}, and \code{mm10}}
}
\value{
Empty \code{\linkS4class{GRN}} object
}
\description{
Initialize a \code{\linkS4class{GRN}} object
}
\examples{
meta.l = list(name = "exampleName", date = "01.03.22")
GRN = initializeGRN(objectMetadata = meta.l, outputFolder = "output", genomeAssembly = "hg38")
}