% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GAPIT_PC_selection.R
\name{gapit_results_in_filepath}
\alias{gapit_results_in_filepath}
\title{Identify trait names for GAPIT results.}
\usage{
gapit_results_in_filepath(path = ".")
}
\arguments{
\item{path}{File path to the csv files that GAPIT has created, a character
string or \code{file.path()}.}
}
\value{
A vector of phenotype names.
}
\description{
Creates a vector of phenotype names from GAPIT results in some directory.
}
\examples{
\dontrun{gapit_results_in_filepath(path = file.path("inst", "extdata"))}
\dontrun{gapit_results_in_filepath(path = "./inst/extdata")}

}
