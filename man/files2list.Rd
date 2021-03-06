% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fcts RSF.R
\name{files2list}
\alias{files2list}
\title{Convert a folder with individual files to object provided by ssf_ind or rsf_ind}
\usage{
files2list(lsfiles, cleanModel = F)
}
\arguments{
\item{lsfiles}{The list of files to be imported, for example using the dir() command.}

\item{cleanModel}{Whether the cleanModel function should be applied. .}
}
\value{
A list of the same format as ssf_ind or rsf_ind.
}
\description{
This function convert a folder of individual files (for example if analysis was ran on a high performance cluster) to the same format provided by ssf_ind or rsf_ind. This facilitate the use of other functions such as aictab_ind and pop_avg.
}
