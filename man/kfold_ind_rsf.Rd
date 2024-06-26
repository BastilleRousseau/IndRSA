% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fcts RSF.R
\name{kfold_ind_rsf}
\alias{kfold_ind_rsf}
\title{Perform kfold cross-validation at the individual level for an RSF .}
\usage{
kfold_ind_rsf(
  m = 1,
  mod_ls,
  ls = ls,
  cutoff = 0,
  k = 5,
  nrepet = 5,
  nbins = 10,
  grph = T
)
}
\arguments{
\item{m}{model number (based on number in list of formula provided to rsf_ind)}

\item{mod_ls}{A list of list of model generated by rsf_ind}

\item{ls}{A list of list of formula fed into rsf_ind}

\item{cutoff}{A cutoff value to exclude individuals with bad fit, default = -1 indicating model that did not converge will be excluded. Values > 0 will exclude based on coefficient}

\item{k}{number of fold (default = 5)}

\item{nrepet}{Number of repetitions (default =10)}

\item{nbins}{Number of bins (default =10)}
}
\value{
A vector of mean kfold score for each individual and an (optional) histogram
}
\description{
Perform kfold cross-validation at the individual level and return histogram and mean kfold accross individual
}
\examples{
data(goats)
ls1<-list()
ls1[[1]]<-as.formula(STATUS~ELEVATION+SLOPE+ET+ASPECT+HLI+TASP)
ls1[[2]]<-as.formula(STATUS~ET+ASPECT+HLI+TASP)
out<-rsf_ind(goats$ID, data=goats, form_ls=ls1)
kfold_ind_rsf(m=1, out, ls=ls1)
}
