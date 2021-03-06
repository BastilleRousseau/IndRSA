% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fcts RSF.R
\name{rf}
\alias{rf}
\title{Extraction of proximity from random forest classification}
\usage{
rf(coef, ntree = 10000, ...)
}
\arguments{
\item{coef}{A matrix of model coefficient (from function ind_coef)}
}
\value{
A proximity matrix
}
\description{
Apply random forest classification
}
\examples{
data(goats)
ls1<-list()
ls1[[1]]<-as.formula(STATUS~ELEVATION+SLOPE+ET+ASPECT+HLI+TASP)
ls1[[2]]<-as.formula(STATUS~ET+ASPECT+HLI+TASP)
out<-rsf_ind(goats$ID, data=goats, form_ls=ls1)
coef<-ind_coef(m=1, out)
prox<-rf(coef)
}
