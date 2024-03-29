% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fcts RSF.R
\name{simu_rev3}
\alias{simu_rev3}
\title{Temporal reversal (with THREE time periods) based on simulated coefficients}
\usage{
simu_rev3(simu, col1, col2, col3)
}
\arguments{
\item{simu}{An array of simulated individual coefficients based on their uncertainties(output of simu_coefs)}

\item{col1}{Column of the variable of interest for the 1st temporal period}

\item{col2}{Column of the variable of interest for the 2nd temporal period}

\item{col3}{Column of the variable of interest for the 3rd temporal period}
}
\value{
A matrix of reversal (one value for each covariates and replicates)
}
\description{
Calculate temporal reversal for each replicate of simulated coefficients
}
\examples{
data(goats)
goats$Season<-c("1", "2", "3")
ls1<-list()
ls1[[1]]<-as.formula(STATUS~(ELEVATION+SLOPE+ET+ASPECT+HLI+TASP):Season)
out<-rsf_ind(goats$ID, data=goats, form_ls=ls1)
coef<-ind_coef(m=1, out)
se<-ind_se(m=1, out)
simu<-simu_coefs(coef, se, n=100)
head(coef)
rev_elevation<-simu_rev3(simu, 3, 4,5) #Calculate specialization for elevation covariate
quantile(rev_elevation) #Show variation around estimate of elevation covariate
mean(rev_elevation) #Calculate average reversal for elevation covariate
}
