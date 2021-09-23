#
# # ## Package creation -
# # install.packages("devtools")
# # library("devtools")
#  devtools::install_github("klutometis/roxygen")
#  library(roxygen2)
#  setwd("C:/Users/Guillaume/OneDrive/Elephant/Analysis/RSF")
#  create("IndRSA")
# # # #
# # # #
# # # # Add Files
# setwd("./IndRSA")
# devtools::document()
# 1#
# # # Create pdf
#   pack <- "IndRSA"
#   path <- find.package(pack)
#   system(paste(shQuote(file.path(R.home("bin"), "R")),"CMD", "Rd2pdf", shQuote(path)))
#


# ####################################
# #Clean model (#Taken from: https://www.r-bloggers.com/trimming-the-fat-from-glm-models-in-r/)
# ###################################
#' Remove elements from glm object to save space
#'
#' Remove elements from glm objects (Taken from: https://www.r-bloggers.com/trimming-the-fat-from-glm-models-in-r/)
#' @param cm A glm object
#' @return A glm object
#' @export
cmodel = function(cm) {
  cm$y = c()
  cm$model = c()
  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
  cm$qr = c()
  cm$linear.predictors = c()
  cm$weights = c()
  cm$prior.weights = c()
  cm$data = c()
  cm
}

#' Convert a folder with individual files to object provided by ssf_ind or rsf_ind
#'
#' This function convert a folder of individual files (for example if analysis was ran on a high performance cluster) to the same format provided by ssf_ind or rsf_ind. This facilitate the use of other functions such as aictab_ind and pop_avg.
#' @param lsfiles The list of files to be imported, for example using the dir() command.
#' @param cleanModel Whether the cleanModel function should be applied. .
#' @return A list of the same format as ssf_ind or rsf_ind.
#' @export
files2list<-function(lsfiles, cleanModel=F) {
  out4<-list()
  pb = txtProgressBar(min = 0, max = length(lsfiles), initial = 0, style=3)
  for (i in 1:length(lsfiles)) {
    out<-get(load(lsfiles[[i]]))
    out2<-unlist(lapply(out, try(logLik)))
    try(out3<-lapply(out, function(x) sqrt(diag(vcov(x)))))
    if (cleanModel) try((out<-lapply(out, cmodel_ssf)))  #Only for ssf need to adjust here for RSF
    out4[[i]]<-list(out, out2, out3)
    setTxtProgressBar(pb,i)
  }
  names(out4)<-lsfiles
  return(out4)}

#' Remove elements from clogit object to save space
#'
#' Remove elements from glm objects
#' @param cm A coxph object
#' @return A coxph object
#' @export
cmodel_ssf = function(cm) {
  cm$linear.predictors = c()
  cm$residuals = c()
  cm$means = c()
  cm$nevent = c()
  cm$y = c()
  cm$formula = c()
  cm
}

#' Apply a list of candidate RSF models to a single individual
#'
#' Apply a list of candidate RSF models to a single individual
#' @param sub A subset of data from a single individual
#' @param form_ls A list of formulas for the different candidate models
#' @param cleamModel Whether the model should be "cleaned" to save memory space
#' @param method Weither typical glm or bias-reduction glm should be fitted (see package brglm)
#' @return A list of glm objects
#' @export
rsf_mod<-function(sub, form_ls, cleanModel=F, method=method) {
  out<-lapply(1:length(form_ls), function(x) brglm::brglm(form_ls[[x]], data=sub, family=binomial, link="logit", method=method, control.glm=glm.control(maxit=1000)))
  out2<-unlist(lapply(out, logLik))
  try(out3<-lapply(out, function(x) sqrt(diag(vcov(x)))))
  if (cleanModel) (out<-lapply(out, cmodel))
  return(list(out, out2, out3))
}

#' Apply a list of candidate SSF models to a single individual
#'
#' Apply a list of candidate SSF models to a single individual
#' @param sub A subset of data from a single individual
#' @param form_ls A list of formulas for the different candidate models
#' @param cleamModel Whether the model should be "cleaned" to save memory space
#' @param method Whether exact or approximate ML should be performed (see package survival)
#' @return A list of coxph objects
#' @export
ssf_mod<-function(sub, form_ls, cleanModel=F, method="approximate") {
  require(survival)
  out<-lapply(1:length(form_ls), function(x) survival::clogit(form_ls[[x]], data=sub, method=method))
  out2<-unlist(lapply(out, try(logLik)))
  try(out3<-lapply(out, function(x) sqrt(diag(vcov(x)))))
  if (cleanModel) try((out<-lapply(out, cmodel_ssf)))
  return(list(out, out2, out3))
}


#' Apply a list of candidate models to multiple individuals
#'
#' Apply rsf_mod to each individual of a dataset
#' @param id A vector indicating the individuals
#' @param data The dataset containing all data
#' @param form_ls A list of formulas for the different candidate models
#' @param cleamModel Whether the model should be "cleaned" to save memory space (default = F)
#' @param method Weither typical glm or bias-reduction glm should be fitted (default="glm.fit) (see package brglm)
#' @return A list of list of glm objects
#' @export
#' @examples
#' data(goats)
#' ls1<-list()
#' ls1[[1]]<-as.formula(STATUS~ELEVATION+SLOPE+ET+ASPECT+HLI+TASP)
#' ls1[[2]]<-as.formula(STATUS~ET+ASPECT+HLI+TASP)
#' out<-rsf_ind(goats$ID, data=goats, form_ls=ls1)
rsf_ind<-function(id,data, form_ls, cleanModel=F, method="glm.fit") { #id is a vector of nrow(data)
  if(length(id) != nrow(data)) (stop("id should be the same length as data"))
  id1<-sort(unique(id))
  out<-pbapply::pblapply(1:length(id1), function(x) rsf_mod(sub=data[id==id1[x],], form_ls=form_ls, method=method, cleanModel=cleanModel))
  names(out)<-id1
  return(out)
}

#' Apply a list of SSF candidate models to multiple individuals
#'
#' Apply ssf_mod to each individual of a dataset
#' @param id A vector indicating the individuals
#' @param data The dataset containing all data
#' @param form_ls A list of formulas for the different candidate models
#' @param cleamModel Whether the model should be "cleaned" to save memory space (default = F)
#' @param method Whether exact or approximate ML should be performed (see package survival)
#' @return A list of of coxph objects
#' @export
#' @examples
ssf_ind<-function(id,data, form_ls, cleanModel=T, method="approximate") { #id is a vector of nrow(data)
  if(length(id) != nrow(data)) (stop("id should be the same length as data"))
  id1<-sort(unique(id))
  out<-pbapply::pblapply(1:length(id1), function(x) try(ssf_mod(sub=data[id==id1[x],], form_ls=form_ls, method=method, cleanModel=cleanModel)))
  names(out)<-id1
  return(out)
}


#' Identify potential individual with bad fits based on coefficients and model convergence
#'
#' Identify potential individual with bad fits based on coefficients and model convergence
#' @param mod_ls A list of list of model generated by rsf_ind
#' @param cutoff Value a coeffiecient may take to indicate bad fit (default=1000)
#' @return A list with first element giving individuals with bad fits based on coefficients and second element containing individuals with bad fit based on convergence
#' @export
#' @examples
#' data(goats)
#' ls1<-list()
#' ls1[[1]]<-as.formula(STATUS~ELEVATION+SLOPE+ET+ASPECT+HLI+TASP)
#' ls1[[2]]<-as.formula(STATUS~ET+ASPECT+HLI+TASP)
#' out<-rsf_ind(goats$ID, data=goats, form_ls=ls1)
#' bad_fit(out, cutoff=20) #None
bad_fit<-function(mod_ls, cutoff=1000) {
  tt<-unlist(lapply(mod_ls, function(x) sum(abs(unlist(lapply(x[[1]], coef)))>cutoff, na.rm=T)))
  tt2<-unlist(lapply(mod_ls, function(x) sum(unlist(lapply(x[[1]], function(y) y$converged))==FALSE, na.rm=T)))
  names(tt)<-names(mod_ls)
  names(tt2)<-names(mod_ls)
  return(list(tt[which(tt>0)],tt2[which(tt2>0)]  ))
}

#' Identify potential individual with bad fits based on coefficients and model convergence for a specific model
#'
#' Identify potential individual with bad fits based on coefficients and model convergence
#' @param mod_ls A list of list of model generated by rsf_ind
#' @param cutoff Value a coeffiecient may take to indicate bad fit (default=1000)
#' @param m model number (based on number in list of formula provided to rsf_ind)
#' @return A list with first element giving individuals with bad fits based on coefficients and second element containing individuals with bad fit based on convergence
#' @export
#' @examples
#' data(goats)
#' ls1<-list()
#' ls1[[1]]<-as.formula(STATUS~ELEVATION+SLOPE+ET+ASPECT+HLI+TASP)
#' ls1[[2]]<-as.formula(STATUS~ET+ASPECT+HLI+TASP)
#' out<-rsf_ind(goats$ID, data=goats, form_ls=ls1)
#' bad_fit(out, cutoff=20, m=1) #None
bad_fit1<-function(mod_ls, cutoff=1000, m=1) {
  tt<-unlist(lapply(mod_ls, function(x) sum(abs(coef(x[[1]][[m]]))>cutoff, na.rm=T)))
  tt2<-unlist(lapply(mod_ls, function(x) sum(x[[1]][[m]]$converged==FALSE, na.rm=T)))
  names(tt)<-names(mod_ls)
  names(tt2)<-names(mod_ls)
  return(list(tt[which(tt>0)],tt2[which(tt2>0)]  ))
}


#' Remove potential individual with bad fits based on coefficients
#'
#' Remove potential individual with bad fits based on coefficients
#' @param mod_ls A list of list of model generated by rsf_ind
#' @param cutoff Value a coeffiecient may take to indicate bad fit (default=1000)
#' @return A list excluding individual with bad fits.
#' @export
rm_bad_fit<-function(mod_ls, cutoff=1000) {
  tt<-unlist(lapply(mod_ls, function(x) sum(abs(unlist(lapply(x[[1]], coef)))>cutoff, na.rm=T)))
  tt2<-which(tt==0)
  return(mod_ls[tt2])
}

#' Remove potential individual with bad fits based on coefficients for a specific model
#'
#' Remove potential individual with bad fits based on coefficients
#' @param mod_ls A list of list of model generated by rsf_ind
#' @param cutoff Value a coeffiecient may take to indicate bad fit (default=1000)
#' @param m model number (based on number in list of formula provided to rsf_ind)
#' @return A list excluding individual with bad fits.
#' @export
rm_bad_fit1<-function(mod_ls, cutoff=1000, m=1) {
  tt<-unlist(lapply(mod_ls, function(x) sum(abs(coef(x[[1]][[m]]))>cutoff, na.rm=T)))
  tt2<-which(tt==0)
  return(mod_ls[tt2])
}

#' Remove potential individual with bad fits based on model convergence
#'
#' Remove potential individual with bad fits based on model convergence
#' @param mod_ls A list of list of model generated by rsf_ind
#' @return A list excluding individual with bad fits.
#' @export
rm_conv_fit<-function(mod_ls) {
  tt<-unlist(lapply(mod_ls, function(x) sum(unlist(lapply(x[[1]], function(y) y$converged))==FALSE, na.rm=T)))
  tt2<-which(tt==0)
  return(mod_ls[tt2])
}

#' Remove potential individual with bad fits based on model convergence for a specific model
#'
#' Remove potential individual with bad fits based on model convergence

#' @param m model number (based on number in list of formula provided to rsf_ind)
#' @return A list excluding individual with bad fits.
#' @export
rm_conv_fit1<-function(mod_ls, m=1) {
  tt<-unlist(lapply(mod_ls, function(x) sum(x[[1]][[m]]$converged==FALSE, na.rm=T)))
  tt2<-which(tt==0)
  return(mod_ls[tt2])
}


#' Perform model selection over all individuals
#'
#' Perform AIC model selection over all individuals by adding up likelihood of individual model (based code partly taken from package AICcmodavg)
#' @param mod_ls A list of list of model generated by rsf_ind
#' @param cutoff A cutoff value to exclude individuals with bad fit, default = -1 indicating model that did not converge will be excluded. Values > 0 will exclude based on coefficient values
#' @return A AIC model selection table
#' @export
#' @examples
#' data(goats)
#' ls1<-list()
#' ls1[[1]]<-as.formula(STATUS~ELEVATION+SLOPE+ET+ASPECT+HLI+TASP)
#' ls1[[2]]<-as.formula(STATUS~ET+ASPECT+HLI+TASP)
#' out<-rsf_ind(goats$ID, data=goats, form_ls=ls1)
#' aictab_ind(out)
aictab_ind<-function(mod_ls, cutoff=0, K=NULL) {
  if(cutoff>0) {mod_ls<-rm_bad_fit(mod_ls, cutoff=cutoff)}
  if(cutoff==-1) {mod_ls<-rm_conv_fit(mod_ls)}
  Aiccf<-function(k, n, loglik) {
    2*k-2*loglik+(2*k*(k+1))/(n-k-1)}

  ll_ls<-function(mod_ls) {
    n_mod<-length(mod_ls[[1]][[1]])
    n_id<-length(mod_ls)
    ll_ls<-unlist(lapply(1:n_mod, function(y) sum(unlist(lapply(mod_ls, function(x) x[[2]][y])))))
    return(ll_ls)
  }

  n_mod<-length(mod_ls[[1]][[1]])
  n_id<-length(mod_ls)
  #h <- function(w) if( any( grepl( "is.na", w) ) ) invokeRestart( "muffleWarning" )
  Results <- NULL
  Results <- data.frame(Modnames = 1:n_mod)
  Results$LL <-ll_ls(mod_ls)
  Results$K <- unlist(lapply(mod_ls[[1]][[1]], function(x) length(coef(x))))
  if(!is.null(K)) { Results$K<-K}
  Results$AICc <- Aiccf(Results$K, n_id, Results$LL)
  Results$Delta_AICc <- Results$AICc - min(Results$AICc)
  Results$ModelLik <- exp(-0.5 * Results$Delta_AICc)
  Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)
  return(Results)
}



#' Extract population average of top model and extract individual coefficients
#'
#' Extract population average of top model and extract individual coefficients. Population average can be calculated based on bootstrap (Prokopenko et al. 2016 JAppEco) or weighted based on standard errors (Murtaugh 2007 Ecology)
#' @param m model number (based on number in list of formula provided to rsf_ind)
#' @param mod_ls A list of list of model generated by rsf_ind
#' @param cutoff A cutoff value to exclude individuals with bad fit, default = -1 indicating model that did not converge will be excluded. Values > 0 will exclude based on coefficient values
#' @param method If = "boot", population average is based on bootstrap, if = "murtaugh" based on standard errors weighting. See Prokopenko et al 2016 or Murtaugh 2007 for details.
#' @param nboot Number of bootstrap iterations, default = 1000. Only applicable if method = "boot".
#' @param id_year Whether id_year (instead of individual) are provided. Individual and year needs to be separated by an underscore for the function to work properly.
#' @return A list containing a table population average with confidence intervals and a table of individual coefficients
#' @export
#' @examples
#' data(goats)
#' ls1<-list()
#' ls1[[1]]<-as.formula(STATUS~ELEVATION+SLOPE+ET+ASPECT+HLI+TASP)
#' ls1[[2]]<-as.formula(STATUS~ET+ASPECT+HLI+TASP)
#' out<-rsf_ind(goats$ID, data=goats, form_ls=ls1)
#' pop_avg(m=1, out, method="murtaugh")
pop_avg<-function(m=1, mod_ls,cutoff=-1, nudge=0.01, method="murtaugh", nboot=1000, id_year=F) {
  if(cutoff>0) {mod_ls<-rm_bad_fit1(mod_ls, cutoff=cutoff, m=m)}
  if(cutoff==-1) {mod_ls<-rm_conv_fit1(mod_ls,  m=m)}
  i<-length(mod_ls)
  add_weights<-function(x, name=name) {
    t<-data.frame(table(x$name))
    t$Freq<-1/t$Freq
    out<-merge(x, t, by.x="name", by.y="Var1")
    return(out)
  }
  co<-lapply(1:i, function(x) data.frame(t(data.frame(coef(mod_ls[[x]][[1]][[m]])))))
   coef2<-plyr::rbind.fill(co)
  rownames(coef2)<-names(mod_ls)
  if(id_year==T) {coef2$name<-matrix(unlist(strsplit(as.character(rownames(coef2)), "_")), ncol=2, byrow=T)[,1]} ##### Needs to be adjusted when not Id year
  if(id_year==F) {coef2$name<-rownames(coef2)}
  coef2$ID<-names(mod_ls)
  coef2<-add_weights(coef2)

  if(method=="boot") {
    mean_ci_boot<-function(boot, lci=0.025, uci=0.975) {
      n<-length(boot[[1]])
      out<-data.frame()
      for (i in 1:n) {
        out<-rbind(out,
                   c(mean(unlist(lapply(boot, function(x) x[i]))),
                     quantile(unlist(lapply(boot, function(x) x[i])), probs=c(lci, uci), na.rm=T)))
      }
      names(out)<-c("Mean", "LCI", "UCI")
      rownames(out)<-names(boot[[1]])
      return(out)
    }

    boot<-list()
    for(i in 1:nboot){
      boot[[i]]<-apply(coef2[sample(nrow(coef2), nrow(coef2), replace=T, prob=coef2$Freq), 2:(ncol(coef2)-2) ],2, median, na.rm=T) #Modify
    }
    pop<-mean_ci_boot(boot)
    pop$Prop<-unlist(lapply(1:nrow(pop), function(x) ifelse(pop[x,1]>0, sum(coef2[,(x+1)]>0, na.rm=T),sum(coef2[,(x+1)]<0, na.rm=T))/length(mod_ls))) #Not useful since all are below zero
  }

  if(method=="murtaugh") {
    se<-lapply(1:i, function(x) data.frame(t(data.frame(mod_ls[[x]][[3]][[m]]))))
    se2<-plyr::rbind.fill(se)
    se2[se2<nudge]<-nudge
    cc<-coef2[,2:(ncol(coef2)-2)]

    ls<-lapply(1:ncol(se2), function(x) try(lm(cc[,x]~1, weights=1/se2[,x]^2)))
    ii<-which(unlist(lapply(ls, class))=="lm")
    pop<-data.frame(matrix(unlist(lapply(ls[ii], function(x) cbind(coef(x), confint(x, method="Wald")))), byrow=3, ncol=3))
    names(pop)<-c("Mean", "LCI", "UCI")
    rownames(pop)<-names(coef2[,2:(ncol(coef2)-2)])[ii]
    #pop$Prop<-unlist(lapply(1:nrow(pop), function(x) ifelse(pop[x,1]>0, sum(coef2[,(x+1)]>0, na.rm=T),sum(coef2[,(x+1)]<0, na.rm=T))/length(mod_ls))) #Not useful since all are below zero
  }
  return(list(pop, coef2))
}

#' Extract individual standard errors
#'
#' Extract individual standard errors
#' @param m model number (based on number in list of formula provided to rsf_ind)
#' @param mod_ls A list of list of model generated by rsf_ind
#' @param cutoff A cutoff value to exclude individuals with bad fit, default = -1 indicating model that did not converge will be excluded. Values > 0 will exclude based on coefficient
#' @param id_year Whether id_year (instead of individual) are provided. Individual and year needs to be separated by an underscore for the function to work properly.
#' @return A table of individual standard errors for each coefficients
#' @export
#' @examples
#' data(goats)
#' ls1<-list()
#' ls1[[1]]<-as.formula(STATUS~ELEVATION+SLOPE+ET+ASPECT+HLI+TASP)
#' ls1[[2]]<-as.formula(STATUS~ET+ASPECT+HLI+TASP)
#' out<-rsf_ind(goats$ID, data=goats, form_ls=ls1)
#' ind_se(m=1, out)
ind_se<-function(m=1, mod_ls,cutoff=0, id_year=F) {
  if(cutoff>0) {mod_ls<-rm_bad_fit1(mod_ls, cutoff=cutoff, m=m)}
  if(cutoff==-1) {mod_ls<-rm_conv_fit1(mod_ls,  m=m)}
  i<-length(mod_ls)
  add_weights<-function(x, name=name) {
    t<-data.frame(table(x$name))
    t$Freq<-1/t$Freq
    out<-merge(x, t, by.x="name", by.y="Var1")
    return(out)
  }
  se<-lapply(1:i, function(x) data.frame(t(data.frame(mod_ls[[x]][[3]][[m]]))))
  se2<-plyr::rbind.fill(se)
  rownames(se2)<-names(mod_ls)
  #se2$name<-matrix(unlist(strsplit(as.character(rownames(se2)), "_")), ncol=2, byrow=T)[,1]
  if(id_year==T) {se2$name<-matrix(unlist(strsplit(as.character(rownames(se2)), "_")), ncol=2, byrow=T)[,1]} ##### Needs to be adjusted when not Id year
  if(id_year==F) {se2$name<-rownames(se2)}

  se2$ID<-names(mod_ls)
  se2<-add_weights(se2)
  return(se2)
}

#' Extract individual coefficients
#'
#' Extract individual coefficients.
#' @param m model number (based on number in list of formula provided to rsf_ind)
#' @param mod_ls A list of list of model generated by rsf_ind
#' @param cutoff A cutoff value to exclude individuals with bad fit, default = -1 indicating model that did not converge will be excluded. Values > 0 will exclude based on coefficient
#' @param id_year Whether id_year (instead of individual) are provided. Individual and year needs to be separated by an underscore for the function to work properly.
#' @return A table of individual coefficients
#' @export
#' @examples
#' data(goats)
#' ls1<-list()
#' ls1[[1]]<-as.formula(STATUS~ELEVATION+SLOPE+ET+ASPECT+HLI+TASP)
#' ls1[[2]]<-as.formula(STATUS~ET+ASPECT+HLI+TASP)
#' out<-rsf_ind(goats$ID, data=goats, form_ls=ls1)
#' ind_coef(m=1, out)
ind_coef<-function(m=1, mod_ls,cutoff=0, id_year=F) {
  if(cutoff>0) {mod_ls<-rm_bad_fit1(mod_ls, cutoff=cutoff, m=m)}
  if(cutoff==-1) {mod_ls<-rm_conv_fit1(mod_ls,  m=m)}
  i<-length(mod_ls)
  add_weights<-function(x, name=name) {
    t<-data.frame(table(x$name))
    t$Freq<-1/t$Freq
    out<-merge(x, t, by.x="name", by.y="Var1")
    return(out)
  }
  co<-lapply(1:i, function(x) data.frame(t(data.frame(coef(mod_ls[[x]][[1]][[m]])))))
  coef2<-plyr::rbind.fill(co)
  rownames(coef2)<-names(mod_ls)
  if(id_year==T) {coef2$name<-matrix(unlist(strsplit(as.character(rownames(coef2)), "_")), ncol=2, byrow=T)[,1]} ##### Needs to be adjusted when not Id year
  if(id_year==F) {coef2$name<-rownames(coef2)}
  coef2$ID<-names(mod_ls)
  coef2<-add_weights(coef2)
  return(coef2)
}

#' Extraction of proximity from random forest classification
#'
#' Apply random forest classification
#' @param coef A matrix of model coefficient (from function ind_coef)
#' @return A proximity matrix
#' @export
#' @examples
#' data(goats)
#' ls1<-list()
#' ls1[[1]]<-as.formula(STATUS~ELEVATION+SLOPE+ET+ASPECT+HLI+TASP)
#' ls1[[2]]<-as.formula(STATUS~ET+ASPECT+HLI+TASP)
#' out<-rsf_ind(goats$ID, data=goats, form_ls=ls1)
#' coef<-ind_coef(m=1, out)
#' prox<-rf(coef)
rf<-function(coef, ntree=10000, ...) {
  out<-randomForest::randomForest(coef[,c(3:(ncol(coef)-2))], proximity=T, importance=T, ntree=ntree)
  prox<-out$proximity
  prox[lower.tri(prox, diag=T)]<-NA
  return(prox)
}


#' Evaluate ratio of used and random locations of individuals in a RSF table
#'
#' Evaluate ratio of used and random locations of individuals in a RSF table
#' @param id A vector of individual for each observation
#' @param value A vector indicating if each observation is used (=1) or random(=0)
#' @return A list indicating the range in ratio, range in random locations, and range in used location.
#' @export
eval_ratio<-function(id, value) {
  tt<-table(id, value)
  prop<-tt[,1]/tt[,2]
  out<-list(quantile(prop), quantile(tt[,1]), quantile(tt[,2]))
  names(out)<-c("Prop", "#Rnd", "#Used")
  return(out)
}


#' Resample a RSF table to keep constant ratio of used/random locations across individuals
#'
#' Resample a RSF table to keep constant ratio of used/random locations across individuals. Resampling is done with replacement.
#' @param data The RSF dataset to resample
#' @param id A vector of individual for each observation
#' @param value A vector indicating if each observation is used (=1) or random(=0)
#' @param ratio The ratio of random:used location (default =3, meaning 3 random locations for each used location)
#' @return A RSF dataset
#' @export
resample_rsf<-function(data, id="Id_Year", value="Value", ratio=3) {
  id1<-sort(unique(data[,id]))
  tt<-function(sub, value=value, ratio=ratio) {
    sub1<-sub[sub[,value]==1,]
    sub0<-sub[sub[,value]==0,]
    sub0a<-sub0[sample(1:nrow(sub0), nrow(sub1)*ratio, replace=T),]
    return(rbind(sub1,sub0a))
  }
  out<-data.frame()
  pb = txtProgressBar(min = 0, max = length(id1), initial = 0, style=3)
  for (i in 1:length(id1)) {
    sub<-data[data[,id]==id1[i],]
    sub2<-tt(sub, value=value, ratio)
    out<-rbind(out, sub2)
    setTxtProgressBar(pb,i)
  }
  return(out)
}

#' Perform kfold cross-validation at the individual level .
#'
#' Perform kfold cross-validation at the individual level and return histogram, mean kfold accros individual and min/max value
#' @param m model number (based on number in list of formula provided to rsf_ind)
#' @param mod_ls A list of list of model generated by rsf_ind
#' @param cutoff A cutoff value to exclude individuals with bad fit, default = -1 indicating model that did not converge will be excluded. Values > 0 will exclude based on coefficient
#' @param k number of fold (default = 5)
#' @param nrepet Number of repetitions (default =10)
#' @param nbins Number of bins (default =10)
#' @param jitter Logical, whether to add some random noise to the predictions (useful when the model is fitted on categorica  variables, which can produces error in the ranking process).
#' @param reproducible Logical, whether to use a fixed seed for each repetition.
#' @return A data frame with the correlations (\code{cor}) and the type of value (\code{type}).
#' @export
#' @examples
#' data(goats)
#' ls1<-list()
#' ls1[[1]]<-as.formula(STATUS~ELEVATION+SLOPE+ET+ASPECT+HLI+TASP)
#' ls1[[2]]<-as.formula(STATUS~ET+ASPECT+HLI+TASP)
#' out<-rsf_ind(goats$ID, data=goats, form_ls=ls1)
#' kfold_ind(m=1, out, ls=ls1)
kfold_ind<-function(m=1, mod_ls, ls=ls, cutoff=0, k=5, nrepet=5, nbins=10, grph=T) {
  if(cutoff>0) {mod_ls<-rm_bad_fit1(mod_ls, cutoff=cutoff, m=m)}
  if(cutoff==-1) {mod_ls<-rm_conv_fit1(mod_ls, m=m)}
  method=mod_ls[[1]][[1]][[1]]$method
  #x<-m
  #form_ls<-ls
  kk<-pbapply::pblapply(mod_ls, function(y) kfoldRSF(y[[1]][[m]], k=k, form_ls=ls, nrepet=nrepet, x=m, method=method, nbins=nbins, random=F))
  mean<-lapply(kk, mean, na.rm=T)
  if (grph) {plot(hist(unlist(mean)))}
  return(unlist(mean))
}

#' Perform kfold cross-validation on a RSF output.
#'
#' Perform kfold cross-validation on a RSF output. Similar to what is recommended in Boyce 2002. Function developped with Mathieu Basille
#' @param mod A RSF model (glm or glmer)
#' @param k number of fold (default = 5)
#' @param nrepet Number of repetitions (default =10)
#' @param nbins Number of bins (default =10)
#' @param jitter Logical, whether to add some random noise to the predictions (useful when the model is fitted on categorical variables, which can produces error in the ranking process).
#' @param reproducible Logical, whether to use a fixed seed for each repetition.
#' @return A data frame with the correlations (\code{cor}) and the type of value (\code{type}).
#' @export
kfoldRSF <- function(mod, k = 5, nrepet = 10, nbins = 10,  jitter = TRUE,
                     random = TRUE, method=method, x=m, form_ls=ls, reproducible = TRUE)
{
  if (!inherits(mod, c("glm", "mer", "glmerMod")))
    stop("Model of class '(g)lm' or '(g)lmer' expected")
  if (inherits(mod, c("glmerMod")))
    require(lme4)
  dt <- model.frame(mod)
  kfold <- rd <- numeric(length = nrepet)
  resp <- as.character(attr(terms(mod), "variables"))[attr(terms(mod),
                                                           "response") + 1]
  for (i in 1:nrepet) {
    dt$sets <- "train"
    if(reproducible)
      set.seed(i)
    dt$sets[sample(which(dt[, resp] == 1), sum(dt[, resp] ==
                                                 1)/k)] <- "test"
    reg <- update(mod, data = subset(dt, sets == "train"))
    cof<-coef(reg)
    cof[is.na(cof)]<-0
    if (inherits(mod, "glm"))
      predall <- exp(as.numeric(model.matrix(terms(reg),
                                             dt) %*% cof))
    else if (inherits(mod, "glmerMod"))
      predall <- exp(as.numeric(model.matrix(terms(reg),
                                             dt) %*% fixef(reg)))
    if (jitter) {
      if(reproducible)
        set.seed(i)
      predall <- jitter(predall)
    }
    quant <- quantile(predall[dt[, resp] == 0], probs = seq(from = 0,
                                                            to = 1, length.out = nbins + 1))
    quant[1] <- -Inf
    quant[length(quant)] <- Inf
    int <- factor(findInterval(predall[dt$sets == "test"],
                               quant), levels = 1:nbins)
    kfold[i] <- cor(1:nbins, table(int), method = "spearman")
    if (random) {
      if (reproducible)
        set.seed(i)
      dt$sets[sample(which(dt[, resp] == 0), sum(dt[, resp] ==
                                                   1)/k)] <- "rd"
      int <- factor(findInterval(predall[dt$sets == "rd"],
                                 quant), levels = 1:nbins)
      rd[i] <- cor(1:nbins, table(int), method = "spearman")
    }
  }
  if (random)
    return(data.frame(kfold = c(kfold, rd), type = rep(c("obs",
                                                         "rand"), each = nrepet)))
  else return(kfold)
}


#' goats - Mountain goats data set
#'
#' GPS collar data of mountain goats (Oreamnos americanus) from Lele and Keim (2006).
#'
#' @format A data frame with 19014 rows and 8 variables
#' @format STATUS a numeric vector, 1: used, 0: available
#' @format ID a numeric vector, individuals
#' @format ELEVATION a numeric vector (m)
#' @format SLOPE a numeric vector (degrees, steep)
#' @format ET a numeric vector, access to escape terrain (distance from steep slopes, m)
#' @format ASPECT a numeric vector (degrees)
#' @format HLI a numeric vector, heat load index (0-1)
#' @format TASP a numeric vector, transformed aspect
"goats"


#############################
#### Functions related to characterizing individual variation - new manuscript
##################################

#' Simulate normally-distributed individual coefficients from RSF/SSF based on their standard errors
#'
#' Simulate individual coefficients based on standard errors (to propagate uncertainty)
#' @param coef A matrix of individual coefficients (output of ind_coef)
#' @param coef A matrix of individual coefficients (output of ind_se)
#' @param n Number of random coefficients to generate for each individual (default=1000)
#' @return An array of individual coefficients
#' @export
#' @examples
#' data(goats)
#' ls1<-list()
#' ls1[[1]]<-as.formula(STATUS~ELEVATION+SLOPE+ET+ASPECT+HLI+TASP)
#' ls1[[2]]<-as.formula(STATUS~ET+ASPECT+HLI+TASP)
#' out<-rsf_ind(goats$ID, data=goats, form_ls=ls1)
#' coef<-ind_coef(m=1, out)
#' se<-ind_se(m1, out)
#' simu<-simu_coefs(coef, se, n=100)
simu_coefs<-function(coef, se, n=1000) {
  out<-array(NA, c(nrow(coef), ncol(coef), n))
  for (i in 1:nrow(coef)) {
    for (j in 2:(ncol(coef)-2)) {
      out[i,j,]<-rnorm(n,coef[i,j], se[i,j])
    }}
  return(out)
}

#' Un-weighted population average based on simulated coefficients
#'
#' Calculate population average for each replicate of simulated coefficients
#' @param simu An array of simulated individual coefficients based on their uncertainties(output of simu_coefs)
#' @return A matrix of population averages (one value for each covariates and replicates)
#' @export
#' @examples
#' data(goats)
#' ls1<-list()
#' ls1[[1]]<-as.formula(STATUS~ELEVATION+SLOPE+ET+ASPECT+HLI+TASP)
#' out<-rsf_ind(goats$ID, data=goats, form_ls=ls1)
#' coef<-ind_coef(m=1, out)
#' se<-ind_se(m=1, out)
#' simu<-simu_coefs(coef, se, n=100)
#' avg<-simu_avg(simu)
#' colnames(avg)<-names(coef)
#' head(avg)
#' apply(avg, 2, quantile, na.rm=T) #Show variation around estimate of each covariate
#' colMeans(avg) #Calculate average population average
simu_avg<-function(simu) {
 tt<-(matrix(unlist(lapply(1:dim(simu)[3], function(x) colMeans(simu[,,x], na.rm=T))), nrow=dim(simu)[3], ncol=dim(simu)[2], byrow=T))
}



#' Specialization based on simulated coefficients
#'
#' Calculate specialization for each replicate of simulated coefficients
#' @param simu An array of simulated individual coefficients based on their uncertainties(output of simu_coefs)
#' @return A matrix of specialization (one value for each covariates and replicates)
#' @export
#' @examples
#' data(goats)
#' ls1<-list()
#' ls1[[1]]<-as.formula(STATUS~ELEVATION+SLOPE+ET+ASPECT+HLI+TASP)
#' out<-rsf_ind(goats$ID, data=goats, form_ls=ls1)
#' coef<-ind_coef(m=1, out)
#' se<-ind_se(m=1, out)
#' simu<-simu_coefs(coef, se, n=100)
#' spe<-simu_spe(simu)
#' colnames(spe)<-names(coef)
#' head(spe)
#' apply(spe, 2, quantile, na.rm=T) #Show variation around estimate of each covariate
#' colMeans(spe) #Calculate average specialization for each covariate
simu_spe<-function(simu) {
  return(matrix(unlist(lapply(1:dim(simu)[3], function(x) apply(simu[,,x], 2, function(y) mean(abs(y), na.rm=T)))), nrow=dim(simu)[3], ncol=dim(simu)[2], byrow=T))
}

#' Individual variation (Heterogeneity) based on simulated coefficients
#'
#' Calculate heterogeneity for each replicate of simulated coefficients
#' @param simu An array of simulated individual coefficients based on their uncertainties(output of simu_coefs)
#' @return A matrix of heterogeneity (one value for each covariates and replicates)
#' @export
#' @examples
#' data(goats)
#' ls1<-list()
#' ls1[[1]]<-as.formula(STATUS~ELEVATION+SLOPE+ET+ASPECT+HLI+TASP)
#' out<-rsf_ind(goats$ID, data=goats, form_ls=ls1)
#' coef<-ind_coef(m=1, out)
#' se<-ind_se(m=1, out)
#' simu<-simu_coefs(coef, se, n=100)
#' sd<-simu_sd(simu)
#' colnames(sd)<-names(coef)
#' head(sd)
#' apply(sd, 2, quantile, na.rm=T) #Show variation around estimate of each covariate
#' colMeans(sd) #Calculate average heterogeneity for each covariate
simu_sd<-function(simu) {
  return(matrix(unlist(lapply(1:dim(simu)[3], function(x) apply(simu[,,x], 2, function(y) sd(y, na.rm=T)))), nrow=dim(simu)[3], ncol=dim(simu)[2], byrow=T))
}



#' Temporal consistency (with TWO time periods) based on simulated coefficients
#'
#' Calculate temporal consistency for each replicate of simulated coefficients
#' @param simu An array of simulated individual coefficients based on their uncertainties(output of simu_coefs)
#' @param col1 Column of the variable of interest for the 1st temporal period
#' @param col2 Column of the variable of interest for the 2nd temporal period
#' @return A matrix of consistency (one value for each covariates and replicates)
#' @export
#' @examples
#' data(goats)
#' goats$Season<-c("1", "2")
#' ls1<-list()
#' ls1[[1]]<-as.formula(STATUS~(ELEVATION+SLOPE+ET+ASPECT+HLI+TASP):Season)
#' out<-rsf_ind(goats$ID, data=goats, form_ls=ls1)
#' coef<-ind_coef(m=1, out)
#' se<-ind_se(m=1, out)
#' simu<-simu_coefs(coef, se, n=100)
#' head(coef)
#' cons_elevation<-simu_cons2(simu, 3, 4) #Calculate specialization for elevation covariate
#' quantile(cons_elevation) #Show variation around estimate of elevation covariate
#' mean(cons_elevation) #Calculate average consistency for elevation covariate
simu_cons2<-function(simu, col1, col2) {
  return(unlist(lapply(1:dim(simu)[3], function(x)   mean(abs(simu[,col1,x]-simu[,col2,x]), na.rm=T)  )))
}

#' Temporal consistency (with THREE time periods) based on simulated coefficients
#'
#' Calculate temporal consistency for each replicate of simulated coefficients
#' @param simu An array of simulated individual coefficients based on their uncertainties(output of simu_coefs)
#' @param col1 Column of the variable of interest for the 1st temporal period
#' @param col2 Column of the variable of interest for the 2nd temporal period
#' @param col3 Column of the variable of interest for the 3rd temporal period
#' @return A matrix of consistency (one value for each covariates and replicates)
#' @export
#' @examples
#' data(goats)
#' goats$Season<-c("1", "2", "3")
#' ls1<-list()
#' ls1[[1]]<-as.formula(STATUS~(ELEVATION+SLOPE+ET+ASPECT+HLI+TASP):Season)
#' out<-rsf_ind(goats$ID, data=goats, form_ls=ls1)
#' coef<-ind_coef(m=1, out)
#' se<-ind_se(m=1, out)
#' simu<-simu_coefs(coef, se, n=100)
#' head(coef)
#' cons_elevation<-simu_cons3(simu, 3, 4,5) #Calculate specialization for elevation covariate
#' quantile(cons_elevation) #Show variation around estimate of elevation covariate
#' mean(cons_elevation) #Calculate average consistency for elevation covariate
simu_cons3<-function(simu, col1, col2, col3) {
  return(unlist(lapply(1:dim(simu)[3], function(x)   mean(abs(simu[,col1,x]-simu[,col2,x])+abs(simu[,col1,x]-simu[,col3,x])+abs(simu[,col3,x]-simu[,col2,x])/3, na.rm=T) )))
}


#' Temporal reversal (with TWO time periods) based on simulated coefficients
#'
#' Calculate temporal reversal for each replicate of simulated coefficients
#' @param simu An array of simulated individual coefficients based on their uncertainties(output of simu_coefs)
#' @param col1 Column of the variable of interest for the 1st temporal period
#' @param col2 Column of the variable of interest for the 2nd temporal period
#' @return A matrix of reversal (one value for each covariates and replicates)
#' @export
#' @examples
#' data(goats)
#' goats$Season<-c("1", "2")
#' ls1<-list()
#' ls1[[1]]<-as.formula(STATUS~(ELEVATION+SLOPE+ET+ASPECT+HLI+TASP):Season)
#' out<-rsf_ind(goats$ID, data=goats, form_ls=ls1)
#' coef<-ind_coef(m=1, out)
#' se<-ind_se(m=1, out)
#' simu<-simu_coefs(coef, se, n=100)
#' head(coef)
#' rev_elevation<-simu_rev2(simu, 3, 4) #Calculate specialization for elevation covariate
#' quantile(rev_elevation) #Show variation around estimate of elevation covariate
#' mean(rev_elevation) #Calculate average reversal for elevation covariate
simu_rev2<-function(simu, col1, col2) {
  return(unlist(lapply(1:dim(simu)[3], function(x) mean(ifelse(sign(simu[,col1,x])!= sign(simu[,col2,x]), 1, 0), na.rm=T))))
}



#' Temporal reversal (with THREE time periods) based on simulated coefficients
#'
#' Calculate temporal reversal for each replicate of simulated coefficients
#' @param simu An array of simulated individual coefficients based on their uncertainties(output of simu_coefs)
#' @param col1 Column of the variable of interest for the 1st temporal period
#' @param col2 Column of the variable of interest for the 2nd temporal period
#' @param col3 Column of the variable of interest for the 3rd temporal period
#' @return A matrix of reversal (one value for each covariates and replicates)
#' @export
#' @examples
#' data(goats)
#' goats$Season<-c("1", "2", "3")
#' ls1<-list()
#' ls1[[1]]<-as.formula(STATUS~(ELEVATION+SLOPE+ET+ASPECT+HLI+TASP):Season)
#' out<-rsf_ind(goats$ID, data=goats, form_ls=ls1)
#' coef<-ind_coef(m=1, out)
#' se<-ind_se(m=1, out)
#' simu<-simu_coefs(coef, se, n=100)
#' head(coef)
#' rev_elevation<-simu_rev3(simu, 3, 4,5) #Calculate specialization for elevation covariate
#' quantile(rev_elevation) #Show variation around estimate of elevation covariate
#' mean(rev_elevation) #Calculate average reversal for elevation covariate
simu_rev3<-function(simu, col1, col2, col3) {
  return(unlist(lapply(1:dim(simu)[3], function(x) mean(c(ifelse(sign(simu[,col1,x])!= sign(simu[,col2,x]), 1, 0), ifelse(sign(simu[,col1,x])== sign(simu[,col3,x]), 1, 0), ifelse(sign(simu[,col3,x])== sign(simu[,col2,x]), 1, 0)),na.rm=T))))
}


cv<-function(x, na.rm=T) {sd(x, na.rm=na.rm)/mean(x,na.rm=na.rm)}

simu_cv<-function(simu) {
  return(matrix(unlist(lapply(1:dim(simu)[3], function(x) apply(simu[,,x], 2, function(y) cv(y, na.rm=T)))), nrow=dim(simu)[3], ncol=dim(simu)[2], byrow=T))
}



#
# rsf_mod2<-function(sub, form_ls, cleanModel=F) {
#   out<-lapply(1:length(form_ls), function(x) glm(form_ls[[x]], data=sub, family=binomial))
#   if (cleanModel) (out<-lapply(out, cleanModel))
#   return(out)
# }
#
#
# rsf_ind2<-function(id,data, form_ls, cleanModel=F) { #id is a vector of nrow(data)
#   if(length(id) != nrow(data)) (stop("id should be the same length as data"))
#   id1<-sort(unique(id))
#   out<-pbapply::pblapply(1:length(id1), function(x) rsf_mod2(sub=data[id==id1[x],], form_ls=form_ls, cleanModel=cleanModel))
#   names(out)<-id1
#   return(out)
# }
# #coef2score
# coef2score<-function(popavg) {
#   out<-data.frame(matrix(unlist(lapply(3:(ncol(popavg)-2), function(x) exp(popavg[,2]+popavg[,x]))), nrow=nrow(popavg), ncol=length(3:(ncol(popavg)-2))))
#   names(out)<-names(popavg[,c(3:(ncol(popavg)-2))])
#   out$ID<-popavg[,"ID"]
#   out$Freq<-popavg[,"Freq"]
#   return(out)
# }
