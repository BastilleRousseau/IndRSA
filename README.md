


## IndRSA  ##

**IndRSA** is an R package for individual-level modelling of resource selection.

This is the development area for the package `IndRSA`, which provides a series of functions to analyze resource selection of animals at the individual level (instead of using mixed-effects models). 

*References*: Bastille-Rousseau, G., and G. Wittemyer. 2019. Leveraging multidimensional heterogeneity in resource selection to define movement tactics of animals. Ecology Letters 22:1417–1427.

For questions: gbr |at| siu.edu

## Installation of the development version  ##

You need to use the package `devtools`from Hadley Wickham. 
    
    library(devtools)
    install_github("BastilleRousseau/IndRSA")


## Getting started ##

The package main functions are `rsf_mod`, `pop_avg`, `aictab_ind` and `kfold_ind`.  For a list of documented functions see the Reference manual. Examples of how to use the main functions are also provided in the vignette. 

## Important note ##
As of now, the functions are implemented to take an "ID_Year" as input for individuals in rsf_mod. The underscore "_" is critical for proper use of later functions. Even if your indivduals are over the same period or over multiple years, adding _XXXX is required. This will be fixed in future update.   
