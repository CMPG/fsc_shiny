# FSC Shiny app

This Shiny application has been developed to facilitate the interpretation of `fastsimcoal` results.

## Running the application

The app can be found online on [the CMPG website](http://cmpg.unibe.ch/shiny/FSC/).

Alternatively, one can run the latest version of the app locally using `runGithub()`:

```R
library(shiny)
runGitHub("fsc_shiny", "CMPG")
```

Or, after cloning this git repository, using `runApp()`:

```R
# First clone the repository with git. If you have cloned it into
# ~/fsc_shiny, first go to that directory, then use runApp().

setwd("~/fsc_shiny")
runApp()
```

## Using the application

FSC takes as input a parameter file (.par) and the observed SFS (.obs) and outputs various files, including one or several files containing the expected SFS. Different files can be uploaded to be visualisedin the app:

* the parameter file (.par);
* observed and expected SFS (depending on the number of populations, one could consider marginal or joint SFS);
* likelihoods of different parameters values (for single or multiple runs).
  
A __single .par file__ can be uploaded to plot the demographic model. All other files should be uploaded as __a single ZIP archive__ containing the working directory used for running fastsimcoal.

Various example files are provided in the `example_files` subdirectory.

## About

This Shiny application has been developed by Alexandre Gouy, based on R scripts written by Laurent Excoffier, Vitor Sousa and Alexandre Gouy from the CMPG lab at University of Bern.