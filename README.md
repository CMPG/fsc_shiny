# FSC Shiny app

This Shiny application has been developed to facilitate the interpretation of `fastsimcoal` results.

## Running the app

It is available on [the CMPG website](http://cmpg.unibe.ch/shiny/FSC/).

Alternatively, one can run the app locally using `runGithub()`:

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

## About the app

This Shiny application has been developed by Alexandre Gouy, based on R scripts written by Laurent Excoffier, Vitor Sousa and Alexandre Gouy from the CMPG lab at University of Bern.