reproducibility-analysis
============

[![License: CC BY-NC-SA 4.0](https://licensebuttons.net/l/by-nc-sa/4.0/80x15.png)](https://creativecommons.org/licenses/by-nc-sa/4.0/)

## About this repository

This repository contains data and code to reproduce results from

> Minocher, et al. "Estimating the reproducibility of social learning research published between 1955 and 2018"

This code was written by Riana Minocher under Creative Commons License [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/). See LICENSE.md for details.

The code was written in R 4.0.3. Statistical models are fit using the Stan MCMC engine via the `rstan` package (2.21.2), which requires a C++ compiler. Installation instructions are available at https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started. The `rethinking` package (2.12) is required to process fitted model outputs - installation instructions at http://xcelab.net/rm/software/.

## How to run

To fit the model and process the output

1. In R, set the working directory to this folder (i.e. "reproducibility-analysis") containing the `README.md` and `analysis.r`.

2. In the R console, execute the line `source("analysis.r")`.

If the required packages are correctly installed, the code will take a few minutes to run and create the manuscript figures and tables in the `output/` folder. 

The `model.stan` file contains the model code, and is called by the `analysis.R` script.

## The support folder

The `support` folder contains files that are supplementary to the main analysis.  

The `bib` folder contains a list of our full sample of publications. 

The script `pps.R` contains code to perform the prior predictive simulation, to verify reasonable behaviour of model priors.

The script `robustness-checks.R` performs a series of alternative analyses on the same data, to verify that our analysis choice does not bias our results. This scripts executes the models `robustness-checks-m1.stan`,  `robustness-checks-m2.stan` and  `robustness-checks-code-availability.stan`. 

The script `sample.R` contains information on our justification of the sample size for the subset.

These scripts produce figures and tables which are stored to the `output` folder. 

The `protocols` folder contains information we circulated amongst research team members on coding data, the template protocol for recording data and the data request template. 

## About the data

The table `anon_database.csv`, within the input folder, contains data on the 560 papers sampled in this study. Each row corresponds to a single publication, with the following data columns:

```

# key = unique identifier for a paper
# year = year of publication of paper
# type = type of data included in study, which is a composite of design (observational or experimental) and species (human or non-human)
# emailed = T/F whether the author was contacted about materials
# reason_no_email = why, if emailed is F
# downloaded = T/F whether (any) materials were obtained online
# reply_received = T/F whether author replied to request
# data_sent = T/F if materials were received from authors
# data_available = T/F whether downloaded OR data_sent
# n_results = number of results coded for this study, if sampled as part of second phase
# data_complete = number of results for which data was complete, or usable
# analysis_clear = number of results for which analysis steps were clear or executable
# reproduced = number of results which were successfully reproduced, i.e. corresponded to published results
# reproduced = number of results which were successfully reproduced, i.e. corresponded to published results
# scripted = 1/0 whether the study we evaluated had available code 


```
