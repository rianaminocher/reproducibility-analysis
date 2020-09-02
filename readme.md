reproducibility-analysis
============

[![License: CC BY-NC-SA 4.0](https://licensebuttons.net/l/by-nc-sa/4.0/80x15.png)](https://creativecommons.org/licenses/by-nc-sa/4.0/)

This repository contains data and analysis code for:

> Minocher, et al. "Reproducibility of social learning research declines exponentially over 63 years of publication"

This code was written by Riana Minocher and Bret Beheim under Creative Commons License [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/). See LICENSE.md for details.

The code was developed on R v3.6.1. This code uses the Stan MCMC engine via the `rstan` package (v2.19.2), which requires a C++ compiler. Installation instructions are available at https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started. We also use the `rethinking` package (v1.93), installation instructions at http://xcelab.net/rm/software/.

## How to run

To run the primary analyses

1. In R, set the working directory to this folder (e.g. "reproducibility-analysis") containing the `README.md` and `analysis.r`.

2. Type `source("analysis.r")` and hit enter.

If the required packages are correctly installed, the code will take a few minutes to run and return the output assets to the `output/` folder. The code is designed to run in sequence on the data described by the glossary below.

## About the data

There are two data frames in this analysis. The `pubs_anon.csv` file contains data from part I of the study, on the availability of materials from all sampled publications (n = 560 papers). Each row corresponds to a single publication with the following data columns:

```

# key = unique identifier for each study
# year = year of publication of study
# data_type = type of data (study design) included in study
# species = study species including human or non-human subjects
# id = unique id for the first or corresponding author
# emailed = T/F whether the author was contacted about materials
# reason_no_email = why, if emailed is F
# downloaded = T/F whether (any) materials were obtained online
# reply_received = T/F whether author replied to request
# undelivered = T/F if email attempt was unsuccessful
# auto = T/F if email attempt returned auto response
# data_sent = T/F if materials were received from authors
# data_available = T/F whether downloaded OR data_sent

```

The second data frame is `results_anon.csv`, where each row corresponds to a single result (n = 111 results), from a publication in our subsample of Phase II reanalyses, with the following data columns:

```

# key = unique identifier for the study this result was drawn from
# reproduced = T/F whether result was successfully reproduced
# failure_category = reason for unsuccessful reproduction
# scripted = T/F if the study materials contained analytical code
# species = study species including human or non-human subjects
# type = type of data (study design) included in study
# year = year of publication of study

```
