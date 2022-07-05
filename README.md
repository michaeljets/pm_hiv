# Simulation example for "Optimizing Treatment for Human Immunodeficiency Virus to Improve Clinical Outcomes using Precision Medicine."

## Author: Michael Jetsupphasuk

This repository gives a simulated example of the methods used in the paper titled "Optimizing Treatment for Human Immunodeficiency Virus to Improve Clinical Outcomes using Precision Medicine." Since the data used in the study is non-public, an example using simulated data is given in lieu of the actual data and code used in the study analysis. 

Run the following R scripts in order to reproduce results:

1. `create_data.R`
2. `get_cv_treatrules.R`
3. `boot_cv.R`
4. `evaluate.R`

Point estimates may take 45 minutes to compute and bootstrapped standard errors may take 19 hours (using 200 bootstraps), depending on computing power.

The R script `RISTfunctions.r` was taken from https://sites.google.com/site/teazrq/software. 