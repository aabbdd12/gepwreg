# gepwreg: Percentile Weights Regression for Stata

This repository hosts the official Stata package for **Percentile Weights Regression (PWR)**, developed by Araar Abdelkrim.

## Description
The `gepwreg` command estimates the local derivative of an outcome variable with respect to covariates at a specific percentile of the outcome (or another ranking variable) distribution. 

### Key Features:
* **Analytical Standard Errors:** Closed-form linearisation-based standard errors (functional delta method) that account for the sampling variability of the estimated kernel weights.
* **MSE-Optimal Bandwidth:** Includes a practical two-step plug-in estimator for optimal bandwidth selection.
* **Complex Survey Design:** Fully supports survey data (strata, PSU) via Taylor linearisation (`vce(svy)`).
* **Generalised Ranking:** Allows ranking observations on any variable $z \neq y$.

## Installation
You can install the latest version directly from this GitHub repository by typing the following command in Stata:

```stata
net install gepwreg, from("[https://raw.githubusercontent.com/Araar-Abdelkrim/gepwreg/main](https://raw.githubusercontent.com/Araar-Abdelkrim/gepwreg/main)") replace