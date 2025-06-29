
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nhgCI

<!-- badges: start -->
<!-- badges: end -->

The **nhgCI** package provides methods for constructing exact and
approximate confidence intervals for the parameters of the Negative
Hypergeometric distribution.

It supports the following methods:

- Analog to Clopper-Pearson

- Conditional Minimal Cardinality (CMC)

- Crow & Gardner (CG)

- Blaker

The package handles cases where either: - The number of successes (`M`)
is unknown, or - The total population size (`N`) is unknown.

## Installation

You can install the development version from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("rroggenk/nhgCI")
```


## Example

Here’s a basic example of using `nhgCI`:

``` r
library(nhgCI)

# Confidence intervals for unknown M
nhgCI_M(N = 30, m = 5, conf_level = 0.95, method = "Analog to Clopper-Pearson")

# Confidence intervals for unknown N
nhgCI_N(M = 10, m = 5, conf_level = 0.95, method = "CMC")
```

## Notes

- For `nhgCI_M()`, confidence intervals are computed for all appropriate
  values of observed failures `x` automatically.
- For `nhgCI_N()`, the number of observed failures `x` is theoretically
  unbounded.  
  You may need to increase the `max_N` parameter if your desired `x`
  value does not appear.

Some computations (especially for large sample sizes or high precision)
may take a while to complete.

[![DOI](https://zenodo.org/badge/974013857.svg)](https://doi.org/10.5281/zenodo.15617252)

