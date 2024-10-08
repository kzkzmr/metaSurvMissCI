
<!-- README.md is generated from README.Rmd. Please edit that file -->

# metaSurvMissCI

<!-- badges: start -->
<!-- badges: end -->

This package include two functions: `impute_se_surv` and `forest_surv`.
`impute_se_surv` impute missing precision information for meta analysis
of survival rates and calculate survival rates and their SEs on the
transformed scale. `forest_surv` provides a forest plot of the
meta-analysis of survival rates.

For more information about each function, type `?impute_se_surv` and
`?forest_surv`.

## Installation

You can install the development version of metaSurvMissCI from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kzkzmr/metaSurvMissCI")
```

## Example

``` r
library(metaSurvMissCI)
library(meta)
```

``` r
data("metadata_chordoma")
impdata <- impute_se_surv(data = metadata_chordoma, St = "PFS5y",
                          LCL = "PFSL5y", UCL = "PFSU5y", n = "n",
                          nt = "n_5yPFS", ne = "ne_PFS", p = "pr_PFS")
meta_imp <- metagen(TE = tr_St, seTE = tr_SE, studlab = Study,
                     data = impdata, method.tau = "REML")
forest_surv(meta_imp, xlim = c(30, 100), estlab = "5-year PFS (%)")
```

<img src="man/figures/README-example-1.png" width="100%" />
