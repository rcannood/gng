# gng: Growing Neural Gas in RCpp

<!-- badges: start -->

[![R-CMD-check](https://github.com/rcannood/gng/workflows/R-CMD-check/badge.svg)](https://github.com/rcannood/gng/actions)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/gng)](https://cran.r-project.org/package=gng)
[![Coverage
Status](https://codecov.io/gh/rcannood/gng/branch/master/graph/badge.svg)](https://codecov.io/gh/rcannood/gng?branch=master)
<!-- badges: end -->

An implementation of the Growing Neural Gas algorithm in Rcpp.

## Example

You can run gng as follows:

    library(gng)
    data(iris)

    x <- as.matrix(iris[,1:4])
    gng_fit <- gng(x)

And visualise it as follows:

    plot_gng(gng_fit, iris[,5], max_size = 0.05, max_size_legend = .15)

![](man/figures/README_plot-1.png)
