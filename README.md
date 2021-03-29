# gng: Growing Neural Gas in RCpp

<!-- badges: start -->

[![R-CMD-check](https://github.com/rcannood/gng/workflows/R-CMD-check/badge.svg)](https://github.com/rcannood/gng/actions)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/gng)](https://cran.r-project.org/package=gng)
[![Coverage
Status](https://codecov.io/gh/rcannood/gng/branch/master/graph/badge.svg)](https://codecov.io/gh/rcannood/gng?branch=master)
<!-- badges: end -->

An implementation of the Growing Neural Gas algorithm in Rcpp.

## Example

Here’s an example of running a GNG on the iris dataset (which, arguably,
doesn’t make much sense).

    library(gng)
    data(iris)

    x <- as.matrix(iris[,1:4])
    gng_fit <- gng(x)

You can visualise the GNG nodes as follows.

    plot_gng(gng_fit, plot_labels = iris[,5])

![](man/figures/README_plot-1.png)

    plot_gng(gng_fit, plot_labels = iris[,5], plot_expression = NULL)

![](man/figures/README_plot-2.png)

    plot_gng(gng_fit)

![](man/figures/README_plot-3.png)
