gng
===

[![Build
Status](https://travis-ci.org/rcannood/gng.svg?branch=master)](https://travis-ci.org/rcannood/gng)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/rcannood/gng?branch=master&svg=true)](https://ci.appveyor.com/project/rcannood/gng)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/gng)](https://cran.r-project.org/package=gng)
[![Coverage
Status](https://codecov.io/gh/rcannood/gng/branch/master/graph/badge.svg)](https://codecov.io/gh/rcannood/gng?branch=master)

An implementation of the Growing Neural Gas algorithm in Rcpp.

Example
-------

You can run gng as follows:

    library(gng)
    library(ggplot2)
    data(iris)
    gng_fit <- gng(as.matrix(iris[,1:4]), max_nodes = 10)

And visualise it as follows:

    autoplot(gng_fit)

![](man/figures/README_plot-1.png)
