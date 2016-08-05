library(ggplot2)
library(cowplot)
library(dplyr)

set.seed(1)

#' Ik reduceer de ginhoux data naar een 10-dimensionale ruimte.
library(SCORPIUS)
data(ginhoux)
dist <- correlation.distance(ginhoux$expression)
X <- reduce.dimensionality(dist, ndim = 5)

gng.fit <- gng(X)

group.name <- ginhoux$sample.info$group.name

g1 <- flowsomlike_plot(gng.fit, group.name, max.size = .075)

space.df <- data.frame(gng.fit$space.proj, group.name)
g2 <- ggplot() + geom_point(aes(GNG_X, GNG_Y, colour = group.name), space.df)

cowplot::plot_grid(g1, g2, nrow = 1)
