#' Produce a FlowSOM-like plot
#'
#' @param gng.fit The GNG produced by the \code{\link{gng}} function
#' @param labels Labels for the training samples
#' @param max.size The maximum size of visualised nodes
#'
#' @export
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr data_frame
#' @importFrom cowplot plot_grid
#' @importFrom ggtree subview
#'
#' @examples
#' # To be added soon!
flowsomlike_plot <- function(gng.fit, labels, max.size = .075) {
  library(ggplot2)

  nodes <- gng.fit$node.proj
  edges <- gng.fit$edges
  node.expr <- gng.fit$node.space

  categories <- if (is.factor(labels)) levels(labels) else sort(unique(labels))
  node.counts <-
    sapply(categories, function(ca) {
      sapply(seq_len(nrow(nodes)), function(i) {
        sum(labels == ca & gng.fit$clustering == i)
      })
    })


  annotation.colours <- list(
    count = setNames(RColorBrewer::brewer.pal(ncol(node.counts), "Set1"), colnames(node.counts)),
    expr = setNames(RColorBrewer::brewer.pal(ncol(node.expr), "Dark2"), colnames(node.expr))
  )

  # scale expression between 0 and 1
  node.expr.sc <- apply(node.expr, 2, function(x) (x - min(x)) / (max(x) - min(x)))

  # attach positions of edges
  colnames(nodes) <- c("X", "Y")
  colnames(edges) <- c("i", "j")
  edges.with.pos <- data.frame(
    edges,
    i = nodes[edges[,1],,drop=F],
    j = nodes[edges[,2],,drop=F]
  )

  make.plot <- function(x.count, x.expr, annotation.colours) {
    count.df <- data_frame(
      celltype = factor(names(x.count), levels = names(x.count)),
      count = x.count,
      percent = x.count / sum(x.count),
      right = cumsum(percent),
      left = right - percent
    )
    expr.df <- data_frame(
      marker = factor(names(x.expr), levels = names(x.expr)),
      expression = x.expr,
      right = as.integer(marker) / length(marker),
      left = right - min(right)
    )
    ggplot() +
      geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 2), fill = "white") +
      geom_rect(aes(xmin = left, xmax = right, ymin = 0, ymax = 1, fill = paste0("count.", celltype)), count.df) +
      geom_rect(aes(xmin = left, xmax = right, ymin = 1, ymax = 1 + expression, fill = paste0("expr.", marker)), expr.df) +
      geom_hline(yintercept = c(0, 1, 2)) +
      scale_fill_manual(values = unlist(annotation.colours)) +
      xlim(0, 1) +
      ylim(0, 2) +
      theme_nothing() +
      coord_polar()
  }
  node.plots <- lapply(seq_len(nrow(node.counts)), function(i) {
    x.count <- node.counts[i,]
    x.expr <- node.expr.sc[i,]
    make.plot(x.count, x.expr, annotation.colours)
  })
  make.legend.plot <- function(annotation.colours) {
    num.count <- length(annotation.colours$count)
    x.count <- setNames(rep(1, num.count), names(annotation.colours$count))
    num.expr <- length(annotation.colours$expr)
    x.expr <- setNames(seq(.5, 1, length.out = num.expr), names(annotation.colours$expr))

    count.df <- data_frame(
      celltype = factor(names(x.count), levels = names(x.count)),
      count = x.count,
      percent = x.count / sum(x.count),
      right = cumsum(percent),
      left = right - percent
    )
    expr.df <- data_frame(
      marker = factor(names(x.expr), levels = names(x.expr)),
      expression = x.expr,
      right = as.integer(marker) / length(marker),
      left = right - min(right)
    )

    g1 <- ggplot() +
      geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 2), fill = "white") +
      geom_rect(aes(xmin = left, xmax = right, ymin = 0, ymax = 1, fill = paste0("count.", celltype)), count.df) +
      geom_text(aes((left + right)/2, 2.5, label = celltype, colour = celltype), count.df) +
      geom_hline(yintercept = c(0, 1, 2)) +
      scale_fill_manual(values = unlist(annotation.colours)) +
      scale_colour_manual(values = annotation.colours$count) +
      xlim(0, 1) +
      ylim(0, 3) +
      theme_nothing() +
      coord_polar()

    g2 <- ggplot() +
      geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 2), fill = "white") +
      geom_rect(aes(xmin = left, xmax = right, ymin = 1, ymax = 1 + expression, fill = paste0("expr.", marker)), expr.df) +
      geom_text(aes((left + right)/2, 2.5, label = marker, colour = marker), expr.df) +
      geom_hline(yintercept = c(0, 1, 2)) +
      scale_fill_manual(values = unlist(annotation.colours)) +
      scale_colour_manual(values = annotation.colours$expr) +
      xlim(0, 1) +
      ylim(0, 3) +
      theme_nothing() +
      coord_polar()

    cowplot::plot_grid(g1, g2, nrow = 1)
  }

  p <- ggplot(data.frame(nodes), aes(X, Y)) +
    geom_segment(aes(x = i.X, xend = j.X, y = i.Y, yend = j.Y), data.frame(edges.with.pos)) +
    theme_nothing()
  for (i in seq_along(node.plots)) {
    x <- nodes[i,"X"]
    y <- nodes[i,"Y"]
    subplot <- node.plots[[i]]
    area <- gng.fit$nodes$number[[i]] / max(gng.fit$nodes$number) * pi/4*max.size^2
    width <- sqrt(area) * 4 / pi

    p <- ggtree::subview(p, subplot, x, y, width, width)
  }
  legends <- make.legend.plot(annotation.colours)
  ggtree::subview(p, legends, .15, .1, .2, .1)
}
