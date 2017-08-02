#' Produce a FlowSOM-like plot
#'
#' @param gng_fit The GNG produced by the \code{\link{gng}} function
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
flowsomlike_plot <- function(gng_fit, labels, max.size = .075) {
  library(ggplot2)

  nodes <- gng_fit$node_proj
  edges <- gng_fit$edges
  node_expr <- gng_fit$node_space

  categories <- if (is.factor(labels)) levels(labels) else sort(unique(labels))
  node_counts <-
    sapply(categories, function(ca) {
      sapply(seq_len(nrow(nodes)), function(i) {
        sum(labels == ca & gng_fit$clustering == i)
      })
    })


  annotation_colours <- list(
    count = setNames(RColorBrewer::brewer.pal(ncol(node_counts), "Set1"), colnames(node_counts)),
    expr = setNames(RColorBrewer::brewer.pal(ncol(node_expr), "Dark2"), colnames(node_expr))
  )

  # scale expression between 0 and 1
  node_expr_sc <- apply(node_expr, 2, function(x) (x - min(x)) / (max(x) - min(x)))

  # attach positions of edges
  colnames(nodes) <- c("X", "Y")
  colnames(edges) <- c("i", "j")
  edges_with_pos <- data.frame(
    edges,
    i = nodes[edges[,1],,drop=F],
    j = nodes[edges[,2],,drop=F]
  )

  make_plot <- function(x_count, x_expr, annotation_colours) {
    count_df <- data_frame(
      celltype = factor(names(x_count), levels = names(x_count)),
      count = x_count,
      percent = x_count / sum(x_count),
      right = cumsum(percent),
      left = right - percent
    )
    expr_df <- data_frame(
      marker = factor(names(x_expr), levels = names(x_expr)),
      expression = x_expr,
      right = as.integer(marker) / length(marker),
      left = right - min(right)
    )
    ggplot() +
      geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 2), fill = "white") +
      geom_rect(aes(xmin = left, xmax = right, ymin = 0, ymax = 1, fill = paste0("count_", celltype)), count_df) +
      geom_rect(aes(xmin = left, xmax = right, ymin = 1, ymax = 1 + expression, fill = paste0("expr_", marker)), expr_df) +
      geom_hline(yintercept = c(0, 1, 2)) +
      scale_fill_manual(values = unlist(annotation_colours)) +
      xlim(0, 1) +
      ylim(0, 2) +
      theme_nothing() +
      coord_polar()
  }
  node_plots <- lapply(seq_len(nrow(node_counts)), function(i) {
    x_count <- node_counts[i,]
    x_expr <- node_expr_sc[i,]
    make_plot(x_count, x_expr, annotation_colours)
  })
  make_legend_plot <- function(annotation_colours) {
    num_count <- length(annotation_colours$count)
    x_count <- setNames(rep(1, num_count), names(annotation_colours$count))
    num_expr <- length(annotation_colours$expr)
    x_expr <- setNames(seq(.5, 1, length.out = num_expr), names(annotation_colours$expr))

    count_df <- data_frame(
      celltype = factor(names(x_count), levels = names(x_count)),
      count = x_count,
      percent = x_count / sum(x_count),
      right = cumsum(percent),
      left = right - percent
    )
    expr_df <- data_frame(
      marker = factor(names(x_expr), levels = names(x_expr)),
      expression = x_expr,
      right = as.integer(marker) / length(marker),
      left = right - min(right)
    )

    g1 <- ggplot() +
      geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 2), fill = "white") +
      geom_rect(aes(xmin = left, xmax = right, ymin = 0, ymax = 1, fill = paste0("count_", celltype)), count_df) +
      geom_text(aes((left + right)/2, 2.5, label = celltype, colour = celltype), count_df) +
      geom_hline(yintercept = c(0, 1, 2)) +
      scale_fill_manual(values = unlist(annotation_colours)) +
      scale_colour_manual(values = annotation_colours$count) +
      xlim(0, 1) +
      ylim(0, 3) +
      theme_nothing() +
      coord_polar()

    g2 <- ggplot() +
      geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 2), fill = "white") +
      geom_rect(aes(xmin = left, xmax = right, ymin = 1, ymax = 1 + expression, fill = paste0("expr_", marker)), expr_df) +
      geom_text(aes((left + right)/2, 2.5, label = marker, colour = marker), expr_df) +
      geom_hline(yintercept = c(0, 1, 2)) +
      scale_fill_manual(values = unlist(annotation_colours)) +
      scale_colour_manual(values = annotation_colours$expr) +
      xlim(0, 1) +
      ylim(0, 3) +
      theme_nothing() +
      coord_polar()

    cowplot::plot_grid(g1, g2, nrow = 1)
  }

  p <- ggplot(data.frame(nodes), aes(X, Y)) +
    geom_segment(aes(x = i.X, xend = j.X, y = i.Y, yend = j.Y), data.frame(edges_with_pos)) +
    theme_nothing()
  for (i in seq_along(node_plots)) {
    x <- nodes[i,"X"]
    y <- nodes[i,"Y"]
    subplot <- node_plots[[i]]
    area <- gng_fit$nodes$number[[i]] / max(gng_fit$nodes$number) * pi/4*max_size^2
    width <- sqrt(area) * 4 / pi

    p <- ggtree::subview(p, subplot, x, y, width, width)
  }
  legends <- make_legend_plot(annotation_colours)
  ggtree::subview(p, legends, .15, .1, .2, .1)
}
