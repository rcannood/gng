#' Produce a FlowSOM-like plot
#'
#' @param gng_fit The GNG produced by the \code{\link{gng}} function
#' @param labels Labels for the training samples. Can also be NULL, if need be.
#' @param max.size The maximum size of visualised nodes
#'
#' @importFrom igraph graph_from_data_frame layout_with_kk
#' @importFrom ggforce geom_circle geom_arc_bar
#' @importFrom RColorBrewer brewer.pal
#' @importFrom reshape2 melt
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 ggplot geom_segment scale_fill_identity coord_equal labs
#'
#' @export
flowsomlike_plot <- function(gng_fit, labels, max.size = .075) {
  nodes <- gng_fit$nodes %>% mutate(name = as.character(name))
  edges <- gng_fit$edges %>% mutate(i = as.character(i), j = as.character(j))
  node_expr <- gng_fit$node_space
  clustering <- gng_fit$clustering

  gr <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes$name)

  # apply dimred to graph
  lay <- igraph::layout_with_kk(gr)
  lay <- apply(lay, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  colnames(lay) <- c("X", "Y")
  rownames(lay) <- nodes$name
  lay_df <- data.frame(nodes, lay, stringsAsFactors = FALSE, row.names = NULL)

  # combine edges with dimred
  gr_df_with_pos <- data.frame(
    edges,
    i = lay[edges$i, , drop = F],
    j = lay[edges$j, , drop = F],
    row.names = NULL
  )

  # make a colour scheme
  annotation_colours <- list(
    expr = set_names(RColorBrewer::brewer.pal(ncol(node_expr), "Dark2"), colnames(node_expr))
  )

  # scale expression between 0 and 1
  node_expr_sc <- apply(node_expr, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  node_expr_df <- node_expr_sc %>%
    reshape2::melt(varnames = c("name", "gene"), value.name = "expr") %>%
    mutate(name = as.character(name)) %>%
    left_join(lay_df %>% mutate(name = as.character(name)), by = "name") %>%
    group_by(name, gene) %>%
    mutate(
      gene_index = match(gene, colnames(node_expr)),
      start = (gene_index - 1) / ncol(node_expr) * 2 * pi,
      end = gene_index / ncol(node_expr) * 2 * pi
    ) %>%
    ungroup() %>%
    mutate(colour = annotation_colours$expr[gene])

  arc_df <- node_expr_df %>% mutate(r0 = 0.5 * max_size, r = (.5 + expr / 2) * max_size) %>% select(node = name, X, Y, start, end, r0, r, colour)

  # if labels are provided
  if (!is.null(labels)) {

    # check how many of each label are in each node
    categories <- if (is.factor(labels)) levels(labels) else sort(unique(labels))
    category_counts <- crossing(category = categories, cluster = seq_len(nrow(nodes))) %>%
      rowwise() %>%
      mutate(number = sum(labels == category, clustering == cluster)) %>%
      ungroup()

    annotation_colours$count <- set_names(RColorBrewer::brewer.pal(length(categories), "Set1"), categories)

    # generate pie df with positioning
    pie_df <- data_frame(label = as.character(labels), name = as.character(clustering)) %>%
      group_by(name, label) %>%
      summarise(n = n()) %>%
      mutate(
        value = n / sum(n) * 2 * pi,
        end = cumsum(value),
        start = end - value
      ) %>%
      ungroup() %>%
      left_join(lay_df, by = "name") %>%
      mutate(colour = annotation_colours$count[label])

    # add to arc df
    arc_df <- bind_rows(
      arc_df,
      pie_df %>% mutate(r0 = 0, r = .5 * max_size) %>% select(node = name, X, Y, start, end, r0, r, colour)
    )
  } else {
    arc_df$r0 <- 0
  }

  # Make a line plot
  ggplot() +
    geom_segment(aes(x = i.X, xend = j.X, y = i.Y, yend = j.Y), gr_df_with_pos) +
    ggforce::geom_circle(aes(x0 = X, y0 = Y, r = max_size), fill = "white", lay_df) +
    ggforce::geom_arc_bar(aes(x0 = X, y0 = Y, r0 = r0, r = r, start = start, end = end, fill = colour), data = arc_df %>% filter(!(start == 0 & end == 2 * pi))) +
    ggforce::geom_circle(aes(x0 = X, y0 = Y, r = r, fill = colour), data = arc_df %>% filter((start == 0 & end == 2 * pi))) +
    scale_fill_identity() +
    cowplot::theme_cowplot() +
    coord_equal() +
    labs(x = "X", y = "Y")
}
