#' Produce a FlowSOM-like plot
#'
#' @param gng_fit The GNG produced by the \code{\link{gng}} function
#' @param plot_labels Labels for the training samples. Is NULL, by default.
#' @param plot_expression Whether or not to plot the expression.
#' @param max_size The maximum size of visualised nodes
#'
#' @importFrom igraph graph_from_data_frame layout_with_kk
#' @importFrom ggforce geom_circle geom_arc_bar
#' @importFrom RColorBrewer brewer.pal
#' @importFrom reshape2 melt
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 ggplot geom_segment scale_fill_identity coord_equal labs aes geom_text
#'
#' @export
#'
#' @examples
#' data(iris)
#' gng_out <- gng(
#'   as.matrix(iris[,1:4]),
#'   max_iter = 20000,
#'   epsilon_b = 0.05,
#'   epsilon_n = 0.001,
#'   age_max = 200,
#'   max_nodes = 30,
#'   lambda = 200,
#'   alpha = 0.5,
#'   beta = 0.99,
#'   assign_cluster = TRUE,
#'   verbose = TRUE,
#'   make_logs = FALSE,
#'   cpp = TRUE
#' )
#' plot_gng(gng_out, iris[,5], max_size = 0.075)
#'
plot_gng <- function(
  gng_fit,
  plot_labels = NULL,
  plot_expression = gng_fit$node_space,
  max_size = .075
) {
  max_size_legend <- .2

  nodes <- gng_fit$nodes %>% mutate(name = as.character(name))
  edges <- gng_fit$edges %>% mutate(i = as.character(i), j = as.character(j))

  gr <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes$name)

  # apply dimred to graph
  lay <- apply(igraph::layout_with_kk(gr), 2, function(x) (x - min(x)) / (max(x) - min(x)))
  colnames(lay) <- c("X", "Y")
  rownames(lay) <- nodes$name
  lay_df <- data.frame(nodes, lay, r = max_size, stringsAsFactors = FALSE, row.names = NULL)

  # combine edges with dimred
  gr_df_with_pos <- data.frame(
    edges,
    i = lay[edges$i, , drop = F],
    j = lay[edges$j, , drop = F],
    row.names = NULL
  )

  # make a colour scheme
  annotation_colours <- list()
  arc_df <- NULL

  do_plot_expression <- !is.null(plot_expression) && !(is.logical(plot_expression) && !plot_expression)
  do_plot_labels <- !is.null(plot_labels) && !(is.logical(plot_labels) && !plot_labels)

  # If the user wants to plot the expression

  if (do_plot_expression) {
    annotation_colours$expr <- set_names(RColorBrewer::brewer.pal(ncol(plot_expression), "Dark2"), colnames(plot_expression))

    # scale expression between 0 and 1
    plot_expression_sc <- apply(plot_expression, 2, function(x) (x - min(x)) / (max(x) - min(x)))
    plot_expression_df <- plot_expression_sc %>%
      reshape2::melt(varnames = c("name", "gene"), value.name = "expr") %>%
      mutate(name = as.character(name)) %>%
      left_join(lay_df %>% mutate(name = as.character(name)), by = "name") %>%
      group_by(name, gene) %>%
      mutate(
        gene_index = match(gene, colnames(plot_expression)),
        start = (gene_index - 1) / ncol(plot_expression) * 2 * pi,
        end = gene_index / ncol(plot_expression) * 2 * pi
      ) %>%
      ungroup() %>%
      mutate(colour = annotation_colours$expr[gene])

    arc_df <-
      plot_expression_df %>%
      mutate(
        r0 = ifelse(is.null(do_plot_labels), 0, 0.5 * max_size),
        r = ifelse(is.null(do_plot_labels), expr, .5 + expr / 2) * max_size,
        plot_label = FALSE
      ) %>%
      select(node = name, X, Y, start, end, r0, r, colour, plot_label) %>%
      bind_rows(arc_df)

    # create legend plot
    num <- length(annotation_colours$expr)
    rads <- seq(0, 2 * pi, length.out = num + 1)

    leg_df <- data_frame(
      node = names(annotation_colours$expr),
      X = 1.4,
      Y = ifelse(is.null(do_plot_labels), .5, 0.75),
      start = rads %>% head(-1),
      end = rads %>% tail(-1),
      r0 = ifelse (is.null(do_plot_labels), 0, 0.5 * max_size_legend),
      r = { if (is.null(do_plot_labels)) seq(.5, 1, length.out = num) else seq(.75, 1, length.out = num) } * max_size_legend,
      colour = annotation_colours$expr,
      plot_label = TRUE
    )

    arc_df <- bind_rows(arc_df, leg_df)
    lay_df <- lay_df %>% add_row(name = "Expression", X = 1.4, Y = ifelse(is.null(do_plot_labels), .5, 0.75), r = max_size_legend)
  }



  # if labels are provided
  if (!is.null(plot_labels)) {
    clustering <- gng_fit$clustering

    # check how many of each label are in each node
    categories <- if (is.factor(plot_labels)) levels(plot_labels) else sort(unique(plot_labels))
    category_counts <-
      crossing(category = categories, cluster = seq_len(nrow(nodes))) %>%
      rowwise() %>%
      mutate(number = sum(categories == category, clustering == cluster)) %>%
      ungroup()

    annotation_colours$count <- set_names(RColorBrewer::brewer.pal(length(categories), "Set1"), categories)

    # generate pie df with positioning
    pie_df <- data_frame(label = as.character(plot_labels), name = as.character(clustering)) %>%
      group_by(name, label) %>%
      summarise(n = n()) %>%
      mutate(
        value = n / sum(n) * 2 * pi,
        end = cumsum(value),
        start = end - value
      ) %>%
      ungroup() %>%
      left_join(lay_df %>% select(-r), by = "name") %>%
      mutate(colour = annotation_colours$count[label])

    # add to arc df
    arc_df <-
      pie_df %>%
      mutate(r0 = 0, r = ifelse(do_plot_expression, .5 * max_size, max_size), plot_label = FALSE) %>%
      select(node = name, X, Y, start, end, r0, r, colour, plot_label) %>%
      bind_rows(arc_df)

    # create legend plot
    num <- length(annotation_colours$count)
    rads <- seq(0, 2 * pi, length.out = num + 1)
    leg_df <- data_frame(
      node = names(annotation_colours$count),
      X = 1.4,
      Y = ifelse(is.null(do_plot_expression), .5, 0.25),
      start = rads %>% head(-1),
      end = rads %>% tail(-1),
      r0 = 0,
      r = ifelse(is.null(do_plot_expression), 1, .5) * max_size_legend,
      colour = annotation_colours$count,
      plot_label = TRUE
    )

    arc_df <- bind_rows(arc_df, leg_df)
    lay_df <- lay_df %>% add_row(name = "Expression", X = 1.4, Y = ifelse(is.null(do_plot_labels), .5, 0.25), r = max_size_legend)
  }

  # Make a line plot
  label_df <- arc_df %>%
    filter(plot_label) %>%
    mutate(
      rad = (start + end) / 2,
      xpos = X + max_size_legend * 1.2 * sin(rad),
      ypos = Y + max_size_legend * 1.2 * cos(rad)
    )
  ggplot() +
    geom_segment(aes(x = i.X, xend = j.X, y = i.Y, yend = j.Y), gr_df_with_pos) +
    ggforce::geom_circle(aes(x0 = X, y0 = Y, r = r), fill = "white", lay_df) +
    ggforce::geom_arc_bar(aes(x0 = X, y0 = Y, r0 = r0, r = r, start = start, end = end, fill = colour), data = arc_df %>% filter(!(start == 0 & end == 2 * pi))) +
    ggforce::geom_circle(aes(x0 = X, y0 = Y, r = r, fill = colour), data = arc_df %>% filter((start == 0 & end == 2 * pi))) +
    geom_text(aes(xpos, ypos, label = node), label_df) +
    scale_fill_identity() +
    cowplot::theme_nothing() +
    coord_equal() +
    labs(x = NULL, y = NULL)
}
