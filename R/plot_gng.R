#' Produce a FlowSOM-like plot
#'
#' @param gng_fit The GNG produced by the \code{\link{gng}} function
#' @param plot_labels Labels for the training samples. Is NULL, by default.
#' @param plot_expression Whether or not to plot the expression.
#' @param max_size The maximum size of visualised nodes
#' @param max_size_legend The maximum size of the legend nodes.
#'
#' @importFrom igraph graph_from_data_frame layout_with_kk
#' @importFrom ggforce geom_circle geom_arc_bar theme_no_axes
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 ggplot geom_segment scale_fill_identity coord_equal labs aes geom_text
#' @importFrom utils head tail
#'
#' @export
#'
#' @importFrom dynutils scale_minmax
#'
#' @examples
#' data(iris)
#' gng_fit <- gng(x = as.matrix(iris[,1:4]))
#' plot_gng(gng_fit, plot_labels = iris[,5], max_size = 0.05)
#'
plot_gng <- function(
  gng_fit,
  plot_labels = NULL,
  plot_expression = gng_fit$node_space,
  max_size = .05,
  max_size_legend = .15
) {
  nodes <- gng_fit$nodes %>% mutate(name = as.character(.data$name))
  edges <- gng_fit$edges %>% mutate(i = as.character(.data$i), j = as.character(.data$j))

  gr <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes$name)

  # apply dimred to graph
  lay <- dynutils::scale_minmax(igraph::layout_with_kk(gr))
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
    plot_expression_sc <- dynutils::scale_minmax(plot_expression)
    plot_expression_df <- plot_expression_sc %>%
      as.data.frame() %>%
      rownames_to_column("name") %>%
      gather("gene", "expr", -.data$name) %>%
      left_join(lay_df %>% mutate(name = as.character(.data$name)), by = "name") %>%
      group_by(.data$name, .data$gene) %>%
      mutate(
        gene_index = match(.data$gene, colnames(plot_expression)),
        start = (.data$gene_index - 1) / ncol(plot_expression) * 2 * pi,
        end = .data$gene_index / ncol(plot_expression) * 2 * pi
      ) %>%
      ungroup() %>%
      mutate(colour = annotation_colours$expr[.data$gene])

    arc_df <-
      plot_expression_df %>%
      mutate(
        r0 = ifelse(is.null(do_plot_labels), 0, 0.5 * max_size),
        r = { if (is.null(do_plot_labels)) .data$expr else .5 + .data$expr / 2 } * max_size,
        plot_label = FALSE
      ) %>%
      select(node = .data$name, .data$X, .data$Y, .data$start, .data$end, .data$r0, .data$r, .data$colour, .data$plot_label) %>%
      bind_rows(arc_df)

    # create legend plot
    num <- length(annotation_colours$expr)
    rads <- seq(0, 2 * pi, length.out = num + 1)

    leg_df <- tibble(
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
  if (do_plot_labels) {
    clustering <- gng_fit$clustering

    # check how many of each label are in each node
    categories <- if (is.factor(plot_labels)) levels(plot_labels) else sort(unique(plot_labels))
    category_counts <-
      crossing(category = categories, cluster = seq_len(nrow(nodes))) %>%
      rowwise() %>%
      mutate(number = sum(categories == .data$category, clustering == .data$cluster)) %>%
      ungroup()

    annotation_colours$count <- set_names(RColorBrewer::brewer.pal(length(categories), "Set1"), categories)

    # generate pie df with positioning
    pie_df <- tibble(
      label = as.character(plot_labels),
      name = as.character(clustering)
    ) %>%
      group_by(.data$name, .data$label) %>%
      summarise(n = n(), .groups = "drop_last") %>%
      mutate(
        value = .data$n / sum(.data$n) * 2 * pi,
        end = cumsum(.data$value),
        start = .data$end - .data$value
      ) %>%
      ungroup() %>%
      left_join(lay_df %>% select(-.data$r), by = "name") %>%
      mutate(colour = annotation_colours$count[.data$label])

    # add to arc df
    arc_df <-
      pie_df %>%
      mutate(r0 = 0, r = ifelse(do_plot_expression, .5 * max_size, max_size), plot_label = FALSE) %>%
      select(node = .data$name, .data$X, .data$Y, .data$start, .data$end, .data$r0, .data$r, .data$colour, .data$plot_label) %>%
      bind_rows(arc_df)

    # create legend plot
    num <- length(annotation_colours$count)
    rads <- seq(0, 2 * pi, length.out = num + 1)
    leg_df <- tibble(
      node = names(annotation_colours$count),
      X = 1.4,
      Y = ifelse(is.null(do_plot_expression), .5, 0.25),
      start = rads %>% head(-1),
      end = rads %>% tail(-1),
      r0 = 0,
      r = ifelse(do_plot_expression, 1, .5) * max_size_legend,
      colour = annotation_colours$count,
      plot_label = TRUE
    )

    arc_df <- bind_rows(arc_df, leg_df)
    lay_df <- lay_df %>% add_row(name = "Expression", X = 1.4, Y = ifelse(is.null(do_plot_labels), .5, 0.25), r = max_size_legend)
  }

  if (!do_plot_labels && !do_plot_expression) {
    arc_df <- lay_df %>%
      rename(node = .data$name) %>%
      mutate(
        start = 0,
        end = 2 * pi,
        r0 = 0,
        colour = "lightgray",
        plot_label = FALSE
      ) %>%
      bind_rows(arc_df)
  }

  # Make a line plot
  label_df <- arc_df %>%
    filter(.data$plot_label) %>%
    mutate(
      rad = (.data$start + .data$end) / 2,
      xpos = .data$X + max_size_legend * 1.2 * sin(.data$rad),
      ypos = .data$Y + max_size_legend * 1.2 * cos(.data$rad)
    )
  ggplot() +
    geom_segment(aes(x = .data$i.X, xend = .data$j.X, y = .data$i.Y, yend = .data$j.Y), gr_df_with_pos) +
    ggforce::geom_circle(aes(x0 = .data$X, y0 = .data$Y, r = .data$r), fill = "white", lay_df) +
    ggforce::geom_arc_bar(
      aes(x0 = .data$X, y0 = .data$Y, r0 = .data$r0, r = .data$r, start = .data$start, end = .data$end, fill = .data$colour),
      data = arc_df %>% filter(!(.data$start == 0 & .data$end == 2 * pi))
    ) +
    ggforce::geom_circle(aes(x0 = .data$X, y0 = .data$Y, r = .data$r, fill = .data$colour), data = arc_df %>% filter((.data$start == 0 & .data$end == 2 * pi))) +
    geom_text(aes(.data$xpos, .data$ypos, label = .data$node), label_df) +
    scale_fill_identity() +
    ggforce::theme_no_axes() +
    coord_equal() +
    labs(x = NULL, y = NULL)
}
