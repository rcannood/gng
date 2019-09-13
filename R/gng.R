#' Growing Neural Gas
#'
#' @references Fritzke, Bernd. "A growing neural gas network learns topologies."
#'   Advances in neural information processing systems 7 (1995): 625-632.
#'
#' @param x The input data. Must be a matrix!
#' @param max_iter The max number of iterations.
#' @param epsilon_b Move the winning node by epsilon_b times the distance.
#' @param epsilon_n Move the neighbours of the winning node by epsilon_n times the distance.
#' @param age_max Remove edges older than age_max.
#' @param max_nodes The maximum number of nodes.
#' @param lambda Insert new nodes every lambda iterations.
#' @param alpha The decay parameter for error when a node is added.
#' @param beta The decay parameter for error in every node every iteration.
#' @param verbose Will output progress if \code{TRUE}.
#'
#' @importFrom RANN nn2
#' @importFrom dynutils calculate_distance
#' @importFrom lmds cmdscale_landmarks
#' @export
#'
#' @examples
#' library(ggplot2)
#' x <- as.matrix(iris[,1:4])
#' gng_out <- gng(x, max_nodes = 10)
#' autoplot(gng_out)
gng <- function(
  x,
  max_iter = 20000,
  epsilon_b = .05,
  epsilon_n = .001,
  age_max = 200,
  max_nodes = 30,
  lambda = 200,
  alpha = .5,
  beta = .99,
  verbose = TRUE
) {
  gng_out <- gng_cpp(
    x = x,
    max_iterations = max_iter,
    epsilon_b = epsilon_b,
    epsilon_n = epsilon_n,
    age_max = age_max,
    max_nodes = max_nodes,
    lambda = lambda,
    alpha = alpha,
    beta = beta,
    verbose = verbose
  )

  # rename edges and nodes
  gng_out$nodes$name <- as.character(gng_out$nodes$name)
  gng_out$edges$i <- gng_out$nodes$name[match(gng_out$edges$i, gng_out$nodes$index)]
  gng_out$edges$j <- gng_out$nodes$name[match(gng_out$edges$j, gng_out$nodes$index)]
  gng_out$clustering <- gng_out$nodes$name[RANN::nn2(gng_out$node_space, x, k = 1)$nn.idx[,1]]
  names(gng_out$clustering) <- rownames(x)

  # project gng and points into 2D space
  xn <- gng_out$node_space
  x2 <- rbind(xn, x)
  dis <- dynutils::calculate_distance(xn, x2, method = "euclidean")
  attr(dis, "landmark_ix") <- seq_len(nrow(xn))
  lmds <- lmds::cmdscale_landmarks(dis, ndim = 2)
  colnames(lmds) <- paste0("comp_", seq_len(ncol(lmds)))
  gng_out$node_proj <- lmds[seq_len(nrow(xn)), , drop = FALSE]
  gng_out$space_proj <- lmds[-seq_len(nrow(xn)), , drop = FALSE]

  # return output
  class(gng_out) <- "gng"
  gng_out
}

#' @export
#' @importFrom ggplot2 ggplot geom_point aes_string geom_segment theme_bw autoplot
autoplot.gng <- function(object, ...) {
  # if (use_projected) {
    node_space <- object$node_proj
    sample_space <- object$space_proj
  # } else {
  #   node_space <- gng_out$node_space
  #   sample_space <- x
  # }

  nodes <-
    object$nodes %>%
    left_join(as.data.frame(node_space) %>% rownames_to_column("name"), by = "name")

  edges <-
    object$edges %>%
    left_join(as.data.frame(node_space) %>% rename_all(~paste0("from_", .)) %>% rownames_to_column("i"), by = "i") %>%
    left_join(as.data.frame(node_space) %>% rename_all(~paste0("to_", .)) %>% rownames_to_column("j"), by = "j")

  points <-
    as.data.frame(sample_space) %>%
    rownames_to_column("sample") %>%
    mutate(cluster = object$clustering)

  ggplot(mapping = aes_string("comp_1", "comp_2")) +
    geom_point(data = points, colour = "darkgray", alpha = .5) +
    geom_segment(aes_string(x = "from_comp_1", xend = "to_comp_1", y = "from_comp_2", yend = "to_comp_2"), edges) +
    geom_point(data = nodes, size = 5) +
    theme_bw()
}

