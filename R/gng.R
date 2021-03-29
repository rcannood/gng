#' Growing Neural Gas
#'
#' @references Fritzke, Bernd. "A growing neural gas network learns topologies." Advances in neural information processing systems 7 (1995): 625-632.
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
#' @param assign_cluster Whether or not to assign each sample to a GNG node.
#' @param verbose Will output progress if `TRUE`.
#' @param cpp Whether or not to use the C++ implementation over the R implementation. The C++ implementation is a lot faster.
#' @param make_logs_at At which iterations to store the GNG for visualisation purposes, only usable when `cpp == FALSE`.
#'
#' @export
#'
#' @importFrom FNN knnx.index
#'
#' @examples
#' data(iris)
#' gng_fit <- gng(x = as.matrix(iris[,1:4]))
#' plot_gng(gng_fit, plot_labels = iris[,5], max_size = 0.05)
#'
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
  assign_cluster = TRUE,
  verbose = TRUE,
  cpp = TRUE,
  make_logs_at = NULL
) {
  if (cpp) {
    if (!is.null(make_logs_at)) stop("make_logs_at is only supported using the R implementaton of GNG")

    o <- gng_cpp(
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
    o$nodes$name <- as.character(o$nodes$name)
    o$edges$i <- o$nodes$name[match(o$edges$i, o$nodes$index)]
    o$edges$j <- o$nodes$name[match(o$edges$j, o$nodes$index)]
  } else {
    o <- gng_r(
      x = x,
      max_iter = max_iter,
      epsilon_b = epsilon_b,
      epsilon_n = epsilon_n,
      age_max = age_max,
      max_nodes = max_nodes,
      lambda = lambda,
      alpha = alpha,
      beta = beta,
      verbose = verbose,
      make_logs_at = make_logs_at
    )
  }
  # Determine assignment of the samples to the nodes
  if (assign_cluster) {
    node_ix <- FNN::knnx.index(o$node_space, x, k = 1)[,1]
    o$clustering <- o$nodes$name[node_ix]
  }
  o
}
