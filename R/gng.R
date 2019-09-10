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
#' @param verbose Will output progress if \code{TRUE}.
#'
#' @importFrom RANN nn2
#' @export
#'
#' @examples
#' ## example with noisy spiral
#' library(stats)
#' n=2000
#' t=runif(n)^.7*10
#' al=.15;bet=.5;
#' x1=bet*exp(al*t)*cos(t)+rnorm(n,0,.1)
#' y1=bet*exp(al*t)*sin(t)+rnorm(n,0,.1)
#'
#' x <- cbind(x1, y1)
#' gng_out <- gng(x)
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
  gng_out$nodes$name <- as.character(gng_out$nodes$name)
  gng_out$edges$i <- gng_out$nodes$name[match(gng_out$edges$i, gng_out$nodes$index)]
  gng_out$edges$j <- gng_out$nodes$name[match(gng_out$edges$j, gng_out$nodes$index)]
  gng_out$clustering <- gng_out$nodes$name[RANN::nn2(gng_out$node_space, x, k = 1)$nn.idx[,1]]
  gng_out
}

