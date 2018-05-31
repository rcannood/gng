#' Growing Neural Gas
#'
#' @references Fritzke, Bernd. "A growing neural gas network learns topologies." Advances in neural information processing systems 7 (1995): 625-632.
#'
#' @param x The input data
#' @param max_iter The max number of iterations
#' @param epsilon_b Move the winning node by epsilon_b times the distance
#' @param epsilon_n Move the neighbours of the winning node by epsilon_n times the distance
#' @param age_max Remove edges older than age_max
#' @param max_nodes The maximum number of nodes
#' @param lambda Insert new nodes every lambda iterations
#' @param alpha The decay parameter for error when a node is added
#' @param beta The decay parameter for error in every node every iteration
#' @param assign_cluster Whether or not to assign each sample to a GNG node
#' @param verbose Will output progress if \code{TRUE}
#' @param make_logs The decay for error of all nodes every iteration
#' @param cpp Whether or not to use the C++ implementation over the R implementation.
#'
#' @export
gng <- function(x, max_iter = 20000, epsilon_b = .05, epsilon_n = .001, age_max = 200, max_nodes = 30, lambda = 200, alpha = .5, beta = .99, assign_cluster = T, verbose = T, make_logs = F, cpp = T) {
  if (cpp) {
    o <- gng_cpp(x, max_iter, epsilon_b, epsilon_n, age_max, max_nodes, lambda, alpha, beta, verbose)
  } else {
    o <- gng_r(x, max_iter, epsilon_b, epsilon_n, age_max, max_nodes, lambda, alpha, beta, verbose, make_logs)
  }
  # Determine assignment of the samples to the nodes
  if (assign_cluster) {
    o$clustering <- apply(x, 1, function(xi) which.min(apply(o$node_space, 1, function(si) mean((xi-si)^2))))
  }
  o
}

#' @importFrom stats runif
gng_r <- function(x, max_iter, epsilon_b, epsilon_n, age_max, max_nodes, lambda, alpha, beta, verbose, make_logs) {
  # Calculate ranges of the dimensionality of the dataset
  ranges <- apply(x, 2, range)

  # Initialise a node position matrix
  S <- matrix(NA, nrow = max_nodes, ncol = ncol(x), dimnames = list(NULL, colnames(x)))

  # A dataframe containing the error of each node
  S_meta <- data.frame(i = seq_len(nrow(S)), error = rep(0, nrow(S)))

  # A list containing the neighbours of each node
  E <- lapply(seq_len(max_nodes), function(i) numeric(0))
  E[[1]] <- 2
  E[[2]] <- 1

  # A matrix containing the ages of each edge
  Age <- matrix(NA, nrow = max_nodes, ncol = max_nodes)

  # 0. start with two units a and b at random positions in R^n
  S[1:2,] <- apply(ranges, 2, function(r) {
    stats::runif(2, r[[1]], r[[2]])
  })
  Age[1,2] <- 0
  Age[2,1] <- 0

  # This variable is the index of a new node
  nextr <- 3

  # I separated out the distance and move functions in case I want to try a different approach
  distance_function <- function(xi, Si) {
    diff <- xi - Si
    sqrt(mean(diff * diff))
  }
  move_function <- function(xi, Si, epsilon) {
    Si + epsilon * (xi - Si)
  }
  sample_input_signal <- function() {
    x[sample.int(nrow(x), 1),]
  }

  if (make_logs) {
    node_loglist <- list(
      data.frame(added = T, node = c(1, 2), iteration = 0)
    )
    nodepos_loglist <- list(
      data.frame(iteration = 0, node = c(1, 2), S[1:2,])
    )
    edge_loglist <- list(
      data.frame(added = T, i = 1, j = 2, iteration = 0, reason = 1, age = 0)
    )
  }

  current.iter <- 0

  # while stopping criterion not met.
  # In the future, this function might check whether nodes have somewhat converged to a steady position.
  while (current.iter < max_iter) {
    current.iter <- current.iter + 1
    if (verbose && current.iter %% 1000 == 0) cat("Iteration ", current.iter, "\n", sep = "")

    # 1. generate input signal
    xi <- sample_input_signal()

    # 2. find nearest unit s1 and second nearest unit s2
    sdist <- apply(S, 1, function(Si) distance_function(xi, Si))
    sord <- order(sdist)
    s1 <- sord[[1]]
    s2 <- sord[[2]]

    # 3. increment the age of all edges eminating from s1
    Age[s1,] <- Age[s1,] + 1
    Age[,s1] <- Age[,s1] + 1

    # 4. add the squared distance between input signal and the nearest unit to the error variable
    S_meta$error[[s1]] <- S_meta$error[[s1]] + sdist[[s1]]

    # 5. move s1 and its direct topological neighbours towards input signal by fractions epsilon_b and epsilon_n respectively
    neighs <- E[[s1]]
    S[s1, ] <- move_function(xi, S[s1,], epsilon_b)
    for (n1 in neighs) {
      S[n1,] <- move_function(xi, S[n1,], epsilon_n)
    }

    if (make_logs) {
      nodepos_loglist[[length(nodepos_loglist)+1]] <- data.frame(iteration = current.iter, node = c(s1, neighs), S[c(s1, neighs), , drop = F])
    }

    # 6. if s1 and s2 are connected by an edge, set the age of this edge to zero. otherwise, create edge of age 0
    if (is.na(Age[[s1, s2]])) {
      E[[s1]] <- c(E[[s1]], s2)
      E[[s2]] <- c(E[[s2]], s1)
    }
    Age[[s1, s2]] <- 0
    Age[[s2, s1]] <- 0

    edge.nods <- unique(neighs, s2)

    if (make_logs) {
      edge_loglist[[length(edge_loglist)+1]] <- data.frame(added = Age[s1, edge.nods] <= age_max, i = s1, j = edge.nods, iteration = current.iter, reason = 2, age = Age[s1, edge.nods])
    }

    # 7. remove edges with an age larger than age_max. if a point has no remaining edge, remove as well
    rem <- which(Age[s1,] > age_max)
    if (length(rem) > 0) {
      Age[s1, rem] <- NA
      Age[rem, s1] <- NA
      E[[s1]] <- setdiff(E[[s1]], rem)
      for (sj in rem) {
        E[[sj]] <- setdiff(E[[sj]], s1)
      }
      # edge_loglist[[length(edge_loglist)+1]] <- data.frame(added = F, i = s1, j = rem, iteration = current.iter, reason = "age")
      s1rem <- c(s1, rem)
      removed.nodes <- s1rem[which(sapply(s1rem, function(i) length(E[s1rem]) == 0))]

      if (make_logs && length(removed.nodes) > 0) {
        node_loglist[[length(node_loglist)+1]] <- data.frame(added = F, node = removed.nodes, iteration = current.iter)
      }
    }

    # 8. If the number of input signals generated so far is an integer multiple of a parameter lambda, insert a new node
    if (current.iter %% lambda == 0) {
      if (nextr <= max_nodes) {
        r <- nextr
        nextr <- nextr + 1

        # 8a. determine node p with the maximum accumulated error
        p <- which.max(S_meta$error)

        # 8b. determine node q neighbouring p with the largest error
        np <- E[[p]]
        q <- np[which.max(S_meta$error[np])]

        # 8c. insert new node r halfway between p and q
        S[r,] <- (S[p,] + S[q,]) / 2

        if (make_logs) {
          nodepos_loglist[[length(nodepos_loglist)+1]] <- data.frame(iteration = current.iter, node = r, S[r, , drop = F])
        }

        # 8d. remove (p, q), add (p, r) and (q, r)
        Age[p, q] <- NA
        Age[q, p] <- NA
        E[[p]] <- c(setdiff(E[[p]], q), r)
        E[[q]] <- c(setdiff(E[[q]], p), r)
        E[[r]] <- c(p, q)
        Age[p, r] <- 0
        Age[r, p] <- 0
        Age[q, r] <- 0
        Age[r, q] <- 0

        # 8e. decrease error of p and q by factor alpha and initialise error of r as the error of p
        Ep <- S_meta$error[[p]]
        Eq <- S_meta$error[[q]]
        S_meta$error[c(p, q, r)] <- c(alpha * Ep, alpha * Eq, alpha * Ep)

        if (make_logs) {
          node_loglist[[length(node_loglist)+1]] <- data.frame(added = T, node = r, iteration = current.iter+.01)
          edge_loglist[[length(edge_loglist)+1]] <- data.frame(added = c(F, T, T), i = c(p, p, q), j = c(q, r, r), iteration = current.iter + .01, reason = 3, age = c(age_max + 1, 0, 0))
        }
      }
    }

    # 9. decrease all variables by multiplying them with a constant beta
    S_meta$error <- S_meta$error * beta
  }

  # Construct GNG output data structures
  nodes <- data.frame(node = seq_len(nrow(S)), S_meta[,2,drop=F])
  node_space <- S
  filt <- !is.na(node_space[,1])
  nodes <- nodes[filt,,drop=F]
  node_space <- node_space[filt,,drop=F]

  edges <- bind_rows(lapply(seq(2, nrow(S)), function(nj) {
    ni <- E[[nj]]
    ni <- ni[ni < nj]
    if (length(ni) > 0) {
      data.frame(i = ni, j = nj)
    } else {
      NULL
    }
  }))

  # Process the logs
  if (make_logs) {
    edge_log <- bind_rows(edge_loglist)
    node_log <- bind_rows(node_loglist)
    nodepos_log <- bind_rows(nodepos_loglist)
    edge_log[,2:3] <- t(apply(edge_log[,2:3], 1, sort))
  } else {
    edge_log <- NULL
    node_log <- NULL
    nodepos_log <- NULL
  }

  list(
    nodes = nodes,
    node_space = node_space,
    edges = edges,
    node_log = node_log,
    nodepos_log = nodepos_log,
    edge_log = edge_log
  )
}

#' Perform FR dimensionality reduction on GNG and project cells to same space
#'
#' @param gng_out The object generated by \code{\link{gng}}
#' @param make_projection Whether or not to also project the cells using a randomForest.
#'
#' @importFrom stats predict na.omit
#' @importFrom randomForestSRC rfsrc
#' @importFrom igraph graph_from_edgelist layout_with_fr
#'
#' @export
gng_project <- function(gng_out, make_projection = TRUE) {
  # Project the nodes to a 2D plane
  igr <- igraph::graph_from_edgelist(as.matrix(gng_out$edges), directed = FALSE)
  node_proj <- igraph::layout_with_fr(igr)
  colnames(node_proj) <- c("GNG_x", "GNG_Y")

  # Apply minmax-scale
  mins <- apply(node_proj, 2, min, na.rm = T)
  maxs <- apply(node_proj, 2, max, na.rm = T)
  scale <- maxs - mins
  scale[scale == 0] <- 1
  node_proj <- t(apply(node_proj, 1, function(x) (x - mins) / scale))

  # Map the positions to the original space x
  if (make_projection) {
    rf <- randomForestSRC::rfsrc(Multivar(GNG_x, GNG_Y) ~ ., data.frame(stats::na.omit(S), node_proj, check.names = F))
    pred <- stats::predict(rf, data.frame(x, check.names = F, stringsAsFactors = F))
    space_proj <- sapply(colnames(node_proj), function(n) pred$regrOutput[[n]]$predicted)
  } else {
    space_proj <- NULL
  }

  list(
    node_proj = node_proj,
    space_proj = space_proj,
    igraph = igr
  )
}

#' FlowSOM wrapper
#'
#' @param data fsom$FlowSOM$data as generated by \code{ReadInput}
#' @param colsToUse ids of the columns used
#' @param nNodes the number of nodes
#' @param nPasses the number of passes of the data
#' @param ... Extra params to pass to \code{\link{gng}}
#'
#' @export
flowSOM_wrapper <- function(data, colsToUse, nNodes = 100, nPasses = 50, ...) {
  requireNamespace("FlowSOM")

  max_iter <- nPasses * nrow(data)
  gng.fit <- gng(scale(data[,colsToUse]), max_nodes = nNodes, max_iter = max_iter, ...)
  fsomObject <- list("MST"=list(),"map"=list())
  fsomObject$MST$graph <- gng.fit$igr
  fsomObject$MST$l <- gng.fit$node_proj
  fsomObject$map$mapping <- as.matrix(gng.fit$clustering,ncol=1)
  fsomObject$map$grid <- cbind(seq(nNodes),rep(1,nNodes))
  fsomObject <- FlowSOM::UpdateNodeSize(fsomObject)

  t <- sqrt(table(fsomObject$map$mapping[, 1]))
  rescaled <- 15 * t / max(t)
  fsomObject$MST$size <- numeric(nrow(fsomObject$map$grid))
  fsomObject$MST$size[as.numeric(names(t))] <- rescaled

  return(fsomObject)
}
