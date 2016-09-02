#' @title Growing Neural Gas
#'
#' @description It is recommended to read the paper by Fritzke Bernd in order to gain a better insight into the given parameters.
#'
#' Step 1, generating the input signal could be improved upon.
#'
#' @references Fritzke, Bernd. "A growing neural gas network learns topologies." Advances in neural information processing systems 7 (1995): 625-632.
#'
#' @param X The input data
#' @param max.iter The max number of iterations
#' @param epsilon.b Move the winning node by epsilon.b times the distance
#' @param epsilon.n Move the neighbours of the winning node by epsilon.n times the distance
#' @param age.max Remove edges older than age.max
#' @param max.nodes The maximum number of nodes
#' @param lambda Insert new nodes every lambda iterations
#' @param alpha The decay parameter for error when a node is added
#' @param beta The decay parameter for error in every node every iteration
#' @param make.projection Project the points from X to a 2D plane
#' @param make.logs The decay for error of all nodes every iteration
#' @param verbose Will output progress if \code{TRUE}
#'
#' @return The GNG (need to be more specific in the future).
#'
#' @importFrom dplyr bind_rows
#' @importFrom igraph graph_from_edgelist layout_with_fr
#' @importFrom randomForestSRC rfsrc
#'
#' @export
#'
#' @examples
#' # To be added soon!
gng <- function(X, max.iter = 20000, epsilon.b = .05, epsilon.n = .001, age.max = 200, max.nodes = 30, lambda = 200, alpha = .5, beta = .99, make.projection = T, make.logs = F, verbose = T) {
  # Calculate ranges of the dimensionality of the dataset
  ranges <- apply(X, 2, range)

  # Initialise a node position matrix
  S <- matrix(NA, nrow = max.nodes, ncol = ncol(X), dimnames = list(NULL, colnames(X)))

  # A dataframe containing the error of each node
  S.meta <- data.frame(i = seq_len(nrow(S)), error = rep(0, nrow(S)))

  # A list containing the neighbours of each node
  E <- lapply(seq_len(max.nodes), function(i) numeric(0))
  E[[1]] <- 2
  E[[2]] <- 1

  # A matrix containing the ages of each edge
  Age <- matrix(NA, nrow = max.nodes, ncol = max.nodes)

  # 0. start with two units a and b at random positions in R^n
  S[1:2,] <- apply(ranges, 2, function(x) {
    runif(2, x[[1]], x[[2]])
  })
  Age[1,2] <- 0
  Age[2,1] <- 0

  # This variable is the index of a new node
  nextr <- 3

  # I separated out the distance and move functions in case I want to try a different approach
  distance.function <- function(Xi, Si) {
    sqrt(mean((Xi - Si)^2))
  }
  move.function <- function(Xi, Si, epsilon) {
    Si + epsilon * (Xi - Si)
  }
  sample.input.signal <- function(X) {
    X[sample.int(nrow(X), 1),]
  }

  if (make.logs) {
    node.loglist <- list(
      data.frame(added = T, node = c(1, 2), iteration = 0)
    )
    nodepos.loglist <- list(
      data.frame(iteration = 0, node = c(1, 2), S[1:2,])
    )
    edge.loglist <- list(
      data.frame(added = T, i = 1, j = 2, iteration = 0, reason = 1, age = 0)
    )
  }

  current.iter <- 0

  # while stopping criterion not met.
  # In the future, this function might check whether nodes have somewhat converged to a steady position.
  while (current.iter < max.iter) {
    current.iter <- current.iter + 1
    if (verbose && current.iter %% 1000 == 0) cat("Iteration ", current.iter, "\n", sep = "")

    # 1. generate input signal
    Xi <- sample.input.signal(X)

    # 2. find nearest unit s1 and second nearest unit s2
    sdist <- apply(S, 1, function(Si) distance.function(Xi, Si))
    sord <- order(sdist)
    s1 <- sord[[1]]
    s2 <- sord[[2]]

    # 3. increment the age of all edges eminating from s1
    Age[s1,] <- Age[s1,] + 1
    Age[,s1] <- Age[,s1] + 1

    # 4. add the squared distance between input signal and the nearest unit to the error variable
    S.meta$error[[s1]] <- S.meta$error[[s1]] + sdist[[s1]]

    # 5. move s1 and its direct topological neighbors towards input signal by franctions epsilon.b and epsilon.n respoctively
    neighs <- E[[s1]]
    S[s1, ] <- move.function(Xi, S[s1,], epsilon.b)
    for (n1 in neighs) {
      S[n1,] <- move.function(Xi, S[n1,], epsilon.n)
    }

    if (make.logs) {
      nodepos.loglist[[length(nodepos.loglist)+1]] <- data.frame(iteration = current.iter, node = c(s1, neighs), S[c(s1, neighs), , drop = F])
    }

    # 6. if s1 and s2 are connected by an edge, set the age of this edge to zero. otherwise, create edge of age 0
    if (is.na(Age[[s1, s2]])) {
      E[[s1]] <- c(E[[s1]], s2)
      E[[s2]] <- c(E[[s2]], s1)
    }
    Age[[s1, s2]] <- 0
    Age[[s2, s1]] <- 0

    edge.nods <- unique(neighs, s2)

    if (make.logs) {
      edge.loglist[[length(edge.loglist)+1]] <- data.frame(added = Age[s1, edge.nods] <= age.max, i = s1, j = edge.nods, iteration = current.iter, reason = 2, age = Age[s1, edge.nods])
    }

    # 7. remove edges with an age larger than age.max. if a point has no remaining edge, remove as well
    rem <- which(Age[s1,] > age.max)
    if (length(rem) > 0) {
      Age[s1, rem] <- NA
      Age[rem, s1] <- NA
      E[[s1]] <- setdiff(E[[s1]], rem)
      for (sj in rem) {
        E[[sj]] <- setdiff(E[[sj]], s1)
      }
      # edge.loglist[[length(edge.loglist)+1]] <- data.frame(added = F, i = s1, j = rem, iteration = current.iter, reason = "age")
      s1rem <- c(s1, rem)
      removed.nodes <- s1rem[which(sapply(s1rem, function(i) length(E[s1rem]) == 0))]

      if (make.logs && length(removed.nodes) > 0) {
        node.loglist[[length(node.loglist)+1]] <- data.frame(added = F, node = removed.nodes, iteration = current.iter)
      }
    }

    # 8. If the number of input signals generated so far is an integer multiple of a parameter lambda, insert a new node
    if (current.iter %% lambda == 0) {
      if (nextr <= max.nodes) {
        r <- nextr
        nextr <- nextr + 1

        # 8a. determine node p with the maximum accumulated error
        p <- which.max(S.meta$error)

        # 8b. determine node q neighbouring p with the largest error
        np <- E[[p]]
        q <- np[which.max(S.meta$error[np])]

        # 8c. insert new node r halfway between p and q
        S[r,] <- (S[p,] + S[q,]) / 2

        if (make.logs) {
          nodepos.loglist[[length(nodepos.loglist)+1]] <- data.frame(iteration = current.iter, node = r, S[r, , drop = F])
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
        Ep <- S.meta$error[[p]]
        Eq <- S.meta$error[[q]]
        S.meta$error[c(p, q, r)] <- c(alpha * Ep, alpha * Eq, alpha * Ep)

        if (make.logs) {
          node.loglist[[length(node.loglist)+1]] <- data.frame(added = T, node = r, iteration = current.iter+.01)
          edge.loglist[[length(edge.loglist)+1]] <- data.frame(added = c(F, T, T), i = c(p, p, q), j = c(q, r, r), iteration = current.iter + .01, reason = 3, age = c(age.max + 1, 0, 0))
        }
      }
    }

    # 9. decrease all variables by multiplying them with a constant beta
    S.meta$error <- S.meta$error * beta
  }

  # Construct GNG output data structures
  nodes <- data.frame(node = seq_len(nrow(S)), S.meta[,2,drop=F])
  node.space <- S
  filt <- !is.na(node.space[,1])
  nodes <- nodes[filt,,drop=F]
  node.space <- node.space[filt,,drop=F]

  edges <- dplyr::bind_rows(lapply(seq(2, nrow(S)), function(nj) {
    ni <- E[[nj]]
    ni <- ni[ni < nj]
    if (length(ni) > 0) {
      data.frame(i = ni, j = nj)
    } else {
      NULL
    }
  }))

  # Determine assignment of the samples to the nodes
  clustering <- apply(X, 1, function(x) which.min(apply(S, 1, function(s) distance.function(x, s))))
  nodes$number <- sapply(seq_len(nrow(nodes)), function(i) sum(clustering == i))

  # Project the nodes to a 2D plane
  igr <- igraph::graph_from_edgelist(as.matrix(edges), directed = F)
  node.proj <- igraph::layout_with_fr(igr)
  colnames(node.proj) <- c("GNG_X", "GNG_Y")

  # Apply minmax-scale
  mins <- apply(node.proj, 2, min, na.rm = T)
  maxs <- apply(node.proj, 2, max, na.rm = T)
  scale <- maxs - mins
  scale[scale == 0] <- 1
  node.proj <- t(apply(node.proj, 1, function(x) (x - mins) / scale))

  # Map the positions to the original space X
  if (make.projection) {
    nodes.df <- data.frame(S, node.proj, check.names = F)
    rf <- randomForestSRC::rfsrc(Multivar(GNG_X, GNG_Y) ~ ., nodes.df)
    pred <- predict(rf, data.frame(X, check.names = F, stringsAsFactors = F))
    space.proj <- sapply(colnames(node.proj), function(n) pred$regrOutput[[n]]$predicted)
  } else {
    space.proj <- NULL
  }

  # Process the logs
  if (make.logs) {
    edge.log <- dplyr::bind_rows(edge.loglist)
    node.log <- dplyr::bind_rows(node.loglist)
    nodepos.log <- dplyr::bind_rows(nodepos.loglist)
    edge.log[,2:3] <- t(apply(edge.log[,2:3], 1, sort))
  } else {
    edge.log <- NULL
    node.log <- NULL
    nodepos.log <- NULL
  }

  list(
    nodes = nodes,
    node.space = node.space,
    edges = edges,
    node.log = node.log,
    nodepos.log = nodepos.log,
    edge.log = edge.log,
    clustering = clustering,
    node.proj = node.proj,
    space.proj = space.proj,
    igraph = igr
  )
}

#' FlowSOM wrapper
#'
#' @param data fsom$FlowSOM$data as generated by \code{ReadInput}
#' @param colsToUse ids of the columns used
#' @param nNodes the number of nodes
#' @param nPasses the number of passes of the data
#'
#' @export
flowSOM_wrapper <- function(data, colsToUse, nNodes = 100, nPasses = 50, ...) {
  max.iter <- nPasses * nrow(data)
  gng.fit <- gng(scale(data[,colsToUse]), max.nodes = nNodes, max.iter = max.iter, ...)
  fsomObject <- list("MST"=list(),"map"=list())
  fsomObject$MST$graph <- gng.fit$igr
  fsomObject$MST$l <- gng.fit$node.proj
  fsomObject$map$mapping <- as.matrix(gng.fit$clustering,ncol=1)
  fsomObject$map$grid <- cbind(seq(nNodes),rep(1,nNodes))
  fsomObject <- FlowSOM::UpdateNodeSize(fsomObject)

  t <- sqrt(table(fsomObject$map$mapping[, 1]))
  rescaled <- 15 * t / max(t)
  fsomObject$MST$size <- numeric(nrow(fsomObject$map$grid))
  fsomObject$MST$size[as.numeric(names(t))] <- rescaled

  return(fsomObject)
}
