# gng_r exists only for educational purposes
#' @importFrom stats runif
gng_r <- function(x, max_iter, epsilon_b, epsilon_n, age_max, max_nodes, lambda, alpha, beta, verbose, make_logs_at) {
  # Calculate ranges of the dimensionality of the dataset
  ranges <- apply(x, 2, range)

  # Initialise a node position matrix
  S <- matrix(NA, nrow = max_nodes, ncol = ncol(x), dimnames = list(NULL, colnames(x)))

  # A dataframe containing the error of each node
  S_meta <- data.frame(i = seq_len(nrow(S)), error = rep(0, nrow(S)))

  # A list containing the neighbours of each node
  E <- map(seq_len(max_nodes), function(i) numeric(0))
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

  current.iter <- 0L

  # create data structures for saving the intermediary gngs
  log <- list()

  if (current.iter %in% make_logs_at) {
    current_log <- list(current.iter = current.iter, S = S, S_meta = S_meta, Age = Age, E = E)
    log[[length(log) + 1]] <- current_log
  }

  # while stopping criterion not met.
  # In the future, this function might check whether nodes have somewhat converged to a steady position.
  while (current.iter <= max_iter) {
    current.iter <- current.iter + 1L
    if (verbose && current.iter %% 1000L == 0L) cat("Iteration ", current.iter, "\n", sep = "")

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
    S[s1, ] <- move_function(xi, S[s1,], epsilon_b)

    neighs <- E[[s1]]
    for (n1 in neighs) {
      S[n1,] <- move_function(xi, S[n1,], epsilon_n)
    }

    # 6. if s1 and s2 are connected by an edge, set the age of this edge to zero. otherwise, create edge of age 0
    if (is.na(Age[[s1, s2]])) {
      E[[s1]] <- c(E[[s1]], s2)
      E[[s2]] <- c(E[[s2]], s1)
    }
    Age[[s1, s2]] <- 0
    Age[[s2, s1]] <- 0

    edge.nods <- unique(neighs, s2)

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
      }
    }

    # 9. decrease all variables by multiplying them with a constant beta
    S_meta$error <- S_meta$error * beta

    # store intermediary gng if so desired
    if (current.iter %in% make_logs_at) {
      current_log <- list(current.iter = current.iter, S = S, S_meta = S_meta, Age = Age, E = E)
      log[[length(log) + 1]] <- current_log
    }
  }

  # Construct GNG output data structures
  nodes <- data.frame(node = seq_len(nrow(S)), S_meta[, 2, drop = FALSE])
  node_space <- S
  filt <- !is.na(node_space[,1])
  nodes <- nodes[filt, , drop = FALSE]
  node_space <- node_space[filt, , drop = FALSE]

  edges <- map_df(seq(2, nrow(S)), function(nj) {
    ni <- E[[nj]]
    ni <- ni[ni < nj]
    if (length(ni) > 0) {
      tibble(i = ni, j = nj)
    } else {
      NULL
    }
  })

  list(
    nodes = nodes,
    node_space = node_space,
    edges = edges,
    log = log
  )
}
