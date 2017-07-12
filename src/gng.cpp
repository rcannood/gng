#include <Rcpp.h>
#include <boost/config.hpp>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/lexical_cast.hpp>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

using namespace Rcpp;

struct node_prop_t {
  int index;
  double error;
  double dist;
  NumericVector position;
};
struct edge_prop_t {
  int age;
};

typedef boost::adjacency_list < boost::listS, boost::vecS, boost::undirectedS, node_prop_t, edge_prop_t > Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::graph_traits<Graph>::vertex_iterator Vertex_iter;
typedef boost::graph_traits<Graph>::edge_iterator Edge_iter;
typedef boost::graph_traits<Graph>::out_edge_iterator Edge_o_iter;
typedef boost::graph_traits<Graph>::adjacency_iterator Adj_iter;

//' @useDynLib GNG
// [[Rcpp::export]]
List gng_cpp(NumericMatrix x, int max_iterations, float epsilon_b, float epsilon_n, int age_max, int max_nodes, int lambda, float alpha, float beta, bool verbose) {
  RNGScope rngScope;

  Progress p(max_iterations, verbose);

  std::pair<Edge_o_iter, Edge_o_iter> eop;
  std::pair<Edge_iter, Edge_iter> ep;
  std::pair<Vertex_iter, Vertex_iter> vp;
  std::pair<Adj_iter, Adj_iter> ap;
  Vertex s0, s1;
  std::pair<Edge, bool> e01;

  // Calculate ranges of the dimensionality of the dataset
  int num_features = x.ncol();
  int num_samples = x.nrow();
  NumericVector mins(num_features), maxs(num_features);

  for (int i = 0; i < num_features; i++) {
    mins(i) = min(x(_,i));
    maxs(i) = max(x(_,i));
  }

  // 0a. Initialise the graph with two nodes
  Graph graph(2);
  int next_r = 0;

  // 0b. The locations are uniformily sampled
  for (vp = vertices(graph); vp.first != vp.second; ++vp.first) {
    Vertex node = *vp.first;
    NumericVector pos(num_features);
    for (int j = 0; j < num_features; j++) {
      pos(j) = R::runif(mins(j), maxs(j));
    }
    graph[node].position = pos;
    graph[node].index = next_r;
    graph[node].error = 0;
    next_r++;
  }

  // 0b. Add an edge between the two nodes
  vp = vertices(graph);
  s0 = *vp.first;
  vp.first++;
  s1 = *vp.first;
  e01 = edge(s0, s1, graph);
  if (!e01.second) {
    e01 = boost::add_edge(s0, s1, graph);
  }
  graph[e01.first].age = 0;

  // While stopping criterion not met.
  // In the future, this function might check whether nodes have converged
  int current_iter = 0;
  while (current_iter < max_iterations) {
    ++current_iter;

    p.increment(); // update progress

    if (current_iter % 10000 == 0 && Progress::check_abort() )
      return NULL;

    // 1. generate input signal
    int i = floor(R::runif(0, num_samples));
    NumericVector x_i = x(i, _);

    // 2a. calculate distances between each node and the input signal
    for (vp = vertices(graph); vp.first != vp.second; ++vp.first) {
      Vertex node = *vp.first;
      NumericVector n_i = graph[node].position;
      graph[node].dist = sqrt(sum(pow(x_i - n_i, 2.0)));
    }

    //  find nearest node s0 and second nearest node s1
    vp = vertices(graph);
    s0 = *vp.first;
    vp.first++;
    s1 = *vp.first;
    vp.first++;
    if (graph[s1].dist < graph[s0].dist) {
      Vertex swap = s1;
      s1 = s0;
      s0 = swap;
    }

    for (; vp.first != vp.second; ++vp.first) {
      Vertex node = *vp.first;
      if (graph[node].dist < graph[s0].dist) {
        s1 = s0;
        s0 = node;
      } else if (graph[node].dist < graph[s1].dist) {
        s1 = node;
      }
    }

    // 3. increment the age of all edges eminating from s0
    for (eop = boost::out_edges(s0, graph); eop.first != eop.second; ++eop.first) {
      graph[*eop.first].age++;
    }

    // 4. add the distance between the input signal and the s0 to the error
    graph[s0].error += graph[s0].dist;

    // 5. move s0 and its direct topological neighbours toward input signal by
    //    fractions epsilon_b and epsilon_n respectively
    graph[s0].position += epsilon_b * (x_i - graph[s0].position);

    for (ap = boost::adjacent_vertices(s0, graph); ap.first != ap.second; ++ap.first) {
      Vertex node = *ap.first;
      graph[node].position += epsilon_n * (x_i - graph[node].position);
    }

    // 6. if s0 and s1 are connected by an edge, set the age of this edge to zero.
    //    otherwise, create a new edge of age 0
    e01 = edge(s0, s1, graph);
    if (!e01.second) {
      e01 = boost::add_edge(s0, s1, graph);
    }
    graph[e01.first].age = 0;

    // 7. remove edges with an age larger than age_max.
    ep = edges(graph);
    for (Edge_iter enext = ep.first; ep.first != ep.second; ep.first = enext) {
      ++enext;
      Edge edge = *ep.first;
      if (graph[edge].age > age_max) {
        Vertex from = boost::source(edge, graph);
        Vertex to = boost::target(edge, graph);

        boost::remove_edge(edge, graph);
        if (boost::degree(from, graph) == 0) {
          boost::remove_vertex(from, graph);
        }
        if (boost::degree(to, graph) == 0) {
          boost::remove_vertex(to, graph);
        }
      }
    }

    // 8. insert a new node
    if (current_iter % lambda == 0 && boost::num_vertices(graph) < max_nodes) {
      // 8a. determine node p with maximum accumulated error
      vp = boost::vertices(graph);
      Vertex p = *vp.first;
      vp.first++;
      for (; vp.first != vp.second; ++vp.first) {
        Vertex node = *vp.first;
        if (graph[p].error < graph[node].error) {
          p = node;
        }
      }

      // 8b. determine node q neighbouring p with the largest error
      ap = adjacent_vertices(p, graph);
      Vertex q = *ap.first;
      ap.first++;
      for (; ap.first != ap.second; ++ap.first) {
        Vertex node = *ap.first;
        if (graph[q].error < graph[node].error) {
          q = node;
        }
      }

      // 8c. insert a new node r halfway between p and q
      Vertex r = boost::add_vertex(graph);
      graph[r].index = next_r;
      next_r++;
      graph[r].position = (graph[p].position + graph[q].position) / 2;
      graph[r].error = 0;

      // 8d. remove (p, q) and add (p, r) and (q, r)
      remove_edge(p, q, graph);
      e01 = add_edge(p, r, graph);
      graph[e01.first].age = 0;
      e01 = add_edge(q, r, graph);
      graph[e01.first].age = 0;

      // 8e. decrease error or p and q by factor alpha
      //     and initialise error of r as the error of p
      graph[p].error *= alpha;
      graph[q].error *= alpha;
      graph[r].error = graph[p].error;
    }

    // 9. decrease all errors by factor beta
    for (vp = boost::vertices(graph); vp.first != vp.second; ++vp.first) {
      Vertex node = *vp.first;
      graph[node].error *= beta;
    }
  }

  // Create nodes output
  int num_nodes = boost::num_vertices(graph);
  CharacterVector nodedf_names(num_nodes);
  IntegerVector nodedf_index(num_nodes);
  NumericVector nodedf_error(num_nodes);
  NumericMatrix node_space (num_nodes, num_features);

  int nodedf_i = 0;
  for (vp = vertices(graph); vp.first != vp.second; ++vp.first) {
    Vertex node = *vp.first;
    int index = graph[node].index+1;
    nodedf_index(nodedf_i) = index;
    nodedf_names(nodedf_i) = "Comp" + boost::lexical_cast<std::string>(index);
    nodedf_error(nodedf_i) = graph[node].error;
    node_space(nodedf_i, _) = graph[node].position;
    nodedf_i++;
  }

  // CharacterVector nodedf_names = CharacterVector("Comp") + nodedf_index;
  DataFrame nodes = DataFrame::create(
    Named("name") = nodedf_names,
    Named("index") = nodedf_index,
    Named("error") = nodedf_error);

  rownames(node_space) = nodedf_names;
  colnames(node_space) = colnames(x);

  // Create edges output
  int num_edges = boost::num_edges(graph);
  IntegerVector edgedf_from(num_edges), edgedf_to(num_edges);
  NumericVector edgedf_age(num_edges);

  int edgedf_i = 0;
  for (ep = edges(graph); ep.first != ep.second; ++ep.first) {
    Edge edge = *ep.first;

    edgedf_from(edgedf_i) = graph[boost::source(edge, graph)].index+1;
    edgedf_to(edgedf_i) = graph[boost::target(edge, graph)].index+1;
    edgedf_age(edgedf_i) = graph[edge].age;

    edgedf_i++;
  }

  DataFrame edges = DataFrame::create(
    Named("i") = edgedf_from,
    Named("j") = edgedf_to,
    Named("age") = edgedf_age);

  // also return clustering

  // Return everything
  return List::create(
    _["nodes"] = nodes,
    _["node_space"] = node_space,
    _["edges"] = edges
  );
}

// x <- matrix(runif(200), ncol = 10)
// Rcpp::sourceCpp('src/gng.cpp')
// gng_cpp(x, max_iterations = 20000, epsilon_b = .05, epsilon_n = .001, age_max = 200, max_nodes = 30, lambda = 200, alpha = .5, beta = .99, verbose = T, weights = .1)
