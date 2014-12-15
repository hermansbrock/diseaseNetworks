#' Network simulation code
#' 
#' This code simulates a network using either the Erdos, Waxman or Keeling network method.
#' 
#' @param n Number of nodes.
#' @param method String, determines how we simulate the network. Either 'erdos', 'waxman' or 'keeling'. No capitals.
#' @param percent Determines the percentage of edges we retain. Used in method='erods'.
#' @param alpha Ratio of long edges to short edges. Used in method='waxman' or 'keeling'.
#' @param delta Probability of an edge between any two nodes. Used in method='waxman' or 'keeling'.
#' @param Fn Number of focal points. Used in method='keeling'.
#' @param f Proportion points are moved towards closest focal point. Used in method='keeling'.
#' @return nodes nx2 matrix, rows represent the (x,y) coordinate of a node.
#' @return edges Matrix, row i represent the existence of an edge between nodes in first two columns. Third column is set to 0 as default for simulating diseases.
#' @return A Adjacency matrix.
#' @return node_states Vector with each node state: 0=susceptible; 1=infective; 2=recovered.
#' @author Brock Hermans <brock.hermans@@adelaide.edu.au>
#' @export
#' @note 9 December 2014
#' @examples
#' net <- network(10,'keeling',1,0.5,0.5,4,0.9)
network <- function(n,method,prob=1,alpha=0.5,delta=0.5,Fn=1,f=1){
  if (method=='erdos'){
    net <- erdos(n,prob)
  }
  else if (method=='waxman'){
    net <- waxman(n,alpha,delta)
  }
  else if (method=='keeling'){
    net <- keeling(n,Fn,f,alpha,delta)
  }
  else {
    stop("Invalid method: Must be 'erdos', 'waxman' or 'keeling'")
  }
  
  edges2 <- matrix(c(0,0,0),ncol=3)
  for (i in 1:nrow(net$edges)){
    if (net$edges[i,4]==1){
      edges2 <- rbind(edges2,c(net$edges[i,5]-1,net$edges[i,6]-1,0))
      edges2 <- rbind(edges2,c(net$edges[i,6]-1,net$edges[i,5]-1,0))
    }
  }
  
  edges2 <- edges2[-1,]
  net$edges <- edges2
  class(net) = "network"
  return(net)
}
print.network = function(x,...){
  cat("First 6 nodes\n\n")
  print(head(x$nodes))
  cat("First 6 edges\n\n")
  print(head(x$edges))
  print(head(x$graph))
}