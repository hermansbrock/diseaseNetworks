#' Plot network objects
#' 
#' Takes a network structure and plots the network. Also plots subgraphs if given
#' 
#' @param x Network object.
#' @param size Integer of size of points in plot.
#' @return plot of network
#' @author Brock Hermans <brock.hermans@@adelaide.edu.au
#' @export
#' @note 29 August 2014
#' @examples
#' network <- networkSimulate(0.1,10)
#' plot(network)
plot.network <- function(x,size=4,...) {
  require(ggplot2)  
  df = x$points
  df$s1[df$s1==0] = "Susceptible"
  df$s1[df$s1==2] = "Recovered"
  df$s1[df$s1==1] = "Infected"
  from <- df[x$edges[,1],]
  to <- df[x$edges[,2],]
  segments = data.frame(x1=from$x,y1=from$y,x2 = to$x,y2=to$y)
  p = ggplot(aes(x,y,col=s1),data=df)
  #
  # If we have subgraphs identified, plot the network with subgraphs identified
  if(!is.null(x$graph)){
    segments$graph = factor(x$graph[x$edges[,1],2])
    p = p + 
      geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2,col=graph),data=segments) +
      labs(col="Subgraphs")
  }
  #
  # If we don't know the subgraphs, do a simple plot
  else {p <- p +
          geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),data=segments,col="black") +
          labs(col="Status")}
  p <- p + geom_point(size=size) + theme_bw()
  return(p)  
}