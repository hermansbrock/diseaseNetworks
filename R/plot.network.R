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
#' network <- network(10,0.5)
#' plot(network)
plot.network <- function(x,size=4,...){
  require(ggplot2)  
  df = x$nodes
  df <- data.frame(cbind(df,x$node_states))
  colnames(df) <- c("x","y","s1")
  df$s1[df$s1==0] = "Susceptible"
  df$s1[df$s1==2] = "Recovered"
  df$s1[df$s1==1] = "Infected"
  x$edges <- x$edges+1
  
  from <- df[x$edges[,1],]
  to <- df[x$edges[,2],]
  segments = data.frame(x1=from$x,y1=from$y,x2 = to$x,y2=to$y)
  p = ggplot(aes(x,y,col=s1),data=df)
  p <- p +
    geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),data=segments,col="black") +
    labs(col="Status")
  p <- p + geom_point(size=size) + theme_bw()+xlim(0,1)+ylim(0,1)
  return(p)  
}


