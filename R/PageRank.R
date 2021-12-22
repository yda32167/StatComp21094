#' This is some descriptio of this function.
#' @title Use R to give a function of PageRank.
#' 
#' @description PageRank is a Google proprietary algorithm used to measure the importance of a particular webpage relative to other webpages in the search engine index.
#' 
#' @details Considering and not considering the situation of the damping system.
#' 
#' @param pages is a page is the adjacency matrix of the webpage.
#'
#' @import knitr
#'
#' @return a vector
#' @export
#' @examples
pagerk <-function(pages){
  
  adjacencyMatrix<-function(pages){
    n<-max(apply(pages,2,max))
    A <- matrix(0,n,n)
    for(i in 1:nrow(pages)) A[pages[i,]$dist,pages[i,]$src]<-1
    A
  }
  
  probabilityMatrix<-function(G){
    cs <- colSums(G)
    cs[cs==0] <- 1
    n <- nrow(G)
    A <- matrix(0,nrow(G),ncol(G))
    for (i in 1:n) A[i,] <- A[i,] + G[i,]/cs
    A
  }

  dProbabilityMatrix<-function(G,d=0.85){
    cs <- colSums(G)
    cs[cs==0] <- 1
    n <- nrow(G)
    delta <- (1-d)/n
    A <- matrix(delta,nrow(G),ncol(G))
    for (i in 1:n) A[i,] <- A[i,] + d*G[i,]/cs
    A
  }
  
  eigenMatrix<-function(G,iter=100){
    iter<-10
    n<-nrow(G)
    x <- rep(1,n)
    for (i in 1:iter) x <- G %*% x
    x/sum(x)
  }
  
  names(pages)<-c("src","dist")
  A<-adjacencyMatrix(pages)
  G<-probabilityMatrix(A)
  q<-eigenMatrix(G,100)
  q
  
  A<-adjacencyMatrix(pages)
  G<-dprobabilityMatrix(A)
  q<-eigenMatrix(G,100)
  q
}