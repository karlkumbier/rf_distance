#' Signed iterative random forest distance (dist_sirf)
#'
#' Compute sirf-based distance between all pairs of observations passed through
#' a fitted random forest.
#'
#' @param fit a fitted random forest of class ranger (probability = FALSE) or
#' randomForest.
#' @param x feature matrix, indicating observations to compute pairwise
#' distances for.
#' @param n.core number of workers to distribute readForest over.
#'
#'@return a numeric matrix indicating sirf-based distance over each pair of
#'observations in x.
#'
#' @export
#'
#' @importFrom iRF readForest
#' @importFrom Matrix t Matrix rowMeans
#' @importFrom parallel mclapply
dist_forest_irf <- function(fit, x, rf=NULL, n.core=1) {

  # Pass observations through fitted RF to generate node membership matrix
  if (is.null(rf)) {
    rf <- readForest(fit, x, n.core=n.core, oob.importance=FALSE)
  }

  # Generate all pairwise combinations of observations to compute distance over
  n <- nrow(x)
  g1 = rep(1:n, times=((1:n) - 1))
  g2 = unlist(lapply(1:(n-1), function(i) 1:i), use.names=FALSE)
  
  # Appears to be faster than sparse matrix operations....
  rf$node.feature <- as.matrix(rf$node.feature) != 0
  
  # Compute distance between all pairs within each tree and average
  ntree <- max(rf$tree.info$tree)
  dist.forest <- numeric(length(g1))
  for (k in 1:ntree) {
    dist.forest <- dist.forest + dist_tree_irf(rf, k, g1, g2)
  }
  

  dist.forest <- dist.forest / ntree
  dist.forest <- 1 - dist.forest

  # Reformat as distance matrix
  dmat <- matrix(0, nrow=n, ncol=n)
  # TODO: replace lower.tri with indices
  dmat[lower.tri(dmat, diag=FALSE)] <- dist.forest
  dmat <- dmat + t(dmat)

  
  return(dmat)

}

#' Signed iterative random forest distance (dist_forest_irf_cpp)
#'
#' Compute sirf-based distance between all pairs of observations passed through
#' a fitted random forest, making calls to CPP distance implementation.
#'
#' @param fit a fitted random forest of class ranger (probability = FALSE) or
#' randomForest.
#' @param x feature matrix, indicating observations to compute pairwise
#' distances for.
#' @param n.core number of workers to distribute readForest over.
#'
#'@return a numeric matrix indicating sirf-based distance over each pair of
#'observations in x.
#'
#' @export
#'
#' @importFrom iRF readForest
#' @importFrom Matrix t Matrix rowMeans
#' @importFrom parallel mclapply
#' 
# TODO: replace lower.tri for optimization
# TODO: subset to nodes based on RF weight
dist_forest_irf_cpp <- function(fit, x, rf=NULL, n.core=1) {
  
  n <- nrow(x)

  # Pass observations through fitted RF to generate node membership matrix
  if (is.null(rf)) {
    rf <- readForest(fit, x, n.core=n.core, oob.importance=FALSE)
  }  
  
  # Read out node feature input
  nf <- as.matrix(rf$node.feature != 0)
  
  # Read out tree indices for each leaf node
  ntree <- fit$num.trees
  trees <- rf$tree.info$tree

  # Get index corresponding to first leaf node for each tree
  tree.id <- sapply(1:ntree, function(k) min(which(trees == k)))
  
  # Convert node.obs into n x ntree node memership matrix
  no.t <- t(rf$node.obs)
  no <- t(matrix(no.t@i, ncol=n))
  n.node.tree <- cumsum(c(0, diff(tree.id)))
  adj <- matrix(rep(n.node.tree, each=n), nrow=n) 
  no <- no - adj  
  
  # Adjust tree indexing for cpp call
  tree.id <- c(tree.id - 2, length(trees) - 1)
  
  # Initialize grid of all pairwise combinations for n observaitons
  grid <- expand.grid(0:(n-1), 0:(n-1)) %>% filter(Var1 > Var2)
  
  # Compute distance
  dist.forest <- forestDist(no, nf, tree.id, grid$Var1, grid$Var2)
  dist.forest <- dist.forest / ntree
  dist.forest <- 1 - dist.forest
  
  # Reformat as distance matrix
  dmat <- matrix(0, nrow=n, ncol=n)
  dmat[lower.tri(dmat, diag=FALSE)] <- dist.forest
  dmat <- dmat + t(dmat)
  return(dmat)
}

dist_tree_irf <- function(rf, k, g1, g2) {
  # Compute pairwise sirf distance for tree-k.
  #
  # Args:
  #   rf: output of iRF::readForest.
  #   k: tree to calculate distances for
  
  # Filter to nodes from given tree
  id.k <- rf$tree.info$tree == k

  # Map observation pair indices to leaf node pair indices
  no.k <- t(rf$node.obs[,id.k])@i + 1

  # Compute intersect over union across leaf nodes  
  nf.t <- t(rf$node.feature[id.k,])
  k.intersect <- (rf$node.feature[id.k,]) %*% nf.t 
  k.union <- sum(id.k) - (!rf$node.feature[id.k,]) %*% !nf.t
  dist.k <- k.intersect / k.union
  pnode <- ncol(dist.k)
  return(dist.k[no.k[g1] + pnode * (no.k[g2] - 1)])
}


