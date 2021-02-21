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
#' @importFrom BBmisc chunk
#' 
# TODO: compare rf to irf with one selected node from node.feature
# TODO: subset to nodes based on RF weight
dist_forest_irf <- function(fit, x, read.forest=NULL, n.core=1) {

  
  # Pass observations through fitted RF to generate node membership matrix
  if (is.null(read.forest)) {
    read.forest <- readForest(fit$rf.list, x, n.core=n.core, oob.importance=FALSE)
  }

  # Generate all pairwise combinations of observations to compute distance over
  grid <- expand.grid(1:nrow(x), 1:nrow(x))
  grid <- filter(grid, Var1 > Var2)

  # Appears to be faster than sparse matrix operations....
  read.forest$node.feature <- as.matrix(read.forest$node.feature) != 0
  
  # Compute distance between all pairs within each tree and average
  ntree <- max(read.forest$tree.info$tree)
  #tree.blocks <- chunk(1:ntree, n.chunks=n.core)
  
  dist.forest <- lapply(1:ntree, dist_tree_irf, 
                        read.forest=read.forest, 
                        grid=grid)
  
  dist.forest <- 1 - Reduce('+', dist.forest) / ntree
  # TODO: check that this is returning correct distances - unit test

  # Reformat as distance matrix
  n <- nrow(read.forest$node.obs)
  dmat <- matrix(0, nrow=n, ncol=n)
  dmat[lower.tri(dmat, diag=FALSE)] <- dist.forest
  dmat <- dmat + t(dmat)
  
  return(dmat)
  
}

dist_trees_irf <- function(read.forest, k, grid=NULL) {
  # Wrapper function to run dist_tree_irf over multiple trees
  out <- lapply(k, dist_tree_irf, read.forest=read.forest, grid=grid)
  return(out)
}

dist_tree_irf <- function(read.forest, k, grid=NULL) {
  # Compute pairwise sirf distance for tree-k.
  #
  # Args:
  #   read.forest: output of iRF::readForest.
  #   k: tree to calculate distances for
  #   grid
  
  if (is.null(grid)) {
    # Generate all pairwise combinations of observations to compute distance over
    grid <- expand.grid(1:nrow(read.forest$node.obs), 1:nrow(read.forest$node.obs))
    grid <- filter(grid, Var1 > Var2)
  }
  
  # Filter to nodes from given tree
  id.k <- read.forest$tree.info$tree == k

  # Map observation pair indices to leaf node pair indices
  no.k <- t(read.forest$node.obs[,id.k])@i + 1

  # Compute intersect over union across leaf nodes  
  nf.t <- t(read.forest$node.feature[id.k,])
  k.intersect <- (read.forest$node.feature[id.k,]) %*% nf.t 
  k.union <- sum(id.k) - (!read.forest$node.feature[id.k,]) %*% !nf.t
  dist.k <- k.intersect / k.union
  pnode <- ncol(dist.k)
  
  return(dist.k[no.k[grid[,1]] + pnode * (no.k[grid[,2]] - 1)])
}

