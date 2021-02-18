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
dist_forest_irf <- function(fit, x, ntree=NULL, n.core=1) {

  # Pass observations through fitted RF to generate node membership matrix
  read.forest <- readForest(fit$rf.list, x, n.core=n.core, oob.importance=FALSE)

  # Generate all pairwise combinations of observations to compute distance over
  grid <- expand.grid(1:nrow(x), 1:nrow(x))
  grid <- filter(grid, Var1 > Var2)

  # Compute distance between all pairs within each tree and average
  if (is.null(ntree)) ntree <- max(read.forest$tree.info$tree)
  dist.forest <- lapply(1:ntree, dist_tree_irf, read.forest=read.forest, grid=grid)
  rf.dist <- Reduce('+', dist.forest) / ntree

  # Reformat as distance matrix
  n <- nrow(read.forest$node.obs)
  dmat <- matrix(0, nrow=n, ncol=n)
  dmat[lower.tri(dmat, diag=FALSE)] <- rf.dist
  dmat <- dmat + t(dmat)
  
  return(dmat)
}

d1 = dist_tree_irf(read.forest, 1, grid)
d2 = dist_tree_irf_2(read.forest, 1, grid)

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
  nf.k <- read.forest$node.feature[id.k,] != 0
    
  # Map observation pair indices to leaf node pair indices
  no.k <- t(read.forest$node.obs[,id.k])@i + 1
  no.ki <- no.k[grid[,1]]
  no.kj <- no.k[grid[,2]]
  
  # Compoute intersect over union for pairwise ints
  nf.intersect <- Matrix::rowMeans(nf.k[no.ki,] & nf.k[no.kj,])
  nf.union <- Matrix::rowMeans(nf.k[no.ki,] | nf.k[no.kj,])

  return(1 - nf.intersect / nf.union)
}