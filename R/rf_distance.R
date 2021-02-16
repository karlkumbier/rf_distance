#' Random forest similarity
#'
#' Compute rf similarirty between all pairs of observations passed through
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
dist_forest_rf <- function(fit, x) {

  # Pass data points through rf to generate leaf node membership
  leaf.nodes <- predict(fit$rf.list, x, type='terminalNodes')$predictions

  # Generate all pairwise combinations of observations to compute distance over
  grid <- expand.grid(1:nrow(x), 1:nrow(x))
  grid <- filter(grid, Var1 > Var2)

  # Compute pairwise distances over each node in a tree
  rf.dist <- apply(leaf.nodes, MAR=2, dist_tree_rf)
  rf.dist <- rowMeans(rf.dist)

  # Average across trees in the rf  
  dmat <- matrix(0, nrow=nrow(x), ncol=nrow(x))
  dmat[lower.tri(dmat, diag=FALSE)] <- rf.dist
  dmat <- dmat + t(dmat)

  return(dmat)
}

dist_tree_rf <- function(leaf.nodes, grid=NULL) {
  
  if (is.null(grid)) {
    # Generate all pairwise combinations of observations to compute distance over
    grid <- expand.grid(1:length(leaf.nodes), 1:length(leaf.nodes))
    grid <- filter(grid, Var1 > Var2)
  }

  dmat <- as.numeric(leaf.nodes[grid[,1]] != leaf.nodes[grid[,2]])
  return(dmat)
}
