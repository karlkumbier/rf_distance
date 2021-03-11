library(iRF)
library(Matrix)
library(parallel)

Rcpp::sourceCpp('~/github/Rpackages/rf_distance/src/dist.cpp')
source('~/github/Rpackages/rf_distance/R/irf_distance.R')

################################################################################
# Test that irf distance can computed on small fit and returns relevant object
################################################################################
n <- 50
p <- 10
x <- matrix(rnorm(n * p), nrow=n)
y <- as.numeric(x[,1] > 0 & x[,2] > 0)

fit <- iRF(x=x, y=as.factor(y), n.iter=1, type='ranger')
rf.dist <- dist_forest_irf_cpp(fit$rf.list, x)

all(dim(rf.dist) == c(n, n))
all(diag(rf.dist) == 0)
all(rf.dist <= 1 & rf.dist >= 0)

################################################################################
# Check that irf distance is computing correct distances for simple test case
################################################################################
# Leaf node distance = 1
nf <- cbind(
  rep(1:0, times=2),
  rep(0:1, times=2)
)

# Two trees
trees <- rep(1:2, each=2)

# Observations overlap in 1 or no trees
no <- rbind(
  c(1, 0, 1, 0),
  c(0, 1, 1, 0),
  c(1, 0, 0, 1),
  c(0, 1, 0, 1),
)

# Format inputs for distance function
fit <- list()
fit$num.trees <- 2

rf <- list()
rf$node.feature <- nf
rf$node.obs <- Matrix(no, sparse=TRUE)
rf$tree.info <- list(tree=trees)

x <- matrix(0, nrow=nrow(no), ncol=1)
rf.dist <- dist_forest_irf_cpp(fit, x, rf)

# TODO: test for indexing issue with grid - i.e. "grid" indices do not match
# lower triangle indices
