#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

arma::mat outer(arma::mat);
arma::mat treeDist(arma::mat, int);

//' Computes distances between pairs of observations using iRF-encoded binary
//' features.
//'
//' @param obsn node membership matrix indicating which leaf node observations
//'   fall in for each tree, rows being observations, cols being trees.
//' @param featn binary matrix indicating iRF decision path encodings, rows
//'   being nodes, cols being features.
//' @param idtree vector indicating rows of featn that correspond to new trees.
//' @param ridx first index of a pair for which to compute distance.
//' @param cidx second index of a pair for which to compute distance.
//' @return Product of v1 and v2
// [[Rcpp::export]]
arma::mat forestDist(arma::umat obsn,
                     arma::mat featn,
                     arma::uvec idtree,
                     arma::uvec ridx,
                     arma::uvec cidx) {
  
  // pre-allocate variables for tree distances
  arma::mat nfk;
  arma::mat distk;
  int pnode;
  arma::uvec nok;
  
  // pre-allocate matrix for pairwise distance 
  arma::mat rfdist(ridx.n_elem, 1, arma::fill::zeros);

  for (int k = 0; k < (idtree.n_elem - 1); k++) {

    // compute leaf node distances for tree k
    nfk = featn.rows(idtree[k] + 1, idtree[k + 1]);
    distk = treeDist(nfk, featn.n_cols);
  
    // map leaf node distances to observation distance in tree k
    pnode = distk.n_cols;
    nok = obsn.col(k);
    arma::uvec indices = nok(ridx) + pnode * nok(cidx);
    rfdist += distk(nok(ridx) + pnode * nok(cidx));
  }
  
  return rfdist;
}

// [[Rcpp::export]]
arma::mat treeDist(arma::mat x, int nk) {
  arma::mat intk = outer(x);
  arma::mat unik = nk - outer(1 - x);
  return intk / unik;
}

arma::mat outer(arma::mat x) {
  return x * x.t();
}
