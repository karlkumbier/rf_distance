// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// forestDist
arma::mat forestDist(arma::umat obsn, arma::mat featn, arma::uvec idtree, arma::uvec ridx, arma::uvec cidx);
RcppExport SEXP _rfdistance_forestDist(SEXP obsnSEXP, SEXP featnSEXP, SEXP idtreeSEXP, SEXP ridxSEXP, SEXP cidxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::umat >::type obsn(obsnSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type featn(featnSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type idtree(idtreeSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ridx(ridxSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type cidx(cidxSEXP);
    rcpp_result_gen = Rcpp::wrap(forestDist(obsn, featn, idtree, ridx, cidx));
    return rcpp_result_gen;
END_RCPP
}
// treeDist
arma::mat treeDist(arma::mat x, int nk);
RcppExport SEXP _rfdistance_treeDist(SEXP xSEXP, SEXP nkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type nk(nkSEXP);
    rcpp_result_gen = Rcpp::wrap(treeDist(x, nk));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rfdistance_forestDist", (DL_FUNC) &_rfdistance_forestDist, 5},
    {"_rfdistance_treeDist", (DL_FUNC) &_rfdistance_treeDist, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_rfdistance(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
