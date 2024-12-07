// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// projected_gradient_descent
arma::vec projected_gradient_descent(arma::mat matrix, double epsilon, int max_iter, double tol);
RcppExport SEXP _SA24204187_projected_gradient_descent(SEXP matrixSEXP, SEXP epsilonSEXP, SEXP max_iterSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type matrix(matrixSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(projected_gradient_descent(matrix, epsilon, max_iter, tol));
    return rcpp_result_gen;
END_RCPP
}
// gibbsSamplingC
NumericMatrix gibbsSamplingC(int N, int thin);
RcppExport SEXP _SA24204187_gibbsSamplingC(SEXP NSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(gibbsSamplingC(N, thin));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SA24204187_projected_gradient_descent", (DL_FUNC) &_SA24204187_projected_gradient_descent, 4},
    {"_SA24204187_gibbsSamplingC", (DL_FUNC) &_SA24204187_gibbsSamplingC, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_SA24204187(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
