// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// scorebed
arma::mat scorebed(const std::string fileName, int N, int P, const arma::mat input, arma::Col<int> extract, arma::Col<int> keep);
RcppExport SEXP _OAPRS_scorebed(SEXP fileNameSEXP, SEXP NSEXP, SEXP PSEXP, SEXP inputSEXP, SEXP extractSEXP, SEXP keepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type fileName(fileNameSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type input(inputSEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type extract(extractSEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type keep(keepSEXP);
    rcpp_result_gen = Rcpp::wrap(scorebed(fileName, N, P, input, extract, keep));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_OAPRS_scorebed", (DL_FUNC) &_OAPRS_scorebed, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_OAPRS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
