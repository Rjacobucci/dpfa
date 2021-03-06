// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rztpois_single
int rztpois_single(double lambda);
RcppExport SEXP _dpfa_rztpois_single(SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(rztpois_single(lambda));
    return rcpp_result_gen;
END_RCPP
}
// rztpois_cpp
arma::rowvec rztpois_cpp(unsigned int n, double lambda);
RcppExport SEXP _dpfa_rztpois_cpp(SEXP nSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(rztpois_cpp(n, lambda));
    return rcpp_result_gen;
END_RCPP
}
// calcC_kn
arma::Cube<int> calcC_kn(arma::Cube<int> ZZip_3D, arma::vec bias_0, arma::Cube<int> W_3D, arma::mat Phi);
RcppExport SEXP _dpfa_calcC_kn(SEXP ZZip_3DSEXP, SEXP bias_0SEXP, SEXP W_3DSEXP, SEXP PhiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Cube<int> >::type ZZip_3D(ZZip_3DSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type bias_0(bias_0SEXP);
    Rcpp::traits::input_parameter< arma::Cube<int> >::type W_3D(W_3DSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Phi(PhiSEXP);
    rcpp_result_gen = Rcpp::wrap(calcC_kn(ZZip_3D, bias_0, W_3D, Phi));
    return rcpp_result_gen;
END_RCPP
}
// calcTheta
arma::mat calcTheta(arma::vec rk, arma::mat ZZip, arma::mat x_kn, double p0);
RcppExport SEXP _dpfa_calcTheta(SEXP rkSEXP, SEXP ZZipSEXP, SEXP x_knSEXP, SEXP p0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type rk(rkSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ZZip(ZZipSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_kn(x_knSEXP);
    Rcpp::traits::input_parameter< double >::type p0(p0SEXP);
    rcpp_result_gen = Rcpp::wrap(calcTheta(rk, ZZip, x_kn, p0));
    return rcpp_result_gen;
END_RCPP
}
// calcW
arma::mat calcW(arma::vec sk, arma::mat ZZip, arma::mat C_k1n, double p0);
RcppExport SEXP _dpfa_calcW(SEXP skSEXP, SEXP ZZipSEXP, SEXP C_k1nSEXP, SEXP p0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type sk(skSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ZZip(ZZipSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C_k1n(C_k1nSEXP);
    Rcpp::traits::input_parameter< double >::type p0(p0SEXP);
    rcpp_result_gen = Rcpp::wrap(calcW(sk, ZZip, C_k1n, p0));
    return rcpp_result_gen;
END_RCPP
}
// mult_cpp
List mult_cpp(arma::mat X_mtn, arma::mat Psi, arma::mat Theta, arma::mat ZZip);
RcppExport SEXP _dpfa_mult_cpp(SEXP X_mtnSEXP, SEXP PsiSEXP, SEXP ThetaSEXP, SEXP ZZipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X_mtn(X_mtnSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ZZip(ZZipSEXP);
    rcpp_result_gen = Rcpp::wrap(mult_cpp(X_mtn, Psi, Theta, ZZip));
    return rcpp_result_gen;
END_RCPP
}
// crt_cpp
arma::rowvec crt_cpp(arma::mat X_kn, arma::vec rk);
RcppExport SEXP _dpfa_crt_cpp(SEXP X_knSEXP, SEXP rkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X_kn(X_knSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type rk(rkSEXP);
    rcpp_result_gen = Rcpp::wrap(crt_cpp(X_kn, rk));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dpfa_rztpois_single", (DL_FUNC) &_dpfa_rztpois_single, 1},
    {"_dpfa_rztpois_cpp", (DL_FUNC) &_dpfa_rztpois_cpp, 2},
    {"_dpfa_calcC_kn", (DL_FUNC) &_dpfa_calcC_kn, 4},
    {"_dpfa_calcTheta", (DL_FUNC) &_dpfa_calcTheta, 4},
    {"_dpfa_calcW", (DL_FUNC) &_dpfa_calcW, 4},
    {"_dpfa_mult_cpp", (DL_FUNC) &_dpfa_mult_cpp, 4},
    {"_dpfa_crt_cpp", (DL_FUNC) &_dpfa_crt_cpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_dpfa(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
