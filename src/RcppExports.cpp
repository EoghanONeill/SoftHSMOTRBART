// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// get_children_cpp
arma::uvec get_children_cpp(arma::mat tree_mat, int parent);
RcppExport SEXP _SoftHSMOTRBART_get_children_cpp(SEXP tree_matSEXP, SEXP parentSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type tree_mat(tree_matSEXP);
    Rcpp::traits::input_parameter< int >::type parent(parentSEXP);
    rcpp_result_gen = Rcpp::wrap(get_children_cpp(tree_mat, parent));
    return rcpp_result_gen;
END_RCPP
}
// phi_app_soft
arma::mat phi_app_soft(arma::mat X_stand, arma::mat anc, double tau);
RcppExport SEXP _SoftHSMOTRBART_phi_app_soft(SEXP X_standSEXP, SEXP ancSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X_stand(X_standSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type anc(ancSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(phi_app_soft(X_stand, anc, tau));
    return rcpp_result_gen;
END_RCPP
}
// phi_app_softHS
arma::mat phi_app_softHS(arma::mat treemat, arma::mat internalmat, arma::mat xmat, arma::mat probmat);
RcppExport SEXP _SoftHSMOTRBART_phi_app_softHS(SEXP treematSEXP, SEXP internalmatSEXP, SEXP xmatSEXP, SEXP probmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type treemat(treematSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type internalmat(internalmatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xmat(xmatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type probmat(probmatSEXP);
    rcpp_result_gen = Rcpp::wrap(phi_app_softHS(treemat, internalmat, xmat, probmat));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SoftHSMOTRBART_get_children_cpp", (DL_FUNC) &_SoftHSMOTRBART_get_children_cpp, 2},
    {"_SoftHSMOTRBART_phi_app_soft", (DL_FUNC) &_SoftHSMOTRBART_phi_app_soft, 3},
    {"_SoftHSMOTRBART_phi_app_softHS", (DL_FUNC) &_SoftHSMOTRBART_phi_app_softHS, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_SoftHSMOTRBART(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}