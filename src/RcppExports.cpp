// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// erdos
List erdos(int n, double prob);
RcppExport SEXP diseaseNetworks_erdos(SEXP nSEXP, SEXP probSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< int >::type n(nSEXP );
        Rcpp::traits::input_parameter< double >::type prob(probSEXP );
        List __result = erdos(n, prob);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// waxman
List waxman(int n, double a, double d);
RcppExport SEXP diseaseNetworks_waxman(SEXP nSEXP, SEXP aSEXP, SEXP dSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< int >::type n(nSEXP );
        Rcpp::traits::input_parameter< double >::type a(aSEXP );
        Rcpp::traits::input_parameter< double >::type d(dSEXP );
        List __result = waxman(n, a, d);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// apply_cluster
NumericMatrix apply_cluster(NumericMatrix nodes, NumericMatrix F, double f);
RcppExport SEXP diseaseNetworks_apply_cluster(SEXP nodesSEXP, SEXP FSEXP, SEXP fSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type nodes(nodesSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type F(FSEXP );
        Rcpp::traits::input_parameter< double >::type f(fSEXP );
        NumericMatrix __result = apply_cluster(nodes, F, f);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// keeling
List keeling(int n, int F, double f, double a, double d);
RcppExport SEXP diseaseNetworks_keeling(SEXP nSEXP, SEXP FSEXP, SEXP fSEXP, SEXP aSEXP, SEXP dSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< int >::type n(nSEXP );
        Rcpp::traits::input_parameter< int >::type F(FSEXP );
        Rcpp::traits::input_parameter< double >::type f(fSEXP );
        Rcpp::traits::input_parameter< double >::type a(aSEXP );
        Rcpp::traits::input_parameter< double >::type d(dSEXP );
        List __result = keeling(n, F, f, a, d);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// sim_epidemic
IntegerVector sim_epidemic(IntegerMatrix edges, IntegerVector node_states, double gamma, double Beta, int trace);
RcppExport SEXP diseaseNetworks_sim_epidemic(SEXP edgesSEXP, SEXP node_statesSEXP, SEXP gammaSEXP, SEXP BetaSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type edges(edgesSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type node_states(node_statesSEXP );
        Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP );
        Rcpp::traits::input_parameter< double >::type Beta(BetaSEXP );
        Rcpp::traits::input_parameter< int >::type trace(traceSEXP );
        IntegerVector __result = sim_epidemic(edges, node_states, gamma, Beta, trace);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// finalSize
NumericVector finalSize(double Beta, int n);
RcppExport SEXP diseaseNetworks_finalSize(SEXP BetaSEXP, SEXP nSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< double >::type Beta(BetaSEXP );
        Rcpp::traits::input_parameter< int >::type n(nSEXP );
        NumericVector __result = finalSize(Beta, n);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// Prod
double Prod(NumericVector a);
RcppExport SEXP diseaseNetworks_Prod(SEXP aSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP );
        double __result = Prod(a);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// myWhich
int myWhich(NumericVector a, double b);
RcppExport SEXP diseaseNetworks_myWhich(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP );
        Rcpp::traits::input_parameter< double >::type b(bSEXP );
        int __result = myWhich(a, b);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// mle
double mle(NumericVector Ne, NumericVector Beta);
RcppExport SEXP diseaseNetworks_mle(SEXP NeSEXP, SEXP BetaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type Ne(NeSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Beta(BetaSEXP );
        double __result = mle(Ne, Beta);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
