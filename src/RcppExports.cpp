// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// eval_cpp
Rcpp::List eval_cpp(const std::string& rhs, const Rcpp::IntegerVector& ances, const Rcpp::NumericMatrix& states, const Rcpp::NumericMatrix& forTime, const Rcpp::RObject& lambdas, const Rcpp::NumericVector& mus, const Rcpp::NumericMatrix& Q, const std::string& method, double atol, double rtol, bool is_complete_tree, size_t num_steps);
RcppExport SEXP _secsse_eval_cpp(SEXP rhsSEXP, SEXP ancesSEXP, SEXP statesSEXP, SEXP forTimeSEXP, SEXP lambdasSEXP, SEXP musSEXP, SEXP QSEXP, SEXP methodSEXP, SEXP atolSEXP, SEXP rtolSEXP, SEXP is_complete_treeSEXP, SEXP num_stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type rhs(rhsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type ances(ancesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type states(statesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type forTime(forTimeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::RObject& >::type lambdas(lambdasSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type mus(musSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type atol(atolSEXP);
    Rcpp::traits::input_parameter< double >::type rtol(rtolSEXP);
    Rcpp::traits::input_parameter< bool >::type is_complete_tree(is_complete_treeSEXP);
    Rcpp::traits::input_parameter< size_t >::type num_steps(num_stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(eval_cpp(rhs, ances, states, forTime, lambdas, mus, Q, method, atol, rtol, is_complete_tree, num_steps));
    return rcpp_result_gen;
END_RCPP
}
// calc_ll_cpp
Rcpp::List calc_ll_cpp(const std::string& rhs, const Rcpp::IntegerVector& ances, const Rcpp::NumericMatrix& states, const Rcpp::NumericMatrix& forTime, const Rcpp::RObject& lambdas, const Rcpp::NumericVector& mus, const Rcpp::NumericMatrix& Q, const std::string& method, double atol, double rtol, bool is_complete_tree, bool see_states);
RcppExport SEXP _secsse_calc_ll_cpp(SEXP rhsSEXP, SEXP ancesSEXP, SEXP statesSEXP, SEXP forTimeSEXP, SEXP lambdasSEXP, SEXP musSEXP, SEXP QSEXP, SEXP methodSEXP, SEXP atolSEXP, SEXP rtolSEXP, SEXP is_complete_treeSEXP, SEXP see_statesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type rhs(rhsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type ances(ancesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type states(statesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type forTime(forTimeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::RObject& >::type lambdas(lambdasSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type mus(musSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type atol(atolSEXP);
    Rcpp::traits::input_parameter< double >::type rtol(rtolSEXP);
    Rcpp::traits::input_parameter< bool >::type is_complete_tree(is_complete_treeSEXP);
    Rcpp::traits::input_parameter< bool >::type see_states(see_statesSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_ll_cpp(rhs, ances, states, forTime, lambdas, mus, Q, method, atol, rtol, is_complete_tree, see_states));
    return rcpp_result_gen;
END_RCPP
}
// ct_condition_cpp
Rcpp::NumericVector ct_condition_cpp(const std::string rhs, const Rcpp::NumericVector& state, const double t, const Rcpp::RObject& lambdas, const Rcpp::NumericVector& mus, const Rcpp::NumericMatrix& Q, const std::string& method, double atol, double rtol);
RcppExport SEXP _secsse_ct_condition_cpp(SEXP rhsSEXP, SEXP stateSEXP, SEXP tSEXP, SEXP lambdasSEXP, SEXP musSEXP, SEXP QSEXP, SEXP methodSEXP, SEXP atolSEXP, SEXP rtolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type rhs(rhsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type state(stateSEXP);
    Rcpp::traits::input_parameter< const double >::type t(tSEXP);
    Rcpp::traits::input_parameter< const Rcpp::RObject& >::type lambdas(lambdasSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type mus(musSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type atol(atolSEXP);
    Rcpp::traits::input_parameter< double >::type rtol(rtolSEXP);
    rcpp_result_gen = Rcpp::wrap(ct_condition_cpp(rhs, state, t, lambdas, mus, Q, method, atol, rtol));
    return rcpp_result_gen;
END_RCPP
}
// secsse_sim_cpp
Rcpp::List secsse_sim_cpp(const std::vector<double>& m_R, const Rcpp::List& lambdas_R, const Rcpp::NumericMatrix& q_R, double max_time, double max_species, bool max_species_extant, double min_species, const std::vector<double>& init_states, std::string condition, int num_concealed_states, bool non_extinction, bool verbose, int max_tries, int seed, const std::vector<double>& conditioning_vec, bool return_tree_size_hist);
RcppExport SEXP _secsse_secsse_sim_cpp(SEXP m_RSEXP, SEXP lambdas_RSEXP, SEXP q_RSEXP, SEXP max_timeSEXP, SEXP max_speciesSEXP, SEXP max_species_extantSEXP, SEXP min_speciesSEXP, SEXP init_statesSEXP, SEXP conditionSEXP, SEXP num_concealed_statesSEXP, SEXP non_extinctionSEXP, SEXP verboseSEXP, SEXP max_triesSEXP, SEXP seedSEXP, SEXP conditioning_vecSEXP, SEXP return_tree_size_histSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type m_R(m_RSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type lambdas_R(lambdas_RSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type q_R(q_RSEXP);
    Rcpp::traits::input_parameter< double >::type max_time(max_timeSEXP);
    Rcpp::traits::input_parameter< double >::type max_species(max_speciesSEXP);
    Rcpp::traits::input_parameter< bool >::type max_species_extant(max_species_extantSEXP);
    Rcpp::traits::input_parameter< double >::type min_species(min_speciesSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type init_states(init_statesSEXP);
    Rcpp::traits::input_parameter< std::string >::type condition(conditionSEXP);
    Rcpp::traits::input_parameter< int >::type num_concealed_states(num_concealed_statesSEXP);
    Rcpp::traits::input_parameter< bool >::type non_extinction(non_extinctionSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type max_tries(max_triesSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type conditioning_vec(conditioning_vecSEXP);
    Rcpp::traits::input_parameter< bool >::type return_tree_size_hist(return_tree_size_histSEXP);
    rcpp_result_gen = Rcpp::wrap(secsse_sim_cpp(m_R, lambdas_R, q_R, max_time, max_species, max_species_extant, min_species, init_states, condition, num_concealed_states, non_extinction, verbose, max_tries, seed, conditioning_vec, return_tree_size_hist));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_secsse_eval_cpp", (DL_FUNC) &_secsse_eval_cpp, 12},
    {"_secsse_calc_ll_cpp", (DL_FUNC) &_secsse_calc_ll_cpp, 12},
    {"_secsse_ct_condition_cpp", (DL_FUNC) &_secsse_ct_condition_cpp, 9},
    {"_secsse_secsse_sim_cpp", (DL_FUNC) &_secsse_secsse_sim_cpp, 16},
    {NULL, NULL, 0}
};

RcppExport void R_init_secsse(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}