#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _xyz2_absolute_covariates(SEXP, SEXP, SEXP, SEXP);
extern SEXP _xyz2_absolute_covariates_pairs(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _xyz2_calculate_residuals(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _xyz2_calculate_xbeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _xyz2_clean_all_effects(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _xyz2_clean_pairs(SEXP);
extern SEXP _xyz2_colsum_index(SEXP, SEXP);
extern SEXP _xyz2_create_lambda_sequence(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _xyz2_equalpairs(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _xyz2_estimate_background_interaction_frequency(SEXP, SEXP, SEXP);
extern SEXP _xyz2_find_strongest_pairs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _xyz2_gaussiglmnet(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _xyz2_interaction_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _xyz2_interaction_search_low_level(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _xyz2_iterate(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _xyz2_naive_interaction_search(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _xyz2_order_vector(SEXP, SEXP);
extern SEXP _xyz2_prod_matrix_vector(SEXP, SEXP);
extern SEXP _xyz2_projected_equal_pairs(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _xyz2_sample_int_replace(SEXP, SEXP);
extern SEXP _xyz2_sample_uniform(SEXP, SEXP);
extern SEXP _xyz2_scale_intr(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _xyz2_scale_main(SEXP, SEXP, SEXP);
extern SEXP _xyz2_scan_intr_effects(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _xyz2_scan_main_effects(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _xyz2_soft_threshold(SEXP, SEXP, SEXP, SEXP);
extern SEXP _xyz2_sort_using_order_intmat(SEXP, SEXP);
extern SEXP _xyz2_sort_using_order_numvec(SEXP, SEXP);
extern SEXP _xyz2_translate_to_binary(SEXP, SEXP);
extern SEXP _xyz2_update_intr_final(SEXP, SEXP);
extern SEXP _xyz2_update_intr_vars(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _xyz2_warm_start(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_xyz2_absolute_covariates",                       (DL_FUNC) &_xyz2_absolute_covariates,                        4},
    {"_xyz2_absolute_covariates_pairs",                 (DL_FUNC) &_xyz2_absolute_covariates_pairs,                  5},
    {"_xyz2_calculate_residuals",                       (DL_FUNC) &_xyz2_calculate_residuals,                       11},
    {"_xyz2_calculate_xbeta",                           (DL_FUNC) &_xyz2_calculate_xbeta,                           11},
    {"_xyz2_clean_all_effects",                         (DL_FUNC) &_xyz2_clean_all_effects,                          5},
    {"_xyz2_clean_pairs",                               (DL_FUNC) &_xyz2_clean_pairs,                                1},
    {"_xyz2_colsum_index",                              (DL_FUNC) &_xyz2_colsum_index,                               2},
    {"_xyz2_create_lambda_sequence",                    (DL_FUNC) &_xyz2_create_lambda_sequence,                     5},
    {"_xyz2_equalpairs",                                (DL_FUNC) &_xyz2_equalpairs,                                 5},
    {"_xyz2_estimate_background_interaction_frequency", (DL_FUNC) &_xyz2_estimate_background_interaction_frequency,  3},
    {"_xyz2_find_strongest_pairs",                      (DL_FUNC) &_xyz2_find_strongest_pairs,                       6},
    {"_xyz2_gaussiglmnet",                              (DL_FUNC) &_xyz2_gaussiglmnet,                              10},
    {"_xyz2_interaction_search",                        (DL_FUNC) &_xyz2_interaction_search,                         8},
    {"_xyz2_interaction_search_low_level",              (DL_FUNC) &_xyz2_interaction_search_low_level,               7},
    {"_xyz2_iterate",                                   (DL_FUNC) &_xyz2_iterate,                                   15},
    {"_xyz2_naive_interaction_search",                  (DL_FUNC) &_xyz2_naive_interaction_search,                   5},
    {"_xyz2_order_vector",                              (DL_FUNC) &_xyz2_order_vector,                               2},
    {"_xyz2_prod_matrix_vector",                        (DL_FUNC) &_xyz2_prod_matrix_vector,                         2},
    {"_xyz2_projected_equal_pairs",                     (DL_FUNC) &_xyz2_projected_equal_pairs,                      5},
    {"_xyz2_sample_int_replace",                        (DL_FUNC) &_xyz2_sample_int_replace,                         2},
    {"_xyz2_sample_uniform",                            (DL_FUNC) &_xyz2_sample_uniform,                             2},
    {"_xyz2_scale_intr",                                (DL_FUNC) &_xyz2_scale_intr,                                 5},
    {"_xyz2_scale_main",                                (DL_FUNC) &_xyz2_scale_main,                                 3},
    {"_xyz2_scan_intr_effects",                         (DL_FUNC) &_xyz2_scan_intr_effects,                         13},
    {"_xyz2_scan_main_effects",                         (DL_FUNC) &_xyz2_scan_main_effects,                         11},
    {"_xyz2_soft_threshold",                            (DL_FUNC) &_xyz2_soft_threshold,                             4},
    {"_xyz2_sort_using_order_intmat",                   (DL_FUNC) &_xyz2_sort_using_order_intmat,                    2},
    {"_xyz2_sort_using_order_numvec",                   (DL_FUNC) &_xyz2_sort_using_order_numvec,                    2},
    {"_xyz2_translate_to_binary",                       (DL_FUNC) &_xyz2_translate_to_binary,                        2},
    {"_xyz2_update_intr_final",                         (DL_FUNC) &_xyz2_update_intr_final,                          2},
    {"_xyz2_update_intr_vars",                          (DL_FUNC) &_xyz2_update_intr_vars,                           5},
    {"_xyz2_warm_start",                                (DL_FUNC) &_xyz2_warm_start,                                 5},
    {NULL, NULL, 0}
};

void R_init_xyz2(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}