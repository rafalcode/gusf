
int strmat_repeats_primitives(STRING *string, int print_stats);
int strmat_repeats_supermax(STRING *string, int min_percent, int min_length);
int strmat_repeats_nonoverlapping(STRING *string, int print_stats);
int strmat_repeats_bigpath(STRING *string, int build_policy,
                          int build_threshold, int print_stats);
int strmat_repeats_tandem(STRING *string, int build_policy,
                          int build_threshold, int print_stats);
int strmat_repeats_vocabulary(STRING *string, int build_policy,
                              int build_threshold, int print_stats, char mode);
int strmat_repeats_linear_occs(STRING *string, int build_policy,
                               int build_threshold, int print_stats, char mode);

