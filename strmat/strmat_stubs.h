
int strmat_naive_match(STRING *pattern, STRING *text, int stats);
int strmat_bmbad_match(STRING *pattern, STRING *text, int stats);
int strmat_bmext_match(STRING *pattern, STRING *text, int stats);
int strmat_bmgood_match(STRING *pattern, STRING *text, int stats);
int strmat_bmextgood_match(STRING *pattern, STRING *text, int stats);
int strmat_kmp_sp_z_match(STRING *pattern, STRING *text, int stats);
int strmat_kmp_spprime_z_match(STRING *pattern, STRING *text, int stats);
int strmat_kmp_sp_orig_match(STRING *pattern, STRING *text, int stats);
int strmat_kmp_spprime_orig_match(STRING *pattern, STRING *text, int stats);
int strmat_ac_match(STRING **pattern_ary, int num_patterns, STRING *text,
                    int stats);
int strmat_bmset_naive_match(STRING **pattern_ary, int num_patterns,
                             STRING *text, int stats);
int strmat_bmset_badonly_match(STRING **patterns, int num_patterns,
                               STRING *text, int stats);
int strmat_bmset_2trees_match(STRING **patterns, int num_patterns,
                              STRING *text, int stats);
int strmat_bmset_1tree_match(STRING **patterns, int num_patterns,
                             STRING *text, int stats);
int strmat_z_build(STRING *str, int stats);
int strmat_z_match(STRING *pat, STRING *text, int stats);
