
#ifndef _REPEATS_LINEAR_OCCS_H_
#define _REPEATS_LINEAR_OCCS_H_

#ifdef STRMAT
#include "stree_strmat.h"
#else
#include "stree.h"
#endif
#include "stree_decomposition.h"
#include "repeats_vocabulary.h"

VOCABULARY_STRUCT *linear_occs_prep(SUFFIX_TREE tree,
                                    DECOMPOSITION_STRUCT *decomposition,
                                    char *string, char *raw_string,
                                    int length);
void linear_occs_free(VOCABULARY_STRUCT *vocabulary);

void linear_occs_find_tandem_repeats(VOCABULARY_STRUCT *v);
void linear_occs_write(VOCABULARY_STRUCT *v, char *type);

#endif
