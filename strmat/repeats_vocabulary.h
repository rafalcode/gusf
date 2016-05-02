
#ifndef _REPEATS_VOCABULARY_H_
#define _REPEATS_VOCABULARY_H_

#ifdef STRMAT
#include "stree_strmat.h"
#else
#include "stree.h"
#endif
#include "stree_decomposition.h"

typedef struct TANDEM {
  int len;
  struct TANDEM *next;
} tandem;

typedef struct {
  char *string, *raw_string;
  int length;

  SUFFIX_TREE tree;
  DECOMPOSITION decomposition;

  tandem *tandem_space, **tandems, **last;
  int next_tandem;
  int *tlens1, *tlens2, *dvector;

  int *PREF, *PREF2, *SUFF;

  unsigned int num_tandem_repeats;
  unsigned int num_primitive_tandem_repeats;
  unsigned int num_tandem_arrays;

  unsigned int num_tandem_repeat_occs;
  unsigned int num_primitive_tandem_repeat_occs;
  unsigned int num_tandem_array_occs;

#ifdef STATS
  unsigned int num_prep;
  unsigned int num_compares_for_tandem_repeats;
  unsigned int num_compares_for_primitive_tandem_repeats;
  unsigned int num_compares_for_tandem_arrays;
#endif

} VOCABULARY_STRUCT;


VOCABULARY_STRUCT *vocabulary_prep(SUFFIX_TREE tree,
                                   DECOMPOSITION_STRUCT *decomposition,
                                   char *string, char *raw_string, int length);
void vocabulary_free(VOCABULARY_STRUCT *vocabulary);

void vocabulary_find_tandem_repeats(VOCABULARY_STRUCT *v);
void vocabulary_find_primitive_tandem_repeats(VOCABULARY_STRUCT *v);
void vocabulary_find_tandem_arrays(VOCABULARY_STRUCT *v);

void vocabulary_write(VOCABULARY_STRUCT *v, char *type);
void vocabulary_count(VOCABULARY_STRUCT *v,
                      unsigned int *num, unsigned int *occ);

#endif
