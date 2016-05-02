
#ifndef _REPEATS_BIGPATH_H_
#define _REPEATS_BIGPATH_H_

#include <limits.h>

typedef struct BIGPATH_ENTRY {
  struct BIGPATH_ENTRY *next, *prev, *old_prev;
  int mark;
} bp_entry;

typedef struct {
  char *string, *raw_string;
  int length;
  char a[CHAR_MAX];
  int alpha_size;

  SUFFIX_TREE tree;

  bp_entry *entries, *entries2, ***list, ***last;
  STREE_NODE *big_child;

  unsigned int num_nonoverlapping_maximal_pairs;

#ifdef STATS
  unsigned int num_prep;
  unsigned int num_steps_for_lists;
  unsigned int num_compares;
#endif

} bp_struct;

bp_struct *bigpath_prep(SUFFIX_TREE tree,
                       char *string, char *raw_string, int length);
void bigpath_free(bp_struct *m);
void bigpath_find(bp_struct *m);

#endif
