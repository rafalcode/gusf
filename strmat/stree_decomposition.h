
#ifndef _STREE_DECOMPOSITION_H_
#define _STREE_DECOMPOSITION_H_

#ifdef STRMAT
#include "stree_strmat.h"
#else
#include "stree.h"
#endif

typedef struct {
  SUFFIX_TREE tree;
  char *string, *raw_string;
  int length;

  int *prev, *block;
  int num_blocks, max_block_length;

#ifdef STATS
  unsigned int num_compares, num_edge_traversals;
#endif

} DECOMPOSITION_STRUCT, *DECOMPOSITION;


DECOMPOSITION_STRUCT *decomposition_prep(SUFFIX_TREE tree, char *string,
                                         char *raw_string, int length);
void decomposition_free(DECOMPOSITION_STRUCT *d);

void decomposition_print(DECOMPOSITION_STRUCT *d);
#define decomposition_get_num_blocks(D) ((D)->num_blocks)
#define decomposition_get_max_block_length(D) ((D)->max_block_length)
#define decomposition_get_block(D,I) ((D)->block[I])
#define decomposition_get_prev(D,I) ((D)->prev[I])

void lempel_ziv(DECOMPOSITION_STRUCT *d);
void lempel_ziv_original(DECOMPOSITION_STRUCT *d);
void lempel_ziv_nonoverlapping(DECOMPOSITION_STRUCT *d);

#endif
