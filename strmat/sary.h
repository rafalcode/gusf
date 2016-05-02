
#ifndef _SARY_H_
#define _SARY_H_

typedef struct {
  char *S;
  int M, copyflag;

  int *Pos, *lcp, *lcp_leaves;

  int num_compares, num_tree_ops, num_lcp_ops;
} SARY_STRUCT;

SARY_STRUCT *sary_qsort_build(char *S, int M, int copyflag);
SARY_STRUCT *sary_zerkle_build(char *S, int M, int copyflag);
SARY_STRUCT *sary_stree_build(char *S, int M, int copyflag);
void sary_free(SARY_STRUCT *sary);

#endif
