
#ifndef _SARY_MATCH_H_
#define _SARY_MATCH_H_

#include "sary.h"

typedef enum { NAIVE_MATCH, MLR_MATCH, LCP_MATCH } SARY_MATCH_TYPE;
typedef struct {
  SARY_MATCH_TYPE type;

  char *T;
  int M, copyflag;

  SARY_STRUCT *sary;

  int i, iprime, N;

  int num_compares, search_depth;
} SARYMAT_STRUCT;

SARYMAT_STRUCT *sary_match_naive_prep(char *T, int M, int copyflag);
char *sary_match_naive_first(SARYMAT_STRUCT *smstruct, char *P, int N);
char *sary_match_naive_next(SARYMAT_STRUCT *smstruct);

SARYMAT_STRUCT *sary_match_mlr_prep(char *T, int M, int copyflag);
char *sary_match_mlr_first(SARYMAT_STRUCT *smstruct, char *P, int N);
char *sary_match_mlr_next(SARYMAT_STRUCT *smstruct);

SARYMAT_STRUCT *sary_match_lcp_prep(char *T, int M, int copyflag);
char *sary_match_lcp_first(SARYMAT_STRUCT *smstruct, char *P, int N);
char *sary_match_lcp_next(SARYMAT_STRUCT *smstruct);


void sary_match_free(SARYMAT_STRUCT *smstruct);

#endif
