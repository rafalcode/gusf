
#ifndef _BMSET_NAIVE_H_
#define _BMSET_NAIVE_H_

#include "bm.h"

typedef struct {
  int num_patterns, patterns_size;
  BM_STRUCT **patterns;
  int *ids;

  char *T;
  int N, initflag, errorflag;
  int startflag, endflag, output;
  char **matches;

  int prep_compares, num_compares;
  int num_shifts;
} BMSET_NAIVE_STRUCT;


BMSET_NAIVE_STRUCT *bmset_naive_alloc(void);
int bmset_naive_add_string(BMSET_NAIVE_STRUCT *node, char *P, int M, int id,
                           int copyflag);
int bmset_naive_del_string(BMSET_NAIVE_STRUCT *node, char *P, int M, int id);
void bmset_naive_search_init(BMSET_NAIVE_STRUCT *node, char *T, int N);
char *bmset_naive_search(BMSET_NAIVE_STRUCT *node, int *length_out,
                         int *id_out);
void bmset_naive_free(BMSET_NAIVE_STRUCT *node);

#endif
