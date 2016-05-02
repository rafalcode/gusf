
#ifndef _NAIVE_H_
#define _NAIVE_H_

typedef struct {
  char *S;
  int M, copyflag;
  int num_compares;
} NAIVE_STRUCT;

NAIVE_STRUCT *naive_prep(char *S, int M, int copyflag);
char *naive_search(NAIVE_STRUCT *node, char *T, int N, int initmatch);
void naive_free(NAIVE_STRUCT *node);

#endif
