
#ifndef _Z_H_
#define _Z_H_

typedef struct {
  char *S;
  int M, copyflag;
  int *Z;
  int prep_compares, num_compares;
} Z_STRUCT;

Z_STRUCT *z_build(char *S, int M, int copyflag);
char *z_search(Z_STRUCT *node, char *T, int N, int initmatch);
void z_free(Z_STRUCT *node);

#endif
