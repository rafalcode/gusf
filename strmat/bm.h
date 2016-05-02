
#ifndef _BM_H_
#define _BM_H_

typedef enum {
  BM_BAD, BM_EXT, BM_GOOD, BM_EXTGOOD
} BMALG_TYPE;

typedef struct {
  BMALG_TYPE type;

  char *P;
  int M, copyflag;
  int *R, *Rnext, *Lprime, *lprime;

  int prep_compares, num_compares;
  int num_shifts, shift_cost;
  int num_init_mismatch;
} BM_STRUCT;


BM_STRUCT *bmbad_prep(char *S, int M, int copyflag);
char *bmbad_search(BM_STRUCT *node, char *T, int N, int initmatch);
BM_STRUCT *bmext_prep(char *S, int M, int copyflag);
char *bmext_search(BM_STRUCT *node, char *T, int N, int initmatch);
BM_STRUCT *bmgood_prep(char *S, int M, int copyflag);
char *bmgood_search(BM_STRUCT *node, char *T, int N, int initmatch);
BM_STRUCT *bmextgood_prep(char *S, int M, int copyflag);
char *bmextgood_search(BM_STRUCT *node, char *T, int N, int initmatch);
void bm_free(BM_STRUCT *node);

#endif
