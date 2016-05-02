
#ifndef _KMP_H_
#define _KMP_H_

typedef struct {
  char *P;
  int M, copyflag;
  int *F;

  int prep_compares, num_compares;
  int num_shifts, total_shifts, num_init_mismatch;
} KMP_STRUCT;


KMP_STRUCT *kmp_sp_z_prep(char *S, int M, int copyflag);
KMP_STRUCT *kmp_spprime_z_prep(char *S, int M, int copyflag);
KMP_STRUCT *kmp_sp_orig_prep(char *S, int M, int copyflag);
KMP_STRUCT *kmp_spprime_orig_prep(char *S, int M, int copyflag);
char *kmp_search(KMP_STRUCT *node, char *T, int N, int initmatch);
void kmp_free(KMP_STRUCT *node);

#endif
