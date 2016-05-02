
/* 
 * kmp.c
 *
 * The implementation of four variants of the Knuth-Morris-Pratt algorithm.
 * The variants are 1) using sp values and the Z-values preprocessing,
 * 2) using sp' values and the Z-values preprocessing, 3) using sp values
 * and the original preprocessing, and 4) using sp' values and the original
 * preprocessing.

 *  Notes:
 *   	8/94  -  Original implementation (Sean Davis)
 *      9/94  -  Redid implementation (James Knight)
 *      8/95  -  Rewrote the code with the new program organization and
 *               added all of the variations.
 *               (James Knight)
 *      3/96  -  Modularized the code (James Knight)
 *      7/96  -  Finished the modularization (James Knight)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef STRMAT
#include "z.h"
#endif
#include "kmp.h"


/*
 *
 * Knuth-Morris-Pratt preprocessing using sp values and Z-values algorithm.
 *
 */
KMP_STRUCT *kmp_sp_z_prep(char *P, int M, int copyflag)
{
  int i, j, *Z, *F, *sp, *spprime;
  char *buf;
#ifdef STRMAT
  Z_STRUCT *zstruct;
#endif
  KMP_STRUCT *node;

  P--;            /* Shift to make sequence be P[1],...,P[M] */

  /*
   * Allocate everything.
   */
  if ((node = malloc(sizeof(KMP_STRUCT))) == NULL)
    return NULL;

  memset(node, 0, sizeof(KMP_STRUCT));
  node->M = M;
  node->copyflag = copyflag;

  if ((F = node->F = malloc((M + 2) * sizeof(int))) == NULL) {
    free(node);
    return NULL;
  }
  memset(F, 0, (M + 2) * sizeof(int));

  if (!copyflag)
    node->P = P;
  else {
    if ((buf = malloc(M + 2)) == NULL) {
      kmp_free(node);
      return NULL;
    }

    buf[0] = buf[M+1] = '\0';
    memcpy(buf + 1, P + 1, M);
    node->P = buf;
  }

  if ((sp = malloc((M + 1) * sizeof(int))) == NULL ||
      (spprime = malloc((M + 1) * sizeof(int))) == NULL) {
    if (sp != NULL)
      free(sp);
    kmp_free(node);
    return NULL;
  }
  memset(sp, 0, (M + 1) * sizeof(int));
  memset(spprime, 0, (M + 1) * sizeof(int));

  /*
   * Compute the Z values for the pattern, either using the Z-values
   * algorithm in "z.c", or by brute force.
   */
  Z = NULL;
#ifdef STRMAT

  if ((zstruct = z_build(P+1, M, 0)) != NULL) {
    Z = zstruct->Z;
    zstruct->Z = NULL;

#ifdef STATS
    node->prep_compares += zstruct->prep_compares;
#endif

    z_free(zstruct);
  }

#else

  if ((Z = malloc((M + 1) * sizeof(int))) != NULL) {
    for (i=2; i <= M; i++) {
      count = 0;
      for (j=i,k=1; j <= M && P[j] == P[k]; j++,k++)
        count++;
      Z[i] = count;

#ifdef STATS
      node->prep_compares += count + 1;
#endif
    }
  }

#endif

  if (Z == NULL) {
    free(spprime);
    free(sp);
    kmp_free(node);
    return NULL;
  }

  /*
   * Compute the spprime values and then the sp values.
   */
  for (j=M; j >= 2; j--) {
    i = j + Z[j] - 1;
    spprime[i] = Z[j];

#ifdef STATS
    node->prep_compares++;
#endif
  }

  sp[M] = spprime[M];
  for (i=M-1; i >= 2; i--) {
    if (spprime[i] >= sp[i+1] - 1)
      sp[i] = spprime[i];
    else
      sp[i] = sp[i+1] - 1;

#ifdef STATS
    node->prep_compares++;
#endif
  }

  /*
   * Compute the F values.
   */
  F[1] = 1;
  for (i=2; i <= M + 1; i++) {
    F[i] = sp[i-1] + 1;

#ifdef STATS
    node->prep_compares++;
#endif
  }

  free(Z);
  free(spprime);
  free(sp);

  return node;
}


/*
 *
 * Knuth-Morris-Pratt preprocessing using spprime values and Z-values alg.
 *
 */
KMP_STRUCT *kmp_spprime_z_prep(char *P, int M, int copyflag)
{
  int i, j, *Z, *F, *spprime;
  char *buf;
#ifdef STRMAT
  Z_STRUCT *zstruct;
#endif
  KMP_STRUCT *node;

  P--;            /* Shift to make sequence be P[1],...,P[M] */

  /*
   * Allocate everything.
   */
  if ((node = malloc(sizeof(KMP_STRUCT))) == NULL)
    return NULL;

  memset(node, 0, sizeof(KMP_STRUCT));
  node->M = M;
  node->copyflag = copyflag;

  if ((F = node->F = malloc((M + 2) * sizeof(int))) == NULL) {
    free(node);
    return NULL;
  }
  memset(F, 0, (M + 2) * sizeof(int));

  if (!copyflag)
    node->P = P;
  else {
    if ((buf = malloc(M + 2)) == NULL) {
      kmp_free(node);
      return NULL;
    }

    buf[0] = buf[M+1] = '\0';
    memcpy(buf + 1, P + 1, M);
    node->P = buf;
  }

  if ((spprime = malloc((M + 1) * sizeof(int))) == NULL) {
    kmp_free(node);
    return NULL;
  }
  memset(spprime, 0, (M + 1) * sizeof(int));

  /*
   * Compute the Z values for the pattern, either using the Z-values
   * algorithm in "z.c", or by brute force.
   */
  Z = NULL;
#ifdef STRMAT

  if ((zstruct = z_build(P+1, M, 0)) != NULL) {
    Z = zstruct->Z;
    zstruct->Z = NULL;

#ifdef STATS
    node->prep_compares += zstruct->prep_compares;
#endif

    z_free(zstruct);
  }

#else

  if ((Z = malloc((M + 1) * sizeof(int))) != NULL) {
    for (i=2; i <= M; i++) {
      count = 0;
      for (j=i,k=1; j <= M && P[j] == P[k]; j++,k++)
        count++;
      Z[i] = count;

#ifdef STATS
      node->prep_compares += count + 1;
#endif
    }
  }

#endif

  if (Z == NULL) {
    free(spprime);
    kmp_free(node);
    return NULL;
  }

  /*
   * Compute the spprime values and then the F values.
   */
  for (j=M; j >= 2; j--) {
    i = j + Z[j] - 1;
    spprime[i] = Z[j];

#ifdef STATS
    node->prep_compares++;
#endif
  }

  F[1] = 1;
  for (i=2; i <= M + 1; i++) {
    F[i] = spprime[i-1] + 1;

#ifdef STATS
    node->prep_compares++;
#endif
  }

  free(Z);
  free(spprime);

  return node;
}


/*
 *
 * Original Knuth-Morris-Pratt preprocessing using sp values.
 *
 */
KMP_STRUCT *kmp_sp_orig_prep(char *P, int M, int copyflag)
{
  int i, v, *F, *sp;
  char x, *buf;
  KMP_STRUCT *node;

  P--;            /* Shift to make sequence be P[1],...,P[M] */

  /*
   * Allocate everything.
   */
  if ((node = malloc(sizeof(KMP_STRUCT))) == NULL)
    return NULL;

  memset(node, 0, sizeof(KMP_STRUCT));
  node->M = M;
  node->copyflag = copyflag;

  if ((F = node->F = malloc((M + 2) * sizeof(int))) == NULL) {
    free(node);
    return NULL;
  }
  memset(F, 0, (M + 2) * sizeof(int));

  if (!copyflag)
    node->P = P;
  else {
    if ((buf = malloc(M + 2)) == NULL) {
      kmp_free(node);
      return NULL;
    }

    buf[0] = buf[M+1] = '\0';
    memcpy(buf + 1, P + 1, M);
    node->P = buf;
  }

  if ((sp = malloc((M + 1) * sizeof(int))) == NULL) {
    kmp_free(node);
    return NULL;
  }
  memset(sp, 0, (M + 1) * sizeof(int));

  /*
   * Compute the sp values and then the F values.
   */
  sp[1] = 0;
  for (i=1; i <= M - 1; i++) {
    x = P[i+1];
    v = sp[i];
    while (v != 0 && P[v+1] != x) {
      v = sp[v];

#ifdef STATS
      node->prep_compares++;
#endif
    }
#ifdef STATS
    node->prep_compares++;
#endif

    if (P[v+1] == x)
      sp[i+1] = v + 1;
    else
      sp[i+1] = 0;

#ifdef STATS
    node->prep_compares++;
#endif
  }

  F[1] = 1;
  for (i=2; i <= M + 1; i++) {
    F[i] = sp[i-1] + 1;

#ifdef STATS
    node->prep_compares++;
#endif
  }

  free(sp);

  return node;
}


/*
 *
 * Original Knuth-Morris-Pratt preprocessing using spprime values.
 *
 */
KMP_STRUCT *kmp_spprime_orig_prep(char *P, int M, int copyflag)
{
  int i, v, *F, *sp, *spprime;
  char x, *buf;
  KMP_STRUCT *node;

  P--;            /* Shift to make sequence be P[1],...,P[M] */

  /*
   * Allocate everything.
   */
  if ((node = malloc(sizeof(KMP_STRUCT))) == NULL)
    return NULL;

  memset(node, 0, sizeof(KMP_STRUCT));
  node->M = M;
  node->copyflag = copyflag;

  if ((F = node->F = malloc((M + 2) * sizeof(int))) == NULL) {
    free(node);
    return NULL;
  }
  memset(F, 0, (M + 2) * sizeof(int));

  if (!copyflag)
    node->P = P;
  else {
    if ((buf = malloc(M + 2)) == NULL) {
      kmp_free(node);
      return NULL;
    }

    buf[0] = buf[M+1] = '\0';
    memcpy(buf + 1, P + 1, M);
    node->P = buf;
  }

  if ((sp = malloc((M + 1) * sizeof(int))) == NULL ||
      (spprime = malloc((M + 1) * sizeof(int))) == NULL) {
    if (sp != NULL)
      free(sp);
    kmp_free(node);
    return NULL;
  }
  memset(sp, 0, (M + 1) * sizeof(int));
  memset(spprime, 0, (M + 1) * sizeof(int));

  /*
   * Compute the sp values and then the spprime values.
   */
  sp[1] = 0;
  for (i=1; i <= M - 1; i++) {
    x = P[i+1];
    v = sp[i];
    while (v != 0 && P[v+1] != x) {
      v = sp[v];

#ifdef STATS
      node->prep_compares++;
#endif
    }
#ifdef STATS
    node->prep_compares++;
#endif

    if (P[v+1] == x)
      sp[i+1] = v + 1;
    else
      sp[i+1] = 0;

#ifdef STATS
    node->prep_compares++;
#endif
  }

  spprime[1] = 0;
  for (i=2; i <= M; i++) {
    v = sp[i];
    if (P[v+1] != P[i+1])
      spprime[i] = v;
    else
      spprime[i] = spprime[v];

#ifdef STATS
    node->prep_compares++;
#endif
  }

  /*
   * Compute the F values.
   */
  F[1] = 1;
  for (i=2; i <= M + 1; i++) {
    F[i] = spprime[i-1] + 1;

#ifdef STATS
    node->prep_compares++;
#endif
  }

  free(spprime);
  free(sp);

  return node;
}


/*
 *
 * Knuth-Morris-Pratt search algorithm.
 *
 */
char *kmp_search(KMP_STRUCT *node, char *T, int N, int initmatch)
{
  int p, c, M, *F;
  char *P;

  T--;            /* Shift to make sequence be T[1],...,T[N] */

  P = node->P;
  M = node->M;
  F = node->F;

  /*
   * Run the Knuth-Morris-Pratt algorithm over the sequence.
   * Stop at the first position to match the complete pattern.
   */
  if (!initmatch) {
    c = 1;
    p = 1;
  }
  else {
    c = M + 1;
    p = F[M+1];

#ifdef STATS
    node->total_shifts += M+1 - F[M+1];
    node->num_shifts++;
#endif
  }

  while (c <= N) {
    /*
     * Compare the suffix of P[p],...,P[M] to T[c],...,T[c+M-p].
     *
     * If stopped at the end of the pattern, return a match.  If
     * stopped by a mismatch at the beginning of the pattern, advance
     * c.  Otherwise, apply the failure function to p.
     */
#ifdef STATS

    while (p <= M && c <= N && P[p] == T[c]) {
      c++;
      p++;

      node->num_compares++;
    }
    node->num_compares++;

    if (p == M + 1)
      return &T[c-M];
    else if (p == 1) {
      c++;
      node->num_init_mismatch++;
    }
    else {
      node->total_shifts += p - F[p];
      p = F[p];
      node->num_shifts++;
    }

#else

    while (p <= M && c <= N && P[p] == T[c]) {
      c++;
      p++;
    }

    if (p == M + 1)
      return &T[c-M];
    else if (p == 1)
      c++;
    else
      p = F[p];

#endif
  }

  return NULL;
}


/*
 * kmp_free
 *
 * Free up the allocated KMP_STRUCT structure.
 *
 * Parameters:   node  -  a KMP_STRUCT structure
 *
 * Returns:  nothing.
 */
void kmp_free(KMP_STRUCT *node)
{
  if (node == NULL)
    return;

  if (node->copyflag && node->P != NULL)
    free(node->P);
  if (node->F != NULL)
    free(node->F);
  free(node);
}


