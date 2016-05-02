/* 
 * bm.c
 *
 * The implementation of four variants of the Boyer-Moore algorithm.  The
 * variants are 1) using only the simple bad character rule, 2) using only
 * the extended bad character rule, 3) using the simple bad character rule
 * and the strong good suffix rule, and 4) using the extended bad character
 * rule and the strong good suffix rule.
 *
 * NOTES:
 *    7/94  -  Original Implementation.  (James Knight)
 *    8/94  -  Removed the DEBUG tests.  (James Knight)
 *    8/94  -  Created new procedures for searching while collecting
 *              statistics.  (James Knight)
 *    8/94  -  Fixed bug which would cause division by zero of 
 *             problem->bm_num_shifts in bm_print_stats.  (James Knight)
 *    8/94  -  Change the print routines to use mprintf.  (James Knight)
 *    9/94  -  Converted it to ANSI format.  (James Knight)
 *    9/94  -  Redid the print routines to conform to print policies.
 *             (James Knight)
 *    6/95  -  Took out the optimized version, and fixed up the structure.
 *             (James Knight)
 *    7/95  -  Rewrote the code with the new program organization and
 *             added all of the variations.
 *             (James Knight)
 *    3/96  -  Modularized the code (James Knight)
 *    7/96  -  Finished the modularization  (James Knight)
 */

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#ifdef STRMAT
#include "z.h"
#endif
#include "bm.h"


/*
 *
 * Boyer-Moore searching using only the bad character rule.
 *
 */
BM_STRUCT *bmbad_prep(char *P, int M, int copyflag)
{
  int i, *R;
  char *buf;
  BM_STRUCT *node;

  P--;            /* Shift to make sequence be P[1],...,P[M] */

  /*
   * Allocate everything.
   */
  if ((node = malloc(sizeof(BM_STRUCT))) == NULL)
    return NULL;

  memset(node, 0, sizeof(BM_STRUCT));
  node->type = BM_BAD;
  node->copyflag = copyflag;
  node->M = M;

  if ((R = node->R = malloc(128 * sizeof(int))) == NULL) {
    bm_free(node);
    return NULL;
  }
  memset(R, 0, 128 * sizeof(int));

  if (!copyflag)
    node->P = P;
  else {
    if ((buf = malloc(M + 2)) == NULL) {
      bm_free(node);
      return NULL;
    }

    buf[0] = buf[M+1] = '\0';
    memcpy(buf + 1, P + 1, M);
    node->P = buf;
  }

  /*
   * Simple bad character rule.
   */
  for (i=1; i <= M; i++) {
    R[(int) P[i]] = i;

#ifdef STATS
    node->prep_compares++;
#endif
  }

  return node;
}


char *bmbad_search(BM_STRUCT *node, char *T, int N, int initmatch)
{
  int i, k, h, M, bshift, *R;
  char *P;

  if (node->type != BM_BAD)
    return NULL;

  T--;            /* Shift to make sequence be T[1],...,T[M] */

  P = node->P;
  M = node->M;
  R = node->R;

  /*
   * Run the Boyer-Moore algorithm over the sequence.  Stop at the first
   * position to match the complete pattern.
   */
  k = M;
  if (initmatch) {
    /*
     * If there's a complete match, just shift by one
     * (there's no bad character to shift by).
     */
    k++;

#ifdef STATS
    node->shift_cost++;
    node->num_shifts++;
#endif
  }

  while (k <= N) {
    /*
     * Compare the pattern to the string ending at k (in T).
     */
#ifdef STATS

    i = M;
    h = k;
    while (i > 0 && P[i] == T[h]) {
      i--;
      h--;
      node->num_compares++;
    }
    node->num_compares++;

#else

    for (i=M,h=k; i > 0 && P[i] == T[h]; i--,h--) 
      ;

#endif

    /*
     * Either return on a complete match to the pattern, or
     * use the bad character rule to perform a shift.
     */
    if (i == 0)
      return &T[k-M+1];
    else {
#ifdef STATS

      if (i == M)
        node->num_init_mismatch++;

      bshift = i - R[(int) T[h]];
      if (bshift < 1)
        bshift = 1;
      node->shift_cost++;

      k += bshift;
      node->num_shifts++;

#else

      bshift = i - R[(int) T[h]];
      if (bshift < 1)
        bshift = 1;

      k += bshift;

#endif
    }
  }

  return NULL;
}

 
/*
 *
 * Boyer-Moore searching using only the extended bad character rule.
 *
 */
BM_STRUCT *bmext_prep(char *P, int M, int copyflag)
{
  int i, *R, *Rnext;
  char *buf;
  BM_STRUCT *node;

  P--;            /* Shift to make sequence be P[1],...,P[M] */

  /*
   * Allocate everything.
   */
  if ((node = malloc(sizeof(BM_STRUCT))) == NULL)
    return NULL;

  memset(node, 0, sizeof(BM_STRUCT));
  node->type = BM_EXT;
  node->copyflag = copyflag;
  node->M = M;

  if ((R = node->R = malloc(128 * sizeof(int))) == NULL) {
    bm_free(node);
    return NULL;
  }
  memset(R, 0, 128 * sizeof(int));

  if ((Rnext = node->Rnext = malloc((M + 1) * sizeof(int))) == NULL) {
    bm_free(node);
    return NULL;
  }
  memset(Rnext, 0, (M + 1) * sizeof(int));

  if (!copyflag)
    node->P = P;
  else {
    if ((buf = malloc(M + 2)) == NULL) {
      bm_free(node);
      return NULL;
    }

    buf[0] = buf[M+1] = '\0';
    memcpy(buf + 1, P + 1, M);
    node->P = buf;
  }

  /*
   * Extended bad character rule.  R holds the heads of the lists and
   * Rnext holds the rest of the values.
   */
  for (i=1; i <= M; i++) {
    Rnext[i] = R[(int) P[i]];
    R[(int) P[i]] = i;

#ifdef STATS
    node->prep_compares++;
#endif
  }

  return node;
}


char *bmext_search(BM_STRUCT *node, char *T, int N, int initmatch)
{
  int i, k, h, M, pos, bshift, *R, *Rnext;
  char *P;

  if (node->type != BM_EXT)
    return NULL;

  T--;            /* Shift to make sequence be T[1],...,T[M] */

  P = node->P;
  M = node->M;
  R = node->R;
  Rnext = node->Rnext;

  /*
   * Run the Boyer-Moore algorithm over the sequence.  Stop at the first
   * position to match the complete pattern.
   */
  if (!initmatch)
    k = M;
  else {
    /*
     * If there's a complete match, just shift by one
     * (there's no bad character to shift by).
     */
    k = M + 1;

#ifdef STATS
    node->shift_cost++;
    node->num_shifts++;
#endif
  }

  while (k <= N) {
    /*
     * Compare the pattern to the string ending at k (in T).
     */
#ifdef STATS

    i = M;
    h = k;
    while (i > 0 && P[i] == T[h]) {
      i--;
      h--;
      node->num_compares++;
    }
    node->num_compares++;

#else

    for (i=M,h=k; i > 0 && P[i] == T[h]; i--,h--) 
      ;

#endif

    /*
     * Either return on a complete match to the pattern, or
     * use the extended bad character rule to perform a shift.
     */
    if (i == 0)
      return &T[k-M+1];
    else {
#ifdef STATS

      if (i == M)
        node->num_init_mismatch++;

      pos = R[(int) T[h]];
      node->shift_cost++;

      while (pos >= i) {
        pos = Rnext[pos];
        node->shift_cost++;
      }
      bshift = i - pos;

      k += bshift;
      node->num_shifts++;

#else

      pos = R[(int) T[h]];
      while (pos >= i)
        pos = Rnext[pos];
      bshift = i - pos;

      k += bshift;

#endif
    }
  }

  return NULL;
}

 
/*
 *
 * Boyer-Moore searching using the bad character rule and the strong good
 * suffix rule.
 *
 */
BM_STRUCT *bmgood_prep(char *P, int M, int copyflag)
{
  int i, j, *R, *Z, *Lprime, *lprime;
  char *buf;
#ifdef STRMAT
  char *Pbar;
  Z_STRUCT *zstruct;
#else
  int k, count;
#endif
  BM_STRUCT *node;

  P--;            /* Shift to make sequence be P[1],...,P[M] */

  /*
   * Allocate everything.
   */
  if ((node = malloc(sizeof(BM_STRUCT))) == NULL)
    return NULL;

  memset(node, 0, sizeof(BM_STRUCT));
  node->type = BM_GOOD;
  node->M = M;
  node->copyflag = copyflag;

  if ((R = node->R = malloc(128 * sizeof(int))) == NULL) {
    bm_free(node);
    return NULL;
  }
  memset(R, 0, 128 * sizeof(int));

  if ((Lprime = node->Lprime = malloc((M + 1) * sizeof(int))) == NULL) {
    bm_free(node);
    return NULL;
  }
  memset(Lprime, 0, (M + 1) * sizeof(int));

  if ((lprime = node->lprime = malloc((M + 1) * sizeof(int))) == NULL) {
    bm_free(node);
    return NULL;
  }
  memset(lprime, 0, (M + 1) * sizeof(int));

  if (!copyflag) 
    node->P = P;
  else {
    if ((buf = malloc(M + 2)) == NULL) {
      bm_free(node);
      return NULL;
    }

    buf[0] = buf[M+1] = '\0';
    memcpy(buf + 1, P + 1, M);
    node->P = buf;
  }

  /*
   * Simple bad character rule.
   */
  for (i=1; i <= M; i++) {
    R[(int) P[i]] = i;

#ifdef STATS
    node->prep_compares++;
#endif
  }

  /*
   * Good Suffix rule.
   *
   * First, compute the Z values for the reversed string, either using the
   * Z-values algorithm in "z.c", if available, or using a brute force
   * computation.
   */
  Z = NULL;
#ifdef STRMAT

  zstruct = NULL;
  if ((Pbar = malloc(M + 1)) != NULL) {
    for (i=0; i < M; i++)
      Pbar[i] = P[M-i];
    Pbar[M] = '\0';

    if ((zstruct = z_build(Pbar, M, 0)) != NULL) {
      Z = zstruct->Z;

#ifdef STATS
      node->prep_compares += zstruct->prep_compares;
#endif
    }

    free(Pbar);
  }

#else

  if ((Z = malloc((M + 1) * sizeof(int))) != NULL) {
    for (i=1; i < M; i++) {
      count = 0;
      for (j=i,k=M; j > 0 && P[j] == P[k]; j--,k--)
        count++;
      Z[i] = count;

#ifdef STATS
      node->prep_compares += count + 1;
#endif
    }
  }

#endif

  if (Z == NULL) {
    bm_free(node);
    return NULL;
  }

  /*
   * Next, compute the L' and l' values from the Z values.
   */
  for (j=1; j < M; j++) {
    i = M - Z[M-j+1] + 1;
    Lprime[i] = j;

#ifdef STATS
    node->prep_compares++;
#endif
  }

  lprime[M] = (P[1] == P[M] ? 1 : 0);
#ifdef STATS
  node->prep_compares++;
#endif

  for (i=M-1,j=2; i >= 2; i--,j++) {
    if (Z[M-j+1] == j)
      lprime[i] = j;
    else
      lprime[i] = lprime[i+1];

#ifdef STATS
    node->prep_compares++;
#endif
  }

#ifdef STRMAT
  z_free(zstruct);
#else
  free(Z);
#endif

  return node;
}


char *bmgood_search(BM_STRUCT *node, char *T, int N, int initmatch)
{
  int i, k, h, M, gshift, bshift, *R, *Lprime, *lprime;
  char *P;

  if (node->type != BM_GOOD)
    return NULL;

  T--;            /* Shift to make sequence be T[1],...,T[N] */

  P = node->P;
  M = node->M;
  R = node->R;
  Lprime = node->Lprime;
  lprime = node->lprime;

  /*
   * Run the Boyer-Moore algorithm over the sequence.  Stop at the first
   * position to match the complete pattern.
   */
  k = M;
  if (initmatch) {
    /*
     * Perform a good shift to shift past the match.
     */
    gshift = M - lprime[2];
    k += gshift;

#ifdef STATS
    node->shift_cost++;
    node->num_shifts++;
#endif
  }

  while (k <= N) {
    /*
     * Compare the pattern to the string ending at k (in T).
     */
#ifdef STATS

    i = M;
    h = k;
    while (i > 0 && P[i] == T[h]) {
      i--;
      h--;
      node->num_compares++;
    }
    node->num_compares++;

#else

    for (i=M,h=k; i > 0 && P[i] == T[h]; i--,h--) 
      ;

#endif

    /*
     * Either return on a complete match to the pattern, or use
     * the bad character and good suffix rules to perform a shift.
     */
    if (i == 0)
      return &T[k-M+1];
    else {
#ifdef STATS

      if (i == M)
        node->num_init_mismatch++;

      bshift = i - R[(int) T[h]];
      if (bshift < 1)
        bshift = 1;
      node->shift_cost++;

      gshift = 0;
      if (i < M) {
        i++;     /* Add one here because bad shift deals with mismatched */
                 /* character, whereas good shift deals with the last    */
                 /* matched character.                                   */
        if (Lprime[i] > 0)
          gshift = M - Lprime[i];
        else
          gshift = M - lprime[i];
        node->shift_cost++;
      }

      k += (bshift > gshift ? bshift : gshift);
      node->num_shifts++;

#else

      bshift = i - R[(int) T[h]];
      if (bshift < 1)
        bshift = 1;

      gshift = 0;
      if (i < M) {
        i++;     /* Add one here because bad shift deals with mismatched */
                 /* character, whereas good shift deals with the last    */
                 /* matched character.                                   */
        if (Lprime[i] > 0)
          gshift = M - Lprime[i];
        else
          gshift = M - lprime[i];
      }

      k += (bshift > gshift ? bshift : gshift);

#endif
    }
  }

  return NULL;
}  


/*
 *
 * Boyer-Moore searching using the extended bad character rule and 
 * the strong good suffix rule.
 *
 */
BM_STRUCT *bmextgood_prep(char *P, int M, int copyflag)
{
  int i, j, *R, *Rnext, *Lprime, *lprime, *Z;
  char *buf;
#ifdef STRMAT
  char *Pbar;
  Z_STRUCT *zstruct;
#else
  int k, count;
#endif  
  BM_STRUCT *node;

  P--;            /* Shift to make sequence be P[1],...,P[M] */

  /*
   * Allocate everything.
   */
  if ((node = malloc(sizeof(BM_STRUCT))) == NULL)
    return NULL;

  memset(node, 0, sizeof(BM_STRUCT));
  node->type = BM_EXTGOOD;
  node->M = M;
  node->copyflag = copyflag;

  if ((R = node->R = malloc(128 * sizeof(int))) == NULL) {
    bm_free(node);
    return NULL;
  }
  memset(R, 0, 128 * sizeof(int));

  if ((Rnext = node->Rnext = malloc((M + 1) * sizeof(int))) == NULL) {
    bm_free(node);
    return NULL;
  }
  memset(Rnext, 0, (M + 1) * sizeof(int));

  if ((Lprime = node->Lprime = malloc((M + 1) * sizeof(int))) == NULL) {
    bm_free(node);
    return NULL;
  }
  memset(Lprime, 0, (M + 1) * sizeof(int));

  if ((lprime = node->lprime = malloc((M + 1) * sizeof(int))) == NULL) {
    bm_free(node);
    return NULL;
  }
  memset(lprime, 0, (M + 1) * sizeof(int));

  if (!copyflag) 
    node->P = P;
  else {
    if ((buf = malloc(M + 2)) == NULL) {
      bm_free(node);
      return NULL;
    }

    buf[0] = buf[M+1] = '\0';
    memcpy(buf + 1, P + 1, M);
    node->P = buf;
  }

  /*
   * Extended bad character rule.  R holds the heads of the lists and
   * Rnext holds the rest of the values.
   */
  for (i=1; i <= M; i++) {
    Rnext[i] = R[(int) P[i]];
    R[(int) P[i]] = i;

#ifdef STATS
    node->prep_compares++;
#endif
  }

  /*
   * Good Suffix rule.
   *
   * First, compute the Z values for the reversed string, either using the
   * Z-values algorithm in "z.c", if available, or using a brute force
   * computation.
   */
  Z = NULL;
#ifdef STRMAT

  zstruct = NULL;
  if ((Pbar = malloc(M + 1)) != NULL) {
    for (i=0; i < M; i++)
      Pbar[i] = P[M-i];
    Pbar[M] = '\0';

    if ((zstruct = z_build(Pbar, M, 0)) != NULL) {
      Z = zstruct->Z;

#ifdef STATS
      node->prep_compares += zstruct->prep_compares;
#endif
    }

    free(Pbar);
  }

#else

  if ((Z = malloc((M + 1) * sizeof(int))) != NULL) {
    for (i=1; i < M; i++) {
      count = 0;
      for (j=i,k=M; j > 0 && P[j] == P[k]; j--,k--)
        count++;
      Z[i] = count;

#ifdef STATS
      node->prep_compares += count + 1;
#endif
    }
  }

#endif

  if (Z == NULL) {
    bm_free(node);
    return NULL;
  }

  /*
   * Next, compute the L' and l' values from the Z values.
   */
  for (j=1; j < M; j++) {
    i = M - Z[M-j+1] + 1;
    Lprime[i] = j;

#ifdef STATS
    node->prep_compares++;
#endif
  }

  lprime[M] = (P[1] == P[M] ? 1 : 0);
#ifdef STATS
  node->prep_compares++;
#endif

  for (i=M-1,j=2; i >= 2; i--,j++) {
    if (Z[M-j+1] == j)
      lprime[i] = j;
    else
      lprime[i] = lprime[i+1];

#ifdef STATS
    node->prep_compares++;
#endif
  }

#ifdef STRMAT
  z_free(zstruct);
#else
  free(Z);
#endif

  return node;
}


char *bmextgood_search(BM_STRUCT *node, char *T, int N, int initmatch)
{
  int i, k, h, M, pos, gshift, bshift, *R, *Rnext, *Lprime, *lprime;
  char *P;

  if (node->type != BM_EXTGOOD)
    return NULL;

  T--;            /* Shift to make sequence be T[1],...,T[N] */

  P = node->P;
  M = node->M;
  R = node->R;
  Rnext = node->Rnext;
  Lprime = node->Lprime;
  lprime = node->lprime;

  /*
   * Run the Boyer-Moore algorithm over the sequence.  Stop at the first
   * position to match the complete pattern.
   */
  k = M;
  if (initmatch) {
    /*
     * Perform a good shift to shift past the match.
     */
    gshift = M - lprime[2];
    k += gshift;

#ifdef STATS
    node->shift_cost++;
    node->num_shifts++;
#endif
  }

  while (k <= N) {
    /*
     * Compare the pattern to the string ending at k (in T).
     */
#ifdef STATS

    i = M;
    h = k;
    while (i > 0 && P[i] == T[h]) {
      i--;
      h--;
      node->num_compares++;
    }
    node->num_compares++;

#else

    for (i=M,h=k; i > 0 && P[i] == T[h]; i--,h--) 
      ;

#endif

    /*
     * Either return on a complete match to the pattern, or
     * use the extended bad character and good suffix rules
     * to perform a shift.
     */
    if (i == 0)
      return &T[k-M+1];
    else {
#ifdef STATS

      if (i == M)
        node->num_init_mismatch++;

      pos = R[(int) T[h]];
      node->shift_cost++;
      while (pos >= i) {
        pos = Rnext[pos];
        node->shift_cost++;
      }
      bshift = i - pos;

      gshift = 0;
      if (i < M) {
        i++;     /* Add one here because bad shift deals with mismatched */
                 /* character, whereas good shift deals with the last    */
                 /* matched character.                                   */
        if (Lprime[i] > 0)
          gshift = M - Lprime[i];
        else
          gshift = M - lprime[i];
        node->shift_cost++;
      }

      k += (bshift > gshift ? bshift : gshift);
      node->num_shifts++;

#else

      pos = R[(int) T[h]];
      while (pos >= i)
        pos = Rnext[pos];
      bshift = i - pos;

      gshift = 0;
      if (i < M) {
        i++;     /* Add one here because bad shift deals with mismatched */
                 /* character, whereas good shift deals with the last    */
                 /* matched character.                                   */
        if (Lprime[i] > 0)
          gshift = M - Lprime[i];
        else
          gshift = M - lprime[i];
      }

      k += (bshift > gshift ? bshift : gshift);

#endif
    }
  }

  return NULL;
}  



/*
 * bm_free
 *
 * Free up the allocated BM_STRUCT structure.
 *
 * Parameters:   node  -  a BM_STRUCT structure
 *
 * Returns:  nothing.
 */
void bm_free(BM_STRUCT *node)
{
  if (node == NULL)
    return;

  if (node->copyflag && node->P != NULL)
    free(node->P);
  if (node->R != NULL)
    free(node->R);
  if (node->Rnext != NULL)
    free(node->Rnext);
  if (node->Lprime != NULL)
    free(node->Lprime);
  if (node->lprime != NULL)
    free(node->lprime);
  
  free(node);
}
