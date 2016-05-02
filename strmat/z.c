/*
 * z.c
 *
 * Contains the prcedures to calculate and print the Z values of a string.
 * Where Z[i] is the length of the longest prefix of S[i..n] which also occurs
 * as a prefix of S.  See Sections 1.4 and 1.5 of Dan Gusfield's text for
 * more information.
 *
 * NOTES:
 *    8/94  -  Original Implementation  (Sean Davis)
 *    9/94  -  Redid Implementation  (James Knight)
 *    7/95  -  Rewrote the code with the new program organization and
 *             combined the z.c and zmatch.c contents.  (James Knight)
 *    3/96  -  Modularized the code  (James Knight)
 *    7/96  -  Finished the modularization  (James Knight)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "z.h"


/*
 * z_build
 *
 * Implementation of the Z values algorithm for a single sequence.
 *
 * Parameters:   S         -  the sequence
 *               M         -  the sequence length
 *               copyflag  -  make a copy of the sequence?
 *
 * Returns:  An initialized Z_STRUCT structure.
 */
Z_STRUCT *z_build(char *S, int M, int copyflag)
{
  int j, k, l, r, beta, kprime, *Z;
  char *buf;
  Z_STRUCT *node;

  S--;            /* Shift to make sequence be S[1],...,S[M] */

  /*
   * Allocate everything.
   */
  if ((node = malloc(sizeof(Z_STRUCT))) == NULL)
    return NULL;
  memset(node, 0, sizeof(Z_STRUCT));
  node->M = M;
  node->copyflag = copyflag;

  if ((Z = node->Z = malloc((M + 1) * sizeof(int))) == NULL) {
    free(node);
    return NULL;
  }

  if (!copyflag)
    node->S = S;
  else {
    if ((buf = malloc(M + 2)) == NULL) {
      free(node->Z);
      free(node);
      return NULL;
    }

    buf[0] = buf[M+1] = '\0';
    memcpy(buf + 1, S + 1, M);
    node->S = buf;
  }

  /*
   * Run the Z Values Algorithm.
   */
  l = r = 0;
  Z[0] = Z[1] = 0;
  for (k=2; k <= M; k++) {
    /*
     * Case 1.
     */
    if (k > r) {
#ifdef STATS
      j = 0;
      while (k + j <= M && S[1+j] == S[k+j]) {
        j++;
        node->prep_compares++;
      }
      node->prep_compares++;
#else
      for (j=0; k + j <= M && S[1+j] == S[k+j]; j++) ;
#endif

      Z[k] = j;
      r = k + Z[k] - 1;
      l = k;
    }

    else {
      beta = r - k + 1;
      kprime = k - l + 1;

      /*
       * Case 2a and case 2b.
       */
      if (Z[kprime] < beta)
        Z[k] = Z[kprime];
      else {
#ifdef STATS
        j = 1;
        while (r + j <= M && S[beta+j] == S[r+j]) {
          j++;
          node->prep_compares++;
        }
        node->prep_compares++;
#else
        for (j=1; r + j <= M && S[beta+j] == S[r+j]; j++) ;
#endif

        Z[k] = r + j - k;
        r = r + j - 1;
        l = k;
      }
    }
  }

  return node;
}


/*
 * z_search
 *
 * Search a sequence text using the Z values algorithm.
 *
 * Parameters:   node      -  a preprocessed pattern
 *               T         -  the sequence
 *               N         -  the sequence length
 *               initmatch -  does the text begin with a match?
 *
 * Returns:  The location of the first match to the pattern, or NULL.
 */
char *z_search(Z_STRUCT *node, char *T, int N, int initmatch)
{
  int j, l, r, k, beta, kprime, M, *Z, ZT;
  char *P;

  T--;            /* Shift to make sequence be T[1],...,T[M] */

  P = node->S;
  M = node->M;
  Z = node->Z;

  /*
   * Run the Z values algorithm over the sequence.  Stop at the first
   * position to match the complete pattern (i.e., Z[k] == M).
   */
  l = 0;
  r = (!initmatch ? 0 : M - 1);
  k = (!initmatch ? 1 : 2);

  for ( ; k + M <= N; k++) {
    /*
     * Case 1.
     */
    if (k > r) {
#ifdef STATS
      j = 0;
      while (1 + j <= M && P[1+j] == T[k+j]) {
        j++;
        node->num_compares++;
      }
      node->num_compares++;
#else
      for (j=0; 1 + j <= M && P[1+j] == T[k+j]; j++) ;
#endif

      ZT = j;
      r = k + j - 1;
      l = k;
    }

    else {
      beta = r - k + 1;
      kprime = k - l + 1;

      /*
       * Case 2a and case 2b.
       */
      if (Z[kprime] < beta)
        ZT = Z[kprime];
      else {
#ifdef STATS
        j = 1;
        while (beta + j <= M && r + j <= N && P[beta+j] == T[r+j]) {
          j++;
          node->num_compares++;
        }
        node->num_compares++;
#else
        for (j=1; beta + j <= M && r + j <= N && P[beta+j] == T[r+j]; j++) ;
#endif

        ZT = r + j - k;
        r = r + j - 1;
        l = k;
      }
    }

    if (ZT == M)
      return &T[k];
  }

  return NULL;
}  


/*
 * z_free
 *
 * Free up the allocated Z_STRUCT structure.
 *
 * Parameters:   node  -  a Z_STRUCT structure
 *
 * Returns:  nothing.
 */
void z_free(Z_STRUCT *node)
{
  if (node == NULL)
    return;

  if (node->Z != NULL)
    free(node->Z);
  if (node->copyflag && node->S != NULL)
    free(node->S);

  free(node);
}

