/*
 * naive.c
 *
 * Performs the simple, naive matching of a pattern to a text.
 *
 * NOTES:
 *    7/95  -  Original Implementation  (James Knight)
 *    3/96  -  Modularized the code  (James Knight)
 *    7/96  -  Finished the modularization  (James Knight)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "naive.h"


/*
 * naive_prep
 *
 * Preprocessing for the naive search algorithm.
 *
 * Parameters:   S         -  the sequence
 *               M         -  the sequence length
 *               copyflag  -  make a copy of the sequence?
 *
 * Returns:  An initialized NAIVE_STRUCT structure.
 */
NAIVE_STRUCT *naive_prep(char *S, int M, int copyflag)
{
  char *buf;
  NAIVE_STRUCT *node;

  S--;            /* Shift to make sequence be S[1],...,S[M] */

  /*
   * Allocate everything.
   */
  if ((node = malloc(sizeof(NAIVE_STRUCT))) == NULL)
    return NULL;
  memset(node, 0, sizeof(NAIVE_STRUCT));
  
  node->M = M;
  node->copyflag = copyflag;

  if (!copyflag)
    node->S = S;
  else {
    if ((buf = malloc(M + 2)) == NULL) {
      free(node);
      return NULL;
    }

    buf[0] = buf[M+1] = '\0';
    memcpy(buf + 1, S + 1, M);
    node->S = buf;
  }

  return node;
}



/*
 * naive_search
 *
 * Search a sequence text using the naive search algorithm.
 *
 * Parameters:   node      -  a preprocessed pattern
 *               T         -  the sequence
 *               N         -  the sequence length
 *               initmatch -  does the text begin with a match?
 *
 * Returns:  The location of the first match to the pattern, or NULL.
 */
char *naive_search(NAIVE_STRUCT *node, char *T, int N, int initmatch)
{
  int j, k, M;
  char *P;

  T--;            /* Shift to make sequence be T[1],...,T[M] */
  P = node->S;
  M = node->M;

  /*
   * Perform the matching.
   */
  k = (!initmatch ? 1 : 2);
  for ( ; k + M <= N; k++) {

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

    if (j == M)
      return &T[k];
  }

  return NULL;
}


/*
 * naive_free
 *
 * Free up the allocated NAIVE_STRUCT structure.
 *
 * Parameters:   node  -  a NAIVE_STRUCT structure
 *
 * Returns:  nothing.
 */
void naive_free(NAIVE_STRUCT *node)
{
  if (node == NULL)
    return;

  if (node->copyflag && node->S != NULL)
    free(node->S);

  free(node);
}


