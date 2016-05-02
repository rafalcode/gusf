
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bmopt.h"


BMOPT_STRUCT *bmopt_prep(char *P, int M)
{
  int i, j, size, count, *buf, *R, *B, *Z, *L, *l;
  char *p, *t;
  BMOPT_STRUCT *node;

  if ((node = malloc(sizeof(BMOPT_STRUCT))) == NULL)
    return NULL;

  if ((node->P = malloc(M+1)) == NULL) {
    free(node);
    return NULL;
  }
  memcpy(node->P, P, M);
  node->P[M] = '\0';
  node->M = M;

  size = (770 + 3 * (M + 1)) * sizeof(int);
  if ((buf = malloc(size)) == NULL) {
    free(node->P);
    free(node);
    return NULL;
  }
  memset(buf, 0, size);

  node->R = R = buf + 129;  buf += 385;
  node->B = B = buf + 129;  buf += 385;
  Z = buf;  buf += M + 1;
  node->L = L = buf;  buf += M + 1;
  l = buf;

  /*
   * Preprocess the keyword.
   *
   * Simple bad character rule.
   */
  for (i=1,p=P; i <= M; i++,p++)
    R[(int) *p] = i;

  for (i=-128; i < 256; i++)
    B[i] = M - R[i];

  /*
   * Good Suffix rule.
   *
   * First compute the Z-values for the reversed string (i.e., Z[j] is
   * the length of the longest suffix of P[1..j] which is also a suffix
   * of the whole string P[1..M]).  Then compute the L' and l' values.
   */
  for (j=1; j < M; j++) {
    for (p=P+j-1,t=P+M-1,count=0; p >= P && *p == *t; p--,t--,count++) ;
    Z[j] = count;
  }

  for (j=1; j < M; j++) {
    i = M - Z[j] + 1;
    L[i] = j;
  }

  l[M] = (P[0] == P[M-1] ? 1 : 0);
  for (i=M-1,j=2; i >= 2; i--,j++) {
    l[i] = (Z[j] == j ? j : l[i+1]);
  }

  /*
   * Combine the L' and l' values, and presubtract them from M to get
   * the correct goodshift values.
   */
  for (i=1; i <= M; i++) {
    if (L[i] > 0)
      L[i] = M - L[i];
    else
      L[i] = M - l[i];
    if (L[i] == 0)
      L[i] = 1;
  }

  return node;
}

char *bmopt_search(BMOPT_STRUCT *node, char *T, char *Tend)
{
  int i, shift, gshift, M, *R, *B, *L;
  char *p, *pend, *t, *t2, *P;

  P = node->P;
  M = node->M;
  R = node->R;
  B = node->B;
  L = node->L;

  pend = P + M;

  t = T + M - 1;
  while (t < Tend) {
    /*
     * Run a small loop handling the cases of initial mismatches.
     */
    t += (shift = B[(int) *t]);
    while (shift) {
      t += B[(int) *t];
      t += B[(int) *t];
      t += (shift = B[(int) *t]);
      if (t >= Tend)
        goto BMOPT_LOOP_END;
    }

    /*
     * Check the possible match.
     */
    p = pend - 2;
    i = M - 1;
    t2 = t - 1;
    while (i && *p == *t2) {
      i--;
      p--;
      t2--;
    }

    if (i) {
      gshift = L[i+1];
      if ((i -= R[(int) *t2]) > gshift)
        t += i;
      else
        t += gshift;
    }
    else
      return t-M+1;
  }
BMOPT_LOOP_END:

  return NULL;
}

void bmopt_free(BMOPT_STRUCT *node)
{
  if (node == NULL)
    return;

  if (node->P != NULL)
    free(node->P);
  if (node->R != NULL)
    free(node->R);
  free(node);
}
