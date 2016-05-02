
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sary_match.h"


/*
 * sary_match_*_prep
 *
 * Preprocessing for the suffix array exact matching algorithm using
 * a naive binary search.
 *
 * Parameters:   T         -  the text
 *               M         -  the text length
 *               copyflag  -  make a copy of the sequence?
 *
 * Returns:  An initialized SARYMAT_STRUCT structure.
 */
static SARYMAT_STRUCT *int_sary_match_prep(char *T, int M, int copyflag,
                                           SARY_MATCH_TYPE type);

SARYMAT_STRUCT *sary_match_naive_prep(char *T, int M, int copyflag)
{  return int_sary_match_prep(T, M, copyflag, NAIVE_MATCH);  }
SARYMAT_STRUCT *sary_match_mlr_prep(char *T, int M, int copyflag)
{  return int_sary_match_prep(T, M, copyflag, MLR_MATCH);  }
SARYMAT_STRUCT *sary_match_lcp_prep(char *T, int M, int copyflag)
{  return int_sary_match_prep(T, M, copyflag, LCP_MATCH);  }

static SARYMAT_STRUCT *int_sary_match_prep(char *T, int M, int copyflag,
                                           SARY_MATCH_TYPE type)
{
  char *buf;
  SARYMAT_STRUCT *node;

  T--;            /* Shift to make sequence be T[1],...,T[M] */

  /*
   * Allocate everything.
   */
  if ((node = malloc(sizeof(SARYMAT_STRUCT))) == NULL)
    return NULL;
  memset(node, 0, sizeof(SARYMAT_STRUCT));

  node->type = type;
  node->copyflag = copyflag;
  node->M = M;

  if (!copyflag)
    node->T = T;
  else {
    if ((buf = malloc(M + 2)) == NULL) {
      free(node);
      return NULL;
    }

    buf[0] = buf[M+1] = '\0';
    memcpy(buf + 1, T + 1, M);
    node->T = buf;
  }

  /*
   * Create the suffix array.
   */
  switch (type) {
  case NAIVE_MATCH:
  case MLR_MATCH:
    if ((node->sary = sary_qsort_build(T+1, M, 0)) == NULL) {
      sary_match_free(node);
      return NULL;
    }
    break;

  case LCP_MATCH:
    if ((node->sary = sary_stree_build(T+1, M, 0)) == NULL) {
      sary_match_free(node);
      return NULL;
    }
    break;
  }

  return node;
}


/*
 * sary_match_naive_first
 *
 * Use the suffix array to search a text for matches to a pattern.
 *
 * Note:  This procedure and sary_match_naive_next do NOT return the
 *        matches in the order they appear in the text, because of the
 *        the way the matches are found (i.e. using the suffix array).
 *
 * Parameters:   node      -  the preprocessed text
 *               P         -  the pattern
 *               N         -  the pattern's length
 *
 * Returns:  The location of a match to the pattern, or NULL.
 */
char *sary_match_naive_first(SARYMAT_STRUCT *smstruct, char *P, int N)
{
  int k, L, R, midpoint, matchpos, pos, M, *Pos;
  char *T;

  P--;            /* Shift to make sequence be P[1],...,P[N] */
  
  T = smstruct->T;
  M = smstruct->M;
  Pos = smstruct->sary->Pos;

  /*
   * Perform a binary search to find the smallest position of a match
   * in the suffix array.
   *
   * In this binary search, the region between L and R (inclusively)
   * is the region yet to be determined (i.e., everything from L-1 down
   * is strictly less than the pattern and everything from R+1 up is 
   * greater than or equal to the pattern).
   */
  L = 1;
  R = M;
  while (L <= R) {
    midpoint = (L + R) / 2;
    pos = Pos[midpoint];

#ifdef STATS
    smstruct->search_depth++;

    for (k=0; 1 + k <= N && pos + k <= M && P[1+k] == T[pos+k]; k++)
      smstruct->num_compares++;
    smstruct->num_compares++;
#else
    for (k=0; 1 + k <= N && pos + k <= M && P[1+k] == T[pos+k]; k++) ;
#endif

    if (1 + k > N && pos + k > M)
      R = midpoint - 1;
    else if (1 + k > N && pos + k <= M)
      R = midpoint - 1;
    else if (1 + k <= N && pos + k > M)
      L = midpoint + 1;
    else if (1 + k <= N && pos + k <= M && P[1+k] > T[pos+k])
      L = midpoint + 1;
    else if (1 + k <= N && pos + k <= M && P[1+k] < T[pos+k])
      R = midpoint - 1;

  }

  /*
   * If a match exists, the first match must be at R+1.
   */
  smstruct->i = R+1;

  if (R == M) {
    smstruct->iprime = R;
    return NULL;
  }

  pos = Pos[R+1];

#ifdef STATS
  for (k=0; 1 + k <= N && pos + k <= M && P[1+k] == T[pos+k]; k++)
    smstruct->num_compares++;
  smstruct->num_compares++;
#else
  for (k=0; 1 + k <= N && pos + k <= M && P[1+k] == T[pos+k]; k++) ;
#endif

  if (k < N) {
    smstruct->iprime = R;
    return NULL;
  }

  /*
   * Now that we know that matches exist, perform a second binary search
   * to find the largest position containing a match.
   *
   * In this binary search, the region between L and R (inclusively)
   * is the region yet to be determined (i.e., everything from L-1 down
   * is less than or equal to the pattern and everything from R+1 up is 
   * strictly greater than the pattern).
   *
   * NOTE:  Remember, a text suffix is "equal" to the pattern here whenever
   *        the complete pattern matches the beginning of the suffix
   *        (regardless of what comes after).  That's why this search
   *        is not symmetrical with the search above.
   */
  L = 1;
  R = M;
  while (L <= R) {
    midpoint = (L + R) / 2;
    pos = Pos[midpoint];

#ifdef STATS
    for (k=0; 1 + k <= N && pos + k <= M && P[1+k] == T[pos+k]; k++)
      smstruct->num_compares++;
    smstruct->num_compares++;
#else
    for (k=0; 1 + k <= N && pos + k <= M && P[1+k] == T[pos+k]; k++) ;
#endif

    if (1 + k > N && pos + k > M)
      L = midpoint + 1;
    else if (1 + k > N && pos + k <= M)
      L = midpoint + 1;
    else if (1 + k <= N && pos + k > M)
      L = midpoint + 1;
    else if (1 + k <= N && pos + k <= M && P[1+k] > T[pos+k])
      L = midpoint + 1;
    else if (1 + k <= N && pos + k <= M && P[1+k] < T[pos+k])
      R = midpoint - 1;
  }

  smstruct->iprime = L - 1;

  matchpos = Pos[smstruct->i++];
  return &T[matchpos];
}


/*
 * sary_match_naive_next
 *
 * Return the next match of a pattern to a text, if another match exists.
 *
 * Parameters:  smstruct  -  A preprocessed and initialized search
 * 
 * Returns:  the location of a match in the text, or NULL.
 */
char *sary_match_naive_next(SARYMAT_STRUCT *smstruct)
{
  int matchpos, *Pos;
  char *T;

  T = smstruct->T;
  Pos = smstruct->sary->Pos;

  if (smstruct->i > smstruct->iprime)
    return NULL;
  else {
    matchpos = Pos[smstruct->i++];
    return &T[matchpos];
  }
}


/*
 * sary_match_mlr_first
 *
 * Use the suffix array to search a text for matches to a pattern.
 *
 * Note:  This procedure and sary_match_mlr_next do NOT return the
 *        matches in the order they appear in the text, because of the
 *        the way the matches are found (i.e. using the suffix array).
 *
 * Parameters:   node      -  the preprocessed text
 *               P         -  the pattern
 *               N         -  the pattern's length
 *
 * Returns:  The location of a match to the pattern, or NULL.
 */
char *sary_match_mlr_first(SARYMAT_STRUCT *smstruct, char *P, int N)
{
  int k, L, R, l, r, mlr, l_init, r_init, midpoint, matchpos, pos, M, *Pos;
  char *T;

  P--;            /* Shift to make sequence be P[1],...,P[N] */
  
  T = smstruct->T;
  M = smstruct->M;
  Pos = smstruct->sary->Pos;

  /*
   * Compute the initial value of l (by finding the longest common
   * prefix between the pattern and Pos[1]).
   *
   * And check to see if the pattern is smaller than all of the suffixes.
   */
#ifdef STATS
  for (k=0; 1 + k <= N && Pos[1] + k <= M && P[1+k] == T[Pos[1]+k]; k++)
    smstruct->num_compares++;
  smstruct->num_compares++;
#else
  for (k=0; 1 + k <= N && Pos[1] + k <= M && P[1+k] == T[Pos[1]+k]; k++) ;
#endif

  if (1 + k <= N && Pos[1] + k <= M && P[1+k] < T[Pos[1]+k]) {
    smstruct->i = 1;
    smstruct->iprime = 0;
    return NULL;
  }

  l_init = k;

  /*
   * Compute the initial value of r (by finding the longest common
   * prefix between the pattern and Pos[M]).
   *
   * And check to see if the pattern is larger than all of the suffixes.
   */
#ifdef STATS
  for (k=0; 1 + k <= N && Pos[M] + k <= M && P[1+k] == T[Pos[M]+k]; k++)
    smstruct->num_compares++;
  smstruct->num_compares++;
#else
  for (k=0; 1 + k <= N && Pos[M] + k <= M && P[1+k] == T[Pos[M]+k]; k++) ;
#endif

  if (1 + k <= N && (Pos[M] + k > M || P[1+k] > T[Pos[M]+k])) {
    smstruct->i = 1;
    smstruct->iprime = 0;
    return NULL;
  }

  r_init = k;

  /*
   * Now, check to see if the pattern equals the first suffix.
   * If so, set L and R so that we skip the binary search (because
   * we've found the starting point).
   *
   * Otherwise, initialize L, R, l and r for the binary search.
   */
  if (1 + l_init > N) {
    L = 0;
    R = 1;
    l = r = 0;
  }
  else {
    L = 1;
    R = M;
    l = l_init;
    r = r_init;
  }

  /*
   * Perform a binary search to find the smallest position of a match
   * in the suffix array.
   *
   * In this binary search, the region between L and R (exclusively)
   * is the region yet to be determined (i.e., everything from L down
   * is strictly less than the pattern and everything from R up is 
   * greater than or equal to the pattern).
   */
  while (L + 1 < R) {
    midpoint = (L + R) / 2;
    pos = Pos[midpoint];

    mlr = (l < r ? l : r);

#ifdef STATS
    smstruct->search_depth++;

    for (k=mlr; 1 + k <= N && pos + k <= M && P[1+k] == T[pos+k]; k++)
      smstruct->num_compares++;
    smstruct->num_compares++;
#else
    for (k=mlr; 1 + k <= N && pos + k <= M && P[1+k] == T[pos+k]; k++) ;
#endif

    if (1 + k > N && pos + k > M)
      R = midpoint;
    else if (1 + k > N && pos + k <= M)
      R = midpoint;
    else if (1 + k <= N && pos + k > M)
      L = midpoint;
    else if (1 + k <= N && pos + k <= M && P[1+k] > T[pos+k])
      L = midpoint;
    else if (1 + k <= N && pos + k <= M && P[1+k] < T[pos+k])
      R = midpoint;

    if (L == midpoint)
      l = k;
    else
      r = k;
  }

  /*
   * If a match exists, the first match must be at R.
   */
  smstruct->i = R;
  pos = Pos[R];

#ifdef STATS
  for (k=0; 1 + k <= N && pos + k <= M && P[1+k] == T[pos+k]; k++)
    smstruct->num_compares++;
  smstruct->num_compares++;
#else
  for (k=0; 1 + k <= N && pos + k <= M && P[1+k] == T[pos+k]; k++) ;
#endif

  if (k < N) {
    smstruct->iprime = R - 1;
    return NULL;
  }

  /*
   * Now that we know that matches exist, perform a second binary search
   * to find the largest position containing a match.
   *
   * In this binary search, the region between L and R (exclusively)
   * is the region yet to be determined (i.e., everything from L down
   * is less than or equal to the pattern and everything from R up is 
   * strictly greater than the pattern).
   *
   * NOTE:  Remember, a text suffix is "equal" to the pattern here whenever
   *        the complete pattern matches the beginning of the suffix
   *        (regardless of what comes after).  That's why this search
   *        is not symmetrical with the search above.
   */
  if (1 + r_init > N) {
    L = M;
    R = M + 1;
    l = r = 0;
  }
  else {
    L = 1;
    l = l_init;
    R = M;
    r = r_init;
  }

  while (L + 1 < R) {
    midpoint = (L + R) / 2;
    pos = Pos[midpoint];

    mlr = (l < r ? l : r);

#ifdef STATS
    for (k=mlr; 1 + k <= N && pos + k <= M && P[1+k] == T[pos+k]; k++)
      smstruct->num_compares++;
    smstruct->num_compares++;
#else
    for (k=mlr; 1 + k <= N && pos + k <= M && P[1+k] == T[pos+k]; k++) ;
#endif

    if (1 + k > N && pos + k > M)
      L = midpoint;
    else if (1 + k > N && pos + k <= M)
      L = midpoint;
    else if (1 + k <= N && pos + k > M)
      L = midpoint;
    else if (1 + k <= N && pos + k <= M && P[1+k] > T[pos+k])
      L = midpoint;
    else if (1 + k <= N && pos + k <= M && P[1+k] < T[pos+k])
      R = midpoint;

    if (L == midpoint)
      l = k;
    else
      r = k;
  }

  smstruct->iprime = L;

  matchpos = Pos[smstruct->i++];
  return &T[matchpos];
}


/*
 * sary_match_naive_next
 *
 * Return the next match of a pattern to a text, if another match exists.
 *
 * Parameters:  smstruct  -  A preprocessed and initialized search
 * 
 * Returns:  the location of a match in the text, or NULL.
 */
char *sary_match_mlr_next(SARYMAT_STRUCT *smstruct)
{
  int matchpos, *Pos;
  char *T;

  T = smstruct->T;
  Pos = smstruct->sary->Pos;

  if (smstruct->i > smstruct->iprime)
    return NULL;
  else {
    matchpos = Pos[smstruct->i++];
    return &T[matchpos];
  }
}


/*
 * sary_match_lcp_first
 *
 * Use the suffix array to search a text for matches to a pattern.
 *
 * Note:  This procedure and sary_match_naive_next do NOT return the
 *        matches in the order they appear in the text, because of the
 *        the way the matches are found (i.e. using the suffix array).
 *
 * Parameters:   node      -  the preprocessed text
 *               P         -  the pattern
 *               N         -  the pattern's length
 *
 * Returns:  The location of a match to the pattern, or NULL.
 */
char *sary_match_lcp_first(SARYMAT_STRUCT *smstruct, char *P, int N)
{
  int k, L, R, l, r, mlr, l_init, r_init, midpoint;
  int M, L_to_M, M_to_R, *Lcp, pos, *Pos, matchpos;
  char *T;

  P--;            /* Shift to make sequence be P[1],...,P[N] */
  
  T = smstruct->T;
  M = smstruct->M;
  Pos = smstruct->sary->Pos;
  Lcp = smstruct->sary->lcp;

  /*
   * Store the pattern length for use in sary_match_lcp_next.
   */
  smstruct->N = N;

  /*
   * Compute the initial value of l (by finding the longest common
   * prefix between the pattern and Pos[1]).
   *
   * And check to see if the pattern is smaller than all of the suffixes.
   */
#ifdef STATS
  for (k=0; 1 + k <= N && Pos[1] + k <= M && P[1+k] == T[Pos[1]+k]; k++)
    smstruct->num_compares++;
  smstruct->num_compares++;
#else
  for (k=0; 1 + k <= N && Pos[1] + k <= M && P[1+k] == T[Pos[1]+k]; k++) ;
#endif

  if (1 + k > N) {
    smstruct->i = 2;
    smstruct->iprime = M;
    return &T[Pos[1]];
  }
  else if (Pos[1] + k <= M && P[1+k] < T[Pos[1]+k]) {
    smstruct->i = 1;
    smstruct->iprime = 0;
    return NULL;
  }

  l_init = k;

  /*
   * Compute the initial value of r (by finding the longest common
   * prefix between the pattern and Pos[M]).
   *
   * And check to see if the pattern is larger than all of the suffixes.
   */
#ifdef STATS
  for (k=0; 1 + k <= N && Pos[M] + k <= M && P[1+k] == T[Pos[M]+k]; k++)
    smstruct->num_compares++;
  smstruct->num_compares++;
#else
  for (k=0; 1 + k <= N && Pos[M] + k <= M && P[1+k] == T[Pos[M]+k]; k++) ;
#endif

  if (1 + k <= N && (Pos[M] + k > M || P[1+k] > T[Pos[M]+k])) {
    smstruct->i = 1;
    smstruct->iprime = 0;
    return NULL;
  }

  r_init = k;

  /*
   * Perform a binary search to find the smallest position of a match
   * in the suffix array.
   *
   * In this binary search, the region between L and R (exclusively)
   * is the region yet to be determined (i.e., everything from L down
   * is strictly less than the pattern and everything from R up is 
   * greater than or equal to the pattern).
   */
  L = 1;
  l = l_init;
  R = M;
  r = r_init;

  L_to_M = 2;
  M_to_R = 3;

  while (L + 1 < R) {
    midpoint = (L + R) / 2;
    pos = Pos[midpoint];

#ifdef STATS
    smstruct->search_depth++;
#endif

    /*
     * The super-accelerant cases.
     */
    if (l > r && Lcp[L_to_M] > l)
      L = midpoint;
    else if (l > r && Lcp[L_to_M] < l) {
      R = midpoint;
      r = Lcp[L_to_M];
    }
    else if (l < r && Lcp[M_to_R] > r)
      R = midpoint;
    else if (l < r && Lcp[M_to_R] < r) {
      L = midpoint;
      l = Lcp[M_to_R];
    }
    else {
      /*
       * All of the other cases (requiring explicit matching):
       *    1)  l == r
       *    2)  l > r && Lcp[L_to_M] == l
       *    3)  l < r && Lcp[M_to_R] == r
       */
      mlr = (l >= r ? l : r);

#ifdef STATS
      for (k=mlr; 1 + k <= N && pos + k <= M && P[1+k] == T[pos+k]; k++)
        smstruct->num_compares++;
      smstruct->num_compares++;
#else
      for (k=mlr; 1 + k <= N && pos + k <= M && P[1+k] == T[pos+k]; k++) ;
#endif

      if (1 + k > N && pos + k > M)
        R = midpoint;
      else if (1 + k > N && pos + k <= M)
        R = midpoint;
      else if (1 + k <= N && pos + k > M)
        L = midpoint;
      else if (1 + k <= N && pos + k <= M && P[1+k] > T[pos+k])
        L = midpoint;
      else if (1 + k <= N && pos + k <= M && P[1+k] < T[pos+k])
        R = midpoint;

      if (L == midpoint)
        l = k;
      else
        r = k;
    }

    /*
     * Find the Lcp values for the next step in the binary search.
     */
    if (L == midpoint) {
      L_to_M = 2 * M_to_R;
      M_to_R = 2 * M_to_R + 1;
    }
    else {
      M_to_R = 2 * L_to_M + 1;
      L_to_M = 2 * L_to_M;
    }
  }

  /*
   * If a match exists, the first match must be at R.
   *
   * It is not necessary here to compute the largest index of a suffix
   * that matches the pattern, because the Lcp(a,a+1) values provide
   * a way to determine whether each next suffix (after the first match)
   * matches the pattern in constant time.  See sary_match_lcp_next below.
   */
  smstruct->i = R;
  pos = Pos[R];

#ifdef STATS
  for (k=0; 1 + k <= N && pos + k <= M && P[1+k] == T[pos+k]; k++)
    smstruct->num_compares++;
  smstruct->num_compares++;
#else
  for (k=0; 1 + k <= N && pos + k <= M && P[1+k] == T[pos+k]; k++) ;
#endif

  if (k < N) {
    smstruct->iprime = R - 1;
    return NULL;
  }

  smstruct->iprime = M;
  matchpos = Pos[smstruct->i++];
  return &T[matchpos];
}


/*
 * sary_match_lcp_next
 *
 * Return the next match of a pattern to a text, if another match exists.
 *
 * Because the Lcp(a,a+1) are already computed, those can be used to
 * determine in constant time whether the next suffix of the suffix
 * array matches the pattern (i.e., if it's common prefix with the
 * previous suffix is as long as the pattern).  This saves us from
 * performing the second binary search above.
 *
 * Parameters:  smstruct  -  A preprocessed and initialized search
 * 
 * Returns:  the location of a match in the text, or NULL.
 */
char *sary_match_lcp_next(SARYMAT_STRUCT *smstruct)
{
  int matchpos, M, *Pos, *leaves;
  char *T;

  T = smstruct->T;
  M = smstruct->M;
  Pos = smstruct->sary->Pos;
  leaves = smstruct->sary->lcp_leaves;

  if (smstruct->i > M)
    return NULL;
  else if (leaves[smstruct->i] < smstruct->N) {
    smstruct->i = M + 1;
    return NULL;
  }
  else {
    matchpos = Pos[smstruct->i++];
    return &T[matchpos];
  }
}


/*
 * sary_match_free
 *
 * Free up an SARYMAT_STRUCT structure.
 *
 * Parameters:  smstruct  -  an SARYMAT_STRUCT structure
 *
 * Returns:  nothing.
 */
void sary_match_free(SARYMAT_STRUCT *smstruct)
{
  if (smstruct == NULL)
    return;

  if (smstruct->sary != NULL)
    sary_free(smstruct->sary);
  if (smstruct->copyflag && smstruct->T)
    free(smstruct->T);

  free(smstruct);
}

