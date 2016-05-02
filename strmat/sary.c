
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sary.h"
#include "sary_zerkle.h"
#ifdef STRMAT
#include "stree_strmat.h"
#include "stree_ukkonen.h"
#else
#include "stree.h"
#endif


/*
 * sary_qsort_build
 *
 * Build a suffix array using a quick sort of the indices into the string.
 *
 * Parameters:  S         -  the input string
 *              M         -  the string's length
 *              copyflag  -  whether to copy the input string
 *
 * Returns:  an initialized SARY_STRUCT structure
 */
static SARY_STRUCT *sarystruct;
static int sarycmp(int *a, int *b);

SARY_STRUCT *sary_qsort_build(char *S, int M, int copyflag)
{
  int i, *Pos;
  char *buf;
  SARY_STRUCT *sary;

  if (S == NULL)
    return NULL;

  S--;            /* Shift to make sequence be S[1],...,S[M] */

  /*
   * Allocate everything.
   */
  if ((sary = malloc(sizeof(SARY_STRUCT))) == NULL)
    return NULL;
  memset(sary, 0, sizeof(SARY_STRUCT));

  sary->M = M;
  sary->copyflag = copyflag;

  if (!copyflag)
    sary->S = S;
  else {
    if ((buf = malloc(M + 2)) == NULL) {
      free(sary);
      return NULL;
    }

    buf[0] = buf[M+1] = '\0';
    memcpy(buf + 1, S + 1, M);
    sary->S = buf;
  }
    
  if ((Pos = sary->Pos = malloc((M + 1) * sizeof(int))) == NULL) {
    sary_free(sary);
    return NULL;
  }

  /*
   * Compute the suffix array.
   */
  for (i=1; i <= M; i++)
    Pos[i] = i;

  sarystruct = sary;
  qsort(Pos+1, M, sizeof(int), (int (*)()) sarycmp);

  return sary;
}


/*
 * sarycmp
 *
 * Compare two suffixes to determine which is lexically smaller.
 *
 * Parameters:  a  -  Address where the first index into the suffix is stored
 *              b  -  Address where the other index is stored
 *
 * Returns:  -1,0,1 if the first suffix is less than, equal to or greater
 *           than the second suffix.
 */
static int sarycmp(int *a, int *b)
{
  int i, j, k, M;
  char *S;

  i = *a;
  j = *b;
  S = sarystruct->S;
  M = sarystruct->M;

  for (k=0; i + k <= M && j + k <= M && S[i+k] == S[j+k]; k++) ;

#ifdef STATS
  sarystruct->num_compares += k + 1;
#endif

  if (i + k > M)
    return (j + k > M ? 0 : -1);
  else if (j + k <= M && S[i+k] > S[j+k])
    return 1;
  else
    return -1;
}


/*
 * sary_zerkle_build
 *
 * Build a suffix array using the Zerkle code (implementing Dan's O(N log N)
 * algorithm).
 *
 * Parameters:  S         -  the input string
 *              M         -  the string's length
 *              copyflag  -  whether to copy the input string
 *
 * Returns:  an initialized SARY_STRUCT structure
 */
SARY_STRUCT *sary_zerkle_build(char *S, int M, int copyflag)
{
  int *Pos;
  char *buf;
  SARY_STRUCT *sary;

  if (S == NULL)
    return NULL;

  S--;            /* Shift to make sequence be S[1],...,S[M] */

  /*
   * Allocate everything.
   */
  if ((sary = malloc(sizeof(SARY_STRUCT))) == NULL)
    return NULL;
  memset(sary, 0, sizeof(SARY_STRUCT));

  sary->M = M;
  sary->copyflag = copyflag;

  if (!copyflag)
    sary->S = S;
  else {
    if ((buf = malloc(M + 2)) == NULL) {
      free(sary);
      return NULL;
    }

    buf[0] = buf[M+1] = '\0';
    memcpy(buf + 1, S + 1, M);
    sary->S = buf;
  }
    
  if ((Pos = sary->Pos = malloc((M + 1) * sizeof(int))) == NULL) {
    sary_free(sary);
    return NULL;
  }

  /*
   * Compute the suffix array.
   */
  zerkle(S + 1, M, Pos);

  return sary;
}


/*
 * sary_stree_build
 *
 * Build a suffix array by first building the suffix tree.  Also,
 * compute the lcp values for doing matching to the suffix array.
 *
 * Parameters:  S         -  the input string
 *              M         -  the string's length
 *              copyflag  -  whether to copy the input string
 *              tree      -  the suffix tree for S
 *
 * Returns:  an initialized SARY_STRUCT structure
 */

static void compute_arrays(SARY_STRUCT *sary, SUFFIX_TREE tree,
                           STREE_NODE node, int current_depth);
static int compute_lcp_values(SARY_STRUCT *sary, int min, int max, int index);

SARY_STRUCT *sary_stree_build(char *S, int M, int copyflag)
{
  int i, flag, lcp_size, midpoint, *Pos, *lcp, *leaves;
  char *buf;
  SARY_STRUCT *sary;
  SUFFIX_TREE tree;

  if (S == NULL)
    return NULL;

  S--;            /* Shift to make sequence be S[1],...,S[M] */

  /*
   * Allocate everything.
   */
  if ((sary = malloc(sizeof(SARY_STRUCT))) == NULL)
    return NULL;
  memset(sary, 0, sizeof(SARY_STRUCT));

  sary->M = M;
  sary->copyflag = copyflag;

  if (!copyflag)
    sary->S = S;
  else {
    if ((buf = malloc(M + 2)) == NULL) {
      free(sary);
      return NULL;
    }

    buf[0] = buf[M+1] = '\0';
    memcpy(buf + 1, S + 1, M);
    sary->S = buf;
  }
    
  if ((Pos = sary->Pos = malloc((M + 1) * sizeof(int))) == NULL) {
    sary_free(sary);
    return NULL;
  }

  for (lcp_size=1; lcp_size < M - 1; lcp_size*=2) ;
  lcp_size = lcp_size * 2 + 1;

  if ((lcp = sary->lcp = malloc(lcp_size * sizeof(int))) == NULL) {
    sary_free(sary);
    return NULL;
  }
  for (i=0; i < lcp_size; i++)
    lcp[i] = -1;

  if ((leaves = sary->lcp_leaves = malloc(M * sizeof(int))) == NULL) {
    sary_free(sary);
    return NULL;
  }

  /*
   * Build the suffix tree for the string.
   */
#ifdef STRMAT
  tree = stree_new_tree(128, 0, SORTED_LIST, 0);
#else
  tree = stree_new_tree(128, 0);
#endif
  if (tree == NULL) {
    sary_free(sary);
    return NULL;
  }

#ifdef STRMAT
  flag = stree_ukkonen_add_string(tree, S+1, S+1, M, 1);
#else
  flag = stree_add_string(tree, S+1, M, 1);
#endif
  if (flag <= 0) {
    stree_delete_tree(tree);
    sary_free(sary);
    return NULL;
  }

#ifdef STATS
  sary->num_compares = tree->num_compares;
  sary->num_tree_ops = tree->child_cost;
#endif

  /*
   * Compute the suffix array and lcp leaf values from the suffix tree.
   * Then, compute the lcp tree values.
   *
   * NOTE:  The lcp tree values are placed in the tree so that
   *        during a binary search, the paths down the tree match
   *        the greater-than/less-than choices made by the possible
   *        searches.
   */
  compute_arrays(sary, tree, stree_get_root(tree), 0);

  midpoint = (1 + M) / 2;
  compute_lcp_values(sary, 1, midpoint, 2);
  compute_lcp_values(sary, midpoint, M, 3);

  stree_delete_tree(tree);

  return sary;
}


static void compute_arrays(SARY_STRUCT *sary, SUFFIX_TREE tree,
                           STREE_NODE node, int current_depth)
{
  static int leafnum, min_depth;
  int i, pos, id, edgelen;
  char *str;
  STREE_NODE child;

  if (node == stree_get_root(tree)) {
    leafnum = 1;
    min_depth = 0;
  }

  if (leafnum > 1 && min_depth > current_depth)
    min_depth = current_depth;

  /*
   * For each leaf, fill in the next value of the suffix array and
   * the next value of the last row of the lcp tree.
   */
  for (i=1; stree_get_leaf(tree, node, i, &str, &pos, &id); i++) {
    sary->Pos[leafnum] = pos + 1;

    if (leafnum > 1)
      sary->lcp_leaves[leafnum] = min_depth;

    min_depth = current_depth;
    leafnum++;

#ifdef STATS
    sary->num_lcp_ops++;
#endif
  }

  /*
   * Recurse on the children.
   */
#ifdef STRMAT
  stree_sort_children(tree, node);
#endif
  child = stree_get_children(tree, node);
  while (child != NULL) {
    edgelen = stree_get_edgelen(tree, child);
    compute_arrays(sary, tree, child, current_depth + edgelen);

    if (leafnum > 1 && min_depth > current_depth)
      min_depth = current_depth;

#ifdef STATS
    sary->num_lcp_ops++;
#endif

    child = stree_get_next(tree, child);
  }
}


static int compute_lcp_values(SARY_STRUCT *sary, int min, int max, int index)
{
  int midpoint, value1, value2;

  if (max - min == 1)
    sary->lcp[index] = sary->lcp_leaves[max];
  else {
    midpoint = (min + max) / 2;
    value1 = compute_lcp_values(sary, min, midpoint, index * 2);
    value2 = compute_lcp_values(sary, midpoint, max, index * 2 + 1);
    sary->lcp[index] = (value1 <= value2 ? value1 : value2);

#ifdef STATS
    sary->num_lcp_ops++;
#endif
  }

  return sary->lcp[index];
}


void sary_free(SARY_STRUCT *sary)
{
  if (sary->lcp_leaves != NULL)
    free(sary->lcp_leaves);
  if (sary->lcp != NULL)
    free(sary->lcp);
  if (sary->Pos != NULL)
    free(sary->Pos);
  if (sary->copyflag && sary->S != NULL)
    free(sary->S);

  free(sary);
}
