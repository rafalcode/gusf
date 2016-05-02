/*
 * bmset_naive.c
 *
 * A "naive" implementation of the Boyer-Moore set matching algorithm which
 * simply runs the Boyer-Moore algorithm separately for each of the patterns.
 *
 * NOTES:
 *    9/95  -  Original Implementation (James Knight)
 *    3/96  -  Modularized the code (James Knight)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bmset_naive.h"



/*
 * bmset_naive_alloc
 *
 * Creates a new BMSET_NAIVE_STRUCT structure and initializes its fields.
 *
 * Parameters:    none.
 *
 * Returns:  A dynamically allocated BMSET_NAIVE_STRUCT structure.
 */
BMSET_NAIVE_STRUCT *bmset_naive_alloc(void)
{
  BMSET_NAIVE_STRUCT *node;

  if ((node = malloc(sizeof(BMSET_NAIVE_STRUCT))) == NULL)
    return NULL;
  memset(node, 0, sizeof(BMSET_NAIVE_STRUCT));

  return node;
}


/*
 * bmset_naive_add_string
 *
 * Adds a string to the BMSET_NAIVE_STRUCT structure's keyword tree.
 *
 * NOTE:  The `id' value given must be unique to any of the strings
 *        added to the tree, and should be a small integer greater than
 *        0.
 *
 *        The best id's to use are to number the strings from 1 to K.
 *
 * Parameters:   node      -  a BMSET_NAIVE_STRUCT structure
 *               P         -  the sequence
 *               M         -  the sequence length
 *               id        -  the sequence identifier
 *               copyflag  -  make a copy of the sequence?
 *
 * Returns:  non-zero on success, zero on error.
 */
int bmset_naive_add_string(BMSET_NAIVE_STRUCT *node, char *P, int M, int id,
                           int copyflag)
{
  int i;
  BM_STRUCT *bmstruct;

  if (node->errorflag)
    return 0;

  /*
   * Allocate space for the new string's information.
   */
  if (node->num_patterns == node->patterns_size) {
    if (node->patterns_size == 0) {
      node->patterns_size = 16;
      node->patterns = malloc(node->patterns_size * sizeof(BM_STRUCT *));
      node->ids = malloc(node->patterns_size * sizeof(int));
      node->matches = malloc(node->patterns_size * sizeof(char *));
    }
    else {
      node->patterns_size += node->patterns_size;
      node->patterns = realloc(node->patterns,
                               node->patterns_size * sizeof(BM_STRUCT *));
      node->ids = realloc(node->ids, node->patterns_size * sizeof(int));
      node->matches = realloc(node->matches,
                              node->patterns_size * sizeof(char *));
    }
    if (node->patterns == NULL || node->ids == NULL || node->matches == NULL) {
      node->errorflag = 1;
      return 0;
    }
  }

  /*
   * Check for duplicates.
   */
  for (i=0; i < node->num_patterns; i++) {
    if (node->ids[i] == id) {
      fprintf(stderr, "Error in naive Boyer-Moore Set preprocessing.  "
              "Duplicate identifiers.\n");
      return 0;
    }
  }

  if ((bmstruct = bmgood_prep(P, M, copyflag)) == NULL) {
    fprintf(stderr, "Error in naive Boyer-Moore Set preprocessing.  "
            "Cannot preprocess pattern %d.\n", id);
    return 0;
  }
#ifdef STATS
  node->prep_compares += bmstruct->prep_compares;
#endif

  node->patterns[node->num_patterns] = bmstruct;
  node->ids[node->num_patterns] = id;
  node->num_patterns++;
  node->initflag = 0;
  
  return 1;
}


/*
 * bmset_naive_del_string
 *
 * Deletes a string from the set of patterns.
 *
 * Parameters:   node  -  a BMSET_NAIVE_STRUCT structure
 *               P     -  the sequence to be deleted
 *               M     -  its length
 *               id    -  its identifier
 *
 * Returns:  non-zero on success, zero on error.
 */
int bmset_naive_del_string(BMSET_NAIVE_STRUCT *node, char *P, int M, int id)
{
  int i;

  if (node->errorflag)
    return 0;

  for (i=0; i < node->num_patterns; i++)
    if (node->ids[i] == id)
      break;

  if (i == node->num_patterns) {
    fprintf(stderr, "Error in naive Boyer-Moore Set preprocessing.  "
            "String to be deleted is not in tree.\n");
    return 0;
  }

  bm_free(node->patterns[i]);
  for (i++; i < node->num_patterns; i++) {
    node->patterns[i-1] = node->patterns[i];
    node->ids[i-1] = node->ids[i];
  }

  node->num_patterns--;
  node->initflag = 0;

  return 1;
}


/*
 * bmset_naive_search_init
 *
 * Initializes the variables used during an Boyer-Moore Set search.
 * See bmset_naive_search for an example of how it should be used.
 *
 * Parameters:  node  -  an BMSET_NAIVE_STRUCT structure
 *              T     -  the sequence to be searched
 *              N     -  the length of the sequence
 *
 * Returns:  nothing.
 */
void bmset_naive_search_init(BMSET_NAIVE_STRUCT *node, char *T, int N)
{
  int i;

  if (node->errorflag)
    return;

  node->T = T;
  node->N = N;
  for (i=0; i < node->num_patterns; i++)
    node->matches[i] = NULL;

  node->initflag = node->startflag = 1;
  node->endflag = 0;
  node->output = -1;

#ifdef STATS
  node->num_compares = node->num_shifts = 0;
#endif
}


/*
 * bmset_naive_search
 *
 * Scans a text to look for the next occurrence of one of the patterns
 * in the text.  An example of how this search should be used is the
 * following:
 *
 *    s = T;
 *    len = N; 
 *    contflag = 0;
 *    bmset_naive_search_init(node, T, N);
 *    while ((s = bmset_naive_search(node, &matchlen, &matchid) != NULL) {
 *      >>> Pattern `matchid' matched from `s' to `s + matchlen - 1'. <<<
 *    }
 *
 * where `node', `T' and `N' are assumed to be initialized appropriately.
 *
 * Parameters:  node           -  a preprocessed BMSET_NAIVE_STRUCT structure
 *              length_out     -  where to store the new match's length
 *              id_out         -  where to store the identifier of the
 *                                pattern that matched
 *
 * Returns:  the left end of the text that matches a pattern, or NULL
 *           if no match occurs.  (It also stores values in `*length_out',
 *           and `*id_out' giving the match's length and pattern identifier.
 */
char *bmset_naive_search(BMSET_NAIVE_STRUCT *node, int *length_out,
                         int *id_out)
{
  int i, N, M, len, pos, match, matchlen, matchpos;
  char *s, *T, **matches;

  if (node->errorflag)
    return NULL;
  else if (!node->initflag) {
    fprintf(stderr, "Error in naive Boyer-Moore Set search.  "
            "bmset_naive_search_init was not called.\n");
    return NULL;
  }
  else if (node->endflag)
    return NULL;

  T = node->T;
  N = node->N;
  matches = node->matches;

  /*
   * If just starting, run an initial Boyer-Moore scan for each of the 
   * patterns.
   */
  if (node->startflag) {
    for (i=0; i < node->num_patterns; i++)
      node->matches[i] = bmgood_search(node->patterns[i], T, N, 0);
    node->startflag = 0;
  }

  /*
   * If a previous match has been output, advance that Boyer-Moore scan.
   */
  if (node->output != -1) {
    i = node->output;
    s = node->matches[i];
    len = N - (s - T);
    node->matches[i] = bmgood_search(node->patterns[i], s, len, 1);
  }

  /*
   * Find the leftmost match, using the following criteria:
   *   1) Smallest right endpoint.
   *   2) Break ties by largest length.
   *
   * The use of the smallest right endpoint (which may not seem intuitive
   * when the left endpoint is being returned) is needed for compatibility
   * with Aho-Corasick, so that both algorithms return their answers in
   * the same order.
   */
#ifdef STATS
  node->num_compares = node->num_shifts = 0;
#endif

  match = -1;
  matchlen = 0;
  matchpos = 0;
  for (i=0; i < node->num_patterns; i++) {
    if (node->matches[i]) {
      M = node->patterns[i]->M;
      pos = (node->matches[i] - T) + M;
      if (match == -1 || pos < matchpos || (pos == matchpos && M > matchlen)) {
        match = i;
        matchlen = M;
        matchpos = pos;
      }
    }

#ifdef STATS
    node->num_compares += node->patterns[i]->num_compares;
    node->num_shifts += node->patterns[i]->num_shifts;
#endif
  }

  if (match == -1) {
    node->endflag = 1;
    return NULL;
  }
  else {
    if (id_out)
      *id_out = node->ids[match];
    if (length_out)
      *length_out = matchlen;

    node->output = match;
    return node->matches[match];
  }
}


/*
 * bmset_naive_free
 *
 * Free up the allocated BMSET_NAIVE_STRUCT structure.
 *
 * Parameters:   node  -  a BMSET_NAIVE_STRUCT structure
 *
 * Returns:  nothing.
 */
void bmset_naive_free(BMSET_NAIVE_STRUCT *node)
{
  int i;

  if (node == NULL)
    return;

  if (node->patterns != NULL) {
    for (i=0; i < node->num_patterns; i++)
      bm_free(node->patterns[i]);
    free(node->patterns);
  }
  if (node->ids != NULL)
    free(node->ids);
  if (node->matches != NULL)
    free(node->matches);

  free(node);
}
