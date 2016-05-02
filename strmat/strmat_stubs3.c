
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "strmat.h"
#include "strmat_match.h"
#include "sary.h"
#include "sary_match.h"



static int print_lcp_values(int *lcp, int min, int max, int index, int depth);



int strmat_sary_qsort(STRING *string, int print_stats)
{
  int i, len, M, *Pos;
  char format[32], buf[36];
  SARY_STRUCT *sary;

  if (string == NULL || string->sequence == NULL || string->length == 0)
    return 0;

  /*
   * Compute the suffix array using the qsort algorithm.
   */
  sary = sary_qsort_build(string->sequence, string->length, 0);
  if (sary == NULL)
    return 0;

  /*
   * Print the statistics.
   */
  if (print_stats) {
    mprintf("Statistics:\n");
    mprintf("   Number of Compares:  %d\n", sary->num_compares);
    mprintf("   Number of Tree Ops:  %d\n", sary->num_tree_ops);
    mprintf("   Number of LCP Ops:   %d\n", sary->num_lcp_ops);
    mprintf("\n");
  }

  /*
   * Print the suffix array values.
   */
  mprintf("The Suffix Array:\n");

  len = my_itoalen(string->length);
  sprintf(format, "  %%%dd:  %%s\n", len);

  buf[30] = buf[31] = buf[32] = '.';
  buf[33] = '\0';

  Pos = sary->Pos;
  M = sary->M;
  for (i=1; i <= M; i++) {
    strncpy(buf, &string->raw_seq[Pos[i] - 1], 30);
    if (mprintf(format, Pos[i], buf) == 0)
      return 0;
  }

  sary_free(sary);

  return 1;
}


int strmat_sary_zerkle(STRING *string, int print_stats)
{
  int i, len, M, *Pos;
  char format[32], buf[36];
  SARY_STRUCT *sary;

  if (string == NULL || string->sequence == NULL || string->length == 0)
    return 0;

  /*
   * Compute the suffix array using the qsort algorithm.
   */
  sary = sary_zerkle_build(string->sequence, string->length, 0);
  if (sary == NULL)
    return 0;

  /*
   * Print the suffix array values.
   */
  mprintf("The Suffix Array:\n");

  len = my_itoalen(string->length);
  sprintf(format, "  %%%dd:  %%s\n", len);

  buf[30] = buf[31] = buf[32] = '.';
  buf[33] = '\0';

  Pos = sary->Pos;
  M = sary->M;
  for (i=1; i <= M; i++) {
    strncpy(buf, &string->raw_seq[Pos[i] - 1], 30);
    if (mprintf(format, Pos[i], buf) == 0)
      return 0;
  }

  sary_free(sary);

  return 1;
}


int strmat_sary_stree(STRING *string, int print_stats)
{
  int i, len, midpoint, M, *Pos;
  char format[32], buf[36];
  SARY_STRUCT *sary;

  if (string == NULL || string->sequence == NULL || string->length == 0)
    return 0;

  /*
   * Build the suffix array from the suffix tree.
   */
  sary = sary_stree_build(string->sequence, string->length, 0);
  if (sary == NULL)
    return 0;

  /*
   * Print the statistics.
   */
  if (print_stats) {
    mprintf("Statistics:\n");
    mprintf("   Number of Compares:  %d\n", sary->num_compares);
    mprintf("   Number of Tree Ops:  %d\n", sary->num_tree_ops);
    mprintf("   Number of LCP Ops:   %d\n", sary->num_lcp_ops);
    mprintf("\n\n");
  }

  /*
   * Print the suffix array values.
   */
  mprintf("The Suffix Array:\n");
  len = my_itoalen(string->length);
  sprintf(format, "  %%%dd)  %%%dd:  %%s\n", len, len);

  buf[30] = buf[31] = buf[32] = '.';
  buf[33] = '\0';

  Pos = sary->Pos;
  M = sary->M;
  for (i=1; i <= M; i++) {
    strncpy(buf, &string->raw_seq[Pos[i] - 1], 30);
    if (mprintf(format, i, Pos[i], buf) == 0)
      return 0;
  }
  mprintf("\n");

  /*
   * Print the lcp values.
   */
  mprintf("The LCP Values:\n");
  midpoint = (1 + M) / 2;
  if (print_lcp_values(sary->lcp, 1, midpoint, 2, 3)) {
    print_lcp_values(sary->lcp, midpoint, M, 3, 3);
    mprintf("\n");
  }

  /*
   * Free everything.
   */
  sary_free(sary);

  return 1;
}


static int print_lcp_values(int *lcp, int min, int max, int index, int depth)
{
  int i, midpoint;

  for (i=0; i < depth; i++)
    mputc(' ');

  if (mprintf("%d,%d  =  %d\n", min, max, lcp[index]) == 0)
    return 0;

  if (max - min == 1)
    return 1;

  midpoint = (max + min) / 2;
  if (!print_lcp_values(lcp, min, midpoint, index * 2, depth + 3))
    return 0;
  else
    return print_lcp_values(lcp, midpoint, max, index * 2 + 1, depth + 3);
}



/*
 * strmat_sary_match_naive
 *
 * Performs the exact matching algorithm for a pattern and text using
 * a suffix array and the naive binary search algorithm.
 *
 * Parameters:   pattern  -  the pattern sequence
 *               text     -  the text sequence
 *               stats    -  should stats be printed
 *
 * Returns:  non-zero on success, zero on an error.
 */

static int int_strmat_sary_match(STRING *pattern, STRING *text, int stats,
                                 SARY_MATCH_TYPE type);

int strmat_sary_match_naive(STRING *pattern, STRING *text, int stats)
{  return int_strmat_sary_match(pattern, text, stats, NAIVE_MATCH);  }
int strmat_sary_match_mlr(STRING *pattern, STRING *text, int stats)
{  return int_strmat_sary_match(pattern, text, stats, MLR_MATCH);  }
int strmat_sary_match_lcp(STRING *pattern, STRING *text, int stats)
{  return int_strmat_sary_match(pattern, text, stats, LCP_MATCH);  }

static int int_strmat_sary_match(STRING *pattern, STRING *text, int stats,
                                 SARY_MATCH_TYPE type)
{
  int M, N, pos, flag, matchcount, num_compares;
  char *s, *P, *T;
  MATCHES matchlist, matchtail, newmatch, back, current, next;
  SARYMAT_STRUCT *smstruct;

  if (pattern == NULL || pattern->sequence == NULL || pattern->length == 0 ||
      text == NULL || text->sequence == NULL || text->length == 0)
    return 0;

  P = pattern->sequence;
  M = pattern->length;
  T = text->sequence;
  N = text->length;

  num_compares = 0;
  matchlist = matchtail = NULL;
  matchcount = 0;

  /*
   * "Preprocess" the pattern.
   */
  smstruct = NULL;
  switch (type) {
  case NAIVE_MATCH:  smstruct = sary_match_naive_prep(T, N, 0);  break;
  case MLR_MATCH:    smstruct = sary_match_mlr_prep(T, N, 0);  break;
  case LCP_MATCH:    smstruct = sary_match_lcp_prep(T, N, 0);  break;
  }
  if (smstruct == NULL)
    return 0;

  /*
   * Perform the matching.
   */
  matchlist = matchtail = NULL;
  matchcount = 0;

  s = NULL;
  switch (type) {
  case NAIVE_MATCH:  s = sary_match_naive_first(smstruct, P, M);  break;
  case MLR_MATCH:    s = sary_match_mlr_first(smstruct, P, M);  break;
  case LCP_MATCH:    s = sary_match_lcp_first(smstruct, P, M);  break;
  }

  while (s != NULL) {
    pos = s - T + 1;

    newmatch = alloc_match();
    if (newmatch == NULL) {
      free_matches(matchlist);
      sary_match_free(smstruct);
      mprintf("Memory Error:  Ran out of memory.\n");
      return 0;
    }
    newmatch->type = ONESEQ_EXACT;
    newmatch->lend = pos;
    newmatch->rend = pos + M - 1;

    if (matchlist == NULL)
      matchlist = matchtail = newmatch;
    else {
      matchtail->next = newmatch;
      matchtail = newmatch;
    }
    matchcount++;

    switch (type) {
    case NAIVE_MATCH:  s = sary_match_naive_next(smstruct);  break;
    case MLR_MATCH:    s = sary_match_mlr_next(smstruct);  break;
    case LCP_MATCH:    s = sary_match_lcp_next(smstruct);  break;
    }
  }

  /*
   * Bubble sort the matches.
   */
  if (matchlist != NULL && matchlist->next != NULL) {
    flag = 1;
    while (flag) {
      flag = 0;
      back = NULL;
      current = matchlist;
      while (current->next != NULL) {
        if (current->next->lend < current->lend) {
          /*
           * Move current->next before current in the list.
           */
          next = current->next;
          current->next = next->next;
          next->next = current;
          if (back == NULL)
            back = matchlist = next;
          else
            back = back->next = next;
          
          flag = 1;
        }
        else {
          back = current;
          current = current->next;
        }
      }
    }
  }
    
  /*
   * Print the statistics and the matches.
   */
  print_matches(text, NULL, 0, matchlist, matchcount);

  if (stats) {
    mprintf("Statistics:\n");
#ifdef STATS
    mprintf("  Preprocessing:\n");
    mprintf("     Text Length:         %d\n", N);
    mprintf("     Number of Compares:  %d\n", smstruct->sary->num_compares);
    mprintf("     Number of Tree Ops:  %d\n", smstruct->sary->num_tree_ops);
    mprintf("     Number of LCP Ops:   %d\n", smstruct->sary->num_lcp_ops);
    mprintf("\n");
    mprintf("  Searching:\n");
    mprintf("     Pattern Length:          %d\n", M);
    mprintf("     Number of Compares:      %d\n", smstruct->num_compares);
    mprintf("     Depth of Binary Search:  %d\n", smstruct->search_depth);
    mprintf("     Number of Output Ops:    %d\n", matchcount);
#else
    mputs("   No statistics available.\n");
#endif
    mputc('\n');
  }

  /*
   * Free everything allocated.
   */
  free_matches(matchlist);
  sary_match_free(smstruct);

  return 1;
}


