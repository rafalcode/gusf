
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "strmat.h"
#include "strmat_match.h"
#include "ac.h"
#include "bm.h"
#include "bmset_naive.h"
#include "bmset.h"
#include "kmp.h"
#include "naive.h"
#include "z.h"

/*
 * strmat_naive_match
 *
 * Performs the exact matching algorithm for a pattern and text using
 * the naive search algorithm.
 *
 * Parameters:   pattern  -  the pattern sequence
 *               text     -  the text sequence
 *               stats    -  should stats be printed
 *
 * Returns:  non-zero on success, zero on an error.
 */
int strmat_naive_match(STRING *pattern, STRING *text, int stats)
{
  int M, N, len, pos, flag, matchcount;
  char *s, *P, *T;
  MATCHES matchlist, matchtail, newmatch;
  NAIVE_STRUCT *nstruct;

  if (pattern == NULL || pattern->sequence == NULL || pattern->length == 0 ||
      text == NULL || text->sequence == NULL || text->length == 0)
    return 0;

  P = pattern->sequence;
  M = pattern->length;
  T = text->sequence;
  N = text->length;

  /*
   * "Preprocess" the pattern.
   */
  if ((nstruct = naive_prep(P, M, 0)) == NULL)
    return 0;

  /*
   * Perform the matching.
   */
  matchlist = matchtail = NULL;
  matchcount = 0;

  flag = 0;
  s = T;
  len = N;
  while ((s = naive_search(nstruct, s, len, flag)) != NULL) {
    pos = s - T + 1;

    newmatch = alloc_match();
    if (newmatch == NULL) {
      free_matches(matchlist);
      naive_free(nstruct);
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

    len = N - pos;
    flag = 1;
  }

  /*
   * Print the statistics and the matches.
   */
  print_matches(text, NULL, 0, matchlist, matchcount);

  if (stats) {
    mprintf("Statistics:\n");
#ifdef STATS
    mprintf("   Text Length:                %d\n", N);
    mprintf("   Number of Comparisons:      %d\n", nstruct->num_compares);
    mprintf("   Avg. Compares per Position: %.2f\n",
            (float) nstruct->num_compares / (float) N);
#else
    mputs("   No statistics available.\n");
#endif
    mputc('\n');
  }

  /*
   * Free everything allocated.
   */
  free_matches(matchlist);
  naive_free(nstruct);

  return 1;
}


/*
 * strmat_bm*_match
 *
 * Performs the exact matching algorithm for a pattern and text using
 * one of the four variations of the Boyer-Moore algorithm.
 *
 * Since the four functions are essentially the same, they've been
 * packaged into one internal procedure that is called by each of 
 * the four functions called by the rest of strmat (i.e., the stubs
 * within the stubs).
 *
 * Parameters:   pattern  -  the pattern sequence
 *               text     -  the text sequence
 *               stats    -  should stats be printed
 *
 * Returns:  non-zero on success, zero on an error.
 */
static int internal_bm_match(STRING *pattern, STRING *text, int stats,
                             BMALG_TYPE flag);

int strmat_bmbad_match(STRING *pattern, STRING *text, int stats)
{  return internal_bm_match(pattern, text, stats, BM_BAD);  }
int strmat_bmext_match(STRING *pattern, STRING *text, int stats)
{  return internal_bm_match(pattern, text, stats, BM_EXT);  }
int strmat_bmgood_match(STRING *pattern, STRING *text, int stats)
{  return internal_bm_match(pattern, text, stats, BM_GOOD);  }
int strmat_bmextgood_match(STRING *pattern, STRING *text, int stats)
{  return internal_bm_match(pattern, text, stats, BM_EXTGOOD);  }

static int internal_bm_match(STRING *pattern, STRING *text, int stats,
                             BMALG_TYPE flag)
{
  int M, N, len, pos, matchcount;
  char *s, *P, *T;
  MATCHES matchlist, matchtail, newmatch;
  BM_STRUCT *bmstruct;

  if (pattern == NULL || pattern->sequence == NULL || pattern->length == 0 ||
      text == NULL || text->sequence == NULL || text->length == 0)
    return 0;

  P = pattern->sequence;
  M = pattern->length;
  T = text->sequence;
  N = text->length;

  matchlist = matchtail = NULL;
  matchcount = 0;

  /*
   * Do the preprocessing.
   */
  bmstruct = NULL;
  switch (flag) {
  case BM_BAD:     bmstruct = bmbad_prep(P, M, 0);  break;
  case BM_EXT:     bmstruct = bmext_prep(P, M, 0);  break;
  case BM_GOOD:    bmstruct = bmgood_prep(P, M, 0);  break;
  case BM_EXTGOOD: bmstruct = bmextgood_prep(P, M, 0);  break;
  }
  if (bmstruct == NULL) {
    fprintf(stderr, "Error in Boyer-Moore preprocessing.  Stopping search.\n");
    return 0;
  }

  /*
   * Perform the matching.
   */
  s = T;
  len = N;

  switch (flag) {
  case BM_BAD:     s = bmbad_search(bmstruct, s, len, 0);  break;
  case BM_EXT:     s = bmext_search(bmstruct, s, len, 0);  break;
  case BM_GOOD:    s = bmgood_search(bmstruct, s, len, 0);  break;
  case BM_EXTGOOD: s = bmextgood_search(bmstruct, s, len, 0);  break;
  }
  while (s != NULL) {
    pos = s - T + 1;

    newmatch = alloc_match();
    if (newmatch == NULL) {
      free_matches(matchlist);
      bm_free(bmstruct);
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

    len = N - pos;
    switch (flag) {
    case BM_BAD:     s = bmbad_search(bmstruct, s, len, 1);  break;
    case BM_EXT:     s = bmext_search(bmstruct, s, len, 1);  break;
    case BM_GOOD:    s = bmgood_search(bmstruct, s, len, 1);  break;
    case BM_EXTGOOD: s = bmextgood_search(bmstruct, s, len, 1);  break;
    }
  }

  /*
   * Print the matches and the statistics.
   */
  print_matches(text, NULL, 0, matchlist, matchcount);

  if (stats) {
    mprintf("Statistics:\n");
#ifdef STATS
    mprintf("   Preprocessing Comparisons:  %d\n", bmstruct->prep_compares);
    mprintf("\n");
    mprintf("   Text Length:                %d\n", N);
    mprintf("   Number of Comparisons:      %d\n", bmstruct->num_compares);
    mprintf("   Avg. Compares per Position: %.2f\n",
            (float) bmstruct->num_compares / (float) N);
    mprintf("\n");
    mprintf("   Number of Init. Mismatches: %d\n",
            bmstruct->num_init_mismatch);
    if (bmstruct->num_shifts > bmstruct->num_init_mismatch)
      mprintf("   Average Length of Matches:  %.2f\n",
              (float) (bmstruct->num_compares - 
                       bmstruct->num_shifts + matchcount) / 
              (float) (bmstruct->num_shifts - bmstruct->num_init_mismatch));
    mprintf("   Number of Shifts:           %d\n", bmstruct->num_shifts);
    if (bmstruct->num_shifts != bmstruct->shift_cost)
      mprintf("   Cost of Computing Shifts:   %d\n", bmstruct->shift_cost);
    if (bmstruct->num_shifts > 0)
      mprintf("   Average Shift Length:       %.2f\n", 
              (float) N / (float) bmstruct->num_shifts);
#else
    mprintf("   No statistics available.\n");
#endif
    mputc('\n');
  }

  /*
   * Free everything allocated.
   */
  free_matches(matchlist);
  bm_free(bmstruct);

  return 1;
}


/*
 * strmat_kmp*_match
 *
 * Performs the exact matching algorithm for a pattern and text using
 * one of the four variations of the Knuth-Morris-Pratt algorithm.
 *
 * Since the four functions are essentially the same, they've been
 * packaged into one internal procedure that is called by each of 
 * the four functions called by the rest of strmat (i.e., the stubs
 * within the stubs).
 *
 * Parameters:   pattern  -  the pattern sequence
 *               text     -  the text sequence
 *               stats    -  should stats be printed
 *
 * Returns:  non-zero on success, zero on an error.
 */
typedef enum { SP_Z, SPPRIME_Z, SP_ORIG, SPPRIME_ORIG } kmp_types;
static int internal_kmp_match(STRING *pattern, STRING *text,
                              int stats, kmp_types flag);

int strmat_kmp_sp_z_match(STRING *pattern, STRING *text, int stats)
{  return internal_kmp_match(pattern, text, stats, SP_Z);  }
int strmat_kmp_spprime_z_match(STRING *pattern, STRING *text, int stats)
{  return internal_kmp_match(pattern, text, stats, SPPRIME_Z);  }
int strmat_kmp_sp_orig_match(STRING *pattern, STRING *text, int stats)
{  return internal_kmp_match(pattern, text, stats, SP_ORIG);  }
int strmat_kmp_spprime_orig_match(STRING *pattern, STRING *text, int stats)
{  return internal_kmp_match(pattern, text, stats, SPPRIME_ORIG);  }

static int internal_kmp_match(STRING *pattern, STRING *text,
                              int stats, kmp_types flag)
{
  int M, N, len, pos, matchcount, matchflag;
  char *s, *P, *T;
  MATCHES matchlist, matchtail, newmatch;
  KMP_STRUCT *kmpstruct;

  if (pattern == NULL || pattern->sequence == NULL || pattern->length == 0 ||
      text == NULL || text->sequence == NULL || text->length == 0)
    return 0;

  P = pattern->sequence;
  M = pattern->length;
  T = text->sequence;
  N = text->length;

  matchlist = matchtail = NULL;
  matchcount = 0;

  /*
   * Do the preprocessing.
   */
  kmpstruct = NULL;
  switch (flag) {
  case SP_Z:          kmpstruct = kmp_sp_z_prep(P, M, 0);  break;
  case SPPRIME_Z:     kmpstruct = kmp_spprime_z_prep(P, M, 0);  break;
  case SP_ORIG:       kmpstruct = kmp_sp_orig_prep(P, M, 0);  break;
  case SPPRIME_ORIG:  kmpstruct = kmp_spprime_orig_prep(P, M, 0);  break;
  }
  if (kmpstruct == NULL) {
    fprintf(stderr, "Error in Knuth-Morris-Pratt preprocessing."
                    "  Stopping search.\n");
    return 0;
  }

  /*
   * Perform the matching.
   */
  s = T;
  len = N;
  matchflag = 0;
  while ((s = kmp_search(kmpstruct, s, len, matchflag)) != NULL) {
    pos = s - T + 1;

    newmatch = alloc_match();
    if (newmatch == NULL) {
      free_matches(matchlist);
      kmp_free(kmpstruct);
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

    len = N - pos;
    matchflag = 1;
  }

  /*
   * Print the statistics and the matches.
   */
  print_matches(text, NULL, 0, matchlist, matchcount);

  if (stats) {
    mprintf("Statistics:\n");
#ifdef STATS
    mprintf("   Preprocessing Comparisons:  %d\n", kmpstruct->prep_compares);
    mprintf("\n");
    mprintf("   Text Length:                %d\n", N);
    mprintf("   Number of Comparisons:      %d\n", kmpstruct->num_compares);
    mprintf("   Avg. Compares per Position: %.2f\n",
            (float) kmpstruct->num_compares / (float) N);
    mprintf("\n");
    mprintf("   Number of Init. Mismatches: %d\n",
            kmpstruct->num_init_mismatch);
    mprintf("   Number of Failures:         %d\n", kmpstruct->num_shifts);
    if (kmpstruct->num_shifts > 0)
      mprintf("   Avg. Failure Distance:      %.2f\n",
              (float) kmpstruct->total_shifts / (float) kmpstruct->num_shifts);
#else
    mprintf("   No statistics available.\n");
#endif
    mputc('\n');
  }

  /*
   * Free everything allocated.
   */
  free_matches(matchlist);
  kmp_free(kmpstruct);

  return 1;
}
  

/*
 * strmat_ac_match
 *
 * Performs the Aho-Corasick set matching algorithm for a set of patterns
 * and a text.
 *
 * Parameters:   pattern_ary   -  the pattern sequences
 *               num_patterns  -  the number of patterns
 *               text          -  the text sequence
 *               stats         -  should stats be printed
 *
 * Returns:  non-zero on success, zero on an error.
 */
int strmat_ac_match(STRING **pattern_ary, int num_patterns, STRING *text,
                    int stats)
{
  int i, M, N, pos, matchcount, matchlen, matchid, total_length;
  char *s, *P, *T;
  MATCHES matchlist, matchtail, newmatch;
  AC_STRUCT *acstruct;

  if (pattern_ary == NULL || num_patterns == 0 || text == NULL ||
      text->sequence == NULL || text->length == 0)
    return 0;

  T = text->sequence;
  N = text->length;

  matchlist = matchtail = NULL;
  matchcount = 0;

  /*
   * Perform the Aho-Corasick preprocessing.
   */
  if ((acstruct = ac_alloc()) == NULL)
    return 0;

  total_length = 0;
  for (i=0; i < num_patterns; i++) {
    P = pattern_ary[i]->sequence;
    M = pattern_ary[i]->length;
    total_length += M;

    if (ac_add_string(acstruct, P, M, i+1) == 0) {
      ac_free(acstruct);
      return 0;
    }
  }
  ac_prep(acstruct);

  /*
   * Perform the matching.
   */
  ac_search_init(acstruct, T, N);
  while ((s = ac_search(acstruct, &matchlen, &matchid)) != NULL) {
    pos = s - T + 1;

    newmatch = alloc_match();
    if (newmatch == NULL) {
      free_matches(matchlist);
      ac_free(acstruct);
      mprintf("Memory Error:  Ran out of memory.\n");
      return 0;
    }
    newmatch->type = SET_EXACT;
    newmatch->lend = pos;
    newmatch->rend = pos + matchlen - 1;
    newmatch->id = matchid;

    if (matchlist == NULL)
      matchlist = matchtail = newmatch;
    else {
      matchtail->next = newmatch;
      matchtail = newmatch;
    }
    matchcount++;
  }

  /*
   * Print the statistics and the matches.
   */
  print_matches(text, NULL, 0, matchlist, matchcount);

  if (stats) {
    mprintf("Statistics:\n");
#ifdef STATS
    mprintf("   Preprocessing:\n");
    mprintf("      Sum of Pattern Sizes:       %d\n", total_length);
    mprintf("      Number of Created Edges:    %d\n",
            acstruct->prep_new_edges);
    mprintf("      Number of Traversed Edges:  %d\n",
            acstruct->prep_old_edges);
    mprintf("      Failure Link Comparisons:   %d\n",
            acstruct->prep_fail_compares);

    mprintf("\n   Searching:\n");
    mprintf("      Text Length:                %d\n", N);
    mprintf("      Number of Compares:         %d\n", acstruct->num_compares);
    mprintf("      Avg. Compares per Position: %.2f\n",
            (float) acstruct->num_compares / (float) N);
    mprintf("\n");
    mprintf("      Tree Edges Traversed:       %d\n",
            acstruct->edges_traversed);
    mprintf("      Fail Links Traversed:       %d\n", acstruct->num_failures);
    mprintf("      Output Link Traversed:      %d\n",
            acstruct->outlinks_traversed);
#else
    mputs("   No statistics available.\n");
#endif
    mputc('\n');
  }

  /*
   * Free everything allocated.
   */
  free_matches(matchlist);
  ac_free(acstruct);

  return 1;
}


/*
 * strmat_bmset_naive_match
 *
 * Performs the exact matching of a set of patterns and a text by using
 * separate Boyer-Moore matchers for each pattern in the set.  The
 * Boyer-Moore matchers are run in parallel, so that the output matches
 * are produced in order.
 *
 * Parameters:   pattern_ary   -  the pattern sequences
 *               num_patterns  -  the number of patterns
 *               text          -  the text sequence
 *               stats         -  should stats be printed
 *
 * Returns:  non-zero on success, zero on an error.
 */
int strmat_bmset_naive_match(STRING **pattern_ary, int num_patterns,
                             STRING *text, int stats)
{
  int i, M, N, pos, matchcount, matchlen, matchid, total_length;
  char *s, *P, *T;
#ifdef STATS
  float avg;
#endif
  MATCHES matchlist, matchtail, newmatch;
  BMSET_NAIVE_STRUCT *bmstruct;

  if (pattern_ary == NULL || num_patterns == 0 || text == NULL ||
      text->sequence == NULL || text->length == 0)
    return 0;

  T = text->sequence;
  N = text->length;

  /*
   * Perform the naive Boyer-Moore Set preprocessing.
   */
  if ((bmstruct = bmset_naive_alloc()) == NULL)
    return 0;

  total_length = 0;
  for (i=0; i < num_patterns; i++) {
    P = pattern_ary[i]->sequence;
    M = pattern_ary[i]->length;
    total_length += M;

    if (bmset_naive_add_string(bmstruct, P, M, i+1, 0) == 0) {
      bmset_naive_free(bmstruct);
      return 0;
    }
  }

  /*
   * Perform the matching.
   */
  matchlist = matchtail = NULL;
  matchcount = 0;

  bmset_naive_search_init(bmstruct, T, N);
  while ((s = bmset_naive_search(bmstruct, &matchlen, &matchid)) != NULL) {
    pos = s - T + 1;

    newmatch = alloc_match();
    if (newmatch == NULL) {
      free_matches(matchlist);
      bmset_naive_free(bmstruct);
      mprintf("Memory Error:  Ran out of memory.\n");
      return 0;
    }
    newmatch->type = SET_EXACT;
    newmatch->lend = pos;
    newmatch->rend = pos + matchlen - 1;
    newmatch->id = matchid;

    if (matchlist == NULL)
      matchlist = matchtail = newmatch;
    else {
      matchtail->next = newmatch;
      matchtail = newmatch;
    }
    matchcount++;
  }

  /*
   * Print the matches and the statistics.
   */
  print_matches(text, NULL, 0, matchlist, matchcount);

  if (stats) {
    mprintf("Statistics:\n");
#ifdef STATS
    mprintf("   Preprocessing:\n");
    mprintf("      Sum of Pattern Sizes:           %d\n", total_length);
    mprintf("      Preprocessing Comparisons:      %d\n",
            bmstruct->prep_compares);

    mprintf("\n   Searching:\n");
    mprintf("      Text Length:                    %d\n", N);
    mprintf("      Number of Compares:             %d\n",
            bmstruct->num_compares);
    mprintf("      Avg. Compares per Position:     %.2f\n",
            (float) bmstruct->num_compares / (float) N);
    mprintf("\n");
    mprintf("      Number of Shifts:               %d\n",
            bmstruct->num_shifts);
    if (bmstruct->num_shifts > 0) {
      avg = 0.0;
      for (i=0; i < num_patterns; i++)
        avg += (float) N / (float) bmstruct->patterns[i]->num_shifts;

      mprintf("      Avg. Shift Length per Pattern:  %.2f\n",
              avg / (float) num_patterns);
    }
#else
    mputs("   No statistics available.\n");
#endif
    mputc('\n');
  }

  /*
   * Free everything allocated.
   */
  free_matches(matchlist);
  bmset_naive_free(bmstruct);

  return 1;
}


/*
 * strmat_bmset_match
 *
 * Performs the Boyer-Moore Set matching algorithm for a set of patterns
 * and a text.
 *
 * Parameters:   pattern_ary   -  the pattern sequences
 *               num_patterns  -  the number of patterns
 *               text          -  the text sequence
 *               stats         -  should stats be printed
 *
 * Returns:  non-zero on success, zero on an error.
 */

static int int_strmat_bmset_match(STRING **patterns, int num_patterns,
                                  STRING *text, int stats, BMSETALG_TYPE type);

int strmat_bmset_badonly_match(STRING **patterns, int num_patterns,
                               STRING *text, int stats)
{  return int_strmat_bmset_match(patterns, num_patterns,
                                 text, stats, BMSET_BADONLY);  }
int strmat_bmset_2trees_match(STRING **patterns, int num_patterns,
                              STRING *text, int stats)
{  return int_strmat_bmset_match(patterns, num_patterns,
                                 text, stats, BMSET_2TREES);  }
int strmat_bmset_1tree_match(STRING **patterns, int num_patterns,
                             STRING *text, int stats)
{  return int_strmat_bmset_match(patterns, num_patterns,
                                 text, stats, BMSET_1TREE);  }

static int int_strmat_bmset_match(STRING **patterns, int num_patterns,
                                  STRING *text, int stats, BMSETALG_TYPE type)
{
  int i, M, N, pos, status, matchcount, matchlen, matchid, total_length;
  char *s, *P, *T;
  MATCHES matchlist, matchtail, newmatch;
  BMSET_STRUCT *bmstruct;

  if (patterns == NULL || num_patterns == 0 || text == NULL ||
      text->sequence == NULL || text->length == 0)
    return 0;

  T = text->sequence;
  N = text->length;

  matchlist = matchtail = NULL;
  matchcount = 0;

  /*
   * Perform the Boyer-Moore set matching preprocessing.
   */
  bmstruct = NULL;
  switch (type) {
  case BMSET_BADONLY:  bmstruct = bmset_badonly_alloc();  break;
  case BMSET_2TREES:   bmstruct = bmset_2trees_alloc();  break;
  case BMSET_1TREE:    bmstruct = bmset_1tree_alloc();  break;
  }
  if (bmstruct == NULL)
    return 0;

  total_length = 0;
  for (i=0; i < num_patterns; i++) {
    P = patterns[i]->sequence;
    M = patterns[i]->length;
    total_length += M;

    status = 0;
    switch (type) {
    case BMSET_BADONLY:  status = bmset_badonly_add_string(bmstruct,
                                                           P, M, i+1);  break;
    case BMSET_2TREES:   status = bmset_2trees_add_string(bmstruct,
                                                          P, M, i+1);  break;
    case BMSET_1TREE:    status = bmset_1tree_add_string(bmstruct,
                                                         P, M, i+1);  break;
    }
    if (status == 0) {
      switch (type) {
      case BMSET_BADONLY:  bmset_badonly_free(bmstruct);  break;
      case BMSET_2TREES:   bmset_2trees_free(bmstruct);  break;
      case BMSET_1TREE:    bmset_1tree_free(bmstruct);  break;
      }
      return 0;
    }
  }

  switch (type) {
  case BMSET_BADONLY:  bmset_badonly_prep(bmstruct);  break;
  case BMSET_2TREES:   bmset_2trees_prep(bmstruct);  break;
  case BMSET_1TREE:    bmset_1tree_prep(bmstruct);  break;
  }

  /*
   * Perform the matching.
   */
  switch (type) {
  case BMSET_BADONLY:  bmset_badonly_search_init(bmstruct, T, N);  break;
  case BMSET_2TREES:   bmset_2trees_search_init(bmstruct, T, N);  break;
  case BMSET_1TREE:    bmset_1tree_search_init(bmstruct, T, N);  break;
  }

  s = NULL;
  switch (type) {
  case BMSET_BADONLY:  s = bmset_badonly_search(bmstruct,
                                                &matchlen, &matchid);  break;
  case BMSET_2TREES:   s = bmset_2trees_search(bmstruct,
                                               &matchlen, &matchid);  break;
  case BMSET_1TREE:    s = bmset_1tree_search(bmstruct,
                                              &matchlen, &matchid);  break;
  }
  while (s != NULL) {
    pos = s - T + 1;

    newmatch = alloc_match();
    if (newmatch == NULL) {
      free_matches(matchlist);
      switch (type) {
      case BMSET_BADONLY:  bmset_badonly_free(bmstruct);  break;
      case BMSET_2TREES:   bmset_2trees_free(bmstruct);  break;
      case BMSET_1TREE:    bmset_1tree_free(bmstruct);  break;
      }
      mprintf("Memory Error:  Ran out of memory.\n");
      return 0;
    }
    newmatch->type = SET_EXACT;
    newmatch->lend = pos;
    newmatch->rend = pos + matchlen - 1;
    newmatch->id = matchid;

    if (matchlist == NULL)
      matchlist = matchtail = newmatch;
    else {
      matchtail->next = newmatch;
      matchtail = newmatch;
    }
    matchcount++;

    switch (type) {
    case BMSET_BADONLY:  s = bmset_badonly_search(bmstruct,
                                                  &matchlen, &matchid);  break;
    case BMSET_2TREES:   s = bmset_2trees_search(bmstruct,
                                                 &matchlen, &matchid);  break;
    case BMSET_1TREE:    s = bmset_1tree_search(bmstruct,
                                                &matchlen, &matchid);  break;
    }
  }

  /*
   * Print the matches and the statistics.
   */
  print_matches(text, NULL, 0, matchlist, matchcount);

  if (stats) {
    mprintf("Statistics:\n");
#ifdef STATS
    mprintf("   Preprocessing:\n");
    mprintf("       Sum of Pattern Sizes:  %d\n", bmstruct->totallen);
    mprintf("       Tree Construction:     %d\n", bmstruct->prep_tree_ops);
    mprintf("       Value Computation:     %d\n", bmstruct->prep_value_ops);
    mprintf("\n");
    mprintf("   Search:\n"); 
    mprintf("       Text Length:                %d\n", N);
    mprintf("       Number of Comparisons:      %d\n", bmstruct->num_compares);
    mprintf("       Avg. Compares per Position: %.2f\n",
            (float) bmstruct->num_compares / (float) N);
    mprintf("\n");
    mprintf("       Num. Tree Edges Traversed:  %d\n", bmstruct->num_edges);
    mprintf("       Cost of Tree Traversal:     %d\n", bmstruct->edge_cost);
    mprintf("\n");
    mprintf("       Number of Init. Mismatches: %d\n",
            bmstruct->num_init_mismatch);
    if (bmstruct->num_shifts > bmstruct->num_init_mismatch)
      mprintf("       Average Length of Matches:  %.2f\n",
              (float) bmstruct->total_matches / 
              (float) (bmstruct->num_shifts - bmstruct->num_init_mismatch));
    mprintf("\n");
    mprintf("       Number of Shifts:           %d\n", bmstruct->num_shifts);
    mprintf("       Cost of Computing Shifts:   %d\n", bmstruct->shift_cost);
    if (bmstruct->num_shifts > 0)
      mprintf("       Average Shift Length:       %.2f\n", 
              (float) N / (float) bmstruct->num_shifts);
#else
    mputs("   No statistics available.\n");
#endif
    mputc('\n');
  }

  /*
   * Free everything allocated.
   */
  free_matches(matchlist);
  switch (type) {
  case BMSET_BADONLY:  bmset_badonly_free(bmstruct);  break;
  case BMSET_2TREES:   bmset_2trees_free(bmstruct);  break;
  case BMSET_1TREE:    bmset_1tree_free(bmstruct);  break;
  }

  return 1;
}


/*
 * strmat_z_build
 *
 * Build the Z values for a string and output those values (and stats).
 *
 * Parameters:   string  -  a sequence
 *               stats   -  should stats be printed
 *
 * Returns:  non-zero on success, zero on an error.
 */
int strmat_z_build(STRING *string, int stats)
{
  int i, len, *Z;
  char format[64], buffer[64];
  Z_STRUCT *zvalues;

  if (string == NULL || string->sequence == NULL || string->length == 0)
    return 0;

  /*
   * Compute the Z values for the sequence.
   */
  zvalues = z_build(string->sequence, string->length, 0);
  if (zvalues == NULL) {
    fprintf(stderr, "Error during Z values construction.\n");
    return 0;
  }

  /*
   * Print the non-zero Z values.
   */
  Z = zvalues->Z;

  sprintf(format, "%d", string->length);
  len = strlen(format);
  sprintf(format, "   Position %%%dd: %%%dd - %%s\n", len, len);

  buffer[30] = '.';
  buffer[31] = '.';
  buffer[32] = '.';
  buffer[33] = '\0';

  mprintf("Z Values (non-zero values only):\n");
  for (i=2; i <= string->length; i++) {
    if (Z[i] == 0)
      continue;

    if (Z[i] < 30) {
      strncpy(buffer, string->raw_seq + i - 1, Z[i]);
      buffer[Z[i]] = '\0';
    }
    else
      strncpy(buffer, string->raw_seq + i - 1, 30);
      
    if (mprintf(format, i, Z[i], buffer) == 0) {
      z_free(zvalues);
      return 1;
    }
  }
  mputc('\n');

  /*
   * Print the statistics for the sequence.
   */
  if (stats) {
    mprintf("Statistics:\n");
#ifdef STATS
    mprintf("   Text Length:           %d\n", string->length);
    mprintf("   Number of Comparisons: %d\n\n", zvalues->prep_compares);
#else
    mprintf("   No statistics available.\n\n");
#endif
  }

  z_free(zvalues);

  return 1;
}



/*
 * strmat_z_match
 *
 * Performs the exact matching algorithm for a pattern and text using
 * the Z values matching algorithm.
 *
 * Parameters:   pattern  -  the pattern sequence
 *               text     -  the text sequence
 *               stats    -  should stats be printed
 *
 * Returns:  non-zero on success, zero on an error.
 */
int strmat_z_match(STRING *pattern, STRING *text, int stats)
{
  int M, N, len, pos, flag, matchcount;
  char *s, *P, *T;
  MATCHES matchlist, matchtail, newmatch;
  Z_STRUCT *zvalues;

  if (pattern == NULL || pattern->sequence == NULL || pattern->length == 0 ||
      text == NULL || text->sequence == NULL || text->length == 0)
    return 0;

  P = pattern->sequence;
  M = pattern->length;
  T = text->sequence;
  N = text->length;

  /*
   * Compute the Z values of the pattern.
   */
  if ((zvalues = z_build(P, M, 0)) == NULL)
    return 0;

  /*
   * Perform the matching.
   */
  matchlist = matchtail = NULL;
  matchcount = 0;
  flag = 0;
  s = T;
  len = N;
  while ((s = z_search(zvalues, s, len, flag)) != NULL) {
    pos = s - T + 1;

    newmatch = alloc_match();
    if (newmatch == NULL) {
      free_matches(matchlist);
      z_free(zvalues);
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

    len = N - pos;
    flag = 1;
  }

  /*
   * Print the statistics and the matches.
   */
  print_matches(text, NULL, 0, matchlist, matchcount);

  if (stats) {
    mprintf("Statistics:\n");
#ifdef STATS
    mprintf("   Pattern Length:             %d\n", M);
    mprintf("   Preprocessing Comparisons:  %d\n", zvalues->prep_compares);
    mprintf("   Text Length:                %d\n", N);
    mprintf("   Number of Comparisons:      %d\n", zvalues->num_compares);
    mprintf("   Avg. Compares per Position: %.2f\n",
            (float) zvalues->num_compares / (float) N);
#else
    mputs("   No statistics available.\n");
#endif
    mputc('\n');
  }

  /*
   * Free everything allocated.
   */
  free_matches(matchlist);
  z_free(zvalues);

  return 1;
}
