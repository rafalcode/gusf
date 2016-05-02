
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "strmat.h"
#include "strmat_alpha.h"
#include "strmat_match.h"
#include "strmat_print.h"
#include "stree_strmat.h"
#include "stree_ukkonen.h"
#include "repeats_primitives.h"
#include "repeats_supermax.h"
#include "repeats_nonoverlapping.h"
#include "repeats_bigpath.h"
#include "repeats_tandem.h"
#include "repeats_vocabulary.h"
#include "repeats_linear_occs.h"
#include "strmat_stubs4.h"



/*
 * strmat_repeats_primitives
 *
 * Find all primitive tandem repeats with Crochemore's algorithm.
 *
 * Parameters:   string       -  the input string
 *               print_stats  -  flag telling whether to print the stats
 *
 * Returns:  non-zero on success, zero on error
 */
int strmat_repeats_primitives(STRING *string, int print_stats)
{
  primitives_struct *pr_struct;

  if (string == NULL || string->sequence == NULL || string->length == 0)
    return 0;

  /*
   * Allocate the memory.
   */
  mprintf("Preprocessing...\n");
  pr_struct = primitives_prep(string->sequence,string->raw_seq,string->length);
  if(pr_struct == NULL)
    return 0;

  /*
   * Report the results.
   */
  mprintf("\nThe following primitive tandem repeats were found:\n\n");
  primitives_find(pr_struct);
  mend(17);

  /*
   * Write summary of results.
   */
  mstart(stdin, stdout, OK, OK, 0, NULL);
  mprintf("\nSummary:\n");
  mprintf("   Primitive Tandem Repeat Occurrences: %d\n",
          pr_struct->num_primitive_tandem_repeat_occs);

  /*
   * Write statistics and free memory
   */
  if (print_stats) {
    mprintf("\nStatistics:\n");
#ifdef STATS
    mprintf("   String Length:          %d\n", string->length);
    mprintf("   Number of Compares:     %d\n", pr_struct->num_compares);
#else
    mprintf("   No statistics available.\n");
#endif
    mprintf("\n");
  }

  primitives_free(pr_struct);

  return 1;
}


/*
 * strmat_repeats_supermax
 *
 * Find the supermaximals for a string.
 *
 * Parameters:  string       -  the string
 *              min_percent  -  min percent for any reported supermaximal
 *              min_length   -  min length for any reported supermaximal
 *
 * Returns:  non-zero on success, zero on error.
 */
int strmat_repeats_supermax(STRING *string, int min_percent, int min_length)
{
  int pos;
  char buffer[64];
  SUPERMAXIMALS list, next;

  if (string == NULL || string->sequence == NULL || string->length == 0)
    return 0;

  /*
   * Compute the supermaximals, print them out and free the space.
   */
  list = supermax_find(string->sequence, string->length, min_percent,
                       min_length);

  mprintf("Supermaximals:\n");
  if (list == NULL)
    mprintf("   none\n");
  else {
    for ( ; list != NULL; list=next) {
      next = list->next;
      pos = list->S - string->sequence;
      if (list->M <= 50) {
        strncpy(buffer, string->raw_seq + pos, list->M);
        buffer[list->M] = '\0';
      }
      else {
        memcpy(buffer, string->raw_seq + pos, 50);
        buffer[50] = buffer[51] = buffer[52] = '.';
        buffer[53] = '\0';
      }

      mprintf("   %s    %d/%d  %d%%\n", buffer, list->num_witness,
              list->num_leaves, list->percent);

      free(list);
    }
  }

  return 1;
}


/*
 * strmat_repeats_nonoverlapping
 *
 * Find all nonoverlapping repeats with the extended version of Crochemore's
 * algorithm.
 *
 * Parameters:   string       -  the input string
 *               print_stats  -  flag telling whether to print the stats
 *
 * Returns:  non-zero on success, zero on error
 */
int strmat_repeats_nonoverlapping(STRING *string, int print_stats)
{
  nonoverlapping_struct *no_struct;

  if (string == NULL || string->sequence == NULL || string->length == 0)
    return 0;

  /*
   * Allocate the memory.
   */
  mprintf("Preprocessing...\n");
  no_struct = nonoverlapping_prep(string->sequence,string->raw_seq,
                                  string->length);
  if(no_struct == NULL)
    return 0;

  /*
   * Report the results.
   */
  mprintf("\nThe following nonoverlapping maximal pairs were found:\n\n");
  nonoverlapping_find(no_struct);
  mend(17);

  /*
   * Write summary of results.
   */
  mstart(stdin, stdout, OK, OK, 0, NULL);
  mprintf("\nSummary:\n");
  mprintf("   Nonoverlapping Maximal Pairs: %d\n",
          no_struct->num_nonoverlapping_maximal_pairs);

  /*
   * Write statistics and free memory
   */
  if (print_stats) {
    mprintf("\nStatistics:\n");
#ifdef STATS
    mprintf("   String Length:          %d\n", string->length);
    mprintf("   Number of Compares:     %d\n", no_struct->num_compares);
#else
    mprintf("   No statistics available.\n");
#endif
    mprintf("\n");
  }

  nonoverlapping_free(no_struct);

  return 1;
}


/*
 * strmat_repeats_bigpath
 *
 * Find all maximal pairs which are a) nonoverlapping
 *                               or b) max-gapped
 * using the algorithm by Pedersen and Stoye.
 *
 * Parameters:   string           -  the input string
 *               build_policy     -  suffix tree build policy
 *               build_threshold  -  threshold used by LIST_THEN_ARRAY
 *               print_stats      -  flag telling whether to print the stats
 *
 * Returns:  non-zero on success, zero on error
 */
int strmat_repeats_bigpath(STRING *string, int build_policy,
                          int build_threshold, int print_stats)
{
  SUFFIX_TREE tree;
  bp_struct *b;

  if (string == NULL || string->sequence == NULL || string->length == 0)
    return 0;

  /*
   * Build the tree (can't use stree_ukkonen_build because of copyflag).
   */
  mprintf("\nBuilding the suffix tree...\n");
  if ((tree = stree_new_tree(string->alpha_size, 0, build_policy,
                             build_threshold)) == NULL)
    return 0;

  if (stree_ukkonen_add_string(tree, string->sequence, string->raw_seq,
                               string->length, 1) <= 0) {
    stree_delete_tree(tree);
    return 0;
  }

  /*
   * Preprocess the suffix tree and initialize repeat_struct
   */
  mprintf("Preprocessing...\n");
  b = bigpath_prep(tree, string->sequence, string->raw_seq, string->length);
  if(b == NULL) {
    stree_delete_tree(tree);
    return 0;
  }

  /*
   * Report the results.
   */
  mprintf("\nThe following nonoverlapping maximal pairs were found:\n\n");
  bigpath_find(b);
  mend(17);

  /*
   * Write summary of results.
   */
  mstart(stdin, stdout, OK, OK, 0, NULL);
  mprintf("\nSummary:\n");
  mprintf("   Nonoverlapping Maximal Pairs: %d\n",
          b->num_nonoverlapping_maximal_pairs);

  /*
   * Write statistics and free memory
   */
  if (print_stats) {
    mprintf("\nStatistics:\n");
#ifdef STATS
    mprintf("   String Length:        %d\n", string->length);
    mprintf("   Number of Compares\n");
    mprintf("     for suffix tree:    %d\n", tree->num_compares);
    mprintf("     for maximal pairs:  %d\n", b->num_compares);
#else
    mprintf("   No statistics available.\n");
#endif
    mprintf("\n");
  }

  bigpath_free(b);
  stree_delete_tree(tree);

  return 1;
}


/*
 * strmat_repeats_tandem
 *
 * Find all (branching) (maximal) occurrences of tandem arrays in a string.
 *
 * Parameters:   string           -  the input string
 *               build_policy     -  suffix tree build policy
 *               build_threshold  -  threshold used by LIST_THEN_ARRAY
 *               print_stats      -  flag telling whether to print the stats
 *
 * Returns:  non-zero on success, zero on error
 */
int strmat_repeats_tandem(STRING *string, int build_policy,
                          int build_threshold, int print_stats)
{
  SUFFIX_TREE tree;
  TANDEM_STRUCT *tandem_struct;

  if (string == NULL || string->sequence == NULL || string->length == 0)
    return 0;

  /*
   * Build the tree (can't use stree_ukkonen_build because of copyflag).
   */
  mprintf("\nBuilding the suffix tree...\n");
  if ((tree = stree_new_tree(string->alpha_size, 0, build_policy,
                             build_threshold)) == NULL)
    return 0;

  if (stree_ukkonen_add_string(tree, string->sequence, string->raw_seq,
                               string->length, 1) <= 0) {
    stree_delete_tree(tree);
    return 0;
  }

  /*
   * Preprocess the suffix tree and initialize repeat_struct
   */
  mprintf("Preprocessing...\n");
  tandem_struct = tandem_prep(tree, string->sequence, string->raw_seq,
                              string->length);
  if(tandem_struct == NULL) {
    stree_delete_tree(tree);
    return 0;
  }

  /*
   * Look up repeats
   */
  mprintf("\nThe following tandem arrays/tandem repeats were found:\n\n");
  tandem_lookup(tandem_struct);
  mend(14);

  /* Write summary of results */
  mstart(stdin, stdout, OK, OK, 0, NULL);
  mprintf("\nSummary:\n");
  mprintf("   Branching Primitive Tandem Repeats:           %d\n",
          tandem_struct->num_branching_primitive_tandem_repeats);
  mprintf("   Non-branching Primitive Tandem Repeats:       %d\n",
          tandem_struct->num_non_branching_primitive_tandem_repeats);
  mprintf("   Right-maximal Primitive Tandem Arrays (k>2):  %d\n",
          tandem_struct->num_right_maximal_primitive_tandem_arrays);
  mprintf("   Branching Non-primitive Tandem Repeats:       %d\n",
          tandem_struct->num_branching_non_primitive_tandem_repeats);
  mprintf("   Non-branching Non-primitive Tandem Repeats:   %d\n",
          tandem_struct->num_non_branching_non_primitive_tandem_repeats);

  /*
   * Write statistics and free memory
   */
  if (print_stats) {
    mprintf("\nStatistics:\n");
#ifdef STATS
    mprintf("   String Length:          %d\n", tandem_struct->length);
    mprintf("   Number of Tree Nodes:   %d\n", stree_get_num_nodes(tree));
    mprintf("   Preprocessing Steps:    %d\n", tandem_struct->num_prep);
    mprintf("   Number of Compares:     %d\n", tandem_struct->num_compares);
#else
    mprintf("   No statistics available.\n");
#endif
    mprintf("\n");
  }

  tandem_free(tandem_struct);
  stree_delete_tree(tree);

  return 1;
}

/*
 * strmat_repeats_vocabulary
 *
 * Find the vocabulary of tandem repeats in a string (plus several statistics).
 *
 * Parameters:   string           -  the input string
 *               build_policy     -  suffix tree build policy
 *               build_threshold  -  threshold used by LIST_THEN_ARRAY
 *               print_stats      -  flag telling whether to print the stats
 *               mode             -  flag telling whether
 *                                     A) most efficient decomposition
 *                                     B) non-overlapping Lempel-Ziv
 *
 * Returns:  non-zero on success, zero on error
 */
int strmat_repeats_vocabulary(STRING *string, int build_policy,
                              int build_threshold, int print_stats,
                              char mode)
{
  SUFFIX_TREE tree;
  DECOMPOSITION_STRUCT *decomp_struct;
  VOCABULARY_STRUCT *vocabulary_struct;

  if (string == NULL || string->sequence == NULL || string->length == 0)
    return 0;

  /*
   * Build the tree (can't use stree_ukkonen_build because of copyflag).
   */
  mprintf("\nBuilding the suffix tree...\n");
  if ((tree = stree_new_tree(string->alpha_size, 0, build_policy,
                             build_threshold)) == NULL)
    return 0;

  if (stree_ukkonen_add_string(tree, string->sequence, string->raw_seq,
                               string->length, 1) <= 0) {
    stree_delete_tree(tree);
    return 0;
  }

  /*
   * Compute Lempel-Ziv decomposition.
   */
  mprintf("\nComputing Lempel-Ziv decomposition...\n");
  decomp_struct = decomposition_prep(tree, string->sequence,
                                     string->raw_seq, string->length);
  if(decomp_struct == NULL) {
    stree_delete_tree(tree);
    return 0;
  }
  if(mode == 'A')
    lempel_ziv(decomp_struct);
  else if(mode == 'B')
    lempel_ziv_nonoverlapping(decomp_struct);
  else {
    fprintf(stderr,"Unknown mode `%c'.\n",mode);
    return 0;
  }

  /*
   * Allocate the memory.
   */
  mprintf("Preprocessing...\n");
  vocabulary_struct = vocabulary_prep(tree, decomp_struct, string->sequence,
                                      string->raw_seq, string->length);
  if(vocabulary_struct == NULL) {
    stree_delete_tree(tree);
    decomposition_free(decomp_struct);
    return 0;
  }

  /*
   * Find vocabulary of tandem repeats.
   */
  mprintf("\nLocating the Tandem Repeats...\n");
  vocabulary_find_tandem_repeats(vocabulary_struct);
  vocabulary_count(vocabulary_struct,
                   &vocabulary_struct->num_tandem_repeats,
                   &vocabulary_struct->num_tandem_repeat_occs);
  mprintf("The following tandem repeats were found:\n\n");
  vocabulary_write(vocabulary_struct,"tandem repeat");
  mend(23);

  /*
   * Find vocabulary of primitive tandem repeats.
   */
  mstart(stdin, stdout, OK, OK, 0, NULL);
  mprintf("\nLocating the Primitive Tandem Repeats...\n");
  vocabulary_find_primitive_tandem_repeats(vocabulary_struct);
  vocabulary_count(vocabulary_struct,
                   &vocabulary_struct->num_primitive_tandem_repeats,
                   &vocabulary_struct->num_primitive_tandem_repeat_occs);
  mprintf("The following primitive tandem repeats were found:\n\n");
  vocabulary_write(vocabulary_struct,"primitive tandem repeat");
  mend(23);

  /*
   * Find vocabulary of tandem arrays.
   */
  mstart(stdin, stdout, OK, OK, 0, NULL);
  mprintf("\nLocating the Tandem Arrays...\n");
  vocabulary_find_tandem_arrays(vocabulary_struct);
  vocabulary_count(vocabulary_struct,
                   &vocabulary_struct->num_tandem_arrays,
                   &vocabulary_struct->num_tandem_array_occs);
  mprintf("The following tandem arrays were found:\n\n");
  vocabulary_write(vocabulary_struct,"tandem array");
  mend(17);

  /*
   * Write summary of results.
   */
  mstart(stdin, stdout, OK, OK, 0, NULL);
  mprintf("\nSummary:\n");
  mprintf("   Tandem Repeats:                       %d\n",
          vocabulary_struct->num_tandem_repeats);
  mprintf("   Tandem Repeat Occurrences:            %d\n",
          vocabulary_struct->num_tandem_repeat_occs);
  mprintf("   Primitive Tandem Repeats:             %d\n",
          vocabulary_struct->num_primitive_tandem_repeats);
  mprintf("   Primitive Tandem Repeat Occurrences:  %d\n",
          vocabulary_struct->num_primitive_tandem_repeat_occs);
  mprintf("   Tandem Arrays:                        %d\n",
          vocabulary_struct->num_tandem_arrays);
  mprintf("   Tandem Array Occurrences:             %d\n",
          vocabulary_struct->num_tandem_array_occs);

  /*
   * Write statistics and free memory
   */
  if (print_stats) {
    mprintf("\nStatistics:\n");
#ifdef STATS
    mprintf("   String Length:                   %d\n", string->length);
    mprintf("   Number of Tree Nodes:            %d\n",
            stree_get_num_nodes(tree));
    mprintf("   Preprocessing Steps:             %d\n",
            vocabulary_struct->num_prep);
    mprintf("   Number of Compares\n");
    mprintf("     for Tandem Repeats:            %d\n",
            vocabulary_struct->num_compares_for_tandem_repeats);
    mprintf("     for Primitive Tandem Repeats:  %d\n",
            vocabulary_struct->num_compares_for_primitive_tandem_repeats);
    mprintf("     for Tandem Arrays:             %d\n",
            vocabulary_struct->num_compares_for_tandem_arrays);
#else
    mprintf("   No statistics available.\n");
#endif
    mprintf("\n");
  }

  vocabulary_free(vocabulary_struct);
  decomposition_free(decomp_struct);
  stree_delete_tree(tree);

  return 1;
}



/*
 * strmat_repeats_linear_occs
 *
 * Find all occurrences of tandem repeats in a string in time O(n + z).
 *
 * Parameters:   string           -  the input string
 *               build_policy     -  suffix tree build policy
 *               build_threshold  -  threshold used by LIST_THEN_ARRAY
 *               print_stats      -  flag telling whether to print the stats
 *               mode             -  flag telling whether
 *                                     A) most efficient decomposition
 *                                     B) non-overlapping Lempel-Ziv
 *
 * Returns:  non-zero on success, zero on error
 */
int strmat_repeats_linear_occs(STRING *string, int build_policy,
                               int build_threshold, int print_stats,
                               char mode)
{
  SUFFIX_TREE tree;
  DECOMPOSITION_STRUCT *decomp_struct;
  VOCABULARY_STRUCT *vocabulary_struct;

  if (string == NULL || string->sequence == NULL || string->length == 0)
    return 0;

  /*
   * Build the tree (can't use stree_ukkonen_build because of copyflag).
   */
  mprintf("\nBuilding the suffix tree...\n");
  if ((tree = stree_new_tree(string->alpha_size, 0, build_policy,
                             build_threshold)) == NULL)
    return 0;

  if (stree_ukkonen_add_string(tree, string->sequence, string->raw_seq,
                               string->length, 1) <= 0) {
    stree_delete_tree(tree);
    return 0;
  }

  /*
   * Compute Lempel-Ziv decomposition.
   */
  mprintf("Computing Lempel-Ziv decomposition...\n");
  decomp_struct = decomposition_prep(tree, string->sequence,
                                     string->raw_seq, string->length);
  if(decomp_struct == NULL) {
    stree_delete_tree(tree);
    return 0;
  }
  if(mode == 'A')
    lempel_ziv(decomp_struct);
  else if(mode == 'B')
    lempel_ziv_nonoverlapping(decomp_struct);
  else {
    fprintf(stderr,"Unknown mode `%c'.\n",mode);
    return 0;
  }


  /*
   * Allocate the memory.
   */
  mprintf("Preprocessing...\n");
  vocabulary_struct = linear_occs_prep(tree, decomp_struct, string->sequence,
                                       string->raw_seq, string->length);
  if(vocabulary_struct == NULL) {
    stree_delete_tree(tree);
    decomposition_free(decomp_struct);
    return 0;
  }

  /*
   * Find the occurrences of tandem repeats.
   */
  mprintf("Locating the Tandem Repeats...\n");
  linear_occs_find_tandem_repeats(vocabulary_struct);

  /*
   * Report the results.
   */
  mprintf("\nThe following tandem repeats were found:\n\n");
  linear_occs_write(vocabulary_struct,"tandem repeat");
  mend(17);

  /*
   * Write summary of results.
   */
  mstart(stdin, stdout, OK, OK, 0, NULL);
  mprintf("\nSummary:\n");
  mprintf("   Tandem Repeat Occurrences: %d\n",
          vocabulary_struct->num_tandem_repeat_occs);

  /*
   * Write statistics and free memory
   */
  if (print_stats) {
    mprintf("\nStatistics:\n");
#ifdef STATS
    mprintf("   String Length:             %d\n", string->length);
    mprintf("   Number of Tree Nodes:      %d\n", stree_get_num_nodes(tree));
    mprintf("   Preprocessing Steps:       %d\n", vocabulary_struct->num_prep);
    mprintf("   Number of Compares:        %d\n",
            vocabulary_struct->num_compares_for_tandem_repeats);
#else
    mprintf("   No statistics available.\n");
#endif
    mprintf("\n");
  }

  linear_occs_free(vocabulary_struct);
  decomposition_free(decomp_struct);
  stree_delete_tree(tree);

  return 1;
}

