/*
 * stree_decomposition.c
 *
 * Find the Ziv-Lempel decomposition and variants of it in O(n) time
 * using the (ready-built) suffix tree.
 *
 * NOTES:
 *    6/98  -  Original implementation of the algorithms (Jens Stoye)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#ifdef STRMAT
#include "stree_strmat.h"
#else
#include "stree.h"
#endif
#include "more.h"
#include "stree_decomposition.h"

/*=============================================================================
 *
 * decomposition_prep
 *
 * Memory allocation and initialization of a decomposition structure.
 *
 * Parameters:  tree        -  a suffix tree
 *              string      -  the string
 *              raw_string  -  the raw string
 *              length      -  the length of the string
 *
 * Returns:  A decomposition structure
 *
 *---------------------------------------------------------------------------*/
DECOMPOSITION_STRUCT *decomposition_prep(SUFFIX_TREE tree, char *string,
                                         char *raw_string, int length)
{
  DECOMPOSITION_STRUCT *d;

  if(tree == NULL)
    return NULL;

  /*
   * Allocate the memory.
   */
  if((d = malloc(sizeof(DECOMPOSITION_STRUCT))) == NULL)
    return NULL;
  memset(d, 0, sizeof(DECOMPOSITION_STRUCT));

  d->tree = tree;
  d->string = string;
  d->raw_string = raw_string;
  d->length = length;

  /* maxmally length blocks (plus end of last block) */
  if((d->prev = malloc((length+1) * sizeof(int))) == NULL) {
    decomposition_free(d);
    return NULL;
  }
  memset(d->prev, 0, (length+1) * sizeof(int));

  if((d->block = malloc((length+1) * sizeof(int))) == NULL) {
    decomposition_free(d);
    return NULL;
  }
  memset(d->block, 0, (length+1) * sizeof(int));

  d->num_blocks = 0;
  d->max_block_length = 0;

#ifdef STATS
  d->num_compares = 0;
  d->num_edge_traversals = 0;
#endif

  return d;
} /* decomposition_prep() */

/*=============================================================================
 *
 * decomposition_free
 *
 * Free memory of decomposition structure.
 *
 * Parameters:  d  -  a decomposition structure
 *
 * Returns:  nothing
 *
 *---------------------------------------------------------------------------*/
void decomposition_free(DECOMPOSITION_STRUCT *d)
{
  if(d->prev!= NULL)
    free(d->prev);
  if(d->block!= NULL)
    free(d->block);
  free(d);

} /* decomposition_free() */

/*=============================================================================
 *
 * decomposition_print
 *
 * Print decomposition.
 *
 * Parameters:  d  -  a decomposition structure
 *
 * Returns:  nothing
 *
 *---------------------------------------------------------------------------*/
void decomposition_print(DECOMPOSITION_STRUCT *d)
{
  int i;

  for(i=0; i<d->num_blocks; i++)
    mprintf("Block %d: start=%d, length=%d (leftmost occurrence: %d)\n",
            i+1,d->block[i]+1,d->block[i+1]-d->block[i],d->prev[i]+1);
  mprintf("End of last block: %d\n",d->block[i]);

} /* decomposition_print() */

/*=============================================================================
 *
 * lempel_ziv
 *
 * Compute the Lempel-Ziv decomposition (f-factorization).
 *
 * CAVEAT: THIS WORKS ONLY IF THE LABEL ATTACHED AT EACH EDGE OF THE
 *         SUFFIX TREE POINTS TO THE LEFTMOST POSSIBLE SUBSTRING, AS FOR
 *         EXAMPLE WHEN THE TREE IS COMPUTED WITH UKKONEN'S ALGORITHM.
 *
 * Parameters:  d  -  a decomposition structure
 *
 * Returns:  nothing
 * 
 *---------------------------------------------------------------------------*/
void lempel_ziv(DECOMPOSITION_STRUCT *d)
{
  int i, j, block_length;
  STREE_NODE node, child;

  d->prev[0] = -1;
  d->block[0] = 0;
  d->block[1] = 1;
  d->num_blocks = 1;
  d->max_block_length = 1;

  for(i=j=1; i<d->length; i=j) {

    node = stree_get_root(d->tree);
    child = stree_find_child(d->tree,node,d->string[j]);

    while(child != NULL && stree_get_edgestr(d->tree,child) < &d->string[j]) {
      j += stree_get_edgelen(d->tree,child);
      node = child;
      if(j < d->length)
        child = stree_find_child(d->tree,node,d->string[j]);
      else
        child = NULL;
#ifdef STATS
      d->num_edge_traversals++;
      d->num_compares++;
#endif
    }

    d->num_blocks++;
    if(node == stree_get_root(d->tree)) {
      d->prev[d->num_blocks-1] = -1;
      d->block[d->num_blocks] = ++j;
    }
    else {
      d->prev[d->num_blocks-1] = stree_get_edgestr(d->tree,node) - d->string
                                 + stree_get_edgelen(d->tree,node) - j + i;
      d->block[d->num_blocks] = j;
    }

#ifdef STATS
    d->num_compares += 2;
#endif

    if((block_length = d->block[d->num_blocks]-d->block[d->num_blocks-1]) >
       d->max_block_length)
      d->max_block_length = block_length;

  } /* for i */

} /* lempel_ziv() */

/*=============================================================================
 *
 * lempel_ziv_nonoverlapping
 *
 * Compute non-overlapping Lempel-Ziv decomposition (as described in the Book,
 * Section 7.17).
 *
 * CAVEAT: THIS WORKS ONLY IF THE LABEL ATTACHED AT EACH EDGE OF THE
 *         SUFFIX TREE POINTS TO THE LEFTMOST POSSIBLE SUBSTRING, AS FOR
 *         EXAMPLE WHEN THE TREE IS COMPUTED WITH UKKONEN'S ALGORITHM.
 *
 * Parameters:  d  -  a decomposition structure
 *
 * Returns:  nothing
 *
 *---------------------------------------------------------------------------*/
void lempel_ziv_nonoverlapping(DECOMPOSITION_STRUCT *d)
{
  int i, j, block_length, edgelen;
  char *edgestr;
  STREE_NODE node, child;

  d->prev[0] = -1;
  d->block[0] = 0;
  d->block[1] = 1;
  d->num_blocks = 1;
  d->max_block_length = 1;

  for(i=j=1; i<d->length; i=j) {

    node = stree_get_root(d->tree);
    child = stree_find_child(d->tree,node,d->string[j]);

    while(child != NULL &&
          (edgestr=stree_get_edgestr(d->tree,child)) + 
          (edgelen=stree_get_edgelen(d->tree,child)) < &d->string[i]) {
      j += edgelen;
      node = child;
      if(j < d->length)
        child = stree_find_child(d->tree,node,d->string[j]);
      else
        child = NULL;
#ifdef STATS
      d->num_edge_traversals++;
      d->num_compares++;
#endif
    }

    d->num_blocks++;
    if(node == stree_get_root(d->tree) && edgestr+edgelen != &d->string[i]) {
      d->prev[d->num_blocks-1] = -1;
      d->block[d->num_blocks] = ++j;
    }
    else {
      if(child != NULL && edgestr < &d->string[i]) {
        d->prev[d->num_blocks-1] = edgestr - d->string - j + i;
        j += &d->string[i] - edgestr;
      }
      else
        d->prev[d->num_blocks-1] = stree_get_edgestr(d->tree,node) - d->string
                                   + stree_get_edgelen(d->tree,node) - j + i;
      d->block[d->num_blocks] = j;
    }

#ifdef STATS
    d->num_compares += 2;
#endif

    if((block_length = d->block[d->num_blocks]-d->block[d->num_blocks-1]) >
       d->max_block_length)
      d->max_block_length = block_length;

  } /* for i */

} /* lempel_ziv_nonoverlapping() */

/****** EOF (stree_decomposition.c) ******************************************/

