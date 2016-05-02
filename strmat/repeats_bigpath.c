/*****************************************************************************/
/*
 * repeats_bigpath.c
 *
 * Implementation of the nonoverlapping maximal pairs algorithm
 * which constructs the leaf-lists of all nodes of a suffix tree in
 * overall O(n log n) time. The nonoverlapping maximal pairs are
 * reported in O(n log n + z) time.
 *
 * NOTES:
 *    8/98  -  Original implementation of the algorithms (Jens Stoye)
 *    2/99  -  Correction of some minor bugs (Jens Stoye)
 */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#ifdef STRMAT
#include "stree_strmat.h"
#else
#include "stree.h"
#endif
#include "more.h"
#include "repeats_bigpath.h"

/*===========================================================================*/
/*
 * bp_append_entry
 *
 * Append entry e to list number c of node with node_id: O(1) time.
 *
 * Parameters:  b        -  a bigpath structure
 *              e        -  an entry
 *              node_id  -  index of the node
 *              c        -  character to the left
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void bp_append_entry(bp_struct *b, bp_entry *e, int node_id, int c)
{
  e->prev = b->last[node_id][c];
  e->next = NULL;
  e->mark = -1;
  if(b->last[node_id][c] != NULL)
    b->last[node_id][c]->next = e;
  if(b->list[node_id][c] == NULL)
    b->list[node_id][c] = e;
  b->last[node_id][c] = e;

#ifdef STATS
  b->num_steps_for_lists++;
#endif

  return;
} /* bp_append_entry() */

/*===========================================================================*/
/*
 * bp_copy_entry
 *
 * If e->mark defined, replace entry e by a copy of itself (e2) and append e
 * to the list with index e->mark; otherwise, only replace e by a copy of
 * itself (e2): O(1) time.
 *
 * Parameters:  b        -  a bigpath structure
 *              e        -  an entry
 *              node_id  -  index of the node
 *              c        -  character to the left
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void bp_copy_entry(bp_struct *b, bp_entry *e, int node_id, int c)
{
  bp_entry *e2;

  e2 = &b->entries2[e-b->entries];

  e2->prev = e->prev;
  if(e->prev != NULL)
    e->prev->next = e2;
  else
    b->list[node_id][c] = e2;
  e2->next = e->next;
  if(e->next != NULL)
    e->next->prev = e2;
  else
    b->last[node_id][c] = e2;

  if(e->mark != -1)
    bp_append_entry(b,e,e->mark,c);

#ifdef STATS
  b->num_steps_for_lists++;
#endif

  return;
} /* bp_copy_entry() */

/*===========================================================================*/
/*
 * bp_remove_entry
 *
 * Save previous left neighbor of entry[pos] and then remove entry2[pos] from
 * its list: O(1) time.
 *
 * Parameters:  b        -  a bigpath structure
 *              pos      -  index of an entry
 *              node_id  -  index of the node
 *              c        -  character to the left
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void bp_remove_entry(bp_struct *b, int pos, int node_id, int c)
{
  bp_entry *e, *e2;

  e = &b->entries[pos];
  e2 = &b->entries2[pos];
  e->old_prev = e2->prev;

  if(e2->prev != NULL)
    e2->prev->next = e2->next;
  else
    b->list[node_id][c] = e2->next;
  if(e2->next != NULL)
    e2->next->prev = e2->prev;
  else
    b->last[node_id][c] = e2->prev;

#ifdef STATS
  b->num_steps_for_lists++;
#endif

  return;
} /* bp_remove_entry() */

/*===========================================================================*/
/*
 * bp_find_big
 *
 * Recursive procedure which finds the big child of each node: total O(n) time.
 *
 * Parameters:  b     -  a bigpath structure
 *              node  -  a suffix tree node
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
int bp_find_big(bp_struct *b, STREE_NODE node)
{
  int big_num, num, num_leaves;
  STREE_NODE big_child, child;

  big_num = -1;
  big_child = NULL;
  num = 0;

  /* depth-first */
  for(child=stree_get_children(b->tree,node);
      child!=NULL;
      child=stree_get_next(b->tree,child)) {
    num_leaves = bp_find_big(b,child);
    if(num_leaves > big_num) {
      big_num = num_leaves;
      big_child = child;
    }
    num += num_leaves;
  }

  b->big_child[stree_get_ident(b->tree,node)] = big_child;

#ifdef STATS
  b->num_prep++;
#endif

  return num + stree_get_num_leaves(b->tree,node);

} /* bp_find_big() */

/*===========================================================================*/
/*
 * bigpath_prep
 *
 * Preprocessing for the bigpath array algorithm: O(1) time.
 *
 * Parameters:  tree        -  a suffix tree
 *              string      -  the string
 *              raw_string  -  the raw string
 *              length      -  the length of the string
 *
 * Returns:  a bigpath structure
 */
/*---------------------------------------------------------------------------*/
bp_struct *bigpath_prep(SUFFIX_TREE tree,
                        char *string, char *raw_string, int length)
{
  bp_struct *b;
  unsigned int num_nodes, num_leaves;
  STREE_NODE root;
  int i, c, alph, root_id;

  if (tree == NULL)
    return NULL;

  /*
   * Allocate the memory.
   */
  if ((b = malloc(sizeof(bp_struct))) == NULL)
    return NULL;
  memset(b, 0, sizeof(bp_struct));

  b->string = string;
  b->raw_string = raw_string;
  b->length = length;

  b->tree = tree;

  num_nodes = stree_get_num_nodes(tree);
  num_leaves = length;

  if ((b->entries = malloc(num_leaves * sizeof(bp_entry))) == NULL) {
    bigpath_free(b);
    return NULL;
  }
  memset(b->entries, 0, num_leaves * sizeof(bp_entry));

  if ((b->entries2 = malloc(num_leaves * sizeof(bp_entry))) == NULL) {
    bigpath_free(b);
    return NULL;
  }
  memset(b->entries2, 0, num_leaves * sizeof(bp_entry));

  if ((b->list = malloc(num_nodes * sizeof(bp_entry**))) == NULL) {
    bigpath_free(b);
    return NULL;
  }
  memset(b->list, 0, num_nodes * sizeof(bp_entry**));

  if ((b->last = malloc(num_nodes * sizeof(bp_entry**))) == NULL) {
    bigpath_free(b);
    return NULL;
  }
  memset(b->last, 0, num_nodes * sizeof(bp_entry**));

  if ((b->big_child = malloc(num_nodes * sizeof(STREE_NODE)))== NULL) {
    bigpath_free(b);
    return NULL;
  }
  memset(b->big_child, 0, num_nodes * sizeof(STREE_NODE));

  b->num_nonoverlapping_maximal_pairs = 0;

#ifdef STATS
  b->num_prep = 0;
  b->num_steps_for_lists = 0;
  b->num_compares = 0;
#endif

  /*
   * Find alphabet size and a pointers
   */
  for(i=0; i<CHAR_MAX; i++)
    b->a[i] = -1;
  alph = 0;
  for(i=0; i<length; i++) { /* including delimiters */
    c = (int)(unsigned char)string[i];
    if(b->a[c] == -1)
      b->a[c] = alph++;
  }
  b->alpha_size = alph;

  /*
   * Allocate memory of the lists for the individual characters.
   */
  for(i=0; i<num_nodes; i++) {
    if ((b->list[i] = malloc((alph+1) * sizeof(bp_entry*))) == NULL) {
      bigpath_free(b);
      return NULL;
    }
    memset(b->list[i], 0, (alph+1) * sizeof(bp_entry*));

    if ((b->last[i] = malloc((alph+1) * sizeof(bp_entry*))) == NULL) {
      bigpath_free(b);
      return NULL;
    }
    memset(b->last[i], 0, (alph+1) * sizeof(bp_entry*));
  }

  /*
   * Create basic lists.
   */
  root = stree_get_root(b->tree);
  root_id = stree_get_ident(b->tree,root);

  bp_append_entry(b,&b->entries[0],root_id,alph);
  for(i=1; i<b->length ; i++) {
    c = b->a[(int)(unsigned char)b->string[i-1]];
    bp_append_entry(b,&b->entries[i],root_id,c);
  }

  /*
   * Find the big children.
   */
  root = stree_get_root(b->tree);
  bp_find_big(b,root);

  return b;

} /* bigpath_prep() */

/*===========================================================================*/
/*
 * bigpath_free
 *
 * Free memory of a bigpath structure.
 *
 * Parameters:  b  -  a bigpath structure
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void bigpath_free(bp_struct *b)
{
  int i, num_nodes;

  num_nodes = stree_get_num_nodes(b->tree);
  for(i=0; i<num_nodes; i++) {
    if(b->list[i] != NULL)
      free(b->list[i]);
    if(b->last[i] != NULL);
      free(b->last[i]);
  }
  if(b->entries != NULL)
    free(b->entries);
  if(b->entries2 != NULL)
    free(b->entries2);
  if(b->list != NULL)
    free(b->list);
  if(b->last != NULL)
    free(b->last);
  if(b->big_child != NULL)
    free(b->big_child);
  free(b);
} /* bigpath_free() */

/*===========================================================================*/
/*
 * bp_mark
 *
 * Recursive procedure marking the entries of a subtree with root_id: total
 * O(size of subtree) time.
 *
 * Parameters:  b        -  a bigpath structure
 *              node     -  a suffix tree node
 *              root_id  -  id of root of the subtree (the small child)
 */
/*---------------------------------------------------------------------------*/
void bp_mark(bp_struct *b, STREE_NODE node, int root_id)
{
  int leavesnum, i, pos, dummy_id;
  STREE_NODE child;
  char *dummy_string;

  /* recurse depth-first */
  for(child=stree_get_children(b->tree,node);
      child!=NULL;
      child=stree_get_next(b->tree,child))
    bp_mark(b,child,root_id);

  /* leaves: set marks */
  leavesnum = stree_get_num_leaves(b->tree,node);
  for(i=1; i<=leavesnum; i++) {
    stree_get_leaf(b->tree,node,i,&dummy_string,&pos,&dummy_id);
    b->entries[pos].mark = root_id;
  }

#ifdef STATS
  b->num_steps_for_lists++;
#endif

  return;
} /* bp_mark() */

/*===========================================================================*/
/*
 * bp_write
 *
 * Write nonoverlapping maximal pair: O(1) time.
 *
 * Parameters:  type        -  repeat type
 *              pos1        -  starting position of first occurrence
 *              pos2        -  starting position of second occurrence
 *              len         -  length of the repeat
 *              raw_string  -  the raw string
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void bp_write(char *type, int pos1, int pos2, int len, char *raw_string)
{
  int i, textlen, restlen;
  char *s, *t, buffer[77];

  sprintf(buffer,"%s (%d,%d,%d): ",type,pos1+1,pos2+1,len);
  buffer[76] = '\0';

  textlen = strlen(buffer);
  restlen = 76-textlen;
  for(i=0,s=&buffer[textlen],t=&raw_string[pos1];
      i<restlen && i<len;
      i++,s++,t++)
    *s = isprint((int)(*t)) ? *t : '#';
  *s = '\0';
  mputs(buffer);
  if(len > restlen)
    mputs("...");
  mputc('\n');
  return;
} /* bp_write() */

/*===========================================================================*/
/*
 * bp_report_entry
 *
 * Report nonoverlapping maximal pair of entry versus list: O(z) time.
 *
 * Parameters:  b       -  a bigpath structure
 *              id      -  index of a leaf list
 *              i       -  index of an entry of another leaf list
 *              c       -  character to the left of entry i
 *              d       -  string-depth of list id
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void bp_report_entry(bp_struct *b, int id, int i, int c, int d)
{
  int cc;
  bp_entry *e;

  if(d > 0) {

    /* report nonoverlapping repeats */
    for(cc=0; cc<=b->alpha_size; cc++)
      if(cc != c) {
        for(e=b->list[id][cc]; e!=NULL && e-b->entries2+d<=i; e=e->next) {
          bp_write("nonoverlapping maximal pair",e-b->entries2,i,d,
                        b->raw_string);
          b->num_nonoverlapping_maximal_pairs++;
#ifdef STATS
          b->num_compares++;
#endif
        }
        for(e=b->last[id][cc]; e!=NULL && i+d<=e-b->entries2; e=e->prev) {
          bp_write("nonoverlapping maximal pair",i,e-b->entries2,d,
                   b->raw_string);
          b->num_nonoverlapping_maximal_pairs++;
#ifdef STATS
          b->num_compares++;
#endif
        }
      }
#ifdef STATS
    b->num_compares += 2;
#endif

  }
  return;
} /* bp_report_entry() */

/*===========================================================================*/
/*
 * bp_report_list
 *
 * Report nonoverlapping maximal pairs of list versus list: O(z) time.
 *
 * Parameters:  b        -  a bigpath structure
 *              id1      -  index of first leaf list
 *              id2      -  index of second leaf list
 *              d        -  string-depth of list1
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void bp_report_list(bp_struct *b, int id1, int id2, int d)
{
  int c;
  bp_entry *e;

  for(c=0; c<=b->alpha_size; c++)
    for(e=b->list[id2][c]; e!=NULL; e=e->next)
      bp_report_entry(b,id1,e-b->entries,c,d);

  return;
} /* bp_report_list() */

/*===========================================================================*/
/*
 * bp_find_rec
 *
 * Find all nonoverlapping and max-gapped maximal pairs: total O(n log n + z)
 * time.
 *
 * Parameters:  b        -  a bigpath structure
 *              root     -  a suffix tree node
 *              depth    -  string-depth of root
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void bp_find_rec(bp_struct *b, STREE_NODE root, int depth)
{
  int c,d, i, root_id,node_id,child_id,big_child_id, leavesnum, pos, dummy_id;
  STREE_NODE node, child, big_child;
  bp_entry *e, *e_next, **tmp;
  char *dummy_string;

  /* marking phase (along big path) */
  node = root;
  while(node != NULL) {
    node_id = stree_get_ident(b->tree,node);
    big_child = b->big_child[node_id];

    for(child=stree_get_children(b->tree,node);
        child!=NULL;
        child=stree_get_next(b->tree,child))
      if(child != big_child) {
        child_id = stree_get_ident(b->tree,child);
        bp_mark(b,child,child_id);
      }
    /* direct leaves are not marked or copied */
    node = big_child;
  }

  /* copying phase (along list); unmarked entries ar not copied */
  root_id = stree_get_ident(b->tree,root);
  for(c=0; c<=b->alpha_size; c++)
    for(e=b->list[root_id][c]; e!=NULL; e=e_next) {
      e_next = e->next;
      bp_copy_entry(b,e,root_id,c);
    }

  /* pruning phase (along big path, incl. reporting) */
  node = root;
  d = depth;
  while(node != NULL) {
    node_id = stree_get_ident(b->tree,node);
    big_child = b->big_child[node_id];

    for(child=stree_get_children(b->tree,node);
        child!=NULL;
        child=stree_get_next(b->tree,child))
      if(child != big_child) {
        child_id = stree_get_ident(b->tree,child);
        for(c=0; c<=b->alpha_size; c++)
          for(e=b->list[child_id][c]; e!=NULL; e=e->next)
            bp_remove_entry(b,e-b->entries,node_id,c);
        bp_report_list(b,node_id,child_id,d);
      }
    leavesnum = stree_get_num_leaves(b->tree,node);
    for(i=1; i<=leavesnum; i++) {
      stree_get_leaf(b->tree,node,i,&dummy_string,&pos,&dummy_id);
      c = pos==0 ? b->alpha_size : b->a[(int)(unsigned char)b->string[pos-1]];
      bp_remove_entry(b,pos,node_id,c);
      bp_report_entry(b,node_id,pos,c,d);
    }

    if(big_child != NULL) {
      big_child_id = stree_get_ident(b->tree,big_child);
      tmp = b->list[big_child_id];              /* swap not really necessary */
      b->list[big_child_id] = b->list[node_id]; /* but avoids memory leakage */
      b->list[node_id] = tmp;
      tmp = b->last[big_child_id];
      b->last[big_child_id] = b->last[node_id];
      b->last[node_id] = tmp;
      d += stree_get_edgelen(b->tree,big_child);
    }
    node = big_child;
  }

  /* recursive calls (along big path) */
  node = root;
  d = depth;
  while(node != NULL) {
    node_id = stree_get_ident(b->tree,node);
    big_child = b->big_child[node_id];

    for(child=stree_get_children(b->tree,node);
        child!=NULL;
        child=stree_get_next(b->tree,child))
      if(child != big_child)
        bp_find_rec(b,child,d+stree_get_edgelen(b->tree,child));

    if(big_child != NULL)
      d += stree_get_edgelen(b->tree,big_child);
    node = big_child;
  }

} /* bp_find_rec() */

/*===========================================================================*/
/*
 * bigpath_find
 *
 * Find all nonoverlapping maximal pairs.
 *
 * Parameters:  b        -  a bigpath structure
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void bigpath_find(bp_struct *b)
{
  bp_find_rec(b,stree_get_root(b->tree),0);
  return;
} /* bigpath_find() */

/****** EOF (repeats_bigpath.c) **********************************************/

