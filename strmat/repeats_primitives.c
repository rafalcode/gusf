/*
 * repeats_primitives.c
 *
 * This file contains an implementation of Crochemore's algorithm (1981) to 
 * find all occurrences of primitive tandem repeats in O(n log n) time.
 *
 * NOTES:
 *    7/98  -  Original implementation of the algorithms (Jens Stoye)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include "strmat.h"
#include "more.h"
#include "repeats_primitives.h"

/*
 *
 * Forward References.
 *
 */
void pr_remove_list(pr_list *l);
void pr_remove_node(pr_node *n);

/*===========================================================================*/
/*
 * pr_append_entry
 *
 * Append entry e to list l: O(1) time.
 *
 * Parameters:  e  -  an entry
 *              l  -  a list
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void pr_append_entry(pr_entry *e, pr_list *l)
{
  e->next = NULL;
  e->prev = l->last;
  e->inList = l;
  if(l->last != NULL)
    l->last->next = e;
  if(l->entries == NULL)
    l->entries = e;
  l->last = e;
  l->len++;
  return;
} /* pr_append_entry() */

/*===========================================================================*/
/*
 * pr_remove_entry
 *
 * Remove entry e from its list: O(1) time.
 *
 * Parameters:  e  -  an entry
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void pr_remove_entry(pr_entry *e)
{
  if(e->inList != NULL) { /* only if entry was not removed earlier */
    if(e->prev != NULL)
      e->prev->next = e->next;
    else
      e->inList->entries = e->next;
    if(e->next != NULL)
      e->next->prev = e->prev;
    else
      e->inList->last = e->prev;
    e->inList->len--;
    if(e->inList->len == 0)
      pr_remove_list(e->inList);
    e->inList = NULL;
  }
  return;
} /* pr_remove_entry() */

/*===========================================================================*/
/*
 * pr_move_entry
 *
 * Remove entry e from its list and append it to list l: O(1) time.
 *
 * Parameters:  e  -  an entry
 *              l  -  a list
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void pr_move_entry(pr_entry *e, pr_list *l)
{
  pr_remove_entry(e);
  pr_append_entry(e,l);
  return;
} /* pr_move_entry() */

/*===========================================================================*/
/*
 * pr_new_list
 *
 * Create new (empty) list: O(1) time.
 *
 * Parameters:  p  -  a primitives structure
 *
 * Returns:  the new list
 */
/*---------------------------------------------------------------------------*/
pr_list *pr_new_list(primitives_struct *p)
{
  pr_list *l_new;

  l_new = &p->lists[p->nextList++];
  l_new->next = l_new->prev = NULL;
  l_new->atNode = NULL;
  l_new->entries = l_new->last = NULL;
  l_new->len = 0;
  return l_new;
} /* pr_new_list() */

/*===========================================================================*/
/*
 * pr_append_list
 *
 * Append list l to the lists at node n: O(1) time.
 *
 * Parameters:  l  -  a list
 *              n  -  a node
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void pr_append_list(pr_list *l, pr_node *n)
{
  l->next = NULL;
  l->prev = n->last;
  l->atNode = n;
  if(n->last != NULL)
    n->last->next = l;
  if(n->lists == NULL)
    n->lists = l;
  n->last = l;
  return;
} /* pr_append_list() */

/*===========================================================================*/
/*
 * pr_remove_list
 *
 * Remove list l from its node: O(1) time.
 *
 * Parameters:  l  -  a list
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void pr_remove_list(pr_list *l)
{
  if(l->prev != NULL)
    l->prev->next = l->next;
  else
    l->atNode->lists = l->next;
  if(l->next != NULL)
    l->next->prev = l->prev;
  else
    l->atNode->last = l->prev;
  l->atNode = NULL;
  return;
} /* pr_remove_list() */

/*===========================================================================*/
/*
 * pr_replace_list
 *
 * Replace list l by a copy of itself (l_new): O(length of l) time.
 *
 * Parameters:  p  -  a primitives structure
 *              l  -  a list
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void pr_replace_list(primitives_struct *p, pr_list *l)
{
  pr_list *l_new;
  pr_entry *e, *e_new;

  /* create new list (manually in p->lists2) */
  l_new = &p->lists2[p->nextList2++];

  /* replace l by new list */
  l_new->next = l->next;
  l_new->prev = l->prev;
  l_new->atNode = l->atNode;
  l_new->entries = l_new->last = NULL;
  l_new->len = 0;
  if(l->prev != NULL)
    l->prev->next = l_new;
  else
    l->atNode->lists = l_new;
  if(l->next != NULL)
    l->next->prev = l_new;
  else
    l->atNode->last = l_new;

  /* remove l from its node */
  l->atNode = NULL;

  /* fill new list with copies of l's entries */
  for(e=l->entries; e!=NULL; e=e->next) {
    e_new = &p->entries2[e-p->entries];
    pr_append_entry(e_new,l_new);
  }

  return;
} /* pr_replace_list() */

/*===========================================================================*/
/*
 * pr_new_node
 *
 * Create new (empty) node: O(1) time.
 *
 * Parameters:  p  -  a primitives structure
 *
 * Returns:  the new node
 */
/*---------------------------------------------------------------------------*/
pr_node *pr_new_node(primitives_struct *p)
{
  pr_node *n_new;

  n_new = &p->nodelist[p->nextNode++];
  n_new->next = n_new->prev = NULL;
  n_new->inNodelist = NULL;
  n_new->lists = n_new->last = n_new->last_source_list = NULL;
  return n_new;
} /* pr_new_node() */

/*===========================================================================*/
/*
 * pr_append_node
 *
 * Append node n to nodes of p: O(1) time.
 *
 * Parameters:  p  -  a primitives structure
 *              n  -  a node
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void pr_append_node(primitives_struct *p, pr_node *n)
{
  n->next = NULL;
  n->prev = p->last;
  n->inNodelist = p;
  if(p->last != NULL)
    p->last->next = n;
  if(p->nodes == NULL)
    p->nodes = n;
  p->last = n;
  return;
} /* pr_append_node() */

/*===========================================================================*/
/*
 * pr_remove_node
 *
 * Remove node n from its nodelist: O(1) time.
 *
 * Parameters:  n  -  a node
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void pr_remove_node(pr_node *n)
{
  if(n->prev != NULL)
    n->prev->next = n->next;
  else
    n->inNodelist->nodes2 = n->next;
  if(n->next != NULL)
    n->next->prev = n->prev;
  else
    n->inNodelist->last2 = n->prev;
  n->inNodelist = NULL;
  return;
} /* pr_remove_node() */

/*===========================================================================*/
/*
 * primitives_prep
 *
 * Create new (empty) primitives structure p_new: O(1) time.
 *
 * Parameters:  string         -  the string
 *              raw_string     -  the raw string
 *              length         -  the length of the string
 *
 * Returns:  a primitives structure
 */
/*---------------------------------------------------------------------------*/
primitives_struct *primitives_prep(char *string, char *raw_string, int length)
{
  primitives_struct *p_new;

  if((p_new = malloc(sizeof(primitives_struct))) == NULL)
    return NULL;
  memset(p_new, 0, sizeof(primitives_struct));

  p_new->string = string;
  p_new->raw_string = raw_string;
  p_new->length = length;

  if((p_new->entries = malloc(length * sizeof(pr_entry))) == NULL) {
    primitives_free(p_new);
    return NULL;
  }
  memset(p_new->entries, 0, length * sizeof(pr_entry));

  if((p_new->entries2 = malloc(length * sizeof(pr_entry))) == NULL) {
    primitives_free(p_new);
    return NULL;
  }
  memset(p_new->entries2, 0, length * sizeof(pr_entry));

  if((p_new->lists = malloc(2*length * sizeof(pr_list))) == NULL) {
    primitives_free(p_new);
    return NULL;
  }
  memset(p_new->lists, 0, 2*length * sizeof(pr_list));

  if((p_new->lists2 = malloc(length * sizeof(pr_list))) == NULL) {
    primitives_free(p_new);
    return NULL;
  }
  memset(p_new->lists2, 0, length * sizeof(pr_list));

  p_new->nextList = p_new->nextList2 = 0;

  if((p_new->nodelist = malloc(length * sizeof(pr_node))) == NULL) {
    primitives_free(p_new);
    return NULL;
  }
  memset(p_new->nodelist, 0, length * sizeof(pr_node));

  if((p_new->nodelist2 = malloc(length * sizeof(pr_node))) == NULL) {
    primitives_free(p_new);
    return NULL;
  }
  memset(p_new->nodelist2, 0, length * sizeof(pr_node));

  p_new->nextNode = 0;
  p_new->nodes = p_new->nodes2 = p_new->last, p_new->last2 = NULL;

  p_new->num_primitive_tandem_repeat_occs = 0;

#ifdef STATS
  p_new->num_compares = 0;
#endif


  return p_new;

} /* primitives_prep() */

/*===========================================================================*/
/*
 * pr_next_level
 *
 * Step forward to next level (i.e. swap nodelists): O(1) time.
 *
 * Parameters:  p  -  a primitives structure
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void pr_next_level(primitives_struct *p)
{
  pr_node *n_tmp;

  /* entries and lists are not swapped: they are moved in the main procedure */

  n_tmp = p->nodelist;
  p->nodelist = p->nodelist2;
  p->nodelist2 = n_tmp;

  p->nodes2 = p->nodes;
  p->nodes = NULL;
  p->last2 = p->last;
  p->last = NULL;
  p->nextNode = p->nextList2 = 0;

  return;
} /* pr_next_level() */

/*===========================================================================*/
/*
 * primitives_free
 *
 * Delete primitives structure p: O(1) time.
 *
 * Parameters:  p  -  a primitives structure
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void primitives_free(primitives_struct *p)
{
  free(p->entries);
  free(p->entries2);
  free(p->lists);
  free(p->lists2);
  free(p->nodelist);
  free(p->nodelist2);
  free(p);
  return;
} /* primitives_free() */

/*===========================================================================*/
/*
 * pr_write
 *
 * Write a primitive tandem repeat.
 *
 * Parameters:  raw_string  -  the raw string
 *              pos         -  starting position of the repeat
 *              len         -  length of the repeat
 *              rep         -  number of repeats
 *              type        -  repeat type
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void pr_write(char *raw_string, int pos, int len, int rep, char *type)
{
  int i, textlen, restlen;
  char *s, *t, buffer[77];

  sprintf(buffer,"%s: (%d,%d,%d) ",type,pos+1,len,rep);
  buffer[76] = '\0';

  textlen = strlen(buffer);
  restlen = 76-textlen;
  for(i=0,s=&buffer[textlen],t=&raw_string[pos];
      i<restlen && i<len*rep;
      i++,s++,t++)
    *s = isprint((int)(*t)) ? *t : '#';
  *s = '\0';
  mputs(buffer);
  if(len > restlen)
    mputs("...");
  mputc('\n');

} /* pr_write() */

/*===========================================================================*/
/*
 * pr_report
 *
 * Report all primitive tandem repeats in this iteration.
 *
 * Parameters:  p          -  a primitives structure
 *              iteration  -  iteration
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void pr_report(primitives_struct *p, int iteration)
{
  pr_node *n;
  pr_list *l;
  pr_entry *e, *previous;

  for(n=p->nodes; n!=NULL; n=n->next)
    for(l=n->lists; l!=NULL; l=l->next)
      if(l->entries != l->last) { /* more than one entry in list l */
        previous = l->entries;
        for(e=l->entries->next; e!=NULL; e=e->next) {
          if(e-previous == iteration) {
            pr_write(p->raw_string,previous-p->entries,iteration,2,
                     "primitive tandem repeat");
            p->num_primitive_tandem_repeat_occs++;
          }
          previous = e;
        }
      }
  return;
} /* pr_report() */

/*===========================================================================*/
/*
 * pr_create_basic_lists
 *
 * Create basic lists: O(CHAR_MAX) space, O(length + CHAR_MAX) time.
 *
 * Parameters:  p  -  a primitives structure
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void pr_create_basic_lists(primitives_struct *p)
{
  int i,c;
  pr_list *occ[CHAR_MAX+1];
  pr_node *n;

  for(i=0; i<=CHAR_MAX; i++)
    occ[i] = NULL;

  n = pr_new_node(p);
  pr_append_node(p,n);
  for(i=0; i<p->length ; i++) {
    c = (int)p->string[i];
    if(occ[c] == NULL)
      pr_append_list((occ[c]=pr_new_list(p)),n);
    pr_append_entry(&p->entries[i],occ[c]);
  }
  return;
} /* pr_create_basic_lists() */

/*===========================================================================*/
/*
 * primitives_find
 *
 * This is Crochemore's O(n log n) time algorithm.
 *
 * Parameters:  p  -  a primitives structure
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void primitives_find(primitives_struct *p)
{
  pr_node *n,*n_next,*nn,*new_n;
  pr_list *l,*l_next,*maxlist;
  pr_entry *e,*ee;
  int i, maxlistlen, pos;

  /* create basic lists */
  pr_create_basic_lists(p);

  for(i=1; i<p->length && p->nodes!=NULL; i++) {

    /* report primitive repeats */
    pr_report(p,i);

    /* step forward to next level */
    pr_next_level(p);

    /* for each node of the previous level */
    for(n=p->nodes2; n!=NULL; n=n_next) {
      n_next = n->next;

      /* find big list */
      maxlistlen = 0;
      for(l=n->lists; l!=NULL; l=l->next)
        if(l->len > maxlistlen) {
          maxlistlen = l->len;
          maxlist = l;
        }

      /* copy small lists/move big list (except singletons) */
      for(l=n->lists; l!=NULL; l=l_next) {
        l_next = l->next;
        if(l == maxlist)
          pr_remove_list(l);
        else
          pr_replace_list(p,l);
        if(l->len == 1)
          l->entries->inList = NULL; /* mark entry removed */
        else {
          new_n = pr_new_node(p);
          pr_append_node(p,new_n);
          pr_append_list(l,new_n);
        }
      } /* for l */
    } /* for n */

    /* using the remaining (small) lists, pull entries out */
    for(n=p->nodes2; n!=NULL; n=n->next)
      for(l=n->lists; l!=NULL; l=l->next)
        for(e=l->entries; e!=NULL; e=e->next) {
          pos = e-p->entries2;
          if(pos != 0) {
            ee = &p->entries[pos-1];
            if(ee->inList != NULL) {
              nn = ee->inList->atNode;
              if(nn->last_source_list != l) {
                pr_append_list(pr_new_list(p),nn);
                nn->last_source_list = l;
              }
              pr_move_entry(ee,nn->last);
            }
          }
#ifdef STATS
          p->num_compares++;
#endif
        }

    /* remove entry length-iteration */
    pr_remove_entry(&p->entries[p->length-i]);

#ifdef STATS
    p->num_compares++;
#endif

  } /* for i */

} /* primitives_find() */

/****** EOF (repeats_primitives.c) *******************************************/

