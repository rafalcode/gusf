/*
 * repeats_nonoverlapping.c
 *
 * This file contains an implementation of an extension of Crochemore's
 * algorithm (1981) which finds all occurrences of nonoverlapping maximal
 * repeats in O(n log n + z) time.
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
#include "repeats_nonoverlapping.h"

/*
 *
 * Forward References.
 *
 */
void no_remove_list(no_list *l);
void no_remove_node(no_node *n);

/*===========================================================================*/
/*
 * no_append_entry
 *
 * Append entry e to list l: O(1) time.
 *
 * Parameters:  e  -  an entry
 *              l  -  a list
 *              d  -  a character class
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void no_append_entry(no_entry *e, no_list *l, int d)
{
  e->next = NULL;
  e->prev = l->last[d];
  e->inList = l;
  if(l->last[d] != NULL)
    l->last[d]->next = e;
  if(l->entries[d] == NULL)
    l->entries[d] = e;
  l->last[d] = e;
  l->len++;
  return;
} /* no_append_entry() */

/*===========================================================================*/
/*
 * no_remove_entry
 *
 * Remove entry e from its list: O(1) time.
 *
 * Parameters:  e  -  an entry
 *              d  -  a character class
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void no_remove_entry(no_entry *e, int d)
{
  if(e->inList != NULL) { /* only if entry was not removed earlier */
    if(e->prev != NULL)
      e->prev->next = e->next;
    else
      e->inList->entries[d] = e->next;
    if(e->next != NULL)
      e->next->prev = e->prev;
    else
      e->inList->last[d] = e->prev;
    e->inList->len--;
    if(e->inList->len == 0)
      no_remove_list(e->inList);
    e->inList = NULL;
  }
  return;
} /* no_remove_entry() */

/*===========================================================================*/
/*
 * no_move_entry
 *
 * Remove entry e from its list and append it to list l: O(1) time.
 *
 * Parameters:  e  -  an entry
 *              l  -  a list
 *              d  -  a character class
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void no_move_entry(no_entry *e, no_list *l, int d)
{
  no_remove_entry(e,d);
  no_append_entry(e,l,d);
  return;
} /* no_move_entry() */

/*===========================================================================*/
/*
 * no_new_list
 *
 * Create new (empty) list: O(1) time.
 *
 * Parameters:  p  -  a nonoverlapping structure
 *
 * Returns:  the new list
 */
/*---------------------------------------------------------------------------*/
no_list *no_new_list(nonoverlapping_struct *p)
{
  int d, number;
  no_list *l_new;

  number = p->nextList++;
  l_new = &p->lists[number];
  l_new->next = l_new->prev = NULL;
  l_new->atNode = NULL;
  l_new->entries = &p->entryList[number*p->alpha_size];
  l_new->last = &p->lastList[number*p->alpha_size];
  for(d=0; d<p->alpha_size; d++)
    l_new->entries[d] = l_new->last[d] = NULL;
  l_new->len = 0;
  return l_new;
} /* no_new_list() */

/*===========================================================================*/
/*
 * no_append_list
 *
 * Append list l to the lists at node n: O(1) time.
 *
 * Parameters:  l  -  a list
 *              n  -  a node
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void no_append_list(no_list *l, no_node *n)
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
} /* no_append_list() */

/*===========================================================================*/
/*
 * no_remove_list
 *
 * Remove list l from its node: O(1) time.
 *
 * Parameters:  l  -  a list
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void no_remove_list(no_list *l)
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
} /* no_remove_list() */

/*===========================================================================*/
/*
 * no_replace_list
 *
 * Replace list l by a copy of itself (l_new): O(length of l) time.
 *
 * Parameters:  p  -  a nonoverlapping structure
 *              l  -  a list
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void no_replace_list(nonoverlapping_struct *p, no_list *l)
{
  int d, number;
  no_list *l_new;
  no_entry *e, *e_new;

  /* create new list (manually in p->lists2) */
  number = p->nextList2++;
  l_new = &p->lists2[number];

  /* replace l by new list */
  l_new->next = l->next;
  l_new->prev = l->prev;
  l_new->atNode = l->atNode;
  l_new->entries = &p->entryList2[number*p->alpha_size];
  l_new->last = &p->lastList2[number*p->alpha_size];
  for(d=0; d<p->alpha_size; d++)
    l_new->entries[d] = l_new->last[d] = NULL;
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
  for(d=0; d<p->alpha_size; d++)
    for(e=l->entries[d]; e!=NULL; e=e->next) {
      e_new = &p->entries2[e-p->entries];
      no_append_entry(e_new,l_new,d);
    }

  return;
} /* no_replace_list() */

/*===========================================================================*/
/*
 * no_new_node
 *
 * Create new (empty) node: O(1) time.
 *
 * Parameters:  p  -  a nonoverlapping structure
 *
 * Returns:  the new node
 */
/*---------------------------------------------------------------------------*/
no_node *no_new_node(nonoverlapping_struct *p)
{
  no_node *n_new;

  n_new = &p->nodelist[p->nextNode++];
  n_new->next = n_new->prev = NULL;
  n_new->inNodelist = NULL;
  n_new->lists = n_new->last = n_new->last_source_list = NULL;
  return n_new;
} /* no_new_node() */

/*===========================================================================*/
/*
 * no_append_node
 *
 * Append node n to nodes of p: O(1) time.
 *
 * Parameters:  p  -  a nonoverlapping structure
 *              n  -  a node
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void no_append_node(nonoverlapping_struct *p, no_node *n)
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
} /* no_append_node() */

/*===========================================================================*/
/*
 * nonoverlapping_prep
 *
 * Create new (empty) nonoverlapping structure p_new: O(1) time.
 *
 * Parameters:  string         -  the string
 *              raw_string     -  the raw string
 *              length         -  the length of the string
 *
 * Returns:  a nonoverlapping structure
 */
/*---------------------------------------------------------------------------*/
nonoverlapping_struct *nonoverlapping_prep(char *string, char *raw_string,
                                           int length)
{
  nonoverlapping_struct *p_new;
  int i, c, len, alph;

  if((p_new = malloc(sizeof(nonoverlapping_struct))) == NULL)
    return NULL;
  memset(p_new, 0, sizeof(nonoverlapping_struct));

  if((p_new->string_space = malloc((length+2) * sizeof(char))) == NULL) {
    nonoverlapping_free(p_new);
    return NULL;
  }
  memset(p_new->string_space, 0, (length+2) * sizeof(char));

  /*
   * Append delimiters at beginning and end of the string.
   */
  p_new->string_space[0] = CHAR_MAX;
  for(i=0; i<length; i++)
    p_new->string_space[i+1] = string[i];
  p_new->string_space[length+1] = CHAR_MAX;
  p_new->string = &p_new->string_space[1];

  p_new->raw_string = raw_string;
  len = p_new->length = length + 1;

  /*
   * Find alphabet size and a pointers
   */
  for(i=0; i<=CHAR_MAX; i++)
    p_new->a[i] = -1;
  alph = 0;
  for(i=-1; i<len; i++) { /* including delimiters */
    c = (int)(unsigned char)p_new->string[i];
    if(p_new->a[c] == -1)
      p_new->a[c] = alph++;
  }
  p_new->alpha_size = alph;

  if((p_new->entries = malloc(len * sizeof(no_entry))) == NULL) {
    nonoverlapping_free(p_new);
    return NULL;
  }
  memset(p_new->entries, 0, len * sizeof(no_entry));

  if((p_new->entries2 = malloc(len * sizeof(no_entry))) == NULL) {
    nonoverlapping_free(p_new);
    return NULL;
  }
  memset(p_new->entries2, 0, len * sizeof(no_entry));

  if((p_new->entryList = malloc(2*len*alph * sizeof(no_entry*))) == NULL) {
    nonoverlapping_free(p_new);
    return NULL;
  }
  memset(p_new->entryList, 0, 2*len*alph * sizeof(no_entry*));

  if((p_new->entryList2 = malloc(2*len*alph * sizeof(no_entry*))) ==NULL) {
    nonoverlapping_free(p_new);
    return NULL;
  }
  memset(p_new->entryList2, 0, 2*len*alph * sizeof(no_entry*));

  if((p_new->lastList = malloc(2*len*alph * sizeof(no_entry*))) == NULL) {
    nonoverlapping_free(p_new);
    return NULL;
  }
  memset(p_new->lastList, 0, 2*len*alph * sizeof(no_entry*));

  if((p_new->lastList2 = malloc(2*len*alph * sizeof(no_entry*))) == NULL) {
    nonoverlapping_free(p_new);
    return NULL;
  }
  memset(p_new->lastList2, 0, 2*len*alph * sizeof(no_entry*));

  if((p_new->lists = malloc(2*len * sizeof(no_list))) == NULL) {
    nonoverlapping_free(p_new);
    return NULL;
  }
  memset(p_new->lists, 0, 2*len * sizeof(no_list));

  if((p_new->lists2 = malloc(len * sizeof(no_list))) == NULL) {
    nonoverlapping_free(p_new);
    return NULL;
  }
  memset(p_new->lists2, 0, len * sizeof(no_list));

  p_new->nextList = p_new->nextList2 = 0;

  if((p_new->nodelist = malloc(len * sizeof(no_node))) == NULL) {
    nonoverlapping_free(p_new);
    return NULL;
  }
  memset(p_new->nodelist, 0, len * sizeof(no_node));

  if((p_new->nodelist2 = malloc(len * sizeof(no_node))) == NULL) {
    nonoverlapping_free(p_new);
    return NULL;
  }
  memset(p_new->nodelist2, 0, len * sizeof(no_node));

  p_new->nextNode = 0;
  p_new->nodes = p_new->nodes2 = p_new->last, p_new->last2 = NULL;

  p_new->num_nonoverlapping_maximal_pairs = 0;

#ifdef STATS
  p_new->num_compares = 0;
#endif

  return p_new;

} /* nonoverlapping_prep() */

/*===========================================================================*/
/*
 * no_next_level
 *
 * Step forward to next level (i.e. swap nodelists): O(1) time.
 *
 * Parameters:  p  -  a nonoverlapping structure
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void no_next_level(nonoverlapping_struct *p)
{
  no_node *n_tmp;

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
} /* no_next_level() */

/*===========================================================================*/
/*
 * nonoverlapping_free
 *
 * Delete nonoverlapping structure p: O(1) time.
 *
 * Parameters:  p  -  a nonoverlapping structure
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void nonoverlapping_free(nonoverlapping_struct *p)
{
  free(p->string_space);
  free(p->entries);
  free(p->entries2);
  free(p->entryList);
  free(p->entryList2);
  free(p->lastList);
  free(p->lastList2);
  free(p->lists);
  free(p->lists2);
  free(p->nodelist);
  free(p->nodelist2);
  free(p);
  return;
} /* nonoverlapping_free() */

/*===========================================================================*/
/*
 * no_write
 *
 * Write a nonoverlapping maximal repeat.
 *
 * Parameters:  raw_string  -  the raw string
 *              pos1        -  starting position of first occurrence
 *              pos2        -  starting position of second occurrence
 *              len         -  length of the repeat
 *              type        -  repeat type
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void no_write(char *raw_string, int pos1, int pos2, int len, char *type)
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
} /* no_write() */

/*===========================================================================*/
/*
 * no_report
 *
 * Report all nonoverlapping maximal repeats in this iteration.
 *
 * Parameters:  p          -  a nonoverlapping structure
 *              iteration  -  iteration
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void no_report(nonoverlapping_struct *p, int iteration)
{
  no_node *n;
  no_list *l,*ll;
  int d,dd;
  no_entry *e,*ee;

  for(n=p->nodes; n!=NULL; n=n->next)
    for(l=n->lists; l!=NULL; l=l->next)
      for(d=0; d<p->alpha_size; d++)
        for(e=l->entries[d]; e!=NULL; e=e->next)
          for(ll=n->lists; ll!=NULL; ll=ll->next)
            if(ll != l)
              for(dd=0; dd<p->alpha_size; dd++)
                if(dd != d)
                  for(ee=ll->last[dd];
                      ee!=NULL && ee-e>=iteration;
                      ee=ee->prev) {
                    no_write(p->raw_string,e-p->entries,ee-p->entries,
                             iteration,"nonoverlapping maximal pair");
                    p->num_nonoverlapping_maximal_pairs++;
                  }
  return;
} /* no_report() */

/*===========================================================================*/
/*
 * no_create_basic_lists
 *
 * Create basic lists: O(alpha_size) space, O(length + alpha_size) time.
 *
 * Parameters:  p  -  a nonoverlapping structure
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void no_create_basic_lists(nonoverlapping_struct *p)
{
  int i,c,d;
  no_list **occ;
  no_node *n;

  occ = malloc(p->alpha_size * sizeof(no_list*));
  for(i=0; i<p->alpha_size; i++)
    occ[i] = NULL;

  n = no_new_node(p);
  no_append_node(p,n);
  for(i=0; i<p->length ; i++) { /* including right delimiter */
    c = p->a[(int)(unsigned char)p->string[i]];
    d = p->a[(int)(unsigned char)p->string[i-1]];
    if(occ[c] == NULL)
      no_append_list((occ[c]=no_new_list(p)),n);
    no_append_entry(&p->entries[i],occ[c],d);
  }

  free(occ);

  return;
} /* no_create_basic_lists() */

/*===========================================================================*/
/*
 * nonoverlapping_find
 *
 * This is the modified version of Crochemore's algorithm.
 *
 * Parameters:  p  -  a nonoverlapping structure
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void nonoverlapping_find(nonoverlapping_struct *p)
{
  no_node *n,*n_next,*nn,*new_n;
  no_list *l,*l_next,*maxlist;
  no_entry *e,*ee;
  int i,d, maxlistlen, pos;

  /* create basic lists */
  no_create_basic_lists(p);

  for(i=1; i<p->length && p->nodes!=NULL; i++) {

    /* step forward to next level */
    no_next_level(p);

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
          no_remove_list(l);
        else
          no_replace_list(p,l);
        if(l->len == 1) { /* don't remove brackets! */
          for(d=0; d<p->alpha_size; d++)
            if(l->entries[d] != NULL)
              l->entries[d]->inList = NULL; /* mark entry removed */
        }
        else {
          new_n = no_new_node(p);
          no_append_node(p,new_n);
          no_append_list(l,new_n);
        }
      } /* for l */
    } /* for n */

    /* using the remaining (small) lists, pull entries out */
    for(n=p->nodes2; n!=NULL; n=n->next)
      for(l=n->lists; l!=NULL; l=l->next)
        for(d=0; d<p->alpha_size; d++)
          for(e=l->entries[d]; e!=NULL; e=e->next) {
            pos = e-p->entries2;
            if(pos != 0) {
              ee = &p->entries[pos-1];
              if(ee->inList != NULL) {
                nn = ee->inList->atNode;
                if(nn->last_source_list != l) {
                  no_append_list(no_new_list(p),nn);
                  nn->last_source_list = l;
                }
                no_move_entry(ee,nn->last,
                              p->a[(int)(unsigned char)p->string[pos-2]]);
              }
            }
#ifdef STATS
            p->num_compares++;
#endif
          }

    /* report nonoverlapping repeats (before removing the last entry) */
    no_report(p,i);

    /* remove entry length-iteration */
    no_remove_entry(&p->entries[p->length-i],
                    p->a[(int)(unsigned char)p->string[p->length-i-1]]);

#ifdef STATS
    p->num_compares++;
#endif

  } /* for i */

} /* nonoverlapping_find() */

/****** EOF (repeats_nonoverlapping.c) ***************************************/

