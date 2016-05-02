/*
 * repeats_vocabulary.c
 *
 * Find the vocabulary of tandem repeats in O(n) time.
 * Find the vocabulary of primitive tandem repeats in O(n) time.
 * Find the vocabulary of primitive tandem arrays in O(n) time.
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
#ifdef STRMAT
#include "stree_strmat.h"
#else
#include "stree.h"
#endif
#include "more.h"
#include "stree_decomposition.h"
#include "repeats_vocabulary.h"

#define MAX(X,Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))
#define MIN3(X,Y,Z) (MIN(X,MIN(Y,Z)))
#define ABS(X) (((X) < 0) ? -(X) : (X))

/*===========================================================================*/
/*
 * vocabulary_prep
 *
 * Preprocessing for the tandem vocabulary algorithm.
 *
 * Parameters:  tree           -  a suffix tree
 *              decomposition  -  a decomposition structure
 *              string         -  the string
 *              raw_string     -  the raw string
 *              length         -  the length of the string
 *
 * Returns:  A vocabulary structure
 */
/*---------------------------------------------------------------------------*/
VOCABULARY_STRUCT *vocabulary_prep(SUFFIX_TREE tree,
                                   DECOMPOSITION_STRUCT *decomposition,
                                   char *string, char *raw_string, int length)
{
  VOCABULARY_STRUCT *vocabulary;
  unsigned int num_nodes;
  int i, max_block_length;

  if(tree==NULL || decomposition==NULL)
    return NULL;

  /*
   * Allocate the memory.
   */
  if((vocabulary = malloc(sizeof(VOCABULARY_STRUCT))) == NULL)
    return NULL;
  memset(vocabulary, 0, sizeof(VOCABULARY_STRUCT));

  vocabulary->string = string;
  vocabulary->raw_string = raw_string;
  vocabulary->length = length;

  vocabulary->tree = tree;
  vocabulary->decomposition = decomposition;

  num_nodes = stree_get_num_nodes(tree);
  max_block_length = decomposition_get_max_block_length(decomposition);

  /* for each possible tandem repeat */
  if((vocabulary->tandem_space = malloc(2*length * sizeof(tandem))) == NULL) {
    vocabulary_free(vocabulary);
    return NULL;
  }
  memset(vocabulary->tandem_space, 0, 2*length * sizeof(tandem));

  vocabulary->next_tandem = 0;

  /* for each position of the text */
  if((vocabulary->tandems = malloc(length * sizeof(tandem*))) == NULL) {
    vocabulary_free(vocabulary);
    return NULL;
  }
  memset(vocabulary->tandems, 0, length * sizeof(tandem*));

  for(i=0; i<length; i++)
    vocabulary->tandems[i] = NULL;

  /* for each node (not only internal) */
  if((vocabulary->tlens1 = malloc(num_nodes * sizeof(int))) == NULL) {
    vocabulary_free(vocabulary);
    return NULL;
  }
  memset(vocabulary->tlens1, 0, num_nodes * sizeof(int));

  if((vocabulary->tlens2 = malloc(num_nodes * sizeof(int))) == NULL) {
    vocabulary_free(vocabulary);
    return NULL;
  }
  memset(vocabulary->tlens2, 0, num_nodes * sizeof(int));

  for(i=0; i<num_nodes; i++)
    vocabulary->tlens1[i] = vocabulary->tlens2[i] = 0;

  /* for each possible depth of the suffix tree (actually 0..depth=length+1) */
  if((vocabulary->dvector = malloc((length+1) * sizeof(int))) == NULL) {
    vocabulary_free(vocabulary);
    return NULL;
  }
  memset(vocabulary->dvector, 0, (length+1) * sizeof(int));

  /* for the longest two blocks of the decomposition */
  if((vocabulary->PREF = malloc(2*max_block_length * sizeof(int))) == NULL) {
    vocabulary_free(vocabulary);
    return NULL;
  }
  memset(vocabulary->PREF, 0, 2*max_block_length * sizeof(int));

  if((vocabulary->PREF2 = malloc(2*max_block_length * sizeof(int))) == NULL) {
    vocabulary_free(vocabulary);
    return NULL;
  }
  memset(vocabulary->PREF2, 0, 2*max_block_length * sizeof(int));

  if((vocabulary->SUFF = malloc(2*max_block_length * sizeof(int))) == NULL) {
    vocabulary_free(vocabulary);
    return NULL;
  }
  memset(vocabulary->SUFF, 0, 2*max_block_length * sizeof(int));

  vocabulary->num_tandem_repeats = 0;
  vocabulary->num_primitive_tandem_repeats = 0;
  vocabulary->num_tandem_arrays = 0;

  vocabulary->num_tandem_repeat_occs = 0;
  vocabulary->num_primitive_tandem_repeat_occs = 0;
  vocabulary->num_tandem_array_occs = 0;

#ifdef STATS
  vocabulary->num_prep = tree->num_compares + decomposition->num_compares;
  vocabulary->num_compares_for_tandem_repeats = 0;
  vocabulary->num_compares_for_primitive_tandem_repeats = 0;
  vocabulary->num_compares_for_tandem_arrays = 0;
#endif

  return vocabulary;

} /* vocabulary_prep() */

/*===========================================================================*/
/*
 * vocabulary_free
 *
 * Free memory of a vocabulary structure.
 *
 * Parameters:  vocabulary  -  a vocabulary structure.
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void vocabulary_free(VOCABULARY_STRUCT *vocabulary)
{
  if(vocabulary->tandem_space != NULL)
    free(vocabulary->tandem_space);
  if(vocabulary->tandems!= NULL)
    free(vocabulary->tandems);
  if(vocabulary->tlens1 != NULL)
    free(vocabulary->tlens1);
  if(vocabulary->tlens2 != NULL)
    free(vocabulary->tlens2);
  if(vocabulary->dvector != NULL)
    free(vocabulary->dvector);
  if(vocabulary->PREF != NULL)
    free(vocabulary->PREF);
  if(vocabulary->PREF2 != NULL)
    free(vocabulary->PREF2);
  if(vocabulary->SUFF != NULL)
    free(vocabulary->SUFF);
  free(vocabulary);

} /* vocabulary_free() */

/*===========================================================================*/
/*
 * vocabulary_tandem_insert
 *
 * Prepend tandem repeat of length len before the tandem list no. pos.
 *
 * Parameters:  v    -  a vocabulary structure
 *              pos  -  starting position of the tandem repeat
 *              len  -  length of the tandem repeat
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void vocabulary_tandem_insert(VOCABULARY_STRUCT *v, int pos, int len)
{
  tandem *new_t;

  new_t = &v->tandem_space[v->next_tandem++];
  new_t->len = len;
  new_t->next = v->tandems[pos];
  v->tandems[pos] = new_t;

} /* vocabulary_tandem_insert() */

/*===========================================================================*/
/*
 * leftreps
 *
 * Find all tandem repeats in uv with midpoint in u (but not between u and v)
 * where u = v->string[pos1..pos2-1] and v = v->string[pos2..pos3-1].
 *
 * Parameters:  voc   -  a vocabulary structure
 *              pos1  -  position of the first character of u
 *              pos2  -  position of the first character of v
 *              pos3  -  position after the last character of v
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void leftreps(VOCABULARY_STRUCT *voc, int pos1, int pos2, int pos3)
{
  int l, r, k, j, p, beta, kprime, ulen, ulast, vlen;
  int *PREF, *PREF2, *SUFF;
  char *u, *v;

  u = &(voc->string[pos1]);
  v = &(voc->string[pos2]);
  ulen = pos2 - pos1;
  ulast = ulen - 1;
  vlen = pos3 - pos2;

  PREF = voc->PREF;
  PREF2 = voc->PREF2;
  SUFF = voc->SUFF;

  /*
   * Compute PREF (the Z algorithm backwards).
   */
  l = r = 0;
  PREF[0] = 0;
  for(k=1; k<ulen; k++) {
    if(k > r) {
      for(j=0; k+j<ulen && u[ulast-j]==u[ulast-k-j]; j++)
        /* nothing */ ;
      PREF[k] = j;
      r = k + j;
      l = k;
    }
    else {
      beta = r - k;
      kprime = k - l;

      if(PREF[kprime] < beta)
        PREF[k] = PREF[kprime];
      else {
        for(j=0; r+j<ulen && u[ulast-r-j]==u[ulast-beta-j]; j++)
          /* nothing */ ;
        PREF[k] = beta + j;
        r += j;
        l = k;
      }
    }
  }

  /*
   * Compute PREF2 (the Z algorithm on the first ulen letters of v).
   */
  l = r = 0;
  PREF2[0] = 0;
  for(k=1; k<MIN(ulen,vlen); k++) {
    if(k > r) {
      for(j=0; k+j<MIN(ulen,vlen) && v[j]==v[k+j]; j++)
        /* nothing */ ;
      PREF2[k] = j;
      r = k + j;
      l = k;
    }
    else {
      beta = r - k;
      kprime = k - l;

      if(PREF2[kprime] < beta)
        PREF2[k] = PREF2[kprime];
      else {
        for(j=0; r+j<MIN(ulen,vlen) && v[r+j]==v[beta+j]; j++)
          /* nothing */ ;
        PREF2[k] = beta + j;
        r += j;
        l = k;
      }
    }
  }

  /*
   * Compute SUFF (the Z algorithm on two strings).
   */
  l = r = -1;
  for(k=0; k<ulen; k++) {
    if(k > r) {
      for(j=0; j<vlen && k+j<ulen && v[j]==u[k+j]; j++)
        /* nothing */ ;
      SUFF[k] = j;
      r = k + j;
      l = k;
    }
    else {
      beta = r - k;
      kprime = k - l;

      if(PREF2[kprime] < beta)
        SUFF[k] = PREF2[kprime];
      else {
        for(j=0; r+j<ulen && beta+j<vlen && u[r+j]==v[beta+j]; j++)
          /* nothing */ ;
        SUFF[k] = beta + j;
        r += j;
        l = k;
      }
    }
  }

#ifdef STATS
  voc->num_compares_for_tandem_repeats += ulen;
#endif

  /*
   * Find the starting position p of the left-most tandem repeat of each run:
   *
   * There is a run of tandem repeats of length l if and only if
   *
   *   SUFF[ulen-l] + PREF[l] >= l
   *   =>  - PREF[l] <= SUFF[ulen-l] - l
   *   =>  pos2 - PREF[l] <= pos2 + SUFF[ulen-l] - l
   *   =>  pos2 - PREF[l] - l <= pos2 + SUFF[ulen-l] - 2*l
   *
   * With the additional constraints
   *
   *  p >= pos2 - 2*l + 1   to ensure that the tandem repeat touches v
   *  p <= pos2 - l - 1     to ensure that the center is in u
   *
   */
  for(l=1; l<ulen; l++)
    if((p=MAX(pos2-PREF[l]-l,pos2-2*l+1))<=MIN(pos2+SUFF[ulen-l]-2*l,pos2-l-1))
      vocabulary_tandem_insert(voc,p,2*l);

} /* leftreps() */

/*===========================================================================*/
/*
 * rightreps
 *
 * Find all tandem repeats in uvw with midpoint in v or between u and v
 * where u = t->string[pos1..pos2-1], v = t->string[pos2..posM-1], and
 * w = t->string[posM..pos3-1].
 *
 * Parameters:  voc   -  a vocabulary structure
 *              pos1  -  position of the first character of u
 *              pos2  -  position of the first character of v
 *              posM  -  position of the first character of w
 *              pos3  -  position after the last character of w
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void rightreps(VOCABULARY_STRUCT *voc, int pos1, int pos2, int posM, int pos3)
{
  int l, r, k, j, p, beta, kprime, ulen, ulast, vlen, vlast;
  int *PREF, *PREF2, *SUFF;
  char *u, *v;

  u = &(voc->string[pos1]);
  v = &(voc->string[pos2]);
  ulen = pos2 - pos1;
  ulast = ulen - 1;
  vlen = pos3 - pos2;
  vlast = vlen - 1;

  PREF = voc->PREF;
  PREF2 = voc->PREF2;
  SUFF = voc->SUFF;

  /*
   * Compute PREF (the Z algorithm).
   */
  l = r = 0;
  PREF[0] = 0;
  for(k=1; k<vlen; k++) {
    if(k > r) {
      for(j=0; k+j<vlen && v[j]==v[k+j]; j++)
        /* nothing */ ;
      PREF[k] = j;
      r = k + j;
      l = k;
    }
    else {
      beta = r - k;
      kprime = k - l;

      if(PREF[kprime] < beta)
        PREF[k] = PREF[kprime];
      else {
        for(j=0; r+j<vlen && v[r+j]==v[beta+j]; j++)
          /* nothing */ ;
        PREF[k] = beta + j;
        r += j;
        l = k;
      }
    }
  }

  /*
   * Compute PREF2 (the Z algorithm backwards on the last vlen letters of u).
   */
  l = r = 0;
  PREF2[0] = 0;
  for(k=1; k<MIN(ulen,vlen); k++) {
    if(k > r) {
      for(j=0; k+j<MIN(ulen,vlen) && u[ulast-j]==u[ulast-k-j]; j++)
        /* nothing */ ;
      PREF2[k] = j;
      r = k + j;
      l = k;
    }
    else {
      beta = r - k;
      kprime = k - l;

      if(PREF2[kprime] < beta)
        PREF2[k] = PREF2[kprime];
      else {
        for(j=0; r+j<MIN(ulen,vlen) && u[ulast-r-j]==u[ulast-beta-j]; j++)
          /* nothing */ ;
        PREF2[k] = beta + j;
        r += j;
        l = k;
      }
    }
  }

  /*
   * Compute SUFF (the Z algorithm in backward direction on two strings).
   */
  l = r = -1;
  for(k=0; k<vlen; k++) {
    if(k > r) {
      for(j=0; j<ulen && k+j<vlen && u[ulast-j]==v[vlast-k-j]; j++)
        /* nothing */ ;
      SUFF[k] = j;
      r = k + j;
      l = k;
    }
    else {
      beta = r - k;
      kprime = k - l;

      if(PREF2[kprime] < beta)
        SUFF[k] = PREF2[kprime];
      else {
        for(j=0;r+j<vlen && beta+j<ulen && v[vlast-r-j]==u[ulast-beta-j];j++)
          /* nothing */ ;
        SUFF[k] = beta + j;
        r += j;
        l = k;
      }
    }
  }

#ifdef STATS
  voc->num_compares_for_tandem_repeats += vlen;
#endif

  /*
   * Find the starting position p of the left-most tandem repeat of each run:
   *
   * There is a run of tandem repeats of length l if and only if
   *
   *   PREF[l] + SUFF[vlen-l] >= l
   *   =>  - SUFF[vlen-l] <= PREF[l] - l
   *   =>  pos2 - SUFF[vlen-l] <= pos2 + PREF[l] - l
   *
   * With the additional constraints
   *
   *  p >= pos2 - l      to ensure that the center is in v or between u and v
   *  p <= pos2 - 1      to ensure that the tandem repeat touches u
   *  p <= posM - l - 1  to ensure that the center is to the left of posM
   *
   */
  for(l=1; l<vlen; l++)
    if((p=MAX(pos2-SUFF[vlen-l],pos2-l))<=MIN3(pos2+PREF[l]-l,pos2-1,posM-l-1))
      vocabulary_tandem_insert(voc,p,2*l);

  /* test additionally for l == vlen */
  if((p=pos2-SUFF[0]) <= MIN(pos2-vlen,posM-vlen-1))
    vocabulary_tandem_insert(voc,p,2*vlen);

} /* rightreps() */

/*===========================================================================*/
/*
 * vocabulary_tloc_insert
 *
 * Mark tandem repeat/tandem array location in the suffix tree.
 *
 * Parameters:  v     -  a vocabulary structure
 *              node  -  the node below the location
 *              len   -  the length of the tandem repeat/tandem array
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void vocabulary_tloc_insert(VOCABULARY_STRUCT *v, STREE_NODE node, int len)
{
  int id;

  id = stree_get_ident(v->tree,node);

  if(v->tlens1[id] == 0)
    v->tlens1[id] = len;
  else if(v->tlens2[id] == 0)
    v->tlens2[id] = len;
  else {
    fprintf(stderr,"ERROR: Three tandem repeats in an edge!!!\n");
    exit(1);
  }
} /* vocabulary_tloc_insert() */

/*===========================================================================*/
/*
 * vocabulary_collect
 *
 * Find locations of a leftmost-covering set of tandem repeats in the tree: 
 * bottom-up keep the list from the leaf with the smallest label.
 *
 * Parameters:  v           -  a vocabulary structure
 *              node        -  a suffix tree node
 *              depth       -  the string-depth of node
 *              curr_pos    -  the pos-label of the currently smallest leaf
 *              curr_tlist  -  pointer into the current list of tandem repeats
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void vocabulary_collect(VOCABULARY_STRUCT *v, STREE_NODE node, int depth,
                        int *curr_pos, tandem **curr_tlist)
{
  STREE_NODE child;
  int i, edgelen, leavesnum, new_pos, dummy_id;
  tandem *new_tlist;
  char *dummy_string;

  *curr_pos = INT_MAX;

  /* process subtrees */
  for(child = stree_get_children(v->tree,node);
      child != NULL;
      child = stree_get_next(v->tree,child)) {

    edgelen = stree_get_edgelen(v->tree,child);

    /* recurse depth-first */
    vocabulary_collect(v,child,depth+edgelen,&new_pos,&new_tlist);

    /* store location of tandems which are inside this edge */
    while(new_tlist!=NULL && new_tlist->len>depth) {
      vocabulary_tloc_insert(v,child,new_tlist->len);
      new_tlist = new_tlist->next;
    }

    /* keep list from smallest leaf as return value */
    if(new_pos < *curr_pos) {
      *curr_pos = new_pos;
      *curr_tlist = new_tlist;
    }
  }

  /* then process direct leaves */
  leavesnum = stree_get_num_leaves(v->tree,node);
  for(i=1; i<=leavesnum; i++) {
    stree_get_leaf(v->tree,node,i,&dummy_string,&new_pos,&dummy_id);
    new_tlist = v->tandems[new_pos];

    /* keep list from smallest leaf as return value */
    if(new_pos < *curr_pos) {
      *curr_pos = new_pos;
      *curr_tlist = new_tlist;
    }
  }

} /* vocabulary_collect() */

/*===========================================================================*/
/*
 * vocabulary_sl_walk
 *
 * Perform suffix-link walk: Find suffix-link location of a location in the
 * suffix tree, and test for continuation (rotation).
 * If successful, continue.
 *
 * Parameters:  v            -  a vocabulary structure
 *              node         -  the node above the current location
 *              depth        -  the string-depth of node
 *              child        -  the node below the current location
 *              tlen         -  the string-depth of the current location
 *              period       -  the period of the tandem repeat/tandem array
 *              num_inserts  -  counter for number of inserts
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void vocabulary_sl_walk(VOCABULARY_STRUCT *v, STREE_NODE node, int depth,
                        STREE_NODE child, int tlen, int period,
                        int *num_inserts)
{
  int offset, edgelen, id;
  char *edgestr;

  /* save old edge label */
  edgestr = stree_get_edgestr(v->tree,child);
  offset = tlen - depth;

  /* follow suffix link */
  node = stree_get_suffix_link(v->tree,node);
  depth--;
  child = stree_find_child(v->tree,node,*edgestr);

  /* canonize (skip/count) */
  while(offset > (edgelen=stree_get_edgelen(v->tree,child))) {
    node = child;
    depth += edgelen;
    offset -= edgelen;
    edgestr += edgelen;
    child = stree_find_child(v->tree,node,*edgestr);
  }

  /* prepare test for continuation */
  id = stree_get_ident(v->tree,child);
  edgestr = stree_get_edgestr(v->tree,child);
  /* edgelen was already set in while() clause */

  /* if necessary, find next edge */
  if(offset == edgelen) {
    node = child;
    depth += edgelen;
    offset -= edgelen;
    edgestr += edgelen;
    child = stree_find_child(v->tree,node,*(edgestr-period));
    if(child != NULL) {
      id = stree_get_ident(v->tree,child);
      edgestr = stree_get_edgestr(v->tree,child);
    }
  }

  /* test if continuation exists */
  if(child!=NULL && *(edgestr+offset)==*(edgestr+offset-period) &&
     ABS(v->tlens1[id])!=tlen && ABS(v->tlens2[id])!=tlen) {
    vocabulary_tloc_insert(v,child,-tlen);
    (*num_inserts)++;
    vocabulary_sl_walk(v,node,depth,child,tlen,period,num_inserts);
  }
} /* vocabulary_sl_walk() */

/*===========================================================================*/
/*
 * vocabulary_mod_sl_walk
 *
 * Perform modified suffix-link walk: Find location in string-depth tlen whose
 * prefix is the right-rotation of the current location.
 * Continue even if unsuccessful, but only steps times.
 *
 * Parameters:  v            -  a vocabulary structure
 *              node         -  the node above the current location
 *              depth        -  the string-depth of node
 *              child        -  the node below the current location
 *              tlen         -  the string-depth of the seeked location
 *              period       -  the period of the tandem repeat/tandem array
 *              offset       -  edge-offset of the current location
 *              steps        -  number of steps to be performed
 *              num_inserts  -  counter for number of inserts
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void vocabulary_mod_sl_walk(VOCABULARY_STRUCT *v, STREE_NODE node, int depth,
                            STREE_NODE child, int tlen, int period, int offset,
                            int steps, int *num_inserts)
{
  STREE_NODE cTmp;
  int o, edgelen;
  char *edgestr;

  if(steps > 1) {

    /* save old edge label */
    edgestr = stree_get_edgestr(v->tree,child);

    /* follow suffix link */
    node = stree_get_suffix_link(v->tree,node);
    depth--;
    child = stree_find_child(v->tree,node,*edgestr);

    /* canonize (skip/count) */
    o = 0;
    while(offset > (edgelen=stree_get_edgelen(v->tree,child))) {
      node = child;
      depth += edgelen;
      offset -= edgelen;
      o += edgelen;
      child = stree_find_child(v->tree,node,*(edgestr+o));
    }

    /* prepare stepping down */
    edgestr = stree_get_edgestr(v->tree,child);

    /* if necessary, find next edge */
    cTmp = child;
    if(offset == edgelen) {
      cTmp = stree_find_child(v->tree,child,*(edgestr+offset-period));
      if(cTmp != NULL) {
        node = child;
        depth += edgelen;
        offset -= edgelen;
        child = cTmp;
        edgelen = stree_get_edgelen(v->tree,child);
        edgestr = stree_get_edgestr(v->tree,child);
      }
    }

    /* walk down (single steps) */
    while(cTmp!=NULL && *(edgestr+offset)==*(edgestr+offset-period) &&
          depth+offset<tlen) {
      offset++;

#ifdef STATS
      v->num_compares_for_tandem_arrays += 2;
#endif

      if(depth+offset == tlen) {
        vocabulary_tloc_insert(v,child,-tlen);
        (*num_inserts)++;
      }

      if(offset == edgelen) {
        cTmp = stree_find_child(v->tree,child,*(edgestr+edgelen-period));
        if(cTmp != NULL) {
          node = child;
          depth += edgelen;
          offset -= edgelen;
          child = cTmp;
          edgelen = stree_get_edgelen(v->tree,child);
          edgestr = stree_get_edgestr(v->tree,child);
        }
      }
    }

    vocabulary_mod_sl_walk(v,node,depth,child,tlen,period,offset,steps-1,
                           num_inserts);

  }
} /* vocabulary_mod_sl_walk() */

/*===========================================================================*/
/*
 * vocabulary_rotate
 *
 * Find rotated tandem repeats using suffix-link walks.
 *
 * Parameters:  v      -  a vocabulary structure
 *              node   -  a suffix tree node
 *              depth  -  the string-depth of node
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void vocabulary_rotate(VOCABULARY_STRUCT *v, STREE_NODE node, int depth)
{
  STREE_NODE child;
  int id, edgelen, tlen, dummy=0;

  /* depth-first */
  for(child=stree_get_children(v->tree,node);
      child!=NULL;
      child=stree_get_next(v->tree,child)) {

    id = stree_get_ident(v->tree,child);
    edgelen = stree_get_edgelen(v->tree,child);

    if((tlen=v->tlens1[id]) > 0)
      vocabulary_sl_walk(v,node,depth,child,tlen,tlen/2,&dummy);
    if((tlen=v->tlens2[id]) > 0)
      vocabulary_sl_walk(v,node,depth,child,tlen,tlen/2,&dummy);

    vocabulary_rotate(v,child,depth+edgelen);
  }

} /* vocabulary_rotate() */

/*===========================================================================*/
/*
 * vocabulary_capitalize
 *
 * Change all tandem repeat/tandem array entries to positive values.
 *
 * Parameters:  v     -  a vocabulary structure
 *              node  -  a suffix tree node
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void vocabulary_capitalize(VOCABULARY_STRUCT *v, STREE_NODE node)
{
  STREE_NODE child;
  int id, tlen;

  /* depth-first */
  for(child=stree_get_children(v->tree,node);
      child!=NULL;
      child=stree_get_next(v->tree,child)) {

    id = stree_get_ident(v->tree,child);

    if((tlen=v->tlens1[id]) < 0)
      v->tlens1[id] = -tlen;
    if((tlen=v->tlens2[id]) < 0)
      v->tlens2[id] = -tlen;

    vocabulary_capitalize(v,child);
  }

} /* vocabulary_capitalize() */

/*===========================================================================*/
/*
 * vocabulary_filter
 *
 * Filter out the non-primitive tandem repeats.
 *
 * Parameters:  v     -  a vocabulary structure
 *              node  -  a suffix tree node
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void vocabulary_filter(VOCABULARY_STRUCT *v, STREE_NODE node)
{
  STREE_NODE child;
  int id, tlen, prim_tlen;

  /* depth-first */
  for(child=stree_get_children(v->tree,node);
      child!=NULL;
      child=stree_get_next(v->tree,child)) {

    id = stree_get_ident(v->tree,child);

    if((tlen=v->tlens1[id]) > 0) {
      if((prim_tlen=v->dvector[tlen]) > 0) {
        if(tlen+prim_tlen <= v->length)
          v->dvector[tlen+prim_tlen] = prim_tlen;
      }
      else
        if(2*tlen <= v->length)
          v->dvector[2*tlen] = tlen;

#ifdef STATS
      v->num_compares_for_primitive_tandem_repeats += 4;
#endif
    }

    if((tlen=v->tlens2[id]) > 0) {
      if((prim_tlen=v->dvector[tlen]) > 0) {
        if(tlen+prim_tlen <= v->length)
          v->dvector[tlen+prim_tlen] = prim_tlen;
      }
      else
        if(2*tlen <= v->length)
          v->dvector[2*tlen] = tlen;

#ifdef STATS
      v->num_compares_for_primitive_tandem_repeats += 4;
#endif
    }

    vocabulary_filter(v,child);

    if((tlen=v->tlens1[id]) > 0) {
      if((prim_tlen=v->dvector[tlen]) > 0) {
        if(tlen+prim_tlen <= v->length) {
          v->dvector[tlen+prim_tlen] = 0;
        }
        v->tlens1[id] = 0;
      }
      else {
        if(2*tlen <= v->length) {
          v->dvector[2*tlen] = 0;
        }
      }
    }

    if((tlen=v->tlens2[id]) > 0) {
      if((prim_tlen=v->dvector[tlen]) > 0) {
        if(tlen+prim_tlen <= v->length) {
          v->dvector[tlen+prim_tlen] = 0;
        }
        v->tlens2[id] = 0;
      }
      else {
        if(2*tlen <= v->length) {
          v->dvector[2*tlen] = 0;
        }
      }
    }
  }
} /* vocabulary_filter() */

/*===========================================================================*/
/*
 * vocabulary_minimize
 *
 * Filter out those tandem repeats which are not in a minimal leftmost
 * covering set.
 *
 * Parameters:  v       -  a vocabulary structure
 *              node    -  the node above the current location
 *              depth   -  the string-depth of node
 *              child   -  the child below the current location
 *              tlen    -  the length of the current tandem repeat
 *              cStart  -  the child below the location where we started
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void vocabulary_minimize(VOCABULARY_STRUCT *v, STREE_NODE node, int depth,
                         STREE_NODE child, int tlen, STREE_NODE cStart)
{
  int offset, period, id, edgelen;
  char *edgestr;

  period = tlen / 2;
  offset = tlen - depth;
  edgestr = stree_get_edgestr(v->tree,child);

  /* follow suffix link */
  node = stree_get_suffix_link(v->tree,node);
  depth--;
  child = stree_find_child(v->tree,node,*edgestr);

  /* canonize (skip/count) */
  while(offset > (edgelen=stree_get_edgelen(v->tree,child))) {
    node = child;
    depth += edgelen;
    offset -= edgelen;
    edgestr += edgelen;
    child = stree_find_child(v->tree,node,*edgestr);
  }

  /* prepare test for continuation */
  id = stree_get_ident(v->tree,child);

  /* if necessary, find next edge */
  if(offset == edgelen) {
    node = child;
    depth += edgelen;
    offset -= edgelen;
    edgestr += edgelen;
    child = stree_find_child(v->tree,node,*(edgestr-period));
    if(child != NULL)
      id = stree_get_ident(v->tree,child);
  }


  if(child!=NULL && child!=cStart && v->tlens1[id]==tlen) {
    v->tlens1[id] = 0;
    vocabulary_minimize(v,node,depth,child,tlen,cStart);
  }
  if(child!=NULL && child!=cStart && v->tlens2[id]==tlen) {
    v->tlens2[id] = 0;
    vocabulary_minimize(v,node,depth,child,tlen,cStart);
  }

#ifdef STATS
  v->num_compares_for_tandem_arrays += 2;
#endif

} /* vocabulary_minimize() */

/*===========================================================================*/
/*
 * vocabulary_minimize_rec
 *
 * Find a minimal leftmost covering set of (primitive) tandem repeats.
 *
 * Parameters:  v      -  a vocabulary structure
 *              node   -  a suffix tree node
 *              depth  -  the string-depth of node
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void vocabulary_minimize_rec(VOCABULARY_STRUCT *v, STREE_NODE node, int depth)
{
  STREE_NODE child;
  int id, edgelen, tlen;

  for(child=stree_get_children(v->tree,node);
      child!=NULL;
      child=stree_get_next(v->tree,child)) {

    id = stree_get_ident(v->tree,child);
    edgelen = stree_get_edgelen(v->tree,child);

    if((tlen=v->tlens1[id]) > 0)
      vocabulary_minimize(v,node,depth,child,tlen,child);
    if((tlen=v->tlens2[id]) > 0)
      vocabulary_minimize(v,node,depth,child,tlen,child);

    vocabulary_minimize_rec(v,child,depth+edgelen);
  }
} /* vocabulary_minimize_rec() */

/*===========================================================================*/
/*
 * vocabulary_find_arrays
 *
 * Starting at the location of a tandem repeat from a minimal leftmost covering
 * set of primitive tandem repeats, walk down the tree and find all tandem
 * arrays as long as they exist. At the end, do a modified suffix-link walk
 * if the previous level was complete.
 *
 * Parameters:  v      -  a vocabulary structure
 *              node   -  the node above the current location
 *              depth  -  the string-depth of node
 *              child  -  the node below the current location
 *              tlen   -  length of the current tandem repeat
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void vocabulary_find_arrays(VOCABULARY_STRUCT *v, STREE_NODE node, int depth,
                            STREE_NODE child, int tlen)
{
  STREE_NODE cTmp;
  int period, offset, edgelen, num_inserts;
  char *edgestr;

  period = tlen / 2;
  offset = tlen - depth;
  edgelen = stree_get_edgelen(v->tree,child);
  edgestr = stree_get_edgestr(v->tree,child);

  /* mark rotated tandem repeats */
  num_inserts = 1;
  vocabulary_sl_walk(v,node,depth,child,tlen,period,&num_inserts);

  /* if necessary, find next edge */
  cTmp = child;
  if(offset == edgelen) {
    cTmp = stree_find_child(v->tree,child,*(edgestr+offset-period));
    if(cTmp != NULL) {
      node = child;
      depth += edgelen;
      offset -= edgelen;
      child = cTmp;
      edgelen = stree_get_edgelen(v->tree,child);
      edgestr = stree_get_edgestr(v->tree,child);
    }
  }

  /* do full rotations as long as possible */
  tlen += period;
  while(cTmp!=NULL && *(edgestr+offset)==*(edgestr+offset-period)) {

#ifdef STATS
    v->num_compares_for_tandem_arrays += 2;
#endif

    /* walk down (single steps) */
    offset++;

    if(depth+offset == tlen) {
      vocabulary_tloc_insert(v,child,-tlen);
      num_inserts = 1;
      vocabulary_mod_sl_walk(v,node,depth,child,tlen,period,offset,period,
                             &num_inserts);
      tlen += period;
    }

    if(offset == edgelen) {
      cTmp = stree_find_child(v->tree,child,*(edgestr+edgelen-period));
      if(cTmp != NULL) {
        node = child;
        depth += edgelen;
        offset -= edgelen;
        child = cTmp;
        edgelen = stree_get_edgelen(v->tree,child);
        edgestr = stree_get_edgestr(v->tree,child);
      }
    }
  }

  /* if necessary, do one additional rotation */
  if(num_inserts==period && depth+offset > tlen-period)
    vocabulary_mod_sl_walk(v,node,depth,child,tlen,period,offset,period,
                           &num_inserts);

} /* vocabulary_find_arrays() */

/*===========================================================================*/
/*
 * vocabulary_find_arrays_rec
 *
 * Given the locations of a minimal leftmost covering set of tandem repeats,
 * find the locations of all primitive tandem arrays.
 *
 * Parameters:  v      -  a vocabulary structure
 *              node   -  a suffix tree node
 *              depth  -  the string-depth of node
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void vocabulary_find_arrays_rec(VOCABULARY_STRUCT *v, STREE_NODE node,
                                int depth)
{
  STREE_NODE child;
  int id, edgelen, tlen;

  for(child=stree_get_children(v->tree,node);
      child!=NULL;
      child=stree_get_next(v->tree,child)) {

    id = stree_get_ident(v->tree,child);
    edgelen = stree_get_edgelen(v->tree,child);

    if((tlen=v->tlens1[id]) > 0)
      vocabulary_find_arrays(v,node,depth,child,tlen);
    if((tlen=v->tlens2[id]) > 0)
      vocabulary_find_arrays(v,node,depth,child,tlen);

    vocabulary_find_arrays_rec(v,child,depth+edgelen);
  }
} /* vocabulary_find_arrays_rec() */

/*===========================================================================*/
/*
 * vocabulary_find_tandem_repeats
 *
 * Find the tandem repeat decoration of the suffix tree.
 *
 * Parameters:  v  -  a vocabulary structure
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void vocabulary_find_tandem_repeats(VOCABULARY_STRUCT *v)
{
  int i, num_blocks, block_i, block_iplus1, block_iplus2, dummy_leaf_label;
  tandem *dummy_tlist;

  /*
   * Find a leftmost covering set using the Lempel-Ziv decomposition.
   */
  num_blocks = decomposition_get_num_blocks(v->decomposition);

  /* leftreps with i=0 is not necessary: can't contain a tandem repeat */
  for(i=1; i<num_blocks-1; i++) {
    block_i = decomposition_get_block(v->decomposition,i);
    block_iplus1 = decomposition_get_block(v->decomposition,i+1);
    block_iplus2 = decomposition_get_block(v->decomposition,i+2);
    leftreps(v,block_i,block_iplus1,block_iplus2);
    rightreps(v,0,block_i,block_iplus1,block_iplus2);
  }

  /* one extra rightreps for the last block */
  if(num_blocks > 1) {
    block_i = decomposition_get_block(v->decomposition,i);
    block_iplus1 = decomposition_get_block(v->decomposition,i+1);
    rightreps(v,0,block_i,block_iplus1,block_iplus1);
  }

  /* Find the locations of a (possibly smaller) leftmost covering set. */
  vocabulary_collect(v,stree_get_root(v->tree),0,
                     &dummy_leaf_label,&dummy_tlist);

  /* Find the locations of the right-rotated tandem repeats. */
  vocabulary_rotate(v,stree_get_root(v->tree),0);
  vocabulary_capitalize(v,stree_get_root(v->tree));

} /* vocabulary_find_tandem_repeats() */

/*===========================================================================*/
/*
 * vocabulary_find_primitive_tandem_repeats
 *
 * Given the tandem repeat decoration of the suffix tree, filter out the
 * non-primitive tandem repeats.
 *
 * Parameters:  d  -  a vocabulary structure
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void vocabulary_find_primitive_tandem_repeats(VOCABULARY_STRUCT *v)
{
  vocabulary_filter(v,stree_get_root(v->tree));

} /* vocabulary_find_primitive_tandem_repeats() */

/*===========================================================================*/
/*
 * vocabulary_find_tandem_arrays
 *
 * Given the (primitive) tandem repeat decoration of the suffix tree, find the
 * locations of all (primitive) tandem arrays.
 *
 * Parameters:  v      -  a vocabulary structure
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void vocabulary_find_tandem_arrays(VOCABULARY_STRUCT *v)
{
  /* Keep a minimal leftmost covering set. */
  vocabulary_minimize_rec(v,stree_get_root(v->tree),0);

  /* Find primitive tandem arrays. */
  vocabulary_find_arrays_rec(v,stree_get_root(v->tree),0);
  vocabulary_capitalize(v,stree_get_root(v->tree));

} /* vocabulary_find_tandem_arrays() */

/*===========================================================================*/
/*
 * vocabulary_write_repeat
 *
 * Write a single tandem repeat/tandem array.
 *
 * Parameters:  type        -  the type of the repeat
 *              raw_string  -  the raw string containing the repeat
 *              pos         -  the starting position of the repeat
 *              len         -  the length of the repeat
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void vocabulary_write_repeat(char *raw_string, int pos, int len, char *type)
{
  int i, textlen, restlen;
  char *s, *t, buffer[77];

  sprintf(buffer,"%s: ",type);
  buffer[76] = '\0';

  textlen = strlen(buffer);
  restlen = 76-textlen;
  for (s=&buffer[textlen],t=&raw_string[pos],i=0;
       i<restlen && i<len;
       s++,t++,i++)
    *s = isprint((int)(*t)) ? *t : '#';
  *s = '\0';
  mputs(buffer);
  if(len > restlen)
    mputs("...");
  mputc('\n');

} /* vocabulary_write_repeat() */

/*===========================================================================*/
/*
 * vocabulary_write_rec
 *
 * Write the vocabulary of repeats from the subtree rooted a node.
 *
 * Parameters:  v      -  a vocabulary structure
 *              node   -  a suffix tree node
 *              depth  -  the string-depth of node
 *              type   -  the type of the repeat
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void vocabulary_write_rec(VOCABULARY_STRUCT *v, STREE_NODE node, int depth,
                          char *type)
{
  STREE_NODE child;
  int id, edgelen, tlen;
  char *edgestr;

  for(child=stree_get_children(v->tree,node);
      child!=NULL;
      child=stree_get_next(v->tree,child)) {

    id = stree_get_ident(v->tree,child);
    edgelen = stree_get_edgelen(v->tree,child);
    edgestr = stree_get_edgestr(v->tree,child);

    if((tlen=v->tlens1[id]) > 0)
      vocabulary_write_repeat(v->raw_string,edgestr-depth-v->string,tlen,type);
    if((tlen=v->tlens2[id]) > 0)
      vocabulary_write_repeat(v->raw_string,edgestr-depth-v->string,tlen,type);

    vocabulary_write_rec(v,child,depth+edgelen,type);
  }

} /* vocabulary_write_rec() */

/*===========================================================================*/
/*
 * vocabulary_write
 *
 * Write the vocabulary of repeats from the decorated suffix tree.
 *
 * Parameters:  v     -  a vocabulary structure
 *              type  -  the type of the repeat
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void vocabulary_write(VOCABULARY_STRUCT *v, char *type)
{
  vocabulary_write_rec(v,stree_get_root(v->tree),0,type);

} /* vocabulary_write() */

/*===========================================================================*/
/*
 * vocabulary_count_leaves
 *
 * Count the number of leaves in the subtree below node.
 *
 * Parameters:  tree  -  a suffix tree structure
 *              node  -  a suffix tree node
 *
 * Returns:  the number of leaves in the subtree below node.
 */
/*---------------------------------------------------------------------------*/
int vocabulary_count_leaves(SUFFIX_TREE tree, STREE_NODE node)
{
  STREE_NODE child;
  int num;

  if(stree_get_num_children == 0)
    num = 1;
  else {
    num = stree_get_num_leaves(tree,node);
    for(child=stree_get_children(tree,node);
        child!=NULL;
        child=stree_get_next(tree,child))
      num += vocabulary_count_leaves(tree,child);
  }
  return num;

} /* vocabulary_count_leaves() */
/*===========================================================================*/
/*
 * vocabulary_count_rec
 *
 * Count size of the vocabulary and number of occurrences in subtree below
 * node.
 *
 * Parameters:  v     -  a vocabulary structure
 *              node  -  a suffix tree node
 *              num   -  return parameter for size of the vocabulary
 *              occ   -  return parameter for number of occurrences
 *
 * Returns:  size of vocabulary in variable num
 *           number of occurrences in variable occ.
 */
/*---------------------------------------------------------------------------*/
void vocabulary_count_rec(VOCABULARY_STRUCT *v, STREE_NODE node,
                          unsigned int *num, unsigned int *occ)
{
  STREE_NODE child;
  int id;

  for(child=stree_get_children(v->tree,node);
      child!=NULL;
      child=stree_get_next(v->tree,child)) {

    id = stree_get_ident(v->tree,child);

    if(v->tlens1[id] > 0) {
      (*num)++;
      *occ += vocabulary_count_leaves(v->tree,child);
    }
    if(v->tlens2[id] > 0) {
      (*num)++;
      *occ += vocabulary_count_leaves(v->tree,child);
    }

    vocabulary_count_rec(v,child,num,occ);
  }

} /* vocabulary_count_rec() */

/*===========================================================================*/
/*
 * vocabulary_count
 *
 * Count size of the vocabulary and number of occurrences.
 * (This is a naive, recursive implementation!)
 *
 * Parameters:  v    -  a vocabulary structure
 *              num  -  return parameter for size of the vocabulary
 *              occ  -  return parameter for number of occurrences
 *
 * Returns:  size of vocabulary in variable num
 *           number of occurrences in variable occ.
 */
/*---------------------------------------------------------------------------*/
void vocabulary_count(VOCABULARY_STRUCT *v,
                      unsigned int *num, unsigned int *occ)
{
  *num = 0;
  *occ = 0;
  vocabulary_count_rec(v,stree_get_root(v->tree),num,occ);

} /* vocabulary_count() */

/*===========================================================================*/

