/*
 * repeats_linear_occs.c
 *
 * Find all occurrences of tandem repeats in O(n + z) time.
 *
 * NOTES:
 *    7/98  -  Original implementation of the algorithm (Jens Stoye)
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
#include "repeats_vocabulary.h"
#include "repeats_linear_occs.h"

#define MAX(X,Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))
#define MIN3(X,Y,Z) (MIN(X,MIN(Y,Z)))
#define ABS(X) (((X) < 0) ? -(X) : (X))

/*===========================================================================*/
/*
 * linear_occs_prep
 *
 * Preprocessing for the linear tandem repeat occurrences algorithm.
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
VOCABULARY_STRUCT *linear_occs_prep(SUFFIX_TREE tree,
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

  /* for each position of the text */
  if((vocabulary->tandems = malloc(length * sizeof(tandem*))) == NULL) {
    linear_occs_free(vocabulary);
    return NULL;
  }
  memset(vocabulary->tandems, 0, length * sizeof(tandem*));

  if((vocabulary->last = malloc(length * sizeof(tandem*))) == NULL) {
    linear_occs_free(vocabulary);
    return NULL;
  }
  memset(vocabulary->last, 0, length * sizeof(tandem*));

  for(i=0; i<length; i++)
    vocabulary->tandems[i] = vocabulary->last[i] = NULL;

  /* for the longest two blocks of the decomposition */
  if((vocabulary->PREF = malloc(2*max_block_length * sizeof(int))) == NULL) {
    linear_occs_free(vocabulary);
    return NULL;
  }
  memset(vocabulary->PREF, 0, 2*max_block_length * sizeof(int));

  if((vocabulary->PREF2 = malloc(2*max_block_length * sizeof(int))) == NULL) {
    linear_occs_free(vocabulary);
    return NULL;
  }
  memset(vocabulary->PREF2, 0, 2*max_block_length * sizeof(int));

  if((vocabulary->SUFF = malloc(2*max_block_length * sizeof(int))) == NULL) {
    linear_occs_free(vocabulary);
    return NULL;
  }
  memset(vocabulary->SUFF, 0, 2*max_block_length * sizeof(int));

  vocabulary->num_tandem_repeat_occs = 0;

#ifdef STATS
  vocabulary->num_prep = tree->num_compares + decomposition->num_compares;
  vocabulary->num_compares_for_tandem_repeats = 0;
#endif

  return vocabulary;

} /* linear_occs_prep() */

/*===========================================================================*/
/*
 * linear_occs_free
 *
 * Free memory of a vocabulary structure (the part used by the linear
 * occurrences algorithm).
 *
 * Parameters:  vocabulary  -  a vocabulary structure.
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void linear_occs_free(VOCABULARY_STRUCT *vocabulary)
{
  int i;
  tandem *t, *tTmp;

  if(vocabulary->tandems!= NULL) {
    for(i=0; i<vocabulary->length; i++)
      for(t=vocabulary->tandems[i]; t!=NULL; t=tTmp) {
        tTmp = t->next;
        free(t);
      }
    free(vocabulary->tandems);
  }
  if(vocabulary->last != NULL)
    free(vocabulary->last);
  if(vocabulary->PREF != NULL)
    free(vocabulary->PREF);
  if(vocabulary->PREF2 != NULL)
    free(vocabulary->PREF2);
  if(vocabulary->SUFF != NULL)
    free(vocabulary->SUFF);
  free(vocabulary);

} /* linear_occs_free() */

/*===========================================================================*/
/*
 * linear_occs_tandem_append
 *
 * Append tandem repeat of length len at the end of tandem list no. pos.
 * Allocate memory locally here.
 *
 * Parameters:  v    -  a vocabulary structure
 *              pos  -  starting position of the tandem repeat
 *              len  -  length of the tandem repeat
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void linear_occs_tandem_append(VOCABULARY_STRUCT *v, int pos, int len)
{
  tandem *new_t;

  if((new_t = malloc(sizeof(tandem))) == NULL) {
    fprintf(stderr,"Can't allocate memory for tandem repeats.\n");
    return;
  }

  new_t->len = len;
  new_t->next = NULL;
  if(v->tandems[pos] == NULL)
    v->tandems[pos] = new_t;
  else
    v->last[pos]->next = new_t;
  v->last[pos] = new_t;

  v->num_tandem_repeat_occs++;

} /* linear_occs_tandem_append() */

/*===========================================================================*/
/*
 * lo_leftreps
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
void lo_leftreps(VOCABULARY_STRUCT *voc, int pos1, int pos2, int pos3)
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
   * Find the starting positions p of the tandem repeats:
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
    for(p=MAX(pos2-PREF[l]-l,pos2-2*l+1);
        p<=MIN(pos2+SUFF[ulen-l]-2*l,pos2-l-1);
        p++)
      linear_occs_tandem_append(voc,p,2*l);

} /* lo_leftreps() */

/*===========================================================================*/
/*
 * lo_rightreps
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
void lo_rightreps(VOCABULARY_STRUCT *voc, int pos1, int pos2, int posM,
                  int pos3)
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
        for(j=0; r+j<vlen && beta+j<ulen && v[vlast-r-j]==u[ulast-beta-j]; j++)
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
   * Find the starting positions p of the tandem repeat:
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
    for(p=MAX(pos2-SUFF[vlen-l],pos2-l);
        p<=MIN3(pos2+PREF[l]-l,pos2-1,posM-l-1);
        p++)
      linear_occs_tandem_append(voc,p,2*l);

  /* test additionally for l == vlen */
  for(p=pos2-SUFF[0]; p<=MIN(pos2-vlen,posM-vlen-1); p++)
    linear_occs_tandem_append(voc,p,2*vlen);

} /* lo_rightreps() */

/*===========================================================================*/
/*
 * vocabulary_tandem_copy
 *
 * Copy tandem repeats from source to destination as long as len <= maxlen.
 *
 * Parameters:  v    -  a vocabulary structure
 *              pos  -  starting position of the tandem repeat
 *              len  -  length of the tandem repeat
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void vocabulary_tandem_copy(VOCABULARY_STRUCT *v, int source, int target,
                            int maxlen)
{
  tandem *tmp_tandems, *tmp_last, *copy_t;

  tmp_tandems = v->tandems[target];
  tmp_last = v->last[target];
  v->tandems[target] = v->last[target] = NULL;

  for(copy_t=v->tandems[source];
      copy_t!=NULL && copy_t->len<=maxlen;
      copy_t=copy_t->next)
    linear_occs_tandem_append(v,target,copy_t->len);

  if(v->last[target] == NULL)
    v->tandems[target] = tmp_tandems;
  else
    v->last[target]->next = tmp_tandems;
  if(tmp_last != NULL)
    v->last[target] = tmp_last;

} /* vocabulary_tandem_copy() */

/*===========================================================================*/
/*
 * linear_occs_find_tandem_repeats
 *
 * Find all z occurrences of tandem repeats in O(n+z) time.
 *
 * Parameters:  v  -  a vocabulary structure
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void linear_occs_find_tandem_repeats(VOCABULARY_STRUCT *v)
{
  int i, p, num_blocks, block_i, block_iplus1, block_iplus2, prev_i;

  num_blocks = decomposition_get_num_blocks(v->decomposition);

  /* leftreps with i=0 is not necessary: can't contain a tandem repeat. */
  for(i=1; i<num_blocks-1; i++) {
    block_i = decomposition_get_block(v->decomposition,i);
    block_iplus1 = decomposition_get_block(v->decomposition,i+1);
    block_iplus2 = decomposition_get_block(v->decomposition,i+2);
    prev_i = decomposition_get_prev(v->decomposition,i);
    lo_leftreps(v,block_i,block_iplus1,block_iplus2);
    lo_rightreps(v,0,block_i,block_iplus1,block_iplus2);

    /* blocks with prev[i] undefined are skipped automatically since all blocks
     * of width 1 are ignored. */
    for(p=0; p<block_iplus1-block_i-1; p++)
      vocabulary_tandem_copy(v,prev_i+p,block_i+p,block_iplus1-block_i-p);
  }

  /* one extra rightreps for the last block */
  if(num_blocks > 1) {
    block_i = decomposition_get_block(v->decomposition,i);
    block_iplus1 = decomposition_get_block(v->decomposition,i+1);
    prev_i = decomposition_get_prev(v->decomposition,i);
    lo_rightreps(v,0,block_i,block_iplus1,block_iplus1);
    for(p=0; p<block_iplus1-block_i-1; p++)
      vocabulary_tandem_copy(v,prev_i+p,block_i+p,block_iplus1-block_i-p);
  }

} /* linear_occs_compute() */

/*===========================================================================*/
/*
 * linear_occs_write_repeat
 *
 * Write a tandem repeat.
 *
 * Parameters:  raw_string  -  the raw string
 *              pos         -  starting position of the repeat
 *              len         -  length of the repeat
 *              type        -  repeat type
 *
 * Returns:  nothing
 */
/*---------------------------------------------------------------------------*/
void linear_occs_write_repeat(char *raw_string, int pos, int len, char *type)
{
  int i, textlen, restlen;
  char *s, *t, buffer[77];

  sprintf(buffer,"%s: (%d,%d,%d) ",type,pos,len/2,2);
  buffer[76] = '\0';

  textlen = strlen(buffer);
  restlen = 76-textlen;
  for (s=&buffer[textlen],t=&raw_string[pos],i=0;
       i<restlen && i<len; i++,s++,t++)
    *s = (isprint((int)(*t)) ? *t : '#');
  *s = '\0';
  mputs(buffer);
  if(len > restlen)
    mputs("...");
  mputc('\n');

} /* linear_occs_write_repeat() */

/*===========================================================================*/
/* -- LINEAR_OCCS_WRITE --
 * Write the occurrences contained in the tandem repeat lists. */
/*---------------------------------------------------------------------------*/
void linear_occs_write(VOCABULARY_STRUCT *v, char *type)
{
  int i;
  tandem *t;

  for(i=0; i<v->length; i++)
    for(t=v->tandems[i]; t!=NULL; t=t->next)
      linear_occs_write_repeat(v->raw_string,i,t->len,type);

} /* linear_occs_write_rec() */

/*===========================================================================*/

