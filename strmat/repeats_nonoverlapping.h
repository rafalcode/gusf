
#ifndef _REPEATS_NONOVERLAPPING_H_
#define _REPEATS_NONOVERLAPPING_H_

#include <limits.h>

typedef struct NONOVERLAPPING_ENTRY {
  struct NONOVERLAPPING_ENTRY *next, *prev;
  struct NONOVERLAPPING_LIST *inList;
} no_entry;

typedef struct NONOVERLAPPING_LIST {
  struct NONOVERLAPPING_LIST *next, *prev;
  struct NONOVERLAPPING_NODE *atNode;
  no_entry **entries, **last;
  int len;
} no_list;

typedef struct NONOVERLAPPING_NODE {
  struct NONOVERLAPPING_NODE *next, *prev;
  struct NONOVERLAPPING_STRUCT *inNodelist;
  no_list *lists, *last, *last_source_list;
} no_node;

typedef struct NONOVERLAPPING_STRUCT {
  char *string_space, *string, *raw_string;
  int length;
  char a[CHAR_MAX+1];
  int alpha_size;

  no_entry *entries,*entries2,**entryList,**entryList2,**lastList,**lastList2;
  no_list *lists, *lists2;
  int nextList, nextList2;
  no_node *nodelist, *nodelist2, *nodes, *last, *nodes2, *last2;
  int nextNode;

  unsigned int num_nonoverlapping_maximal_pairs;

#ifdef STATS
  unsigned int num_compares;
#endif

} nonoverlapping_struct;

nonoverlapping_struct *nonoverlapping_prep(char *string, char *raw_string,
                                           int length);
void nonoverlapping_free(nonoverlapping_struct *nonoverlapping);
void nonoverlapping_find(nonoverlapping_struct *n);

#endif
