
#ifndef _REPEATS_PRIMITIVES_H_
#define _REPEATS_PRIMITIVES_H_

typedef struct PRIMITIVES_ENTRY {
  struct PRIMITIVES_ENTRY *next, *prev;
  struct PRIMITIVES_LIST *inList;
} pr_entry;

typedef struct PRIMITIVES_LIST {
  struct PRIMITIVES_LIST *next, *prev;
  struct PRIMITIVES_NODE *atNode;
  pr_entry *entries, *last;
  int len;
} pr_list;

typedef struct PRIMITIVES_NODE {
  struct PRIMITIVES_NODE *next, *prev;
  struct PRIMITIVES_STRUCT *inNodelist;
  pr_list *lists, *last, *last_source_list;
} pr_node;

typedef struct PRIMITIVES_STRUCT {
  char *string, *raw_string;
  int length;

  pr_entry *entries, *entries2;
  pr_list *lists, *lists2;
  int nextList, nextList2;
  pr_node *nodelist, *nodelist2, *nodes, *last, *nodes2, *last2;
  int nextNode;

  unsigned int num_primitive_tandem_repeat_occs;

#ifdef STATS
  unsigned int num_compares;
#endif

} primitives_struct;

primitives_struct *primitives_prep(char *string, char *raw_string, int length);
void primitives_free(primitives_struct *primitives);
void primitives_find(primitives_struct *p);

#endif
