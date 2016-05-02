
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
#include "stree_weiner.h"
#include "stree_lca.h"
#include "stree_decomposition.h"
#include "strmat_stubs2.h"


/*
 * strmat_ukkonen_build
 *
 * Performs the Ukkonen suffix tree construction algorithm for
 * one or more strings.
 *
 * Parameters:   strings          -  the input strings
 *               num_strings      -  the number of input strings
 *               build_policy     -  suffix tree build policy
 *               build_threshold  -  threshold used by LIST_THEN_ARRAY
 *               print_stats      -  flag telling whether to print the stats
 *               print_tree       -  flag telling whether to print the tree
 *
 * Returns:  non-zero on success, zero on error
 */
int strmat_ukkonen_build(STRING **strings, int num_strings, int build_policy,
                         int build_threshold, int print_stats, int print_tree)
{
  int i, max_length, total_length;
  SUFFIX_TREE tree;

  if (strings == NULL)
    return 0;

  /*
   * Find the total and maximum length of the input strings.
   */
  max_length = -1;
  total_length = 0;
  for (i=0; i < num_strings; i++) {
    total_length += strings[i]->length;
    if (max_length == -1 || strings[i]->length > max_length)
      max_length = strings[i]->length;
  }

  /*
   * Build the tree, then print the output.
   */
  tree = stree_gen_ukkonen_build(strings, num_strings, build_policy,
                                 build_threshold);
  if (tree == NULL)
    return 0;

  if (print_stats) {
    mprintf("\nStatistics:\n");
#ifdef STATS
    mprintf("   Sum of Sequence Sizes:       %d\n", total_length);
    mprintf("   Number of Tree Nodes:        %d\n", stree_get_num_nodes(tree));
    mprintf("   Size of Optimized Tree:      %d\n", tree->tree_size);
    mprintf("   Bytes per Character:         %.2f\n",
            (float) tree->tree_size / (float) total_length);
    mprintf("\n");
    mprintf("   Number of Comparisons:       %d\n", tree->num_compares);
    mprintf("   Cost of Constructing Edges:  %d\n", tree->creation_cost);
    mprintf("   Number of Edges Traversed:   %d\n", tree->edges_traversed);
    mprintf("   Cost of Traversing Edges:    %d\n", tree->child_cost);
    mprintf("   Number of Links Traversed:   %d\n", tree->links_traversed);
#else
    mprintf("   No statistics available.\n");
#endif

    if (mputc('\n') == 0) {
      stree_delete_tree(tree);
      return 1;
    }
  }

  if (print_tree) {
    mprintf("Suffix Tree:\n");
    if (max_length < 40)
      small_print_tree(tree, stree_get_root(tree), 0, (num_strings > 1));
    else
      large_print_tree(tree, stree_get_root(tree), (num_strings > 1));
    mputc('\n');
  }

  stree_delete_tree(tree);

  return 1;
}



/*
 * strmat_weiner_build
 *
 * Performs the Weiner suffix tree construction algorithm for
 * one or more strings.
 *
 * Parameters:   strings          -  the input strings
 *               num_strings      -  the number of input strings
 *               build_policy     -  suffix tree build policy
 *               build_threshold  -  threshold used by LIST_THEN_ARRAY
 *               print_stats      -  flag telling whether to print the stats
 *               print_tree       -  flag telling whether to print the tree
 *
 * Returns:  non-zero on success, zero on error
 */
int strmat_weiner_build(STRING **strings, int num_strings, int build_policy,
                        int build_threshold, int print_stats, int print_tree)
{
  int i, max_length, total_length, num_nodes, size;
  SUFFIX_TREE tree;

  if (strings == NULL)
    return 0;

  /*
   * Find the total and maximum length of the input strings.
   */
  max_length = -1;
  total_length = 0;
  for (i=0; i < num_strings; i++) {
    total_length += strings[i]->length;
    if (max_length == -1 || strings[i]->length > max_length)
      max_length = strings[i]->length;
  }

  /*
   * Build the tree, then print the output.
   */
  tree = stree_gen_weiner_build(strings, num_strings, build_policy,
                                build_threshold);
  if (tree == NULL)
    return 0;

  if (print_stats) {
    num_nodes = stree_get_num_nodes(tree);
    size = num_nodes * strings[0]->alpha_size;

    mprintf("\nStatistics:\n");
#ifdef STATS
    mprintf("   Sum of Sequence Sizes:        %d\n", total_length);
    mprintf("   Number of Tree Nodes:         %d\n", num_nodes);
    mprintf("   Size of Optimized Tree:       %d\n", tree->tree_size);
    mprintf("   Bytes per Character:          %.2f\n",
            (float) tree->tree_size / (float) total_length);
    mprintf("   Size of Unoptimized Vectors:  %d\n",
            size * sizeof(STREE_NODE) + size / 8);
    mprintf("\n");
    mprintf("   Number of Comparisons:        %d\n", tree->num_compares);
    mprintf("   Cost of Constructing Edges:   %d\n", tree->creation_cost);
    mprintf("   Number of Edges Traversed:    %d\n", tree->edges_traversed);
    mprintf("   Cost of Traversing Edges:     %d\n", tree->child_cost);
    mprintf("   Number of Links Traversed:    %d\n", tree->links_traversed);
#else
    mprintf("   No statistics available.\n");
#endif

    if (mputc('\n') == 0) {
      stree_delete_tree(tree);
      return 1;
    }
  }

  if (print_tree) {
    mprintf("Suffix Tree:\n");
    if (max_length < 40)
      small_print_tree(tree, stree_get_root(tree), 0, (num_strings > 1));
    else
      large_print_tree(tree, stree_get_root(tree), (num_strings > 1));
    mputc('\n');
  }

  stree_delete_tree(tree);

  return 1;
}


/*
 * strmat_stree_match
 *
 * Perform an exact matching of a pattern to a suffix tree of one or
 * more strings.
 *
 *
 * Parameters:   pattern          -  the input pattern
 *               strings          -  the input strings
 *               num_strings      -  the number of input strings
 *               build_policy     -  suffix tree build policy
 *               build_threshold  -  threshold used by LIST_THEN_ARRAY
 *               print_stats      -  flag telling whether to print the stats
 *
 * Returns:  non-zero on success, zero on error
 */
static MATCHES matchlist;
static int matchcount, matcherror, patlen;

static int add_match(SUFFIX_TREE tree, STREE_NODE node)
{
  int i, pos, id;
  char *seq;
  MATCHES newmatch;
  
  for (i=1; stree_get_leaf(tree, node, i, &seq, &pos, &id); i++) {
    newmatch = alloc_match();
    if (newmatch == NULL) {
      free_matches(matchlist);
      matchlist = NULL;
      matcherror = 1;
      return 0;
    }

    /*
     * Shift positions by 1 here (from 0..N-1 to 1..N).
     */
    newmatch->type = TEXT_SET_EXACT;
    newmatch->lend = pos + 1;
    newmatch->rend = pos + patlen;
    newmatch->textid = id;

    newmatch->next = matchlist;
    matchlist = newmatch;
    matchcount++;
  }

  return 1;
}


int strmat_stree_match(STRING *pattern, STRING **strings, int num_strings,
                       int build_policy, int build_threshold, int print_stats)
{
  int flag, pos, matchlen;
#ifdef STATS
  int num_compares, edges_traversed, child_cost;
#endif
  MATCHES back, current, next;
  STREE_NODE node;
  SUFFIX_TREE tree;

  if (pattern == NULL || strings == NULL)
    return 0;

  /*
   * Build the suffix tree.
   */
  mprintf("Building the tree...\n\n");
  tree = stree_gen_ukkonen_build(strings, num_strings, build_policy,
                                 build_threshold);
  if (tree == NULL)
    return 0;

  stree_reset_stats(tree);

  /*
   * Match the pattern string to a path in the suffix tree.
   */
  matchlen = stree_match(tree, pattern->sequence, pattern->length,
                         &node, &pos);
  if (matchlen < 0) {
    stree_delete_tree(tree);
    return 0;
  }

#ifdef STATS
  num_compares = tree->num_compares;
  edges_traversed = tree->edges_traversed;
  child_cost = tree->child_cost;

  stree_reset_stats(tree);
#endif

  /*
   * Traverse the subtree, finding the matches.
   */
  matchlist = NULL;
  matchcount = matcherror = 0;
  patlen = pattern->length;

  if (matchlen == pattern->length) {
    stree_traverse_subtree(tree, node, add_match, (int (*)()) NULL);
    if (matcherror) {
      stree_delete_tree(tree);
      return 0;
    }
      
    /*
     * Bubble sort the matches.
     */
    flag = 1;
    while (flag) {
      flag = 0;
      back = NULL;
      current = matchlist;
      while (current->next != NULL) {
        if (current->next->textid < current->textid ||
            (current->next->textid == current->textid &&
             current->next->lend < current->lend)) {
          /*
           * Move current->next before current in the list.
           */
          next = current->next;
          current->next = next->next;
          next->next = current;
          if (back == NULL)
            back = matchlist = next;
          else
            back = back->next = next;
          
          flag = 1;
        }
        else {
          back = current;
          current = current->next;
        }
      }
    }
  }

  /*
   * Print the matches and the statistics.
   */
  print_matches(NULL, strings, num_strings, matchlist, matchcount);

  if (print_stats) {
    mprintf("Statistics:\n");
#ifdef STATS
    mprintf("   Matching:\n");
    mprintf("      Pattern Length:          %d\n", pattern->length);
    mprintf("      Number of Comparisons:   %d\n", num_compares);
    mprintf("      Number Edges Traversed:  %d\n", edges_traversed);
    mprintf("      Cost of Edge Traversal:  %d\n", child_cost);
    mprintf("\n");
    mprintf("   Subtree Traversal:\n");
    mprintf("      Number of Matches:       %d\n", matchcount);
    mprintf("      Number Edges Traversed:  %d\n", tree->edges_traversed);
    mprintf("      Cost of Edge Traversal:  %d\n", tree->child_cost);
#else
    mprintf("   No statistics available.\n");
#endif
    mputc('\n');
  }

  /*
   * Free everything allocated.
   */
  free_matches(matchlist);
  stree_delete_tree(tree);

  return 1;
}


/*
 * strmat_stree_lca
 *
 * Performs the constant time LCA algorithm on a suffix tree for a
 * string.
 *
 * Parameters:   strings          -  the input strings
 *               num_strings      -  the number of input strings
 *               build_policy     -  suffix tree build policy
 *               build_threshold  -  threshold used by LIST_THEN_ARRAY
 *               print_stats      -  flag telling whether to print the stats
 *
 * Returns:  non-zero on success, zero on error
 */
static void compute_nodemap(SUFFIX_TREE tree, STREE_NODE node,
                            STREE_NODE *nodemap);
int int_strmat_stree_lca(STRING **strings, int num_strings, int build_policy,
                         int build_threshold, int print_stats, LCA_TYPE type);

int strmat_stree_lca(STRING **strings, int num_strings, int build_policy,
                     int build_threshold, int print_stats)
{  return int_strmat_stree_lca(strings, num_strings, build_policy,
                               build_threshold, print_stats, LCA_LINEAR);  }
int strmat_stree_naive_lca(STRING **strings, int num_strings, int build_policy,
                           int build_threshold, int print_stats)
{  return int_strmat_stree_lca(strings, num_strings, build_policy,
                               build_threshold, print_stats, LCA_NAIVE);  }

int int_strmat_stree_lca(STRING **strings, int num_strings, int build_policy,
                         int build_threshold, int print_stats, LCA_TYPE type)
{
  int i, num_nodes, num1, num2, len, num_lcas, max_length;
  char *s, *line, buffer[64];
  STREE_NODE x, y, z, *nodemap;
  SUFFIX_TREE tree;
  LCA_STRUCT *lcastruct;

  if (strings == NULL)
    return 0;

  /*
   * Build the tree.
   */
  printf("Building the suffix tree...\n");
  tree = stree_gen_ukkonen_build(strings, num_strings, build_policy,
                                 build_threshold);
  if (tree == NULL)
    return 0;

  num_nodes = stree_get_num_nodes(tree);
  max_length = -1;
  for (i=0; i < num_strings; i++)
    if (max_length == -1 || strings[i]->length > max_length)
      max_length = strings[i]->length;

  /*
   * Preprocess the suffix tree.
   */
  printf("Preprocessing...\n");
  lcastruct = NULL;
  switch (type) {
  case LCA_NAIVE:   lcastruct = lca_naive_prep(tree);  break;
  case LCA_LINEAR:  lcastruct = lca_prep(tree);  break;
  case LCA_NLOGN:   return 0;
  }
  if (lcastruct == NULL) {
    stree_delete_tree(tree);
    return 0;
  }

  /*
   * Build the map of suffix tree nodes.
   */
  nodemap = malloc(num_nodes * sizeof(STREE_NODE));
  if (nodemap == NULL) {
    lca_free(lcastruct);
    stree_delete_tree(tree);
    return 0;
  }

  compute_nodemap(tree, stree_get_root(tree), nodemap);

  /*
   * Query the user to enter nodes, and compute the LCA of those nodes.
   */
  printf("\n");
  printf("Commands (0-%d 0-%d - Find LCA of two nodes (identify by number),\n",
         num_nodes-1, num_nodes-1);
  printf("          ! - print suffix tree, Ctl-D - quit)\n");

  num_lcas = 0;
  while (1) {
    printf("Enter nodes: ");

    if ((line = my_getline(stdin, NULL)) == NULL) {
      printf("\n\n");
      break;
    }

    if (line[0] == '\0')
      continue;

    if (line[0] == '!') {
      mstart(stdin, stdout, OK, OK, 5, NULL);
      mputc('\n');
      mprintf("Suffix Tree:\n");
      if (max_length < 40)
        small_print_tree(tree, stree_get_root(tree), 0, (num_strings > 1));
      else
        large_print_tree(tree, stree_get_root(tree), (num_strings > 1));
      mputc('\n');
      mend(2);
    }
    else if (sscanf(line, "%d %d", &num1, &num2) == 2 && num1 >= 0 &&
             num1 < num_nodes && num2 >= 0 && num2 < num_nodes) {
      x = nodemap[num1];
      y = nodemap[num2];

      z = NULL;
      switch (type) {
      case LCA_NAIVE:   z = lca_naive_lookup(lcastruct, x, y);  break;
      case LCA_LINEAR:  z = lca_lookup(lcastruct, x, y);  break;
      case LCA_NLOGN:  return 0;
      }
      num_lcas++;

      if (x == stree_get_root(tree))
        printf("   Node %d:  (root)\n", stree_get_ident(tree, x));
      else {
        len = stree_get_labellen(tree, x);
        stree_get_label(tree, x, buffer, 50, 0);
        for (s=buffer,i=0; *s && i < 50; s++,i++)
          if (!isprint((int)(*s)))
            *s = '#';
        if (len > 50) {
          buffer[50] = buffer[51] = buffer[52] = '.';
          buffer[53] = '\0';
        }
        else if (len == 50)
          buffer[50] = '\0';

        printf("   Node %d:  %s\n", stree_get_ident(tree, x), buffer);
      }

      if (y == stree_get_root(tree))
        printf("   Node %d:  (root)\n", stree_get_ident(tree, y));
      else {
        len = stree_get_labellen(tree, y);
        stree_get_label(tree, y, buffer, 50, 0);
        for (s=buffer,i=0; *s && i < 50; s++,i++)
          if (!isprint((int)(*s)))
            *s = '#';
        if (len > 50) {
          buffer[50] = buffer[51] = buffer[52] = '.';
          buffer[53] = '\0';
        }
        else if (len == 50)
          buffer[50] = '\0';

        printf("   Node %d:  %s\n", stree_get_ident(tree, y), buffer);
      }

      if (z == stree_get_root(tree))
        printf("   LCA Node %d:  (root)\n", stree_get_ident(tree, z));
      else {
        len = stree_get_labellen(tree, z);
        stree_get_label(tree, z, buffer, 50, 0);
        for (s=buffer,i=0; *s && i < 50; s++,i++)
          if (!isprint((int)(*s)))
            *s = '#';
        if (len > 50) {
          buffer[50] = buffer[51] = buffer[52] = '.';
          buffer[53] = '\0';
        }
        else if (len == 50)
          buffer[50] = '\0';

        printf("   LCA Node %d:  %s\n", stree_get_ident(tree, z), buffer);
      }

      putchar('\n');
    }
    else
      printf("  Invalid input line.  Please reenter.\n\n");
  }

  if (print_stats) {
    printf("\nStatistics:\n");
#ifdef STATS
    printf("   Preprocessing Steps:    %d\n", lcastruct->num_prep);
    printf("\n");
    printf("   Number LCA's Computed:  %d\n", num_lcas);
    printf("   LCA Compute Steps:      %d\n", lcastruct->num_compares);
#else
    printf("   No statistics available.\n");
#endif

    putchar('\n');
  }

  free(nodemap);
  switch (type) {
  case LCA_NAIVE:   lca_naive_free(lcastruct);  break;
  case LCA_LINEAR:  lca_free(lcastruct);  break;
  case LCA_NLOGN:  return 0;
  }
  stree_delete_tree(tree);

  return 1;
}


/*
 * compute_nodemap
 *
 * Compute the mapping from identifiers to suffix tree nodes and store
 * that mapping in "nodemap"
 *
 * Parameters:  tree     -  a suffix tree
 *              node     -  a suffix tree node
 *              nodemap  -  the nodemap being computed
 *
 * Returns: nothing
 */
static void compute_nodemap(SUFFIX_TREE tree, STREE_NODE node, STREE_NODE *map)
{
  STREE_NODE child;

  map[stree_get_ident(tree, node)] = node;

  child = stree_get_children(tree, node);
  while (child != NULL) {
    compute_nodemap(tree, child, map);
    child = stree_get_next(tree, child);
  }
}


/*
 * strmat_stree_walkaround
 *
 * Interactively walk around the nodes of a suffix tree.
 *
 * Parameters:   strings          -  the input strings
 *               num_strings      -  the number of input strings
 *               build_policy     -  suffix tree build policy
 *               build_threshold  -  threshold used by LIST_THEN_ARRAY
 *
 * Returns:  non-zero on success, zero on error
 */
static void print_stree_node(SUFFIX_TREE tree, STREE_NODE node,
                             int gen_stree_flag, int mend_num_lines);

int strmat_stree_walkaround(STRING **strings, int num_strings,
                            int build_policy, int build_threshold)
{
  int alphabet;
  char mapch, *choice;
  STREE_NODE node, child;
  SUFFIX_TREE tree;

  if (strings == NULL)
    return 0;

  /*
   * Build the tree.
   */
  printf("Building the suffix tree...\n");
  tree = stree_gen_ukkonen_build(strings, num_strings, build_policy,
                                 build_threshold);
  if (tree == NULL)
    return 0;

  alphabet = strings[0]->alphabet;

  /*
   * Main interactive loop.
   */
  node = stree_get_root(tree);
  while (1) {
    /*
     * Print the details of the current node.
     */
    printf("\n\n");
    print_stree_node(tree, node, (num_strings > 1), 5);

    /*
     * Ask the user where to move in the tree.
     */
    printf("\n");
    printf("Commands (d%% - move down to a child, u - move up to parent,\n");
    printf("          l - move across suffix link, Ctl-D - quit)\n");
    printf("Enter Move: ");

    if ((choice = my_getline(stdin, NULL)) == NULL)
      break;

    if (choice[0] == '\0')
      continue;

    /*
     * Execute the move command.
     */
    switch (toupper(choice[0])) {
    case 'D':
      if (choice[1] == '\0') {
        printf("\nYou must specify the first character "
               "on an edge to a child.\n");
        continue;
      }

      mapch = mapchar(alphabet, choice[1]);

      if ((child = stree_find_child(tree, node, mapch)) == NULL) {
        printf("\nNo child's edge begins with '%c'.\n", choice[1]);
        continue;
      }

      node = child;
      break;

    case 'U':
      if (node == stree_get_root(tree))
        printf("\nYou cannot move up from the root.\n");
      else
        node = stree_get_parent(tree, node);
      break;

    case 'L':
      if (node == stree_get_root(tree))
        printf("\nThe root has no suffix link.\n");
      else
        node = stree_get_suffix_link(tree, node);
      break;

    default:
      printf("\nThat is not a choice.\n");
    }
  }

  stree_delete_tree(tree);

  return 1;
}


/*
 * print_stree_node
 *
 * Print the details about a node of a suffix tree.
 *
 * Parameters:  tree            -  a suffix tree
 *              node            -  a tree node
 *              gen_stree_flag  -  is this a generalized suffix tree?
 *              mend_num_lines  -  how many lines will appear after this
 *
 * Returns: nothing
 */
static void print_stree_node(SUFFIX_TREE tree, STREE_NODE node,
                             int gen_stree_flag, int mend_num_lines)
{
  int i, j, index, ident, idwidth, labellen, edgelen;
  int leafnum, pos;
  char *edgestr, label[36];
  char *str;
  STREE_NODE child;

  mstart(stdin, stdout, OK, OK, 5, NULL);

  /*
   * Get the node and edge information.
   */
  labellen = stree_get_labellen(tree, node);
  stree_get_label(tree, node, label, 30, 1);
  for (i=0; i < 30 && label[i]; i++)
    if (!isprint((int)label[i]))
      label[i] = '#';
  label[i] = '\0';

  edgelen = stree_get_edgelen(tree, node);

  ident = stree_get_ident(tree, node);
  idwidth = my_itoalen(ident);

  /*
   * Print the node info: ident, label, edge label, leaves & suffix link.
   */
  if (node == stree_get_root(tree))
    mprintf("Current node is Node %d, the Root\n", ident);
  else {
    mprintf("Current node is Node %d, labeled `%s%s'\n", ident,
            (labellen > 30 ? "..." : ""), label);
    mprintf("     Leaves:  ");
    leafnum = 1;
    for (;stree_get_leaf(tree, node, leafnum, &str, &pos, &index); leafnum++) {
      if (leafnum > 1)
        mputs(", ");

      if (gen_stree_flag)
        mprintf("%d:%d", index, pos + 1);
      else
        mprintf("%d", pos + 1);
    }
    if (leafnum == 1)
      mputs("(none)");
    mputc('\n');
    mprintf("       Edge:  %s%s\n", (edgelen > 30 ? "..." : ""),
            (edgelen > 30 ? label : (label + (labellen - edgelen))));
    mprintf("     Parent:  Node %d\n",
            stree_get_ident(tree, stree_get_parent(tree, node)));
    mprintf("  Suf. Link:  Node %d\n",
            stree_get_ident(tree, stree_get_suffix_link(tree, node)));
  }

  /*
   * Print the outgoing edges.
   */
  if (stree_get_num_children(tree, node) == 0)
    mprintf("   Children:\n       (none)\n");
  else {
    mprintf("   Children:\n");
    child = stree_get_children(tree, node);
    while (child != NULL) {
      edgestr = stree_get_rawedgestr(tree, child);
      edgelen = stree_get_edgelen(tree, child);

      for (j=0; j < 30 && j < edgelen; j++)
        label[j] = (isprint((int)edgestr[j]) ? edgestr[j] : '#');
      label[j] = '\0';

      if (edgelen > 30) {
        label[30] = label[31] = label[32] = '.';
        label[33] = '\0';
      }

      mprintf("       %s  ->  Node %d", label,
              stree_get_ident(tree, child));

      leafnum = 1;
      while (stree_get_leaf(tree, child, leafnum, &str, &pos, &index)) {
        if (leafnum == 1)
          mprintf("    (Leaf #");
        else
          mprintf(", ");

        if (gen_stree_flag)
          mprintf("%d:%d", index, pos + 1);
        else
          mprintf("%d", pos + 1);

        leafnum++;
      }
      if (leafnum > 1)
        mputc(')');

      mputc('\n');

      child = stree_get_next(tree, child);
    }
  }

  mend(mend_num_lines);
}


/*
 * strmat_stree_lempel_ziv
 *
 * Compute Lempel-Ziv decomposition of a string.
 *
 * Parameters:   string           -  the input string
 *               build_policy     -  suffix tree build policy
 *               build_threshold  -  threshold used by LIST_THEN_ARRAY
 *               print_stats      -  flag telling whether to print the stats
 *               mode             -  flag telling whether
 *                                     A) most efficient decomposition
 *                                     B) original Lempel-Ziv
 *                                     C) non-overlapping Lempel-Ziv
 *
 * Returns:  non-zero on success, zero on error
 */
int strmat_stree_lempel_ziv(STRING *string, int build_policy,
                            int build_threshold, int print_stats,
                            char mode)
{
  SUFFIX_TREE tree;
  DECOMPOSITION_STRUCT *decomp_struct;

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
  decomp_struct = decomposition_prep(tree, string->sequence,
                                     string->raw_seq, string->length);
  if(decomp_struct == NULL) {
    stree_delete_tree(tree);
    return 0;
  }

  /*
   * Compute the decomposition
   */
  mprintf("Computing the decomposition...\n");
  if(mode == 'A')
    lempel_ziv(decomp_struct);
  else if(mode == 'B')
    lempel_ziv_nonoverlapping(decomp_struct);
  else {
    fprintf(stderr,"Unknown mode `%d'.\n",mode);
    return 0;
  }

  /*
   * Print Lempel-Ziv decomposition
   */
  mprintf("\nThe %sLempel-Ziv decomposition is:\n\n",
          mode=='A' ? "" : mode=='B' ? "original " : "non-overlapping ");
  decomposition_print(decomp_struct);
  mprintf("\n");
  mend(14);

  /* Write summary of results */
  mstart(stdin,stdout,OK,OK,0,NULL);
  mprintf("\nSummary:\n");
  mprintf("   Number of Blocks:              %d\n",decomp_struct->num_blocks);
  mprintf("   Maximal Block Length:          %d\n",
          decomp_struct->max_block_length);
  mprintf("   Average Block Length:          %.1f\n",
          (float)decomp_struct->length / decomp_struct->num_blocks);

  /*
   * Write statistics and free memory
   */
  if (print_stats) {
    mprintf("\nStatistics:\n");
#ifdef STATS
    mprintf("   String Length:                 %d\n",decomp_struct->length);
    mprintf("   Suffix Tree\n");
    mprintf("      Number of Tree Nodes:       %d\n",
            stree_get_num_nodes(tree));
    mprintf("      Number of Compares:         %d\n",tree->num_compares);
    mprintf("   Decomposition\n");
    mprintf("      Number of Compares:         %d\n",
            decomp_struct->num_compares);
    mprintf("      Number of Edge Traversals:  %d\n",
            decomp_struct->num_edge_traversals);
#else
    mprintf("   No statistics available.\n");
#endif
    mprintf("\n");
  }

  decomposition_free(decomp_struct);
  stree_delete_tree(tree);

  return 1;
}



