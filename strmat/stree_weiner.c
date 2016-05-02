/*
 * stree_ukkonen.c
 *
 * The implementation of Ukkonen's suffix tree algorithm, for use with
 * strmat's suffix tree implementation.
 *
 * NOTES:
 *    7/94  -  Initial implementation of Weiner's algorithm (Joyce Lau)   
 *    9/95  -  Reimplemented suffix trees, optimized the data structure,
 *             created streeopt.[ch]   (James Knight)
 *    4/96  -  Modularized the code  (James Knight)
 *    7/96  -  Finished the modularization  (James Knight)
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include "strmat.h"
#include "stree_strmat.h"
#include "stree_weiner.h"


/*
 *
 * Forward references.
 *
 */
static int reset_links();
static int grow_links(int newid);
static int I(int id, char ch);
static STREE_NODE L(int id, char ch);
static int set_I(int id, char ch);
static int set_L(int id, char ch, STREE_NODE linkptr);
static int copy_links(int dest, int src);
static void free_links(void);


/*
 * stree_weiner_add_string
 *
 * Use Weiner's suffix tree construction algorithm to add a string
 * to a suffix tree.
 *
 * Parameters:  tree  -  a suffix tree
 *              S     -  the string to add
 *              Sraw  -  the raw version of the string
 *              M     -  the string length
 *              strid -  the string identifier
 *
 * Returns:  non-zero on success, zero on error.
 */
int stree_weiner_add_string(SUFFIX_TREE tree, char *S, char *Sraw,
                            int M, int strid)
{
  int i, t, l, id, headlen, edgelen, edgepos;
  char c;
  STREE_NODE root, node, previous, v, vprime, vdblprime, w;
  STREE_LEAF leaf;

  id = int_stree_insert_string(tree, S, Sraw, M, strid);
  if (id == -1)
    return 0;

  root = stree_get_root(tree);

  /*
   * Add the suffix consisting of only the last character of S to the tree.
   */
  i = M - 1;
  node = stree_find_child(tree, root, S[i]);
  if (node == NULL) {
    if ((leaf = int_stree_new_leaf(tree, id, i, i)) == NULL)
      return 0;

    if (int_stree_connect(tree, root, (STREE_NODE) leaf) == NULL) {
      int_stree_free_leaf(tree, leaf);
      return 0;
    }

    leaf->id = tree->num_nodes;
    tree->num_nodes++;
    node = (STREE_NODE) leaf;

    if (!set_I(root->id, S[i]))
      return 0;
  }
  else {
    if (stree_get_edgelen(tree, node) > 1) {
      if ((w = int_stree_edge_split(tree, node, 1)) == NULL)
        return 0;
      w->id = tree->num_nodes - 1;
      if (!copy_links(w->id, node->id))
        return 0;

      node = w;
    }

    if (int_stree_isaleaf(tree, node) &&
        (node = int_stree_convert_leafnode(tree, node)) == NULL)
      return 0;
    
    if (!int_stree_add_intleaf(tree, node, id, i))
      return 0;

    if (!set_L(root->id, S[i], node))
      return 0;
    if (node->suffix_link == NULL)
      node->suffix_link = root;
  }

  /*
   * Run Weiner's suffix tree construction algorithm.
   *
   * This implementation differs from the book description in 
   * several ways:
   *    1) All of the offsets into the sequence, and the phases of
   *       the algorithm, use the C array indices of 0 to M-1, not 1 to M.
   *    2) The constructed suffix tree only has suffix links in
   *       the internal nodes of the tree (to save space).  However,
   *       the stree_get_suffix_link function will return the suffix links
   *       even for leaves (it computes the leaves' suffix links on the
   *       fly).
   *    3) The Indicator and Link vectors are actually implemented using
   *       a link list data structure, rather than the bitmap/pointer-array
   *       structure described in the book.  The original implement used
   *       the book's data structures and Weiner's algorithm ran out of
   *       memory for any decently sized problem.  The tree size statistics
   *       given with the running of this program, however, will report
   *       the space used as if it were using the book's data structures.
   */
  for (i--; i >= 0; i--) {
    v = NULL;

    /*
     * Step 1.
     */
    while (I(node->id, S[i]) == 0 && node != root) {
      if (!set_I(node->id, S[i]))
        return 0;

      node = stree_get_parent(tree, node);

#ifdef STATS
      tree->num_compares++;
      tree->edges_traversed++;
#endif
    }
#ifdef STATS
    tree->num_compares++;
#endif

    /*
     * Step 2.
     */
    if (node == root && I(node->id, S[i]) == 0) {
      headlen = 0;

      if (!set_I(node->id, S[i]))
        return 0;
    }
    else {
      /*
       * Step 3.
       */
      v = node;

      t = 0;
      previous = NULL;
      while (L(node->id, S[i]) == NULL && node != root) {
        t += stree_get_edgelen(tree, node);
        previous = node;
        node = stree_get_parent(tree, node);

#ifdef STATS
        tree->num_compares++;
        tree->edges_traversed++;
#endif
      }
#ifdef STATS
      tree->num_compares++;
#endif

      if (node == root && L(node->id, S[i]) == NULL) {
        /*
         * Case 3a.
         */
        node = stree_find_child(tree, node, S[i]);
        headlen = t + 1;

#ifdef STATS
        tree->edges_traversed++;
#endif
      }
      else {
        /*
         * Case 3b.
         */
        vprime = node;
        vdblprime = L(node->id, S[i]);
        l = t;

#ifdef STATS
        tree->links_traversed++;
#endif

        if (l == 0) {
          node = vdblprime;
          headlen = stree_get_edgelen(tree, vdblprime);
        }
        else {
          c = stree_getch(tree, previous);
          node = stree_find_child(tree, vdblprime, c);
          headlen = l;

#ifdef STATS
          tree->edges_traversed++;
#endif
        }
      }
    }
          
    /*
     * Step 4.
     */
    if (headlen == stree_get_edgelen(tree, node))
      w = node;
    else {
      if ((w = int_stree_edge_split(tree, node, headlen)) == NULL)
        return 0;
      w->id = tree->num_nodes - 1;
      if (!copy_links(w->id, node->id))
        return 0;
    }

    /*
     * NOTE:  Since Head(i) is the label marking the path from the root to w,
     *        its length is the length of that label.  And, so the first
     *        character of the leaf edge is i + |Head(i)|, and the length
     *        of the edge to the leaf is M minus the edge's first character
     *        position.
     */
    headlen = stree_get_labellen(tree, w);
    edgepos = i + headlen;
    edgelen = M - edgepos;

    if (edgelen == 0) {
      if (int_stree_isaleaf(tree, w) &&
          (w = int_stree_convert_leafnode(tree, w)) == NULL)
        return 0;

      if (!int_stree_add_intleaf(tree, w, id, i))
        return 0;

      node = w;
    }
    else {
      if ((leaf = int_stree_new_leaf(tree, id, edgepos, i)) == NULL ||
          (w = int_stree_connect(tree, w, (STREE_NODE) leaf)) == NULL) {
        if (leaf != NULL)
          int_stree_free_leaf(tree, leaf);
        return 0;
      }
      leaf->id = tree->num_nodes;
      tree->num_nodes++;

      node = (STREE_NODE) leaf;
    }

    if (v != NULL) {
      if (!set_L(v->id, S[i], w))
        return 0;
    
      if (w->suffix_link == NULL)
        w->suffix_link = v;
      assert(w->suffix_link == v);
    }
  }

  return 1;
}


/*
 * The links data structure and get/set procedures.
 *
 * Note that this implementaion is not the naive
 * bit-vector/array-of-pointers-vector implementation of the L and I
 * vectors.  An earlier version of this algorithm used that
 * implementation, and it was using way too much memory (i.e., 50MB for
 * a string of 200,000 characters).

 * This implemenation uses a linked list of the set L and I values
 * for each node in the tree.  The linked list nodes give the list of
 * characters for which the I value is set, and subset of nodes whose
 * `link' field is non-NULL give the characters for which the L value
 * is set.  (So, when just the I value is set for a character x, there
 * is a node in the linked list for x, but its link field is NULL.)
 */

typedef struct listnode {
  char ch;
  STREE_NODE link;
  struct listnode *next;
} LINKNODE, *LINKS;

static LINKS *links = NULL;
static int links_size = 0;


/*
 * reset_links
 *
 * Initialize the links data structure for the next suffix tree build.
 *
 * Parameters:  none
 *
 * Returns:  non-zero on success, zero on error.
 */
static int reset_links()
{
  int i;

  links_size = 1024;
  if ((links = malloc(links_size * sizeof(LINKS))) == NULL)
    return 0;

  for (i=0; i < 1024; i++)
    links[i] = NULL;

  return 1;
}


/*
 * grow_links
 *
 * The links data structure is an array of linked lists (using the
 * identifiers associated with each suffix tree node).  This function
 * grows the links data structure when new suffix tree nodes are added
 * during the suffix tree construction.
 *
 * Parameter:  newid  -  a (possibly) new suffix tree node identifier.
 *
 * Returns:  non-zero on success in growing the structure, zero on error.
 */
static int grow_links(int newid)
{
  int i, newsize;

  if (newid < links_size)
    return 1;

  newsize = links_size + links_size;
  if (newid >= newsize)
    newsize = newid + 1;

  if ((links = realloc(links, newsize * sizeof(LINKS))) == NULL)
    return 0;

  for (i=links_size; i < newsize; i++)
    links[i] = NULL;
  links_size = newsize;

  return 1;
}


/*
 * I
 *
 * This function returns the I value for a suffix tree node and a character.
 *
 * Parameters:  id  -  a suffix tree node identifier.
 *              ch  -  a character.
 *
 * Returns:  The appropriate I value.
 */
static int I(int id, char ch)
{
  LINKS node;

  if (id >= links_size)
    return 0;

  for (node=links[id]; node != NULL; node=node->next)
    if (node->ch == ch)
      return 1;

  return 0;
}


/*
 * L
 *
 * This function returns the L value for a suffix tree node and a character.
 *
 * Parameters:  id  -  a suffix tree node identifier.
 *              ch  -  a character.
 *
 * Returns:  The appropriate L value.
 */
static STREE_NODE L(int id, char ch)
{
  LINKS node;

  if (id >= links_size)
    return NULL;

  for (node=links[id]; node != NULL; node=node->next)
    if (node->ch == ch)
      return node->link;

  return NULL;
}  


/*
 * set_I
 *
 * This function sets the I value for a suffix tree node and a character
 * to 1.
 *
 * Parameters:  id  -  a suffix tree node identifier.
 *              ch  -  a character.
 *
 * Returns:  non-zero on success, zero on error.
 */
static int set_I(int id, char ch)
{
  LINKS node, newlink;

  if (id >= links_size && !grow_links(id))
    return 0;

  for (node=links[id]; node != NULL; node=node->next)
    if (node->ch == ch)
      return 0;

  if ((newlink = malloc(sizeof(LINKNODE))) == NULL)
    return 0;

  newlink->ch = ch;
  newlink->link = NULL;
  newlink->next = links[id];
  links[id] = newlink;

  return 1;
}


/*
 * set_L
 *
 * This function sets the L value for a suffix tree node and a character
 * to a suffix tree node.
 *
 * Parameters:  id      -  a suffix tree node identifier.
 *              ch      -  a character.
 *              linkptr -  node to set the link pointer to
 *
 * Returns:  non-zero on success, zero on error.
 */
static int set_L(int id, char ch, STREE_NODE linkptr)
{
  LINKS node, newlink;

  if (id >= links_size && !grow_links(id))
    return 0;

  for (node=links[id]; node != NULL; node=node->next) {
    if (node->ch == ch) {
      if (node->link != NULL && node->link != linkptr)
        return 0;
      else {
        node->link = linkptr;
        return 1;
      }
    }
  }

  if ((newlink = malloc(sizeof(LINKNODE))) == NULL)
    return 0;

  newlink->ch = ch;
  newlink->link = linkptr;
  newlink->next = links[id];
  links[id] = newlink;

  return 1;
}


/*
 * copy_links
 *
 * This function copies the currently set I and L value of the `src'
 * suffix tree node to the `dest' suffix tree node.
 *
 * Parameters:  dest  -  the suffix tree identifier for the destination
 *              src   -  the suffix tree identifier for the source
 *
 * Returns:  non-zero on success, zero on error.
 */
static int copy_links(int dest, int src)
{
  LINKS node, newlink;

  assert(dest >= 0 && src >= 0 && src < links_size);

  if (dest >= links_size && !grow_links(dest))
    return 0;

  for (node=links[src]; node != NULL; node=node->next) {
    if ((newlink = malloc(sizeof(LINKNODE))) == NULL)
      return 0;

    newlink->ch = node->ch;
    newlink->link = NULL;
    newlink->next = links[dest];
    links[dest] = newlink;
  }

  return 1;
}


/*
 * free_links
 *
 * This function frees the links data structure.
 *
 * Parameters:  none
 * 
 * Returns: none
 */
static void free_links(void)
{
  int i;
  LINKS node, next;

  for (i=0; i < links_size; i++) {
    if (links[i] != NULL) {
      for (node=links[i]; node != NULL; node=next) {
        next = node->next;
        free(node);
      }
    }
  }

  free(links);
  links = NULL;
  links_size = 0;
}


/*
 * 
 * Construction shell functions for use in strmat.
 *
 */

/*
 * stree_weiner_build
 *
 * Build a suffix tree for a single string.
 *
 * Parameters:  string           -  the string
 *              build_policy     -  what type of build policy to use
 *              build_threshold  -  for LIST_THEN_ARRAY, when to move from
 *                                  linked list to complete array
 *
 * Returns:  the suffix tree, or NULL on an error.
 */
SUFFIX_TREE stree_weiner_build(STRING *string, int build_policy,
                                int build_threshold)
{
  SUFFIX_TREE tree;

  if (string == NULL || string->sequence == NULL || string->length == 0)
    return NULL;
  
  tree = stree_new_tree(string->alpha_size, 1, build_policy, build_threshold);
  if (tree == NULL)
    return NULL;

  reset_links();

  if (stree_weiner_add_string(tree, string->sequence,
                              string->raw_seq, string->length, 1) < 1) {
    stree_delete_tree(tree);
    free_links();
    return NULL;
  }

  free_links();

  return tree;
}


/*
 * stree_gen_weiner_build
 *
 * Build a generalized suffix tree for multiple strings.
 *
 * Parameters:  strings          -  the strings
 *              num_strings      -  the number of strings
 *              build_policy     -  what type of build policy to use
 *              build_threshold  -  for LIST_THEN_ARRAY, when to move from
 *                                  linked list to complete array
 *
 * Returns:  the suffix tree, or NULL on an error.
 */
SUFFIX_TREE stree_gen_weiner_build(STRING **strings, int num_strings,
                                    int build_policy, int build_threshold)
{
  int i;
  SUFFIX_TREE tree;

  if (strings == NULL || num_strings == 0)
    return NULL;

  tree = stree_new_tree(strings[0]->alpha_size, 0,
                        build_policy, build_threshold);
  if (tree == NULL)
    return NULL;

  reset_links();

  for (i=0; i < num_strings; i++) {
    if (stree_weiner_add_string(tree, strings[i]->sequence,
                                strings[i]->raw_seq,
                                strings[i]->length, i+1) < 1) {
      stree_delete_tree(tree);
      free_links();
      return NULL;
    }
  }

  free_links();

  return tree;
}
