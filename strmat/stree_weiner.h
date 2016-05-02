
#ifndef _STREE_WEINER_H_
#define _STREE_WEINER_H_

int stree_weiner_add_string(SUFFIX_TREE tree, char *S, char *Sraw,
                            int M, int strid);

SUFFIX_TREE stree_weiner_build(STRING *string, int build_policy,
                               int build_threshold);
SUFFIX_TREE stree_gen_weiner_build(STRING **strings, int num_strings,
                                   int build_policy, int build_threshold);

#endif
