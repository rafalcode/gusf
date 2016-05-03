/*
 * strmat.c
 *
 *  GENERAL DESCRIPTION
 *  This file contains the main body of the user interface.
 *  All other functions are accessed through the functions
 *  (directly or indirectly) here.    
 *
 *  ADDING ALGORITHMS
 *  All patterns and texts for particular algorithms are
 *  returned by the get_string function.  The get_string
 *  function returns a *STRING after querying the user.
 *  The get_string function returns NULL if the
 *  user chooses to return the the previous menu, so
 *  calls made to get_string should check for NULL, and
 *  cease execution of the algorithm.
 *
 * NOTES:
 *    8/94  -  Original Implementation (as menu.c) (Gene Fodor)
 *    3/96  -  Combined main.c and menu.c, changed its name to strmat.c
 *             and altered the code as the algorithm code was modularized
 *             (James Knight)
 *    7/96  -  Finished the modularization (James Knight)
 *    2/00  -  Changed type of function main() to int (Jens Stoye)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "strmat.h"
#include "strmat_alpha.h"
#include "strmat_print.h"
#include "strmat_seqary.h"
#include "strmat_stubs.h"
#include "strmat_stubs2.h"
#include "strmat_stubs3.h"
#include "strmat_stubs4.h"
#include "strmat_util.h"

#include "sary_match.h"
#include "stree_ukkonen.h"


#define ON 1
#define OFF 0

static char *choice;
static int ch_len;
static int stree_build_policy = SORTED_LIST;
static int stree_build_threshold = 10;
static int stree_print_flag = ON;
static int stats_flag = ON;

// FILE *fpout = stdout;
#define fpout stdout


/*
 * Forward prototypes.
 */

void util_menu(void);
void basic_alg_menu(void);
void z_alg_menu(void);
void suf_tree_menu(void);
void suf_ary_menu(void);
void repeats_menu(void);
void set_display_options(void);



/**********************************************************************
 *  Function main()                                              
 *                                                                   
 *  Parameter:                                                       
 * 
 *  This is the main function.  It allocates space for the initial data
 *  structures and displays the main menu. You can choose an algorithm
 *  sub-menu, or the String Utilities Menu.
 *                                                                   
 **********************************************************************/
int main(void)
{
  /*
   * Initialize the array of sequences.
   */
  create_seq_array();

  /*
   * The main loop.
   */
  while (1) {
    printf("\n**          Main Menu        **\n");
    printf("\n");
    printf("1)  Basic Search Algorithms\n");
    printf("2)  Z-value Algorithms\n");
    printf("3)  Suffix Tree Algorithms\n");
    printf("4)  Suffix Array Algorithms\n");
    printf("5)  Repeat Algorithms\n");
    printf("*)  String Utilities\n");
    printf("0)  Exit Program\n");
    printf("\nEnter Selection: ");

    while ((choice = my_getline(stdin, &ch_len)) == NULL) ;
    switch (choice[0]) {
    case '0':
      return 0;

    case '1':
      basic_alg_menu(); 
      break;

    case '2':
      z_alg_menu();
      break;

    case '3':
      suf_tree_menu();
      break;
      
    case '4':
      suf_ary_menu();
      break;

    case '5':
      repeats_menu();
      break;

    case '*':
      util_menu();
      break;

    default:
      printf("\nThat is not a choice.\n");
      break;
    }
  }
} /* main() */


/**********************************************************************
 *  Function  util_menu()
 *                                                                    
 *  Parameter:                                                       
 *   
 *                                                              
 *  This function is the backbone of the interface.  All file input/ouput
 *  happens here.  This is also where the user can key in their own string
 *  if desired.  Here the user can view, delete and list the available 
 *  sequences.
 *  
 *                                                                   
 **********************************************************************/
void util_menu()  
{
  int num_seqs, num_lines;

  while (1) {
    num_seqs = get_num_sequences();

    num_lines = 14;
    printf("\n**   String Utilites Menu    **\n\n");
    printf("1)  Read formatted file\n");
    printf("2)  Read unformatted file\n");
    printf("3)  Create new sequence\n");
    if (num_seqs == 0)
      printf("4)  List sequences (currently available: None)\n");
    else
      printf("4)  List sequences (currently available: 1 - %d)\n", num_seqs);
    printf("5)  Print sequence\n");
    printf("6)  Save sequences\n");
    printf("7)  Delete sequences\n");
    printf("8)  Set output options\n");
    printf("0)  Exit\n");
    printf("\nEnter Selection: ");

    while ((choice = my_getline(stdin, &ch_len)) == NULL) ;
    switch (choice[0]) {
    case '0':
      return;

    case '1':
      fread_formatted();
      break;

    case '2':
      fread_unformatted();
      break;

    case '3':
      type_in_seq();
      break;

    case '4':
      list_sequences(num_lines);
      break;

    case '5':
      print_seq(num_lines);
      break;    

    case '6':
      fwrite_formatted();      
      break;

    case '7':
      delete_seq();
      break;

    case '8':
      set_display_options();
      break;

    default:
      printf("\nThat is not a choice.\n");
    }
  }
}


/**********************************************************************
 *  Function  basic_alg_menu();
 *                                                                    
 *  Parameter:                                                       
 *   
 *                                                              
 *  This function prompts user to choose a basic algorithm to execute
 *  such as Boyer-Moore.  The utilities menu is also available from
 *  this menu.
 *                                                                   
 **********************************************************************/
void basic_alg_menu()
{
  int i, status, num_lines, alpha_size, num_patterns;
  char ch;
  STRING *text, *pattern, **patterns;

  alpha_size = 0;
  while (1)  {
    num_lines = 22;
    printf("\n**   Basic Search Algorithm Menu    **\n\n");
    printf("1)  Naive Algorithm\n");
    printf("2)  Boyer-Moore Variations\n");
    printf("     a) Bad character rule\n");
    printf("     b) Extended bad character rule\n");
    printf("     c) Good suffix & bad character rules\n");
    printf("     d) Good suffix & extended bad character rules\n");
    printf("3)  Knuth-Morris-Pratt (original preprocessing)\n");
    printf("     a) using sp values\n");
    printf("     b) using sp' values\n");
    printf("4)  Aho-Corasick Set Matching\n");
    printf("5)  Boyer-Moore Set Matching\n");
    printf("     a) Bad character rule only\n");
    printf("     b) Good suffix rule using keyword and suffix trees\n");
    printf("     c) Good suffix rule using suffix tree only\n");
    printf("     d) using original Boyer-Moore (1c) on each pattern\n");
    printf("*)  String Utilites\n");
    printf("0)  Exit\n");
    printf("\nEnter Selection: ");
  
    while ((choice = my_getline(stdin, &ch_len)) == NULL) ;
    switch (choice[0]) {
    case '0':
      return;

    case '1':
      if (!(pattern = get_string("pattern")) || !(text = get_string("text")))
        continue;

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mprintf("\nThe pattern:\n");
      terse_print_string(pattern);
      mprintf("\nThe text:\n");
      terse_print_string(text);
      mputc('\n');
        
      status = map_sequences(text, pattern, NULL, 0);
      if (status != -1) {
        mprintf ("Executing naive search algorithm...\n\n");
        strmat_naive_match(pattern, text, stats_flag);
        unmap_sequences(text, pattern, NULL, 0);
      }
      mend(num_lines);
      putchar('\n');
      break;

    case '2':
      ch = choice[1];
      if (toupper(ch) != 'A' && toupper(ch) != 'B' &&
          toupper(ch) != 'C' && toupper(ch) != 'D') {
        printf("\nYou must specify the Boyer-Moore variation"
               " (as in '2a' or '2c').\n");
        continue;
      }

      if (!(pattern = get_string("pattern")) || !(text = get_string("text")))
        continue;

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mprintf("\nThe pattern:\n");
      terse_print_string(pattern);
      mprintf("\nThe text:\n");
      terse_print_string(text);
      mputc('\n');

      status = map_sequences(text, pattern, NULL, 0);
      if (status != -1) {
        mprintf("Executing Boyer-Moore algorithm...\n\n");
        switch (toupper(ch)) {
        case 'A':  strmat_bmbad_match(pattern, text, stats_flag);  break;
        case 'B':  strmat_bmext_match(pattern, text, stats_flag);  break;
        case 'C':  strmat_bmgood_match(pattern, text, stats_flag);  break;
        case 'D':  strmat_bmextgood_match(pattern, text, stats_flag);  break;
        }
        unmap_sequences(text, pattern, NULL, 0);
      }
      mend(num_lines);
      putchar('\n');
      break;

    case '3':
      ch = choice[1];
      if (toupper(ch) != 'A' && toupper(ch) != 'B') {
        printf("\nYou must specify the KMP variation (as in '3a' or '3b').\n");
        continue;
      }

      if (!(pattern = get_string("pattern")) || !(text = get_string("text")))
        continue;

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mprintf("\nThe pattern:\n");
      terse_print_string(pattern);
      mprintf("\nThe text:\n");
      terse_print_string(text);
      mputc('\n');
        
      status = map_sequences(text, pattern, NULL, 0);
      if (status != -1) {
        if (toupper(ch) == 'A') {
          mprintf ("Executing KMP with sp values...\n\n");
          strmat_kmp_sp_orig_match(pattern, text, stats_flag);
        }
        else {
          mprintf ("Executing KMP with sp' values...\n\n");
          strmat_kmp_spprime_orig_match(pattern, text, stats_flag);
        }
        unmap_sequences(text, pattern, NULL, 0);
      }
      mend(num_lines);
      putchar('\n');
      break;

    case '4':
      if (!(patterns = get_string_ary("list of patterns", &num_patterns)))
        continue;
      if (!(text = get_string("text"))) {
        free(patterns);
        continue;
      }

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mprintf("\nThe patterns:\n");
      for (i=0; i < num_patterns; i++) {
        mprintf("%2d)", i + 1);
        terse_print_string(patterns[i]);
      }
      mprintf("\nThe text:\n");
      terse_print_string(text);
      mputc('\n');

      status = map_sequences(text, NULL, patterns, num_patterns);
      if (status != -1) {
        mprintf("Executing Aho-Corasick algorithm...\n\n");
        strmat_ac_match(patterns, num_patterns, text, stats_flag);
        unmap_sequences(text, NULL, patterns, num_patterns);
      }
      mend(num_lines);
      putchar('\n');

      free(patterns);
      break;

    case '5':
      ch = toupper(choice[1]);
      if (ch != 'A' && ch != 'B' && ch != 'C' && ch != 'D') {
        printf("\nYou must specify the set Boyer-Moore variation"
               " (as in '5a' or '5c').\n");
        continue;
      }

      if (!(patterns = get_string_ary("list of patterns", &num_patterns)))
        continue;
      if (!(text = get_string("text"))) {
        free(patterns);
        continue;
      }

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mprintf("\nThe patterns:\n");
      for (i=0; i < num_patterns; i++) {
        mprintf("%2d)", i + 1);
        terse_print_string(patterns[i]);
      }
      mprintf("\nThe text:\n");
      terse_print_string(text);
      mputc('\n');

      status = map_sequences(text, NULL, patterns, num_patterns);
      if (status != -1) {
        mprintf("Executing Boyer-Moore set matching algorithm...\n\n");
        switch (toupper(ch)) {
        case 'A':  strmat_bmset_badonly_match(patterns, num_patterns,
                                              text, stats_flag);       break;
        case 'B':  strmat_bmset_2trees_match(patterns, num_patterns,
                                             text, stats_flag);        break;
        case 'C':  strmat_bmset_1tree_match(patterns, num_patterns,
                                            text, stats_flag);         break;
        case 'D':  strmat_bmset_naive_match(patterns, num_patterns,
                                            text, stats_flag);         break;
        }
        unmap_sequences(text, NULL, patterns, num_patterns);
      }
      mend(num_lines);
      putchar('\n');

      free(patterns);
      break;

    case '*':
      util_menu();
      break;

    default:
      printf("\nThat is not a choice.\n");
    }
  }
}





/**********************************************************************
 *  Function  z_alg_menu();
 *                                                                    
 *  Parameter:                                                       
 *   
 *                                                              
 *  This function prompts user to create and utilize the Z values 
 *  for a string.  The utilities menu is also available from this menu.
 *                                                                   
 **********************************************************************/
void z_alg_menu()
{
  int status, num_lines;
  char ch;
  STRING *text, *pattern;

  while (1) {
    num_lines = 12;
    printf("\n**   Z-value Algorithm Menu    **\n\n");
    printf("1)  Build Z values for a sequence\n");
    printf("2)  Exact matching using Z values\n");
    printf("3)  Knuth-Morris-Pratt  (Z-values preprocessing)\n");
    printf("     a) using sp values\n");
    printf("     b) using sp' values\n");
    printf("*)  String Utilites\n");
    printf("0)  Exit\n");
    printf("\nEnter Selection: ");
  
    while ((choice = my_getline(stdin, &ch_len)) == NULL) ;
    switch (choice[0]) {
    case '0':
      return;

    case '1':
      if (!(text = get_string("text")))
        continue;

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mputs("The string:\n");
      terse_print_string(text);

      status = map_sequences(text, NULL, NULL, 0);
      if (status != -1) {
        mputs("Building Z values...\n\n");
        strmat_z_build(text, stats_flag);
        unmap_sequences(text, NULL, NULL, 0);
      }
      mend(num_lines);
      putchar('\n');
      break;

    case '2':
      if (!(pattern = get_string("pattern")) || !(text = get_string("text")))
        continue;

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mputs("\nThe pattern:\n");
      terse_print_string(pattern);
      mputs("\nThe text:\n");
      terse_print_string(text);
      mputc('\n');

      status = map_sequences(text, pattern, NULL, 0);
      if (status != -1) {
        mputs("Executing exact matching with Z values algorithm...\n\n");
        strmat_z_match(pattern, text, stats_flag);
        unmap_sequences(text, pattern, NULL, 0);
      }
      mend(num_lines);
      putchar('\n');
      break;

    case '3':
      ch = choice[1];
      if (toupper(ch) != 'A' && toupper(ch) != 'B') {
        printf("\nYou must specify the KMP variation (as in '3a' or '3b').\n");
        continue;
      }

      if (!(pattern = get_string("pattern")) || !(text = get_string("text")))
        continue;

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mprintf("\nThe pattern:\n");
      terse_print_string(pattern);
      mprintf("\nThe text:\n");
      terse_print_string(text);
      mputc('\n');
        
      status = map_sequences(text, pattern, NULL, 0);
      if (status != -1) {
        if (toupper(ch) == 'A') {
          mprintf ("Executing KMP with sp values...\n\n");
          strmat_kmp_sp_z_match(pattern, text, stats_flag);
        }
        else {
          mprintf ("Executing KMP with sp' values...\n\n");
          strmat_kmp_spprime_z_match(pattern, text, stats_flag);
        }
        unmap_sequences(text, pattern, NULL, 0);
      }
      mend(num_lines);
      putchar('\n');
      break;

    case '*':
      util_menu();
      break;

    default:
      printf("\nThat is not a choice.\n");
    }
  }
}

/**********************************************************************
 *  Function  suf_tree_menu()
 *                                                                    
 *  Parameter:                                                       
 *   
 *                                                              
 *  This function prompts user to choose an algorithm to create 
 *  a suffix tree and use it.  
 *  The utilities menu is also available from this menu.
 *                                                                   
 **********************************************************************/
void suf_tree_menu()
{
  int i, status, num_lines, num_strings;
  char ch;
  STRING *pattern, **strings, *text;

  while (1) {
    num_lines = 18;

    printf("\n**   Suffix Tree Menu    **\n\n");
    printf("1)  Build a suffix tree using Ukkonen's algorithm\n");
    printf("2)  Build a suffix tree using Weiner's algorithm\n");
    printf("3)  Exact matching using a suffix tree for the text\n");
    printf("4)  Walk around a suffix tree\n");
    printf("5)  Compute the LCA values for a suffix tree\n");
    printf("     a) using the naive LCA algorithm\n");
    printf("     b) using the constant time LCA algorithm\n");
    printf("6)  Compute Lempel-Ziv decomposition\n");
    printf("     a) original version (f-factorization)\n");
    printf("     b) nonoverlapping blocks (as in the book)\n");
    printf("8)  Set suffix tree build policy (current: ");
    switch(stree_build_policy)  {
    case LINKED_LIST:      printf("linked list)\n");  break;
    case SORTED_LIST:      printf("sorted list)\n");  break;
    case LIST_THEN_ARRAY:  printf("list then array, threshold %d)\n",
                                  stree_build_threshold);  break;
    case COMPLETE_ARRAY:   printf("complete array)\n");  break;
    }
    printf("9)  Suffix tree print toggle (current: %s)\n",
           (stree_print_flag == ON ? "on" : "off"));
    printf("*)  String Utilites\n");
    printf("0)  Exit\n");
    printf("\nEnter Selection: ");
 
    while ((choice = my_getline(stdin, &ch_len)) == NULL) ;

    switch (choice[0]) {
    case '0':
      return;

    case '1':
      strings = get_string_ary("list of sequences", &num_strings);
      if (strings == NULL)
        continue;

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mprintf("\nThe sequences:\n");
      for (i=0; i < num_strings; i++) {
        mprintf("%2d)", i + 1);
        terse_print_string(strings[i]);
      }
      mputc('\n');

      status = map_sequences(NULL, NULL, strings, num_strings);
      if (status != -1) {
        mprintf("Executing Ukkonen's Algorithm...\n\n");
        strmat_ukkonen_build(strings, num_strings, stree_build_policy,
                             stree_build_threshold, stats_flag,
                             stree_print_flag);
        unmap_sequences(NULL, NULL, strings, num_strings);
      }
      mend(num_lines);
      putchar('\n');

      free(strings);
      break;

    case '2':
      strings = get_string_ary("list of sequences", &num_strings);
      if (strings == NULL)
        continue;

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mprintf("\nThe sequences:\n");
      for (i=0; i < num_strings; i++) {
        mprintf("%2d)", i + 1);
        terse_print_string(strings[i]);
      }
      mputc('\n');

      status = map_sequences(NULL, NULL, strings, num_strings);
      if (status != -1) {
        mprintf("Executing Weiner's Algorithm...\n\n");
        strmat_weiner_build(strings, num_strings, stree_build_policy,
                            stree_build_threshold, stats_flag,
                            stree_print_flag);
        unmap_sequences(NULL, NULL, strings, num_strings);
      }
      mend(num_lines);
      putchar('\n');

      free(strings);
      break;

    case '3':
      if (!(pattern = get_string("pattern")) ||
          !(strings = get_string_ary("list of sequences", &num_strings)))
        continue;

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mprintf("\nThe pattern:\n");
      terse_print_string(pattern);
      mprintf("\nThe texts:\n");
      for (i=0; i < num_strings; i++) {
        mprintf("%2d)", i + 1);
        terse_print_string(strings[i]);
      }
      mputc('\n');

      status = map_sequences(NULL, pattern, strings, num_strings);
      if (status != -1) {
        mprintf("Executing exact matching with a suffix tree...\n\n");
        strmat_stree_match(pattern, strings, num_strings, stree_build_policy,
                           stree_build_threshold, stats_flag);
        unmap_sequences(NULL, pattern, strings, num_strings);
      }
      mend(num_lines);
      putchar('\n');

      free(strings);
      break;
        
    case '4':
      strings = get_string_ary("list of sequences", &num_strings);
      if (strings == NULL)
        continue;

      status = map_sequences(NULL, NULL, strings, num_strings);
      if (status != -1) {
        strmat_stree_walkaround(strings, num_strings, stree_build_policy,
                                stree_build_threshold);
        unmap_sequences(NULL, NULL, strings, num_strings);
      }
      putchar('\n');

      free(strings);
      break;

    case '5':
      ch = toupper(choice[1]);
      if (ch != 'A' && ch != 'B') {
        printf("\nYou must specify which type of LCA algorithm to use "
               "(as in '3a' or '3b').\n");
        continue;
      }

      strings = get_string_ary("list of sequences", &num_strings);
      if (strings == NULL)
        continue;

      status = map_sequences(NULL, NULL, strings, num_strings);
      if (status != -1) {
        if (ch == 'A')
          strmat_stree_naive_lca(strings, num_strings, stree_build_policy,
                                 stree_build_threshold, stats_flag);
        else
          strmat_stree_lca(strings, num_strings, stree_build_policy,
                           stree_build_threshold, stats_flag);
        unmap_sequences(NULL, NULL, strings, num_strings);
      }
      putchar('\n');

      free(strings);
      break;

    case '6':
      ch = toupper(choice[1]);
      if (ch!='A' && ch!='B') {
        printf("\nYou must specify which type of decomposition to compute "
               "(as in '6a' or '6b').\n");
        continue;
      }

      if (!(text = get_string("string")))
        continue;

      mstart(stdin, fpout, OK, OK, 0, NULL);
      mprintf("\nThe string:\n");
      terse_print_string(text);
      mputc('\n');

      status = map_sequences(text, NULL, NULL, 0);
      if(status != -1) {
        strmat_stree_lempel_ziv(text, stree_build_policy,
                                stree_build_threshold, stats_flag,ch);
        unmap_sequences(text, NULL, NULL, 0);
      }
      mend(num_lines);
      putchar('\n');
      break;

    case '8':
      choice = "0";
      while (choice[0] != '1' && choice[0] != '2' &&
             choice[0] != '3' && choice[0] != '4') {
        printf("\n**  Suffix Tree Build Policies **\n");
        printf("\n(1 - linked list, 2 - sorted list, 3 - linked list/array,"
               " 4 - complete array)\n");
        printf("Enter Build Policy [%d]: ",
               (stree_build_policy == LINKED_LIST ? 1
                  : (stree_build_policy == SORTED_LIST ? 2
                       : (stree_build_policy == LIST_THEN_ARRAY ? 3 : 4))));

        if ((choice = my_getline(stdin, &ch_len)) == NULL || choice[0] == '\0')
          break;
      
        switch (choice[0]) {
        case '1':
          stree_build_policy = LINKED_LIST;
          break;

        case '2':
          stree_build_policy = SORTED_LIST;
          break;
        case '3':
          stree_build_policy = LIST_THEN_ARRAY;
          break;

        case '4':
          stree_build_policy = COMPLETE_ARRAY;
          break;

        default:
          printf("\nThat is not a choice.\n");
        }
      }
      if (stree_build_policy == LIST_THEN_ARRAY) {
        printf("\nEnter Build Threshold [%d]: ", stree_build_threshold);
        if ((choice = my_getline(stdin, &ch_len)) != NULL)
          sscanf(choice, "%d", &stree_build_threshold);
      }
      putchar('\n');
      break;

    case '9':
      if (stree_print_flag == ON)
        stree_print_flag = OFF;
      else
        stree_print_flag = ON;
      break;

    case '*':
      util_menu();
      break;
   
    default:
      printf("\nThat is not a choice.\n");
    }
  }
}

/**********************************************************************
 *  Function  suf_ary_menu()
 *                                                                    
 *  Parameter:                                                       
 *   
 *                                                              
 *  This function prompts user to choose an algorithm to create 
 *  a suffix array and use it.  
 *  The utilities menu is also available from this menu.
 *                                                                   
 **********************************************************************/
void suf_ary_menu()
{
  int status, num_lines;
  STRING *spt, *pattern, *text;

  while (1)  {
    num_lines = 13;
    printf("\n**   Suffix Array Menu    **\n\n");
    printf("1)  Build suffix array using quick sort\n");
    printf("2)  Build suffix array (Zerkle's version)\n");
    printf("3)  Build suffix array from a suffix tree\n");
    printf("4)  Exact matching using suffix array and naive algorithm\n");
    printf("5)  Exact matching using suffix array and mlr accelerant\n");
    printf("6)  Exact matching using suffix array and lcp super-accelerant\n");
    printf("*)  String Utilites\n");
    printf("0)  Exit\n");
    printf("\nEnter Selection: ");
  
    while ((choice = my_getline(stdin, &ch_len)) == NULL) ;
    switch (choice[0]) {
    case '0':
      return;

    case '1':
      if (!(spt = get_string("sequence")))
        continue;

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mprintf("\nThe sequence:\n");
      terse_print_string(spt);
      mputc('\n');

      status = map_sequences(spt, NULL, NULL, 0);
      if (status != -1) {
        mprintf("Building suffix array using qsort...\n\n");
        strmat_sary_qsort(spt, stats_flag);
        unmap_sequences(spt, NULL, NULL, 0);
      }
      mend(num_lines);
      putchar('\n');
      break;
        
    case '2':
      if (!(spt = get_string("sequence")))
        continue;

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mprintf("\nThe sequence:\n");
      terse_print_string(spt);
      mputc('\n');

      status = map_sequences(spt, NULL, NULL, 0);
      if (status != -1) {
        mprintf("Executing Zerkle's algorithm...\n\n");
        strmat_sary_zerkle(spt, stats_flag);
        unmap_sequences(spt, NULL, NULL, 0);
      }
      mend(num_lines);
      putchar('\n');
      break;
      
    case '3':
      if (!(spt = get_string("sequence")))
        continue;

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mprintf("\nThe sequence:\n");
      terse_print_string(spt);
      mputc('\n');

      status = map_sequences(spt, NULL, NULL, 0);
      if (status != -1) {
        mprintf("Building suffix array from suffix tree...\n\n");
        strmat_sary_stree(spt, stats_flag);
        unmap_sequences(spt, NULL, NULL, 0);
      }
      mend(num_lines);
      putchar('\n');
      break;

    case '4':
      if (!(pattern = get_string("pattern")) || !(text = get_string("text")))
        continue;

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mprintf("\nThe pattern:\n");
      terse_print_string(pattern);
      mprintf("\nThe text:\n");
      terse_print_string(text);    
      mputc('\n');

      status = map_sequences(text, pattern, NULL, 0);
      if (status != -1) {
        mprintf("Executing exact matching using suffix array...\n\n");
        strmat_sary_match_naive(pattern, text, stats_flag);
        unmap_sequences(text, pattern, NULL, 0);
      }
      mend(num_lines);
      putchar('\n');
      break;  

    case '5':
      if (!(pattern = get_string("pattern")) || !(text = get_string("text")))
        continue;

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mprintf("\nThe pattern:\n");
      terse_print_string(pattern);
      mprintf("\nThe text:\n");
      terse_print_string(text);    
      mputc('\n');

      status = map_sequences(text, pattern, NULL, 0);
      if (status != -1) {
        mprintf("Executing exact matching using suffix array...\n\n");
        strmat_sary_match_mlr(pattern, text, stats_flag);
        unmap_sequences(text, pattern, NULL, 0);
      }
      mend(num_lines);
      putchar('\n');
      break;  

    case '6':
      if (!(pattern = get_string("pattern")) || !(text = get_string("text")))
        continue;

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mprintf("\nThe pattern:\n");
      terse_print_string(pattern);
      mprintf("\nThe text:\n");
      terse_print_string(text);    
      mputc('\n');

      status = map_sequences(text, pattern, NULL, 0);
      if (status != -1) {
        mprintf("Executing exact matching using suffix array...\n\n");
        strmat_sary_match_lcp(pattern, text, stats_flag);
        unmap_sequences(text, pattern, NULL, 0);
      }
      mend(num_lines);
      putchar('\n');
      break;  

    case '*':
      util_menu();
      break;

    default:
      printf("\nThat is not a choice.\n");
    }
  }
}


/**********************************************************************
 *  Function  repeats_menu()
 *                                                                    
 *  Parameter:                                                       
 *   
 *                                                              
 *  This function prompts user to choose an algorithm to compute
 *  repeats and similar things.
 *  The utilities menu is also available from this menu.
 *                                                                   
 **********************************************************************/
void repeats_menu()
{
  static int smax_percent = 0;
  static int smax_minlen = 0;
  int status, num_lines;
  char ch;
  STRING *text;

  while (1) {
    num_lines = 20;

    printf("\n**   Repeats Menu    **\n\n");
    printf("1)  Find primitive tandem repeats (Crochemore's algorithm)\n");
    printf("2)  Find supermaximals and near supermaximals of a string\n");
    printf("3)  Find nonoverlapping maximals of a string"
           " (Crochemore variant)\n");
    printf("4)  Find nonoverlapping maximals of a string"
           " (big path algorithm)\n");
    printf("5)  Find tandem repeats/tandem arrays using the suffix tree\n");
    printf("6)  Find vocabulary of tandem repeats (and more) using\n");
    printf("     a) Ziv-Lempel decomposition\n");
    printf("     b) nonoverlapping blocks decomposition (as in the book)\n");
    printf("7)  Find occurrences in linear time (without suffix tree) using\n");
    printf("     a) Ziv-Lempel decomposition\n");
    printf("     b) nonoverlapping blocks decomposition (as in the book)\n");
    printf("*)  String Utilites\n");
    printf("0)  Exit\n");
    printf("\nEnter Selection: ");
 
    while ((choice = my_getline(stdin, &ch_len)) == NULL) ;

    switch (choice[0]) {
    case '0':
      return;

    case '1':
      if (!(text = get_string("string")))
        continue;

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mprintf("\nThe string:\n");
      terse_print_string(text);
      mputc('\n');

      status = map_sequences(text, NULL, NULL, 0);
      if (status != -1) {
        strmat_repeats_primitives(text, stats_flag);
        unmap_sequences(text, NULL, NULL, 0);
      }
      mend(num_lines);
      putchar('\n');
      break;

    case '2':
      if (!(text = get_string("text")))
        continue;

      smax_percent = get_bounded("Percent Supermaximal", 0, 100, smax_percent);
      printf("\n");
      if (smax_percent == -1)
        continue;

      smax_minlen = get_bounded("Supermax. Minimum Length", 0, text->length,
                                smax_minlen);
      printf("\n");
      if (smax_minlen == -1)
        continue;
                                         
      mstart(stdin, fpout, OK, OK, 5, NULL);
      mprintf("\nThe text:\n");
      terse_print_string(text);
      mputc('\n');

      status = map_sequences(text, NULL, NULL, 0);
      if (status != -1) {
        mprintf("Finding the supermaximals...\n\n");
        strmat_repeats_supermax(text, smax_percent, smax_minlen);
        unmap_sequences(text, NULL, NULL, 0);
      }
      mend(num_lines);
      putchar('\n');
      break;

    case '3':
      if (!(text = get_string("string")))
        continue;

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mprintf("\nThe string:\n");
      terse_print_string(text);
      mputc('\n');

      status = map_sequences(text, NULL, NULL, 0);
      if (status != -1) {
        strmat_repeats_nonoverlapping(text, stats_flag);
        unmap_sequences(text, NULL, NULL, 0);
      }
      mend(num_lines);
      putchar('\n');
      break;

    case '4':
      if (!(text = get_string("string")))
        continue;

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mprintf("\nThe string:\n");
      terse_print_string(text);
      mputc('\n');

      status = map_sequences(text, NULL, NULL, 0);
      if (status != -1) {
        strmat_repeats_bigpath(text, stree_build_policy, stree_build_threshold,
                               stats_flag);
        unmap_sequences(text, NULL, NULL, 0);
      }
      mend(num_lines);
      putchar('\n');
      break;

    case '5':
      if (!(text = get_string("string")))
        continue;

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mprintf("\nThe string:\n");
      terse_print_string(text);
      mputc('\n');

      status = map_sequences(text, NULL, NULL, 0);
      if (status != -1) {
        strmat_repeats_tandem(text, stree_build_policy, stree_build_threshold,
                              stats_flag);
        unmap_sequences(text, NULL, NULL, 0);
      }
      mend(num_lines);
      putchar('\n');
      break;

    case '6':
      ch = toupper(choice[1]);
      if (ch!='A' && ch!='B') {
        printf("\nYou must specify which type of decomposition to use"
               "(as in '6a' or '6b').\n");
        continue;
      }

      if (!(text = get_string("string")))
        continue;

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mprintf("\nThe string:\n");
      terse_print_string(text);
      mputc('\n');

      status = map_sequences(text, NULL, NULL, 0);
      if (status != -1) {
        strmat_repeats_vocabulary(text, stree_build_policy,
                                  stree_build_threshold, stats_flag,ch);
        unmap_sequences(text, NULL, NULL, 0);
      }
      mend(num_lines);
      putchar('\n');
      break;

    case '7':
      ch = toupper(choice[1]);
      if (ch!='A' && ch!='B') {
        printf("\nYou must specify which type of decomposition to use"
               "(as in '7a' or '7b').\n");
        continue;
      }

      if (!(text = get_string("string")))
        continue;

      mstart(stdin, fpout, OK, OK, 5, NULL);
      mprintf("\nThe string:\n");
      terse_print_string(text);
      mputc('\n');

      status = map_sequences(text, NULL, NULL, 0);
      if (status != -1) {
        strmat_repeats_linear_occs(text, stree_build_policy,
                                   stree_build_threshold, stats_flag,ch);
        unmap_sequences(text, NULL, NULL, 0);
      }
      mend(num_lines);
      putchar('\n');
      break;

    case '*':
      util_menu();
      break;
   
    default:
      printf("\nThat is not a choice.\n");
    }
  }
}

/********************************************************************
*  Function set_display options
*
*  Parameters none.
*
*  Prints a strings information using the more interface
*  the output is sized to the screen and is displayed with
*  row and col labels.
*
*********************************************************************/
void set_display_options(void)
{
  FILE *newfpout;

  while (1) {
    printf("\nOptions (1 - redirect output to file, 2 - reset to screen,\n");
    printf("         3 - turn stats %s, 0 - Exit)\n",
           (stats_flag ? "off" : "on"));
    printf("Enter Selection: ");

    while ((choice = my_getline(stdin, &ch_len)) == NULL) ;
    switch (choice[0]) {
    case '0':
      return;

    case '1':
      choice[0] = '\0';
      while (choice != NULL && choice[0] == '\0') {
        printf("\nEnter file name (Ctl-D to cancel): ");
        choice = my_getline(stdin, &ch_len);
      }
      if (choice == NULL) {
        printf("\n\n");
        continue;
      }

      printf("\nRedirecting output to %s...", choice);
      if ((newfpout = fopen(choice, "w")) == NULL)
        fprintf(stderr,"\nError:  could not open file %s for output.\n\n",
                choice);
      else {
        if (fpout != stdout)
          fclose(fpout);
        fpout = newfpout;
        printf("done.\n\n");
      }
      break;

    case '2':
      if (fpout != stdout) {
        fclose(fpout);
        fpout = stdout;
      }
      putchar('\n');
      break;

    case '3':
      stats_flag = (stats_flag ? OFF : ON);
      putchar('\n');
      break;

    default:
      printf("\nThat is not a choice.\n");
    }
  }
}




int my_itoalen(int num)
{
  int i;

  for (i=1; num >= 10; i++)
    num /= 10;
  return i;
}


char *my_getline(FILE *fp, int *len_out)
{
  static int bufsize = 0;
  static char *buffer = NULL;
  int size, len;
  char *s;

  if (!buffer) {
    bufsize = 20;
    if ((buffer = malloc(bufsize)) == NULL) {
      fprintf(stderr, "Memory Error:  Ran out of memory.\n");
      return NULL;
    }
  }  

  s = buffer + bufsize - 2;
  *s = '\0';
  if (fgets(buffer, bufsize, fp) == NULL)
    return NULL;
  else if (!*s || *s == '\n') {
    len = strlen(buffer);
    if (buffer[len-1] == '\n')
      buffer[--len] = '\0';
    if (len_out) *len_out = len;
    return buffer;
  }

  while (1) {
    size = bufsize - 1;
    bufsize += bufsize;
    if ((buffer = realloc(buffer, bufsize)) == NULL) {
      fprintf(stderr, "Memory Error:  Ran out of memory.\n");
      return NULL;
    }

    s = buffer + bufsize - 2;
    *s = '\0';
    if (fgets(buffer + size, bufsize - size, fp) == NULL) {
      len = size;
      if (buffer[len-1] == '\n')
        buffer[--len] = '\0';
      if (len_out) *len_out = len;
      return buffer;
    }
    else if (!*s || *s == '\n') {
      len = size + strlen(buffer + size);
      if (buffer[len-1] == '\n')
        buffer[--len] = '\0';
      if (len_out) *len_out = len;
      return buffer;
    }
  }
}

