/**********************************************************************
 * SUF.C
 *
 * Revision 1.0, December 13 1992
 *
 * This is an implementation of the increment-by-one suffix array
 * algorithm described by Dan Gusfield in an unpublished paper.
 * Contact zerkle@cs.ucdavis.edu or gusfield@cs.ucdavis.edu for more
 * information.
 *
 * This generates a suffix array in n log(n) time and linear space.
 * Along the way, it generates Hgt values, which can be used to generate
 * Lcp values.  See the paper for more details on this.
 *
 * Basically, it functions as follows:
 *
 * setup()       reads the input and sets up the various lists and arrays
 *
 * main()        repeatedly calles make_pass() until the suffix array has
 *               been properly generated
 *
 * make_pass()   completes a pass (as described in the paper) by processing
 *               the parents which were acted-upon classes in the last
 *               pass.
 *
 * proc_parent() process the children of this parent, but skip the
 *               largest one.
 *
 * proc_fwd()    process a child before its largest sibling
 *
 * proc_back()   process a child that is after its largest sibling class
 *
 *
 * This program can be compiled to generate debugging information.  If
 * so, a great deal of output will be generated, detailing every single
 * action of the algorithm and reporting the state after each pass.
 * Simply compile with the DEBUG macro set.  If using gcc, this would
 * be "gcc -DDEBUG -O -o suf suf.c".  Otherwise, the only output is
 * the suffix array.
 *
 * The input string is taken on stdin (keyboard or redirection), and
 * must consist only of alphabetical characters.  The puncutation mark
 * "!" is added as an "end of string" character.  This extra suffix is
 * always in position 0 and can be ignored.
 *
 * Each character of input requires 58 bytes of storage (assuming that
 * characters use one byte and integers use 4).  This can be reduced
 * to 57 by compiling with debugging turned off.
 *
 * MAX is the size of the largest possible input.  Raise this number if
 * you want to work with bigger strings.
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "sary_zerkle.h"

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif


/**********************************************************************
 **********************************************************************
 **********                                                  **********
 **********              Structure Definitions               **********
 **********                                                  **********
 **********************************************************************
 **********************************************************************/
/* This structure holds the input */
struct in_struct
  {
    int sufnum;
    char data;
  };

/* This structure holds the info needed for any particular class */
struct class_struct
  {
    int start;		/* Where in Pos the class starts */
    int size;		/* Size of the class */
    int NAL;		/* Next and Last available location */
    int LAL;
    int newstart;	/* Start and size of this class in next pass */
    int newsize;	/* they change as pieces are split off */
    unsigned actedclass:1;    /* True if acted on by this processed class */
    unsigned actedpass:1;     /* True if acted on in this pass */
  };

/* This structure holds the info needed for a parent/acted upon class */
struct parent_struct
  {
    int start;		/* Where in Pos the class starts */
    int size;		/* Size of the class */
  };


/**********************************************************************
 **********************************************************************
 **********                                                  **********
 **********              Function Prototypes                 **********
 **********                                                  **********
 **********************************************************************
 **********************************************************************/

static void startup(char *s, int n, int *posarray);
static int cmp_str(struct in_struct *first,struct in_struct *second);
static void proc_fwd(int classnum);
static void allocspace(void);
static void make_pass(void);
static void proc_parent(int parentnum);
static void setupclasses(void);

#ifdef DEBUG
static void debugoutput();
static void pclass(int classnum);
#endif


/**********************************************************************
 **********************************************************************
 **********                                                  **********
 **********                Global Variables                  **********
 **********                                                  **********
 **********************************************************************
 **********************************************************************/

static int pass;		/* pass number through the data */

static int N;			/* size of input */

static int *Hgt;		/* Hgt array as in paper */

static int *Pos;		/* Pos array as in paper */

#ifdef DEBUG
static char *originput;	/* original input */
#endif

static int *which;		/* which class each suffix belongs to */

static int *where;		/* Where in Pos[] each suffix is */

static struct class_struct *classes;	/* start & size of each class */

/* active parents (acted on last pass) */
static struct parent_struct *parentsd;	/* data */
static int *parents;			/* list of names of parents */
static int numparents;			/* number of them */

/* list of classes acted on this pass, (parents in next pass) */
static int *actedon;			/* list of names */
static int numactedon;			/* number of them */

/* list of classes acted by this class, */
static int *cl_actedon;		/* list of names */
static int numcl_actedon;		/* number of them */

/* name of the next class that will be created */
static int nextname;

/* This indicates what a class should look like as it is acted upon */
static int *move;


/**********************************************************************
 **********************************************************************
 **********                                                  **********
 **********                 Startup Functions                **********
 **********                                                  **********
 **********************************************************************
 **********************************************************************/

/**********************************************************************
 * STARTUP()
 *
 * This function gets the sequence, allocates memory, and sets up all the
 * arrays to the condition they need to be before the first pass.
 **********************************************************************/
void startup(char *s, int n, int *posarray)
{
  struct in_struct *input; /* input data */
  int i;

  N = n;
  Pos = posarray;

  /* Allocate space for input */
  if ((input=(struct in_struct *)
             calloc((size_t)(N+1),sizeof(struct in_struct)))==NULL)
    {
      fprintf(stderr,"Not enough memory for input array\n");
      exit(1);
    } /* if */

#ifdef DEBUG
  puts("Reading input...");
#endif

  /* Get Sequence */
  i=0;
  for (i=0; i < N; i++) {
    input[i].sufnum = i;
    input[i].data = *s++;
  }

  /* Stick on an end-of-string character */
  input[N].sufnum=N;  /* attach the suffix number to this item */
  input[N].data='\0';
  N++;

#ifdef DEBUG
  puts("Allocating space...");
#endif

  /* Allocate space for the wide variety of arrays needed for this mess */
  allocspace();

#ifdef DEBUG
  puts("Setting up originput, Pos, Hgt, and where...");
#endif

#ifdef DEBUG
  /* Copy original data over */
  for(i=0;i<N;i++)
    originput[i]=input[i].data;
#endif

  /* Sort input */
  qsort(input, (size_t)N, sizeof(struct in_struct), (int (*)())cmp_str);

  /* copy input into Pos array and record divisions in Hgt array */
  for(i=0;i<N;i++)
    Pos[i]=input[i].sufnum;

  /* Set up Hgt array */
  Hgt[0]=0;
  for(i=1;i<N;i++)
    if (input[i].data==input[i-1].data)
      Hgt[i]=(-1);
    else
      Hgt[i]=0;  /* Start of new class */

  /* Set up "where" array */
  for (i=0;i<N;i++)
    where[Pos[i]]=i;

#ifdef DEBUG
  puts("Setting up classes...");
#endif

  setupclasses();

  /* Free input array */
  free(input);

#ifdef DEBUG
  puts("Done setting up");
#endif

} /* startup() */


/*************************************************************************
 * SETUPCLASSES()
 *
 * This sets up the various lists to contain the regular and parent classes.
 *************************************************************************/
void setupclasses(void)
{
  int i;

  /* Set up original parent class (all the suffixes) */
  /* Give it the name "0" and make it the whole Pos array */
  parentsd[0].start=0;
  parentsd[0].size=N;

  /* Put said parent class in the list */
  parents[0]=0;
  numparents=1;

  /* Start up the first (regular) class */
  classes[0].start=classes[0].newstart=0;
  nextname=0;
  which[N]=0;

  /* Set up the rest of the classes */
  for (i=1;i<N;i++)
    {
      if (Hgt[i]==0)
        {
	  /* Finish up the previous class ... */
          classes[nextname].size=i-classes[nextname].start;
          classes[nextname].newsize=classes[nextname].size;
	  classes[nextname].NAL=classes[nextname].start;
	  classes[nextname].LAL=i-1;
	  classes[nextname].actedclass=FALSE;
	  classes[nextname].actedpass=FALSE;

	  /* And start in on the next class */
	  nextname++;
	  classes[nextname].start=classes[nextname].newstart=i;
        } /* if */

      /* Indicate which class this suffix is in */
      which[Pos[i]]=nextname;

#ifdef DEBUG
      printf("%d gets %d\n",Pos[i],nextname);
#endif

    } /* for */

  /* Finish up the very last class ... */
  classes[nextname].size=i-classes[nextname].start;
  classes[nextname].newsize=classes[nextname].size;
  classes[nextname].newstart=classes[nextname].start;
  classes[nextname].NAL=classes[nextname].start;
  classes[nextname].LAL=i-1;
  classes[nextname].actedclass=FALSE;
  classes[nextname].actedpass=FALSE;
  nextname++;
} /* setupclasses() */

/*************************************************************************
 * ALLOCSPACE()
 *
 * This allocates all the arrays that will be needed for this stunt.
 *************************************************************************/
void allocspace(void)
{
  /* Allocate space for Hgt array */
  if ((Hgt=(int *) calloc((size_t)N,sizeof(int)))==NULL)
    {
      fprintf(stderr,"Not enough memory for Hgt array\n");
      exit(1);
    } /* if */

#ifdef DEBUG
  /* Allocate space for data array (to contain input in original order) */
  if ((originput=(char *) calloc((size_t)N,sizeof(char)))==NULL)
    {
      fprintf(stderr,"Not enough memory for original data array\n");
      exit(1);
    } /* if */
#endif

  /* Allocate space for "where" array */
  if ((where=(int *) calloc((size_t)N,sizeof(int)))==NULL)
    {
      fprintf(stderr,"Not enough memory for where array\n");
      exit(1);
    } /* if */

  /* Allocate space for "which" array */
  if ((which=(int *) calloc((size_t)N,sizeof(int)))==NULL)
    {
      fprintf(stderr,"Not enough memory for where array\n");
      exit(1);
    } /* if */

  /* Allocate space for parent data */
  if ((parentsd=(struct parent_struct *)
             calloc((size_t)N,sizeof(struct parent_struct)))==NULL)
    {
      fprintf(stderr,"Not enough memory for parent data\n");
      exit(1);
    } /* if */

  /* Allocate space for parent list */
  if ((parents=(int *) calloc((size_t)N,sizeof(int)))==NULL)
    {
      fprintf(stderr,"Not enough memory for parents list\n");
      exit(1);
    } /* if */

  /* Allocate space for "acted on by this pass" list */
  if ((actedon=(int *) calloc((size_t)N,sizeof(int)))==NULL)
    {
      fprintf(stderr,"Not enough memory for acted on list\n");
      exit(1);
    } /* if */

  /* Allocate space for "acted on by this class" list */
  if ((cl_actedon=(int *) calloc((size_t)N,sizeof(int)))==NULL)
    {
      fprintf(stderr,"Not enough memory for class acted on list\n");
      exit(1);
    } /* if */

  /* Allocate space for the move array */
  if ((move=(int *) calloc((size_t)N,sizeof(int)))==NULL)
    {
      fprintf(stderr,"Not enough memory for class acted on list\n");
      exit(1);
    } /* if */

  /* Allocate space for classes data */
  if ((classes=(struct class_struct *)
             calloc((size_t)N,sizeof(struct class_struct)))==NULL)
    {
      fprintf(stderr,"Not enough memory for classes array\n");
      exit(1);
    } /* if */


} /* allocspace() */

/*************************************************************************
 * CMP_STR()
 *
 * This compares two structs and indicates which is greater.
 * Used by the qsort() function in setup().
 *************************************************************************/
int cmp_str(struct in_struct *first,struct in_struct *second)
{
  return first->data - second->data;
} /* cmp_str() */



/**********************************************************************
 **********************************************************************
 **********                                                  **********
 **********              Processing Functions                **********
 **********                                                  **********
 **********************************************************************
 **********************************************************************/


/*************************************************************************
 * PROC_FWD()
 *
 * This processes a class forward.  Any swaps recorded will be with the
 * NAL of the acted upon class.  This is used on the children of a parent
 * class before the largest child.  After that, PROC_BACK() is used.
 *************************************************************************/
void proc_fwd(int classnum)
{
  int j;	/* to step through the elements of the class */
  int start;	/* start of class is Pos */
  int size;	/* size of class in Pos */
  int affsuf;	/* the suffix affected by this suffix */
  int affclass; /* the class affected by this suffix */

  /* Abbreviations for this affected class and newly created class */
  struct class_struct *thisclass, *newclass;

#ifdef DEBUG
  printf("Processing class %d forward\n",classnum);
#endif

  /* Determine start and size of class */
  start=classes[classnum].start;
  size=classes[classnum].size;

  /* Clear the list of classes acted upon by this class to empty */
  numcl_actedon=0;

  /* Go through each element of this class */
  for (j=start;j<size+start;j++)
    /* Don't do anything for suffix 0 */
    if (Pos[j])
      {

#ifdef DEBUG
        printf("Examining suffix %d\n",Pos[j]);
#endif

        /* Determine the number of the suffix at this location in Pos */
        affsuf=Pos[j]-1;

        /* Figure out which class it is in */
        affclass=which[affsuf];

        /* Don't do anything if it is a singleton class */
        if (classes[affclass].size!=1)
          {

            /* Record the needed swap */
	    /* In the next pass, the affected suffix will go in the */
	    /* position it is here placed in the move[] array */
	    move[classes[affclass].NAL]=affsuf;

#ifdef DEBUG
            printf("Affected suffix %d moves to %d in\n",
	           affsuf, classes[affclass].NAL);
            pclass(affclass);
#endif

	    /* Point NAL to the next location in the class */
	    classes[affclass].NAL++;

	    /* Indicate that this class has been "acted on" by... */
	    /* ...this processed class */
	    if (!classes[affclass].actedclass)
	      {
	        /* Set the "acted on" flag */
		classes[affclass].actedclass=TRUE;

		/* put it in the list */
                cl_actedon[numcl_actedon]=affclass;
                numcl_actedon++;

	        /* ...and acted on in this pass */
	        if (!classes[affclass].actedpass)
	          {
	            /* Set the "acted on" flag */
		    classes[affclass].actedpass=TRUE;

		    /* put it in the list */
                    actedon[numactedon]=affclass;
                    numactedon++;

	            /* ...and acted on in this pass */
	          } /* if (not yet acted on by this pass) */

	      } /* if (not yet acted on by this class) */

	  } /* if (it's not a singleton class) */

      } /* if (suffix is not suffix 0) */
  /* End of FOR loop going through all suffixes in this class */

  /* Update arrays indicating where all classes, suffices, and so on go */
  /* Pos is not updated until the whole pass is finished, though. */

  /* Indicate creation of new classes */

  /* If NAL got all the way to the end of this class, it means that */
  /* whatever is left of it (after previous splits) will be its own */
  /* class in the next pass.  Carry over the name from this class to */
  /* the new class in the next segment.  This program is set up so that */
  /* such a class need only be touched to remove the flag indicating that */
  /* it was acted on by the last processed class */

  /* If not, this class must be split at the position of NAL */
  /* Give the first half a new name, and let the remainder keep the */
  /* old name.  This new class with the old name will have different */
  /* start and size, so indicate those in the newstart and newsize   */
  /* members of the structure.  These will be changed again if the class */
  /* is split again. */

  for (j=0;j<numcl_actedon;j++)
    {

      thisclass=&(classes[cl_actedon[j]]);

#ifdef DEBUG
      printf("Considering splitting class %d, NAL %d, LAL %d\n",
             cl_actedon[j],thisclass->NAL,thisclass->LAL);
#endif

      /* Do the rest only if splitting to make a class with a new name */
      if (thisclass->LAL!=thisclass->NAL-1)

        {

#ifdef DEBUG
          printf("Actually splitting class %d at %d to make class %d\n",
	         cl_actedon[j],thisclass->NAL,nextname);
#endif

          Hgt[thisclass->newstart]=pass;

	  /* Create the new class */
	  newclass=&(classes[nextname]);

	  newclass->start=newclass->newstart=thisclass->newstart;
	  newclass->NAL=newclass->newstart;
	  newclass->newsize=thisclass->NAL-thisclass->newstart;
	  newclass->size=newclass->newsize;
	  newclass->LAL=newclass->start+newclass->size-1;
	  newclass->actedclass=newclass->actedpass=FALSE;

	  /* Modify the existing class */
	  thisclass->newstart=thisclass->NAL;
	  thisclass->newsize=thisclass->LAL-thisclass->NAL+1;

#ifdef DEBUG
	  pclass(cl_actedon[j]);
	  pclass(nextname);
#endif

	  nextname++;
        } /* if */

      /* Get it ready to be possibly acted on by the next class */
      /* Do this whether splitting or not */
      thisclass->actedclass=FALSE;

    } /* for */


} /* proc_fwd() */


/*************************************************************************
 * PROC_BACK()
 *
 * This processes a class backward.  Any swaps recorded will be with the
 * LAL of the acted upon class.  This is used on the children of a parent
 * class before the largest child.  Before this, proc_fwd is used.
 *
 * Note that suffixes are still examined left-to-right (order doesn't
 * matter in processing a class).  However, the acted-on classes will
 * have suffixes swapped towards the right (LAL), instead of left (NAL).
 *************************************************************************/
void proc_back(int classnum)
{
  int j;	/* to step through the elements of the class */
  int start;	/* start of class is Pos */
  int size;	/* size of class in Pos */
  int affsuf;	/* the suffix affected by this suffix */
  int affclass; /* the class affected by this suffix */

  /* Abbreviations for this affected class and newly created class */
  struct class_struct *thisclass, *newclass;

#ifdef DEBUG
  printf("Processing class %d backward\n",classnum);
#endif

  /* Determine start and size of class */
  start=classes[classnum].start;
  size=classes[classnum].size;

  /* Clear the list of classes acted upon by this class to empty */
  numcl_actedon=0;

  /* Go through each element of this class */
  for (j=start;j<size+start;j++)
    /* Don't do anything for suffix 0 */
    if (Pos[j])
      {

#ifdef DEBUG
        printf("Examining suffix %d\n",Pos[j]);
#endif

        /* Determine the number of the suffix at this location in Pos */
        affsuf=Pos[j]-1;

        /* Figure out which class it is in */
        affclass=which[affsuf];

        /* Don't do anything if it is a singleton class */
        if (classes[affclass].size!=1)
          {

            /* Record the needed swap */
	    /* In the next pass, the affected suffix will go in the */
	    /* Position indicated by the move[] array */
	    move[classes[affclass].LAL]=affsuf;

#ifdef DEBUG
            printf("Affected suffix %d moves to %d in\n",
	           affsuf, classes[affclass].LAL);
            pclass(affclass);
#endif

	    /* Point LAL to the next location in the class */
	    classes[affclass].LAL--;

	    /* Indicate that this class has been "acted on" by... */
	    /* ...this processed class */
	    if (!classes[affclass].actedclass)
	      {
	        /* Set the "acted on" flag */
		classes[affclass].actedclass=TRUE;

		/* put it in the list */
                cl_actedon[numcl_actedon]=affclass;
                numcl_actedon++;

	        /* ...and acted on in this pass */
	        if (!classes[affclass].actedpass)
	          {
	            /* Set the "acted on" flag */
		    classes[affclass].actedpass=TRUE;

		    /* put it in the list */
                    actedon[numactedon]=affclass;
                    numactedon++;

	            /* ...and acted on in this pass */
	          } /* if (not yet acted on by this pass) */

	      } /* if (not yet acted on by this class) */

	  } /* if (it's not a singleton class) */

      } /* if (suffix is not suffix 0) */
  /* End of FOR loop going through all suffixes in this class */

  /* Update arrays indicating where all classes, suffices, and so on go */
  /* Pos is not updated until the whole pass is finished, though. */

  /* Indicate creation of new classes */

  /* If LAL got all the way to the beginning of this class, it means that */
  /* whatever is left of it (after previous splits) will be its own */
  /* class in the next pass.  Carry over the name from this class to */
  /* the new class in the next segment.  This program is set up so that */
  /* such a class need only be touched to remove the flag indicating that */
  /* it was acted on by the last processed class */

  /* If not, this class must be split at the position of LAL */
  /* Give the right half a new name, and let the remainder keep the */
  /* old name.  This new class with the old name will have different */
  /* start and size, so indicate those in the newstart and newsize   */
  /* members of the structure.  These will be changed again if the class */
  /* is split again. */

  for (j=0;j<numcl_actedon;j++)
    {

      thisclass=&(classes[cl_actedon[j]]);

#ifdef DEBUG
      printf("Considering splitting class %d, NAL %d, LAL %d\n",
             cl_actedon[j],thisclass->NAL,thisclass->LAL);
#endif

      /* Do the rest only if splitting to make a class with a new name */
      if (thisclass->LAL!=thisclass->NAL-1)

        {

#ifdef DEBUG
          printf("Actually splitting class %d at %d to make class %d\n",
	         cl_actedon[j],thisclass->LAL+1,nextname);
          puts("Before: ");
	  pclass(cl_actedon[j]);
	  puts("After:");
#endif

          Hgt[thisclass->LAL+1]=pass;

	  /* Create the new class */
	  newclass=&(classes[nextname]);

	  newclass->start=newclass->newstart=thisclass->LAL+1;
	  newclass->NAL=newclass->newstart;
	  newclass->newsize=thisclass->newstart+thisclass->newsize
	                    -thisclass->LAL-1;
	  newclass->size=newclass->newsize;
	  newclass->LAL=newclass->start+newclass->size-1;
	  newclass->actedclass=newclass->actedpass=FALSE;

	  /* Modify the existing class */
	  thisclass->newsize=thisclass->LAL-thisclass->NAL+1;

#ifdef DEBUG
	  pclass(cl_actedon[j]);
	  pclass(nextname);
#endif

	  nextname++;
        } /* if */

      /* Get it ready to be possibly acted on by the next class */
      /* Do this whether splitting or not */
      thisclass->actedclass=FALSE;

    } /* for */


} /* proc_fwd() */


/*************************************************************************
 * MAKE_PASS()
 *
 * This processes all the needed classes at this stage.  Basically, it
 * initializes some lists, processes all the parents who were acted on
 * in the last pass, then swaps the suffixes around in the Pos[], where[]
 * and which[] arrays.
 *************************************************************************/
void make_pass(void)
{
  int i;
  int j;
  int classstart,classend,classsize,classNAL,classLAL;
  int firstnew;
  int *listtmp;
  struct class_struct *thisclass;

#ifdef DEBUG
  printf("Making pass %d\n",pass);
#endif

  /* Get the name of what the next new class will be */
  firstnew=nextname;

  /* Initialize list of "acted upon" classes to empty */
  numactedon=0;

  /* Run through all the parents that were acted on in the last pass */
  for (i=0;i<numparents;i++)
    proc_parent(parents[i]);  /* Process all children but largest */

  /* Move all the suffixes into their new positions */
  for (i=0;i<numactedon;i++)  /* for each acted-on class */
    {
      classstart=classes[actedon[i]].start;
      classsize=classes[actedon[i]].size;
      classNAL=classes[actedon[i]].NAL;
      classLAL=classes[actedon[i]].LAL;
      classend=classstart+classsize-1;

      /* Move suffixes swapped with NAL to the beginning of the class */
      for (j=classstart;j<classNAL;j++)
        {
	  /* Move out the current resident */
	  Pos[where[move[j]]]=Pos[j];
	  where[Pos[j]]=where[move[j]];

	  /* and move in the new one */
	  Pos[j]=move[j];
	  where[Pos[j]]=j;
	} /* for */

      /* Now move suffixes swapped with LAL to the end */
      for (j=classend;j>classLAL;j--)
        {
	  /* Move out the current resident */
	  Pos[where[move[j]]]=Pos[j];
	  where[Pos[j]]=where[move[j]];

	  /* and move in the new one */
	  Pos[j]=move[j];
	  where[Pos[j]]=j;
	} /* for */

    } /* for */

  /* Update which[] for all suffixes that are elements of new classes */
  for (i=firstnew;i<nextname;i++)
    for (j=classes[i].start;j<classes[i].start+classes[i].size;j++)
      which[Pos[j]]=i;

#ifdef DEBUG
  printf("Total %d classes acted on in this pass\n",numactedon);
#endif

  /* Update all the "acted on" classes, and record them as parent classes */
  for (i=0;i<numactedon;i++)
    {
      thisclass=&classes[actedon[i]];

      parentsd[actedon[i]].start=thisclass->start;
      parentsd[actedon[i]].size=thisclass->size;

      thisclass->start=thisclass->NAL=thisclass->newstart;
      thisclass->size=thisclass->newsize;
      thisclass->LAL=thisclass->start+thisclass->size-1;
      thisclass->actedpass=FALSE;
    } /* for */

  /* make "acted on" list into parent list for next pass */
  listtmp=actedon;
  actedon=parents;
  parents=listtmp;
  numparents=numactedon;
} /* make_pass() */

/*************************************************************************
 * PROC_PARENT()
 *
 * This processes all the children of the given parent, except the largest
 * one.  It goes left->right until the largest one, then right->left
 * (from the other side) until the largest.
 *************************************************************************/
void proc_parent(int parentnum)
{
  int largechild;	/* Largest child of this parent */
  int largesize;	/* size of said child */
  int pstart, psize;	/* copies of this parent's start and size */
  int cclass;		/* name of a child class */
  int p;		/* Position of a child's element in Pos */

  /* Get copies, for convenience */
  pstart=parentsd[parentnum].start;
  psize=parentsd[parentnum].size;

#ifdef DEBUG
  printf("Processing parent %d, starts: %d, size %d\n",
          parentnum,pstart,psize);
#endif

  /* Find out which child is largest */
  largechild=(-1);
  largesize=0;
  p=pstart;

  /* Look at each of the children */
  while (p<pstart+psize)
    {
      /* find class that starts at that position */
      cclass=which[Pos[p]];

#ifdef DEBUG
      printf("child %d has size %d\n",cclass,classes[cclass].size);
#endif

      /* See if it's bigger than the biggest so far */
      if (classes[cclass].size>largesize)
        {
	  largechild=cclass;
	  largesize=classes[cclass].size;
	} /* if */

      /* Skip to next child class */
      p+=classes[cclass].size;

    } /* while */

#ifdef DEBUG
  printf("Largest child of parent %d is %d\n",parentnum,largechild);
#endif

  /* process all the children forward to the largest class */
  p=pstart;
  while (which[Pos[p]]!=largechild)
    {
      /* process this class */
      proc_fwd(which[Pos[p]]);

      /* skip to next class */
      p+=classes[which[Pos[p]]].size;
    } /* while */

  /* Now, process children backward to largest class */
  /* Find the last class */
  cclass=which[Pos[pstart+psize-1]];
  while (cclass!=largechild)
    {
      proc_back(cclass);

      /* go to previous class */
      cclass=which[Pos[classes[cclass].start-1]];
    } /* while */


} /* proc_parent() */

/*************************************************************************
 *************************************************************************
 **********                                                     **********
 **********                     Main Function                   **********
 **********                                                     **********
 *************************************************************************
 *************************************************************************/


int zerkle(char *seq, int n, int *posarray)
{
  int i;

  startup(seq, n, posarray);

#ifdef DEBUG
  putchar('\n');

  puts("Original input");
  for (i=0;i<N;i++)
    printf("  %c",originput[i]);

  putchar('\n');

  puts("POS array");
  for (i=0;i<N;i++)
    printf("%3d",Pos[i]);

  putchar('\n');

  puts("First character of suffixes indicated by Pos");
  for (i=0;i<N;i++)
    printf("%3c",originput[Pos[i]]);

  putchar('\n');

  puts("Whichclass array:");
  for (i=0;i<N;i++)
    printf("%3d",which[i]);

  putchar('\n');

  puts("Hgt array");
  for (i=0;i<N;i++)
    printf("%3d",Hgt[i]);

  putchar('\n');
  putchar('\n');

  puts("Class positions:");
  for(i=0;i<N;i++)
    printf("%3d",which[Pos[i]]);

  putchar('\n');

  for (i=0;i<nextname;i++)
	pclass(i);
#endif /* DEBUG */

  /* Find the suffix array */
  pass=1;
  while (nextname<N)
    {
      make_pass();
#ifdef DEBUG
      printf("Result of Pass %d:\n",pass);
      debugoutput();
#endif
      pass++;
    }

  /*
   * Pos[0] is the position of the end-of-string character, so
   * set the suffix_array to the array starting at Pos[1].
   *
   * The suffices are computed as 0..N-1 (plus end-of-string N), so
   * we must shift them to 1..N (and ignore the end-of-string).
   */
  for (i=1; i <= N; i++)
    Pos[i]++;

  /* Output the Pos array 
  for (i=0;i<N;i++)
    printf(" %2d",Pos[i]);

  putchar('\n');
  */

  free(Hgt);
#ifdef DEBUG
  free(originput);
#endif
  free(where);
  free(which);
  free(parentsd);
  free(parents);
  free(actedon);
  free(cl_actedon);
  free(move);
  free(classes);

  return 1;
} /* sarray_zerkle_build() */


/*************************************************************************
 *************************************************************************
 **********                                                     **********
 **********                  Debugging Functions                **********
 **********                                                     **********
 *************************************************************************
 *************************************************************************/
#ifdef DEBUG
void debugoutput(void)
{
  int i;

  printf("     ");
  for (i=0;i<N;i++)
    printf("%3d",i);
  putchar('\n');

  printf("INPUT");
  for (i=0;i<N;i++)
    printf("  %c",originput[i]);
  putchar('\n');

  printf("POS  ");
  for (i=0;i<N;i++)
    printf("%3d",Pos[i]);
  putchar('\n');

  printf("WHICH");
  for (i=0;i<N;i++)
    printf("%3d",which[i]);
  putchar('\n');

  printf("CLASS");
  for(i=0;i<N;i++)
    printf("%3d",which[Pos[i]]);
  putchar('\n');

}

void pclass(int classnum)
{
  struct class_struct *clpt;

  clpt=&classes[classnum];

  printf("Class %d: st%d si%d NAL%d LAL%d nst %d nsi%d ac%d ap%d\n",
         classnum, clpt->start, clpt->size, clpt->NAL, clpt->LAL,
         clpt->newstart, clpt->newsize, clpt->actedclass, clpt->actedpass);
} /* pclass() */
#endif
