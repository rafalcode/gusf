#
# Makefile
#
# NOTES:
#
#    8/94  -  Created the Makefile and added the Boyer-Moore files.
#             (James Knight)
#    8/94  -  Added more.c and its appropriate dependencies.  (James Knight)
#    8/94  -  Added the main.c, menu.c, util.c, seq_ary.c,
#             file_access.c, and James Knights' my_getline.c
#             Also added Sean Davis' kmp_sp.c, kmp.c, and z.c.
#             (Gene Fodor) 
#    9/94  -  Added stree.c.  (James Knight)
#    9/94  -  Added stree_ukkonen.c, stree_weiner.c, stree_print.c,
#             stree_match.c, sary_qsort.c, sary_zerkle.c, ac.c, 
#             and stree_to_sary.c
#             (Gene Fodor)
#    8/95  -  Reorganized things and corrected all of the dependencies.
#             (James Knight)
#    3/96  -  Updated to reflect the reorganized files  (James Knight)
#    7/96  -  Updated again to reflect the finished reorganization
#             (James Knight)
#   10/97  -  Added the inclusion of SAMPLE_SEQUENCES in tarfile (Jens Stoye)
#    1/98  -  Added tandem.[ch] and its dependencies (Jens Stoye)
#    6/98  -  Added stree_decomposition.[ch] (Jens Stoye)
#    7/98  -  Renamed stree_supermax.[ch] to repeats_supermax.[ch].
#             Renamed tandem.[ch] to repeats_tandem.[ch];
#             Added repeats_primitives.[ch], repeats_nonoverlapping.[ch],
#             repeats_vocabulary.[ch], repeats_linear_occs.[ch],
#             and strmat_stubs4.[ch] (Jens Stoye)
#    8/98  -  Added repeats_maxgap.[ch] (Jens Stoye)
#    2/99  -  Renamed repeats_maxgap.[ch] to repeats_bigpath.[ch];
#             Removed some small bugs in various modules (Jens Stoye)
#

#
# The flags for using the gcc compiler (with flags for either
# debugging and error checking the program or optimizing the code).
#

CC=gcc

#CFLAGS= -g -Wall -Wshadow -pedantic -DSTRMAT -DSTATS
CFLAGS= -O3 -DSTRMAT -DSTATS


#
# The source files, object files, libraries and executable name.
#
SRCFILES= strmat.c \
          ac.c bm.c bmset.c bmset_naive.c kmp.c more.c naive.c \
          sary.c sary_match.c sary_zerkle.c \
          stree_strmat.c stree_ukkonen.c stree_weiner.c stree_lca.c \
          stree_decomposition.c \
          repeats_primitives.c repeats_supermax.c repeats_nonoverlapping.c \
          repeats_bigpath.c repeats_tandem.c repeats_vocabulary.c \
          repeats_linear_occs.c \
          strmat_alpha.c strmat_fileio.c strmat_match.c \
          strmat_print.c strmat_seqary.c strmat_stubs.c strmat_stubs2.c \
          strmat_stubs3.c strmat_stubs4.c strmat_util.c z.c 

OBJFILES= strmat.o \
          ac.o bm.o bmset.o bmset_naive.o kmp.o more.o naive.o \
          sary.o sary_match.o sary_zerkle.o \
          stree_strmat.o stree_ukkonen.o stree_weiner.o stree_lca.o \
          stree_decomposition.o \
          repeats_primitives.o repeats_supermax.o repeats_nonoverlapping.o \
          repeats_bigpath.o repeats_tandem.o repeats_vocabulary.o \
          repeats_linear_occs.o \
          strmat_alpha.o strmat_fileio.o strmat_match.o \
          strmat_print.o strmat_seqary.o strmat_stubs.o strmat_stubs2.o \
          strmat_stubs3.o strmat_stubs4.o strmat_util.o z.o 

LIBS= 

EXECFILE= strmat

#
# The make rule for the executable
#
$(EXECFILE) : $(OBJFILES)
	$(CC) $(CFLAGS) -o $(EXECFILE) $(OBJFILES) $(LIBS)

#
# The dependencies for each of the *.o files.
#

ac.o: ac.h
bm.o: bm.h z.h
bmset_naive.o: bm.h bmset_naive.h
bmset.o: ac.h stree_strmat.h stree_ukkonen.h bmset.h
kmp.o: kmp.h z.h
naive.o: naive.h
z.o: z.h

stree_strmat.o: stree_strmat.h
stree_ukkonen.o: strmat.h stree_strmat.h stree_ukkonen.h
stree_weiner.o: strmat.h stree_strmat.h stree_weiner.h

stree_lca.o: stree_strmat.h stree_lca.h
stree_decomposition.o: stree_strmat.h more.h stree_decomposition.h

repeats_primitives.o: stree_strmat.h more.h repeats_primitives.h
repeats_supermax.o: stree_strmat.h repeats_supermax.h
repeats_nonoverlapping.o: stree_strmat.h more.h repeats_nonoverlapping.h
repeats_bigpath.o: stree_strmat.h more.h repeats_bigpath.h
repeats_tandem.o: stree_strmat.h more.h repeats_tandem.h
repeats_vocabulary.o: stree_strmat.h more.h repeats_vocabulary.h
repeats_linear_occs.o: stree_strmat.h more.h repeats_vocabulary.h \
                       repeats_linear_occs.h

sary.o: sary_zerkle.h sary.h
sary_match.o: strmat.h stree_strmat.h stree_ukkonen.h sary_match.h
sary_zerkle.o: sary_zerkle.h

more.o : more.h
strmat.o: strmat_alpha.h strmat_seqary.h strmat_util.h \
          strmat_print.h \
          strmat_stubs.h strmat_stubs2.h strmat_stubs3.h strmat_stubs4.h \
          sary_match.h stree_ukkonen.h \
          strmat.h
strmat_alpha.o: strmat.h strmat_alpha.h
strmat_fileio.o: strmat.h strmat_alpha.h strmat_fileio.h
strmat_match.o: strmat.h strmat_match.h
strmat_print.o: strmat.h strmat_alpha.h stree_strmat.h strmat_print.h
strmat_seqary.o : strmat.h strmat_seqary.h
strmat_stubs.o: strmat.h strmat_match.h naive.h bm.h bmset_naive.h ac.h \
                kmp.h z.h bmset.h strmat_stubs.h
strmat_stubs2.o: strmat.h strmat_match.h stree_ukkonen.h stree_weiner.h \
                 stree_lca.h stree_decomposition.h strmat_stubs2.h
strmat_stubs3.o: strmat.h stree_ukkonen.h sary.h sary_match.h strmat_stubs3.h
strmat_stubs4.o: strmat.h strmat_match.h stree_ukkonen.h \
                 repeats_primitives.h repeats_supermax.h \
                 repeats_nonoverlapping.h repeats_bigpath.h repeats_tandem.h \
                 repeats_vocabulary.h repeats_linear_occs.h strmat_stubs4.h
strmat_util.o : strmat.h strmat_seqary.h strmat_fileio.h strmat_alpha.h \
                strmat_util.h


#
# Other Standard make rules
#

lint : 
	lint $(SRCFILES) | more

clean:
	rm -f *.o *~

remove: clean
	rm -f $(EXECFILE)

tarfile:
	mkdir strmat.dir
	cp -rf *.c *.h Makefile README FILES SAMPLE_SEQUENCES doc strmat.dir
	tar cvf strmat.tar strmat.dir
	gzip strmat.tar
	/bin/rm -r strmat.dir

