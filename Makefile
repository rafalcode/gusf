CC=gcc
CFLAGS=-g -Wall
EXECUTABLES=p8 a.out

# small tiny shapes for testing on eigbirds.
p8: p8.c
	${CC} ${CFLAGS} -o $@ $^


# this prog calculates the Z for each position in the string.
# the first value of course is trivial and, and technical not part of the Gusfield procedure.
# Z=0 is also an uninteresting value, and so may be left out. We're therefore interested in positions
# where Z has a value.
z0: z0.c
	${CC} ${CFLAGS} -o $@ $^

.PHONY: clean

clean:
	rm -f ${EXECUTABLES}
