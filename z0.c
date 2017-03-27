/* Calculates the size of the zbox forThe first of a series of Dan Gusfield code snippets */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{
	int i;
	char *s="aabaabcaxaabaabcy";

	unsigned slen=strlen(s);
	unsigned *zl=calloc(slen, sizeof(unsigned));
	char *sp1, *sp2; // two ptrs to travel along s.

	for(i=0;i<slen;++i) {
		sp1=s+i;
		sp2=s;
		while(*sp1) {
			if(*sp1++ == *sp2++)
				zl[i]++;
			else
				break;
		}
	}

	printf("The target string was \"%s\"\n", s); 
	for(i=0;i<slen;++i) 
		printf("%u ", zl[i]);
	printf("\n"); 
	printf("Note first value is trivial, really should be left out\n"); 

	free(zl);

	return 0;
}
