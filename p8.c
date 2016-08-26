/* The first of a series of Dan Gusfield code snippets */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define IOF 1 /* index offset, 0 or 1 */

int main(int argc, char *argv[])
{
   char *s="aabaabcaxaabaabcy";
   printf("%s\n", s); 

   printf("Pstring length is: %zu\n", strlen(s));
   int i=10;
   int zi=7;
   printf("Z_%d is %d, str \"%.*s\"\n", i, zi, zi, s+i-IOF);
   int j=15;
   int zj=2;
   printf("Z_%d is %d, str \"%.*s\"\n", j, zj, zj, s+j-IOF);

   return 0;
}
