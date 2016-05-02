
typedef struct {
  char *P;
  int M;
  int *R, *B, *L;
} BMOPT_STRUCT;

BMOPT_STRUCT *bmopt_prep(char *P, int M);
char *bmopt_search(BMOPT_STRUCT *node, char *T, char *Tend);
void bmopt_free(BMOPT_STRUCT *node);


