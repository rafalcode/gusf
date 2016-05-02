
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "strmat.h"
#include "strmat_alpha.h"
#include "strmat_fileio.h"



#define TOTAL_ALPHABETS 5
char *alpha_names[] = {"UNKNOWN", "DNA", "RNA", "PROTEIN", "ASCII"};
#define TOTAL_DB_NUMBER 3
char *db_names[] = {"Text", "GenBank", "UNKNOWN"};

#define ASCII_MAX_NUM_CHARS 79
#define NONASCII_NUM_CHARS 60


int find_next_sequence(FILE *fp)
{
  int len, fpos;
  char *s, *line;

  while ((line = my_getline(fp, &len)) != NULL) {
    if (strncmp(line, "TYPE:", 5) == 0) {
      fpos = ftell(fp);
      fseek(fp, fpos - len - 1, 0);
      return OK;
    }
    else {
      for (s=line; *s; s++)
        if (!isspace((int)(*s)))
          return FORMAT_ERROR;
    }
  }
  if (feof(fp))
    return ERROR;
  else
    return FORMAT_ERROR;
}


int read_sequence(FILE *fp, STRING *sptr)
{
  int i, len, other_chars, seqlen;
  char ch, newline_ch, *s, *line, buffer[32];

  /*
   * Skip blank lines.
   */
  while ((line = my_getline(fp, &len)) != NULL && line[0] == '\0') ;
  if (line == NULL)
    return ERROR;

  /*
   * Parse the TYPE line.
   */
  if (sscanf(line, "TYPE: %s", buffer) == 0)
    return FORMAT_ERROR;

  for (i=0; i < TOTAL_DB_NUMBER; i++)
    if (strcmp(buffer, db_names[i]) == 0)
      break;
  if (i == TOTAL_DB_NUMBER)
    return FORMAT_ERROR;
  else
    sptr->db_type = i; 

  /*
   * Parse the IDENT line.
   */
  if ((line = my_getline(fp, &len)) == NULL)
    return PREMATURE_EOF;
  else if (strncmp(line, "IDENT:", 6) != 0) {
    if (strncmp(line, "INDENT:", 7) == 0) {
      fprintf(stderr, "\nWarning:  Please rewrite this input file."
              "  It's file format is outdated.\n");
      sscanf(line, "INDENT: %s", sptr->ident);
    }
    else 
      return FORMAT_ERROR;
  }
  else
    sscanf(line, "IDENT: %s", sptr->ident);

  /*
   * Parse the TITLE line.
   */
  if ((line = my_getline(fp, &len)) == NULL)
    return PREMATURE_EOF;
  else if (strncmp(line, "TITLE:", 6) != 0)
    return FORMAT_ERROR;

  for (s=line+6; *s && isspace((int)(*s)); s++) ;
  strncpy(sptr->title, s, TITLE_LENGTH);
  sptr->title[TITLE_LENGTH] = '\0';

  /*
   * Parse the ALPHABET line.
   */
  if ((line = my_getline(fp, &len)) == NULL)
    return PREMATURE_EOF;
  else if (sscanf(line, "ALPHABET: %s", buffer) == 0)
    return FORMAT_ERROR;

  for (i=0; i < TOTAL_ALPHABETS; i++)
    if (strcmp(buffer, alpha_names[i]) == 0)
      break;
  if (i == TOTAL_ALPHABETS)
    return FORMAT_ERROR;
  else
    sptr->raw_alpha = i;

  newline_ch = rawmapchar(sptr->raw_alpha, '\n');
  if (newline_ch == -1)
    return ERROR;

  /*
   * Parse the LENGTH line.
   */
  if ((line = my_getline(fp, &len)) == NULL)
    return PREMATURE_EOF;
  else if (sscanf(line, "LENGTH: %d", &sptr->length) == 0)
    return FORMAT_ERROR;

  /*
   * Parse the SEQUENCE line.
   */
  if ((line = my_getline(fp, &len)) == NULL)
    return PREMATURE_EOF;
  else if (strcmp(line, "SEQUENCE:") != 0)
    return FORMAT_ERROR;

  /*
   * Allocate the sequence buffers.
   */
  if ((sptr->raw_seq = malloc(sptr->length + 1)) == NULL) {
    fprintf(stderr, "\nRan out of memory.  Unable to store new sequence.\n");
    return ERROR;
  }
  if ((sptr->sequence = malloc(sptr->length + 1)) == NULL) {
    fprintf(stderr, "\nRan out of memory.  Unable to store new sequence.\n");
    return ERROR;
  }

  /*
   * Read the sequence.
   */
  other_chars = 0;
  seqlen = 0;
  while ((line = my_getline(fp, &len)) != NULL && strcmp(line, "//") != 0) {
    for (s=line; *s; s++) {
      ch = rawmapchar(sptr->raw_alpha, *s);
      if (ch == -1)
        return ERROR;
      else if (ch != '\0') {
        if (seqlen == sptr->length) {
          fprintf(stderr, 
                  "\nError: Sequence %s contains too many characters.\n",
                  sptr->title);
          return ERROR;
        }
        sptr->raw_seq[seqlen++] = ch;
      }
      else if (!isspace((int)(*s)))
        other_chars = 1;
    }

    if (newline_ch != '\0' && len < ASCII_MAX_NUM_CHARS &&
        seqlen < sptr->length)
      sptr->raw_seq[seqlen++] = newline_ch;
  }
  if (line == NULL)
    return PREMATURE_EOF;

  if (seqlen < sptr->length) {
    fprintf(stderr, "\nError:  Sequence %s contains too few characters.\n",
            sptr->title);
    return ERROR;
  }

  if (other_chars)
    return MISMATCH;
  else
    return OK;
}


int write_sequence(FILE *fp, STRING *sptr)
{
  int i, j, pos;
  char ch, newline_ch;

  if (sptr->db_type < 0 || sptr->db_type >= TOTAL_DB_NUMBER ||
      sptr->raw_alpha < 0 || sptr->raw_alpha >= TOTAL_ALPHABETS ||
      (newline_ch = rawmapchar(sptr->raw_alpha, '\n')) == -1)
    return FORMAT_ERROR;

  /*
   * Print the header.
   */
  fprintf(fp, "TYPE:  %s\n", db_names[sptr->db_type]);
  fprintf(fp, "IDENT:  %s\n", sptr->ident);
  fprintf(fp, "TITLE:  %s\n", sptr->title);
  fprintf(fp, "ALPHABET:  %s\n", alpha_names[sptr->raw_alpha]);
  fprintf(fp, "LENGTH:  %d\n", sptr->length);
  fprintf(fp, "SEQUENCE:\n");

  /*
   * Print the sequence.
   */
    
  if (newline_ch == '\0') {
    for (i=0; i + NONASCII_NUM_CHARS <= sptr->length; ) {
      for (j=0; j < NONASCII_NUM_CHARS; j++)
        fputc(sptr->raw_seq[i++], fp);
      fputc('\n', fp);
    }
    if (i < sptr->length) {
      while (i < sptr->length)
        fputc(sptr->raw_seq[i++], fp);
      fputc('\n', fp);
    }
  }
  else {
    pos = 0;
    for (i=0; i < sptr->length; i++) {
      ch = sptr->raw_seq[i];
      fputc(ch, fp);
      if (ch == '\n')
        pos = 0;
      else if (++pos == ASCII_MAX_NUM_CHARS) {
        fputc('\n', fp);
        pos = 0;
      }
    }
    fputc('\n', fp);
  }

  fprintf(fp, "//\n");

  return OK;
}

