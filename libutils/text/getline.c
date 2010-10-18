
#include "utils.h"

/* Arbitrary initial size for a line, the width of a 'standard' terminal */
#define INITIAL_BUF_SIZE 80

/*
 * Return a pointer to the next line from an input stream, allocating
 * memory as necessary.  The caller is reponsible for any memory
 * allocated.  NULL is returned on EOF, else a pointer to the line.
 *
 * Note that a line is considered to end either with a newline or a
 * carriage return.  Thus both DOS and Unix style line endings should
 * be handled properly.  Unless you embed random carriage returns or
 * newlines in your oppositely formatted text files, in which case you
 * lose I guess.
 */
char *utils_getline(FILE *fp)
{
  char  *nextline;
  int    ch;
  size_t sz = INITIAL_BUF_SIZE, i = 0;

  assert(fp != NULL);

  nextline = (char *)MALLOC(INITIAL_BUF_SIZE);

  while ((ch = fgetc(fp)) != EOF)
  {
    if ((ch == '\n') || (ch == '\r'))
      break;
    else
      nextline[i] = ch;

    /* Check if a realloc is needed. */
    if (++i == sz)
    {
      sz *= 2;
      nextline = REALLOC(nextline, sz);
    }
  }

  nextline[i] = '\0';

  /* If we got at least one character before EOF, return the line. */
  if (i > 0)
    return (nextline);
  else
    return (NULL);
}
