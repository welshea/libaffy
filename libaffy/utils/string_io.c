#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define MEM_OVERHEAD 1.01    /* speed hack -- overallocate to avoid reallocs */


/* we can't just pass strcmp to qsort, it needs a wrapper */
int compare_string(const void *f1, const void *f2)
{
  char **cf1 = (char **)f1;
  char **cf2 = (char **)f2;

  return (strcmp(*cf1, *cf2));
}


/* realloc input string, store new max array length (including NULL) */
/* handles \r\n \n \r, including mixes of EOL characters within same file */
/* strips EOL from end of string */
char * fgets_strip_realloc(char **return_string, int *return_max_length,
			   FILE *infile)
{
    char c;
    char *string = *return_string;
    int length = 0;
    int total_length;
    int max_length = *return_max_length;
    char old_c = '\0';
    int anything_flag = 0;

    while((c = fgetc(infile)) != EOF)
    {
    	anything_flag = 1;
    
    	/* EOL: \n or \r\n */
    	if (c == '\n')
    	{
    	    /* MSDOS, get rid of the previously stored \r */
    	    if (old_c == '\r')
    	    {
    	    	string[length - 1] = '\0';
    	    }

    	    old_c = c;
    	    
    	    break;
    	}
    	/* EOL: \r */
    	/* may be a Mac text line, back up a character */
    	else if (old_c == '\r')
    	{
    	    /* fseek(infile, -1 * sizeof(char), SEEK_CUR); */
    	    ungetc(c, infile);

    	    break;
    	}
    	
    	old_c = c;
    
    	length++;
    	total_length = length + 1;

    	if (total_length > max_length)
    	{
    	    max_length = MEM_OVERHEAD * total_length;
    	    string = (char *) realloc(string, max_length * sizeof(char));
    	}

    	string[length-1] = c;
    }
    
    /* check for dangling \r from reading in Mac lines */
    if (length && string[length-1] == '\r')
    {
    	string[length-1] = '\0';
    }
    
    if (length == 0)
    {
    	if (!anything_flag)
    	{
    	    return NULL;
    	}
    	
    	if (1 > max_length)
    	{
    	    string = (char *) realloc(string, sizeof(char));
    	    max_length = 1;
    	}
    }

    string[length] = '\0';
    
    *return_string = string;
    *return_max_length = max_length;

    return string;
}


/* replace tabs with NULLs, fill array of field pointers,
 * return number of fields
 * WARNING -- clobbers tabs in original input string
 */
int split_tabs(char *string, char ***fields, int *return_max_field)
{
    char *cptr, *sptr;
    int count = 0;
    int max_field = *return_max_field;
    
    sptr = string;
    for (cptr = string; *cptr; cptr++)
    {
    	if (*cptr == '\t')
    	{
    	    count++;

    	    if (count > max_field)
    	    {
    	    	max_field = MEM_OVERHEAD * count;
    	    	*fields = realloc(*fields, max_field * sizeof(char *));
    	    }

    	    (*fields)[count-1] = sptr;
    	    sptr = cptr + 1;
    	    *cptr = '\0';
    	}
    }
    
    /* final field */
    count++;
    if (count > max_field)
    {
        max_field = MEM_OVERHEAD * count;
        *fields = realloc(*fields, max_field * sizeof(char *));
    }
    (*fields)[count-1] = sptr;
    
    *return_max_field = max_field;
    
    return count;
}
