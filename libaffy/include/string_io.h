/*
 * 04/06/11: Added fgets_strip_realloc() and split_tabs() (EAW)
 * 06/01/18: Added compare_string() function (EAW)
 */

extern char * fgets_strip_realloc(char **return_string, int *return_max_length,
                                  FILE *infile);
extern int split_tabs(char *string, char ***fields, int *return_max_field);
extern int compare_string(const void *f1, const void *f2);
