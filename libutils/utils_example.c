
/*
  utils_example

  A simple program to demonstrate (and verify) that the utils header
  file works correctly.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "argp.h"

int debug_level = 2;
bool test_arguments = false;

const char *argp_program_version="utils_example v3.0";
const char *argp_program_bug_address="<Eric.Welsh@moffitt.org>";
static struct argp_option options[] = {
  {"verbose", 'v', 0, 0, "Produce verbose output" },
  {"test-arguments", 't', 0, 0, "Test arguments" },
  {0}
};
static error_t
parse_opt(int key, char *arg, struct argp_state *state);
static struct argp argp={options, parse_opt,
	0,
	"libutils example program - demonstrate the use of libutils"
};

int main(int argc, char **argv)
{
  bool bvalue = true;
  char buf[MAXBUF];
  int n = 100;
  int i1 = 10, i2 = 20, imax, imin;
  double d1 = 10.1, d2 = 10.2, dmax, dmin;
  char *p;
  FILE *fp;

  argp_parse(&argp, argc, argv, 0, 0, 0);

  printf("------------------------------------\n"
         "Example and test program for utils.h\n");

  if ( test_arguments == true ) {
     printf("Argument processing ok.\n");
  }

  if (bvalue == true)
    printf("Boolean variable ok.\n");
  printf("MAXBUF ok.\n");

  sprintf(buf, "foo");
  if (STREQ(buf, "foo"))
    printf("STREQ ok.\n");

  printf("Printing string %s and int %d using various methods.\n", buf, n);
  info(("Using info for printing (%s,%d)\n", buf, n));
  warn(("Using warn for printing (%s,%d)\n", buf, n));
  debug((1, "Using debug for printing (%s,%d)\n", buf, n));
  if (fork() == 0)
    die(("Die function for printing (%s,%d)\n", buf, n));
  else
    wait(NULL);

  printf("Checking max and min...\n");
  imax = max_macro(i1, i2);
  imin = min_macro(i1, i2);
  printf("For %d and %d, min is %d, max is %d\n", i1, i2, imin, imax);
  dmax = max_macro(d1, d2);
  dmin = min_macro(d1, d2);
  printf("For %f and %f, min is %f, max is %f\n", d1, d2, dmin, dmax);

  printf("\nChecking alloc functions (p is NULL)...\n");
  p = NULL;
  p = MALLOC(1 * sizeof(char));
  printf(" p after malloc is %p\n", p);

  p = NULL;
  p = CALLOC(1, sizeof(char));
  printf(" p after calloc is %p\n", p);

  p = REALLOC(p, 100 * sizeof(char));
  printf(" p after realloc is %p\n", p);

  printf("\nFOPEN check on noexistent file:\n");
  if (fork() == 0)
    fp = FOPEN("foobar.nowhere", "r");
  else
    wait(NULL);

  sprintf(buf, "This is a test line\n");
  printf("\nThis line has a new line in it: \n%s", buf);
  chomp(buf);
  printf("After chomping: %s.\n", buf);

  printf("\n\nChecks complete.\n"
         "------------------------------------\n");

  return 0;
}

static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
   switch (key) {
	case 't':
	   test_arguments=true;
	   break;
    	case 'v':
	   info(("Verbose mode enabled\n"));
	   break;
	default:
	   return ARGP_ERR_UNKNOWN;
   }
   return 0;
} 
