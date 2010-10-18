/**
 * 
 * This code is designed to be an interface to existing C code that
 * implements the RMA program. It is designed to be easily accessible
 * from within the C++ environment, namely wxWidgets.
 * 
 * One of the principal data structures needed for executing RMA is the
 * RMA_FLAGS structure, which holds all of the important settings. These
 * are implemented below in java bean style to make C/C++ interface simpler.
 * 
 */

#include <stdlib.h>
#include "affy.h"
#include "affy_rma.h"
#include "utils.h"

/** The flags that hold RMA settings */
static AFFY_RMA_FLAGS f;
/** The output file to write rma'd expressions to */
static char *output_file="exprs-rma.txt";
 
/**
 * Execute the rma function using the set RMA_FLAGS values. Results are
 * written to the output_file.
 */
void ri_callRMA(char **files) 
{
  AFFY_CHIPSET *results;
  AFFY_ERROR *err;
  char buf[1024];
  err=affy_get_default_error();
  info("Starting rma (wd=%s)\n",getcwd(buf,1024));
  results=affy_rma(files, &f,err);
  info("Writing output to %s..",output_file);
  affy_write_expressions(results,output_file,err);
  info("done.\n");
}

 
/**
 * Initialization routine.
 * NOTE: This must be called prior to using the bean-type routines for
 * RMA_Flags.
 */
void ri_init() 
{
	affy_rma_set_defaults(&f);
}


/**
 * Return background correction flag.
 */
int ri_getBackground()
{
	return f.use_background_correction;
}
/**
 * Set background correction flag.
 */
void ri_setBackground(int b)
{
	f.use_background_correction=b;
}


int ri_isAFFXProbeNormalization() 
{
	return f.normalize_affx_probes;
}

void ri_setAFFXProbeNormalization(int b) 
{
	f.normalize_affx_probes=b;	
}
/**
 * Return one of "Quantile" "None" or "Mean"
 */
char *ri_getNormalization()
{
	if ( ! f.use_normalization ) return "None";
	if ( f.use_mean_normalization ) return "Mean";

	return "Quantile";
}
/**
 * Set normalization based on text input
 */
void ri_setNormalization(const char *text)
{
	if ( STREQ("None",text) ) {
		f.use_normalization=false;
	} else if ( STREQ("Mean",text) ) {
		f.use_normalization=true;
		f.use_mean_normalization=true;
	} else {
		f.use_normalization=true;
		f.use_mean_normalization=false;
	}
}

void ri_setCDFDirectory(const char *text)
{
	f.cdf_directory=strdup(text);
}
char *ri_getCDFDirectory()
{
	return f.cdf_directory;
}

/**
 * Return output file name
 */
char *ri_getOutputFile()
{
	return output_file;
}
/**
 * Set output file name
 */
void ri_setOutputFile(const char *f)
{
	output_file=strdup(f);
}
