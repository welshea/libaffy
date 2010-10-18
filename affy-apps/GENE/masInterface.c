/**
 * 
 * This code is designed to be an interface to existing C code that
 * implements the MAS program. It is designed to be easily accessible
 * from within the C++ environment, namely wxWidgets.
 * 
 * One of the principal data structures needed for executing MAS is the
 * MAS_FLAGS structure, which holds all of the important settings. These
 * are implemented below in java bean style to make C/C++ interface simpler.
 * 
 */

#include <stdlib.h>
#include "utils.h"
#include "affy.h"
#include "affy_mas5.h"

/** The flags that hold MAS settings */
static AFFY_MAS5_FLAGS f;
/** The output file to write mas'd expressions to */
static char *output_file="exprs-mas.txt";
 
/**
 * Execute the mas function using the set MAS_FLAGS values. Results are
 * written to the output_file.
 */
void mi_callMAS(char **files) 
{
  AFFY_CHIPSET *results;
  char buf[1024];
  AFFY_ERROR *err;

  err=affy_get_default_error();
  info("Starting mas (wd=%s)\n",getcwd(buf,1024));
  results=affy_mas5(files, &f,err);
  info("Writing output to %s..",output_file);
  affy_write_expressions(results,output_file,err);
  info("done.\n");
}

 
/**
 * Initialization routine.
 * NOTE: This must be called prior to using the bean-type routines for
 * MAS_Flags.
 */
void mi_init() 
{
  affy_mas5_set_defaults(&f);
}


/**
 * Return background correction flag.
 */
int mi_getBackground()
{
  return f.use_background_correction;
}
/**
 * Set background correction flag.
 */
void mi_setBackground(int b)
{
  f.use_background_correction=b;
}

int mi_getQuantileNormalization()
{
  return f.use_quantile_normalization;
}
void mi_setQuantileNormalization(int b) 
{
  f.use_quantile_normalization=b;
}

int mi_getBioconductorCompatability() 
{
   return f.bioconductor_compatability;
}
void mi_setBioconductorCompatability(int b)
{
   f.bioconductor_compatability=b;
}
int mi_getMeanNormalization()
{
  return f.use_mean_normalization;
}

void mi_setMeanNormalization(int b) 
{
	f.use_mean_normalization=b;
}

int mi_getMeanNormalizationValue()
{
	return f.mean_normalization_target_mean;
}
void mi_setMeanNormalizationValue(int val)
{
	f.mean_normalization_target_mean=val;
}

int mi_getScaleProbesets()
{
	return f.use_probeset_scaling;
}
void mi_setScaleProbesets(int b)
{
	f.use_probeset_scaling=b;
}
int mi_getScaleProbesetsValue()
{
	return f.scale_target;
}
void mi_setScaleProbesetsValue(int val)
{
	f.scale_target=val;
}


void mi_setCDFDirectory(const char *text)
{
	f.cdf_directory=strdup(text);
}
char *mi_getCDFDirectory()
{
	return f.cdf_directory;
}

/**
 * Return output file name
 */
char *mi_getOutputFile()
{
	return output_file;
}
/**
 * Set output file name
 */
void mi_setOutputFile(const char *f)
{
	output_file=strdup(f);
}
