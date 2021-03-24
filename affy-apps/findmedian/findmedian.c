
/**************************************************************************
 *
 * Filename: findmedian.c
 * 
 * Purpose:  Utility to find the nearest chip to a population of chips
 *
 * Creation: 14 March, 2011
 *
 * Author:   Eric A. Welsh
 *
 *
 * Update History
 * --------------
 * 03/14/11: File creation (EAW)
 * 05/13/13: added support for --salvage
 * 05/13/13: runs correctly on a single chip now
 * 06/13/13: fixed --spreadsheet regression introduced by --salvage code
 * 09/12/13: fixed --probeset regression introduced by --salvage code
 * 09/12/13: added support for --log2 --nolog2 --unlog2 options
 * 02/14/14: divide distances by normalized number of points in comparison
 * 03/06/14: floor CEL file data at 1 to prevent NaNs
 * 03/06/14: --ignore-weak flag is now default, use --include-weak to include weak
 * 05/15/14: include number of non-weak points when normalizing the number
 *           of points (see 2/14/14 changes)
 * 10/24/17: fixed memory allocation segfault in --spreadsheet mode when #cols > #rows
 * 01/10/18: fixed --meancenter so it actually performs mean centering,
 *             rather than doing nothing.
 *           changed missing data handling so that --meancenter and --pearson
 *             work correctly with --ignore-weak
 *           changed --pearson distance from [1-r] to [sqrt(0.5*(1-r))],
 *             which satisfies triangle inequality
 * 05/15/18: fixed a memory pool/free-on-exit issue, no obvious user impact
 * 08/29/18: added -x probeset exclusion list for --spreadsheet mode
 * 09/11/18: extend -x probeset exclusion list to CEL files as well
 * 09/14/18: add -s spikein support for flag compatability with iron_generic
 * 10/22/18: fix segfaults due to not checking exclude flags before excluding
 * 10/23/18: add global brightness metrics, add a header line to the output
 * 03/13/19: skip AFFX/control probes/probesets
 * 04/08/19: fix potential sqrt(-#) metric pearson distance roundoff error
 * 05/22/19: fix --probesets crash with new control probeset checking
 * 05/22/19: added --ignore-chip-mismatch support
 * 08/13/19: fixed printing of avg_rmsd to not print the pointer address
 * 04/01/20: fixed segfault with --include-weak introduced 10/23/18
 * 05/22/20: add -o option to write to a file, rather than STDOUT
 * 07/10/20: fopen() as "w" rather than "wt", as "wt" is not standard syntax
 * 08/03/20: add output buffering and #include <stdlib.h> in attempt to make
 *           file output happy when compiled with emcc (Emscripten), neither
 *           should matter or make any difference...
 * 08/12/20: move Average RMSD line, print additional median sample line;
 *           should make it easier to parse median sample, while still
 *           maintaining backwards comptability (original line is still the
 *           2nd to last line).  It only breaks anything that required the
 *           Average RMSD line to be the last line (it is now 3rd to last).
 * 08/12/20: disable searching current working directory for CEL files (EAW)
 * 08/12/20: pass flags to affy_create_chipset() (EAW)
 * 03/24/21: corrected --ignore-weak description (ignore <= 0, not <= 1)
 * 03/24/21: add progress indicator for missing data pre-scan
 *
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "affy.h"
#include "affy_mas5.h"
#include "affy_rma.h"
#include "argp.h"

#define MIN_VALUE 1
#define MISSING -FLT_MAX
#define SEARCH_WORKING_DIR 0

AFFY_COMBINED_FLAGS   flags;
char                 *directory     = ".";
char                 *cdf_directory = ".";
char                **filelist;
char                 *outfile_name = NULL;
int                  *mempool = NULL;
int                   debug_level = 2;
int                   opt_distance_rmsd_flag = 1;
int                   opt_distance_pearson_flag = 0;
int                   opt_mean_center_flag = 0;
int                   opt_spreadsheet_flag = 0;
int                   opt_probeset_flag = 0;
int                   opt_ignore_weak = 1;
int                   opt_log2   = 1;
int                   opt_unlog2 = 0;
int                   opt_nolog2 = 0;

/* Administrative options */
const char *argp_program_version     = affy_version;
const char *argp_program_bug_address = "<Eric.Welsh@moffitt.org>";
static struct argp_option options[]  = 
{
  { "dir",          'd', "DIRECTORY", 0, "Use DIRECTORY as working directory"  },
  { "cdf",          'c', "CDFDIR",    0, "Use CDFDIR as location for CDF file" },
  { "probesets",    's',           0, 0, "Use probesets instead of probes" },
  { "rmsd",         'r',           0, 0, "Use RMSD distances only [default]" },
  { "output",       'o', "OUTPUTFILE",0, "Output OUTPUTFILE instead of STDOUT" },
  { "pearson",      'p',           0, 0, "Use Pearson correlation distances only" },
  { "geomean",      'g',           0, 0, "Use geometric mean of RMSD and Pearson distances" },
  { "meancenter",   'm',           0, 0, "Mean-center log intensities prior to distance calculations" },
  { "unnormalized", 'u',           0, 0, "Do not perform any pre-normalization [default]" },
  { "spreadsheet",  't',           0, 0, "Input is a single spreadsheet of (un-logged) intensities" },
  { "salvage",       24,           0, 0, "Attempt to salvage corrupt CEL files (may still result in corrupt data!)" },
  { "log2",          25,           0, 0, "Take log2 of input data prior to dist calcs [default]" },
  { "nolog2",        26,           0, 0, "Do not transform data prior to dist calcs" },
#if 0
  { "unlog2",        27,           0, 0, "Exponentiate input data prior to dist calcs" },
#endif
  { "ignore-weak",   28,           0, 0, "Ignore weak intensities (<= 0 un-logged) [default]" },
  { "include-weak",  29,           0, 0, "Include weak intensities (<= 0 un-logged)" },
  { "ignore-chip-mismatch", 137,   0, 0, "Do not abort when multiple chips types are detected" },
  { "probeset-exclusions",'x',"EXCLUSIONSFILE",0,"Do not load probes in EXCLUSIONSFILE from spreadsheet" },
  { "probeset-spikeins",'S',"SPIKEINSSFILE",0,"Do not load probes in SPIKEINSFILE from spreadsheet" },
  { NULL }
};

static error_t parse_opt(int key, char *arg, struct argp_state *state);
static struct argp argp = { options, 
                            parse_opt, 
                            0, 
                            "findmedian - Pairwise normalization median sample finder"};


double calculate_pearson_r_float_skip_missing(float *array1, float *array2,
                                              int n)
{
  double x_avg = 0, y_avg = 0;
  double sum_xy_diff = 0, sum_x2_diff = 0, sum_y2_diff = 0;
  double x_diff, y_diff;
  double temp;
  int i;
  int count = 0;

  for (i = 0; i < n; i++)
  {
    if (array1[i] == MISSING || array2[i] == MISSING)
      continue;
  
    x_avg += array1[i];
    y_avg += array2[i];
    
    count++;
  }
  x_avg /= count;
  y_avg /= count;

  for (i = 0; i < n; i++)
  {
    if (array1[i] == MISSING || array2[i] == MISSING)
      continue;

    x_diff = array1[i] - x_avg;
    sum_x2_diff += x_diff * x_diff;

    y_diff = array2[i] - y_avg;
    sum_y2_diff += y_diff * y_diff;

    sum_xy_diff += x_diff * y_diff;
  }

  if (sum_x2_diff && sum_y2_diff)
  {
    sum_x2_diff = sqrt(sum_x2_diff);
    sum_y2_diff = sqrt(sum_y2_diff);
    temp = sum_x2_diff * sum_y2_diff;
  
    if (temp)
    {
      temp = sum_xy_diff / temp;
      
      /* round off errors can lead to slightly greater than 1 */
      if (temp >  1.0)
        return    1.0;
      if (temp < -1.0)
        return   -1.0;
      
      return temp;
    }
  }

  return 0;
}


affy_uint32 fill_all_points(char *filename, float **all_points,
                            int num_all_points,
                            char **sample_names,
                            AFFY_CDFFILE *cdf,
                            AFFY_ERROR *err)
{
  FILE *data_file;
  int max_string_len = 0;
  char *string       = NULL;
  char **fields      = NULL;
  char *sptr;
  int num_fields = 0;
  int max_field  = 0;

  int numprobes = 0;
  int numchips = 0;
  int ok_chip_flag, ok_probe_flag;
  int i;

  double log2;
  double value;

  log2 = log(2.0);
  
  data_file = fopen(filename, "rb");
  if (!data_file)
    AFFY_HANDLE_ERROR("can not open data file", AFFY_ERROR_NOTFOUND, err, 0);
  
  /* read header line */
  fgets_strip_realloc(&string, &max_string_len, data_file);
  num_fields = split_tabs(string, &fields, &max_field);
  
  /* first field is probe, the rest are samples, skip empty columns */
  for (i = 1; i < num_fields; i++)
  {
    ok_chip_flag = 0;

    /* skip empty columns */
    for (sptr = fields[i]; sptr; sptr++)
    {
      if (!isspace(*sptr))
      {
        ok_chip_flag = 1;
        break;
      }
    }
    
    /* save sample names */
    if (ok_chip_flag)
    {
      sample_names[numchips] = h_strdup(fields[i]);
      if (sample_names[numchips] == NULL)
      {
        AFFY_HANDLE_ERROR("strdup failed", AFFY_ERROR_OUTOFMEM, err, 0);
      }
      hattach(sample_names[numchips], mempool);
      
      numchips++;
    }
  }

  /* assume only a single line per probe */
  /* multiple identical probes will be treated as separate probes */
  /* skip blank probes */
  while(fgets_strip_realloc(&string, &max_string_len, data_file))
  {
    num_fields = split_tabs(string, &fields, &max_field);

    ok_probe_flag = 0;

    /* skip empty columns */
    for (sptr = fields[0]; sptr; sptr++)
    {
      if (!isspace(*sptr))
      {
        ok_probe_flag = 1;
        break;
      }
    }
    
    /* HACK -- skip probes we've decided to exclude */
    /* easier to skip them here, than skip them in the distance functions */
    /* since this is only used in findmedian, so it is OK for now... */
    if (affy_is_control_string(fields[0]))
    {
        ok_probe_flag = 0;
    }
    if (flags.use_exclusions && cdf->exclusions)
    {
        if (bsearch(&fields[0], &cdf->exclusions[0],
            cdf->numexclusions, sizeof(char *), compare_string))
        {
            ok_probe_flag = 0;
        }
    }
    /* we want to exclude spikeins too */
    if (flags.use_spikeins && cdf->spikeins)
    {
        if (bsearch(&fields[0], &cdf->spikeins[0],
            cdf->numspikeins, sizeof(char *), compare_string))
        {
            ok_probe_flag = 0;
        }
    }
    
    if (ok_probe_flag)
    {
      numchips = 0;

      for (i = 1; i < num_fields; i++)
      {
        ok_chip_flag = 0;

        /* skip empty columns */
        for (sptr = fields[i]; sptr; sptr++)
        {
          if (!isspace(*sptr))
          {
            ok_chip_flag = 1;
            break;
          }
        }

        /* store intensity data */
        if (ok_chip_flag)
        {
          value = atof(fields[i]);
          
          if (opt_ignore_weak && value <= 0)
          {
            value = MISSING;
          }
          else
          {
            if (opt_log2)
            {
              if (value < MIN_VALUE)
                value = MIN_VALUE;

              value = log(value) / log2;
            }
            else if (opt_unlog2)
              value = pow(2, value);
          }

          all_points[numchips++][numprobes] = value;
        }
      }

      numprobes++;
    }
  }
  
  fclose(data_file);

  if (string)
    free(string);
  if (fields)
    free(fields);
  
  return numprobes;
}


double normalize_mean_center(float *array1, int n)
{
    double avg = 0;
    int i, count;
    
    count = 0;
    for (i = 0; i < n; i++)
    {
       if (opt_ignore_weak && array1[i] == MISSING)
           continue;

       avg += array1[i];
       count++;
    }
    if (count)
       avg /= count;
    
    /* normalize */
    for (i = 0; i < n; i++)
    {
      if (opt_ignore_weak && array1[i] == MISSING)
          continue;

      array1[i] = array1[i] - avg;
    }
    
    return 0;
}


double normalize_trim_mean(float *array1, int n)
{
    double *array2 = NULL;
    double avg = 0;
    int i, n2;
    int start, end, width;
    
    array2 = (double *) malloc(n * sizeof(double));

    n2 = 0;
    for (i = 0; i < n; i++)
    {
       if (opt_ignore_weak && array1[i] == MISSING)
          continue;
       
       array2[n2++] = array1[i];
    }
    
    if (n2)
       qsort(array2, n2, sizeof(double), dcompare);
    
    width = 0.05 * n2;
    start = width;
    end = n2 - width - 1;
    
    if (end <= start)
    {
       start = 0;
       end = n2 - 1;
    }

    width = end - start + 1;
    
    
    for (i = start; i <= end; i++)
    {
       avg += array2[i];
    }
    if (n2)
       avg /= width;
    
    /* normalize */
    for (i = 0; i < n; i++)
    {
      if (opt_ignore_weak && array1[i] == MISSING)
          continue;

      array1[i] = array1[i] - avg;
    }
    
    return 0;
}


int pairgen_find_median_chip_distance(float **all_points,
                     affy_uint32 num_probes, unsigned int max_chips,
                     int method_flag,
                     char **sample_names, double *rmsd,
                     double *avg_rmsd, AFFY_ERROR *err)
{
  FILE        *outfile = NULL;
  char        *buffer_out = NULL;
  affy_int32  i, j, point_idx;
  double      *distances = NULL;
  double     **distance_rows = NULL;
  double      *distances_pearson = NULL;
  double     **distance_rows_pearson = NULL;
  double     **count_rows = NULL;
  double      *counts = NULL;
  double      *means_sample = NULL;
  int         *counts_sample_non_weak = NULL;
  double       diff, sum, average;
  double       best_score = 9E99;
  float       *fptr1, *fptr2;
  int          best_chip = -1;
  int          max_chips_squared = max_chips * max_chips;
  int          count;
  int          max_count = 0;

  /* open file for writing */    
  if (outfile_name)
  {
    /* open as text, so it will translate EOL automatically */
    outfile = fopen(outfile_name, "w");

    if (!outfile)
    {
      fprintf(stderr, "ERROR -- can't open output file %s\n",
              outfile_name);

      return 1;
    }
  }
  else
  {
    outfile = stdout;
  }

  /* we shouldn't need buffering for efficiency, but Emscripten is breaking
   * on file output *only* in this program, and *only* when compiled with
   * Emscripten, and this is one of the last differences between this and
   * other working programs I have left to try to magically make Emscripten
   * happy.
   */
  buffer_out = (char *) malloc(1048576 * sizeof(char));
  setvbuf(outfile, buffer_out, _IOFBF, 1048576);

  if (opt_ignore_weak)
    fprintf(stderr, "Ignoring points with weak values\n");

  means_sample = h_subcalloc(mempool, max_chips, sizeof(double));
  if (means_sample == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed",
                           AFFY_ERROR_OUTOFMEM,
                           err,
                           cleanup);

  if (opt_ignore_weak)
  {
    counts_sample_non_weak = h_subcalloc(mempool, max_chips, sizeof(int));
    if (counts_sample_non_weak == NULL)
      AFFY_HANDLE_ERROR_GOTO("calloc failed",
                              AFFY_ERROR_OUTOFMEM,
                              err,
                              cleanup);

    counts = h_subcalloc(mempool, max_chips_squared, sizeof(double));
    if (counts == NULL)
      AFFY_HANDLE_ERROR_GOTO("calloc failed",
                              AFFY_ERROR_OUTOFMEM,
                              err,
                              cleanup);

    count_rows = h_subcalloc(mempool, max_chips, sizeof(double *));
    if (count_rows == NULL)
      AFFY_HANDLE_ERROR_GOTO("calloc failed",
                             AFFY_ERROR_OUTOFMEM,
                             err,
                             cleanup);

    for (i = 0; i < max_chips; i++)
      count_rows[i]  = counts + i * max_chips;
    memset(counts, 0, max_chips_squared * sizeof(double));


    fprintf(stderr, "Pre-scan pair-wise missing data ");

    /* count non-weak points for each reference chip */
    for (i = 0; i < max_chips; i++)
    {
      fptr1 = all_points[i];

      fprintf(stderr, ".");

      count = 0;
      for (point_idx = 0; point_idx < num_probes; point_idx++)
      {
        /* skip points with missing values and/or super-weak intensities */
        if (opt_ignore_weak)
        {
          if (fptr1[point_idx] != MISSING)
          {
            means_sample[i] += fptr1[point_idx];
            count++;
          }
        }
        else
        {
          means_sample[i] += fptr1[point_idx];
          count++;
        }
      }
      counts_sample_non_weak[i] = count;
      
      if (count)
          means_sample[i] /= count;

      for (j = i+1; j < max_chips; j++)
      {
        fptr2 = all_points[j];

        count = 0;
        for (point_idx = 0; point_idx < num_probes; point_idx++)
        {
          if (opt_ignore_weak)
          {
            /* skip points with missing values and/or super-weak intensities */
            if (fptr1[point_idx] == MISSING || fptr2[point_idx] == MISSING)
              continue;
          }

          count++;
        }

        /* avoid potential divide by zero problems later */
        if (count < 1)
          count = 1;

        count_rows[i][j] = count;
        count_rows[j][i] = count;

        if (count > max_count)
          max_count = count;
      }
    }
    fprintf(stderr, "\n");
  }
  else
  {
    /* calculate the log2 mean of all points, treating missing as zero */
    for (i = 0; i < max_chips; i++)
    {
      fptr1 = all_points[i];

      count = 0;
      for (point_idx = 0; point_idx < num_probes; point_idx++)
      {
          means_sample[i] += fptr1[point_idx];
          count++;
      }

      if (count)
          means_sample[i] /= count;
    }
  }

  /* RMSD distances */
  if (method_flag == 'r' || method_flag == 'g')
  {
    distances = h_subcalloc(mempool, max_chips_squared, sizeof(double));
    if (distances == NULL)
      AFFY_HANDLE_ERROR_GOTO("calloc failed",
                             AFFY_ERROR_OUTOFMEM,
                             err,
                             cleanup);

    distance_rows = h_subcalloc(mempool, max_chips, sizeof(double *));
    if (distance_rows == NULL)
      AFFY_HANDLE_ERROR_GOTO("calloc failed",
                             AFFY_ERROR_OUTOFMEM,
                             err,
                             cleanup);

    /* initialize row pointers */
    for (i = 0; i < max_chips; i++)
      distance_rows[i] = distances + i * max_chips;
    memset(distances, 0, max_chips_squared * sizeof(double));


    fprintf(stderr, "Finding median sample in RMSD space: ");

    for (i = 0; i < max_chips; i++)
    {
      fptr1 = all_points[i];

      fprintf(stderr, ".");
      
      for (j = i+1; j < max_chips; j++)
      {
        fptr2 = all_points[j];

        sum = 0.0;

        count = 0;
        for (point_idx = 0; point_idx < num_probes; point_idx++)
        {
          if (opt_ignore_weak)
          {
            /* skip points with missing values and/or super-weak intensities */
            if (fptr1[point_idx] == MISSING || fptr2[point_idx] == MISSING)
              continue;
          }
        
#if 0
          /* divide difference by average */
          diff = 2 * (fptr1[point_idx] - fptr2[point_idx]) /
                     (fptr1[point_idx] + fptr2[point_idx]);
#else
          /* we're already in log-ratio space, no need to normalize by avg */
          diff = fptr1[point_idx] - fptr2[point_idx];
#endif
          sum += diff * diff;
          count++;
        }

        if (count)
          sum = sqrt(sum / count);
    
        distance_rows[i][j] = sum;
        distance_rows[j][i] = sum;
      }
    }
    fprintf(stderr, "\n");
  }

  /* Pearson distances */
  if (method_flag == 'p' || method_flag == 'g')
  {
    distances_pearson = h_subcalloc(mempool, max_chips_squared, sizeof(double));
    if (distances_pearson == NULL)
      AFFY_HANDLE_ERROR_GOTO("calloc failed",
                             AFFY_ERROR_OUTOFMEM,
                             err,
                             cleanup);

    distance_rows_pearson = h_subcalloc(mempool, max_chips, sizeof(double *));
    if (distance_rows_pearson == NULL)
      AFFY_HANDLE_ERROR_GOTO("calloc failed",
                             AFFY_ERROR_OUTOFMEM,
                             err,
                             cleanup);

    /* initialize distances_pearson[] row pointers */
    for (i = 0; i < max_chips; i++)
      distance_rows_pearson[i] = distances_pearson + i * max_chips;

    memset(distances_pearson, 0, max_chips_squared * sizeof(double));

    fprintf(stderr, "Finding median sample in Pearson space: ");

    for (i = 0; i < max_chips; i++)
    {
      fprintf(stderr, ".");
    
      fptr1 = all_points[i];
  
      for (j = i+1; j < max_chips; j++)
      {
        fptr2 = all_points[j];

        if (opt_ignore_weak)
        {
          sum = calculate_pearson_r_float_skip_missing(fptr1, fptr2,
                                                       num_probes);
        }
        else
        {
          sum = calculate_pearson_r_float(fptr1, fptr2, num_probes);
        }
        
        /* transform r into a metric distance [sqrt(0.5 * (1-r))] */
        /* see Stijn van Dongen and Anton J. Enright 2012
         * "Metric distance derived from cosine similarity and Pearson
         *  and Spearman correlations"
         */
        /* deal with potential nan-inducing roundoff error */
        sum = 1.0 - sum;
        if (sum < 0.0)
          sum = 0.0;
        sum = sqrt(0.5 * sum);

        distance_rows_pearson[i][j] = sum;
        distance_rows_pearson[j][i] = sum;
      }
    }
    fprintf(stderr, "\n");
  }

  /* point distances at pearson distance matrix */
  if (method_flag == 'p')
  {
    distance_rows = distance_rows_pearson;
  }
  /* take geometric means of RMSD and Pearson distances */
  else if (method_flag == 'g')
  {
    fprintf(stderr, "Finding median using geometric mean of RMSD and Pearson distances\n");

    for (i = 0; i < max_chips; i++)
    {
      for (j = i + 1; j < max_chips; j++)
      {
        distance_rows[i][j] = distance_rows[j][i] =
                              sqrt(distance_rows[i][j] *
                                   distance_rows_pearson[i][j]);
      }
    }
  }
  
  
  /* normalize counts */
  if (opt_ignore_weak)
  {
    for (i = 0; i < max_chips; i++)
    {
      for (j = 0; j < max_chips; j++)
      {
        if (i == j)
          continue;

        /* Conditional probability of #overlap given #query.
         *
         * This is assymetric, but that's what we want in this case.
         * We are searching for candidate reference samples.
         * Say that we have a reference with 100% present, and the query
         * is only 25% present.  We don't want to use some similarity metric
         * that combines those two #present together, since we want to know
         * if the query is maximally covered or not.  #overlap / #query will
         * be 1.0, since it is fully covered.  The other way around,
         * #overlap / #target, would be 0.25.  This is the scoring behavior
         * we want, since the first case indicates a good reference sample
         * and the 2nd case indicates a poorer reference sample.
         */
        distance_rows[i][j] /= count_rows[i][j] /
                               counts_sample_non_weak[j];
      }
    }
  }


  /* print header line */
  fprintf(outfile, "%s",   "Score");
  fprintf(outfile, "\t%s", "SampleIndex");
  fprintf(outfile, "\t%s", "MeanDistance");
  fprintf(outfile, "\t%s", "SampleName");
  fprintf(outfile, "\t%s", "MeanLog2Abundance");
  fprintf(outfile, "\n");

  /* calculate means of distances to other chips */
  average = 0;
  for (i = 0; i < max_chips; i++)
  {
    sum = 0;
    for (j = 0; j < max_chips; j++)
    {
      if (i == j)
        continue;

      sum += distance_rows[i][j];
    }
    
    if (max_chips > 1)
      sum = sum / (max_chips - 1);
    
    average += sum;
    
    if (sum < best_score)
    {
      best_score = sum;
      best_chip = i;
    }

    fprintf(outfile, "Score\t%d\t%f\t%s\t%f\n",
            i,
            sum,
            sample_names[i],
            means_sample[i]);
  }
  
  average /= max_chips;
  
  *rmsd = best_score;
  *avg_rmsd = average;


  fprintf(outfile, "Average RMSD:\t%f\n", *avg_rmsd);
  fprintf(outfile, "Median CEL:\t%d\t%f\t%s\t%f\n",
          best_chip, *rmsd, sample_names[best_chip], means_sample[best_chip]);
  fprintf(outfile, "%s\n", sample_names[best_chip]);

  fclose(outfile);
  free(buffer_out);

  return best_chip;


cleanup:
  if (outfile) fclose(outfile);
  if (buffer_out)  free(buffer_out);
  *rmsd = 999999;
  h_free(mempool);

  return -1;
}


int main(int argc, char **argv)
{
  float           **all_points = NULL;     /* save memory with doubles */
  char            **sample_names = NULL;
  AFFY_CHIPSET     *cs = NULL;
  AFFY_CHIPSET     *temp = NULL;
  AFFY_CDFFILE     *cdf;
  AFFY_CELFILE     *cel;
  int               status = EXIT_FAILURE, i;
  unsigned int      max_chips;
  affy_uint32       num_all_points;
  affy_uint32       num_points = 0;
  affy_uint32       num_probesets;
  affy_int32        point_idx;
  char             *chip_type, **p;
  LIBUTILS_PB_STATE pbs;
  AFFY_ERROR       *err = NULL;
  
  double            log2;
  double            rmsd, avg_rmsd;
  double            value;
  int temp_idx;
  
  log2 = log(2.0);

  pb_init(&pbs);

  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
  {
    fprintf(stderr, "malloc failed: out of memory\n");
    exit(EXIT_FAILURE);
  }

  err = affy_get_default_error();

  affy_rma_set_defaults(&flags);
  affy_mas5_set_defaults(&flags);

  flags.bioconductor_compatability = false;
  flags.use_background_correction  = true;
  flags.bg_rma                     = true;
  flags.bg_mas5                    = false;
  flags.use_tukey_biweight         = true;
  flags.use_median_polish          = false;

  argp_parse(&argp, argc, argv, 0, 0, 0);

  /* If files is NULL, open all CEL files in the current working directory */
  if (filelist == NULL && SEARCH_WORKING_DIR)
    filelist = affy_list_files(directory, ".cel", err);

  /* Give up if we have no files to operate on */
  if ((filelist == NULL) || (filelist[0] == NULL))
  {
    if (SEARCH_WORKING_DIR)
    {
      fprintf(stderr, 
              "no CEL files specified or found in current working directory, exiting\n");
    }
    else
    {
      fprintf(stderr, 
              "no input files specified, exiting\n");
    }

    goto cleanup;
  }

  /* Input is a spreadsheet of intensities */
  if (opt_spreadsheet_flag)
  {
    get_generic_spreadsheet_bounds(filelist[0], &num_all_points,
                                   &max_chips, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    info("NumSamples:\t%d\tNumProbes:\t%d", max_chips, num_all_points);

    cdf = NULL;
    if (flags.use_exclusions || flags.use_spikeins)
    {
        fprintf(stderr, "Exclusion probeset filename:         %s\n",
                flags.exclusions_filename);

        /* create a chipset and cdf to store the exclusions */
        cs = create_blank_generic_chipset(max_chips, num_all_points, err);
        AFFY_CHECK_ERROR_GOTO(err, cleanup);

        /* load in probesets to exclude from IRON training */
        if (flags.use_exclusions)
            affy_load_exclusions_file(flags.exclusions_filename,
                                      cs->cdf, mempool, err);
        if (flags.use_spikeins)
            affy_load_spikeins_file(flags.spikeins_filename,
                                    cs->cdf, mempool, err);
        cdf = cs->cdf;
    }

    /* allocate chip names */
    sample_names = h_subcalloc(mempool, max_chips, sizeof(char **));
    if (sample_names == NULL)
      AFFY_HANDLE_ERROR_GOTO("calloc failed",
                             AFFY_ERROR_OUTOFMEM,
                             err,
                             cleanup);

    /* Set up the 2D array of chips/probes */
    all_points = h_subcalloc(mempool, max_chips, sizeof(float *));
    if (all_points == NULL)
    {
      fprintf(stderr, "calloc failed: out of memory\n");
      goto cleanup;
    }

    all_points[0] = h_subcalloc(mempool, max_chips * num_all_points,
                                sizeof(float));
    if (all_points[0] == NULL)
    {
      fprintf(stderr, "calloc failed: out of memory\n");
      goto cleanup;
    }

    for (point_idx = 1; point_idx < max_chips; point_idx++)
      all_points[point_idx] = all_points[point_idx - 1] + num_all_points;
    
    num_points = fill_all_points(filelist[0], all_points, num_all_points,
                                 sample_names, cdf, err);
    
    info("Finished reading %u samples, %d variables", max_chips, num_points);
  }

  /* Input is separate CEL files, probe-level */
  else if (opt_probeset_flag == 0)
  {
    /* Count up chips */
    for (p = filelist, max_chips = 0; *p != NULL; p++)
      max_chips++;

    /* Get the array type. */
    chip_type = affy_get_cdf_name_from_cel(filelist[0], err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    hattach(chip_type, mempool);

    /* Create temp chipset */
    cs = affy_create_chipset(1, chip_type, cdf_directory, &flags, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    temp = affy_clone_chipset(cs, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    cs = affy_resize_chipset(cs, max_chips, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    cdf = cs->cdf;

    hattach(cs, mempool);

    /* load in probesets to exclude */
    if (flags.use_exclusions)
        affy_load_exclusions_file(flags.exclusions_filename,
                                cdf, mempool, err);
    if (flags.use_spikeins)
        affy_load_spikeins_file(flags.spikeins_filename,
                                cdf, mempool, err);

    /* count number of PM's we will keep */
    num_all_points = 0;
    memset(cdf->seen_xy[0], 0, cdf->numrows*cdf->numcols*sizeof(affy_uint8));
    for (temp_idx = 0; temp_idx < cdf->numprobes; temp_idx++)
    {
        int x, y;
      
        x = cdf->probe[temp_idx]->pm.x;
        y = cdf->probe[temp_idx]->pm.y;
        
        if (cdf->seen_xy[x][y])
          continue;
        cdf->seen_xy[x][y] = 1;

        /* skip non-probeset probes */
        if ((cdf->cell_type[x][y] == AFFY_UNDEFINED_LOCATION) ||
            (cdf->cell_type[x][y] == AFFY_QC_LOCATION))
        {
          continue;
        }
        
        if (affy_is_control_string(cdf->probe[temp_idx]->ps->name))
          continue;

        /* skip excluded probesets */
        if (flags.use_exclusions && cdf->exclusions)
        {
          if (bsearch(&cdf->probe[temp_idx]->ps->name, &cdf->exclusions[0],
              cdf->numexclusions, sizeof(char *), compare_string))
          {
            continue;
          }
        }
        if (flags.use_spikeins && cdf->spikeins)
        {
          if (bsearch(&cdf->probe[temp_idx]->ps->name, &cdf->spikeins[0],
              cdf->numspikeins, sizeof(char *), compare_string))
          {
            continue;
          }
        }
        
        num_all_points++;
    }

    /* Set up the 2D array of chips/probes */
    all_points = h_subcalloc(mempool, num_all_points, sizeof(float *));
    if (all_points == NULL)
    {
      fprintf(stderr, "calloc failed: out of memory\n");
      goto cleanup;
    }

    all_points[0] = h_subcalloc(mempool, max_chips * num_all_points,
                                sizeof(float));
    if (all_points[0] == NULL)
    {
      fprintf(stderr, "calloc failed: out of memory\n");
      goto cleanup;
    }

    for (point_idx = 1; point_idx < max_chips; point_idx++)
      all_points[point_idx] = all_points[point_idx - 1] + num_all_points;

    for (i = 0; i < max_chips; i++)
    {
      /* Load each chip, abort on failure */
      affy_load_chipset_single(cs, filelist[i],
                               flags.ignore_chip_mismatch, err);
      if (err->type != AFFY_ERROR_NONE)
        AFFY_CHECK_ERROR_GOTO(err, cleanup);

      /* Temp chipset now contains the most recently loaded chip */
      temp->chip[0] = cs->chip[(cs->num_chips) - 1];
      temp->num_chips = 1;
      
      cel = temp->chip[0]->cel;
      
      if (cel->corrupt_flag && flags.salvage_corrupt == false)
        AFFY_HANDLE_ERROR_GOTO("corrupt CEL file",
                               AFFY_ERROR_BADFORMAT,
                               err,
                               cleanup);

      point_idx = 0;
      memset(cdf->seen_xy[0], 0, cdf->numrows*cdf->numcols*sizeof(affy_uint8));
      for (temp_idx = 0; temp_idx < cdf->numprobes; temp_idx++)
      {
        int x, y;
      
        x = cdf->probe[temp_idx]->pm.x;
        y = cdf->probe[temp_idx]->pm.y;

        if (cdf->seen_xy[x][y])
          continue;
        cdf->seen_xy[x][y] = 1;

        /* skip non-probeset probes */
        if ((cdf->cell_type[x][y] == AFFY_UNDEFINED_LOCATION) ||
            (cdf->cell_type[x][y] == AFFY_QC_LOCATION))
        {
          continue;
        }

        if (affy_is_control_string(cdf->probe[temp_idx]->ps->name))
          continue;

        /* skip excluded probesets */
        if (flags.use_exclusions && cdf->exclusions)
        {
          if (bsearch(&cdf->probe[temp_idx]->ps->name, &cdf->exclusions[0],
              cdf->numexclusions, sizeof(char *), compare_string))
          {
            continue;
          }
        }
        if (flags.use_spikeins && cdf->spikeins)
        {
          if (bsearch(&cdf->probe[temp_idx]->ps->name, &cdf->spikeins[0],
              cdf->numspikeins, sizeof(char *), compare_string))
          {
            continue;
          }
        }
        
        value = cel->data[x][y].value;

        if (opt_ignore_weak && value <= 0)
        {
          value = MISSING;
        }
        else
        {
          if (opt_log2)
          {
            if (value < MIN_VALUE)
              value = MIN_VALUE;

            value = log(value) / log2;
          }
          else if (opt_unlog2)
            value = pow(2, value);
        }

        all_points[i][point_idx++] = value;
      }

      num_points = point_idx;

      affy_mostly_free_chip(temp->chip[0]);
    }

    /* allocate chip names */
    sample_names = h_subcalloc(mempool, max_chips, sizeof(char **));
    if (sample_names == NULL)
      AFFY_HANDLE_ERROR_GOTO("calloc failed",
                             AFFY_ERROR_OUTOFMEM,
                             err,
                             cleanup);

    /* fill sample_names */
    for (i = 0; i < max_chips; i++)
      sample_names[i] = filelist[i];

    info("Finished reading %u CEL files", max_chips);
  }

  /* Input is separate CEL files, probeset-level */
  else if (opt_probeset_flag == 1)
  {
    /* Count up chips */
    for (p = filelist, max_chips = 0; *p != NULL; p++)
      max_chips++;

    /* Get the array type. */
    chip_type = affy_get_cdf_name_from_cel(filelist[0], err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    hattach(chip_type, mempool);

    /* Create temp chipset */
    cs = affy_create_chipset(1, chip_type, cdf_directory, &flags, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    temp = affy_clone_chipset(cs, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    cs = affy_resize_chipset(cs, max_chips, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    cdf = cs->cdf;

    hattach(cs, mempool);

    num_probesets = cs->cdf->numprobesets;

    /* load in probesets to exclude */
    if (flags.use_exclusions)
        affy_load_exclusions_file(flags.exclusions_filename,
                                  cs->cdf, mempool, err);
    if (flags.use_spikeins)
        affy_load_spikeins_file(flags.spikeins_filename,
                                cs->cdf, mempool, err);

    /* Set up the 2D array of chips/probesets */
    all_points = h_subcalloc(mempool, num_probesets, sizeof(float *));
    if (all_points == NULL)
    {
      fprintf(stderr, "calloc failed: out of memory\n");
      goto cleanup;
    }

    all_points[0] = h_subcalloc(mempool, max_chips * num_probesets,
                                sizeof(float));
    if (all_points[0] == NULL)
    {
      fprintf(stderr, "calloc failed: out of memory\n");
      goto cleanup;
    }

    for (point_idx = 1; point_idx < max_chips; point_idx++)
      all_points[point_idx] = all_points[point_idx - 1] + num_probesets;

    for (i = 0; i < max_chips; i++)
    {
      /* Load each chip, abort on failure */
      affy_load_chipset_single(cs, filelist[i],
                               flags.ignore_chip_mismatch, err);
      if (err->type != AFFY_ERROR_NONE)
        AFFY_CHECK_ERROR_GOTO(err, cleanup);

      /* Temp chipset now contains the most recently loaded chip */
      temp->chip[0] = cs->chip[(cs->num_chips) - 1];
      temp->num_chips = 1;
      
      /* Tukey's Biweight */
      affy_mas5_signal(temp, &flags, err);
      AFFY_CHECK_ERROR_GOTO(err, cleanup);

      num_all_points = 0;
      for (point_idx = 0; point_idx < num_probesets; point_idx++)
      {
        if (affy_is_control_string(cdf->probeset[point_idx].name))
          continue;

        /* skip excluded probesets */
        if (flags.use_exclusions && cdf->exclusions)
        {
          if (bsearch(&cdf->probeset[point_idx].name, &cdf->exclusions[0],
              cdf->numexclusions, sizeof(char *), compare_string))
          {
            continue;
          }
        }
        if (flags.use_spikeins && cdf->spikeins)
        {
          if (bsearch(&cdf->probeset[point_idx].name, &cdf->spikeins[0],
              cdf->numspikeins, sizeof(char *), compare_string))
          {
            continue;
          }
        }

        value = temp->chip[0]->probe_set[point_idx];
        
        if (opt_ignore_weak && value <= 0)
        {
          value = MISSING;
        }
        else
        {
          if (opt_log2)
          {
            if (value < MIN_VALUE)
              value = MIN_VALUE;

            value = log(value) / log2;
          }
          else if (opt_unlog2)
            value = pow(2, value);
        }

        all_points[i][num_all_points++] = value;
      }

      num_points = num_all_points;

      affy_mostly_free_chip(temp->chip[0]);
    }

    /* allocate chip names */
    sample_names = h_subcalloc(mempool, max_chips, sizeof(char **));
    if (sample_names == NULL)
      AFFY_HANDLE_ERROR_GOTO("calloc failed",
                             AFFY_ERROR_OUTOFMEM,
                             err,
                             cleanup);

    /* fill sample_names */
    for (i = 0; i < max_chips; i++)
      sample_names[i] = filelist[i];

    info("Finished reading %u CEL files", max_chips);
  }
  

  if (opt_mean_center_flag)
  {
/*    fprintf(stderr, "Normalizing probes to mean-centered unit variance\n");
*/
    fprintf(stderr, "Mean-centering vectors\n");
    
    for (i = 0; i < max_chips; i++)
    {
      normalize_mean_center(all_points[i], num_points);
/*      normalize_trim_mean(all_points[i], num_points); */
    }
  }


  if (opt_distance_rmsd_flag && opt_distance_pearson_flag)
  {
    temp_idx = pairgen_find_median_chip_distance(all_points,
                                  num_points, max_chips,
                                  'g',
                                  sample_names,
                                  &rmsd, &avg_rmsd,
                                  err);
  }
  else if (opt_distance_rmsd_flag)
  {
    temp_idx = pairgen_find_median_chip_distance(all_points,
                                  num_points, max_chips,
                                  'r',
                                  sample_names,
                                  &rmsd, &avg_rmsd,
                                  err);
  }
  else if (opt_distance_pearson_flag)
  {
    temp_idx = pairgen_find_median_chip_distance(all_points,
                                  num_points, max_chips,
                                  'p',
                                  sample_names,
                                  &rmsd, &avg_rmsd,
                                  err);
  }

  if (opt_spreadsheet_flag == 0)
    print_corrupt_chips_to_stderr(cs);

  /* everything OK so far, exit 0 */
  status = 0;

cleanup:  
  if (cs)
    affy_free_chipset(cs);
  if (temp)
    affy_free_chipset(temp);
  h_free(mempool);
  pb_cleanup(&pbs);

  if (err)
    free(err);

  exit(status);
}

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
  switch (key)
  {
    case 24:
      flags.salvage_corrupt = true;
      break;
    case 25:
      opt_log2   = 1;
      opt_nolog2 = 0;
      opt_unlog2 = 0;
      break;
    case 26:
      opt_log2   = 0;
      opt_nolog2 = 1;
      opt_unlog2 = 0;
      break;
    case 27:
      opt_log2   = 0;
      opt_nolog2 = 0;
      opt_unlog2 = 1;
      break;
    case 28:
      opt_ignore_weak = 1;
      break;
    case 29:
      opt_ignore_weak = 0;
      break;
    case 137:
      flags.ignore_chip_mismatch = true;
      break;

    case 'd':
      directory = h_strdup(arg);
      hattach(directory, mempool);
      break;
    case 'c':
      cdf_directory = h_strdup(arg);
      hattach(cdf_directory, mempool);
      break;
    case 'o':
      outfile_name = h_strdup(arg);
      hattach(outfile_name, mempool);
      break;
    case 'p':
      opt_distance_rmsd_flag = 0;
      opt_distance_pearson_flag = 1;
      break;
    case 'r':
      opt_distance_rmsd_flag = 1;
      opt_distance_pearson_flag = 0;
      break;
    case 'g':
      opt_distance_rmsd_flag = 1;
      opt_distance_pearson_flag = 1;
      break;
    case 'm':
      opt_mean_center_flag = 1;
      break;
    case 's':
      opt_probeset_flag = 1;
      break;
    case 't':
      opt_spreadsheet_flag = 1;
      break;
    case 'u':
      opt_mean_center_flag = 0;
      break;
    case 'S':
      flags.use_spikeins = true;
      if (arg != NULL)
      {
        flags.spikeins_filename = h_strdup(arg);
        hattach(flags.spikeins_filename, mempool);
      }
      break;
    case 'x':
      flags.use_exclusions = true;
      if (arg != NULL)
      {
        flags.exclusions_filename = h_strdup(arg);
        hattach(flags.exclusions_filename, mempool);
      }
      break;
    case ARGP_KEY_ARG: /* Do not try to get filenames one at a time, use all */
      return ARGP_ERR_UNKNOWN;
    case ARGP_KEY_ARGS:
      filelist = state->argv + state->next;
      break;
    default:
      return ARGP_ERR_UNKNOWN;
  }

  return (0);
}
