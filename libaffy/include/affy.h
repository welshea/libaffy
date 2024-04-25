/**************************************************************************
 *
 * Filename:  affy.h
 * 
 * Purpose:   Top-level header file for libaffy.
 *
 * Creation:  04/04/2005
 *
 * Author:    Steven Eschrich
 *
 * Copyright: Copyright (C) 2007, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 04/04/05: Import from old libaffy (AMH)
 * 04/08/05: DAT pixeldata now a 2D array of unsigned ints (AMH) 
 * 04/19/05: Add AFFY_CELL and AFFY_PIXREGION types (AMH)
 * 09/27/07: Add data structures for the new Calvin file format (AMH)
 * 10/05/07: Update binary I/O functions (AMH)
 * 03/05/08: New error handling system (AMH)
 * 09/16/10: Add support for P/A calls (EAW)
 * 09/20/10: Pooled memory allocator (AMH)
 * 11/16/10: Chipset affinity array and t_value for Median Polish reuse (EAW)
 * 11/19/10: Added affy_mean() function prototype, pass median flags (EAW)
 * 03/11/11: Added affy_mean_geometric_floor_1() function prototype (EAW)
 * 04/06/11: Added string_io.h (EAW)
 * 04/07/11: Added generic chipset/cdf support (EAW)
 * 04/08/11: Extended AFFY_POINT y to int32 to support generic chips (EAW)
 * 01/10/13: renamed AFFY_PAIRWISE_LINEAR to AFFY_PAIRWISE_GLOBAL_SCALING (EAW)
 * 01/10/13: added new AFFY_PAIRWISE_LINEAR_SCALING (EAW)
 * 05/09/13: added AFFY_POINT16, since calvin I/O code requires both X and Y
 *           to be 16-bit when reading in data from files, and the earlier
 *           change to AFFY_POINT to make Y 32-bit messed things up (EAW)
 * 05/10/13: added affy_mostly_free_cel_file() , affy_mostly_free_chip() (EAW)
 * 10/22/13: various: numrows, numcols 16->32 bit (EAW)
 * 03/06/14: added no_mm_flag to AFFY_CDFFILE (EAW)
 * 03/10/14: #ifdef out CEL qc fields to save memory (EAW)
 * 06/01/18: added support for probeset exclusions during IRON training (EAW)
 * 09/11/18: pass mempool to probeset exclusions loader (EAW)
 * 09/14/18: added spikein related exclusions (EAW)
 * 03/13/19: added affy_is_control_probeset(), affy_is_control_string() (EAW)
 * 03/14/19: changed some int mempool to void mempool (EAW)
 * 05/22/19: added --ignore-chip-mismatch support (EAW)
 * 08/13/19: removed "unsigned" from affy_is_control_string() (EAW)
 * 10/15/19: added affy_floor_probeset_non_zero_to_one() (EAW)
 * 08/12/20: pass flags to more functions (EAW)
 * 01/10/24: pass flags to affy_mean_normalization() (EAW)
 * 04/25/24: added affy_median_normalization() (EAW)
 *
 **************************************************************************/

#ifndef _AFFY_H_
#define _AFFY_H_

/* All system-level includes needed by libaffy are included here. */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <limits.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <stdarg.h>

#ifdef __cplusplus
extern "C"
{
#endif

#include "utils.h"
#include "halloc.h"
#include "affy_basetypes.h"
#include "affy_version.h"
#include "affy_apps_common.h"
#include "string_io.h"

  /* Types of cell locations on an Affy chip. */
#define AFFY_UNDEFINED_LOCATION 0
#define AFFY_QC_LOCATION        1
#define AFFY_NORMAL_LOCATION    2

  /* Binary-format magic numbers. */
#define AFFY_DAT_FILEMAGIC        0xFC
#define AFFY_CDF_BINARYFILE_MAGIC 67
#define AFFY_CEL_BINARYFILE_MAGIC 64
#define AFFY_CALVIN_FILEMAGIC     59

  /* Everybody's favorite number */
#define AFFY_PI 3.14159265358979323846

  /* Options for affy_write_probe_values() */
#define AFFY_USE_PM 1
 
  /* Options for affy_write_expressions() */
#define AFFY_WRITE_EXPR_DEFAULT 0
#define AFFY_WRITE_EXPR_PA      1
#define AFFY_WRITE_EXPR_UNLOG   2
#define AFFY_WRITE_EXPR_LOG     4

  /* Options for affy_pairwise_normalization() */
#define AFFY_PAIRWISE_DEFAULT           0
#define AFFY_PAIRWISE_PM_ONLY           1
#define AFFY_PAIRWISE_GLOBAL_SCALING    2
#define AFFY_PAIRWISE_LINEAR_SCALING    3

  /**************************************************************************/

  /* First, some primitive types used to compose the higher-level, more
     abstract things. */

  /* 
   * Generic error codes used within the library.
   */
  typedef enum
    {
      AFFY_ERROR_NONE         = 0,
      AFFY_ERROR_NOTFOUND     = 1,
      AFFY_ERROR_SYSPERM      = 2,
      AFFY_ERROR_NOTREADY     = 3,
      AFFY_ERROR_LIMITREACHED = 4,
      AFFY_ERROR_IO           = 5,
      AFFY_ERROR_WRONGTYPE    = 6,
      AFFY_ERROR_OUTOFMEM     = 7,
      AFFY_ERROR_BADPARAM     = 8,
      AFFY_ERROR_BADFORMAT    = 9,
      AFFY_ERROR_NOTSUPP      = 10,
      AFFY_ERROR_UNKNOWN      = 99,
      AFFY_ERROR_USER         = 100   /* extended error codes start here */
    } AFFY_ERROR_TYPE;

  /*
   * Error/exception handling block.
   */
  typedef struct affy_error_s
  {
    AFFY_ERROR_TYPE type;           /* Type of error             */
    time_t          timestamp;      /* Time of occurrence        */
    const char     *descr;          /* Extended description      */
    const char     *module;         /* Module/filename of origin */
    int             location;       /* Location/line # of error  */
    void          (*handler)(struct affy_error_s *err);
  } AFFY_ERROR;

  /*
   * Calvin I/O context structure analogous to AFFY_TEXTIO.
   *
   * Users of the library should not rely on the contents of this
   * structure, it is for private use of the Calvin I/O layer.
   */
  typedef struct affy_calvinio_s
  {
    FILE       *fp;
  
    affy_uint8  file_version;
    affy_uint32 first_datagroup;
    affy_uint32 num_datagroups;
  } AFFY_CALVINIO;

  /*
   * A simple structure for the text I/O functions to bundle a filestream
   * pointer along with some internal information.
   *
   * Users of the library should not rely on the contents of this
   * structure, it is for private use of the text I/O layer.
   */
  typedef struct affy_textio_s
  {
    FILE *fp;
    char *buf;
    int   max_buf_len;
    bool  skip_read;
  } AFFY_TEXTIO;

  /* Forward reference for AFFY_PIXREGION */
  typedef struct affy_cell_s AFFY_CELL;

  /* 
   * An AFFY_PIXREGION defines a subset or window into a grid of pixel
   * values.  **data is a 2D array of the dimensions given.
   */
  typedef struct affy_pixregion_s
  {
    affy_uint32    numrows;
    affy_uint32    numcols;
    /* 
     * If the celfile is loaded, *cell is a backpointer to the
     * parent cell struct, so we can zero out the cache pointer when
     * free'ing.
     */
    AFFY_CELL     *cell;
    unsigned int **data;          /* Must hold at least 16 bits; C89 says so. */
  } AFFY_PIXREGION;

  /*
   * A point somewhere on a cel or pixel map.
   */
  typedef struct affy_point_s
  {
    affy_int16 x;
    affy_int32 y;
  } AFFY_POINT;

  /*
   * A point somewhere on a cel or pixel map.
   * Only used in Calvin CEL file I/O, where both x and y must be 16-bit
   */
  typedef struct affy_point_s16
  {
    affy_int16 x;
    affy_int16 y;
  } AFFY_POINT16;

  /*
   * An individual cell obtained from a CEL file.  A cell consists of some
   * number of pixels which can be accessed directly if the DAT file is
   * available.  
   */
  struct affy_cell_s
  {
    double          value;        /* Intensity (mean pixel intensity).     */
#ifdef STORE_CEL_QC
    double          stddev;       /* Std. deviation of pixel intensity.    */
    affy_int16      numpixels;    /* Number of pixels composing this cell. */
    AFFY_PIXREGION *pixels;       /* The actual pixels making up the cell. */
#endif
  };

  /*
   * A probe. The basic unit of information, it consists of both a 
   * PM (perfect match) and MM (mismatch) value.
   */
  typedef struct affy_probe_s
  {
    int                     index; /* Unique id for probe                  */
    AFFY_POINT              mm;    /* Location of Mismatch value           */
    AFFY_POINT              pm;    /* Location of Perfect match value      */
    struct affy_probeset_s *ps;    /* Pointer back to parent probe set     */
  } AFFY_PROBE;

  /* 
   * A ProbeSet consists of a set of individual probes. Each probe set
   * will eventually consist of a single expression value computed from
   * the various pm/mm probes.
   */
  typedef struct affy_probeset_s
  {
    int         index;            /* Unique id for probe set                */
    char       *name;             /* A text description of probe set        */
    int         numprobes;        /* Total number of probes                 */
    AFFY_PROBE *probe;            /* The probes themselves                  */
  } AFFY_PROBESET;

  /*
   * CDF file definitions. The CDF contains the meta-information about
   * a particular microarray chip, so there are many definitions in this
   * section. They link all of the various higher level details together.
   *
   * The CDF file consists of mappings from x,y coordinates into probes,
   * which are grouped by probe sets.
   */
  typedef struct affy_cdffile_s
  {
    char          *array_type;    /* Name of file/chip type                  */
    affy_uint32    numrows;       /* Row/col dimensions.                     */
    affy_uint32    numcols;
    affy_int32     numprobes;     /* Count of total probes                   */
    affy_int32     numexclusions; /* #probesets to exclude from IRON fit     */
    affy_int32     numspikeins;   /* #probesets to leave unnormalized        */
    affy_int32     numprobesets;  /* Number of probe sets                    */
    affy_int32     numqcunits;    /* Num quality control units               */
    affy_uint8   **cell_type;     /* What kind of cell (normal or QC)        */
    affy_uint8   **seen_xy;       /* Have we already parsed this x,y coord?  */
#ifdef STORE_XY_REF
    AFFY_PROBE  ***xy_ref;        /* row,col map to individual probe         */
#endif
    AFFY_PROBESET *probeset;      /* An array of probe sets, made of probes  */
    AFFY_PROBE   **probe;         /* A linear array of probes                */
    char         **exclusions;    /* An array of probe/probeset name strings */
    char         **spikeins;      /* An array of probe/probeset name strings */
    affy_int8      no_mm_flag;    /* set if CDF file is missing MM probes    */
    affy_int8      dupe_probes_flag;  /* probes shared between probesets   */
  } AFFY_CDFFILE;

  /*
   * Given the above definitions, a CEL file can be a simple entity: a
   * matrix of data. The CDF provides all other necessary mapping
   * (probes->probeset).
   */
  typedef struct affy_celfile_s
  {
    char        *filename;
    affy_int32   numrows;
    affy_int32   numcols;
    affy_uint32  nummasks;
    affy_uint32  numoutliers;
    AFFY_CELL  **data;
    bitstr_t   **mask;
    bitstr_t   **outlier;
    char         corrupt_flag;
  } AFFY_CELFILE;

  /*
   * DAT file definition; a binary representation of individual
   * pixel intensity values.  The raw data for the higher level files.
   */
  typedef struct affy_datfile_s
  {
    char          *experiment_name;
    affy_uint16    pixel_width;
    affy_uint16    pixel_height;
    affy_uint16    scanspeed;
    double         temperature;
    double         laser_power;
    char           timestamp[19];           /* field width 18, plus NUL */
    affy_uint32    numpixels;
    affy_uint32    minpixel;
    affy_uint32    maxpixel;
    affy_uint32    numsamples_dc_offset;
    AFFY_POINT     grid_ul;
    AFFY_POINT     grid_ur;
    AFFY_POINT     grid_ll;
    AFFY_POINT     grid_lr;
    affy_uint16    cellmargin;
    char          *scannerid;
    char          *probe_array_type;
    double         meanpixel;
    double         std_dev_pixel;
    double         avg_dc_offset;
    double         std_dev_dc_offset;
    AFFY_PIXREGION pixels;
  } AFFY_DATFILE;

  /*
   * We need a higher abstraction than a celfile, to pull everything 
   * together and provide a convenient storage place for expressions, etc.
   *
   * This is a CHIP, or more generally (and usefully) a collection of one
   * or more CHIPS.
   *
   * A single chip consists of pointers to the meta-information and
   * actual information (CEL file). The summary information is stored
   * in a single array probe_set.
   */
  typedef struct affy_chip_s
  {
    char         *filename;       /* Filename of chip                   */
    AFFY_CDFFILE *cdf;            /* The chip description file          */
    AFFY_CELFILE *cel;            /* For the CEL file information       */
    AFFY_DATFILE *dat;            /* Underlying pixel data (if present) */

    /* The true information - the probe sets */
    int     numprobesets;
    double *probe_set;

    /* P/A call p-values */
    double *probe_set_call_pvalue;

    /* A convenience ptr: not globally used (RMA) */
    double *pm;
  } AFFY_CHIP;

  /*
   * And finally, a chipset is a group of the same type of chips. 
   */
  typedef struct affy_chipset_s
  {
    /* Some common parameters needed, from the CDF file. */
    unsigned int  max_chips;
    unsigned int  num_chips;       /* Index of the next empty spot. */
    affy_uint32   numrows;
    affy_uint32   numcols;
    /* Chip/array type designator (forms the filestem of the CDF file) */
    char         *array_type;
    AFFY_CDFFILE *cdf;
    AFFY_CHIP   **chip;
    /* probeset affinities and t_value for reuse of Median Polish parameters */
    double      **affinities;
    double       *t_values;
    char          mp_allocated_flag;
    char          mp_populated_flag;
  } AFFY_CHIPSET;

  /* 
   * These definitions attempt to model the internal structure of 
   * the new Affymetrix "Calvin" format, which is a self-describing,
   * binary file containing name/value pairs along with some metadata.
   */

  typedef enum 
    {
      AFFY_CALVIN_BYTE = 0,
      AFFY_CALVIN_UBYTE,
      AFFY_CALVIN_SHORT,
      AFFY_CALVIN_USHORT,
      AFFY_CALVIN_INT,
      AFFY_CALVIN_UINT,
      AFFY_CALVIN_FLOAT,
      AFFY_CALVIN_DOUBLE,
      AFFY_CALVIN_STRING,
      AFFY_CALVIN_WSTRING,
      AFFY_CALVIN_UNKNOWN
    } AFFY_CALVIN_DATA_TYPE;

  typedef union affy_calvin_data_s
  {
    affy_int8    byte_val;
    affy_uint8   ubyte_val;
    affy_int16   short_val;
    affy_uint16  ushort_val;
    affy_int32   int_val;
    affy_uint32  uint_val;
    affy_float32 float_val;
    affy_float64 double_val;
    char        *string_val;
  } AFFY_CALVIN_DATA;

  typedef struct affy_calvin_param_s
  {
    char                  *name;
    AFFY_CALVIN_DATA       value;
    AFFY_CALVIN_DATA_TYPE  type;
  } AFFY_CALVIN_PARAM;

  typedef struct affy_calvin_column_s
  {
    char                 *name;
    AFFY_CALVIN_DATA_TYPE type;
    affy_int32            size;
  } AFFY_CALVIN_COLUMN;

  typedef struct affy_calvin_dataset_s
  {
    char       *name;
    affy_uint32 num_params;
    affy_uint32 num_cols;
    affy_uint32 num_rows;

    /* Internal use only */
    affy_uint32 cols_read;
    affy_uint32 rows_read;

    AFFY_CALVIN_COLUMN  *columns;
    AFFY_CALVIN_PARAM   *params;
    AFFY_CALVIN_DATA   **data;
  } AFFY_CALVIN_DATASET;

  typedef struct affy_calvin_dataset_io_s
  {
    AFFY_CALVIN_DATASET *metadata;
    affy_uint32          initial_offset;
    affy_uint32          row_length;
    AFFY_CALVINIO       *calvin_io;
  } AFFY_CALVIN_DATASET_IO;

  typedef struct affy_calvin_column_mapping_s
  {
    const char *name;
    size_t      offset;
  } AFFY_CALVIN_COLUMN_MAPPING;

  typedef struct affy_calvin_datagroup_s
  {
    char       *name;
    affy_uint32 num_datasets;

    AFFY_CALVIN_DATASET **datasets;
  } AFFY_CALVIN_DATAGROUP;

  typedef struct affy_calvin_dataheader_s
  {
    char *type_identifier;
    char *file_identifier;
    char *timestamp;   /* XXX keep as string or convert? */
    char *locale;

    affy_uint32 num_params;
    affy_uint32 num_parent_headers;

    struct affy_calvin_dataheader_s **parent_headers;
    AFFY_CALVIN_PARAM               *params;
  } AFFY_CALVIN_DATAHEADER;

  typedef struct affy_calvin_fileheader_s
  {
    affy_uint8 file_version;     /* File format version (== 1 normally) */
    affy_int32 num_datagroups;   /* Number of data groups */
  } AFFY_CALVIN_FILEHEADER;

  /* 
   * This simple top-level struct ties together all the elements of a
   * Calvin file.
   */
  typedef struct affy_calvin_container_s
  {
    AFFY_CALVIN_FILEHEADER  *file_header;   /* Container info */
    AFFY_CALVIN_DATAHEADER  *data_header;   /* Extended container metadata */
    AFFY_CALVIN_DATAGROUP  **data_groups;   /* Array of data groups */
  } AFFY_CALVIN_CONTAINER;


  /**************************************************************************/

  /** Function prototypes. **/

  /* Memory management. */
  typedef int affy_mempool;

  affy_mempool *affy_pool_create(void);
  void          affy_pool_destroy(affy_mempool *pool);
  void         *affy_pool_alloc(affy_mempool *pool, size_t len);
  void          affy_pool_free(void *mem);
  void          affy_pool_attach(affy_mempool *child, affy_mempool *parent);

  /* File I/O stuff. */
  AFFY_CDFFILE          *affy_load_cdf_file(char *chip_type,
                                            char *dir,
                                            AFFY_COMBINED_FLAGS *f,
                                            AFFY_ERROR *err);
  AFFY_CDFFILE          *affy_load_cdf_file_byname(char *filename,
                                                   char *chip_type,
                                                   AFFY_ERROR *err);
  void                   affy_load_binary_cdf_file(FILE *fp,
						   AFFY_CDFFILE *cdf,
                                                   LIBUTILS_PB_STATE *pbs,
						   AFFY_ERROR *err);
  void                   affy_load_binary_cel_file(FILE *fp,
						   AFFY_CELFILE *cf,
                                                   LIBUTILS_PB_STATE *pbs,
						   AFFY_ERROR *err);
  void                   affy_write_binary_cel_file(FILE *fp,
                                                    AFFY_CHIP *cp,
                                                    AFFY_ERROR *err);
  void                   affy_load_calvin_cel_file(FILE *fp,
						   AFFY_CELFILE *cf,
                                                   LIBUTILS_PB_STATE *pbs,
						   AFFY_ERROR *err);
  void                   affy_load_text_cel_file(FILE *fp, 
                                                 AFFY_CELFILE *cf, 
                                                 LIBUTILS_PB_STATE *pbs,
                                                 AFFY_ERROR *err);
  void                   affy_load_calvin_dat_file(FILE *fp,
						   AFFY_DATFILE *df,
                                                   LIBUTILS_PB_STATE *pbs,
						   AFFY_ERROR *err);
  void                   affy_load_binary_dat_file(FILE *fp,
						   AFFY_DATFILE *df,
                                                   LIBUTILS_PB_STATE *pbs,
						   AFFY_ERROR *err);
  void                   affy_load_text_cdf_file(FILE *fp, 
                                                 AFFY_CDFFILE *cdf, 
                                                 LIBUTILS_PB_STATE *pbs,
                                                 AFFY_ERROR *err);
  AFFY_CDFFILE          *create_blank_generic_cdf(unsigned int max_chips,
                                                  unsigned int numprobes,
                                                  AFFY_ERROR *err);
  void                   affy_load_chipset_single(AFFY_CHIPSET *cs, 
                                                  char *pathname,
                                                  bool ignore_chip_mismatch,
                                                  AFFY_ERROR *err);
  void                   affy_load_chipset(AFFY_CHIPSET *cs, 
                                           char **filelist,
                                           bool ignore_chip_mismatch);
  AFFY_CHIPSET          *affy_create_chipset(unsigned int max_chips,
                                             char *chip_type,
                                             char *cdf_hint,
                                             AFFY_COMBINED_FLAGS *f,
                                             AFFY_ERROR *err);
  AFFY_CHIPSET          *create_blank_generic_chipset(unsigned int max_chips,
                                                      unsigned int numprobes,
                                                      AFFY_ERROR *err);
  AFFY_CHIPSET          *affy_resize_chipset(AFFY_CHIPSET *cs,
                                             unsigned int max_chips,
                                             AFFY_ERROR *err);
  AFFY_CHIPSET          *affy_clone_chipset(AFFY_CHIPSET *cur_chip,
                                            AFFY_ERROR *err);
  AFFY_CHIPSET          *affy_clone_chipset_one_chip(AFFY_CHIPSET *cur_chip,
                                                     int chip_idx,
                                                     AFFY_ERROR *err);
  AFFY_CHIP             *affy_clone_chip(AFFY_CHIP *cur_chip, AFFY_ERROR *err);
  AFFY_CHIP             *affy_load_chip(char *filename, AFFY_ERROR *err);
  AFFY_CELFILE          *affy_load_cel_file(char *filename, AFFY_ERROR *err);
  AFFY_DATFILE          *affy_load_dat_file(char *filename, AFFY_ERROR *err);
  void                   affy_free_cel_file(AFFY_CELFILE *cf);
  void                   affy_mostly_free_cel_file(AFFY_CELFILE *cf);
  void                   affy_free_dat_file(AFFY_DATFILE *df);
  void                   affy_free_cdf_file(AFFY_CDFFILE *cdf);
  void                   affy_free_chip(AFFY_CHIP *ch);
  void                   affy_mostly_free_chip(AFFY_CHIP *ch);
  void                   affy_free_chipset(AFFY_CHIPSET *cs);
  void                   affy_free_calvin_datagroup(AFFY_CALVIN_DATAGROUP *dg);
  void                   affy_free_calvin_fileheader(AFFY_CALVIN_FILEHEADER *fh);
  void                   affy_free_calvin_dataset(AFFY_CALVIN_DATASET *ds);
  void                   affy_free_calvin_column(AFFY_CALVIN_COLUMN *col);
  void                   affy_free_calvin_dataheader(AFFY_CALVIN_DATAHEADER *dh);
  void                   affy_free_calvin_container(AFFY_CALVIN_CONTAINER *cc);
  char                 **affy_list_files(char *directory, 
					 char *extension, 
					 AFFY_ERROR *err);
  void                   print_corrupt_chips_to_stderr(AFFY_CHIPSET *cs);

  void                   affy_load_exclusions_file(char *filename,
                                                   AFFY_CDFFILE *cdf,
                                                   void *mempool,
                                                   AFFY_ERROR *err);
  void                   affy_load_spikeins_file(char *filename,
                                                 AFFY_CDFFILE *cdf,
                                                 void *mempool,
                                                 AFFY_ERROR *err);

  /* generic spreadsheet functions */
  void get_generic_spreadsheet_bounds(char *filename,
                                      affy_uint32 *return_max_rows,
                                      affy_uint32 *return_max_cols,
                                      AFFY_ERROR *err);

  /* Calvin-related functions. */
  AFFY_CALVINIO *affy_calvinio_init(FILE *fp, AFFY_ERROR *err);
  void           affy_calvinio_free(AFFY_CALVINIO *cio);
  AFFY_CALVIN_CONTAINER *affy_calvin_load_container(AFFY_CALVINIO *cio, 
                                                    AFFY_ERROR *err);
  AFFY_CALVIN_DATASET_IO *affy_calvin_prepare_dataset(AFFY_CALVINIO *cio,
                                                      affy_uint32 dg_index,
                                                      affy_uint32 ds_index,
                                                      AFFY_ERROR *err);
  void affy_calvin_close_dataset(AFFY_CALVIN_DATASET_IO *dio);
  void affy_calvin_read_data(AFFY_CALVINIO *cio,
                             AFFY_CALVIN_DATA *dest,
                             AFFY_CALVIN_DATA_TYPE type,
                             AFFY_ERROR *err);
  affy_int32 affy_calvin_find_datagroup_index(AFFY_CALVINIO *cio,
                                              const char *datagroup_name,
                                              AFFY_ERROR *err);
  affy_int32 affy_calvin_find_dataset_index(AFFY_CALVINIO *cio,
                                            affy_uint32 dg_index,
                                            const char *dataset_name,
                                            AFFY_ERROR *err);
  affy_int32 affy_calvin_find_column_index(AFFY_CALVIN_DATASET_IO *dio,
                                           const char *column_name,
                                           AFFY_ERROR *err);
  AFFY_CALVIN_DATAGROUP *affy_calvin_get_datagroup_metadata(AFFY_CALVINIO *cio,
                                                            affy_uint32 dg_index,
                                                            AFFY_ERROR *err);
  AFFY_CALVIN_DATASET *affy_calvin_get_dataset_metadata(AFFY_CALVINIO *cio,
                                                        affy_uint32 dg_index,
                                                        affy_uint32 ds_index,
                                                        AFFY_ERROR *err);
  AFFY_CALVIN_FILEHEADER *affy_calvin_get_file_metadata(AFFY_CALVINIO *cio,
                                                        AFFY_ERROR *err);
  AFFY_CALVIN_DATAHEADER *affy_calvin_get_dataheader(AFFY_CALVINIO *cio,
                                                     AFFY_ERROR *err);
  AFFY_CALVIN_PARAM *affy_calvin_find_param(AFFY_CALVIN_PARAM *params,
                                            affy_uint32 num_params,
                                            const char *param_name);
  void affy_calvin_read_dataset_col(AFFY_CALVIN_DATASET_IO *dio,
                                    LIBUTILS_PB_STATE *pbs,
                                    affy_uint32 col_index,
                                    void *dest,
                                    AFFY_ERROR *err);
  void affy_calvin_read_dataset_rows(AFFY_CALVIN_DATASET_IO *dio,
                                     LIBUTILS_PB_STATE *pbs,                  
                                     affy_uint32 start_row,
                                     affy_uint32 num_rows,
                                     void *base,
                                     size_t base_sz,
                                     const AFFY_CALVIN_COLUMN_MAPPING *offsets,
                                     AFFY_ERROR *err);

  /* interface TBD */
  /* void affy_calvin_read_dataset(); */

  /* Text handling functions. */
  char        *affy_textio_get_next_line(AFFY_TEXTIO *tf);
  void         affy_textio_unget_next_line(AFFY_TEXTIO *tf);
  void         affy_textio_reset_next_line(AFFY_TEXTIO *tf);
  void         affy_textio_skip_to_next_header(AFFY_TEXTIO *tf);
  AFFY_TEXTIO *affy_textio_init(FILE *fp, AFFY_ERROR *err);
  void         affy_textio_free(AFFY_TEXTIO *tf);

  /* Binary handling functions. */

  /* Unadorned readXX() functions to perform no swapping. */
  int   affy_read8(FILE *input, void *buf);
  int   affy_read16(FILE *input, void *buf);
  int   affy_read32(FILE *input, void *buf);
  int   affy_read64(FILE *input, void *buf);
  int   affy_readchars(FILE *input, char *buf, size_t numbytes);
  int   affy_write8(FILE *output, void *buf);
  int   affy_write16(FILE *output, void *buf);
  int   affy_write32(FILE *output, void *buf);
  int   affy_write64(FILE *output, void *buf);
  int   affy_writechars(FILE *output, char *buf);

  /* Functions which explicitly handle endian swapping. */
  int   affy_read16_le(FILE *input, void *buf);
  int   affy_read32_le(FILE *input, void *buf);
  int   affy_read64_le(FILE *input, void *buf);
  int   affy_read16_be(FILE *input, void *buf);
  int   affy_read32_be(FILE *input, void *buf);
  int   affy_read64_be(FILE *input, void *buf);
  int   affy_write16_le(FILE *output, void *buf);
  int   affy_write32_le(FILE *output, void *buf);
  int   affy_write64_le(FILE *output, void *buf);
  int   affy_write16_be(FILE *output, void *buf);
  int   affy_write32_be(FILE *output, void *buf);
  int   affy_write64_be(FILE *output, void *buf);

  /* "Vectorized" reading function for convenience. */
  int   affy_readmulti(FILE *fp, const char *fmt, ...);
        
  /* Utility functions */
  void            affy_mean_normalization(AFFY_CHIPSET *d, double target_mean,
                                          AFFY_COMBINED_FLAGS *f);
  void            affy_median_normalization(AFFY_CHIPSET *d, double target_mean,
                                          AFFY_COMBINED_FLAGS *f);
  char           *affy_get_cdf_name(const char *buf, AFFY_ERROR *err);
  char           *affy_get_cdf_name_from_cel(const char *filename, 
                                             AFFY_ERROR *err);
  double        **affy_matrix_from_cel(AFFY_CELFILE *cf, AFFY_ERROR *err);
  const char     *affy_strerror(AFFY_ERROR_TYPE err);
  AFFY_ERROR     *affy_get_default_error(void);
  void            affy_clone_error(AFFY_ERROR *e1, AFFY_ERROR *e2);
  
  void            affy_dump_dat_hdr(AFFY_DATFILE *df);
  void            affy_dump_calvin_container(AFFY_CALVIN_CONTAINER *cc);
  void            affy_print_calvin_param(AFFY_CALVIN_PARAM *cp);
  void            affy_print_calvin_value(AFFY_CALVIN_DATA data, 
                                          AFFY_CALVIN_DATA_TYPE type);
  void            affy_write_expressions(AFFY_CHIPSET *c, 
                                         char         *filename, 
                                         unsigned int  opts,
                                         AFFY_ERROR   *err);
  void            affy_write_probe_values(AFFY_CHIPSET *c, 
                                          char *filename, 
                                          int opts,
                                          AFFY_ERROR *err);
  void            affy_write_expressions_gct(AFFY_CHIPSET *c, 
                                             char *filename, 
                                             AFFY_ERROR *err);
  bool            affy_ismasked(AFFY_CHIP *chip, int x, int y);
  bool            affy_isoutlier(AFFY_CHIP *chip, int x, int y);
  bool            affy_is_control_probe(AFFY_PROBE *p_probe);
  bool            affy_is_control_probeset(AFFY_PROBESET *ps);
  bool            affy_is_control_string(char *string);
  bool            affy_isqc(AFFY_CHIP *chip, int x, int y);
  bool            affy_isundefined(AFFY_CHIP *chip, int x, int y);
  AFFY_PIXREGION *affy_pixels_from_cell(AFFY_CHIP *cp, 
                                        int x, 
                                        int y, 
                                        AFFY_ERROR *err);
  AFFY_POINT      affy_cell_to_pixel(AFFY_CHIP *cp, int x, int y);

  void            affy_pixregion2tiff(AFFY_PIXREGION *p, 
                                      char *filename, 
                                      AFFY_ERROR *err);
  void            affy_pixregion2text(AFFY_PIXREGION *p, 
                                      char *filename, 
                                      AFFY_ERROR *err);
  void            affy_write_pixel_region(AFFY_PIXREGION *pr,
                                          char *filename, 
                                          AFFY_ERROR *err);
  void            affy_free_pixregion(AFFY_PIXREGION *pr);
  void            print_flags(AFFY_COMBINED_FLAGS *f, char *output_file_name);

  /* Statistical functions (from util). */
  double affy_median_save(double *x, int length, AFFY_COMBINED_FLAGS *f,
                          AFFY_ERROR *err);
  double affy_median(double *x, int length, AFFY_COMBINED_FLAGS *f);
  double affy_mean(double *x, int length);
  double affy_mean_geometric_floor_1(double *x, int length);
  void   affy_get_row_median(double **z, 
                             double *rdelta, 
                             int startrow, 
                             int startcol, 
                             int numrows, 
                             int numcolumns,
                             AFFY_COMBINED_FLAGS *f,
                             AFFY_ERROR *err);
  void   affy_get_column_median(double **z, 
                                double *cdelta, 
                                int startrow, 
                                int startcol, 
                                int numrows, 
                                int numcolumns,
                                AFFY_COMBINED_FLAGS *f,
                                AFFY_ERROR *err);
  int    affy_median_sort(const void *p1, const void *p2);
  int    affy_qnorm_compare(const void *p1, const void *p2);
  void   affy_rank_order(double *rank, double **x, int n);
  double calculate_pearson_r_float(float *array1, float *array2, int n);
  double calculate_pearson_r_float_skip_weak(float *array1, float *array2,
                                             int n);
  double calculate_pearson_r_double(double *array1, double *array2, int n);
  void   affy_quantile_normalization(AFFY_CHIPSET *d, affy_uint8 pm_only,
                                     AFFY_ERROR *err);
  void   affy_quantile_normalization_probeset(AFFY_CHIPSET *d,
                                              AFFY_ERROR *err);
  void   affy_pairwise_normalization(AFFY_CHIPSET *d,
                                     AFFY_CHIP *model_chip,
                                     unsigned int opts,
                                     AFFY_COMBINED_FLAGS *f,
                                     AFFY_ERROR *err);
  void   affy_pairwise_normalization_probeset(AFFY_CHIPSET *d,
                                              AFFY_CHIP *model_chip,
                                              int unlog_flag,
                                              AFFY_COMBINED_FLAGS *f,
                                              AFFY_ERROR *err);
  void   affy_floor_probe(AFFY_CHIPSET *cs,
                          double floor_value,
                          AFFY_ERROR *err);
  void   affy_floor_probeset(AFFY_CHIPSET *cs,
                             double floor_value,
                             AFFY_ERROR *err);
  void   affy_floor_probeset_to_min_non_zero(AFFY_CHIPSET *cs,
                                             AFFY_ERROR *err);
  void   affy_floor_probeset_non_zero_to_one(AFFY_CHIPSET *cs,
                                             AFFY_ERROR *err);

  void   affy_unlog_probeset(AFFY_CHIPSET *cs, AFFY_ERROR *err);
  void   fill_normalization_scales(char *filestem,
                                   double *signals1,
                                   double *signals2,
                                   double *signals2_scales,
                                   char *mask_array,
                                   int num_spots,
                                   double rank_frac_cutoff,
                                   double rank_frac_cutoff2,
                                   int condense_training_flag,
                                   AFFY_COMBINED_FLAGS *f,
                                   double *return_training_frac,
                                   double *return_rmsd,
                                   AFFY_ERROR *err);

  void   affy_pnorm_both(double x, 
                         double *cum, 
                         double *ccum, 
                         int i_tail,
                         int log_p);
  double affy_pnorm5(double x, 
                     double mu, 
                     double sigma, 
                     int lower_tail,
                     int log_p);
  void   affy_kernel_density(double *x, 
                             int nx, 
                             double *weights, 
                             double *dy,
                             double *dx, 
                             int N,
                             AFFY_ERROR *err);
  double affy_max_density(double *x, int n, AFFY_ERROR *err);
  double affy_trunc(double x);

  /*
   * These error handling macros are adapted from the GNU Scientific
   * Library.
   *
   * Macro arguments should be:
   *
   *     descr:   string containing a reason/description of the fault
   *     errvar:  name of an AFFY_ERROR pointer
   *     errtype: type of error that occurred (see AFFY_ERROR_TYPE)
   *     retval: when errors are propagated back up, and the surrounding 
   *             function returns a value, this is the value to return
   */
#define AFFY_HANDLE_ERROR(desc, errtype, errvar, retval)        \
  do {                                                          \
    assert(errvar != NULL);                                     \
    errvar->type      = errtype;                                \
    errvar->timestamp = time(NULL);                             \
    errvar->descr     = desc;                                   \
    errvar->module    = __FILE__;                               \
    errvar->location  = __LINE__;                               \
    if (errvar->handler != NULL)                                \
      (errvar->handler)(errvar);                                \
    return (retval);                                            \
  } while (0)

#define AFFY_HANDLE_ERROR_VOID(desc, errtype, errvar)   \
  do {                                                  \
    assert(errvar != NULL);                             \
    errvar->type      = errtype;                        \
    errvar->timestamp = time(NULL);                     \
    errvar->descr     = desc;                           \
    errvar->module    = __FILE__;                       \
    errvar->location  = __LINE__;                       \
    if (errvar->handler != NULL)                        \
      (errvar->handler)(errvar);                        \
    return;                                             \
  } while (0)

#define AFFY_HANDLE_ERROR_VOID_ZERO(desc, errtype, errvar)   \
  do {                                                       \
    assert(errvar != NULL);                                  \
    errvar->type      = errtype;                             \
    errvar->timestamp = time(NULL);                          \
    errvar->descr     = desc;                                \
    errvar->module    = __FILE__;                            \
    errvar->location  = __LINE__;                            \
    if (errvar->handler != NULL)                             \
      (errvar->handler)(errvar);                             \
    return 0;                                                  \
  } while (0)

#define AFFY_HANDLE_ERROR_GOTO(desc, errtype, errvar, label)    \
  do {                                                          \
    assert(errvar != NULL);                                     \
    errvar->type      = errtype;                                \
    errvar->timestamp = time(NULL);                             \
    errvar->descr     = desc;                                   \
    errvar->module    = __FILE__;                               \
    errvar->location  = __LINE__;                               \
    if (errvar->handler != NULL)                                \
      (errvar->handler)(errvar);	                        \
    goto label;                                                 \
  } while (0)

#define AFFY_CHECK_ERROR_VOID(errvar)           \
  if (errvar->type != AFFY_ERROR_NONE) return;

#define AFFY_CHECK_ERROR_VOID_ZERO(errvar)           \
  if (errvar->type != AFFY_ERROR_NONE) return 0;

#define AFFY_CHECK_ERROR(errvar, retval)                \
  if (errvar->type != AFFY_ERROR_NONE) return (retval);
       
#define AFFY_CHECK_ERROR_GOTO(errvar, label)            \
  if (errvar->type != AFFY_ERROR_NONE) goto label;  

#ifdef __cplusplus
};
#endif

#endif
