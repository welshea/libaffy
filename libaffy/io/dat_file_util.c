
/**************************************************************************
 *
 * Filename:  dat_file_util.c
 * 
 * Purpose:   Binary DAT file utility routines.
 *
 * Creation:  04/04/05
 *
 * Author:    Andrew M. Hoerter
 *
 * Copyright: Copyright (C) 2007, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 04/04/05: Creation (AMH)
 * 04/08/05: Reorder & clarify DAT header dump output, fix pixeldata to be
 *           2D (AMH)
 * 04/19/05: Use an AFFY_PIXREGION now (AMH)
 * 03/11/08: New error handling scheme (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 *
 **************************************************************************/

#include <affy.h>

/* 
 * affy_dump_dat_hdr(): Print a DAT header in human-readable form.
 *
 * Inputs: *dat contains an initialized DAT file header.
 * Outputs: None.
 * Side effects: Print header data to stdout.
 */
void affy_dump_dat_hdr(AFFY_DATFILE *dat)
{
  assert(dat != NULL);

  /* Dump out all the known information. */
  printf("DAT header information: \n");
  printf("\tExperiment name: %s\n", dat->experiment_name);
  printf("\tProbe array type: %s\n", dat->probe_array_type);
  printf("\tPixels per line: %" AFFY_PRNu32 "\n", dat->pixels.numcols);
  printf("\tNumber of lines: %" AFFY_PRNu32 "\n", dat->pixels.numrows);

/*
  printf("\tPixels per row (CLS=): %" AFFY_PRNu32 "\n", dat->pixels_per_row);
  printf("\tNumber of rows (RWS=): %" AFFY_PRNu32 "\n", dat->numrows);
 */
  printf("\tPixel width: %" AFFY_PRNu16 "\n", dat->pixel_width);
  printf("\tPixel height: %" AFFY_PRNu16 "\n", dat->pixel_height);
  printf("\tTotal number of pixels: %" AFFY_PRNu32 "\n", dat->numpixels);
  printf("\tMinimum pixel intensity: %" AFFY_PRNu32 "\n", dat->minpixel);
  printf("\tMaximum pixel intensity: %" AFFY_PRNu32 "\n", dat->maxpixel);
  printf("\tMean pixel intensity: %f\n", dat->meanpixel);
  printf("\tStandard deviation of pixel intensity: %f\n", dat->std_dev_pixel);
  printf("\tScan speed: %" AFFY_PRNu16 "\n", dat->scanspeed);
  printf("\tTemperature in degrees C: %f\n", dat->temperature);
  printf("\tLaser power reading: %f\n", dat->laser_power);
  printf("\tTime of scan: %s\n", dat->timestamp);
  printf("\tCell margin: %" AFFY_PRNu16 "\n", dat->cellmargin);
  printf("\tScanner ID: %s\n", dat->scannerid);
  printf("\tUpper-left grid coordinates: %" AFFY_PRNd16 " %" AFFY_PRNd32
         "\n", dat->grid_ul.x, dat->grid_ul.y);
  printf("\tUpper-right grid coordinates: %" AFFY_PRNd16 " %" AFFY_PRNd32
         "\n", dat->grid_ur.x, dat->grid_ur.y);
  printf("\tLower-left grid coordinates: %" AFFY_PRNd16 " %" AFFY_PRNd32
         "\n", dat->grid_ll.x, dat->grid_ll.y);
  printf("\tLower-right grid coordinates: %" AFFY_PRNd16 " %" AFFY_PRNd32
         "\n", dat->grid_lr.x, dat->grid_lr.y);
  printf("\t# DC offset samples: %" AFFY_PRNu32 "\n",
         dat->numsamples_dc_offset);
  printf("\tAverage DC offset: %f\n", dat->avg_dc_offset);
  printf("\tStandard deviation of DC offset: %f\n",
         dat->std_dev_dc_offset);
}

/* 
 * affy_free_dat_file(): Free a DAT file header and any associated storage.
 *
 * Inputs: *dat contains an initialized DAT file header.
 * Outputs: None.
 * Side effects: Associated memory is freed.
 */
void affy_free_dat_file(AFFY_DATFILE *dat)
{
  h_free(dat);
}

