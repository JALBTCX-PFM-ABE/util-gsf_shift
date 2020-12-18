
/*********************************************************************************************

    This is public domain software that was developed by or for the U.S. Naval Oceanographic
    Office and/or the U.S. Army Corps of Engineers.

    This is a work of the U.S. Government. In accordance with 17 USC 105, copyright protection
    is not available for any work of the U.S. Government.

    Neither the United States Government, nor any employees of the United States Government,
    nor the author, makes any warranty, express or implied, without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, or assumes any liability or
    responsibility for the accuracy, completeness, or usefulness of any information,
    apparatus, product, or process disclosed, or represents that its use would not infringe
    privately-owned rights. Reference herein to any specific commercial products, process,
    or service by trade name, trademark, manufacturer, or otherwise, does not necessarily
    constitute or imply its endorsement, recommendation, or favoring by the United States
    Government. The views and opinions of authors expressed herein do not necessarily state
    or reflect those of the United States Government, and shall not be used for advertising
    or product endorsement purposes.

*********************************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>

#include "nvutility.h"

#include "gsf.h"

#include "version.h"


extern int32_t gsfError;


void usage ()
{
    fprintf (stderr, 
        "\nUsage: gsf_shift[-s SHIFT -z DEPTH_CORRECTOR] <GSF_FILE> \n");
    fprintf (stderr, "\nWhere:\n\n");
    fprintf (stderr, "\t-s and -z are mutually exclusive:\n\n");
    fprintf (stderr, "\t-s = Shift the data along the track by SHIFT seconds (time)\n");
    fprintf (stderr, 
             "\t-z = Correct the depth by adding DEPTH_CORRECTOR to the depth (stored in GSF depth corrector field)\n\n");
    fflush (stderr);
}



int32_t main (int32_t argc, char **argv)
{
    gsfDataID           gsfID;
    gsfRecords          rec;
    float               offset = 0.0, corrector = 0.0, min_depth, max_depth;
    double              dist, az, lat, lon, mult, off;
    int32_t             hnd, outhnd, percent = 0, old_percent = -1, ret, i, ioffset, eof, current;
    uint8_t             option_s = NVFalse, option_z = NVFalse;
    char                comment[16384], infile[512], outfile[512], c;
    uint8_t             cFlag;
    extern char         *optarg;
    extern int          optind;


    int32_t write_history (int32_t, char **, char *, char *, int32_t);
    void newgp (double, double, double, double, double *, double *);


    printf("\n\n %s \n\n",VERSION);


    while ((c = getopt (argc, argv, "s:z:")) != EOF)
      {
	switch (c)
          {
          case 's':
            sscanf (optarg, "%f", &offset);
            option_s = NVTrue;
            break;

          case 'z':
            sscanf (optarg, "%f", &corrector);
            option_z = NVTrue;
            break;

          default:
            usage ();
            exit (-1);
            break;
          }
      }


    /*  There's no point in running this without an option.  */

    if ((option_s && option_z) || !(option_s || option_z))
      {
        usage ();
        exit (-1);
      }


    /* Make sure we got the mandatory file name argument.  */

    if (optind >= argc)
      {
        usage ();
        exit (-1);
      }


    strcpy (infile, argv[optind]);

    sprintf (outfile, "gsf_shift_%d.tmp", getpid ());


    /* Create the output file */

    ret = gsfOpen (outfile, GSF_CREATE, &outhnd);
    if (ret)
      {
	fprintf (stderr, "Unable to create temporary output file: %s\n", outfile);
	fprintf (stderr, "File: %s not processed\n", infile);
        fflush (stderr);
	return (-1);
      }


    /* Try to open the specified gsf file */

    ret = gsfOpen (infile, GSF_UPDATE, &hnd);
    if (ret)
      {
	fprintf (stderr, "Could not open GSF file %s\n", infile);
	gsfPrintError (stderr);
        fflush (stderr);
	return (-1);
    }


    printf("File : %s\n\n", infile);


    /* This is the main processing loop, where we read each record, apply the correctors, and write each record. 
       Note that all pings in the file are processed, even if the ignore ping flag has been set.  This is done so
       that if later the edited data are viewed, they will have had the same corrections applied. All beams, except
       those with the ignore beam flag set but no reason specified are processed. (Ignore beam with no reason mask
       indicates that there was no detection made by the sonar.)  */

    eof = 0;
    while (!eof)
      {
	ret = gsfRead (hnd, GSF_NEXT_RECORD, &gsfID, &rec, NULL, 0);
	if (ret < 0)
          {
	    if (gsfError != GSF_READ_TO_END_OF_FILE)
              {
		fprintf (stderr, "Error: %d, reading gsf file: %s", gsfError, infile);
                fflush (stderr);
		return (-1);
	    }
	    eof = 1;
	    continue;
	}

	if (gsfID.recordID == GSF_RECORD_SWATH_BATHYMETRY_PING)
          {
            if (option_s)
              {
                dist = (rec.mb_ping.speed * 1852.0 / 3600.0) * fabsf (offset);

                if (offset < 0.0)
                  {
                    az = rec.mb_ping.heading + 180.0;
                    if (az > 360.0) az -= 360.0;
                  }
                else
                  {
                    az = rec.mb_ping.heading;
                  }

                lat = rec.mb_ping.latitude;
                lon = rec.mb_ping.longitude;

                newgp (lat, lon, az, dist, &rec.mb_ping.latitude, &rec.mb_ping.longitude);
              }
            else
              {
                gsfGetScaleFactor (hnd, GSF_SWATH_BATHY_SUBRECORD_DEPTH_ARRAY, &cFlag, &mult, &off);


                min_depth = 99999.9;
                max_depth = -99999.9;

                if (!(rec.mb_ping.ping_flags & GSF_IGNORE_PING))
                  {
                    for (i = 0 ; i < rec.mb_ping.number_beams ; i++)
                      {
                        if (rec.mb_ping.beam_flags != (uint8_t *) NULL)
                          {
                            /* If NV_GSF_IGNORE_BEAM is set, and nothing else is set then there was no detection made by the 
                               sonar, correct all beams but a "NULL" beam.  It makes sense to correct edited beams in case 
                               they are un-edited later.  */

                            if (!check_flag (rec.mb_ping.beam_flags[i], NV_GSF_IGNORE_NULL_BEAM))
                              {
                                /* Make sure this data has a depth field */

                                if (rec.mb_ping.depth != (double *) NULL)
                                  {
                                    rec.mb_ping.depth[i] += corrector;

                                    if (rec.mb_ping.depth[i] < min_depth) min_depth = rec.mb_ping.depth[i];
                                    if (rec.mb_ping.depth[i] > max_depth) max_depth = rec.mb_ping.depth[i];
                                  }


                                /* Correct the nominal depth field if it is present */

                                if (rec.mb_ping.nominal_depth != (double *) NULL) 
                                  rec.mb_ping.nominal_depth[i] += corrector;
                              }
                          }
                      }
                    rec.mb_ping.depth_corrector += corrector;
                  }


                /*  Before we write the record we need to scan it for negative values (drying heights) that may have 
                    been created when we applied the depth corrector.  We also want to make sure that we haven't pushed
                    a value over the precision limit of the scale factors.  This requires a change to the scale factors.
                    We will use a scale factor of 0.005 m for depths of 326 meters or less, 0.01 m for depths of 654 
                    meters or less, 0.1 m for depths to 6552 meters, and 1.0 m for anything over that.  JCD  */

                if (min_depth != 99999.9)
                  {
                    if (min_depth < 0.0) 
                      {
                        ioffset = -((int32_t) min_depth + 1);
                      }
                    else
                      {
                        ioffset = 0;
                      }

                    if (max_depth < 327.0) 
                      {
                        mult = 0.005;
                      }
                    else if (max_depth < 655.0)
                      {
                        mult = 0.01;
                      }
                    else if (max_depth < 6553.0)
                      {
                        mult = 0.1;
                      }
                    else
                      {
                        mult = 1.0;
                      }

                    i = gsfLoadScaleFactor (&rec.mb_ping.scaleFactors, GSF_SWATH_BATHY_SUBRECORD_DEPTH_ARRAY, cFlag, 
                                            mult, ioffset);
                    if (i) 
                      {
                        fprintf (stderr, "%s %s %d - Error loading scale factors: %d\n", __FILE__, __FUNCTION__, __LINE__, i);
                        fflush (stderr);
                      }
                  }
              }
          }


	/* Write the current record to the output file */

	ret = gsfWrite (outhnd, &gsfID, &rec);
	if (ret < 0)
          {
	    fprintf(stderr, "Error: %d, writting file: %s", gsfError, outfile);
            fflush (stderr);
          }


	/* Print a percent complete message */

	current = gsfPercent (hnd);
	if (old_percent != current)
          {
	    fprintf (stdout," %02d%% complete\r", current);
	    fflush (stdout);
	    old_percent = current;
          }
      }


    percent = 100;
    printf ("%3d%% processed    \n", percent);
        
    gsfClose(hnd);
    printf("\n");


    /*  Open the file non-indexed so that we can write a history record.  */

    if (gsfOpen (outfile, GSF_UPDATE, &hnd))
    {
        gsfPrintError (stderr);
        exit (-1);
    }


    /*  Write a history record describing the kludge.  */

    if (option_s)
      {
        sprintf (comment, "Hacked with a %f second shifted position (here be dragons)\n", offset);
      }
    else
      {
        sprintf (comment, 
                 "Added a depth corrector of %f meters to the depths (in addition to any preexisting corrector)\n", 
                 corrector);
      }


    ret = write_history (argc, argv, comment, outfile, hnd);
    if (ret) 
      {
        fprintf(stderr, "Error: %d - writing gsf history record\n", ret);
        fflush (stderr);
      }

    gsfClose (hnd);

    ret = rename (outfile, infile);
    if (ret)
      {
	fprintf (stderr, "Error updating file: %s\n", infile);
	fprintf (stderr, "Results remain in file: %s\n", outfile);
        fflush (stderr);
	return (-1);
      }




    /*  We need to remove the index file because changing scale factors can cause relocation of records.  */

    if (option_z)
      {
        strcpy (comment, outfile);
        comment[strlen (comment) - 3] = 'n';
        remove (comment);

        strcpy (outfile, infile);
        outfile[strlen (outfile) - 3] = 'n';
        remove (outfile);
      }


    return (0);
}
