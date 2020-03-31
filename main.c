#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "fftw3.h"
#include "libmseed.h"

#include "fft.h"
#include "standard_deviation.h"

static void
usage ()
{
  printf ("Usage: ./ms2ampmax <mseedfile>\n");
  printf ("\nOutput format: \n");
  printf ("an image of spectogram\n");
}

static void
dumpdata (double *data, uint64_t totalSamples, FILE *fptr)
{
  uint64_t i;
  for (i = 0; i < totalSamples; i++)
  {
    fprintf (fptr, "%lf\n", data[i]);
  }
}

static void demean (double *data, uint64_t totalSamples, double demean);
static void plotFrequencyResponse (const char *inputfile, double maxFrequency);

int
main (int argc, char **argv)
{
  char *mseedfile    = NULL;
  MS3TraceList *mstl = NULL;
  MS3TraceID *tid    = NULL;
  MS3TraceSeg *seg   = NULL;
  double sampleRate;
  uint32_t flags = 0;
  int8_t verbose = 0;
  int rv;

  const char *fftoutputFile = "fftoutput.txt";

#ifdef DUMPDATA
  FILE *fptr;
  fptr = fopen ("dumpdata.txt", "w");
  if (fptr == NULL)
  {
    printf ("Error!");
    return -1;
  }
#endif
  FILE *output = fopen (fftoutputFile, "w");
  if (output == NULL)
  {
    printf ("Error opening!");
    return -1;
  }

  /* Simple argument parsing */
  if (argc != 2)
  {
    usage ();
    return -1;
  }
  mseedfile = argv[1];

  /* Set flags to validate CRCs while reading */
  flags |= MSF_VALIDATECRC;
  /* Unpack the data */
  flags |= MSF_RECORDLIST;

  /* Read all miniSEED into a trace list */
  rv = ms3_readtracelist (&mstl, mseedfile, NULL, 0, flags, verbose);
  if (rv != MS_NOERROR)
  {
    ms_log (2, "Cannot read miniSEED from file: %s\n", ms_errorstr (rv));
    return -1;
  }

  /* Traverse trace list structures */
  tid = mstl->traces;
  while (tid)
  {
    /* Allocate the data array of every trace */
    double *data;
    uint64_t totalSamples = 0;
    /* Count how many samples of this trace */
    seg = tid->first;
    while (seg)
    {
      totalSamples += seg->samplecnt;
      seg = seg->next;
    }
#ifdef DEBUG
    printf ("estimated samples of this trace: %" PRId64 "\n", totalSamples);
#endif

    /* Get the data of this trace */
    seg            = tid->first;
    uint64_t index = 0;
    data           = (double *)malloc (sizeof (double) * totalSamples);
    if (data == NULL)
    {
      ms_log (2, "something wrong when mallocing data array\n");
      return -1;
    }
    while (seg)
    {
      /* Get sampling rate */
      sampleRate = seg->samprate;

      /* Unpack and get the data */
      if (seg->recordlist && seg->recordlist->first)
      {
        void *sptr;
        size_t i;
        uint8_t samplesize;
        char sampletype;
        /* Determine sample size and type based on encoding of first record */
        ms_encoding_sizetype (seg->recordlist->first->msr->encoding, &samplesize, &sampletype);

        /* Unpack data samples using record list.
         * No data buffer is supplied, so it will be allocated and assigned to the segment.
         * Alternatively, a user-specified data buffer can be provided here. */
        int64_t unpacked = mstl3_unpack_recordlist (tid, seg, NULL, 0, verbose);

        if (unpacked != seg->samplecnt)
        {
          ms_log (2, "Cannot unpack samples for %s\n", tid->sid);
        }
        else
        {
          if (sampletype == 'a')
          {
            printf ("%*s",
                    (seg->numsamples > INT_MAX) ? INT_MAX : (int)seg->numsamples,
                    (char *)seg->datasamples);
          }
          else
          {
            for (i = 0; i < seg->numsamples; i++)
            {
              sptr = (char *)seg->datasamples + (i * samplesize);
              if (sampletype == 'i')
              {
                data[index] = (double)(*(int32_t *)sptr);
              }
              else if (sampletype == 'f')
              {
                data[index] = (double)(*(float *)sptr);
              }
              else if (sampletype == 'd')
              {
                data[index] = (double)(*(double *)sptr);
              }

              index++;
            }
          }
        }
      }

      seg = seg->next;
    }

    /* FFT and print */
    double mean, SD;
    getMeanAndSD (data, totalSamples, &mean, &SD);
    demean (data, totalSamples, mean);
#ifdef DUMPDATA
    dumpdata (data, totalSamples, fptr);
#endif
    fftToFileHalf (data, totalSamples, sampleRate, output);

    //plotFrequencyResponse (fftoutputFile, sampleRate / 2.);

    free (data);
    tid = tid->next;
  }

  /* Make sure everything is cleaned up */
  if (mstl)
    mstl3_free (&mstl, 0);
#ifdef DUMPDATA
  fclose (fptr);
#endif
  fclose (output);
  return 0;
}

static void
demean (double *data, uint64_t totalSamples, double mean)
{
  uint64_t i;
  for (i = 0; i < totalSamples; i++)
  {
    data[i] -= mean;
  }
}

/* Plot Frequency Response using `gnuplot` */
static void
plotFrequencyResponse (const char *inputfile, double maxFrequency)
{
  char config[] =
      "set title 'Frequency Domain'\n"
      "set ylabel 'Amplitude'\n"
      "set xlabel 'Frequency [Hz]'\n";

  /* Set title, xlabel, and ylabel */
  FILE *plot;
  plot = popen ("gnuplot -p", "w");
  fprintf (plot, config);

  /* Read inputfile (fftoutput.txt) and plot half of the sample points */
  fprintf (plot, "inputfile = '%s'\n", inputfile);
  fprintf (plot, "stats inputfile nooutput\n");
  fprintf (plot, "set xrange[:%lf]\n", maxFrequency);

  fprintf (plot, "plot \
            inputfile using ($1):(abs($2)/STATS_records)\n");

  pclose (plot);
}
