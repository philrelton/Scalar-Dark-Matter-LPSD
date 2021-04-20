/********************************************************************************
    lpsd.c
			  
    2003, 2004 by Michael Troebs, mt@lzh.de and Gerhard Heinzel, ghh@mpq.mpg.de

    calculate spectra from time series using discrete Fourier 
    transforms at frequencies equally spaced on a logarithmic axis
    
    lpsd does everything except user interface and data output
    
 ********************************************************************************/
#define SINCOS
#define FAST 1


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <fftw3.h>
#include "config.h"
#include "ask.h"
#include "IO.h"
#include "genwin.h"
#include "debug.h"
#include "lpsd.h"
#include "misc.h"
#include "errors.h"

/*
20.03.2004: http://www.caddr.com/macho/archives/iolanguage/2003-9/549.html
*/
#ifndef __linux__
#include <windows.h>
struct timezone
{
  int tz_minuteswest;
  int tz_dsttime;
};
void
gettimeofday (struct timeval *tv, struct timezone *tz
	      __attribute__ ((unused)))
{
  long int count = GetTickCount ();
  tv->tv_sec = (int) (count / 1000);
  tv->tv_usec = (count % 1000) * 1000;
}

#else
#include <sys/time.h>		/* gettimeofday, timeval */
#endif

#ifdef __linux__
extern double round (double x);
/*
gcc (GCC) 3.3.1 (SuSE Linux) gives warning: implicit declaration of function `round'
without this line - math.h defines function round, but its declaration seems to be
missing in math.h 
*/

#else
/*
int round (double x) {
  return ((int) (floor (x + 0.5)));
}
*/
#endif

/********************************************************************************
 * 	global variables						   	
 ********************************************************************************/
static int nread;
static double winsum;
static double winsum2;
static double nenbw;		/* normalized equivalent noise bandwidth */
static double *dwin;		/* pointer to window function for FFT */

/********************************************************************************
 * 	functions								
 ********************************************************************************/

/* 
	copies nfft values from data to segm
	if drift removal is selected, a linear regression of data is 
	performed and the subtracted values are copied to segm
*/
static void
remove_drift (double *segm, double *data, int nfft, int LR)
{
  int i;
  long double sx, sy, stt, sty, xm, t;
double a,b;
  if (LR == 2)
    {				/* subtract straight line through first and last point */
      a = data[0];
      b = data[nfft - 1] - data[0] / (double) (nfft - 1.0);
      for (i = 0; i < nfft; i++)
	{
	  segm[i] = data[i] - (a + b * i);
	}
    }
  else if (LR == 1)
    {				/* linear regression */

      sx = sy = 0;
      for (i = 0; i < nfft; i++)
	{
	  sx += i;
	  sy += data[i];
	}
      xm = sx / nfft;
      stt = sty = 0;
      for (i = 0; i < nfft; i++)
	{
	  t = i - xm;
	  stt += t * t;
	  sty += t * data[i];
	}
      b = sty / stt;
      a = (sy - sx * b) / nfft;
      for (i = 0; i < nfft; i++)
	{
	  segm[i] = data[i] - (a + b * i);
	}
    }
  else if (LR == 0)
    {				/* copy data */
      for (i = 0; i < nfft; i++)
	{
	  segm[i] = data[i];
	}
    }
}

static void
remove_drift2 (double *a, double *b, double *data, int nfft, int LR)
{
  int i;
  long double sx, sy, stt, sty, xm, ndbl;

  if (LR == 2)
    {				/* subtract straight line through first and last point */
      *a = data[0];
      *b = data[nfft - 1] - data[0] / (double) (nfft - 1.0);
    }
  else if (LR == 1)
    {				/* linear regression */
      ndbl = (long double) nfft;
      xm = (ndbl - 1.0L) / 2.0L;
      sx = ndbl * xm;
      stt = (ndbl * ndbl - 1.0L) * ndbl / 12.0L;
      sy = sty = 0.L;
      for (i = 0; i < nfft; i++)
	{
	  sy += data[i];
	  sty += (i - xm) * data[i];
	}
      *b = sty / stt;
      *a = (sy - sx * *b) / nfft;
    }
  else if (LR == 0)
    {				/* copy data */
      *a = 1.0;
      *b = 0.0;
    }
}


/********************************************************************************
 *	calculates DFT 
 *		
 *	Parameters
 *		nfft	dimension of fft
 *		bin	bin to be calculated
 *		rslt	array for DFT as spectral density and spectrum
 *			and variance
 *			rslt[0]=PSD, rslt[1]=variance(PSD) 
 *			rslt[2]=PS rslt[3]=variance(PS)
 ********************************************************************************/
static void
getDFT (int nfft, double bin, double fsamp, double ovlp, int LR, double *rslt,
	int *avg)
{
  double *dwincs;		/* pointer to array containing window function*cos,window function*sin */
  int i, j;
  double dft_re, dft_im;	/* real and imaginary part of DFT */
  int start;			/* first index in data array */
  double *data;			/* start address of data */
  double dft2;			/* sum of real part squared and imag part squared */
  int nsum;			/* number of summands */
  double *segm;			/* contains data of one segment without drift */

  double west_q, west_r, west_temp;
  double west_sumw;		/* temp variable for West's averaging */
  double west_m, west_t;

  /* calculate window function */
  dwincs = (double *) xmalloc (2 * nfft * sizeof (double));
  assert (dwincs != 0);

  makewinsincos (nfft, bin, dwincs, &winsum, &winsum2, &nenbw);

  data = get_data ();
  assert (data != 0);

  segm = (double *) xmalloc (nfft * sizeof (double));
  assert (segm != 0);

  /* remove drift from first data segment */
  remove_drift (&segm[0], &data[0], nfft, LR);

  start = 0;
  dft2 = 0.;
  nsum = 1;

  /* calculate first DFT */
  dft_re = dft_im = 0;
  for (i = 0, j = 0; i < nfft; i++, j += 2)
    {
      dft_re += dwincs[j] * segm[i];
      dft_im += dwincs[j + 1] * segm[i];
    }
  dft2 = dft_re * dft_re + dft_im * dft_im;
  west_sumw = 1.;
  west_m = dft2;
  west_t = 0.;

  start += nfft * (1.0 - (double) (ovlp / 100.));	/* go to next segment */
  /* process other segments if available */
  while (start + nfft < nread)
    {
      remove_drift (&segm[0], &data[start], nfft, LR);

      /* calculate DFT */
      dft_re = dft_im = 0;
      for (i = 0, j = 0; i < nfft; i++, j += 2)
	{
	  dft_re += dwincs[j] * segm[i];
	  dft_im += dwincs[j + 1] * segm[i];
	}
      dft2 = dft_re * dft_re + dft_im * dft_im;

      west_q = dft2 - west_m;
      west_temp = west_sumw + 1.;
      west_r = west_q / west_temp;
      west_m += west_r;
      west_t += west_r * west_sumw * west_q;
      west_sumw = west_temp;

      nsum++;
      start += nfft * (1.0 - (double) (ovlp / 100.));	/* go to next segment */
    }

  /* return result */
  rslt[0] = west_m;
  /* if only one DFT has been computed, then stddev equals DFT 
     otherwise, divide variance by n-1, then take root
   */
  if (nsum > 2)
    rslt[1] = sqrt (west_t / ((double) nsum - 1.));
  else
    rslt[1] = rslt[0];

  rslt[2] = rslt[0];
  rslt[3] = rslt[1];
  rslt[0] *= 2. / (fsamp * winsum2);	/* power spectral density */
  rslt[1] *= 2. / (fsamp * winsum2);	/* variance of power spectral density */
  rslt[2] *= 2. / (winsum * winsum);	/* power spectrum */
  rslt[3] *= 2. / (winsum * winsum);	/* variance of power spectrum */

  *avg = nsum;

  /* clean up */
  xfree (segm);
  xfree (dwincs);
}

static void
getDFT2 (int nfft, double bin, double fsamp, double ovlp, int LR,
	 double *rslt, int *avg)
{
  double *dwincs;		/* pointer to array containing window function*cos,window function*sin */
  int i;
  double dft_re, dft_im;	/* real and imaginary part of DFT */
  int start;			/* first index in data array */
  double *data;			/* start address of data */
  double dft2;			/* sum of real part squared and imag part squared */
  int nsum;			/* number of summands */
  double y;			/* time series detrended with window */
  double *winp, *datp;

  double total;		/* Running sum of DFTs */

  /* calculate window function */
  dwincs = (double *) xmalloc (2 * nfft * sizeof (double));
  assert (dwincs != 0);

  makewinsincos (nfft, bin, dwincs, &winsum, &winsum2, &nenbw);

  data = get_data ();
  assert (data != 0);

  start = 0;
  dft2 = 0.;
  nsum = 1;

  /* calculate first DFT */
  dft_re = dft_im = 0;
  datp = data;
  winp = dwincs;
  for (i = 0; i < nfft; i++)
    {
      y = *(datp++);
      dft_re += *(winp++) * y;
      dft_im += *(winp++) * y;
    }
  dft2 = dft_re * dft_re + dft_im * dft_im;
  total = dft2;

  start += nfft * (1.0 - (double) (ovlp / 100.));	/* go to next segment */
  /* process other segments if available */
  while (start + nfft < nread)
    {
      /* calculate DFT */
      dft_re = dft_im = 0;
      datp = data + start;
      winp = dwincs;
      for (i = 0; i < nfft; i++)
	{
	  y = *(datp++);
	  dft_re += *(winp++) * y;
	  dft_im += *(winp++) * y;
	}
      dft2 = dft_re * dft_re + dft_im * dft_im;

      total += dft2;
      nsum++;
      start += nfft * (1.0 - (double) (ovlp / 100.));	/* go to next segment */
    }

  /* return result */
  rslt[0] = total / nsum;
  
  /* This sets the variance to zero. This is not true, but we are not using the variance. */
  rslt[1] = 0;
  
  rslt[2] = rslt[0];
  rslt[3] = rslt[1];
  rslt[0] *= 2. / (fsamp * winsum2);	/* power spectral density */
  rslt[1] *= 2. / (fsamp * winsum2);	/* variance of power spectral density */
  rslt[2] *= 2. / (winsum * winsum);	/* power spectrum */
  rslt[3] *= 2. / (winsum * winsum);	/* variance of power spectrum */

  *avg = nsum;

  /* clean up */
  xfree (dwincs);
}


/*
	calculates paramaters for DFTs
	
	input
		nread		number of data
		fsamp		sampling frequency
		ndiv		desired number of entries in spectrum
		sollavg		desired number of averages
	output
		ndiv		actual number of entries in spectrum
		fspec		frequencies in spectrum
		bins		bins for DFTs
		nffts		dimensions for DFTs
 ********************************************************************************
 	Naming convention	source code	publication
				i		j
				fresc		r_{min}
				fresb		r_{avg}
				fresa		r'
				fres		r''
				ndft		L(j)
				bin		m(j)
 ********************************************************************************/
static void
calc_params (tCFG * cfg, tDATA * data)
{
  double fres, f;
  int i, i0, ndft;
  double bin;
  double navg;
  double ovfact, xov, g;

  ovfact = 1. / (1. - (*cfg).ovlp / 100.);
  xov = (1. - (*cfg).ovlp / 100.);
  g = log ((*cfg).fmax / (*cfg).fmin);

  i = (*cfg).nspec * (*cfg).iter;
  i0 = i;
  f = (*cfg).fmin * exp (i * g / ((*cfg).Jdes - 1.));
  while (f <= (*cfg).fmax && i / (*cfg).nspec < (*cfg).iter + 1)
   {
      fres = f * (exp (g / ((*cfg).Jdes - 1.)) - 1);
      ndft = round ((*cfg).fsamp / fres);
      bin = (f / fres);
      navg = ((double) ((nread - ndft)) * ovfact) / ndft + 1;
      (*data).fspec[i - i0] = f;
      (*data).nffts[i - i0] = ndft;
      (*data).bins[i - i0] = bin;
      i++;
      f = (*cfg).fmin * exp (i * g / ((*cfg).Jdes - 1.));
  }
  (*cfg).nspec = i - i0;
  (*cfg).fmin = (*data).fspec[0];
  (*cfg).fmax = (*data).fspec[(*cfg).nspec - 1];
}

void
calculate_lpsd (tCFG * cfg, tDATA * data)
{
  int k;			/* 0..nspec */
  int k_start = 0;		/* N. lines in save file. Post fail start point */
  char ch;			/* For scanning through checkpointing file */
  int Nsave = (*cfg).nspec / 100; /* Frequency of data checkpointing */
  int j; 			/* Iteration variables for checkpointing data */
  FILE * file1;			/* Output file, temp for checkpointing */
  double rslt[4];		/* rslt[0]=PSD, rslt[1]=variance(PSD) rslt[2]=PS rslt[3]=variance(PS) */
  double progress;

  struct timeval tv;
  double start, now, print;

  /* Check output file for saved checkpoint */
  file1 = fopen((*cfg).ofn, "r");
  if (file1){
      while((ch=fgetc(file1)) != EOF){
          if(ch == '\n'){
              k_start++;
          }
      }
  fclose(file1);
  printf("Backup collected. Starting from k = %i\n", k_start);
  }
  else{
      printf("No backup file. Starting from fmin\n");
      k_start = 0;
  }
  printf ("Checkpointing every %i iterations\n", Nsave);
  printf ("Computing output:  00.0%%");
  fflush (stdout);
  gettimeofday (&tv, NULL);
  start = tv.tv_sec + tv.tv_usec / 1e6;
  now = start;
  print = start;
  
  /* Start calculation of LPSD fromi saved checkpoint or zero */
  for (k = k_start; k < (*cfg).nspec; k++)
    {
    if (FAST)
      getDFT2 ((*data).nffts[k], (*data).bins[k], (*cfg).fsamp, (*cfg).ovlp,
	      (*cfg).LR, &rslt[0], &(*data).avg[k]);
else      getDFT ((*data).nffts[k], (*data).bins[k], (*cfg).fsamp, (*cfg).ovlp,
	       (*cfg).LR, &rslt[0], &(*data).avg[k]);
      (*data).psd[k] = rslt[0];
      (*data).varpsd[k] = rslt[1];
      (*data).ps[k] = rslt[2];
      (*data).varps[k] = rslt[3];
      gettimeofday (&tv, NULL);
      now = tv.tv_sec + tv.tv_usec / 1e6;
      if (now - print > PSTEP)
	{
	  print = now;
	  progress = (100 * ((double) k)) / ((double) ((*cfg).nspec));
	  printf ("\b\b\b\b\b\b%5.1f%%", progress);
	  fflush (stdout);
	}

      /* If k is a multiple of Nsave then write data to backup file */
      if(k % Nsave  == 0 && k != k_start){
          file1 = fopen((*cfg).ofn, "a");
          for(j=k-Nsave; j<k; j++){
		fprintf(file1, "%e	", (*data).psd[j]);
		fprintf(file1, "%e	", (*data).ps[j]);
		fprintf(file1, "%d	", (*data).avg[j]);
		fprintf(file1, "\n");
          }
          fclose(file1);
      }
      else if(k == (*cfg).nspec - 1){
          file1 = fopen((*cfg).ofn, "a");
          for(j=Nsave*(k/Nsave); j<(*cfg).nspec; j++){
		fprintf(file1, "%e	", (*data).psd[j]);
		fprintf(file1, "%e	", (*data).ps[j]);
		fprintf(file1, "%d	", (*data).avg[j]);
		fprintf(file1, "\n");
          }
          fclose(file1);
      }
    }
  /* finish */
  printf ("\b\b\b\b\b\b  100%%\n");
  fflush (stdout);
  gettimeofday (&tv, NULL);
  printf ("Duration (s)=%5.3f\n\n", tv.tv_sec - start + tv.tv_usec / 1e6);
}

void
calculate_fftw (tCFG * cfg, tDATA * data)
{
  int nfft;			/* dimension of DFT */
  FILE *wfp;
  fftw_plan plan;
  double *rawdata;		/* start address of data */
  double *out;
  double *segm;			/* contains data of one segment without drift */
  int i, j;
  double d;
  int start;
  double *west_sumw;
  double west_q, west_r, west_temp;
  int navg;
  double *fft_ps, *fft_varps;

  struct timeval tv;
  double stt;

  gettimeofday (&tv, NULL);
  stt = tv.tv_sec + tv.tv_usec / 1e6;

  nfft = (*cfg).nfft;

  dwin = (double *) xmalloc (nfft * sizeof (double));
  segm = (double *) xmalloc (nfft * sizeof (double));
  out = (double *) xmalloc (nfft * sizeof (double));
  west_sumw = (double *) xmalloc ((nfft / 2 + 1) * sizeof (double));
  fft_ps = (double *) xmalloc ((nfft / 2 + 1) * sizeof (double));
  fft_varps = (double *) xmalloc ((nfft / 2 + 1) * sizeof (double));

  /* calculate window function */
  makewin (nfft, 0, dwin, &winsum, &winsum2, &nenbw);

  /* import fftw "wisdom" */
  if ((wfp = fopen ((*cfg).wfn, "r")) == NULL)
    message1 ("Cannot open '%s'", (*cfg).wfn);
  else
    {
      if (fftw_import_wisdom_from_file (wfp) == 0)
	message ("Error importing wisdom");
      fclose (wfp);
    }
  /* plan DFT */
  printf ("Planning...");
  fflush (stdout);

  plan = fftw_plan_r2r_1d (nfft, segm, out, FFTW_R2HC, FFTW_ESTIMATE);
  printf ("done.\n");
  fflush (stdout);

  rawdata = get_data ();
  assert (rawdata != 0);

  printf ("Computing output\n");

  /* remove drift from first data segment */
  remove_drift (&segm[0], &rawdata[0], nfft, (*cfg).LR);
  /* multiply data with window function */
  for (i = 0; i < nfft; i++)
    segm[i] = segm[i] * dwin[i];

  fftw_execute (plan);

  d = 2 * (out[0] * out[0]);
  (*data).fft_ps[0] = d;
  west_sumw[0] = 1.;
  (*data).fft_varps[0] = 0;
  for (j = 1; j < nfft / 2 + 1; j++)
    {
      d = 2 * (out[j] * out[j] + out[nfft - j] * out[nfft - j]);
      (*data).fft_ps[j] = d;
      west_sumw[j] = 1.;
      (*data).fft_varps[j] = 0;
    }
  navg = 1;
  start = nfft * (1.0 - (double) ((*cfg).ovlp / 100.));

  /* remaining segments */
  while (start + nfft < nread)
    {

      printf (".");
      fflush (stdout);
      if (navg % 75 == 0)
	printf ("\n");

      navg++;
      remove_drift (&segm[0], &rawdata[start], nfft, (*cfg).LR);

      /* multiply data with window function */
      for (i = 0; i < nfft; i++)
	segm[i] = segm[i] * dwin[i];

      fftw_execute (plan);

      d = 2 * (out[0] * out[0]);
      west_q = d - (*data).fft_ps[0];
      west_temp = west_sumw[0] + 1;
      west_r = west_q / west_temp;
      (*data).fft_ps[0] += west_r;
      (*data).fft_varps[0] += west_r * west_sumw[0] * west_q;
      west_sumw[0] = west_temp;

      for (j = 1; j < nfft / 2 + 1; j++)
	{
	  d = 2 * (out[j] * out[j] + out[nfft - j] * out[nfft - j]);
	  west_q = d - (*data).fft_ps[j];
	  west_temp = west_sumw[j] + 1;
	  west_r = west_q / west_temp;
	  (*data).fft_ps[j] += west_r;
	  (*data).fft_varps[j] += west_r * west_sumw[j] * west_q;
	  west_sumw[j] = west_temp;
	}
      start += nfft * (1.0 - (double) ((*cfg).ovlp / 100.));	/* go to next segment */
    }

  if (navg > 1)
    {
      for (i = 0; i < nfft / 2 + 1; i++)
	{
	  (*data).fft_varps[i] =
	    sqrt ((*data).fft_varps[i] / ((double) navg - 1));
	}
    }
  else
    {
      for (i = 0; i < nfft / 2 + 1; i++)
	{
	  (*data).fft_varps[i] = (*data).fft_ps[i];
	}
    }
  /* normalizations and additional information */
  j = 0;
  for (i = 0; i < nfft / 2 + 1; i++)
    {
      if (((*cfg).fres * i >= (*cfg).fmin) &&
	  ((*cfg).fres * i <= (*cfg).fmax) && ((*cfg).sbin <= i))
	{
	  (*data).fspec[j] = (*cfg).fres * i;
	  (*data).ps[j] = (*data).fft_ps[i] / (winsum * winsum);
	  (*data).varps[j] = (*data).fft_varps[i] / (winsum * winsum);
	  (*data).psd[j] = (*data).fft_ps[i] / ((*cfg).fsamp * winsum2);
	  (*data).varpsd[j] = (*data).fft_varps[i] / ((*cfg).fsamp * winsum2);
	  (*data).avg[j] = navg;
	  (*data).nffts[j] = nfft;
	  (*data).bins[j] = (double) i;
	  j++;
	}
    }
  printf ("done.\n");

  gettimeofday (&tv, NULL);
  printf ("Duration (s)=%5.3f\n\n", tv.tv_sec - stt + tv.tv_usec / 1e6);

  /* write wisdom to file */
  if ((wfp = fopen ((*cfg).wfn, "w")) == NULL)
    message1 ("Cannot open '%s'", (*cfg).wfn);
  else
    {
      fftw_export_wisdom_to_file (wfp);
      fclose (wfp);
    }
  /* clean up */
  fftw_destroy_plan (plan);

  /* forget wisdom, free memory */
  fftw_forget_wisdom ();
  xfree (fft_ps);
  xfree (fft_varps);
  xfree (west_sumw);
  xfree (dwin);
  xfree (segm);
  xfree (out);
}

/*
	works on cfg, data structures of the calling program
*/
void
calculateSpectrum (tCFG * cfg, tDATA * data)
{
  nread = (*data).nread;

  /* read data file into memory */
  /* and subtract mean data value */
  printf ("\nReading data, subtracting mean...\n");
  nread = floor (((*cfg).tmax - (*cfg).tmin) * (*cfg).fsamp + 1);
  read_file ((*cfg).ifn, (*cfg).ulsb, (*data).mean,
	     (int) ((*cfg).tmin * (*cfg).fsamp), (*data).nread,
	     (*data).comma);

  if ((*cfg).METHOD == 0)
    {
      calc_params (cfg, data);
      calculate_lpsd (cfg, data);
    }
  else if ((*cfg).METHOD == 1)
    {
      calculate_fftw (cfg, data);
    }
}
