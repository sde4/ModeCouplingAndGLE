/*
 * This function performs integration of Langevin equation
 */

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <omp.h>

#include "struct_def.h"
#include "array_def.h"

/* Definitions for routines */
double ForceSys(sys_var s, run_param r, int j);

void	IntegrateSys(sys_var *s, run_param r, gsl_rng * gr)
{

  int     j, k;
  float   *xi, *thet;
  float   fr, *gam, m, *sig;
  float   *q_t, *qdot_t, *c_t, *f_t, q_tp1, qdot_tp1, f_tp1;
  float   teng;
  FILE    *outfp;
  char    outfname[40];

  /* Allocating memory */
  xi 		= CreateVector(r.nmodes);
  thet 		= CreateVector(r.nmodes);
  gam 		= CreateVector(r.nmodes);
  sig 		= CreateVector(r.nmodes);
  q_t 		= CreateVector(r.nmodes);
  qdot_t 	= CreateVector(r.nmodes);
  c_t 		= CreateVector(r.nmodes);
  f_t 		= CreateVector(r.nmodes);

  #pragma omp parallel for
  for (j=0; j<r.nmodes; j++)
  {
    gam[j]    = gsl_vector_get(s->gamvec, j);
    sig[j]    = gsl_vector_get(s->sigvec, j);

    xi[j]     = gsl_ran_gaussian(gr, 1.0);
    thet[j]   = gsl_ran_gaussian(gr, 1.0);
    //printf("%f\t%f\n", xi, thet);
  }
  
  #pragma omp parallel for private(j,q_tp1)
  for (j=0; j<r.nmodes; j++)
  {
    q_t[j]    = gsl_vector_get(s->qvec, j);
    qdot_t[j] = gsl_vector_get(s->qdotvec, j);
    f_t[j]    = ForceSys(*s, r, j);
    //printf("%f\n", f_t);

    c_t[j]    = r.dt*r.dt/2.0 * (f_t[j] - gam[j]*qdot_t[j] ) + 
                sig[j]*pow(r.dt, 1.5)*( 0.5*xi[j] + 1.0/(2.0*pow(3.0, 0.5))*thet[j]);

    q_tp1     = q_t[j] + r.dt*qdot_t[j] + c_t[j];
    gsl_vector_set (s->qvec, j, q_tp1);
  }
    
  #pragma omp parallel for private(j,qdot_tp1,f_tp1)
  for (j=0; j<r.nmodes; j++)
  { 
    f_tp1     = ForceSys(*s, r, j);
    
    qdot_tp1  = qdot_t[j] + r.dt/2.0*( f_tp1 + f_t[j] ) - 
                r.dt*gam[j]*qdot_t[j] + sig[j]*pow(r.dt, 0.5)*xi[j] - gam[j]*c_t[j];

    gsl_vector_set (s->qdotvec, j, qdot_tp1);
    //printf("Thread %d is doing iteration %d.\n", omp_get_thread_num( ), j);
  }
  // Free up memory
  free(xi);
  free(thet);
  free(gam);
  free(sig);
  free(q_t);
  free(qdot_t);
  free(c_t);
  free(f_t);

  // Writing in a file !!
  if ((r.step)%r.nfreq==0) 
  {
    printf("# STEP %10d OF %10d IS %2.2f \%\n", r.step, r.nsteps, (float) r.step*100.0/r.nsteps);
    // outfp = fopen("modenormaldist.txt", "a");
    // fprintf(outfp, "%f\t%f\n", xi, thet);
    // fclose(outfp);
    
    for (j=0; j<9; j++)
    {
      sprintf(outfname, "modedisp.%04d.txt", j+1);
      outfp = fopen(outfname, "a");
      fprintf(outfp, "%5.5e\n", gsl_vector_get(s->qvec, j));
      fclose(outfp);

      sprintf(outfname, "modevelc.%04d.txt", j+1);
      outfp = fopen(outfname, "a");
      fprintf(outfp, "%5.5e\n", gsl_vector_get(s->qdotvec, j));
      fclose(outfp);

      m         = gsl_vector_get(s->mvec, j);
      fr        = gsl_vector_get(s->frvec, j);
      teng      = 0.5*m*gsl_vector_get(s->qdotvec, j)*gsl_vector_get(s->qdotvec, j) +
      0.5*m*pow(2*PI*fr, 2.0)*gsl_vector_get(s->qvec, j)*gsl_vector_get(s->qvec, j);

      sprintf(outfname, "modeteng.%04d.txt", j+1);
      outfp = fopen(outfname, "a");
      fprintf(outfp, "%5.5e\n", teng);
      fclose(outfp);
    
      // sprintf(outfname, "modenormaldist.1.txt", r.step);
      // outfp = fopen(outfname, "a");
      // fprintf(outfp, "%f\t%f\n", xi, thet);
      // fclose(outfp); 
    }

  }

  return;
}
