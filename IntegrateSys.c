/*
 * This function performs integration of Langevin equation
 */

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <omp.h>

#include "struct_def.h"

/* Definitions for routines */
double ForceSys(sys_var s, run_param r, int j);

void	IntegrateSys(sys_var *s, run_param r)
{

  int     j, k;
  double  xi, thet;
  double  fr, gam, m, sig;
  double  q_t, qdot_t, f_t, c_t, q_tp1, f_tp1, qdot_tp1;
  double  teng;
  FILE    *outfp;
  char    outfname[40];

  const   gsl_rng_type * gT;
  gsl_rng * gr;
  gsl_rng_env_setup();
  gT    = gsl_rng_default;
  gr    = gsl_rng_alloc(gT);

  for (j=0; j<r.nmodes; j++)
  {
    gam       = gsl_vector_get(s->gamvec, j);
    sig       = gsl_vector_get(s->sigvec, j);

    xi        = gsl_ran_gaussian(gr, 1.0);
    thet      = gsl_ran_gaussian(gr, 1.0);

    q_t       = gsl_vector_get(s->qvec, j);
    qdot_t    = gsl_vector_get(s->qdotvec, j);
    f_t       = ForceSys(*s, r, j);
    //printf("%f\n", f_t);

    c_t       = r.dt*r.dt/2.0 * (f_t - gam*qdot_t ) + 
    sig*pow(r.dt, 1.5)*( 0.5*xi + 1.0/(2.0*pow(3.0, 0.5))*thet);

    q_tp1     = q_t + r.dt*qdot_t + c_t;
    gsl_vector_set (s->qvec, j, q_tp1);
    
    f_tp1     = ForceSys(*s, r, j);
    
    qdot_tp1  = qdot_t + r.dt/2.0*( f_tp1 + f_t ) - 
    r.dt*gam*qdot_t + sig*pow(r.dt, 0.5)*xi - gam*c_t;

    gsl_vector_set (s->qdotvec, j, qdot_tp1);

    // Writing in a file !!
    if ((r.step+1)%r.nfreq==0) 
    {
      sprintf(outfname, "modedisp.%04d.txt", j);
      outfp = fopen(outfname, "a");
      fprintf(outfp, "%5.5e\n", gsl_vector_get(s->qvec, j));
      fclose(outfp);

      sprintf(outfname, "modevelc.%04d.txt", j);
      outfp = fopen(outfname, "a");
      fprintf(outfp, "%5.5e\n", gsl_vector_get(s->qdotvec, j));
      fclose(outfp);

      m         = gsl_vector_get(s->mvec, j);
      fr        = gsl_vector_get(s->frvec, j);
      teng      = 0.5*m*gsl_vector_get(s->qdotvec, j)*gsl_vector_get(s->qdotvec, j) +
      0.5*m*pow(2*PI*fr, 2.0)*gsl_vector_get(s->qvec, j)*gsl_vector_get(s->qvec, j);

      sprintf(outfname, "modeteng.%04d.txt", j);
      outfp = fopen(outfname, "a");
      fprintf(outfp, "%5.5e\n", teng);
      fclose(outfp);

    }
  }
  if ((r.step)%r.nfreq==0) 
  {
    printf("# STEP %10d OF %10d\n", r.step, r.nsteps);
  }

  gsl_rng_free(gr);
  return;
}
