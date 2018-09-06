/*
 * This function performs integration of Langevin equation
 */

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <omp.h>

#include "struct_def.h"
#include "array_def.h"

/* Definitions for routines */
double ForceSys(sys_var , run_param , int );

void	IntegrateSys(sys_var *s, run_param r, gsl_rng * gr)
{

  int     j, k, sind;
  float   *xi, *thet;
  float   fr, *gam, m, *sig;
  float   *q_t, *qdot_t, *c_t, *f_t, q_tp1, qdot_tp1, f_tp1;
  float   teng;
  FILE    *outfp;
  char    outfname[40];

  /* Allocating memory */
  xi 		  = CreateVector(r.nmodes);
  thet    = CreateVector(r.nmodes);
  gam 		= CreateVector(r.nmodes);
  sig 		= CreateVector(r.nmodes);
  q_t 		= CreateVector(r.nmodes);
  qdot_t 	= CreateVector(r.nmodes);
  c_t 		= CreateVector(r.nmodes);
  f_t 		= CreateVector(r.nmodes);

  #pragma omp parallel for private(j,sind)
  for (j=0; j<s->NmodeIRcou; j++)
  {
    sind           = (int) gsl_matrix_get(s->IRs_pqrcountmat, j, 0);
    gam[sind-1]    = gsl_vector_get(s->gamvec, sind-1);
    sig[sind-1]    = gsl_vector_get(s->sigvec, sind-1);

    xi[sind-1]     = gsl_ran_gaussian(gr, 1.0);
    thet[sind-1]   = gsl_ran_gaussian(gr, 1.0);
    //printf("%f\t%f\n", xi, thet);
  }
  
  #pragma omp parallel for private(j,sind,q_tp1)
  for (j=0; j<s->NmodeIRcou; j++)
  {
    sind      = (int) gsl_matrix_get(s->IRs_pqrcountmat, j, 0);
    q_t[sind-1]    = gsl_vector_get(s->qvec, sind-1);
    qdot_t[sind-1] = gsl_vector_get(s->qdotvec, sind-1);
    f_t[sind-1]    = ForceSys(*s, r, j);
    //printf("%f\n", f_t[sind-1]);

    c_t[sind-1]    = r.dt*r.dt/2.0 * (f_t[sind-1] - gam[sind-1]*qdot_t[sind-1] ) + 
                sig[sind-1]*pow(r.dt, 1.5)*( 0.5*xi[sind-1] + 1.0/(2.0*pow(3.0, 0.5))*thet[sind-1]);

    q_tp1     	   = q_t[sind-1] + r.dt*qdot_t[sind-1] + c_t[sind-1];
    gsl_vector_set (s->qvec, sind-1, q_tp1);
  }
    
  #pragma omp parallel for private(j,sind,qdot_tp1,f_tp1)
  for (j=0; j<s->NmodeIRcou; j++)
  { 
    sind      = (int) gsl_matrix_get(s->IRs_pqrcountmat, j, 0);
    f_tp1     = ForceSys(*s, r, j);
    
    qdot_tp1  = qdot_t[sind-1] + r.dt/2.0*( f_tp1 + f_t[sind-1] ) - 
                r.dt*gam[sind-1]*qdot_t[sind-1] + sig[sind-1]*pow(r.dt, 0.5)*xi[sind-1] - gam[sind-1]*c_t[sind-1];

    gsl_vector_set (s->qdotvec, sind-1, qdot_tp1);
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
  double systeng = 0.0;
  if ((r.step)%r.nfreq==0) 
  {
    printf("# STEP %10d OF %10d IS %2.2f \%\n", r.step, r.nsteps, (float) r.step*100.0/r.nsteps);
    // outfp = fopen("modenormaldist.txt", "a");
    // fprintf(outfp, "%f\t%f\n", xi, thet);
    // fclose(outfp);
    
    for (j=0; j<s->NmodeIRcou; j++)
    {
      sind      = (int) gsl_matrix_get(s->IRs_pqrcountmat, j, 0);
      sprintf(outfname, "modedisp.%04d.txt", sind);
      outfp = fopen(outfname, "a");
      fprintf(outfp, "%5.5e\n", gsl_vector_get(s->qvec, sind-1));
      fclose(outfp);

      sprintf(outfname, "modevelc.%04d.txt", sind);
      outfp = fopen(outfname, "a");
      fprintf(outfp, "%5.5e\n", gsl_vector_get(s->qdotvec, sind-1));
      fclose(outfp);

      m         = gsl_vector_get(s->mvec, sind-1);
      fr        = gsl_vector_get(s->frvec, sind-1);
      teng      = 0.5*m*gsl_vector_get(s->qdotvec, sind-1)*gsl_vector_get(s->qdotvec, sind-1) +
      0.5*m*pow(2*PI*fr, 2.0)*gsl_vector_get(s->qvec, sind-1)*gsl_vector_get(s->qvec, sind-1);
      systeng   += teng;

      sprintf(outfname, "modeteng.%04d.txt", sind);
      outfp = fopen(outfname, "a");
      fprintf(outfp, "%5.5e\n", teng);
      fclose(outfp);      

      // sprintf(outfname, "modenormaldist.1.txt", r.step);
      // outfp = fopen(outfname, "a");
      // fprintf(outfp, "%f\t%f\n", xi, thet);
      // fclose(outfp); 
    }
    sprintf(outfname, "modetotteng.txt");
    outfp = fopen(outfname, "a");
    fprintf(outfp, "%5.5e\n", systeng);
    fclose(outfp);
  }

  return;
}
