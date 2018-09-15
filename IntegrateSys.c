/*
 * This function performs integration of Langevin equation
 */

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <omp.h>

#include "struct_def.h"
#include "array_def.h"

/* Definitions for routines */
for_var ForceSys(sys_var , run_param , int );

void	IntegrateSys(sys_var *s, run_param r, gsl_rng * gr, state_var sv)
{

  int     j, k, sind;
  float   *xi, *thet;
  float   fr, *gam, m, *sig, sigma;
  float   q_t, *qdot_t, *c_t, *f_t, *q_tp1, qdot_tp1, f_tp1;
  float   teng;
  for_var f;
  FILE    *outfp, *outfp1;
  char    outfname[40];
  


  /* Allocating memory */
  xi 		= CreateVector(r.nmodes);
  thet    	= CreateVector(r.nmodes);
  gam 		= CreateVector(r.nmodes);
  sig 		= CreateVector(r.nmodes);
  q_tp1 	= CreateVector(r.nmodes);
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
  
 #pragma omp parallel for private(j,sind,q_t,f)
  for (j=0; j<s->NmodeIRcou; j++)
  {
    sind      	   = (int) gsl_matrix_get(s->IRs_pqrcountmat, j, 0);
    f    	   = ForceSys(*s, r, j);
    f_t[sind-1]    = f.f1+f.f3;
    //printf("%f\n", f_t[sind-1]);
    qdot_t[sind-1] = gsl_vector_get(s->qdotvec, sind-1);
    c_t[sind-1]    = r.dt*r.dt/2.0 * (f_t[sind-1] - gam[sind-1]*qdot_t[sind-1] ) + 
                	sig[sind-1]*pow(r.dt, 1.5)*( 0.5*xi[sind-1] + 1.0/(2.0*pow(3.0, 0.5))*thet[sind-1]);
    q_t	    	   = gsl_vector_get(s->qvec, sind-1);
    q_tp1[sind-1]  = q_t + r.dt*qdot_t[sind-1] + c_t[sind-1];

    /*
    // randomizing the displacements
    if (sind!=1){
    fr   ``	   = gsl_vector_get(s->frvec, j);
    m 		   = gsl_vector_get(s->mvec, j);
    sigma 	   = pow(kB * sv.T / (m * pow(2 * PI * fr, 2.0)), 0.5);
    q_tp1 	   = gsl_ran_gaussian(gr, sigma);
    }
    */ 
  } 
  #pragma omp parallel for private(j,sind)
  for (j=0; j<s->NmodeIRcou; j++)
  {
    sind      	   = (int) gsl_matrix_get(s->IRs_pqrcountmat, j, 0);
    gsl_vector_set (s->qvec, sind-1, q_tp1[sind-1]);
  } 
  #pragma omp parallel for private(j,sind,qdot_tp1,f,f_tp1)
  for (j=0; j<s->NmodeIRcou; j++)
  { 
    sind      = (int) gsl_matrix_get(s->IRs_pqrcountmat, j, 0);
    f         = ForceSys(*s, r, j);
    f_tp1     = f.f1+f.f3;
    gsl_vector_set (s->fvec, sind-1, f_tp1);                                //// writing just the nonlinear part
    gsl_vector_set (s->enonvec, sind-1, f.ep4);
    
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
  free(q_tp1);
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
    
    // trajectory file !!
    outfp1 = fopen("modetraj.txt", "a");
    fprintf(outfp1, "%d\n", r.step);
    double systeng = 0.0;
    for (j=0; j<s->NmodeIRcou; j++){
    //for (j=0; j<4; j++){
      sind      =  (int) gsl_matrix_get(s->IRs_pqrcountmat, j, 0);
      
      // sprintf(outfname, "modedisp.%04d.txt", sind);
      // outfp = fopen(outfname, "a");
      // fprintf(outfp, "%5.5e\n", gsl_vector_get(s->qvec, sind-1));
      // fclose(outfp);

      // sprintf(outfname, "modevelc.%04d.txt", sind);
      // outfp = fopen(outfname, "a");
      // fprintf(outfp, "%5.5e\n", gsl_vector_get(s->qdotvec, sind-1));
      // fclose(outfp);

      // sprintf(outfname, "modeforc.%04d.txt", sind);
      // outfp = fopen(outfname, "a");
      // f = ForceSys(*s, r, i);
      // gsl_vector_set(s->fvec, sind-1, f.f3);          // storing just the nonlinear part
      // fprintf(outfp, "%5.5e\n", gsl_vector_get(s->fvec, sind-1));
      // fclose(outfp);

      q_t      = gsl_vector_get(s->qvec, sind-1);
      qdot_tp1 = gsl_vector_get(s->qdotvec, sind-1);
      f_tp1    = gsl_vector_get (s->fvec, sind-1)/q_t;
      // f_tp1    = pow(-gsl_vector_get (s->fvec, sind-1)/q_t, 0.5)/(2*PI);
      m        = gsl_vector_get(s->mvec, sind-1);
      fr       = gsl_vector_get(s->frvec, sind-1);
      teng     = 0.5*m*qdot_tp1*qdot_tp1 + 0.5*m*2*PI*fr*2*PI*fr*q_t*q_t;
      systeng += teng +  m*gsl_vector_get(s->enonvec, sind-1);

      // sprintf(outfname, "modeteng.%04d.txt", sind);
      // outfp = fopen(outfname, "a");
      // fprintf(outfp, "%5.5e\n", teng);
      // fclose(outfp);
      fprintf(outfp1, "%d %5.5e %5.5e %5.5e %5.5e\n", sind, q_t, qdot_tp1, f_tp1, teng);


      // sprintf(outfname, "modenormaldist.1.txt", r.step);
      // outfp = fopen(outfname, "a");
      // fprintf(outfp, "%f\t%f\n", xi, thet);
      // fclose(outfp); 
    }
    outfp = fopen("modetotteng.txt", "a");
    fprintf(outfp, "%5.5e\n", systeng);
    fclose(outfp);
    fclose(outfp1);
  }

  return;
}
