/*
 * This function initializes the system of Langevin equation
 */

#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h> 
#include <omp.h>
#include "struct_def.h"

/* Definitions for routines */
for_var ForceSys(sys_var , run_param , int );

void StateRead(sys_var* s, run_param r) {
  printf("\n# 3. State Initialization ...\n");

  int i, j;
  int sind;
  double fr;
  double m;
  double q_t, qdot_t, f_t, teng, fdat;
  for_var f;
  FILE *infp, *outfp, *outfp1;
  char outfname[40];

  printf("#    o Input file required: 1. in.restart\n");
  // check for the required input files

  if ((infp=fopen("read.restart", "r")) == NULL) 
  {
     printf("#      ERROR: File 1 not found! \n");
     exit(1);
  }
  
  /**************************/
  // State initializations !!
  /**************************/
  infp = fopen("read.restart", "r");
  fscanf(infp, "%d", &s->step);
  s->istep = s->step;
  for (i=0; i<s->NmodeIRcou; i++){
    fscanf(infp, "%d %lf %lf %lf %lf", &sind, &q_t, &qdot_t, &fdat, &fdat);
    gsl_vector_set(s->qvec, i, q_t);
    gsl_vector_set(s->qdotvec, i, qdot_t);

  }
  fclose(infp);

  // Initial perturbation energy to one of the modes
  m 	 = gsl_vector_get(s->mvec, r.pertmodind-1);
  // qdot_t = gsl_vector_get(s->qdotvec, r.pertmodind-1);
  // qdot_t += pow(2*r.pertEval/m, 0.5);
  qdot_t = pow(2*r.pertEval/m, 0.5);
  gsl_vector_set(s->qdotvec, r.pertmodind-1, qdot_t);


  /**************************/
  // Writing files !!
  /**************************/
  // trajectory file !!
  sprintf(outfname, "modetraj.%d.txt", s->statetyp);
  outfp1 = fopen(outfname, "w");
  fprintf(outfp1, "%d\n", s->step);
  double systeng = 0.0;
  for (i=0; i<s->NmodeIRcou; i++){
    sind      =  (int) gsl_matrix_get(s->IRs_pqrcountmat, i, 0);
    
    // sprintf(outfname, "modedisp.%04d.txt", sind);
    // outfp = fopen(outfname, "w");
    // fprintf(outfp, "%5.5e\n", gsl_vector_get(s->qvec, sind-1));
    // fclose(outfp);

    // sprintf(outfname, "modevelc.%04d.txt", sind);
    // outfp = fopen(outfname, "w");
    // fprintf(outfp, "%5.5e\n", gsl_vector_get(s->qdotvec, sind-1));
    // fclose(outfp);

    // sprintf(outfname, "modeforc.%04d.txt", sind);
    // outfp = fopen(outfname, "w");
    // f = ForceSys(*s, r, i);
    // gsl_vector_set(s->fvec, sind-1, f.f3);          // storing just the nonlinear part
    // fprintf(outfp, "%5.5e\n", gsl_vector_get(s->fvec, sind-1));
    // fclose(outfp);

    q_t        = gsl_vector_get(s->qvec, sind-1);
    qdot_t     = gsl_vector_get(s->qdotvec, sind-1);
    f        = ForceSys(*s, r, i);
    gsl_vector_set (s->fvec, sind-1, f.f1+f.f3);                                // writing just the nonlinear part
    f_t      = gsl_vector_get (s->fvec, sind-1); ///q_t;
    // f_t      = pow(-gsl_vector_get (s->fvec, sind-1)/q_t, 0.5)/(2*PI);
    gsl_vector_set (s->enonvec, sind-1, f.ep4);
    m        = gsl_vector_get(s->mvec, sind-1);
    fr       = gsl_vector_get(s->frvec, sind-1);
    teng     = 0.5*m*qdot_t*qdot_t + 0.5*m*2*PI*fr*2*PI*fr*q_t*q_t + m*gsl_vector_get(s->enonvec, sind-1);
    systeng += teng;

    // sprintf(outfname, "modeteng.%04d.txt", sind);
    // outfp = fopen(outfname, "w");
    // fprintf(outfp, "%5.5e\n", teng);
    // fclose(outfp);
    fprintf(outfp1, "%d %5.5e %5.5e %5.5e %5.5e\n", sind, q_t, qdot_t, f_t, teng);
  }
  fclose(outfp1);
  sprintf(outfname, "modetotteng.%d.txt", s->statetyp);
  outfp = fopen(outfname, "w");
  fprintf(outfp, "%5.5e\n", systeng);
  fclose(outfp);

  return;
}
