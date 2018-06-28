/*
 * This function performs integration of Langevin equation
 */

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <omp.h>

#include "struct_def.h"

void	SysInit(sys_var *s, run_param r, mat_const mc, state_var sv, sys_const sc)
{
  printf("\n# 2. System Initialization ...\n");

  int i, j, k, cou1, cou2;
  cou1=0;
  for (i=0; i<sc.nmax; i++)
  {
    for (j=0; j<sc.mmax; j++)
    {
      gsl_matrix_set(s->modindmat, cou1, 0, i+1);
      gsl_matrix_set(s->modindmat, cou1, 1, j+1);

      cou1++;
    }
  }

  int     mord, nord;
  double  N, lamh, fr;
  double  gam;
  double  m;
  double  sig;
  double  sigma, q, qdot, teng;
  FILE    *outfp;
  char    outfname[40];

  const gsl_rng_type * gT;
  gsl_rng * gr;
  gsl_rng_env_setup();
  gT    = gsl_rng_default;
  gr    = gsl_rng_alloc(gT);

  for (i=0; i<r.nmodes; i++)
  {
    mord      = gsl_matrix_get(s->modindmat, i, 0);
    nord      = gsl_matrix_get(s->modindmat, i, 1);

    /**************************/
    // frequency calculation !!
    /**************************/
    N         = pow( (mord*mord+nord*nord)/2.0, 0.5);
    lamh      = sc.L/N;
    fr        = 1/lamh*pow(mc.Et*sv.e_pre/(2*mc.rho), 0.5); 
    gsl_vector_set(s->frvec, i, fr);


    /**************************/
    // mass calculation !!
    /**************************/    
    m         = mc.rho*sc.L*sc.L;
    gsl_vector_set(s->mvec, i, m);
    


    /**************************/    
    // coupling term calculation !!
    /**************************/    
    for (j=0; j<r.nmodes; j++)
    {
      gsl_matrix_set(s->alphamat, i, j, 0);
    }



    /**************************/
    // friction calculation !!
    /**************************/    
    gam       = pow(PI, 4.0)*kB*sv.T*mc.DEt*( 9.0*pow(mord, 4.0) + 9.0*pow(nord, 4.0) + 2.0*pow(mord, 2.0)*pow(nord, 2.0))/(128.0*pow(sc.L, 2.0)*pow(m, 2.0)*pow(2*PI*fr, 3.0));
    gsl_vector_set(s->gamvec, i, gam);




    /**************************/
    // noise calculation !!
    /**************************/
    sig       = pow(2*kB*sv.T*gam/m, 0.5);
    gsl_vector_set(s->sigvec, i, sig);



    /**************************/
    // displacement and velocity initialization !!
    /**************************/
    sigma     = pow(kB*sv.T/(m*pow(2*PI*fr, 2.0)), 0.5);
    q         = gsl_ran_gaussian(gr, sigma);
    gsl_vector_set(s->qvec, i, q);

    sigma     = pow(kB*sv.T/m, 0.5);
    qdot      = gsl_ran_gaussian(gr, sigma);
    gsl_vector_set(s->qdotvec, i, qdot);


    /**************************/
    // Writing in a file !!
    /**************************/
    sprintf(outfname, "modedisp.%04d.txt", i);
    outfp = fopen(outfname, "a");
    fprintf(outfp, "%5.5e\n", gsl_vector_get(s->qvec, i));
    fclose(outfp);

    sprintf(outfname, "modevelc.%04d.txt", i);
    outfp = fopen(outfname, "a");
    fprintf(outfp, "%5.5e\n", gsl_vector_get(s->qdotvec, i));
    fclose(outfp);

    teng      = 0.5*m*gsl_vector_get(s->qdotvec, i)*gsl_vector_get(s->qdotvec, i) +
    0.5*m*pow(2*PI*fr, 2.0)*gsl_vector_get(s->qvec, i)*gsl_vector_get(s->qvec, i);

    sprintf(outfname, "modeteng.%04d.txt", i);
    outfp = fopen(outfname, "a");
    fprintf(outfp, "%5.5e\n", teng);
    fclose(outfp);


  }
  gsl_rng_free(gr);
  return;
}
