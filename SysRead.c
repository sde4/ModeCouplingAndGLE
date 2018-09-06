/*
 * This function initializes the system of Langevin equation
 */

#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h> 
#include <omp.h>
#include "struct_def.h"

/* Definitions for routines */

void SysRead(sys_var* s, run_param r, mat_const mc, state_var sv, sys_const sc, disc_const dc, gsl_rng* gr) {
  printf("\n# 2. System Initialization ...\n");

  int 	i, j, cou1, cou2;
  int 	idat;
  int   sind;
  double fdat;
  double N, Nj, lamh, fr, fri;
  double gam;
  double m;
  double sig;
  double sigma, q, qdot, teng;
  FILE * infp, * outfp;
  char outfname[40];

  /**************************/
  // Reading files !!
  /**************************/

  printf("#    o Input file required: 1. frvec.dat 2. IRs_pqrcombmatsorted.dat\n");
  printf("#                           3. IRs_pqrcountmat.dat 4. alphamat.dat\n");
  printf("#                           5. mvec.dat\n");
  // check for the required input files

  if ((infp=fopen("frvec.dat", "r")) == NULL) 
  {
     printf("#      ERROR: File 1 not found! \n");
     exit(1);
  }
  if ((infp=fopen("IRs_pqrcombmatsorted.dat", "r")) == NULL) 
  {
     printf("#      ERROR: File 2 not found! \n");
     exit(1);
  }
  if ((infp=fopen("IRs_pqrcountmat.dat", "r")) == NULL) 
  {
     printf("#      ERROR: File 3 not found! \n");
     exit(1);
  }
  if ((infp=fopen("alphamat.dat", "r")) == NULL) 
  {
     printf("#      ERROR: File 4 not found! \n");
     exit(1);
  }
  if ((infp=fopen("mvec.dat", "r")) == NULL) 
  {
     printf("#      ERROR: File 5 not found! \n");
     exit(1);
  }
  
  // frequency file !!
  infp = fopen("frvec.dat", "r");
  for (i = 0; i < r.nmodes; i++) {
    fscanf(infp, "%lf", &fdat);
    gsl_vector_set(s->frvec, i, fdat);
  }
  fclose(infp);

  // IR comb sorted file !!
  infp = fopen("IRs_pqrcombmatsorted.dat", "r");
  cou1 = 0;
  while(!feof(infp))
  {
    idat = fgetc(infp);
    if(idat == '\n')
    {
      cou1++;
    }
  }
  s->IRcou = cou1;
  s->IRs_pqrcombmatsorted    = gsl_matrix_alloc(s->IRcou, 4);
  fclose(infp);
  infp = fopen("IRs_pqrcombmatsorted.dat", "r");
  for (i = 0; i < s->IRcou; i++) {
    fscanf(infp, "%d", &idat);
    gsl_matrix_set(s->IRs_pqrcombmatsorted, i, 0, idat);
    fscanf(infp, "%d", &idat);
    gsl_matrix_set(s->IRs_pqrcombmatsorted, i, 1, idat);
    fscanf(infp, "%d", &idat);
    gsl_matrix_set(s->IRs_pqrcombmatsorted, i, 2, idat);
    fscanf(infp, "%d", &idat);
    gsl_matrix_set(s->IRs_pqrcombmatsorted, i, 3, idat);
  }
  fclose(infp);

  // IR count file !!
  infp = fopen("IRs_pqrcountmat.dat", "r");
  cou1 = 0;
  while(!feof(infp))
  {
    idat = fgetc(infp);
    if(idat == '\n')
    {
      cou1++;
    }
  } 
  fclose(infp);
  infp = fopen("IRs_pqrcountmat.dat", "r");
  s->NmodeIRcou = cou1;
  for (i = 0; i < s->NmodeIRcou; i++) {
    fscanf(infp, "%d", &idat);
    gsl_matrix_set(s->IRs_pqrcountmat, i, 0, idat);
    fscanf(infp, "%d", &idat);
    gsl_matrix_set(s->IRs_pqrcountmat, i, 1, idat);
    fscanf(infp, "%d", &idat);
    gsl_matrix_set(s->IRs_pqrcountmat, i, 2, idat);
  }
  fclose(infp);
  
  // coupling file !!
  s->alphamat  	  	  = gsl_matrix_alloc(s->IRcou, 2);
  infp = fopen("alphamat.dat", "r");
  for (i = 0; i < s->IRcou; i++) {
    for (j = 0; j < 2; j++) {
      fscanf(infp, "%lf", &fdat);
      gsl_matrix_set(s->alphamat, i, j, fdat);
    }
  }
  fclose(infp);
  
  // mass file !!
  infp = fopen("mvec.dat", "r");
  for (i = 0; i < r.nmodes; i++) {
    fscanf(infp, "%lf", &fdat);
    gsl_vector_set(s->mvec, i, fdat);
  }
  fclose(infp);

  /**************************/
  // Rest of the initializations !!
  /**************************/
  for (i = 0; i < r.nmodes; i++) {
    fri = gsl_vector_get(s->frvec, i);
    m = gsl_vector_get(s->mvec, i);

    gam = mc.gam;
    gsl_vector_set(s->gamvec, i, gam);

    sig = pow(2 * kB * sv.T * gam / m, 0.5);
    gsl_vector_set(s->sigvec, i, sig);

    /**************************/
    // displacement and velocity 
    // initialization !!
    /**************************/
    sigma = pow(kB * sv.T / (m * pow(2 * PI * fri, 2.0)), 0.5);
    q         = gsl_ran_gaussian(gr, sigma);
    // q = 0.0;
    // q = pow(kB * sv.T / (m * pow(2 * PI * fri, 2.0)), 0.5); 		// 1/2 total energy is stored as KE
    gsl_vector_set(s->qvec, i, q);

    sigma = pow(kB * sv.T / m, 0.5);
    qdot      = gsl_ran_gaussian(gr, sigma);
    // qdot = pow(kB * sv.T / m, 0.5); 		// 1/2 total energy is stored as KE
    gsl_vector_set(s->qdotvec, i, qdot);
  }
  // Initial perturbation of 20 KBT to mode 1
  // fri = gsl_vector_get(s->frvec, 0);
  // m = gsl_vector_get(s->mvec, 0);
  gsl_vector_set(s->gamvec, 0, 0.0);
  gsl_vector_set(s->sigvec, 0, 0.0);
  // gsl_vector_set(s->qvec, 0, pow(10 * kB * sv.T / (m * pow(2 * PI * fri, 2.0)), 0.5));
  // gsl_vector_set(s->qdotvec, 0, pow(10 * kB * sv.T / m, 0.5));


  /**************************/
  // Writing files !!
  /**************************/
  
  // friction file !!
  outfp = fopen("gamvec.txt", "w");
  for (i = 0; i < r.nmodes - 1; i++) {
    fprintf(outfp, "%5.5e\n", gsl_vector_get(s->gamvec, i));
  }
  fprintf(outfp, "%5.5e", gsl_vector_get(s->gamvec, i));
  fclose(outfp);

  // noise file !!
  outfp = fopen("noisevec.txt", "w");
  for (i = 0; i < r.nmodes - 1; i++) {
    fprintf(outfp, "%5.5e\n", gsl_vector_get(s->sigvec, i));
  }
  fprintf(outfp, "%5.5e", gsl_vector_get(s->sigvec, i));
  fclose(outfp);

  // displacement and velocity files !!
  double systeng = 0.0;
  for (i=0; i<s->NmodeIRcou; i++){
    sind      =  (int) gsl_matrix_get(s->IRs_pqrcountmat, i, 0);
    
    sprintf(outfname, "modedisp.%04d.txt", sind);
    outfp = fopen(outfname, "a");
    fprintf(outfp, "%5.5e\n", gsl_vector_get(s->qvec, sind-1));
    fclose(outfp);

    sprintf(outfname, "modevelc.%04d.txt", sind);
    outfp = fopen(outfname, "a");
    fprintf(outfp, "%5.5e\n", gsl_vector_get(s->qdotvec, sind-1));
    fclose(outfp);

    m = gsl_vector_get(s->mvec, sind-1);
    fr = gsl_vector_get(s->frvec, sind-1);
    teng = 0.5*m*gsl_vector_get(s->qdotvec, sind-1)*gsl_vector_get(s->qdotvec, sind-1) +
      0.5*m*pow(2*PI*fr, 2.0)*gsl_vector_get(s->qvec, sind-1)*gsl_vector_get(s->qvec, sind-1);
    systeng += teng;

    sprintf(outfname, "modeteng.%04d.txt", sind);
    outfp = fopen(outfname, "a");
    fprintf(outfp, "%5.5e\n", teng);
    fclose(outfp);
  }
  sprintf(outfname, "modetotteng.txt");
  outfp = fopen(outfname, "a");
  fprintf(outfp, "%5.5e\n", systeng);
  fclose(outfp);

  return;
}
