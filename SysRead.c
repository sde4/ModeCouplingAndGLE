/*
 * This function initializes the system of Langevin equation
 */

#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h> 
#include <omp.h>
#include "struct_def.h"

/* Definitions for routines */
for_var ForceSys(sys_var , run_param , int );

void SysRead(sys_var* s, run_param r, mat_const mc, state_var sv, sys_const sc, disc_const dc, gsl_rng* gr) {
  printf("\n# 2. System Initialization ...\n");

  int i, j, k, l, cou1, cou2, cou3, cou4, cou5;
  int idat;
  int sind;
  double fdat;
  double N, Nj, lamh, fr, fri;
  double gam;
  double m;
  double sig;
  double sigma, q_t, qdot_t, f_t, teng;
  for_var f;
  FILE *infp, *outfp, *outfp1;
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
 
  // Count no. of modes in each symmetry group !! 
  cou1 = 0; cou2 = 0; cou3 = 0; cou4 = 0; cou5 = 0;
  for (j = 0; j < sc.nmax; j++) {
    for (i = 0; i < sc.mmax; i++) {
      // Mode grouping based on symmetry!!
      // 4 groups: SS, SA, AS, AA
      if (((i + 1) % 2 == 1) && ((j + 1) % 2 == 1)) {
        cou2++;
      }
      if (((i + 1) % 2 == 1) && ((j + 1) % 2 == 0)) {
        cou3++;
      }
      if (((i + 1) % 2 == 0) && ((j + 1) % 2 == 1)) {
        cou4++;
      }
      if (((i + 1) % 2 == 0) && ((j + 1) % 2 == 0)) {
        cou5++;
      }
      cou1++;
    }
  }
  s->SScou = cou2;
  s->SAcou = cou3;
  s->AScou = cou4;
  s->AAcou = cou5;
  
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
  s->IRs_pqrcombmatsorted    = gsl_matrix_alloc(s->IRcou, 6);
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
    fscanf(infp, "%d", &idat);
    gsl_matrix_set(s->IRs_pqrcombmatsorted, i, 4, idat);
    fscanf(infp, "%d", &idat);
    gsl_matrix_set(s->IRs_pqrcombmatsorted, i, 5, idat);
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
    q_t         = gsl_ran_gaussian(gr, sigma);
    // q_t = 1E-4;
    // q_t = pow(kB * sv.T / (m * pow(2 * PI * fri, 2.0)), 0.5); 		// 1/2 total energy is stored as PE
    gsl_vector_set(s->qvec, i, q_t);

    // sigma = pow(kB * sv.T / m, 0.5);
    // qdot_t      = gsl_ran_gaussian(gr, sigma);
    qdot_t = 0.0;
    // qdot_t = pow(kB * sv.T / m, 0.5); 					// 1/2 total energy is stored as KE
    gsl_vector_set(s->qdotvec, i, qdot_t);

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

  // trajectory file !!
  outfp1 = fopen("modetraj.txt", "w");
  fprintf(outfp1, "%d\n", 0);
  double systeng = 0.0;
  for (i=0; i<s->NmodeIRcou; i++){
  //for (i=0; i<4; i++){
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
    f_t      = gsl_vector_get (s->fvec, sind-1)/q_t;
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
  outfp = fopen("modetotteng.txt", "w");
  fprintf(outfp, "%5.5e\n", systeng);
  fclose(outfp);
  fclose(outfp1);

  return;
}
