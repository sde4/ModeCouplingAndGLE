/*
 * This function initializes the system of Langevin equation
 */

#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h> 
#include <omp.h>
#include "struct_def.h"

/* Definitions for routines */
void ModeCombinations4NonZeroCouplingConstants(sys_var* , run_param );
void ModeCombinations4InternalResonance(sys_var* );
void ThirdOrderCouplingCalculation(sys_var* , mat_const , disc_const );

void SysInit(sys_var* s, run_param r, mat_const mc, state_var sv, sys_const sc, disc_const dc) {
  printf("\n# 2. System Initialization ...\n");

  int i, j, k, l, cou1, cou2, cou3, cou4, cou5;
  int mord, nord, pord, qord;
  int sind, pind, qind, rind;
  double N, Nj, lamx, lamy, lamh, R, fr;
  double fri, frj, frk, frl;
  double gam;
  double m;
  double sig;
  FILE *outfp;
  char outfname[40];

  cou1 = 0; cou2 = 0; cou3 = 0; cou4 = 0; cou5 = 0;
  for (k=1; k<=sc.max; k++) {
    for (j=1; j<=k; j++) {
      for (i=1; i<=k; i++) {
        if (j==k || i==k) {
          gsl_matrix_set(s->modindmat, cou1, 0, i);
          gsl_matrix_set(s->modindmat, cou1, 1, j);
          gsl_vector_set(s->modindvec, cou1, cou1 + 1);

          // Mode grouping based on symmetry!!
          // 4 groups: SS, SA, AS, AA
          if ((i%2 == 1) && (j%2 == 1)) {
            gsl_vector_set(s->SSmodindvec, cou2, cou1 + 1);
            cou2++;
          }
          if ((i%2 == 1) && (j%2 == 0)) {
            gsl_vector_set(s->SAmodindvec, cou3, cou1 + 1);
            cou3++;
          }
          if ((i%2 == 0) && (j%2 == 1)) {
            gsl_vector_set(s->ASmodindvec, cou4, cou1 + 1);
            cou4++;
          }
          if ((i%2 == 0) && (j%2 == 0)) {
            gsl_vector_set(s->AAmodindvec, cou5, cou1 + 1);
            cou5++;
          }
          cou1++;
        }
      }
    }
  }
  s->SScou = cou2;
  s->SAcou = cou3;
  s->AScou = cou4;
  s->AAcou = cou5;

  // File write !!
  // modindmat file !!
  outfp = fopen("modindmat.dat", "w");
  for (i = 0; i < r.nmodes - 1; i++) {
    fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->modindmat, i, 0));
    fprintf(outfp, "%d\n", (int) gsl_matrix_get(s->modindmat, i, 1));
  }
  fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->modindmat, i, 0));
  fprintf(outfp, "%d\n", (int) gsl_matrix_get(s->modindmat, i, 1));
  fclose(outfp);

  // modindvec file !!
  outfp = fopen("modindvec.dat", "w");
  for (i = 0; i < r.nmodes - 1; i++) {
    fprintf(outfp, "%d\n", (int) gsl_vector_get(s->modindvec, i));
  }
  fprintf(outfp, "%d", (int) gsl_vector_get(s->modindvec, i));
  fclose(outfp);

  // SS modindvec file !!
  outfp = fopen("SSmodindvec.dat", "w");
  for (i = 0; i < s->SScou-1; i++) {
    fprintf(outfp, "%d\n", (int) gsl_vector_get(s->SSmodindvec, i));
  }
  fprintf(outfp, "%d", (int) gsl_vector_get(s->SSmodindvec, i));
  fclose(outfp);

  // SA modindvec file !!
  outfp = fopen("SAmodindvec.dat", "w");
  for (i = 0; i < s->SAcou-1; i++) {
    fprintf(outfp, "%d\n", (int) gsl_vector_get(s->SAmodindvec, i));
  }
  fprintf(outfp, "%d", (int) gsl_vector_get(s->SAmodindvec, i));
  fclose(outfp);

  // AS modindvec file !!
  outfp = fopen("ASmodindvec.dat", "w");
  for (i = 0; i < s->AScou-1; i++) {
    fprintf(outfp, "%d\n", (int) gsl_vector_get(s->ASmodindvec, i));
  }
  fprintf(outfp, "%d", (int) gsl_vector_get(s->ASmodindvec, i));
  fclose(outfp);

  // AA modindvec file !!
  outfp = fopen("AAmodindvec.dat", "w");
  for (i = 0; i < s->AAcou-1; i++) {
    fprintf(outfp, "%d\n", (int) gsl_vector_get(s->AAmodindvec, i));
  }
  fprintf(outfp, "%d", (int) gsl_vector_get(s->AAmodindvec, i));
  fclose(outfp);

  /**************************/
  // frequency calculation !!
  /**************************/
  for (i = 0; i < r.nmodes; i++) {
    mord = gsl_matrix_get(s->modindmat, i, 0);
    nord = gsl_matrix_get(s->modindmat, i, 1);
    lamx = sc.Lx/mord;
    lamy = sc.Ly/nord;
    lamh = pow(2.0/( pow(lamx, -2.0) + pow(lamy, -2.0) ), 0.5);
    R    = 2*PI*PI*mc.kb/(mc.Et*sv.e_pre*lamh*lamh);
    fri  = 1/(2*lamh)*pow(2*mc.Et*sv.e_pre*(1+R)/(mc.rho), 0.5);
    gsl_vector_set(s->frvec, i, fri);
  }

  // File write !!
  // frequency file !!
  outfp = fopen("frvec.dat", "w");
  for (i = 0; i < r.nmodes - 1; i++) {
    fprintf(outfp, "%f\n", gsl_vector_get(s->frvec, i));
  }
  fprintf(outfp, "%f", gsl_vector_get(s->frvec, i));
  fclose(outfp);

  /**************************/
  // Mode combinations leading
  // to nonzero coupling 
  // constants
  /**************************/
  ModeCombinations4NonZeroCouplingConstants(s, r);

  /**************************/
  // Mode combinations leading 
  // to Internal resonances!
  /**************************/
  ModeCombinations4InternalResonance(s);

  /**************************/
  // Third order coupling 
  // constant calculation!
  /**************************/
  ThirdOrderCouplingCalculation(s, mc, dc);

  /**************************/
  // Rest of the initializations !!
  /**************************/
  for (i = 0; i < r.nmodes; i++) {
    fri = gsl_vector_get(s->frvec, i);

    /**************************/
    // mass calculation !!
    /**************************/
    m = mc.rho * sc.Lx * sc.Ly;
    gsl_vector_set(s->mvec, i, m);

    /**************************/
    // friction calculation !!
    /**************************/
    gam = mc.gam; //1E3*pow(PI, 4.0)*kB*sv.T*mc.DEt*( 9.0*pow(mord, 4.0) + 9.0*pow(nord, 4.0) + 2.0*pow(mord, 2.0)*pow(nord, 2.0))/(128.0*pow(sc.L, 2.0)*pow(m, 2.0)*pow(2*PI*fri, 3.0));
    gsl_vector_set(s->gamvec, i, gam);

    /**************************/
    // noise calculation !!
    /**************************/
    sig = pow(2 * kB * sv.T * gam / m, 0.5);
    gsl_vector_set(s->sigvec, i, sig);
  }
  // Initial perturbation of 20 KBT to mode 1
  // gsl_vector_set(s->gamvec, 0, 0.0);
  // gsl_vector_set(s->sigvec, 0, 0.0);
  
  // File write !!
  // mass file !!
  outfp = fopen("mvec.dat", "w");
  for (i = 0; i < r.nmodes - 1; i++) {
    fprintf(outfp, "%f\n", gsl_vector_get(s->mvec, i));
  }
  fprintf(outfp, "%f", gsl_vector_get(s->mvec, i));
  fclose(outfp);


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

  return;
}
