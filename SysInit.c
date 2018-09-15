/*
 * This function initializes the system of Langevin equation
 */

#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h> 
#include <omp.h>
#include "struct_def.h"

/* Definitions for routines */
void ModeCombinations4NonZeroCouplingConstants(sys_var* , run_param );
void ModeCombinations4InternalResonance(sys_var* , double );
void ThirdOrderCouplingCalculation(sys_var* , mat_const , disc_const );
for_var ForceSys(sys_var , run_param , int );

void SysInit(sys_var* s, run_param r, mat_const mc, state_var sv, sys_const sc, disc_const dc, gsl_rng* gr) {
  printf("\n# 2. System Initialization ...\n");

  int i, j, k, l, cou1, cou2, cou3, cou4, cou5;
  int mord, nord, pord, qord;
  int sind, pind, qind, rind;
  double N, Nj, lamh, fr;
  double fri, frj, frk, frl;
  double tol = 5E-4; // detuning parameter for IRs
  double gam;
  double m;
  double sig;
  double sigma, q_t, qdot_t, f_t, teng;
  for_var f;
  FILE *outfp, *outfp1;
  char outfname[40];

  cou1 = 0; cou2 = 0; cou3 = 0; cou4 = 0; cou5 = 0;
  for (j = 0; j < sc.nmax; j++) {
    for (i = 0; i < sc.mmax; i++) {
      gsl_matrix_set(s->modindmat, cou1, 0, i + 1);
      gsl_matrix_set(s->modindmat, cou1, 1, j + 1);
      gsl_vector_set(s->modindvec, cou1, cou1 + 1);

      // Mode grouping based on symmetry!!
      // 4 groups: SS, SA, AS, AA
      if (((i + 1) % 2 == 1) && ((j + 1) % 2 == 1)) {
        gsl_vector_set(s->SSmodindvec, cou2, cou1 + 1);
        cou2++;
      }
      if (((i + 1) % 2 == 1) && ((j + 1) % 2 == 0)) {
        gsl_vector_set(s->SAmodindvec, cou3, cou1 + 1);
        cou3++;
      }
      if (((i + 1) % 2 == 0) && ((j + 1) % 2 == 1)) {
        gsl_vector_set(s->ASmodindvec, cou4, cou1 + 1);
        cou4++;
      }
      if (((i + 1) % 2 == 0) && ((j + 1) % 2 == 0)) {
        gsl_vector_set(s->AAmodindvec, cou5, cou1 + 1);
        cou5++;
      }
      cou1++;
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
    fri = 1 / (2 * PI) * pow(pow(mord * PI / sc.Lx, 2.0) + pow(nord * PI / sc.Ly, 2.0), 0.5) * pow(2 * mc.Et * sv.e_pre / (mc.rho), 0.5);
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
  ModeCombinations4InternalResonance(s, tol);

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

    /**************************/
    // displacement and velocity initialization !!
    /**************************/
    // sigma = pow(kB * sv.T / (m * pow(2 * PI * fri, 2.0)), 0.5);
    // q_t         = gsl_ran_gaussian(gr, sigma);
    q_t = 1E-4;
    // q_t = pow(kB * sv.T / (m * pow(2 * PI * fri, 2.0)), 0.5); 			// 1/2 total energy is stored as PE
    gsl_vector_set(s->qvec, i, q_t);

    sigma = pow(kB * sv.T / m, 0.5);
    qdot_t      = gsl_ran_gaussian(gr, sigma);
    // qdot_t = pow(kB * sv.T / m, 0.5); 						// 1/2 total energy is stored as KE
    gsl_vector_set(s->qdotvec, i, qdot_t);

  }
  // Initial perturbation of 20 KBT to mode 1
  gsl_vector_set(s->gamvec, 0, 0.0);
  gsl_vector_set(s->sigvec, 0, 0.0);
  // gsl_vector_set(s->qvec, 0, pow(10 * kB * sv.T / (m * pow(2 * PI * fri, 2.0)), 0.5));
  // gsl_vector_set(s->qdotvec, 0, pow(10 * kB * sv.T / m, 0.5));

  
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
