/*
 * This function initializes the system of Langevin equation
 */

#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h> 
#include <omp.h>
#include "struct_def.h"

/* Definitions for routines */
void ModeCombinations4NonZeroCouplingConstants(sys_var * s, run_param r);
void ModeCombinations4InternalResonance(sys_var * s, run_param r, double tol);

void SysInit(sys_var * s, run_param r, mat_const mc, state_var sv, sys_const sc, gsl_rng * gr) {
  printf("\n# 2. System Initialization ...\n");

  int i, j, k, l, cou1, cou2, cou3, cou4, cou5;
  cou1 = 0;
  cou2 = 0;
  cou3 = 0;
  cou4 = 0;
  cou5 = 0;
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

  int mord, nord, pord, qord;
  int sind, pind, qind, rind;
  double N, Nj, lamh, fr;
  double fri, frj, frk, frl;
  double tol = 5E-4; // detuning parameter for IRs
  double gam;
  double m;
  double alpha;
  double sig;
  double sigma, q, qdot, teng;
  FILE * outfp;
  char outfname[40];

  for (i = 0; i < r.nmodes; i++) {
    mord = gsl_matrix_get(s->modindmat, i, 0);
    nord = gsl_matrix_get(s->modindmat, i, 1);

    /**************************/
    // frequency calculation !!
    /**************************/
    fri = 1 / (2 * PI) * pow(pow(mord * PI / sc.Lx, 2.0) + pow(nord * PI / sc.Ly, 2.0), 0.5) * pow(2 * mc.Et * sv.e_pre / (mc.rho), 0.5);
    gsl_vector_set(s->frvec, i, fri);
  }

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
  ModeCombinations4InternalResonance(s, r, tol);

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
    // coupling term calculation !!
    /**************************/
    for (j = 0; j < r.nmodes; j++) {
      pord = gsl_matrix_get(s->modindmat, j, 0);
      qord = gsl_matrix_get(s->modindmat, j, 1);
      alpha = mc.alpha; //2.0*pow(PI, 4.0)*mc.Et*(mord*mord*pord*pord + nord*nord*qord*qord)/(64*sc.L*sc.L);
      gsl_matrix_set(s->alphamat, i, j, alpha);
    }

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
    sigma = pow(kB * sv.T / (m * pow(2 * PI * fri, 2.0)), 0.5);
    // q         = gsl_ran_gaussian(gr, sigma);
    q = 0.0;
    gsl_vector_set(s->qvec, i, q);

    sigma = pow(kB * sv.T / m, 0.5);
    // qdot      = gsl_ran_gaussian(gr, sigma);
    qdot = pow(kB * sv.T / m, 0.5); // 1/2 total energy is stored as KE
    gsl_vector_set(s->qdotvec, i, qdot);
  }
  // Initial perturbation of 20 KBT to mode 1
  gsl_vector_set(s->gamvec, 0, 0.0);
  gsl_vector_set(s->sigvec, 0, 0.0);
  gsl_vector_set(s->qvec, 0, 0.0);
  gsl_vector_set(s->qdotvec, 0, pow(20 * kB * sv.T / m, 0.5));

  /**************************/
  // Writing in a file !!
  /**************************/

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
  for (i = 0; i < r.nmodes / 4 - 1; i++) {
    fprintf(outfp, "%d\n", (int) gsl_vector_get(s->SSmodindvec, i));
  }
  fprintf(outfp, "%d", (int) gsl_vector_get(s->SSmodindvec, i));
  fclose(outfp);

  // SA modindvec file !!
  outfp = fopen("SAmodindvec.dat", "w");
  for (i = 0; i < r.nmodes / 4 - 1; i++) {
    fprintf(outfp, "%d\n", (int) gsl_vector_get(s->SAmodindvec, i));
  }
  fprintf(outfp, "%d", (int) gsl_vector_get(s->SAmodindvec, i));
  fclose(outfp);

  // AS modindvec file !!
  outfp = fopen("ASmodindvec.dat", "w");
  for (i = 0; i < r.nmodes / 4 - 1; i++) {
    fprintf(outfp, "%d\n", (int) gsl_vector_get(s->ASmodindvec, i));
  }
  fprintf(outfp, "%d", (int) gsl_vector_get(s->ASmodindvec, i));
  fclose(outfp);

  // AA modindvec file !!
  outfp = fopen("AAmodindvec.dat", "w");
  for (i = 0; i < r.nmodes / 4 - 1; i++) {
    fprintf(outfp, "%d\n", (int) gsl_vector_get(s->AAmodindvec, i));
  }
  fprintf(outfp, "%d", (int) gsl_vector_get(s->AAmodindvec, i));
  fclose(outfp);

  // NZ comb file !!
  outfp = fopen("NZs_pqrcombmat.dat", "w");
  for (i = 0; i < s->NZcou - 1; i++) {
    fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->NZs_pqrcombmat, i, 0));
    fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->NZs_pqrcombmat, i, 1));
    fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->NZs_pqrcombmat, i, 2));
    fprintf(outfp, "%d\n", (int) gsl_matrix_get(s->NZs_pqrcombmat, i, 3));
  }
  fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->NZs_pqrcombmat, i, 0));
  fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->NZs_pqrcombmat, i, 1));
  fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->NZs_pqrcombmat, i, 2));
  fprintf(outfp, "%d\n", (int) gsl_matrix_get(s->NZs_pqrcombmat, i, 3));
  fclose(outfp);

  // frequency file !!
  outfp = fopen("frvec.dat", "w");
  for (i = 0; i < r.nmodes - 1; i++) {
    fprintf(outfp, "%f\n", gsl_vector_get(s->frvec, i));
  }
  fprintf(outfp, "%f", gsl_vector_get(s->frvec, i));
  fclose(outfp);

  // IR comb file !!
  outfp = fopen("IRs_pqrcombmat.dat", "w");
  for (i = 0; i < s->IRcou - 1; i++) {
    fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 0));
    fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 1));
    fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 2));
    fprintf(outfp, "%d\n", (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 3));
  }
  fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 0));
  fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 1));
  fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 2));
  fprintf(outfp, "%d\n", (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 3));
  fclose(outfp);

  // IR comb sorted file !!
  outfp = fopen("IRs_pqrcombmatsorted.dat", "w");
  for (i = 0; i < s->IRcou - 1; i++) {
    fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 0));
    fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 1));
    fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 2));
    fprintf(outfp, "%d\n", (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 3));
  }
  fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 0));
  fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 1));
  fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 2));
  fprintf(outfp, "%d\n", (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 3));
  fclose(outfp);

  // IR count file !!
  outfp = fopen("IRs_pqrcountmat.dat", "w");
  for (i = 0; i < s->NmodeIRcou - 1; i++) {
    fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcountmat, i, 0));
    fprintf(outfp, "%d\n", (int) gsl_matrix_get(s->IRs_pqrcountmat, i, 1));
  }
  fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcountmat, i, 0));
  fprintf(outfp, "%d\n", (int) gsl_matrix_get(s->IRs_pqrcountmat, i, 1));
  fclose(outfp);

  // mass file !!
  outfp = fopen("mvec.dat", "w");
  for (i = 0; i < r.nmodes - 1; i++) {
    fprintf(outfp, "%f\n", gsl_vector_get(s->mvec, i));
  }
  fprintf(outfp, "%f", gsl_vector_get(s->mvec, i));
  fclose(outfp);

  // coupling file !!
  outfp = fopen("coupmat.dat", "w");
  for (i = 0; i < r.nmodes - 1; i++) {
    for (j = 0; j < r.nmodes - 1; j++) {
      fprintf(outfp, "%5.5e\t", gsl_matrix_get(s->alphamat, i, j));
    }
    fprintf(outfp, "%5.5e\n", gsl_matrix_get(s->alphamat, i, j));
  }
  for (j = 0; j < r.nmodes - 1; j++) {
    fprintf(outfp, "%5.5e\t", gsl_matrix_get(s->alphamat, i, j));
  }
  fprintf(outfp, "%5.5e\n", gsl_matrix_get(s->alphamat, i, j));
  fclose(outfp);

  // friction file !!
  outfp = fopen("gamvec.dat", "w");
  for (i = 0; i < r.nmodes - 1; i++) {
    fprintf(outfp, "%5.5e\n", gsl_vector_get(s->gamvec, i));
  }
  fprintf(outfp, "%5.5e", gsl_vector_get(s->gamvec, i));
  fclose(outfp);

  // noise file !!
  outfp = fopen("noisevec.dat", "w");
  for (i = 0; i < r.nmodes - 1; i++) {
    fprintf(outfp, "%5.5e\n", gsl_vector_get(s->sigvec, i));
  }
  fprintf(outfp, "%5.5e", gsl_vector_get(s->sigvec, i));
  fclose(outfp);

  // displacement and velocity files !!
  double systeng = 0.0;
  for (i = 0; i < r.nmodes; i++) {
    sprintf(outfname, "modedisp.%04d.txt", i + 1);
    outfp = fopen(outfname, "a");
    fprintf(outfp, "%5.5e\n", gsl_vector_get(s->qvec, i));
    fclose(outfp);

    sprintf(outfname, "modevelc.%04d.txt", i + 1);
    outfp = fopen(outfname, "a");
    fprintf(outfp, "%5.5e\n", gsl_vector_get(s->qdotvec, i));
    fclose(outfp);

    m = gsl_vector_get(s->mvec, i);
    fr = gsl_vector_get(s->frvec, i);
    teng = 0.5 * m * gsl_vector_get(s->qdotvec, i) * gsl_vector_get(s->qdotvec, i) +
      0.5 * m * pow(2 * PI * fr, 2.0) * gsl_vector_get(s->qvec, i) * gsl_vector_get(s->qvec, i);
    systeng += teng;

    sprintf(outfname, "modeteng.%04d.txt", i + 1);
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
