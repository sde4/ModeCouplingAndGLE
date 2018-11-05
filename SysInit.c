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
  double m, gam, fri, Ath;
  double sig;
  double N, Nj, R, ilamhx, ilamhy, ilamh, lamh;
  double frj, frk, frl;
  double gamavg, DissViscElas, DissOscillator;
  FILE *outfp;
  char outfname[40];

  cou1 = 0; cou2 = 0; cou3 = 0; cou4 = 0; cou5 = 0;
  for (k=1; k<=sc.max; k++) {
    for (j=1; j<=k; j++) {
      for (i=1; i<=j; i++) {    // neglecting symmetric pairs i.e. only one of (m,n) and (n,m) will be chosen assuming the membrane is square 
        if (j==k || i==k) {
          mord = i;
          nord = j;
          gsl_matrix_set(s->modindmat, cou1, 0, mord);
          gsl_matrix_set(s->modindmat, cou1, 1, nord);
          gsl_vector_set(s->modindvec, cou1, cou1 + 1);

          // Mode grouping based on symmetry!!
          // 4 groups: SS, SA, AS, AA
          if ((mord%2 == 1) && (nord%2 == 1)) {
            gsl_vector_set(s->SSmodindvec, cou2, cou1 + 1);
            cou2++;
          }
          if ((mord%2 == 1) && (nord%2 == 0)) {
            gsl_vector_set(s->SAmodindvec, cou3, cou1 + 1);
            cou3++;
          }
          if ((mord%2 == 0) && (nord%2 == 1)) {
            gsl_vector_set(s->ASmodindvec, cou4, cou1 + 1);
            cou4++;
          }
          if ((mord%2 == 0) && (nord%2 == 0)) {
            gsl_vector_set(s->AAmodindvec, cou5, cou1 + 1);
            cou5++;
          }
          cou1++;
        }
      }
    }
  }
  // for (k=1; k<=sc.max; k++) {
  //   for (j=1; j<=k; j++) {
  //     for (i=1; i<=j; i++) {    // neglecting symmetric pairs i.e. only one of (m,n) and (n,m) will be chosen assuming the membrane is square 
  //       if (j==k || i==k) {
  //         mord = 25*i;
  //         nord = 25*j;
  //         gsl_matrix_set(s->modindmat, cou1, 0, mord);
  //         gsl_matrix_set(s->modindmat, cou1, 1, nord);
  //         gsl_vector_set(s->modindvec, cou1, cou1 + 1);

  //         // Mode grouping based on symmetry!!
  //         // 4 groups: SS, SA, AS, AA
  //         if ((mord%2 == 1) && (nord%2 == 1)) {
  //           gsl_vector_set(s->SSmodindvec, cou2, cou1 + 1);
  //           cou2++;
  //         }
  //         if ((mord%2 == 1) && (nord%2 == 0)) {
  //           gsl_vector_set(s->SAmodindvec, cou3, cou1 + 1);
  //           cou3++;
  //         }
  //         if ((mord%2 == 0) && (nord%2 == 1)) {
  //           gsl_vector_set(s->ASmodindvec, cou4, cou1 + 1);
  //           cou4++;
  //         }
  //         if ((mord%2 == 0) && (nord%2 == 0)) {
  //           gsl_vector_set(s->AAmodindvec, cou5, cou1 + 1);
  //           cou5++;
  //         }
  //         cou1++;
  //       }
  //     }
  //   }
  // }
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
  // mass calculation !!
  /**************************/
  m = mc.rho * sc.Lx * sc.Ly/4;   // meff=m/4 https://arxiv.org/pdf/1305.0557.pdf

  /**************************/
  // frequency calculation !!
  /**************************/
  for (i = 0; i < r.nmodes; i++) {
    mord    = gsl_matrix_get(s->modindmat, i, 0);
    nord    = gsl_matrix_get(s->modindmat, i, 1);
    ilamhx  = mord/sc.Lx;
    ilamhy  = nord/sc.Ly;
    ilamh   = pow((ilamhx*ilamhx + ilamhy*ilamhy)/2.0, 0.5);        // 2/lamh^2 = 1/lamx^2 + 1/lamy^2;
    lamh    = 1.0/ilamh;
    R       = 2*PI*PI*mc.kb/(mc.Et*sv.e_pre*lamh*lamh);

    fri     = 1/(2*lamh)*pow(2*mc.Et*sv.e_pre*(1+R)/(mc.rho), 0.5);
    gsl_vector_set(s->frvec, i, fri);

    Ath     = pow(kB * sv.T / (m * pow(2 * PI * fri, 2.0)), 0.5);
    gsl_vector_set(s->Athvec, i, Ath);
  }
  // File write !!
  // frequency file !!
  outfp = fopen("frvec.dat", "w");
  for (i = 0; i < r.nmodes - 1; i++) {
    fprintf(outfp, "%10.10f\n", gsl_vector_get(s->frvec, i));
  }
  fprintf(outfp, "%10.10f", gsl_vector_get(s->frvec, i));
  fclose(outfp);

  // Thermal amplitude file !!
  outfp = fopen("Athvec.dat", "w");
  for (i = 0; i < r.nmodes - 1; i++) {
    fprintf(outfp, "%5.5e\n", gsl_vector_get(s->Athvec, i));
  }
  fprintf(outfp, "%5.5e", gsl_vector_get(s->Athvec, i));
  fclose(outfp);


  /**************************/
  // friction calculation !!
  /**************************/
  gamavg = 0;
  if (mc.gam == -1.0){
    for (j=1; j<=round(pow(r.nmodecut, 0.5)); j++) {
      for (i=1; i<=round(pow(r.nmodecut, 0.5)); i++) { 
        mord    = i;
        nord    = j;
        ilamhx  = mord/sc.Lx;
        ilamhy  = nord/sc.Ly;
        ilamh   = pow((ilamhx*ilamhx + ilamhy*ilamhy)/2.0, 0.5);        // 2/lamh^2 = 1/lamx^2 + 1/lamy^2;
        lamh    = 1.0/ilamh;
        R       = 2*PI*PI*mc.kb/(mc.Et*sv.e_pre*lamh*lamh);

        fri     = 1/(2*lamh)*pow(2*mc.Et*sv.e_pre*(1+R)/(mc.rho), 0.5);
        Ath     = pow(kB * sv.T / (m * pow(2 * PI * fri, 2.0)), 0.5);

        DissViscElas = mc.DEt*pow(PI,5.0)*pow(Ath, 4.0)*( 9*pow(sc.Ly,4.0)*pow(mord,4.0) + 2*pow(sc.Lx,2.0)*pow(sc.Ly,2.0)*pow(mord,2.0)*pow(nord,2.0) + 9*pow(sc.Lx,4.0)*pow(nord,4.0) )/(256*pow(sc.Lx*sc.Ly,3.0)) * ((2*PI*fri*mc.tausig)/(1 + (2*PI*fri*mc.tausig)*(2*PI*fri*mc.tausig)) );
        DissOscillator = 1/2.0*PI*m*(2*PI*fri)*pow(Ath, 2.0);
        gamavg  += DissViscElas/DissOscillator;   		     // intrinsic
      }
    }
    gamavg = gamavg/round(r.nmodecut);
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
    Ath = gsl_vector_get(s->Athvec, i);

    /**************************/
    // mass assignment !!
    /**************************/
    gsl_vector_set(s->mvec, i, m);

    /**************************/
    // friction assignment !!
    /**************************/
    if (mc.gam > 0.0 ){
      DissViscElas = mc.DEt*pow(PI,5.0)*pow(Ath, 4.0)*( 9*pow(sc.Ly,4.0)*pow(mord,4.0) + 2*pow(sc.Lx,2.0)*pow(sc.Ly,2.0)*pow(mord,2.0)*pow(nord,2.0) + 9*pow(sc.Lx,4.0)*pow(nord,4.0) )/(256*pow(sc.Lx*sc.Ly,3.0)) * ((2*PI*fri*mc.tausig)/(1 + (2*PI*fri*mc.tausig)*(2*PI*fri*mc.tausig)) );
      DissOscillator = 1/2.0*PI*m*(2*PI*fri)*pow(Ath, 2.0);
      gam = mc.gam + DissViscElas/DissOscillator;   	// intrinsic + user specified
    }
    else if (mc.gam == -1.0){
      gam = gamavg;   				// intrinsic for device scale obtained from constitutive reln
    }
    else if (mc.gam == -2.0){
      // DissViscElas = mc.DEt*pow(PI,5.0)*pow(Ath, 4.0)*( 9*pow(sc.Ly,4.0)*pow(mord,4.0) + 2*pow(sc.Lx,2.0)*pow(sc.Ly,2.0)*pow(mord,2.0)*pow(nord,2.0) + 9*pow(sc.Lx,4.0)*pow(nord,4.0) )/(256*pow(sc.Lx*sc.Ly,3.0)) * ((2*PI*fri*mc.tausig)/(1 + (2*PI*fri*mc.tausig)*(2*PI*fri*mc.tausig)) );
      // DissOscillator = 1/2.0*PI*m*(2*PI*fri)*pow(Ath, 2.0);
      // gam = DissViscElas/DissOscillator;   				// intrinsic for device scale obtained from constitutive reln
      gam = 1000./pow(10.0, 0.63834*pow(ilamh*1E-4, -0.30952));         // intrinsic for MD scale obtained from fitting MD data parameterized for prestrain = 0.001, fitted parameters time-ns, length -Ang
    }
    else if (mc.gam == 0.0){
      gam = 0.0;   					// zero
    }
    gsl_vector_set(s->gamvec, i, gam);

    /**************************/
    // noise assignment !!
    /**************************/
    sig = pow(2 * kB * sv.T * gam / m, 0.5);
    gsl_vector_set(s->sigvec, i, sig);
  }
  
  // File write !!
  // mass file !!
  outfp = fopen("mvec.txt", "w");
  for (i = 0; i < r.nmodes - 1; i++) {
    fprintf(outfp, "%5.5e\n", gsl_vector_get(s->mvec, i));
  }
  fprintf(outfp, "%5.5e", gsl_vector_get(s->mvec, i));
  fclose(outfp);

  // friction file !!
  sprintf(outfname, "gamvec.%d.txt", s->statetyp);
  outfp = fopen(outfname, "w");
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
