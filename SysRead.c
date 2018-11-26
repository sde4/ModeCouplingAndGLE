/*
 * This function initializes the system of Langevin equation
 */

#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h> 
#include <omp.h>
#include "struct_def.h"

void SysRead(sys_var* s, run_param r, mat_const mc, state_var sv, sys_const sc, disc_const dc) {
  printf("\n# 2. System Initialization ...\n");

  int i, j, k, l, cou1, cou2, cou3, cou4, cou5;
  int mord, nord, pord, qord;
  int idat;
  int sind;
  double m, gam, fri, Ath;
  double sig;
  double N, Nj, R, ilamhx, ilamhy, ilamh, lamh;
  double gamavg, DissViscElas, DissOscillator;
  double fdat;
  for_var f;
  FILE *infp, *outfp;
  char outfname[40];

  /**************************/
  // Reading files !!
  /**************************/

  printf("#    o Input file required: 1. frvec.dat 2. IRs_pqrcombmatsorted.dat\n");
  printf("#                           3. IRs_pqrcountmat.dat 4. alphamatsorted.dat\n");
  printf("#                           5. Athvec.dat 6. frq2vec.dat\n");
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
  if ((infp=fopen("alphamatsorted.dat", "r")) == NULL) 
  {
     printf("#      ERROR: File 4 not found! \n");
     exit(1);
  }
  if ((infp=fopen("Athvec.dat", "r")) == NULL) 
  {
     printf("#      ERROR: File 5 not found! \n");
     exit(1);
  }
  if ((infp=fopen("frq2vec.dat", "r")) == NULL) 
  {
     printf("#      ERROR: File 5 not found! \n");
     exit(1);
  }
 
  // Count no. of modes in each symmetry group !! 
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
  //for (k=1; k<=sc.max; k++) {
  //  for (j=1; j<=k; j++) {
  //    for (i=1; i<=j; i++) {    // neglecting symmetric pairs i.e. only one of (m,n) and (n,m) will be chosen assuming the membrane is square 
  //      if (j==k || i==k) {
  //        mord = 25*i;
  //        nord = 25*j;
  //        gsl_matrix_set(s->modindmat, cou1, 0, mord);
  //        gsl_matrix_set(s->modindmat, cou1, 1, nord);
  //        gsl_vector_set(s->modindvec, cou1, cou1 + 1);

  //        // Mode grouping based on symmetry!!
  //        // 4 groups: SS, SA, AS, AA
  //        if ((mord%2 == 1) && (nord%2 == 1)) {
  //          gsl_vector_set(s->SSmodindvec, cou2, cou1 + 1);
  //          cou2++;
  //        }
  //        if ((mord%2 == 1) && (nord%2 == 0)) {
  //          gsl_vector_set(s->SAmodindvec, cou3, cou1 + 1);
  //          cou3++;
  //        }
  //        if ((mord%2 == 0) && (nord%2 == 1)) {
  //          gsl_vector_set(s->ASmodindvec, cou4, cou1 + 1);
  //          cou4++;
  //        }
  //        if ((mord%2 == 0) && (nord%2 == 0)) {
  //          gsl_vector_set(s->AAmodindvec, cou5, cou1 + 1);
  //          cou5++;
  //        }
  //        cou1++;
  //      }
  //    }
  //  }
  //}
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

  infp = fopen("frq2vec.dat", "r");
  for (i = 0; i < r.nmodes; i++) {
    fscanf(infp, "%lf", &fdat);
    gsl_vector_set(s->frq2vec, i, fdat);
  }
  fclose(infp);

  // Athvec file !!
  infp = fopen("Athvec.dat", "r");
  for (i = 0; i < r.nmodes; i++) {
    fscanf(infp, "%lf", &fdat);
    gsl_vector_set(s->Athvec, i, fdat);
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
  s->alphamatsorted  	  	  = gsl_matrix_alloc(s->IRcou, 2);
  infp = fopen("alphamatsorted.dat", "r");
  for (i = 0; i < s->IRcou; i++) {
    for (j = 0; j < 2; j++) {
      fscanf(infp, "%lf", &fdat);
      gsl_matrix_set(s->alphamatsorted, i, j, fdat);
    }
  }
  fclose(infp);
  
  /**************************/
  // mass calculation !!
  /**************************/
  m = mc.rho * sc.Lx * sc.Ly/4;   // meff=m/4 https://arxiv.org/pdf/1305.0557.pdf

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
  if (mc.gam == -2.0){
    for (j=1; j<=round(pow(r.nmodecut, 0.5)); j++) {
      for (i=1; i<=round(pow(r.nmodecut, 0.5)); i++) { 
        mord    = i;
        nord    = j;
        ilamhx  = mord/sc.Lx;
        ilamhy  = nord/sc.Ly;
        ilamh   = pow((ilamhx*ilamhx + ilamhy*ilamhy)/2.0, 0.5);        // 2/lamh^2 = 1/lamx^2 + 1/lamy^2;
        lamh    = 1.0/ilamh;

        gam     = 1000.*1E2/pow(10.0, 0.63834*pow(ilamh*1E-6, -0.30952))/2.0;       // intrinsic for MD scale obtained from fitting MD data parameterized for prestrain = 0.001, fitted parameters time-ns, length -Ang, Remember tau_d = 2*tau_e thats why factor 2 in the beginning
        gamavg  += gam;   		     // intrinsic
      }
    }
    gamavg = gamavg/round(r.nmodecut);
  }

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
      // gam = 1000.*1E2/pow(10.0, 0.63834*pow(ilamh*1E-6, -0.30952))/2.0;       // intrinsic for MD scale obtained from fitting MD data parameterized for prestrain = 0.001, fitted parameters time-ns, length -Ang, Remember tau_d = 2*tau_e thats why factor 2 in the beginning
      gam = gamavg;
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
  // Setting gam 0 for initially perturbed mode
  if (r.pertEval != 0.0 && (mc.gam==-1.0 ||mc.gam==-2.0) ){
    m 	 = gsl_vector_get(s->mvec, r.pertmodind-1);
    fri  = gsl_vector_get(s->frvec, r.pertmodind-1);
    Ath  = pow(2*r.pertEval/(m * pow(2 * PI * fri, 2.0)), 0.5);
    gsl_vector_set(s->Athvec, r.pertmodind-1, Ath);

    mord = gsl_matrix_get(s->modindmat, r.pertmodind-1, 0);
    nord = gsl_matrix_get(s->modindmat, r.pertmodind-1, 1);

    DissViscElas = mc.DEt*pow(PI,5.0)*pow(Ath, 4.0)*( 9*pow(sc.Ly,4.0)*pow(mord,4.0) + 2*pow(sc.Lx,2.0)*pow(sc.Ly,2.0)*pow(mord,2.0)*pow(nord,2.0) + 9*pow(sc.Lx,4.0)*pow(nord,4.0) )/(256*pow(sc.Lx*sc.Ly,3.0)) * ((2*PI*fri*mc.tausig)/(1 + (2*PI*fri*mc.tausig)*(2*PI*fri*mc.tausig)) );
    DissOscillator = 1/2.0*PI*m*(2*PI*fri)*pow(Ath, 2.0);
    
    // gam = DissViscElas/DissOscillator;   		// intrinsic
    // sig = pow(2 * kB * sv.T * gam / m, 0.5);
    gam = 0.0;
    sig = 0.0;
    gsl_vector_set(s->gamvec, r.pertmodind-1, gam);
    gsl_vector_set(s->sigvec, r.pertmodind-1, sig);
  }

  /**************************/
  // Writing files !!
  /**************************/
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
