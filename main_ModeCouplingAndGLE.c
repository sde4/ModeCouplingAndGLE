#include "struct_def.h"

/* Definitions for routines */
void SysInit(sys_var *s, run_param r, mat_const mc, state_var sv, sys_const sc);
void IntegrateSys(sys_var *s, run_param r);

int 	main(int argc, char* argv[])
{
  printf("\n \n");
  printf("##################################################################### \n");
  printf("#                @Subhadeep De  <sde4@illinois.edu>                 # \n"); 
  printf("##################################################################### \n\n");
  printf("# Description: This program implements GLE and modecoupling for squar \n");
  printf("# e grahene resonators \n\n");

  //printf("# Phonon analyzer started ...\n\n");
  
  // if (argc!=4) 
  // {
  //  printf("ERROR: specify three arguments: 1. nmodes 2. run_id \n");
  //  exit(1);
  // }

  printf("# 1. Main started ...\n");
  run_param r;
  sys_var s;
  mat_const mc;
  state_var sv;
  sys_const sc;

  // system constants !!

/*  int nmodes    = atoi(argv[1]);
  double sides  = atoi(argv[2]);
  int run_id    = atoi(argv[3]);
*/  

  sc.mmax             = 3;
  sc.nmax             = 3;
  sc.L                = 5E4;
  sc.run_id           = 1;

  // Material constants !!
  mc.Et               = 340*Npm_eVpA2;                        // eV/A^2
  mc.DEt              = 40*Npm_eVpA2;                         // eV/A^2
  mc.rho              = 7.4E-7*kgpm2_amupA2*amuA2pns2_eV;     // eV/(A^4/ns^2)

  // State variables !!
  sv.T                = 300;                                  // K
  sv.e_pre            = 1E-4;

  // Run Parameters !!
  r.dt                = 0.001;                                // ns
  r.runtime           = 2E6;                                  // ns
  r.nfreq             = 100;
  r.nmodes            = sc.mmax*sc.nmax;   

  printf("# No. of modes: %d  Run id: %d \n", r.nmodes, sc.run_id);

  /* System initialization */
  s.modindmat         = gsl_matrix_alloc(r.nmodes, 2);
  s.frvec             = gsl_vector_alloc(r.nmodes);
  s.gamvec            = gsl_vector_alloc(r.nmodes);
  s.mvec              = gsl_vector_alloc(r.nmodes);
  s.alphamat          = gsl_matrix_alloc(r.nmodes, r.nmodes);
  s.qvec              = gsl_vector_alloc(r.nmodes);
  s.qdotvec           = gsl_vector_alloc(r.nmodes);

  s.sigvec            = gsl_vector_alloc(r.nmodes);
 
  // System Initialization
  SysInit(&s, r, mc, sv, sc);

  // Integrate
  int     i;
  r.nsteps            = (int) r.runtime/r.dt;

  printf("\n# 2. Integration started ...\n");
  for (i=1; i<=r.nsteps; i++)
  {
    r.step = i;

    IntegrateSys(&s, r);

  }
}

