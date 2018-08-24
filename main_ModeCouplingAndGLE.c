#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "struct_def.h"

/* Definitions for routines */
void SysInit(sys_var *s, run_param r, mat_const mc, state_var sv, sys_const sc, gsl_rng * gr);
void IntegrateSys(sys_var *s, run_param r, gsl_rng * gr);

int 	main(int argc, char* argv[])
{
  printf("\n \n");
  printf("##################################################################### \n");
  printf("#                @Subhadeep De  <sde4@illinois.edu>                 # \n"); 
  printf("##################################################################### \n\n");
  printf("# Description: This program implements GLE and modecoupling for squar \n");
  printf("# e grahene resonators \n\n");

  printf("# Langevin dynamics simulator for coupled modes started ...\n\n");
  
  if (argc!=6) 
  {
   printf("ERROR: specify four arguments: 1. mmax 2. nmax 3. gam 4. alpha 5. seed \n");
   exit(1);
  }

  printf("# 1. Main started ...\n");
  run_param r;
  sys_var s;
  mat_const mc;
  state_var sv;
  sys_const sc;
  int seed 	      = atoi(argv[5]); 

  // system constants !!

/*  int nmodes    = atoi(argv[1]);
  double sides  = atoi(argv[2]);
  int run_id    = atoi(argv[3]);
*/  

  sc.mmax             = atoi(argv[1]);
  sc.nmax             = atoi(argv[2]);
  sc.Lx               = 2E4;
  sc.Ly               = 1E4;
  sc.run_id           = 1;

  // Material constants !!
  mc.Et               = 340*Npm_eVpA2;                        // eV/A^2
  mc.DEt              = 40*Npm_eVpA2;                         // eV/A^2
  mc.rho              = 7.4E-7*kgpm2_amupA2*amuA2pns2_eV;     // eV/(A^4/ns^2)
  mc.gam 	            = atof(argv[3]);			                  // 1/ns
  mc.alpha            = atof(argv[4]);

  // State variables !!
  sv.T                = 300;                                  // K
  sv.e_pre            = 1E-4;

  // Run Parameters !!
  r.dt                = 0.02;                                 // ns
  r.runtime           = 2E6;                                  // ns
  r.nfreq             = 50;
  r.nmodes            = sc.mmax*sc.nmax;   

  printf("# No. of modes: %d  gam: %f alpha: %f seed: %d \n", r.nmodes, mc.gam, mc.alpha, seed);

  /* System initialization */
  s.modindmat         = gsl_matrix_alloc(r.nmodes, 2);
  s.modindvec         = gsl_vector_alloc(r.nmodes);
  s.SSmodindvec       = gsl_vector_alloc(r.nmodes/4);         // check the dimensions
  s.SAmodindvec       = gsl_vector_alloc(r.nmodes/4);
  s.ASmodindvec       = gsl_vector_alloc(r.nmodes/4);
  s.AAmodindvec       = gsl_vector_alloc(r.nmodes/4);
  s.NZs_pqrcombmat    = gsl_matrix_alloc(r.nmodes/4*r.nmodes/4*r.nmodes/4*(6+6*3)*4, 4);
  s.IRs_pqrcountmat   = gsl_matrix_alloc(r.nmodes, 2);
  s.frvec             = gsl_vector_alloc(r.nmodes);
  s.gamvec            = gsl_vector_alloc(r.nmodes);
  s.mvec              = gsl_vector_alloc(r.nmodes);
  s.alphamat          = gsl_matrix_alloc(r.nmodes, r.nmodes);
  s.qvec              = gsl_vector_alloc(r.nmodes);
  s.qdotvec           = gsl_vector_alloc(r.nmodes);
  s.sigvec            = gsl_vector_alloc(r.nmodes);


  // Random number initialization
  const gsl_rng_type * gT;
  gsl_rng * gr;
  gsl_rng_env_setup();
  gT    = gsl_rng_default;
  gr    = gsl_rng_alloc(gT);

  // System Initialization
  SysInit(&s, r, mc, sv, sc, gr);

  // Integrate
  int     i;
  r.nsteps            = (int) r.runtime/r.dt;

  printf("\n# 3. Integration started ...\n");
  for (i=1; i<=r.nsteps; i++)
  {
    r.step = i;

    IntegrateSys(&s, r, gr);

  }
  gsl_rng_free(gr);
}

