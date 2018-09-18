#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "struct_def.h"

/* Definitions for routines */
void SysInit(sys_var*, run_param , mat_const , state_var , sys_const , disc_const , gsl_rng* );
void SysRead(sys_var*, run_param , mat_const , state_var , sys_const , disc_const , gsl_rng* );
void IntegrateSys(sys_var*, run_param , gsl_rng* , state_var);

int 	main(int argc, char* argv[])
{
  printf("\n \n");
  printf("##################################################################### \n");
  printf("#                @Subhadeep De  <sde4@illinois.edu>                 # \n"); 
  printf("##################################################################### \n\n");
  printf("# Description: This program implements GLE and modecoupling for squar \n");
  printf("# e grahene resonators \n\n");

  printf("# Langevin dynamics simulator for coupled modes started ...\n\n");
  
  if (argc!=9) 
  {
   printf("ERROR: specify four arguments: 1. mmax  2. nmax  3. timestep (ns)  4. temp (K)  5. gam  6. tol  7. Init(0)/Read(1)  8. seed \n");
   exit(1);
  }

  printf("# 1. Main started ...\n");
  run_param r;
  sys_var s;
  mat_const mc;
  state_var sv;
  sys_const sc;
  disc_const dc;
  int seed 	      = atoi(argv[8]); 

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
  mc.gam 	      = atof(argv[5]);			                  // 1/ns
  // mc.alpha            = atof(argv[4]);

  // State variables !!
  sv.T                = atof(argv[4]);                        // K
  sv.e_pre            = 1E-4;

  // Run Parameters !!
  r.dt                = atof(argv[3]);                        // ns
  r.runtime           = 1E6;                                  // ns
  r.nfreq             = 20;
  r.nmodes            = sc.mmax*sc.nmax;   

  // Discretization Constants !!
  dc.Nfx              = 51;                                  // set by mathematica wrapper
  dc.Nfy              = 26;                                  // set by mathematica wrapper
  dc.Lx 	      = sc.Lx/(1E4);			     // um - using different unit of length for discretization
  dc.Ly 	      = sc.Ly/(1E4);			     // um - using different unit of length for discretization
  dc.hx 	      = dc.Lx/(dc.Nfx-1.0);		     // um - using different unit of length for discretization
  dc.hy 	      = dc.Ly/(dc.Nfy-1.0); 	 	     // um - using different unit of length for discretization
  dc.gamunitconv      = 1.0/pow(1E4, 4.0);		     // Unit of gam is L^-4 - um^-4: convert to A^-4


  /* System initialization */
  s.modindmat         = gsl_matrix_alloc(r.nmodes, 2);
  s.modindvec         = gsl_vector_alloc(r.nmodes);
  s.SSmodindvec       = gsl_vector_alloc(r.nmodes);         // check the dimensions
  s.SAmodindvec       = gsl_vector_alloc(r.nmodes);
  s.ASmodindvec       = gsl_vector_alloc(r.nmodes);
  s.AAmodindvec       = gsl_vector_alloc(r.nmodes);
  s.IRs_pqrcountmat   = gsl_matrix_alloc(r.nmodes, 3);
  s.frvec             = gsl_vector_alloc(r.nmodes);
  s.gamvec            = gsl_vector_alloc(r.nmodes);
  s.mvec              = gsl_vector_alloc(r.nmodes);
  s.qvec              = gsl_vector_alloc(r.nmodes);
  s.qdotvec           = gsl_vector_alloc(r.nmodes);
  s.fvec              = gsl_vector_alloc(r.nmodes);
  s.enonvec           = gsl_vector_alloc(r.nmodes);
  s.sigvec            = gsl_vector_alloc(r.nmodes);
  s.tol 	      = atof(argv[6]);			   // detuning parameter for IRs

  const char *mode[2];
  mode[0] = "Init";
  mode[1] = "Read";
  printf("# No. of modes: %d\n# dumpstep: 	  %lf\n# temp: 	  %lf\n# gam: 	  %lf\n# tol: 	  %lf\n# mode: 	  %s\n# seed: 	  %d\n",r.nmodes, r.dt*r.nfreq, sv.T, mc.gam, s.tol, mode[atoi(argv[7])], seed);

  // Random number initialization
  const gsl_rng_type * gT;
  gsl_rng * gr;
  gsl_rng_env_setup();
  gT    = gsl_rng_default;
  gr    = gsl_rng_alloc(gT);

  // System Initialization
  if (atoi(argv[7])==0)
    SysInit(&s, r, mc, sv, sc, dc, gr);
  else 
    SysRead(&s, r, mc, sv, sc, dc, gr);

  // Integrate
  int     i;
  r.nsteps            = (int) r.runtime/r.dt;

  printf("\n# 3. Integration started ...\n");
  for (i=1; i<=r.nsteps; i++)
  {
    r.step = i;

    IntegrateSys(&s, r, gr, sv);

  }
  gsl_rng_free(gr);
}

