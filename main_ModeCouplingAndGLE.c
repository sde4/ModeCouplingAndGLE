#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "struct_def.h"
#include <time.h>

/* Definitions for routines */
void SysInit(sys_var*, run_param , mat_const , state_var , sys_const , disc_const );
void SysRead(sys_var*, run_param , mat_const , state_var , sys_const , disc_const );
void StateInit(sys_var*, run_param , state_var , gsl_rng* );
void StateRead(sys_var*, run_param );
void IntegrateSys(sys_var*, run_param , gsl_rng* , state_var );
void StateStore(sys_var* ,run_param );

int 	main(int argc, char* argv[])
{
  struct timespec begin, end;
  clock_gettime(CLOCK_MONOTONIC_RAW, &begin);

  printf("\n \n");
  printf("##################################################################### \n");
  printf("#                @Subhadeep De  <sde4@illinois.edu>                 # \n"); 
  printf("##################################################################### \n\n");
  printf("# Description: This program implements GLE and modecoupling for squar \n");
  printf("# e grahene resonators \n\n");

  printf("# Langevin dynamics simulator for coupled modes started ...\n\n");
  
  int gargc = 11;
  if ((argc-1)!=gargc) 
  {
   printf("ERROR: no. of arguments given %d, needed %d\n", (argc-1), gargc);
   printf("ERROR: specify 9 arguments: 1. max 2. timestep (ns)  3. runtime (ns)  4. temp (K)  5. gam  6. tol  7. sys-init(0)/read(1)  8. state-init(0)/read(1)  9. pertmodeind  10. pertEscal  11. seed\n");
   exit(1);
  }

  printf("# 1. Main started ...\n");
  run_param r;
  sys_var s;
  mat_const mc;
  state_var sv;
  sys_const sc;
  disc_const dc;
  int seed 	      = atoi(argv[11]); 
  // system constants !!

/*  int nmodes    = atoi(argv[1]);
  double sides  = atoi(argv[2]);
  int run_id    = atoi(argv[3]);
*/  

  sc.max             = atoi(argv[1]);
  // sc.nmax             = atoi(argv[2]);
  // sc.Lx               = 0.04E4; // A
  // sc.Ly               = 0.04E4; // A
  // sc.Lx               = 10; // um
  // sc.Ly               = 10; // um
  sc.Lx               = 0.005; // 100um
  sc.Ly               = 0.005; // 100um
  sc.run_id           = 1;

  // State variables !!
  sv.T                = atof(argv[4]);                       				// K
  sv.e_pre            = 1E-6;

  // Material constants !!
  mc.kb               = 1.5;                        	      				 	// eV
  // mc.Et               = 340*Npm_eVpA2;                       			 	// eV/A^2
  // mc.DEt              = 40*Npm_eVpA2;                        			 	// eV/A^2
  // mc.rho              = 7.4E-7*kgpm2_amupA2*amuA2pns2_eV;     				// eV/(A^4/ns^2)
  // mc.Et               = 331.0*Npm_eVpum2;                        				// eV/um^2
  // mc.DEt              = 3.09754e+01*Npm_eVpum2;						// eV/um^2
  // mc.tausig           = 4.35758e+02/1.0E3;                         				// ns
  // mc.rho              = 7.4E-7*kgpm2_amupum2*amuum2pns2_eV;    				// eV/(um^4/ns^2)
  mc.Et               = 331.0*Npm_eVp100um2;                        			// eV/100um^2
  mc.DEt              = 3.09754e+01*Npm_eVp100um2;					// eV/100um^2
  mc.tausig           = 4.35758e+02/1.0E3;                         			// ns
  mc.rho              = 7.4E-7*kgpm2_amup100um2*amu100um2pns2_eV;    			// eV/(100um^4/ns^2)
  mc.gam 	      = atof(argv[5]);			                  		// 1/ns
  // mc.alpha            = atof(argv[4]);

  // Run Parameters !!
  double R, lamhcut, nmodecut;
  r.dt                = atof(argv[2]);                        				// ns
  r.runtime           = atof(argv[3]);                                  		// ns
  r.nfreq             = 20;
  r.nmodes            = sc.max*(sc.max+1)/2;						// Change if needed !!!!!!!!!!!!!!!!!!!   
  r.pertmodind        = (int) atoi(argv[9]);   
  r.pertEval          = kB*sv.T*atof(argv[10]);  	      				// eV 
  R		      = 5.0;								// assuming stretching energy << bending energy for R>10
  lamhcut	      = pow(2*PI*PI*mc.kb/(mc.Et*sv.e_pre*R), 0.5);
  nmodecut            = sc.Lx*sc.Ly/(lamhcut*lamhcut);
  r.modefact          = nmodecut/r.nmodes;

  // Discretization Constants !!
  dc.Nfx              = 51;                                  				// set by mathematica wrapper
  dc.Nfy              = 51;                                  				// set by mathematica wrapper
  // dc.Lx 	      = sc.Lx/(1E4);			     				// um - using different unit of length for discretization
  // dc.Ly 	      = sc.Ly/(1E4);			     				// um - using different unit of length for discretization
  // dc.Lx 	      = sc.Lx;			     					// um - using same unit of length for discretization
  // dc.Ly 	      = sc.Ly;			     					// um - using same unit of length for discretization
  // dc.hx 	      = dc.Lx/(dc.Nfx-1.0);		     				// um - using different unit of length for discretization
  // dc.hy 	      = dc.Ly/(dc.Nfy-1.0); 	 	     				// um - using different unit of length for discretization
  dc.Lx 	      = sc.Lx;			     					// 100um - using same unit of length for discretization
  dc.Ly 	      = sc.Ly;			     					// 100um - using same unit of length for discretization
  dc.hx 	      = dc.Lx/(dc.Nfx-1.0);		     				// 100um - using different unit of length for discretization
  dc.hy 	      = dc.Ly/(dc.Nfy-1.0); 	 	     				// 100um - using different unit of length for discretization
  // dc.gamunitconv      = 1.0/pow(1E4, 4.0);		     				// Unit of gam is L^-4 - um^-4: convert to A^-4
  // dc.gamunitconv      = 1.0;		    	     					// Unit of gam is L^-4 - um^-4: convert to um^-4
  dc.gamunitconv      = 1.0;		    	     					// Unit of gam is L^-4 - 100um^-4: convert to 100um^-4


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
  s.Athvec            = gsl_vector_alloc(r.nmodes);
  s.mvec              = gsl_vector_alloc(r.nmodes);
  s.qvec              = gsl_vector_alloc(r.nmodes);
  s.qdotvec           = gsl_vector_alloc(r.nmodes);
  s.fvec              = gsl_vector_alloc(r.nmodes);
  s.enonvec           = gsl_vector_alloc(r.nmodes);
  s.sigvec            = gsl_vector_alloc(r.nmodes);
  s.tol 	      = atof(argv[6]);			   // detuning parameter for IRs
  s.systyp 	      = atoi(argv[7]);
  s.statetyp 	      = atoi(argv[8]);

  const char *sysmode[3], *statemode[3];
  sysmode[0] = "Init";
  sysmode[1] = "Read";
  statemode[0] = "Init";
  statemode[1] = "Read";
  printf("#    o No. of modes:\t%d\n#    o Mode scal fact:\t%lf\n#    o timestep:\t%2.2e ns\n#    o dumpstep:\t%2.2e ns\n#    o runtime:\t\t%2.2e ns\n#    o temp:\t\t%2.2e K\n#    o gam:\t\t%2.2e 1/ns\n#    o tol:\t\t%2.2e\n#    o sysmode:\t\t%s\n#    o statemode:\t%s\n#    o seed:\t\t%d\n", r.nmodes, r.modefact, r.dt, r.dt*r.nfreq, r.runtime, sv.T, mc.gam, s.tol, sysmode[s.systyp], statemode[s.statetyp], seed);

  // Random number initialization
  const gsl_rng_type * gT;
  gsl_rng * gr;
  gsl_rng_env_setup();
  gT    = gsl_rng_default;
  gr    = gsl_rng_alloc(gT);

  // System Initialization
  if (s.systyp==0)
    SysInit(&s, r, mc, sv, sc, dc);
  else
    SysRead(&s, r, mc, sv, sc, dc);

  if (s.statetyp==0)
    StateInit(&s, r, sv, gr);
  else
    StateRead(&s, r);

  // Integrate
  int     i;
  r.nsteps            = (int) (r.runtime/r.dt + s.step);

  printf("\n# 4. Integration started ...\n");
  for (i=s.istep; i<r.nsteps; i++)
  {
    //r.step = i;

    s.step+=1;
    IntegrateSys(&s, r, gr, sv);

  }
  StateStore(&s, r);
  gsl_rng_free(gr);
 
  clock_gettime(CLOCK_MONOTONIC_RAW, &end); 
  printf("# Total time = %f seconds\n",
            (end.tv_nsec - begin.tv_nsec) / 1000000000.0 +
            (end.tv_sec  - begin.tv_sec));
}

