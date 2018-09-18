#ifndef strct_def
#define strct_def

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#define PI            3.14159265358979323846
#define kB            8.6173303E-5            // eV/K
#define Npm_eVpA2     0.0624150912
#define kgpm2_amupA2  6.022E6
#define amuA2pns2_eV  1.03642696E-10

typedef struct {
  double        dt;                            // ns
  double        runtime;                       // ns
  int           nfreq;
  int           nmodes;
  int           step;
  int           nsteps;
  }             run_param;

typedef struct {
  gsl_matrix    *modindmat;
  gsl_vector    *modindvec;
  gsl_vector    *SSmodindvec;
  gsl_vector    *SAmodindvec;
  gsl_vector    *ASmodindvec;
  gsl_vector    *AAmodindvec;
  gsl_vector    *frvec;
  gsl_matrix    *IRs_pqrcombmatsorted;
  gsl_matrix    *IRs_pqrcountmat;
  gsl_matrix    *alphamatsorted;
  gsl_vector    *detuningsorted;
  gsl_vector    *gamvec;
  gsl_vector    *mvec;
  gsl_vector    *qvec;
  gsl_vector    *qdotvec;
  gsl_vector    *fvec;
  gsl_vector    *enonvec;
  gsl_vector    *sigvec;
  int           SScou;
  int           SAcou;
  int           AScou;
  int           AAcou;
  int           NZcou;
  int           IRcou;
  int           NmodeIRcou;
  double 	tol;
  }             sys_var;

typedef struct {
  double        Et;                           // eV/A^2
  double        DEt;                          // eV/A^2
  double        rho;                          // eV/(A^4/ns^2)
  double        gam;                          // 1/ns
  double        alpha;                        
  }             mat_const;

typedef struct {
  double        T;                            // K
  double        e_pre;
  }             state_var;

typedef struct {
  int           mmax;
  int           nmax;
  double        Lx;
  double        Ly;
  int           run_id;
  }             sys_const;

typedef struct {
  int           Nfx;
  int           Nfy;
  double        Lx;
  double        Ly;
  double        hx;
  double        hy;
  double 	gamunitconv;
  }             disc_const;
  
typedef struct {
  double        f1;
  double        f3;                       
  double        ep4;                       
  }             for_var;

#endif
