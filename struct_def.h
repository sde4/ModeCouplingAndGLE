#ifndef strct_def
#define strct_def

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#define PI            3.1415926535
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
  gsl_vector    *frvec;
  gsl_vector    *gamvec;
  gsl_vector    *mvec;
  gsl_matrix    *alphamat;
  gsl_vector    *qvec;
  gsl_vector    *qdotvec;
  gsl_vector    *sigvec;
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


#endif
