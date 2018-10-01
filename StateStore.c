/*
 * This function initializes the system of Langevin equation
 */

#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h> 
#include <omp.h>
#include "struct_def.h"

void StateStore(sys_var* s, run_param r) {
  printf("\n# 3. Storing system state ...\n");

  int i, j, k, l, cou1, cou2, cou3, cou4, cou5;
  int idat;
  int sind;
  double fdat;
  double q_t, qdot_t;
  FILE *outfp;
  char outfname[40];

  // trajectory file !!
  outfp = fopen("read.restart", "w");
  fprintf(outfp, "%d\n", s->step);
  // printf("%d\n", s->step);
  for (i=0; i<s->NmodeIRcou; i++){
    sind      =  (int) gsl_matrix_get(s->IRs_pqrcountmat, i, 0);
    
    q_t        = gsl_vector_get(s->qvec, sind-1);
    qdot_t     = gsl_vector_get(s->qdotvec, sind-1);
    fprintf(outfp, "%d %5.5e %5.5e %5.5e %5.5e\n", sind, q_t, qdot_t, 0.0, 0.0);
  }
  fclose(outfp);

  return;
}
