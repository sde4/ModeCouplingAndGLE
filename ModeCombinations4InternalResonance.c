/*
 * This function finds the mode combinations which lead to non-zero coupling constant
 */

#include <omp.h>
#include "struct_def.h"

void ModeCombinations4InternalResonance(sys_var * s, run_param r, double tol) {
  printf("## 2.2 Mode combinations for Internal resonances ...\n");

  /**************************/
  // Mode combinations leading 
  // to Internal resonances!
  /**************************/

  int i, j, k, l;
  int cou1, cou2, cou3, cou4, cou5;
  int sind, pind, qind, rind;
  double frs, frp, frq, frr;


  cou1 = 0; // counts the total number of IR combinations
  s->IRs_pqrcombmat    = gsl_matrix_alloc(s->NZcou, 4);
  for (i = 0; i < s->NZcou; i++) {

    sind = gsl_matrix_get(s->NZs_pqrcombmat, i, 0);
    pind = gsl_matrix_get(s->NZs_pqrcombmat, i, 1);
    qind = gsl_matrix_get(s->NZs_pqrcombmat, i, 2);
    rind = gsl_matrix_get(s->NZs_pqrcombmat, i, 3);

    frs = gsl_vector_get(s->frvec, sind - 1);
    frp = gsl_vector_get(s->frvec, pind - 1);
    frq = gsl_vector_get(s->frvec, qind - 1);
    frr = gsl_vector_get(s->frvec, rind - 1);

    if (fabs(1 - fabs((frp + frq + frr) / frs)) < tol) {
      gsl_matrix_set(s->IRs_pqrcombmat, cou1, 0, sind); // mode number starts from 1
      gsl_matrix_set(s->IRs_pqrcombmat, cou1, 1, pind);
      gsl_matrix_set(s->IRs_pqrcombmat, cou1, 2, qind);
      gsl_matrix_set(s->IRs_pqrcombmat, cou1, 3, rind);
      cou1++;
    } else if (fabs(1 - fabs((frp + frq - frr) / frs)) < tol) {
      gsl_matrix_set(s->IRs_pqrcombmat, cou1, 0, sind);
      gsl_matrix_set(s->IRs_pqrcombmat, cou1, 1, pind);
      gsl_matrix_set(s->IRs_pqrcombmat, cou1, 2, qind);
      gsl_matrix_set(s->IRs_pqrcombmat, cou1, 3, rind);
      cou1++;
    } else if (fabs(1 - fabs((frp - frq + frr) / frs)) < tol) {
      gsl_matrix_set(s->IRs_pqrcombmat, cou1, 0, sind);
      gsl_matrix_set(s->IRs_pqrcombmat, cou1, 1, pind);
      gsl_matrix_set(s->IRs_pqrcombmat, cou1, 2, qind);
      gsl_matrix_set(s->IRs_pqrcombmat, cou1, 3, rind);
      cou1++;
    } else if (fabs(1 - fabs((frp - frq - frr) / frs)) < tol) {
      gsl_matrix_set(s->IRs_pqrcombmat, cou1, 0, sind);
      gsl_matrix_set(s->IRs_pqrcombmat, cou1, 1, pind);
      gsl_matrix_set(s->IRs_pqrcombmat, cou1, 2, qind);
      gsl_matrix_set(s->IRs_pqrcombmat, cou1, 3, rind);
      cou1++;
    }
  }
  s->IRcou = cou1;
}
