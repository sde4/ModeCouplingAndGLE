/*
 * This function finds the mode combinations which lead to non-zero coupling constant
 */

#include <omp.h>
#include "struct_def.h"

void ModeCombinations4InternalResonance(sys_var * s) {
  printf("#    - Mode combinations for Internal resonances ...\n");

  /**************************/
  // Mode combinations leading 
  // to Internal resonances!
  /**************************/

  int i, j, k, l;
  int cou1, cou2, cou3, cou4, cou5;
  int sind, tsind, pind, qind, rind;
  int idat1, idat2;
  double frs, frp, frq, frr;
  double fdat1, fdat2;
  FILE *infp, *infp1, *outfp, *outfp1, *outfp2;

  infp = fopen("NZs_pqrcombmat.dat", "r");
  gsl_matrix *IRs_pqrcombmat    = gsl_matrix_alloc(1, 6);
  outfp1 = fopen("IRs_pqrcombmat.dat", "w");
  outfp2 = fopen("detuning.dat", "w");

  cou1 = 0; // counts the total number of IR combinations
  for (i = 0; i < s->NZcou; i++) {  
    // sind  = (int) gsl_matrix_get(IRs_pqrcombmat, i, 0);
    // pind  = (int) gsl_matrix_get(IRs_pqrcombmat, i, 1);
    // qind  = (int) gsl_matrix_get(IRs_pqrcombmat, i, 2);
    // rind  = (int) gsl_matrix_get(IRs_pqrcombmat, i, 3);
    // idat1 = (int) gsl_matrix_get(IRs_pqrcombmat, i, 4);
    // idat2 = (int) gsl_matrix_get(IRs_pqrcombmat, i, 5);
    fscanf(infp, "%d %d %d %d %d %d", &sind, &pind, &qind, &rind, &idat1, &idat2);

    /*
    gsl_matrix_set(IRs_pqrcombmat, i, 0, sind);
    gsl_matrix_set(IRs_pqrcombmat, i, 1, pind);
    gsl_matrix_set(IRs_pqrcombmat, i, 2, qind);
    gsl_matrix_set(IRs_pqrcombmat, i, 3, rind);
    cou1++;
    */

    frs = gsl_vector_get(s->frvec, sind - 1);
    frp = gsl_vector_get(s->frvec, pind - 1);
    frq = gsl_vector_get(s->frvec, qind - 1);
    frr = gsl_vector_get(s->frvec, rind - 1);

    if (fabs( frs-fabs(frp+frq+frr) )/(frs+frp+frq+frr) < s->tol) {
      gsl_matrix_set(IRs_pqrcombmat, 0, 0, sind); // mode number starts from 1
      gsl_matrix_set(IRs_pqrcombmat, 0, 1, pind);
      gsl_matrix_set(IRs_pqrcombmat, 0, 2, qind);
      gsl_matrix_set(IRs_pqrcombmat, 0, 3, rind);
      gsl_matrix_set(IRs_pqrcombmat, 0, 4, idat1);
      gsl_matrix_set(IRs_pqrcombmat, 0, 5, idat2);
      fprintf(outfp1, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(IRs_pqrcombmat, 0, 0), (int) gsl_matrix_get(IRs_pqrcombmat, 0, 1), (int) gsl_matrix_get(IRs_pqrcombmat, 0, 2), (int) gsl_matrix_get(IRs_pqrcombmat, 0, 3), (int) gsl_matrix_get(IRs_pqrcombmat, 0, 4), (int) gsl_matrix_get(IRs_pqrcombmat, 0, 5));
      fprintf(outfp2, "%5.5e\n", fabs( frs-fabs(frp+frq+frr) )/(frs+frp+frq+frr));
      cou1++;
    } else if (fabs( frs-fabs(frp+frq-frr) )/(frs+frp+frq+frr) < s->tol) {
      gsl_matrix_set(IRs_pqrcombmat, 0, 0, sind);
      gsl_matrix_set(IRs_pqrcombmat, 0, 1, pind);
      gsl_matrix_set(IRs_pqrcombmat, 0, 2, qind);
      gsl_matrix_set(IRs_pqrcombmat, 0, 3, rind);
      gsl_matrix_set(IRs_pqrcombmat, 0, 4, idat1);
      gsl_matrix_set(IRs_pqrcombmat, 0, 5, idat2);
      fprintf(outfp1, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(IRs_pqrcombmat, 0, 0), (int) gsl_matrix_get(IRs_pqrcombmat, 0, 1), (int) gsl_matrix_get(IRs_pqrcombmat, 0, 2), (int) gsl_matrix_get(IRs_pqrcombmat, 0, 3), (int) gsl_matrix_get(IRs_pqrcombmat, 0, 4), (int) gsl_matrix_get(IRs_pqrcombmat, 0, 5));
      fprintf(outfp2, "%5.5e\n", fabs( frs-fabs(frp+frq-frr) )/(frs+frp+frq+frr));
      cou1++;
    } else if (fabs( frs-fabs(frp-frq+frr) )/(frs+frp+frq+frr) < s->tol) {
      gsl_matrix_set(IRs_pqrcombmat, 0, 0, sind);
      gsl_matrix_set(IRs_pqrcombmat, 0, 1, pind);
      gsl_matrix_set(IRs_pqrcombmat, 0, 2, qind);
      gsl_matrix_set(IRs_pqrcombmat, 0, 3, rind);
      gsl_matrix_set(IRs_pqrcombmat, 0, 4, idat1);
      gsl_matrix_set(IRs_pqrcombmat, 0, 5, idat2);
      fprintf(outfp1, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(IRs_pqrcombmat, 0, 0), (int) gsl_matrix_get(IRs_pqrcombmat, 0, 1), (int) gsl_matrix_get(IRs_pqrcombmat, 0, 2), (int) gsl_matrix_get(IRs_pqrcombmat, 0, 3), (int) gsl_matrix_get(IRs_pqrcombmat, 0, 4), (int) gsl_matrix_get(IRs_pqrcombmat, 0, 5));
      fprintf(outfp2, "%5.5e\n", fabs( frs-fabs(frp-frq+frr) )/(frs+frp+frq+frr));
      cou1++;
    } else if (fabs( frs-fabs(frp-frq-frr) )/(frs+frp+frq+frr) < s->tol ) {
      gsl_matrix_set(IRs_pqrcombmat, 0, 0, sind);
      gsl_matrix_set(IRs_pqrcombmat, 0, 1, pind);
      gsl_matrix_set(IRs_pqrcombmat, 0, 2, qind);
      gsl_matrix_set(IRs_pqrcombmat, 0, 3, rind);
      gsl_matrix_set(IRs_pqrcombmat, 0, 4, idat1);
      gsl_matrix_set(IRs_pqrcombmat, 0, 5, idat2);
      fprintf(outfp1, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(IRs_pqrcombmat, 0, 0), (int) gsl_matrix_get(IRs_pqrcombmat, 0, 1), (int) gsl_matrix_get(IRs_pqrcombmat, 0, 2), (int) gsl_matrix_get(IRs_pqrcombmat, 0, 3), (int) gsl_matrix_get(IRs_pqrcombmat, 0, 4), (int) gsl_matrix_get(IRs_pqrcombmat, 0, 5));
      fprintf(outfp2, "%5.5e\n", fabs( frs-fabs(frp-frq-frr) )/(frs+frp+frq+frr));
      cou1++;
    }
  }
  fclose(infp);
  fclose(outfp1);
  fclose(outfp2);

  s->IRcou = cou1;
  gsl_matrix_free(IRs_pqrcombmat);
   
  return;
}
