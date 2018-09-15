/*
 * This function finds the mode combinations which lead to non-zero coupling constant
 */

#include <omp.h>
#include "struct_def.h"

void ModeCombinations4InternalResonance(sys_var * s, double tol) {
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
  FILE *infp, *infp1, *outfp, *outfp1;

  infp = fopen("NZs_pqrcombmat.dat", "r");
  // s->IRs_pqrcombmat    = gsl_matrix_alloc(s->NZcou, 6);
  s->IRs_pqrcombmat    = gsl_matrix_alloc(1, 6);
  outfp = fopen("IRs_pqrcombmat.dat", "w");
  outfp1 = fopen("detuning.dat", "w");

  cou1 = 0; // counts the total number of IR combinations
  for (i = 0; i < s->NZcou; i++) {  
    // sind  = (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 0);
    // pind  = (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 1);
    // qind  = (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 2);
    // rind  = (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 3);
    // idat1 = (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 4);
    // idat2 = (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 5);
    fscanf(infp, "%d %d %d %d %d %d", &sind, &pind, &qind, &rind, &idat1, &idat2);

    /*
    gsl_matrix_set(s->IRs_pqrcombmat, i, 0, sind);
    gsl_matrix_set(s->IRs_pqrcombmat, i, 1, pind);
    gsl_matrix_set(s->IRs_pqrcombmat, i, 2, qind);
    gsl_matrix_set(s->IRs_pqrcombmat, i, 3, rind);
    cou1++;
    */

    frs = gsl_vector_get(s->frvec, sind - 1);
    frp = gsl_vector_get(s->frvec, pind - 1);
    frq = gsl_vector_get(s->frvec, qind - 1);
    frr = gsl_vector_get(s->frvec, rind - 1);

    if (fabs(1 - fabs((frp + frq + frr) / frs)) < tol) {
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 0, sind); // mode number starts from 1
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 1, pind);
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 2, qind);
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 3, rind);
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 4, idat1);
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 5, idat2);
      fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 0), (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 1), (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 2), (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 3), (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 4), (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 5));
      fprintf(outfp1, "%5.5e\n", fabs(1 - fabs((frp + frq + frr) / frs)));
      cou1++;
    } else if (fabs(1 - fabs((frp + frq - frr) / frs)) < tol) {
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 0, sind);
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 1, pind);
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 2, qind);
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 3, rind);
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 4, idat1);
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 5, idat2);
      fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 0), (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 1), (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 2), (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 3), (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 4), (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 5));
      fprintf(outfp1, "%5.5e\n", fabs(1 - fabs((frp + frq - frr) / frs)));
      cou1++;
    } else if (fabs(1 - fabs((frp - frq + frr) / frs)) < tol) {
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 0, sind);
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 1, pind);
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 2, qind);
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 3, rind);
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 4, idat1);
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 5, idat2);
      fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 0), (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 1), (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 2), (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 3), (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 4), (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 5));
      fprintf(outfp1, "%5.5e\n", fabs(1 - fabs((frp - frq + frr) / frs)));
      cou1++;
    } else if (fabs(1 - fabs((frp - frq - frr) / frs)) < tol) {
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 0, sind);
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 1, pind);
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 2, qind);
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 3, rind);
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 4, idat1);
      gsl_matrix_set(s->IRs_pqrcombmat, 0, 5, idat2);
      fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 0), (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 1), (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 2), (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 3), (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 4), (int) gsl_matrix_get(s->IRs_pqrcombmat, 0, 5));
      fprintf(outfp1, "%5.5e\n", fabs(1 - fabs((frp - frq - frr) / frs)));
      cou1++;
    }
  }
  fclose(infp);
  fclose(outfp);
  fclose(outfp1);

  s->IRcou = cou1;
  s->IRs_pqrcombmatsorted    = gsl_matrix_alloc(s->IRcou, 6);
  gsl_vector *detuningsorted             = gsl_vector_alloc(s->IRcou);
  cou1=0; // counts the total number of IR combinations
  cou3=0; // counts the number of distinct modes participating in IR
  sind=1;
  while(cou1<s->IRcou){
    cou2=0; // counts the number of IR combinations for mode sind
    infp = fopen("IRs_pqrcombmat.dat", "r");
    infp1 = fopen("detuning.dat", "r");
    for (j = 0; j < s->IRcou; j++){
      fscanf(infp, "%d %d %d %d %d %d", &tsind, &pind, &qind, &rind, &idat1, &idat2);
      fscanf(infp1, "%lf", &fdat1);
      // if (gsl_matrix_get(s->IRs_pqrcombmat, j, 0)== sind){
      if (tsind == sind){
        // pind  = (int) gsl_matrix_get(s->IRs_pqrcombmat, j, 1);
	// qind  = (int) gsl_matrix_get(s->IRs_pqrcombmat, j, 2);
	// rind  = (int) gsl_matrix_get(s->IRs_pqrcombmat, j, 3);
	// idat1 = (int) gsl_matrix_get(s->IRs_pqrcombmat, j, 4);
	// idat2 = (int) gsl_matrix_get(s->IRs_pqrcombmat, j, 5);

	gsl_matrix_set(s->IRs_pqrcombmatsorted, cou1, 0, sind);
	gsl_matrix_set(s->IRs_pqrcombmatsorted, cou1, 1, pind);
	gsl_matrix_set(s->IRs_pqrcombmatsorted, cou1, 2, qind);
	gsl_matrix_set(s->IRs_pqrcombmatsorted, cou1, 3, rind);	
        gsl_matrix_set(s->IRs_pqrcombmatsorted, cou1, 4, idat1);
        gsl_matrix_set(s->IRs_pqrcombmatsorted, cou1, 5, idat2);
        gsl_vector_set(detuningsorted, cou1, fdat1);
	cou2++;
	cou1++;
      }
    }
    fclose(infp);
    fclose(infp1);
    if (cou2>0){
      gsl_matrix_set(s->IRs_pqrcountmat, cou3, 0, sind);
      gsl_matrix_set(s->IRs_pqrcountmat, cou3, 1, cou2);
      gsl_matrix_set(s->IRs_pqrcountmat, cou3, 2, cou1);
      cou3++;
    }
    sind++;
  }
  s->NmodeIRcou = cou3;
  
  // File write !!
  // IR comb file !!
  // outfp = fopen("IRs_pqrcombmat.dat", "w");
  // for (i = 0; i < s->IRcou - 1; i++) {
  //   fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 0));
  //   fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 1));
  //   fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 2));
  //   fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 3));
  //   fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 4));
  //   fprintf(outfp, "%d\n", (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 5));
  // }
  // fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 0));
  // fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 1));
  // fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 2));
  // fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 3));
  // fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 4));
  // fprintf(outfp, "%d\n", (int) gsl_matrix_get(s->IRs_pqrcombmat, i, 5));
  // fclose(outfp);

  // IR comb sorted file !!
  outfp = fopen("IRs_pqrcombmatsorted.dat", "w");
  for (i = 0; i < s->IRcou - 1; i++) {
    fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 0));
    fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 1));
    fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 2));
    fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 3));
    fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 4));
    fprintf(outfp, "%d\n", (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 5));
  }
  fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 0));
  fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 1));
  fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 2));
  fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 3));
  fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 4));
  fprintf(outfp, "%d\n", (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 5));
  fclose(outfp);

  // detuning sorted file !!
  outfp = fopen("detuningsorted.dat", "w");
  for (i = 0; i < s->IRcou - 1; i++) {
    fprintf(outfp, "%5.5e\n", gsl_vector_get(detuningsorted, i));
  }
  fprintf(outfp, "%5.5e", gsl_vector_get(detuningsorted, i));
  fclose(outfp);

  // IR count file !!
  outfp = fopen("IRs_pqrcountmat.dat", "w");
  for (i = 0; i < s->NmodeIRcou - 1; i++) {
    fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcountmat, i, 0));
    fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcountmat, i, 1));
    fprintf(outfp, "%d\n", (int) gsl_matrix_get(s->IRs_pqrcountmat, i, 2));
  }
  fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcountmat, i, 0));
  fprintf(outfp, "%d\t", (int) gsl_matrix_get(s->IRs_pqrcountmat, i, 1));
  fprintf(outfp, "%d\n", (int) gsl_matrix_get(s->IRs_pqrcountmat, i, 2));
  fclose(outfp);
 
  gsl_vector_free(detuningsorted);
  
  return;
}
