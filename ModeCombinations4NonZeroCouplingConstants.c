/*
 * This function finds the mode combinations which lead to non-zero coupling constant
 */

#include <omp.h>
#include "struct_def.h"

void ModeCombinations4NonZeroCouplingConstants(sys_var * s, run_param r) {
  printf("#    - Mode combinations for non-zero coupling constants ...\n");

  /**************************/
  // Mode combinations leading
  // to nonzero coupling 
  // constants
  /**************************/


  int i, j, k, l;
  int cou1, cou2, cou3, cou4, cou5;
  int sind, pind, qind, rind;
  FILE * outfp;
  
  // s->NZs_pqrcombmat    = gsl_matrix_alloc(s->SScou * s->SAcou * s->AScou * s->AAcou * (12+6*12+6*7+6*4), 6);
  gsl_matrix *NZs_pqrcombmat    = gsl_matrix_alloc(1, 6);
  // File write !!
  // NZ comb file !!
  outfp = fopen("NZs_pqrcombmat.dat", "w");


  cou1 = 0;
  cou2 = 0;   // ids such that all with same ids have same value

/**/
  // modes coming from 4  distinct groups (implies s, p, q, r are all distinct)
  // s		p	q	r
  // SS		SA	AS	AA  	24(12*2) combinations
  for (i=0; i<s->SScou; i++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, i);
    for (j=0; j<s->SAcou; j++) {
      pind = (int) gsl_vector_get(s->SAmodindvec, j);
      for (k=0; k<s->AScou; k++) {
        qind = (int) gsl_vector_get(s->ASmodindvec, k);
        for (l=0; l<s->AAcou; l++) {
          rind = (int) gsl_vector_get(s->AAmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;

	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
	  cou2++;

          // s p r q | s r p q
          // p s q r | p q s r
          // q p r s | q r p s
          // r s q p | r q s p
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;

	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;

	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;          
          cou2++;

          // s q r p | s r q p
          // p q r s | p r q s
          // q s p r | q p s r
          // r s p q | r p s q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;

	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and s, p, q, r all distinct
  // s		p	q	r
  // SS		SS	SA	SA 	24(12*2) combinations
  for (i=0; i<s->SScou; i++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, i);
    for (j=i+1; j<s->SScou; j++) {
      pind = (int) gsl_vector_get(s->SSmodindvec, j);
      for (k=0; k<s->SAcou; k++) {
        qind = (int) gsl_vector_get(s->SAmodindvec, k);
        for (l=k+1; l<s->SAcou; l++) {
          rind = (int) gsl_vector_get(s->SAmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
    	  cou2++;
  
          // s p r q | s r p q
          // p s q r | p q s r
          // q p r s | q r p s
          // r s q p | r q s p
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
   
   	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;          
          cou2++;
  
          // s q r p | s r q p
          // p q r s | p r q s
          // q s p r | q p s r
          // r s p q | r p s q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and s, p, q, r all distinct
  // s		p	q	r
  // SS		SS	AS	AS 	24(12*2) combinations
  for (i=0; i<s->SScou; i++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, i);
    for (j=i+1; j<s->SScou; j++) {
      pind = (int) gsl_vector_get(s->SSmodindvec, j);
      for (k=0; k<s->AScou; k++) {
        qind = (int) gsl_vector_get(s->ASmodindvec, k);
        for (l=k+1; l<s->AScou; l++) {
          rind = (int) gsl_vector_get(s->ASmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
    	  cou2++;
  
          // s p r q | s r p q
          // p s q r | p q s r
          // q p r s | q r p s
          // r s q p | r q s p
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
   
   	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;          
          cou2++;
  
          // s q r p | s r q p
          // p q r s | p r q s
          // q s p r | q p s r
          // r s p q | r p s q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and s, p, q, r all distinct
  // s		p	q	r
  // SS		SS	AA	AA 	24(12*2) combinations
  for (i=0; i<s->SScou; i++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, i);
    for (j=i+1; j<s->SScou; j++) {
      pind = (int) gsl_vector_get(s->SSmodindvec, j);
      for (k=0; k<s->AAcou; k++) {
        qind = (int) gsl_vector_get(s->AAmodindvec, k);
        for (l=k+1; l<s->AAcou; l++) {
          rind = (int) gsl_vector_get(s->AAmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
    	  cou2++;
  
          // s p r q | s r p q
          // p s q r | p q s r
          // q p r s | q r p s
          // r s q p | r q s p
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
   
   	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;          
          cou2++;
  
          // s q r p | s r q p
          // p q r s | p r q s
          // q s p r | q p s r
          // r s p q | r p s q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and s, p, q, r all distinct
  // s		p	q	r
  // SA		SA	AS	AS 	24(12*2) combinations
  for (i=0; i<s->SAcou; i++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, i);
    for (j=i+1; j<s->SAcou; j++) {
      pind = (int) gsl_vector_get(s->SAmodindvec, j);
      for (k=0; k<s->AScou; k++) {
        qind = (int) gsl_vector_get(s->ASmodindvec, k);
        for (l=k+1; l<s->AScou; l++) {
          rind = (int) gsl_vector_get(s->ASmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
    	  cou2++;
  
          // s p r q | s r p q
          // p s q r | p q s r
          // q p r s | q r p s
          // r s q p | r q s p
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
   
   	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;          
          cou2++;
  
          // s q r p | s r q p
          // p q r s | p r q s
          // q s p r | q p s r
          // r s p q | r p s q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and s, p, q, r all distinct
  // s		p	q	r
  // SA		SA	AA	AA 	24(12*2) combinations
  for (i=0; i<s->SAcou; i++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, i);
    for (j=i+1; j<s->SAcou; j++) {
      pind = (int) gsl_vector_get(s->SAmodindvec, j);
      for (k=0; k<s->AAcou; k++) {
        qind = (int) gsl_vector_get(s->AAmodindvec, k);
        for (l=k+1; l<s->AAcou; l++) {
          rind = (int) gsl_vector_get(s->AAmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
    	  cou2++;
  
          // s p r q | s r p q
          // p s q r | p q s r
          // q p r s | q r p s
          // r s q p | r q s p
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
   
   	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;          
          cou2++;
  
          // s q r p | s r q p
          // p q r s | p r q s
          // q s p r | q p s r
          // r s p q | r p s q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and s, p, q, r all distinct
  // s		p	q	r
  // AS		AS	AA	AA 	24(12*2) combinations
  for (i=0; i<s->AScou; i++) {
    sind = (int) gsl_vector_get(s->ASmodindvec, i);
    for (j=i+1; j<s->AScou; j++) {
      pind = (int) gsl_vector_get(s->ASmodindvec, j);
      for (k=0; k<s->AAcou; k++) {
        qind = (int) gsl_vector_get(s->AAmodindvec, k);
        for (l=k+1; l<s->AAcou; l++) {
          rind = (int) gsl_vector_get(s->AAmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
    	  cou2++;
  
          // s p r q | s r p q
          // p s q r | p q s r
          // q p r s | q r p s
          // r s q p | r q s p
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
   
   	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;          
          cou2++;
  
          // s q r p | s r q p
          // p q r s | p r q s
          // q s p r | q p s r
          // r s p q | r p s q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from the same group and s, p, q, r all distinct
  // s		p	q	r
  // SS		SS	SS	SS 	24(12*2) combinations
  for (i=0; i<s->SScou; i++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, i);
    for (j=i+1; j<s->SScou; j++) {
      pind = (int) gsl_vector_get(s->SSmodindvec, j);
      for (k=j+1; k<s->SScou; k++) {
        qind = (int) gsl_vector_get(s->SSmodindvec, k);
        for (l=k+1; l<s->SScou; l++) {
          rind = (int) gsl_vector_get(s->SSmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
    	  cou2++;
  
          // s p r q | s r p q
          // p s q r | p q s r
          // q p r s | q r p s
          // r s q p | r q s p
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
   
   	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;          
          cou2++;
  
          // s q r p | s r q p
          // p q r s | p r q s
          // q s p r | q p s r
          // r s p q | r p s q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from the same group and s, p, q, r all distinct
  // s		p	q	r
  // SA		SA	SA	SA 	24(12*2) combinations
  for (i=0; i<s->SAcou; i++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, i);
    for (j=i+1; j<s->SAcou; j++) {
      pind = (int) gsl_vector_get(s->SAmodindvec, j);
      for (k=j+1; k<s->SAcou; k++) {
        qind = (int) gsl_vector_get(s->SAmodindvec, k);
        for (l=k+1; l<s->SAcou; l++) {
          rind = (int) gsl_vector_get(s->SAmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
    	  cou2++;
  
          // s p r q | s r p q
          // p s q r | p q s r
          // q p r s | q r p s
          // r s q p | r q s p
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
   
   	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;          
          cou2++;
  
          // s q r p | s r q p
          // p q r s | p r q s
          // q s p r | q p s r
          // r s p q | r p s q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from the same group and s, p, q, r all distinct
  // s		p	q	r
  // AS		AS	AS	AS 	24(12*2) combinations
  for (i=0; i<s->AScou; i++) {
    sind = (int) gsl_vector_get(s->ASmodindvec, i);
    for (j=i+1; j<s->AScou; j++) {
      pind = (int) gsl_vector_get(s->ASmodindvec, j);
      for (k=j+1; k<s->AScou; k++) {
        qind = (int) gsl_vector_get(s->ASmodindvec, k);
        for (l=k+1; l<s->AScou; l++) {
          rind = (int) gsl_vector_get(s->ASmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
    	  cou2++;
  
          // s p r q | s r p q
          // p s q r | p q s r
          // q p r s | q r p s
          // r s q p | r q s p
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
   
   	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;          
          cou2++;
  
          // s q r p | s r q p
          // p q r s | p r q s
          // q s p r | q p s r
          // r s p q | r p s q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from the same group and s, p, q, r all distinct
  // s		p	q	r
  // AA		AA	AA	AA 	24(12*2) combinations
  for (i=0; i<s->AAcou; i++) {
    sind = (int) gsl_vector_get(s->AAmodindvec, i);
    for (j=i+1; j<s->AAcou; j++) {
      pind = (int) gsl_vector_get(s->AAmodindvec, j);
      for (k=j+1; k<s->AAcou; k++) {
        qind = (int) gsl_vector_get(s->AAmodindvec, k);
        for (l=k+1; l<s->AAcou; l++) {
          rind = (int) gsl_vector_get(s->AAmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
    	  cou2++;
  
          // s p r q | s r p q
          // p s q r | p q s r
          // q p r s | q r p s
          // r s q p | r q s p
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
   
   	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;          
          cou2++;
  
          // s q r p | s r q p
          // p q r s | p r q s
          // q s p r | q p s r
          // r s p q | r p s q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and p=s, r!=q
  // s		p	q	r
  // SS		SS	SA	SA 	4C2*2!=12 combinations
  for (i=0; i<s->SScou; i++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->SSmodindvec, j);
      for (k=0; k<s->SAcou; k++) {
        qind = (int) gsl_vector_get(s->SAmodindvec, k);
        for (l=k+1; l<s->SAcou; l++) {
          rind = (int) gsl_vector_get(s->SAmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
    	  cou2++;
  
          // s q r p | s r q p
          // q s p r
          // r s p q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and p=s, r!=q
  // s		p	q	r
  // SS		SS	AS	AS 	4C2*2!=12 combinations
  for (i=0; i<s->SScou; i++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->SSmodindvec, j);
      for (k=0; k<s->AScou; k++) {
        qind = (int) gsl_vector_get(s->ASmodindvec, k);
        for (l=k+1; l<s->AScou; l++) {
          rind = (int) gsl_vector_get(s->ASmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
    	  cou2++;
  
          // s q r p | s r q p
          // q s p r
          // r s p q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and p=s, r!=q
  // s		p	q	r
  // SS		SS	AA	AA 	4C2*2!=12 combinations
  for (i=0; i<s->SScou; i++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->SSmodindvec, j);
      for (k=0; k<s->AAcou; k++) {
        qind = (int) gsl_vector_get(s->AAmodindvec, k);
        for (l=k+1; l<s->AAcou; l++) {
          rind = (int) gsl_vector_get(s->AAmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
    	  cou2++;
  
          // s q r p | s r q p
          // q s p r
          // r s p q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and p=s, r!=q
  // s		p	q	r
  // SA		SA	AS	AS 	4C2*2!=12 combinations
  for (i=0; i<s->SAcou; i++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->SAmodindvec, j);
      for (k=0; k<s->AScou; k++) {
        qind = (int) gsl_vector_get(s->ASmodindvec, k);
        for (l=k+1; l<s->AScou; l++) {
          rind = (int) gsl_vector_get(s->ASmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
    	  cou2++;
  
          // s q r p | s r q p
          // q s p r
          // r s p q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and p=s, r!=q
  // s		p	q	r
  // SA		SA	AA	AA 	4C2*2!=12 combinations
  for (i=0; i<s->SAcou; i++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->SAmodindvec, j);
      for (k=0; k<s->AAcou; k++) {
        qind = (int) gsl_vector_get(s->AAmodindvec, k);
        for (l=k+1; l<s->AAcou; l++) {
          rind = (int) gsl_vector_get(s->AAmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
    	  cou2++;
  
          // s q r p | s r q p
          // q s p r
          // r s p q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and p=s, r!=q
  // s		p	q	r
  // AS		AS	AA	AA 	4C2*2!=12 combinations
  for (i=0; i<s->AScou; i++) {
    sind = (int) gsl_vector_get(s->ASmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->ASmodindvec, j);
      for (k=0; k<s->AAcou; k++) {
        qind = (int) gsl_vector_get(s->AAmodindvec, k);
        for (l=k+1; l<s->AAcou; l++) {
          rind = (int) gsl_vector_get(s->AAmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
    	  cou2++;
  
          // s q r p | s r q p
          // q s p r
          // r s p q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and p=s, r!=q
  // s		p	q	r
  // SA		SA	SS	SS 	4C2*2!=12 combinations
  for (i=0; i<s->SAcou; i++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->SAmodindvec, j);
      for (k=0; k<s->SScou; k++) {
        qind = (int) gsl_vector_get(s->SSmodindvec, k);
        for (l=k+1; l<s->SScou; l++) {
          rind = (int) gsl_vector_get(s->SSmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
    	  cou2++;
  
          // s q r p | s r q p
          // q s p r
          // r s p q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and p=s, r!=q
  // s		p	q	r
  // AS		AS	SS	SS 	4C2*2!=12 combinations
  for (i=0; i<s->AScou; i++) {
    sind = (int) gsl_vector_get(s->ASmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->ASmodindvec, j);
      for (k=0; k<s->SScou; k++) {
        qind = (int) gsl_vector_get(s->SSmodindvec, k);
        for (l=k+1; l<s->SScou; l++) {
          rind = (int) gsl_vector_get(s->SSmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
    	  cou2++;
  
          // s q r p | s r q p
          // q s p r
          // r s p q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and p=s, r!=q
  // s		p	q	r
  // AA		AA	SS	SS 	4C2*2!=12 combinations
  for (i=0; i<s->AAcou; i++) {
    sind = (int) gsl_vector_get(s->AAmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->AAmodindvec, j);
      for (k=0; k<s->SScou; k++) {
        qind = (int) gsl_vector_get(s->SSmodindvec, k);
        for (l=k+1; l<s->SScou; l++) {
          rind = (int) gsl_vector_get(s->SSmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
    	  cou2++;
  
          // s q r p | s r q p
          // q s p r
          // r s p q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and p=s, r!=q
  // s		p	q	r
  // AS		AS	SA	SA 	4C2*2!=12 combinations
  for (i=0; i<s->AScou; i++) {
    sind = (int) gsl_vector_get(s->ASmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->ASmodindvec, j);
      for (k=0; k<s->SAcou; k++) {
        qind = (int) gsl_vector_get(s->SAmodindvec, k);
        for (l=k+1; l<s->SAcou; l++) {
          rind = (int) gsl_vector_get(s->SAmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
    	  cou2++;
  
          // s q r p | s r q p
          // q s p r
          // r s p q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and p=s, r!=q
  // s		p	q	r
  // AA		AA	SA	SA 	4C2*2!=12 combinations
  for (i=0; i<s->AAcou; i++) {
    sind = (int) gsl_vector_get(s->AAmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->AAmodindvec, j);
      for (k=0; k<s->SAcou; k++) {
        qind = (int) gsl_vector_get(s->SAmodindvec, k);
        for (l=k+1; l<s->SAcou; l++) {
          rind = (int) gsl_vector_get(s->SAmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
    	  cou2++;
  
          // s q r p | s r q p
          // q s p r
          // r s p q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and p=s, r!=q
  // s		p	q	r
  // AA		AA	AS	AS 	4C2*2!=12 combinations
  for (i=0; i<s->AAcou; i++) {
    sind = (int) gsl_vector_get(s->AAmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->AAmodindvec, j);
      for (k=0; k<s->AScou; k++) {
        qind = (int) gsl_vector_get(s->ASmodindvec, k);
        for (l=k+1; l<s->AScou; l++) {
          rind = (int) gsl_vector_get(s->ASmodindvec, l);

          // s p q r | s q p r
          // p s r q | p r s q
          // q s r p | q r s p
          // r p q s | r q p s
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
    	  cou2++;
  
          // s q r p | s r q p
          // q s p r
          // r s p q
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from the same group and p=s, r!=q
  // s		p	q	r
  // SS		SS	SS	SS 	4C2*2!=12 combinations
  for (i=0; i<s->SScou; i++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->SSmodindvec, j);
      for (k=0; k<s->SScou; k++) {
        qind = (int) gsl_vector_get(s->SSmodindvec, k);
        for (l=k+1; l<s->SScou; l++) {
          rind = (int) gsl_vector_get(s->SSmodindvec, l);
	  if (qind!=sind && rind!=sind){

            // s p q r | s q p r
            // p s r q | p r s q
            // q s r p | q r s p
            // r p q s | r q p s
            gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
            
  	    gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
            
  	    gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
  
  	    gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
    	    cou2++;
  
            // s q r p | s r q p
            // q s p r
            // r s p q
            gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
  
  	    gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
            
  	    gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
            cou2++;
	  }
        }
      }
    }
  }


  // mode pairs coming from the same group and p=s, r!=q
  // s		p	q	r
  // SA		SA	SA	SA 	4C2*2!=12 combinations
  for (i=0; i<s->SAcou; i++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->SAmodindvec, j);
      for (k=0; k<s->SAcou; k++) {
        qind = (int) gsl_vector_get(s->SAmodindvec, k);
        for (l=k+1; l<s->SAcou; l++) {
          rind = (int) gsl_vector_get(s->SAmodindvec, l);
	  if (qind!=sind && rind!=sind){

            // s p q r | s q p r
            // p s r q | p r s q
            // q s r p | q r s p
            // r p q s | r q p s
            gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
            
  	    gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
            
  	    gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
  
  	    gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
    	    cou2++;
  
            // s q r p | s r q p
            // q s p r
            // r s p q
            gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
  
  	    gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
            
  	    gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
            cou2++;
	  }
        }
      }
    }
  }


  // mode pairs coming from the same group and p=s, r!=q
  // s		p	q	r
  // AS		AS	AS	AS 	4C2*2!=12 combinations
  for (i=0; i<s->AScou; i++) {
    sind = (int) gsl_vector_get(s->ASmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->ASmodindvec, j);
      for (k=0; k<s->AScou; k++) {
        qind = (int) gsl_vector_get(s->ASmodindvec, k);
        for (l=k+1; l<s->AScou; l++) {
          rind = (int) gsl_vector_get(s->ASmodindvec, l);
	  if (qind!=sind && rind!=sind){

            // s p q r | s q p r
            // p s r q | p r s q
            // q s r p | q r s p
            // r p q s | r q p s
            gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
            
  	    gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
            
  	    gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
  
  	    gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
    	    cou2++;
  
            // s q r p | s r q p
            // q s p r
            // r s p q
            gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
  
  	    gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
            
  	    gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
            cou2++;
	  }
        }
      }
    }
  }


  // mode pairs coming from the same group and p=s, r!=q
  // s		p	q	r
  // AA		AA	AA	AA 	4C2*2!=12 combinations
  for (i=0; i<s->AAcou; i++) {
    sind = (int) gsl_vector_get(s->AAmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->AAmodindvec, j);
      for (k=0; k<s->AAcou; k++) {
        qind = (int) gsl_vector_get(s->AAmodindvec, k);
        for (l=k+1; l<s->AAcou; l++) {
          rind = (int) gsl_vector_get(s->AAmodindvec, l);
	  if (qind!=sind && rind!=sind){

            // s p q r | s q p r
            // p s r q | p r s q
            // q s r p | q r s p
            // r p q s | r q p s
            gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
            
  	    gsl_matrix_set(NZs_pqrcombmat, 0, 0, pind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
            
  	    gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
  
  	    gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, sind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
    	    cou2++;
  
            // s q r p | s r q p
            // q s p r
            // r s p q
            gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
  
  	    gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
            
  	    gsl_matrix_set(NZs_pqrcombmat, 0, 0, rind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
            gsl_matrix_set(NZs_pqrcombmat, 0, 3, qind);
            gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
            gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	    fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
            cou1++;
            cou2++;
	  }
        }
      }
    }
  }


/**/
  // mode pairs coming from 2  distinct groups and p=s, r=q
  // s		p	q	r
  // SS		SS	SA	SA 	4C2=6 combinations
  for (i=0; i<s->SScou; i++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->SSmodindvec, j);
      for (k=0; k<s->SAcou; k++) {
        qind = (int) gsl_vector_get(s->SAmodindvec, k);
        for (l=k; l<k+1; l++) {
          rind = (int) gsl_vector_get(s->SAmodindvec, l);

          // s p q r | s q p r
          // q s r p | q r s p
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++; 
    	  cou2++;
  
          // s q r p
          // q s p r
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and p=s, r=q
  // s		p	q	r
  // SS		SS	AS	AS 	4C2=6 combinations
  for (i=0; i<s->SScou; i++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->SSmodindvec, j);
      for (k=0; k<s->AScou; k++) {
        qind = (int) gsl_vector_get(s->ASmodindvec, k);
        for (l=k; l<k+1; l++) {
          rind = (int) gsl_vector_get(s->ASmodindvec, l);

          // s p q r | s q p r
          // q s r p | q r s p
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++; 
    	  cou2++;
  
          // s q r p
          // q s p r
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and p=s, r=q
  // s		p	q	r
  // SS		SS	AA	AA 	4C2=6 combinations
  for (i=0; i<s->SScou; i++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->SSmodindvec, j);
      for (k=0; k<s->AAcou; k++) {
        qind = (int) gsl_vector_get(s->AAmodindvec, k);
        for (l=k; l<k+1; l++) {
          rind = (int) gsl_vector_get(s->AAmodindvec, l);

          // s p q r | s q p r
          // q s r p | q r s p
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++; 
    	  cou2++;
  
          // s q r p
          // q s p r
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and p=s, r=q
  // s		p	q	r
  // SA		SA	AS	AS 	4C2=6 combinations
  for (i=0; i<s->SAcou; i++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->SAmodindvec, j);
      for (k=0; k<s->AScou; k++) {
        qind = (int) gsl_vector_get(s->ASmodindvec, k);
        for (l=k; l<k+1; l++) {
          rind = (int) gsl_vector_get(s->ASmodindvec, l);

          // s p q r | s q p r
          // q s r p | q r s p
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++; 
    	  cou2++;
  
          // s q r p
          // q s p r
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and p=s, r=q
  // s		p	q	r
  // SA		SA	AA	AA 	4C2=6 combinations
  for (i=0; i<s->SAcou; i++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->SAmodindvec, j);
      for (k=0; k<s->AAcou; k++) {
        qind = (int) gsl_vector_get(s->AAmodindvec, k);
        for (l=k; l<k+1; l++) {
          rind = (int) gsl_vector_get(s->AAmodindvec, l);

          // s p q r | s q p r
          // q s r p | q r s p
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++; 
    	  cou2++;
  
          // s q r p
          // q s p r
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from 2  distinct groups and p=s, r=q
  // s		p	q	r
  // AS		AS	AA	AA 	4C2=6 combinations
  for (i=0; i<s->AScou; i++) {
    sind = (int) gsl_vector_get(s->ASmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->ASmodindvec, j);
      for (k=0; k<s->AAcou; k++) {
        qind = (int) gsl_vector_get(s->AAmodindvec, k);
        for (l=k; l<k+1; l++) {
          rind = (int) gsl_vector_get(s->AAmodindvec, l);

          // s p q r | s q p r
          // q s r p | q r s p
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++; 
    	  cou2++;
  
          // s q r p
          // q s p r
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from the same group and p=s, r=q
  // s		p	q	r
  // SS		SS	SS	SS 	4C2=6 combinations
  for (i=0; i<s->SScou; i++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->SSmodindvec, j);
      for (k=i+1; k<s->SScou; k++) {
        qind = (int) gsl_vector_get(s->SSmodindvec, k);
        for (l=k; l<k+1; l++) {
          rind = (int) gsl_vector_get(s->SSmodindvec, l);

          // s p q r | s q p r
          // q s r p | q r s p
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++; 
    	  cou2++;
  
          // s q r p
          // q s p r
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from the same group and p=s, r=q
  // s		p	q	r
  // SA		SA	SA	SA 	4C2=6 combinations
  for (i=0; i<s->SAcou; i++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->SAmodindvec, j);
      for (k=i+1; k<s->SAcou; k++) {
        qind = (int) gsl_vector_get(s->SAmodindvec, k);
        for (l=k; l<k+1; l++) {
          rind = (int) gsl_vector_get(s->SAmodindvec, l);

          // s p q r | s q p r
          // q s r p | q r s p
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++; 
    	  cou2++;
  
          // s q r p
          // q s p r
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from the same group and p=s, r=q
  // s		p	q	r
  // AS		AS	AS	AS 	4C2=6 combinations
  for (i=0; i<s->AScou; i++) {
    sind = (int) gsl_vector_get(s->ASmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->ASmodindvec, j);
      for (k=i+1; k<s->AScou; k++) {
        qind = (int) gsl_vector_get(s->ASmodindvec, k);
        for (l=k; l<k+1; l++) {
          rind = (int) gsl_vector_get(s->ASmodindvec, l);

          // s p q r | s q p r
          // q s r p | q r s p
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++; 
    	  cou2++;
  
          // s q r p
          // q s p r
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from the same group and p=s, r=q
  // s		p	q	r
  // AA		AA	AA	AA 	4C2=6 combinations
  for (i=0; i<s->AAcou; i++) {
    sind = (int) gsl_vector_get(s->AAmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->AAmodindvec, j);
      for (k=i+1; k<s->AAcou; k++) {
        qind = (int) gsl_vector_get(s->AAmodindvec, k);
        for (l=k; l<k+1; l++) {
          rind = (int) gsl_vector_get(s->AAmodindvec, l);

          // s p q r | s q p r
          // q s r p | q r s p
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 2);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++; 
    	  cou2++;
  
          // s q r p
          // q s p r
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, rind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, pind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
  
  	  gsl_matrix_set(NZs_pqrcombmat, 0, 0, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from the same group and p=s=r=q
  // s		p	q	r
  // SS		SS	SS	SS 	4C2=6 combinations
  for (i=0; i<s->SScou; i++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->SSmodindvec, j);
      for (k=i; k<i+1; k++) {
        qind = (int) gsl_vector_get(s->SSmodindvec, k);
        for (l=i; l<i+1; l++) {
          rind = (int) gsl_vector_get(s->SSmodindvec, l);

          // s p q r
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from the same group and p=s=r=q
  // s		p	q	r
  // SA		SA	SA	SA 	4C2=6 combinations
  for (i=0; i<s->SAcou; i++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->SAmodindvec, j);
      for (k=i; k<i+1; k++) {
        qind = (int) gsl_vector_get(s->SAmodindvec, k);
        for (l=i; l<i+1; l++) {
          rind = (int) gsl_vector_get(s->SAmodindvec, l);

          // s p q r
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from the same group and p=s=r=q
  // s		p	q	r
  // AS		AS	AS	AS 	4C2=6 combinations
  for (i=0; i<s->AScou; i++) {
    sind = (int) gsl_vector_get(s->ASmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->ASmodindvec, j);
      for (k=i; k<i+1; k++) {
        qind = (int) gsl_vector_get(s->ASmodindvec, k);
        for (l=i; l<i+1; l++) {
          rind = (int) gsl_vector_get(s->ASmodindvec, l);

          // s p q r
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }


  // mode pairs coming from the same group and p=s=r=q
  // s		p	q	r
  // AA		AA	AA	AA 	4C2=6 combinations
  for (i=0; i<s->AAcou; i++) {
    sind = (int) gsl_vector_get(s->AAmodindvec, i);
    for (j=i; j<i+1; j++) {
      pind = (int) gsl_vector_get(s->AAmodindvec, j);
      for (k=i; k<i+1; k++) {
        qind = (int) gsl_vector_get(s->AAmodindvec, k);
        for (l=i; l<i+1; l++) {
          rind = (int) gsl_vector_get(s->AAmodindvec, l);
  
          // s p q r
          gsl_matrix_set(NZs_pqrcombmat, 0, 0, sind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 1, pind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 2, qind); 
          gsl_matrix_set(NZs_pqrcombmat, 0, 3, rind);
          gsl_matrix_set(NZs_pqrcombmat, 0, 4, cou2);
          gsl_matrix_set(NZs_pqrcombmat, 0, 5, 1);
	  fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, 0, 0), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 1), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 2), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 3), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 4), (int) gsl_matrix_get(NZs_pqrcombmat, 0, 5));
          cou1++;
          cou2++;
        }
      }
    }
  }
/**/

  s->NZcou = cou1;
  fclose(outfp);
  gsl_matrix_free(NZs_pqrcombmat);

  // File write !!
  // NZ comb file !!
  // outfp = fopen("NZs_pqrcombmat.dat", "w");
  // for (i = 0; i < s->NZcou - 1; i++) {
  //   fprintf(outfp, "%d\t", (int) gsl_matrix_get(NZs_pqrcombmat, i, 0));
  //   fprintf(outfp, "%d\t", (int) gsl_matrix_get(NZs_pqrcombmat, i, 1));
  //   fprintf(outfp, "%d\t", (int) gsl_matrix_get(NZs_pqrcombmat, i, 2));
  //   fprintf(outfp, "%d\t", (int) gsl_matrix_get(NZs_pqrcombmat, i, 3));
  //   fprintf(outfp, "%d\t", (int) gsl_matrix_get(NZs_pqrcombmat, i, 4));
  //   fprintf(outfp, "%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, i, 5));
  // }
  // fprintf(outfp, "%d\t", (int) gsl_matrix_get(NZs_pqrcombmat, i, 0));
  // fprintf(outfp, "%d\t", (int) gsl_matrix_get(NZs_pqrcombmat, i, 1));
  // fprintf(outfp, "%d\t", (int) gsl_matrix_get(NZs_pqrcombmat, i, 2));
  // fprintf(outfp, "%d\t", (int) gsl_matrix_get(NZs_pqrcombmat, i, 3));
  // fprintf(outfp, "%d\t", (int) gsl_matrix_get(NZs_pqrcombmat, i, 4));
  // fprintf(outfp, "%d\n", (int) gsl_matrix_get(NZs_pqrcombmat, i, 5));
  // fclose(outfp);
  
  return;
}
