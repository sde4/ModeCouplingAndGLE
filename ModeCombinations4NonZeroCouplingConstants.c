/*
 * This function finds the mode combinations which lead to non-zero coupling constant
 */

#include <omp.h>
#include "struct_def.h"

void ModeCombinations4NonZeroCouplingConstants(sys_var * s, run_param r) {
  printf("##   -Mode combinations for non-zero coupling constants ...\n");

  /**************************/
  // Mode combinations leading
  // to nonzero coupling 
  // constants
  /**************************/


  int i, j, k, l;
  int cou1, cou2, cou3, cou4, cou5;
  int sind, pind, qind, rind;

  // s=1
  sind = 1;
  cou1 = 0;
  // modes coming from 4 distinct group
  // s		p	q	r
  // SS		SA	AS	AA
  // 		SA	AA	AS
  // 		AA	SA	AS
  // 		AA	AS	SA
  // 		AS	AA	SA
  // 		AS	SA	AA
  for (j = 0; j < r.nmodes / 4; j++) {
    pind = (int) gsl_vector_get(s->SAmodindvec, j);
    for (k = 0; k < r.nmodes / 4; k++) {
      qind = (int) gsl_vector_get(s->ASmodindvec, k);
      for (l = 0; l < r.nmodes / 4; l++) {
        rind = (int) gsl_vector_get(s->AAmodindvec, l);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
        cou1++;
      }
    }
  }

  for (j = 0; j < r.nmodes / 4; j++) {
    pind = (int) gsl_vector_get(s->SAmodindvec, j);
    for (k = 0; k < r.nmodes / 4; k++) {
      qind = (int) gsl_vector_get(s->AAmodindvec, k);
      for (l = 0; l < r.nmodes / 4; l++) {
        rind = (int) gsl_vector_get(s->ASmodindvec, l);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
        cou1++;
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    pind = (int) gsl_vector_get(s->AAmodindvec, j);
    for (k = 0; k < r.nmodes / 4; k++) {
      qind = (int) gsl_vector_get(s->SAmodindvec, k);
      for (l = 0; l < r.nmodes / 4; l++) {
        rind = (int) gsl_vector_get(s->ASmodindvec, l);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
        cou1++;
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    pind = (int) gsl_vector_get(s->AAmodindvec, j);
    for (k = 0; k < r.nmodes / 4; k++) {
      qind = (int) gsl_vector_get(s->ASmodindvec, k);
      for (l = 0; l < r.nmodes / 4; l++) {
        rind = (int) gsl_vector_get(s->SAmodindvec, l);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
        cou1++;
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    pind = (int) gsl_vector_get(s->ASmodindvec, j);
    for (k = 0; k < r.nmodes / 4; k++) {
      qind = (int) gsl_vector_get(s->AAmodindvec, k);
      for (l = 0; l < r.nmodes / 4; l++) {
        rind = (int) gsl_vector_get(s->SAmodindvec, l);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
        cou1++;
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    pind = (int) gsl_vector_get(s->ASmodindvec, j);
    for (k = 0; k < r.nmodes / 4; k++) {
      qind = (int) gsl_vector_get(s->SAmodindvec, k);
      for (l = 0; l < r.nmodes / 4; l++) {
        rind = (int) gsl_vector_get(s->AAmodindvec, l);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
        cou1++;
      }
    }
  }
  // mode pairs coming from 2 distinct group
  // s		p	q	r
  // SS		SS	SA	SA
  // 			AS	AS
  // 			AA	AA
  for (j = 0; j < r.nmodes / 4; j++) {
    pind = (int) gsl_vector_get(s->SSmodindvec, j);
    for (k = 0; k < r.nmodes / 4; k++) {
      qind = (int) gsl_vector_get(s->SAmodindvec, k);
      for (l = 0; l < r.nmodes / 4; l++) {
        rind = (int) gsl_vector_get(s->SAmodindvec, l);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
        cou1++;
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    pind = (int) gsl_vector_get(s->SSmodindvec, j);
    for (k = 0; k < r.nmodes / 4; k++) {
      qind = (int) gsl_vector_get(s->ASmodindvec, k);
      for (l = 0; l < r.nmodes / 4; l++) {
        rind = (int) gsl_vector_get(s->ASmodindvec, l);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
        cou1++;
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    pind = (int) gsl_vector_get(s->SSmodindvec, j);
    for (k = 0; k < r.nmodes / 4; k++) {
      qind = (int) gsl_vector_get(s->AAmodindvec, k);
      for (l = 0; l < r.nmodes / 4; l++) {
        rind = (int) gsl_vector_get(s->AAmodindvec, l);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
        cou1++;
      }
    }
  }
  // s		p	q	r
  // SS		SA	SS	SA
  // 		AS		AS
  // 		AA		AA
  for (j = 0; j < r.nmodes / 4; j++) {
    pind = (int) gsl_vector_get(s->SAmodindvec, j);
    for (k = 0; k < r.nmodes / 4; k++) {
      qind = (int) gsl_vector_get(s->SSmodindvec, k);
      for (l = 0; l < r.nmodes / 4; l++) {
        rind = (int) gsl_vector_get(s->SAmodindvec, l);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
        cou1++;
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    pind = (int) gsl_vector_get(s->ASmodindvec, j);
    for (k = 0; k < r.nmodes / 4; k++) {
      qind = (int) gsl_vector_get(s->SSmodindvec, k);
      for (l = 0; l < r.nmodes / 4; l++) {
        rind = (int) gsl_vector_get(s->ASmodindvec, l);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
        cou1++;
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    pind = (int) gsl_vector_get(s->AAmodindvec, j);
    for (k = 0; k < r.nmodes / 4; k++) {
      qind = (int) gsl_vector_get(s->SSmodindvec, k);
      for (l = 0; l < r.nmodes / 4; l++) {
        rind = (int) gsl_vector_get(s->AAmodindvec, l);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
        cou1++;
      }
    }
  }
  // s		p	q	r
  // SS		SA	SA	SS
  // 		AS	AS
  // 		AA	AA
  for (j = 0; j < r.nmodes / 4; j++) {
    pind = (int) gsl_vector_get(s->SAmodindvec, j);
    for (k = 0; k < r.nmodes / 4; k++) {
      qind = (int) gsl_vector_get(s->SAmodindvec, k);
      for (l = 0; l < r.nmodes / 4; l++) {
        rind = (int) gsl_vector_get(s->SSmodindvec, l);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
        cou1++;
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    pind = (int) gsl_vector_get(s->ASmodindvec, j);
    for (k = 0; k < r.nmodes / 4; k++) {
      qind = (int) gsl_vector_get(s->ASmodindvec, k);
      for (l = 0; l < r.nmodes / 4; l++) {
        rind = (int) gsl_vector_get(s->SSmodindvec, l);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
        cou1++;
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    pind = (int) gsl_vector_get(s->AAmodindvec, j);
    for (k = 0; k < r.nmodes / 4; k++) {
      qind = (int) gsl_vector_get(s->AAmodindvec, k);
      for (l = 0; l < r.nmodes / 4; l++) {
        rind = (int) gsl_vector_get(s->SSmodindvec, l);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
        gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
        cou1++;
      }
    }
  }

  // p=1
  pind = 1;
  // modes coming from 4 distinct group
  // p		s	q	r
  // SS		SA	AS	AA
  // 		SA	AA	AS
  // 		AA	SA	AS
  // 		AA	AS	SA
  // 		AS	AA	SA
  // 		AS	SA	AA
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        qind = (int) gsl_vector_get(s->ASmodindvec, k);
        for (l = 0; l < r.nmodes / 4; l++) {
          rind = (int) gsl_vector_get(s->AAmodindvec, l);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
          cou1++;
        }
      }
    }
  }

  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        qind = (int) gsl_vector_get(s->AAmodindvec, k);
        for (l = 0; l < r.nmodes / 4; l++) {
          rind = (int) gsl_vector_get(s->ASmodindvec, l);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
          cou1++;
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->AAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        qind = (int) gsl_vector_get(s->SAmodindvec, k);
        for (l = 0; l < r.nmodes / 4; l++) {
          rind = (int) gsl_vector_get(s->ASmodindvec, l);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
          cou1++;
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->AAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        qind = (int) gsl_vector_get(s->ASmodindvec, k);
        for (l = 0; l < r.nmodes / 4; l++) {
          rind = (int) gsl_vector_get(s->SAmodindvec, l);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
          cou1++;
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->ASmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        qind = (int) gsl_vector_get(s->AAmodindvec, k);
        for (l = 0; l < r.nmodes / 4; l++) {
          rind = (int) gsl_vector_get(s->SAmodindvec, l);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
          cou1++;
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->ASmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        qind = (int) gsl_vector_get(s->SAmodindvec, k);
        for (l = 0; l < r.nmodes / 4; l++) {
          rind = (int) gsl_vector_get(s->AAmodindvec, l);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
          cou1++;
        }
      }
    }
  }
  // mode pairs coming from 2 distinct group
  // p		s	q	r
  // SS		SS	SA	SA
  // 			AS	AS
  // 			AA	AA
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        qind = (int) gsl_vector_get(s->SAmodindvec, k);
        for (l = 0; l < r.nmodes / 4; l++) {
          rind = (int) gsl_vector_get(s->SAmodindvec, l);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
          cou1++;
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        qind = (int) gsl_vector_get(s->ASmodindvec, k);
        for (l = 0; l < r.nmodes / 4; l++) {
          rind = (int) gsl_vector_get(s->ASmodindvec, l);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
          cou1++;
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        qind = (int) gsl_vector_get(s->AAmodindvec, k);
        for (l = 0; l < r.nmodes / 4; l++) {
          rind = (int) gsl_vector_get(s->AAmodindvec, l);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
          cou1++;
        }
      }
    }
  }
  // p		s	q	r
  // SS		SA	SS	SA
  // 		AS		AS
  // 		AA		AA
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        qind = (int) gsl_vector_get(s->SSmodindvec, k);
        for (l = 0; l < r.nmodes / 4; l++) {
          rind = (int) gsl_vector_get(s->SAmodindvec, l);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
          cou1++;
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->ASmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        qind = (int) gsl_vector_get(s->SSmodindvec, k);
        for (l = 0; l < r.nmodes / 4; l++) {
          rind = (int) gsl_vector_get(s->ASmodindvec, l);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
          cou1++;
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->AAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        qind = (int) gsl_vector_get(s->SSmodindvec, k);
        for (l = 0; l < r.nmodes / 4; l++) {
          rind = (int) gsl_vector_get(s->AAmodindvec, l);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
          cou1++;
        }
      }
    }
  }
  // p		s	q	r
  // SS		SA	SA	SS
  // 		AS	AS
  // 		AA	AA
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        qind = (int) gsl_vector_get(s->SAmodindvec, k);
        for (l = 0; l < r.nmodes / 4; l++) {
          rind = (int) gsl_vector_get(s->SSmodindvec, l);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
          cou1++;
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->ASmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        qind = (int) gsl_vector_get(s->ASmodindvec, k);
        for (l = 0; l < r.nmodes / 4; l++) {
          rind = (int) gsl_vector_get(s->SSmodindvec, l);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
          cou1++;
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->AAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        qind = (int) gsl_vector_get(s->AAmodindvec, k);
        for (l = 0; l < r.nmodes / 4; l++) {
          rind = (int) gsl_vector_get(s->SSmodindvec, l);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
          gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
          cou1++;
        }
      }
    }
  }
  // q=1
  qind = 1;
  // modes coming from 4 distinct group
  // q		s	p	r
  // SS		SA	AS	AA
  // 		SA	AA	AS
  // 		AA	SA	AS
  // 		AA	AS	SA
  // 		AS	AA	SA
  // 		AS	SA	AA
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->ASmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            rind = (int) gsl_vector_get(s->AAmodindvec, l);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
            cou1++;
          }
        }
      }
    }
  }

  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->AAmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            rind = (int) gsl_vector_get(s->ASmodindvec, l);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
            cou1++;
          }
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->AAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->SAmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            rind = (int) gsl_vector_get(s->ASmodindvec, l);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
            cou1++;
          }
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->AAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->ASmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            rind = (int) gsl_vector_get(s->SAmodindvec, l);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
            cou1++;
          }
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->ASmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->AAmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            rind = (int) gsl_vector_get(s->SAmodindvec, l);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
            cou1++;
          }
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->ASmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->SAmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            rind = (int) gsl_vector_get(s->AAmodindvec, l);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
            cou1++;
          }
        }
      }
    }
  }
  // mode pairs coming from 2 distinct group
  // q		s	p	r
  // SS		SS	SA	SA
  // 			AS	AS
  // 			AA	AA
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->SAmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            rind = (int) gsl_vector_get(s->SAmodindvec, l);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
            cou1++;
          }
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->ASmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            rind = (int) gsl_vector_get(s->ASmodindvec, l);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
            cou1++;
          }
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->AAmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            rind = (int) gsl_vector_get(s->AAmodindvec, l);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
            cou1++;
          }
        }
      }
    }
  }
  // q		s	p	r
  // SS		SA	SS	SA
  // 		AS		AS
  // 		AA		AA
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->SSmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            rind = (int) gsl_vector_get(s->SAmodindvec, l);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
            cou1++;
          }
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->ASmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->SSmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            rind = (int) gsl_vector_get(s->ASmodindvec, l);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
            cou1++;
          }
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->AAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->SSmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            rind = (int) gsl_vector_get(s->AAmodindvec, l);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
            cou1++;
          }
        }
      }
    }
  }
  // q		s	p	r
  // SS		SA	SA	SS
  // 		AS	AS
  // 		AA	AA
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->SAmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            rind = (int) gsl_vector_get(s->SSmodindvec, l);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
            cou1++;
          }
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->ASmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->ASmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            rind = (int) gsl_vector_get(s->SSmodindvec, l);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
            cou1++;
          }
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->AAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->AAmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            rind = (int) gsl_vector_get(s->SSmodindvec, l);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
            gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
            cou1++;
          }
        }
      }
    }
  }
  // r=1
  rind = 1;
  // modes coming from 4 distinct group
  // r		s	p	q
  // SS		SA	AS	AA
  // 		SA	AA	AS
  // 		AA	SA	AS
  // 		AA	AS	SA
  // 		AS	AA	SA
  // 		AS	SA	AA
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->ASmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            qind = (int) gsl_vector_get(s->AAmodindvec, l);
            if (qind != 1) {
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
              cou1++;
            }
          }
        }
      }
    }
  }

  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->AAmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            qind = (int) gsl_vector_get(s->ASmodindvec, l);
            if (qind != 1) {
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
              cou1++;
            }
          }
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->AAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->SAmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            qind = (int) gsl_vector_get(s->ASmodindvec, l);
            if (qind != 1) {
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
              cou1++;
            }
          }
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->AAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->ASmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            qind = (int) gsl_vector_get(s->SAmodindvec, l);
            if (qind != 1) {
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
              cou1++;
            }
          }
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->ASmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->AAmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            qind = (int) gsl_vector_get(s->SAmodindvec, l);
            if (qind != 1) {
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
              cou1++;
            }
          }
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->ASmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->SAmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            qind = (int) gsl_vector_get(s->AAmodindvec, l);
            if (qind != 1) {
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
              cou1++;
            }
          }
        }
      }
    }
  }
  // mode pairs coming from 2 distinct group
  // r		s	p	q
  // SS		SS	SA	SA
  // 			AS	AS
  // 			AA	AA
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->SAmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            qind = (int) gsl_vector_get(s->SAmodindvec, l);
            if (qind != 1) {
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
              cou1++;
            }
          }
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->ASmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            qind = (int) gsl_vector_get(s->ASmodindvec, l);
            if (qind != 1) {
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
              cou1++;
            }
          }
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->SSmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->AAmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            qind = (int) gsl_vector_get(s->AAmodindvec, l);
            if (qind != 1) {
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
              cou1++;
            }
          }
        }
      }
    }
  }
  // r		s	p	q
  // SS		SA	SS	SA
  // 		AS		AS
  // 		AA		AA
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->SSmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            qind = (int) gsl_vector_get(s->SAmodindvec, l);
            if (qind != 1) {
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
              cou1++;
            }
          }
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->ASmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->SSmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            qind = (int) gsl_vector_get(s->ASmodindvec, l);
            if (qind != 1) {
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
              cou1++;
            }
          }
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->AAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->SSmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            qind = (int) gsl_vector_get(s->AAmodindvec, l);
            if (qind != 1) {
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
              cou1++;
            }
          }
        }
      }
    }
  }
  // r		s	p	q
  // SS		SA	SA	SS
  // 		AS	AS
  // 		AA	AA
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->SAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->SAmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            qind = (int) gsl_vector_get(s->SSmodindvec, l);
            if (qind != 1) {
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
              cou1++;
            }
          }
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->ASmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->ASmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            qind = (int) gsl_vector_get(s->SSmodindvec, l);
            if (qind != 1) {
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
              cou1++;
            }
          }
        }
      }
    }
  }
  for (j = 0; j < r.nmodes / 4; j++) {
    sind = (int) gsl_vector_get(s->AAmodindvec, j);
    if (sind != 1) {
      for (k = 0; k < r.nmodes / 4; k++) {
        pind = (int) gsl_vector_get(s->AAmodindvec, k);
        if (pind != 1) {
          for (l = 0; l < r.nmodes / 4; l++) {
            qind = (int) gsl_vector_get(s->SSmodindvec, l);
            if (qind != 1) {
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 0, sind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 1, pind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 2, qind);
              gsl_matrix_set(s->NZs_pqrcombmat, cou1, 3, rind);
              cou1++;
            }
          }
        }
      }
    }
  }
  s->NZcou = cou1;
}
