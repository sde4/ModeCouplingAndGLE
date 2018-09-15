/*
 * This function calculates the system force
 */

#include "struct_def.h"

for_var	ForceSys(sys_var s, run_param r, int j)
{
  int 	sind, pind, qind, rind;
  int 	cou2, lo, hi;
  double fr;
  double q_s, q_p, q_q, q_r;
  double alpha, mult;
  for_var f;

  sind  = (int) gsl_matrix_get(s.IRs_pqrcountmat, j, 0);
  fr 	= gsl_vector_get(s.frvec, sind-1);
  q_s 	= gsl_vector_get(s.qvec, sind-1);

  f.f1 	= -pow(2*PI*fr, 2.0)*q_s;
	
  lo 	= gsl_matrix_get(s.IRs_pqrcountmat, j, 2)- gsl_matrix_get(s.IRs_pqrcountmat, j, 1);
  hi 	= gsl_matrix_get(s.IRs_pqrcountmat, j, 2);

  f.f3 	= 0.0;
  f.ep4 = 0.0;
  for (cou2=lo; cou2<hi; cou2++)
  {
    pind  = (int) gsl_matrix_get(s.IRs_pqrcombmatsorted, cou2, 1);
    qind  = (int) gsl_matrix_get(s.IRs_pqrcombmatsorted, cou2, 2);
    rind  = (int) gsl_matrix_get(s.IRs_pqrcombmatsorted, cou2, 3);
    mult  = 	  gsl_matrix_get(s.IRs_pqrcombmatsorted, cou2, 5);

    // if ( !((sind==pind && qind==rind) || (sind==qind && pind==rind) || (sind==rind && qind==pind)) ){
    // if ( !((sind==pind) || (sind==qind) || (sind==rind)) ){
    // if ( !((sind!=pind) && (sind!=qind) && (sind!=rind)) ){
    // if ( !((sind!=pind) && (sind!=qind) && (sind!=rind) && (pind!=qind) && (pind!=rind) && (qind!=rind)) ){
    // if ( ((sind==pind) && (qind==rind)) || ((sind==qind) && (pind==rind)) || ((sind==rind) && (pind==qind)) ){
    // if ( (sind==pind) || (qind==rind) || (sind==qind) || (pind==rind) || (sind==rind) || (pind==qind) ){
    // printf("%d\t%d\t%d\t%d\n", sind, pind, qind, rind);
    q_p   = gsl_vector_get(s.qvec, pind-1);
    q_q   = gsl_vector_get(s.qvec, qind-1);
    q_r   = gsl_vector_get(s.qvec, rind-1);
		
    alpha = gsl_matrix_get(s.alphamat, cou2, 0)*mult;
    f.f3  = f.f3 - alpha*q_p*q_q*q_r;
    f.ep4 = f.ep4 + alpha*q_p*q_q*q_r*q_s/4.0;
    // }
  }

  return f;
}
