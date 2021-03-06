/*
 * This function calculates the system force
 */

#include "struct_def.h"

for_var	ForceSys(sys_var s, run_param r, int j)
{
  int 	sind, pind, qind, rind;
  int 	cou2, lo, hi;
  double fr, frq2;
  double q_s, q_p, q_q, q_r;
  double alpha, mult, f3ssss, ep4ssss;
  for_var f;

  sind  = (int) gsl_matrix_get(s.IRs_pqrcountmat, j, 0);
  fr 	= gsl_vector_get(s.frvec, sind-1);
  frq2 	= gsl_vector_get(s.frq2vec, sind-1);
  q_s 	= gsl_vector_get(s.qvec, sind-1);

  f.f1 	= -pow(2*PI*fr, 2.0)*q_s;
	
  lo 	= gsl_matrix_get(s.IRs_pqrcountmat, j, 2)- gsl_matrix_get(s.IRs_pqrcountmat, j, 1);
  hi 	= gsl_matrix_get(s.IRs_pqrcountmat, j, 2);

  f.f3 	= -pow(2*PI*frq2, 2.0)*q_s*q_s*q_s;
  f.ep4 = 0.5*pow(2*PI*frq2, 2.0)*q_s*q_s*q_s*q_s;
  f3ssss = 0.0;
  ep4ssss = 0.0;
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
    // if ( sind==1 ){
    // printf("%d\t%d\t%d\t%d\n", sind, pind, qind, rind);

    q_p   = gsl_vector_get(s.qvec, pind-1);
    q_q   = gsl_vector_get(s.qvec, qind-1);
    q_r   = gsl_vector_get(s.qvec, rind-1);
		
    alpha = gsl_matrix_get(s.alphamatsorted, cou2, 0)*mult;

    if ( sind==pind && pind==qind && qind==rind ){
      f3ssss  = -alpha*q_p*q_q*q_r;
      ep4ssss = alpha*q_p*q_q*q_r*q_s/4.0;
    }
    else
    {
      f.f3  = f.f3 - alpha*q_p*q_q*q_r;
      f.ep4 = f.ep4 + alpha*q_p*q_q*q_r*q_s/4.0;
    }
  }
  f.f3 = f.f3*r.modefact + f3ssss;
  f.ep4 = f.ep4*r.modefact + ep4ssss;

  return f;
}
