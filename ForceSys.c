/*
 * This function calculates the system force
 */

#include "struct_def.h"

double	ForceSys(sys_var s, run_param r, int j)
{
  int 	sind, pind, qind, rind;
  int 	cou2, lo, hi;
  double fr;
  double q_s, q_p, q_q, q_r;
  double f_t, f2_t, alpha;

  sind  = (int) gsl_matrix_get(s.IRs_pqrcountmat, j, 0);
  fr 	= gsl_vector_get(s.frvec, sind-1);
  q_s 	= gsl_vector_get(s.qvec, sind-1);

  f_t 	= -pow(2*PI*fr, 2.0)*q_s;
	
  lo 	= gsl_matrix_get(s.IRs_pqrcountmat, j, 2)- gsl_matrix_get(s.IRs_pqrcountmat, j, 1);
  hi 	= gsl_matrix_get(s.IRs_pqrcountmat, j, 2);
  f2_t 	= 0.0;
  
  for (cou2=lo; cou2<hi; cou2++)
  {
    pind  = (int) gsl_matrix_get(s.IRs_pqrcombmatsorted, cou2, 1);
    qind  = (int) gsl_matrix_get(s.IRs_pqrcombmatsorted, cou2, 2);
    rind  = (int) gsl_matrix_get(s.IRs_pqrcombmatsorted, cou2, 3);

    q_p   = gsl_vector_get(s.qvec, pind-1);
    q_q   = gsl_vector_get(s.qvec, qind-1);
    q_r   = gsl_vector_get(s.qvec, rind-1);
		
    alpha = gsl_matrix_get(s.alphamat, cou2, 0);
    f2_t  = f2_t - alpha*q_p*q_q*q_r;
  }

  return f_t+f2_t;
}
