/*
 * This function calculates the system force
 */

#include "struct_def.h"

double	ForceSys(sys_var s, run_param r, int j)
{
	int 	k;
	double 	fr;
	double 	q_j, q_k;
	double 	f_t, f2_t, alpha;

	fr 		= gsl_vector_get(s.frvec, j);
	q_j 		= gsl_vector_get(s.qvec, j);

	f_t 		= -pow(2*PI*fr, 2.0)*q_j;
	f2_t 		= 0.0;
	for (k=0; k<r.nmodes; k++)
	{
		q_k 		= gsl_vector_get(s.qvec, k);
		alpha 		= gsl_matrix_get(s.alphamat, j, k);
		f2_t 		= f2_t - alpha*q_k*q_k*q_j;
	}

	return f_t+f2_t;
}
