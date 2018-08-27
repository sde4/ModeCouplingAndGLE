/*
 * This function finds the mode combinations which lead to non-zero coupling constant
 */

#include <omp.h>
#include "struct_def.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

/* Definitions for routines */
double GalerkinProj(gsl_vector* , gsl_vector*, gsl_vector* , disc_const );
void ModeshapesDisp(gsl_vector* , int , int , disc_const );

void ThirdOrderCouplingCalculation(sys_var* s, mat_const mc, disc_const dc) {
  printf("#    - Third order mode coupling calculations ...\n");

  /**************************/
  // Third order mode
  // coupling!
  /**************************/

  int i, j, k, l, n;
  int cou1, cou2, cou3, cou4;
  int sind, pind, qind, rind;
  FILE *infp, *outfp;

  double dat, gams_pqr, Hn_pq, Hn_rs, xeta_n4, norm2D;
  int mords, nords, mordp, nordp, mordq, nordq, mordr, nordr;

  gsl_matrix * Fmat 	  = gsl_matrix_alloc (dc.Nfx*dc.Nfy, dc.Nfx*dc.Nfy);
  gsl_vector * Feval 	  = gsl_vector_alloc (dc.Nfx*dc.Nfy);
  gsl_matrix * Fevec 	  = gsl_matrix_alloc (dc.Nfx*dc.Nfy, dc.Nfx*dc.Nfy);
  s->alphamat  	  	  = gsl_matrix_alloc(s->IRcou, 2);
  gsl_vector * phi_s 	  = gsl_vector_alloc(dc.Nfx*dc.Nfy);
  gsl_vector * phi_p 	  = gsl_vector_alloc(dc.Nfx*dc.Nfy);
  gsl_vector * phi_q 	  = gsl_vector_alloc(dc.Nfx*dc.Nfy);
  gsl_vector * phi_r 	  = gsl_vector_alloc(dc.Nfx*dc.Nfy);
  gsl_vector * psi_n 	  = gsl_vector_alloc(dc.Nfx*dc.Nfy);

  // First compute the eigen vectors and eigen values of the 
  // FD discretized BH operator.
  // Discretization is done using Mathematica
  // Fmat import !!

  printf("#         o Input file required: 1. in.Fmat\n");
  // check for the required input files

  if ((infp=fopen("in.Fmat", "r")) == NULL) 
  {
     printf("#           ERROR: File not found! \n");
     exit(1);
  }
  infp = fopen("in.Fmat", "r");
  for (i = 0; i < dc.Nfx*dc.Nfy; i++) {
    for (j = 0; j < dc.Nfx*dc.Nfy; j++) {
      fscanf(infp, "%lf", &dat);
      gsl_matrix_set(Fmat, i, j, dat);
    }
  }
  fclose(infp);
  // outfp = fopen("Fmat.check", "w");
  // for (i = 0; i < dc.Nfx*dc.Nfy; i++) {
  //   for (j = 0; j < dc.Nfx*dc.Nfy; j++) {
  //     fprintf(infp, "%5.5e\t", gsl_matrix_get(Fmat, i, j));
  //   }
  //   fprintf(infp, "\n");
  // }
  // fclose(outfp);

  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (dc.Nfx*dc.Nfy);
  gsl_eigen_symmv (Fmat, Feval, Fevec, w);
  gsl_eigen_symmv_free (w);
  gsl_eigen_symmv_sort (Feval, Fevec, GSL_EIGEN_SORT_ABS_ASC);

  printf("#          Mode groups and couling constants ...\n");
  for (i = 0; i < s->IRcou; i++){
    sind = gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 0);
    pind = gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 1);
    qind = gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 2);
    rind = gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 3);
    
    mords = gsl_matrix_get(s->modindmat, sind-1, 0); 
    nords = gsl_matrix_get(s->modindmat, sind-1, 1);
    ModeshapesDisp(phi_s, mords, nords, dc);
    norm2D  = gsl_blas_dnrm2(phi_s)*pow(dc.hx*dc.hy, 0.5);
    gsl_vector_scale (phi_s, 1.0/norm2D);
 
    mordp = gsl_matrix_get(s->modindmat, pind-1, 0); 
    nordp = gsl_matrix_get(s->modindmat, pind-1, 1);
    ModeshapesDisp(phi_p, mordp, nordp, dc);
    norm2D  = gsl_blas_dnrm2(phi_p)*pow(dc.hx*dc.hy, 0.5);
    gsl_vector_scale (phi_p, 1.0/norm2D);

    mordq = gsl_matrix_get(s->modindmat, qind-1, 0); 
    nordq = gsl_matrix_get(s->modindmat, qind-1, 1);
    ModeshapesDisp(phi_q, mordq, nordq, dc);
    norm2D  = gsl_blas_dnrm2(phi_q)*pow(dc.hx*dc.hy, 0.5);
    gsl_vector_scale (phi_q, 1.0/norm2D);

    mordr = gsl_matrix_get(s->modindmat, rind-1, 0); 
    nordr = gsl_matrix_get(s->modindmat, rind-1, 1);
    ModeshapesDisp(phi_r, mordr, nordr, dc);
    norm2D  = gsl_blas_dnrm2(phi_r)*pow(dc.hx*dc.hy, 0.5);
    gsl_vector_scale (phi_r, 1.0/norm2D);

    gams_pqr = 0.0;
    for (n= 3; n < dc.Nfx*dc.Nfy; n++){				//excluding first 3 with zero eigen values
      gsl_matrix_get_col (psi_n, Fevec, n);
      norm2D  = gsl_blas_dnrm2(psi_n)*pow(dc.hx*dc.hy, 0.5);
      gsl_vector_scale (psi_n, 1.0/norm2D);
      xeta_n4 = gsl_vector_get(Feval, n);
       
      Hn_pq   = GalerkinProj(psi_n, phi_p, phi_q, dc);
      Hn_rs   = GalerkinProj(psi_n, phi_r, phi_s, dc);
      
      gams_pqr += Hn_pq*Hn_rs/(2.0*xeta_n4); 			// the unit of gamma here is um^-4 and needs to be converted to A^-4
    }
    gsl_matrix_set(s->alphamat, i, 0, gams_pqr*mc.Et/mc.rho*dc.gamunitconv);
    gsl_matrix_set(s->alphamat, i, 1, gams_pqr*pow(dc.Lx*dc.Ly, 3.0));
    printf("#          %03d(%02dx%02d)\t%03d(%02dx%02d)\t%03d(%02dx%02d)\t%03d(%02dx%02d)\t%f\t%5.5e\n", 
		  sind, mords, nords, pind, mordp, nordp, qind, mordq, nordq, rind, mordr, nordr,
		  gams_pqr*pow(dc.Lx*dc.Ly, 3.0), gams_pqr*mc.Et/mc.rho*dc.gamunitconv);
  }
  
  return;
}
